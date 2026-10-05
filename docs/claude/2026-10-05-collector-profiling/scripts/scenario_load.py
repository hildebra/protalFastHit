#!/usr/bin/env python3
"""The work of the default scenarios' drawn reads (Ultima, PacBio, Nanopore) in a GTDB r226 build, as the collector
does it now and with one pass over a sample's genomes: genome reads (decompress + parse), templates, reads, and how
much of that holds the GIL, per collection (training: 3 samples per scenario, test: 2).

The drawing is modelled after long_read_templates (vectorised): each round a sample's (or chunk's) missing bases as
reads of genomes drawn by abundance x length, the setup's gamma lengths, starts uniform over the contigs, cut where
a contig ends; a genome is read once per round it gets a read in. Communities: the scenarios' species per sample
scaled to the r226 v11 pool (scenarios.fit_species: 5,377 species the training database has, 2,621 it lacks),
x1.4 genomes per species (strains 0.3,0.1), lognormal abundances (sigma 1.5 and 2.0 alternating) or the host
scenario's power law, genome lengths and contigs like r226's (median 3.4 Mb, ~78 contigs). Unit costs (laptop core,
unit_costs.py) are given as arguments.

usage: scenario_load.py --scripts PROTAL/scripts [--known 5377 --lacking 2621]
"""

import argparse
import json
import math
import sys

import numpy as np

# Unit costs measured by unit_costs.py (Python 3.12, Core Ultra 7 258V, one core): seconds.
COST = {
    "genome_zlib": 0.0138, "genome_parse": 0.0072,  # parse incl. file read: holds the GIL; zlib releases it
    "template_read": {"ultima": 3.7e-6, "hifi": 26.9e-6, "qshmm": 16.6e-6},  # per read, GIL held
    "host_read": 2.5e-6,  # a host template (Host.draw), per read, GIL held
    "reads_per_mb": {"ultima": 0.223, "hifi": 0.112, "qshmm": 0.591},  # hifi_reads.py (numpy) or pbsim3 (a child)
    "rename_per_mb": 0.031,  # pbsim3's reads renamed and gzipped by Python, GIL held
    "host_fragment": 2.46e-6,  # scenarios.host_pe_chunk's loop, per host read pair, GIL held
    "art_pair": 54e-6,  # ART per 150 bp pair (2026-10-03 report), a child
}
GENOMES_PER_SPECIES = 1.4


def genome_structure(rng):
    length = int(np.clip(rng.lognormal(np.log(3.4e6), 0.35), 1.2e6, 9e6))
    contigs = int(np.clip(rng.lognormal(np.log(78), 0.8), 1, 600))
    cuts = np.sort(rng.choice(np.arange(1, length), contigs - 1, replace=False)) if contigs > 1 else np.array([], int)
    sizes = np.diff(np.concatenate(([0], cuts, [length])))
    return length, sizes[sizes >= 100]


def community(rng, genomes, abundance):
    structures = [genome_structure(rng) for _ in range(genomes)]
    lengths = np.array([s[0] for s in structures], dtype=float)
    if abundance.startswith("powerlaw"):
        alpha = float(abundance.split(":")[1])
        ab = 1.0 / np.arange(1, genomes + 1) ** alpha
    else:
        ab = rng.lognormal(0.0, float(abundance), genomes)
    return ab * lengths, [s[1] for s in structures]


def draw_rounds(rng, weights, contigs, bases, mean, sd):
    """[genomes read in each round] for one sample or chunk of `bases`."""
    cum = np.cumsum(weights)
    total = cum[-1]
    shape, scale = (mean / sd) ** 2, sd * sd / mean
    starts_of = [np.cumsum(c - 100 + 1) for c in contigs]
    got, rounds = 0, []
    while got < bases:
        need = bases - got
        n = int(need / mean * 1.05) + 16
        lengths = np.clip(np.rint(rng.gamma(shape, scale, n)), 100, 1_000_000).astype(np.int64)
        keep = np.searchsorted(np.cumsum(lengths), need) + 1  # until the planned bases reach the need
        lengths = lengths[:keep]
        g = np.searchsorted(cum, rng.random(len(lengths)) * total, side="right")
        g = np.minimum(g, len(weights) - 1)
        order = np.argsort(g, kind="stable")
        gs, ls = g[order], lengths[order]
        bounds = np.flatnonzero(np.diff(gs)) + 1
        made = 0
        for idx, part in zip(np.split(gs, bounds), np.split(ls, bounds)):
            k = int(idx[0])
            starts = starts_of[k]
            at = rng.integers(0, starts[-1], len(part))
            c = np.searchsorted(starts, at, side="right")
            start = at - np.where(c > 0, starts[np.maximum(c - 1, 0)], 0)
            made += int(np.minimum(part, contigs[k][c] - start).sum())
        rounds.append(np.unique(gs))
        got += made
    return rounds


def unit_load(rng, unit, genomes, abundances, chunk, host_share):
    """Genome reads (now: per chunk and round; one pass: per sample and round), reads and bases of a drawn unit's
    samples."""
    setup = unit["setup"]
    bases = unit["bases"]
    k = max(1, -(-bases // chunk)) if chunk and bases > chunk else 1
    out = {"samples": unit["samples"], "chunks": k * unit["samples"], "parses_now": 0, "parses_onepass": 0, "reads": 0,
           "host_reads": 0, "bases": bases * unit["samples"]}
    for s in range(unit["samples"]):
        weights, contigs = community(rng, genomes, abundances[s % len(abundances)])
        community_bases = bases * (1 - host_share)
        per_chunk = [draw_rounds(rng, weights, contigs, community_bases / k, setup["length_mean"], setup["length_sd"])
                     for _ in range(k)]
        out["parses_now"] += sum(len(r) for rounds in per_chunk for r in rounds)
        depth = max(len(rounds) for rounds in per_chunk)
        for r in range(depth):
            out["parses_onepass"] += len(set().union(*(set(rounds[r].tolist()) for rounds in per_chunk if r < len(rounds))))
        out["reads"] += community_bases / setup["length_mean"]
        out["host_reads"] += bases * host_share / setup["length_mean"]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--known", type=int, default=5377)
    ap.add_argument("--lacking", type=int, default=2621)
    ap.add_argument("--chunk", type=int, default=250_000_000)
    ap.add_argument("--json", help="write the table here too")
    a = ap.parse_args()
    sys.path.insert(0, a.scripts)
    import collect_training_data as ctd
    import scenarios

    rng = np.random.default_rng(5)
    defs = scenarios.definitions()
    rows = []
    for collection, samples in (("training", 3), ("test", 2)):
        class Opts:
            pass
        opts = Opts()
        opts.scenarios = ",".join(f"{n}:{samples}" for n in defs)
        opts.scenario_samples, opts.scenario_file = samples, None
        opts.read_types = ["pe", "se", "pb", "ont"]
        opts.pb_setup, opts.ont_setup = "hifi:15000:3000:3", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36"
        points, units = ctd.scenario_units(opts)
        for point in points:
            d, _ = scenarios.fit_species(point["definition"], a.known, a.lacking)
            lo, hi = scenarios.species_bounds(d["species"])
            genomes = round((lo + hi) / 2 * GENOMES_PER_SPECIES)
            abundances = d["abundance"].split(":")[1].split(",") if d["abundance"].startswith("lognormal") else [d["abundance"]]
            host = d["host_share"]
            for unit in units:
                if unit.get("scenario") != point["scenario"]:
                    continue
                if unit["type"] == "pe":
                    pairs = int(point["read_pairs"])
                    rows.append({"collection": collection, "scenario": point["scenario"], "unit": unit["name"],
                                 "type": "pe", "samples": unit["samples"], "genomes": genomes,
                                 "host_pairs": point["host_pairs"] * unit["samples"], "pairs": pairs * unit["samples"]})
                    continue
                load = unit_load(rng, unit, genomes, abundances, a.chunk, unit.get("host_share", 0))
                rows.append({"collection": collection, "scenario": point["scenario"], "unit": unit["name"],
                             "type": unit["type"], "method": unit["setup"]["method"], "genomes": genomes, **load})

    def costs(row, parses):
        method = row["method"]
        mb = row["bases"] / 1e6
        gil = parses * COST["genome_parse"] + row["reads"] * COST["template_read"][method] + \
            row["host_reads"] * COST["host_read"] + (mb * COST["rename_per_mb"] if method == "qshmm" else 0)
        free = parses * COST["genome_zlib"]
        reads = mb * COST["reads_per_mb"][method]  # numpy (hifi, ultima) or pbsim3's own process
        return gil, free, reads

    print(f"{'collection':10} {'unit':32} {'samples':>7} {'chunks':>6} {'genomes':>7} {'genome reads now':>16} "
          f"{'one pass':>9} {'GIL-held CPU h now':>18} {'one pass':>9} {'all CPU h now':>13} {'one pass':>9}")
    totals = {}
    for row in rows:
        if row["type"] == "pe":
            continue
        g_now, f_now, r = costs(row, row["parses_now"])
        g_one, f_one, _ = costs(row, row["parses_onepass"])
        row.update(gil_now_h=g_now / 3600, all_now_h=(g_now + f_now + r) / 3600, gil_onepass_h=g_one / 3600,
                   all_onepass_h=(g_one + f_one + r) / 3600, reads_h=r / 3600)
        t = totals.setdefault(row["collection"], {"parses_now": 0, "parses_onepass": 0, "gil_now_h": 0, "gil_onepass_h": 0,
                                                  "all_now_h": 0, "all_onepass_h": 0, "reads_h": 0, "chunks": 0,
                                                  "bases": 0})
        for key in t:
            t[key] += row[key]
        print(f"{row['collection']:10} {row['unit']:32} {row['samples']:7d} {row['chunks']:6d} {row['genomes']:7d} "
              f"{row['parses_now']:16,d} {row['parses_onepass']:9,d} {g_now / 3600:18.2f} {g_one / 3600:9.2f} "
              f"{(g_now + f_now + r) / 3600:13.2f} {(g_one + f_one + r) / 3600:9.2f}")
    for row in rows:
        if row["type"] == "pe":
            gil = row["host_pairs"] * COST["host_fragment"]
            row.update(host_gil_h=gil / 3600, host_art_h=row["host_pairs"] * COST["art_pair"] / 3600)
            print(f"{row['collection']:10} {row['unit']:32} {row['samples']:7d} {'':6} {row['genomes']:7d} "
                  f"paired-end, {row['pairs']:,d} pairs" + (f"; host pairs {row['host_pairs']:,d}: GIL-held "
                                                            f"{gil / 3600:.2f} h, ART {row['host_pairs'] * COST['art_pair'] / 3600:.2f} h"
                                                            if row["host_pairs"] else ""))
            t = totals[row["collection"]]
            t["host_gil_h"] = t.get("host_gil_h", 0) + gil / 3600
    for collection, t in totals.items():
        print(f"{collection}: {t['chunks']} chunks, {t['bases'] / 1e9:.1f} Gb drawn; genome reads {t['parses_now']:,} now, "
              f"{t['parses_onepass']:,} in one pass; GIL-held CPU {t['gil_now_h']:.2f} h now (+{t.get('host_gil_h', 0):.2f} h "
              f"host fragments), {t['gil_onepass_h']:.2f} h in one pass; all CPU {t['all_now_h']:.2f} h now, "
              f"{t['all_onepass_h']:.2f} h in one pass (of it reads {t['reads_h']:.2f} h)")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({"rows": rows, "totals": totals, "cost": COST}, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
