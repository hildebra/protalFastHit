#!/usr/bin/env python3
"""Write a simulate_metagenomes --from_manifest manifest for one phylo-audit scenario.

Samples are named after the tip they carry (t01..t12): each sample holds tip tNN of every focal
species (Mockella alpha, Calidella fervens, Testella one) plus background species from the world.

Scenarios:
  depth      all focal strains at --depth x (baseline = 20)
  mixed      as depth, plus sample 'mix' holding two strains per focal species at --mix_major /
             (1-mix_major) of --depth, major = --mix_a, minor = --mix_b
  congener   as depth, but in --congener_samples the Mockella congeners (beta, gamma) are at
             --congener_depth x instead of the background depth
"""
import argparse
import csv
import os
import random

PHYLO = os.path.expanduser("~/audit5/phylo")
WORLD = os.path.expanduser("~/audit5/world/gtdb_r226/simulation/genomes.tsv")
FOCAL = {"Malpha": "Mockella alpha", "Cferv": "Calidella fervens", "Tone": "Testella one"}
BACKGROUND = {"Mockella beta": 10.0, "Mockella gamma": 10.0, "Fakibacter gamma": 10.0}
READ_LEN = 150


def load_table(path, has_header=True):
    rows = []
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        if has_header:
            next(rd)
        for r in rd:
            if r:
                rows.append({"accession": r[0], "taxonomy": r[1], "fasta": r[2], "length": int(r[3])})
    return rows


def species_of(tax):
    return tax.split("s__")[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=["depth", "mixed", "congener"], default="depth")
    ap.add_argument("--depth", type=float, default=20.0)
    ap.add_argument("--bg_depth", type=float, default=None, help="override background depth")
    ap.add_argument("--mix_a", default="t03")
    ap.add_argument("--mix_b", default="t10")
    ap.add_argument("--mix_major", type=float, default=0.7)
    ap.add_argument("--congener_samples", default="t02,t05,t08,t11")
    ap.add_argument("--congener_depth", type=float, default=100.0)
    ap.add_argument("--ntips", type=int, default=12)
    ap.add_argument("--tip_depths", default=None,
                    help="comma list of focal-strain depths for t01..tNN (uneven-depth scenario)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rnd = random.Random(args.seed)

    world = load_table(WORLD)
    by_species = {}
    for g in world:
        by_species.setdefault(species_of(g["taxonomy"]), []).append(g)
    strains = {}
    for short in FOCAL:
        for g in load_table(os.path.join(PHYLO, "strains", short, "genomes.tsv"), has_header=False):
            tip = g["accession"].rsplit("_", 1)[-1]
            strains[(short, tip)] = g

    tips = [f"t{i+1:02d}" for i in range(args.ntips)]
    samples = []  # (sample, [(genome, depth)])
    cong = set(args.congener_samples.split(",")) if args.scenario == "congener" else set()
    tip_depth = {t: args.depth for t in tips}
    if args.tip_depths:
        for t, d in zip(tips, args.tip_depths.split(",")):
            tip_depth[t] = float(d)
    for tip in tips:
        members = [(strains[(s, tip)], tip_depth[tip]) for s in FOCAL]
        for sp, d in BACKGROUND.items():
            g = rnd.choice(by_species[sp])
            dd = args.bg_depth if args.bg_depth is not None else d
            if tip in cong and sp.startswith("Mockella"):
                dd = args.congener_depth
            members.append((g, dd))
        samples.append((tip, members))
    if args.scenario == "mixed":
        members = []
        for s in FOCAL:
            members.append((strains[(s, args.mix_a)], args.depth * args.mix_major))
            members.append((strains[(s, args.mix_b)], args.depth * (1 - args.mix_major)))
        for sp, d in BACKGROUND.items():
            members.append((rnd.choice(by_species[sp]), args.bg_depth if args.bg_depth is not None else d))
        samples.append(("mix", members))

    with open(args.out, "w") as fh:
        fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\t"
                 "relative_abundance\tfasta_path\n")
        for name, members in samples:
            rps = [max(1, round(d * g["length"] / (2 * READ_LEN))) for g, d in members]
            tot = sum(rps)
            for (g, d), rp in zip(members, rps):
                vc = rp * 2 * READ_LEN / g["length"]
                fh.write(f"{name}\t{g['accession']}\t{species_of(g['taxonomy'])}\t{g['taxonomy']}\t"
                         f"{g['length']}\t{rp}\t{vc:.4f}\t{rp/tot:.6f}\t{g['fasta']}\n")


if __name__ == "__main__":
    main()
