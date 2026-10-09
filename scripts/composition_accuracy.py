#!/usr/bin/env python3
"""How accurate the profiles' composition is on simulated samples: the share of the reads the called species explain,
the unknown share that ends a profile ("?"), the average genome size, the species' genome sizes and depths, and the
species the rest would be (protal's Composition.h, docs/running.md#what-the-called-species-explain-the-unknown-share).

For each sample of a collection (collect_training_data.py -o; build_gtdb_database.py: OUT/work/training and
OUT/work/test), the composition is computed again as protal computes it, from the depth and genome size of every taxon
in the sample's <profile>.log and the reads in its <profile>.composition, over the species a model calls: by default those of
--calls (the trainer's PREFIX.calls.tsv.gz: the final model's calls on the test set, the calls with species held out
on the training samples), else those the profile called. The truth is the simulator's manifest of the sample's
community (a long-read or Ultima sample's: the paired-end community it replays, drawn by bases):
  explained          the share of the reads from the species the database has (the called ones: explained_by_calls)
  unknown            the share of the cells of the species the database lacks (not called: unknown_given_calls), as
                     protal's "?" counts genomes; the host's reads, which protal counts as unknown genomes of the
                     called species' average size, are left out: samples with a host are summarised apart
  average genome     of all cells, and of the cells of the species called and present
  missing species    the species present and not called, against protal's MissingSpeciesAt{Median,Lowest}Depth
  per species        a present species' GenomeSize over the length of its genomes simulated, and its depth over the
                     depth its reads give (paired-end: the manifest's vertical coverage times the sample's
                     FragmentBaseShare; single-end: half the pairs' coverage; drawn reads: their bases by the genomes'
                     weights)
A paired-end or single-end sample's host share is its reads beyond the community's; a drawn sample's that of the
paired-end sample whose community it replays (drawn by bases at the same share).

Writes --out (one line per sample) and prints a summary per set and read type (--summary writes it too): the median
and the 10-90% range of each estimate's error against its truth. --species-out writes every sample's species present
or called, protal's depth and genome size beside the truth, so that a bias can be traced to the species it comes from
off the cluster, where the profiles and manifests are not.

  python3 scripts/composition_accuracy.py --calls OUT/model_logs/trained_model.calls.tsv.gz --training OUT/work/training \\
      --test OUT/work/test --heldout OUT/model_logs/heldout_species.txt --read-type pe \\
      --out OUT/model_logs/composition_accuracy.tsv --species-out OUT/model_logs/composition_species.tsv.gz
"""
import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collect_training_data as collect  # noqa: E402

SETS = ("training", "test")
NAN = float("nan")
HOST_FREE = 0.01  # a sample with less of its reads from the host counts as without one
COLUMNS = ["set", "read_type", "point", "scenario", "sample", "host_share", "species", "absent_species", "called",
           "true_called", "missing_true", "scanned_fragments", "true_explained", "true_explained_by_calls", "explained",
           "true_unknown", "true_unknown_given_calls", "unknown", "true_ags", "true_ags_called", "ags", "ags_p10", "ags_p90",
           "missing_at_median", "missing_at_lowest", "depth_ratio", "size_ratio"]
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)  # Composition.h's kQuantiles
# --species-out: per sample, every species present or called. genomes: those its reads came from (the manifest's; an
# in-silico strain's name starts with insilico_); true_read_pairs: the paired-end community's (a drawn sample's reads
# are drawn from it by bases); depth and genome_size: protal's (VCov, GenomeSize; NA where the profile lacks the taxon).
SPECIES_COLUMNS = ["set", "read_type", "point", "scenario", "sample", "species", "taxid", "genomes", "present",
                   "in_database", "called", "true_read_pairs", "true_cell_share", "true_genome_length", "true_depth",
                   "depth", "genome_size"]


def number(text):
    try:
        return float(text) if text not in (None, "", "NA", "nan") else NAN
    except ValueError:
        return NAN


def read_profile_log(path):
    """{taxid: (name, predicted, depth, genome size)} of a <profile>.log (only the columns needed: the gene columns of a
    GTDB-sized database make its lines long)."""
    taxa = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "GenomeSize" not in header:
            return None  # a profile of protal before 2026-10-08
        at = {name: header.index(name) for name in ("Predicted", "TaxID", "Name", "VCov", "GenomeSize")}
        last = max(at.values())
        for line in fh:
            fields = line.split("\t", last + 1)
            if len(fields) <= last:
                continue
            taxa[fields[at["TaxID"]]] = (fields[at["Name"]], fields[at["Predicted"]] == "1", number(fields[at["VCov"]]),
                                         number(fields[at["GenomeSize"]]))
    return taxa


def read_composition(path):
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else None


def weighted_quantile(pairs, total, q):
    """Composition.h's WeightedQuantile: the smallest size whose cumulative weight reaches q of the total."""
    cumulative = 0.0
    for size, weight in pairs:
        cumulative += weight
        if cumulative >= q * total * (1 - 1e-12):
            return size
    return pairs[-1][0] if pairs else NAN


def compose(species, scanned_fragments, scanned_bases, share, has_sizes):
    """The composition of called species [(depth, genome size)] as Composition.h's Compute computes it."""
    depths = sorted(d for d, _ in species if d > 0)
    sized = sorted((s, d) for d, s in species if d > 0 and s > 0)
    ge = sum(depths)
    sized_depth = sum(d for _, d in sized)
    sized_bases = sum(s * d for s, d in sized)
    r = {"ags": NAN, "quantiles": [NAN] * 5, "explained": NAN, "unknown": NAN, "missing_at_median": NAN,
         "missing_at_lowest": NAN}
    if sized_depth > 0:
        r["ags"] = sized_bases / sized_depth
        r["quantiles"] = [weighted_quantile(sized, sized_depth, q) for q in QUANTILES]
        attributed = sized_bases + (ge - sized_depth) * r["ags"]
    else:
        attributed = 0.0 if not depths else NAN
    share = share if 0 < share <= 1 else 1.0
    if not scanned_bases > 0 or math.isnan(attributed):
        return r
    total = scanned_bases * share
    r["explained"] = attributed / total
    if not has_sizes:
        return r
    if not depths:
        r["unknown"] = 1.0
        return r
    unknown = max(0.0, total - attributed) / r["ags"]
    r["unknown"] = unknown / (ge + unknown)
    median = depths[len(depths) // 2] if len(depths) % 2 else (depths[len(depths) // 2 - 1] + depths[len(depths) // 2]) / 2
    r["missing_at_median"], r["missing_at_lowest"] = unknown / median, unknown / depths[0]
    return r


def read_heldout(path):
    """The species of heldout_species.txt (the first column), or an empty set."""
    if not path:
        return set()
    with open(path) as fh:
        return {line.split("\t")[0].strip() for line in fh if line.strip() and not line.startswith("#")}


def sample_files(collection, point, sample):
    """(the sample's <profile>.log, the manifest rows of its community, whether its reads are drawn, the community's
    paired-end <profile>.composition for a drawn sample) as the collector lays them out; None if missing."""
    points = os.path.join(collection, "points")
    drawn = os.path.join(points, point, "sim", "samples.tsv")
    if os.path.isfile(drawn):  # long reads and Ultima reads: drawn from a paired-end community sample's genomes
        with open(drawn) as fh:
            next(fh)
            community = next((f[3] for f in (line.rstrip("\n").split("\t") for line in fh) if f[0] == sample), None)
        log = os.path.join(points, point, "protal", "profiles", sample + ".profile.log")
        if community is None or not os.path.isfile(log):
            return None
        base = community.rsplit("_s_", 1)[0]
        rows = [r for r in collect.manifest_rows(os.path.join(points, base, "sim")) if r["sample"] == community]
        return log, rows, True, os.path.join(points, base, "protal", "profiles", community + ".profile.composition")
    if point.endswith("_se") and sample.endswith("_se"):  # a paired-end point's first reads, as single-end reads
        base = point[: -len("_se")]
        log = os.path.join(points, base, "protal_se", "profiles", sample + ".profile.log")
        sim = os.path.join(points, base, "sim")
        if not os.path.isfile(log) or not os.path.isfile(os.path.join(sim, "manifest.tsv")):
            return None
        return log, [r for r in collect.manifest_rows(sim) if r["sample"] == sample[: -len("_se")]], False, None
    sim = os.path.join(points, point, "sim")
    meta = os.path.join(sim, "protal.meta")
    if not os.path.isfile(meta):
        return None
    _, rows, _ = collect.map_rows(meta)
    profile = next((r["PROFILE"] for r in rows if r["SAMPLEID"] == sample), None)
    if profile is None or not os.path.isfile(profile + ".log"):
        return None
    return profile + ".log", [r for r in collect.manifest_rows(sim) if r["sample"] == sample], False, None


def samples_from_calls(path, read_type):
    """{(set, point, sample): (scenario, {called taxids})} of a read type's rows of a PREFIX.calls.tsv.gz."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("meta_read_type", read_type) != read_type or row.get("set") not in SETS:
                continue
            key = (row["set"], row["meta_design"], row["meta_sample"])
            entry = out.setdefault(key, (row.get("meta_scenario") or "", set()))
            if row.get("call") == "1":
                entry[1].add(str(row["taxon"]))
    return out


def samples_of_collection(collection, set_name, read_type):
    """{(set, point, sample): ("", None)} of a collection's paired-end (or, for se, single-end) samples: their profiles'
    own calls. Drawn samples need --calls (their read type is in its rows)."""
    out = {}
    for meta in sorted(glob.glob(os.path.join(collection, "points", "*", "sim", "protal.meta"))):
        point = os.path.basename(os.path.dirname(os.path.dirname(meta)))
        _, rows, _ = collect.map_rows(meta)
        for row in rows:
            if read_type == "pe":
                out[(set_name, point, row["SAMPLEID"])] = ("", None)
            elif read_type == "se":
                out[(set_name, point + "_se", row["SAMPLEID"] + "_se")] = ("", None)
    return out


def evaluate(job):
    """One sample's line (COLUMNS) and its present species' (depth ratio, size ratio), or (None, why)."""
    (set_name, point, sample), scenario, called_ids, collection, heldout, read_type = job
    files = sample_files(collection, point, sample)
    if files is None:
        return None, "files missing"
    log, rows, drawn, community_composition = files
    taxa = read_profile_log(log)
    composition = read_composition(log[: -len(".log")] + ".composition")
    if taxa is None or composition is None or not rows:
        return None, "no composition (a profile of protal before 2026-10-08)"
    scanned_fragments, scanned_bases = number(composition["ScannedFragments"]), number(composition["ScannedBases"])
    share = number(composition["FragmentBaseShare"])
    if math.isnan(scanned_bases):
        return None, "its SAM does not say how many reads were scanned"
    if called_ids is None:
        called_ids = {t for t, v in taxa.items() if v[1]}
    has_sizes = any(v[3] > 0 for v in taxa.values()) or composition.get("UnknownShare", "NA") != "NA"
    called = {t: taxa[t] for t in called_ids if t in taxa}
    est = compose([(v[2], v[3]) for v in called.values()], scanned_fragments, scanned_bases, share, has_sizes)

    # The truth: per species, its pairs, cells (relative abundance), cells times length, coverage.
    truth = collections.defaultdict(lambda: [0, 0.0, 0.0, 0.0])
    for r in rows:
        t = truth[r["taxonomy"].split(";")[-1]]
        t[0] += int(r["read_pairs"])
        t[1] += float(r["relative_abundance"])
        t[2] += float(r["relative_abundance"]) * float(r["genome_length"])
        t[3] += float(r["vertical_coverage"])
    cells = sum(t[1] for t in truth.values())
    weight = sum(t[2] for t in truth.values())  # a drawn sample's bases, by genome
    if drawn:
        pe = read_composition(community_composition) if community_composition else None
        pairs = sum(t[0] for t in truth.values())
        host = max(0.0, 1 - pairs / number(pe["ScannedFragments"])) if pe and number(pe["ScannedFragments"]) > 0 else NAN
        read_share = {s: (1 - host) * t[2] / weight for s, t in truth.items()}
        depth = {s: scanned_bases * (1 - host) * t[1] / weight for s, t in truth.items()}
    else:
        pairs = sum(t[0] for t in truth.values())
        host = max(0.0, 1 - pairs / scanned_fragments) if scanned_fragments > 0 else NAN
        read_share = {s: t[0] / scanned_fragments for s, t in truth.items()}
        depth = {s: t[3] * (share if read_type == "pe" else 0.5) for s, t in truth.items()}
    present = set(truth)
    in_db = present - heldout
    names = {t: v[0] for t, v in called.items()}
    true_called = {n for n in names.values() if n in present}
    ags_called_cells = sum(truth[s][1] for s in true_called)
    ratios = []
    for t, (name, _, d, size) in called.items():
        if name in present and truth[name][1] > 0:
            length = truth[name][2] / truth[name][1]
            ratios.append((d / depth[name] if depth[name] > 0 else NAN, size / length if size > 0 else NAN))
    median = lambda values: statistics.median(v for v in values if not math.isnan(v)) if any(not math.isnan(v) for v in values) else NAN
    line = {"set": set_name, "read_type": read_type, "point": point, "scenario": scenario, "sample": sample,
            "host_share": host, "species": len(present), "absent_species": len(present & heldout), "called": len(called),
            "true_called": len(true_called), "missing_true": len(present - true_called),
            "scanned_fragments": scanned_fragments,
            "true_explained": sum(read_share[s] for s in in_db), "true_explained_by_calls": sum(read_share[s] for s in true_called),
            "explained": est["explained"],
            "true_unknown": sum(truth[s][1] for s in present - in_db) / cells,
            "true_unknown_given_calls": sum(truth[s][1] for s in present - true_called) / cells, "unknown": est["unknown"],
            "true_ags": weight / cells, "true_ags_called": sum(truth[s][2] for s in true_called) / ags_called_cells
            if ags_called_cells > 0 else NAN, "ags": est["ags"], "ags_p10": est["quantiles"][0], "ags_p90": est["quantiles"][4],
            "missing_at_median": est["missing_at_median"], "missing_at_lowest": est["missing_at_lowest"],
            "depth_ratio": median([r[0] for r in ratios]), "size_ratio": median([r[1] for r in ratios])}
    # Per species (--species-out): every species present or called, with protal's depth and genome size where the
    # profile has the taxon (called or not) and the truth where it is present.
    by_name = {}
    for t, (name, _, d, size) in taxa.items():
        by_name.setdefault(name, (t, d, size))
    genomes_of = collections.defaultdict(list)
    for r in rows:
        genomes_of[r["taxonomy"].split(";")[-1]].append(r.get("genome", ""))
    called_names = set(names.values())
    species_rows = []
    for name in sorted(present | called_names):
        here = name in present
        taxid, d, size = by_name.get(name, ("", NAN, NAN))
        t = truth[name] if here else None
        species_rows.append({
            "set": set_name, "read_type": read_type, "point": point, "scenario": scenario, "sample": sample,
            "species": name, "taxid": taxid, "genomes": ",".join(genomes_of.get(name, [])), "present": int(here),
            "in_database": int(name not in heldout), "called": int(name in called_names),
            "true_read_pairs": t[0] if here else NAN, "true_cell_share": t[1] / cells if here and cells > 0 else NAN,
            "true_genome_length": t[2] / t[1] if here and t[1] > 0 else NAN, "true_depth": depth[name] if here else NAN,
            "depth": d, "genome_size": size if size > 0 else NAN})
    return (line, ratios, species_rows), None


def spread(values):
    """'median (10-90%: low to high, n)' of the finite values, or '-'."""
    values = sorted(v for v in values if not math.isnan(v))
    if not values:
        return "-"
    pick = lambda q: values[min(len(values) - 1, int(q * len(values)))]
    return f"{statistics.median(values):+.4f} (10-90%: {pick(0.1):+.4f} to {pick(0.9):+.4f}, n={len(values)})"


def summary(lines, ratios_of, label):
    """The summary of a set's samples of a read type: the errors' medians and ranges."""
    out = [f"{label}: {len(lines)} samples"]
    for which, chosen in (("without a host", [l for l in lines if l["host_share"] < HOST_FREE]),
                          ("with a host's reads", [l for l in lines if l["host_share"] >= HOST_FREE]),
                          ("host share not known", [l for l in lines if math.isnan(l["host_share"])])):
        if not chosen:
            continue
        diff = lambda a, b: [l[a] - l[b] for l in chosen]
        rel = lambda a, b: [l[a] / l[b] - 1 if l[b] > 0 else NAN for l in chosen]
        ratios = [r for l in chosen for r in ratios_of[id(l)]]
        missing = [l["missing_at_median"] / l["missing_true"] for l in chosen if l["missing_true"] > 0]
        within = [l["missing_at_median"] <= l["missing_true"] * 2 and l["missing_true"] <= l["missing_at_lowest"]
                  for l in chosen if l["missing_true"] > 0 and not math.isnan(l["missing_at_lowest"])]
        out += [f"  {which} ({len(chosen)}):",
                f"    explained share - its truth given the calls: {spread(diff('explained', 'true_explained_by_calls'))}",
                f"    explained share - the share of the database's species: {spread(diff('explained', 'true_explained'))}",
                f"    unknown share (?) - its truth given the calls: {spread(diff('unknown', 'true_unknown_given_calls'))}",
                f"    unknown share (?) - the cells of the species the database lacks: {spread(diff('unknown', 'true_unknown'))}",
                f"    average genome size / the called present species' - 1: {spread(rel('ags', 'true_ags_called'))}",
                f"    average genome size / all cells' - 1: {spread(rel('ags', 'true_ags'))}",
                f"    species' depth / their reads' - 1: {spread([r[0] - 1 for r in ratios])}",
                f"    species' genome size / their genomes' - 1: {spread([r[1] - 1 for r in ratios])}",
                f"    species missing at the median depth / those present and not called: {spread(missing)}"
                + (f"; truth within [estimate at the median / 2, at the lowest] in {sum(within)} of {len(within)}" if within else "")]
    return out


def brief(lines):
    """One line of the medians over the test set's samples without a host (the training set's without a test set), for
    the build's console."""
    chosen = [l for l in lines if l["set"] == "test"] or lines
    chosen = [l for l in chosen if l["host_share"] < HOST_FREE]
    if not chosen:
        return "no sample without a host to judge"
    med = lambda values: statistics.median(v for v in values if not math.isnan(v)) if any(not math.isnan(v) for v in values) else NAN
    explained = med([l["explained"] - l["true_explained_by_calls"] for l in chosen])
    unknown = med([l["unknown"] - l["true_unknown_given_calls"] for l in chosen])
    ags = med([l["ags"] / l["true_ags_called"] - 1 if l["true_ags_called"] > 0 else NAN for l in chosen])
    depth = med([l["depth_ratio"] - 1 for l in chosen])
    show = lambda v, f: "-" if math.isnan(v) else f.format(v)
    return (f"median errors over {len(chosen)} {chosen[0]['set']} samples without a host: explained share "
            f"{show(explained, '{:+.3f}')}, unknown share (?) {show(unknown, '{:+.3f}')}, average genome size "
            f"{show(ags * 100, '{:+.1f}%')}, species' depth {show(depth * 100, '{:+.1f}%')}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--training", help="the training samples' collection (collect_training_data.py -o)")
    ap.add_argument("--test", help="the test samples' collection")
    ap.add_argument("--calls", help="the trainer's PREFIX.calls.tsv.gz: the calls to compose (default: the profiles' own; "
                                    "then paired-end and single-end samples only)")
    ap.add_argument("--heldout", help="heldout_species.txt: the species the training database lacks")
    ap.add_argument("--read-type", default="pe", choices=list(collect.READ_TYPES))
    ap.add_argument("--out", required=True, help="the table, one line per sample")
    ap.add_argument("--summary", help="also write the summary here")
    ap.add_argument("--species-out", help="also a table of every sample's species present or called (SPECIES_COLUMNS; "
                                          ".gz: gzipped): protal's depth and genome size against the truth")
    ap.add_argument("--threads", type=int, default=1, help="samples read at once")
    args = ap.parse_args(argv)
    collections_ = {"training": args.training, "test": args.test}
    heldout = read_heldout(args.heldout)
    if args.calls:
        samples = {k: v for k, v in samples_from_calls(args.calls, args.read_type).items() if collections_.get(k[0])}
    else:
        if args.read_type not in ("pe", "se"):
            sys.exit("--read-type pb or ont needs --calls (it names the samples of a read type)")
        samples = {}
        for name, folder in collections_.items():
            if folder:
                samples.update(samples_of_collection(folder, name, args.read_type))
    jobs = [(key, scenario, called, collections_[key[0]], heldout, args.read_type)
            for key, (scenario, called) in sorted(samples.items())]
    if args.threads > 1 and len(jobs) > 1:
        with concurrent.futures.ProcessPoolExecutor(min(args.threads, len(jobs))) as pool:
            results = list(pool.map(evaluate, jobs, chunksize=4))
    else:
        results = [evaluate(job) for job in jobs]
    lines, ratios_of, species_rows, skipped = [], {}, [], collections.Counter()
    for result, why in results:
        if result is None:
            skipped[why] += 1
            continue
        line, ratios, rows = result
        lines.append(line)
        ratios_of[id(line)] = ratios
        species_rows.extend(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    def cell(v):
        return "NA" if isinstance(v, float) and math.isnan(v) else (f"{v:.6g}" if isinstance(v, float) else v)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for line in lines:
            w.writerow([cell(line[c]) for c in COLUMNS])
    if args.species_out:
        opener = gzip.open if args.species_out.endswith(".gz") else open
        with opener(args.species_out, "wt", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(SPECIES_COLUMNS)
            for row in species_rows:
                w.writerow([cell(row[c]) for c in SPECIES_COLUMNS])
    text = [f"Composition of the {args.read_type} samples against their truth ("
            + ("the model's calls, " + os.path.basename(args.calls) if args.calls else "the profiles' own calls") + f"): {args.out}"]
    for name in SETS:
        chosen = [l for l in lines if l["set"] == name]
        if chosen:
            text += summary(chosen, ratios_of, f"{name} set")
    if skipped:
        text.append("skipped: " + ", ".join(f"{n} ({why})" for why, n in sorted(skipped.items())))
    print("\n".join(text))
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("\n".join(text) + "\n")
    return lines


if __name__ == "__main__":
    main()
