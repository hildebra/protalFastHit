#!/usr/bin/env python3
"""score_stress.py [STRESS_DIR] - abundance under each depth-identity rule of run.sh.

Truth (design/design.tsv): a species' depth is the sum of its genomes' depths; the species present are those
of the sample, less the held-out ones against db_missing. Profiles are of --knob 0 (every taxon with reads),
so the scores are of abundance, not of detection:
  bray      Bray-Curtis between the true relative abundances of the species present and the reported ones
            renormalised over them (a species present but not reported counts as 0)
  leak      the share of all reported abundance on taxa that are not present (relatives' reads: on held-out
            species' congeners that are not in the sample, or misplaced reads)
  by kind   median log2(reported / true) of the species present, both renormalised over them: negative =
            undercounted. Kinds: the representative alone; one strain alone, by its distance from the
            reference; the representative with a distant minor strain (close_major) and the reverse
            (far_major); against db_missing also the species present whose congener is held out and in the
            sample (its reads land on them).
"""
import collections
import csv
import glob
import math
import os
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/stress")

design = collections.defaultdict(list)
with open(os.path.join(B, "design", "design.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        design[row["sample"]].append(row)
with open(os.path.join(B, "heldout.txt")) as fh:
    heldout = {line.strip().removeprefix("s__") for line in fh if line.strip()}
genus_of = {}
with open(os.path.join(B, "gtdb", "simulation", "genomes.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        genus_of[row["gtdb_taxonomy"].split(";s__")[-1]] = row["gtdb_taxonomy"].split(";")[5]


def reported(path):
    result = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            at = next((i for i, x in enumerate(f) if x.startswith("d__")), None)
            if at is None or at + 1 >= len(f):
                continue
            try:
                result[f[at].split(";")[-1].removeprefix("s__")] = float(f[at + 1])
            except ValueError:
                pass
    return result


def kind_of(rows, species):
    first = next(r for r in rows if r["species"] == species)
    if first["kind"] != "alone_strain":
        return first["kind"]
    d = float(first["distance"])
    return "alone_strain <2%" if d < 0.02 else "alone_strain 2-4%" if d < 0.04 else "alone_strain >=4%"


def score(rule, db, sample):
    paths = glob.glob(os.path.join(B, "prof", rule, db, sample, "*.profile"))
    if not paths:
        return None
    got = reported(paths[0])
    rows = design[sample]
    depth = collections.Counter()
    for r in rows:
        depth[r["species"]] += float(r["depth"])
    present = set(depth) - (heldout if db == "missing" else set())
    in_sample_heldout = set(depth) & heldout if db == "missing" else set()
    true = {s: depth[s] / sum(depth[x] for x in present) for s in present}
    shown = sum(got.get(s, 0) for s in present)
    est = {s: got.get(s, 0) / shown if shown else 0 for s in present}
    bray = 0.5 * sum(abs(true[s] - est[s]) for s in present)
    leak = 1 - shown / sum(got.values()) if got else 0
    kinds = collections.defaultdict(list)
    for s in present:
        if est[s] <= 0:
            continue
        err = math.log2(est[s] / true[s])
        kinds[kind_of(rows, s)].append(err)
        if any(genus_of[h] == genus_of[s] for h in in_sample_heldout):
            kinds["congener of a missing species"].append(err)
    return bray, leak, kinds


def main():
    rules = sorted(os.listdir(os.path.join(B, "prof")))
    order = ["none", "top004", "top008", "top012", "gm004", "gm006", "gm008", "gq20m004", "gq20m006", "fl006f008",
             "sc004k2"]
    rules = [r for r in order if r in rules] + [r for r in rules if r not in order]
    samples = sorted(design)
    print("| rule | " + " | ".join(f"{db}, {setup}: Bray / leak" for db in ("full", "missing") for setup in ("2x100", "2x150")) + " |")
    print("|---" * 5 + "|")
    kinds_by_rule = {}
    for rule in rules:
        cells = []
        kinds_by_rule[rule] = {"full": collections.defaultdict(list), "missing": collections.defaultdict(list)}
        for db in ("full", "missing"):
            for setup in ("s100", "s150"):
                scores = [score(rule, db, s) for s in samples if s.startswith(setup)]
                scores = [x for x in scores if x]
                if not scores:
                    cells.append("-")
                    continue
                for _, _, k in scores:
                    for name, errs in k.items():
                        kinds_by_rule[rule][db][name] += errs
                cells.append(f"{statistics.mean(x[0] for x in scores):.4f} / {statistics.mean(x[1] for x in scores):.4f}")
        print(f"| {rule} | " + " | ".join(cells) + " |")
    for db in ("full", "missing"):
        names = sorted({n for r in rules for n in kinds_by_rule[r][db]})
        print(f"\nmedian log2(reported / true) by kind, db_{db} (both read lengths; n = species x samples)")
        print("| rule | " + " | ".join(f"{n} (n={len(kinds_by_rule[rules[0]][db][n])})" for n in names) + " |")
        print("|---" * (len(names) + 1) + "|")
        for rule in rules:
            k = kinds_by_rule[rule][db]
            print(f"| {rule} | " + " | ".join(f"{statistics.median(k[n]):+.3f}" if k[n] else "-" for n in names) + " |")


if __name__ == "__main__":
    main()
