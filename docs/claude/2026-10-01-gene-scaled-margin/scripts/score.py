#!/usr/bin/env python3
"""score.py WORLD_DIR [PROF_SUBDIR] - abundance under each depth rule of run_rules.sh (profiles in
WORLD_DIR/PROF_SUBDIR, default prof7, one folder per rule).

Truth (design/design.tsv): a species' depth is the sum of its genomes' depths; the species present are those
of the sample, less the held-out ones against db_missing. The profiles are of --knob 0 (every taxon with reads),
so the scores are of abundance, not of detection:
  bray      Bray-Curtis between the true relative abundances of the species present and the reported ones
            renormalised over them (a species present but not reported counts as 0)
  leak      the share of all reported abundance on taxa that are not present
  by kind   median log2(reported / true) of the species present, both renormalised over them (negative:
            undercounted), by what the sample holds of the species (alone_ref: the representative alone;
            alone_strain by its distance from the reference; close_major: the representative with a distant
            minor strain; far_major the reverse) and by its congeners in the sample: the closest one's distance
            at the markers (the two species' divergences from their genus' ancestor, at a gene of factor 1),
            and against db_missing whether a held-out congener is in the sample, and one with 3 times the
            species' depth or more ("outnumbered"). mean |log2| is given beside the median for these.
"""
import collections
import csv
import glob
import math
import os
import statistics
import sys

B = os.path.expanduser(sys.argv[1])
PROF = sys.argv[2] if len(sys.argv) > 2 else "prof7"
ORDER = ["none", "top004", "top008", "gtop3", "gtop0", "gtop2", "gtop4", "gtop3m10", "gm008", "gmedc33", "gmedc32",
         "gsplit13s", "gsplit08s", "gsplit07s", "gsplit08u", "gcong05a", "gcong05d"]

design = collections.defaultdict(list)
with open(os.path.join(B, "design", "design.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        design[row["sample"]].append(row)
with open(os.path.join(B, "heldout.txt")) as fh:
    heldout = {line.strip().removeprefix("s__") for line in fh if line.strip()}
genus_of, species_div = {}, {}
with open(os.path.join(B, "gtdb", "simulation", "genomes.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        genus_of[row["gtdb_taxonomy"].split(";s__")[-1]] = row["gtdb_taxonomy"].split(";")[5]
with open(os.path.join(B, "gtdb", "simulation", "divergence.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        species_div[row["species"].removeprefix("s__")] = float(row["species_divergence"])


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


def congener_kinds(s, depth, db):
    """The congener categories of species s in a sample with species depths `depth`."""
    kinds = []
    congeners = [c for c in depth if c != s and genus_of[c] == genus_of[s]]
    if species_div and congeners:
        closest = min(species_div[s] + species_div[c] for c in congeners)
        kinds.append("closest congener <4%" if closest < 0.04 else "closest congener 4-7%" if closest < 0.07
                     else "closest congener >=7%")
    elif not congeners:
        kinds.append("no congener")
    if db == "missing":
        missing = [c for c in congeners if c in heldout]
        if missing:
            kinds.append("held-out congener present")
            if any(depth[c] >= 3 * depth[s] for c in missing):
                kinds.append("outnumbered by a held-out congener")
    return kinds


def score(rule, db, sample):
    paths = glob.glob(os.path.join(B, PROF, rule, db, sample, "*.profile"))
    if not paths:
        return None
    got = reported(paths[0])
    rows = design[sample]
    depth = collections.Counter()
    for r in rows:
        depth[r["species"]] += float(r["depth"])
    present = set(depth) - (heldout if db == "missing" else set())
    true = {s: depth[s] / sum(depth[x] for x in present) for s in present}
    shown = sum(got.get(s, 0) for s in present)
    est = {s: got.get(s, 0) / shown if shown else 0 for s in present}
    bray = 0.5 * sum(abs(true[s] - est[s]) for s in present)
    leak = 1 - shown / sum(got.values()) if got else 0
    kinds, congeners = collections.defaultdict(list), collections.defaultdict(list)
    for s in present:
        if est[s] <= 0:
            continue
        err = math.log2(est[s] / true[s])
        kinds[kind_of(rows, s)].append(err)
        for k in congener_kinds(s, depth, db):
            congeners[k].append(err)
    return bray, leak, kinds, congeners


def table(title, rules, by_rule, db, with_abs):
    names = sorted({n for r in rules for n in by_rule[r][db]})
    if not names:
        return
    print(f"\n{title}, db_{db} (both read lengths; n = species x samples)")
    print("| rule | " + " | ".join(f"{n} (n={len(by_rule[rules[0]][db][n])})" for n in names) + " |")
    print("|---" * (len(names) + 1) + "|")
    for rule in rules:
        k = by_rule[rule][db]
        cell = (lambda v: f"{statistics.median(v):+.3f} / {statistics.mean(abs(x) for x in v):.3f}") if with_abs else \
               (lambda v: f"{statistics.median(v):+.3f}")
        print(f"| {rule} | " + " | ".join(cell(k[n]) if k[n] else "-" for n in names) + " |")


def main():
    found = sorted(os.listdir(os.path.join(B, PROF)))
    rules = [r for r in ORDER if r in found] + [r for r in found if r not in ORDER]
    samples = sorted(design)
    print("| rule | " + " | ".join(f"{db}, {setup}" for db in ("full", "missing") for setup in ("2x100", "2x150")) +
          " | sum | leak, full / missing |")
    print("|---" * 7 + "|")
    kinds_by_rule, congeners_by_rule = {}, {}
    for rule in rules:
        cells, total, leaks = [], 0.0, {}
        kinds_by_rule[rule] = {db: collections.defaultdict(list) for db in ("full", "missing")}
        congeners_by_rule[rule] = {db: collections.defaultdict(list) for db in ("full", "missing")}
        for db in ("full", "missing"):
            leak = []
            for setup in ("s100", "s150"):
                scores = [x for x in (score(rule, db, s) for s in samples if s.startswith(setup)) if x]
                if not scores:
                    cells.append("-")
                    continue
                for _, _, k, c in scores:
                    for name, errs in k.items():
                        kinds_by_rule[rule][db][name] += errs
                    for name, errs in c.items():
                        congeners_by_rule[rule][db][name] += errs
                bray = statistics.mean(x[0] for x in scores)
                total += bray
                leak += [x[1] for x in scores]
                cells.append(f"{bray:.4f}")
            leaks[db] = f"{statistics.mean(leak):.4f}" if leak else "-"
        print(f"| {rule} | " + " | ".join(cells) + f" | {total:.3f} | {leaks['full']} / {leaks['missing']} |")
    for db in ("full", "missing"):
        table("median log2(reported / true) by kind", rules, kinds_by_rule, db, False)
    for db in ("full", "missing"):
        table("median log2(reported / true) / mean |log2| by congeners in the sample", rules, congeners_by_rule, db, True)


if __name__ == "__main__":
    main()
