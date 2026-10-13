#!/usr/bin/env python3
"""probe_summary.py - what a sparse first stage would keep, from a protal run with --stage1_probe.

Reads the run's misc/<sample>_stage1_probe.tsv (every rule: c, T, m, the reads kept, recall on the reads with records),
misc/<sample>_stage1_probe_taxa.tsv (per taxon of the best records: its reads with records and those each rule keeps,
m = 1 and 2) and misc/<sample>_first_lookups_probe.tsv (short reads' mates by the best score of their first k lookups),
and the sample's profile (the called species), and writes a Markdown summary:

  1. the share of reads kept and the recall, by c and T, for m = 1 and 2;
  2. the cheapest rule (fewest reads kept) at a recall of at least 99%, 99.5% and 99.9%;
  3. for those rules, the called species that lose half or more of their reads with records;
  4. short reads: giving up a mate whose best score among its first k lookups is at most s.

  python3 probe_summary.py --misc OUT/misc --sample pe1 --profile OUT/pe1.profile > pe1_probe.md

--pick stride (protal since 2026-10-13) summarises the rules for a stride pick of the cores (every c-th in read order)
instead: misc/<sample>_stage1_probe_stride.tsv and the taxa table's s<c>_t<T>_m<m> columns
(docs/claude/2026-10-13-stride-pick compares the two).

docs/claude/2026-10-12-stage1-probe; the tallies: src/Core/Stage1Probe.h.
"""

import argparse
import csv
import os
import sys

LEVELS = [1, 2, 4, 8, 16, 32]
SCORES = [12, 13, 14, 15, 16]
RECALLS = [0.99, 0.995, 0.999]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def called_species(profile):
    """The lineages of the species a profile calls (its second column), without the '?' line."""
    called = set()
    with open(profile) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or line.startswith("#") or fields[0] == "?":
                continue
            called.add(fields[1])
    return called


def pct(value):
    return f"{100 * value:.2f}%"


def rules_table(rules, min_hits):
    rows = {(int(r["c"]), int(r["min_score"])): r for r in rules if int(r["min_hits"]) == min_hits}
    out = [f"| c \\ T (m = {min_hits}) | " + " | ".join(str(t) for t in SCORES) + " |",
           "|---|" + "---|" * len(SCORES)]
    for c in LEVELS:
        cells = []
        for t in SCORES:
            r = rows.get((c, t))
            cells.append(f"{pct(float(r['share_kept']))} / {pct(float(r['recall']))}" if r else "")
        out.append(f"| {c} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def cheapest(rules, recall):
    """The rule keeping the fewest reads at a recall of at least `recall` (ties: the sparser index)."""
    fitting = [r for r in rules if float(r["recall"]) >= recall]
    if not fitting:
        return None
    return min(fitting, key=lambda r: (float(r["share_kept"]), -int(r["c"])))


def species_losses(taxa, called, rule, limit=15, prefix="c"):
    """The called species of which `rule` (m 1 or 2) keeps half or fewer of the reads with records (prefix "s": the
    stride pick's columns)."""
    column = f"{prefix}{rule['c']}_t{rule['min_score']}_m{rule['min_hits']}"
    if not taxa or column not in taxa[0]:
        return None, []
    rows = [t for t in taxa if t["name"] in called]
    losing = []
    for t in rows:
        aligned, kept = int(t["aligned"]), int(t[column])
        if aligned and kept <= aligned / 2:
            losing.append((kept / aligned, aligned, kept, t["name"]))
    losing.sort()
    return len(rows), losing[:limit] if limit else losing


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--misc", required=True, help="the run's misc folder")
    parser.add_argument("--sample", required=True, help="the sample's name (its misc files' prefix)")
    parser.add_argument("--profile", help="the sample's profile (for the called species)")
    parser.add_argument("--pick", choices=("hash", "stride"), default="hash",
                        help="the rules of the hash sample of the cores (the default) or of the stride pick (protal since 2026-10-13)")
    args = parser.parse_args()

    suffix, prefix = ("_stride", "s") if args.pick == "stride" else ("", "c")
    rules = read_tsv(os.path.join(args.misc, f"{args.sample}_stage1_probe{suffix}.tsv"))
    taxa_path = os.path.join(args.misc, f"{args.sample}_stage1_probe_taxa.tsv")
    taxa = read_tsv(taxa_path) if os.path.exists(taxa_path) else []
    first_path = os.path.join(args.misc, f"{args.sample}_first_lookups_probe.tsv")
    first = read_tsv(first_path) if os.path.exists(first_path) else []
    called = called_species(args.profile) if args.profile else set()

    reads, aligned = int(rules[0]["reads"]), int(rules[0]["aligned"])
    print(f"# Stage-1 probe of {args.sample}" + (" (stride pick)" if args.pick == "stride" else "") + "\n")
    print(f"{reads:,} reads (pairs for paired-end reads), {aligned:,} of them with records ({pct(aligned / reads if reads else 0)}); "
          f"{len(called)} species called.\n")
    sample_of = "every c-th core in read order" if args.pick == "stride" else "the sparse index"
    print(f"## Reads kept / recall, by {sample_of} (1/c of the cores) and the minimum best score T\n")
    for m in (1, 2):
        print(rules_table(rules, m) + "\n")

    print("## The cheapest rule at a recall\n")
    print("Of all rules, and of those with m = 1 or 2 (the taxa table's, which give the called species' losses).\n")
    print("| recall at least | rules | c | T | m | reads kept | recall | called species | of them losing half or more of their reads |")
    print("|---|---|---|---|---|---|---|---|---|")
    chosen = []
    for target in RECALLS:
        for label, pool in (("all", rules), ("m <= 2", [r for r in rules if int(r["min_hits"]) <= 2])):
            r = cheapest(pool, target)
            if r is None:
                print(f"| {pct(target)} | {label} | none | | | | | | |")
                continue
            if int(r["min_hits"]) <= 2 and r not in chosen:
                chosen.append(r)
            total, losing = species_losses(taxa, called, r, limit=0, prefix=prefix) if int(r["min_hits"]) <= 2 else (None, [])
            print(f"| {pct(target)} | {label} | {r['c']} | {r['min_score']} | {r['min_hits']} | {pct(float(r['share_kept']))} | "
                  f"{pct(float(r['recall']))} | {total if total is not None else 'n/a'} | {len(losing) if total is not None else 'n/a'} |")
    print()
    for r in chosen:
        total, losing = species_losses(taxa, called, r, prefix=prefix)
        if not losing:
            continue
        print(f"Called species losing half or more of their reads at c={r['c']} T={r['min_score']} m={r['min_hits']} "
              f"(kept share, reads with records, kept, lineage):\n")
        for share, n, kept, name in losing:
            print(f"- {pct(share)} of {n} ({kept}): {name}")
        print()

    if first:
        print("## Short reads: giving up a mate whose first k lookups score at most s\n")
        print("Shares of all mates given up / of the mates that align lost / of the mates of aligned fragments lost.\n")
        print("| k \\ s | " + " | ".join(str(s) for s in (9, 10, 11, 12, 13)) + " |")
        print("|---|" + "---|" * 5)
        by_k = {}
        for row in first:
            by_k.setdefault(row["first_lookups"], []).append(row)
        for k, rows in by_k.items():
            cells = []
            for s in (9, 10, 11, 12, 13):
                low = [r for r in rows if int(r["best_score"]) <= s]
                def total(rs, cols):
                    return sum(int(r[c]) for r in rs for c in cols)
                cols_all = ["mates_fragment_unaligned", "mates_other_mate_aligned", "mates_aligned"]
                all_mates, low_mates = total(rows, cols_all), total(low, cols_all)
                aligned_mates, low_aligned = total(rows, ["mates_aligned"]), total(low, ["mates_aligned"])
                frag_mates = total(rows, ["mates_other_mate_aligned", "mates_aligned"])
                low_frag = total(low, ["mates_other_mate_aligned", "mates_aligned"])
                share = lambda a, b: pct(a / b) if b else "n/a"
                cells.append(f"{share(low_mates, all_mates)} / {share(low_aligned, aligned_mates)} / {share(low_frag, frag_mates)}")
            print(f"| {k} | " + " | ".join(cells) + " |")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
