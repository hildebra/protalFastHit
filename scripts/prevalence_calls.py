#!/usr/bin/env python3
"""Species calls adjusted by their prevalence across the samples of a run (a postprocessing step).

The samples of a study share an environment: a species present in most of them is likelier present in the one at
hand than the model's training prior says, and one seen nowhere else less so (DADA2's pseudo-pooling, in the terms of
protal's scores). For every taxon the script takes its probabilities in the run's other samples (never the sample
itself, so a call cannot feed on itself), smooths their mean towards the run's base rate with a few pseudo-samples,
and moves the sample's probability by the odds ratio of that prevalence to the base rate:

    odds(p') = odds(p) * odds(prevalence) / odds(base rate)

where the base rate is the model's training prior if the run's database model states it (--prior), else the mean
probability over every taxon with records in the run. The odds ratio is capped (--max-odds-ratio, 4: a prevalence
can at most quadruple or quarter the odds). By default (--direction up) only the boosts of prevalent species are
applied; --direction both applies the cuts of rare ones too, the full update; down only the cuts. A taxon with
records in one sample only keeps its probability. Calls are then made at --knob (0.5, protal's default without a
knob curve), and the taxa that change are listed.

The adjustment is a Bayesian update under the assumption that the run's samples share species: a species present
in most of them is likelier in the one at hand, one seen nowhere else less so. Where the samples do not share
species the prior is wrong and the full update cuts recall: on 42 simulated test samples drawn independently,
--direction both took 24% of the calls away (23% with the cap), nearly all of them true, while the default changed
nothing there (docs/claude/2026-10-03-false-positive-fixes). So use --direction both only on the samples of one study,
check `--check`'s numbers (how many calls flipped, how far the probabilities moved), and validate on a study with a
truth before trusting it. The default boosts prevalent species without cutting rare ones, at the price of boosting a
recurring artefact too (scripts/recurrent_calls.py finds those).

Reads every <sample>.profile.log of a protal run (protal -o: every taxon the sample's reads hit, with Predicted,
Probability, Lineage, Name and TaxID). Writes prevalence_calls.tsv (sample, taxid, species, hits, p, prevalence,
p_adjusted, called, called_adjusted) and, with --profiles DIR, one <sample>.profile per sample in protal's format
(representative genome, lineage, abundance) with the adjusted calls, abundances renormalised over them; where the
sample's own <sample>.profile beside its .profile.log ends with the unknown share ("?", protal since 2026-10-08), the
adjusted calls share the rest and the "?" line is kept.

Usage: prevalence_calls.py OUTPUT_DIR [more dirs or .profile.log files] [-o prevalence_calls.tsv] [--knob 0.5]
       [--prior P] [--pseudo-samples 2] [--profiles DIR] [--check]
"""
import argparse
import csv
import glob
import math
import os
import sys
from collections import defaultdict


def total_hits(summary):
    key = "Total Hits:"
    at = summary.find(key)
    if at < 0:
        return None
    digits = ""
    for c in summary[at + len(key):].strip():
        if c.isdigit():
            digits += c
        else:
            break
    return int(digits) if digits else None


def read_profile_log(path):
    """[(taxid, name, lineage, rep_genome, abundance, hits, called, probability)] of a .profile.log."""
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                taxid = row.get("TaxID", "").strip()
                if not taxid:
                    continue
                rows.append((taxid, row.get("Name", "").strip(), row.get("Lineage", "").strip(), row.get("RepGenome", "").strip(),
                             float(row.get("Abundance", "nan") or "nan"), total_hits(row.get("Summary", "") or ""),
                             str(row.get("Predicted", "")).strip().lower() in ("1", "true", "yes"), float(row.get("Probability", "nan"))))
            except (TypeError, ValueError):
                continue
    return rows


def unknown_share(log_path):
    """The unknown share ("?" line) of the <sample>.profile beside a <sample>.profile.log, or None without one."""
    path = log_path[: -len(".log")]
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[1] == "?":
                try:
                    return float(fields[2])
                except ValueError:
                    return None
    return None


def profile_logs(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "**", "*.profile.log"), recursive=True))
        elif p.endswith(".profile.log"):
            files.append(p)
        else:
            sys.exit(f"{p}: not a folder nor a .profile.log")
    if not files:
        sys.exit("no .profile.log files found")
    return files


def odds(p):
    p = min(max(p, 1e-9), 1 - 1e-9)
    return p / (1 - p)


def from_odds(o):
    return o / (1 + o)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="protal output folders (searched for *.profile.log) or .profile.log files")
    ap.add_argument("-o", "--output", default="prevalence_calls.tsv")
    ap.add_argument("--knob", type=float, default=0.5, help="the threshold for the adjusted calls (default 0.5)")
    ap.add_argument("--prior", type=float, default=None,
                    help="the model's base rate, the share of present taxa among candidates in its training (the model header's "
                         "protal_prior if it has one); default: the mean probability over the run's taxa with records")
    ap.add_argument("--pseudo-samples", type=float, default=2.0,
                    help="pseudo-samples at the base rate added to a taxon's prevalence, so that one or two samples move it little (default 2)")
    ap.add_argument("--max-odds-ratio", type=float, default=4.0,
                    help="the most a taxon's prevalence may multiply or divide its odds by (default 4; 0: no cap)")
    ap.add_argument("--direction", choices=("both", "up", "down"), default="up",
                    help="apply only the boosts of prevalent taxa (up, default), the boosts and the cuts (both: the full update, "
                         "for samples that share species) or only the cuts of rare ones (down)")
    ap.add_argument("--profiles", default=None, help="write <sample>.profile files with the adjusted calls into this folder")
    ap.add_argument("--check", action="store_true", help="also print how far the probabilities moved and how many calls flipped")
    args = ap.parse_args()

    files = profile_logs(args.paths)
    samples, unknown = {}, {}
    for path in files:
        sample = os.path.basename(path)[: -len(".profile.log")]
        samples[sample] = read_profile_log(path)
        unknown[sample] = unknown_share(path)
    if len(samples) < 2:
        sys.exit("prevalence needs two samples or more")
    # Every probability by taxon and sample.
    p_by_taxon = defaultdict(dict)
    for sample, rows in samples.items():
        for taxid, _, _, _, _, _, _, p in rows:
            if not math.isnan(p):
                p_by_taxon[taxid][sample] = p
    all_p = [p for d in p_by_taxon.values() for p in d.values()]
    base = args.prior if args.prior is not None else (sum(all_p) / len(all_p) if all_p else 0.5)
    base = min(max(base, 1e-6), 1 - 1e-6)
    n_samples = len(samples)

    def prevalence(taxid, sample):
        """The taxon's mean probability over the other samples (absent from a sample: 0), smoothed towards the base rate."""
        others = [p for s, p in p_by_taxon[taxid].items() if s != sample]
        missing = n_samples - 1 - len(others)  # samples without a record of the taxon: as good as 0
        total = sum(others) + base * args.pseudo_samples
        return total / (len(others) + missing + args.pseudo_samples)

    out_rows, flips_on, flips_off, moves = [], 0, 0, []
    adjusted = defaultdict(list)
    for sample, rows in sorted(samples.items()):
        for taxid, name, lineage, rep, abundance, hits, called, p in rows:
            if math.isnan(p):
                continue
            if len(p_by_taxon[taxid]) <= 1:
                prev, p_adj = float("nan"), p
            else:
                prev = prevalence(taxid, sample)
                ratio = odds(prev) / odds(base)
                if args.max_odds_ratio > 0:
                    ratio = min(max(ratio, 1 / args.max_odds_ratio), args.max_odds_ratio)
                if (args.direction == "up" and ratio < 1) or (args.direction == "down" and ratio > 1):
                    ratio = 1
                p_adj = from_odds(odds(p) * ratio)
            called_adj = p_adj >= args.knob
            flips_on += called_adj and not called
            flips_off += called and not called_adj
            moves.append(abs(p_adj - p))
            out_rows.append(dict(sample=sample, taxid=taxid, species=name, hits="" if hits is None else hits, p=f"{p:.6g}",
                                 prevalence="" if math.isnan(prev) else f"{prev:.4g}", p_adjusted=f"{p_adj:.6g}",
                                 called=int(called), called_adjusted=int(called_adj)))
            if called_adj:
                adjusted[sample].append((rep, lineage, abundance))
    with open(args.output, "w", newline="") as fh:
        fields = ["sample", "taxid", "species", "hits", "p", "prevalence", "p_adjusted", "called", "called_adjusted"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    if args.profiles:
        os.makedirs(args.profiles, exist_ok=True)
        for sample in samples:
            taxa = adjusted.get(sample, [])
            total = sum(a for _, _, a in taxa if not math.isnan(a)) or 1.0
            # The unknown share of the sample's profile stays; the adjusted calls share the rest.
            rest = 1.0 - unknown[sample] if unknown[sample] is not None else 1.0
            with open(os.path.join(args.profiles, sample + ".profile"), "w") as fh:
                for rep, lineage, abundance in taxa:
                    fh.write(f"{rep}\t{lineage}\t{(abundance / total * rest) if not math.isnan(abundance) else 0:.6g}\n")
                if unknown[sample] is not None:
                    fh.write(f"?\t?\t{unknown[sample]:.6g}\n")
    before = sum(r["called"] for r in out_rows)
    after = sum(r["called_adjusted"] for r in out_rows)
    print(f"{n_samples} samples, {len(out_rows)} taxa with records, base rate {base:.4f} ({'given' if args.prior is not None else 'the run mean'}); "
          f"calls {before} -> {after}: {flips_on} gained, {flips_off} lost; {args.output}")
    if args.check:
        moves.sort()
        q = lambda f: moves[min(len(moves) - 1, int(f * len(moves)))] if moves else float("nan")
        print(f"probability moved by median {q(0.5):.4f}, 90% {q(0.9):.4f}, max {q(1.0):.4f}; "
              f"taxa in one sample only: {sum(1 for d in p_by_taxon.values() if len(d) <= 1)} of {len(p_by_taxon)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
