#!/usr/bin/env python3
"""Per-sample alignment work and time of two GTDB builds' collection runs, matched by sample.

Reads OUTDIR/logs/protal_runs_{training,test}.log of two builds (A, B) and compares, sample by
sample, what protal did (k-mers looked up, flex cells scanned, candidate alignments tried, aligned
from anchors, alignments made, records written) and how long it took (the "Aligning reads took"
timer and the "Profiling sample ... took" parts). Where the work counts agree, a time ratio B/A is
the machine (and code speed), not the samples.

Usage: run_log_work.py BUILD_A_OUTDIR BUILD_B_OUTDIR
"""
import argparse
import os
import re
import statistics
from collections import defaultdict

ALIGN_HDR = re.compile(r"^Align the .* reads of sample (\S+)")
ALIGN_TOOK = re.compile(r"^Aligning reads took (.*)$")
SUMMARY = re.compile(r"^Sample (\S+): (\d+) (?:read pairs|reads), (\d+) (?:mates )?with an anchor; (\d+) candidate "
                     r"alignments tried: (\d+) refused by the k-mer screen, (\d+) by the indel bound, (\d+) aligned "
                     r"from the anchor's exact matches and (\d+) as whole windows; (\d+) alignments made, (\d+) records")
SEEDING = re.compile(r"^Sample (\S+) seeding: (\d+) k-mers looked up, (\d+) in the index, (\d+) blocks scanned with "
                     r"(\d+) flex cells")
PROFILE = re.compile(r"^Profiling sample (\S+) took ([\d.]+)s: (.*)$")
PART = re.compile(r"([a-zA-Z][a-zA-Z ]+?) ([\d.]+)s(?:,|$)")


def seconds(text):
    """'1m 3s 20ms', '17s 499ms', '30ms' -> seconds."""
    total = 0.0
    for value, unit in re.findall(r"(\d+)(ms|s|m|h)\b", text):
        total += int(value) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600}[unit]
    return total


def read_type(sample):
    if sample.endswith("_se") or "_se_" in sample:
        return "se"
    if "ont_b" in sample:
        return "ont"
    if "pb_b" in sample:
        return "pb"
    return "pe"


def parse(outdir):
    samples = defaultdict(dict)
    for collection in ("training", "test"):
        path = os.path.join(outdir, "logs", f"protal_runs_{collection}.log")
        current = None
        with open(path, errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\n")
                m = ALIGN_HDR.match(line)
                if m:
                    current = (collection, m.group(1))
                    samples[current]["collection"] = collection
                    continue
                m = ALIGN_TOOK.match(line)
                if m and current:
                    samples[current]["align_s"] = seconds(m.group(1))
                    continue
                m = SUMMARY.match(line)
                if m:
                    s = samples[(collection, m.group(1))]
                    keys = ("reads", "anchored", "tried", "screened", "indel", "from_anchor", "windows", "made",
                            "records")
                    for key, value in zip(keys, m.groups()[1:]):
                        s[key] = int(value)
                    continue
                m = SEEDING.match(line)
                if m:
                    s = samples[(collection, m.group(1))]
                    for key, value in zip(("kmers", "in_index", "blocks", "flex"), m.groups()[1:]):
                        s[key] = int(value)
                    continue
                m = PROFILE.match(line)
                if m:
                    s = samples[(collection, m.group(1))]
                    s["profile_s"] = float(m.group(2))
                    for name, value in PART.findall(m.group(3)):
                        s["prof " + name.strip()] = float(value)
    return samples


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("a")
    parser.add_argument("b")
    args = parser.parse_args()
    a, b = parse(args.a), parse(args.b)
    common = sorted(set(a) & set(b))
    print(f"samples: A {len(a)}, B {len(b)}, both {len(common)}")

    work_keys = ("reads", "kmers", "flex", "tried", "screened", "indel", "from_anchor", "windows", "made", "records")
    print("\n## Work: samples whose count is identical in both builds, and B/A of the sums")
    print("type  " + " ".join(f"{k:>12}" for k in work_keys))
    by_type = defaultdict(list)
    for name in common:
        by_type[read_type(name[1])].append(name)
    for t in ("pe", "se", "pb", "ont"):
        names = by_type[t]
        same = []
        ratio = []
        for k in work_keys:
            pairs = [(a[n].get(k), b[n].get(k)) for n in names if a[n].get(k) is not None and b[n].get(k) is not None]
            same.append(sum(1 for x, y in pairs if x == y) / len(pairs) if pairs else float("nan"))
            sa, sb = sum(x for x, _ in pairs), sum(y for _, y in pairs)
            ratio.append(sb / sa if sa else float("nan"))
        print(f"{t:<5} " + " ".join(f"{s:>11.1%} " for s in same) + "  (identical)")
        print(f"{'':<5} " + " ".join(f"{r:>12.4f}" for r in ratio) + "  (B/A)")

    print("\n## Time, B/A: sums and the median of per-sample ratios (samples over 1 s in A)")
    time_keys = ["align_s", "profile_s"] + sorted({k for n in common for k in a[n] if k.startswith("prof ")})
    for t in ("pe", "se", "pb", "ont"):
        names = by_type[t]
        print(f"{t}: {len(names)} samples")
        for k in time_keys:
            pairs = [(a[n][k], b[n][k]) for n in names if k in a[n] and k in b[n]]
            sa, sb = sum(x for x, _ in pairs), sum(y for _, y in pairs)
            big = [y / x for x, y in pairs if x >= 1.0]
            med = statistics.median(big) if big else float("nan")
            print(f"  {k:<42} A {sa:9.1f} s  B {sb:9.1f} s  B/A {sb / sa if sa else float('nan'):6.3f}  "
                  f"per-sample median {med:6.3f} (n={len(big)})")

    print("\n## Aligning time per unit of work (B/A), samples with identical work counts")
    for t in ("pe", "se", "pb", "ont"):
        names = [n for n in by_type[t]
                 if all(a[n].get(k) == b[n].get(k) for k in ("kmers", "tried", "made")) and a[n].get("align_s", 0) >= 1]
        if not names:
            print(f"{t}: no sample with identical work and >= 1 s")
            continue
        sa = sum(a[n]["align_s"] for n in names)
        sb = sum(b[n]["align_s"] for n in names)
        print(f"{t}: {len(names)} samples, aligning {sa:.1f} -> {sb:.1f} s, B/A {sb / sa:.3f}")


if __name__ == "__main__":
    main()
