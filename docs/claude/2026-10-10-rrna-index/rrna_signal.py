#!/usr/bin/env python3
"""Would 16S reads of a lineage the database lacks still be placed? A held-out test on GTDB's 16S sequences.

    python3 rrna_signal.py ~/GTDB/r226 --sample 1000 --workers 4 > signal_output.txt    # ~15 minutes on 4 cores

The pool is one full-length 16S per GTDB species (MIMAG length, from genomic_files_all/ssu_all: the representative's
own longest hit where it is full length, else the species' longest), oriented by the 515F primer site. For a random
sample of pool species S and each level (species, genus, family, order), the clade of S at that level is held out (all
its pool sequences removed) and S's 16S is compared with what is left:

- full length: the identity of S's 16S (100 bp trimmed off each end) to its nearest relative inside the parent clade
  (same genus for a held-out species, same family for a held-out genus, ...), among the 8 parent-clade sequences that
  share the most 12-mers with it;
- reads: error-free 150 bp windows of S's 16S every 50 bp, each aligned to its candidates: the 5 sequences sharing the
  most 12-mers with it in the parent clade and in each of the two clades above (held-out clade excluded), and the 10
  best overall. A read is alignable when its best identity is >= 0.90 (protal's --max_score_ani default); the LCA of its
  best hits (fewest edits: strict; within 1 edit: loose) is then right to the parent and stops there (A), right to the
  parent with all hits in one sibling clade (B: a classifier would name a known sibling), right only above the parent
  (C), or in a wrong clade (D).

Identity = 1 - edits / pattern length, edits from Myers' bit-parallel algorithm (the pattern aligned end to end, free
ends in the text). Estimates are optimistic where the candidate pre-selection misses a closer outside hit (it then
under-counts misplaced reads) and ignore sequencing errors, intragenomic 16S copies and the species without a 16S.
"""

import argparse
import collections
import gzip
import multiprocessing
import os
import random
import re
import sys

import numpy as np
import pandas as pd

K = 12
FULL16 = {"Bacteria": 1200, "Archaea": 900}
LEVELS = [("species", 6), ("genus", 5), ("family", 4), ("order", 3)]  # held-out rank index; the parent is one up
READ, STEP = 150, 50
TOP_FULL, TOP_READ, TOP_GLOBAL = 8, 5, 10
RANK_NAMES = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PRIMER_F = re.compile(r"GTG[CT]CAGC[AC]GCCGCGGTAA")
PRIMER_R = re.compile(r"TTACCGCGGC[GT]GCTG[AG]CAC")
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate("ACGT"):
    LUT[ord(b)] = i
COMP = str.maketrans("ACGTN", "TGCAN")

# set in main() before the worker processes fork
SEQS, TAX, ALL_K, ALL_I, N = None, None, None, None, 0


def edits_in(pattern, text):
    """Fewest edits that turn pattern into a substring of text (Myers 1999, bit-parallel, Python integers)."""
    m = len(pattern)
    peq = {}
    for i, c in enumerate(pattern):
        if c in "ACGT":
            peq[c] = peq.get(c, 0) | (1 << i)
    mask = (1 << m) - 1
    high = 1 << (m - 1)
    pv, mv, score, best = mask, 0, m, m
    for c in text:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = (ph << 1) & mask
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
        if score < best:
            best = score
    return best


def kmers(seq):
    c = LUT[np.frombuffer(seq.encode("ascii"), np.uint8)]
    n = len(c) - K + 1
    if n <= 0:
        return np.empty(0, np.uint32)
    v = np.zeros(n, np.uint32)
    bad = np.zeros(n, bool)
    for t in range(K):
        ct = c[t:t + n]
        v = (v << np.uint32(2)) | (ct & 3).astype(np.uint32)
        bad |= ct == 255
    return np.unique(v[~bad])


def shared(qk):
    """12-mers each pool sequence shares with the query's."""
    left = np.searchsorted(ALL_K, qk, "left")
    right = np.searchsorted(ALL_K, qk, "right")
    lens = right - left
    total = int(lens.sum())
    if total == 0:
        return np.zeros(N, np.int64)
    pos = np.repeat(left - (np.cumsum(lens) - lens), lens) + np.arange(total)
    return np.bincount(ALL_I[pos], minlength=N)


def top(cands, counts, n):
    if len(cands) <= n:
        return cands
    part = np.argpartition(-counts[cands], n - 1)[:n]
    return cands[part]


def load_pool(release):
    meta = []
    for dom in ("bac120", "ar53"):
        path = os.path.join(release, f"{dom}_metadata_r{release_number(release)}.tsv.gz")
        meta.append(pd.read_csv(path, sep="\t", usecols=["accession", "gtdb_genome_representative", "gtdb_taxonomy"],
                                low_memory=False))
    meta = pd.concat(meta, ignore_index=True)
    species_of = dict(zip(meta["accession"], meta["gtdb_genome_representative"]))
    taxonomy = dict(zip(meta["gtdb_genome_representative"], meta["gtdb_taxonomy"]))
    best = {}  # species -> (own representative's, length, sequence)
    path = os.path.join(release, "genomic_files_all", f"ssu_all_r{release_number(release)}.fna.gz")
    acc, seq = None, []

    def keep():
        if acc is None:
            return
        s = "".join(seq).upper()
        sp = species_of.get(acc)
        if sp is None:
            return
        dom = taxonomy[sp].split(";")[0][3:]
        if not FULL16[dom] <= len(s) <= 2000:
            return
        key = (acc == sp, len(s))
        if sp not in best or key > best[sp][:2]:
            best[sp] = (key[0], key[1], s)

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                keep()
                acc, seq = line[1:line.index("~")], []
            else:
                seq.append(line.strip())
        keep()
    species = sorted(best)
    seqs, flipped, unoriented = [], 0, 0
    for sp in species:
        s = best[sp][2]
        if not PRIMER_F.search(s):
            if PRIMER_R.search(s):
                s = s.translate(COMP)[::-1]
                flipped += 1
            else:
                unoriented += 1
        seqs.append(s)
    print(f"# pool: {len(species):,} species with a full-length 16S ({sum(best[sp][0] for sp in species):,} the "
          f"representative's own); reverse-complemented {flipped:,}, without the 515F site {unoriented:,}", file=sys.stderr)
    return species, seqs, [taxonomy[sp].split(";") for sp in species]


def release_number(release):
    m = re.search(r"r(\d+)", os.path.basename(os.path.normpath(release)))
    return m.group(1) if m else sys.exit(f"cannot tell the release number from {release}")


def placement(s, hits):
    """(deepest rank the hits' LCA gets right, the LCA's rank): the LCA is wrong below the first when the second is
    deeper. Rank indices 0 (domain) to 6 (species); -1 when even the domain is wrong."""
    hits = np.asarray(hits)
    lca = right = -1
    for r in range(7):
        names = TAX[r][hits]
        if (names != names[0]).any():
            break
        lca = r
        if names[0] == TAX[r][s] and right == r - 1:
            right = r
    return right, lca


def run_one(s):
    q = SEQS[s]
    counts = shared(kmers(q))
    pattern = q[100:-100]
    windows = [q[w:w + READ] for w in range(0, len(q) - READ + 1, STEP)]
    wcounts = [shared(kmers(r)) for r in windows]
    out = {}
    for level, h in LEVELS:
        kept = TAX[h] != TAX[h][s]
        inside = np.flatnonzero(kept & (TAX[h - 1] == TAX[h - 1][s]))
        if len(inside) == 0:
            out[level] = None
            continue
        full = min(edits_in(pattern, SEQS[c]) for c in top(inside, counts, TOP_FULL))
        # each read's candidates: the best by shared 12-mers in its parent clade and in the 2 clades above, and overall
        groups = [inside] + [np.flatnonzero(kept & (TAX[a] == TAX[a][s])) for a in range(h - 2, max(h - 4, 0), -1)]
        everything = np.flatnonzero(kept)
        reads = []
        for r, wc in zip(windows, wcounts):
            cands = set(top(everything, wc, TOP_GLOBAL).tolist())
            for g in groups:
                cands.update(top(g, wc, TOP_READ).tolist())
            eds = {c: edits_in(r, SEQS[c]) for c in cands}
            best = min(eds.values())
            strict = placement(s, [c for c, e in eds.items() if e == best])
            loose = placement(s, [c for c, e in eds.items() if e <= best + 1])
            reads.append((best, strict, loose))
        out[level] = (1 - full / len(pattern), reads)
    return s, out


def summarise(results, n_pool):
    print(f"pool species {n_pool:,}; sampled {len(results):,}")
    print()
    print("## Full-length identity of a held-out clade's 16S to its nearest relative in the parent clade")
    print("level held out | parent | with a relative | none in parent | p10 | p25 | median | p75 | p90 | >= 98.7% | "
          ">= 94.5% | >= 86.5%")
    parents = {"species": "genus", "genus": "family", "family": "order", "order": "class"}
    for level, _ in LEVELS:
        vals = [r[level] for r in results if r[level] is not None]
        none = sum(1 for r in results if r[level] is None)
        ids = np.array([v[0] for v in vals])
        if not len(ids):
            continue
        q = np.percentile(ids, [10, 25, 50, 75, 90]) * 100
        print(f"{level} | {parents[level]} | {len(ids)} | {none} | " + " | ".join(f"{x:.1f}" for x in q) + " | " +
              " | ".join(f"{100 * (ids >= t).mean():.0f}%" for t in (0.987, 0.945, 0.865)))
    print()
    print("## 150 bp reads of the held-out clade (error-free windows, every 50 bp): identity to the best remaining 16S")
    print("level held out | reads | identity p10 | p25 | median | >= 0.97 | >= 0.90 (alignable) | >= 0.85")
    for level, _ in LEVELS:
        reads = [x for r in results if r[level] is not None for x in r[level][1]]
        if not reads:
            continue
        ident = 1 - np.array([x[0] for x in reads]) / READ
        print(f"{level} | {len(reads):,} | {np.percentile(ident, 10):.3f} | {np.percentile(ident, 25):.3f} | "
              f"{np.median(ident):.3f} | {100 * (ident >= 0.97).mean():.0f}% | {100 * (ident >= 0.90).mean():.0f}% | "
              f"{100 * (ident >= 0.85).mean():.0f}%")
    print()
    print("## Where a best-hit LCA puts the alignable reads (identity >= 0.90)")
    print("A = right to the parent and stops there (an unknown member of the parent: the ideal); B = right to the parent but")
    print("all best hits in one sibling clade (the LCA names a known sibling: the novel clade looks like a known one);")
    print("C = right only above the parent (shallower, not wrong); D = the best hits lie in a wrong clade (misplaced).")
    print("strict: the hits with the fewest edits; loose: within 1 edit of those. 'B identical' = B reads at identity 1.")
    print("level held out | parent | LCA | alignable reads | A | B | B identical | C | D")
    for level, h in LEVELS:
        reads = [x for r in results if r[level] is not None for x in r[level][1] if 1 - x[0] / READ >= 0.90]
        if not reads:
            continue
        for which, idx in (("strict", 1), ("loose", 2)):
            cat = collections.Counter()
            for x in reads:
                right, lca = x[idx]
                if right >= h - 1:
                    cat["A" if lca == h - 1 else "B"] += 1
                    if lca > h - 1 and x[0] == 0:
                        cat["Bid"] += 1
                elif lca == right:
                    cat["C"] += 1
                else:
                    cat["D"] += 1
            n = len(reads)
            print(f"{level} | {RANK_NAMES[h - 1]} | {which} | {n:,} | " +
                  " | ".join(f"{100 * cat[k] / n:.0f}%" for k in ("A", "B", "Bid", "C", "D")))
    print()
    print("## Per held-out clade: share of its reads (all windows) placed right to the parent (A + B, loose LCA), and misplaced")
    print("level held out | clades | A+B p25 | median | p75 | clades with >= 3 such reads per 16S | D median | clades with D > A+B")
    for level, h in LEVELS:
        good, bad = [], []
        for r in results:
            if r[level] is None:
                continue
            reads = r[level][1]
            ok = [x for x in reads if 1 - x[0] / READ >= 0.90]
            good.append(sum(1 for x in ok if x[2][0] >= h - 1))
            bad.append(sum(1 for x in ok if x[2][0] < h - 1 and x[2][1] > x[2][0]))
            good[-1] /= len(reads)
            bad[-1] /= len(reads)
        if good:
            good, bad = np.array(good), np.array(bad)
            n_windows = (1500 - READ) // STEP + 1
            print(f"{level} | {len(good)} | {100 * np.percentile(good, 25):.0f}% | {100 * np.median(good):.0f}% | "
                  f"{100 * np.percentile(good, 75):.0f}% | {100 * (good * n_windows >= 3).mean():.0f}% | "
                  f"{100 * np.median(bad):.0f}% | {100 * (bad > good).mean():.0f}%")


def main():
    global SEQS, TAX, ALL_K, ALL_I, N
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("release", help="GTDB release folder (e.g. ~/GTDB/r226)")
    p.add_argument("--sample", type=int, default=1000, help="species to hold out (default 1000)")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    assert edits_in("ACGT", "TTACGTTT") == 0 and edits_in("ACGT", "TTAGGTTT") == 1 and edits_in("ACGT", "AC") == 2
    release = os.path.expanduser(args.release)
    species, SEQS, lineages = load_pool(release)
    N = len(SEQS)
    TAX = []
    for rank in range(7):
        TAX.append(pd.factorize(pd.Series([";".join(t[:rank + 1]) for t in lineages]))[0])
    sizes = [0] * N
    ks = []
    for i, s in enumerate(SEQS):
        k = kmers(s)
        ks.append(k)
        sizes[i] = len(k)
    ALL_I = np.repeat(np.arange(N, dtype=np.uint32), sizes)
    ALL_K = np.concatenate(ks)
    order = np.argsort(ALL_K, kind="stable")
    ALL_K, ALL_I = ALL_K[order], ALL_I[order]
    del ks, order
    print(f"# index: {len(ALL_K):,} 12-mers", file=sys.stderr)
    sample = random.Random(args.seed).sample(range(N), min(args.sample, N))
    results = [None] * len(sample)
    with multiprocessing.get_context("fork").Pool(args.workers) as pool:
        for done, (s, out) in enumerate(pool.imap_unordered(run_one, sample, chunksize=4), 1):
            results[sample.index(s)] = out
            if done % 50 == 0:
                print(f"# {done}/{len(sample)}", file=sys.stderr, flush=True)
    summarise(results, N)


if __name__ == "__main__":
    main()
