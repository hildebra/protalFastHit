#!/usr/bin/env python3
"""Score protal's strain MSA cells against the simulated truth.

usage: msa_truth.py <raw.msa.fna> <raw.partition.txt> <truth_markers.tsv> [--msa X --colmap Y]
                    [--mix mix=t03] [--json]

Without --msa the raw MSA itself is scored; with --msa/--colmap (qcmsa_trace output) the filtered
MSA is scored, its columns mapped back to raw columns.

Per sample row, over reference-position columns (insertion columns are counted separately):
  valid       ACGT or IUPAC cells
  missing     '-'/N where the true strain has a base
  correct     base == truth
  wrong       base != truth (split into false_snp: truth==ref, missed_snp: called ref but truth!=ref,
              wrong_alt: both differ)
  iupac       ambiguity cells (iupac_ok: code contains the true base)
  true_snps   positions where the strain differs from the reference (in covered genes)
  snp_called  true SNP written correctly; snp_missing: written '-'/N; snp_as_ref: written as ref
Column summary (all rows): true variable columns and, for the filtered MSA, what site cleanup
removed: true singletons vs true parsimony-informative vs truly constant (error-only) columns.
"""
import argparse
import json
import re
from collections import Counter, defaultdict

IUPAC = {"R": "AG", "Y": "CT", "W": "AT", "S": "CG", "M": "AC", "K": "GT",
         "B": "CGT", "H": "ACT", "D": "AGT", "V": "ACG"}
MISSING = set("-Nn.")


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].strip())
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


def parse_partition(path):
    out = []
    for line in open(path):
        m = re.search(r"gene(\d+)\s*=\s*(\d+)-(\d+)", line)
        if m:
            out.append((m.group(1), int(m.group(2)) - 1, int(m.group(3)) - 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("partition")
    ap.add_argument("truth")
    ap.add_argument("--msa", default=None)
    ap.add_argument("--colmap", default=None)
    ap.add_argument("--mix", default=None, help="NAME=TIP: score row NAME against TIP")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    rnames, rseqs = read_fasta(a.raw)
    ref = rseqs[0]
    part = parse_partition(a.partition)
    # raw column -> (gene, refpos) or (gene, None) for insertion columns
    colinfo = {}
    glen = {}
    for g, s, e in part:
        p = 0
        for c in range(s, e + 1):
            if ref[c] != "-":
                colinfo[c] = (g, p)
                p += 1
            else:
                colinfo[c] = (g, None)
        glen[g] = p
    truth = defaultdict(dict)
    with open(a.truth) as fh:
        next(fh)
        for line in fh:
            tip, gid, marker, seq, nins = line.rstrip("\n").split("\t")
            truth[tip][gid] = seq
    # sanity: gene lengths
    anytip = next(iter(truth))
    bad = [g for g in glen if g in truth[anytip] and len(truth[anytip][g]) != glen[g]]
    if bad:
        raise SystemExit(f"gene length mismatch for genes {bad[:5]}")

    if a.msa:
        names, seqs = read_fasta(a.msa)
        cols = []
        with open(a.colmap) as fh:
            next(fh)
            for line in fh:
                k, raw_col, g = line.split("\t")
                cols.append(int(raw_col))
    else:
        names, seqs = rnames, rseqs
        cols = list(range(len(ref)))
    mixmap = {}
    if a.mix:
        n, t = a.mix.split("=")
        mixmap[n] = t

    per = {}
    col_true = defaultdict(list)  # analysed col index -> list of true bases of rows
    for name, seq in zip(names, seqs):
        if name.endswith("_reference"):
            continue
        tip = mixmap.get(name, name)
        if tip not in truth:
            continue
        c = Counter()
        for k, rc in enumerate(cols):
            g, p = colinfo[rc]
            ch = seq[k]
            if p is None:
                c["ins_cols_base" if ch not in MISSING else "ins_cols_gap"] += 1
                continue
            t = truth[tip].get(g)
            if t is None:
                continue
            t = t[p]
            r = ref[rc]
            col_true[k].append(t)
            if t != r and t != "-":
                c["true_snps"] += 1
            if ch in MISSING:
                if t == "-":
                    c["del_ok"] += 1
                else:
                    c["missing"] += 1
                    if t != r:
                        c["snp_missing"] += 1
                continue
            c["valid"] += 1
            if ch in IUPAC:
                c["iupac"] += 1
                if t in IUPAC[ch]:
                    c["iupac_ok"] += 1
                continue
            if ch == t:
                c["correct"] += 1
                if t != r:
                    c["snp_called"] += 1
            else:
                c["wrong"] += 1
                if t == "-":
                    c["base_at_true_deletion"] += 1
                elif t == r:
                    c["false_snp"] += 1
                elif ch == r:
                    c["missed_snp"] += 1
                    c["snp_as_ref"] += 1
                else:
                    c["wrong_alt"] += 1
        per[name] = dict(c)

    # column summary over analysed columns: true pattern among the scored rows
    colsum = Counter()
    for k, rc in enumerate(cols):
        g, p = colinfo[rc]
        if p is None:
            colsum["ins_cols"] += 1
            continue
        tb = [x for x in col_true.get(k, []) if x != "-"]
        r = ref[rc]
        cnt = Counter(tb + [r])
        if len(cnt) <= 1:
            colsum["true_constant"] += 1
        else:
            maj = max(cnt.values())
            minor = sum(cnt.values()) - maj
            colsum["true_singleton" if minor == 1 else "true_informative"] += 1
    res = {"rows": per, "columns": dict(colsum), "n_cols": len(cols)}
    if a.msa:
        # what site cleanup / gene removal dropped, from the raw MSA
        kept = set(cols)
        raw_true = defaultdict(list)
        raw_rows = [(n, s) for n, s in zip(rnames, rseqs) if not n.endswith("_reference")
                    and mixmap.get(n, n) in truth and n in set(names)]
        dropped = Counter()
        kept_genes = {colinfo[c][0] for c in cols}
        for rc, (g, p) in colinfo.items():
            if rc in kept or p is None:
                continue
            where = "in_kept_gene" if g in kept_genes else "in_dropped_gene"
            tb = []
            for n, s in raw_rows:
                t = truth[mixmap.get(n, n)].get(g)
                if t is not None and t[p] != "-":
                    tb.append(t[p])
            cnt = Counter(tb + [ref[rc]])
            if len(cnt) <= 1:
                dropped[where + ":true_constant"] += 1
            else:
                maj = max(cnt.values())
                minor = sum(cnt.values()) - maj
                dropped[where + (":true_singleton" if minor == 1 else ":true_informative")] += 1
        res["dropped_columns"] = dict(dropped)
    if a.json:
        print(json.dumps(res))
    else:
        keys = ["valid", "missing", "correct", "wrong", "false_snp", "missed_snp", "wrong_alt", "iupac",
                "iupac_ok", "true_snps", "snp_called", "snp_missing", "snp_as_ref"]
        print("row\t" + "\t".join(keys))
        for n, c in per.items():
            print(n + "\t" + "\t".join(str(c.get(k, 0)) for k in keys))
        print("columns", json.dumps(res["columns"]))
        if "dropped_columns" in res:
            print("dropped", json.dumps(res["dropped_columns"]))


if __name__ == "__main__":
    main()
