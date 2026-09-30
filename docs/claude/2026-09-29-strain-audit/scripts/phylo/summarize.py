#!/usr/bin/env python3
"""Tabulate results*.jsonl files written by analyze.py.

usage: summarize.py <results.jsonl>... [--cols default|long] [--warnings]
"""
import argparse
import json


def f(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--warnings", action="store_true")
    ap.add_argument("--detail", action="store_true")
    a = ap.parse_args()
    hdr = ["run", "sp", "variant", "seqs", "cols", "pinf", "single", "RF", "nRF", "RFref", "qrt",
           "patR", "patRatio", "termR", "intR", "supT", "minT", "T95", "F", "F95", "miss%", "fSNP",
           "iupac", "snpRec", "status"]
    print("\t".join(hdr))
    for fn in a.files:
        for line in open(fn):
            r = json.loads(line)
            t = r.get("truth", {}).get("cells", {}) if isinstance(r.get("truth"), dict) else {}
            valid = t.get("valid", 0)
            miss = t.get("missing", 0)
            rec = None
            if t.get("true_snps"):
                rec = t.get("snp_called", 0) / t["true_snps"]
            row = [r.get("run", "") + (":" + r["tag"] if r.get("tag") else ""), r["species"], r["variant"],
                   f(r.get("n_seqs")), f(r.get("iq_cols")), f(r.get("iq_pinf")), f(r.get("iq_single")),
                   (f"{r['rf']}/{r['max_rf']}" if "rf" in r else "-"), f(r.get("nrf"), 2),
                   (f"{r['rf_ref']}/{r['max_rf_ref']}" if "rf_ref" in r else "-"), f(r.get("quartet"), 3),
                   f(r.get("pat_pearson")), f(r.get("pat_ratio"), 2), f(r.get("term_ratio"), 2),
                   f(r.get("int_ratio"), 2), f(r.get("true_split_support_mean"), 1),
                   f(r.get("true_split_support_min"), 0),
                   (f"{r['n_true_splits_ge95']}/{r['n_true_splits']}" if "n_true_splits" in r else "-"),
                   f(r.get("n_false_splits")), f(r.get("n_false_ge95")),
                   (f"{100*miss/(valid+miss):.2f}" if valid + miss else "-"),
                   f(t.get("false_snp", 0) + t.get("wrong_alt", 0)) if t else "-",
                   f(t.get("iupac", 0)) if t else "-", f(rec, 3), r.get("status", "")]
            print("\t".join(row))
            if a.warnings and r.get("warnings"):
                for w in r["warnings"]:
                    print("    WARN " + w[:160])
            if a.detail:
                if r.get("false_splits"):
                    print("    false:", r["false_splits"])
                if r.get("missed_true_splits"):
                    print("    missed:", r["missed_true_splits"])
                if r.get("mix_sister") is not None:
                    print("    mix sister:", r["mix_sister"], "term", r.get("mix_term_len"),
                          "median term", r.get("median_term_len"))
                q = r.get("qcmsa")
                if q:
                    print("    qcmsa:", {k: q.get(k) for k in ("samples_in", "samples_kept", "genes_in", "genes_kept",
                                                            "genes_filtered_coverage", "genes_filtered_mrate2",
                                                            "sites_in", "sites_kept",
                                                            "sites_removed_low_parsimony", "outlier_cells",
                                                            "coverage_gap_filled_cells", "samples_removed")})
                if r.get("truth", {}).get("dropped_columns"):
                    print("    dropped:", r["truth"]["dropped_columns"])


if __name__ == "__main__":
    main()
