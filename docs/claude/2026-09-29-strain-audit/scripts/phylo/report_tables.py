#!/usr/bin/env python3
"""Compact markdown tables for the report from runs/*/results*.jsonl.
usage: report_tables.py <spec>...   spec = run[:tag]/species/variant,variant,...  (species comma list ok)"""
import json
import os
import sys

PHYLO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(run, tag):
    fn = os.path.join(PHYLO, "runs", run, f"results{'_' + tag if tag else ''}.jsonl")
    out = {}
    if os.path.exists(fn):
        for line in open(fn):
            r = json.loads(line)
            out[(r["species"], r["variant"])] = r
    return out


def fmt(x, nd=2):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


print("| run | sp | variant | tips | cols | genes | RF/max | quartet | patr r | pat ratio | term ratio | int ratio "
      "| true splits >=95 | false (>=95) | N/- % | false SNP | SNP recall |")
print("|" + "---|" * 17)
for spec in sys.argv[1:]:
    runtag, sps, vs = spec.split("/")
    run, _, tag = runtag.partition(":")
    res = load(run, tag)
    for sp in sps.split(","):
        for v in vs.split(","):
            r = res.get((sp, v))
            if r is None:
                continue
            t = (r.get("truth") or {}).get("cells", {}) if isinstance(r.get("truth"), dict) else {}
            q = r.get("qcmsa") or {}
            genes = f"{q.get('genes_kept')}/{q.get('genes_in')}" if q else "all"
            miss = (100 * t.get("missing", 0) / (t.get("valid", 0) + t.get("missing", 0))) if t.get("valid") else None
            rec = t.get("snp_called", 0) / t["true_snps"] if t.get("true_snps") else None
            if r.get("status") != "ok":
                print(f"| {runtag} | {sp} | {v} | {r.get('n_seqs','-')} | - | {genes} | {r.get('status')} |" + " |" * 10)
                continue
            n_tips = r.get("n_tips")
            rf = f"{r['rf']}/{r['max_rf']}" if "rf" in r else "-"
            t95 = f"{r['n_true_splits_ge95']}/{r['n_true_splits']}" if "n_true_splits" in r else "-"
            fs = f"{r.get('n_false_splits')} ({r.get('n_false_ge95')})" if "n_false_splits" in r else "-"
            fsnp = t.get("false_snp", 0) + t.get("wrong_alt", 0) if t else None
            print(f"| {runtag} | {sp} | {v} | {n_tips} | {r.get('iq_cols')} | {genes} | {rf} | {fmt(r.get('quartet'),3)} "
                  f"| {fmt(r.get('pat_pearson'),3)} | {fmt(r.get('pat_ratio'))} | {fmt(r.get('term_ratio'))} "
                  f"| {fmt(r.get('int_ratio'))} | {t95} | {fs} | {fmt(miss,1)} | {fmt(fsnp)} | {fmt(rec,2)} |")
