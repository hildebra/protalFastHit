"""How many of a build's false negatives had reads of their own on the taxon: own_best_on_taxon of
model_logs/error_reads/<type>/taxa.tsv.gz (the taxon's own fragments whose best alignment is on it, any MAPQ), by
read type and set, and the F1 if the false negatives without 3 such reads were not counted as present (an evaluation
that requires the reads to show the species).

usage: own_reads.py <build dir> <read type> [<read type> ...]
"""
import sys
import pandas as pd

D, types = sys.argv[1], sys.argv[2:]
for T in types:
    SUF = "" if T == "pe" else "_" + T
    calls = pd.read_csv(f"{D}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t", usecols=["truth", "set", "call"])
    err = pd.read_csv(f"{D}/model_logs/error_reads/{T}/taxa.tsv.gz", sep="\t")
    for s in ("training", "test"):
        x = calls[calls["set"] == s]
        tp = int(((x["truth"] == 1) & (x["call"] == 1)).sum())
        fp = int(((x["truth"] == 0) & (x["call"] == 1)).sum())
        fn = int(((x["truth"] == 1) & (x["call"] == 0)).sum())
        f = err[(err["set"] == s) & (err["error"] == "FN")]
        own = pd.to_numeric(f["own_best_on_taxon"], errors="coerce").fillna(0)
        n0, n12 = int((own == 0).sum()), int(((own >= 1) & (own <= 2)).sum())
        n3 = len(f) - n0 - n12
        f1 = 2 * tp / (2 * tp + fp + fn)
        f1_own = 2 * tp / (2 * tp + fp + fn - n0 - n12)
        print(f"{T} {s:8s}: FN {fn}: no own read on the taxon {n0} ({n0 / len(f):.0%}), 1-2 {n12} ({n12 / len(f):.0%}), "
              f"3+ {n3} ({n3 / len(f):.0%}); F1 {f1:.4f}, counting as present only the missed taxa with 3+ own reads "
              f"{f1_own:.4f} ({f1_own - f1:+.4f})")
