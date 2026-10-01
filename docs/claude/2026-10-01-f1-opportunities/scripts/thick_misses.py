#!/usr/bin/env python3
"""thick_misses.py [READ_TYPE] - the present taxa of 0.7.1's test set (refitted model, seed 1) missed although they have more
than 100 fragments, beside the medians of the present taxa called with as many: which features set them apart.
Writes results/thick_misses_<rt>.md."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["V2"] = os.path.expanduser(os.environ.get("V071", "~/bench071/V071"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-alignment-features"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
import exp_lib as X  # noqa: E402
from model_features import NORMALIZED_FEATURES as BASE  # noqa: E402

rt = sys.argv[1] if len(sys.argv) > 1 else "pe"
train, test = X.load("training", rt), X.load("test", rt)
test["p"] = X.test_predict(train, test, BASE, 1)
thick = test[(test["truth"] == 1) & (test["fragments"] > 100)]
missed, called = thick[thick["p"] < 0.5], thick[thick["p"] >= 0.5]
cols = ["p", "fragments", "depth", "identity", "top_identity", "low_identity_share", "hit_gene_fraction", "gene_presence_ratio",
        "depth_cv", "uniqueness", "lu_per_kb", "lsu_per_kb", "variant_sites_per_kb", "multiallelic_sites_per_kb", "RAF1",
        "low_mapq_share", "congener_fit_share", "other_genus_fit_share", "meta_rep_genome", "meta_novel_congener"]
cols = [c for c in cols if c in test]
lines = [f"# {rt}: {len(missed)} present taxa with more than 100 fragments missed, {len(called)} called", "",
         "| " + " | ".join(["taxon"] + cols) + " |", "|" + "---|" * (len(cols) + 1)]
for _, r in missed.sort_values("p").iterrows():
    lines.append("| " + " | ".join([r["taxon_name"][3:]] + [f"{r[c]:.3g}" for c in cols]) + " |")
lines.append("| median of those called | " + " | ".join(f"{called[c].median():.3g}" for c in cols) + " |")
out = os.path.join(HERE, "..", "results", f"thick_misses_{rt}.md")
open(out, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
