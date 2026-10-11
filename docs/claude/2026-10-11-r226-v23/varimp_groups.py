"""The models' split-gain importance (model_logs/trained_model*.varimp.tsv) summed by feature group (model_features.py),
for several builds side by side.

usage: varimp_groups.py <scripts dir> <local dir> <build> [<build> ...]
"""
import sys
import pandas as pd

sys.path.insert(0, sys.argv[1])
from model_features import DEFAULT_FEATURE_SET, FEATURE_GROUPS  # noqa: E402

L, BUILDS = sys.argv[2], sys.argv[3:]
group_of = {}  # a column in several groups goes to the first of the default set (relative_* to distance)
for g in DEFAULT_FEATURE_SET.split("+") + list(FEATURE_GROUPS):
    cols = FEATURE_GROUPS[g]
    for c in cols:
        group_of.setdefault(c, g)

rows = {}
for b in BUILDS:
    for t in ("pe", "se", "pb", "ont"):
        suf = "" if t == "pe" else "_" + t
        v = pd.read_csv(f"{L}/{b}/protal0.7.9_r226_{b}/model_logs/trained_model{suf}.varimp.tsv", sep="\t")
        v["group"] = v["feature"].map(group_of).fillna("?")
        rows[f"{b} {t}"] = v.groupby("group")["importance"].sum() / v["importance"].sum()
tab = pd.DataFrame(rows).fillna(0).sort_values(f"{BUILDS[-1]} pe", ascending=False)
pd.set_option("display.width", 250)
print("share of the split gain by group")
print((100 * tab).round(2).to_string())
print()
print("the 25 features of most gain in the last build, pe and pb, and their groups")
for t in ("pe", "pb"):
    suf = "" if t == "pe" else "_" + t
    v = pd.read_csv(f"{L}/{BUILDS[-1]}/protal0.7.9_r226_{BUILDS[-1]}/model_logs/trained_model{suf}.varimp.tsv", sep="\t")
    v["share %"] = (100 * v["importance"] / v["importance"].sum()).round(2)
    v["group"] = v["feature"].map(group_of)
    v["rank"] = range(1, len(v) + 1)
    print(t)
    print(v.head(25)[["rank", "feature", "group", "share %"]].to_string(index=False))
