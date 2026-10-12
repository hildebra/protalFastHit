"""The weighted site shift's effect on the species' depth, paired: the small real world's builds with and without the
shift (docs/claude/2026-10-11-r226-v23/weighted_shift_runs.sh; the same samples in both arms), each present species
called in both arms, depth with / depth without - 1, and each arm's depth / true depth - 1
(model_logs/composition_species*.tsv.gz).

usage: shift_depth_pairs.py <off build> <on build> [<off build> <on build> ...]
"""
import sys
import pandas as pd

KEY = ["set", "sample", "species"]
for off, on in zip(sys.argv[1::2], sys.argv[2::2]):
    for suf in ("", "_pb"):
        arms = {}
        for name, d in (("off", off), ("on", on)):
            c = pd.read_csv(f"{d}/model_logs/composition_species{suf}.tsv.gz", sep="\t",
                            usecols=KEY + ["present", "called", "true_depth", "depth"])
            c = c[(c["present"] == 1) & (c["called"] == 1) & c["depth"].notna() & (c["true_depth"] > 0)]
            arms[name] = c.set_index(KEY)
        j = arms["off"][["true_depth", "depth"]].join(arms["on"][["depth"]], rsuffix="_on", how="inner")
        ratio = j["depth_on"] / j["depth"] - 1
        err_off = j["depth"] / j["true_depth"] - 1
        err_on = j["depth_on"] / j["true_depth"] - 1
        print(f"{on} {'pe' if suf == '' else 'pb'}: {len(j)} species called in both arms; depth on / off - 1: median "
              f"{ratio.median():+.4f}, mean {ratio.mean():+.4f}, lower in {(ratio < -1e-9).mean():.1%}, higher in "
              f"{(ratio > 1e-9).mean():.1%}; depth / true - 1, median: off {err_off.median():+.4f}, on {err_on.median():+.4f}")
