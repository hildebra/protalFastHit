#!/usr/bin/env python3
"""congener_diag.py - where the gene-median rule (gm006) loses to top - 0.08 (top008) against db_missing at
2x150: the species with the largest extra Bray-Curtis, their held-out congeners in the sample, and the share
of the extra error on species with such a congener. Reads ~/stress as score_stress.py does."""
import collections, csv, glob, os, sys
sys.argv = [sys.argv[0], os.path.expanduser("~/stress")]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_stress as S
rows_by = S.design
contrib = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0, ""])
for sample in [s for s in sorted(rows_by) if s.startswith("s150")]:
    rows = rows_by[sample]
    depth = collections.Counter()
    for r in rows: depth[r["species"]] += float(r["depth"])
    present = set(depth) - S.heldout
    held_in = {h: depth[h] for h in set(depth) & S.heldout}
    true = {s: depth[s] / sum(depth[x] for x in present) for s in present}
    est = {}
    for rule in ("top008", "gm006", "none"):
        got = S.reported(glob.glob(os.path.join(S.B, "prof", rule, "missing", sample, "*.profile"))[0])
        shown = sum(got.get(s, 0) for s in present)
        est[rule] = {s: got.get(s, 0) / shown for s in present}
    for s in present:
        cong = [(h, d) for h, d in held_in.items() if S.genus_of[h] == S.genus_of[s]]
        a, b = abs(true[s] - est["top008"][s]) / 2, abs(true[s] - est["gm006"][s]) / 2
        key = (sample, s)
        contrib[key] = [true[s], est["top008"][s], est["gm006"][s], est["none"][s], depth[s], cong, b - a]
top = sorted(contrib.items(), key=lambda kv: -kv[1][6])[:10]
print("largest extra Bray-Curtis of gm006 over top008 (db_missing, 2x150)")
for (sample, s), (t, a, b, n, d, cong, diff) in top:
    print(f"{sample} {s:28s} true {t:.4f} top008 {a:.4f} gm006 {b:.4f} none {n:.4f} own depth {d:5.1f}x  held-out congeners in sample: "
          + ", ".join(f"{h} {x:.1f}x" for h, x in cong) + f"  (+{diff:.4f})")
tot = sum(v[6] for v in contrib.values())
withc = sum(v[6] for v in contrib.values() if v[5])
print(f"total extra Bray-Curtis summed over the 6 samples {tot:.4f}, of it on species with a held-out congener in the sample {withc:.4f}")
