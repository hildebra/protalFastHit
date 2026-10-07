#!/usr/bin/env python3
"""Were the FN species' own reads that went to a congener cut off by --align_top? A read's ZN tag counts the taxa its
seeds could not tell apart (an anchor at least 0.8 as long as its longest), also those beyond align_top (3, and the
anchors as long as the third) that were never aligned. ZN above 3 means taxa as good as the tried ones went untried.
Per fate of the FN species' own fragments (sensitivity.own_fates) and for the FP taxa's counted fragments: ZN's
distribution, from the SAMs (the fragment's largest ZN over its records).

    python3 crowding.py --records ~/v15/err_analysis --error-records ~/v15/v15_error_records > crowding.txt
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

import fragments
import sensitivity

ZN = re.compile(rb"\tZN:i:(\d+)")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--error-records", required=True, help="the unpacked v15_error_records (its SAMs)")
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def crowding(error_records, rt):
    """{(<set>:<sample>, qname): the fragment's largest ZN}."""
    out = {}
    for path in sorted(glob.glob(os.path.join(error_records, rt, "*", "*", "*.sam.zst"))):
        which = os.path.basename(os.path.dirname(os.path.dirname(path)))
        sample = which + ":" + os.path.basename(path)[:-len(".sam.zst")]
        text = subprocess.run(["zstd", "-dcq", path], check=True, capture_output=True).stdout
        for line in text.split(b"\n"):
            if not line or line.startswith(b"@"):
                continue
            m = ZN.search(line)
            if m:
                key = (sample, line.split(b"\t", 1)[0].decode())
                n = int(m.group(1))
                if n > out.get(key, -1):
                    out[key] = n
    return out


def summary(zn):
    zn = zn.dropna()
    return {"fragments": len(zn), "with ZN": len(zn), "ZN median": zn.median(), "ZN > 3": round((zn > 3).mean(), 3),
            "ZN > 5": round((zn > 5).mean(), 3), "ZN > 10": round((zn > 10).mean(), 3)}


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        taxa = fragments.load(opts.records, rt, "taxa")
        r = fragments.load(opts.records, rt, "records")[["sample", "qname", "taxon", "zf"]]
        fn = taxa[taxa["error"] == "FN"].copy()
        fn["taxid"] = fn["taxid"].astype(str)
        own = sensitivity.own_fates(f, r, fn)
        own = own[own["fate"] != "nowhere, seeded on other taxa only"]
        zn = crowding(opts.error_records, rt)
        rows = []
        for fate, g in own.groupby("fate"):
            v = pd.Series([zn.get((s, q), np.nan) for s, q in zip(g["sample"], g["qname"])], dtype=float)
            rows.append({"group": "FN own: " + fate, "all": len(g), **summary(v)})
        w = f[f["counted"] & (f["role"] == "FP")]
        v = pd.Series([zn.get((s, q), np.nan) for s, q in zip(w["sample"], w["qname"])], dtype=float)
        rows.append({"group": "FP taxa's counted fragments", "all": len(w), **summary(v)})
        out = pd.DataFrame(rows).drop(columns="fragments")
        print(f"\n# {rt}\n")
        print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
