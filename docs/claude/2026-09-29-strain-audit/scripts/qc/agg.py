#!/usr/bin/env python3
"""Aggregate ablate.json of several runs: what each default filter removes (default vs that filter off)."""
import sys, json, os
from collections import defaultdict

OFF = {"gene-min-samples 3": "gene-min-samples 0", "gene-min-hcov 0.3": "gene-min-hcov 0",
       "gene-min-mean-depth 1.0": "gene-min-mean-depth 0", "reapply-hcov 1000": "reapply-hcov 0",
       "MRate2 (Tukey+cells)": "mrate2 off", "min-parsimony-samples 2": "min-parsimony 0"}
for path in sys.argv[1:]:
    res = json.load(open(path))
    by = defaultdict(dict)
    for r in res:
        by[r["config"]][r["species"]] = r
    run = os.path.basename(os.path.dirname(path))
    d = by["default"]
    n_sp = len(d)
    tot = lambda cfg, k: sum(by[cfg][s][k] for s in by[cfg])
    print(f"== {run}: {n_sp} species, sample rows {tot('default','samples_in')}, genes {tot('default','genes_in')}, "
          f"var sites {tot('default','var_in')}; default keeps: species with MSA "
          f"{sum(1 for s in d if d[s]['cols_out'])}, samples {tot('default','samples_out')}, genes {tot('default','genes_out')}, "
          f"var {tot('default','var_out')}; all-off keeps var {tot('all off','var_out')}")
    for f, cfg in OFF.items():
        o = by[cfg]
        sp_lost = sum(1 for s in d if not d[s]["cols_out"] and o[s]["cols_out"])
        dsam = tot(cfg, "samples_out") - tot("default", "samples_out")
        dgen = tot(cfg, "genes_out") - tot("default", "genes_out")
        dvar = tot(cfg, "var_out") - tot("default", "var_out")
        dk = [(d[s]["dist_kept"], o[s]["dist_kept"]) for s in d if d[s]["dist_kept"] is not None and o[s]["dist_kept"] is not None]
        dks = f"{sum(a for a, _ in dk)/len(dk):.2f} vs {sum(b for _, b in dk)/len(dk):.2f}" if dk else "-"
        print(f"   {f:26} species lost {sp_lost}, samples {dsam}, genes {dgen}, variable sites {dvar}; "
              f"mean pairwise distance kept default vs off: {dks}")
