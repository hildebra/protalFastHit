"""The (species, gene) copies a sample's kept records touch (the sum of present_genes over its taxa: a lower bound on
the ancestry and column-weight cache entries its profiling builds, as records the filters drop touch copies too), per
sample, by read type and kind of sample, from collect_training_data.py tables.

usage: copies_per_sample.py <training table> [<training table> ...]
"""
import sys
import pandas as pd

for path in sys.argv[1:]:
    t = pd.read_csv(path, sep="\t", usecols=["meta_sample", "present_genes", "fragments"])
    kind = t["meta_sample"].str.startswith("sc_").map({True: "scenario", False: "design"})
    g = t.groupby([kind, "meta_sample"]).agg(taxa=("present_genes", "size"), copies=("present_genes", "sum"),
                                             fragments=("fragments", "sum"))
    print(path)
    print(g.groupby(level=0).describe(percentiles=[0.5, 0.9]).T.round(0).to_string())
    print()
