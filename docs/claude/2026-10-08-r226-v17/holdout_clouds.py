#!/usr/bin/env python3
"""Does the training database's hold-out respect the clouds of near relatives around each species? From the build's
taxonomy and held-out list: for the species held out alone, how many congeners the training database keeps (the
places their reads can land as false positives); for the kept species, how many congeners are held out (whose reads
they take); and, where a build's tables give the nearest kept congener's distance (db_nearest_congener) and the FP
rows' identity, how near the held-out sisters are: a held-out species within a strain's distance of a kept one
teaches the model that a cloud of reads at that identity is absent.

    python3 holdout_clouds.py --build local/v17 > holdout_clouds.txt
"""
import argparse
import collections
import csv
import os
import sys

import numpy as np


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True, help="the build's local folder: internal_taxonomy.dmp, heldout_species.txt")
    return p.parse_args(argv)


def read_taxonomy(path):
    parent, rank, name = {}, {}, {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            parent[row[0]], rank[row[0]], name[row[0]] = row[1], row[4], row[3]
    return parent, rank, name


def ancestor(parent, rank, taxid, want):
    for _ in range(64):
        if rank.get(taxid) == want:
            return taxid
        if parent.get(taxid, taxid) == taxid:
            return None
        taxid = parent[taxid]
    return None


def main(argv=None):
    opts = parse_args(argv)
    parent, rank, name = read_taxonomy(os.path.join(opts.build, "internal_taxonomy.dmp"))
    held = {}
    with open(os.path.join(opts.build, "heldout_species.txt")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].strip():
                held[f[0].strip()] = f[1] if len(f) > 1 and f[1] else "species"
    species = {t for t, r in rank.items() if r == "species"}
    genus_of = {t: ancestor(parent, rank, t, "genus") for t in species}
    family_of = {t: ancestor(parent, rank, t, "family") for t in species}
    kept = {t for t in species if name[t] not in held}
    held_ids = {t: held[name[t]] for t in species if name[t] in held}
    by_genus = collections.defaultdict(lambda: [0, 0])  # kept, held out
    by_family = collections.defaultdict(lambda: [0, 0])
    for t in species:
        by_genus[genus_of[t]][0 if t in kept else 1] += 1
        by_family[family_of[t]][0 if t in kept else 1] += 1
    print(f"{len(species)} species in the taxonomy, {len(kept)} kept, {len(held_ids)} held out "
          f"({collections.Counter(held_ids.values()).most_common()})")
    print(f"genera: {len(by_genus)}, with 2 or more species {sum(1 for k, h in by_genus.values() if k + h >= 2)}; "
          f"species whose genus has another species: {sum(1 for t in species if sum(by_genus[genus_of[t]]) >= 2) / len(species):.3f}")

    # Held out alone: the kept congeners (the places its reads can land as a false positive at genus level)
    alone = [t for t, r in held_ids.items() if r == "species"]
    kept_cong = np.array([by_genus[genus_of[t]][0] for t in alone])
    kept_fam = np.array([by_family[family_of[t]][0] - by_genus[genus_of[t]][0] for t in alone])
    print(f"\nheld out alone ({len(alone)}): kept congeners none {np.mean(kept_cong == 0):.3f}, 1 {np.mean(kept_cong == 1):.3f}, "
          f"2-4 {np.mean((kept_cong >= 2) & (kept_cong <= 4)):.3f}, 5 or more {np.mean(kept_cong >= 5):.3f}; of those without, "
          f"kept species of other genera in the family: none {np.mean(kept_fam[kept_cong == 0] == 0):.3f}")
    # Held out in a clade: by construction no kept species at the clade's rank
    for r in ("genus", "family", "order", "class", "phylum"):
        members = [t for t, rr in held_ids.items() if rr == r]
        if members:
            kc = np.array([by_genus[genus_of[t]][0] for t in members])
            print(f"held out as a {r} clade ({len(members)}): kept congeners none {np.mean(kc == 0):.3f} (the clade's "
                  f"congeners are held out with them{'' if r == 'genus' else ', and their genus too'})")

    # Kept: the held-out congeners (whose reads they take)
    held_cong = np.array([by_genus[genus_of[t]][1] for t in kept])
    print(f"\nkept ({len(kept)}): held-out congeners none {np.mean(held_cong == 0):.3f}, 1 {np.mean(held_cong == 1):.3f}, "
          f"2 or more {np.mean(held_cong >= 2):.3f}; kept species whose genus is otherwise all held out: "
          f"{sum(1 for t in kept if by_genus[genus_of[t]][0] == 1 and by_genus[genus_of[t]][1] >= 1)}")
    sizes = np.array([sum(by_genus[genus_of[t]]) for t in kept])
    for lo, hi in ((2, 2), (3, 5), (6, 20), (21, 10 ** 6)):
        sel = (sizes >= lo) & (sizes <= hi)
        if sel.any():
            print(f"  kept species in genera of {lo}-{hi if hi < 10 ** 6 else 'more'} species ({int(sel.sum())}): "
                  f"with a held-out congener {np.mean(held_cong[sel] >= 1):.3f}")

    # The build's tables, if here: the FP rows' identity by class, and the nearest kept congener's distance
    tables = [os.path.join(opts.build, w, "training_data.tsv") for w in ("training", "test")]
    tables = [p for p in tables if os.path.isfile(p)]
    if not tables:
        return 0
    try:
        import pandas as pd
    except ImportError:
        return 0
    cols = ["meta_sample", "taxon", "truth", "meta_novel_level", "meta_relative_rank", "meta_rep_genome", "identity",
            "fragments", "db_nearest_congener", "db_congeners_02", "db_congeners_05"]
    t = pd.concat([pd.read_csv(p, sep="\t", usecols=lambda c: c in cols) for p in tables], ignore_index=True)
    if "db_nearest_congener" not in t:
        return 0
    novel = (t["truth"] == 0) & (t["meta_novel_level"].astype(str) == "species") & (t["meta_relative_rank"].astype(str) == "genus")
    present = t["truth"] == 1
    print(f"\n{len(t)} rows of the build's tables: absent beside a held-out congener {int(novel.sum())}, present {int(present.sum())}")
    for label, sel in (("absent beside a held-out congener", novel & (t["fragments"] >= 3)),
                       ("present, another genome than the representative", present & (t["meta_rep_genome"] == 0) & (t["fragments"] >= 3)),
                       ("present, the representative", present & (t["meta_rep_genome"] == 1) & (t["fragments"] >= 3))):
        x = t.loc[sel, "identity"]
        q = x.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(4).tolist()
        print(f"  {label} (3+ fragments, {len(x)}): identity 10/25/50/75/90% {q}; at or above 0.985: {np.mean(x >= 0.985):.3f}, "
              f"0.97-0.985: {np.mean((x >= 0.97) & (x < 0.985)):.3f}")
    d = t.loc[t["db_nearest_congener"] >= 0, "db_nearest_congener"]
    print(f"  db_nearest_congener (the kept congener nearest the taxon's reference, by the training database): "
          f"10/25/50/75/90% {d.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(4).tolist()}; within 0.02: {np.mean(d <= 0.02):.3f}, "
          f"within 0.05: {np.mean(d <= 0.05):.3f}; -1 (none within 0.15 or no table): {np.mean(t['db_nearest_congener'] < 0):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
