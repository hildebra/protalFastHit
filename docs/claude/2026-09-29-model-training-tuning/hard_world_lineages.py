#!/usr/bin/env python3
"""Write the lineages of a harder mini world: 16 genera of 4 species, two genera per family,
12 bacterial and 4 archaeal genera. The 4th species of every genus is left out of the database
(heldout.txt). Genera 1-8 (6 bacterial, 2 archaeal) feed the training samples; the test samples
draw from all 16, so half their species were never seen in training."""
import sys

out = sys.argv[1]
bac = ["Proxima", "Vicinia", "Cognata", "Affinia", "Similia", "Paria", "Gemella", "Fratria",
       "Sororia", "Consobrina", "Propinqua", "Germana"]
arc = ["Thermoproxa", "Thermovicina", "Halocognata", "Haloaffina"]
species = ["alpha", "beta", "gamma", "novus"]
lineages, heldout, train_genera = [], [], []
for i, genus in enumerate(bac + arc):
    domain = "Archaea" if genus in arc else "Bacteria"
    fam = i // 2  # two genera per family, families in pairs per order
    order = fam // 2
    lin = f"d__{domain};p__Phyl{order % 3}ota;c__Class{order}ia;o__Ord{order}ales;f__Fam{fam}aceae;g__{genus}"
    for sp in species:
        lineages.append(f"{lin};s__{genus} {sp}")
    heldout.append(f"s__{genus} novus")
    if (domain == "Bacteria" and i < 6) or genus in arc[:2]:
        train_genera.append(genus)
open(f"{out}/lineages.txt", "w").write("\n".join(lineages) + "\n")
open(f"{out}/heldout.txt", "w").write("\n".join(heldout) + "\n")
open(f"{out}/train_genera.txt", "w").write("\n".join(train_genera) + "\n")
print(len(lineages), "species;", len(heldout), "held out; training genera:", ", ".join(train_genera))
