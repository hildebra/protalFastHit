#!/usr/bin/env python3
"""congener_lineages.py [SEED] - lineages for simulate_gtdb_release.py --lineages with many congeners: 12 genera
of 10 species each (in 4 families of 3 genera) and 40 species in 25 small genera of 1-3 species, all bacteria,
every name unique. Writes one GTDB lineage per line to stdout (160 species).
"""
import random
import sys

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
SYLLABLES = ["ba", "co", "li", "ter", "mo", "na", "ru", "vi", "do", "xa", "pe", "sal", "fi", "gor", "lu", "ne",
             "thi", "ca", "ro", "mi", "zu", "bel", "cor", "dam", "fla", "gen", "hal", "ix", "jo", "ka", "lem"]
used = set()


def stem():
    while True:
        s = "".join(rng.choice(SYLLABLES) for _ in range(rng.choice((2, 3))))
        if s not in used:
            used.add(s)
            return s


def name(ending):
    return stem().capitalize() + ending


phylum, klass, order = "p__" + name("ota"), "c__" + name("ia"), "o__" + name("ales")
families = ["f__" + name("aceae") for _ in range(8)]
sizes = [10] * 12
small = []
while sum(small) < 40:
    small.append(min(rng.choice((1, 1, 2, 3)), 40 - sum(small)))
lines = []
for i, size in enumerate(sizes + small):
    family = families[i // 3] if i < 12 else families[4 + i % 4]
    genus = name(rng.choice(["bacter", "coccus", "monas", "ella", "spira", "vibrio"]))
    epithets = set()
    while len(epithets) < size:
        epithets.add(stem() + rng.choice(["us", "is", "a", "um", "ensis", "ii"]))
    for epithet in sorted(epithets):
        lines.append(";".join(["d__Bacteria", phylum, klass, order, family, f"g__{genus}", f"s__{genus} {epithet}"]))
sys.stdout.write("".join(line + "\n" for line in lines))
