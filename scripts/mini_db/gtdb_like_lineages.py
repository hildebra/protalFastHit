#!/usr/bin/env python3
"""Write lineages shaped like GTDB's for simulate_gtdb_release.py --lineages.

A world for tuning the presence model's training: many species, most genera with one species and a
few with many (sizes drawn from a truncated power law, as in GTDB), genera grouped into families,
orders, classes and phyla the same way, and a share of archaea. Every name is unique, across ranks
and domains, as in GTDB.

  python3 scripts/mini_db/gtdb_like_lineages.py --species 900 --archaea 0.08 --seed 1 > lineages.txt
  python3 scripts/mini_db/simulate_gtdb_release.py --outdir world --lineages lineages.txt \\
      --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06

simulate_gtdb_release.py numbers species with three digits, so at most 999.
"""

import argparse
import random
import sys

SYLLABLES = ["ba", "co", "li", "ter", "mo", "na", "ru", "vi", "do", "xa", "pe", "sal", "fi", "gor", "lu", "ne",
             "thi", "ca", "ro", "mi", "zu", "bel", "cor", "dam", "fla", "gen", "hal", "ix", "jo", "ka", "lem",
             "nor", "pra", "qui", "ser", "tul", "ur", "ven", "wo", "yal"]
ENDINGS = {"g": ["bacter", "coccus", "monas", "ella", "spira", "vibrio", "plasma", "ococcus", "illus", "ia"],
           "f": ["aceae"], "o": ["ales"], "c": ["ia"], "p": ["ota"]}


class Names:
    def __init__(self, rng):
        self.rng, self.used = rng, set()

    def stem(self, parts=(2, 3)):
        while True:
            s = "".join(self.rng.choice(SYLLABLES) for _ in range(self.rng.choice(parts)))
            if s not in self.used:
                self.used.add(s)
                return s

    def taxon(self, rank):
        return self.stem().capitalize() + self.rng.choice(ENDINGS[rank])


def power_law(rng, alpha, cap):
    """An integer from 1..cap with P(k) proportional to k^-alpha."""
    weights = [k ** -alpha for k in range(1, cap + 1)]
    return rng.choices(range(1, cap + 1), weights)[0]


def group(rng, items, alpha, cap):
    """items split into consecutive groups of power-law sizes."""
    out, i = [], 0
    while i < len(items):
        size = power_law(rng, alpha, cap)
        out.append(items[i:i + size])
        i += size
    return out


def lineages(n_species, archaea, seed):
    rng = random.Random(seed)
    names = Names(rng)
    out = []
    n_arch = max(1, round(archaea * n_species)) if archaea > 0 else 0
    for domain, n in (("Bacteria", n_species - n_arch), ("Archaea", n_arch)):
        if n == 0:
            continue
        genera, left = [], n
        while left > 0:  # species per genus: most 1, a few up to 40
            size = min(left, power_law(rng, 1.9, 40))
            genera.append(size)
            left -= size
        # genera -> families -> orders -> classes -> phyla, each in power-law groups
        levels = [list(range(len(genera)))]
        for rank, alpha, cap in (("f", 1.8, 15), ("o", 1.8, 10), ("c", 1.8, 8), ("p", 1.8, 6)):
            levels.append(group(rng, list(range(len(levels[-1]))), alpha, cap))
        # name every node top-down
        path_of = {}
        for p_index, phylum in enumerate(levels[4]):
            p = "p__" + names.taxon("p")
            for c_index in phylum:
                c = "c__" + names.taxon("c")
                for o_index in levels[3][c_index]:
                    o = "o__" + names.taxon("o")
                    for f_index in levels[2][o_index]:
                        f = "f__" + names.taxon("f")
                        for g_index in levels[1][f_index]:
                            path_of[g_index] = (p, c, o, f)
        for g_index, size in enumerate(genera):
            genus = names.taxon("g")
            epithets = set()
            while len(epithets) < size:
                epithets.add(names.stem((2, 3)) + rng.choice(["us", "is", "a", "um", "ensis", "ii"]))
            for epithet in sorted(epithets):
                out.append(";".join([f"d__{domain}", *path_of[g_index], f"g__{genus}", f"s__{genus} {epithet}"]))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--species", type=int, default=900)
    p.add_argument("--archaea", type=float, default=0.08, help="share of archaeal species (default 0.08)")
    p.add_argument("--seed", type=int, default=1)
    opts = p.parse_args(argv)
    if not 1 <= opts.species <= 999:
        sys.exit("--species must be 1..999 (simulate_gtdb_release.py numbers species with three digits)")
    lines = lineages(opts.species, opts.archaea, opts.seed)
    sys.stdout.write("".join(line + "\n" for line in lines))


if __name__ == "__main__":
    main()
