"""protal's ancestry sites on real GTDB marker genes: how a species' real strains and its real congeners (each left out
in turn, as a novel species the database lacks) fall at the species' derived sites.

For each species T of the chosen genera (select_taxa.py, extract_subset.py), per marker gene: T's representative
copy compared with each congener's representative copy (AncestrySites.h CompareCopies, through ancestry_oracle.py),
and T's sites from those comparisons (Consensus: where >= 0.9 of >= 3 congeners compared carry another base, one
dissenter tolerated from six; with fewer congeners the nearest congener's differences, unpolarised, as protal does
without column weights). Then, against T's representative copy:
  - each strain of T (the species' other genomes): at each site it covers, T's base, the congeners' base, or another;
  - each congener C of T, with T's sites recomputed without C (leave one out: C is the novel species the database
    lacks, and the reads of C are what a false positive on T is made of): the same.
Per (T, query) pair, pooled over genes: the sites covered, the share with T's base ("agreement", what the
ancestry_agreement feature measures), the share with the congeners' base, the identity to T's copies on the compared
columns; and the sites fixed among T's strains (no strain of T other than the query carries another base there).
Per species and gene: the sites per kb, and how many sites a 150-base read covers.

usage: real_ancestry.py <subset dir> <taxa.tsv> <out prefix> [processes]
Writes <out>.pairs.tsv, <out>.species.tsv and <out>.summary.txt (also on stdout).
"""
import bisect
import collections
import glob
import multiprocessing
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
import ancestry_oracle as oracle  # noqa: E402

SUB, TAXA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
PROCS = int(sys.argv[4]) if len(sys.argv) > 4 else 4
READ = 150

species_of, genus_of, role_of = {}, {}, {}
with open(TAXA) as f:
    next(f)
    for line in f:
        acc, role, sp, gen, fam = line.rstrip("\n").split("\t")
        species_of[acc], genus_of[sp], role_of[acc] = sp, gen, role


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts)
    return seqs


# Per gene: the representative copy of each species, and each genome's copy.
genes = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{SUB}/genomic_files_all/*/fna/*.fna"))
rep_copy, genome_copy = {}, {}
for g in genes:
    all_seqs = read_fasta(glob.glob(f"{SUB}/genomic_files_all/*/fna/{g}.fna")[0])
    genome_copy[g] = all_seqs
    rep_copy[g] = {species_of[a]: s for a, s in all_seqs.items() if role_of.get(a) == "representative"}
strains_of = collections.defaultdict(list)
for a, r in role_of.items():
    if r == "strain":
        strains_of[species_of[a]].append(a)
congeners_of = collections.defaultdict(list)
for sp, gen in genus_of.items():
    congeners_of[gen].append(sp)


def classify(sites, own, own_kmers, other):
    """At T's sites, the query's state: (covered, T's base, the congeners' base); the query's identity to T."""
    c = oracle.compare_copies(own, own_kmers, other)
    if c is None:
        return None
    diff = dict(zip(c.positions, c.bases))
    covered = agree = cons = 0
    per_site = []
    for pos, base in zip(sites.positions, sites.bases):
        if pos >= len(c.covered) or not c.covered[pos]:
            per_site.append(None)
            continue
        covered += 1
        q = diff.get(pos)
        if q is None:
            agree += 1
            per_site.append(True)
        else:
            cons += q == base
            per_site.append(False)
    return covered, agree, cons, c.identity, c.compared, per_site


def species_task(sp):
    """All of one species' pairs and per-gene site counts."""
    others = [c for c in congeners_of[genus_of[sp]] if c != sp]
    pairs = collections.defaultdict(lambda: [0, 0, 0, 0, 0.0, 0, 0, 0])  # covered agree cons genes ident_sum fixed_cov fixed_agree compared
    gene_rows = []
    for g in genes:
        own = rep_copy[g].get(sp)
        if not own:
            continue
        kmers = oracle.unique_kmers(own)
        comps = {}
        for c in others:
            if c in rep_copy[g]:
                s = oracle.compare_copies(own, kmers, rep_copy[g][c])
                if s is not None:
                    s.congener = c
                    comps[c] = s
        if not comps:
            continue
        ordered = [comps[c] for c in others if c in comps]
        sites = oracle.consensus(len(own), ordered)
        # Per 150-base window starting every 50 bases: the sites it covers.
        windows = [sites.count(b, b + READ) for b in range(0, max(1, len(own) - READ + 1), 50)]
        gene_rows.append((sp, g, len(own), len(ordered), len(sites.positions), sites.compared,
                          statistics.mean(windows) if windows else 0, sum(w == 0 for w in windows) / max(len(windows), 1)))
        # The strains: their states at T's sites (with all congeners), and the per-site states for the fixed sites.
        strain_states = {}
        for a in strains_of[sp]:
            q = genome_copy[g].get(a)
            if not q:
                continue
            r = classify(sites, own, kmers, q)
            if r is None:
                continue
            strain_states[a] = r
        for a, (cov, ag, co, ident, cmp_, per_site) in strain_states.items():
            p = pairs[(a, "strain")]
            p[0] += cov; p[1] += ag; p[2] += co; p[3] += 1; p[4] += ident; p[7] += cmp_
            # Fixed: no other strain of T carries another base at the site (and some other strain covers it).
            for i, state in enumerate(per_site):
                if state is None:
                    continue
                others_cov = [s[5][i] for b, s in strain_states.items() if b != a and s[5][i] is not None]
                if others_cov and all(others_cov):
                    p[5] += 1
                    p[6] += state
        # The congeners, each left out of T's sites in turn.
        for c in comps:
            rest = [s for s in ordered if s.congener != c]
            if not rest:
                continue
            loo = oracle.consensus(len(own), rest)
            r = classify(loo, own, kmers, rep_copy[g][c])
            if r is None:
                continue
            cov, ag, co, ident, cmp_, per_site = r
            p = pairs[(c, "congener")]
            p[0] += cov; p[1] += ag; p[2] += co; p[3] += 1; p[4] += ident; p[7] += cmp_
            # The fixed sites of T (by its strains) among the left-out sites: how C falls there.
            for i, state in enumerate(per_site):
                if state is None:
                    continue
                st = [s[5] for s in strain_states.values()]
                # strain states were taken at the full-congener sites; match by position
                pos = loo.positions[i]
                k = bisect.bisect_left(sites.positions, pos)
                if k < len(sites.positions) and sites.positions[k] == pos:
                    cov_states = [x[k] for x in st if x[k] is not None]
                    if cov_states and all(cov_states):
                        p[5] += 1
                        p[6] += state
    out = []
    for (q, kind), p in pairs.items():
        out.append((sp, q, kind, p[3], p[0], p[1], p[2], p[4] / p[3], p[5], p[6], p[7]))
    return out, gene_rows


def quantiles(xs, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    if not xs:
        return "n=0"
    xs = sorted(xs)
    return f"n={len(xs)}: " + " ".join(f"{int(q * 100)}%={xs[min(len(xs) - 1, int(q * len(xs)))]:.3f}" for q in qs) \
        + f", mean {statistics.mean(xs):.3f}"


def auc(pos, neg):
    if not pos or not neg:
        return float("nan")
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank, i, rp = 1, 0, 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (rank + rank + (j - i) - 1) / 2
        rp += r * sum(1 for k in range(i, j) if allv[k][1])
        rank += j - i
        i = j
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


if __name__ == "__main__":
    species = sorted(genus_of)
    with multiprocessing.Pool(PROCS) as pool:
        results = pool.map(species_task, species, chunksize=1)
    pairs = [r for res, _ in results for r in res]
    gene_rows = [r for _, res in results for r in res]
    with open(f"{OUT}.pairs.tsv", "w") as o:
        o.write("species\tquery\tkind\tgenes\tsites_covered\tagree\tcongener_base\tidentity\tfixed_covered\tfixed_agree\tcompared\n")
        for r in pairs:
            o.write("\t".join(str(round(x, 5)) if isinstance(x, float) else str(x) for x in r) + "\n")
    with open(f"{OUT}.species.tsv", "w") as o:
        o.write("species\tgene\tlength\tcongeners\tsites\tcompared\tsites_per_read\tread_share_without_site\n")
        for r in gene_rows:
            o.write("\t".join(str(round(x, 5)) if isinstance(x, float) else str(x) for x in r) + "\n")

    lines = []
    say = lines.append
    say(f"species {len(species)}, genera {len(congeners_of)}, strains {sum(len(v) for v in strains_of.values())}, genes {len(genes)}")
    by_sp = collections.defaultdict(list)
    for r in gene_rows:
        by_sp[r[0]].append(r)
    per_kb = [1000 * sum(r[4] for r in rs) / max(1, sum(r[5] for r in rs)) for rs in by_sp.values()]
    say("sites per kb compared, per species: " + quantiles(per_kb))
    say("sites per 150-base read, per species (mean over its genes' windows): "
        + quantiles([statistics.mean(r[6] for r in rs) for rs in by_sp.values()]))
    say("share of 150-base reads without a site, per species: "
        + quantiles([statistics.mean(r[7] for r in rs) for rs in by_sp.values()]))
    for kind in ("strain", "congener"):
        rows = [r for r in pairs if r[2] == kind and r[4] >= 10]
        say(f"-- {kind}s (pairs with >= 10 sites covered; congeners left out of the sites)")
        say("   identity to the representative: " + quantiles([r[7] for r in rows]))
        say("   sites covered: " + quantiles([r[4] for r in rows]))
        say("   agreement (T's base): " + quantiles([r[5] / r[4] for r in rows]))
        say("   the congeners' base: " + quantiles([r[6] / r[4] for r in rows]))
        fx = [r for r in rows if r[8] >= 5]
        say(f"   agreement at the fixed sites (>= 5 covered, {len(fx)} pairs): " + quantiles([r[9] / r[8] for r in fx]))
    s = [r[5] / r[4] for r in pairs if r[2] == "strain" and r[4] >= 10]
    c = [r[5] / r[4] for r in pairs if r[2] == "congener" and r[4] >= 10]
    say(f"AUC of agreement, strains over left-out congeners: {auc(s, c):.3f}")
    sf = [r[9] / r[8] for r in pairs if r[2] == "strain" and r[8] >= 5]
    cf = [r[9] / r[8] for r in pairs if r[2] == "congener" and r[8] >= 5]
    say(f"AUC of the fixed-site agreement: {auc(sf, cf):.3f}")
    # By identity band: the hard cases are congeners as close as strains.
    say("-- by identity to T's copies: agreement of strains | of left-out congeners (pairs)")
    for lo, hi in ((0.0, 0.85), (0.85, 0.9), (0.9, 0.95), (0.95, 0.97), (0.97, 0.99), (0.99, 1.01)):
        a = [r[5] / r[4] for r in pairs if r[2] == "strain" and r[4] >= 10 and lo <= r[7] < hi]
        b = [r[5] / r[4] for r in pairs if r[2] == "congener" and r[4] >= 10 and lo <= r[7] < hi]
        med = lambda xs: f"{statistics.median(xs):.3f}" if xs else "-"
        say(f"   [{lo:.2f}, {hi:.2f}): {med(a)} ({len(a)}) | {med(b)} ({len(b)})")
    text = "\n".join(lines)
    print(text)
    with open(f"{OUT}.summary.txt", "w") as o:
        o.write(text + "\n")
