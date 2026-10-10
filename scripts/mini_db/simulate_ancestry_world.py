#!/usr/bin/env python3
"""simulate_ancestry_world.py - a synthetic GTDB-like release grown on explicit trees, with its samples, for the
true-positive test of protal's ancestry sites, strain alleles and polymorphic sites
(scripts/ancestry_truth_test.py; docs/claude/2026-10-09-ancestry-true-positive-test/README.md).

simulate_gtdb_release.py evolves every species of a genus from the genus node and every genome of a species from the
species node: stars, in which a novel congener disagrees at every ancestry site and two strains share nothing but the
representative's own mutations. Here each test genus is a tree with known branch lengths (substitutions per site):

    genus root (height H)
     ├── T (h_T): S's split from its nearest database relative, the "sister" (a twin S' in twin genera)
     │    ├── sister (database)
     │    └── S's stem, length L = h_T - c_S; the novel species N_b leave it at fraction b of it
     │         └── S's crown (c_S): a coalescent genealogy of --tips genomes: the representative, the allele genomes
     │              (in the full reference only) and the query strains (never in the release)
     └── the other database congeners (a coalescent of their own)

The marker genes (substitutions only, codon-aware, no stop codon) evolve down that tree, so every site's history is
known. A query strain Q of S carries the congeners' base at the sites of the representative's own lineage below
MRCA(Q, rep) and S's base at S's stem; a novel species N_b carries S's base at b of the stem and the ancestral base
elsewhere. Roles in S's genealogy (simulation/roles.tsv):
  rep        the representative (database)
  allele     an allele genome (--allele_genomes per genus): the release's other genome of S, so in the full reference
             and strain_alleles.tsv, never simulated except as Q_leak
  Q_near     a query strain whose sister tip is an allele genome (recently diverged from a known strain)
  Q_lone     a query strain whose sister tip is no allele genome (recently diverged from an unknown one)
  Q_deep     a query strain of the crown's other side, in a clade without allele genomes (a long-diverged lineage)
  Q_rep      the query strain nearest the representative
  Q_leak     (in the samples) an allele genome itself: the leaky positive control
  decoy      a tip of Q_deep's clade, drawn as they are, labelled a species of its own (class D): nothing but its
             label tells it from Q_deep, so any feature that separates the two is a leak (--no_decoy: none)
  pop        the rest of the genealogy: donors of within-species recombination, never simulated
Novel species (never in the release): N_<b> for each --novel_b; N_1 leaves at S's crown, a lineage as old as S's
diversity but outside it (a deep strain shares its side of the crown with S's other genomes, N_1 with none). Derived
genomes:
  Q_imp<p>   a Q_near (else Q_lone, Q_deep) with segments of the sister's genes in a share p of its genes
             (--import_shares): congeneric recombination into a strain
  N_imp<p>   N_0.5 with segments of S's genomes in a share p of its genes: congeneric recombination into a novel species
  Q_ils      a Q_deep with the ancestral base at --ils_share of S's stem sites (retained ancestral polymorphism)
  P_f<f>     planted, in genera with >= 3 congeners and 1-4 allele genomes, no twin: S's representative carrying a share f
             of the farthest allele genome's differences, filled with new ones to --planted_divergence
  A_a<a>     planted likewise: the congeners' base at a share 1 - a of the consensus ancestry sites
Recombination regimes (--regimes, per genus): within-species imports (each S genome, per gene, from another S
genome) and congeneric ones (every genome of the genus, per gene, from another species), segments of --segment bases
on average.

Factors, crossed over the test genera (--replicates of each cell): database congeners (--congeners; 2 leaves the
ancestry sites to the nearest congener's differences, protal's kMinCongeners is 3), a twin (--twins: the sister 0.02-
0.03 away instead of 0.05-0.08), allele genomes (--allele_genomes) and the regime (--recombination). Each genus also
draws its species' width (crown height c_S) and stem length, so lineages diverged recently and long ago are both
there. Background genera (--background_genera) of a few species fill the samples.

Writes into --outdir a release in simulate_gtdb_release.py's layout (taxonomy, metadata, marker genes of the
representatives and of all genomes, genome FASTAs, simulation/genomes.tsv) holding the database's genomes only, and:
  simulation/genera.tsv        each test genus: its factors, heights, target species and roles
  simulation/roles.tsv         every genome of a test genus: role, genus, b/f/a, MRCA heights, and its differences
                               from the target's representative by class (truth for the analyses)
  simulation/all_genomes.tsv   every genome (release and query): accession, gtdb_taxonomy, fasta_path, genome_length,
                               role, genus, in_release
  simulation/truth_sites.tsv.gz  per test genus and marker: S's stem sites, the representative's lineage, protal's
                               ancestry sites (its rule on the database's copies, by ancestry_oracle.py's port; the
                               planted A_a genomes and the cons_*/fixed_* counts of roles.tsv use them), each allele
                               genome's differences
  simulation/imports.tsv       every imported segment (genome, marker, start, end, donor)
  samples/design.tsv           every sample's focus genome per test genus (set, sample, genus, genome, role, class,
                               depth) and its background genomes
  samples/pe/manifest.tsv      paired-end samples for simulate_metagenomes --from_manifest (exact read pairs per genome)
  samples/hifi/long_samples.tsv, long_genomes.tsv   HiFi samples for simulate_metagenomes --long_samples
  samples/exact/samples.tsv, communities/   error-free paired-end samples for simulate_reads.py --error_rate 0 (the
                               planted genomes, Q_leak and the oracle's)
  samples/<set>/truth/<sample>.tsv   protal --profile_truth files: the species present (a Q's target, the background;
                               a novel species' own lineage, which the database lacks)

  python3 scripts/mini_db/simulate_ancestry_world.py --outdir world --seed 1
  python3 scripts/mini_db/simulate_ancestry_world.py --outdir small --quick        (two genera, a few samples)
"""

import argparse
import csv
import gzip
import itertools
import math
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import ancestry_oracle as oracle  # noqa: E402
import simulate_gtdb_release as sgr  # noqa: E402
from gtdb_like_lineages import Names  # noqa: E402

STOP = sgr.STOP
BASES = "ACGT"
MIN_CONGENERS = oracle.MIN_CONGENERS  # protal's ancestry::kMinCongeners
REGIMES = "none:0:0,low:0.1:0.03,high:0.3:0.15"
Q_CLASS = ("Q_near", "Q_lone", "Q_deep", "Q_rep", "Q_imp", "Q_ils")
N_CLASS = ("N_", "N_imp")


def role_class(role):
    """'Q' (the target is present), 'N' (absent: a novel species' reads land on it), 'P' (planted or leaky: present,
    analysed apart), 'D' (the decoy: absent by its label, a deep strain by its history) or '' (in the release, or
    never simulated)."""
    if role.startswith(("P_", "A_")) or role == "Q_leak":
        return "P"
    if role == "decoy":
        return "D"
    if role.startswith("Q_"):
        return "Q"
    if role.startswith("N_"):
        return "N"
    return ""


def floats(text):
    return [float(x) for x in str(text).split(",") if x.strip()]


def ints(text):
    return [int(x) for x in str(text).split(",") if x.strip()]


def label(x):
    """0.5 -> '0.5', 1.0 -> '1', 0.95 -> '0.95': b, f and a in role names."""
    return f"{x:g}"


# ---- sequences ----------------------------------------------------------------------------------------------------

def substitute(rng, seq, d, rates=None):
    """Bernoulli(d) substitutions per site, the first and last codon kept, none making a stop codon:
    simulate_gtdb_release.mutate without its codon indels. With `rates` (a multiplier per site, --site_rates: the
    same at every branch of the world, so that a column keeps its rate class down the whole tree) site i changes with
    probability d x rates[i], at most 1."""
    if d <= 0:
        return seq
    s = list(seq)
    if rates is None:
        log_q = math.log(1 - d)
        pos = 3 + int(math.log(1 - rng.random()) / log_q)
        while pos < len(s) - 3:
            c = pos - pos % 3
            alts = [b for b in BASES if b != s[pos] and "".join(s[c:pos] + [b] + s[pos + 1:c + 3]) not in STOP]
            if alts:
                s[pos] = rng.choice(alts)
            pos += 1 + int(math.log(1 - rng.random()) / log_q)
        return "".join(s)
    for pos in range(3, len(s) - 3):
        if rng.random() >= d * rates[pos]:
            continue
        c = pos - pos % 3
        alts = [b for b in BASES if b != s[pos] and "".join(s[c:pos] + [b] + s[pos + 1:c + 3]) not in STOP]
        if alts:
            s[pos] = rng.choice(alts)
    return "".join(s)


def set_base(s, pos, base):
    """s[pos] = base (s a list) unless that makes its codon a stop; whether it was set."""
    c = pos - pos % 3
    if "".join(s[c:pos] + [base] + s[pos + 1:c + 3]) in STOP:
        return False
    s[pos] = base
    return True


def import_segment(rng, recipient, donor, mean):
    """recipient with a segment of donor (the same gene, the same coordinates: substitutions only) of about `mean`
    bases at a uniform place, the first and last codon kept; a boundary codon that would be a stop keeps the
    recipient's bases. -> (sequence, (start, end))."""
    n = len(recipient)
    length = min(n - 6, 30 + int(rng.expovariate(1.0 / mean)))
    if length <= 0:
        return recipient, None
    start = rng.randint(3, n - 3 - length)
    end = start + length
    s = list(recipient)
    s[start:end] = donor[start:end]
    for c in {start - start % 3, (end - 1) - (end - 1) % 3}:
        if "".join(s[c:c + 3]) in STOP:
            s[c:c + 3] = recipient[c:c + 3]
    return "".join(s), (start, end)


def as_array(seq):
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


def consensus_sites(rep, congeners):
    """protal's ancestry sites of the representative's copy (AncestrySites.h: the copies compared along their shared
    12-mers, then the consensus over the congeners, 0.9 of those compared at a position, the nearest's differences
    where fewer than MIN_CONGENERS were), by ancestry_oracle.py's port. The congeners are the database's
    ones of the genus, nearest first, as species_neighbours.tsv lists them (all within its 0.15 here). -> (mask,
    bases) as arrays (bases as ASCII)."""
    mask = np.zeros(len(rep), dtype=bool)
    bases = np.zeros(len(rep), dtype=np.uint8)
    own_kmers = oracle.unique_kmers(rep)
    comparisons = []
    for other in sorted(congeners, key=lambda c: -sum(a == b for a, b in zip(rep, c))):
        c = oracle.compare_copies(rep, own_kmers, other)
        if c is not None:
            comparisons.append(c)
    if comparisons:
        sites = oracle.consensus(len(rep), comparisons)
        for p, b in zip(sites.positions, sites.bases):
            mask[p] = True
            bases[p] = ord(BASES[b])
    return mask, bases


# ---- trees --------------------------------------------------------------------------------------------------------

class Tree:
    """Nodes with a parent and a height (substitutions per site above the present); a branch is the height
    difference. Parents come before their children."""

    def __init__(self):
        self.name, self.parent, self.height = [], [], []

    def add(self, name, parent, height):
        self.name.append(name)
        self.parent.append(parent)
        self.height.append(height)
        return len(self.name) - 1

    def children(self):
        out = [[] for _ in self.name]
        for i, p in enumerate(self.parent):
            if p is not None:
                out[p].append(i)
        return out

    def ancestors(self, i):
        out = []
        while i is not None:
            out.append(i)
            i = self.parent[i]
        return out

    def mrca_height(self, a, b):
        seen = set(self.ancestors(a))
        for x in self.ancestors(b):
            if x in seen:
                return self.height[x]
        return None


def coalescent(tree, rng, root, n, prefix):
    """n tips (height 0) under `root` by Kingman's coalescent, scaled so that its last merge is `root`. -> the tips."""
    if n == 1:
        return [tree.add(f"{prefix}00", root, 0.0)]
    items, merges, parent, t, next_id = list(range(n)), [], {}, 0.0, n
    while len(items) > 1:
        k = len(items)
        t += rng.expovariate(k * (k - 1) / 2)
        a, b = rng.sample(items, 2)
        items.remove(a)
        items.remove(b)
        parent[a] = parent[b] = next_id
        merges.append((next_id, t))
        items.append(next_id)
        next_id += 1
    scale = tree.height[root] / t
    node_of = {merges[-1][0]: root}
    for mid, tm in reversed(merges[:-1]):
        node_of[mid] = tree.add(f"{prefix}i{mid}", node_of[parent[mid]], tm * scale)
    return [tree.add(f"{prefix}{i:02d}", node_of[parent[i]], 0.0) for i in range(n)]


def evolve(tree, rng, root_seqs, site_rates=None):
    """{marker: [sequence of each node]}: every marker down the tree from the root's sequence (site_rates: per marker
    the sites' rate multipliers, or None)."""
    out = {}
    for m, seq in root_seqs.items():
        rates = site_rates.get(m) if site_rates else None
        seqs = [None] * len(tree.name)
        for i, p in enumerate(tree.parent):
            seqs[i] = seq if p is None else substitute(rng, seqs[p], tree.height[p] - tree.height[i], rates)
        out[m] = seqs
    return out


# ---- the world ----------------------------------------------------------------------------------------------------

class World:
    def __init__(self, args):
        self.args = args
        self.rng = random.Random(args.seed)
        self.names = Names(random.Random(args.seed + 1))
        lengths, by_set = sgr.read_markers(args.markers)
        self.markers = list(by_set["bac120"])
        self.markers_by_set = by_set
        self.root_seq = {m: sgr.random_cds(self.rng, lengths[m]) for m in self.markers}
        # Among-site rate variation (--site_rates ALPHA): a gamma(alpha, mean 1) multiplier per site of each marker,
        # the same on every branch, as real columns keep their rate class across the tree; 0: every site alike.
        self.site_rates = {}
        if args.site_rates > 0:
            for m in self.markers:
                self.site_rates[m] = [self.rng.gammavariate(args.site_rates, 1.0 / args.site_rates) for _ in range(len(self.root_seq[m]))]
        self.node_seq = {}            # (lineage prefix) -> {marker: sequence}
        self.release = []             # the database's genomes (dicts)
        self.query = []               # the other genomes, written as they are made: (accession, lineage, path, length, role, genus)
        self.roles = []               # rows of roles.tsv
        self.genera = []              # rows of genera.tsv
        self.imports = []             # rows of imports.tsv
        self.truth_lines = []         # lines of truth_sites.tsv.gz
        self.focus = {}               # genus id -> {role: [accession]} of the genomes the samples may take
        self.species_index = 0
        self.regimes = {}
        for item in args.regimes.split(","):
            name, within, congeneric = item.split(":")
            self.regimes[name] = (float(within), float(congeneric))
        unknown = [r for r in args.recombination.split(",") if r not in self.regimes]
        if unknown:
            sys.exit(f"--recombination names regimes --regimes lacks: {', '.join(unknown)}")

    # Lineage nodes above the genera, as simulate_gtdb_release evolves them (BRANCH per rank), cached by prefix.
    def lineage_seqs(self, ranks):
        seqs = self.root_seq
        for depth, rank in enumerate(ranks):
            key = ";".join(ranks[:depth + 1])
            if key not in self.node_seq:
                self.node_seq[key] = {m: substitute(self.rng, s, sgr.BRANCH[rank[0]], self.site_rates.get(m)) for m, s in seqs.items()}
            seqs = self.node_seq[key]
        return seqs

    def next_species(self):
        self.species_index += 1
        if self.species_index > 999:
            sys.exit("more than 999 database species: fewer genera, congeners or background genera")
        return self.species_index

    def genome(self, acc, lineage, genes, layout, is_rep, rep_acc, strain):
        """A genome dict as the release writers take it: the genes in the genus's layout between its species'
        spacers, in --contigs pieces."""
        order, strand, spacers, breaks = layout
        contigs, cur = [], [spacers[0]]
        for i, m in enumerate(order):
            seq = genes[m] if strand[m] == "+" else sgr.revcomp(genes[m])
            cur += [seq, spacers[i + 1]]
            if i in breaks:
                contigs.append("".join(cur))
                cur = []
        contigs.append("".join(cur))
        name = lineage.split(";")[-1][3:]
        return {"accession": acc, "gtdb_acc": ("RS_" if is_rep else "GB_") + acc, "rep_gtdb_acc": "RS_" + rep_acc,
                "is_rep": is_rep, "lineage": lineage, "name": name, "strain": strain, "genes": genes,
                "contigs": contigs}

    def layout(self, rng, gc):
        """A genus's gene order and strands, its spacers (the --background_length, shared by its genomes) and its
        contig breaks."""
        order = list(self.markers)
        rng.shuffle(order)
        strand = {m: rng.choice("+-") for m in order}
        spacer = self.args.background_length // (len(order) + 1)
        spacers = [sgr.random_dna(rng, spacer, gc) for _ in range(len(order) + 1)]
        breaks = set(rng.sample(range(len(order) - 1), max(0, self.args.contigs - 1)))
        return order, strand, spacers, breaks

    def write_query(self, g, role, genus_id):
        path = os.path.join(self.args.outdir, "simulation", "query_genomes", f"{g['accession']}_genomic.fna.gz")
        length = sgr.write_genome(path, {"accession": g["accession"], "species": {"name": g["name"]},
                                         "strain": g["strain"], "contigs": g["contigs"]})
        self.query.append((g["accession"], g["lineage"], os.path.abspath(path), length, role, genus_id))

    # -- a test genus --

    def test_genus(self, genus_id, n_congeners, twin, n_alleles, regime, family):
        a, rng = self.args, random.Random(f"{self.args.seed}:genus:{genus_id}")
        genus = self.names.taxon("g")
        prefix = family + [f"g__{genus}"]
        epithets = iter(sorted({self.names.stem((2, 3)) + rng.choice(["us", "is", "a", "um", "ensis", "ii"])
                                for _ in range(64)}))
        lineage = {}

        def species(key):
            lineage[key] = ";".join(prefix + [f"s__{genus} {next(epithets)}"])
            return lineage[key]

        # Heights (substitutions per site): the genus root H, S's split T from its sister, S's crown.
        tree = Tree()
        H = rng.uniform(0.05, 0.06)
        if twin:
            h_T = rng.uniform(0.010, 0.015)
            c_S = rng.uniform(0.003, h_T - 0.004)
        else:
            h_T = rng.uniform(0.025, 0.04)
            c_S = rng.uniform(0.005, 0.02)
        L = h_T - c_S
        root = tree.add("root", None, H)
        T = tree.add("T", root, h_T)
        sister = tree.add("sister", T, 0.0)
        others = []
        if n_congeners - 1 == 1:
            others = [tree.add("other00", root, 0.0)]
        elif n_congeners > 2:
            O = tree.add("O", root, rng.uniform(0.01, H - 0.005))
            others = coalescent(tree, rng, O, n_congeners - 1, "other")
        # S's stem, with the novel species' attachment points; N_1 at the crown.
        attach, prev = {}, T
        for b in sorted(set(a.novel_b)):
            if b <= 0:
                attach[b] = T
            elif b < 1:
                prev = tree.add(f"stem{label(b)}", prev, h_T - b * L)
                attach[b] = prev
        crown = tree.add("crown", prev, c_S)
        novel = {b: tree.add(f"N_{label(b)}", attach.get(b, crown), 0.0) for b in sorted(set(a.novel_b))}
        tips = coalescent(tree, rng, crown, a.tips, "s")

        roles = self.choose_roles(tree, rng, crown, tips, n_alleles)
        genus_seqs = self.lineage_seqs(prefix)
        evolved = evolve(tree, rng, genus_seqs, self.site_rates)
        seq = {i: {m: evolved[m][i] for m in self.markers} for i in tips + [sister, T, crown] + others + list(novel.values())}
        del evolved

        # Recombination: within S, then between species, from the genomes as they were before any import.
        within, congeneric = self.regimes[regime]
        species_of = {i: "S" for i in tips}
        species_of[sister] = "sister"
        for i in others:
            species_of[i] = f"other{i}"
        for b, i in novel.items():
            species_of[i] = f"N_{label(b)}"
        before = {i: dict(seq[i]) for i in species_of}
        by_species = {}
        for i, s in species_of.items():
            by_species.setdefault(s, []).append(i)
        imported = {i: {} for i in species_of}  # node -> {marker: [(start, end, donor node)]}
        for i in sorted(species_of):
            for m in self.markers:
                if species_of[i] == "S" and within > 0 and rng.random() < within:
                    donor = rng.choice([j for j in tips if j != i])
                    seq[i][m], span = import_segment(rng, seq[i][m], before[donor][m], a.segment)
                    if span:
                        imported[i].setdefault(m, []).append((*span, donor))
                if congeneric > 0 and rng.random() < congeneric:
                    donor_species = rng.choice(sorted(s for s in by_species if s != species_of[i]))
                    donor = rng.choice(by_species[donor_species])
                    seq[i][m], span = import_segment(rng, seq[i][m], before[donor][m], a.segment)
                    if span:
                        imported[i].setdefault(m, []).append((*span, donor))

        rep = roles["rep"][0]
        name_of = {i: tree.name[i] for i in seq}
        stem_mask = {m: as_array(seq[crown][m]) != as_array(seq[T][m]) for m in self.markers}

        # Derived genomes: congeneric imports into a strain and into a novel species, retained ancestral states.
        derived = {}  # role -> (sequences, base node, imports {marker: [(start, end, donor name)]})
        base_q = (roles.get("Q_near") or roles.get("Q_lone") or roles.get("Q_deep"))[0]
        n_genes = len(self.markers)
        for p in a.import_shares:
            chosen = set(rng.sample(self.markers, round(p * n_genes)))
            genes, spans = dict(seq[base_q]), {}
            for m in sorted(chosen):
                genes[m], span = import_segment(rng, genes[m], seq[sister][m], a.segment)
                if span:
                    spans.setdefault(m, []).append((*span, sister))
            derived[f"Q_imp{label(p)}"] = (genes, base_q, spans)
            if 0.5 in novel:
                genes, spans = dict(seq[novel[0.5]]), {}
                for m in sorted(chosen):
                    donor = rng.choice(tips)
                    genes[m], span = import_segment(rng, genes[m], seq[donor][m], a.segment)
                    if span:
                        spans.setdefault(m, []).append((*span, donor))
                derived[f"N_imp{label(p)}"] = (genes, novel[0.5], spans)
        if a.ils_share > 0 and roles.get("Q_deep"):
            q = roles["Q_deep"][0]
            genes = {}
            for m in self.markers:
                s, stem_t = list(seq[q][m]), seq[T][m]
                sites = [p for p in np.flatnonzero(stem_mask[m]) if s[p] == seq[crown][m][p]]
                for p in rng.sample(sites, round(a.ils_share * len(sites))):
                    set_base(s, int(p), stem_t[p])
                genes[m] = "".join(s)
            derived["Q_ils"] = (genes, q, {})

        # The database's copies: S's representative, its allele genomes, the congeners' representatives.
        db_congeners = [sister] + others
        cons = {}
        for m in self.markers:
            cons[m] = consensus_sites(seq[rep][m], [seq[c][m] for c in db_congeners])
        alleles = roles.get("allele", [])

        # Planted genomes (genera with the consensus, 1-4 allele genomes, no twin).
        planted = {}
        if a.planted and n_congeners >= MIN_CONGENERS and 1 <= n_alleles <= 4 and not twin and alleles:
            far = max(alleles, key=lambda i: sum(int((as_array(seq[i][m]) != as_array(seq[rep][m])).sum())
                                                 for m in self.markers))
            for f in a.planted:
                planted[f"P_f{label(f)}"] = self.plant(rng, seq, rep, alleles, cons, far=far, f=f)
            for share in a.planted:
                planted[f"A_a{label(share)}"] = self.plant(rng, seq, rep, alleles, cons, agree=share)

        # Species and accessions. The database: S, the sister, the other congeners.
        s_index = self.next_species()
        s_lineage = species("S")
        rep_acc = f"GCF_999{s_index:03d}001.1"
        gc = 0.38 + rng.random() * 0.26
        layout = self.layout(rng, gc)
        acc_of = {rep: rep_acc}
        self.release.append(self.genome(rep_acc, s_lineage, seq[rep], layout, True, rep_acc, "SIM-rep"))
        for k, i in enumerate(alleles, 2):
            acc_of[i] = f"GCA_999{s_index:03d}{k:03d}.1"
            self.release.append(self.genome(acc_of[i], s_lineage, seq[i], layout, False, rep_acc, f"SIM-{name_of[i]}"))
        for key, i in [("sister", sister)] + [(f"other{k}", i) for k, i in enumerate(others)]:
            idx = self.next_species()
            acc_of[i] = f"GCF_999{idx:03d}001.1"
            self.release.append(self.genome(acc_of[i], species(key), seq[i], layout, True, acc_of[i], "SIM-rep"))
        # The query genomes (never in the release): accessions in the 998 range.
        counter = itertools.count(1)
        queries = []  # (role, node or None, genes, lineage, imports)
        for role in ("Q_near", "Q_lone", "Q_deep", "Q_rep"):
            for i in roles.get(role, []):
                queries.append((role, i, seq[i], s_lineage, imported[i]))
        for b, i in novel.items():
            queries.append((f"N_{label(b)}", i, seq[i], species(f"N_{label(b)}"), imported[i]))
        for i in roles.get("decoy", []):
            queries.append(("decoy", i, seq[i], species("decoy"), imported[i]))
        for role, (genes, base, spans) in derived.items():
            line = lineage.get(f"N_{label(0.5)}") if role.startswith("N_") else s_lineage
            merged = {m: list(imported[base].get(m, [])) for m in self.markers}
            for m, v in spans.items():
                merged[m] += v
            queries.append((role, None, genes, line, {m: v for m, v in merged.items() if v}))
        for role, genes in planted.items():
            queries.append((role, None, genes, s_lineage, {}))
        focus = {}
        for role, node, genes, line, imp in queries:
            acc = f"GCA_998{genus_id:03d}{next(counter):03d}.1"
            if node is not None:
                acc_of[node] = acc
            g = self.genome(acc, line, genes, layout, False, rep_acc, f"SIM-{role}")
            self.write_query(g, role, genus_id)
            focus.setdefault(role, []).append(acc)
            self.genome_truth(genus_id, acc, role, node, genes, imp, tree, seq, rep, crown, T, stem_mask, cons,
                              alleles, novel, roles, name_of, acc_of)
        if alleles:
            focus["Q_leak"] = [acc_of[alleles[0]]]
        for i in alleles:
            self.genome_truth(genus_id, acc_of[i], "allele", i, seq[i], imported[i], tree, seq, rep, crown, T,
                              stem_mask, cons, alleles, novel, roles, name_of, acc_of)
        self.focus[genus_id] = focus
        for i in species_of:
            for m, spans in imported[i].items():
                for start, end, donor in spans:
                    self.imports.append((acc_of.get(i, name_of[i]), m, start, end, acc_of.get(donor, name_of[donor])))
        for role, (genes, base, spans) in derived.items():
            for m, v in spans.items():
                for start, end, donor in v:
                    self.imports.append((focus[role][0], m, start, end, acc_of.get(donor, name_of[donor])))

        # Truth sites, per marker: S's stem, the representative's lineage, the consensus sites, the allele genomes.
        for m in self.markers:
            r, cr, t = as_array(seq[rep][m]), as_array(seq[crown][m]), as_array(seq[T][m])
            lines = [("stem", "", np.flatnonzero(cr != t), t), ("rep_lineage", "", np.flatnonzero(r != cr), cr),
                     ("consensus", "", np.flatnonzero(cons[m][0]), cons[m][1])]
            for i in alleles:
                x = as_array(seq[i][m])
                lines.append(("allele", acc_of[i], np.flatnonzero(x != r), x))
            for kind, acc, pos, bases in lines:
                text = ",".join(f"{p}{chr(bases[p])}" for p in pos)
                self.truth_lines.append(f"{genus_id}\t{m}\t{kind}\t{acc}\t{text}\n")

        self.genera.append({
            "genus": genus_id, "genus_name": genus, "target": s_lineage, "target_rep": rep_acc,
            "congeners": n_congeners, "twin": int(twin), "allele_genomes": n_alleles, "alleles_made": len(alleles),
            "regime": regime, "within_rate": within, "congeneric_rate": congeneric, "H": round(H, 5),
            "h_T": round(h_T, 5), "c_S": round(c_S, 5), "L": round(L, 5), "planted": int(bool(planted)),
            "sister": acc_of[sister], "congener_reps": ",".join(acc_of[i] for i in db_congeners),
            "roles": ";".join(f"{r}:{len(v)}" for r, v in sorted(focus.items()))})

    def choose_roles(self, tree, rng, crown, tips, n_alleles):
        """The roles of S's genealogy's tips (see the module text). -> {role: [node]}."""
        children = tree.children()
        under = {}

        def tips_under(i):
            if i not in under:
                under[i] = [i] if not children[i] else [t for c in children[i] for t in tips_under(c)]
            return under[i]

        rep = rng.choice(tips)
        # The crown's two sides of the genealogy (N_1, the decoy, leaves the crown too).
        side_a, side_b = [c for c in children[crown] if not tree.name[c].startswith("N_")]
        if rep in tips_under(side_b):
            side_a, side_b = side_b, side_a
        # Q_deep: a clade of the crown's other side kept free of allele genomes (2-4 tips, the deepest such).
        candidates = [i for i in [side_b] + [d for d in range(len(tree.name)) if side_b in tree.ancestors(d)[1:]]
                      if 2 <= len(tips_under(i)) <= 4]
        x_node = max(candidates, key=lambda i: (tree.height[i], -i)) if candidates else side_b
        x = tips_under(x_node)
        deep = rng.sample(x, min(2, len(x)))
        # The decoy: a tip of the same clade, drawn as the deep strains are, labelled a species of its own (with two
        # tips in the clade, one deep strain and the decoy).
        decoy = []
        if self.args.decoy and len(x) >= 2:
            if len(x) == 2:
                deep = deep[:1]
            decoy = [rng.choice([t for t in x if t not in deep])]
        used = {rep, *x}
        near_rep = min((t for t in tips if t not in used), key=lambda t: (tree.mrca_height(t, rep), t))
        used.add(near_rep)
        genealogy = set(tips)
        cherries = [i for i in range(len(tree.name)) if len(children[i]) == 2
                    and all(c in genealogy for c in children[i]) and not set(children[i]) & used]
        rng.shuffle(cherries)
        n_near = min(3, n_alleles)
        roles = {"rep": [rep], "Q_deep": deep, "decoy": decoy, "Q_rep": [near_rep], "Q_near": [], "Q_lone": [],
                 "allele": [], "pop": []}
        for k, c in enumerate(cherries[:3]):
            a, q = rng.sample(children[c], 2)
            roles["Q_near" if k < n_near else "Q_lone"].append(q)
            if k < n_near:
                roles["allele"].append(a)
            used |= {a, q}
        free = [t for t in tips if t not in used]
        rng.shuffle(free)
        roles["allele"] += free[:max(0, n_alleles - len(roles["allele"]))]
        roles["pop"] = [t for t in tips if t not in used and t not in roles["allele"]]
        return {r: v for r, v in roles.items() if v}

    def plant(self, rng, seq, rep, alleles, cons, far=None, f=None, agree=None):
        """A planted genome from S's representative: with `far`, a share f of that allele genome's differences; with
        `agree`, the consensus base at a share 1 - agree of the consensus sites; then new differences, away from every
        allele genome's differences and the consensus sites, up to --planted_divergence."""
        genes = {}
        for m in self.markers:
            r = as_array(seq[rep][m])
            s = list(seq[rep][m])
            mask, bases = cons[m]
            poly = np.zeros(len(r), dtype=bool)
            for i in alleles:
                poly |= as_array(seq[i][m]) != r
            done = 0
            if far is not None:
                x = seq[far][m]
                diffs = [int(p) for p in np.flatnonzero(as_array(x) != r)]
                for p in rng.sample(diffs, round(f * len(diffs))):
                    done += set_base(s, p, x[p])
            else:
                sites = [int(p) for p in np.flatnonzero(mask)]
                for p in rng.sample(sites, round((1 - agree) * len(sites))):
                    done += set_base(s, p, chr(bases[p]))
            target = round(self.args.planted_divergence * len(r))
            free = [int(p) for p in np.flatnonzero(~poly & ~mask) if 3 <= p < len(r) - 3]
            rng.shuffle(free)
            for p in free:
                if done >= target:
                    break
                done += set_base(s, p, rng.choice([b for b in BASES if b != s[p]]))
            genes[m] = "".join(s)
        return genes

    def genome_truth(self, genus_id, acc, role, node, genes, imports, tree, seq, rep, crown, T, stem_mask, cons,
                     alleles, novel, roles, name_of, acc_of):
        """roles.tsv's row of a genome: its differences from the target's representative by class, summed over the
        markers, and where it sits in the trees."""
        c = dict.fromkeys(["length", "diffs", "at_rep_lineage", "at_rep_lineage_poly", "at_stem", "known", "other",
                           "imported_diffs",
                           "cons_sites", "cons_agree", "cons_congener", "fixed_sites", "fixed_agree", "fixed_congener",
                           "poly_sites", "poly_known", "poly_novel", "explainable"], 0)
        genes_imported = genes_congeneric = 0
        # the genome's own species: S's genealogy for S's genomes; a novel species is a species of one genome
        own_species = set() if role.startswith("N_") else {i for v in roles.values() for i in v}
        for m in self.markers:
            g, r = as_array(genes[m]), as_array(seq[rep][m])
            cr, t = as_array(seq[crown][m]), as_array(seq[T][m])
            diff = g != r
            rep_lineage = r != cr
            at_r = diff & rep_lineage & (g == cr)
            at_stem = diff & stem_mask[m] & (g == t)
            poly = np.zeros(len(r), dtype=bool)
            same = np.zeros(len(r), dtype=bool)
            for i in alleles:
                if acc_of.get(i) == acc:
                    continue  # an allele genome is not its own allele
                x = as_array(seq[i][m])
                poly |= x != r
                same |= (x != r) & (x == g)
            known = diff & ~at_r & ~at_stem & same
            mask, bases = cons[m]
            fixed = mask & ~poly
            imported = np.zeros(len(r), dtype=bool)
            for start, end, _ in imports.get(m, []):
                imported[start:end] = True
            genes_imported += bool(imports.get(m))
            genes_congeneric += any(donor not in own_species for _, _, donor in imports.get(m, []))
            c["length"] += len(r)
            c["diffs"] += int(diff.sum())
            c["at_rep_lineage"] += int(at_r.sum())
            c["at_rep_lineage_poly"] += int((at_r & poly).sum())  # those the allele genomes make polymorphic
            c["at_stem"] += int(at_stem.sum())
            c["known"] += int(known.sum())
            c["other"] += int((diff & ~at_r & ~at_stem & ~known).sum())
            c["imported_diffs"] += int((diff & imported).sum())
            c["cons_sites"] += int(mask.sum())
            c["cons_agree"] += int((mask & ~diff).sum())
            c["cons_congener"] += int((mask & (g == bases)).sum())
            c["fixed_sites"] += int(fixed.sum())
            c["fixed_agree"] += int((fixed & ~diff).sum())
            c["fixed_congener"] += int((fixed & (g == bases)).sum())
            c["poly_sites"] += int(poly.sum())
            c["poly_known"] += int((poly & diff & same).sum())
            c["poly_novel"] += int((poly & diff & ~same).sum())
            c["explainable"] += int((diff & same).sum())
        b = f = share = ""
        if role.startswith("N_") and not role.startswith("N_imp"):
            b = role[2:]
        if role.startswith("N_imp"):
            b = "0.5"
        if role.startswith("P_f"):
            f = role[3:]
        if role.startswith("A_a"):
            share = role[3:]
        mrca_rep = nearest_allele = nearest_acc = ""
        if role.startswith(("Q_imp", "Q_ils")):  # where the genome it was derived from sits
            node = (roles.get("Q_deep") if role == "Q_ils" else (roles.get("Q_near") or roles.get("Q_lone") or
                                                                  roles.get("Q_deep")))[0]
        if role.startswith("N_imp") and 0.5 in novel:
            node = novel[0.5]
        if node is not None and node in novel.values():
            mrca_rep = round(tree.height[tree.parent[node]], 6)
        elif node is not None and node != rep:
            mrca_rep = round(tree.mrca_height(node, rep), 6)
            others = [i for i in alleles if i != node]
            if others:
                height, nearest = min((tree.mrca_height(node, i), acc_of[i]) for i in others)
                nearest_allele, nearest_acc = round(height, 6), nearest
        self.roles.append({"accession": acc, "genus": genus_id, "role": role, "class": role_class(role), "b": b,
                           "f": f, "a": share, "mrca_rep": mrca_rep, "mrca_nearest_allele": nearest_allele,
                           "nearest_allele": nearest_acc, "genes_imported": genes_imported,
                           "genes_congeneric": genes_congeneric, "distance": round(c["diffs"] / c["length"], 6), **c})

    # -- background genera --

    def background_genus(self, k, family):
        rng = random.Random(f"{self.args.seed}:background:{k}")
        genus = self.names.taxon("g")
        prefix = family + [f"g__{genus}"]
        genus_seqs = self.lineage_seqs(prefix)
        layout = self.layout(rng, 0.38 + rng.random() * 0.26)
        out = []
        epithets = sorted({self.names.stem((2, 3)) + "us" for _ in range(self.args.background_species * 3)})
        for epithet in epithets[:self.args.background_species]:
            idx = self.next_species()
            acc = f"GCF_999{idx:03d}001.1"
            genes = {m: substitute(rng, s, sgr.BRANCH["s"], self.site_rates.get(m)) for m, s in genus_seqs.items()}
            g = self.genome(acc, ";".join(prefix + [f"s__{genus} {epithet}"]), genes, layout, True, acc, "SIM-rep")
            self.release.append(g)
            out.append(acc)
        return out

    def build(self):
        a = self.args
        names = self.names
        phylum = ["d__Bacteria", "p__" + names.taxon("p"), "c__" + names.taxon("c")]
        cells = list(itertools.product(a.congeners, a.twins, a.allele_genomes, a.recombination.split(",")))
        genus_id = 0
        orders, families = {}, {}
        for rep_index in range(a.replicates):
            for n_c, twin, n_a, regime in cells:
                genus_id += 1
                if a.genera_per_family <= 1:
                    order = orders.setdefault(genus_id % 4, "o__" + names.taxon("o"))
                    family = phylum + [order, "f__" + names.taxon("f")]
                else:
                    # Families of --genera_per_family test genera, nested in four orders.
                    family_index = (genus_id - 1) // a.genera_per_family
                    order = orders.setdefault(family_index % 4, "o__" + names.taxon("o"))
                    family = families.setdefault(family_index, phylum + [order, "f__" + names.taxon("f")])
                self.test_genus(genus_id, n_c, bool(twin), n_a, regime, family)
                sys.stderr.write(f"genus {genus_id}: {n_c} congeners, twin {twin}, {n_a} allele genomes, {regime}\n")
        self.background = []
        bg_order = "o__" + names.taxon("o")
        for k in range(a.background_genera):
            self.background += self.background_genus(k, phylum + [bg_order, "f__" + names.taxon("f")])
        self.n_test_genera = genus_id


# ---- the release --------------------------------------------------------------------------------------------------

def write_metadata(path, genomes):
    cols = ["accession", "ambiguous_bases", "checkm2_completeness", "checkm2_contamination", "contig_count",
            "gc_percentage", "genome_size", "gtdb_genome_representative", "gtdb_representative", "gtdb_taxonomy",
            "ncbi_organism_name"]
    with sgr.gzip_text(path) as fh:
        fh.write("\t".join(cols) + "\n")
        for g in genomes:
            seq = "".join(g["contigs"])
            gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
            fh.write("\t".join([g["gtdb_acc"], "0", "100.00", "0.00", str(len(g["contigs"])), f"{gc:.2f}",
                                str(len(seq)), g["rep_gtdb_acc"], "t" if g["is_rep"] else "f", g["lineage"],
                                f"{g['name']} strain {g['strain']}"]) + "\n")


def write_release(world, out, release):
    rel = f"r{release}"
    with open(os.path.join(out, "VERSION.txt"), "w", newline="\n") as fh:
        fh.write(f"v{release}.0\nSynthetic GTDB-like release grown on trees by simulate_ancestry_world.py (seed "
                 f"{world.args.seed}, {len(world.release)} genomes)\n")
    genomes = [{**g, "species": {"lineage": g["lineage"], "name": g["name"]}} for g in world.release]
    for mset in ("bac120", "ar53"):
        chosen = genomes if mset == "bac120" else []
        sgr.write_taxonomy(os.path.join(out, f"{mset}_taxonomy_{rel}.tsv"), chosen)
        write_metadata(os.path.join(out, f"{mset}_metadata_{rel}.tsv.gz"), chosen)
        sgr.write_marker_files(os.path.join(out, "genomic_files_reps", f"{mset}_marker_genes_reps_{rel}"), mset,
                               world.markers_by_set[mset], [g for g in chosen if g["is_rep"]])
        sgr.write_marker_files(os.path.join(out, "genomic_files_all", f"{mset}_marker_genes_all_{rel}"), mset,
                               world.markers_by_set[mset], chosen)
    sim = os.path.join(out, "simulation")
    reps_db = os.path.join(out, "genomic_files_reps", f"gtdb_genomes_reps_{rel}", "database")
    rows = []
    with open(os.path.join(sim, "genomes.tsv"), "w", newline="\n") as gt:
        gt.write("accession\tgtdb_taxonomy\tfasta_path\tgenome_length\tgtdb_representative\n")
        for g in genomes:
            fasta = (sgr.genome_path(reps_db, g["accession"]) if g["is_rep"] else
                     os.path.join(sim, "genomes_nonreps", f"{g['accession']}_genomic.fna.gz"))
            length = sgr.write_genome(fasta, g)
            gt.write(f"{g['accession']}\t{g['lineage']}\t{os.path.abspath(fasta)}\t{length}\t"
                     f"{'t' if g['is_rep'] else 'f'}\n")
            rows.append((g["accession"], g["lineage"], os.path.abspath(fasta), length,
                         "rep" if g["is_rep"] else "allele", ""))
    return rows


def write_tables(world, out, release_rows):
    sim = os.path.join(out, "simulation")
    genus_of = {r["accession"]: r["genus"] for r in world.roles}
    for g in world.genera:
        genus_of[g["target_rep"]] = g["genus"]
        for acc in g["congener_reps"].split(","):
            genus_of[acc] = g["genus"]
    role_of = {r["accession"]: r["role"] for r in world.roles}
    for g in world.genera:
        role_of[g["target_rep"]] = "rep"
        role_of[g["sister"]] = "twin" if g["twin"] else "congener"
        for acc in g["congener_reps"].split(",")[1:]:
            role_of[acc] = "congener"
    with open(os.path.join(sim, "all_genomes.tsv"), "w", newline="\n") as fh:
        fh.write("accession\tgtdb_taxonomy\tfasta_path\tgenome_length\trole\tgenus\tin_release\n")
        for acc, lineage, path, length, role, _ in release_rows:
            fh.write(f"{acc}\t{lineage}\t{path}\t{length}\t{role_of.get(acc, 'background')}\t"
                     f"{genus_of.get(acc, '')}\t1\n")
        for acc, lineage, path, length, role, genus in world.query:
            fh.write(f"{acc}\t{lineage}\t{path}\t{length}\t{role}\t{genus}\t0\n")
    with open(os.path.join(sim, "genera.tsv"), "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(world.genera[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(world.genera)
    with open(os.path.join(sim, "roles.tsv"), "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(world.roles[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(world.roles)
    with open(os.path.join(sim, "imports.tsv"), "w", newline="\n") as fh:
        fh.write("accession\tmarker\tstart\tend\tdonor\n")
        for row in world.imports:
            fh.write("\t".join(map(str, row)) + "\n")
    with sgr.gzip_text(os.path.join(sim, "truth_sites.tsv.gz")) as fh:
        fh.write("genus\tmarker\tkind\tgenome\tsites (0-based position and base: the ancestral base for stem, the "
                 "crown's for rep_lineage, the congeners' for consensus, the allele genome's for allele)\n")
        fh.writelines(world.truth_lines)


# ---- the samples --------------------------------------------------------------------------------------------------

def balanced(rng, items, n):
    """n draws cycling through `items` in shuffled rounds, so that each comes up as evenly as n allows."""
    if not items:
        sys.exit("a test genus has no genome of a class the samples need")
    out = []
    while len(out) < n:
        round_ = list(items)
        rng.shuffle(round_)
        out += round_
    return out[:n]


def write_samples(world, out):
    a = world.args
    rng = random.Random(f"{a.seed}:samples")
    lengths, lineage_of = {}, {}
    with open(os.path.join(out, "simulation", "all_genomes.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lengths[row["accession"]] = int(row["genome_length"])
            lineage_of[row["accession"]] = (row["gtdb_taxonomy"], row["fasta_path"])
    # Every genome has every marker at its root length (substitutions only).
    marker_bases = sum(len(s) for s in world.root_seq.values())
    targets = {g["genus"]: g["target"] for g in world.genera}
    design = []  # (set, sample, genus, genome, role, class, depth)

    def focus_lists(genus, classes):
        focus = world.focus[genus]
        return [(acc, role) for role, accs in sorted(focus.items()) if role_class(role) in classes for acc in accs]

    def draws(genus, n, depths, with_planted=False):
        q = [(acc, role, d) for acc, role in focus_lists(genus, ("Q",)) for d in depths]
        nn = [(acc, role, d) for acc, role in focus_lists(genus, ("N", "D")) for d in depths]
        if with_planted:
            p = [(acc, role, d) for acc, role in focus_lists(genus, ("P",)) for d in depths]
            if any(r.startswith(("P_", "A_")) for _, r, _ in p):
                return balanced(rng, p, n)
            # without planted genomes: the leaky control and the oracle's strains and novel species
            return balanced(rng, p + q[:len(depths) * 2] + nn[:len(depths) * 2], n)
        half = n // 2
        out = balanced(rng, q, n - half) + balanced(rng, nn, half)
        rng.shuffle(out)
        return out

    def present(genus, acc, role):
        return targets[genus] if role_class(role) in ("Q", "P") else lineage_of[acc][0]

    def communities(set_name, n, depths, with_planted=False, background=True):
        per_genus = {g["genus"]: draws(g["genus"], n, depths, with_planted) for g in world.genera}
        samples = []
        for s in range(n):
            name = f"{set_name}{s + 1:03d}"
            members = []
            for genus, picks in sorted(per_genus.items()):
                acc, role, depth = picks[s]
                members.append((acc, role, depth, genus, present(genus, acc, role)))
                design.append((set_name, name, genus, acc, role, role_class(role), depth))
            if background and world.background:
                for acc in rng.sample(world.background, min(a.background_per_sample, len(world.background))):
                    depth = rng.choice(depths)
                    members.append((acc, "background", depth, "", lineage_of[acc][0]))
                    design.append((set_name, name, "", acc, "background", "", depth))
            samples.append((name, members))
        return samples

    def write_truth(folder, name, members):
        os.makedirs(os.path.join(folder, "truth"), exist_ok=True)
        path = os.path.join(folder, "truth", f"{name}.tsv")
        with open(path, "w", newline="\n") as fh:
            for line in sorted({m[4] for m in members}):
                fh.write(line + "\n")
        return os.path.abspath(path)

    def pairs_for(acc, fragments):
        """Read pairs that put about `fragments` on the genome's marker genes."""
        return max(1, math.ceil(fragments * lengths[acc] / marker_bases))

    # Paired-end samples: a manifest for simulate_metagenomes --from_manifest.
    folder = os.path.join(out, "samples", "pe")
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "manifest.tsv"), "w", newline="\n") as fh:
        fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\t"
                 "relative_abundance\tfastq_r1\tfastq_r2\tfasta_path\n")
        for name, members in communities("pe", a.samples, a.depths):
            write_truth(folder, name, members)
            rows = [(acc, pairs_for(acc, depth)) for acc, _, depth, _, _ in members]
            total = sum(p * 300 / lengths[acc] for acc, p in rows)
            for acc, pairs in rows:
                lineage, path = lineage_of[acc]
                cov = pairs * 300 / lengths[acc]
                fh.write(f"{name}\t{acc}\t{lineage.split(';')[-1]}\t{lineage}\t{lengths[acc]}\t{pairs}\t{cov:.6g}\t"
                         f"{cov / total:.6g}\t\t\t{path}\n")
    # HiFi samples: simulate_metagenomes --long_samples / --long_genomes (weights: coverage x length).
    folder = os.path.join(out, "samples", "hifi")
    os.makedirs(os.path.join(folder, "reads"), exist_ok=True)
    with open(os.path.join(folder, "long_samples.tsv"), "w", newline="\n") as fs, \
         open(os.path.join(folder, "long_genomes.tsv"), "w", newline="\n") as fg:
        fs.write("sample\tout\tbases\tseed\n")
        fg.write("sample\tgenome\tfasta\tweight\thost\n")
        for k, (name, members) in enumerate(communities("hifi", a.long_samples, a.long_depths)):
            write_truth(folder, name, members)
            bases = 0
            for acc, _, cov, _, _ in members:
                weight = cov * lengths[acc]
                bases += weight
                fg.write(f"{name}\t{acc}\t{lineage_of[acc][1]}\t{weight:.6g}\t0\n")
            fs.write(f"{name}\t{os.path.abspath(os.path.join(folder, 'reads', name + '.fq.gz'))}\t{int(bases)}\t"
                     f"{a.seed * 1000 + k + 1}\n")
    # Error-free paired-end samples: simulate_reads.py communities (abundance x length = the pairs wanted).
    folder = os.path.join(out, "samples", "exact")
    os.makedirs(os.path.join(folder, "communities"), exist_ok=True)
    with open(os.path.join(folder, "samples.tsv"), "w", newline="\n") as fh:
        fh.write("sample\tcommunity\tpairs\tseed\ttruth\n")
        for k, (name, members) in enumerate(communities("exact", a.exact_samples, [a.exact_depth], with_planted=True,
                                                        background=False)):
            truth = write_truth(folder, name, members)
            path = os.path.join(folder, "communities", f"{name}.tsv")
            pairs = 0
            with open(path, "w", newline="\n") as cf:
                cf.write("accession\trelative_abundance\n")
                for acc, _, depth, _, _ in members:
                    n = pairs_for(acc, depth)
                    pairs += n
                    cf.write(f"{acc}\t{n / lengths[acc]:.9g}\n")
            fh.write(f"{name}\t{os.path.abspath(path)}\t{pairs}\t{a.seed * 1000 + k + 1}\t{truth}\n")
    with open(os.path.join(out, "samples", "design.tsv"), "w", newline="\n") as fh:
        fh.write("set\tsample\tgenus\tgenome\trole\tclass\tdepth\n")
        for row in design:
            fh.write("\t".join(map(str, row)) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--release", default="226")
    p.add_argument("--markers", default=os.path.join(HERE, "markers_r226.tsv"))
    p.add_argument("--congeners", type=ints, default=[2, 6], help="database congeners of the target (default 2,6)")
    p.add_argument("--twins", type=ints, default=[0, 1], help="0: the sister 0.05-0.08 away, 1: a twin 0.02-0.03 away")
    p.add_argument("--allele_genomes", type=ints, default=[0, 2, 8], help="allele genomes per target (default 0,2,8)")
    p.add_argument("--recombination", default="none,low,high", help="regimes of --regimes, crossed with the rest")
    p.add_argument("--regimes", default=REGIMES,
                   help=f"name:within:congeneric, the per-gene probability of an import (default {REGIMES})")
    p.add_argument("--replicates", type=int, default=1, help="genera per cell of the factors (default 1)")
    p.add_argument("--site_rates", type=float, default=0.0, metavar="ALPHA",
                   help="among-site rate variation: a gamma(ALPHA, mean 1) rate multiplier per site of each marker, the "
                        "same on every branch (default 0: every site alike; 0.5-1 is realistic). For the column weights "
                        "(protal --build --column_weights), which see nothing without it")
    p.add_argument("--genera_per_family", type=int, default=1,
                   help="test genera per family (default 1: a family each, as before 2026-10-10); with more, the genera "
                        "of a family share its node's sequences (simulate_gtdb_release's branch per rank) and the column "
                        "weights' family references come from the other genera")
    p.add_argument("--tips", type=int, default=24, help="genomes in the target's genealogy (default 24)")
    p.add_argument("--novel_b", type=floats, default=[0, 0.5, 0.8, 0.95, 1.0],
                   help="where the novel species leave the target's stem, as the share of it they keep (default "
                        "0,0.5,0.8,0.95,1; 1: at the crown, a decoy)")
    p.add_argument("--no_decoy", dest="decoy", action="store_false",
                   help="no decoy (a tip of the deep strains' clade labelled a species of its own: no feature should "
                        "tell it from them); for a training world")
    p.add_argument("--import_shares", type=floats, default=[0.2, 0.4],
                   help="genes with a congener's segment in Q_imp and N_imp genomes (default 0.2,0.4)")
    p.add_argument("--segment", type=float, default=400, help="mean length of an imported segment (default 400)")
    p.add_argument("--ils_share", type=float, default=0.1,
                   help="share of the stem sites a Q_ils genome has the ancestral base at (default 0.1; 0: none)")
    p.add_argument("--planted", type=floats, default=[0, 0.25, 0.5, 0.75, 1.0],
                   help="f and a of the planted genomes (default 0,0.25,0.5,0.75,1; empty: none)")
    p.add_argument("--planted_divergence", type=float, default=0.02)
    p.add_argument("--background_genera", type=int, default=8)
    p.add_argument("--background_species", type=int, default=3)
    p.add_argument("--background_length", type=int, default=40000, help="non-marker DNA per genome (default 40000)")
    p.add_argument("--contigs", type=int, default=3)
    p.add_argument("--samples", type=int, default=96, help="paired-end samples (default 96)")
    p.add_argument("--long_samples", type=int, default=48, help="HiFi samples (default 48)")
    p.add_argument("--exact_samples", type=int, default=24, help="error-free paired-end samples (default 24)")
    p.add_argument("--depths", type=floats, default=[3, 10, 30, 100],
                   help="marker fragments per focus genome (default 3,10,30,100)")
    p.add_argument("--long_depths", type=floats, default=[0.5, 1, 2, 4],
                   help="HiFi vertical coverage per focus genome (default 0.5,1,2,4)")
    p.add_argument("--exact_depth", type=float, default=100)
    p.add_argument("--background_per_sample", type=int, default=10)
    p.add_argument("--quick", action="store_true",
                   help="a small world for the tests: 2 genera (6 congeners, 2 allele genomes, no twin, regimes none "
                        "and high), 12 tips, 1 background genus, 8+4+6 samples")
    a = p.parse_args(argv)
    if a.quick:
        a.congeners, a.twins, a.allele_genomes, a.recombination = [6], [0], [2], "none,high"
        a.tips, a.background_genera, a.samples, a.long_samples, a.exact_samples = 12, 1, 8, 4, 6
        a.depths, a.long_depths = [10, 100], [1, 2]
    if a.tips < 8:
        sys.exit("--tips must be at least 8")
    os.makedirs(os.path.join(a.outdir, "simulation"), exist_ok=True)
    world = World(a)
    world.build()
    release_rows = write_release(world, a.outdir, a.release)
    write_tables(world, a.outdir, release_rows)
    write_samples(world, a.outdir)
    n_query = len(world.query)
    sys.stderr.write(f"Wrote {a.outdir}: {world.n_test_genera} test genera and {a.background_genera} background "
                     f"genera, {len(world.release)} release genomes, {n_query} query genomes; samples "
                     f"{a.samples} pe, {a.long_samples} HiFi, {a.exact_samples} error-free\n")


if __name__ == "__main__":
    main()
