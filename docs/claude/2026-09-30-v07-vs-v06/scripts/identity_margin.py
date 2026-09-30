#!/usr/bin/env python3
"""identity_margin.py - what --depth_identity_margin keeps of a taxon's reads, by where the reads come from.

Usage: identity_margin.py SAM MANIFEST SAMPLE TAXONOMY DIVERGENCE [MARGIN ...]

Reads the SAM protal wrote (plain or .zst), takes every primary alignment, and gives each read protal's
identity (Strain.h AlignmentIdentity: M / (M + X + I + D), soft clips left out) and its aligned reference
length. A taxon's threshold is, as in Profiler.h, the 98th percentile of its reads' identities (by aligned
bases) minus the margin. The simulated read names start with their genome's accession, so every read is
either the taxon's own (a genome of its species: the representative, which is the database's reference,
or another strain) or a relative's (a genome of another species). For each margin: the share of aligned
bases kept, for own reads by the genome's expected difference from the reference (its divergence from the
species ancestor plus the representative's; 0 for the representative), and for relatives' reads.
"""
import collections
import csv
import re
import subprocess
import sys

sam, manifest, sample, taxonomy, divergence = sys.argv[1:6]
margins = [float(m) for m in sys.argv[6:]] or [0.04, 0.08]

species_of_genome, rep_of_species = {}, {}
with open(manifest) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["sample"] == sample:
            species_of_genome[row["genome"]] = row["species"]
div = {}
with open(divergence) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        div[row["accession"]] = float(row["strain_divergence"])
        if row["accession"].startswith("GCF_"):
            rep_of_species[row["species"]] = row["accession"]
name_of_taxid = {}
with open(taxonomy) as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[4] == "species":
            name_of_taxid[f[0]] = f[3].removeprefix("s__")


def identity(cigar):
    matches = diff = ref = 0
    for n, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar):
        n = int(n)
        if op == "M":
            matches += n
        if op in "XID":
            diff += n
        if op in "MXD":
            ref += n
    cols = matches + diff
    return (matches / cols if cols else 0.0), ref


reads = collections.defaultdict(list)  # taxid -> [(identity, length, source genome)]
opener = subprocess.Popen(["zstd", "-dc", sam], stdout=subprocess.PIPE, text=True) if sam.endswith(".zst") else None
fh = opener.stdout if opener else open(sam)
for line in fh:
    if line.startswith("@"):
        continue
    f = line.split("\t", 6)
    flag = int(f[1])
    if flag & 0x904 or f[2] == "*":  # secondary, supplementary, unmapped
        continue
    genome = re.match(r"(GC[AF]_\d{9}\.\d+)", f[0]).group(1)
    ident, length = identity(f[5])
    reads[f[2].split("_")[0]].append((ident, length, genome))

own = collections.defaultdict(lambda: collections.defaultdict(float))       # margin -> bin -> kept bases
own_all = collections.defaultdict(float)
rel = collections.defaultdict(float)
rel_all = 0.0
bins = [(-1, 0.0005, "reference itself"), (0.0005, 0.01, "0-1%"), (0.01, 0.02, "1-2%"), (0.02, 0.03, "2-3%"),
        (0.03, 0.045, "3-4.5%")]
for taxid, rs in reads.items():
    species = name_of_taxid.get(taxid)
    rs_sorted = sorted(rs)
    total = sum(r[1] for r in rs)
    cumulative, top = 0, rs_sorted[-1][0]
    for ident, length, _ in rs_sorted:
        cumulative += length
        if cumulative >= 0.98 * total:
            top = ident
            break
    for ident, length, genome in rs:
        src = species_of_genome.get(genome)
        if src == species:
            d = 0.0 if genome == rep_of_species.get(species) else div[genome] + div[rep_of_species[species]]
            b = next(label for lo, hi, label in bins if lo < d <= hi or (label == "3-4.5%" and d > hi))
            own_all[b] += length
            for m in margins:
                if ident >= top - m:
                    own[m][b] += length
        else:
            rel_all += length
            for m in margins:
                if ident >= top - m:
                    rel[m] += length

print(f"{sample}: kept share of aligned bases, by margin")
print("reads of".ljust(34) + "".join(f"margin {m:<8}" for m in margins) + "bases")
for _, _, label in bins:
    if own_all[label]:
        print(f"own species, {label} from ref".ljust(34) + "".join(f"{own[m][label] / own_all[label]:<15.3f}" for m in margins)
              + f"{own_all[label]:.0f}")
if rel_all:
    print("relatives (other species)".ljust(34) + "".join(f"{rel[m] / rel_all:<15.3f}" for m in margins) + f"{rel_all:.0f}")
