#!/usr/bin/env python3
"""The long reads a species lost to its congeners in a protal SAM, by why. A read whose best records (its consensus)
are of another taxon, and that has records or ZA candidates of the species, is
  locus       the species has no candidate within 5 edits on any of its genes, or is worse on all of them (the read's
              DNA is nearer the congener: e.g. a locus the species' reference does not share)
  cycle       the species fits at least as well on most genes where it has a candidate (ZA edits <= 0) and better on
              one at least, yet the read went to the congener: the consensus took another taxon (ReadConsensus.h: every
              taxon beaten by one on their common genes, then the one with candidates on most genes)
  mixed       the rest
and per class the reads, records and genes, and the congeners they went to.

usage: lost_reads.py SAM(.zst) GENOME2TIID SPECIES_TAXID
"""
import collections
import subprocess
import sys

sam, genome2tiid, species = sys.argv[1], sys.argv[2], sys.argv[3]
names = {}
with open(genome2tiid) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            names[f[1]] = f[3].split(";")[-1]

reads = collections.defaultdict(list)  # read -> [(taxid, gene, zr, {za taxid: edits})]
sam_text = subprocess.Popen(["zstd", "-dc", sam], stdout=subprocess.PIPE, text=True).stdout
for line in sam_text:
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    flag = int(f[1])
    if flag & (4 | 256):
        continue
    taxid, gene = f[2].split("_", 1)
    tags = dict(t.split(":", 2)[0::2] for t in f[11:])
    za = {}
    if tags.get("ZA", "*") != "*":
        for item in tags["ZA"].split(","):
            t, edits = item.split(":")
            za[t] = int(edits)
    reads[f[0]].append((taxid, gene, tags.get("ZR", "0"), za))

classes = collections.defaultdict(lambda: {"reads": 0, "records": 0, "to": collections.Counter(),
                                           "genes": collections.Counter()})
kept = 0
for name, recs in reads.items():
    taxa = collections.Counter(t for t, _, zr, _ in recs if zr != "2")
    if not taxa:
        continue
    consensus = taxa.most_common(1)[0][0]
    if consensus == species:
        kept += 1
        continue
    # The species' standing on each gene of the read: its edits against the record written (ZA), or absent.
    standing = [za.get(species) for t, _, zr, za in recs if t != species]
    if all(s is None for s in standing) and all(t != species for t, _, _, _ in recs):
        continue  # a read with no sign of the species: not lost by it
    present = [s for s in standing if s is not None]
    if not present or all(s > 0 for s in present):
        kind = "locus"
    elif sum(s <= 0 for s in present) * 2 >= len(present) and any(s < 0 for s in present):
        kind = "cycle"
    else:
        kind = "mixed"
    c = classes[kind]
    c["reads"] += 1
    c["records"] += len(recs)
    c["to"][names.get(consensus, consensus)] += 1
    for t, gene, _, _ in recs:
        c["genes"][gene] += 1

print(f"# {sam}: species {species} {names.get(species, '?')}; reads with its consensus: {kept}")
print("class\treads\trecords\twent_to\ttop_genes")
for kind in ("locus", "cycle", "mixed"):
    c = classes[kind]
    print(f"{kind}\t{c['reads']}\t{c['records']}\t{dict(c['to'].most_common(3))}\t"
          f"{', '.join(g for g, _ in c['genes'].most_common(10))}")
