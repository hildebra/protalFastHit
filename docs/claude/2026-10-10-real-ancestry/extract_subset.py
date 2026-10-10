"""Stream GTDB r226's packed marker genes and keep the records of the chosen genomes (select_taxa.py's taxa.tsv): the
representatives from the reps archive, every chosen genome from the all-genomes archive, nucleotides only, laid out
as the server's extracted archives (gtdb_to_protal_db.py reads them), with the chosen genomes' taxonomy and metadata
rows. The archives stay packed.

usage: extract_subset.py <GTDB release dir> <taxa.tsv> <out dir> [marker set: bac120, ar53]
"""
import gzip
import os
import sys
import tarfile

REL, TAXA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
R = "r226"
SET = sys.argv[4] if len(sys.argv) > 4 else "bac120"
chosen, reps = set(), set()
with open(TAXA) as f:
    next(f)
    for line in f:
        acc, role = line.split("\t")[:2]
        chosen.add(acc)
        if role == "representative":
            reps.add(acc)

for kind, keep in (("reps", reps), ("all", chosen)):
    src = f"{REL}/genomic_files_{kind}/{SET}_marker_genes_{kind}_{R}.tar.gz"
    base = f"{OUT}/genomic_files_{kind}"
    kept = files = 0
    with tarfile.open(src, "r|gz") as tar:
        for member in tar:
            if not member.isfile() or "/fna/" not in member.name:
                continue
            files += 1
            target = os.path.join(base, member.name)
            os.makedirs(os.path.dirname(target), exist_ok=True)
            data = tar.extractfile(member).read().decode()
            out, take = [], False
            for line in data.splitlines():
                if line.startswith(">"):
                    take = line[1:].split()[0] in keep
                if take:
                    out.append(line)
                    kept += line.startswith(">")
            with open(target, "w") as o:
                o.write("\n".join(out) + ("\n" if out else ""))
    print(f"{kind}: {files} marker files, {kept} records of {len(keep)} genomes")

for name in (f"{SET}_taxonomy_{R}.tsv.gz", f"{SET}_metadata_{R}.tsv.gz"):
    with gzip.open(f"{REL}/{name}", "rt") as f, gzip.open(f"{OUT}/{name}", "wt") as o:
        for i, line in enumerate(f):
            if ("_metadata_" in name and i == 0) or line.split("\t", 1)[0] in chosen:
                o.write(line)
print("taxonomy and metadata rows written")
