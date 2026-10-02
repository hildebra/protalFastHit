#!/usr/bin/env python3
"""foreign_world.py SCRIPTS_DIR DB OUT [PAIRS] - samples in which a database species' gene gets another genome's reads.

On the operon world (~/opw/world/simulation; build ~/opw/b_gn2, whose database lacks the world's 135 unknown species):
PAIRS (default 8) triples of a recipient O, an unknown species' representative, a donor T, a database species'
representative of another phylum, and a marker A of both, such that each end of A in O's genome faces another marker
within 3 kb, and those pairings are unlikely in T's clades by DB's gene_neighbours.tsv (the run's table, species lines
included; the verdict as protal's Table::Assess), while T's own neighbours of A are expected: a gene transferred
recently from T into O, among O's own genes. O' is O with its gene A replaced by T's, 1% of its bases substituted,
in O's orientation. Sample hgt<i> holds T at 10x, O' at 30x and 20 database species at 1-5x; ctl<i> the same without
O'. Reads: paired-end 2x150 (art_illumina HSXt, 350 +- 50 bp), PacBio HiFi (hifi_reads.py, 15 +- 3 kb) and Nanopore
(pbsim3 QSHMM-ONT-HQ, 8 +- 6 kb, 97%) as the collector draws them (collect_training_data.long_read_sample).
Writes OUT/: genomes/, reads/<sample>_{R1,R2,pb,ont}.fq.gz, triples.tsv, samples.tsv (sample, genome, coverage).
"""
import collections
import csv
import gzip
import os
import random
import shutil
import subprocess
import sys

SCRIPTS, DB, OUT = (os.path.expanduser(a) for a in sys.argv[1:4])
PAIRS = int(sys.argv[4]) if len(sys.argv) > 4 else 8
W = os.path.expanduser("~/opw/world/simulation")
GENE_IDS = os.path.expanduser("~/opw/b_gn2/protal_db/gene2geneid.tsv")
UNKNOWN = os.path.expanduser("~/opw/unknown.txt")
ENV = os.path.expanduser("~/micromamba/envs/protal-db-build")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(SCRIPTS, "mini_db"))
import collect_training_data as collect  # noqa: E402
import gene_neighbours as gn  # noqa: E402

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
MAX_GAP = 3000
SEED = 2026


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    out, name, seq = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append([name, "".join(seq)])
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        out.append([name, "".join(seq)])
    return out


def write_fasta(path, contigs):
    with open(path, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def load_table(path):
    counts = collections.defaultdict(dict)
    informative = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("clade"):
                continue
            c, g, e, p, pe, n, inf = (int(x) for x in line.split("\t")[:7])
            counts[(c, g, e)][(p, pe)] = [None] * n
            informative[(c, g, e)] = inf
    return counts, informative


def main():
    rng = random.Random(SEED)
    genomes = {r["accession"]: r for r in csv.DictReader(open(os.path.join(W, "genomes.tsv")), delimiter="\t")}
    markers = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(W, "marker_positions.tsv")), delimiter="\t"):
        markers[r["accession"]].append(r)
    unknown = {line.strip() for line in open(UNKNOWN) if line.strip()}
    gene_id = dict(line.rstrip("\n").split("\t")[:2] for line in open(GENE_IDS) if line.strip())
    nodes, reps = gn.read_taxonomy(os.path.join(DB, "internal_taxonomy.dmp"))
    counts, informative = load_table(os.path.join(DB, "gene_neighbours.tsv"))
    clades = {c for c, _, _ in informative}

    def chain(taxid):
        out, t, seen = [], taxid, set()
        while t in nodes and t not in seen:
            seen.add(t)
            if t in clades:
                out.append(t)
            t = nodes[t][0]
        return out

    def verdict(taxid, gene, end, partner):
        c = chain(taxid)
        share, populated = gn.smoothed_share(counts, informative, c, gene, end, partner)
        if share is None:
            return "unknown"
        if not populated:
            seen = any(partner in counts.get((k, gene, end), {}) for k in c)
            return "expected" if seen else "unknown"
        return "expected" if share >= 0.2 else "unlikely" if share <= 0.05 else "rare"

    def ends(acc, marker):
        """The partners (marker, its facing end) of each end of `marker` in genome acc within MAX_GAP: {5|3: ...}."""
        rows = sorted((r for r in markers[acc]), key=lambda r: (r["contig"], int(r["start"])))
        out = {}
        for i, r in enumerate(rows):
            if r["marker"] != marker:
                continue
            left_end, right_end = (5, 3) if r["strand"] == "+" else (3, 5)
            if i > 0 and rows[i - 1]["contig"] == r["contig"] and int(r["start"]) - int(rows[i - 1]["end"]) <= MAX_GAP:
                p = rows[i - 1]
                out[left_end] = (p["marker"], 3 if p["strand"] == "+" else 5)
            if i + 1 < len(rows) and rows[i + 1]["contig"] == r["contig"] and int(rows[i + 1]["start"]) - int(r["end"]) <= MAX_GAP:
                n = rows[i + 1]
                out[right_end] = (n["marker"], 5 if n["strand"] == "+" else 3)
        return out

    species_of = lambda acc: genomes[acc]["gtdb_taxonomy"].split(";")[-1]
    phylum_of = lambda acc: genomes[acc]["gtdb_taxonomy"].split(";")[1]
    domain_of = lambda acc: genomes[acc]["gtdb_taxonomy"].split(";")[0]
    rep_accs = [a for a, r in genomes.items() if r["gtdb_representative"] == "t"]
    recipients = [a for a in rep_accs if species_of(a) in unknown]
    donors = [a for a in rep_accs if species_of(a) not in unknown and a in reps]
    rng.shuffle(recipients)
    triples, used = [], set()
    for o in recipients:
        if len(triples) == PAIRS:
            break
        o_markers = {r["marker"] for r in markers[o]}
        candidates = [t for t in donors if phylum_of(t) != phylum_of(o) and domain_of(t) == domain_of(o) and t not in used]
        rng.shuffle(candidates)
        found = None
        for t in candidates[:200]:
            taxid = reps[t]
            shared = sorted(o_markers & {r["marker"] for r in markers[t]})
            rng.shuffle(shared)
            for a in shared:
                if a not in gene_id:
                    continue
                o_ends, t_ends = ends(o, a), ends(t, a)
                if len(o_ends) < 2:
                    continue
                g = int(gene_id[a])
                o_verdicts = [verdict(taxid, g, e, (int(gene_id[p]), pe)) for e, (p, pe) in o_ends.items() if p in gene_id]
                t_verdicts = [verdict(taxid, g, e, (int(gene_id[p]), pe)) for e, (p, pe) in t_ends.items() if p in gene_id]
                if len(o_verdicts) == 2 and all(v == "unlikely" for v in o_verdicts) and t_verdicts and \
                        all(v == "expected" for v in t_verdicts):
                    found = (t, a, o_ends, t_ends)
                    break
            if found:
                break
        if not found:
            continue
        t, a, o_ends, t_ends = found
        used.update({o, t})
        triples.append(dict(o=o, o_species=species_of(o), t=t, t_species=species_of(t), t_taxid=reps[t], marker=a,
                            gene=gene_id[a], o_partners=";".join(f"{e}:{p}/{pe}" for e, (p, pe) in sorted(o_ends.items())),
                            t_partners=";".join(f"{e}:{p}/{pe}" for e, (p, pe) in sorted(t_ends.items()))))
        print(f"triple {len(triples)}: {species_of(o)} receives {a} of {species_of(t)}", flush=True)

    shutil.rmtree(OUT, ignore_errors=True)
    for d in ("genomes", "reads", "tmp"):
        os.makedirs(os.path.join(OUT, d))
    samples = []
    background_pool = [d for d in donors if d not in used]
    for i, tr in enumerate(triples, 1):
        # O' : O with T's gene A, 1% substituted, in O's orientation.
        contigs = read_fasta(genomes[tr["o"]]["fasta_path"])
        t_contigs = dict(read_fasta(genomes[tr["t"]]["fasta_path"]))
        mo = next(r for r in markers[tr["o"]] if r["marker"] == tr["marker"])
        mt = next(r for r in markers[tr["t"]] if r["marker"] == tr["marker"])
        gene = t_contigs[mt["contig"]][int(mt["start"]) - 1:int(mt["end"])]
        if mt["strand"] == "-":
            gene = gene.translate(COMPLEMENT)[::-1]
        gene = "".join(rng.choice([b for b in "ACGT" if b != c]) if c in "ACGT" and rng.random() < 0.01 else c for c in gene)
        if mo["strand"] == "-":
            gene = gene.translate(COMPLEMENT)[::-1]
        for contig in contigs:
            if contig[0] == mo["contig"]:
                contig[1] = contig[1][:int(mo["start"]) - 1] + gene + contig[1][int(mo["end"]):]
        o_path = os.path.join(OUT, "genomes", f"recipient{i}.fna")
        write_fasta(o_path, contigs)
        t_path = os.path.join(OUT, "genomes", f"donor{i}.fna")
        write_fasta(t_path, read_fasta(genomes[tr["t"]]["fasta_path"]))
        background = rng.sample(background_pool, 20)
        bg = []
        for b in background:
            path = os.path.join(OUT, "genomes", f"{b}.fna")
            if not os.path.exists(path):
                write_fasta(path, read_fasta(genomes[b]["fasta_path"]))
            bg.append((path, rng.uniform(1, 5)))
        for kind in ("hgt", "ctl"):
            members = [(t_path, 10.0)] + ([(o_path, 30.0)] if kind == "hgt" else []) + bg
            samples.append((f"{kind}{i}", members))
    with open(os.path.join(OUT, "triples.tsv"), "w") as fh:
        cols = list(triples[0])
        fh.write("sample\t" + "\t".join(cols) + "\n")
        for i, tr in enumerate(triples, 1):
            fh.write(f"hgt{i}\t" + "\t".join(str(tr[c]) for c in cols) + "\n")
    with open(os.path.join(OUT, "samples.tsv"), "w") as fh:
        fh.write("sample\tgenome\tcoverage\n")
        for name, members in samples:
            fh.writelines(f"{name}\t{g}\t{c:.3f}\n" for g, c in members)

    art = os.path.join(ENV, "bin", "art_illumina")
    for j, (name, members) in enumerate(samples):
        tmp = os.path.join(OUT, "tmp", name)
        os.makedirs(tmp, exist_ok=True)
        r1, r2 = [], []
        for k, (genome, cov) in enumerate(members):
            prefix = os.path.join(tmp, f"pe{k}_")
            subprocess.run([art, "-ss", "HSXt", "-i", genome, "-p", "-l", "150", "-f", f"{cov:.3f}", "-m", "350", "-s", "50",
                            "-na", "-rs", str(SEED + 100 * j + k), "-o", prefix], check=True, stdout=subprocess.DEVNULL)
            r1.append(prefix + "1.fq")
            r2.append(prefix + "2.fq")
        for parts, mate in ((r1, "R1"), (r2, "R2")):
            with gzip.open(os.path.join(OUT, "reads", f"{name}_{mate}.fq.gz"), "wt", compresslevel=3) as out:
                for part in parts:
                    with open(part) as fh:
                        shutil.copyfileobj(fh, out)
        weighted = [{"genome": os.path.basename(g), "fasta": g,
                     "weight": cov * sum(len(s) for _, s in read_fasta(g))} for g, cov in members]
        for kind, setup in (("pb", "hifi:15000:3000:3"), ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36")):
            error = collect.long_read_sample({
                "sample": name, "tmp": os.path.join(tmp, kind), "out": os.path.join(OUT, "reads", f"{name}_{kind}.fq.gz"),
                "seed": SEED + 10 * j + (kind == "ont"), "setup": collect.parse_long_setup(setup),
                "bases": int(sum(w["weight"] for w in weighted)), "genomes": weighted,
                "pbsim": os.path.join(ENV, "bin", "pbsim"), "model": os.path.join(ENV, "data", "QSHMM-ONT-HQ.model")})
            if error:
                sys.exit(error)
        shutil.rmtree(tmp)
        print(f"{name}: reads written", flush=True)


if __name__ == "__main__":
    main()
