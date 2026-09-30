#!/usr/bin/env python3
"""test_mini_db.py - checks for simulate_gtdb_release.py, gtdb_to_protal_db.py
and simulate_reads.py.

Runs the scripts into a temporary directory and checks the invariants protal
relies on: reference.map byte offsets, a well-formed internal taxonomy whose
leaves are the reference taxids, byte-identical output for a fixed seed, that
sequence similarity follows the taxonomy, and that simulated reads come from
the right genomes in the right proportions. Does not need a protal binary.

  python3 -m unittest scripts/mini_db/test_mini_db.py
"""

import functools
import gzip
import hashlib
import http.server
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import threading
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = os.path.join(HERE, "simulate_gtdb_release.py")
CONVERT = os.path.join(HERE, "gtdb_to_protal_db.py")
SIMULATE_READS = os.path.join(HERE, "simulate_reads.py")
LINEAGES = os.path.join(HERE, "gtdb_like_lineages.py")
DOWNLOAD = os.path.join(HERE, "..", "download_gtdb.py")

# A stand-in for NCBI's `datasets`: `download genome accession` writes the list into the zip,
# `rehydrate` copies the synthetic release's genomes; accessions in $FAKE_SUPPRESSED fail a request.
FAKE_DATASETS = r'''#!/usr/bin/env python3
import csv, gzip, os, shutil, sys, zipfile
args = sys.argv[1:]
if args[:3] == ["download", "genome", "accession"]:
    accessions = open(args[args.index("--inputfile") + 1]).read().split()
    if set(accessions) & set(os.environ.get("FAKE_SUPPRESSED", "").split(",")):
        sys.exit("Error: some accessions are not valid")
    with zipfile.ZipFile(args[args.index("--filename") + 1], "w") as z:
        z.writestr("ncbi_dataset/fetch.txt", "\n".join(accessions))
elif args[0] == "rehydrate":
    folder = args[args.index("--directory") + 1]
    paths = {r["accession"]: r["fasta_path"] for r in csv.DictReader(open(os.environ["FAKE_TABLE"]), delimiter="\t")}
    for a in open(os.path.join(folder, "ncbi_dataset", "fetch.txt")).read().split():
        os.makedirs(os.path.join(folder, "ncbi_dataset", "data", a), exist_ok=True)
        with gzip.open(paths[a]) as fin, open(os.path.join(folder, "ncbi_dataset", "data", a, a + "_ASM1v1_genomic.fna"), "wb") as fout:
            shutil.copyfileobj(fin, fout)
'''


def gtdb_mirror(release, root):
    """release, a synthetic GTDB release, laid out as GTDB's server has it under root/release226/226.0."""
    base = os.path.join(root, "release226", "226.0")
    os.makedirs(os.path.join(base, "genomic_files_reps"))
    os.makedirs(os.path.join(base, "genomic_files_all"))
    for mset in ("bac120", "ar53"):
        with open(os.path.join(release, f"{mset}_taxonomy_r226.tsv"), "rb") as fin, \
                gzip.open(os.path.join(base, f"{mset}_taxonomy_r226.tsv.gz"), "wb") as fout:
            fout.write(fin.read())
        with open(os.path.join(release, f"{mset}_metadata_r226.tsv.gz"), "rb") as fin, \
                open(os.path.join(base, f"{mset}_metadata_r226.tsv.gz"), "wb") as fout:
            fout.write(fin.read())
        for kind in ("reps", "all"):
            name = f"{mset}_marker_genes_{kind}_r226"
            with tarfile.open(os.path.join(base, f"genomic_files_{kind}", name + ".tar.gz"), "w:gz") as tar:
                tar.add(os.path.join(release, f"genomic_files_{kind}", name), arcname=name)
    with open(os.path.join(base, "VERSION.txt"), "w") as fh:
        fh.write("v226.0\n")
    with open(os.path.join(base, "MD5SUM.txt"), "w") as fh:
        for folder, _dirs, names in os.walk(base):
            for name in sorted(names):
                if name != "MD5SUM.txt":
                    path = os.path.join(folder, name)
                    with open(path, "rb") as f:
                        fh.write(f"{hashlib.md5(f.read()).hexdigest()}  ./{os.path.relpath(path, base)}\n")


def run(*args):
    subprocess.run([sys.executable, *args], check=True, capture_output=True, text=True)


def read_fasta(path):
    records, header = {}, None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()[0]
                records[header] = ""
            else:
                records[header] += line
    return records


def kmers(seq, k=21):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


class MiniDbTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.gtdb = os.path.join(cls.tmp.name, "gtdb")
        cls.db = os.path.join(cls.tmp.name, "db")
        run(SIMULATE, "--outdir", cls.gtdb, "--genome_length", "20000")
        run(CONVERT, "--gtdb", cls.gtdb, "--outdir", cls.db)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def taxonomy(self):
        with open(os.path.join(self.db, "internal_taxonomy.dmp")) as fh:
            next(fh)
            return {int(r[0]): r for r in (l.rstrip("\n").split("\t") for l in fh)}

    def test_reference_map_offsets(self):
        with open(os.path.join(self.db, "reference.fna"), "rb") as fh:
            data = fh.read()
        with open(os.path.join(self.db, "reference.map")) as fh:
            rows = [list(map(int, l.split("\t"))) for l in fh]
        self.assertEqual(len(rows), data.count(b">"))
        for taxid, geneid, start, end in rows:
            header_start = data.rindex(b">", 0, start)
            self.assertEqual(data[header_start:start], f">{taxid}_{geneid}\n".encode())
            self.assertEqual(data[end:end + 1], b"\n")
            self.assertRegex(data[start:end].decode(), r"^[ACGT]+$")

    def test_reference_order(self):
        def map_rows(db):
            with open(os.path.join(db, "reference.map")) as fh:
                return [tuple(map(int, l.split("\t")[:2])) for l in fh]
        rows = map_rows(self.db)  # default: by gene, then taxid
        self.assertEqual(rows, sorted(rows, key=lambda r: (r[1], r[0])))
        genome_db = os.path.join(self.tmp.name, "db_genome_order")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", genome_db, "--order", "genome")
        genome_rows = map_rows(genome_db)
        self.assertEqual(genome_rows, sorted(genome_rows))
        self.assertEqual(sorted(rows), sorted(genome_rows), "the same genes, only the order differs")

    def test_taxonomy_tree(self):
        tax = self.taxonomy()
        roots = [i for i, r in tax.items() if int(r[1]) == i]
        self.assertEqual(len(roots), 1)
        names = [r[3] for r in tax.values()]
        self.assertEqual(len(names), len(set(names)), "names must be unique (string_to_id)")
        for i, r in tax.items():
            if i != roots[0]:
                self.assertIn(int(r[1]), tax)
                self.assertEqual(int(r[5]), int(tax[int(r[1])][5]) + 1, f"level of {r[3]}")
        with open(os.path.join(self.db, "reference.map")) as fh:
            ref_taxids = {int(l.split("\t")[0]) for l in fh}
        parents = {int(r[1]) for i, r in tax.items() if i != roots[0]}
        for t in ref_taxids:
            self.assertEqual(tax[t][4], "species")
            self.assertTrue(tax[t][3].startswith("s__"))
            self.assertRegex(tax[t][6], r"^GCF_999\d{6}\.1$")
            self.assertNotIn(t, parents, "reference taxids must be leaves")

    def test_full_reference_covers_all_genomes(self):
        with open(os.path.join(self.gtdb, "simulation", "genomes.tsv")) as fh:
            n_genomes = sum(1 for _ in fh) - 1
        with open(os.path.join(self.db, "genome2tiid.tsv")) as fh:
            self.assertEqual(sum(1 for _ in fh), n_genomes)
        with open(os.path.join(self.db, "full_reference.fna")) as f, \
             open(os.path.join(self.db, "reference.fna")) as r:
            self.assertGreater(f.read().count(">"), 2 * r.read().count(">"))

    def test_deterministic(self):
        other = os.path.join(self.tmp.name, "gtdb_again")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000")
        genome = os.path.join("genomic_files_reps", "gtdb_genomes_reps_r226", "database", "GCF", "999",
                              "001", "001", "GCF_999001001.1_genomic.fna.gz")
        for f in (os.path.join("genomic_files_all", "bac120_marker_genes_all_r226", "fna", "bac120_TIGR02013.fna"),
                  "bac120_metadata_r226.tsv.gz", genome):
            with open(os.path.join(self.gtdb, f), "rb") as a, open(os.path.join(other, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f} differs between runs with the same seed")

    def test_parallel_conversion(self):
        # The marker files are read by -t workers; the database must not depend on how many.
        for order in ("gene", "genome"):
            one, four = (os.path.join(self.tmp.name, f"db_{order}_t{t}") for t in (1, 4))
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", one, "--order", order, "-t", "1")
            run(CONVERT, "--gtdb", self.gtdb, "--outdir", four, "--order", order, "-t", "4")
            for f in ("reference.fna", "reference.map", "full_reference.fna", "internal_taxonomy.dmp"):
                with open(os.path.join(one, f), "rb") as a, open(os.path.join(four, f), "rb") as b:
                    self.assertEqual(a.read(), b.read(), f"{order} order, {f}")
            self.assertFalse(os.path.exists(os.path.join(four, ".convert_tmp")))

    def test_exclude_species(self):
        excluded = os.path.join(self.tmp.name, "excluded.txt")
        with open(excluded, "w") as fh:
            fh.write("# held out\nMockella beta\n")
        direct, copied = os.path.join(self.tmp.name, "db_direct"), os.path.join(self.tmp.name, "db_copied")
        run(CONVERT, "--gtdb", self.gtdb, "--outdir", direct, "--exclude_species", excluded)
        run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", copied)
        taxids = {r[3]: r[0] for r in self.taxonomy().values() if r[4] == "species"}
        for f in ("reference.fna", "reference.map", "full_reference.fna", "internal_taxonomy.dmp"):
            with open(os.path.join(direct, f), "rb") as a, open(os.path.join(copied, f), "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{f}: --from_db and --gtdb differ")
        with open(os.path.join(self.db, "internal_taxonomy.dmp"), "rb") as a, \
                open(os.path.join(direct, "internal_taxonomy.dmp"), "rb") as b:
            self.assertEqual(a.read(), b.read(), "the taxonomy must keep the species left out")
        with open(os.path.join(direct, "reference.fna")) as fh:
            kept = {line[1:].split("_")[0] for line in fh if line.startswith(">")}
        self.assertNotIn(taxids["s__Mockella beta"], kept)
        self.assertEqual(len(kept), len(taxids) - 1)
        # reference.map still points at each sequence line
        with open(os.path.join(direct, "reference.fna"), "rb") as fna, open(os.path.join(direct, "reference.map")) as fmap:
            data = fna.read()
            for line in fmap:
                tid, gid, start, end = line.split()
                header_end = int(start)
                self.assertEqual(data[data.rfind(b">", 0, header_end):header_end], f">{tid}_{gid}\n".encode())
                self.assertEqual(data[int(end):int(end) + 1], b"\n")
        with self.assertRaises(subprocess.CalledProcessError):
            with open(excluded, "w") as fh:
                fh.write("s__Nonexistent species\n")
            run(CONVERT, "--from_db", self.db, "--exclude_species", excluded, "--outdir", os.path.join(self.tmp.name, "x"))

    def test_download_gtdb(self):
        mirror = os.path.join(self.tmp.name, "mirror")
        gtdb_mirror(self.gtdb, mirror)
        class Quiet(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=mirror))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        datasets = os.path.join(self.tmp.name, "datasets")
        with open(datasets, "w") as fh:
            fh.write(FAKE_DATASETS)
        os.chmod(datasets, 0o755)
        env = dict(os.environ, FAKE_TABLE=os.path.join(self.gtdb, "simulation", "genomes.tsv"),
                   FAKE_SUPPRESSED="GCA_999002003.1")  # the strain of Mockella beta picked, which NCBI no longer has
        out = os.path.join(self.tmp.name, "inputs")
        command = [sys.executable, DOWNLOAD, "-o", out, "--mirror", f"http://127.0.0.1:{server.server_port}",
                   "--datasets", datasets, "--species", "3", "--per_species", "1", "--rep_only_species", "0",
                   "--batch", "4", "-t", "2"]
        first = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)

        release = os.path.join(out, "release")
        for path in ("bac120_taxonomy_r226.tsv.gz", "bac120_metadata_r226.tsv.gz", "VERSION.txt",
                     "genomic_files_reps/bac120_marker_genes_reps_r226", "genomic_files_all/bac120_marker_genes_all_r226"):
            self.assertTrue(os.path.exists(os.path.join(release, path)), path)
        self.assertFalse(os.path.exists(os.path.join(release, "genomic_files_reps", "bac120_marker_genes_reps_r226.tar.gz")),
                         "archives are removed once extracted")
        genomes = sorted(os.listdir(os.path.join(out, "genomes")))
        # 3 species with a strain each and their 3 representatives, less the suppressed strain
        self.assertEqual(len(genomes), 5, genomes)
        self.assertEqual(sum(g.startswith("GCF_") for g in genomes), 3)
        with gzip.open(os.path.join(out, "genomes", genomes[0]), "rt") as fh:
            self.assertTrue(fh.readline().startswith(">"))
        with open(os.path.join(out, "missing.txt")) as fh:
            self.assertEqual(fh.read().split(), ["GCA_999002003.1"])
        with open(os.path.join(out, "simulation_species.txt")) as fh:
            self.assertEqual(len(fh.read().split("\n")) - 1, 3)
        with open(os.path.join(out, "download.json")) as fh:
            state = json.load(fh)
        self.assertEqual(state["release"]["version"], "226.0")
        self.assertEqual((state["genomes"]["delivered"], state["genomes"]["missing"]), (5, 1))

        # A rerun downloads nothing it has, and the release converts as downloaded.
        again = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
        self.assertIn("marker_genes_reps_r226.tar.gz: already there", again.stdout)
        self.assertIn("5 already there, 1 to download", again.stdout)
        run(CONVERT, "--gtdb", release, "--outdir", os.path.join(self.tmp.name, "db_downloaded"))
        with open(os.path.join(self.db, "reference.fna"), "rb") as a, \
                open(os.path.join(self.tmp.name, "db_downloaded", "reference.fna"), "rb") as b:
            self.assertEqual(a.read(), b.read())

        missing = subprocess.run([sys.executable, DOWNLOAD, "-o", out + "_x", "--release", "999", "--mirror",
                                  f"http://127.0.0.1:{server.server_port}", "--no_genomes"], capture_output=True, text=True)
        self.assertNotEqual(missing.returncode, 0)

    def test_lineages_and_clade_holdout(self):
        # The training scripts read lineages from internal_taxonomy.dmp; build_gtdb_database.py holds out whole
        # clades and then single species.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import build_gtdb_database as build
        import collect_training_data as collect
        import lineages
        taxonomy = os.path.join(self.db, "internal_taxonomy.dmp")
        by_id, ids = lineages.from_taxonomy(taxonomy)
        with open(os.path.join(self.db, "genome2tiid.tsv")) as fh:
            for accession, taxid, _rep, lineage in (line.rstrip("\n").split("\t") for line in fh):
                self.assertEqual(by_id[taxid], lineages.from_string(lineage), accession)
        table = os.path.join(self.gtdb, "simulation", "genomes.tsv")
        pool = build.pool_species(table)
        species = {lin["species"]: lin for lin in by_id.values() if "species" in lin}
        families = {}
        for name, lin in species.items():
            families.setdefault(lin["family"], set()).add(name)
        eligible = {f for f, members in families.items() if len(members & pool) >= 2}
        self.assertTrue(eligible, "the synthetic release needs a family with two species")
        chosen = build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3)
        self.assertEqual(chosen, build.choose_holdout(table, taxonomy, 0.5, {"family": 1}, 1.0, 3))
        held_family = {clade for rank, clade in chosen.values() if rank == "family"}
        self.assertEqual(len(held_family), 1)
        family = held_family.pop()
        self.assertIn(family, eligible)
        self.assertEqual({s for s, (rank, _) in chosen.items() if rank == "family"}, families[family])
        rest = sorted(pool - families[family])
        self.assertEqual(sum(rank == "species" for rank, _ in chosen.values()),
                         sum(round(0.5 * len([s for s in rest if species[s]["domain"] == d]))
                             for d in {species[s]["domain"] for s in rest}))
        # a clade holding more than the share, or inside one taken before, is not drawn
        self.assertEqual(build.choose_holdout(table, taxonomy, 0, {"family": 1}, 0.0, 3), {})
        phylum_first = build.choose_holdout(table, taxonomy, 0, {"phylum": 1, "family": 5}, 1.0, 3)
        phylum = next(c for r, c in phylum_first.values() if r == "phylum")
        for s, (rank, clade) in phylum_first.items():
            if rank == "family":
                self.assertNotEqual(species[s]["phylum"], phylum)
        # heldout_species.txt round trip, and the collector's reading of it
        path = os.path.join(self.tmp.name, "heldout.txt")
        with open(path, "w") as fh:
            fh.write("".join(f"{s}\t{r}\t{c}\n" for s, (r, c) in sorted(chosen.items())) + "s__Alone one\n")
        self.assertEqual(build.read_holdout(path), {**chosen, "s__Alone one": ("species", "s__Alone one")})
        self.assertEqual(collect.read_novel(path), build.read_holdout(path))
        self.assertEqual(build.parse_clades("phylum:2,class:4"), {"phylum": 2, "class": 4})
        self.assertEqual(build.parse_clades("none"), {})

    def test_relation_to_novel_species(self):
        # An absent taxon is put down to a species the database lacks when that species is at least as close to
        # it as every present one.
        sys.path.insert(0, os.path.join(HERE, ".."))
        import collect_training_data as collect
        import lineages
        lin = lambda text: lineages.from_string(text)
        taxon = lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G a")
        in_sample = {"s__G b": lin("d__B;p__P;c__C;o__O;f__F;g__G;s__G b"),
                     "s__H c": lin("d__B;p__P;c__C;o__O2;f__F2;g__H;s__H c")}
        self.assertEqual(collect.relation(taxon, in_sample, {}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__G b": ("species", "s__G b")}), ("genus", "species"))
        self.assertEqual(collect.relation(taxon, in_sample, {"s__H c": ("class", "c__C")}), ("genus", ""))
        self.assertEqual(collect.relation(taxon, {"s__H c": in_sample["s__H c"]}, {"s__H c": ("order", "o__O2")}),
                         ("class", "order"))
        self.assertEqual(collect.relation(taxon, {"s__X": lin("d__A;p__Q;s__X")}, {"s__X": ("phylum", "p__Q")}), ("none", ""))

    def test_gtdb_like_lineages(self):
        text = subprocess.run([sys.executable, LINEAGES, "--species", "300", "--archaea", "0.1", "--seed", "3"],
                              check=True, capture_output=True, text=True).stdout
        lineages = [l.split(";") for l in text.splitlines()]
        self.assertEqual(len(lineages), 300)
        self.assertTrue(all([r[:3] for r in l] == ["d__", "p__", "c__", "o__", "f__", "g__", "s__"] for l in lineages))
        self.assertEqual(sum(l[0] == "d__Archaea" for l in lineages), 30)
        self.assertEqual(len({l[6] for l in lineages}), 300)
        parent = {}
        for l in lineages:  # every name has one parent: names are unique across the tree
            for rank in range(1, 7):
                self.assertEqual(parent.setdefault(l[rank], l[rank - 1]), l[rank - 1], l[rank])
        genera = {}
        for l in lineages:
            genera[l[5]] = genera.get(l[5], 0) + 1
        self.assertGreater(sum(1 for n in genera.values() if n == 1), len(genera) / 3, "most genera have one species")
        self.assertGreater(max(genera.values()), 5)
        path = os.path.join(self.tmp.name, "lineages.txt")
        with open(path, "w") as fh:
            fh.write("\n".join(";".join(l) for l in lineages[:20]) + "\n")
        run(SIMULATE, "--outdir", os.path.join(self.tmp.name, "gtdb_like"), "--lineages", path,
            "--genome_length", "5000", "--genomes_per_species", "1")

    def test_divergence_ranges(self):
        def divergence(root):
            with open(os.path.join(root, "simulation", "divergence.tsv")) as fh:
                next(fh)
                return [(float(r[2]), float(r[3])) for r in (l.rstrip("\n").split("\t") for l in fh)]

        self.assertEqual(set(divergence(self.gtdb)), {(0.035, 0.005)})
        other = os.path.join(self.tmp.name, "gtdb_ranges")
        run(SIMULATE, "--outdir", other, "--genome_length", "20000",
            "--strain_divergence", "0.002-0.02", "--species_divergence", "0.01-0.04")
        rates = divergence(other)
        self.assertTrue(all(0.01 <= s <= 0.04 and 0.002 <= g <= 0.02 for s, g in rates))
        self.assertGreater(len({s for s, _ in rates}), 1)
        self.assertEqual(len({g for _, g in rates}), len(rates))

    def simulate_reads(self, prefix, community, pairs=400, error_rate="0"):
        path = os.path.join(self.tmp.name, "community.tsv")
        with open(path, "w") as fh:
            fh.write("# comment\naccession\trelative_abundance\n")
            fh.writelines(f"{acc}\t{ab}\n" for acc, ab in community.items())
        out = os.path.join(self.tmp.name, prefix)
        run(SIMULATE_READS, "--genomes", os.path.join(self.gtdb, "simulation", "genomes.tsv"),
            "--community", path, "--out_prefix", out, "--pairs", str(pairs), "--error_rate", error_rate)
        return out

    def test_reads(self):
        community = {"GCA_999001002.1": 0.75, "GCA_999003003.1": 0.25}
        out = self.simulate_reads("reads", community)
        with open(out + ".truth.tsv") as fh:
            next(fh)
            truth = {f[0]: f for f in (l.rstrip("\n").split("\t") for l in fh)}
        self.assertEqual(set(truth), set(community))
        self.assertEqual(truth["GCA_999003003.1"][1], "s__Fakibacter gamma")
        pairs = {acc: int(f[5]) for acc, f in truth.items()}
        self.assertEqual(sum(pairs.values()), 400)
        # Read pairs follow abundance x genome length.
        weight = {acc: community[acc] * int(f[4]) for acc, f in truth.items()}
        for acc in community:
            self.assertAlmostEqual(pairs[acc], 400 * weight[acc] / sum(weight.values()), delta=1)

        with open(out + "_R1.fq") as f1, open(out + "_R2.fq") as f2:
            r1, r2 = f1.read().split("\n"), f2.read().split("\n")
        self.assertEqual(r1[0::4], r2[0::4], "mates share a read name")
        genomes = {}
        for acc in community:
            fasta = os.path.join(self.gtdb, "simulation", "genomes_nonreps", f"{acc}_genomic.fna.gz")
            genomes[acc] = "|".join(read_fasta(fasta).values())
        comp = str.maketrans("ACGT", "TGCA")
        for name, s1, s2 in zip(r1[0::4], r1[1::4], r2[1::4]):
            acc = name[1:].rsplit("-", 1)[0]
            self.assertIn(acc, community)
            # Error-free: both mates are genome substrings in FR orientation.
            g = genomes[acc]
            self.assertTrue(s1 in g or s1.translate(comp)[::-1] in g, name)
            self.assertTrue(s2 in g or s2.translate(comp)[::-1] in g, name)

        again = self.simulate_reads("reads_again", community)
        for suffix in ("_R1.fq", "_R2.fq", ".truth.tsv"):
            with open(out + suffix, "rb") as a, open(again + suffix, "rb") as b:
                self.assertEqual(a.read(), b.read(), f"{suffix} differs between runs with the same seed")

    def test_similarity_follows_taxonomy(self):
        seqs = read_fasta(os.path.join(self.gtdb, "genomic_files_all", "bac120_marker_genes_all_r226",
                                       "fna", "bac120_TIGR02013.fna"))

        def shared(a, b):
            ka, kb = kmers(seqs[a]), kmers(seqs[b])
            return len(ka & kb) / len(ka | kb)

        strain = shared("RS_GCF_999001001.1", "GB_GCA_999001002.1")      # alpha vs alpha
        congeneric = shared("RS_GCF_999001001.1", "RS_GCF_999002001.1")  # alpha vs beta
        phylum = shared("RS_GCF_999001001.1", "RS_GCF_999003001.1")      # alpha vs gamma
        self.assertGreater(strain, congeneric)
        self.assertGreater(congeneric, phylum)

    def test_no_internal_stops(self):
        for rec, seq in read_fasta(os.path.join(self.gtdb, "genomic_files_reps", "bac120_marker_genes_reps_r226",
                                                "faa", "bac120_TIGR02013.faa")).items():
            self.assertNotIn("*", seq, rec)


if __name__ == "__main__":
    unittest.main()
