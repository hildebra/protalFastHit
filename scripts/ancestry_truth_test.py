#!/usr/bin/env python3
"""ancestry_truth_test.py - the true-positive test of protal's ancestry sites, strain alleles and polymorphic sites
(docs/claude/2026-10-09-ancestry-true-positive-test/README.md): two worlds grown on trees (simulate_ancestry_world.py), their
databases, reads and protal runs, and the four levels of evaluation.

  run       everything into --outdir, resumable (a step that finished is not run again):
              per world (test, train): the world, the conversion, the gene neighbours, protal --build, --unpack_db;
              the test world also gets a copy of its database with the strain alleles shuffled across species
              (--add_tables) and error-free reads;
              the reads: paired-end (simulate_metagenomes --from_manifest), HiFi (--long_samples), error-free
              (simulate_reads.py --error_rate 0);
              protal runs (placeholder models: the features do not need a trained one): "default" (pe, se, pb, the
              error-free pe), "noallele" (--no_allele_scores: pe, error-free), "shuffled" (the shuffled database: pe,
              error-free); the train world's "default" (pe, se, pb);
              the tables (the training dumps joined, with meta_atp_* columns from the world's truth);
              the presence models (machine_learning_cmdline.py, --test-file the test world's) of five feature sets;
              then evaluate.
  evaluate  the levels on an --outdir that run filled:
              L0 the database's tables against the trees (strain_alleles.tsv against the allele genomes; the ancestry
                 sites against the stem and the representative's lineage);
              L1 every read of the focus genomes, by default against --no_allele_scores: which species it went to,
                 and which the allele scores moved;
              L2 the features: protal's against the oracle's (ancestry_oracle.py, the same counts recomputed from the
                 SAMs), the planted calibration lines, separation of the target's rows with its strain present (Q)
                 from those with a novel species' reads on it (N): AUC overall, within identity bands and by stratum,
                 the gain over identity and depth (grouped by genus), and sensitivity at a fixed specificity by how
                 long ago a strain's lineage diverged and how much of it came from a congener;
              L3 the models' calls on the test world, by feature set, and the controls (no allele scores, shuffled
                 alleles);
              then the pre-registered predictions H1-H10 with their verdicts. Written to OUT/report/: summary.md and a
              TSV per table.

  python3 scripts/ancestry_truth_test.py run --outdir OUT --protal build/protal --simulate build/simulate_metagenomes \\
      --train-python ~/micromamba/envs/protal-db-build/bin/python -t 4 [--quick]
  python3 scripts/ancestry_truth_test.py evaluate --outdir OUT

--train-python (default $PROTAL_TRAIN_PYTHON, else this Python) needs scikit-learn and pandas; without them the models
(L3) are left out. --quick makes two small worlds (simulate_ancestry_world.py --quick), for the tests.
"""

import argparse
import collections
import concurrent.futures
import csv
import gzip
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ancestry_oracle as oracle  # noqa: E402

MINI = os.path.join(HERE, "mini_db")
GENERATOR = os.path.join(MINI, "simulate_ancestry_world.py")
CONVERT = os.path.join(MINI, "gtdb_to_protal_db.py")
GENE_NEIGHBOURS = os.path.join(MINI, "gene_neighbours.py")
SIMULATE_READS = os.path.join(MINI, "simulate_reads.py")
PLACEHOLDERS = os.path.join(HERE, "placeholder_models.py")
TRAINER = os.path.join(HERE, "machine_learning_cmdline.py")

WORLD_SEEDS = {"test": 11, "train": 23}
HIFI_SETUP = "hifi:8000:2000:3"
READ_TYPES = {"pe": "pe", "se": "se", "hifi": "pb", "exact": "pe"}
BASE_SET = ("normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+"
            "gaps+untried")
FEATURE_SETS = {"base": BASE_SET, "ancestry": BASE_SET + "+ancestry", "alleles": BASE_SET + "+alleles+polymorphic",
                "ancestry+alleles": BASE_SET + "+ancestry+alleles", "weights": BASE_SET + "+weights",
                "without weights": BASE_SET + "+ancestry+alleles+polymorphic",
                "default": BASE_SET + "+ancestry+alleles+polymorphic+weights"}
ANCESTRY = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share",
            "ancestry_indel_sites_per_record", "ancestry_indel_congener_share"]
ALLELES = ["allele_explained_share", "allele_identity_gain"]
POLYMORPHIC = ["polymorphic_known_share", "polymorphic_novel_share", "ancestry_fixed_gain", "ancestry_fixed_agreement"]
WEIGHTS = ["conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted", "nonsynonymous_share",
           "nonsynonymous_conserved_rate"]
DIAGNOSTIC = ["allele_copy_share", "allele_sites_per_kb", "ancestry_fixed_share", "column_weight_coverage"]
REFERENCE = ["identity", "top_identity", "lu_per_kb", "uniqueness"]
DERIVED = "fixed_agreement_approx"   # ancestry_agreement + ancestry_fixed_gain: what a model can rebuild
CEILING = "oracle_fixed_agreement"   # the oracle's unweighted agreement at the fixed sites
FPR = 0.05


def log(text):
    print(f"[{time.strftime('%H:%M:%S')}] {text}", flush=True)


def run_command(command, logfile, env=None):
    """Runs command (its output to logfile) in the log's folder: protal writes misc/cpu.tsv where it starts."""
    with open(logfile, "a") as fh:
        fh.write("$ " + " ".join(map(str, command)) + "\n")
        fh.flush()
        result = subprocess.run(list(map(str, command)), stdout=fh, stderr=subprocess.STDOUT, env=env,
                                cwd=os.path.dirname(os.path.abspath(logfile)))
    if result.returncode != 0:
        tail = open(logfile).read().splitlines()[-25:]
        sys.exit(f"failed ({result.returncode}): {' '.join(map(str, command))}\n  " + "\n  ".join(tail))


class Paths:
    def __init__(self, out, name):
        self.root = os.path.join(out, name)
        self.world = os.path.join(self.root, "world")
        self.db = os.path.join(self.root, "db")
        self.unpacked = os.path.join(self.db, "unpacked")
        self.shuffled = os.path.join(self.root, "db_shuffled")
        self.reads = os.path.join(self.root, "reads")
        self.runs = os.path.join(self.root, "runs")
        self.tables = os.path.join(self.root, "tables")
        self.logs = os.path.join(self.root, "logs")
        self.done = os.path.join(self.root, ".done")

    def table(self, run, kind):
        return os.path.join(self.tables, f"{run}_{kind}.tsv")


def step(paths, name, fn):
    """Runs fn() unless the step finished before (a marker in .done)."""
    os.makedirs(paths.done, exist_ok=True)
    os.makedirs(paths.logs, exist_ok=True)
    marker = os.path.join(paths.done, name)
    if os.path.exists(marker):
        return
    started = time.time()
    log(f"{os.path.basename(paths.root)}: {name}")
    fn()
    open(marker, "w").write(f"{time.time() - started:.1f}\n")
    log(f"{os.path.basename(paths.root)}: {name} done in {time.time() - started:.0f} s")


# ---- run ----------------------------------------------------------------------------------------------------------

def shuffle_alleles(source, target, seed):
    """strain_alleles.tsv with each gene's rows' alleles rotated across the species that have the gene: every copy
    keeps a table entry of the same size class, none its own (the genes have one length in these worlds)."""
    rows = collections.defaultdict(list)
    header = None
    with open(source) as fh:
        for line in fh:
            if line.startswith(("#", "taxid")):
                header = header or line
                continue
            taxid, gene, alleles = line.rstrip("\n").split("\t")
            rows[int(gene)].append((int(taxid), alleles))
    rng = random.Random(seed)
    out = []
    for gene, items in rows.items():
        alleles = [a for _, a in items]
        if len(items) > 1:
            k = rng.randrange(1, len(items))
            alleles = alleles[k:] + alleles[:k]
        out += [(taxid, gene, a) for (taxid, _), a in zip(items, alleles)]
    out.sort()
    with open(target, "w", newline="\n") as fh:
        fh.write(header or "taxid\tgene\talleles\n")
        for taxid, gene, alleles in out:
            fh.write(f"{taxid}\t{gene}\t{alleles}\n")


def write_map(path, rows):
    columns = ["SAMPLEID", "FIRST", "SECOND", "SAM", "PREFIX", "PROFILE", "PROFILE_TRUTH", "READ_TYPE", "UNMAPPED_READS"]
    with open(path, "w", newline="\n") as fh:
        fh.write(f"#OUTPUT_DIR\t{os.path.dirname(path)}\n#" + "\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "count" if c == "UNMAPPED_READS" else "-") for c in columns) + "\n")


def read_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


def sample_reads(paths, kind):
    """[(sample, first, second, truth)] of a read set."""
    samples = os.path.join(paths.world, "samples")
    out = []
    if kind in ("pe", "se"):
        for row in read_tsv(os.path.join(paths.reads, "pe", "manifest.tsv")):
            if out and out[-1][0] == row["sample"]:
                continue
            truth = os.path.join(samples, "pe", "truth", row["sample"] + ".tsv")
            out.append((row["sample"], row["fastq_r1"], row["fastq_r2"], truth))
    elif kind == "hifi":
        for row in read_tsv(os.path.join(samples, "hifi", "long_samples.tsv")):
            out.append((row["sample"], row["out"], "-", os.path.join(samples, "hifi", "truth", row["sample"] + ".tsv")))
    elif kind == "exact":
        for row in read_tsv(os.path.join(samples, "exact", "samples.tsv")):
            prefix = os.path.join(paths.reads, "exact", row["sample"])
            out.append((row["sample"], prefix + "_R1.fq", prefix + "_R2.fq", row["truth"]))
    return out


def profile_rows(paths, run, kinds):
    rows = []
    folder = os.path.join(paths.runs, run)
    for kind in kinds:
        for sample, first, second, truth in sample_reads(paths, kind):
            name = f"{sample}_se" if kind == "se" else sample
            rows.append({"SAMPLEID": name, "FIRST": first, "SECOND": "-" if kind in ("se", "hifi") else second,
                         "SAM": os.path.join(folder, "alignments", name + ".sam.gz"),
                         "PREFIX": os.path.join(folder, "misc", name),
                         "PROFILE": os.path.join(folder, "profiles", name + ".profile"),
                         "PROFILE_TRUTH": truth, "READ_TYPE": READ_TYPES[kind]})
    return rows


def run_world(opts, name):
    paths = Paths(opts.outdir, name)
    test = name == "test"
    os.makedirs(paths.root, exist_ok=True)
    t = str(opts.threads)

    def world():
        command = [sys.executable, GENERATOR, "--outdir", paths.world, "--seed", WORLD_SEEDS[name]]
        if opts.quick:
            command.append("--quick")
        if not test:
            command += ["--exact_samples", "0", "--no_decoy"]  # a decoy's label is noise to a model
        command += opts.world_args
        shutil.rmtree(paths.world, ignore_errors=True)
        run_command(command, os.path.join(paths.logs, "world.log"))

    def convert():
        shutil.rmtree(paths.db, ignore_errors=True)
        run_command([sys.executable, CONVERT, "--gtdb", paths.world, "--outdir", paths.db, "-t", t],
                    os.path.join(paths.logs, "convert.log"))

    def neighbours():
        run_command([sys.executable, GENE_NEIGHBOURS, "--db", paths.db, "--genome_table",
                     os.path.join(paths.world, "simulation", "genomes.tsv")], os.path.join(paths.logs, "neighbours.log"))

    def build():
        run_command([opts.protal, "--build", "--no_profile", "-t", t, "--db", paths.db, "--reference",
                     os.path.join(paths.db, "reference.fna"), "--full_reference",
                     os.path.join(paths.db, "full_reference.fna"), "--strain_alleles", opts.strain_alleles,
                     "--allele_genome_share", "1"], os.path.join(paths.logs, "build.log"))

    def unpack():
        shutil.rmtree(paths.unpacked, ignore_errors=True)
        run_command([opts.protal, "--unpack_db", "--db", os.path.join(paths.db, "database.protal"), "--unpack_dir",
                     paths.unpacked, "-t", t], os.path.join(paths.logs, "unpack.log"))

    def shuffled():
        shutil.rmtree(paths.shuffled, ignore_errors=True)
        os.makedirs(os.path.join(paths.shuffled, "tables"))
        shutil.copy(os.path.join(paths.db, "database.protal"), paths.shuffled)
        table = os.path.join(paths.shuffled, "tables", "strain_alleles.tsv")
        shuffle_alleles(os.path.join(paths.unpacked, "strain_alleles.tsv"), table, WORLD_SEEDS[name])
        run_command([opts.protal, "--add_tables", table, "--db", paths.shuffled, "-t", t],
                    os.path.join(paths.logs, "shuffled.log"))

    def reads_pe():
        out = os.path.join(paths.reads, "pe")
        shutil.rmtree(out, ignore_errors=True)
        run_command([opts.simulate, "--from_manifest", os.path.join(paths.world, "samples", "pe", "manifest.tsv"),
                     "--output_dir", out, "-t", t, "--seed", WORLD_SEEDS[name]], os.path.join(paths.logs, "reads_pe.log"))

    def reads_hifi():
        samples = os.path.join(paths.world, "samples", "hifi")
        run_command([opts.simulate, "--long_samples", os.path.join(samples, "long_samples.tsv"), "--long_genomes",
                     os.path.join(samples, "long_genomes.tsv"), "--long_setup", HIFI_SETUP, "-t", t],
                    os.path.join(paths.logs, "reads_hifi.log"))

    def reads_exact():
        out = os.path.join(paths.reads, "exact")
        os.makedirs(out, exist_ok=True)
        rows = read_tsv(os.path.join(paths.world, "samples", "exact", "samples.tsv"))
        genomes = os.path.join(paths.world, "simulation", "all_genomes.tsv")
        logfile = os.path.join(paths.logs, "reads_exact.log")

        def one(row):
            command = [sys.executable, SIMULATE_READS, "--genomes", genomes, "--community", row["community"],
                       "--out_prefix", os.path.join(out, row["sample"]), "--pairs", row["pairs"], "--seed", row["seed"],
                       "--error_rate", "0"]
            result = subprocess.run(command, capture_output=True, text=True)
            return command, result

        with concurrent.futures.ThreadPoolExecutor(max(1, opts.threads)) as pool:
            for command, result in pool.map(one, rows):
                with open(logfile, "a") as fh:
                    fh.write("$ " + " ".join(command) + "\n" + result.stdout + result.stderr)
                if result.returncode != 0:
                    sys.exit(f"simulate_reads.py failed: {result.stderr[-2000:]}")

    def profile(run, db, kinds, extra=()):
        def go():
            folder = os.path.join(paths.runs, run)
            shutil.rmtree(folder, ignore_errors=True)
            for sub in ("alignments", "profiles", "misc"):
                os.makedirs(os.path.join(folder, sub))
            map_file = os.path.join(folder, "samples.map")
            write_map(map_file, profile_rows(paths, run, kinds))
            models = os.path.join(opts.outdir, "placeholder_models")
            run_command([opts.protal, "--db", db, "--map", map_file, "-t", t, "--no_strains", "--no_qcmsa",
                         "--model", os.path.join(models, "model_pe.xml"), "--model_se",
                         os.path.join(models, "model_se.xml"), "--model_pb", os.path.join(models, "model_PB.xml"),
                         *extra], os.path.join(folder, "protal.log"))
        return go

    def tables():
        os.makedirs(paths.tables, exist_ok=True)
        data = WorldData(paths)
        for run, kinds in runs_of(test).items():
            for kind in kinds:
                write_table(data, paths, run, kind)

    step(paths, "world", world)
    step(paths, "convert", convert)
    step(paths, "gene_neighbours", neighbours)
    step(paths, "build", build)
    step(paths, "unpack", unpack)
    step(paths, "reads_pe", reads_pe)
    step(paths, "reads_hifi", reads_hifi)
    if test:
        step(paths, "reads_exact", reads_exact)
        step(paths, "shuffled_db", shuffled)
    for run, kinds in runs_of(test).items():
        db = paths.shuffled if run == "shuffled" else paths.db
        extra = ["--no_allele_scores"] if run == "noallele" else []
        step(paths, f"profile_{run}", profile(run, db, kinds, extra))
    step(paths, "tables", tables)
    return paths


def runs_of(test):
    if test:
        return {"default": ["pe", "se", "hifi", "exact"], "noallele": ["pe", "exact"], "shuffled": ["pe", "exact"]}
    return {"default": ["pe", "se", "hifi"]}


def have_trainer(python):
    result = subprocess.run([python, "-c", "import sklearn, pandas"], capture_output=True)
    return result.returncode == 0


def train_models(opts):
    test, train = Paths(opts.outdir, "test"), Paths(opts.outdir, "train")
    folder = os.path.join(opts.outdir, "models")
    os.makedirs(folder, exist_ok=True)
    if not have_trainer(opts.train_python):
        log(f"no scikit-learn/pandas in {opts.train_python}: the models (L3) are left out")
        return
    sets = ["base", "default"] if opts.quick else opts.model_sets
    jobs = []
    for kind in (["pe"] if opts.quick else opts.model_types):
        # every set for paired-end reads; the other read types base and default unless --model-sets-all
        for name in sets if kind == "pe" or opts.model_sets_all else [s for s in sets if s in ("base", "default")]:
            jobs.append((kind, name, "default"))
    jobs += [("pe", "default", "noallele"), ("pe", "default", "shuffled")]
    for kind, name, run in jobs:
        prefix = os.path.join(folder, f"{kind}_{name}" + ("" if run == "default" else f"_on_{run}"))
        if os.path.exists(prefix + ".calls.tsv.gz"):
            continue
        log(f"model {kind} {name} (tested on {run})")
        run_command([opts.train_python, TRAINER, "--truth-file", train.table("default", kind), "--output-prefix",
                     prefix, "--features", FEATURE_SETS[name], "--evaluation", "basic", "--folds", "3", "--seed", "1",
                     "--threads", str(opts.threads), "--test-file", test.table(run, kind)],
                    os.path.join(folder, "training.log"))


def command_run(opts):
    os.makedirs(opts.outdir, exist_ok=True)
    models = os.path.join(opts.outdir, "placeholder_models")
    if not os.path.exists(os.path.join(models, "model_PB.xml")):
        run_command([sys.executable, PLACEHOLDERS, "-o", models, "--read_types", "pe,se,pb"],
                    os.path.join(opts.outdir, "placeholders.log"))
    for name in ("test", "train"):
        run_world(opts, name)
    train_models(opts)
    command_evaluate(opts)


# ---- the worlds' truth --------------------------------------------------------------------------------------------

class WorldData:
    """A world's truth and its database's names: genera, roles, the samples' design, taxids."""

    def __init__(self, paths):
        sim = os.path.join(paths.world, "simulation")
        self.paths = paths
        self.genera = {int(r["genus"]): r for r in read_tsv(os.path.join(sim, "genera.tsv"))}
        self.roles = {r["accession"]: r for r in read_tsv(os.path.join(sim, "roles.tsv"))}
        self.genomes = {r["accession"]: r for r in read_tsv(os.path.join(sim, "all_genomes.tsv"))}
        self.design = collections.defaultdict(dict)  # (set, sample) -> {genus: row}
        self.members = collections.defaultdict(list)  # (set, sample) -> [genome accession] in order
        for r in read_tsv(os.path.join(paths.world, "samples", "design.tsv")):
            self.members[(r["set"], r["sample"])].append(r["genome"])
            if r["genus"]:
                self.design[(r["set"], r["sample"])][int(r["genus"])] = r
        self.taxid_of, self.rep_of = {}, {}
        with open(os.path.join(paths.db, "genome2tiid.tsv")) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 3:
                    self.taxid_of[row[0]] = int(row[1])
                    self.rep_of[int(row[1])] = row[2]
        self.gene_ids = {}
        with open(os.path.join(paths.db, "gene2geneid.tsv")) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2 and row[1].isdigit():
                    self.gene_ids[row[0]] = int(row[1])
        self.target_taxid, self.genus_of_taxid, self.sister_taxid = {}, {}, {}
        for g, row in self.genera.items():
            self.target_taxid[g] = self.taxid_of[row["target_rep"]]
            self.genus_of_taxid[self.taxid_of[row["target_rep"]]] = g
            for acc in row["congener_reps"].split(","):
                self.genus_of_taxid[self.taxid_of[acc]] = g
            self.sister_taxid[g] = self.taxid_of[row["sister"]]
        self.targets = set(self.target_taxid.values())

    def origin(self, kind, sample, qname):
        """The accession a read came from: simulate_metagenomes' <contig>-<n>, simulate_reads.py's <accession>-<n>,
        the long reads' g<i>x_<n> (the sample's i-th genome)."""
        if kind == "hifi":
            index = int(qname[1:qname.index("x_")])
            return self.members[("hifi", sample)][index]
        name = qname.rsplit("-", 1)[0]
        return name.rsplit("_contig", 1)[0]


META = ["meta_sample", "meta_read_pairs", "meta_atp_set", "meta_atp_genus", "meta_atp_target", "meta_atp_focus",
        "meta_atp_role", "meta_atp_class", "meta_atp_depth", "meta_atp_b", "meta_atp_congeners", "meta_atp_twin",
        "meta_atp_alleles", "meta_atp_regime", "meta_atp_c_S", "meta_atp_L", "meta_atp_distance", "meta_atp_mrca_rep",
        "meta_atp_mrca_allele", "meta_atp_genes_imported", "meta_atp_genes_congeneric", "meta_atp_lineage_covered",
        "meta_atp_truth_agreement", "meta_atp_truth_fixed"]


def row_meta(data, kind, sample, taxid):
    """The meta_atp_* columns of a taxon's row in a sample."""
    design_set = "pe" if kind == "se" else kind
    base = sample[:-3] if kind == "se" else sample
    meta = dict.fromkeys(META, "")
    meta["meta_sample"], meta["meta_atp_set"] = sample, kind
    genus = data.genus_of_taxid.get(taxid)
    if genus is None:
        return meta
    g = data.genera[genus]
    meta.update({"meta_atp_genus": genus, "meta_atp_target": int(taxid == data.target_taxid[genus]),
                 "meta_atp_congeners": g["congeners"], "meta_atp_twin": g["twin"],
                 "meta_atp_alleles": g["allele_genomes"], "meta_atp_regime": g["regime"], "meta_atp_c_S": g["c_S"],
                 "meta_atp_L": g["L"]})
    focus = data.design.get((design_set, base), {}).get(genus)
    if focus is None:
        return meta
    r = data.roles.get(focus["genome"], {})
    meta.update({"meta_atp_focus": focus["genome"], "meta_atp_role": focus["role"],
                 "meta_atp_class": focus["class"], "meta_atp_depth": focus["depth"], "meta_atp_b": r.get("b", ""),
                 "meta_atp_distance": r.get("distance", ""), "meta_atp_mrca_rep": r.get("mrca_rep", ""),
                 "meta_atp_mrca_allele": r.get("mrca_nearest_allele", ""),
                 "meta_atp_genes_imported": r.get("genes_imported", ""),
                 "meta_atp_genes_congeneric": r.get("genes_congeneric", "")})
    if r:
        cons, fixed = int(r["cons_sites"]), int(r["fixed_sites"])
        meta["meta_atp_truth_agreement"] = f"{int(r['cons_agree']) / cons:.6g}" if cons else ""
        meta["meta_atp_truth_fixed"] = f"{int(r['fixed_agree']) / fixed:.6g}" if fixed else ""
        # of the representative's own derived sites the genome lacks, the share the allele genomes make polymorphic
        lineage = int(r["at_rep_lineage"])
        meta["meta_atp_lineage_covered"] = f"{int(r['at_rep_lineage_poly']) / lineage:.6g}" if lineage else ""
    return meta


def write_table(data, paths, run, kind):
    """The run's training dumps of a read set joined, with the meta columns."""
    rows, header = [], None
    for sample, _, _, _ in sample_reads(paths, kind):
        name = f"{sample}_se" if kind == "se" else sample
        dump = os.path.join(paths.runs, run, "profiles", name + ".profile.truth_annotated")
        if not os.path.exists(dump):
            sys.exit(f"{dump} is missing: did protal profile {name}?")
        with open(dump) as fh:
            reader = csv.reader(fh, delimiter="\t")
            columns = next(reader)
            header = header or columns
            if columns != header:
                sys.exit(f"{dump}: other columns than the run's first dump")
            for values in reader:
                record = dict(zip(columns, values))
                meta = row_meta(data, kind, name, int(record["taxon"]))
                meta["meta_read_pairs"] = 0  # the trainer's depth table: one design point
                rows.append([meta[c] for c in META] + values)
    with open(paths.table(run, kind), "w", newline="\n") as fh:
        fh.write("\t".join(META + (header or [])) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


# ---- statistics ---------------------------------------------------------------------------------------------------

def ranks(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    out = np.empty(len(values))
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        out[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return out


def auc(pos, neg):
    """P(a positive scores above a negative), ties half: nan without both."""
    pos, neg = np.asarray(pos, dtype=float), np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = ranks(np.concatenate([pos, neg]))
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def stratified_auc(pos, neg, pos_strata, neg_strata):
    """The AUC within strata, weighted by their pairs."""
    total = weight = 0.0
    for s in set(pos_strata) & set(neg_strata):
        p = [v for v, t in zip(pos, pos_strata) if t == s]
        n = [v for v, t in zip(neg, neg_strata) if t == s]
        if p and n:
            total += auc(p, n) * len(p) * len(n)
            weight += len(p) * len(n)
    return total / weight if weight else float("nan")


def logistic(X, y, l2=1e-2, iterations=60):
    """Ridge logistic regression (IRLS) on standardised columns -> a function scoring rows."""
    mean, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1
    Z = np.column_stack([np.ones(len(X)), (X - mean) / sd])
    w = np.zeros(Z.shape[1])
    penalty = l2 * np.eye(Z.shape[1])
    penalty[0, 0] = 0
    for _ in range(iterations):
        p = 1 / (1 + np.exp(-np.clip(Z @ w, -30, 30)))
        gradient = Z.T @ (p - y) + penalty @ w
        hessian = (Z * (p * (1 - p))[:, None]).T @ Z + penalty + 1e-9 * np.eye(Z.shape[1])
        step_ = np.linalg.solve(hessian, gradient)
        w -= step_
        if np.abs(step_).max() < 1e-8:
            break
    return lambda A: np.column_stack([np.ones(len(A)), (A - mean) / sd]) @ w


def grouped_cv_scores(X, y, groups, folds=5):
    """Scores of each row by a logistic model fitted without its group's rows (genera in `folds` folds)."""
    unique = sorted(set(groups))
    fold_of = {g: i % folds for i, g in enumerate(unique)}
    fold = np.array([fold_of[g] for g in groups])
    out = np.zeros(len(y))
    for k in range(min(folds, len(unique))):
        train, test = fold != k, fold == k
        if not test.any() or len(set(y[train])) < 2:
            continue
        out[test] = logistic(X[train], y[train])(X[test])
    return out


def bootstrap_difference(score_a, score_b, y, groups, n=200, seed=1):
    """AUC(b) - AUC(a) and its 95% interval over resamples of the groups (genera)."""
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups)
    unique = np.unique(groups)
    index = {g: np.flatnonzero(groups == g) for g in unique}
    full = auc(score_b[y == 1], score_b[y == 0]) - auc(score_a[y == 1], score_a[y == 0])
    diffs = []
    for _ in range(n):
        rows = np.concatenate([index[g] for g in rng.choice(unique, len(unique))])
        yy = y[rows]
        if yy.min() == yy.max():
            continue
        a, b = score_a[rows], score_b[rows]
        diffs.append(auc(b[yy == 1], b[yy == 0]) - auc(a[yy == 1], a[yy == 0]))
    if not diffs:
        return full, float("nan"), float("nan")
    return full, float(np.quantile(diffs, 0.025)), float(np.quantile(diffs, 0.975))


def tpr_at_fpr(pos, neg, fpr=FPR):
    """The share of positives above the threshold that (1 - fpr) of the negatives do not exceed, and the threshold."""
    if not len(pos) or not len(neg):
        return float("nan"), float("nan")
    threshold = float(np.quantile(np.asarray(neg, dtype=float), 1 - fpr))
    return float(np.mean(np.asarray(pos, dtype=float) > threshold)), threshold


def fnum(x, digits=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "-"
    return f"{x:.{digits}f}" if isinstance(x, float) else str(x)


def value(row, column):
    try:
        return float(row[column])
    except (KeyError, ValueError, TypeError):
        return float("nan")


# ---- evaluate -----------------------------------------------------------------------------------------------------

class Report:
    def __init__(self, folder):
        self.folder = folder
        os.makedirs(folder, exist_ok=True)
        self.lines = []
        self.verdicts = []

    def section(self, title):
        self.lines += ["", f"## {title}", ""]

    def text(self, text):
        self.lines.append(text)

    def table(self, name, rows, columns=None, show=40):
        """rows (dicts) to OUT/report/<name>.tsv and, the first `show`, into the summary."""
        if not rows:
            self.lines.append(f"({name}: no rows)")
            return
        columns = columns or list(rows[0])
        with open(os.path.join(self.folder, name + ".tsv"), "w", newline="\n") as fh:
            fh.write("\t".join(columns) + "\n")
            for r in rows:
                fh.write("\t".join(fnum(r.get(c)) if isinstance(r.get(c), float) else str(r.get(c, ""))
                                   for c in columns) + "\n")
        self.lines.append("| " + " | ".join(columns) + " |")
        self.lines.append("|" + "---|" * len(columns))
        for r in rows[:show]:
            self.lines.append("| " + " | ".join(fnum(r.get(c)) if isinstance(r.get(c), float) else str(r.get(c, ""))
                                                for c in columns) + " |")
        if len(rows) > show:
            self.lines.append(f"({len(rows) - show} more rows in {name}.tsv)")

    def verdict(self, key, prediction, outcome, status):
        self.verdicts.append({"hypothesis": key, "prediction": prediction, "outcome": outcome, "verdict": status})

    def write(self, header):
        with open(os.path.join(self.folder, "summary.md"), "w", newline="\n") as fh:
            fh.write(header + "\n")
            fh.write("\n## Verdicts\n\n| hypothesis | prediction | outcome | verdict |\n|---|---|---|---|\n")
            for v in self.verdicts:
                fh.write(f"| {v['hypothesis']} | {v['prediction']} | {v['outcome']} | {v['verdict']} |\n")
            fh.write("\n".join(self.lines) + "\n")
        with open(os.path.join(self.folder, "verdicts.tsv"), "w", newline="\n") as fh:
            fh.write("hypothesis\tprediction\toutcome\tverdict\n")
            for v in self.verdicts:
                fh.write(f"{v['hypothesis']}\t{v['prediction']}\t{v['outcome']}\t{v['verdict']}\n")


def load_table(path):
    rows = read_tsv(path)
    for r in rows:
        agreement, gain = value(r, "ancestry_agreement"), value(r, "ancestry_fixed_gain")
        r[DERIVED] = str(agreement + gain) if agreement >= 0 and gain > -1 else str(agreement)
    return rows


def target_rows(rows):
    """The target species' rows of the focus genomes' classes Q, N and D (the decoy, compared with Q_deep alone; P,
    planted and leaky, apart)."""
    return [r for r in rows if r.get("meta_atp_target") == "1" and r.get("meta_atp_class") in ("Q", "N", "D")]


def evaluate_l0(report, data):
    report.section("L0: the database's tables against the trees")
    # The strain alleles: each stored allele against the allele genomes' true differences in its range.
    truth = collections.defaultdict(dict)  # (genus, marker) -> {kind/genome: {pos: base}}
    with gzip.open(os.path.join(data.paths.world, "simulation", "truth_sites.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            genus, marker, kind, genome, sites = line.rstrip("\n").split("\t")
            key = genome if kind == "allele" else kind
            truth[(int(genus), marker)][key] = {int(s[:-1]): s[-1] for s in sites.split(",") if s}
    stored = oracle.read_alleles(os.path.join(data.paths.unpacked, "strain_alleles.tsv"))
    rows = []
    stored_genome = collections.defaultdict(set)  # (genus, marker) -> allele genomes protal stored
    for genus, g in sorted(data.genera.items()):
        taxid = data.target_taxid[genus]
        copies = with_alleles = alleles = exact = 0
        for marker, gid in data.gene_ids.items():
            sites = truth.get((genus, marker))
            if sites is None:
                continue
            copies += 1
            table = stored.get((taxid, gid), [])
            with_alleles += bool(table)
            for allele in table:
                alleles += 1
                edits = {e.pos: "ACGT"[e.base] for e in allele.edits if e.kind == oracle.SUB}
                for genome, diffs in sites.items():
                    if not genome.startswith("GC"):
                        continue
                    inside = {p: b for p, b in diffs.items() if allele.begin <= p < allele.end}
                    if inside == edits:
                        exact += 1
                        stored_genome[(genus, marker)].add(genome)
                        break
        made = int(g["alleles_made"])
        rows.append({"genus": genus, "congeners": g["congeners"], "twin": g["twin"], "allele_genomes": made,
                     "regime": g["regime"], "copies": copies, "copies_with_alleles": with_alleles,
                     "alleles_per_copy": alleles / copies if copies else float("nan"),
                     "share_exact": exact / alleles if alleles else float("nan")})
    report.text("Per test genus: copies of the target's genes with stored alleles, alleles per copy, and the share of "
                "stored alleles whose edits are exactly an allele genome's differences in the allele's range.")
    report.table("l0_alleles", rows)
    with_genomes = [r for r in rows if int(r["allele_genomes"]) > 0]
    exact = [r["share_exact"] for r in with_genomes if not math.isnan(r["share_exact"])]
    # Q strains whose nearest allele genome protal stored, per gene.
    sister_rows = []
    for acc, r in data.roles.items():
        if r["class"] != "Q" or not r["nearest_allele"]:
            continue
        genus = int(r["genus"])
        markers = [m for m in data.gene_ids if (genus, m) in truth]
        share = sum(r["nearest_allele"] in stored_genome[(genus, m)] for m in markers) / max(1, len(markers))
        sister_rows.append({"genome": acc, "genus": genus, "role": r["role"],
                            "allele_genomes": data.genera[genus]["alleles_made"],
                            "twin": data.genera[genus]["twin"], "mrca_nearest_allele": r["mrca_nearest_allele"],
                            "share_genes_nearest_stored": share})
    report.text("")
    report.text("Per strain with an allele genome in its species: the share of genes in which its nearest allele "
                "genome is among the stored alleles (K = 4, chosen for the sample's coverage since 2026-10-09 evening, "
                "farthest first before; the build rejects alleles as far as the nearest congener's copy).")
    report.table("l0_nearest_allele_stored", sorted(sister_rows, key=lambda x: (x["genus"], x["role"])))
    build_log = os.path.join(data.paths.logs, "build.log")
    if os.path.exists(build_log):
        lines = [line.strip() for line in open(build_log) if "Strain alleles:" in line]
        if lines:
            report.text("")
            report.text("The build's line: " + lines[-1])
    status = "PASS" if exact and min(exact) >= 0.95 else "FAIL" if exact else "N/A"
    report.verdict("H1a table", "stored alleles are the allele genomes' differences (>= 0.95 exact)",
                   f"min share exact {fnum(min(exact) if exact else float('nan'))}", status)

    # The ancestry sites: the oracle's (protal's rule on the database's references) against the history.
    db = oracle.Database(data.paths.unpacked)
    site_rows = []
    for genus, g in sorted(data.genera.items()):
        taxid = data.target_taxid[genus]
        n = in_history = stem = lineage = stem_total = lineage_total = 0
        why = collections.Counter()  # history sites missed and other sites taken, by why
        for marker, gid in data.gene_ids.items():
            sites = truth.get((genus, marker))
            if sites is None:
                continue
            s = db.sites(taxid, gid)
            positions = set(s.positions)
            history = set(sites.get("stem", {})) | set(sites.get("rep_lineage", {}))
            n += len(positions)
            in_history += len(positions & history)
            stem += len(positions & set(sites.get("stem", {})))
            lineage += len(positions & set(sites.get("rep_lineage", {})))
            stem_total += len(sites.get("stem", {}))
            lineage_total += len(sites.get("rep_lineage", {}))
            compared = getattr(s, "compared_at", None) or []
            for p in history - positions:
                k = compared[p] if p < len(compared) else 0
                why["missed: not compared" if k == 0 else "missed: < 3 congeners, the nearest agrees"
                    if k < oracle.MIN_CONGENERS else "missed: a congener carries another base"] += 1
            for p in positions - history:
                k = compared[p] if p < len(compared) else 0
                why["other: the nearest's own (< 3 congeners)" if k < oracle.MIN_CONGENERS else
                    "other: consensus"] += 1
        site_rows.append({"genus": genus, "congeners": g["congeners"], "twin": g["twin"], "regime": g["regime"],
                          "sites": n, "share_in_history": in_history / n if n else float("nan"),
                          "stem_recall": stem / stem_total if stem_total else float("nan"),
                          "rep_lineage_recall": lineage / lineage_total if lineage_total else float("nan"),
                          **{k: v for k, v in sorted(why.items())}})
    report.text("")
    report.text("The ancestry sites (the oracle's port of protal's rule, from the database's references; the oracle "
                "equals protal, L2a) against the trees: the share on S's stem or the representative's lineage, the "
                "recall of each, and why history sites are missed or other sites taken (by the congeners compared "
                "at the position).")
    reasons = sorted({k for r in site_rows for k in r if k.startswith(("missed", "other"))})
    report.table("l0_sites", site_rows, ["genus", "congeners", "twin", "regime", "sites", "share_in_history",
                                         "stem_recall", "rep_lineage_recall"] + reasons)
    consensus = [r for r in site_rows if int(r["congeners"]) >= oracle.MIN_CONGENERS and r["regime"] == "none"]
    if consensus:
        precision = min(r["share_in_history"] for r in consensus)
        recall = min(min(r["stem_recall"], r["rep_lineage_recall"]) for r in consensus)
        status = "PASS" if precision >= 0.98 and recall >= 0.95 else "FAIL"
        report.verdict("H3 sites", "consensus sites = stem + representative's lineage (no recombination, >= 3 "
                       "congeners): >= 0.98 of them, recall >= 0.95",
                       f"share {fnum(precision)}, recall {fnum(recall)}", status)
    fallback = [r for r in site_rows if int(r["congeners"]) < oracle.MIN_CONGENERS and r["regime"] == "none"]
    if fallback:
        share = float(np.mean([r["share_in_history"] for r in fallback]))
        report.verdict("H3 fallback", "with 2 congeners the nearest congener's own derived states enter (share in "
                       "history clearly below 1)", f"mean share {fnum(share)}", "PASS" if share < 0.8 else "FAIL")
    return db


def evaluate_l1(report, data):
    report.section("L1: the focus genomes' reads, with and without the allele scores")
    paths = data.paths
    counts = collections.Counter()
    moves = collections.Counter()
    for sample, _, _, _ in sample_reads(paths, "pe"):
        focus = {r["genome"]: (genus, r) for genus, r in data.design[("pe", sample)].items()}
        assigned = {}
        for run in ("default", "noallele"):
            path = os.path.join(paths.runs, run, "alignments", sample + ".sam.gz")
            for qname, mate, rname, _, mapq, _, _, _ in oracle.best_records(path):
                acc = data.origin("pe", sample, qname)
                if acc not in focus:
                    continue
                assigned.setdefault((qname, mate), {})[run] = int(rname.split("_")[0])
        for (qname, mate), by_run in assigned.items():
            acc = data.origin("pe", sample, qname)
            genus, r = focus[acc]
            target, sister = data.target_taxid[genus], data.sister_taxid[genus]

            def where(taxid):
                if taxid is None:
                    return "none"
                return "target" if taxid == target else "sister" if taxid == sister else \
                    "congener" if data.genus_of_taxid.get(taxid) == genus else "other"
            twin = data.genera[genus]["twin"]
            key = (r["class"], twin)
            a, b = by_run.get("default"), by_run.get("noallele")
            counts[key + ("reads",)] += 1
            counts[key + ("target_default",)] += a == target
            counts[key + ("target_noallele",)] += b == target
            if a != b:
                moves[key + (where(b), where(a))] += 1
    rows = []
    for cls in ("Q", "N"):
        for twin in ("0", "1"):
            key = (cls, twin)
            n = counts[key + ("reads",)]
            if not n:
                continue
            moved = sum(v for k, v in moves.items() if k[:2] == key)
            to_target = sum(v for k, v in moves.items() if k[:2] == key and k[3] == "target")
            rows.append({"class": cls, "twin": twin, "reads": n,
                         "on_target_noallele": counts[key + ("target_noallele",)] / n,
                         "on_target_default": counts[key + ("target_default",)] / n,
                         "moved": moved, "moved_to_target": to_target,
                         "share_moves_to_target": to_target / moved if moved else float("nan")})
    report.text("The paired-end focus reads (best record of each mate): the share on the target without and with the "
                "allele scores, and the reads the scores moved (Q: the target is right; N: no species is).")
    report.table("l1_reads", rows)
    move_rows = [{"class": k[0], "twin": k[1], "from": k[2], "to": k[3], "reads": v} for k, v in sorted(moves.items())]
    report.table("l1_moves", move_rows)
    unsure = settled = 0
    log_file = os.path.join(paths.runs, "default", "protal.log")
    pattern = re.compile(r"(\d+) reads with candidates of two species within \d+ mismatches took those shifts, "
                         r"(\d+) of them to another species")
    if os.path.exists(log_file):
        for line in open(log_file):
            m = pattern.search(line)
            if m:
                unsure += int(m.group(1))
                settled += int(m.group(2))
    report.text(f"The protal log's 'strain alleles:' lines (every sample of the default run): {unsure} unsure reads "
                f"took the shifts, {settled} of them to another species.")
    n_rows = [r for r in rows if r["class"] == "N"]
    if n_rows:
        reads = sum(r["reads"] for r in n_rows)
        drawn = sum(r["moved_to_target"] for r in n_rows)
        report.verdict("H7 novel reads", "the allele scores draw some of a novel species' reads onto the target (it "
                       "carries the alleles' ancestral base at the representative's own sites)",
                       f"{drawn} of {reads} reads ({fnum(drawn / reads if reads else float('nan'))})", "INFO")
    q_moves = [r for r in rows if r["class"] == "Q" and r["moved"]]
    if q_moves:
        precision = sum(r["moved_to_target"] for r in q_moves) / sum(r["moved"] for r in q_moves)
        report.verdict("H7 settling", "the strains' reads the allele scores move go to their species (>= 0.9)",
                       f"{fnum(precision)} of {sum(r['moved'] for r in q_moves)} moved", "PASS" if precision >= 0.9 else "FAIL")
    else:
        report.verdict("H7 settling", "the strains' reads the allele scores move go to their species",
                       "no strain read moved", "INFO")


def evaluate_oracle(report, data, db):
    """protal's features against the oracle's, sample by sample (error-free and paired-end runs)."""
    report.section("L2a: protal's features against the oracle's recomputation")
    features = ANCESTRY + ALLELES + POLYMORPHIC + WEIGHTS + DIAGNOSTIC
    rows, ceiling = [], {}
    for run, kind in (("default", "exact"), ("default", "pe"), ("shuffled", "exact"), ("noallele", "exact")):
        run_db = db if run != "shuffled" else shuffled_database(data)
        diffs = collections.defaultdict(list)
        for sample, _, _, _ in sample_reads(data.paths, kind):
            dump = os.path.join(data.paths.runs, run, "profiles", sample + ".profile.truth_annotated")
            sam = os.path.join(data.paths.runs, run, "alignments", sample + ".sam.gz")
            if not os.path.exists(dump) or not os.path.exists(sam):
                continue
            evidence = oracle.profile_sam(sam, run_db, data.targets)
            for r in read_tsv(dump):
                taxid = int(r["taxon"])
                if taxid not in data.targets or taxid not in evidence:
                    continue
                mine = evidence[taxid].features()
                if run == "default":
                    ceiling[(kind, sample, taxid)] = mine[CEILING]
                for f in features:
                    diffs[f].append(abs(value(r, f) - mine[f]))
        for f in features:
            d = diffs[f]
            if d:
                rows.append({"run": run, "set": kind, "feature": f, "rows": len(d), "max_abs_diff": max(d),
                             "share_within_1e-6": float(np.mean(np.array(d) <= 1e-6))})
    report.text("Every target row: |protal - oracle|. The oracle ports AncestrySites.h, StrainAlleles.h and the "
                "profiler's sums; a feature off here computes something else than its definition.")
    report.table("l2_oracle", rows)
    exact_rows = [r for r in rows if r["set"] == "exact"]
    worst = min((r["share_within_1e-6"] for r in exact_rows), default=float("nan"))
    report.verdict("H1b/H2b oracle", "protal's ancestry, allele and polymorphic features equal the oracle's "
                   "(error-free reads, every row within 1e-6)", f"worst share within 1e-6: {fnum(worst)}",
                   "PASS" if worst == 1.0 else "FAIL" if not math.isnan(worst) else "N/A")
    return ceiling


def shuffled_database(data):
    folder = os.path.join(data.paths.shuffled, "tables")
    db = oracle.Database(data.paths.unpacked)
    db.alleles = oracle.read_alleles(os.path.join(folder, "strain_alleles.tsv"))
    db.has_alleles = bool(db.alleles)
    return db


def evaluate_calibration(report, data):
    report.section("L2b: the planted genomes and the leaky control (error-free reads)")
    rows = [r for r in load_table(data.paths.table("default", "exact")) if r.get("meta_atp_class") == "P"
            and r.get("meta_atp_target") == "1"]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r["meta_atp_role"]].append(r)
    out = []
    for role, rs in sorted(groups.items()):
        truth = data.roles.get(rs[0]["meta_atp_focus"], {})
        cons = int(truth.get("cons_sites", 0) or 0)
        out.append({"role": role, "rows": len(rs),
                    "truth_agreement": int(truth["cons_agree"]) / cons if cons else float("nan"),
                    **{f: float(np.mean([value(r, f) for r in rs])) for f in
                       ["ancestry_agreement", "allele_explained_share", "polymorphic_known_share",
                        "polymorphic_novel_share", "ancestry_fixed_gain", "identity"]}})
    report.text("Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; "
                "A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome "
                "itself). truth_agreement: the genome's agreement at the consensus sites by position.")
    report.table("l2_calibration", out)
    leak = [o for o in out if o["role"] == "Q_leak"]
    if leak:
        explained, novel = leak[0]["allele_explained_share"], leak[0]["polymorphic_novel_share"]
        report.verdict("H1 leak", "an allele genome's own reads: allele_explained_share >= 0.95, "
                       "polymorphic_novel_share <= 0.02", f"{fnum(explained)}, {fnum(novel)}",
                       "PASS" if explained >= 0.95 and novel <= 0.02 else "FAIL")
    planted_a = [o for o in out if o["role"].startswith("A_a")]
    if planted_a:
        worst = max(abs(o["ancestry_agreement"] - float(o["role"][3:])) for o in planted_a)
        report.verdict("H2 calibration", "planted A_a: ancestry_agreement = a (+-0.03)", f"worst |diff| {fnum(worst)}",
                       "PASS" if worst <= 0.03 else "FAIL")
    planted_f = sorted((float(o["role"][3:]), o["allele_explained_share"]) for o in out if o["role"].startswith("P_f"))
    if planted_f:
        monotone = all(b[1] >= a[1] - 0.02 for a, b in zip(planted_f, planted_f[1:]))
        report.verdict("H1 calibration", "planted P_f: allele_explained_share rises with f (a step near 0.5: a read is "
                       "explained once it shares more than half of the allele's edits in its span)",
                       ", ".join(f"f {f:g}: {fnum(v)}" for f, v in planted_f), "PASS" if monotone else "FAIL")


def separation_rows(rows, features, strata_key=None):
    """AUC (Q high) of each feature, Q rows against N rows, optionally within strata."""
    q = [r for r in rows if r["meta_atp_class"] == "Q"]
    n = [r for r in rows if r["meta_atp_class"] == "N"]
    out = {}
    for f in features:
        qv, nv = [value(r, f) for r in q], [value(r, f) for r in n]
        keep_q = [i for i, v in enumerate(qv) if not math.isnan(v)]
        keep_n = [i for i, v in enumerate(nv) if not math.isnan(v)]
        if strata_key:
            out[f] = stratified_auc([qv[i] for i in keep_q], [nv[i] for i in keep_n],
                                    [strata_key(q[i]) for i in keep_q], [strata_key(n[i]) for i in keep_n])
        else:
            out[f] = auc([qv[i] for i in keep_q], [nv[i] for i in keep_n])
    return out, len(q), len(n)


def identity_band(r):
    v = value(r, "identity")
    return -1 if math.isnan(v) else int(v / 0.005)


def evaluate_separation(report, data, kind, ceiling):
    rows = target_rows(load_table(data.paths.table("default", kind)))
    for r in rows:
        r[CEILING] = str(ceiling.get((kind, r["meta_sample"], int(r["taxon"])), float("nan")))
    features = REFERENCE + ANCESTRY[:3] + ALLELES + POLYMORPHIC + WEIGHTS + [DERIVED, CEILING]
    report.section(f"L2c: separation, {kind} (target rows: Q = its strain present, N = a novel species' reads)")
    comparisons = [("all Q vs all N", lambda r: True, lambda r: True)]
    for b in ("0", "0.5", "0.8", "0.95", "1"):
        comparisons.append((f"Q_deep vs N_{b}", lambda r: r["meta_atp_role"] == "Q_deep",
                            lambda r, b=b: r["meta_atp_role"] == f"N_{b}"))
    comparisons += [
        ("Q_near vs N_0.95", lambda r: r["meta_atp_role"] == "Q_near", lambda r: r["meta_atp_role"] == "N_0.95"),
        ("Q_lone vs N_0.95", lambda r: r["meta_atp_role"] == "Q_lone", lambda r: r["meta_atp_role"] == "N_0.95"),
        ("Q_rep vs N_0", lambda r: r["meta_atp_role"] == "Q_rep", lambda r: r["meta_atp_role"] == "N_0"),
        ("Q_imp vs N (b <= 0.95)", lambda r: r["meta_atp_role"].startswith("Q_imp"),
         lambda r: r["meta_atp_role"] in ("N_0", "N_0.5", "N_0.8", "N_0.95")),
        ("Q (no import) vs N_imp", lambda r: not r["meta_atp_role"].startswith("Q_imp"),
         lambda r: r["meta_atp_role"].startswith("N_imp")),
        ("Q_ils vs N_0.95", lambda r: r["meta_atp_role"] == "Q_ils", lambda r: r["meta_atp_role"] == "N_0.95"),
        # the decoy: a tip of Q_deep's clade with another label, so nothing but the label tells them apart
        ("decoy: Q_deep vs decoy", lambda r: r["meta_atp_role"] == "Q_deep", lambda r: r["meta_atp_class"] == "D"),
    ]
    out = []
    for label, q_filter, n_filter in comparisons:
        negative = "D" if label.startswith("decoy") else "N"
        subset = [r for r in rows if r["meta_atp_class"] == "Q" and q_filter(r)] + \
                 [{**r, "meta_atp_class": "N"} for r in rows if r["meta_atp_class"] == negative and n_filter(r)]
        aucs, nq, nn = separation_rows(subset, features)
        out.append({"comparison": label, "Q": nq, "N": nn, **aucs})
    for factor, key in (("allele genomes", "meta_atp_alleles"), ("congeners", "meta_atp_congeners"),
                        ("twin", "meta_atp_twin"), ("regime", "meta_atp_regime"), ("depth", "meta_atp_depth")):
        for level in sorted({r[key] for r in rows}, key=lambda x: (len(x), x)):
            subset = [r for r in rows if r[key] == level]
            aucs, nq, nn = separation_rows(subset, features)
            out.append({"comparison": f"{factor} {level}", "Q": nq, "N": nn, **aucs})
    banded, nq, nn = separation_rows(rows, features, identity_band)
    out.append({"comparison": "all, within identity bands of 0.005", "Q": nq, "N": nn, **banded})
    report.text("AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement"
                " + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies "
                "at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.")
    report.table(f"l2_auc_{kind}", out, ["comparison", "Q", "N"] + features)
    # The pairs of section 2 by the factors that decide them.
    pairs = [("Q_deep vs N_0", "N_0"), ("Q_deep vs N_0.5", "N_0.5"), ("Q_deep vs N_0.8", "N_0.8"),
             ("Q_deep vs N_0.95", "N_0.95")]
    strata = [("congeners 2", lambda r: r["meta_atp_congeners"] == "2"),
              ("congeners 6", lambda r: r["meta_atp_congeners"] == "6"),
              ("no twin", lambda r: r["meta_atp_twin"] == "0"), ("twin", lambda r: r["meta_atp_twin"] == "1"),
              ("no allele genomes", lambda r: r["meta_atp_alleles"] == "0"),
              ("allele genomes", lambda r: r["meta_atp_alleles"] != "0"),
              ("6 congeners, allele genomes, no recombination", lambda r: r["meta_atp_congeners"] == "6" and
               r["meta_atp_alleles"] != "0" and r["meta_atp_regime"] == "none"),
              ("recombination none", lambda r: r["meta_atp_regime"] == "none"),
              ("recombination high", lambda r: r["meta_atp_regime"] == "high")]
    key_features = ["identity", "lu_per_kb", "ancestry_agreement", "allele_explained_share",
                    "polymorphic_known_share", "ancestry_fixed_agreement", DERIVED, CEILING]
    by_stratum = []
    for label, n_role in pairs:
        for stratum, keep in strata:
            subset = [r for r in rows if keep(r) and (r["meta_atp_role"] == "Q_deep" or r["meta_atp_role"] == n_role)]
            aucs, nq, nn = separation_rows(subset, key_features)
            by_stratum.append({"pair": label, "stratum": stratum, "Q": nq, "N": nn, **aucs})
    report.text("")
    report.text("The deep strains against the novel species by where they left S's stem, within the strata that "
                "decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).")
    report.table(f"l2_pairs_{kind}", by_stratum, ["pair", "stratum", "Q", "N"] + key_features, show=60)
    return rows, out, by_stratum


def increment_rows(rows, kind):
    """Grouped-CV logistic models: identity and depth, then with each group; AUC gain with a genus bootstrap."""
    y = np.array([1.0 if r["meta_atp_class"] == "Q" else 0.0 for r in rows])
    groups = [r["meta_atp_genus"] for r in rows]
    base_cols = ["identity", "top_identity", "fragments"]
    sets = {"+ancestry": ANCESTRY[:3], "+alleles": ALLELES, "+polymorphic": POLYMORPHIC, "+weights": WEIGHTS,
            "+all three": ANCESTRY[:3] + ALLELES + POLYMORPHIC, "+ancestry, polymorphic": ANCESTRY[:3] + POLYMORPHIC,
            "+all four": ANCESTRY[:3] + ALLELES + POLYMORPHIC + WEIGHTS,
            "+oracle fixed agreement": [CEILING]}

    def matrix(cols):
        X = np.array([[value(r, c) for c in cols] for r in rows])
        X[np.isnan(X)] = -1
        return X
    base = grouped_cv_scores(matrix(base_cols), y, groups)
    out = []
    for name, cols in sets.items():
        scores = grouped_cv_scores(matrix(base_cols + cols), y, groups)
        d, lo, hi = bootstrap_difference(base, scores, y, groups)
        out.append({"set": kind, "model": "identity, depth " + name, "auc_base": auc(base[y == 1], base[y == 0]),
                    "auc": auc(scores[y == 1], scores[y == 0]), "gain": d, "gain_lo": lo, "gain_hi": hi})
    return out


def bins_q(r):
    """The bins of a Q row for the sensitivity tables: how long ago its lineage diverged (from the representative,
    from the nearest allele genome), how much of the representative's own derived sites the alleles cover, how much
    of it came from a congener, and the genus's factors."""
    out = [("role", r["meta_atp_role"].rstrip("0123456789.") if r["meta_atp_role"].startswith("Q_imp")
            else r["meta_atp_role"])]
    d = value(r, "meta_atp_distance")
    out.append(("distance to rep", "<0.005" if d < 0.005 else "0.005-0.015" if d < 0.015 else
                "0.015-0.03" if d < 0.03 else ">=0.03"))
    m = value(r, "meta_atp_mrca_allele")
    out.append(("MRCA with nearest allele genome", "no allele genome" if math.isnan(m) else
                "recent (<0.002)" if m < 0.002 else "0.002-0.01" if m < 0.01 else "long ago (>=0.01)"))
    c = value(r, "meta_atp_lineage_covered")
    out.append(("rep-lineage sites it lacks that the alleles cover", "no such site" if math.isnan(c) else
                "none (<0.1)" if c < 0.1 else "partly" if c < 0.9 else "all (>=0.9)"))
    imp = value(r, "meta_atp_genes_congeneric")
    share = 0 if math.isnan(imp) else imp / 120
    out.append(("genes with a congener's segment", "none" if share == 0 else "<=0.1" if share <= 0.1 else
                "0.1-0.3" if share <= 0.3 else ">0.3"))
    out.append(("regime", r["meta_atp_regime"]))
    out.append(("twin genus", r["meta_atp_twin"]))
    out.append(("depth", r["meta_atp_depth"]))
    return out


def bins_n(r):
    """The bins of an N row: where it left S's stem (its role), its imports, the genus's factors."""
    imp = value(r, "meta_atp_genes_congeneric")
    return [("role", r["meta_atp_role"]), ("regime", r["meta_atp_regime"]), ("twin genus", r["meta_atp_twin"]),
            ("depth", r["meta_atp_depth"]),
            ("genes with a segment of S or a congener", "none" if not imp or math.isnan(imp) else
             "<=0.3" if imp / 120 <= 0.3 else ">0.3")]


def sensitivity_rows(rows, kind):
    """TPR of the Q rows by how long ago their lineage diverged and how much came from a congener, FPR of the N rows
    by b and imports, at the threshold that 5% of all N rows exceed; per feature (oriented by its overall AUC)."""
    features = ["identity", "ancestry_agreement", "allele_explained_share", "polymorphic_known_share", DERIVED, CEILING]
    q = [r for r in rows if r["meta_atp_class"] == "Q"]
    n = [r for r in rows if r["meta_atp_class"] == "N"]
    cells = collections.defaultdict(dict)  # (class, by, level) -> {feature: rate}, and the rows
    for f in features:
        qv = np.array([value(r, f) for r in q])
        nv = np.array([value(r, f) for r in n])
        if np.isnan(qv).all() or not len(nv):
            continue
        qv[np.isnan(qv)] = -1
        nv[np.isnan(nv)] = -1
        sign = 1 if auc(qv, nv) >= 0.5 else -1
        _, threshold = tpr_at_fpr(sign * qv, sign * nv)
        counts = collections.defaultdict(lambda: [0, 0])
        for cls, rs, vs, binner in (("Q", q, sign * qv, bins_q), ("N", n, sign * nv, bins_n)):
            for r, v in zip(rs, vs):
                for by, level in binner(r):
                    counts[(cls, by, level)][0] += v > threshold
                    counts[(cls, by, level)][1] += 1
        for key, (hit, total) in counts.items():
            cells[key][f] = hit / total if total else float("nan")
            cells[key]["rows"] = total
    order = {"Q": 0, "N": 1}
    return [{"set": kind, "class": cls, "by": by, "level": level, **values}
            for (cls, by, level), values in sorted(cells.items(), key=lambda kv: (order[kv[0][0]], kv[0][1], kv[0][2]))]


def evaluate_models(report, opts, data):
    folder = os.path.join(opts.outdir, "models")
    report.section("L3: the presence models (trained on the train world, called on the test world)")
    if not os.path.isdir(folder) or not any(f.endswith(".calls.tsv.gz") for f in os.listdir(folder)):
        report.text("No models (no scikit-learn in --train-python, or not trained yet).")
        return
    out = []
    by_genus = {}
    for name in sorted(os.listdir(folder)):
        if not name.endswith(".calls.tsv.gz"):
            continue
        label = name[:-len(".calls.tsv.gz")]
        rows = [r for r in read_tsv(os.path.join(folder, name)) if r.get("set") == "test"]
        target = [r for r in rows if r.get("meta_atp_target") == "1" and r.get("meta_atp_class") in ("Q", "N")]
        strata = {"target rows": lambda r: True,
                  "hard: Q_deep vs N_0.8/0.95": lambda r: r["meta_atp_role"] in ("Q_deep", "N_0.8", "N_0.95"),
                  "twin genera": lambda r: r["meta_atp_twin"] == "1",
                  "no allele genomes": lambda r: r["meta_atp_alleles"] == "0",
                  "recombination high": lambda r: r["meta_atp_regime"] == "high"}
        for stratum, keep in strata.items():
            rs = [r for r in target if keep(r)]
            q = [value(r, "p") for r in rs if r["meta_atp_class"] == "Q"]
            n = [value(r, "p") for r in rs if r["meta_atp_class"] == "N"]
            calls_q = [r["call"] == "1" for r in rs if r["meta_atp_class"] == "Q"]
            calls_n = [r["call"] == "1" for r in rs if r["meta_atp_class"] == "N"]
            tpr, _ = tpr_at_fpr(q, n)
            out.append({"model": label, "stratum": stratum, "Q": len(q), "N": len(n), "auc": auc(q, n),
                        "tpr_at_5pct_fpr": tpr,
                        "sensitivity_at_knob": float(np.mean(calls_q)) if calls_q else float("nan"),
                        "fpr_at_knob": float(np.mean(calls_n)) if calls_n else float("nan")})
        by_genus[label] = target
        everything = [r for r in rows if r.get("meta_atp_class") != "D"]  # the decoy's label is noise
        y = np.array([value(r, "truth") for r in everything])
        call = np.array([r["call"] == "1" for r in everything])
        tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
        out.append({"model": label, "stratum": "all test rows", "Q": int((y == 1).sum()), "N": int((y == 0).sum()),
                    "auc": auc([value(r, "p") for r, t in zip(everything, y) if t == 1],
                               [value(r, "p") for r, t in zip(everything, y) if t == 0]),
                    "tpr_at_5pct_fpr": float("nan"),
                    "sensitivity_at_knob": tp / (tp + fn) if tp + fn else float("nan"),
                    "fpr_at_knob": float(np.mean(call[y == 0])) if (y == 0).any() else float("nan"),
                    "f1_at_knob": 2 * tp / (2 * tp + fp + fn) if tp else float("nan")})
    report.text("Models by read type and feature set (base: the default set without ancestry, alleles and polymorphic; "
                "_on_noallele / _on_shuffled: the default model on the test world's runs without allele scores and with "
                "the shuffled table). Q and N of 'all test rows': present and absent taxa of every row.")
    report.table("l3_models", out, ["model", "stratum", "Q", "N", "auc", "tpr_at_5pct_fpr", "sensitivity_at_knob",
                                    "fpr_at_knob", "f1_at_knob"], show=200)
    sensitivity, columns = model_sensitivity(by_genus)
    if sensitivity:
        report.text("")
        report.text("The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: "
                    "sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N "
                    "rows exceed.")
        report.table("l3_sensitivity", sensitivity, columns, show=200)
    for kind in ("pe", "se", "hifi"):
        base, full = by_genus.get(f"{kind}_base"), by_genus.get(f"{kind}_default")
        if not base or not full:
            continue
        key = {(r["meta_sample"], r["taxon"]): r for r in base}
        pairs = [(key[(r["meta_sample"], r["taxon"])], r) for r in full if (r["meta_sample"], r["taxon"]) in key]
        hard = [(a, b) for a, b in pairs if a["meta_atp_role"] in ("Q_deep", "N_0.8", "N_0.95")]
        for label, subset in (("target rows", pairs), ("hard", hard)):
            if not subset:
                continue
            y = np.array([1.0 if a["meta_atp_class"] == "Q" else 0.0 for a, _ in subset])
            sa = np.array([value(a, "p") for a, _ in subset])
            sb = np.array([value(b, "p") for _, b in subset])
            d, lo, hi = bootstrap_difference(sa, sb, y, [a["meta_atp_genus"] for a, _ in subset])
            status = "PASS" if lo > 0 else "FAIL" if hi < 0 else "INCONCLUSIVE"
            report.verdict(f"H8 calls {kind} {label}", "the default set's model separates better than base (AUC gain, "
                           "95% CI by genus above 0)", f"{fnum(d)} [{fnum(lo)}, {fnum(hi)}]", status)
    default_pe, shuffled = by_genus.get("pe_default"), by_genus.get("pe_default_on_shuffled")
    if default_pe and shuffled:
        a = auc([value(r, "p") for r in default_pe if r["meta_atp_class"] == "Q"],
                [value(r, "p") for r in default_pe if r["meta_atp_class"] == "N"])
        b = auc([value(r, "p") for r in shuffled if r["meta_atp_class"] == "Q"],
                [value(r, "p") for r in shuffled if r["meta_atp_class"] == "N"])
        report.verdict("H9 shuffled", "shuffled alleles cost the default model separation (the alleles carry "
                       "information)", f"AUC {fnum(a)} -> {fnum(b)}", "PASS" if b < a else "FAIL")


def model_sensitivity(by_model):
    """The models' calls by bin (bins_q, bins_n), base against the other feature sets, per read type."""
    out, names = [], ["base", "ancestry", "alleles", "default"]
    for kind in ("pe", "se", "hifi"):
        models = [m for m in names if f"{kind}_{m}" in by_model]
        if not models:
            continue
        cells = collections.defaultdict(dict)
        for m in models:
            rows = by_model[f"{kind}_{m}"]
            q = [r for r in rows if r["meta_atp_class"] == "Q"]
            n = [r for r in rows if r["meta_atp_class"] == "N"]
            _, threshold = tpr_at_fpr([value(r, "p") for r in q], [value(r, "p") for r in n])
            for cls, rs, binner in (("Q", q, bins_q), ("N", n, bins_n)):
                counts = collections.defaultdict(lambda: [0, 0, 0])
                for r in rs:
                    for key in binner(r):
                        c = counts[key]
                        c[0] += r["call"] == "1"
                        c[1] += value(r, "p") > threshold
                        c[2] += 1
                for (by, level), (called, above, total) in counts.items():
                    cell = cells[(cls, by, level)]
                    cell["rows"] = total
                    cell[f"{m} called"] = called / total
                    cell[f"{m} at 5% FPR"] = above / total
        order = {"Q": 0, "N": 1}
        for (cls, by, level), cell in sorted(cells.items(), key=lambda kv: (order[kv[0][0]], kv[0][1], kv[0][2])):
            out.append({"set": kind, "class": cls, "by": by, "level": level, **cell})
    columns = ["set", "class", "by", "level", "rows"] + [f"{m} called" for m in names] + \
              [f"{m} at 5% FPR" for m in names]
    return out, columns


def command_evaluate(opts):
    test = Paths(opts.outdir, "test")
    data = WorldData(test)
    report = Report(os.path.join(opts.outdir, "report"))
    db = evaluate_l0(report, data)
    evaluate_l1(report, data)
    ceiling = evaluate_oracle(report, data, db)
    evaluate_calibration(report, data)
    increments, sensitivity = [], []
    for kind in ("pe", "se", "hifi"):
        if not os.path.exists(test.table("default", kind)):
            continue
        rows, aucs, by_stratum = evaluate_separation(report, data, kind, ceiling)
        if kind == "pe":
            check_predictions(report, aucs, by_stratum, rows)
        qn = [r for r in rows if r["meta_atp_class"] in ("Q", "N")]
        increments += increment_rows(qn, kind)
        increments += increment_rows([r for r in qn if r["meta_atp_alleles"] != "0"], kind + ", genera with alleles")
        sensitivity += sensitivity_rows(qn, kind)
    report.section("L2d: the gain of each group over identity and depth (logistic models, genera held out)")
    report.text("AUC of grouped-CV logistic models on the target rows; gain over identity, top identity and fragments "
                "with a 95% interval over resampled genera.")
    report.table("l2_increment", increments)
    for r in increments:
        if r["set"] == "pe, genera with alleles" and r["model"].endswith("+oracle fixed agreement"):
            report.verdict("H5 ceiling, genera with alleles", "the true fixed-site agreement (the oracle's) over "
                           "identity and depth, where the fixed sites exist; compare with '+all three'",
                           f"{fnum(r['gain'])} [{fnum(r['gain_lo'])}, {fnum(r['gain_hi'])}]", "INFO")
        if r["set"] == "pe" and r["model"].endswith("+all three"):
            status = "PASS" if r["gain_lo"] > 0 else "FAIL"
            report.verdict("H4/H5 gain", "the three groups add separation over identity and depth (pe, CI above 0)",
                           f"{fnum(r['gain'])} [{fnum(r['gain_lo'])}, {fnum(r['gain_hi'])}]", status)
        if r["set"] == "pe" and r["model"].endswith("+oracle fixed agreement"):
            report.verdict("H5 ceiling", "the true fixed-site agreement alone (the oracle's) adds separation; compare "
                           "with '+all three' for the gap", f"{fnum(r['gain'])} [{fnum(r['gain_lo'])}, "
                           f"{fnum(r['gain_hi'])}]", "INFO")
    report.section("Sensitivity: by divergence of the lineage and by recombination with congeners")
    report.text("Per feature, the threshold 5% of all N rows exceed (the feature oriented by its AUC); rate_above: the "
                "share of Q rows above it (sensitivity) or of N rows (false-positive rate) in each bin. 'MRCA with "
                "nearest allele genome' is how recently the strain shared an ancestor with a strain the database "
                "knows; 'genes with a congener's segment' how much of it came from a congener.")
    report.table("sensitivity", sensitivity, ["set", "class", "by", "level", "rows", "identity", "ancestry_agreement",
                                              "allele_explained_share", "polymorphic_known_share", DERIVED, CEILING],
                 show=200)
    evaluate_models(report, opts, data)
    header = (f"# Ancestry true-positive test: {opts.outdir}\n\nWorlds: test (seed {WORLD_SEEDS['test']}), train "
              f"(seed {WORLD_SEEDS['train']}); {len(data.genera)} test genera. See "
              "docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.")
    report.write(header)
    log(f"report: {os.path.join(opts.outdir, 'report', 'summary.md')}")


def check_predictions(report, aucs, by_stratum, rows):
    """The separation predictions (H4, H5, H6, H9's decoy) from the pe AUC tables and rows."""
    by = {r["comparison"]: r for r in aucs}
    pairs = {(r["pair"], r["stratum"]): r for r in by_stratum}
    clean = [pairs.get((f"Q_deep vs N_{b}", "6 congeners, allele genomes, no recombination"), {}).get(
        "ancestry_agreement", float("nan")) for b in ("0", "0.5")]
    if not any(math.isnan(v) for v in clean):
        report.verdict("H4 ancestry, clean strata", "the same with >= 3 congeners and no recombination (the "
                       "consensus applies)", f"b 0: {fnum(clean[0])}, b 0.5: {fnum(clean[1])}", "INFO")
    # H5 as registered: deep strains whose representative-lineage sites the alleles cover (>= 0.8) against the late
    # novel species (b 0.8, 0.95) of genera with alleles: does ancestry_fixed_gain add to ancestry_agreement?
    subset = [r for r in rows if (r["meta_atp_role"] == "Q_deep" and value(r, "meta_atp_lineage_covered") >= 0.8) or
              (r["meta_atp_role"] in ("N_0.8", "N_0.95") and r["meta_atp_alleles"] != "0")]
    if len({r["meta_atp_class"] for r in subset}) == 2:
        y = np.array([1.0 if r["meta_atp_class"] == "Q" else 0.0 for r in subset])
        groups = [r["meta_atp_genus"] for r in subset]

        def scores(cols):
            X = np.array([[value(r, c) for c in cols] for r in subset])
            X[np.isnan(X)] = -1
            return grouped_cv_scores(X, y, groups)
        base, plus = scores(["ancestry_agreement"]), scores(["ancestry_agreement", "ancestry_fixed_gain"])
        d, lo, hi = bootstrap_difference(base, plus, y, groups)
        ceiling = auc([value(r, CEILING) for r in subset if r["meta_atp_class"] == "Q"],
                      [value(r, CEILING) for r in subset if r["meta_atp_class"] == "N"])
        direct = auc([value(r, "ancestry_fixed_agreement") for r in subset if r["meta_atp_class"] == "Q"],
                     [value(r, "ancestry_fixed_agreement") for r in subset if r["meta_atp_class"] == "N"])
        report.verdict("H5 fixed sites", "deep strains with their representative-lineage sites covered (>= 0.8) "
                       "against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to "
                       "ancestry_agreement (CI above 0)",
                       f"AUC {fnum(auc(base[y == 1], base[y == 0]))} -> {fnum(auc(plus[y == 1], plus[y == 0]))}, "
                       f"gain {fnum(d)} [{fnum(lo)}, {fnum(hi)}]; ancestry_fixed_agreement alone {fnum(direct)}, "
                       f"the true fixed-site agreement alone {fnum(ceiling)} "
                       f"({int(y.sum())} Q, {int((1 - y).sum())} N rows)", "PASS" if lo > 0 else "FAIL")

    def get(comparison, feature):
        return by.get(comparison, {}).get(feature, float("nan"))
    a0, a5 = get("Q_deep vs N_0", "ancestry_agreement"), get("Q_deep vs N_0.5", "ancestry_agreement")
    if not math.isnan(a0):
        report.verdict("H4 ancestry", "ancestry_agreement separates deep strains from novel species that left the stem "
                       "early (b <= 0.5): AUC >= 0.9", f"b 0: {fnum(a0)}, b 0.5: {fnum(a5)}",
                       "PASS" if min(a0, a5) >= 0.9 else "FAIL")
    near, lone = get("Q_near vs N_0.95", "allele_explained_share"), get("Q_lone vs N_0.95", "allele_explained_share")
    if not math.isnan(near):
        report.verdict("H6 confound", "allele_explained_share favours strains with a stored relative (Q_near) over "
                       "late twins, but not lone strains", f"Q_near {fnum(near)}, Q_lone {fnum(lone)}",
                       "PASS" if near > 0.5 and (math.isnan(lone) or lone < near) else "INFO")
    decoy = by.get("decoy: Q_deep vs decoy", {})
    values = {k: v for k, v in decoy.items() if k not in ("comparison", "Q", "N") and not math.isnan(v)}
    if values:
        worst = max(values, key=lambda k: abs(values[k] - 0.5))
        n = min(decoy.get("Q", 0), decoy.get("N", 0))
        report.verdict("H9 decoy", "the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no "
                       "feature separates (AUC 0.5 +- 0.1)", f"largest |AUC - 0.5|: {worst} {fnum(values[worst])} "
                       f"({decoy.get('Q')} Q, {decoy.get('N')} decoy rows)",
                       "PASS" if abs(values[worst] - 0.5) <= 0.1 else "FAIL" if n >= 10 else "TOO FEW ROWS")
    n1 = by.get("Q_deep vs N_1", {})
    if n1:
        report.verdict("H9b crown lineage", "N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry "
                       "features cannot tell them apart; the allele features can where Q_deep shares its side of the "
                       "crown with allele genomes", f"ancestry_agreement {fnum(n1.get('ancestry_agreement'))}, "
                       f"allele_explained_share {fnum(n1.get('allele_explained_share'))}", "INFO")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="everything, resumable")
    e = sub.add_parser("evaluate", help="the evaluation of a filled --outdir")
    for q in (r, e):
        q.add_argument("--outdir", required=True)
        q.add_argument("-t", "--threads", type=int, default=4)
        q.add_argument("--quick", action="store_true", help="small worlds (simulate_ancestry_world.py --quick)")
    r.add_argument("--protal", default=os.environ.get("PROTAL", "build/protal"))
    r.add_argument("--simulate", default=os.environ.get("SIMULATE", "build/simulate_metagenomes"))
    r.add_argument("--train-python", default=os.environ.get("PROTAL_TRAIN_PYTHON", sys.executable))
    r.add_argument("--strain-alleles", dest="strain_alleles", default="4", help="protal --build --strain_alleles")
    r.add_argument("--world-args", dest="world_args", default="",
                   help="more options for simulate_ancestry_world.py, in one string")
    r.add_argument("--model-sets", dest="model_sets", default=",".join(FEATURE_SETS),
                   help=f"feature sets of the models (default all: {', '.join(FEATURE_SETS)}); se and HiFi get base and "
                        "default of them unless --model-sets-all")
    r.add_argument("--model-sets-all", dest="model_sets_all", action="store_true")
    r.add_argument("--model-types", dest="model_types", default="pe,se,hifi", help="read types with models")
    opts = p.parse_args(argv)
    if opts.command == "run":
        opts.world_args = opts.world_args.split()
        opts.model_sets = [s for s in opts.model_sets.split(",") if s]
        opts.model_types = [s for s in opts.model_types.split(",") if s]
        unknown = [s for s in opts.model_sets if s not in FEATURE_SETS] + \
                  [t for t in opts.model_types if t not in ("pe", "se", "hifi")]
        if unknown:
            sys.exit(f"unknown model sets or read types: {', '.join(unknown)}")
        for tool in (opts.protal, opts.simulate):
            if not (os.path.isfile(tool) and os.access(tool, os.X_OK)):
                sys.exit(f"not an executable: {tool}")
        opts.protal, opts.simulate = os.path.abspath(opts.protal), os.path.abspath(opts.simulate)
        opts.outdir = os.path.abspath(opts.outdir)
        command_run(opts)
    else:
        opts.outdir = os.path.abspath(opts.outdir)
        command_evaluate(opts)


if __name__ == "__main__":
    main()
