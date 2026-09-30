#!/usr/bin/env python3
"""Simulate metagenomes and profile them with protal, to train its presence models.

For every design point, a read setup (length, ART profile, fragment size) and a sequencing depth,
simulate_metagenomes draws random communities from a genome table and protal profiles them
against a database, knowing the true species. The training dumps of all samples
(<profile>.truth_annotated: every taxon protal saw, its features and whether it was present) are
joined into one table, with meta_* columns saying where each row comes from (meta_domain from
the genome table's lineages; with --novel_species and --taxonomy, which rows share a genus with a
species the database lacks, and which present species were simulated from another genome than
the database's reference). Train on it with

    python3 scripts/random_forest_cmdline.py --truth-file OUT/training_data.tsv --output-prefix OUT/model

Read types (--read_types): pe, the paired-end samples above; se, the same samples' first reads
alone, profiled as single-end reads; pb and ont, long reads of the same communities simulated with
pbsim3 (--pb_setup, --ont_setup), --long_read_bases per sample, one design point each. Every read
type gets its table: training_data.tsv (pe), training_data_se.tsv, training_data_pb.tsv,
training_data_ont.tsv; protal profiles each sample with that read type's model and settings.

A genome table with species the database lacks gives the negatives that matter most: a relative
the database has picks up their reads. When whole genera, families, orders, classes or phyla are
missing (--novel_species with ranks, --novel_clades), the relatives are distant: meta_novel_level
marks the absent taxa closest to such species, meta_neighbour_rank how close a present taxon's
nearest other species in the sample is. Archaea (--archaea) have fewer marker genes than bacteria
and need to be in the training data, too. Design points are simulated in parallel (--jobs; ART and
pbsim3 simulate one genome at a time), then the samples of all read types are profiled in one
protal run, which loads the database once. Points already simulated or profiled are skipped, so a
run can be resumed.

usage: collect_training_data.py --db DB --genome_table genomes.tsv -o OUT [options]
"""

import argparse
import collections
import concurrent.futures
import csv
import glob
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lineages  # noqa: E402

READ_TYPES = ("pe", "se", "pb", "ont")
LONG_READ_TYPES = ("pb", "ont")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="protal database")
    p.add_argument("--genome_table", required=True,
                   help="simulate_metagenomes genome table (accession, GTDB taxonomy, FASTA path)")
    p.add_argument("-o", "--out", required=True, help="output directory")
    p.add_argument("--protal", default="protal", help="protal binary (default: protal on PATH)")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes binary")
    p.add_argument("--samples", type=int, default=4, help="samples per design point (default: 4)")
    p.add_argument("--read_pairs", default="1000,5000,20000,100000,500000",
                   help="comma-separated read pairs per sample, one design point each")
    p.add_argument("--read_setups", default="100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50",
                   help="comma-separated LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, one design point each "
                        "(HSXt: HiSeq X, the closest of ART's profiles to NovaSeq; file=R1.txt+R2.txt: quality "
                        "profiles art_profiler_illumina made from real reads, e.g. NovaSeq)")
    p.add_argument("--species_per_sample", default="5-30", help="species per sample, N or MIN-MAX")
    p.add_argument("--strains_per_species", default="",
                   help="probabilities of a second, third, ... strain of a species in a sample, e.g. 0.3,0.1 "
                        "(simulate_metagenomes --strains_per_species; default: one strain each)")
    p.add_argument("--abundance", default="",
                   help="abundance model: lognormal:SIGMA, powerlaw:ALPHA or negbin:R:P (default: the simulator's, "
                        "Poisson-lognormal with sigma 1.3)")
    p.add_argument("--archaea", type=int, default=0, help="archaeal species per sample (default: 0)")
    p.add_argument("--congeners", type=int, default=0,
                   help="species of one genus in every sample of a design point, the genus drawn per point among "
                        "those with that many species (default: 0). Species are otherwise drawn uniformly, so among "
                        "many genera relatives hardly ever share a sample, while in real samples they often do")
    p.add_argument("--read_types", default="pe",
                   help="comma-separated read types to collect: pe, se, pb, ont (default pe)")
    p.add_argument("--long_read_bases", default="300000,1500000,6000000,30000000,150000000",
                   help="bases per long-read sample (pb, ont), one design point each; point i replays the "
                        "communities of paired-end point i (default: about the paired-end points' bases)")
    p.add_argument("--pb_setup", default="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of PacBio reads (default: HiFi-like "
                        "reads from the Sequel error model)")
    p.add_argument("--ont_setup", default="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of Nanopore reads")
    p.add_argument("--pbsim", default="pbsim", help="pbsim3 binary (default: pbsim on PATH)")
    p.add_argument("--pbsim_models", help="folder of pbsim3's .model files (default: found next to the binary)")
    p.add_argument("-t", "--threads", type=int, default=4, help="threads of the protal run (default 4)")
    p.add_argument("--jobs", type=int, default=0,
                   help="design points (and long-read genomes) simulated at a time (default: --threads)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--novel_species",
                   help="species the database lacks, one per line, optionally with the rank they were held out at "
                        "and the clade (heldout_species.txt of build_gtdb_database.py): meta_novel_species counts "
                        "them in each sample, meta_novel_congener marks the taxa of their genera, which their reads "
                        "land on, and meta_novel_level (with --taxonomy) the absent taxa whose closest species in the "
                        "sample is one of them, by the rank it was held out at")
    p.add_argument("--novel_clades", type=int, default=0,
                   help="species of held-out clades (--novel_species with ranks above species) in every sample, per "
                        "rank: each design point takes one clade of each rank, in turn, so that every clade is used "
                        "before one is used again (default 0: species are drawn uniformly, and the few of held-out "
                        "clades hardly appear)")
    p.add_argument("--taxonomy", help="internal_taxonomy.dmp of the database: meta_rep_genome says whether a present "
                                      "species was simulated from its representative genome (the database's "
                                      "reference, 1) or from another strain (0)")
    opts = p.parse_args(argv)
    opts.read_types = [t.strip() for t in opts.read_types.split(",") if t.strip()]
    unknown = [t for t in opts.read_types if t not in READ_TYPES]
    if unknown:
        p.error(f"--read_types: unknown {', '.join(unknown)} (pe, se, pb, ont)")
    return opts


META_COLUMNS = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_domain",
                "meta_novel_species", "meta_novel_congener", "meta_rep_genome", "meta_novel_levels",
                "meta_novel_level", "meta_relative_rank", "meta_neighbour_rank", "meta_read_type"]
# meta_read_pairs: the design point's depth, read pairs (pe; reads for se) or bases (pb, ont).
# meta_novel_levels: the sample's species the database lacks, by the rank they were held out at
# ("species:2,family:1"). meta_relative_rank: the deepest rank the taxon shares with a species simulated in
# the sample ("species" for the simulated species themselves, "none" for no shared domain). meta_novel_level:
# for an absent taxon whose closest simulated species (at least as close as any present one) is one the
# database lacks, the rank that species was held out at: its reads are the likely source of the taxon's.
# meta_neighbour_rank: for a present taxon, the deepest rank it shares with another species simulated in the
# sample (a congener's reads fit it nearly as well, so it may be missed).
TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
MAP_COLUMNS = ["SAMPLEID", "FIRST", "SECOND", "SAM", "PREFIX", "PROFILE", "PROFILE_TRUTH", "READ_TYPE"]


def large_genera(genome_table, size):
    """Genera with at least `size` species in the genome table, sorted."""
    species = collections.defaultdict(set)
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                ranks = lineage.split(";")
                species[ranks[5][3:]].add(ranks[-1])
    return sorted(g for g, s in species.items() if len(s) >= size)


def genus_of(species):
    """The genus of a GTDB species name: s__Bacillus_A cereus -> Bacillus_A."""
    return species[3:].split(" ")[0] if species.startswith("s__") else ""


def read_novel(path):
    """species -> (rank it was held out at, clade) from a list of species, optionally with rank and clade."""
    novel = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if not fields[0] or fields[0].startswith("#"):
                continue
            species = fields[0] if fields[0].startswith("s__") else "s__" + fields[0]
            rank = fields[1] if len(fields) > 1 and fields[1] else "species"
            novel[species] = (rank, fields[2] if len(fields) > 2 and fields[2] else species)
    return novel


def species_lineages(genome_table):
    """species name -> lineage ({rank: name}) from the lineages in the genome table."""
    out = {}
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                out[lineage.split(";")[-1].strip()] = lineages.from_string(lineage)
    return out


def novel_clades(novel, genome_table, seed):
    """{rank: [clade, ...]} of the held-out clades (ranks above species) the genome table can simulate, each
    rank's in a random order (design point i takes clade i modulo their number)."""
    simulated = collections.defaultdict(set)
    for species in species_lineages(genome_table):
        if species in novel and novel[species][0] != "species":
            simulated[novel[species][0]].add(novel[species][1])
    return {rank: random.Random(f"{seed}:{rank}").sample(sorted(clades), len(clades)) for rank, clades in simulated.items()}


def relation(taxon_lineage, in_sample, novel):
    """(meta_relative_rank, meta_novel_level) of a taxon with this lineage among the species simulated in a
    sample ({species: lineage})."""
    depth = {r: i for i, r in enumerate(lineages.RANKS)}
    best_any, best_present, best_novel, level = -1, -1, -1, ""
    for species, lineage in in_sample.items():
        shared = lineages.shared_rank(taxon_lineage, lineage)
        d = depth[shared] if shared else -1
        best_any = max(best_any, d)
        if species in novel:
            if d > best_novel:
                best_novel, level = d, novel[species][0]
        else:
            best_present = max(best_present, d)
    rank = lineages.RANKS[best_any] if best_any >= 0 else "none"
    return rank, (level if best_novel >= 0 and best_novel >= best_present else "")


def representatives(taxonomy):
    """species name -> accession of its representative genome, from internal_taxonomy.dmp."""
    with open(taxonomy) as fh:
        header = next(fh).rstrip("\n").split("\t")
        name, rank, rep = header.index("name"), header.index("rank"), header.index("rep_genome")
        return {f[name]: f[rep] for f in (line.rstrip("\n").split("\t") for line in fh) if f[rank] == "species"}


def manifest_rows(point_dir):
    """The rows of a paired-end point's manifest (one per sample and genome), as dicts."""
    with open(os.path.join(point_dir, "sim", "manifest.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def simulated_genomes(point_dir):
    """sample -> {species name: [genome, ...]} from the simulator's manifest."""
    genomes = {}
    for row in manifest_rows(point_dir):
        species = row["taxonomy"].split(";")[-1]
        genomes.setdefault(row["sample"], {}).setdefault(species, []).append(row["genome"])
    return genomes


def species_domains(genome_table):
    """GTDB species name (s__Genus species, as protal names taxa) -> domain, from the lineages in the table."""
    domains = {}
    with open(genome_table) as fh:
        for line in fh:
            for field in line.rstrip("\n").split("\t"):
                if field.startswith("d__") and ";s__" in field:
                    ranks = field.split(";")
                    domains[ranks[-1].strip()] = ranks[0][3:]
                    break
    return domains


def run(command, log):
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"{command[0]} failed with exit code {rc}; see {log}")


# ---- design points ----------------------------------------------------------------------------------------
# A unit is what one table row set comes from: a design point of one read type. pe and se units share a
# paired-end point (se profiles its first reads); pb and ont units have points of their own, whose samples
# replay the communities of a paired-end point.

def design_points(opts):
    """Paired-end design points: read setups x depths. A setup's ART profile is a built-in one (-ss) or
    file=R1.txt+R2.txt (or file=P.txt for both reads): quality profiles art_profiler_illumina made from real
    reads, which ART uses instead of the built-in one."""
    setups = [s.split(":") for s in opts.read_setups.split(",")]
    if any(len(s) != 4 for s in setups):
        sys.exit(f"--read_setups {opts.read_setups!r}: expected LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, comma-separated")
    depths = opts.read_pairs.split(",")
    for what, values in (("--read_setups", [":".join(s) for s in setups]), ("--read_pairs", depths)):
        twice = sorted(v for v, n in collections.Counter(values).items() if n > 1)
        if twice:  # their points would share a folder
            sys.exit(f"{what} lists {', '.join(twice)} more than once")
    # A point's name tells its setup from the others of its read length (profile, then fragment size).
    lengths = collections.Counter(s[0] for s in setups)
    profiles = collections.Counter((s[0], s[1]) for s in setups)
    points = []
    for i, (length, profile, fragment_mean, fragment_sd) in enumerate(setups):
        tag = "" if lengths[length] == 1 else "_" + (f"custom{i}" if profile.startswith("file=") else profile)
        if profiles[(length, profile)] > 1 and not profile.startswith("file="):
            tag += f"_f{fragment_mean}-{fragment_sd}"
        for pairs in depths:
            points.append({"name": f"rl{length}{tag}_p{pairs}", "read_length": length, "sequencer": profile,
                           "fragment_mean": fragment_mean, "fragment_sd": fragment_sd, "read_pairs": pairs})
    return points


def art_profile_args(profile):
    """simulate_metagenomes options for a read setup's ART profile (see design_points)."""
    if not profile.startswith("file="):
        return ["--sequencer", profile]
    files = profile[len("file="):].split("+")
    for f in files:
        if not os.path.isfile(f):
            sys.exit(f"ART quality profile {f} not found (read setup {profile})")
    return ["--sequencer", "HS25", "--extra_art_args", f"-1 {files[0]} -2 {files[-1]}"]


def parse_long_setup(text):
    """pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL]; the last, for qshmm only, is
    pbsim3's --difference-ratio (it recommends 39/24/36 for ONT, 22/45/33 for Sequel)."""
    parts = text.split(":")
    if len(parts) not in (5, 6) or parts[0] not in ("qshmm", "errhmm") or (len(parts) == 6 and parts[0] != "qshmm"):
        sys.exit(f"long-read setup {text!r}: expected METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN, METHOD "
                 "qshmm or errhmm, and for qshmm optionally :SUB/INS/DEL")
    return {"method": parts[0], "model": parts[1], "length_mean": int(float(parts[2])),
            "length_sd": int(float(parts[3])), "accuracy": float(parts[4]),
            "ratio": parts[5].replace("/", ":") if len(parts) == 6 else ""}


def units_of(opts):
    """The units to collect, in table order."""
    pe_points = design_points(opts)
    units = []
    for read_type in opts.read_types:
        if read_type in ("pe", "se"):
            units += [{"type": read_type, "point": p, "name": p["name"] + ("_se" if read_type == "se" else "")}
                      for p in pe_points]
        else:
            setup = parse_long_setup(opts.pb_setup if read_type == "pb" else opts.ont_setup)
            long_depths = opts.long_read_bases.split(",")
            if len(set(long_depths)) != len(long_depths):  # their points would share a folder
                sys.exit(f"--long_read_bases {opts.long_read_bases!r} lists a depth more than once")
            for i, bases in enumerate(long_depths):
                name = f"{read_type}_b{bases}"
                units.append({"type": read_type, "name": name, "setup": setup, "bases": int(float(bases)),
                              "community": pe_points[i % len(pe_points)],
                              "point": {"name": name, "read_length": str(setup["length_mean"]), "read_pairs": bases}})
    return pe_points, units


def point_dirs(point, opts):
    base = os.path.join(opts.out, "points", point["name"])
    return base, os.path.join(base, "sim"), os.path.join(base, "protal")


def profile_dir(unit, opts):
    """Where protal writes a unit's SAMs, profiles and dumps."""
    base = point_dirs(unit["point"], opts)[0]
    return os.path.join(base, "protal_se" if unit["type"] == "se" else "protal")


def dumps_of(unit, opts):
    return sorted(glob.glob(os.path.join(profile_dir(unit, opts), "**", "*.truth_annotated"), recursive=True))


# ---- simulation -------------------------------------------------------------------------------------------

def abundance_args(text):
    """simulate_metagenomes options of an --abundance model."""
    if not text:
        return []
    parts = text.split(":")
    if parts[0] == "lognormal" and len(parts) == 2:
        return ["--distribution", "poisson_lognormal", "--pln_sigma", parts[1]]
    if parts[0] == "powerlaw" and len(parts) == 2:
        return ["--distribution", "power_law", "--alpha", parts[1]]
    if parts[0] == "negbin" and len(parts) == 3:
        return ["--distribution", "negative_binomial", "--nb_r", parts[1], "--nb_p", parts[2]]
    sys.exit(f"--abundance {text!r}: expected lognormal:SIGMA, powerlaw:ALPHA or negbin:R:P")


# ---- what a point was made from ---------------------------------------------------------------------------
# A rerun reuses a point's samples and dumps only if they were made from the same inputs: each point's folder
# holds the key of what simulated it (simulated.json) and each profile folder that of what profiled it
# (profiled.json; profiling.json while a run profiles it, so that a stopped run resumes). Another database or
# protal (a training database rebuilt with other species held out), seed, genome table or design makes them
# again instead of mixing old dumps into the table.

def identity(path):
    """A file as part of a key: its real path, size and modification time (a rebuild changes them)."""
    real = os.path.realpath(shutil.which(path) or path)
    st = os.stat(real)
    return [real, st.st_size, st.st_mtime_ns]


def db_identity(db):
    if os.path.isfile(db):
        return [identity(db)]
    return [identity(os.path.join(db, n)) for n in sorted(os.listdir(db)) if os.path.isfile(os.path.join(db, n))]


def content_hash(path):
    with open(path, "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()


def same_key(path, key):
    try:
        with open(path) as fh:
            return json.load(fh) == json.loads(json.dumps(key))
    except (OSError, ValueError):
        return False


def write_key(path, key):
    with open(path + ".partial", "w") as fh:
        json.dump(key, fh, indent=1)
    os.replace(path + ".partial", path)


def simulation_key(point, index, opts, clades):
    """What a paired-end point's samples are made from: the simulator's command (but its threads), the genome
    table's content and the simulator binary."""
    command, _ = simulation_command(point, index, opts, 1, clades)
    at = command.index("-t")
    return {"command": command[:at] + command[at + 2:], "genome_table": content_hash(opts.genome_table),
            "simulator": identity(opts.simulator)}


def simulation_command(point, index, opts, threads, clades):
    """The simulator's command for a design point, and None, or why there is none."""
    base, sim, profiles = point_dirs(point, opts)
    command = [opts.simulator, "--genome_table", opts.genome_table, "-o", sim, "-n", str(opts.samples),
               "--sample_prefix", point["name"] + "_s", "--total_read_pairs", point["read_pairs"],
               "--species_per_sample", opts.species_per_sample, "--read_length", point["read_length"],
               *art_profile_args(point["sequencer"]), "--fragment_mean", point["fragment_mean"],
               "--fragment_stdev", point["fragment_sd"], "--seed", str(opts.seed + index),
               "-t", str(threads), "--protal_metafile", profiles, *abundance_args(opts.abundance)]
    if opts.strains_per_species:
        command += ["--strains_per_species", opts.strains_per_species]
    # One --taxon for all demands: the simulator reads only the last.
    taxa = [f"d__Archaea:{opts.archaea}"] if opts.archaea > 0 else []
    if opts.novel_clades > 0:
        taxa += [f"{clades[rank][index % len(clades[rank])]}:{opts.novel_clades}" for rank in sorted(clades)]
    if taxa:
        command += ["--taxon", ",".join(taxa)]
    if opts.congeners > 0:
        genera = large_genera(opts.genome_table, opts.congeners)
        if not genera:
            return command, f"no genus in {opts.genome_table} has {opts.congeners} species (--congeners)"
        command += ["--genus", "g__" + random.Random(opts.seed * 1000 + index).choice(genera) + f":{opts.congeners}"]
    if taxa or opts.congeners > 0:
        command += ["--pick_random_demand_if_fail"]
    return command, None


def simulate(point, index, opts, threads, clades, key):
    """Simulates the samples of a design point; None, or why it failed. clades: {rank: [held-out clade, ...]},
    of which one per rank goes into every sample (--novel_clades)."""
    base, sim, _ = point_dirs(point, opts)
    command, error = simulation_command(point, index, opts, threads, clades)
    if error:
        return error
    os.makedirs(sim, exist_ok=True)
    log = os.path.join(base, "simulate.log")
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        shutil.rmtree(sim, ignore_errors=True)  # no protal.meta: simulated again on a rerun
        return f"{point['name']}: {opts.simulator} failed with exit code {rc}; see {log}"
    write_key(os.path.join(base, "simulated.json"), key)
    return None


def pbsim_model(opts, name):
    """The path of a pbsim3 model: a file, or a name found in --pbsim_models or in pbsim3's data folder."""
    if os.path.isfile(name):
        return name
    folders = [opts.pbsim_models] if opts.pbsim_models else []
    exe = shutil.which(opts.pbsim)
    if exe:
        prefix = os.path.dirname(os.path.dirname(os.path.realpath(exe)))
        folders += sorted(glob.glob(os.path.join(prefix, "share", "pbsim*", "data"))) + \
            sorted(glob.glob(os.path.join(prefix, "share", "pbsim*"))) + [os.path.join(prefix, "data")]
    for folder in folders:
        for candidate in (name, name + ".model"):
            if os.path.isfile(os.path.join(folder, candidate)):
                return os.path.join(folder, candidate)
    sys.exit(f"pbsim3 model {name} not found in {', '.join(folders) or 'no folder'} (--pbsim_models)")


def long_read_genome(task):
    """pbsim3 reads of one genome of a long-read sample (task: dict); the FASTQ files it wrote."""
    tmp = task["tmp"]
    os.makedirs(tmp, exist_ok=True)
    fasta = task["fasta"]
    if fasta.endswith(".gz"):
        plain = os.path.join(tmp, "genome.fna")
        with gzip.open(fasta, "rb") as fin, open(plain, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        fasta = plain
    setup = task["setup"]
    prefix = os.path.join(tmp, "r")
    command = [task["pbsim"], "--strategy", "wgs", "--method", setup["method"], f"--{setup['method']}", task["model"],
               "--genome", fasta, "--depth", f"{task['depth']:.6g}", "--length-mean", str(setup["length_mean"]),
               "--length-sd", str(setup["length_sd"]), "--accuracy-mean", str(setup["accuracy"]),
               "--seed", str(task["seed"]), "--prefix", prefix, "--id-prefix", task["id_prefix"]]
    if setup.get("ratio"):
        command += ["--difference-ratio", setup["ratio"]]
    with open(os.path.join(tmp, "pbsim.log"), "w") as log:
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        return None, f"pbsim failed ({rc}) for {task['genome']}; see {os.path.join(tmp, 'pbsim.log')}"
    return sorted(glob.glob(prefix + "_*.fastq*") + glob.glob(prefix + "_*.fq*")), None


def simulate_long(unit, index, opts, jobs):
    """Long reads (pb, ont) of a design point with pbsim3: the communities of its paired-end point (from its
    manifest), each genome given its share of the sample's bases, by relative abundance times length. Writes
    sim/samples.tsv (sample, reads, truth, community sample) last, so an interrupted point is simulated again."""
    base, sim, _ = point_dirs(unit["point"], opts)
    community_dir = point_dirs(unit["community"], opts)[0]
    _, pe_rows, _ = map_rows(os.path.join(community_dir, "sim", "protal.meta"))
    truth = {row["SAMPLEID"]: row["PROFILE_TRUTH"] for row in pe_rows}
    by_sample = collections.OrderedDict()
    for row in manifest_rows(community_dir):
        by_sample.setdefault(row["sample"], []).append(row)
    model = pbsim_model(opts, unit["setup"]["model"])
    tasks, samples = [], []
    for s, (community, genomes) in enumerate(by_sample.items()):
        sample = f"{unit['name']}_s_{community.rsplit('_', 1)[-1]}"
        weight = [float(g["relative_abundance"]) * float(g["genome_length"]) for g in genomes]
        total = sum(weight) or 1.0
        samples.append((sample, community, len(genomes)))
        for g, genome in enumerate(genomes):
            bases = unit["bases"] * weight[g] / total
            tasks.append({"sample": sample, "genome": genome["genome"], "fasta": genome["fasta_path"],
                          "depth": bases / max(1.0, float(genome["genome_length"])), "setup": unit["setup"],
                          "model": model, "pbsim": opts.pbsim, "seed": (opts.seed * 1000003 + index * 1009 + s) * 101 + g,
                          "id_prefix": f"g{g}x", "tmp": os.path.join(sim, "tmp", sample, str(g))})
    with concurrent.futures.ThreadPoolExecutor(max(1, jobs)) as executor:
        results = list(executor.map(long_read_genome, tasks))
    failures = [error for _, error in results if error]
    if failures:
        return "\n".join(failures[:5])
    reads = os.path.join(sim, "reads")
    os.makedirs(reads, exist_ok=True)
    files = collections.defaultdict(list)
    for task, (fastqs, _) in zip(tasks, results):
        files[task["sample"]] += fastqs
    rows = []
    for sample, community, _ in samples:
        out = os.path.join(reads, sample + ".fq.gz")
        with gzip.open(out, "wb", compresslevel=1) as fout:
            for path in files[sample]:
                with (gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")) as fin:
                    shutil.copyfileobj(fin, fout, 1 << 22)
        rows.append((sample, out, truth[community], community))
    shutil.rmtree(os.path.join(sim, "tmp"), ignore_errors=True)
    with open(os.path.join(sim, "samples.tsv.partial"), "w") as fh:
        fh.write("sample\treads\ttruth\tcommunity\n" + "".join("\t".join(r) + "\n" for r in rows))
    os.replace(os.path.join(sim, "samples.tsv.partial"), os.path.join(sim, "samples.tsv"))
    return None


def simulated(unit, opts):
    """Whether a unit's reads are there: the paired-end map (pe, se), or the long-read samples table."""
    sim = point_dirs(unit["point"], opts)[1]
    return os.path.isfile(os.path.join(sim, "samples.tsv" if unit["type"] in LONG_READ_TYPES else "protal.meta"))


# ---- profiling --------------------------------------------------------------------------------------------

def map_rows(meta):
    """The samples of a simulator map (protal.meta) with every path made absolute, as protal resolves them."""
    dirs, header, rows = {}, None, []
    with open(meta) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#SAMPLEID":
                header = [f.lstrip("#") for f in fields]
            elif fields[0].startswith("#") and len(fields) > 1:
                dirs[fields[0]] = fields[1]
            elif line.strip() and header:
                rows.append(dict(zip(header, fields)))
    out = dirs["#OUTPUT_DIR"]
    where = {"FIRST": dirs.get("#INPUT_DIR", ""), "SECOND": dirs.get("#INPUT_DIR", ""), "PREFIX": out,
             "SAM": os.path.join(out, dirs.get("#SAM_OUTPUT_DIR", "alignments")),
             "PROFILE": os.path.join(out, dirs.get("#PROFILE_OUTPUT_DIR", "profiles"))}
    for row in rows:
        for column, folder in where.items():
            if column in row and row[column] not in ("", "-"):
                row[column] = os.path.join(folder, row[column])
    return header, rows, [where["SAM"], where["PROFILE"]]


def unit_map_rows(unit, opts):
    """The rows of a unit in the combined map (MAP_COLUMNS, absolute paths) and the folders they write to."""
    if unit["type"] in ("pe", "se"):
        _, rows, dirs = map_rows(os.path.join(point_dirs(unit["point"], opts)[1], "protal.meta"))
        if unit["type"] == "pe":
            return [{**{c: row.get(c, "-") for c in MAP_COLUMNS}, "READ_TYPE": "pe"} for row in rows], dirs
        out = profile_dir(unit, opts)
        dirs = [os.path.join(out, "alignments"), os.path.join(out, "profiles")]
        return [{"SAMPLEID": row["SAMPLEID"] + "_se", "FIRST": row["FIRST"], "SECOND": "-",
                 "SAM": os.path.join(dirs[0], row["SAMPLEID"] + "_se.sam.gz"),
                 "PREFIX": os.path.join(out, row["SAMPLEID"] + "_se"),
                 "PROFILE": os.path.join(dirs[1], row["SAMPLEID"] + "_se.profile"),
                 "PROFILE_TRUTH": row["PROFILE_TRUTH"], "READ_TYPE": "se"} for row in rows], dirs
    out = profile_dir(unit, opts)
    dirs = [os.path.join(out, "alignments"), os.path.join(out, "profiles")]
    rows = []
    with open(os.path.join(point_dirs(unit["point"], opts)[1], "samples.tsv")) as fh:
        next(fh)
        for sample, reads, truth, _ in (line.rstrip("\n").split("\t") for line in fh if line.strip()):
            rows.append({"SAMPLEID": sample, "FIRST": reads, "SECOND": "-",
                         "SAM": os.path.join(dirs[0], sample + ".sam.gz"), "PREFIX": os.path.join(out, sample),
                         "PROFILE": os.path.join(dirs[1], sample + ".profile"), "PROFILE_TRUTH": truth,
                         "READ_TYPE": unit["type"]})
    return rows, dirs


def profile(units, opts):
    """Profiles the samples of all units in one protal run: the database is loaded once, and every sample is
    profiled as its READ_TYPE says."""
    folder = os.path.join(opts.out, "profile_all")
    os.makedirs(folder, exist_ok=True)
    rows = []
    for unit in units:
        unit_rows, dirs = unit_map_rows(unit, opts)
        for d in dirs:
            os.makedirs(d, exist_ok=True)
        rows += unit_rows
    combined = os.path.join(folder, "samples.map")
    with open(combined, "w") as fh:
        fh.write(f"#OUTPUT_DIR\t{folder}\n#" + "\t".join(MAP_COLUMNS) + "\n")
        fh.writelines("\t".join(row[c] for c in MAP_COLUMNS) + "\n" for row in rows)
    kinds = collections.Counter(row["READ_TYPE"] for row in rows)
    print(f"profiling {len(rows)} samples ({', '.join(f'{n} {t}' for t, n in kinds.items())}) of {len(units)} design "
          "points in one protal run", flush=True)
    run([opts.protal, "--db", opts.db, "--map", combined, "-t", str(opts.threads), "--no_strains", "--no_qcmsa"],
        os.path.join(folder, "protal.log"))


# ---- the tables -------------------------------------------------------------------------------------------

def community_of(unit, sample, opts):
    """The paired-end point and sample whose community a unit's sample holds."""
    if unit["type"] == "pe":
        return unit["point"], sample
    if unit["type"] == "se":
        return unit["point"], sample[:-len("_se")]
    with open(os.path.join(point_dirs(unit["point"], opts)[1], "samples.tsv")) as fh:
        next(fh)
        for name, _, _, community in (line.rstrip("\n").split("\t") for line in fh if line.strip()):
            if name == sample:
                return unit["community"], community
    sys.exit(f"{sample}: not in the samples of {unit['name']}")


def write_table(read_type, units, opts, context):
    """Joins the dumps of a read type's units into its table."""
    domains, novel, reps, db_lineages, sim_lineages = context
    header, rows = None, 0
    table = os.path.join(opts.out, TABLES[read_type])
    genomes_of = {}
    with open(table + ".partial", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for unit in units:
            point = unit["point"]
            present = absent = strains = congeners = 0
            attributed = collections.Counter()
            for dump in dumps_of(unit, opts):
                sample = os.path.basename(dump).split(".profile")[0]
                community_point, community = community_of(unit, sample, opts)
                if community_point["name"] not in genomes_of:
                    genomes_of[community_point["name"]] = simulated_genomes(point_dirs(community_point, opts)[0])
                in_sample = genomes_of[community_point["name"]].get(community, {})
                novel_here = [s for s in in_sample if s in novel]
                novel_genera = {genus_of(s) for s in novel_here}
                levels = collections.Counter(novel[s][0] for s in novel_here)
                novel_levels = ",".join(f"{rank}:{levels[rank]}" for rank in ("species", *reversed(lineages.RANKS[1:-1]))
                                        if levels.get(rank))
                sample_lineages = {s: sim_lineages.get(s, {}) for s in in_sample}
                with open(dump, newline="") as fh:
                    reader = csv.reader(fh, delimiter="\t")
                    dump_header = next(reader)
                    if header is None:
                        header = dump_header
                        writer.writerow(META_COLUMNS + header)
                    elif dump_header != header:
                        sys.exit(f"{dump} has other columns than the dumps before it (another protal version?)")
                    name, truth = header.index("taxon_name"), header.index("truth")
                    for row in reader:
                        taxon = row[name]
                        is_present = row[truth].lower() in ("1", "true")
                        congener = genus_of(taxon) in novel_genera
                        rep = ""
                        if is_present and reps.get(taxon) and taxon in in_sample:
                            rep = "1" if all(g == reps[taxon] for g in in_sample[taxon]) else "0"
                        relative, level, neighbour = "", "", ""
                        if is_present:
                            relative = "species"
                            others = {s: lin for s, lin in sample_lineages.items() if s != taxon}
                            neighbour = relation(db_lineages.get(taxon) or sim_lineages.get(taxon, {}), others, {})[0]
                        elif taxon in db_lineages:
                            relative, level = relation(db_lineages[taxon], sample_lineages, novel)
                        writer.writerow([unit["name"], sample, point["read_length"], point["read_pairs"],
                                         domains.get(taxon, "unknown"), len(novel_here), int(congener), rep,
                                         novel_levels, level, relative, neighbour, read_type] + row)
                        rows += 1
                        present += is_present
                        absent += not is_present
                        strains += rep == "0"
                        congeners += congener and not is_present
                        attributed[level] += bool(level)
            by_level = ", ".join(f"{n} to {rank}" for rank, n in sorted(attributed.items()) if rank)
            print(f"{unit['name']}: {present} present ({strains} from other genomes than the representative) and "
                  f"{absent} absent taxa ({congeners} congeners of species the database lacks; closest to a species "
                  f"it lacks, by the rank held out: {by_level or 'none'}) in {len(dumps_of(unit, opts))} samples",
                  flush=True)
    os.replace(table + ".partial", table)
    print(f"{rows} taxa in {table}", flush=True)


def main(argv=None):
    opts = parse_args(argv)
    os.makedirs(opts.out, exist_ok=True)
    domains = species_domains(opts.genome_table)
    novel = read_novel(opts.novel_species) if opts.novel_species else {}
    reps = representatives(opts.taxonomy) if opts.taxonomy else {}
    # Lineages of the database's taxa (by name) and of the simulated species, for meta_relative_rank.
    db_lineages = {}
    if opts.taxonomy:
        by_id, _ = lineages.from_taxonomy(opts.taxonomy)
        db_lineages = {lin[max(lin, key=lineages.RANKS.index)]: lin for lin in by_id.values() if lin}
    sim_lineages = species_lineages(opts.genome_table)
    clades = novel_clades(novel, opts.genome_table, opts.seed) if opts.novel_clades > 0 else {}
    if clades:
        print("held-out clades in every sample, one per rank and design point: "
              + ", ".join(f"{len(c)} {rank}" for rank, c in sorted(clades.items())), flush=True)
    pe_points, units = units_of(opts)
    jobs = max(1, opts.jobs or opts.threads)

    # Paired-end simulation, which every read type needs (se reads it, pb and ont replay its communities): ART
    # simulates one genome at a time, so design points run in parallel.
    needed = {u["point"]["name"] for u in units if u["type"] in ("pe", "se")} | \
             {u["community"]["name"] for u in units if u["type"] in LONG_READ_TYPES}
    keys = {p["name"]: simulation_key(p, i, opts, clades) for i, p in enumerate(pe_points) if p["name"] in needed}
    pending = []
    for i, p in enumerate(pe_points):
        base, sim, _ = point_dirs(p, opts)
        if p["name"] not in needed:
            continue
        if os.path.isfile(os.path.join(sim, "protal.meta")):
            if same_key(os.path.join(base, "simulated.json"), keys[p["name"]]):
                continue
            print(f"{p['name']} was simulated from other inputs (or by an older collector): simulating it again", flush=True)
            shutil.rmtree(base)
        pending.append((i, p))
    if pending:
        workers = max(1, min(jobs, len(pending)))
        threads = max(1, opts.threads // workers)
        print(f"simulating {len(pending)} paired-end design points, {workers} at a time", flush=True)
        with concurrent.futures.ThreadPoolExecutor(workers) as executor:
            failures = [f for f in executor.map(lambda ip: simulate(ip[1], ip[0], opts, threads, clades, keys[ip[1]["name"]]),
                                                pending) if f]
        if failures:
            sys.exit("\n".join(failures))
    # Long reads: pbsim3 simulates one genome at a time, so the genomes of a point run in parallel.
    for i, unit in enumerate(u for u in units if u["type"] in LONG_READ_TYPES):
        base = point_dirs(unit["point"], opts)[0]
        keys[unit["name"]] = {"community": keys[unit["community"]["name"]], "setup": unit["setup"], "bases": unit["bases"],
                              "index": i, "seed": opts.seed, "pbsim": identity(opts.pbsim),
                              "model": identity(pbsim_model(opts, unit["setup"]["model"]))}
        if simulated(unit, opts) and not same_key(os.path.join(base, "simulated.json"), keys[unit["name"]]):
            print(f"{unit['name']} was simulated from other inputs (or by an older collector): simulating it again", flush=True)
            shutil.rmtree(base)
        if not simulated(unit, opts):
            print(f"simulating {unit['name']} with pbsim3", flush=True)
            failure = simulate_long(unit, i, opts, jobs)
            if failure:
                sys.exit(f"{unit['name']}: {failure}")
            write_key(os.path.join(base, "simulated.json"), keys[unit["name"]])
    # Profiling: every unit not yet profiled against this database with this protal, in one protal run.
    db, protal = db_identity(opts.db), identity(opts.protal)
    unprofiled = []
    for unit in units:
        folder = profile_dir(unit, opts)
        key = {"simulated": keys[unit["point"]["name"] if unit["type"] in ("pe", "se") else unit["name"]],
               "db": db, "protal": protal, "read_type": unit["type"]}
        if same_key(os.path.join(folder, "profiled.json"), key) and len(dumps_of(unit, opts)) >= opts.samples:
            continue
        if not same_key(os.path.join(folder, "profiling.json"), key):  # else a stopped run's: go on with it
            if dumps_of(unit, opts):
                print(f"{unit['name']} was profiled against another database or with another protal (or by an older "
                      "collector): profiling it again", flush=True)
            shutil.rmtree(folder, ignore_errors=True)
            os.makedirs(folder)
            write_key(os.path.join(folder, "profiling.json"), key)
        unprofiled.append(unit)
    if unprofiled:
        profile(unprofiled, opts)
        for unit in unprofiled:
            folder = profile_dir(unit, opts)
            os.replace(os.path.join(folder, "profiling.json"), os.path.join(folder, "profiled.json"))
    for unit in units:
        if len(dumps_of(unit, opts)) != opts.samples:
            sys.exit(f"{unit['name']}: expected {opts.samples} training dumps in {profile_dir(unit, opts)}, "
                     f"found {len(dumps_of(unit, opts))}")

    context = (domains, novel, reps, db_lineages, sim_lineages)
    for read_type in opts.read_types:
        write_table(read_type, [u for u in units if u["type"] == read_type], opts, context)


if __name__ == "__main__":
    main()
