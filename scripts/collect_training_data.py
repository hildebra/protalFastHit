#!/usr/bin/env python3
"""Simulate metagenomes and profile them with protal, to train its presence models.

For every design point, a read setup (length, ART profile, fragment size) and a sequencing depth,
simulate_metagenomes draws random communities from a genome table and protal profiles them
against a database, knowing the true species. The training dumps of all samples
(<profile>.truth_annotated: every taxon protal saw, its features and whether it was present) are
joined into one table, with meta_* columns saying where each row comes from (meta_domain from
the genome table's lineages; with --novel_species and --taxonomy, which rows share a genus with a
species the database lacks, and which present species were simulated from another genome than
the database's reference, and meta_insilico_strain which of them from an in-silico strain of
insilico_strains.py). Train on it with

    python3 scripts/random_forest_cmdline.py --truth-file OUT/training_data.tsv --output-prefix OUT/model

Read types (--read_types): pe, the paired-end samples above; se, the same samples' first reads
alone, profiled as single-end reads; pb and ont, long reads of the same communities (--pb_setup,
--ont_setup), --long_read_bases per sample, one design point each. The collector draws a long-read
sample's reads (each read's genome by abundance times length, its length from the setup's gamma
distribution, its start uniform, cut where its contig ends); PacBio HiFi reads are made of them by
hifi_reads.py (errors mostly in homopolymers, calibrated qualities, Q50 for 5 kb reads to Q30 for 25 kb
and Q20 for 50 kb: pbsim3 simulates no HiFi reads, and its one-pass reads have quality 0 throughout),
Nanopore reads by pbsim3 with its quality model (--strategy templ), in one run per sample. Every read type gets its table:
training_data.tsv (pe), training_data_se.tsv, training_data_pb.tsv, training_data_ont.tsv; protal
profiles each sample with that read type's model and settings.

A genome table with species the database lacks gives the negatives that matter most: a relative
the database has picks up their reads. When whole genera, families, orders, classes or phyla are
missing (--novel_species with ranks, --novel_clades), the relatives are distant: meta_novel_level
marks the absent taxa closest to such species, meta_neighbour_rank how close a present taxon's
nearest other species in the sample is. Archaea (--archaea) have fewer marker genes than bacteria
and need to be in the training data, too. All simulations share --jobs cores in one queue: the
paired-end design points with threads in proportion to their work, the long-read samples (deep ones
in chunks, --long_read_chunk) on the rest, the longest first, as soon as a design run (the same
communities without reads, in seconds) has given their communities. Then the samples of all read
types are profiled in one protal run, which loads the database once (--prepare_profiling and
--also_profile put two collections into one run). Points already simulated or profiled are skipped,
so a run can be resumed; --simulate_only stops before profiling, so that the simulations can run
while the database is built.

Scenarios (--scenarios, scenarios.py): samples like those of one kind of study (gut, soil, shallow soil,
host-dominated samples, or one defined in --scenario_file), simulated and profiled beside the design's
and joined into the same tables, with meta_scenario naming the scenario. A scenario's communities are
drawn from a genome table of its own that gives the share of species the database lacks, and each of
its read types (Illumina at a mean base quality, Ultima single-end reads, PacBio, Nanopore; a host
genome's share of the reads, --host_genome) sequences the same communities. --samples 0 collects the
scenarios alone.

usage: collect_training_data.py --db DB --genome_table genomes.tsv -o OUT [options]
"""

import argparse
import bisect
import collections
import concurrent.futures
import csv
import glob
import gzip
import hashlib
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lineages  # noqa: E402
import scenarios  # noqa: E402

READ_TYPES = ("pe", "se", "pb", "ont")
LONG_READ_TYPES = ("pb", "ont")
# The folders of the scenarios' points start with this (check_model_parity.py and trace_relatives.py leave them out:
# their samples are deep, and the design's points judge what those check).
SCENARIO_PREFIX = "sc_"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="protal database")
    p.add_argument("--genome_table", required=True,
                   help="simulate_metagenomes genome table (accession, GTDB taxonomy, FASTA path)")
    p.add_argument("-o", "--out", required=True, help="output directory")
    p.add_argument("--protal", default="protal", help="protal binary (default: protal on PATH)")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes binary")
    p.add_argument("--samples", type=int, default=4,
                   help="samples per design point (default: 4; 0: no design points, only the scenarios)")
    p.add_argument("--read_pairs", default="1000,5000,20000,100000,500000",
                   help="comma-separated read pairs per sample, one design point each; DEPTH:SAMPLES gives a point other "
                        "samples than --samples (e.g. 10000000:2 for deep samples)")
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
                        "Poisson-lognormal with sigma 1.3); lognormal:S1,S2,... gives a design point's samples the "
                        "sigmas in turn, so that a model does not learn one sigma's prior")
    p.add_argument("--archaea", type=int, default=0, help="archaeal species per sample (default: 0)")
    p.add_argument("--congeners", default="0", type=congener_spec,
                   help="relatives that share a sample. SHARE:MIN-MAX (e.g. 0.25:2-5): about SHARE of each sample's "
                        "species come in groups of MIN to MAX species of one genus, the genera drawn per sample "
                        "(simulate_metagenomes --congener_groups); N: N species of one genus in every sample of a "
                        "design point, the genus drawn per point among those with that many species; 0 (default): "
                        "none. Species are otherwise drawn uniformly, so among many genera relatives hardly ever "
                        "share a sample, while in real samples they often do, at very different abundances")
    p.add_argument("--read_types", default="pe",
                   help="comma-separated read types to collect: pe, se, pb, ont (default pe)")
    p.add_argument("--long_read_bases", default="300000,1500000,6000000,30000000,150000000",
                   help="bases per long-read sample (pb, ont), one design point each, DEPTH:SAMPLES for other samples "
                        "than --long_read_samples; point i replays the communities of the paired-end points of the i-th "
                        "--read_pairs depth, of every read setup in turn (default: about the paired-end points' bases)")
    p.add_argument("--long_read_samples", type=int, default=0,
                   help="samples per long-read design point (default: --samples); at most the samples of the paired-end "
                        "points whose communities it replays")
    p.add_argument("--pb_setup", default="hifi:15000:3000:3",
                   help="PacBio reads: hifi:LENGTH_MEAN:LENGTH_SD:Q_SD, HiFi reads by hifi_reads.py, their quality by "
                        "their length and Q_SD around it (default), or a pbsim3 setup as --ont_setup's")
    p.add_argument("--ont_setup", default="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of Nanopore reads")
    p.add_argument("--pbsim", default="pbsim", help="pbsim3 binary (default: pbsim on PATH)")
    p.add_argument("--pbsim_models", help="folder of pbsim3's .model files (default: found next to the binary)")
    p.add_argument("-t", "--threads", type=int, default=4, help="threads of the protal run (default 4)")
    p.add_argument("--jobs", type=int, default=0,
                   help="cores the simulations share (default: --threads): the paired-end design points take threads in "
                        "proportion to their work, and the long-read samples, or chunks of them, the rest, the longest first")
    p.add_argument("--long_read_chunk", type=int, default=LONG_READ_CHUNK,
                   help=f"long-read samples of more bases are simulated in chunks of at most this many bases side by "
                        f"side, their reads joined (default {LONG_READ_CHUNK}; 0: one run per sample)")
    p.add_argument("--prepare_profiling", action="store_true",
                   help="simulate as needed, then write the map of the samples to profile (OUT/profile_all/samples.map, "
                        "their folders in samples.map.units) and stop, for a run of another collection to profile them "
                        "with its own (--also_profile)")
    p.add_argument("--also_profile", action="append", default=[],
                   help="the samples.map of another collection (--prepare_profiling), profiled in this collection's protal "
                        "run against --db, so that the database is loaded once; repeatable")
    p.add_argument("--simulate_only", action="store_true",
                   help="simulate the design points and stop: a later run without it profiles them (the simulations "
                        "need no database, so they can run while it is built; --db is not read)")
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
    p.add_argument("--scenarios", default="",
                   help="scenarios to collect besides the design, NAME[:SAMPLES] comma-separated (gut, soil, soil_shallow, "
                        "host, those of --scenario_file; all: every one), each with --scenario_samples samples unless it "
                        "gives its own (scenarios.py)")
    p.add_argument("--scenario_samples", type=int, default=3, help="samples per scenario (default 3)")
    p.add_argument("--scenario_file", help="JSON of scenarios by name, which add to or change the presets")
    p.add_argument("--host_genome", help="FASTA (gzipped or not) of the host genome of scenarios with a host share, "
                                         "e.g. the human genome download_gtdb.py fetches")
    opts = p.parse_args(argv)
    opts.read_types = [t.strip() for t in opts.read_types.split(",") if t.strip()]
    unknown = [t for t in opts.read_types if t not in READ_TYPES]
    if unknown:
        p.error(f"--read_types: unknown {', '.join(unknown)} (pe, se, pb, ont)")
    if opts.samples < 0:
        p.error("--samples cannot be negative")
    try:
        scenarios.selection(opts.scenarios, opts.scenario_samples, scenarios.definitions(opts.scenario_file))
    except scenarios.ScenarioError as e:
        p.error(str(e))
    return opts


META_COLUMNS = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_domain",
                "meta_novel_species", "meta_novel_congener", "meta_rep_genome", "meta_novel_levels",
                "meta_novel_level", "meta_relative_rank", "meta_neighbour_rank", "meta_read_type", "meta_insilico_strain",
                "meta_scenario"]
# meta_read_pairs: the design point's depth, read pairs (pe; reads for se) or bases (pb, ont).
# meta_novel_levels: the sample's species the database lacks, by the rank they were held out at
# ("species:2,family:1"). meta_relative_rank: the deepest rank the taxon shares with a species simulated in
# the sample ("species" for the simulated species themselves, "none" for no shared domain). meta_novel_level:
# for an absent taxon whose closest simulated species (at least as close as any present one) is one the
# database lacks, the rank that species was held out at: its reads are the likely source of the taxon's.
# meta_neighbour_rank: for a present taxon, the deepest rank it shares with another species simulated in the
# sample (a congener's reads fit it nearly as well, so it may be missed).
# meta_insilico_strain: for a present taxon, 1 if a genome it was simulated from is an in-silico strain (a mutated copy
# of a one-genome species' representative, insilico_strains.py; its meta_rep_genome is 0), else 0.
# meta_scenario: the scenario the sample is of (--scenarios), empty for the design's samples.
INSILICO_PREFIX = "insilico_"  # insilico_strains.py's PREFIX: the names of its strains
TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
MAP_COLUMNS = ["SAMPLEID", "FIRST", "SECOND", "SAM", "PREFIX", "PROFILE", "PROFILE_TRUTH", "READ_TYPE"]


def congener_spec(text):
    """--congeners: ("groups", "SHARE:MIN-MAX") for congener groups in every sample, ("genus", N) for N species of one
    genus per design point, or ("genus", 0) for none."""
    text = str(text).strip()
    if ":" in text:
        share, _, sizes = text.partition(":")
        lo, _, hi = sizes.partition("-")
        try:
            ok = 0 <= float(share) <= 1 and lo.isdigit() and hi.isdigit() and 2 <= int(lo) <= int(hi)
        except ValueError:
            ok = False
        if not ok:
            raise argparse.ArgumentTypeError(f"{text!r}: expected SHARE:MIN-MAX (a share 0-1, 2 <= MIN <= MAX), or N")
        return ("groups", f"{float(share):g}:{int(lo)}-{int(hi)}")
    if not text.isdigit():
        raise argparse.ArgumentTypeError(f"{text!r}: expected SHARE:MIN-MAX (e.g. 0.25:2-5) or a number of species")
    return ("genus", int(text))


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


def lineage_prefixes(lineage):
    """The prefixes of a lineage from the domain down, each a tuple of names, as far as it has every rank: two lineages
    share rank r (shared_rank) when they share the prefix of r."""
    out, names = [], []
    for r in lineages.RANKS:
        if r not in lineage:
            break
        names.append(lineage[r])
        out.append(tuple(names))
    return out


class SampleRelations:
    """relation() for every taxon of one sample at once: the sample's simulated species ({species: lineage}, in
    order) by their lineages' prefixes, counted, the database's and the novel ones apart, with the held-out rank of
    the first novel species of each prefix. A taxon then costs its own seven prefixes instead of a comparison with
    every species of the sample (10,000 in a soil scenario)."""

    def __init__(self, in_sample, novel):
        self.count, self.present, self.first_novel, self.own = collections.Counter(), set(), {}, {}
        for species, lineage in in_sample.items():
            prefixes = lineage_prefixes(lineage)
            self.own[species] = set(prefixes)
            for prefix in prefixes:
                self.count[prefix] += 1
                if species in novel:
                    self.first_novel.setdefault(prefix, novel[species][0])
                else:
                    self.present.add(prefix)

    def relation(self, lineage):
        """relation(lineage, in_sample, novel): (meta_relative_rank, meta_novel_level)."""
        best_any = best_present = best_novel = -1
        level = ""
        for depth, prefix in enumerate(lineage_prefixes(lineage)):
            if self.count[prefix]:
                best_any = depth
            if prefix in self.present:
                best_present = depth
            if prefix in self.first_novel:
                best_novel, level = depth, self.first_novel[prefix]
        rank = lineages.RANKS[best_any] if best_any >= 0 else "none"
        return rank, (level if best_novel >= 0 and best_novel >= best_present else "")

    def neighbour(self, lineage, species):
        """relation(lineage, the sample's species but `species`, {})[0]: the deepest rank shared with another one."""
        own = self.own.get(species, set())
        best = -1
        for depth, prefix in enumerate(lineage_prefixes(lineage)):
            if self.count[prefix] - (prefix in own) > 0:
                best = depth
        return lineages.RANKS[best] if best >= 0 else "none"


def representatives(taxonomy):
    """species name -> accession of its representative genome, from internal_taxonomy.dmp."""
    with open(taxonomy) as fh:
        header = next(fh).rstrip("\n").split("\t")
        name, rank, rep = header.index("name"), header.index("rank"), header.index("rep_genome")
        return {f[name]: f[rep] for f in (line.rstrip("\n").split("\t") for line in fh) if f[rank] == "species"}


def manifest_rows(folder):
    """The rows of the manifest.tsv in a paired-end point's sim folder (or its design folder; one row per sample and
    genome), as dicts."""
    with open(os.path.join(folder, "manifest.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def simulated_genomes(point_dir):
    """sample -> {species name: [genome, ...]} from the simulator's manifest."""
    genomes = {}
    for row in manifest_rows(os.path.join(point_dir, "sim")):
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


def clock(seconds):
    """Seconds as h:mm:ss."""
    seconds = int(round(seconds))
    return f"{seconds // 3600}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"


def run_protal(command, log, samples):
    """Runs protal on `samples` samples, saying from its log every minute how far it is (when that changed):
    the sample it aligns, then how many it has profiled. Only what the log gained is read."""
    started = time.time()
    with open(log, "w") as fh:
        process = subprocess.Popen(command, stdout=fh, stderr=subprocess.STDOUT)
    offset, rest, aligning, profiled, told = 0, b"", 0, 0, (0, 0)
    try:
        while True:
            try:
                rc = process.wait(timeout=60)
                break
            except subprocess.TimeoutExpired:
                pass
            with open(log, "rb") as fh:
                fh.seek(offset)
                chunk = fh.read()
            offset += len(chunk)
            lines = (rest + chunk).split(b"\n")
            rest = lines.pop()  # a line protal is still writing
            aligning += sum(line.startswith((b"Align the ", b"Skip ")) for line in lines)
            profiled += sum(line.startswith(b"Write truth to: ") for line in lines)
            if (aligning, profiled) != told:
                told = (aligning, profiled)
                state = f"{profiled} of {samples} samples profiled" if profiled else f"aligning sample {aligning} of {samples}"
                print(f"protal, {clock(time.time() - started)} in: {state}", flush=True)
    except BaseException:  # stopped (Ctrl-C, an error): not without protal
        process.kill()
        process.wait()
        raise
    if rc != 0:
        sys.exit(f"{command[0]} failed with exit code {rc}; see {log}")
    print(f"protal profiled {samples} samples in {clock(time.time() - started)}", flush=True)


# ---- design points ----------------------------------------------------------------------------------------
# A unit is what one table row set comes from: a design point of one read type. pe and se units share a
# paired-end point (se profiles its first reads); pb and ont units have points of their own, whose samples
# replay the communities of a paired-end point.

def depths_and_samples(text, what):
    """A comma-separated list of DEPTH or DEPTH:SAMPLES (--read_pairs, --long_read_bases) -> [(depth as given,
    samples, or None where a depth gives none)]: deep points can have fewer samples than the others."""
    out = []
    for item in text.split(","):
        depth, _, n = item.strip().partition(":")
        try:
            if float(depth) <= 0 or (n and int(n) <= 0):
                raise ValueError
        except ValueError:
            sys.exit(f"{what} {text!r}: expected DEPTH or DEPTH:SAMPLES, comma-separated, both positive")
        out.append((depth, int(n) if n else None))
    return out


def design_points(opts):
    """Paired-end design points: read setups x depths, each depth with its samples (--read_pairs DEPTH[:SAMPLES],
    --samples by default). A setup's ART profile is a built-in one (-ss) or file=R1.txt+R2.txt (or file=P.txt for
    both reads): quality profiles art_profiler_illumina made from real reads, which ART uses instead of the
    built-in one."""
    setups = [s.split(":") for s in opts.read_setups.split(",")]
    if any(len(s) != 4 for s in setups):
        sys.exit(f"--read_setups {opts.read_setups!r}: expected LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, comma-separated")
    depth_samples = [(d, n or opts.samples) for d, n in depths_and_samples(opts.read_pairs, "--read_pairs")]
    depths = [d for d, _ in depth_samples]
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
        for depth_index, (pairs, samples) in enumerate(depth_samples):
            points.append({"name": f"rl{length}{tag}_p{pairs}", "read_length": length, "sequencer": profile,
                           "fragment_mean": fragment_mean, "fragment_sd": fragment_sd, "read_pairs": pairs,
                           "samples": samples, "depth_index": depth_index})
    return points


def art_profile_args(profile, extra=""):
    """simulate_metagenomes options for a read setup's ART profile (see design_points), and `extra` ART options
    (a scenario's quality shifts)."""
    if not profile.startswith("file="):
        return ["--sequencer", profile] + (["--extra_art_args", extra] if extra else [])
    files = profile[len("file="):].split("+")
    for f in files:
        if not os.path.isfile(f):
            sys.exit(f"ART quality profile {f} not found (read setup {profile})")
    return ["--sequencer", "HS25", "--extra_art_args", f"-1 {files[0]} -2 {files[-1]}" + (f" {extra}" if extra else "")]


def art_options(profile):
    """art_illumina's own options for a read setup's ART profile: -ss PROFILE, or -1 R1 -2 R2 of a file= profile."""
    if not profile.startswith("file="):
        return ["-ss", profile]
    files = profile[len("file="):].split("+")
    return ["-1", files[0], "-2", files[-1]]


PBSIM_METHODS = ("qshmm", "errhmm")


def parse_long_setup(text):
    """hifi:LENGTH_MEAN:LENGTH_SD:Q_SD (hifi_reads.py: the SD of the reads' quality around their length's, Phred), or
    pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL]; the last, for qshmm only, is
    pbsim3's --difference-ratio (it recommends 39/24/36 for ONT, 22/45/33 for Sequel). And
    ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD, the Ultima Genomics single-end reads of scenarios (hifi_reads.py's flow
    model: a mean base quality of Q_MEAN per read, SD Q_SD)."""
    parts = text.split(":")
    if parts[0] == "ultima":
        try:
            if len(parts) != 5:
                raise ValueError
            length_mean, length_sd, q_mean, q_sd = (float(v) for v in parts[1:])
        except ValueError:
            sys.exit(f"read setup {text!r}: expected ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD")
        return {"method": "ultima", "model": None, "length_mean": int(length_mean), "length_sd": int(length_sd),
                "q_mean": q_mean, "q_sd": q_sd, "ratio": ""}
    if parts[0] == "hifi":
        try:
            if len(parts) != 4:
                raise ValueError
            length_mean, length_sd, q_sd = (float(v) for v in parts[1:])
        except ValueError:
            sys.exit(f"long-read setup {text!r}: expected hifi:LENGTH_MEAN:LENGTH_SD:Q_SD")
        return {"method": "hifi", "model": None, "length_mean": int(length_mean), "length_sd": int(length_sd),
                "q_sd": q_sd, "ratio": ""}
    if len(parts) not in (5, 6) or parts[0] not in ("qshmm", "errhmm") or (len(parts) == 6 and parts[0] != "qshmm"):
        sys.exit(f"long-read setup {text!r}: expected hifi:LENGTH_MEAN:LENGTH_SD:Q_SD, or "
                 "METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN with METHOD qshmm or errhmm, and for qshmm "
                 "optionally :SUB/INS/DEL")
    return {"method": parts[0], "model": parts[1], "length_mean": int(float(parts[2])),
            "length_sd": int(float(parts[3])), "accuracy": float(parts[4]),
            "ratio": parts[5].replace("/", ":") if len(parts) == 6 else ""}


def drawn(unit):
    """Whether a unit's reads are drawn from templates of its communities' genomes (long reads, a scenario's Ultima
    reads: the unit has a setup and a point of its own), not simulated by ART for a paired-end point (pe, and se of its
    first reads)."""
    return "setup" in unit


def scenario_units(opts):
    """The points and units of --scenarios, after the design's: per scenario a community point, sc_<name>_pe_p<pairs>
    when it has paired-end reads (simulated at the community's part of the read pairs, the host's added later), else
    sc_<name>_community (its communities without reads, simulate_metagenomes --test); its pe unit; and its units of
    reads drawn from those communities (se: Ultima reads, pb, ont), each the scenario's samples. Read types not
    collected (--read_types) are left out. -> (points, units)."""
    if not getattr(opts, "scenarios", ""):
        return [], []
    defs = scenarios.definitions(getattr(opts, "scenario_file", None))
    points, units = [], []
    for name, samples in scenarios.selection(opts.scenarios, opts.scenario_samples, defs):
        d = defs[name]
        reads = {r["type"]: r for r in d["reads"] if r["type"] in opts.read_types}
        if not reads:
            continue
        pe, host = reads.get("pe"), d["host_share"]
        common = {"samples": samples, "depth_index": None, "scenario": name, "definition": d}
        if pe:
            community = max(1, round(pe["depth"] * (1 - host)))
            point = {**common, "name": f"{SCENARIO_PREFIX}{name}_pe_p{pe['depth']}", "read_length": str(pe["length"]),
                     "sequencer": pe["profile"], "fragment_mean": str(pe["fragment_mean"]),
                     "fragment_sd": str(pe["fragment_sd"]), "read_pairs": str(pe["depth"]), "quality": pe.get("quality"),
                     "reads": True, "community_pairs": str(community), "host_pairs": pe["depth"] - community}
        else:  # the read pairs only spread over the genomes, every one with one at least
            point = {**common, "name": f"{SCENARIO_PREFIX}{name}_community", "read_length": "150", "sequencer": "HS25",
                     "fragment_mean": "350", "fragment_sd": "50", "read_pairs": "0", "quality": None, "reads": False,
                     "community_pairs": str(100 * scenarios.species_bounds(d["species"])[1]), "host_pairs": 0}
        points.append(point)
        if pe:
            units.append({"type": "pe", "point": point, "name": point["name"], "samples": samples, "scenario": name})
        for kind in ("se", "pb", "ont"):
            if kind not in reads:
                continue
            given = reads[kind]
            setup = parse_long_setup(given.get("setup") or (opts.pb_setup if kind == "pb" else opts.ont_setup))
            if kind == "se":  # Ultima reads: the depth in reads, drawn as bases of the setup's mean length
                unit_name, bases = f"{SCENARIO_PREFIX}{name}_se_ultima_r{given['depth']}", given["depth"] * setup["length_mean"]
            else:
                unit_name, bases = f"{SCENARIO_PREFIX}{name}_{kind}_b{given['depth']}", given["depth"]
            units.append({"type": kind, "name": unit_name, "setup": setup, "bases": bases, "samples": samples,
                          "communities": [point], "scenario": name, "host_share": host,
                          # its samples' seeds do not depend on the design's long-read points
                          "seed_index": 1_000_000 + scenarios.seed_of(0, f"{name}:{kind}") % 1_000_000,
                          "point": {"name": unit_name, "read_length": str(setup["length_mean"]),
                                    "read_pairs": str(given["depth"])}})
    return points, units


def units_of(opts):
    """The units to collect, in table order. A unit has its samples (a paired-end point's for pe and se); a long-read
    point replays communities of the paired-end points of the depth at its own place in the lists (modulo their
    number), of every read setup in turn, so that it can have more samples than one paired-end point
    (--long_read_samples, DEPTH:SAMPLES in --long_read_bases). The scenarios' points and units follow
    (scenario_units)."""
    pe_points = design_points(opts) if opts.samples > 0 else []
    units = []
    for read_type in opts.read_types if pe_points else ():
        if read_type in ("pe", "se"):
            units += [{"type": read_type, "point": p, "name": p["name"] + ("_se" if read_type == "se" else ""),
                       "samples": p["samples"]} for p in pe_points]
        else:
            setup = parse_long_setup(opts.pb_setup if read_type == "pb" else opts.ont_setup)
            long_depths = depths_and_samples(opts.long_read_bases, "--long_read_bases")
            if len({d for d, _ in long_depths}) != len(long_depths):  # their points would share a folder
                sys.exit(f"--long_read_bases {opts.long_read_bases!r} lists a depth more than once")
            n_depths = max(p["depth_index"] for p in pe_points) + 1
            for i, (bases, given) in enumerate(long_depths):
                name = f"{read_type}_b{bases}"
                communities = [p for p in pe_points if p["depth_index"] == i % n_depths]
                available = sum(p["samples"] for p in communities)
                # --long_read_samples (or --samples) where the communities allow it; DEPTH:SAMPLES exactly.
                samples = given or min(opts.long_read_samples or opts.samples, available)
                if samples > available:
                    sys.exit(f"{name}: {samples} long-read samples, but the paired-end points it replays the communities "
                             f"of ({', '.join(p['name'] for p in communities)}) have {available} samples: give fewer "
                             "(DEPTH:SAMPLES in --long_read_bases) or more paired-end samples")
                units.append({"type": read_type, "name": name, "setup": setup, "bases": int(float(bases)),
                              "samples": samples, "communities": communities,
                              "point": {"name": name, "read_length": str(setup["length_mean"]), "read_pairs": bases}})
    scenario_points, more = scenario_units(opts)
    return pe_points + scenario_points, units + more


def unit_communities(unit, opts, source=None):
    """The community samples a long-read unit replays, in order: [(paired-end point, its sample name)], the first
    unit["samples"] of its community points' samples, point after point (from their manifests, in their sim
    folders or source(point))."""
    source = source or (lambda point: point_dirs(point, opts)[1])
    out = []
    for point in unit["communities"]:
        for sample in dict.fromkeys(row["sample"] for row in manifest_rows(source(point))):
            if len(out) < unit["samples"]:
                out.append((point, sample))
    return out


def point_dirs(point, opts):
    base = os.path.join(opts.out, "points", point["name"])
    return base, os.path.join(base, "sim"), os.path.join(base, "protal")


def profile_dir(unit, opts):
    """Where protal writes a unit's SAMs, profiles and dumps."""
    base = point_dirs(unit["point"], opts)[0]
    return os.path.join(base, "protal_se" if unit["type"] == "se" and not drawn(unit) else "protal")


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
    key = {"command": command[:at] + command[at + 2:], "genome_table": content_hash(point_table(point, opts)),
           "simulator": identity(opts.simulator)}
    if point.get("host_pairs"):  # the host's reads, added to the community's (host_pe_jobs)
        key["host"] = {"genome": point["host_key"], "pairs": point["host_pairs"], "chunk": scenarios.HOST_PAIRS_CHUNK,
                       "art": identity("art_illumina")}
    return key


def scenario_table_path(name, opts):
    """Where a scenario's genome table is (scenarios.scenario_table)."""
    return os.path.join(opts.out, "scenarios", name, "genomes.tsv")


def point_table(point, opts):
    """The genome table a paired-end point draws its communities from: a scenario's own, or --genome_table."""
    return scenario_table_path(point["scenario"], opts) if point.get("scenario") else opts.genome_table


def scenario_command(point, opts, threads):
    """The simulator's command for a scenario's community point (scenario_units): its own genome table, species,
    abundances, strains, congeners and seed; the community's part of the read pairs, at the quality shifts
    prepare_scenarios found (point["art_shift"]); no reads (--test) for a point only of communities."""
    d = point["definition"]
    _, sim, profiles = point_dirs(point, opts)
    table = point_table(point, opts)
    command = [opts.simulator, "--genome_table", table, "-o", sim, "-n", str(point["samples"]),
               "--sample_prefix", point["name"] + "_s", "--total_read_pairs", point["community_pairs"],
               "--species_per_sample", d["species"], "--read_length", point["read_length"],
               *art_profile_args(point["sequencer"], point.get("art_shift", "")), "--fragment_mean", point["fragment_mean"],
               "--fragment_stdev", point["fragment_sd"], "--seed", str(scenarios.seed_of(opts.seed, point["scenario"])),
               "-t", str(threads), "--protal_metafile", profiles, *abundance_args(d["abundance"])]
    if d["strains"]:
        command += ["--strains_per_species", d["strains"]]
    kind, congeners = congener_spec(d["congeners"])
    if kind == "groups":
        command += ["--congener_groups", congeners]
    elif congeners > 0:
        genera = large_genera(table, congeners) if os.path.isfile(table) else []
        if not genera:
            return command, f"no genus in {table} has {congeners} species (scenario {point['scenario']}'s congeners)"
        command += ["--genus", "g__" + random.Random(scenarios.seed_of(opts.seed, point["scenario"])).choice(genera) +
                    f":{congeners}", "--pick_random_demand_if_fail"]
    if not point["reads"]:
        command += ["--test"]
    return command, None


def simulation_command(point, index, opts, threads, clades):
    """The simulator's command for a design point, and None, or why there is none."""
    if point.get("scenario"):
        return scenario_command(point, opts, threads)
    base, sim, profiles = point_dirs(point, opts)
    command = [opts.simulator, "--genome_table", opts.genome_table, "-o", sim, "-n", str(point["samples"]),
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
    kind, congeners = opts.congeners
    if kind == "genus" and congeners > 0:
        genera = large_genera(opts.genome_table, congeners)
        if not genera:
            return command, f"no genus in {opts.genome_table} has {congeners} species (--congeners)"
        command += ["--genus", "g__" + random.Random(opts.seed * 1000 + index).choice(genera) + f":{congeners}"]
    if taxa or (kind == "genus" and congeners > 0):
        command += ["--pick_random_demand_if_fail"]
    if kind == "groups":
        command += ["--congener_groups", congeners]
    return command, None


def simulate(point, index, opts, threads, clades, key):
    """Simulates the samples of a design point; None, or why it failed. clades: {rank: [held-out clade, ...]},
    of which one per rank goes into every sample (--novel_clades). The key of a point with host reads is written once
    they are added (host_pe_jobs)."""
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
    if not point.get("host_pairs"):
        write_key(os.path.join(base, "simulated.json"), key)
    return None


def host_pe_jobs(point, opts, key, started):
    """The jobs that add a scenario point's host read pairs to its samples' reads, once the community's are simulated:
    point["host_pairs"] per sample in chunks of scenarios.HOST_PAIRS_CHUNK side by side (scenarios.host_pe_chunk: ART
    in amplicon mode on fragments of the host, the point's profile and quality shifts), then each sample's chunks
    appended to its read files (gzip members one after the other), and the point's key written (simulated.json), so
    that a point stopped before is simulated again."""
    base, sim, _ = point_dirs(point, opts)
    _, rows, _ = map_rows(os.path.join(sim, "protal.meta"))
    art = shutil.which("art_illumina") or "art_illumina"
    art_args = art_options(point["sequencer"]) + point.get("art_shift", "").split()
    jobs, joins = [], []
    for s, row in enumerate(rows):
        tmp = os.path.join(sim, "tmp_host", row["SAMPLEID"])
        k = -(-point["host_pairs"] // scenarios.HOST_PAIRS_CHUNK)
        tasks = [{"sample": row["SAMPLEID"], "host": opts.host_folder, "chunk": c + 1, "art": art, "art_args": art_args,
                  "pairs": point["host_pairs"] // k + (1 if c < point["host_pairs"] % k else 0),
                  "length": int(point["read_length"]), "fragment_mean": float(point["fragment_mean"]),
                  "fragment_sd": float(point["fragment_sd"]), "seed": scenarios.seed_of(opts.seed, point["name"]) + s * 1009 + c,
                  "tmp": os.path.join(tmp, f"c{c + 1}"), "r1": os.path.join(tmp, f"c{c + 1}_R1.fq.gz"),
                  "r2": os.path.join(tmp, f"c{c + 1}_R2.fq.gz")} for c in range(k)]
        names = [f"host:{row['SAMPLEID']}:{t['chunk']}" for t in tasks]
        jobs += [{"name": n, "run": lambda t=t: (scenarios.host_pe_chunk(t), []), "priority": 1.5e9}
                 for n, t in zip(names, tasks)]

        def join(row=row, tasks=tasks, tmp=tmp):
            for column, read in (("FIRST", "r1"), ("SECOND", "r2")):
                with open(row[column], "ab") as out:
                    for t in tasks:
                        with open(t[read], "rb") as fh:
                            shutil.copyfileobj(fh, out, 16 << 20)
            shutil.rmtree(tmp, ignore_errors=True)
            return None, []
        joins.append(f"hostjoin:{row['SAMPLEID']}")
        jobs.append({"name": joins[-1], "run": join, "after": names, "priority": 2e9})

    def finish():
        shutil.rmtree(os.path.join(sim, "tmp_host"), ignore_errors=True)
        write_key(os.path.join(base, "simulated.json"), key)
        print(f"{point['name']}: {point['host_pairs']} host read pairs added to each of its {len(rows)} samples, "
              f"{clock(time.time() - started)} in all", flush=True)
        return None, []
    jobs.append({"name": f"host:{point['name']}", "run": finish, "after": joins, "priority": 2e9})
    return jobs


def design(point, index, opts, clades):
    """A paired-end point's communities without its reads (simulate_metagenomes --test: the same draws, so the same
    manifest.tsv and protal.meta as the real run), in its design folder, in seconds, for the long reads that replay
    them. -> (the folder, None or why it failed)."""
    base, sim, _ = point_dirs(point, opts)
    folder = os.path.join(base, "design")
    shutil.rmtree(folder, ignore_errors=True)
    command, error = simulation_command(point, index, opts, 1, clades)
    if error:
        return folder, error
    command[command.index("-o") + 1] = folder
    os.makedirs(folder)
    log = os.path.join(base, "design.log")
    with open(log, "w") as fh:
        rc = subprocess.run(command + ([] if "--test" in command else ["--test"]), stdout=fh,
                            stderr=subprocess.STDOUT).returncode
    return folder, None if rc == 0 else f"{point['name']}: {opts.simulator} --test failed with exit code {rc}; see {log}"


def same_design(point, folder, opts):
    """None if a paired-end point's simulated communities are those of its design run (which the long reads
    replayed), else why not."""
    def communities(rows):
        return [{k: v for k, v in row.items() if k not in ("fastq_r1", "fastq_r2")} for row in rows]
    sim = point_dirs(point, opts)[1]
    if communities(manifest_rows(sim)) != communities(manifest_rows(folder)):
        return (f"{point['name']}: its samples' communities (sim/manifest.tsv) differ from its design run's "
                f"({folder}/manifest.tsv), which its long reads replay")
    return None


def pe_threads(pending, slots):
    """Threads of each paired-end point [(index, point)]: all `slots` between them in proportion to their cost
    (samples x read pairs x read length), at least one each, at most 16 per sample (simulate_metagenomes runs a
    sample's genomes on the threads beyond its samples). -> {name: threads}."""
    if not pending:
        return {}
    cost = {p["name"]: pe_point_seconds(p, 1) for _, p in pending}
    cap = {p["name"]: 16 * p["samples"] for _, p in pending}
    total = sum(cost.values()) or 1.0
    threads = {n: max(1, min(cap[n], int(slots * c / total))) for n, c in cost.items()}
    while sum(threads.values()) < slots:  # the rest to the point that takes longest per thread
        open_ = [n for n in threads if threads[n] < cap[n]]
        if not open_:
            break
        name = max(open_, key=lambda n: cost[n] / threads[n])
        threads[name] += 1
    return threads


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


# How the long reads are made, part of a long-read point's key: points simulated otherwise (by pbsim3 per
# genome, which cut every contig's last read to its quota) are simulated again.
LONG_READS = "reads drawn by the collector, one pbsim3 --strategy templ run per sample"


def hifi_model(method="hifi"):
    """hifi_reads.MODEL (FLOW_MODEL for Ultima reads), which says how hifi_reads.py makes its reads."""
    import hifi_reads
    return hifi_reads.FLOW_MODEL if method == "ultima" else hifi_reads.MODEL
PBSIM_MIN_LENGTH = 100  # pbsim3's --length-min: its shortest read (and the shortest sequence it takes)
PBSIM_MAX_LENGTH = 1000000  # its --length-max
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def read_contigs(fasta):
    """The sequences of a FASTA (plain or gzipped) as upper-case bytes, in file order."""
    with open(fasta, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return [b"".join(record.partition(b"\n")[2].split()).upper() for record in (b"\n" + data).split(b"\n>")[1:]]


def read_length(rng, mean, sd):
    """A long read's length: from a gamma distribution of the setup's mean and SD between PBSIM_MIN_LENGTH and
    PBSIM_MAX_LENGTH, as pbsim3 draws them (its read-length table is that gamma density over that range); the
    mean if the SD is 0."""
    if sd <= 0:
        return min(max(int(mean), PBSIM_MIN_LENGTH), PBSIM_MAX_LENGTH)
    shape, scale = (mean / sd) ** 2, sd * sd / mean
    while True:
        length = int(round(rng.gammavariate(shape, scale)))
        if PBSIM_MIN_LENGTH <= length <= PBSIM_MAX_LENGTH:
            return length


def long_read_templates(task):
    """The reads of a long-read sample, drawn as they arise in sequencing, until their bases reach the sample's:
    a read's genome by relative abundance times genome length (task["genomes"]: fasta, weight), its length by
    read_length, its start uniform over the genome's contigs of PBSIM_MIN_LENGTH bases or more (a read ends
    where its contig does), either strand. Written to task["templates"] as FASTA, the reads of a genome
    together (each genome is read once per round of drawing: the cuts at contigs' ends leave a few bases to
    draw again), named g<genome>x_<n>, n = 1, 2, ... (in a chunk of a sample, every task["name_step"]-th from
    task["name_offset"] + 1, so that the chunks' names do not meet). A genome with "host" (a prepared host genome's
    folder, scenarios.Host) is read by memory map instead, its reads drawn at random places. -> (the names in file
    order, None), or (None, why it failed)."""
    rng = random.Random(task["seed"])
    genomes, setup = task["genomes"], task["setup"]
    step, offset = task.get("name_step", 1), task.get("name_offset", 0)  # chunk `offset` of `step` (long_read_chunks)
    cumulative, total = [], 0.0
    for genome in genomes:
        total += genome["weight"]
        cumulative.append(total)
    if total <= 0:
        return None, "no genome with reads to simulate (relative abundances and lengths are 0)"
    names, bases = [], 0
    with open(task["templates"], "wb") as out:
        while bases < task["bases"]:
            planned, need = collections.defaultdict(list), task["bases"] - bases
            while need > 0:
                g = bisect.bisect_right(cumulative, rng.random() * total)
                length = read_length(rng, setup["length_mean"], setup["length_sd"])
                planned[g].append(length)
                need -= length
            for g in sorted(planned):
                if genomes[g].get("host"):  # a host genome (gigabases): drawn by memory map, not read whole
                    host = scenarios.Host.of(genomes[g]["host"])
                    for length in planned[g]:
                        seq = host.draw(rng, length)
                        if rng.random() < 0.5:
                            seq = seq.translate(COMPLEMENT)[::-1]
                        names.append(f"g{g}x_{len(names) * step + offset + 1}")
                        out.write(b">" + names[-1].encode() + b"\n" + seq + b"\n")
                        bases += len(seq)
                    continue
                contigs = [c for c in read_contigs(genomes[g]["fasta"]) if len(c) >= PBSIM_MIN_LENGTH]
                if not contigs:
                    return None, f"{genomes[g]['genome']}: no sequence of {PBSIM_MIN_LENGTH} bases or more in {genomes[g]['fasta']}"
                starts = list(itertools.accumulate(len(c) - PBSIM_MIN_LENGTH + 1 for c in contigs))
                for length in planned[g]:
                    at = rng.randrange(starts[-1])
                    k = bisect.bisect_right(starts, at)
                    start = at - (starts[k - 1] if k else 0)
                    seq = contigs[k][start:start + length]
                    if rng.random() < 0.5:
                        seq = seq.translate(COMPLEMENT)[::-1]
                    names.append(f"g{g}x_{len(names) * step + offset + 1}")
                    out.write(b">" + names[-1].encode() + b"\n" + seq + b"\n")
                    bases += len(seq)
    return names, None


def last_line(path, limit=300):
    """The last non-empty line of a log, at most `limit` characters (an empty string if there is none). Only
    the end of the log is read: it may be long, and still being written."""
    try:
        with open(path, "rb") as fh:
            fh.seek(max(0, fh.seek(0, os.SEEK_END) - 16384))
            tail = fh.read().decode(errors="replace")
    except OSError:
        return ""
    return next((line.strip() for line in reversed(tail.replace("\r", "\n").split("\n")) if line.strip()), "")[:limit]


def long_read_sample(task):
    """One long-read sample (task: dict): its reads drawn as templates (long_read_templates), then a read of
    each, named after its template (g<genome>x_<n>), into task["out"]: by hifi_reads.py for a hifi setup (its flow
    model for an ultima one), else by one pbsim3 run with --strategy templ, which makes one read of each template,
    with the model's errors and qualities, and names it <id prefix>_<n> after the template's place n in the file.
    -> None, or why it failed."""
    tmp = task["tmp"]
    os.makedirs(tmp, exist_ok=True)
    task = {**task, "templates": os.path.join(tmp, "templates.fa")}
    names, error = long_read_templates(task)
    if error:
        return f"{task['sample']}: {error}"
    setup = task["setup"]
    if setup["method"] in ("hifi", "ultima"):
        import hifi_reads  # numpy: only long-read collections need it
        reads = hifi_reads.simulate(task["templates"], task["out"] + ".partial", setup["q_sd"], task["seed"],
                                    setup.get("q_mean"))
        if reads != len(names):
            return f"{task['sample']}: hifi_reads.py made {reads} reads of {len(names)} templates"
        os.replace(task["out"] + ".partial", task["out"])
        shutil.rmtree(tmp, ignore_errors=True)  # its templates: ~6 GB for a 6 Gb sample, not kept for the point's others
        return None
    prefix = os.path.join(tmp, "r")
    command = [task["pbsim"], "--strategy", "templ", "--method", setup["method"], f"--{setup['method']}", task["model"],
               "--template", task["templates"], "--accuracy-mean", str(setup["accuracy"]), "--seed", str(task["seed"]),
               "--prefix", prefix, "--id-prefix", "r"]
    if setup.get("ratio"):
        command += ["--difference-ratio", setup["ratio"]]
    log_path = os.path.join(tmp, "pbsim.log")
    with open(log_path, "w") as log:
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        return f"{task['sample']}: pbsim failed ({rc}): {last_line(log_path)}; see {log_path}"
    fastqs = sorted(set(glob.glob(prefix + ".fq*") + glob.glob(prefix + ".fastq*")))
    if len(fastqs) != 1:
        return f"{task['sample']}: pbsim wrote {len(fastqs)} FASTQ files ({prefix}.fq.gz expected); see {log_path}"
    reads, lines = 0, 0
    with (gzip.open(fastqs[0], "rb") if fastqs[0].endswith(".gz") else open(fastqs[0], "rb")) as fin, \
            gzip.open(task["out"] + ".partial", "wb", compresslevel=1) as fout:
        for line in fin:
            if lines % 4 == 0:
                reads += 1
                if reads > len(names):
                    break
                line = b"@" + names[reads - 1].encode() + b"\n"
            fout.write(line)
            lines += 1
    if reads != len(names) or lines % 4:
        return f"{task['sample']}: pbsim made {reads} reads of {len(names)} templates; see {log_path}"
    os.replace(task["out"] + ".partial", task["out"])
    shutil.rmtree(tmp, ignore_errors=True)  # its templates and pbsim3's reads, not kept for the point's other samples
    return None


# Long-read samples above this many bases are simulated in chunks side by side (--long_read_chunk): a 6 Gb Nanopore
# sample is ~40 min of one pbsim3 run, which the rest of the collection would otherwise wait for.
LONG_READ_CHUNK = 250_000_000


def long_read_chunks(task, chunk):
    """A long-read sample's task as the tasks of its chunks: [task] if it has `chunk` bases or fewer (or chunk is 0),
    else k = ceil(bases / chunk) tasks of a k-th of its bases each, with seeds of their own, the read names of chunk
    c every k-th from c + 1 (long_read_templates), each into its own file in the sample's tmp folder; join_chunks
    then writes the sample's reads."""
    bases = task["bases"]
    if not chunk or bases <= chunk:
        return [task]
    k = -(-bases // chunk)
    part = bases // k
    return [{**task, "bases": part if c < k - 1 else bases - part * (k - 1), "seed": task["seed"] * 1009 + c + 1,
             "out": os.path.join(task["tmp"], f"chunk{c + 1}.fq.gz"), "tmp": os.path.join(task["tmp"], f"c{c + 1}"),
             "name_step": k, "name_offset": c} for c in range(k)]


def join_chunks(task, chunks):
    """The chunks' reads (gzip files, one after the other: a gzip file of several members) as the sample's reads.
    -> None, or why it failed."""
    if len(chunks) == 1:
        return None
    try:
        with open(task["out"] + ".partial", "wb") as out:
            for chunk in chunks:
                with open(chunk["out"], "rb") as fh:
                    shutil.copyfileobj(fh, out, 16 << 20)
        os.replace(task["out"] + ".partial", task["out"])
    except OSError as exc:
        return f"{task['sample']}: joining its chunks: {exc}"
    shutil.rmtree(task["tmp"], ignore_errors=True)
    return None


def long_read_seconds(task):
    """A rough estimate of a long-read task's time on one core, to start the longest first (measured 2026-10-03:
    ~10 s of templates, pbsim3 ~0.4 s and hifi_reads.py ~0.15 s per Mb; short Ultima reads add the drawing of each,
    a guess of 5 us)."""
    setup = task["setup"]
    seconds = 10 + task["bases"] / 1e6 * (0.15 if setup["method"] in ("hifi", "ultima") else 0.45)
    return seconds + (task["bases"] / max(1, setup["length_mean"]) * 5e-6 if setup["method"] == "ultima" else 0)


def pe_point_seconds(point, threads):
    """A rough estimate of a paired-end point's time on `threads` threads (measured 2026-10-03: ~85 ms per genome
    of a sample, ~200 genomes, and ~96 s per million 150 bp pairs)."""
    per_sample = 20 + float(point["read_pairs"]) * int(point["read_length"]) / 150 * 96e-6
    return point["samples"] * per_sample / max(1, threads)


class Scheduler:
    """Runs jobs on `slots` slots (cores), the ready job of highest priority first; a job is ready once the jobs it
    comes after are done. A job: name, run: () -> (None or why it failed, [new jobs as add's keywords]), need: slots
    (at most all), after: names, priority. A job's new jobs join the queue. Once a job has failed no more start;
    the running ones end. run() -> {name: why it failed}."""

    def __init__(self, slots):
        self.slots = max(1, slots)
        self.pending = []

    def add(self, name, run, need=1, after=(), priority=0.0):
        self.pending.append({"name": name, "run": run, "need": max(1, min(need, self.slots)), "after": set(after),
                             "priority": priority})

    def run(self):
        done, failures, running, free = set(), {}, {}, self.slots
        with concurrent.futures.ThreadPoolExecutor(self.slots) as executor:
            while self.pending or running:
                if not failures:
                    ready = sorted((j for j in self.pending if j["after"] <= done), key=lambda j: -j["priority"])
                    for job in ready:  # the first that fits; smaller ones fill what is left
                        if job["need"] <= free or not running:
                            self.pending.remove(job)
                            running[executor.submit(job["run"])] = job
                            free -= job["need"]
                if not running:
                    break  # nothing runs and nothing can start: a job failed, or one waits for a job that never came
                finished, _ = concurrent.futures.wait(running, return_when=concurrent.futures.FIRST_COMPLETED)
                for future in finished:
                    job = running.pop(future)
                    free += job["need"]
                    try:
                        error, new = future.result()
                    except Exception as exc:  # a genome that cannot be read, say: the job fails, and says why
                        error, new = f"{type(exc).__name__}: {exc}", []
                    if error:
                        failures[job["name"]] = error
                        continue
                    done.add(job["name"])
                    for j in new:
                        self.add(**j)
        if not failures:
            for job in self.pending:
                failures[job["name"]] = f"waits for {', '.join(sorted(job['after'] - done))}, which never ran"
        return failures


def long_unit_jobs(index, unit, opts, keys, started, counter, source=None):
    """The jobs of a long-read unit's samples, made once its community points' manifests are there: each sample
    replays the community of a sample of the unit's paired-end points (unit_communities), each genome weighted by
    relative abundance times length (long_read_sample), in chunks (long_read_chunks); once all its samples are
    done, its sim/samples.tsv (sample, reads, truth, community sample) is written, and keys[name] to its
    simulated.json, so that an interrupted point is simulated again. source(point): the folder whose manifest.tsv
    and protal.meta describe the point's communities (default: its sim folder; a design run's before its reads
    are there). The truth files named are the sim folder's either way. A scenario's unit with a host share draws
    that share of its bases from the host genome (opts.host_folder), its weight that share of all."""
    source = source or (lambda point: point_dirs(point, opts)[1])
    index = unit.get("seed_index", index)
    host = unit.get("host_share", 0)
    sim = point_dirs(unit["point"], opts)[1]
    truth, genomes_of = {}, collections.defaultdict(list)
    for point in unit["communities"]:
        folder, sim_folder = source(point), point_dirs(point, opts)[1]
        _, pe_rows, _ = map_rows(os.path.join(folder, "protal.meta"))
        for row in pe_rows:
            truth[row["SAMPLEID"]] = row["PROFILE_TRUTH"] if folder == sim_folder else \
                os.path.join(sim_folder, os.path.relpath(row["PROFILE_TRUTH"], folder))
        for row in manifest_rows(folder):
            genomes_of[row["sample"]].append(row)
    model = pbsim_model(opts, unit["setup"]["model"]) if unit["setup"]["method"] in PBSIM_METHODS else None
    shutil.rmtree(os.path.join(sim, "tmp"), ignore_errors=True)
    os.makedirs(os.path.join(sim, "reads"), exist_ok=True)
    rows, jobs, chunk = [], [], getattr(opts, "long_read_chunk", LONG_READ_CHUNK)
    need = 2 if unit["setup"]["method"] in PBSIM_METHODS else 1  # pbsim3 keeps about two cores busy
    for s, (_, community) in enumerate(unit_communities(unit, opts, source)):
        sample = f"{unit['name']}_s_{s + 1}"
        out = os.path.join(sim, "reads", sample + ".fq.gz")
        rows.append((sample, out, truth[community], community))
        task = {"sample": sample, "out": out, "bases": unit["bases"], "setup": unit["setup"], "model": model,
                "pbsim": opts.pbsim, "seed": (opts.seed * 1000003 + index * 1009 + s) * 101,
                "tmp": os.path.join(sim, "tmp", sample),
                "genomes": [{"genome": g["genome"], "fasta": g["fasta_path"],
                             "weight": float(g["relative_abundance"]) * float(g["genome_length"])}
                            for g in genomes_of[community]]}
        if host > 0:
            weight = sum(g["weight"] for g in task["genomes"]) * host / (1 - host)
            task["genomes"].append({"genome": "host", "fasta": "", "host": opts.host_folder, "weight": weight})
        chunks = long_read_chunks(task, chunk)
        names = [f"long:{sample}:{c + 1}" for c in range(len(chunks))]
        jobs += [{"name": name, "run": lambda part=part: (long_read_sample(part), []), "need": need,
                  "priority": long_read_seconds(part)} for name, part in zip(names, chunks)]
        jobs.append({"name": f"join:{sample}", "run": lambda task=task, chunks=chunks: (join_chunks(task, chunks), []),
                     "after": names, "priority": 2e9})

    def finish():
        shutil.rmtree(os.path.join(sim, "tmp"), ignore_errors=True)
        with open(os.path.join(sim, "samples.tsv.partial"), "w") as fh:
            fh.write("sample\treads\ttruth\tcommunity\n" + "".join("\t".join(r) + "\n" for r in rows))
        os.replace(os.path.join(sim, "samples.tsv.partial"), os.path.join(sim, "samples.tsv"))
        if keys is not None:
            write_key(os.path.join(point_dirs(unit["point"], opts)[0], "simulated.json"), keys[unit["name"]])
        counter["long"] += 1
        print(f"{unit['name']} simulated ({len(rows)} samples): {counter['long']} of {counter['long_total']} long-read "
              f"design points, {clock(time.time() - started)} in all", flush=True)
        return None, []

    jobs.append({"name": f"long:{unit['name']}", "run": finish, "after": [f"join:{r[0]}" for r in rows], "priority": 2e9})
    return jobs


def simulate_long(points, opts, jobs, keys=None):
    """Long reads (pb, ont) of design points [(index, unit)] whose paired-end points are simulated
    (long_unit_jobs), on `jobs` slots, the longest samples (or chunks of them) first. -> {job: why it failed}."""
    scheduler, started = Scheduler(jobs), time.time()
    counter = {"long": 0, "long_total": len(points)}
    for index, unit in points:
        scheduler.add(f"units:{unit['name']}", lambda index=index, unit=unit: (
            None, long_unit_jobs(index, unit, opts, keys, started, counter)), priority=3e9)
    return scheduler.run()


def simulated(unit, opts):
    """Whether a unit's reads are there: the paired-end map (pe, se), or the samples table of drawn reads."""
    sim = point_dirs(unit["point"], opts)[1]
    return os.path.isfile(os.path.join(sim, "samples.tsv" if drawn(unit) else "protal.meta"))


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
    if not drawn(unit):
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


def write_map(path, rows):
    """A protal map of rows (MAP_COLUMNS, absolute paths)."""
    with open(path, "w") as fh:
        fh.write(f"#OUTPUT_DIR\t{os.path.dirname(path)}\n#" + "\t".join(MAP_COLUMNS) + "\n")
        fh.writelines("\t".join(row[c] for c in MAP_COLUMNS) + "\n" for row in rows)


def read_map(path):
    """The rows of a map write_map wrote, as dicts."""
    with open(path) as fh:
        return [dict(zip(MAP_COLUMNS, line.rstrip("\n").split("\t"))) for line in fh
                if line.strip() and not line.startswith("#")]


def profile_map(units, opts, path):
    """The rows of all units' samples, written as the map `path` (their output folders made). -> the rows."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = []
    for unit in units:
        unit_rows, dirs = unit_map_rows(unit, opts)
        for d in dirs:
            os.makedirs(d, exist_ok=True)
        rows += unit_rows
    write_map(path, rows)
    return rows


def profile(units, opts, extra=()):
    """Profiles the samples of all units, and the `extra` map rows (another collection's, --also_profile), in one
    protal run: the database is loaded once, and every sample is profiled as its READ_TYPE says."""
    folder = os.path.join(opts.out, "profile_all")
    combined = os.path.join(folder, "samples.map")
    rows = profile_map(units, opts, combined)
    if extra:
        rows += list(extra)
        write_map(combined, rows)
    kinds = collections.Counter(row["READ_TYPE"] for row in rows)
    print(f"profiling {len(rows)} samples ({', '.join(f'{n} {t}' for t, n in kinds.items())}) of {len(units)} design "
          "points" + (f" and {len(extra)} samples of another collection" if extra else "") + " in one protal run",
          flush=True)
    run_protal([opts.protal, "--db", opts.db, "--map", combined, "-t", str(opts.threads), "--no_strains", "--no_qcmsa"],
               os.path.join(folder, "protal.log"), len(rows))


# ---- the tables -------------------------------------------------------------------------------------------

def community_of(unit, sample, opts):
    """The paired-end point and sample whose community a unit's sample holds."""
    if unit["type"] == "pe":
        return unit["point"], sample
    if unit["type"] == "se" and not drawn(unit):
        return unit["point"], sample[:-len("_se")]
    with open(os.path.join(point_dirs(unit["point"], opts)[1], "samples.tsv")) as fh:
        next(fh)
        for name, _, _, community in (line.rstrip("\n").split("\t") for line in fh if line.strip()):
            if name == sample:
                # The simulator names a point's samples <point>_s_<n>.
                point = next((p for p in unit["communities"] if community.rsplit("_s_", 1)[0] == p["name"]), None)
                if point is None:
                    sys.exit(f"{sample}: its community {community} is of none of {unit['name']}'s paired-end points")
                return point, community
    sys.exit(f"{sample}: not in the samples of {unit['name']}")


def write_table(read_type, units, opts, context):
    """Joins the dumps of a read type's units into its table."""
    domains, novel, reps, db_lineages, sim_lineages = context
    header, rows, totals = None, 0, collections.Counter()
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
                relations = SampleRelations(sample_lineages, novel)
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
                        insilico = ""
                        if is_present and reps.get(taxon) and taxon in in_sample:
                            rep = "1" if all(g == reps[taxon] for g in in_sample[taxon]) else "0"
                        if is_present and taxon in in_sample:
                            insilico = str(int(any(g.startswith(INSILICO_PREFIX) for g in in_sample[taxon])))
                        relative, level, neighbour = "", "", ""
                        if is_present:
                            relative = "species"
                            neighbour = relations.neighbour(db_lineages.get(taxon) or sim_lineages.get(taxon, {}), taxon)
                        elif taxon in db_lineages:
                            relative, level = relations.relation(db_lineages[taxon])
                        writer.writerow([unit["name"], sample, point["read_length"], point["read_pairs"],
                                         domains.get(taxon, "unknown"), len(novel_here), int(congener), rep,
                                         novel_levels, level, relative, neighbour, read_type, insilico,
                                         unit.get("scenario", "")] + row)
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
            totals.update(present=present, absent=absent, samples=len(dumps_of(unit, opts)))
    os.replace(table + ".partial", table)
    print(f"{rows} taxa in {table}: {totals['present']} present, {totals['absent']} absent, from {totals['samples']} "
          f"sample{'s' if totals['samples'] != 1 else ''}", flush=True)


def prepare_scenarios(points, units, opts, novel):
    """What the scenarios' simulations need before their keys can be made: each scenario's genome table
    (scenarios.scenario_table, OUT/scenarios/<name>/genomes.tsv), the ART quality shifts of its paired-end reads
    (scenarios.art_shifts, point["art_shift"]) and, for a host share, the host genome as plain sequence
    (opts.host_folder, OUT/host; point["host_key"]). Stops with why a scenario cannot be simulated."""
    opts.host_folder = None
    scenario_points = [p for p in points if p.get("scenario")]
    if not scenario_points:
        return
    pool = scenarios.species_of_table(opts.genome_table)
    known, lacking = sum(s not in novel for s in pool), sum(s in novel for s in pool)
    for point in scenario_points:
        name, d = point["scenario"], point["definition"]
        if d["novel_share"] > 0 and not novel:
            sys.exit(f"scenario {name}: {d['novel_share']:.0%} of its species lack from the database, which needs the "
                     "species the database lacks (--novel_species)")
        try:
            # A scenario larger than the genome table holds at its share is scaled down to it (scenarios.fit_species);
            # the point's definition is that of its units too.
            d, scaled = scenarios.fit_species(d, known, lacking)
            point["definition"] = d
            note = scenarios.scenario_table(opts.genome_table, novel, name, d, opts.seed, scenario_table_path(name, opts))
        except scenarios.ScenarioError as e:
            sys.exit(f"scenario {name}: {e}")
        if scaled:
            print(f"scenario {name}: {scaled}", flush=True)
        if point["reads"] and point.get("quality") is not None:
            if not shutil.which("art_illumina"):
                sys.exit(f"scenario {name}: its Illumina reads' quality is set with art_illumina, which is not on PATH")
            try:
                qs1, qs2, means = scenarios.art_shifts("art_illumina", art_options(point["sequencer"]), point["read_length"],
                                                       point["fragment_mean"], point["fragment_sd"], point["quality"],
                                                       os.path.join(opts.out, "scenarios", "art_quality.json"))
            except scenarios.ScenarioError as e:
                sys.exit(f"scenario {name}: {e}")
            point["art_shift"] = f"-qs {qs1} -qs2 {qs2}"
            note += (f"; Illumina reads at Q{point['quality']:g} (ART {point['sequencer']}: Q{means[0]:.1f} and "
                     f"Q{means[1]:.1f}, shifted by {qs1:+d} and {qs2:+d})")
        if d["host_share"] > 0:
            note += f"; {d['host_share']:.0%} of the reads from the host"
        print(f"scenario {name} ({point['samples']} samples): {note}", flush=True)
    if any(p["definition"]["host_share"] > 0 for p in scenario_points):
        named = ", ".join(sorted({p["scenario"] for p in scenario_points if p["definition"]["host_share"] > 0}))
        if not opts.host_genome:
            sys.exit(f"scenario {named}: its host's reads need the host genome (--host_genome)")
        try:
            opts.host_folder = scenarios.prepare_host(opts.host_genome, os.path.join(opts.out, "host"))
        except (OSError, scenarios.ScenarioError) as e:
            sys.exit(f"scenario {named}: host genome {opts.host_genome}: {e}")
        host = scenarios.host_identity(opts.host_folder)
        print(f"host genome {opts.host_genome}: {host['bases']} bases in {opts.host_folder}", flush=True)
        for point in scenario_points:
            if point.get("host_pairs"):
                point["host_key"] = host


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
    slots = max(1, opts.jobs or opts.threads)
    prepare_scenarios(pe_points, units, opts, novel)

    # Every read type needs paired-end points: se reads them, pb and ont replay their communities. All simulations
    # run in one queue on `slots` cores (Scheduler): the paired-end points with all the threads between them, and the
    # long-read samples (or chunks of them) as cores come free, the longest first. A long-read point needs only the
    # communities of its paired-end points: a design run (simulate_metagenomes --test, the same communities without
    # reads) gives them in seconds, so long reads need not wait for the paired-end reads.
    needed = {u["point"]["name"] for u in units if not drawn(u)} | \
             {p["name"] for u in units if drawn(u) for p in u["communities"]}
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
    long_units = [u for u in units if drawn(u)]
    long_pending = []
    for i, unit in enumerate(long_units):
        base = point_dirs(unit["point"], opts)[0]
        pbsim = unit["setup"]["method"] in PBSIM_METHODS
        keys[unit["name"]] = {"communities": [keys[p["name"]] for p in unit["communities"]], "samples": unit["samples"],
                              "setup": unit["setup"], "bases": unit["bases"],
                              "index": unit.get("seed_index", i), "seed": opts.seed,
                              "pbsim": identity(opts.pbsim) if pbsim else None,
                              "model": identity(pbsim_model(opts, unit["setup"]["model"])) if pbsim else
                              hifi_model(unit["setup"]["method"]), "reads": LONG_READS}
        if opts.long_read_chunk and unit["bases"] > opts.long_read_chunk:  # chunks make other reads
            keys[unit["name"]]["chunk"] = opts.long_read_chunk
        if unit.get("host_share"):
            keys[unit["name"]]["host"] = {"genome": scenarios.host_identity(opts.host_folder), "share": unit["host_share"]}
        if simulated(unit, opts) and not same_key(os.path.join(base, "simulated.json"), keys[unit["name"]]):
            print(f"{unit['name']} was simulated from other inputs (or by an older collector): simulating it again", flush=True)
            shutil.rmtree(base)
        if not simulated(unit, opts):
            long_pending.append((i, unit))
    if pending or long_pending:
        scheduler, started = Scheduler(slots), time.time()
        counter = {"pe": 0, "long": 0, "long_total": len(long_pending)}
        # The paired-end points share the cores with the long reads by their estimated work, so that both end about
        # together (the long reads take the cores the paired-end points leave, and those they free).
        pe_work = sum(pe_point_seconds(p, 1) for _, p in pending)
        long_work = sum(u["samples"] * long_read_seconds(u) * (2 if u["setup"]["method"] in PBSIM_METHODS else 1)
                        for _, u in long_pending)
        threads_of = pe_threads(pending, max(1, round(slots * pe_work / ((pe_work + long_work) or 1))))
        if pending:
            more = [f"{p['name']} {threads_of[p['name']]}" for _, p in pending if threads_of[p["name"]] > 1]
            print(f"simulating {len(pending)} paired-end design points"
                  + (f" (threads: {', '.join(more)}, the others 1)" if more else ""), flush=True)
        if long_pending:
            print(f"simulating {len(long_pending)} long-read design points "
                  f"({sum(u['samples'] for _, u in long_pending)} samples) as cores come free, the longest samples first"
                  + (f", those above {opts.long_read_chunk} bases in chunks" if opts.long_read_chunk else ""), flush=True)
        designed = {}  # paired-end point name -> its design folder, for the long reads until its reads are there
        for i, p in pending:
            def simulate_point(i=i, p=p):
                began = time.time()
                failure = simulate(p, i, opts, threads_of[p["name"]], clades, keys[p["name"]])
                if not failure and p["name"] in designed:
                    failure = same_design(p, designed[p["name"]], opts)
                counter["pe"] += 1
                print(f"{p['name']} {'failed' if failure else 'simulated'} ({p['samples']} samples) in "
                      f"{clock(time.time() - began)}: {counter['pe']} of {len(pending)} paired-end design points, "
                      f"{clock(time.time() - started)} in all", flush=True)
                if failure or not p.get("host_pairs"):
                    return failure, []
                return None, host_pe_jobs(p, opts, keys[p["name"]], started)  # then its host's reads
            # Before any long-read sample (all start at once): the deepest could otherwise wait for long reads.
            scheduler.add(f"pe:{p['name']}", simulate_point, need=threads_of[p["name"]],
                          priority=1e9 + pe_point_seconds(p, threads_of[p["name"]]))
        pending_names = {p["name"] for _, p in pending}
        for i, p in pending:
            if any(p in u["communities"] for _, u in long_pending):
                def design_point(i=i, p=p):
                    folder, error = design(p, i, opts, clades)
                    designed[p["name"]] = folder
                    return error, []
                scheduler.add(f"design:{p['name']}", design_point, priority=4e9)
        for i, unit in long_pending:
            after = [f"design:{p['name']}" for p in unit["communities"] if p["name"] in pending_names]
            source = lambda point: designed.get(point["name"]) or point_dirs(point, opts)[1]
            scheduler.add(f"units:{unit['name']}", lambda i=i, unit=unit, source=source: (
                None, long_unit_jobs(i, unit, opts, keys, started, counter, source)), after=after, priority=3e9)
        failures = scheduler.run()
        if failures:
            sys.exit("\n".join(f"{name}: {why}" for name, why in list(failures.items())[:10]))
        for folder in designed.values():
            shutil.rmtree(folder, ignore_errors=True)
    if opts.simulate_only:
        print("simulated every design point; profiling left to a run without --simulate_only", flush=True)
        return
    # Profiling: every unit not yet profiled against this database with this protal, in one protal run.
    db, protal = db_identity(opts.db), identity(opts.protal)
    unprofiled = []
    for unit in units:
        folder = profile_dir(unit, opts)
        key = {"simulated": keys[unit["name"] if drawn(unit) else unit["point"]["name"]],
               "db": db, "protal": protal, "read_type": unit["type"]}
        if same_key(os.path.join(folder, "profiled.json"), key) and len(dumps_of(unit, opts)) >= unit["samples"]:
            continue
        if not same_key(os.path.join(folder, "profiling.json"), key):  # else a stopped run's: go on with it
            if dumps_of(unit, opts):
                print(f"{unit['name']} was profiled against another database or with another protal (or by an older "
                      "collector): profiling it again", flush=True)
            shutil.rmtree(folder, ignore_errors=True)
            os.makedirs(folder)
            write_key(os.path.join(folder, "profiling.json"), key)
        unprofiled.append(unit)
    if opts.prepare_profiling:
        combined = os.path.join(opts.out, "profile_all", "samples.map")
        for stale in (combined, combined + ".units"):
            if os.path.exists(stale):
                os.remove(stale)
        if unprofiled:
            rows = profile_map(unprofiled, opts, combined)
            with open(combined + ".units", "w") as fh:
                fh.writelines(profile_dir(unit, opts) + "\n" for unit in unprofiled)
            print(f"{len(rows)} samples of {len(unprofiled)} design points to profile, mapped in {combined} for another "
                  "collection's protal run (--also_profile)", flush=True)
        else:
            print("every design point is profiled", flush=True)
        return
    others = [(path, read_map(path)) for path in opts.also_profile]
    if unprofiled or any(rows for _, rows in others):
        profile(unprofiled, opts, [row for _, rows in others for row in rows])
        for unit in unprofiled:
            folder = profile_dir(unit, opts)
            os.replace(os.path.join(folder, "profiling.json"), os.path.join(folder, "profiled.json"))
        for path, _ in others:  # the other collection's run finds them profiled
            with open(path + ".units") as fh:
                for folder in (line.rstrip("\n") for line in fh if line.strip()):
                    os.replace(os.path.join(folder, "profiling.json"), os.path.join(folder, "profiled.json"))
            os.remove(path + ".units")
    for unit in units:
        if len(dumps_of(unit, opts)) != unit["samples"]:
            sys.exit(f"{unit['name']}: expected {unit['samples']} training dumps in {profile_dir(unit, opts)}, "
                     f"found {len(dumps_of(unit, opts))}")

    context = (domains, novel, reps, db_lineages, sim_lineages)
    for read_type in opts.read_types:
        write_table(read_type, [u for u in units if u["type"] == read_type], opts, context)


if __name__ == "__main__":
    main()
