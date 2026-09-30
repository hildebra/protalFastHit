#!/usr/bin/env python3
"""Simulate metagenomes and profile them with protal, to train its presence model.

For every design point, a read setup (length, ART profile, fragment size) and a sequencing depth,
simulate_metagenomes draws random communities from a genome table and protal profiles them
against a database, knowing the true species. The training dumps of all samples
(<profile>.truth_annotated: every taxon protal saw, its features and whether it was present) are
joined into one table, with meta_* columns saying where each row comes from (meta_domain from
the genome table's lineages; with --novel_species and --taxonomy, which rows share a genus with a
species the database lacks, and which present species were simulated from another genome than
the database's reference). Train on it with

    python3 scripts/random_forest_cmdline.py --truth-file OUT/training_data.tsv --output-prefix OUT/model

A genome table with species the database lacks gives the negatives that matter most: a relative
the database has picks up their reads. When whole genera, families, orders, classes or phyla are
missing (--novel_species with ranks, --novel_clades), the relatives are distant: meta_novel_level
marks the absent taxa closest to such species, meta_neighbour_rank how close a present taxon's
nearest other species in the sample is. Archaea (--archaea) have fewer marker genes than bacteria
and need to be in the training data, too. Design points are simulated in parallel (--jobs; ART
simulates one genome at a time), then all their samples are profiled in one protal run, which
loads the database once. Points already simulated or profiled are skipped, so a run can be resumed.

usage: collect_training_data.py --db DB --genome_table genomes.tsv -o OUT [options]
"""

import argparse
import collections
import concurrent.futures
import csv
import glob
import os
import random
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lineages  # noqa: E402


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
    p.add_argument("--read_setups", default="100:HS20:300:40,150:HS25:350:50,250:MSv3:550:50",
                   help="comma-separated LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, one design point each")
    p.add_argument("--species_per_sample", default="5-30", help="species per sample, N or MIN-MAX")
    p.add_argument("--archaea", type=int, default=0, help="archaeal species per sample (default: 0)")
    p.add_argument("--congeners", type=int, default=0,
                   help="species of one genus in every sample of a design point, the genus drawn per point among "
                        "those with that many species (default: 0). Species are otherwise drawn uniformly, so among "
                        "many genera relatives hardly ever share a sample, while in real samples they often do")
    p.add_argument("-t", "--threads", type=int, default=4, help="threads of the protal run (default 4)")
    p.add_argument("--jobs", type=int, default=0,
                   help="design points simulated at a time (default: --threads; ART is single-threaded)")
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
    return p.parse_args(argv)


META_COLUMNS = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_domain",
                "meta_novel_species", "meta_novel_congener", "meta_rep_genome", "meta_novel_levels",
                "meta_novel_level", "meta_relative_rank", "meta_neighbour_rank"]
# meta_novel_levels: the sample's species the database lacks, by the rank they were held out at
# ("species:2,family:1"). meta_relative_rank: the deepest rank the taxon shares with a species simulated in
# the sample ("species" for the simulated species themselves, "none" for no shared domain). meta_novel_level:
# for an absent taxon whose closest simulated species (at least as close as any present one) is one the
# database lacks, the rank that species was held out at: its reads are the likely source of the taxon's.
# meta_neighbour_rank: for a present taxon, the deepest rank it shares with another species simulated in the
# sample (a congener's reads fit it nearly as well, so it may be missed).


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


def simulated_genomes(point_dir):
    """sample -> {species name: [genome, ...]} from the simulator's manifest."""
    genomes = {}
    with open(os.path.join(point_dir, "sim", "manifest.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
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


def design_points(opts):
    points = []
    for setup in opts.read_setups.split(","):
        length, profile, fragment_mean, fragment_sd = setup.split(":")
        for pairs in opts.read_pairs.split(","):
            points.append({"name": f"rl{length}_p{pairs}", "read_length": length, "sequencer": profile,
                           "fragment_mean": fragment_mean, "fragment_sd": fragment_sd, "read_pairs": pairs})
    return points


def point_dirs(point, opts):
    base = os.path.join(opts.out, "points", point["name"])
    return base, os.path.join(base, "sim"), os.path.join(base, "protal")


def dumps_of(point, opts):
    return sorted(glob.glob(os.path.join(point_dirs(point, opts)[2], "**", "*.truth_annotated"), recursive=True))


def simulate(point, index, opts, threads, clades):
    """Simulates the samples of a design point; None, or why it failed. clades: {rank: [held-out clade, ...]},
    of which one per rank goes into every sample (--novel_clades)."""
    base, sim, profiles = point_dirs(point, opts)
    command = [opts.simulator, "--genome_table", opts.genome_table, "-o", sim, "-n", str(opts.samples),
               "--sample_prefix", point["name"] + "_s", "--total_read_pairs", point["read_pairs"],
               "--species_per_sample", opts.species_per_sample, "--read_length", point["read_length"],
               "--sequencer", point["sequencer"], "--fragment_mean", point["fragment_mean"],
               "--fragment_stdev", point["fragment_sd"], "--seed", str(opts.seed + index),
               "-t", str(threads), "--protal_metafile", profiles]
    # One --taxon for all demands: the simulator reads only the last.
    taxa = [f"d__Archaea:{opts.archaea}"] if opts.archaea > 0 else []
    if opts.novel_clades > 0:
        taxa += [f"{clades[rank][index % len(clades[rank])]}:{opts.novel_clades}" for rank in sorted(clades)]
    if taxa:
        command += ["--taxon", ",".join(taxa)]
    if opts.congeners > 0:
        genera = large_genera(opts.genome_table, opts.congeners)
        if not genera:
            return f"no genus in {opts.genome_table} has {opts.congeners} species (--congeners)"
        command += ["--genus", "g__" + random.Random(opts.seed * 1000 + index).choice(genera) + f":{opts.congeners}"]
    if taxa or opts.congeners > 0:
        command += ["--pick_random_demand_if_fail"]
    os.makedirs(sim, exist_ok=True)
    log = os.path.join(base, "simulate.log")
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        shutil.rmtree(sim, ignore_errors=True)  # no protal.meta: simulated again on a rerun
        return f"{point['name']}: {opts.simulator} failed with exit code {rc}; see {log}"
    return None


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


def profile(points, opts):
    """Profiles the samples of all points in one protal run: the database is loaded once."""
    folder = os.path.join(opts.out, "profile_all")
    os.makedirs(folder, exist_ok=True)
    header, rows = None, []
    for point in points:
        point_header, point_rows, dirs = map_rows(os.path.join(point_dirs(point, opts)[1], "protal.meta"))
        if header is None:
            header = point_header
        elif point_header != header:
            sys.exit(f"{point['name']}: its protal.meta has other columns than the others")
        for d in dirs:
            os.makedirs(d, exist_ok=True)
        rows += point_rows
    combined = os.path.join(folder, "samples.map")
    with open(combined, "w") as fh:
        fh.write(f"#OUTPUT_DIR\t{folder}\n#" + "\t".join(header) + "\n")
        fh.writelines("\t".join(row[c] for c in header) + "\n" for row in rows)
    print(f"profiling {len(rows)} samples of {len(points)} design points in one protal run", flush=True)
    run([opts.protal, "--db", opts.db, "--map", combined, "-t", str(opts.threads), "--no_strains", "--no_qcmsa"],
        os.path.join(folder, "protal.log"))


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
    points = design_points(opts)

    # Simulation: ART simulates one genome at a time, so design points run in parallel.
    pending = [(i, p) for i, p in enumerate(points) if not os.path.isfile(os.path.join(point_dirs(p, opts)[1], "protal.meta"))]
    if pending:
        jobs = max(1, min(opts.jobs or opts.threads, len(pending)))
        threads = max(1, opts.threads // jobs)
        print(f"simulating {len(pending)} design points, {jobs} at a time", flush=True)
        with concurrent.futures.ThreadPoolExecutor(jobs) as executor:
            failures = [f for f in executor.map(lambda ip: simulate(ip[1], ip[0], opts, threads, clades), pending) if f]
        if failures:
            sys.exit("\n".join(failures))
    # Profiling: every point not yet profiled, in one protal run.
    unprofiled = [p for p in points if len(dumps_of(p, opts)) < opts.samples]
    if unprofiled:
        profile(unprofiled, opts)
    for point in points:
        if len(dumps_of(point, opts)) != opts.samples:
            sys.exit(f"{point['name']}: expected {opts.samples} training dumps in {point_dirs(point, opts)[2]}, "
                     f"found {len(dumps_of(point, opts))}")

    header, rows = None, 0
    table = os.path.join(opts.out, "training_data.tsv")
    with open(table + ".partial", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for point in points:
            present = absent = strains = congeners = 0
            attributed = collections.Counter()
            dumps = dumps_of(point, opts)
            genomes = simulated_genomes(os.path.join(opts.out, "points", point["name"]))
            for dump in dumps:
                sample = os.path.basename(dump).split(".profile")[0]
                in_sample = genomes.get(sample, {})
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
                        writer.writerow([point["name"], sample, point["read_length"], point["read_pairs"],
                                         domains.get(taxon, "unknown"), len(novel_here), int(congener), rep,
                                         novel_levels, level, relative, neighbour] + row)
                        rows += 1
                        present += is_present
                        absent += not is_present
                        strains += rep == "0"
                        congeners += congener and not is_present
                        attributed[level] += bool(level)
            by_level = ", ".join(f"{n} to {rank}" for rank, n in sorted(attributed.items()) if rank)
            print(f"{point['name']}: {present} present ({strains} from other genomes than the representative) and "
                  f"{absent} absent taxa ({congeners} congeners of species the database lacks; closest to a species "
                  f"it lacks, by the rank held out: {by_level or 'none'}) in {opts.samples} samples", flush=True)
    os.replace(table + ".partial", table)
    print(f"{rows} taxa in {table}")


if __name__ == "__main__":
    main()
