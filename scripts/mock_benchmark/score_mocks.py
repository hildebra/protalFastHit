#!/usr/bin/env python3
"""Run protal (or sylph) on the mock-community samples that download_mocks.py fetched, and score the profiles against
the truth tables (truth/, prepare_tables.py) and against each other (samples.tsv's groups).

  run    runs a profiler on each set (or each chunk of a set: --chunks N, the chunk from --chunk or SLURM's
         $SLURM_ARRAY_TASK_ID) into a run folder of its own, OUT/<run>/<set>[/chunk_<i>]/, with the profiles in
         profiles/<sample>.profile. A sample whose profile the run folder has is not run again (--rerun: all).
         OUT/<run>/versions.tsv notes per sample profiled the tool, its version, the database, when and where.
           --tool protal  (default) a protal map per set (map.tsv) and one protal run on it (protal.log, misc/), with
                          --no_strains unless --strains; --protal-args adds options (quoted). The run folder is
                          <version>_<commit> of the binary's --version (e.g. 0.8.3_1dcc325, "-modified" when built
                          with uncommitted changes), unless --run-name names it. The SAMs go to alignments/, or with
                          --scratch DIR (protal's --scratch_dir) into protal's folder in DIR, each removed once
                          profiled, so only the profiles and logs stay (with --strains they stay in alignments/, for
                          strain MSAs). Older versions, by their --version: before 0.8.3 (no --scratch_dir, one file
                          per map cell) the SAMs go to a folder of their own in --scratch, removed after the run, and
                          a sample of several files is joined first; before 0.7.0 (0.6.0a says v0.6.0: run folder
                          0.6.0_nocommit) the samples that are not paired-end are left out, and the SAMs are plain
           --tool sylph   sylph profile of the set's paired samples and of its single-read samples (sylph/*.tsv,
                          sylph.log; a sample of several files is joined first, in --scratch if given), against
                          --sylph-db. Each genome is given its species by accession in the GTDB taxonomy of --gtdb
                          (a folder of a release's taxonomy or metadata files, not the truth; the truth's release
                          first, so a database of another release is scored in the truth's names; then, comma-separated,
                          the database's release for the genomes the truth's lacks; a genome neither has is kept
                          under its accession), and its Taxonomic_abundance (a share of cells) is summed per species
                          into a profile in protal's format, without "?".
                          sylph runs once per sample and database (since 2026-10-14): the run folder is the one OUT
                          has for the database, whatever sylph is installed now (or none, if every sample has its
                          profile), so that its profiles are the comparison of every later protal run; only samples
                          without a profile there are run, with the version installed (versions.tsv says which).
                          Without one, and with --rerun (all samples again), it is sylph_<version>_<database name>
                          of the sylph installed; --run-name names another
  score  scores the profiles of each run folder (--runs, default all of them; the sylph run folders are scored with
         the runs named, for comparison, unless --no-sylph) and writes OUT/<run>/scores/:
           samples.tsv     per sample: the species called, the profile's unknown share ("?") and the one the truth
                           expects, and by the truth's kind
                             abundance  true and false positives and negatives, precision, recall, F1, recall of the
                                        species at >= 0.1 % and >= 1 % of the placed cells, the share of the profile on
                                        false positives, and the L1 distance and Bray-Curtis dissimilarity of the
                                        profile and the truth (both over the placed species, adding up to 1)
                             closed     the same, without the distances (members' shares unknown)
                             dropout    the strains left out that were called, the ones kept that were missed
                             negative   the species called (all false)
           species.tsv     per sample and species of the truth or the profile: truth share, profile share, status
           replicates.tsv  per pair of samples of a group: Bray-Curtis of their profiles, Jaccard of their calls
           groups.tsv      per group: samples, mean and largest Bray-Curtis, mean and smallest Jaccard
           summary.tsv     per set, read type and truth kind: medians of the above
         and OUT/scores/summary_by_run.tsv and samples_by_run.tsv, all runs' rows with columns run, tool and version
         (from versions.tsv, else from the run folder's default name; "-" for neither), each set's rows of all runs
         together (printed too).
         Shares compare cells: protal's abundances are shares of genomes, the truth's genome copies (prepare_tables.py).

Usage on a SLURM cluster (protal at r226 needs ~40 GB of memory):
  sbatch -c 32 --mem 64G --wrap "python3 scripts/mock_benchmark/score_mocks.py run --reads /scratch/mocks \\
      --out /path/results --protal /path/protal --db /path/protal_db -t 32 --scratch \\$TMPDIR \\
      --sets d6331_dna,meslier,d6331_hifi"
  sbatch -c 16 --mem 48G --array 0-9 --wrap "python3 scripts/mock_benchmark/score_mocks.py run \\
      --reads /scratch/mocks --db /path/protal_db -t 16 --sets hcom2_invitro --chunks 10"
  sbatch -c 32 --mem 64G --wrap "python3 scripts/mock_benchmark/score_mocks.py run --tool sylph \\
      --reads /scratch/mocks --sylph-db gtdb-r232-c200-dbv2.syl2db --gtdb /path/GTDB/r226,/path/GTDB/r232 -t 32 \\
      --sets d6331_dna,meslier"
  python3 scripts/mock_benchmark/score_mocks.py score --reads /scratch/mocks
Standard library only.
"""
import argparse
import collections
import csv
import glob
import gzip
import itertools
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
NAN = float("nan")
THRESHOLDS = (1e-3, 1e-2)
SAMPLE_COLUMNS = ["set", "sample", "read_type", "truth", "kind", "group", "reads", "bases", "cut", "called", "unknown",
                  "unknown_expected", "tp", "fp", "fn", "precision", "recall", "f1", "recall_0.1pct", "recall_1pct",
                  "fp_share", "l1", "bray_curtis", "absent_called", "present_missed", "false_species"]
SUMMARY_COLUMNS = ["set", "read_type", "kind", "samples", "called", "tp", "fp", "fn", "precision", "recall", "f1",
                   "recall_1pct", "fp_share", "l1", "unknown", "dropouts_called"]
NOT_RUNS = ("reads", "scores")  # folders of OUT that are not run folders


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path + ".part", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(v) for k, v in r.items()})
    os.replace(path + ".part", path)


def fmt(value):
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.6g}"
    return value


def chosen(rows, sets):
    if not sets:
        return rows
    known = set(r["set"] for r in rows)
    unknown = [s for s in sets if s not in known]
    if unknown:
        sys.exit(f"unknown set(s): {', '.join(unknown)}")
    return [r for r in rows if r["set"] in sets]


# ---------------------------------------------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------------------------------------------

def profiles_of(folder):
    """{sample: path of its .profile} anywhere under a run folder (profiles/, or the set's folder itself, where protal
    writes a map's <PREFIX>.profile without a PROFILE column)."""
    found = {}
    for path in glob.glob(os.path.join(folder, "**", "*.profile"), recursive=True):
        found[os.path.basename(path)[:-len(".profile")]] = path
    return found


def version_lines(command):
    try:
        result = subprocess.run([command, "--version"], capture_output=True, text=True, timeout=300)
    except OSError as err:
        sys.exit(f"{command} --version failed ({err}); --protal or --sylph names the binary")
    return (result.stdout + result.stderr).strip().splitlines()


def protal_version(protal):
    """<version>_<commit> from protal's --version ("protal v0.8.3 (commit <40 hex>[, with uncommitted changes])")."""
    lines = version_lines(protal)
    found = re.fullmatch(r"\S+ v([0-9.]+)(?: \(commit ([0-9a-f]{40})(, with uncommitted changes)?\))?",
                         lines[-1] if lines else "")
    if not found:
        sys.exit(f"{protal} --version says {lines[-1] if lines else 'nothing'!r}, no version")
    version, commit, modified = found.groups()
    return f"{version}_{commit[:7] if commit else 'nocommit'}{'-modified' if modified else ''}"


def sylph_version(sylph):
    """The version sylph --version names ("sylph 0.8.1": 0.8.1)."""
    lines = version_lines(sylph)
    found = re.search(r"(\d+\.\d+[\w.-]*)", lines[-1] if lines else "")
    if not found:
        sys.exit(f"{sylph} --version says {lines[-1] if lines else 'nothing'!r}, no version")
    return found.group(1)


def database_name(database):
    """A sylph database's name: its file name without .syldb, .syl2db or .sylsp."""
    return re.sub(r"\.(syldb|syl2db|sylsp)$", "", os.path.basename(database.rstrip("/")))


# Each run folder's versions.tsv: per sample profiled, the tool and its version (protal: <version>_<commit>, as the run
# folder's default name; sylph: its version), the database, when and where. Appended by every run that writes profiles.
VERSION_COLUMNS = ["sample", "set", "tool", "version", "database", "date", "host"]


def read_versions(folder):
    """{sample: its last row of the run folder's versions.tsv}; empty without the file (a run of an older script)."""
    path = os.path.join(folder, "versions.tsv")
    # (a header written twice, by chunks of a SLURM array starting together, is no sample)
    return {r["sample"]: r for r in read_tsv(path) if r["sample"] != "sample"} if os.path.exists(path) else {}


def note_versions(folder, rows):
    path = os.path.join(folder, "versions.tsv")
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=VERSION_COLUMNS, delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)


def tool_and_version(name, versions):
    """(tool, versions joined by ",") of a run folder from its versions.tsv rows, else from its default name
    (<version>_<commit> for protal, sylph_<version>_<database> for sylph); "-" where neither says."""
    if versions:
        tools = sorted({r["tool"] for r in versions.values()})
        found = sorted({r["version"] for r in versions.values()})
        return ",".join(tools), ",".join(found)
    found = re.fullmatch(r"sylph_([^_]+)_.+", name)
    if found:
        return "sylph", found.group(1)
    if re.fullmatch(r"[0-9.]+_[0-9a-f]{7}(-modified)?|[0-9.]+_nocommit(-modified)?", name):
        return "protal", name
    return "-", "-"


def sylph_folders(out, database=None):
    """The sylph run folders of OUT (by their versions.tsv, else by the default name sylph_<version>_<database>), those of
    `database` alone if given."""
    found = []
    if not os.path.isdir(out):
        return found
    for name in sorted(os.listdir(out)):
        folder = os.path.join(out, name)
        if name in NOT_RUNS or not os.path.isdir(folder):
            continue
        versions = read_versions(folder)
        if versions:
            rows = [r for r in versions.values() if r["tool"] == "sylph"]
            if rows and (database is None or any(r["database"] == database for r in rows)):
                found.append(folder)
            continue
        named = re.fullmatch(r"sylph_[^_]+_(.+)", name)
        if named and (database is None or named.group(1) == database):
            found.append(folder)
    return found


def run_command(command, folder, log_name):
    start = time.time()
    with open(os.path.join(folder, log_name), "a") as log:
        log.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(shlex.quote(c) for c in command)}\n")
        log.flush()
        code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    return code, time.time() - start


# What older protal versions lack of the maps and options used here: before 0.7.0 (0.6.0a, whose --version says
# v0.6.0) paired-end reads only, and plain SAMs whatever their name; before 0.8.3 a map cell of several files (a
# sample's lanes) and --scratch_dir.
PAIRED_ONLY_BEFORE = (0, 7, 0)
LANES_AND_SCRATCH_FROM = (0, 8, 3)


def release_of(version):
    """(major, minor, patch) of protal's <version>_<commit> (protal_version), None for None."""
    return tuple(int(x) for x in version.split("_")[0].split(".")) if version else None


def work_folder(args, folder, name):
    """A folder for a run's intermediate files, to be removed after it: a new one in --scratch, else <folder>/<name>."""
    if args.scratch:
        os.makedirs(args.scratch, exist_ok=True)
        return tempfile.mkdtemp(prefix=f"score_mocks_{name}_", dir=args.scratch)
    path = os.path.join(folder, name)
    os.makedirs(path, exist_ok=True)
    return path


def run_protal(args, folder, todo, release):
    """One protal run of a set's map. `release` is the binary's (major, minor, patch); None (unknown: --dry-run with
    --run-name) is taken as the current protal's."""
    current = release is None or release >= LANES_AND_SCRATCH_FROM
    map_path = os.path.join(folder, "map.tsv")
    # SAM and PROFILE named, so that they go to alignments/ and profiles/ (without them protal writes
    # <PREFIX>.sam.zst and <PREFIX>.profile into #OUTPUT_DIR itself). With --scratch and without --strains no SAM
    # named: protal names it and writes it into its folder in --scratch_dir, and removes it once profiled. A protal
    # without --scratch_dir writes the SAMs into a folder of their own in --scratch (#SAM_OUTPUT_DIR), removed after the
    # run, and reads a sample of several files joined into one.
    keep_sams = not (current and args.scratch and not args.strains)
    sam_ending = ".sam" if release and release < PAIRED_ONLY_BEFORE else ".sam.zst"
    temporary = []  # folders removed after the run
    try:
        sam_folder = None
        if not current and args.scratch and not args.dry_run:
            sam_folder = work_folder(args, folder, "sams")
            temporary.append(sam_folder)
        lanes = None
        with open(map_path, "w") as fh:
            fh.write(f"#OUTPUT_DIR\t{folder}\n" + (f"#SAM_OUTPUT_DIR\t{sam_folder}\n" if sam_folder else ""))
            columns = ["SAMPLEID", "PREFIX", "FIRST", "SECOND", "READ_TYPE"] + (["SAM"] if keep_sams else [])
            fh.write("#" + "\t".join(columns + ["PROFILE"]) + "\n")
            for s in todo:
                first, second = s["first"], s["second"]
                if not current and "," in first + second:
                    if args.dry_run:  # the map only: no files joined
                        first, second = first.split(",")[0], second.split(",")[0]
                    else:
                        if lanes is None:
                            lanes = work_folder(args, folder, "joined")
                            temporary.append(lanes)
                        first = joined(first, os.path.join(lanes, s["sample"] + "_1.fq.gz"))
                        second = second if second == "-" else joined(second,
                                                                     os.path.join(lanes, s["sample"] + "_2.fq.gz"))
                row = [s["sample"], s["sample"], first, second, s["read_type"]]
                row += ([s["sample"] + sam_ending] if keep_sams else []) + [s["sample"] + ".profile"]
                fh.write("\t".join(row) + "\n")
        command = [args.protal, "--db", args.db, "--map", map_path, "-t", str(args.threads)]
        if args.scratch and current:
            command += ["--scratch_dir", args.scratch]
        if not args.strains:
            command.append("--no_strains")
        command.append("--force")  # a SAM without its profile may be a killed run's
        command += shlex.split(args.protal_args)
        print(f"  {' '.join(shlex.quote(c) for c in command)}", file=sys.stderr)
        if args.dry_run:
            return 0
        code, seconds = run_command(command, folder, "protal.log")
        print(f"  protal exit {code} after {seconds:.0f} s (log: {folder}/protal.log)", file=sys.stderr)
        return code
    finally:
        for path in temporary:
            shutil.rmtree(path, ignore_errors=True)


def joined(files, path):
    """One file of a sample's comma-separated files (gzip members joined), or the file itself if it is one."""
    files = files.split(",")
    if len(files) == 1:
        return files[0]
    with open(path + ".part", "wb") as out:
        for f in files:
            with open(f, "rb") as fh:
                shutil.copyfileobj(fh, out, 1 << 22)
    os.replace(path + ".part", path)
    return path


def load_taxonomy(folders):
    """{accession without version, GCA_ for GCF_: GTDB lineage} from releases' taxonomy (or metadata) files; the first
    folder that has a genome names it (the truth's release first, a newer one for the genomes it added)."""
    lineages = {}
    for folder in folders:
        files = sorted(glob.glob(os.path.join(folder, "*_taxonomy_r*.tsv*")))
        metadata = not files
        if metadata:
            files = sorted(glob.glob(os.path.join(folder, "*_metadata_r*.tsv*")))
        if not files:
            sys.exit(f"--gtdb {folder}: no bac120/ar53_taxonomy_rNNN.tsv(.gz) or _metadata_rNNN.tsv(.gz)")
        for path in files:
            with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
                column = 1
                if metadata:
                    column = fh.readline().rstrip("\n").split("\t").index("gtdb_taxonomy")
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) > column and f[column].startswith("d__"):
                        lineages.setdefault(accession_key(f[0]), f[column])
    return lineages


def accession_key(text):
    found = re.search(r"GC[AF]_(\d{9})", text)
    return f"GCA_{found.group(1)}" if found else None


def sylph_profiles(tsvs, sample_of_file, lineages, profiles):
    """Write each sample's profile from sylph's tables: the genomes' Taxonomic_abundance summed per GTDB species.
    Returns the genomes the taxonomy lacks (named by their accession in the profile)."""
    species = {sample: {} for sample in sample_of_file.values()}
    unknown_genomes = set()
    for tsv in tsvs:
        for r in read_tsv(tsv):
            sample = sample_of_file.get(os.path.abspath(r["Sample_file"])) or sample_of_file.get(r["Sample_file"])
            if sample is None:
                continue
            key = accession_key(r["Genome_file"])
            lineage = lineages.get(key)
            if lineage is None:
                unknown_genomes.add(r["Genome_file"])
                name = re.search(r"GC[AF]_\d{9}\.\d+", r["Genome_file"])
                lineage = f"s__{name.group(0) if name else os.path.basename(r['Genome_file'])}"
            entry = species[sample].setdefault(lineage.split(";")[-1], [0.0, key or "-", lineage])
            entry[0] += float(r["Taxonomic_abundance"]) / 100
    os.makedirs(profiles, exist_ok=True)
    for sample, rows in species.items():
        with open(os.path.join(profiles, sample + ".profile.part"), "w") as fh:
            for share, key, lineage in sorted(rows.values(), key=lambda v: -v[0]):
                fh.write(f"{key}\t{lineage}\t{share:.6g}\n")
        os.replace(os.path.join(profiles, sample + ".profile.part"), os.path.join(profiles, sample + ".profile"))
    return unknown_genomes


def run_sylph(args, folder, todo, lineages):
    merged = work_folder(args, folder, "joined")
    paired, single = [], []
    for s in todo:
        if args.dry_run:  # the commands only: no files joined
            first, second = s["first"].split(",")[0], s["second"].split(",")[0]
        else:
            first = joined(s["first"], os.path.join(merged, s["sample"] + "_1.fq.gz"))
            second = s["second"] if s["second"] == "-" else joined(s["second"],
                                                                     os.path.join(merged, s["sample"] + "_2.fq.gz"))
        if second != "-":
            paired.append((s["sample"], first, second))
        else:
            single.append((s["sample"], first))
    os.makedirs(os.path.join(folder, "sylph"), exist_ok=True)
    tsvs, code = [], 0
    sample_of_file = {}  # of the samples sylph profiled: a failed call's samples get no profile
    for kind, members in (("paired", paired), ("single", single)):
        if not members:
            continue
        tsv = os.path.join(folder, "sylph", f"{kind}.tsv")
        command = [args.sylph, "profile", args.sylph_db]
        if kind == "paired":
            command += ["-1"] + [m[1] for m in members] + ["-2"] + [m[2] for m in members]
        else:
            command += [m[1] for m in members]
        command += ["-t", str(args.threads), "-o", tsv] + shlex.split(args.sylph_args)
        print(f"  {' '.join(shlex.quote(c) for c in command[:4])} ... ({len(members)} {kind} samples)",
              file=sys.stderr)
        if args.dry_run:
            continue
        result, seconds = run_command(command, folder, "sylph.log")
        print(f"  sylph exit {result} after {seconds:.0f} s (log: {folder}/sylph.log)", file=sys.stderr)
        code = code or result
        if result == 0:
            tsvs.append(tsv)
            for m in members:
                sample_of_file[os.path.abspath(m[1])] = m[0]
    if tsvs:
        unknown = sylph_profiles(tsvs, sample_of_file, lineages, os.path.join(folder, "profiles"))
        if unknown:
            print(f"  {len(unknown)} genome(s) of sylph's database not in the --gtdb taxonomy, e.g. "
                  f"{sorted(unknown)[0]}: kept under their accession", file=sys.stderr)
    shutil.rmtree(merged, ignore_errors=True)
    return code


def sylph_run_name(args, out):
    """sylph's run folder for --sylph-db: the one OUT already has for that database (the most profiles), whatever sylph
    is installed now, so that each sample is profiled once and its profile reused by every later run and score; else,
    and with --rerun, sylph_<installed version>_<database>."""
    database = database_name(args.sylph_db)
    existing = [] if args.rerun else sylph_folders(out, database)
    if existing:
        folder = max(existing, key=lambda f: (len(profiles_of(f)), os.path.getmtime(f)))
        name = os.path.basename(folder)
        print(f"sylph on {database}: the run folder {name} (sylph {tool_and_version(name, read_versions(folder))[1]}, "
              f"{len(profiles_of(folder))} profiles) is reused; only samples without a profile there are run",
              file=sys.stderr)
        return name
    return f"sylph_{sylph_version(args.sylph)}_{database}"


def run(args):
    downloaded_path = os.path.join(args.reads, "reads", "downloaded.tsv")
    if not os.path.exists(downloaded_path):
        sys.exit(f"{downloaded_path} not found: run download_mocks.py first")
    samples = chosen(read_tsv(downloaded_path), [s for s in args.sets.split(",") if s])
    if args.tool == "protal" and not args.db:
        sys.exit("--tool protal needs --db")
    lineages = None
    if args.tool == "sylph":
        if not args.sylph_db or not args.gtdb:
            sys.exit("--tool sylph needs --sylph-db and --gtdb (a folder of GTDB r226's bac120/ar53_taxonomy_r226.tsv"
                     ".gz or _metadata_, the truth's release; the truth tables themselves come with the script)")
        lineages = load_taxonomy([os.path.expanduser(f) for f in args.gtdb.split(",") if f])
    out = os.path.abspath(args.out or args.reads)
    if args.scratch:
        args.scratch = os.path.abspath(os.path.expanduser(args.scratch))
    version = None  # the installed tool's, asked for once a sample is to be run
    if args.run_name:
        name = args.run_name
    elif args.tool == "protal":
        version = protal_version(args.protal)
        name = version
    else:
        name = sylph_run_name(args, out)
    if name in NOT_RUNS or "/" in name:
        sys.exit(f"--run-name {name!r} is not a folder name a run can have")
    run_folder = os.path.join(out, name)
    print(f"run folder {run_folder}", file=sys.stderr)
    done = profiles_of(run_folder)
    held = read_versions(run_folder)
    chunk = args.chunk if args.chunk is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    if args.chunks > 1 and not 0 <= chunk < args.chunks:
        sys.exit(f"--chunk {chunk} is not in 0..{args.chunks - 1}")
    failed = []
    for set_name in dict.fromkeys(s["set"] for s in samples):
        members = [s for s in samples if s["set"] == set_name]
        if args.chunks > 1:
            members = members[chunk::args.chunks]
        todo = [s for s in members if args.rerun or s["sample"] not in done]
        if not todo:
            print(f"{set_name}: all {len(members)} samples have profiles", file=sys.stderr)
            continue
        if version is None and not args.dry_run:
            version = protal_version(args.protal) if args.tool == "protal" else sylph_version(args.sylph)
            before = sorted({r["version"] for r in held.values() if r["tool"] == args.tool} - {version})
            if before:
                print(f"note: {args.tool} {version} profiles the new samples; the run folder's others are of "
                      f"{', '.join(before)} (versions.tsv says which; --rerun runs them all again)", file=sys.stderr)
        release = release_of(version) if args.tool == "protal" else None
        if release and release < PAIRED_ONLY_BEFORE:
            unpaired = [s for s in todo if s["second"] == "-"]
            if unpaired:
                kinds = ", ".join(sorted({s["read_type"] for s in unpaired}))
                print(f"{set_name}: protal {version} reads paired-end reads only: {len(unpaired)} {kinds} sample(s) "
                      "left out", file=sys.stderr)
                todo = [s for s in todo if s["second"] != "-"]
            if not todo:
                continue
        folder = os.path.join(run_folder, set_name)
        if args.chunks > 1:
            folder = os.path.join(folder, f"chunk_{chunk}")
        os.makedirs(folder, exist_ok=True)
        print(f"{set_name}: {len(todo)} samples", file=sys.stderr)
        code = (run_protal(args, folder, todo, release) if args.tool == "protal"
                else run_sylph(args, folder, todo, lineages))
        if args.dry_run:
            continue
        written = [s for s in todo if os.path.exists(os.path.join(folder, "profiles", s["sample"] + ".profile"))]
        database = os.path.abspath(args.db) if args.tool == "protal" else database_name(args.sylph_db)
        note_versions(run_folder, [{"sample": s["sample"], "set": set_name, "tool": args.tool, "version": version,
                                    "database": database, "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                                    "host": platform.node()} for s in written])
        print(f"{set_name}: {len(written)} of {len(todo)} profiles written", file=sys.stderr)
        if code or len(written) < len(todo):
            failed.append(set_name)
    if failed:
        sys.exit(f"failed or incomplete: {', '.join(failed)}")


# ---------------------------------------------------------------------------------------------------------------
# score
# ---------------------------------------------------------------------------------------------------------------

def read_profile(path):
    """({GTDB species: share}, unknown share or NaN) of a protal .profile."""
    species = {}
    unknown = NAN
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith("#"):
                continue
            try:
                share = float(f[2])
            except ValueError:
                continue
            if f[0] == "?" and f[1] == "?":
                unknown = share
                continue
            name = f[1].split(";")[-1]
            species[name] = species.get(name, 0.0) + share
    return species, unknown


def read_truth(path):
    """(kind, {species: share, NaN for present without a share, 0 for must be absent}, share of '?')."""
    kind = None
    rows = {}
    unknown = 0.0
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    for f in lines:
        if f[0] == "#kind":
            kind = f[1]
    body = [f for f in lines if not f[0].startswith("#")]
    for f in body[1:]:
        share = NAN if f[1] == "NA" else float(f[1])
        if f[0] == "?":
            unknown = share
        else:
            rows[f[0]] = share
    if kind not in ("abundance", "closed", "dropout", "negative"):
        sys.exit(f"{path}: no #kind line of abundance, closed, dropout or negative")
    return kind, rows, unknown


def normalised(shares):
    total = sum(v for v in shares.values() if v == v and v > 0)
    return {k: v / total for k, v in shares.items() if v == v and v > 0} if total > 0 else {}


def bray_curtis(a, b):
    """Bray-Curtis dissimilarity of two compositions adding up to 1 (half their L1 distance)."""
    if not a and not b:
        return 0.0
    if not a or not b:
        return 1.0
    return sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b)) / 2


def ratio(n, d):
    return n / d if d else NAN


def score_sample(profile, unknown, kind, truth, truth_unknown):
    """The scores of one profile against its truth (columns of SAMPLE_COLUMNS), and its species rows."""
    called = set(profile)
    p = normalised(profile)
    s = {"called": len(called), "unknown": unknown, "unknown_expected": truth_unknown if kind else NAN}
    species_rows = []
    if kind in ("abundance", "closed"):
        members = {k for k, v in truth.items() if v != v or v > 0}
        tp, fp, fn = called & members, called - members, members - called
        s.update(tp=len(tp), fp=len(fp), fn=len(fn), precision=ratio(len(tp), len(called)),
                 recall=ratio(len(tp), len(members)), fp_share=sum(p.get(k, 0.0) for k in fp),
                 false_species=",".join(sorted(fp, key=lambda k: -p.get(k, 0.0))[:10]))
        s["f1"] = ratio(2 * s["tp"], 2 * s["tp"] + s["fp"] + s["fn"])
        if kind == "abundance":
            t = normalised(truth)
            for threshold, column in zip(THRESHOLDS, ("recall_0.1pct", "recall_1pct")):
                big = {k for k, v in t.items() if v >= threshold}
                s[column] = ratio(len(big & called), len(big))
            s["l1"] = sum(abs(p.get(k, 0.0) - t.get(k, 0.0)) for k in set(p) | set(t))
            s["bray_curtis"] = bray_curtis(p, t)
        for k in sorted(members | called):
            status = "TP" if k in tp else ("FP" if k in fp else "FN")
            species_rows.append({"species": k, "truth_share": truth.get(k, 0.0) if kind == "abundance"
                                 else ("member" if k in members else "-"), "profile_share": p.get(k, 0.0),
                                 "status": status})
    elif kind == "dropout":
        absent = {k for k, v in truth.items() if v == 0}
        present = {k for k, v in truth.items() if v != v or v > 0}
        s["absent_called"] = ",".join(sorted(absent & called)) or "-"
        s["present_missed"] = ",".join(sorted(present - called)) or "-"
        for k in sorted(absent | present):
            status = ("absent, called" if k in called else "absent") if k in absent else (
                "present, called" if k in called else "present, missed")
            species_rows.append({"species": k, "truth_share": 0.0 if k in absent else "kept",
                                 "profile_share": p.get(k, 0.0), "status": status})
    elif kind == "negative":
        s.update(fp=len(called), false_species=",".join(sorted(called, key=lambda k: -p.get(k, 0.0))[:10]))
        for k in sorted(called):
            species_rows.append({"species": k, "truth_share": 0.0, "profile_share": p.get(k, 0.0), "status": "FP"})
    return s, species_rows


def median(values):
    values = [v for v in values if isinstance(v, (int, float)) and v == v]
    return statistics.median(values) if values else NAN



def run_folders(out, names):
    """The run folders to score: those named, or every folder of OUT with profiles."""
    if names:
        folders = [os.path.join(out, n) for n in names]
        missing = [f for f in folders if not os.path.isdir(f)]
        if missing:
            sys.exit(f"no run folder {', '.join(missing)}")
        return folders
    return [os.path.join(out, n) for n in sorted(os.listdir(out))
            if n not in NOT_RUNS and os.path.isdir(os.path.join(out, n)) and profiles_of(os.path.join(out, n))]


def score_run(folder, samples, downloaded, truth_folder, truths):
    """Score one run folder's profiles; writes its scores/ and returns (sample rows, summary rows)."""
    profiles = profiles_of(folder)
    sample_rows, species_rows, compositions = [], [], {}
    for s in samples:
        if s["sample"] not in profiles:
            continue
        profile, unknown = read_profile(profiles[s["sample"]])
        kind, truth, truth_unknown = None, {}, NAN
        if s["truth"] != "-":
            if s["truth"] not in truths:
                truths[s["truth"]] = read_truth(os.path.join(truth_folder, s["truth"] + ".tsv"))
            kind, truth, truth_unknown = truths[s["truth"]]
        scores, rows = score_sample(profile, unknown, kind, truth, truth_unknown)
        d = downloaded.get(s["sample"], {})
        row = {k: NAN for k in SAMPLE_COLUMNS}
        row.update({"set": s["set"], "sample": s["sample"], "read_type": s["read_type"], "truth": s["truth"],
                    "kind": kind or "-", "group": s["group"], "reads": d.get("reads", "NA"),
                    "bases": d.get("bases", "NA"), "cut": d.get("cut", "NA"), "absent_called": "-",
                    "present_missed": "-", "false_species": "-"})
        row.update(scores)
        sample_rows.append(row)
        for r in rows:
            species_rows.append({"set": s["set"], "sample": s["sample"], **r})
        compositions[s["sample"]] = (s, normalised(profile))
    if not sample_rows:
        return [], []

    pairs = []
    by_group = collections.defaultdict(list)
    for name, (s, _) in compositions.items():
        if s["group"] != "-":
            by_group[(s["set"], s["group"])].append(name)
    group_rows = []
    for (set_name, group), names in sorted(by_group.items()):
        if len(names) < 2:
            continue
        bcs, jaccards = [], []
        for a, b in itertools.combinations(sorted(names), 2):
            pa, pb = compositions[a][1], compositions[b][1]
            bc = bray_curtis(pa, pb)
            union = set(pa) | set(pb)
            jaccard = ratio(len(set(pa) & set(pb)), len(union)) if union else 1.0
            bcs.append(bc)
            jaccards.append(jaccard)
            pairs.append({"set": set_name, "group": group, "sample_a": a, "sample_b": b,
                          "read_type_a": compositions[a][0]["read_type"],
                          "read_type_b": compositions[b][0]["read_type"], "bray_curtis": bc, "jaccard": jaccard})
        group_rows.append({"set": set_name, "group": group, "samples": len(names), "pairs": len(bcs),
                           "mean_bray_curtis": statistics.fmean(bcs), "max_bray_curtis": max(bcs),
                           "mean_jaccard": statistics.fmean(jaccards), "min_jaccard": min(jaccards)})

    summary = []
    keys = collections.OrderedDict()
    for r in sample_rows:
        keys.setdefault((r["set"], r["read_type"], r["kind"]), []).append(r)
    for (set_name, read_type, kind), rows in keys.items():
        summary.append({"set": set_name, "read_type": read_type, "kind": kind, "samples": len(rows),
                        "called": median([r["called"] for r in rows]), "tp": median([r["tp"] for r in rows]),
                        "fp": median([r["fp"] for r in rows]), "fn": median([r["fn"] for r in rows]),
                        "precision": median([r["precision"] for r in rows]),
                        "recall": median([r["recall"] for r in rows]), "f1": median([r["f1"] for r in rows]),
                        "recall_1pct": median([r["recall_1pct"] for r in rows]),
                        "fp_share": median([r["fp_share"] for r in rows]), "l1": median([r["l1"] for r in rows]),
                        "unknown": median([r["unknown"] for r in rows]),
                        "dropouts_called": sum(r["absent_called"] not in ("-", "NA") for r in rows)})

    scores = os.path.join(folder, "scores")
    os.makedirs(scores, exist_ok=True)
    write_tsv(os.path.join(scores, "samples.tsv"), SAMPLE_COLUMNS, sample_rows)
    write_tsv(os.path.join(scores, "species.tsv"), ["set", "sample", "species", "truth_share", "profile_share",
                                                    "status"], species_rows)
    write_tsv(os.path.join(scores, "replicates.tsv"), ["set", "group", "sample_a", "sample_b", "read_type_a",
                                                       "read_type_b", "bray_curtis", "jaccard"], pairs)
    write_tsv(os.path.join(scores, "groups.tsv"), ["set", "group", "samples", "pairs", "mean_bray_curtis",
                                                   "max_bray_curtis", "mean_jaccard", "min_jaccard"], group_rows)
    write_tsv(os.path.join(scores, "summary.tsv"), SUMMARY_COLUMNS, summary)
    missing = len(samples) - len(sample_rows)
    print(f"{os.path.basename(folder)}: {len(sample_rows)} samples scored"
          f"{f' ({missing} without a profile)' if missing else ''}, {len(pairs)} replicate pairs; tables in {scores}",
          file=sys.stderr)
    return sample_rows, summary


def score(args):
    out = os.path.abspath(args.out or args.reads)
    samples = chosen(read_tsv(args.samples), [s for s in args.sets.split(",") if s])
    downloaded = {}
    path = os.path.join(args.reads, "reads", "downloaded.tsv")
    if os.path.exists(path):
        downloaded = {r["sample"]: r for r in read_tsv(path)}
    names = [n for n in args.runs.split(",") if n]
    folders = run_folders(out, names)
    if names and not args.no_sylph:
        # sylph's runs, made once (run --tool sylph), are the comparison of every later run
        folders += [f for f in sylph_folders(out) if f not in folders and profiles_of(f)]
    truths = {}
    all_samples, all_summary = [], []
    for folder in folders:
        name = os.path.basename(folder)
        versions = read_versions(folder)
        tool, version = tool_and_version(name, versions)
        sample_rows, summary = score_run(folder, samples, downloaded, args.truth, truths)
        for r in sample_rows:
            own = versions.get(r["sample"])
            all_samples.append({"run": name, "tool": own["tool"] if own else tool,
                                "version": own["version"] if own else version, **r})
        for r in summary:
            found = sorted({versions[s["sample"]]["version"] for s in sample_rows
                            if s["sample"] in versions and (s["set"], s["read_type"], s["kind"]) ==
                            (r["set"], r["read_type"], r["kind"])})
            all_summary.append({"run": name, "tool": tool, "version": ",".join(found) or version, **r})
    if not all_samples:
        sys.exit(f"no profiles of the chosen samples in the run folders of {out}")
    # Side by side: each set's rows of all runs together, the runs in the order scored.
    order = {os.path.basename(f): i for i, f in enumerate(folders)}
    first = {}
    for i, r in enumerate(all_summary):
        first.setdefault((r["set"], r["read_type"], r["kind"]), i)
    all_summary.sort(key=lambda r: (first[(r["set"], r["read_type"], r["kind"])], order[r["run"]]))
    scores = os.path.join(out, "scores")
    os.makedirs(scores, exist_ok=True)
    write_tsv(os.path.join(scores, "samples_by_run.tsv"), ["run", "tool", "version"] + SAMPLE_COLUMNS, all_samples)
    write_tsv(os.path.join(scores, "summary_by_run.tsv"), ["run", "tool", "version"] + SUMMARY_COLUMNS, all_summary)
    print("\t".join(["run", "tool", "version"] + SUMMARY_COLUMNS))
    for r in all_summary:
        print("\t".join(str(fmt(r[c])) for c in ["run", "tool", "version"] + SUMMARY_COLUMNS))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ("run", "score"):
        p = sub.add_parser(name)
        p.add_argument("--reads", required=True, help="download_mocks.py's output folder (-o)")
        p.add_argument("--out", default="", help="where the run folders and scores/ go (default: --reads)")
        p.add_argument("--sets", default="", help="comma-separated sets (default: all)")
        p.add_argument("--samples", default=os.path.join(HERE, "samples.tsv"))
        if name == "run":
            p.add_argument("--tool", choices=("protal", "sylph"), default="protal")
            p.add_argument("--run-name", default="", help="the run folder's name (default: from the tool's "
                           "--version: <version>_<commit> for protal; for sylph the run folder of --sylph-db that OUT "
                           "has, else sylph_<version>_<database>)")
            p.add_argument("-t", "--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
            p.add_argument("--chunks", type=int, default=1, help="split each set into N chunks (SLURM arrays)")
            p.add_argument("--chunk", type=int, default=None, help="the chunk to run (default $SLURM_ARRAY_TASK_ID)")
            p.add_argument("--rerun", action="store_true", help="run samples that have profiles again (sylph: into "
                           "sylph_<installed version>_<database>, not the run folder it has)")
            p.add_argument("--dry-run", "--map-only", dest="dry_run", action="store_true",
                           help="write protal's maps and print the commands, run nothing")
            p.add_argument("--scratch", default="", help="a folder on the node's own disk (e.g. $TMPDIR) for "
                           "intermediate files: protal's SAMs (its --scratch_dir; removed once profiled, kept in "
                           "alignments/ only with --strains; before protal 0.8.3 in a folder there removed after the "
                           "run) and joined read files (default: none, SAMs kept in alignments/)")
            g = p.add_argument_group("protal")
            g.add_argument("--db", default="", help="the protal database")
            g.add_argument("--protal", default="protal", help="the protal binary (default: protal on PATH)")
            g.add_argument("--strains", action="store_true", help="write strain MSAs (default --no_strains)")
            g.add_argument("--protal-args", default="", help="more protal options, quoted")
            g = p.add_argument_group("sylph")
            g.add_argument("--sylph-db", default="", help="sylph's database (.syldb)")
            g.add_argument("--sylph", default="sylph", help="the sylph binary (default: sylph on PATH)")
            g.add_argument("--gtdb", default="", help="folder of a GTDB release's bac120/ar53_taxonomy_rNNN.tsv(.gz) "
                           "(or _metadata_): the truth tables' release (r226) first; comma-separated, later folders "
                           "name the genomes the earlier ones lack (e.g. r226,r232 for an r232 database)")
            g.add_argument("--sylph-args", default="", help="more sylph profile options, quoted")
        else:
            p.add_argument("--runs", default="", help="comma-separated run folders (default: all with profiles); "
                           "sylph's run folders are scored with them, for comparison (--no-sylph: not)")
            p.add_argument("--no-sylph", action="store_true", help="with --runs: the runs named alone")
            p.add_argument("--truth", default=os.path.join(HERE, "truth"), help="folder of the truth tables")
    args = ap.parse_args()
    if args.command == "run":
        run(args)
    else:
        score(args)


if __name__ == "__main__":
    main()
