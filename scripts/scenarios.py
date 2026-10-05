#!/usr/bin/env python3
"""Scenarios: training and test samples like those of real studies, for protal's presence models.

The collector's design (collect_training_data.py) spans depths, community sizes and read setups so that a model
learns them all; it says little about how a model does on one kind of study. A scenario is such a kind: a community
(species per sample, the share of them the training database lacks, abundances, strains, congeners, and a host
genome's share of the reads) sequenced by several technologies, each at its depth, every technology sequencing the
same communities. collect_training_data.py --scenarios simulates and profiles a scenario's samples beside the
design's, with meta_scenario naming the scenario in every row of the tables. build_gtdb_database.py --scenarios puts
some in the training data (hold-in) and others, of another seed, in the test set (hold-out), and the trainer reports
F1, false positive and false negative rates per scenario and read type on both (random_forest_cmdline.py, section
"Scenarios").

The presets (PRESETS; --scenario_file adds or changes them):

    gut           ~400 species, 5% of them lacking from the database; Illumina PE 150 bp at Q35, 20M read pairs;
                  PacBio HiFi and Nanopore reads of the same bases (6 Gb)
    soil          ~10,000 species, 60% lacking; Ultima Genomics single-end 300 bp reads at Q25, 20M reads; Illumina
                  PE 150 bp (Q35), PacBio and Nanopore reads of the same bases (6 Gb)
    soil_shallow  the soil communities at 5M Illumina read pairs, and PacBio and Nanopore reads of the same bases
                  (1.5 Gb)
    host          90% of the reads human, the rest 2-50 bacterial and archaeal species of power-law abundances
                  (alpha 1: a rank-abundance line of slope -1 on log-log axes), 5% of them lacking; Illumina PE 150
                  (Q35) at 10M read pairs, Ultima at 10M reads, PacBio and Nanopore at the same bases (3 Gb)

How the parts are made:
- The share of species the database lacks: a scenario draws its species from a genome table of its own
  (OUT/scenarios/<name>/genomes.tsv, scenario_table): every species of the collection's table on one side of the
  split, the database's and the held-out ones (--novel_species), and a random part of the other side, so that a
  species drawn uniformly from it (as simulate_metagenomes draws them) is one the database lacks with the scenario's
  share. A sample's share then varies around it (hypergeometrically). The table must have more species than a sample
  takes: 10,000 species at 60% need 6,000 held-out species and 4,000 others, more than the default download's pool
  (8,000 species) has. A scenario that does not fit is scaled down to what the table holds (fit_species: the largest
  sample a TABLE_MARGIN-th of the table), and the collector and the build say so; download more species for the full
  size (download_gtdb.py --rep_only_species).
- Illumina reads at a quality: ART's built-in profiles have their own mean base quality (HiSeq X TruSeq: Q40.2 for
  the first reads, Q37.9 for the second). ART shifts every quality, and the errors with it, by -qs and -qs2; the
  shifts that give the target mean (art_shifts) come from a short ART run on a random sequence, kept in
  OUT/scenarios/art_quality.json.
- Ultima Genomics reads: single-end, their length from a gamma distribution (mean 300, SD 40), their errors mostly
  homopolymer length errors, made by hifi_reads.py's flow model (mean base quality per read, homopolymers from two
  bases) from templates the collector draws like those of long reads.
- A host: its genome (download_gtdb.py fetches the human one, T2T-CHM13v2.0, gzipped as NCBI serves it) is written
  once as plain sequence (OUT/host/host.seq, read by memory map; Host) and gives host_share of a sample's reads (pe,
  se: read pairs and reads) or bases (pb, ont): the community is simulated at the rest of the depth, and the host's
  paired-end reads are made by ART in amplicon mode from fragments drawn from the genome (host_pe_chunk), its long
  and Ultima reads from templates drawn among the community's by their share of the bases. A host's reads are in no
  truth file: they reach the profile only through spurious alignments.
"""

import gzip
import hashlib
import json
import math
import mmap
import os
import random
import re
import shutil
import subprocess
import tempfile
import threading
import bisect

import compressed

READ_TYPES = ("pe", "se", "pb", "ont")
# The defaults of a scenario's reads, by type: Illumina paired-end reads (ART), Ultima Genomics single-end reads
# (hifi_reads.py's flow model, ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD); long reads take the collection's setups
# (--pb_setup, --ont_setup) unless a scenario gives its own.
ILLUMINA = {"type": "pe", "length": 150, "profile": "HSXt", "fragment_mean": 350, "fragment_sd": 50, "quality": 35}
ULTIMA = {"type": "se", "setup": "ultima:300:40:25:2"}

PRESETS = {
    "gut": {
        "description": "human gut: ~400 species, 5% of them lacking from the database; Illumina PE 150 bp Q35 at 20M "
                       "read pairs, PacBio HiFi and Nanopore at the same bases (6 Gb)",
        "species": "350-450", "novel_share": 0.05, "abundance": "lognormal:1.5,2.0", "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ILLUMINA, "depth": 20_000_000}, {"type": "pb", "depth": 6_000_000_000},
                  {"type": "ont", "depth": 6_000_000_000}]},
    "soil": {
        "description": "soil: ~10,000 species, 60% of them lacking from the database; Ultima Genomics SE 300 bp Q25 at "
                       "20M reads, Illumina PE 150 bp Q35, PacBio HiFi and Nanopore at the same bases (6 Gb)",
        "species": "9000-11000", "novel_share": 0.6, "abundance": "lognormal:1.5,2.0", "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ULTIMA, "depth": 20_000_000}, {**ILLUMINA, "depth": 20_000_000},
                  {"type": "pb", "depth": 6_000_000_000}, {"type": "ont", "depth": 6_000_000_000}]},
    "soil_shallow": {
        "description": "shallow soil: the soil communities at 5M Illumina PE 150 bp Q35 read pairs, PacBio HiFi and "
                       "Nanopore at the same bases (1.5 Gb)",
        "species": "9000-11000", "novel_share": 0.6, "abundance": "lognormal:1.5,2.0", "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ILLUMINA, "depth": 5_000_000}, {"type": "pb", "depth": 1_500_000_000},
                  {"type": "ont", "depth": 1_500_000_000}]},
    "host": {
        "description": "high host contamination: 90% human reads, 2-50 bacterial and archaeal species of power-law "
                       "abundances (5% lacking from the database); Illumina PE 150 bp Q35 at 10M read pairs, Ultima "
                       "SE 300 bp Q25 at 10M reads, PacBio HiFi and Nanopore at the same bases (3 Gb)",
        "species": "2-50", "novel_share": 0.05, "abundance": "powerlaw:1.0", "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.9,
        "reads": [{**ILLUMINA, "depth": 10_000_000}, {**ULTIMA, "depth": 10_000_000},
                  {"type": "pb", "depth": 3_000_000_000}, {"type": "ont", "depth": 3_000_000_000}]},
}
FIELDS = ("description", "species", "novel_share", "abundance", "strains", "congeners", "host_share", "reads")
NAME = re.compile(r"[a-z][a-z0-9_]*")
# A scenario's table should hold this many times the species of its largest sample, or its samples share most of
# their species: it is said so (scenario_table), not refused.
TABLE_MARGIN = 1.5
# Host paired-end reads are made in chunks of this many pairs side by side (host_pe_chunk).
HOST_PAIRS_CHUNK = 1_000_000


class ScenarioError(ValueError):
    """A scenario that cannot be simulated as defined, with why."""


def species_bounds(text):
    """(MIN, MAX) of a scenario's species per sample, N or MIN-MAX."""
    lo, _, hi = str(text).partition("-")
    try:
        lo, hi = int(lo), int(hi or lo)
    except ValueError:
        raise ScenarioError(f"species {text!r}: expected N or MIN-MAX") from None
    if not 1 <= lo <= hi:
        raise ScenarioError(f"species {text!r}: expected 1 <= MIN <= MAX")
    return lo, hi


def check_definition(name, d):
    """A scenario's definition, checked and with its reads completed from the defaults (ILLUMINA, ULTIMA); a
    ScenarioError says what is wrong."""
    if not NAME.fullmatch(name):
        raise ScenarioError(f"scenario name {name!r}: lower-case letters, digits and _, starting with a letter")
    unknown = sorted(set(d) - set(FIELDS))
    missing = [f for f in FIELDS if f not in d and f != "description"]
    if unknown or missing:
        raise ScenarioError(f"scenario {name}: " + "; ".join(([f"unknown fields {', '.join(unknown)}"] if unknown else []) +
                                                           ([f"missing {', '.join(missing)}"] if missing else [])))
    out = {"description": str(d.get("description", "")), "species": str(d["species"]),
           "abundance": str(d["abundance"] or ""), "strains": str(d["strains"] or ""),
           "congeners": str(d["congeners"] or "0")}
    species_bounds(out["species"])
    for field in ("novel_share", "host_share"):
        try:
            out[field] = float(d[field])
        except (TypeError, ValueError):
            raise ScenarioError(f"scenario {name}: {field} {d[field]!r} is not a number") from None
    if not 0 <= out["novel_share"] <= 1:
        raise ScenarioError(f"scenario {name}: novel_share {out['novel_share']} is not between 0 and 1")
    if not 0 <= out["host_share"] < 1:
        raise ScenarioError(f"scenario {name}: host_share {out['host_share']} is not at least 0 and below 1")
    if not isinstance(d["reads"], list) or not d["reads"]:
        raise ScenarioError(f"scenario {name}: reads is a list of the read types' entries")
    reads, seen = [], set()
    for entry in d["reads"]:
        kind = entry.get("type") if isinstance(entry, dict) else None
        if kind not in READ_TYPES:
            raise ScenarioError(f"scenario {name}: a reads entry needs a type of {', '.join(READ_TYPES)}, got {entry!r}")
        if kind in seen:
            raise ScenarioError(f"scenario {name}: two reads entries of type {kind}")
        seen.add(kind)
        entry = {**(ILLUMINA if kind == "pe" else ULTIMA if kind == "se" else {}), **entry}
        try:
            entry["depth"] = int(float(entry["depth"]))
        except (KeyError, TypeError, ValueError):
            raise ScenarioError(f"scenario {name}: the {kind} reads need a depth (read pairs, reads or bases)") from None
        if entry["depth"] <= 0:
            raise ScenarioError(f"scenario {name}: the {kind} reads' depth must be positive")
        if kind == "pe":
            for field in ("length", "fragment_mean", "fragment_sd"):
                entry[field] = int(entry[field])
            if entry.get("quality") is not None:
                entry["quality"] = float(entry["quality"])
        elif kind == "se" and not str(entry.get("setup", "")).startswith("ultima:"):
            raise ScenarioError(f"scenario {name}: single-end reads are Ultima reads (setup ultima:LENGTH_MEAN:LENGTH_SD:"
                                f"Q_MEAN:Q_SD), got {entry.get('setup')!r}")
        reads.append(entry)
    out["reads"] = reads
    return out


def definitions(path=None):
    """The scenarios by name: the presets, and those of a JSON file ({name: definition}; a definition of a preset's
    name changes only the fields it gives, e.g. {"soil": {"species": "4000-5000"}})."""
    found = {name: dict(d) for name, d in PRESETS.items()}
    if path:
        try:
            with open(path) as fh:
                given = json.load(fh)
        except (OSError, ValueError) as e:
            raise ScenarioError(f"--scenario_file {path}: {e}") from None
        if not isinstance(given, dict):
            raise ScenarioError(f"--scenario_file {path}: expected an object of scenarios by name")
        for name, d in given.items():
            if not isinstance(d, dict):
                raise ScenarioError(f"--scenario_file {path}: scenario {name} is not an object")
            found[name] = {**found.get(name, {}), **d}
    return {name: check_definition(name, d) for name, d in found.items()}


def selection(text, samples, defs, keep_zero=False):
    """[(name, samples)] of --scenarios NAME[:SAMPLES],... ("all": every scenario defined), each with `samples`
    unless it gives its own, those of 0 samples left out (unless keep_zero); [] for "" or "none"."""
    text = (text or "").strip()
    if text.lower() in ("", "none"):
        return []
    out = []
    for item in text.split(","):
        name, _, n = item.strip().partition(":")
        names = list(defs) if name == "all" else [name]
        for one in names:
            if one not in defs:
                raise ScenarioError(f"--scenarios: unknown scenario {one!r} (defined: {', '.join(defs)})")
            try:
                count = int(n) if n else samples
            except ValueError:
                raise ScenarioError(f"--scenarios {item!r}: expected NAME or NAME:SAMPLES") from None
            if count < 0:
                raise ScenarioError(f"--scenarios {item!r}: samples cannot be negative")
            if any(o == one for o, _ in out):
                raise ScenarioError(f"--scenarios names {one} twice")
            out.append((one, count))
    return [(name, n) for name, n in out if n > 0 or keep_zero]


def seed_of(seed, name):
    """A scenario's own seed: its samples do not change when the design or the other scenarios do."""
    return seed * 1_000_003 + int(hashlib.sha1(name.encode()).hexdigest()[:7], 16)


# ---- the scenario's genome table --------------------------------------------------------------------------

def lineage_species(line):
    """The species (s__Genus species) of a genome table line, from its GTDB lineage; None if it has none."""
    lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
    return lineage.split(";")[-1].strip() if lineage else None


def table_split(known, novel, share):
    """How many of the `known` species (the database has them) and of the `novel` ones (it lacks them) a scenario's
    table takes, so that a species drawn from it is novel with probability `share`: all of the side that is short
    and as many of the other as the share asks. -> (known, novel)."""
    if share <= 0:
        return known, 0
    if share >= 1:
        return 0, novel
    if novel < share * (known + novel):
        return min(known, round(novel * (1 - share) / share)), novel
    return known, min(novel, round(known * share / (1 - share)))


def fit_species(definition, known, novel):
    """A scenario's definition with its species per sample scaled down to what a table of `known` and `novel` species
    holds at its share (table_split), when its largest sample would not fit: the largest then a TABLE_MARGIN-th of
    the table (so that samples do not all take the same species), the smallest in proportion. -> (definition, None),
    or (the scaled definition, what was scaled and why); a ScenarioError if the table would be empty."""
    lo, hi = species_bounds(definition["species"])
    k, n = table_split(known, novel, definition["novel_share"])
    if k + n >= hi:
        return definition, None
    top = int((k + n) / TABLE_MARGIN)
    if top < 1:
        raise ScenarioError(f"{definition['novel_share']:.0%} of its species lacking from the database: the genome table "
                            f"has {novel} species the training database lacks and {known} it has, no table at that share")
    bottom = max(1, min(top, round(lo * top / hi)))
    scaled = {**definition, "species": f"{bottom}-{top}" if bottom < top else str(top)}
    return scaled, (f"scaled from {definition['species']} to {scaled['species']} species per sample: at "
                    f"{definition['novel_share']:.0%} lacking from the database the genome table's {novel} species the "
                    f"training database lacks and {known} it has make a table of {k + n} (for the full size: simulate "
                    "from more species, download_gtdb.py --rep_only_species, or hold more out, --holdout)")


def table_plan(known, novel, share, largest):
    """table_split, and a ScenarioError if the table would hold fewer species than `largest`, a sample's most."""
    k, n = table_split(known, novel, share)
    if k + n < largest:
        need_novel, need_known = math.ceil(largest * share), math.ceil(largest * (1 - share))
        raise ScenarioError(
            f"{largest} species per sample, {share:.0%} of them lacking from the database, need at least "
            f"{need_novel} species the training database leaves out and {need_known} it has; the genome table has "
            f"{novel} and {known}, so the scenario's table would hold {k + n}. Simulate from more species "
            f"(download_gtdb.py --rep_only_species; build_gtdb_database.py --holdout leaves more out), or define the "
            f"scenario with fewer species (--scenario_file)")
    return k, n


def species_of_table(genome_table):
    """{species: [its lines]} of a genome table, in file order."""
    found = {}
    with open(genome_table) as fh:
        for line in fh:
            species = lineage_species(line)
            if species:
                found.setdefault(species, []).append(line if line.endswith("\n") else line + "\n")
    return found


def scenario_table(genome_table, novel, name, definition, seed, out):
    """Writes a scenario's genome table to `out` (only if its content changes, so that what was made from it stays
    valid): every genome of the species table_plan takes, chosen at random (seed_of) on the side of the split that
    has more than the share needs. novel: the species the database lacks. -> what it holds, in a few words; a
    ScenarioError if it cannot hold a sample."""
    species = species_of_table(genome_table)
    known = sorted(s for s in species if s not in novel)
    lacking = sorted(s for s in species if s in novel)
    largest = species_bounds(definition["species"])[1]
    try:
        k, n = table_plan(len(known), len(lacking), definition["novel_share"], largest)
    except ScenarioError as e:
        raise ScenarioError(f"scenario {name}: {e}") from None
    rng = random.Random(f"{seed_of(seed, name)}:table")
    chosen = set(rng.sample(known, k)) | set(rng.sample(lacking, n))
    text = "".join(line for s, lines in species.items() if s in chosen for line in lines)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    try:
        with open(out) as fh:
            same = fh.read() == text
    except OSError:
        same = False
    if not same:
        with open(out + ".partial", "w") as fh:
            fh.write(text)
        os.replace(out + ".partial", out)
    note = (f"{k + n} species ({n} the database lacks, {n / max(1, k + n):.1%}; {k} it has) for samples of "
            f"{definition['species']}")
    if k + n < TABLE_MARGIN * largest:
        note += f"; the samples share most of their species (a table of {k + n} for {largest})"
    return note


# ---- Illumina reads at a mean quality ---------------------------------------------------------------------

def mean_qualities(fastqs):
    """The mean base quality (Phred, +33) of each FASTQ."""
    out = []
    for path in fastqs:
        total = count = 0
        with open(path, "rb") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 3:
                    line = line.rstrip(b"\r\n")
                    total += sum(line) - 33 * len(line)
                    count += len(line)
        out.append(total / max(1, count))
    return out


def art_shifts(art, profile_args, length, fragment_mean, fragment_sd, target, cache):
    """ART's -qs and -qs2 (whole Phred steps) that bring its reads of this setup to a mean base quality of `target`:
    the profile's means measured once on a random sequence (300 kb, 5x, seed 1) and kept in the JSON file `cache`,
    with ART's identity. ART shifts every quality and draws the errors from the shifted ones, so the mean moves by
    the shift (but where it clips at Q0 or Q93). profile_args: ART's options of the profile (-ss HSXt, or -1 R1 -2 R2).
    -> (qs, qs2, (mean of the first reads, of the second) at no shift)."""
    exe = shutil.which(art) or art
    st = os.stat(exe)
    key = f"{os.path.realpath(exe)}|{st.st_size}|{st.st_mtime_ns}|{' '.join(profile_args)}|{length}|{fragment_mean}|{fragment_sd}"
    try:
        with open(cache) as fh:
            known = json.load(fh)
    except (OSError, ValueError):
        known = {}
    if key not in known:
        with tempfile.TemporaryDirectory() as tmp:
            rng = random.Random(1)
            with open(os.path.join(tmp, "g.fa"), "w") as fh:
                fh.write(">random\n" + "".join(rng.choice("ACGT") for _ in range(300_000)) + "\n")
            command = [exe, "-q", *profile_args, "-i", os.path.join(tmp, "g.fa"), "-p", "-l", str(length), "-f", "5",
                       "-m", str(fragment_mean), "-s", str(fragment_sd), "-na", "-rs", "1", "-o", os.path.join(tmp, "o")]
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise ScenarioError(f"{' '.join(command)} failed ({result.returncode}): {result.stderr.strip()[-300:]}")
            known[key] = mean_qualities([os.path.join(tmp, "o1.fq"), os.path.join(tmp, "o2.fq")])
        os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
        with open(cache + ".partial", "w") as fh:
            json.dump(known, fh, indent=1)
        os.replace(cache + ".partial", cache)
    first, second = known[key]
    return round(target - first), round(target - second), (first, second)


# ---- a host genome ----------------------------------------------------------------------------------------

def prepare_host(fasta, folder):
    """The host genome as one plain sequence file (folder/host.seq: every contig of 1 kb or more, upper case, one
    after the other) and its index (folder/host.json: contigs, offsets, lengths, and the FASTA it came from), written
    once per FASTA (its path, size and time). -> the folder."""
    st = os.stat(fasta)
    source = [os.path.realpath(fasta), st.st_size, st.st_mtime_ns]
    index = os.path.join(folder, "host.json")
    try:
        with open(index) as fh:
            if json.load(fh)["source"] == source and os.path.isfile(os.path.join(folder, "host.seq")):
                return folder
    except (OSError, ValueError, KeyError):
        pass
    os.makedirs(folder, exist_ok=True)
    contigs, offset = [], 0
    with open(fasta, "rb") as fh:
        opener = gzip.open if fh.read(2) == b"\x1f\x8b" else open
    with opener(fasta, "rb") as fin, open(os.path.join(folder, "host.seq.partial"), "wb") as out:
        name, start = None, 0

        def close(end):
            nonlocal offset
            if name is not None and end - start < 1000:  # a short contig: dropped (its bytes overwritten)
                out.seek(start)
                out.truncate()
                offset = start
            elif name is not None:
                contigs.append([name, start, end - start])
        for line in fin:
            if line.startswith(b">"):
                close(offset)
                name, start = line[1:].split()[0].decode() if len(line) > 2 else f"contig{len(contigs) + 1}", offset
                continue
            seq = line.strip().upper()
            out.write(seq)
            offset += len(seq)
        close(offset)
    if not contigs:
        raise ScenarioError(f"host genome {fasta}: no sequence of 1 kb or more")
    os.replace(os.path.join(folder, "host.seq.partial"), os.path.join(folder, "host.seq"))
    with open(index + ".partial", "w") as fh:
        json.dump({"source": source, "contigs": contigs, "bases": offset}, fh)
    os.replace(index + ".partial", index)
    return folder


def host_identity(folder):
    """What a host's reads are made from, for a design point's key: the FASTA (path, size, time) and its bases."""
    with open(os.path.join(folder, "host.json")) as fh:
        index = json.load(fh)
    return {"source": index["source"], "bases": index["bases"]}


class Host:
    """A prepared host genome (prepare_host), read by memory map: templates drawn at random positions, a contig by
    its length, a start uniform within it. One per folder and process (Host.of), shared by threads."""
    _open, _lock = {}, threading.Lock()

    def __init__(self, folder):
        with open(os.path.join(folder, "host.json")) as fh:
            index = json.load(fh)
        self.contigs = index["contigs"]
        self.ends = []
        total = 0
        for _, _, length in self.contigs:
            total += length
            self.ends.append(total)
        self.bases = total
        with open(os.path.join(folder, "host.seq"), "rb") as fh:  # the map keeps a descriptor of its own
            self.seq = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)

    @classmethod
    def of(cls, folder):
        with cls._lock:
            if folder not in cls._open:
                cls._open[folder] = cls(folder)
            return cls._open[folder]

    def draw(self, rng, length, min_length=100, max_n=0.1):
        """A template of `length` bases (shorter at a contig's end, but at least min_length or `length`: pbsim3 takes
        no shorter one, as the collector's templates of genomes have 100 bases at least) from a random place, on the
        forward strand; one too short or with more than max_n of N is drawn again (up to 50 times)."""
        least = min(length, min_length)
        for _ in range(50):
            at = rng.randrange(self.bases)
            k = bisect.bisect_right(self.ends, at)
            _, offset, contig_length = self.contigs[k]
            start = at - (self.ends[k] - contig_length)
            seq = self.seq[offset + start:offset + min(contig_length, start + length)]
            if len(seq) >= least and seq.count(b"N") <= max_n * len(seq):
                return seq
        return seq


COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def host_pe_chunk(task):
    """`task["pairs"]` host read pairs into task["r1"] and task["r2"] (zstd or gzip by their names, as the sample's
    reads they are appended to; compressed.open_write): fragments drawn from the host
    (fragment length normal, of the setup's mean and SD, at least the read length + 1; either strand), each read by
    ART in amplicon mode (-amp -p -c 1: one pair from the two ends of each fragment) with the setup's profile and
    quality shifts (task["art_args"]), named h<chunk>_<n>. -> None, or why it failed."""
    host = Host.of(task["host"])
    rng = random.Random(task["seed"])
    tmp = task["tmp"]
    os.makedirs(tmp, exist_ok=True)
    fragments = os.path.join(tmp, "fragments.fa")
    length = task["length"]
    with open(fragments, "wb") as fh:
        for i in range(task["pairs"]):
            size = max(length + 1, int(round(rng.gauss(task["fragment_mean"], task["fragment_sd"]))))
            seq = host.draw(rng, size, min_length=length + 1)  # longer than a read: ART reads both its ends
            if rng.random() < 0.5:
                seq = seq.translate(COMPLEMENT)[::-1]
            fh.write(b">h%d_%d\n%s\n" % (task["chunk"], i + 1, seq))
    prefix = os.path.join(tmp, "r")
    command = [task["art"], "-q", "-amp", "-p", *task["art_args"], "-i", fragments, "-l", str(length), "-c", "1",
               "-na", "-rs", str(task["seed"] % 2_000_000_000), "-o", prefix]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        return f"{task['sample']}: host reads: {' '.join(command)} failed ({result.returncode}): {result.stderr.strip()[-300:]}"
    made = []
    for source, dest in ((prefix + "1.fq", task["r1"]), (prefix + "2.fq", task["r2"])):
        with open(source, "rb") as fin, compressed.open_write(dest + ".partial") as fout:
            shutil.copyfileobj(fin, fout, 16 << 20)
        with open(source, "rb") as fh:
            made.append(sum(1 for _ in fh) // 4)
        os.replace(dest + ".partial", dest)
    shutil.rmtree(tmp, ignore_errors=True)
    if made != [task["pairs"]] * 2:
        return f"{task['sample']}: host reads: ART made {made[0]} and {made[1]} reads of {task['pairs']} fragments"
    return None
