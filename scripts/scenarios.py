#!/usr/bin/env python3
"""Scenarios: training and test samples like those of real studies, for protal's presence models.

The collector's design (collect_training_data.py) spans depths, community sizes and read setups so that a model
learns them all; it says little about how a model does on one kind of study. A scenario is such a kind: a community
(species per sample, the share of them the training database lacks, abundances, strains, congeners, and a host
genome's share of the reads) sequenced by several technologies, each at its depth, every technology sequencing the
same communities. collect_training_data.py --scenarios simulates and profiles a scenario's samples beside the
design's, with meta_scenario naming the scenario in every row of the tables. build_gtdb_database.py --scenarios puts
some in the training data (hold-in) and others, of another seed, in the test set (hold-out), and the trainer reports
F1, false positive and false negative rates per scenario and read type on both (machine_learning_cmdline.py, section
"Scenarios").

The presets (PRESETS; --scenario_file adds or changes them):

    gut           150-1,000 species, 5% of them lacking from the database; Illumina PE 150 bp at Q35, 20M read pairs;
                  PacBio HiFi and Nanopore reads of the same bases (6 Gb)
    moderate      1,000-5,000 species, 30% lacking (freshwater, marine, sludge: between the gut and the soil); Illumina
                  PE 150 bp (Q35) at 10M read pairs, Ultima at 10M reads, PacBio and Nanopore at the same bases (3 Gb)
    soil          3,000-11,000 species, 60% lacking; Ultima Genomics single-end 300 bp reads at Q25, 20M reads; Illumina
                  PE 150 bp (Q35), PacBio and Nanopore reads of the same bases (6 Gb)
    soil_shallow  soil communities at 5M Illumina read pairs, and PacBio and Nanopore reads of the same bases (1.5 Gb)
    host          90% of the reads human, the rest 2-50 bacterial and archaeal species of power-law abundances
                  (alpha 1: a rank-abundance line of slope -1 on log-log axes), 5% of them lacking; Illumina PE 150
                  (Q35) at 10M read pairs, Ultima at 10M reads, PacBio and Nanopore at the same bases (3 Gb)
The communities' evenness varies too: lognormal abundances of sigma 1.0, 1.5, 2.0 and 2.5, the samples' in turn.

These depths are each scenario's typical ones: a sample's own depth is drawn around them (depth_range, below).

How the parts are made:
- The depths: each sample's depth is the preset's times a factor from depth_range's LOW to HIGH (default DEPTH_RANGE,
  1/8 to 2: soil from 2.5M to 40M read pairs, shallow soil from 0.6M to 10M; depth_spread S, as before 2026-10-07, is
  1/S to S; 1 every sample at the preset's), log-uniform and stratified, so that a scenario's few samples spread over
  the range (depth_factors). Sample s of every technology has the same factor: the technologies still read the same
  communities at the same bases. With one depth per scenario, a sample's depth (the trainer's sample_log_fragments)
  told its samples apart and a model could learn each sample's own offset, which cross-validation by species does not
  see (r226 v13: gradient boosting fell from 0.937 on the shallow-soil hold-in samples to 0.881 on the hold-out ones;
  docs/claude/2026-10-06-r226-v13-soil). From half to twice the depth (2026-10-06 to -07) the samples were still
  narrow clusters per scenario, and a sample at the edge of its scenario's was scored like the design's samples of its
  depth (r226 v15: the pe shallow-soil sample held out whole at 0.903 against 0.929 with species held out;
  docs/claude/2026-10-07-r226-v15); below 1x the samples are cheaper, so the wider range costs about what 0.5-2x did.
- The species per sample: each sample's own count, log-uniform and stratified from the scenario's MIN to MAX like the
  depths but drawn apart from them (species_counts; simulate_metagenomes --species_per_sample N1,N2,...), so that a
  community's size does not follow its depth and few samples still span the range. The design's samples have at most
  ~4,600 taxa with reads and the soils' (9,000-11,000 species until 2026-10-07) 9,800 and more: nothing lay between.
- The share of species the database lacks: a scenario draws its species from a genome table of its own
  (OUT/scenarios/<name>/genomes.tsv, scenario_table): every species of the collection's table on one side of the
  split, the database's and the held-out ones (--novel_species), and a random part of the other side, so that a
  species drawn uniformly from it (as simulate_metagenomes draws them) is one the database lacks with the scenario's
  share. A sample's share then varies around it (hypergeometrically). The table must have more species than a sample
  takes: 10,000 species at 60% need 6,000 held-out species and 4,000 others, which the default download's 25,000
  species give (since 2026-10-05; 8,000 before). A scenario that does not fit is scaled down to what the table holds
  (fit_species: the largest sample a TABLE_MARGIN-th of the table), and the collector and the build say so; download
  more species for the full size (download_gtdb.py --species, --rep_only_species).
- Illumina reads at a quality: simulate_metagenomes's instrument profiles (IlluminaSimulator.h) have their own mean
  base quality by cycle; --mean_quality shifts each read's curve so that its qualities average the scenario's, and
  the errors follow the qualities.
- Ultima Genomics reads: single-end, their length from a gamma distribution (mean 300, SD 40), their errors mostly
  homopolymer length errors, made by hifi_reads.py's flow model (mean base quality per read, homopolymers from two
  bases) from templates the collector draws like those of long reads.
- A host: its genome (download_gtdb.py fetches the human one, T2T-CHM13v2.0, gzipped as NCBI serves it) is written
  once as plain sequence (OUT/host/host.seq, read by memory map in simulate_metagenomes) and gives host_share of a
  sample's reads (pe, se: read pairs and reads) or bases (pb, ont): the community is simulated at the rest of the
  depth, and simulate_metagenomes makes the host's paired-end reads (--host_pairs, after the community's, from
  fragments drawn from the genome), its long and Ultima reads from templates drawn among the community's by their
  share of the bases. A host's reads are in no truth file: they reach the profile only through spurious alignments.
"""

import gzip
import hashlib
import json
import math
import os
import random
import re


READ_TYPES = ("pe", "se", "pb", "ont")
# The defaults of a scenario's reads, by type: Illumina paired-end reads (an instrument profile of simulate_metagenomes
# at a mean base quality), Ultima Genomics single-end reads
# (hifi_reads.py's flow model, ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD); long reads take the collection's setups
# (--pb_setup, --ont_setup) unless a scenario gives its own.
ILLUMINA = {"type": "pe", "length": 150, "profile": "HSXt", "fragment_mean": 350, "fragment_sd": 50, "quality": 35}
ULTIMA = {"type": "se", "setup": "ultima:300:40:25:2"}

# The communities' evenness: the samples' lognormal sigmas in turn (simulate_metagenomes --distribution lognormal
# --pln_sigma).
EVENNESS = "lognormal:1.0,1.5,2.0,2.5"

PRESETS = {
    "gut": {
        "description": "human gut: 150-1,000 species, 5% of them lacking from the database; Illumina PE 150 bp Q35 at "
                       "20M read pairs, PacBio HiFi and Nanopore at the same bases (6 Gb)",
        "species": "150-1000", "novel_share": 0.05, "abundance": EVENNESS, "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ILLUMINA, "depth": 20_000_000}, {"type": "pb", "depth": 6_000_000_000},
                  {"type": "ont", "depth": 6_000_000_000}]},
    "moderate": {
        "description": "a moderately complex community (freshwater, marine, sludge): 1,000-5,000 species, 30% of them "
                       "lacking from the database; Illumina PE 150 bp Q35 at 10M read pairs, Ultima Genomics SE 300 bp "
                       "Q25 at 10M reads, PacBio HiFi and Nanopore at the same bases (3 Gb)",
        "species": "1000-5000", "novel_share": 0.3, "abundance": EVENNESS, "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ILLUMINA, "depth": 10_000_000}, {**ULTIMA, "depth": 10_000_000},
                  {"type": "pb", "depth": 3_000_000_000}, {"type": "ont", "depth": 3_000_000_000}]},
    "soil": {
        "description": "soil: 3,000-11,000 species, 60% of them lacking from the database; Ultima Genomics SE 300 bp Q25 "
                       "at 20M reads, Illumina PE 150 bp Q35, PacBio HiFi and Nanopore at the same bases (6 Gb)",
        "species": "3000-11000", "novel_share": 0.6, "abundance": EVENNESS, "strains": "0.3,0.1",
        "congeners": "0.25:2-5", "host_share": 0.0,
        "reads": [{**ULTIMA, "depth": 20_000_000}, {**ILLUMINA, "depth": 20_000_000},
                  {"type": "pb", "depth": 6_000_000_000}, {"type": "ont", "depth": 6_000_000_000}]},
    "soil_shallow": {
        "description": "shallow soil: soil communities at 5M Illumina PE 150 bp Q35 read pairs, PacBio HiFi and "
                       "Nanopore at the same bases (1.5 Gb)",
        "species": "3000-11000", "novel_share": 0.6, "abundance": EVENNESS, "strains": "0.3,0.1",
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
FIELDS = ("description", "species", "novel_share", "abundance", "strains", "congeners", "host_share", "reads",
          "depth_range", "depth_spread")
OPTIONAL_FIELDS = ("description", "depth_range", "depth_spread")
# A sample's depth is its scenario's times a factor from LOW to HIGH (depth_factors); depth_spread S is 1/S to S.
DEPTH_RANGE = (0.125, 2.0)
NAME = re.compile(r"[a-z][a-z0-9_]*")
# A scenario's table should hold this many times the species of its largest sample, or its samples share most of
# their species: it is said so (scenario_table), not refused.
TABLE_MARGIN = 1.5


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


def depth_range_of(name, d):
    """(LOW, HIGH) of a definition's sample depth factors: its depth_range ("LOW-HIGH" or [LOW, HIGH]), its depth_spread
    S (1/S to S), or DEPTH_RANGE; a ScenarioError for both or a range that is not 0 < LOW <= HIGH, HIGH/LOW <= 100."""
    if "depth_range" in d and "depth_spread" in d:
        raise ScenarioError(f"scenario {name}: give depth_range or depth_spread, not both")
    if "depth_spread" in d:
        try:
            spread = float(d["depth_spread"])
        except (TypeError, ValueError):
            raise ScenarioError(f"scenario {name}: depth_spread {d['depth_spread']!r} is not a number") from None
        if not 1 <= spread <= 10:
            raise ScenarioError(f"scenario {name}: depth_spread {spread} is not between 1 (every sample at the "
                                "scenario's depth) and 10")
        return 1 / spread, spread
    given = d.get("depth_range", DEPTH_RANGE)
    try:
        low, high = (float(x) for x in (given.split("-") if isinstance(given, str) else given))
    except (TypeError, ValueError):
        raise ScenarioError(f"scenario {name}: depth_range {given!r}: expected LOW-HIGH, e.g. 0.125-2") from None
    if not 0 < low <= high or high / low > 100:
        raise ScenarioError(f"scenario {name}: depth_range {given!r}: expected 0 < LOW <= HIGH, HIGH at most 100 LOW")
    return low, high


def check_definition(name, d):
    """A scenario's definition, checked and with its reads completed from the defaults (ILLUMINA, ULTIMA); a
    ScenarioError says what is wrong."""
    if not NAME.fullmatch(name):
        raise ScenarioError(f"scenario name {name!r}: lower-case letters, digits and _, starting with a letter")
    unknown = sorted(set(d) - set(FIELDS))
    missing = [f for f in FIELDS if f not in d and f not in OPTIONAL_FIELDS]
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
    out["depth_range"] = depth_range_of(name, d)
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


def stratified(rng, samples, low, high):
    """`samples` values from low to high, log-uniform and stratified: one in each of `samples` equal parts of the log
    range, at a random place within it, the parts in random order."""
    parts = list(range(samples))
    rng.shuffle(parts)
    return [low * (high / low) ** ((part + rng.random()) / samples) for part in parts]


def depth_factors(seed, name, samples, depth_range):
    """The depth factor of each of a scenario's `samples` (its depth times the factor is the sample's), from LOW to HIGH
    of depth_range (a number S: 1/S to S), log-uniform and stratified (stratified), so that few samples still spread
    over the range; drawn from the collection's seed and the scenario's name, so that a rerun (and simulate_metagenomes
    --test) gets the same, and another seed (the build's test set) others. All LOW when LOW is HIGH."""
    low, high = (1 / depth_range, depth_range) if isinstance(depth_range, (int, float)) else depth_range
    if samples <= 0:
        return []
    if high <= low:
        return [float(low)] * samples
    return stratified(random.Random(seed_of(seed, f"{name}:depth")), samples, float(low), float(high))


def species_counts(seed, name, samples, species):
    """The species of each of a scenario's `samples` (simulate_metagenomes --species_per_sample N1,N2,...), from MIN to
    MAX of its species (species_bounds), log-uniform and stratified like its depths, drawn apart from them (their own
    random stream), so that a community's size does not follow its depth; None when MIN is MAX or for one sample (the
    simulator draws it from the range)."""
    low, high = species_bounds(species)
    if samples <= 1 or low >= high:
        return None
    values = stratified(random.Random(seed_of(seed, f"{name}:species")), samples, low, high)
    return [min(high, max(low, round(v))) for v in values]


def sample_depths(depth, factors):
    """A depth (read pairs, reads or bases) times each sample's factor, at least 1."""
    return [max(1, round(depth * f)) for f in factors]


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
