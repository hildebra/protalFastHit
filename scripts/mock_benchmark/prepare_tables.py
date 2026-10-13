#!/usr/bin/env python3
"""Build the tables of the benchmark on public mock communities (docs/claude/2026-10-13-public-mock-communities):

  samples   samples.tsv: the samples of the report's starting sets, from ENA's run metadata. One row per sample (a
            library; Ion Proton re-runs and a library split over two runs are one sample): set, sample, runs (comma-
            separated), read_type (pe, se, pb, ont), truth (a table in truth/, or - for none), group (samples to
            compare with each other, - for none) and a note. download_mocks.py fetches the runs, score_mocks.py runs
            protal on them and scores the profiles.
  truth     truth/<name>.tsv: what each mock holds, as GTDB species of the database's release. Header lines
            #name, #kind, #source, #basis, then species, share, members. kind is
              abundance  shares known (genome copies: protal's abundances are shares of cells)
              closed     the members known, not their shares (share NA): any other species is a false positive
              dropout    only the rows matter: share NA must be there, share 0 must not (strains left out by design)
              negative   no species expected
            Shares add up to 1 with the "?" row: members protal cannot place (fungi, strains the release lacks), as
            protal's own "?" line. A DNA-mass design is turned into genome copies by genome size.

GTDB species come from the release's metadata (bac120/ar53_metadata_rNNN.tsv.gz): a member's genome accession first,
else its culture-collection number in the genomes' strain identifiers (majority species), else a named strain (organism
name and strain), else a species name the release has. Every assignment and its evidence is printed, so a newer
release can be checked line by line. Sources: the vendors' manuals (Zymo D6331 v1.2.0, D6322; ATCC MSA-1003), Amos et
al. 2020 Table 1 (NIBSC), Meslier et al. 2022 Table S1, Treichel et al. 2026 Supplementary Tables 1-2, Jin et al. 2023
Supplementary Data 1 (hCom2); the vendor and NIBSC numbers are copied below, the others downloaded. Standard library
only.

Usage:
  python3 scripts/mock_benchmark/prepare_tables.py samples -o scripts/mock_benchmark/samples.tsv
  python3 scripts/mock_benchmark/prepare_tables.py truth --gtdb ~/GTDB/r226 -o scripts/mock_benchmark/truth
"""
import argparse
import collections
import csv
import gzip
import io
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = ("run_accession,sample_alias,sample_title,library_name,instrument_model,library_layout,read_count,"
              "base_count,fastq_ftp,submitted_ftp")
MESLIER_S1 = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-022-01762-z/MediaObjects/"
              "41597_2022_1762_MOESM1_ESM.xlsx")
TREICHEL_ESM = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-026-02334-2/MediaObjects/"
                "41564_2026_2334_MOESM{}_ESM.xlsx")  # 4: compositions, 5: genomes
JIN_ESM = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-39121-0/MediaObjects/"
           "41467_2023_39121_MOESM4_ESM.zip")  # supplementaryData/SupplementaryDataS1.csv: the hCom2 strains
SAMPLE_COLUMNS = ["set", "sample", "runs", "read_type", "truth", "group", "note"]


def fetch(url, tries=4):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                return response.read()
        except OSError as err:
            if attempt == tries - 1:
                raise
            print(f"  retry {url}: {err}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def ena_runs(project):
    """ENA's runs of a project as dicts (ENA_FIELDS), sorted by run accession."""
    query = urllib.parse.urlencode({"accession": project, "result": "read_run", "fields": ENA_FIELDS, "format": "tsv"})
    text = fetch(f"{ENA_FILEREPORT}?{query}").decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        sys.exit(f"ENA returned no runs for {project}")
    return sorted(rows, key=lambda r: r["run_accession"])


# ---------------------------------------------------------------------------------------------------------------
# samples.tsv
# ---------------------------------------------------------------------------------------------------------------

def read_type(run):
    model = run["instrument_model"].lower()
    if any(k in model for k in ("promethion", "gridion", "minion")):
        return "ont"
    if any(k in model for k in ("sequel", "revio", "vega", "rs ii")):
        return "pb"
    if "ion torrent" in model or "ug 100" in model:
        return "se"
    return "pe" if run["library_layout"] == "PAIRED" else "se"


def has_reads(run):
    """Whether a run has FASTQ to fetch: ENA's own, or submitted FASTQ (raw-signal runs have neither)."""
    if run["fastq_ftp"]:
        return True
    return any(re.search(r"\.(fastq|fq)(\.gz)?$", f) for f in run["submitted_ftp"].split(";") if f)


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def microbench_sample(run):
    """MicroBench's sample (ZymoGut, ZymoHMW, ...): its alias reads "ZymoGut metagenome", its title "ZymoGut"."""
    return (run["sample_alias"] or run["sample_title"]).split(" ")[0]


def row(set_name, sample, runs, rtype, truth, group, note=""):
    return {"set": set_name, "sample": safe(sample), "runs": ",".join(runs), "read_type": rtype, "truth": truth,
            "group": group or "-", "note": note.replace("\t", " ") or "-"}


def samples_d6331_dohlman():
    """Zymo D6331 DNA alone and spiked into human DNA (Dohlman et al. 2026), and the controls."""
    out = []
    for r in ena_runs("PRJNA1454942"):
        title, acc = r["sample_title"], r["run_accession"]
        if "DNA only, replicate" in title and "D6331" in title:
            out.append(row("d6331_dna", f"d6331_dna_{acc}", [acc], "pe", "d6331", "d6331_dna", title))
            continue
        m = re.search(r"microbe:host ratio 1:(?:10\^)?(\d+)", title)
        if m:
            exp = 0 if m.group(1) == "1" and "1:1," in title else int(m.group(1))
            out.append(row("d6331_host", f"d6331_host_{acc}", [acc], "pe", "d6331", f"host_1e{exp}", title))
        elif "Human DNA only" in title:
            out.append(row("d6331_host", f"d6331_host_{acc}", [acc], "pe", "negative", "host_only", title))
        elif "extraction negative" in title:
            out.append(row("d6331_host", f"d6331_host_{acc}", [acc], "pe", "negative", "extraction_negative", title))
        else:
            sys.exit(f"PRJNA1454942: unexpected sample {acc}: {title}")
    return out


def samples_meslier():
    """Meslier et al. 2022: three mocks on up to seven platforms; a library's runs make one sample."""
    libraries = collections.OrderedDict()
    for r in ena_runs("PRJEB52977"):
        libraries.setdefault((r["sample_alias"], r["library_name"]), []).append(r)
    out = []
    for (alias, library), runs in libraries.items():
        mock = int(re.search(r"MOCK_00(\d)", alias).group(1))
        platform = library.replace(alias, "").strip("_")
        out.append(row("meslier", f"meslier_mock{mock}_{platform}", [r["run_accession"] for r in runs],
                       read_type(runs[0]), f"meslier_mock{mock}", f"meslier_mock{mock}",
                       f"{runs[0]['instrument_model']}, {len(runs)} run(s) of library {library}"))
    return out


def samples_treichel():
    """Treichel et al. 2026: Mock-even-70, Mock-stag-24, Mock-stag-70 and ten 24-strain variants (NovaSeq; one ONT
    run)."""
    out = []
    for r in ena_runs("PRJEB83573"):
        alias, acc = r["sample_alias"], r["run_accession"]
        if not has_reads(r):
            continue
        if alias.startswith("Mock-even"):
            truth, group = "treichel_even70", "even70_bg" if "bgDNA" in alias else "even70"
        elif alias.startswith("Multi-coverage-binning_Mock"):
            n = int(re.search(r"Mock(\d+)$", alias).group(1))
            truth, group = f"treichel_mcb_mock{n}", f"mcb_mock{n}"
        elif alias.startswith("Mock-stag") and "large" in alias:
            truth, group = "treichel_stag70", "stag70"
        elif alias.startswith("Mock-stag"):
            truth, group = "treichel_stag24", "stag24_bg" if "bgDNA" in alias else "stag24"
        else:
            sys.exit(f"PRJEB83573: unexpected sample {acc}: {alias}")
        out.append(row("treichel", f"treichel_{acc}", [acc], read_type(r), truth, group, alias))
    return out


def samples_long_reads():
    out = []
    for r in ena_runs("PRJNA680590"):  # PacBio's three libraries of one D6331 DNA
        out.append(row("d6331_hifi", f"d6331_hifi_{r['run_accession']}", [r["run_accession"]], "pb", "d6331",
                       "d6331_hifi", r["library_name"] or r["sample_title"]))
    for r in ena_runs("PRJNA627939"):  # ATCC MSA-1003 HiFi; SRR11469785 (PRJNA597914) holds the same data
        if r["run_accession"] == "SRR11606871":
            out.append(row("atcc1003_hifi", "atcc1003_hifi_SRR11606871", ["SRR11606871"], "pb", "atcc_msa1003", "-",
                           "Sequel II HiFi, Hon et al. 2020"))
    for r in ena_runs("PRJNA1181840"):
        out.append(row("d6331_ont", f"d6331_ont_{r['run_accession']}", [r["run_accession"]], "ont", "d6331", "-",
                       "PromethION LSK114, Chen et al. 2025"))
    microbench = ena_runs("PRJEB85558")
    for r in microbench:  # D6331: one flow cell's reads called fast/hac/sup (v5.0.0) and sup v5.2.0
        if microbench_sample(r) == "ZymoGut" and has_reads(r) and (
                r["library_name"].endswith(".fastq.gz") or "bm5.2.0_sup" in r["library_name"]):
            out.append(row("d6331_ont", f"d6331_ont_{r['run_accession']}", [r["run_accession"]], "ont", "d6331",
                           "d6331_microbench", r["library_name"]))
    for r in ena_runs("PRJEB65462") + [r for r in microbench if microbench_sample(r) == "ZymoHMW"]:
        if has_reads(r) and "mod4mC" not in r["library_name"]:  # the modified-base BAMs repeat the sup calls
            note = r["library_name"]
            submitted = os.path.basename(r["submitted_ftp"].split(";")[0])
            if submitted and submitted != note:
                note += " " + submitted
            out.append(row("d6322_basecalls", f"d6322_{r['run_accession']}", [r["run_accession"]], "ont",
                           "zymo_d6322", "d6322_PAS01578", note))
    for r in ena_runs("PRJNA934154"):  # SRR23563655 holds fast, hac and sup calls of one run in one file
        if r["run_accession"] in ("SRR24893244", "SRR24893245", "SRR24893246"):
            out.append(row("d6322_basecalls", f"d6322_{r['run_accession']}", [r["run_accession"]], "ont",
                           "zymo_d6322", "d6322_gridion_5khz", r["library_name"]))
    return out


def samples_strain_mixes():
    """Afolayan et al. 2024: 2-3 strains of E. coli or K. pneumoniae, each mix on R9.4.1 and R10.4.1."""
    out = []
    for r in ena_runs("PRJEB82665"):
        alias = r["sample_alias"]
        group = re.sub(r"_(R9|R10)$", "", alias)
        group = re.sub(r"_\d+-\d+$", "", group)
        out.append(row("strain_mixes", f"strain_mixes_{r['run_accession']}", [r["run_accession"]], "ont",
                       "ecoli_kpneumoniae", group, f"{alias} ({r['library_name']})"))
    return out


def samples_nibsc():
    """MHRA/NIBSC whole-cell reagent in 22 labs x 44 methods x 5 replicates, negatives, and the DNA reagents."""
    out = []
    for r in ena_runs("PRJNA1337664"):
        title, acc = r["sample_title"], r["run_accession"]
        m = re.match(r"Lab (\d+) Method (M\d+) (Biological replicate \d+|Negative control)", title)
        if m and m.group(3).startswith("Biological"):
            truth, group = "nibsc_wc_gut", f"lab{m.group(1)}_{m.group(2)}"
        elif m:
            truth, group = "negative", "-"
        elif title.startswith("DNA Mix RR"):
            truth, group = "nibsc_gutmix", "dna_mix"
        elif title.startswith("DNA HiLo RR"):
            truth, group = "nibsc_guthilo", "dna_hilo"
        else:
            sys.exit(f"PRJNA1337664: unexpected sample {acc}: {title}")
        out.append(row("nibsc_multilab", f"nibsc_{acc}", [acc], "pe", truth, group, title))
    return out


def samples_hcom():
    out = []
    for r in ena_runs("PRJNA885585"):  # 117-strain hCom2-like community in vitro (Jin et al. 2023)
        library = r["library_name"]
        group = re.sub(r"(R\d+)[abc]$", r"\1", library)  # a culture's three libraries
        group = re.sub(r"_(Flex|KAPA)$", "", group)  # an inoculum's three libraries
        out.append(row("hcom2_invitro", f"hcom2_{r['run_accession']}", [r["run_accession"]], "pe", "hcom2_closed",
                       group, library))
    arms = {"hCom1a": "hcom1a_full", "ΔCh": "hcom1a_minus_ch", "ΔCs": "hcom1a_minus_cs",
            "ΔChΔCs": "hcom1a_minus_chcs"}
    for r in ena_runs("PRJNA978640"):  # hCom1a with designed dropouts (Wang et al. 2023)
        arm, kind = r["sample_title"].split("_")[0], r["sample_title"].split("_")[1]
        if arm not in arms:
            sys.exit(f"PRJNA978640: unexpected sample {r['run_accession']}: {r['sample_title']}")
        out.append(row("hcom1a_dropouts", f"hcom1a_{r['run_accession']}", [r["run_accession"]], "pe", arms[arm],
                       f"{arms[arm]}_{kind.lower()}", f"{r['sample_title']} ({r['sample_alias']})"))
    return out


def samples_plant():
    """Plant SynComs of 177-988 strains, mostly missing from GTDB (Selten et al. 2026): no truth, groups only."""
    out = []
    for r in ena_runs("PRJNA1131994"):
        alias = r["sample_alias"]
        parts = alias.split("_")
        group = f"input_{parts[1]}" if parts[0] == "Input" else "_".join(parts[:2])
        out.append(row("plant_syncoms", f"plant_{r['run_accession']}", [r["run_accession"]], "pe", "-", group,
                       f"{alias} ({r['library_name']})"))
    for r in ena_runs("PRJNA1449013"):
        alias, library = r["sample_alias"], r["library_name"]
        group = "input_AtSC_7x200" if re.fullmatch(r"S\d", alias or "") else "AtSC_7x200"
        out.append(row("plant_syncoms", f"plant_{r['run_accession']}", [r["run_accession"]], "pe", "-", group,
                       f"{alias or '-'} ({library})"))
    return out


def write_samples(path):
    rows = (samples_d6331_dohlman() + samples_meslier() + samples_treichel() + samples_long_reads()
            + samples_strain_mixes() + samples_nibsc() + samples_hcom() + samples_plant())
    names = collections.Counter(r["sample"] for r in rows)
    duplicates = [n for n, c in names.items() if c > 1]
    if duplicates:
        sys.exit(f"duplicate sample names: {duplicates[:5]}")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    per_set = collections.Counter(r["set"] for r in rows)
    for s, n in per_set.items():
        print(f"{s}\t{n} samples", file=sys.stderr)
    print(f"{len(rows)} samples written to {path}", file=sys.stderr)


# ---------------------------------------------------------------------------------------------------------------
# truth tables
# ---------------------------------------------------------------------------------------------------------------

# Zymo D6331 Gut Microbiome Standard, instruction manual v1.2.0: genome copy %. Genomes: Zymo's 2023 NCBI deposits
# (E. faecalis IP101412 and the yeasts are in Zymo's D6331.refseq.zip only).
D6331 = [
    ("Faecalibacterium prausnitzii AP34BHI", 14.77, "GCF_028743395.1"),
    ("Veillonella rogosae AC2811 AN NA 2", 19.94, "GCF_028743475.1"),
    ("Roseburia hominis OB EAV1 11 DCM", 12.43, "GCF_028743455.1"),
    ("Bacteroides fragilis OB EAV1 11 D6 FAA", 8.33, "GCF_028743275.1"),
    ("Prevotella corporis OB21 FMU 4", 6.26, "GCF_028743775.1"),
    ("Bifidobacterium adolescentis LMG 10502", 8.83, "GCF_028743295.1"),
    ("Fusobacterium nucleatum 2/1/50A", 7.53, "GCF_028743415.1"),
    ("Limosilactobacillus fermentum B-1840", 9.68, "GCF_030770375.1"),
    ("Clostridioides difficile P4D3A1-1", 1.10, "GCF_028743315.1"),
    ("Akkermansia muciniphila OB21 FAA NB 28", 1.62, "GCF_028743255.1"),
    ("Methanobrevibacter smithii DSM 861", 0.17, "GCF_028743435.1"),
    ("Salmonella enterica B-4212", 0.007, "GCF_028743635.1"),
    ("Enterococcus faecalis IP101412 AER FAA 2", 0.0011, "species:Enterococcus faecalis"),
    ("Clostridium perfringens OB21 TSA 19", 0.00009, "GCF_028743735.1"),
    ("Escherichia coli JM109", 1.82, "GCF_028743375.1"),
    ("Escherichia coli B-3008", 1.82, "GCF_028743355.1"),
    ("Escherichia coli B-2207", 1.64, "GCF_028743335.1"),
    ("Escherichia coli B-766", 1.66, "GCF_028743755.1"),
    ("Escherichia coli B-1109", 1.77, "GCF_028743555.1"),
    ("Candida albicans IHEM 3108", 0.31, "unplaced:fungus"),
    ("Saccharomyces cerevisiae Y-567", 0.32, "unplaced:fungus"),
]
# Zymo D6322 HMW DNA Standard: genome copy %; the strains of D6300 lots since ZRC190633 (NCBI CP1179xx, CP118021).
D6322 = [
    ("Pseudomonas aeruginosa B-3509", 7.8, "GCF_028743595.1"),
    ("Escherichia coli B-1109", 10.9, "GCF_028743555.1"),
    ("Salmonella enterica B-4212", 11.2, "GCF_028743635.1"),
    ("Enterococcus faecalis B-537", 18.8, "GCF_028743535.1"),
    ("Staphylococcus aureus B-41012", 19.6, "GCF_028743615.1"),
    ("Listeria monocytogenes B-33116", 17.8, "GCF_028743575.1"),
    ("Bacillus subtilis B-354 (NCBI: B. spizizenii)", 13.2, "GCF_028743795.1"),
    ("Saccharomyces cerevisiae Y-567", 0.63, "unplaced:fungus"),
]
# ATCC MSA-1003 20 Strain Staggered Mix Genomic Material: the product sheet's % (taken as genomic DNA, as by Portik et
# al. 2022: turned into genome copies by genome size).
ATCC_MSA1003 = [
    ("Escherichia coli ATCC 700926 (K-12 MG1655)", 18, "named:Escherichia coli|MG1655"),
    ("Porphyromonas gingivalis ATCC 33277", 18, "strain:ATCC 33277"),
    ("Cereibacter sphaeroides ATCC 17029", 18, "strain:ATCC 17029"),
    ("Staphylococcus epidermidis ATCC 12228", 18, "strain:ATCC 12228"),
    ("Streptococcus mutans ATCC 700610 (UA159)", 18, "named:Streptococcus mutans|UA159"),
    ("Bacillus pacificus ATCC 10987", 1.8, "strain:ATCC 10987"),
    ("Clostridium beijerinckii ATCC 35702", 1.8, "strain:ATCC 35702"),
    ("Pseudomonas paraeruginosa ATCC 9027", 1.8, "strain:ATCC 9027"),
    ("Staphylococcus aureus ATCC BAA-1556 (USA300 FPR3757)", 1.8, "named:Staphylococcus aureus|FPR3757"),
    ("Streptococcus agalactiae ATCC BAA-611 (2603 V/R)", 1.8, "named:Streptococcus agalactiae|2603V/R"),
    ("Acinetobacter baumannii ATCC 17978", 0.18, "strain:ATCC 17978"),
    ("Cutibacterium acnes ATCC 11828", 0.18, "strain:ATCC 11828"),
    ("Helicobacter pylori ATCC 700392 (26695)", 0.18, "named:Helicobacter pylori|26695"),
    ("Lactobacillus gasseri ATCC 33323", 0.18, "strain:ATCC 33323"),
    ("Neisseria meningitidis ATCC BAA-335 (MC58)", 0.18, "named:Neisseria meningitidis|MC58"),
    ("Phocaeicola vulgatus ATCC 8482", 0.02, "strain:ATCC 8482"),
    ("Bifidobacterium adolescentis ATCC 15703", 0.02, "strain:ATCC 15703"),
    ("Deinococcus radiodurans ATCC BAA-816 (R1)", 0.02, "named:Deinococcus radiodurans|R1"),
    ("Enterococcus faecalis ATCC 47077 (OG1RF)", 0.02, "strain:ATCC 47077"),
    # ATCC 17982's genome is NCBI's Schaalia dentiphila ATCC 17982 (the type of that species), not S. odontolytica
    ("Schaalia odontolytica ATCC 17982 (NCBI: S. dentiphila)", 0.02, "GCF_000154225.1"),
]
# NIBSC Gut-Mix-RR and Gut-HiLo-RR (DNA): genome copy % of Amos et al. 2020, Table 1; WC-Gut-RR: the same strains as
# whole cells at equal numbers (Amos et al. 2022). (name, Gut-Mix %, Gut-HiLo %, how to find the GTDB species)
NIBSC = [
    ("Akkermansia muciniphila DSM 22959", 6.37, 0.18, "strain:DSM 22959"),
    ("Alistipes finegoldii DSM 17242", 4.54, 1.30, "strain:DSM 17242"),
    ("Anaerostipes hadrus DSM 3319", 6.11, 1.75, "strain:DSM 3319"),
    ("Bacteroides thetaiotaomicron DSM 2079", 2.69, 7.72, "strain:DSM 2079"),
    ("Bacteroides uniformis DSM 6597", 3.66, 1.05, "GCF_000154205"),
    ("Bifidobacterium longum subsp. infantis DSM 20088", 6.00, 17.20, "strain:DSM 20088"),
    ("Bifidobacterium longum subsp. longum DSM 20219", 6.92, 19.82, "GCF_900104835"),
    ("Blautia wexlerae DSM 19850", 3.77, 0.11, "GCF_000484655"),
    ("Clostridium butyricum DSM 10702", 3.70, 10.59, "GCF_000409755"),
    ("Collinsella aerofaciens DSM 13712", 6.95, 1.99, "GCF_902501475"),
    ("Escherichia coli DSM 1103 (ATCC 25922)", 3.26, 9.33, "named:Escherichia coli|ATCC25922"),
    ("Eubacterium hallii DSM 3353", 5.16, 1.48, "GCF_000173975"),
    ("Faecalibacterium prausnitzii DSM 17677 (A2-165, now F. duncaniae)", 5.49, 0.16,
     "named:Faecalibacterium duncaniae|A2-165"),
    ("Lactobacillus gasseri DSM 20077 (ATCC 33323)", 8.97, 0.26, "strain:ATCC 33323"),
    ("Parabacteroides distasonis DSM 20701 (ATCC 8503)", 3.52, 10.10, "named:Parabacteroides distasonis|ATCC8503"),
    ("Prevotella copri DSM 18205", 4.83, 13.84, "GCF_000157935"),
    ("Prevotella melaninogenica DSM 7089 (ATCC 25845)", 5.34, 1.53, "named:Prevotella melaninogenica|ATCC25845"),
    ("Roseburia hominis DSM 16839 (A2-183)", 4.72, 1.35, "named:Roseburia hominis|A2-183"),
    ("Roseburia intestinalis DSM 14610 (L1-82)", 3.88, 0.11, "named:Roseburia intestinalis|L1-82"),
    ("Ruminococcus gauvreauii DSM 19829", 4.13, 0.12, "GCF_000425525"),
]
# hCom1a's dropout strains (Wang et al. 2023)
# Members whose genome the release lacks or whose source name it renamed: how to find them instead.
OVERRIDES = {
    "GCA_000009225.1": "named:Pseudomonas|SBW25",  # Meslier: P. fluorescens SBW25, its other assemblies
    # Treichel: no assembly; the source's GTDB-Tk name, B. massiliensis, is not in r226
    "Blautia fusiformis CLA-AA-11": "species:Blautia_A fusiformis",
}
HCOM1A_DROPOUTS = [("Clostridium hylemonae DSM 15053", "strain:DSM 15053"),
                   ("Clostridium scindens ATCC 35704", "strain:ATCC 35704")]


def base_accession(acc):
    return acc.split(".")[0].replace("GCF_", "GCA_")


class Gtdb:
    """The species, genome sizes and strain identifiers of a GTDB release's genomes, read in one pass."""

    def __init__(self, folder, strain_ids=(), named=()):
        files = sorted(os.path.join(folder, f) for f in os.listdir(folder)
                       if re.fullmatch(r"(bac120|ar53)_metadata_r\d+\.tsv\.gz", f))
        if not files:
            sys.exit(f"{folder}: no bac120/ar53_metadata_rNNN.tsv.gz")
        self.release = re.search(r"_r(\d+)\.tsv", files[0]).group(1)
        self.accession = {}  # base accession -> (species, genome size)
        self.sizes = collections.defaultdict(list)
        patterns = {s: re.compile(r"(^|[^0-9A-Za-z])" + re.escape(s).replace(r"\ ", r"[ _-]?") + r"($|[^0-9])")
                    for s in strain_ids}
        self.strain_hits = {s: collections.Counter() for s in strain_ids}
        self.named_hits = {n: collections.Counter() for n in named}
        for path in files:
            with gzip.open(path, "rt") as fh:
                head = fh.readline().rstrip("\n").split("\t")
                ix = {k: i for i, k in enumerate(head)}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    species = f[ix["gtdb_taxonomy"]].split(";")[-1]
                    try:
                        size = float(f[ix["genome_size"]])
                    except ValueError:
                        size = float("nan")
                    self.accession[base_accession(f[0][3:])] = (species, size)
                    self.sizes[species].append(size)
                    ids = f[ix["ncbi_organism_name"]] + " " + f[ix["ncbi_strain_identifiers"]]
                    for s, p in patterns.items():
                        if p.search(ids):
                            self.strain_hits[s][species] += 1
                    compact = ids.replace(" ", "")
                    for n in named:
                        organism, strain = n.split("|")
                        if f[ix["ncbi_organism_name"]].startswith(organism) and strain.replace(" ", "") in compact:
                            self.named_hits[n][species] += 1
        self.species_size = {s: sorted(v)[len(v) // 2] for s, v in self.sizes.items()}

    def lookup(self, how, fallback_species=None):
        """(GTDB species or None, evidence). how: an accession, strain:ID, named:Organism|strain, species:Name."""
        if how.startswith(("GCA_", "GCF_")):
            hit = self.accession.get(base_accession(how))
            if hit:
                return hit[0], f"genome {how}"
        elif how.startswith("strain:"):
            hits = self.strain_hits[how[7:]]
            if hits:
                return hits.most_common(1)[0][0], f"{how[7:]} in {sum(hits.values())} genome(s): {dict(hits)}"
        elif how.startswith("named:"):
            hits = self.named_hits[how[6:]]
            if hits:
                return hits.most_common(1)[0][0], f"{how[6:]} in {sum(hits.values())} genome(s): {dict(hits)}"
        elif how.startswith("species:"):
            name = "s__" + how[8:]
            if name in self.species_size:
                return name, "species name"
        elif how.startswith("unplaced:"):
            return None, how[9:]
        if fallback_species:
            name = fallback_species if fallback_species.startswith("s__") else "s__" + fallback_species
            if name in self.species_size:
                return name, f"{how} not in r{self.release}; the source's species name"
        return None, f"{how} not in r{self.release}"

    def genome_size(self, how, species):
        hit = self.accession.get(base_accession(how)) if how.startswith(("GCA_", "GCF_")) else None
        if hit and hit[1] == hit[1]:
            return hit[1]
        return self.species_size.get(species, float("nan"))


def queries(entries):
    strain = [h[7:] for _, h in entries if h.startswith("strain:")]
    named = [h[6:] for _, h in entries if h.startswith("named:")]
    return strain, named


def write_truth(folder, name, kind, source, basis, members):
    """members: [(member label, GTDB species or None, weight or None or 0)]; weights become shares adding to 1."""
    by_species = collections.OrderedDict()
    for label, species, weight in members:
        key = species or "?"
        entry = by_species.setdefault(key, [0.0, [], False])
        if weight is None:
            entry[2] = True
        else:
            entry[0] += weight
        entry[1].append(label)
    total = sum(w for w, _, na in by_species.values() if not na)
    path = os.path.join(folder, f"{name}.tsv")
    with open(path, "w") as fh:
        fh.write(f"#name\t{name}\n#kind\t{kind}\n#source\t{source}\n#basis\t{basis}\n")
        fh.write("species\tshare\tmembers\n")
        rows = sorted(by_species.items(), key=lambda kv: (kv[0] == "?", -kv[1][0], kv[0]))
        for species, (weight, labels, na) in rows:
            if na:
                share = "NA"
            elif weight == 0:
                share = "0"
            else:
                share = f"{weight / total:.6g}"
            fh.write(f"{species}\t{share}\t{'; '.join(labels)}\n")
    print(f"  {path}: {len(by_species)} rows", file=sys.stderr)


def place(gtdb, label, how, fallback=None):
    species, evidence = gtdb.lookup(how, fallback)
    print(f"    {label}\t{species or '?'}\t{evidence}", file=sys.stderr)
    return species


def read_xlsx(data):
    """The sheets of an .xlsx file: {sheet name: [rows of strings]}."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    z = zipfile.ZipFile(io.BytesIO(data))
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])))
    targets = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    sheets = {}
    for sheet in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", ns):
        target = targets[sheet.get("{%s}id" % ns["r"])].lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        rows = []
        for r in ET.fromstring(z.read(target)).iter("{%s}row" % ns["m"]):
            cells = {}
            for c in r.findall("m:c", ns):
                v = c.find("m:v", ns)
                if c.get("t") == "s" and v is not None:
                    value = shared[int(v.text)]
                elif c.get("t") == "inlineStr":
                    value = "".join(t.text or "" for t in c.iter("{%s}t" % ns["m"]))
                else:
                    value = v.text if v is not None else ""
                letters = re.match(r"[A-Z]+", c.get("r")).group(0)
                col = 0
                for ch in letters:
                    col = col * 26 + ord(ch) - 64
                cells[col - 1] = value.strip()
            if cells:
                rows.append([cells.get(j, "") for j in range(max(cells) + 1)])
        sheets[sheet.get("name")] = rows
    return sheets


def number(text):
    try:
        return float(text)
    except ValueError:
        return None


def write_truths(gtdb_folder, out):
    os.makedirs(out, exist_ok=True)
    meslier = read_xlsx(fetch(MESLIER_S1))
    treichel_mocks = read_xlsx(fetch(TREICHEL_ESM.format(4)))
    treichel_genomes = read_xlsx(fetch(TREICHEL_ESM.format(5)))
    with zipfile.ZipFile(io.BytesIO(fetch(JIN_ESM))) as z:
        jin = list(csv.DictReader(io.StringIO(z.read("supplementaryData/SupplementaryDataS1.csv").decode())))

    static = ([(n, h) for n, _, h in D6331] + [(n, h) for n, _, h in D6322] + [(n, h) for n, _, h in ATCC_MSA1003]
              + [(n, h) for n, _, _, h in NIBSC] + HCOM1A_DROPOUTS + list(OVERRIDES.items()))
    strain_ids, named = queries(static)
    print(f"reading GTDB metadata in {gtdb_folder}", file=sys.stderr)
    gtdb = Gtdb(gtdb_folder, strain_ids, named)
    release = f"GTDB r{gtdb.release}"

    print("d6331", file=sys.stderr)
    write_truth(out, "d6331", "abundance", "ZymoBIOMICS Gut Microbiome Standard D6331, manual v1.2.0", "genome copies",
                [(n, place(gtdb, n, h), w) for n, w, h in D6331])
    print("zymo_d6322", file=sys.stderr)
    write_truth(out, "zymo_d6322", "abundance", "ZymoBIOMICS HMW DNA Standard D6322, manual", "genome copies",
                [(n, place(gtdb, n, h), w) for n, w, h in D6322])
    print("atcc_msa1003", file=sys.stderr)
    members = []
    for n, w, h in ATCC_MSA1003:
        species = place(gtdb, n, h)
        members.append((f"{n} ({w} %)", species, w / gtdb.genome_size(h, species)))
    write_truth(out, "atcc_msa1003", "abundance", "ATCC MSA-1003 product sheet (%, taken as genomic DNA)",
                "genome copies from DNA mass by genome size (the species' median in " + release + ")", members)

    print("nibsc", file=sys.stderr)
    nibsc = [(n, place(gtdb, n, h), mix, hilo) for n, mix, hilo, h in NIBSC]
    write_truth(out, "nibsc_gutmix", "abundance", "NIBSC Gut-Mix-RR (WHO 20/302), Amos et al. 2020 Table 1",
                "genome copies", [(n, s, mix) for n, s, mix, _ in nibsc])
    write_truth(out, "nibsc_guthilo", "abundance", "NIBSC Gut-HiLo-RR (WHO 20/304), Amos et al. 2020 Table 1",
                "genome copies", [(n, s, hilo) for n, s, _, hilo in nibsc])
    write_truth(out, "nibsc_wc_gut", "abundance", "NIBSC WC-Gut-RR (whole cells, equal numbers), Amos et al. 2022",
                "cells", [(n, s, 1.0) for n, s, _, _ in nibsc])

    print("meslier", file=sys.stderr)
    rows = meslier["S1-Composition"]
    head = rows[0]
    col = {name: i for i, name in enumerate(head)}
    accession_col = next(i for i, h in enumerate(head) if h.startswith("Genbank Assembly"))
    mock_cols = [next(i for i, h in enumerate(head) if h.startswith(f"MOCK{k}")) for k in (1, 2, 3)]
    placed = []
    for r in rows[1:]:
        if len(r) <= mock_cols[2] or not r[accession_col].startswith("GCA_"):
            continue
        name = f"{r[0].strip()} {r[col['Strain Name']].strip()} ({r[accession_col]})"
        how = OVERRIDES.get(r[accession_col], r[accession_col])
        placed.append((name, place(gtdb, name, how, r[0].strip()), [number(r[c]) for c in mock_cols]))
    for k in range(3):
        write_truth(out, f"meslier_mock{k + 1}", "abundance", "Meslier et al. 2022, Table S1 (% genomes)",
                    "genome copies", [(n, s, shares[k]) for n, s, shares in placed if shares[k]])

    print("treichel", file=sys.stderr)
    genomes = {}
    for r in treichel_genomes["sMG_supptable2"][1:]:
        if len(r) >= 6 and r[5]:  # an assembly, or (one strain) only its ENA sample
            species = r[4].split(";")[-1] if ";" in r[4] else None
            key = re.sub(r"[\s_]+", " ", r[0]).strip()
            genomes[key] = (r[5], species, number(r[2]))
    treichel_placed = {}

    def treichel_member(label):
        key = re.sub(r"[\s_]+", " ", label).strip()
        if key not in genomes:
            sys.exit(f"Treichel: {label} has no genome in Supplementary Table 2")
        if key not in treichel_placed:
            acc, species, size = genomes[key]
            treichel_placed[key] = (place(gtdb, f"{key} ({acc})", OVERRIDES.get(key, acc), species), size)
        return treichel_placed[key]

    source = "Treichel et al. 2026, Supplementary Tables 1-2 (DNA ratio; genome sizes)"
    basis = "genome copies from DNA mass by genome size"
    for sheet, name in (("Mock-even-70", "treichel_even70"), ("Mock-stag-24", "treichel_stag24"),
                        ("Mock-stag-70", "treichel_stag70")):
        members = []
        for r in treichel_mocks[sheet][2:]:
            ratio = number(r[0]) if r else None
            if ratio is None or len(r) < 2 or not r[1]:
                continue
            species, size = treichel_member(r[1])
            members.append((r[1].replace("_", " "), species, ratio / size))
        write_truth(out, name, "abundance", source, basis, members)
    mcb = treichel_mocks["Mock-stag-24_multi_coverage_bin"]
    for k in range(1, 11):
        members = []
        for r in mcb[2:]:
            if len(r) > k and number(r[k]) is not None:
                species, size = treichel_member(r[0])
                members.append((r[0].replace("_", " "), species, number(r[k]) / size))
        write_truth(out, f"treichel_mcb_mock{k}", "abundance", source + ", multi-coverage-binning mocks", basis,
                    members)

    print("hcom2_closed", file=sys.stderr)
    members = []
    for r in jin:
        ani = number(r.get("ANI vs. NCBI genome", "") or "")
        acc = r.get("Closest NCBI genome", "")
        how = acc if acc.startswith(("GCA_", "GCF_")) and ani is not None and ani >= 95 else "none"
        members.append((r["Strain"], place(gtdb, r["Strain"], how, r.get("species") or None), None))
    write_truth(out, "hcom2_closed", "closed", "Jin et al. 2023, Supplementary Data 1 (the 125 genomes of hCom2's "
                "NinjaMap database; 117 strains grew in vitro)", "membership", members)

    print("hcom1a", file=sys.stderr)
    drop = [(n, place(gtdb, n, h)) for n, h in HCOM1A_DROPOUTS]
    for name, absent in (("hcom1a_full", ()), ("hcom1a_minus_ch", (0,)), ("hcom1a_minus_cs", (1,)),
                         ("hcom1a_minus_chcs", (0, 1))):
        write_truth(out, name, "dropout", "Wang et al. 2023: hCom1a (107 strains) and its arms without C. hylemonae "
                    "and/or C. scindens; only these two species are judged", "dropout",
                    [(n, s, 0 if i in absent else None) for i, (n, s) in enumerate(drop)])

    write_truth(out, "ecoli_kpneumoniae", "closed", "Afolayan et al. 2024: mixes of 2-3 strains of E. coli or "
                "K. pneumoniae (one species per mix)", "membership",
                [("Escherichia coli strains", place(gtdb, "E. coli", "species:Escherichia coli"), None),
                 ("Klebsiella pneumoniae strains", place(gtdb, "K. pneumoniae", "species:Klebsiella pneumoniae"),
                  None)])
    write_truth(out, "negative", "negative", "controls without a community (human DNA, extraction blanks, "
                "negative controls)", "none", [])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("samples", help="write samples.tsv from ENA's run metadata")
    s.add_argument("-o", "--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "samples.tsv"))
    t = sub.add_parser("truth", help="write the truth tables")
    t.add_argument("--gtdb", required=True, help="folder of the release's bac120/ar53_metadata_rNNN.tsv.gz")
    t.add_argument("-o", "--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "truth"))
    args = ap.parse_args()
    if args.command == "samples":
        write_samples(args.out)
    else:
        write_truths(os.path.expanduser(args.gtdb), args.out)


if __name__ == "__main__":
    main()
