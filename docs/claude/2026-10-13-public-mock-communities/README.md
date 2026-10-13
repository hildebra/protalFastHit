# Public mock-community metagenomes for an independent benchmark of protal

2026-10-13, branch `audit-fixes` at `9c27756`. A survey only: protal was not run, and nothing was
downloaded except archive metadata.

**Sources**

- The ENA portal API, queried on 2026-10-13.
- The papers, data releases and vendor manuals cited below.
- The local GTDB r226 metadata (WSL `~/GTDB/r226/{bac120,ar53}_metadata_r226.tsv.gz`), to check where a mock's genomes sit in GTDB.

**How this was checked**

- Four research agents searched separate areas: commercial standards, academic mocks, defined synthetic communities, and long reads. Their findings were then cross-checked.
- Every run recommended here is in [`runs.tsv`](runs.tsv). [`fetch_runs.sh`](fetch_runs.sh) re-read each run from ENA (study, instrument, layout, reads, bases), and all 163 were found with data.
- [`replicate_groups.tsv`](replicate_groups.tsv) lists 699 runs of the commercial and reference standards, grouped into replicate sets. An agent built it from ENA run reports; a sample of its rows was re-checked against ENA.
- Compositions, truth tables and chemistries come from the papers and manuals as read during the survey. Where a fact could not be confirmed, the text says "unverified".

## What a mock can and cannot test

**What it tests.** These are real reads with real errors, library and extraction biases, contaminants and depth. That makes them independent of protal's read simulators, which generate all of its training and test data.

**What it does not test.** Novelty: nearly every member is a type or reference strain in GTDB. For example, all 18 prokaryote genomes of Zymo D6331 (Zymo's 2023 deposits) are in r226 as non-representatives (table in [Truth](#truth-units-names-and-known-problems)). A mock therefore measures the case where every species is in the database. The held-out-clade evaluations of the builds remain the test for missing species and congeners.

**Truth units.** protal reports shares of cells (genomes), not of DNA. Zymo, NIBSC, NBRC and Meslier publish genome-copy truth; other mocks need converting (see [Truth](#truth-units-names-and-known-problems)).

**What protal cannot place.** Fungi (the Zymo yeasts) have no marker genes in GTDB. Either drop them and renormalise, or count them against the `?` share of `.profile` (see [Truth](#truth-units-names-and-known-problems)).

## Where to start

One set per question, all with public truth:

| Question | Dataset | Runs | Why |
|---|---|---|---|
| Species, pe, deep | Zymo D6331 DNA, Dohlman et al. 2026 | SRR38158100, SRR38158090 (306 / 276 Gb NovaSeq X) | 18 bacterial strains (14 species, five of them *E. coli*), 1 archaeon and 2 yeasts. Staggered to 1e-4 %. Deep enough to subsample to any depth. Two replicates (whether separate libraries is not stated). |
| Host background and false positives | same study | 1:1 to 1:10^8 microbe:human, 3 replicates each; human-only SRR38158104/103/102 | Detection limit against host reads, and a sample where every call is false. |
| Complex community, many platforms | Meslier et al. 2022, mocks 1-3 | HiSeq, DNBSEQ G400/T7, Ion S5, MinION R9, Sequel II (`meslier_*` in runs.tsv) | 64-87 strains over 29 phyla with archaea. Genome-copy truth with a GenBank accession per strain. The same DNA on 5-7 platforms. |
| Gut-like community, modern pe and ONT | Treichel et al. 2026 (Clavel lab) | 94 NovaSeq runs + ERR16002257 (R10.4.1 sup v5.2, 117 Gb) | 24-70 strains with *E. coli* and *P. vulgatus* strain sets. 13 compositions, 9 depths, 2 library facilities, with and without mouse DNA. |
| HiFi | Zymo D6331, PacBio | SRR13128014 (+013 low input, +012 amplified) | The standard set of long-read profiling papers; same truth as the deep Illumina. |
| HiFi, staggered | ATCC MSA-1003 | SRR11606871 (59 Gb; SRR11469785 is the same data) | 20 strains at 18 / 1.8 / 0.18 / 0.02 %, as DNA (no extraction bias). |
| ONT R10.4.1 | Zymo D6331 | SRR31222153 (128 Gb), MicroBench ERR14251410 (46.5 Gb, also as hac/fast) | Same truth again, on current chemistry. |
| Basecaller effect, ONT | Zymo D6322, one flow cell (PAS01578) | fast/hac/sup, simplex/duplex, Dorado models v4.2-v5.2 (group `D6322_PAS01578`) | Identical molecules; only the basecall differs. |
| Strains | D6331 (5 *E. coli*); Afolayan et al. 2-3-strain mixes | runs above; PRJEB82665 | Known ratios of conspecific strains (see [Strain-level tests](#strain-level-tests)). |
| Reproducibility | NIBSC/MHRA whole-cell reagent, 22 labs | PRJNA1337664 (270 runs) | 44 extraction methods × 5 replicates of one material. |
| Closed-set false positives at scale | hCom2-like 117-strain community in vitro; hCom1a with designed dropouts | PRJNA885585 (276 runs); PRJNA978640 (28 runs) | Every call outside the member list is false, and dropped strains must be absent. Abundances are measured, not designed (see [SynComs](#defined-synthetic-communities-syncoms)). |
| Members missing from GTDB | Plant SynComs of 177-988 strains | PRJNA1131994, PRJNA1449013 | Real reads of a defined community that r226 mostly lacks. |

## Mocks by read type

Run counts and bases are from ENA. For pe, "Gb" means both mates.

### Illumina and other paired-end short reads (`pe`)

| Mock | Members, truth | Runs | Notes |
|---|---|---|---|
| **Zymo D6331** gut standard (cells) | 21 strains: 14 bacterial species with 5 *E. coli* strains, *M. smithii* and 2 yeasts. Genome-copy % in the manual. | DNA: Dohlman 2026 PRJNA1454942 (above). Cells: HKU PRJNA804004 (12 NovaSeq runs, 5 extraction set-ups); Litichevskiy et al. 2025 (26 positive controls of one study, PRJNA1054518); Auckland SRR33710519/518 | Dohlman's runs include the human reads. |
| **Zymo D6300** (cells), **D6305/D6306** (DNA), **D6310** (log) | 8 bacteria + 2 yeasts. Genome copy 6.1-21.6 % (even); 94.8 % down to 1e-4 % (log). | Olm 2021 SRR12324251-253 (triplicate); Nicholls 2019 ERR2984773 (even), ERR2935805 (log); Yang 2020 PRJEB38036 (BGISEQ, 6 extraction protocols × 6); Auckland D6300 vs D6305 (SRR33710516/517 vs SRR33710466/467) | Zymo replaced five strains from lot ZRC190633 on. Species are unchanged; for strain-level scoring, match the lot to Zymo's reference bundle. |
| **Zymo D6322** HMW DNA | The D6300 bacteria minus *L. fermentum*, plus *S. cerevisiae* | ERR7255689 (MiSeq; same DNA as the ONT runs of Sereika 2022) | Sereika et al. say Zymo's D6322 references are not identical to the strains sequenced. |
| **Zymo D6332** oral | 12 bacteria (4 at 16 %, 4 at 8 %, 4 at 1 % gDNA, 4 *Streptococcus* species) | SRR36842452 (NovaSeq X) | Hard congeners. Per-strain truth not checked. |
| **ATCC MSA-1003 / 1002** (DNA) | 20 strains: staggered (MSA-1003) or 5 % each (MSA-1002) | MSA-1003: SRR8359173; NIBSC MiSeq × 5 (SRR11487948-952). MSA-1002: ERR15084348 (Visci et al. 2026, same DNA on ONT and HiFi) | Whether the percentages are DNA mass or genome copies is not stated on the product sheets. |
| **NIBSC/WHO Gut-Mix-RR, Gut-HiLo-RR** (DNA) | 20 DSM gut strains, 19 species (*B. longum* subsp. *infantis* + *longum*). Genome-copy % in Amos et al. 2020, Table 1 (Mix 2.7-9.0 %, HiLo 0.11-19.8 %). | Amos 2020 MiSeq × 5 each (SRR11487931-941); Anwar et al. 2025, 23 labs (PRJNA1238537), e.g. SRR34838295 (Mix, 227 Gb), SRR34838304/305 (Mix/HiLo, 57/63 Gb) | Two subspecies of one species. |
| **NIBSC WC-Gut-RR** (cells) | The same 20 strains, equal cells (flow cytometry) | Amos 2022 PRJNA855455 (128 runs); MHRA multi-lab PRJNA1337664 (270 runs, 22 labs) | Extraction bias by design. |
| **NBRC** mocks (NITE, Japan) | DNA mock: 20 strains, near-equal genome copies. Cell mock: 18 strains. Closed genomes. | Tourlousse et al. 2022 PRJNA747117: 20 DNA-mock + 36 cell-mock runs from 3 labs | Two *B. longum* strains. |
| **Mori et al. 2023** | 18 gut type strains (one per genus); equal cells, plus a DNA mix | PRJDB10817: 6 extraction kits + DNA mix, each extract on HiSeq 2500 and NovaSeq | No conspecifics. |
| **Treichel et al. 2026** | 70 gut isolates (HiBC/miBC): Mock-even-70, Mock-stag-24 (+10 variants), Mock-stag-70. DNA mass truth. | PRJEB83573: 94 NovaSeq PE150 runs | 4 *E. coli* and 4 *P. vulgatus* strains. Truth and run descriptions in the GitHub repo ClavelLab/Benchmarking-shallow-Metagenomics. |
| **Meslier et al. 2022** | 3 mocks of 71/87/64 strains (91 genomes, archaea included). Table S1 gives % genomes and GenBank accessions. | PRJEB52977: HiSeq 3000, DNBSEQ-G400, DNBSEQ-T7 (75 Gb) | Conspecific pairs: *Shewanella baltica*, *Methanococcus maripaludis*. References at forge.inrae.fr/metagenopolis/benchmark_mock. |
| **HC227** (Goussarov et al. 2024) | 227 strains, 197 species, even DNA mass | ERR5321934 (NovaSeq, 83 Gb) | About 30 conspecific strains. Isolate reads in PRJEB32402; no assemblies deposited. |
| **MBARC-26** (Singer 2016) | 23 bacteria + 3 haloarchaea, finished genomes; DNA molarity | SRR3656745 (52 Gb); PRJNA308761 (6 library variants, up to 144 Gb) | *Thermobacillus composti* contaminant. |
| **BMock12** (Sevim 2019) | 12 marine/lake bacteria; DNA molarity | SRR8073716 (64 Gb) | *Halomonas* pair at ~99 % ANI. *M. coxensis* is absent from the reads. |
| **Shakya 2013** | 48 bacteria + 16 archaea; qPCR truth | SRR606249 (11 Gb) | Classic, but at least 4 contaminant species are reported. Meslier's mocks reuse much of its DNA on current platforms. |
| **HMP mock B** (BEI HM-276D/277D) | 20 bacteria; truth in 16S operon copies | Jones 2015 PRJNA298489 (4 library kits × 2 runs) | Divide by operon counts to get genomes. |
| **GIS20** (Bertrand 2019) | 20 species, 0.1-30 % DNA | ERR3200809 (HiSeq 42 Gb) | Composition table not checked. |
| Element AVITI | "PositiveControl_ZymoCommunityStandard" | SRR30875866 | The only Element mock run found; catalogue number not given. |

### Single-end short reads (`se`)

- **Ion S5:** Meslier mocks 1-3, ERR9765777-779.
- **GAII 75 bp:** HMP mock B even and staggered, SRR172902/903.
- **Any paired run:** its first mates can be run as `se`. That gives the same truth and the same library as the `pe` run, so the read type is the only thing that changes.

### Ultima

There is no Ultima mock with a stated product. ENA lists 2,887 Ultima runs. The only mocks among them are three UG 100 runs:

- SRR31481588 (`mock_dna`), 45 Gb
- SRR31481586 (`mock_gut`), 55 Gb
- SRR31481591 (`mock_cellular`), 44 Gb

They are in PRJNA1189566, the stool metagenomes of Wu et al. 2025, *Nat Med*, doi:10.1038/s41591-025-03719-2, a study with Ultima Genomics among its authors. The paper's full text does not name the standards.

The labels suggest a Zymo DNA standard, D6331 and D6300. A subsample is enough to tell them apart: *M. smithii* and five *E. coli* strains mean D6331, and only the 8 D6300 species means D6300 or D6305. Once the product is known, the three runs give Ultima a truth. protal runs Ultima reads as `se`.

### PacBio HiFi (`pb`)

| Mock | Runs | Notes |
|---|---|---|
| Zymo D6331 | SRR13128014 (18 Gb, standard input), SRR13128013 (26 Gb, low input), SRR13128012 (21 Gb, ultra-low input with amplification) | Same DNA, three library protocols. |
| Zymo D6331, current instruments | `downloads.pacbcloud.com/public/revio/2025Q2/microbes/metagenomics/`: Revio and Vega BAMs of library bc2002 (SPK 3.0), plus an amplified Ampli-Fi library | Not in SRA; base counts unverified. |
| ATCC MSA-1003 | SRR11606871 (59 Gb; Hon et al. 2020); SRR9328980 / SRR9202034 (20 Gb each, one library) | SRR11469785 is the same data as SRR11606871 (identical read and base counts). |
| Zymo D6300 / D6305 | ERR12085099 / ERR12085100 (INRAE, 4-6 Gb) | Cells against DNA of one composition. |
| Meslier mock 1 | ERR9765783 (5.4 Gb) | Size-selected, which shifted abundances. |
| HC227 | ERR6178439-442 (Sequel subreads) | Run `ccs` first: ~5.8 Gb HiFi in all, ~7× per genome. |
| ATCC MSA-1002 | ERR15084350 (0.8 Gb) | Too shallow except as a cross-platform point. |

Older RS II CLR data (MBARC-26, BMock12, GIS20) are of little use for HiFi.

### Oxford Nanopore (`ont`)

**R10.4 and R10.4.1:**

| Mock | Runs | Chemistry, basecaller |
|---|---|---|
| Zymo D6331 | SRR31222153 (128 Gb; Chen et al. 2025, *Nat Commun*) | LSK114, PromethION; model unverified |
| Zymo D6331 | MicroBench (Aalborg, PRJEB85558): ERR14251410 / 409 / 408 = sup / hac / fast of the same 9.2 M reads; ERR15459768 = sup v5.2 | Dorado 0.8.2 v5.0.0; Dorado 1.0.0 v5.2.0 |
| Zymo D6331 | ERR15316007 (12 Gb; Benoit et al. 2026); SRR26387398 (3.6 Gb; Chen et al. 2024, *Genome Biol*) | sup v5; Guppy 6.5.7 sup |
| Zymo D6332 | ERR14251413 (MicroBench, 7.5 Gb), also as hac/fast | Dorado v5.0.0 |
| Zymo D6322 | One flow cell, PAS01578, deposited many times. PRJEB65462: fast/hac/sup v4.3.0, sup v4.2.0, duplex. PRJEB85558: v5.0.0 and v5.2.0, simplex and duplex. | R10.4.1 5 kHz, "high duplex" flow cell |
| Zymo D6322 | SRR24893244 / 245 / 246 = sup / hac / fast of the same 2.51 M reads (GridION 5 kHz) | Dorado 0.3.0 v4.2.0. Avoid SRR23563655: it holds the fast, hac and sup calls of one 4 kHz run in one file. |
| Zymo D6322 | ERR7287988 (R10.4, Kit 12), with ERR7255742 (R9.4.1) and ERR7255689 (MiSeq) of the same DNA | Sereika et al. 2022, Guppy 5.0.16 sup |
| Zymo DNA standard | ERR5396170 (R10.4 Q20 early access) | Bonito 0.3.5; the anchor of Portik et al. 2022 |
| Treichel Mock-stag-70 | ERR16002257 (117 Gb) | LSK114, Dorado 1.0.0 sup v5.2.0 |
| *E. coli* / *K. pneumoniae* 2-3-strain mixes | PRJEB82665 (Afolayan et al. 2024, bioRxiv): 20 mixes, each on R9.4.1 (ERR13964249-268) and R10.4.1 (ERR13964269-288); 10 *K. pneumoniae* ratio mixes (ERR15465276-285) | 0.25-2.6 Gb each |
| ATCC MSA-1002 | ERR15084349 (2.8 Gb) | LSK114 |

**R9.4.1** (old chemistry, but deep or complex):

- Nicholls 2019, Zymo D6300 even and D6310 log:
  - Even: GridION ERR3152364, PromethION ERR3152365 (152 Gb).
  - Log: GridION ERR3152366, PromethION ERR2906230 (148 Gb). The log's flip-flop re-call, ERR3152367, has 0 bases in ENA.
- HC227: ERR6109821 (75 Gb).
- HMP mock B: SRR11700265 / 264 (36 / 72 Gb).
- Meslier mocks 1-3: ERR9765780-782.
- Liu et al. 2022, D6331: SRR17913200 (Kit 9, 28 Gb) and SRR17913199 (Kit 12, 7 Gb). Their flow cells are reported inconsistently, and the deposited reads are already filtered.
- BMock12: SRR8073713 (size-selected) and SRR8351023 (not size-selected), a small library-prep comparison.

**Not in ENA:**

- Loman lab R10 and R10.3 D6300 runs, at lomanlab.github.io/mockcommunity/r10.html.
- ONT's own open-data bucket, which has no defined shotgun mock: `zymo_fecal_2025.05` is the D6323 stool reference without a composition, and `zymo_16s_2025.09` is amplicons.

## Strain-level tests

| Dataset | Strains | Reads |
|---|---|---|
| Zymo D6331 | 5 *E. coli* (JM109, B-3008, B-2207, B-766, B-1109) at 1.6-1.8 % genome copies each | pe, HiFi, ONT R10.4.1: all of the above |
| Afolayan et al. 2024 | 2-3 *E. coli* or *K. pneumoniae* strains at 1:1 to 10:25:65; sequence types the same or different | ONT R9 and R10.4.1 of each mix. Isolate assemblies: accessions not found |
| Treichel et al. 2026 | 4 *E. coli* (two at 99.97 % ANI) and 4 *Phocaeicola vulgatus* | pe; ONT for Mock-stag-70 |
| van Dijk et al. 2022 (StrainGE) | 4 *E. coli* at 80/15/4.9/0.1 % of the *E. coli*, ~1 % of a human-DNA sample | SRR13355226 (HiSeq X). The strain names differ between the paper and the BioProject |
| HC227 | ~30 strains sharing a species with another member | pe, ONT R9, PacBio subreads |
| Meslier 2022 | *S. baltica* OS185/OS223, *M. maripaludis* S2/C5 | all platforms |
| NIBSC, NBRC | two *B. longum* (subspecies or strains) | pe |
| BMock12 | *Halomonas* HL-4/HL-93 at ~99 % ANI | pe, old long reads |

## Technical replicates and duplicates

Replicates let protal be compared with itself where the truth is uncertain: cell mocks, lot changes, contaminants. Profiles of the same material should differ only as much as the material and process do.

[`replicate_groups.tsv`](replicate_groups.tsv) lists the commercial and reference standards' replicate sets (699 runs in 148 groups, with what differs in each). `runs.tsv` groups the curated runs in its `group` column.

| Level (what differs) | Datasets |
|---|---|
| **Nothing but the basecall**: identical molecules | D6322 flow cell PAS01578 (Dorado models v4.2-v5.2 × fast/hac/sup × simplex/duplex, PRJEB65462 + PRJEB85558). D6331 MicroBench sup/hac/fast. D6322 SRR24893244/245/246. Nicholls ERR2906227-230 (first calls) against ERR3152364-366 (flip-flop). |
| **The read type**: first mates as `se` against pairs as `pe` | Any paired run. |
| **The run**: one library sequenced twice | ATCC MSA-1003 HiFi SRR9328980 / SRR9202034. Meslier mock 3 ERR9765748/749 (one library over two runs; merge). Ion Proton runs of Meslier (9 re-runs per mock). Jones 2015 (each kit's library run twice). |
| **The platform**: one DNA on several | Meslier (HiSeq, G400, T7, Ion, MinION, Sequel II). Sereika D6322 (MiSeq, R9.4.1, R10.4). HC227 (NovaSeq, R9.4.1, Sequel). Mori 2023 (HiSeq and NovaSeq of each extract). BMock12, GIS20, MSA-1002 (Visci 2026). PacBio D6331 library bc2002 on Revio and Vega. |
| **The library**: one DNA, several libraries | NIBSC 2020 (5 libraries per reagent). NBRC (6-12 per lab). D6331 HiFi (standard, low, ultra-low input). Jones 2015 (Nextera, TruSeq, Kapa, Kapa PCR-free). MBARC-26 PRJNA308761 (6 variants). Treichel (two facilities, two kits). hCom2-like in vitro: one inoculum as three libraries (SRR21819256, SRR29248670, SRR29248672), and 90 cultures each with three libraries (a/b/c). Schirmer 2016 (4 preps × 3 platforms). |
| **The depth**: subsamples of one library | Treichel's 9 depth steps (subsampled with seqtk: not independent). Dohlman (~1 billion pairs per run: subsample at will). OMM12 PRJNA831628 (one DNA at six depths, 0.25-9.5 Gb). |
| **The biological replicate**: one designed community | OMM12 day 0 in 12 mice (PRJNA838765). hCom1a arms, 3 mice each. Kleiner 2017 (3 mixes per design). Saak 2023 (cultures in triplicate). |
| **The extraction**: one material, several extractions | Yang 2020 (6 protocols × 6). NIBSC 2022 (8 kits × 5). MHRA multi-lab (44 methods × 5). HKU D6331 (5 set-ups). Mori 2023 (6 kits). Auckland and INRAE (cells against DNA of one composition). |
| **The lab** | MHRA whole-cell (22 labs). NIBSC DNA reagents (23 labs). NBRC (3 labs). |
| **The batch** | Litichevskiy 2025: 26 D6331 positive controls across extraction and library batches. |
| **The host background** | Dohlman D6331 at 1:1 to 1:10^8 in human DNA, each in triplicate. Treichel with and without germ-free mouse DNA. |

**The same data under two accessions** (count them once):

- SRR11469785 = SRR11606871 (ATCC MSA-1003 HiFi).
- SRR12829162 and SRR12829163 (Zymo DNA standard) have identical read and base counts.
- PAS01578 is in both PRJEB65462 and PRJEB85558.
- The MicroBench D6331 barcode is also in PRJEB74343.
- PRJNA564477 holds subsamples of the Nicholls runs.
- SRR25211204, "ShortPE_mock" labelled as Sequel II, is short reads made from the D6331 HiFi data, not a sequencing run.

## Truth: units, names and known problems

**Units.** protal's abundances are shares of genomes (cells). Treat each mock's truth as follows:

- **Already genome copies:** the Zymo "Genome Copy" column (Zymo recommends it for shotgun data; the "Cell Number" column is nearly the same), NIBSC (Amos 2020, Table 1), the NBRC DNA mock, and Meslier Table S1 ("% genomes").
- **DNA mass:** divide by genome size. This covers Treichel, HC227, GIS20, possibly ATCC (the basis is not stated), and the Zymo "Genomic DNA" column.
- **DNA molarity** (BMock12, MBARC-26): close to genome copies.
- **16S operon copies** (HMP mock B): divide by each strain's operon count.
- **Cell mocks:** the designed cell counts include the extraction's bias. Score them against the matching DNA mock, or against replicates.

**Yeasts.** Two Zymo standards contain yeasts:

| Standard | Yeast share of genome copies | protal's `?` share should be |
|---|---|---|
| D6300 | 0.94 % | close to 0.94 % |
| D6331 | 0.63 % | close to 0.63 % |

Drop the yeasts from the species truth. Their shares also give a small check of protal's `?` line ([unknown share](../../running.md)): the "should be" column assumes nothing else in the sample is unplaced.

**GTDB r226 names.** These were checked in the local metadata for D6331. All 18 prokaryote genomes Zymo deposited in 2023 (GCF_028743255.1 ... GCF_028743775.1, GCF_030770375.1) are in r226 as non-representatives. Four vendor names map to other GTDB r226 species:

| Vendor name | GTDB r226 species |
|---|---|
| *F. prausnitzii* | *Faecalibacterium prausnitzii_E* |
| *F. nucleatum* 2/1/50A | *Fusobacterium animalis* |
| *M. smithii* | *Methanocatella smithii* |
| *C. perfringens* | *Sarcina perfringens* |

Vendor names also differ elsewhere:

- *P. corporis* is listed as "*P. copri*" on Zymo's website; the manual and NCBI say *P. corporis*.
- ATCC lists *B. pacificus*, *P. paraeruginosa* and *P. vulgatus*.

Map every mock through genome accessions, not names. NCBI taxids alone are not enough either: r226 gives most species placeholder taxids ([NCBI taxids](../2026-10-11-ncbi-taxids.md)).

**Known problems:**

- **Contaminants:** Shakya (≥4 species, found by Mash Screen), MBARC-26 (*T. composti*), and ATCC (*Lysobacter* in the *C. acnes* component).
- **Dropouts:** BMock12 *M. coxensis*.
- **Strain labels:** the D6300 *S. enterica* is "likely Choleraesuis", not LT2 (Nicholls 2019).
- **Detection limit:** *C. perfringens* in D6331 is at 9e-5 % of genome copies.
- **Marker-gene counts:** *Cereibacter sphaeroides* (ATCC) has two chromosomes and duplicated *rpsG*/*rplK*.
- **Distorted libraries:** amplified (SRR13128012, SRR23563655) or size-selected (Meslier HiFi, BMock12) libraries skew abundances.

## Defined synthetic communities (SynComs)

A SynCom is a defined set of sequenced isolates, grown together in vitro, in a gnotobiotic animal or on a plant. Its abundances are measured, not designed, and its truth is mostly the authors' own mapping. What it does offer is a closed set: any call outside the member list is a false positive. That holds even with host DNA, and for communities larger than any mock. All runs below are paired-end short reads; no community of this kind has public long reads.

**Human gut, large:**

| Community | Data | Notes |
|---|---|---|
| hCom2-like, 117 strains, in vitro (Jin et al. 2023, *Nat Commun*) | PRJNA885585: 276 NovaSeq runs, 985 Gb, ~3.6 Gb each. 5 conditions × 6 passages × 3 cultures × 3 libraries, plus inocula. | No host DNA. Several pairs at >99 % ANI (*B. dorei* ×3, *Veillonella*, *Subdoligranulum*). One inoculum library set three ways: SRR21819256, KAPA SRR29248670, Flex SRR29248672. 90 technical triplets of one culture each (e.g. SRR21819180-182). |
| hCom1a, 107 strains, with designed dropouts (Wang et al. 2023, *Cell*) | PRJNA978640: 28 runs, 345 Gb. Each arm's inoculum (SRR24799814, 817, 818, 819) and 3 mice per arm. | Arms: full, minus *C. hylemonae*, minus *C. scindens*, minus both. Strains deliberately absent from otherwise identical communities, so a direct FP/FN test. Every sample also holds a *Paenibacillus barengoltzii* spike-in (10^5 cells). |
| hCom2 in gnotobiotic mice, 5 generations (Jin et al. 2025, *Microbiome*) | PRJNA1119053: 69 runs, 656 Gb | Host-read filtering kept a median 91.6 % of reads. The original Cheng et al. 2022 hCom1/hCom2 metagenomes were never deposited. |
| Genomes for the hCom sets | Zenodo 10.5281/zenodo.7872423 (`hCom2.tar.gz`, the NinjaMap reference); isolate reads PRJNA746600, PRJNA885826, PRJNA1119029 | The survey found all 46 PRJNA746600 assemblies in GTDB. |
| 18 human gut strains in mice (Zimmermann-Kogadeeva et al. 2022, bioRxiv) | PRJEB33195: 10 DNA runs (e.g. ERR10107056) | Culture-collection strains. Member list in the ENA study description. |
| DefCom1/2, 24 and 46 species from 4 donors (Fishbein et al. 2025, *Cell Host Microbe*) | PRJNA1067813: only the ~90 `DefCom*` runs are the closed set | The other ~350 runs are humanised or donor samples. |
| 32 human gut species in vitro, with drugs (Wuyts et al. 2023, *Mol Syst Biol*) | PRJEB46619: 68 WGS runs, ~1.1 Gb each | |

**Mouse:**

| Community | Data | Notes |
|---|---|---|
| OMM12 (Oligo-MM12), antibiotic series | PRJNA838765: 112 runs, 3.8-10.3 Gb. 12 mice at day 0 (e.g. SRR25305803, 828), one 16-point series per arm. | All 12 genomes are in r226 (checked locally). *L. reuteri* I49 is `s__Limosilactobacillus sp011958095` there. Link to Münch et al. 2023 inferred from the BioProject. Antibiotics push members below detection. |
| OMM12, other sets | PRJNA831628: one DNA at six depths, 0.25-9.5 Gb, plus Hi-C. PRJNA1196393 (DSS colitis). PRJNA1528472 (with *C. difficile*, 0.4-0.6 Gb). | |
| Altered Schaedler Flora, 8 strains (Proctor et al. 2022, *mSystems*) | PRJNA838940: SRR19241958/959, one BioSample | Truth is *groEL* qPCR, which disagrees with shotgun for two members. Three members are GTDB placeholders. |

The survey found no public shotgun data for GM15, PedsCom, MDb-MM, SIHUMI/SIHUMIx in vivo, OMM19.1, Com20/21 or the Venturelli-lab communities: they are 16S, qPCR or RNA only.

**Plant:**

- **Selten et al. 2026, *Nat Microbiol*:**
  - PRJNA1131994: 442 runs, ~0.4 Gb each. Arabidopsis 447-strain, barley 364, Lotus 177, and all three together (988 genomes). Includes the inocula (e.g. SRR29710026, SRR29710013) and plants without a SynCom.
  - PRJNA1449013: 88 runs, ~2 Gb each. Seven 200-strain Arabidopsis SynComs; their inocula are SRR37935102-108.
  - Most members are **not in GTDB r226**: three spot-checked AtSC genomes are absent from the local metadata.
  - That turns these into a real-reads novelty test: a defined community whose members are mostly missing from the database, with the plant's own DNA alongside.
- **PRJNA1444504:** a Lotus SynCom with family dropouts (78 runs, 2.75 Gb).

The best-known plant SynComs are amplicon-only: Bai 2015, Durán 2018, Finkel 2020 (185 members), Wippel 2021 and the Vorholt dropouts.

**Environmental and food:**

| Community | Data | Notes |
|---|---|---|
| Schirmer et al. 2016: 59 genomes (49 bacteria, 10 archaea) from Shakya DNA, balanced and unbalanced | PRJEB12287: 47 runs on GAII, HiSeq 2500 and MiSeq, 4 library preps | Which run is which design is not in ENA. |
| Kleiner et al. 2017: 32-member cell mock (bacteria, an archaeon, an alga, phages) | PRJEB19901: 9 NextSeq 2×75 runs. Equal cells, equal protein and uneven designs × 3 mixes (ERR1877474-482). | Truth from cell counts and protein. |
| Saak et al. 2023: 16-member cheese rind (14 bacteria, 2 fungi) | PRJNA852571: 23 NovaSeq runs, 8-16 Gb. Input SRR19846572; time series in triplicate; Actinobacteria-dropout and *G. geotrichum*-dropout arms. | Four *Psychrobacter* and two *Staphylococcus*. |
| Ziegler et al. 2025 (bioRxiv): *Prochlorococcus* MED4 + 5 heterotrophs | PRJNA1333906: 21 NextSeq runs, 3 replicates × 7 days | Two *Thalassospira* at 99.99 % ANI. *Thermus* spike-in. |

Excluded: SOLiD data (Dugat-Bony 2015), communities only partly defined (Hot Lake UCC), and the many marine, bioreactor and food SynComs followed only by amplicons or RNA.

## Gaps

- **Ultima:** no mock with a stated product (see [Ultima](#ultima)).
- **Element:** one run.
- **ONT R10.4.1:** none for the even Zymo D6300/D6305 standards. The Q20 early-access run and the Loman R10/R10.3 releases are the closest.
- **NIST Mosaic Standards Challenge:** Mix A and Mix B, with ddPCR genome-copy truth, are on data.nist.gov (doi:10.18434/mds2-2830), not in INSDC. Which kit sample is which mix is unverified.
- **High-complexity mocks:** apart from HC227 (R9, subreads), none exists on current long-read chemistry. Treichel's Mock-stag-70 (70 strains) is the largest on R10.4.1.

## Files

- [`runs_in.tsv`](runs_in.tsv): the curated runs (dataset, mock, read type, group, run, note).
- [`fetch_runs.sh`](fetch_runs.sh): reads each run of `runs_in.tsv` from the ENA portal API and writes `runs.tsv`. Its `check` column flags runs missing in ENA or without bases. Run it with `sh fetch_runs.sh` (curl, awk).
- [`runs.tsv`](runs.tsv): its output of 2026-10-13. All 163 runs are present, with bases.
- [`replicate_groups.tsv`](replicate_groups.tsv): 699 runs of commercial and reference standards. Columns: standard, project, run, instrument, layout, spots, Gb, group, replicate, what differs, note. Built from ENA run reports by a survey agent; spot-checked.

## Follow-up: scripts to download and score the starting sets (same day)

`scripts/mock_benchmark/` downloads the sets of [Where to start](#where-to-start) and scores protal on them. Its usage
is in [development.md](../../development.md#benchmark-on-public-mock-communities). Nothing was run with protal yet.

| File | What it holds or does |
|---|---|
| `samples.tsv` | 1,330 samples (1,355 runs) in 13 sets, from ENA's metadata by `prepare_tables.py samples` |
| `truth/` | 29 truth tables, by `prepare_tables.py truth --gtdb ~/GTDB/r226` |
| `download_mocks.py` | ENA downloads, cut while streaming |
| `score_mocks.py` | `run` writes protal maps and runs protal; `score` scores the profiles |
| `test_mock_benchmark.py` | 10 tests, 5 s |

**What is in `samples.tsv`.**
- **Replicate groups:** basecall series, platforms of one DNA, labs × methods, and technical triplets.
- **Treated as one sample:** a library's re-runs, e.g. Meslier's nine Ion Proton runs per mock.
- **Left out:**
  - raw-signal (pod5) runs;
  - modified-base BAMs that repeat a sup call;
  - SRR23563655, which holds three basecalls of one run in one file;
  - SRR11469785, which is the same data as SRR11606871.
- **Two runs without an ENA FASTQ** (ERR15995300, ERR16002257) are fetched from their submitted FASTQ.

**Truth tables, by member list:**

| Mocks | Source of the member list |
|---|---|
| Zymo D6331, D6322 | Zymo's NCBI genomes |
| ATCC MSA-1003 | ATCC strain numbers in the genomes' strain identifiers, or named strains (MG1655, UA159, FPR3757, 2603V/R, 26695, MC58, R1) |
| NIBSC DNA reagents and WC-Gut-RR | Amos et al. 2020, Table 1 |
| Meslier mocks 1-3 | Table S1 (% genomes, GenBank accessions) |
| Treichel's 13 compositions | Supplementary Tables 1-2 |
| hCom2 closed set | Jin et al. 2023, Supplementary Data 1 |
| hCom1a dropout arms | Wang et al. 2023 |
| Afolayan strain mixes (closed set) | Afolayan et al. 2024 |
| Controls (negative) | none expected |

**Truth tables, by how the shares are set:**
- Treichel's DNA masses and ATCC's percentages (taken as DNA mass) are turned into genome copies by genome size.

**Names that a name-based truth would have got wrong.** Every assignment is printed with its evidence. These differ from the vendors' names:

| Strain | GTDB r226 species | Note |
|---|---|---|
| Zymo *L. monocytogenes* | `s__Listeria monocytogenes_B` | |
| Zymo *B. subtilis* B-354 | `s__Bacillus spizizenii` | NCBI says the same |
| NIBSC *F. prausnitzii* DSM 17677 (A2-165) | `s__Faecalibacterium duncaniae` | |
| *B. longum* subsp. *infantis* | `s__Bifidobacterium infantis` | a separate species from *B. longum* |
| ATCC 17982, sold as *Schaalia odontolytica* | `s__Pauljensenia dentiphila` | its genome GCF_000154225.1 is NCBI's *S. dentiphila* ATCC 17982 |
| *C. hylemonae* | `s__Extibacter hylemonae` | |
| *C. scindens* | `s__Clostridium_AP scindens` | |
| *P. fluorescens* SBW25 | `s__Pseudomonas_E petroselini` | its Meslier genome is not in r226; found through its other assemblies |

**Unplaced members:**
- Treichel's *Blautia fusiformis* CLA-AA-11 has no assembly. It is placed as `s__Blautia_A fusiformis`; the authors' GTDB-Tk name, *B. massiliensis*, is not in r226.
- One hCom2 member stays unplaced: VE202-27 (*Eisenbergiella*), whose genome r226 lacks.

**Checked so far:**
- **Unit tests:** pass.
- **Real downloads of a few thousand reads each** (WSL, 2026-10-13), each cut as intended with mates in step, and from it the maps of `score_mocks.py run --map-only`:
  - Dohlman pe, cut from 306 Gb;
  - Meslier HiSeq;
  - Treichel pe and ONT, from submitted FASTQ;
  - PacBio HiFi;
  - MicroBench ONT.

**What the default limits keep** (`download_mocks.py --dry-run`):

| Set | Samples | GB at ENA | GB kept |
|---|---|---|---|
| d6331_dna | 2 | 236 | 6 |
| d6331_host | 33 | 931 | 42 |
| meslier | 17 | 187 | 46 |
| treichel | 94 | 509 | 107 |
| d6331_hifi | 3 | 46 | 45 |
| atcc1003_hifi | 1 | 45 | 18 |
| d6331_ont | 5 | 291 | 107 |
| d6322_basecalls | 21 | 563 | 388 |
| strain_mixes | 50 | 53 | 53 |
| nibsc_multilab | 270 | 196 | 196 |
| hcom2_invitro | 276 | 359 | 360 |
| hcom1a_dropouts | 28 | 135 | 74 |
| plant_syncoms | 530 | 210 | 205 |
| all | 1,330 | | ~1,650 |

`--max-bases 5e9` brings the D6322 basecall series to about 100 GB.

## Next steps

Step 1 is done; steps 2 and 3 have scripts (see the follow-up above) and wait for the cluster.

1. **Truth tables:** built for every set but the plant SynComs (`scripts/mock_benchmark/truth/`).
2. **Downloads:** a first set on the cluster: Dohlman D6331, Zymo D6331 HiFi and R10.4.1, Meslier mocks, Treichel, ATCC MSA-1003 HiFi (`download_mocks.py --sets d6331_dna,d6331_hifi,d6331_ont,meslier,treichel,atcc1003_hifi`, ~330 GB at the default limits).
3. **Scoring:** `score_mocks.py run` per set on the cluster, then `score_mocks.py score`: species precision, recall and abundance error per read type against the truth; replicate sets by Bray-Curtis and Jaccard, including the D6322 basecall series.
