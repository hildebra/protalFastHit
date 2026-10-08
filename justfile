set shell := ["bash", "-cu"]

# ---- strain test (M1-M5) parameters -- override on the command line, e.g.
#   just strain-test preset=sensitive
# ----------------------------------------------------------------------------
build_dir   := "build"
protal      := build_dir / "protal"
strain_db   := env_var_or_default("PROTAL_DB_PATH", "/home/fritscher/data/db/tool/protal/protal-db-r226-0.6.0a")
strain_input:= "/home/fritscher/non-git/simulate_metagenomes_test/test2"   # simulated dataset (input)
strain_out  := justfile_directory() / "strain_test_out"
strain_variant := "test1"                 # output subfolder: test1 (full filtering) / test2 (raw)
strain_run  := strain_out / strain_variant
strain_threads := "8"
preset      := "default"
# Extra protal flags. test2 (raw) disables ALL protal gene/sample filtering so the
# MSA keeps every observed gene and detected sample (only M1/M2 SNP filtering stays),
# leaving gene/sample filtering entirely to qcmsa and letting users re-filter.
strain_protal_filter_args := ""

# Delete all build trees
clear:
    rm -rf cmake-build-* {{build_dir}}

# Full M1-M5 strain test: run protal (default settings, with the qcmsa post-filter)
# on the pre-existing dataset alignments, then build the HTML QC report.
strain-test: strain-protal strain-dbcounts strain-report
    @echo "Report: {{strain_run}}/report/report.html"

# "Raw" variant -> strain_test_out/test2: protal filters ONLY SNPs (M1/M2) and does
# not even apply the per-sequence hcov floor. Gene/sample coverage filtering always
# lives in qcmsa. The unfiltered <species>.raw.msa.fna can be re-filtered with other
# thresholds.
strain-test-raw:
    just strain_variant=test2 \
         strain_protal_filter_args="--msa_min_hcov 0" \
         strain-test

# Count marker genes per species in the DB genome (the true gene denominator,
# revealing how many markers were lost to abundance before M3). Cached as a TSV.
strain-dbcounts:
    PROTAL_DB_PATH="{{strain_db}}" python3 scripts/strain_test/db_gene_counts.py \
        --db "{{strain_db}}" \
        --strains {{strain_run}}/strains \
        --out {{strain_run}}/db_gene_counts.tsv

# Simulated input in the layout strain-protal reads: simulate_metagenomes writes the reads to
# strain_input/output/reads, protal aligns them into strain_input/protal/alignments, and the map
# is kept as strain_input/protal_map.tsv. strain_genomes is a simulate_metagenomes --genome_table.
strain_genomes := ""
strain_sim_args := "-n 8 --total_read_pairs 40000 --species_per_sample 3 --seed 31"
strain-input: baseline simulate
    test -n "{{strain_genomes}}" || { echo "strain-input: set strain_genomes=<simulate_metagenomes --genome_table>" >&2; exit 2; }
    {{build_dir}}/simulate_metagenomes --genome_table "{{strain_genomes}}" -o "{{strain_input}}/output" \
        --protal_metafile "{{strain_input}}/protal" -t {{strain_threads}} {{strain_sim_args}}
    PROTAL_DB_PATH="{{strain_db}}" {{protal}} --map "{{strain_input}}/output/protal.meta" -t {{strain_threads}} --no_profile
    cp "{{strain_input}}/output/protal.meta" "{{strain_input}}/protal_map.tsv"

# Build a map that reuses the existing alignments + reads and writes strain
# results into the repo-local run dir, then run protal with the qcmsa post-filter (M5).
# The input map's directory variables are dropped: an input #SAM_OUTPUT_DIR (as
# protal_map_utils generate writes) would otherwise point protal into the run dir,
# where it finds no SAM files and aligns every sample again.
strain-protal:
    mkdir -p {{strain_run}}/strains {{strain_run}}/misc {{strain_run}}/profiles
    # Assemble the map: base OUTPUT_DIR + reuse existing SAMs/reads, local strain output.
    {{ '{' }} \
      printf '#OUTPUT_DIR\t%s\n'          "{{strain_run}}"; \
      printf '#INPUT_DIR\t%s\n'           "{{strain_input}}/output/reads"; \
      printf '#SAM_OUTPUT_DIR\t%s\n'      "{{strain_input}}/protal/alignments"; \
      printf '#PROFILE_OUTPUT_DIR\t%s\n'  "{{strain_run}}/profiles"; \
      printf '#STRAIN_OUTPUT_DIR\t%s\n'   "{{strain_run}}/strains"; \
      printf '#MISC_OUTPUT_DIR\t%s\n'     "{{strain_run}}/misc"; \
      grep -vE '^#(INPUT|OUTPUT|SAM_OUTPUT|PROFILE_OUTPUT|STRAIN_OUTPUT|MISC_OUTPUT)_DIR[[:space:]]' "{{strain_input}}/protal_map.tsv"; \
    {{ '}' }} > {{strain_run}}/strain_test_map.tsv
    # Write the full log to a file (avoids a pipeline whose trailing grep can fail
    # the recipe under `set -o pipefail`, which lmod's BASH_ENV enables), show a
    # progress-bar-stripped tail, then exit with protal's own status so a failed run
    # stops `just strain-test` instead of producing a report from stale outputs.
    PROTAL_DB_PATH="{{strain_db}}" {{protal}} \
        --map {{strain_run}}/strain_test_map.tsv \
        -t {{strain_threads}} \
        --qcmsa_args "--preset {{preset}}" {{strain_protal_filter_args}} \
        > {{strain_run}}/protal_run.log 2>&1; rc=$?; \
    tr '\r' '\n' < {{strain_run}}/protal_run.log | grep -vE '^\[=*>* *\] *[0-9]+ %' | tail -40; \
    exit $rc

# Build a per-species ML tree from each strain MSA of the last protal run with IQ-TREE.
# strain_tree_input selects the MSA: "filtered" = qcmsa output <sp>.msa.fna (default);
# "raw" = protal native <sp>.raw.msa.fna. The species come from strains/species.tsv, which
# protal writes each run, so MSAs an earlier run left behind are not used. Uses the `iqtree`
# conda env by default; override the launcher with strain_iqtree=... . Skips gracefully if
# IQ-TREE is unavailable. strain_iqtree_model is the substitution model: after qcmsa
# --discard-constant it needs +ASC (e.g. "GTR+G+ASC").
strain_tree_input := "filtered"
strain_iqtree := ""
strain_iqtree_model := "GTR+G"
strain_iqtree_seed := "1"
strain-trees:
    #!/usr/bin/env bash
    set -uo pipefail
    # Resolve an IQ-TREE launcher: explicit override > PATH > `iqtree` conda env.
    iq="{{strain_iqtree}}"
    if [ -z "$iq" ]; then
        iq=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
    fi
    if [ -z "$iq" ] && command -v conda >/dev/null 2>&1; then
        for b in iqtree3 iqtree2 iqtree; do
            if conda run -n iqtree "$b" --version >/dev/null 2>&1; then
                iq="conda run -n iqtree $b"; break
            fi
        done
    fi
    if [ -z "$iq" ]; then
        echo "[trees] IQ-TREE not found (PATH or 'iqtree' conda env); skipping."
        exit 0
    fi
    echo "[trees] using: $iq"
    mkdir -p "{{strain_run}}/trees"
    list="{{strain_run}}/strains/species.tsv"
    if [ ! -f "$list" ]; then
        echo "[trees] $list not found: run protal (a version that writes it) first." >&2
        exit 1
    fi
    case "{{strain_tree_input}}" in raw) ext=".raw.msa.fna";; *) ext=".msa.fna";; esac
    built=0
    while IFS=$'\t' read -r sp taxid samples raw filtered; do
        [ "$sp" = species ] && continue
        msa="{{strain_run}}/strains/$sp$ext"
        if [ ! -f "$msa" ]; then
            echo "[trees] $sp: no $sp$ext (see the qcmsa output in the protal log) - skipping"
            continue
        fi
        nseq=$(grep -c '^>' "$msa")
        if [ "$nseq" -lt 4 ]; then
            echo "[trees] $sp: only $nseq sequences (<4, the reference row included) - skipping"
            continue
        fi
        echo "[trees] $sp ($nseq taxa) -> {{strain_run}}/trees/$sp.treefile"
        # Unpartitioned ML tree with 1000 ultrafast bootstraps (partitioned trees came out the
        # same in the 2026-09-29 strain audit). IQ-TREE reads protal's IUPAC codes as ambiguities.
        $iq -s "$msa" -m {{strain_iqtree_model}} -B 1000 -T {{strain_threads}} --seed {{strain_iqtree_seed}} \
            --seqtype DNA --prefix "{{strain_run}}/trees/$sp" -redo \
            > "{{strain_run}}/trees/$sp.iqtree.stdout.log" 2>&1 \
          && built=$((built+1)) \
          || echo "[trees] $sp: IQ-TREE failed (see {{strain_run}}/trees/$sp.iqtree.stdout.log)"
    done < "$list"
    echo "[trees] built $built tree(s) in {{strain_run}}/trees (input={{strain_tree_input}})"

# Re-filter a raw run's MSAs with qcmsa: coverage gating (the share of each gene a row
# writes, the gene's mean depth) PLUS the usual MRate2 + site cleanup. Lets
# you re-filter strain_test_out/test2 (the raw run) with any thresholds without
# re-running protal. Outputs into <run>/refiltered/.
refilter_hcov        := "0.3"
refilter_depth       := "0"    # >0: minimum mean depth of a gene (all its positions)
refilter_min_samples := "1"    # a gene needs more than this many samples passing coverage
refilter_reapply_hcov := "1000" # as protal passes its --msa_min_hcov
refilter_sample_abs  := "0"    # >0: remove a sample multi-allelic in >= N genes (catches conspecific/mixed strains)
refilter_gene_abs    := "0"    # >0: remove a gene multi-allelic in >= N samples
strain-refilter:
    #!/usr/bin/env bash
    set -uo pipefail
    mkdir -p "{{strain_run}}/refiltered"
    shopt -s nullglob
    n=0
    for msa in "{{strain_run}}/strains/"*.raw.msa.fna; do
        sp=$(basename "$msa" .raw.msa.fna)
        part="{{strain_run}}/strains/$sp.raw.partition.txt"
        meta="{{strain_run}}/strains/$sp.meta.tsv"
        [ -f "$part" ] && [ -f "$meta" ] || continue
        python3 scripts/qcmsa.py "$msa" "$part" "$meta" \
            --prefix "{{strain_run}}/refiltered/$sp" --preset {{preset}} \
            --gene-min-hcov {{refilter_hcov}} \
            --gene-min-mean-depth {{refilter_depth}} \
            --gene-min-samples {{refilter_min_samples}} \
            --reapply-hcov {{refilter_reapply_hcov}} \
            --sample-abs-min-bad {{refilter_sample_abs}} \
            --gene-abs-min-bad {{refilter_gene_abs}} \
            > "{{strain_run}}/refiltered/$sp.qcmsa.log" 2>&1 \
          && { echo "[refilter] $sp"; n=$((n+1)); } \
          || echo "[refilter] $sp: failed (see {{strain_run}}/refiltered/$sp.qcmsa.log)"
    done
    echo "[refilter] re-filtered $n species (hcov>={{refilter_hcov}}, depth>={{refilter_depth}}, >{{refilter_min_samples}} samples) -> {{strain_run}}/refiltered"

# (Re)build the self-contained HTML QC report from an existing strain run.
strain-report:
    python3 scripts/strain_test/strain_report.py \
        --strains {{strain_run}}/strains \
        --out {{strain_run}}/report
    @echo "Open {{strain_run}}/report/report.html"

# Remove the strain test outputs.
strain-clean:
    rm -rf {{strain_run}}

# ---- mini DB: sparse synthetic GTDB release (3 species) -> protal DB -----------
# Writes {{mini_db_dir}}/gtdb_r226 (GTDB-layout release) and
# {{mini_db_dir}}/protal_db (use with --db), holding the single-file database.protal
# (~1 MB; raw, index.prx would be ~3 GB regardless of size). PROTAL_BUILD_ARGS=--no_bundle
# keeps separate compressed files, --no_compress builds it raw.
mini_db_dir := "data/mini_db"
mini-db: baseline
    PROTAL={{protal}} bash scripts/mini_db/build_mini_db.sh {{mini_db_dir}}

# Checks of the scripts: the mini DB generator and converter, the GTDB downloads, the collector and scenarios, the
# gene neighbours, the GTDB build script (scripts/mini_db/test_*.py), the in-silico strains, trace_relatives.py,
# error_reads.py, the profile scripts and the strain test's reports. Needs numpy; git, cmake and the zstd CLI for a
# few tests. The GTDB build end to end (test_gtdb_pipeline.py, a few minutes) and protal's builds of gene
# subsets run with PROTAL and SIMULATE set (e.g. PROTAL=$PWD/build/protal SIMULATE=$PWD/build/simulate_metagenomes
# just mini-db-test) and a Python with scikit-learn (PROTAL_TRAIN_PYTHON, default python3). A test whose prerequisite
# is missing is skipped; PROTAL_TESTS_REQUIRED=1 makes it fail instead (scripts/prerequisites.py).
mini-db-test:
    python3 -m unittest scripts/mini_db/test_*.py
    python3 -m unittest scripts/test_insilico_strains.py scripts/test_trace_relatives.py scripts/test_error_reads.py \
        scripts/test_composition_accuracy.py scripts/test_profile_scripts.py scripts/test_strain_scripts.py \
        scripts/test_foreign_rates.py

# Checks that the presence model's PMML export scores as scikit-learn does, of the trainer and of the rules it shares
# with protal (tests/data/golden_model_rules.tsv); needs numpy, pandas, joblib and scikit-learn (skipped without them,
# failed with PROTAL_TESTS_REQUIRED=1), no Java. A few minutes, on one thread.
model-test:
    python3 -m unittest -v scripts/test_model_pmml.py

# C++ unit tests (GoogleTest; needs libgtest-dev).
test:
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON
    cmake --build {{build_dir}} --target protal_tests -- -j$(nproc)
    ctest --test-dir {{build_dir}} --output-on-failure

# End-to-end tests of protal and simulate_metagenomes on a freshly built mini DB (tests/e2e/test_protal_e2e.py; about
# 1.5 minutes on 4 cores). Needs Linux, numpy and the zstd CLI. A test whose prerequisite is missing is skipped;
# PROTAL_TESTS_REQUIRED=1 makes it fail instead.
e2e: mini-db simulate
    PROTAL_TEST_DB={{mini_db_dir}}/protal_db PROTAL={{protal}} SIMULATE={{build_dir}}/simulate_metagenomes \
        python3 -m unittest -v tests/e2e/test_protal_e2e.py

# Profiling accuracy on the mini DB: build it (reused while unchanged), simulate reads
# from a known mock community, profile them and check against the truth (exit 1 if a
# check fails; about 20 s with the build). See examples/mini_db/README.md.
example: baseline
    PROTAL={{protal}} bash examples/mini_db/run.sh

# The ISA flags are per-target (isa_baseline in CMakeLists.txt), so the
# dynamic and static binaries both come out of one tree -- no need for
# separate cmake-build-* dirs.
# Configure the build tree
configure:
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE=Release

# protal: one binary for every x86-64 CPU (x86-64 baseline; the hot functions also
# for x86-64-v3, chosen at run time)
baseline: configure
    cmake --build {{build_dir}} --target protal -- -j$(nproc)

# simulate_metagenomes build
simulate: configure
    cmake --build {{build_dir}} --target simulate_metagenomes -- -j$(nproc)

# Static build: fully static protal + simulate_metagenomes (isa_baseline =>
# -march=x86-64; the x86-64-v3 copies of the hot functions are chosen at run time as in `protal`).
# zstd and ISA-L (and nasm, without one on the PATH) are downloaded and built in the build tree for them
# (lib/static-deps.cmake), so no libzstd.a or libisal.a is needed; `just static_fetch_deps=OFF static` links the system's.
static_fetch_deps := "ON"
static:
    cmake -S . -B {{build_dir}} -DCMAKE_BUILD_TYPE=Release -DPROTAL_STATIC_FETCH_DEPS={{static_fetch_deps}}
    cmake --build {{build_dir}} --target protal_static simulate_metagenomes_static -- -j$(nproc)

# Build all binaries that `just install` ships
build-all: configure
    cmake --build {{build_dir}} --target protal simulate_metagenomes -- -j$(nproc)

# Always rebuilds first so the installed binaries match the working tree (an
# out-of-date build dir used to be installed silently).
# Install protal, protal_map_utils, protal_profile_utils, qcmsa and simulate_metagenomes into prefix/bin.
# It also removes protal_baseline and protal_avx2, the two builds that earlier installs ran through a launcher.
install prefix="$HOME/.local": build-all
    mkdir -p {{prefix}}/bin
    rm -f {{prefix}}/bin/protal_baseline {{prefix}}/bin/protal_avx2
    cp {{build_dir}}/protal                         {{prefix}}/bin/protal
    cp {{build_dir}}/simulate_metagenomes           {{prefix}}/bin/simulate_metagenomes
    cp scripts/protal_map_utils                     {{prefix}}/bin/protal_map_utils
    cp scripts/protal_profile_utils                 {{prefix}}/bin/protal_profile_utils
    cp scripts/qcmsa.py                             {{prefix}}/bin/qcmsa
    chmod +x {{prefix}}/bin/protal {{prefix}}/bin/simulate_metagenomes {{prefix}}/bin/protal_map_utils {{prefix}}/bin/protal_profile_utils {{prefix}}/bin/qcmsa
    @echo "Installed to {{prefix}}/bin: $({{prefix}}/bin/protal --version)"
