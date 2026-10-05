# Reviewer report: GTDB build and training scripts (round 5, strain-fixes @ 39a8585)

Saved from the reviewer's hand-back (it could not write findings.md). Scripts in ../scripts/gtdb/; runs in
WSL ~/audit6/gtdb/ (b1-b7, es, ei, ej, k, dl/). Python ~/protal-train/bin/python (scikit-learn 1.9.0, numpy 1.26.4,
pandas 2.1.4); system python3 lacks scikit-learn/joblib. Offline: fake GTDB mirror on 127.0.0.1 (test_mini_db.gtdb_mirror
+ a Range server with fault injection), fake NCBI datasets. WSL restarted twice (other sessions); the first killed run b1.

Test world: gtdb_like_lineages.py --species 120 --archaea 0.1 --seed 1; simulate_gtdb_release.py --genomes_per_species 3
--genome_length 60000 --strain_divergence 0.002-0.015 --species_divergence 0.015-0.04 --seed 3; download_gtdb.py --species 30
--per_species 2 --rep_only_species 20 --batch 8 (110 genomes, 50 species); build_gtdb_database.py --inputs dl/inputs -t 2
--samples 2 --read-pairs 1000,20000 --read-setups 100:HS20:300:40,150:HS25:350:50 --species-per-sample 8-12 --archaea 1
--holdout-max-share 0.1 --holdout-clades phylum:1,class:1,family:1,genus:2 (4-7 min).

| id | sev | file:line | finding | evidence |
|---|---|---|---|---|
| H1 | High | mini_db/gtdb_to_protal_db.py:292-294, build_gtdb_database.py:434 | reruns blocked after any failure or kill once the background protal_db build has started writing: it leaves unique_kmers.tsv, index.prx.zst[.partial], database.protal.partial in protal_db; exclude_from_db copies every regular file of protal_db (except reference.fna/.map, full_reference.fna) into training_db, including the full DB's unique_kmers.tsv; the training_db build then exits 8, on every rerun (at GTDB scale after the converter reran). After a completed run the copy brings database.protal and build_metadata.tsv (overwritten, wasted I/O) | b1 killed during packing; rerun -> "Command failed (8)"; log "Invalid unique k-mer file .../training_db/unique_kmers.tsv, line 1: gene 30_46 is not in reference.map"; training_db held the copied files (unique_kmers.tsv 610388 B = protal_db's). early.sh s + fix_rerun.sh -> same |
| M1 | Medium | build_gtdb_database.py:254-280, 441-443, 479-480 | a background build failure is noticed only in finish(), after collection, training and parity | bgfail.sh: protal_db build exits 137 at 15:58:11; script builds training_db, collects, trains, parity, then "Command failed (137)" at 15:59:02 |
| M2 | Medium | build_gtdb_database.py:63-70, 254-279 | SIGTERM or an uncaught exception leaves both child builds running; a rerun then runs two builds per folder, deleting each other's files | term.sh: after kill -TERM both builds still ran; rerun showed 4 builds; its training build failed "reference.fna does not exist" |
| M3 | Medium | build_gtdb_database.py:397-451 | no resume: every rerun repeats converter and both index builds (~1.5-2 h each at GTDB scale); only the collector skips work | reruns.sh a: converter again, "Built training_db in 24 s", "Built protal_db in 111 s", collection 0 s; rebuilt training_db/database.protal same md5 |
| M4 | Medium | collect_training_data.py:321-333 | collector does not check which DB/options produced its design points: a rerun with another --seed/--holdout*/genome table/protal reuses dumps profiled against the old training DB, recomputes meta_novel_* from the new holdout, trains on the mix; only parity catches it, after training | E-d (--seed 2): 47 species left out (was 45), collection 0 s, trained, parity fails "gene 18_78 of the SAM header (@SQ) is not in the database". E-k: collector with --db b2/protal_db reused every dump, exit 0 |
| M5 | Medium | build_gtdb_database.py:322-323, 442, 454, 464 | missing tools found only when first used (protal, simulator, art_illumina, scikit-learn/joblib) | early.sh j: FileNotFoundError /nonexistent/protal after conversion; early.sh s: after "Built training_db"; early.sh i: ModuleNotFoundError joblib after the builds and the whole collection |
| M6 | Medium | download_gtdb.py:150-169 | a complete .part gets HTTP 416, retried 3x, exit - on every rerun; fix is deleting it by hand, message does not say so; a kill between transfer end and os.replace (MD5 of a multi-GB file) leaves such a file | dl.sh D2: 416 x3, "could not download ... with MD5", exit 1, .part left |
| M7 | Medium | download_gtdb.py:157-168 | no resume within a run: a dropped connection returns short data without raising; MD5 fails, .part deleted, next attempt from byte 0; three drops and exit | dl.sh D3 (server closes after 5000 B): "checksum mismatch, downloading again" x3, exit 1; 3 GETs, all range=None |
| L1 | Low | collect_training_data.py:209 | design points named rl<len>_p<pairs> only: setups sharing a read length (or a repeated depth) collide and are simulated into one dir concurrently (defaults distinct) | probes.py P1: 4 points -> 1 distinct name |
| L2 | Low | download_gtdb.py:349-358, 404-406 | with NCBI down, failed batches halved to single accessions (~2n-1 requests); run exits 0, says inputs ready | ncbi_down.sh: 213 calls, genomes 0 of 110, exit 0 |
| L3 | Low | build_gtdb_database.py:185-186, random_forest_cmdline.py:355-358 | a held-out clade can take all of one domain's simulated species (cap is a share of all species); the trainer's archaea warning needs >= 1 archaeal row | b2: c__Neurbelia (12 species) = every archaeon held out; domain table shows only Bacteria, no warning |
| L4 | Low | build_gtdb_database.py:387-389 | --inputs with --genome-table always rejected, message blames options not given | early.sh f |
| L5 | Low | build_gtdb_database.py:441-451 | two concurrent builds each get -t, as do collection ART jobs and protal beside the background build: peak 2 x -t threads (docs mention memory only) | protal_calls.log |
| L6 | Low (docs) | building-a-database.md:45,51-52,57,118; model-training.md:3,14,207; random_forest_cmdline.py:27 | docs say model.xml; converter writes model_pe.xml | b2 protal_db: model_pe/se/PB/ONT.xml, no model.xml |
| L7 | Low (docs) | download_gtdb.py:17,65 vs building-a-database.md:146 | rep-genome archive 127 GB in the script, 137 GB in the docs | text |
| L8 | Low | check_model_parity.py:210; collect_training_data.py:349; lineages.py:34-41 | latent: parity hard-codes alignments/; a sample prefix containing ".profile" cut short; a taxonomy cycle -> RecursionError | probes.py P3, P2 |
| L9 | Low | collector's protal.meta | training/ cannot be moved: absolute paths | E-d on a copy |
| L10 | Low (extrapolated) | gtdb_to_protal_db.py:237-256, 474-488 | disk at r226 undocumented: full_reference.fna in both protal_db and training_db, .convert_tmp another full copy until joined; ~3.0x the reference per copy -> ~85 GB per copy at r226: ~170 GB kept + ~85 GB transient, besides the release | b2 sizes |
| L11 | Low | tests | no end-to-end test of build_gtdb_database.py; H1, M1-M4, L3 never exercised | justfile, docs/development.md |

Verified to hold: end to end exit 0; database.protal's model_pe.xml identical to trained_model.xml, placeholders marked;
model_logs complete; PMML equals scikit-learn 1.9.0 exactly (trainer: max difference 0 on 219 rows; probes.py P4: 0.0 on
23k-63k rows, 4 forest settings, values at/next to every threshold; protal parity 4 samples, 100 taxa); tree_.value
fractions in 1.9 handled, single-leaf trees work, no division by zero on write_forest output, a layout change would trip
the trainer's own check; no label leakage (truth, prediction, probability, taxon, taxon_name, total_hits, meta_* excluded);
holdout: nested clades not redrawn, training_db drops held-out species (8835 kept, 4513 left out) keeping taxids,
deterministic (same md5 b1/b2/b6), max(lin, key=RANKS.index) safe; stop() kills the background build on a foreground
failure; collector resume after kill -9 of protal works (SAMs .partial, never reused; resumed table differs only in the last
bit of 4 features); simulator writes protal.meta last; download_gtdb.py: rerun downloads nothing it has, a server ignoring
Range handled, GNU tar 1.35 refuses .. members and symlink writes, failed extraction stops the run.
