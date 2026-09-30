#!/usr/bin/env python3
"""Download what building a protal database and training its model need, on a node with internet.

Run it once per GTDB release on a node that can reach the internet; the folder it writes serves
every later build of that release (build_gtdb_database.py --inputs), and a rerun downloads only
what is missing or does not match its checksum.

    python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226             # GTDB r226
    python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r220 --release 220
    python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir OUT ...

Two parts:
- GTDB's files (data.gtdb.ecogenomic.org, --mirror): taxonomy, metadata, and the marker genes of
  the representatives and of all genomes, checked against the release's MD5SUM.txt and extracted.
  Releases 207 and later (bac120 and ar53 marker sets); the newest point release (e.g. 214.1) unless
  --release names one. With --rep_genomes gtdb also GTDB's archive of all representative genomes
  (127 GB for r226).
- Genomes to simulate training samples from, from NCBI (the `datasets` CLI): GTDB distributes the
  representatives' genomes only, which are the database's own references; real samples hold other
  strains. From the metadata, per domain in GTDB's proportions: --species species with
  non-representative genomes that pass the quality filters (CheckM2 completeness and
  contamination), up to --per_species of those each, drawn at random; and --rep_only_species more
  species simulated from their representative only. With --rep_genomes ncbi (the default) the
  representatives of these species come from NCBI too, so the 127 GB archive is not needed. The
  simulator draws a species uniformly, then one of its genomes: a species with k strains is
  simulated from a strain k/(k+1) of the time.

Written to OUT/:
  release/                GTDB's files as build_gtdb_database.py --gtdb reads them
  genomes/                <accession>.fna.gz of the NCBI genomes
  genomes.tsv             accession, species, role (representative or strain), lineage, CheckM2
                          completeness and contamination, of the genomes in genomes/
  simulation_species.txt  the species to simulate from (build_gtdb_database.py --simulate-species)
  missing.txt             accessions NCBI did not deliver
  download.json           the release, the options, the counts and the checksums of what is there
"""

import argparse
import collections
import concurrent.futures
import datetime
import gzip
import hashlib
import http.client
import json
import os
import random
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request

MIRROR = "https://data.gtdb.ecogenomic.org/releases"
MARKER_SETS = ("bac120", "ar53")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-o", "--out", required=True, help="folder for this release's inputs, reused by later builds")
    p.add_argument("--release", default="226",
                   help="GTDB release: 226 (default), another from 207 on (e.g. 220), or a point release (214.1)")
    p.add_argument("--mirror", default=MIRROR, help=f"base URL of GTDB's releases (default {MIRROR})")
    p.add_argument("--rep_genomes", choices=["ncbi", "gtdb"], default="ncbi",
                   help="representatives' genomes of the simulated species from NCBI (default), or GTDB's archive "
                        "of all representative genomes (127 GB for r226)")
    p.add_argument("--no_genomes", action="store_true", help="GTDB's files only, no genomes from NCBI")
    p.add_argument("--dry_run", action="store_true", help="list GTDB's files to download, with their sizes, and stop")
    p.add_argument("--keep_archives", action="store_true",
                   help="keep the downloaded .tar.gz archives after extracting them (default: removed)")
    p.add_argument("--species", type=int, default=6000, help="species to download strains of (default 6000)")
    p.add_argument("--per_species", type=int, default=2, help="strains per species at most (default 2)")
    p.add_argument("--rep_only_species", type=int, default=2000,
                   help="further species, simulated from their representative only (default 2000: with the other "
                        "defaults, a simulated species is another strain about 50%% of the time)")
    p.add_argument("--min_completeness", type=float, default=90.0, help="CheckM2 completeness, %% (default 90)")
    p.add_argument("--max_contamination", type=float, default=5.0, help="CheckM2 contamination, %% (default 5)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--datasets", default="datasets", help="NCBI datasets binary (default: on PATH)")
    p.add_argument("--batch", type=int, default=500, help="genomes per NCBI request (default 500)")
    p.add_argument("-t", "--threads", type=int, default=8, help="parallel downloads and compression (default 8)")
    return p.parse_args(argv)


# ---- GTDB's files -----------------------------------------------------------------------------------------

def fetch_text(url):
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            return response.read().decode()
    except OSError as e:
        sys.exit(f"cannot read {url}: {e}")


def resolve_release(opts):
    """(release number as in file names, e.g. '226'; version folder, e.g. '226.0')."""
    text = opts.release.lower().removeprefix("r").removeprefix("release")
    number = text.split(".")[0]
    if not number.isdigit():
        sys.exit(f"--release {opts.release}: expected a number such as 226 or 214.1")
    if int(number) < 207:
        sys.exit(f"GTDB r{number} predates the ar53 marker set (it has ar122); protal databases need release 207 or later")
    if "." in text:
        return number, text
    listing = fetch_text(f"{opts.mirror}/release{number}/")
    versions = sorted({v for v in re.findall(rf'href="({number}\.\d+)/?"', listing)}, key=lambda v: int(v.split(".")[1]))
    if not versions:
        sys.exit(f"no GTDB release {number} at {opts.mirror}/release{number}/")
    return number, versions[-1]


def release_files(md5sums, number, rep_genomes):
    """The files of the release to download: [(path in the release, md5, extract?)]."""
    wanted = []

    def one_of(*names):
        for name in names:
            if name in md5sums:
                return name
        return None

    for mset in MARKER_SETS:
        taxonomy = one_of(f"{mset}_taxonomy_r{number}.tsv.gz", f"{mset}_taxonomy_r{number}.tsv")
        metadata = one_of(f"{mset}_metadata_r{number}.tsv.gz", f"{mset}_metadata_r{number}.tar.gz")
        reps = one_of(f"genomic_files_reps/{mset}_marker_genes_reps_r{number}.tar.gz")
        every = one_of(f"genomic_files_all/{mset}_marker_genes_all_r{number}.tar.gz")
        missing = [what for what, name in (("taxonomy", taxonomy), ("metadata", metadata),
                                           ("representatives' marker genes", reps), ("all marker genes", every)) if not name]
        if missing:
            sys.exit(f"GTDB r{number} has no {mset} " + ", ".join(missing) + (
                " (releases before 207 use the ar122 marker set; protal reads bac120 and ar53)" if mset == "ar53" else ""))
        wanted += [(taxonomy, False), (metadata, metadata.endswith(".tar.gz")), (reps, True), (every, True)]
    if rep_genomes == "gtdb":
        genomes = one_of(f"genomic_files_reps/gtdb_genomes_reps_r{number}.tar.gz")
        if not genomes:
            sys.exit(f"GTDB r{number} has no genomic_files_reps/gtdb_genomes_reps_r{number}.tar.gz")
        wanted.append((genomes, True))
    if "VERSION.txt" in md5sums:
        wanted.append(("VERSION.txt", False))
    return [(name, md5sums[name], extract) for name, extract in wanted]


def md5_of(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_size(response, start):
    """The whole file's size by a response to a request from byte start: from Content-Range (bytes a-b/total)
    or Content-Length; None if it does not say."""
    total = re.match(r"bytes \d+-\d+/(\d+)", response.headers.get("Content-Range", ""))
    if total:
        return int(total.group(1))
    length = response.headers.get("Content-Length")
    return start + int(length) if length and length.isdigit() else None


def download(url, dest, md5, attempts=3):
    """dest from url, resuming a partial download (dest.part, of this run or an earlier one), until its MD5 is
    md5. A connection that drops returns short data without an error: the next request goes on from there.
    `attempts` requests in a row that bring nothing new give up; a .part that turns out complete (the server
    has nothing after its end, 416) is checked and kept."""
    part, name = dest + ".part", os.path.basename(dest)
    failed = 0
    while failed < attempts:
        start = os.path.getsize(part) if os.path.exists(part) else 0
        request = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
        size = None
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                resumed = start and response.status == 206
                size = expected_size(response, start if resumed else 0)
                with open(part, "ab" if resumed else "wb") as fh:
                    shutil.copyfileobj(response, fh, 1 << 22)
        except urllib.error.HTTPError as e:
            if e.code != 416 or not start:
                failed += 1
                print(f"  {name}: {e} (attempt {failed} of {attempts})", flush=True)
                continue
            size = start  # nothing after the .part's end: it is whole, or longer than the file
        except (OSError, http.client.HTTPException) as e:
            failed += 1 if (os.path.getsize(part) if os.path.exists(part) else 0) == start else 0
            print(f"  {name}: {e}; going on from byte {os.path.getsize(part) if os.path.exists(part) else 0}", flush=True)
            continue
        have = os.path.getsize(part)
        if size is not None and have < size:
            failed += 1 if have == start else 0
            print(f"  {name}: {have} of {size} bytes (the connection dropped); going on", flush=True)
            continue
        if md5_of(part) == md5:
            os.replace(part, dest)
            return
        failed += 1
        print(f"  {name}: checksum mismatch, downloading again (attempt {failed} of {attempts})", flush=True)
        os.remove(part)
    sys.exit(f"could not download {url} with MD5 {md5}" +
             (f"; {part} holds what came, which a rerun goes on from" if os.path.exists(part) else ""))


def extracted_path(archive):
    """What the converter reads from an archive: bac120_marker_genes_reps_r226.tar.gz -> the folder
    bac120_marker_genes_reps_r226, bac120_metadata_r207.tar.gz -> bac120_metadata_r207.tsv."""
    base = archive[:-len(".tar.gz")]
    return base + ".tsv" if "_metadata_" in os.path.basename(base) else base


def extract(archive, into):
    command = ["tar", "-xf", archive, "-C", into]
    command[1:1] = ["-I", "pigz"] if shutil.which("pigz") else ["-z"]
    if subprocess.run(command).returncode != 0:
        sys.exit(f"extracting {archive} failed")
    if not os.path.exists(extracted_path(archive)):
        sys.exit(f"{archive} did not extract to {extracted_path(archive)}, which the converter reads")


def get_release(opts, state):
    number, version = resolve_release(opts)
    base = f"{opts.mirror}/release{number}/{version}"
    release_dir = os.path.join(opts.out, "release")
    os.makedirs(release_dir, exist_ok=True)
    md5sums = {}
    for line in fetch_text(f"{base}/MD5SUM.txt").splitlines():
        parts = line.split()
        if len(parts) == 2:
            md5sums[parts[1].removeprefix("./")] = parts[0]
    files = state.setdefault("files", {})
    print(f"GTDB r{number} ({version}) from {base}", flush=True)
    if opts.dry_run:
        total = 0
        for name, _md5, _unpack in release_files(md5sums, number, opts.rep_genomes):
            with urllib.request.urlopen(urllib.request.Request(f"{base}/{name}", method="HEAD"), timeout=120) as r:
                size = int(r.headers.get("Content-Length", 0))
            total += size
            print(f"  {name}: {size / 1e9:.2f} GB", flush=True)
        print(f"  total {total / 1e9:.1f} GB from GTDB; genomes from NCBI come on top (about 4 MB each)")
        sys.exit(0)
    for name, md5, unpack in release_files(md5sums, number, opts.rep_genomes):
        dest = os.path.join(release_dir, name)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        done = files.get(name, {})
        there = os.path.exists(extracted_path(dest)) if unpack else os.path.exists(dest)
        if done.get("md5") == md5 and there:
            print(f"  {name}: already there", flush=True)
            continue
        if not (os.path.exists(dest) and md5_of(dest) == md5):
            print(f"  {name}: downloading", flush=True)
            download(f"{base}/{name}", dest, md5)
        entry = {"md5": md5, "bytes": os.path.getsize(dest)}
        if unpack:
            print(f"  {name}: extracting", flush=True)
            extract(dest, os.path.dirname(dest))
            entry["extracted"] = True
            if not opts.keep_archives:
                os.remove(dest)
        files[name] = entry
        save_state(opts, state)
    state["release"] = {"number": number, "version": version, "url": base}
    return number


# ---- genomes from NCBI ------------------------------------------------------------------------------------

def read_metadata(gtdb, release):
    """Every genome: (NCBI accession, species, lineage, is representative, completeness, contamination)."""
    genomes = []
    for mset in MARKER_SETS:
        path = next((os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}") for ext in (".tsv", ".tsv.gz")
                     if os.path.exists(os.path.join(gtdb, f"{mset}_metadata_r{release}{ext}"))), None)
        if path is None:
            continue
        with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}
            completeness = col.get("checkm2_completeness", col.get("checkm_completeness"))
            contamination = col.get("checkm2_contamination", col.get("checkm_contamination"))
            for name in ("accession", "gtdb_taxonomy", "gtdb_representative"):
                if name not in col:
                    sys.exit(f"{path} has no {name} column")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                lineage = f[col["gtdb_taxonomy"]]
                accession = f[col["accession"]]
                genomes.append((accession[3:] if accession[:3] in ("RS_", "GB_") else accession, lineage.split(";")[-1],
                                lineage, f[col["gtdb_representative"]] == "t",
                                float(f[completeness]) if completeness is not None and f[completeness] not in ("", "none") else 100.0,
                                float(f[contamination]) if contamination is not None and f[contamination] not in ("", "none") else 0.0))
    if not genomes:
        sys.exit(f"no bac120/ar53_metadata_r{release}.tsv[.gz] in {gtdb}")
    return genomes


def pick(genomes, opts):
    """(strains [(accession, species, lineage, completeness, contamination)], species with strains,
    species simulated from their representative only, {species: domain}, {species: representative row})."""
    rng = random.Random(opts.seed)
    strains = collections.defaultdict(list)
    domain, lineage, representative = {}, {}, {}
    for acc, species, lin, is_rep, comp, cont in genomes:
        domain[species], lineage[species] = lin.split(";")[0], lin
        if is_rep:
            representative[species] = (acc, species, lin, comp, cont)
        elif comp >= opts.min_completeness and cont <= opts.max_contamination:
            strains[species].append((acc, comp, cont))
    with_strains = collections.defaultdict(list)
    without = collections.defaultdict(list)
    for species in sorted(domain):
        (with_strains if strains.get(species) else without)[domain[species]].append(species)
    chosen, rep_only = [], []
    for d in sorted(set(domain.values())):
        share = sum(1 for s in domain if domain[s] == d) / len(domain)
        picked_here = set(rng.sample(with_strains[d], min(len(with_strains[d]), round(opts.species * share))))
        chosen += sorted(picked_here)
        rest = [s for s in with_strains[d] if s not in picked_here] + without[d]
        rep_only += rng.sample(rest, min(len(rest), round(opts.rep_only_species * share)))
    picked = []
    for species in sorted(chosen):
        options = sorted(strains[species])
        for acc, comp, cont in rng.sample(options, min(opts.per_species, len(options))):
            picked.append((acc, species, lineage[species], comp, cont))
    return picked, sorted(chosen), sorted(rep_only), domain, representative


def gzip_into(source, dest):
    with open(source, "rb") as fin, gzip.open(dest + ".part", "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, 1 << 22)
    os.replace(dest + ".part", dest)


def ncbi_batch(opts, accessions, work):
    """Downloads accessions with datasets into work: {accession: FASTA path} of those delivered, and why the
    request failed ("" if it did not)."""
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    listing = os.path.join(work, "accessions.txt")
    with open(listing, "w") as fh:
        fh.writelines(a + "\n" for a in accessions)
    zip_path, unpacked = os.path.join(work, "batch.zip"), os.path.join(work, "batch")
    steps = ([opts.datasets, "download", "genome", "accession", "--inputfile", listing, "--include", "genome",
              "--dehydrated", "--filename", zip_path],
             ["unzip", "-o", "-q", zip_path, "-d", unpacked],
             [opts.datasets, "rehydrate", "--directory", unpacked, "--max-workers", str(max(1, min(opts.threads, 30)))])
    log_path = os.path.join(work, "log.txt")
    with open(log_path, "w") as log:
        for command in steps:
            rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
            if rc != 0:
                log.close()
                with open(log_path, errors="replace") as fh:
                    tail = " | ".join(line.strip() for line in fh.read().strip().splitlines()[-3:])
                return {}, f"{os.path.basename(command[0])} {command[1]} failed ({rc}): {tail}"
    found = {}
    for root, _dirs, names in os.walk(unpacked):
        for name in names:
            m = re.match(r"(GC[AF]_\d{9}\.\d+)_.*_genomic\.fna$", name)
            if m and m.group(1) in accessions:
                found[m.group(1)] = os.path.join(root, name)
    return found, ""


def get_genomes(opts, state, release):
    gtdb = os.path.join(opts.out, "release")
    genomes = read_metadata(gtdb, release)
    strains, chosen, rep_only, domain, representative = pick(genomes, opts)
    pool = sorted(chosen + rep_only)
    wanted = [(acc, species, "strain", lin, comp, cont) for acc, species, lin, comp, cont in strains]
    if opts.rep_genomes == "ncbi":
        wanted += [(representative[s][0], s, "representative", *representative[s][2:]) for s in pool if s in representative]
    folder = os.path.join(opts.out, "genomes")
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(opts.out, "simulation_species.txt"), "w") as fh:
        fh.writelines(s + "\n" for s in pool)
    have = {name[:-len(".fna.gz")] for name in os.listdir(folder) if name.endswith(".fna.gz")}
    todo = [w[0] for w in wanted if w[0] not in have]
    print(f"genomes: {len(wanted)} wanted ({len(strains)} strains of {len(chosen)} species"
          + (f", the representatives of {len(pool)} species" if opts.rep_genomes == "ncbi" else "")
          + f"), {len(wanted) - len(todo)} already there, {len(todo)} to download", flush=True)
    if todo and not shutil.which(opts.datasets):
        sys.exit(f"the NCBI datasets CLI ({opts.datasets}) is not on PATH: install it (conda: ncbi-datasets-cli; see "
                 "envs/protal-db-build.yaml)")
    missing = []
    work = os.path.join(opts.out, "ncbi_batch")
    batches = [todo[i:i + opts.batch] for i in range(0, len(todo), opts.batch)]
    # A request that fails is halved, down to single accessions, to find the accessions NCBI refuses; one
    # such accession fails at most ~log2(--batch) requests in a row. More failures in a row mean that NCBI
    # (or datasets) is not working, which halving would only ask about ~2n times.
    in_a_row, most_in_a_row, error = 0, max(1, opts.batch).bit_length() + 4, ""
    with concurrent.futures.ProcessPoolExecutor(max(1, opts.threads)) as pool_executor:
        while batches:
            batch = batches.pop(0)
            found, why = ncbi_batch(opts, batch, work)
            if why:
                in_a_row, error = in_a_row + 1, why
                if in_a_row >= most_in_a_row:
                    sys.exit(f"{in_a_row} NCBI requests failed in a row, the last: {error}. Is NCBI reachable from "
                             f"here, and does {opts.datasets} work? A rerun downloads what is still missing")
            else:
                in_a_row = 0
            if not found and len(batch) > 1:  # a failed request: halve it, down to single accessions
                batches[:0] = [batch[:len(batch) // 2], batch[len(batch) // 2:]]
                continue
            list(pool_executor.map(gzip_into, found.values(), [os.path.join(folder, a + ".fna.gz") for a in found]))
            missing += [a for a in batch if a not in found]
            print(f"  {len(found)} of {len(batch)} delivered; {len(batches)} requests left", flush=True)
    shutil.rmtree(work, ignore_errors=True)
    have = {name[:-len(".fna.gz")] for name in os.listdir(folder) if name.endswith(".fna.gz")}
    if wanted and not any(w[0] in have for w in wanted):
        sys.exit(f"NCBI delivered none of the {len(wanted)} genomes" + (f"; the last failure: {error}" if error else ""))
    with open(os.path.join(opts.out, "genomes.tsv"), "w") as fh:
        fh.write("accession\tspecies\trole\tlineage\tcheckm2_completeness\tcheckm2_contamination\n")
        fh.writelines("\t".join(map(str, w)) + "\n" for w in wanted if w[0] in have)
    with open(os.path.join(opts.out, "missing.txt"), "w") as fh:
        fh.writelines(a + "\n" for a in sorted(w[0] for w in wanted if w[0] not in have))
    delivered = [w for w in wanted if w[0] in have]
    per_species = collections.Counter(w[1] for w in delivered if w[2] == "strain")
    other = sum(k / (k + 1) for k in per_species.values()) / max(1, len(pool))
    state["genomes"] = {"wanted": len(wanted), "delivered": len(delivered), "missing": len(wanted) - len(delivered),
                        "strains": sum(per_species.values()), "species_with_strains": len(per_species),
                        "simulation_species": len(pool), "other_strain_share": round(other, 3),
                        "by_domain": dict(collections.Counter(domain[s] for s in pool))}
    print(f"genomes: {len(delivered)} of {len(wanted)} in {folder} ({len(wanted) - len(delivered)} missing, see "
          f"missing.txt); a simulated species is another strain than the representative {100 * other:.0f}% of the time",
          flush=True)


# ---- main -------------------------------------------------------------------------------------------------

def load_state(opts):
    path = os.path.join(opts.out, "download.json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save_state(opts, state):
    path = os.path.join(opts.out, "download.json")
    with open(path + ".part", "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(path + ".part", path)


def main(argv=None):
    opts = parse_args(argv)
    os.makedirs(opts.out, exist_ok=True)
    state = load_state(opts)
    release = get_release(opts, state)
    state["options"] = {k: v for k, v in vars(opts).items() if k not in ("out", "datasets", "threads")}
    save_state(opts, state)
    if not opts.no_genomes:
        get_genomes(opts, state, release)
    state["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    save_state(opts, state)
    print(f"Inputs for GTDB r{release}: {opts.out} (build_gtdb_database.py --inputs {opts.out})", flush=True)


if __name__ == "__main__":
    main()
