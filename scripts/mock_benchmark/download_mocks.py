#!/usr/bin/env python3
"""Download the samples of the benchmark on public mock communities (samples.tsv, prepare_tables.py) from ENA, for
score_mocks.py. Run it where the internet is reachable (a login or transfer node); it resumes where it stopped.

Each run's FASTQ comes from ENA (its own FASTQ, else the submitted FASTQ), into OUT/reads/<set>/<sample>/
(<run>_1.fq.gz and <run>_2.fq.gz, or <run>.fq.gz). Deep runs are cut while they stream: a sample keeps its first
--max-pairs read pairs (single-end reads: reads) or, for long reads, its first --max-bases bases, so a 300 Gb run costs
the download of what is kept. A run within the limit is copied as it is and checked against ENA's MD5. A finished run
leaves <run>.json (files, reads, bases, whether cut); OUT/reads/downloaded.tsv lists the samples whose runs are all
there (set, sample, read_type, first, second, reads, bases, cut), which score_mocks.py reads.

  --dry-run     what the chosen sets hold and what they would take on disk with the limits, without downloading
  --sets        sets of samples.tsv (comma-separated; default all): d6331_dna, d6331_host, meslier, treichel,
                d6331_hifi, atcc1003_hifi, d6331_ont, d6322_basecalls, strain_mixes, nibsc_multilab, hcom2_invitro,
                hcom1a_dropouts, plant_syncoms
  --per-set N   only the first N samples of each set (a quick look at a large set); --only REGEX: samples by name

The first reads of a run stand in for a random subsample: an Illumina run's order follows the flow cell's tiles, a
long-read run's the time of sequencing, neither the organisms. Standard library only.

Usage:
  python3 scripts/mock_benchmark/download_mocks.py -o /scratch/mocks --dry-run
  python3 scripts/mock_benchmark/download_mocks.py -o /scratch/mocks --sets d6331_dna,meslier,d6331_hifi --jobs 8
"""
import argparse
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import threading
import time
import urllib.parse
import urllib.request
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS = ["run_accession", "library_layout", "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
              "submitted_ftp", "submitted_md5", "submitted_bytes"]
LONG = ("pb", "ont")
CHUNK = 1 << 22  # 4 MiB
RETRY_WAIT = 10  # seconds before the first retry of a request (then 20, 30)
DOWNLOADED_COLUMNS = ["set", "sample", "read_type", "first", "second", "reads", "bases", "cut"]
print_lock = threading.Lock()


def log(message):
    with print_lock:
        print(message, file=sys.stderr, flush=True)


def read_samples(path, sets, per_set, only=""):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    known = list(dict.fromkeys(r["set"] for r in rows))
    if sets:
        unknown = [s for s in sets if s not in known]
        if unknown:
            sys.exit(f"unknown set(s) {', '.join(unknown)}; samples.tsv has {', '.join(known)}")
        rows = [r for r in rows if r["set"] in sets]
    if only:
        rows = [r for r in rows if re.search(only, r["sample"])]
    if per_set:
        taken = collections.Counter()
        kept = []
        for r in rows:
            if taken[r["set"]] < per_set:
                kept.append(r)
                taken[r["set"]] += 1
        rows = kept
    return rows


def ena_metadata(runs, cache):
    """ENA's file metadata of the runs ({run: dict}), cached in a TSV beside the reads."""
    known = {}
    if os.path.exists(cache):
        with open(cache) as fh:
            known = {r["run_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    missing = [r for r in runs if r not in known]
    for i in range(0, len(missing), 200):
        batch = missing[i:i + 200]
        data = urllib.parse.urlencode({"result": "read_run", "includeAccessions": ",".join(batch),
                                       "fields": ",".join(ENA_FIELDS), "format": "tsv", "limit": 0}).encode()
        text = request(ENA_SEARCH, data=data).read().decode()
        for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
            known[r["run_accession"]] = r
        log(f"ENA metadata: {min(i + 200, len(missing))} of {len(missing)} runs")
    absent = [r for r in runs if r not in known]
    if absent:
        sys.exit(f"ENA does not know run(s) {', '.join(absent[:10])}")
    tmp = cache + ".part"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ENA_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in sorted(known.values(), key=lambda r: r["run_accession"]):
            w.writerow(r)
    os.replace(tmp, cache)
    return known


def request(url, data=None, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "protal-mock-benchmark"})
            return urllib.request.urlopen(req, timeout=300)
        except OSError as err:
            if attempt == tries - 1:
                raise
            log(f"  retry {url}: {err}")
            time.sleep(RETRY_WAIT * (attempt + 1))


def run_files(meta, read_type):
    """[(url, md5, bytes)] of the first mates (or single reads) and of the second mates ([] for single reads)."""
    def files(field):
        urls = [u for u in meta[field + "_ftp"].split(";") if u]
        md5s = meta[field + "_md5"].split(";") + [""] * len(urls)
        sizes = meta[field + "_bytes"].split(";") + [""] * len(urls)
        out = []
        for u, m, b in zip(urls, md5s, sizes):
            if re.search(r"\.(fastq|fq)(\.gz)?$", u):
                out.append(("https://" + u if "://" not in u else u, m, int(b) if b.isdigit() else 0))
        return out

    for field in ("fastq", "submitted"):
        found = files(field)
        if not found:
            continue
        if read_type == "pe":
            first = [f for f in found if re.search(r"(_|_R)1(_\d+)?\.(fastq|fq)(\.gz)?$", f[0])]
            second = [f for f in found if re.search(r"(_|_R)2(_\d+)?\.(fastq|fq)(\.gz)?$", f[0])]
            if len(first) == 1 and len(second) == 1:
                return first, second
        else:
            single = [f for f in found if not re.search(r"(_|_R)2(_\d+)?\.(fastq|fq)(\.gz)?$", f[0])]
            if single:
                return single[:1], []
    raise ValueError(f"{meta['run_accession']}: no {read_type} FASTQ at ENA (fastq_ftp '{meta['fastq_ftp']}', "
                     f"submitted '{meta['submitted_ftp']}')")


class Limit:
    """What a sample may still take: read pairs (or reads), or bases for long reads; None for no limit."""

    def __init__(self, reads=None, bases=None):
        self.reads, self.bases = reads, bases

    def fits(self, reads, bases):
        return (self.reads is None or reads <= self.reads) and (self.bases is None or bases <= self.bases)

    def take(self, reads, bases):
        if self.reads is not None:
            self.reads = max(0, self.reads - reads)
        if self.bases is not None:
            self.bases = max(0, self.bases - bases)

    def spent(self):
        return self.reads == 0 or self.bases == 0


def copy_whole(url, md5, path):
    """Copy a file as it is, checking ENA's MD5; returns its size."""
    digest = hashlib.md5()
    size = 0
    with request(url) as response, open(path + ".part", "wb") as out:
        while True:
            block = response.read(CHUNK)
            if not block:
                break
            digest.update(block)
            out.write(block)
            size += len(block)
    if md5 and digest.hexdigest() != md5:
        os.remove(path + ".part")
        raise IOError(f"{url}: MD5 {digest.hexdigest()} is not ENA's {md5}")
    os.replace(path + ".part", path)
    return size


def copy_records(url, path, max_reads=None, max_bases=None):
    """Stream a gzip (or plain) FASTQ and keep its first records: up to max_reads, or until max_bases bases are kept
    (the record that passes the limit is kept). Returns (reads, bases, whether the file had more)."""
    reads = bases = 0
    more = False
    with request(url) as response:
        source = gzip.GzipFile(fileobj=response) if url.endswith(".gz") else response
        with gzip.open(path + ".part", "wb", compresslevel=1) as out:
            rest = b""
            line_no = 0  # lines of the current record seen
            length = 0  # its bases
            done = False
            ended = False
            while not done and not ended:
                block = source.read(CHUNK)
                if not block:  # the end: a last line without its newline gets one
                    ended = True
                    if not rest:
                        break
                    block = b"\n"
                data = rest + block
                cut = data.rfind(b"\n") + 1
                data, rest = data[:cut], data[cut:]
                lines = data.split(b"\n")[:-1]
                keep = 0  # lines of data to write
                for i, line in enumerate(lines):
                    if line_no == 1:
                        bases += len(line)
                        length = len(line)
                    elif line_no == 3 and len(line) != length:
                        raise IOError(f"{url}: record {reads + 1} has {len(line)} qualities for {length} bases")
                    line_no += 1
                    if line_no == 4:
                        line_no = 0
                        reads += 1
                        keep = i + 1
                        if (max_reads is not None and reads >= max_reads) or (
                                max_bases is not None and bases >= max_bases):
                            done = True
                            break
                if done:
                    out.write(b"\n".join(lines[:keep]) + b"\n")
                    more = bool(lines[keep:]) or bool(rest) or bool(source.read(1))
                else:
                    out.write(data)
            if not done and line_no:
                raise IOError(f"{url}: the file ends inside a record")
    os.replace(path + ".part", path)
    return reads, bases, more


def fetch_sample(sample, metas, reads_dir, max_pairs, max_bases):
    """Fetch a sample's runs; returns its downloaded.tsv row, or None if a run failed."""
    folder = os.path.join(reads_dir, sample["set"], sample["sample"])
    os.makedirs(folder, exist_ok=True)
    long_reads = sample["read_type"] in LONG
    limit = Limit(bases=max_bases) if long_reads else Limit(reads=max_pairs)
    firsts, seconds = [], []
    total_reads = total_bases = 0
    cut_any = False
    for run in sample["runs"].split(","):
        done_path = os.path.join(folder, run + ".json")
        if os.path.exists(done_path):
            with open(done_path) as fh:
                done = json.load(fh)
            if all(os.path.exists(f) for f in done["files"]):
                firsts.append(done["files"][0])
                if len(done["files"]) > 1:
                    seconds.append(done["files"][1])
                total_reads += done["reads"]
                total_bases += done["bases"]
                cut_any = cut_any or done["cut"]
                limit.take(done["reads"], done["bases"])
                continue
        if limit.spent():
            break
        meta = metas[run]
        try:
            first, second = run_files(meta, sample["read_type"])
        except ValueError as err:
            log(f"{sample['sample']}: {err}")
            return None
        run_reads = int(meta["read_count"] or 0)
        run_bases = int(meta["base_count"] or 0)
        paired = bool(second)
        names = [os.path.join(folder, f"{run}_1.fq.gz"), os.path.join(folder, f"{run}_2.fq.gz")] if paired else [
            os.path.join(folder, f"{run}.fq.gz")]
        sources = first + second
        whole = run_reads > 0 and limit.fits(run_reads, run_bases) and all(u.endswith(".gz") for u, _, _ in sources)
        try:
            if whole:
                for (url, md5, _), path in zip(sources, names):
                    copy_whole(url, md5, path)
                kept_reads, kept_bases, cut = run_reads, run_bases, False
            else:
                counts = []
                for (url, _, _), path in zip(sources, names):
                    if long_reads:
                        counts.append(copy_records(url, path, max_bases=limit.bases))
                    else:
                        counts.append(copy_records(url, path, max_reads=limit.reads))
                if paired and counts[0][0] != counts[1][0]:
                    raise IOError(f"{run}: the mates' files kept {counts[0][0]} and {counts[1][0]} reads")
                kept_reads = counts[0][0]
                kept_bases = sum(c[1] for c in counts)
                cut = any(c[2] for c in counts)
        except (OSError, ValueError, EOFError, zlib.error) as err:  # EOFError, zlib.error: a damaged gzip stream
            log(f"{sample['sample']}: run {run} failed: {err}")
            for path in names:
                for p in (path, path + ".part"):
                    if os.path.exists(p):
                        os.remove(p)
            return None
        with open(done_path, "w") as fh:
            json.dump({"files": names, "reads": kept_reads, "bases": kept_bases, "cut": cut}, fh)
        firsts.append(names[0])
        if paired:
            seconds.append(names[1])
        total_reads += kept_reads
        total_bases += kept_bases
        cut_any = cut_any or cut
        limit.take(kept_reads, kept_bases)
        log(f"{sample['sample']}: {run} {kept_reads:,} {'pairs' if paired else 'reads'}, {kept_bases / 1e9:.2f} Gb"
            f"{' (cut)' if cut else ''}")
    if not firsts:
        return None
    return {"set": sample["set"], "sample": sample["sample"], "read_type": sample["read_type"],
            "first": ",".join(firsts), "second": ",".join(seconds) if seconds else "-", "reads": total_reads,
            "bases": total_bases, "cut": int(cut_any)}


def estimate(samples, metas, max_pairs, max_bases):
    """Per set: samples, runs, ENA's compressed bytes, and the bytes the limits would keep."""
    table = collections.OrderedDict()
    for s in samples:
        entry = table.setdefault(s["set"], [0, 0, 0, 0.0, 0.0])
        entry[0] += 1
        limit = Limit(bases=max_bases) if s["read_type"] in LONG else Limit(reads=max_pairs)
        for run in s["runs"].split(","):
            meta = metas[run]
            try:
                first, second = run_files(meta, s["read_type"])
            except ValueError:
                continue
            size = sum(b for _, _, b in first + second)
            reads, bases = int(meta["read_count"] or 0), int(meta["base_count"] or 0)
            share = 1.0
            if s["read_type"] in LONG and limit.bases is not None and bases:
                share = min(1.0, limit.bases / bases)
            elif limit.reads is not None and reads:
                share = min(1.0, limit.reads / reads)
            limit.take(reads * share, bases * share)
            entry[1] += 1
            entry[2] += bases
            entry[3] += size
            entry[4] += size * share * (1.15 if share < 1 else 1.0)  # gzip -1 is ~15% larger than ENA's files
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True, help="output folder (reads go to OUT/reads)")
    ap.add_argument("--samples", default=os.path.join(HERE, "samples.tsv"))
    ap.add_argument("--sets", default="", help="comma-separated sets (default: all)")
    ap.add_argument("--per-set", type=int, default=0, help="only the first N samples of each set")
    ap.add_argument("--only", default="", help="only the samples whose name matches this regular expression")
    ap.add_argument("--max-pairs", type=float, default=20e6,
                    help="read pairs (single-end: reads) kept per short-read sample (default 20M; 0: all)")
    ap.add_argument("--max-bases", type=float, default=20e9,
                    help="bases kept per long-read sample (default 20 Gb; 0: all)")
    ap.add_argument("--jobs", type=int, default=4, help="samples fetched at once (default 4)")
    ap.add_argument("--dry-run", action="store_true", help="print what the sets would take, download nothing")
    args = ap.parse_args()

    max_pairs = int(args.max_pairs) or None
    max_bases = int(args.max_bases) or None
    sets = [s for s in args.sets.split(",") if s]
    samples = read_samples(args.samples, sets, args.per_set, args.only)
    reads_dir = os.path.join(os.path.abspath(args.out), "reads")
    os.makedirs(reads_dir, exist_ok=True)
    runs = sorted({r for s in samples for r in s["runs"].split(",")})
    log(f"{len(samples)} samples, {len(runs)} runs")
    metas = ena_metadata(runs, os.path.join(reads_dir, "ena_runs.tsv"))

    table = estimate(samples, metas, max_pairs, max_bases)
    print("set\tsamples\truns\tGb sequenced\tENA GB\tGB kept (estimate)")
    for name, (n, r, bases, size, kept) in table.items():
        print(f"{name}\t{n}\t{r}\t{bases / 1e9:.1f}\t{size / 1e9:.1f}\t{kept / 1e9:.1f}")
    need = sum(e[4] for e in table.values())
    free = shutil.disk_usage(reads_dir).free
    print(f"total\t{len(samples)}\t{len(runs)}\t\t\t{need / 1e9:.1f}  (free: {free / 1e9:.0f} GB)")
    if args.dry_run:
        return
    if need > free:
        log("warning: the estimate exceeds the free space")

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch_sample, s, metas, reads_dir, max_pairs, max_bases): s for s in samples}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                rows.append(result)
    path = os.path.join(reads_dir, "downloaded.tsv")
    previous = {}
    if os.path.exists(path):
        with open(path) as fh:
            previous = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    previous.update({r["sample"]: r for r in rows})
    with open(path + ".part", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=DOWNLOADED_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sorted(previous.values(), key=lambda r: (r["set"], r["sample"])):
            w.writerow(r)
    os.replace(path + ".part", path)
    failed = len(samples) - len(rows)
    log(f"{len(rows)} samples complete{f', {failed} failed (run again to retry)' if failed else ''}; listed in {path}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
