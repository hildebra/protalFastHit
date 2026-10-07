#!/usr/bin/env python3
"""What build_gtdb_database.py's two collections hold at their defaults, from the collector's own planning functions
(units_of, streamed_simulations, pe_bytes, sample_bytes, pe_point_seconds, long_read_seconds): per collection and
--stream-above threshold, the streamed simulations (each a protal run of its own), the reads written to the disk, the
reads protal profiles, and the simulators' work on one core by the collector's cost constants.

Usage: plan_model.py [--scripts DIR] [--v15]
--scripts: the folder of collect_training_data.py and scenarios.py (default: this checkout's scripts/).
--v15: the r226 v15 build's design instead (4 scenarios of 6 + 3 samples, depth factors 1/2 to 2, training strains
0.3,0.1), to check the byte estimates against v15's logged "GB of reads" per protal run.
No files are read or written besides the imports; it runs in about a second."""
import argparse
import importlib
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "scripts"))

# build_gtdb_database.py's defaults (collect_command), training and test collection.
TRAINING = ["--samples", "12", "--read_pairs",
            "1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1",
            "--species_per_sample", "20-200", "--long_read_samples", "36", "--strains_per_species", "0.5,0.2",
            "--long_read_bases", "300000,1500000,6000000,30000000,150000000,1500000000:4,6000000000:2", "--seed", "1"]
TEST = ["--samples", "4", "--read_pairs", "500,2000,10000,50000,200000,1000000,5000000:2",
        "--species_per_sample", "10-300", "--long_read_samples", "8", "--strains_per_species", "0.5,0.2",
        "--long_read_bases", "150000,1000000,5000000,25000000,250000000,3000000000:2", "--seed", "1001"]
SCENARIOS = ("gut", "moderate", "soil", "soil_shallow", "host")
V15_SCENARIOS = ("gut", "soil", "soil_shallow", "host")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scripts", default=SCRIPTS)
    p.add_argument("--v15", action="store_true")
    p.add_argument("--free", default="",
                   help="GB free on the samples' disk once the training database is built, comma-separated: what "
                        "build_gtdb_database.py's --stream-above auto chooses (stream_threshold) with a 50 GB genome store, "
                        "the SAMs and profiles at KEPT_SHARE of the reads and --keep-free 30")
    a = p.parse_args()
    sys.path.insert(0, a.scripts)
    ctd = importlib.import_module("collect_training_data")
    scenarios = importlib.import_module("scenarios")

    extra, names, hold_in, hold_out = [], SCENARIOS, 10, 4
    if a.v15:
        names, hold_in, hold_out = V15_SCENARIOS, 6, 3
        defs = {n: {**scenarios.PRESETS[n], "depth_spread": 2} for n in names}
        fh = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(defs, fh)
        fh.close()
        extra = ["--scenario_file", fh.name]

    rows, sizes = [], []
    for what, design, n in (("training", TRAINING, hold_in), ("test", TEST, hold_out)):
        design = list(design)
        if a.v15 and what == "training":
            design[design.index("--strains_per_species") + 1] = "0.3,0.1"
        opts = ctd.parse_args(["--db", "x", "--genome_table", "x", "-o", "x", "--read_types", "pe,se,pb,ont",
                               "--scenarios", ",".join(f"{s}:{n}" for s in names), *design, *extra])
        points, units = ctd.units_of(opts)
        if hasattr(ctd, "simulation_bytes"):  # the collector since the change this report proposed
            sizes += ctd.simulation_bytes(units).values()
        pe_points = {u["point"]["name"]: u["point"] for u in units if not ctd.drawn(u)}
        drawn = [u for u in units if ctd.drawn(u)]

        def written(name):  # bytes a simulation writes: a paired-end point's both reads, a drawn unit's reads
            if name in pe_points:
                pt = pe_points[name]
                return sum(c + h for c, h in zip(ctd.community_pairs_of(pt), ctd.host_pairs_of(pt))) * 2 * \
                    int(pt["read_length"]) * ctd.PE_BYTES
            u = next(u for u in drawn if u["name"] == name)
            return sum(ctd.bases_of(u)) * ctd.DRAWN_BYTES

        sims = sorted(set(pe_points) | {u["name"] for u in drawn})
        total = sum(written(s) for s in sims)
        pe_work = sum(ctd.pe_point_seconds(pt, 1, opts) for pt in pe_points.values() if pt.get("reads", True))
        long_work = sum(sum(ctd.long_read_seconds({**u, "bases": b}) for b in ctd.bases_of(u)) for u in drawn)
        print(f"{what}: {len(sims)} simulations ({len(pe_points)} paired-end points, {len(drawn)} drawn units), "
              f"{sum(u['samples'] for u in units)} unit samples; reads {total / 1e9:.0f} GB "
              f"(pe points {sum(written(s) for s in pe_points) / 1e9:.0f} GB, drawn {sum(written(u['name']) for u in drawn) / 1e9:.0f} GB); "
              f"simulators by the collector's constants: pe {pe_work / 3600:.1f} core-h, drawn {long_work / 3600:.1f} core-h")
        for above in (2, 4, 8, 16, 0):
            opts.stream_above = above
            streamed = ctd.streamed_simulations(units, opts) if above else set()
            on_disk = sum(written(s) for s in sims if s not in streamed)
            long_runs = sum(1 for s in streamed if s not in pe_points)
            largest = max((written(s) for s in sims if s not in streamed), default=0)
            print(f"  --stream-above {above or 'none':>4}: {len(streamed):2d} streamed protal runs ({len(streamed) - long_runs} "
                  f"paired-end with their se, {long_runs} long-read or Ultima), {sum(written(s) for s in streamed) / 1e9:5.0f} GB "
                  f"streamed, {on_disk / 1e9:5.0f} GB written (largest simulation {largest / 1e9:.0f} GB)")
            rows.append((what, above, len(streamed), on_disk))
        if what == "training":
            for s in sims:
                print(f"    {s}: {written(s) / 1e9:.1f} GB")
    if a.free:
        import build_gtdb_database as build
        reads = sum(r for _, r in sizes)
        fixed = 50e9 + build.KEPT_SHARE * reads + 30e9
        for free in (float(f) * 1e9 for f in a.free.split(",")):
            threshold = build.stream_threshold(sizes, free - fixed)
            if threshold is None:
                print(f"--free {free / 1e9:g}: nothing fits ({fixed / 1e9:.0f} GB needed besides the reads)")
                continue
            streamed = [r for largest, r in sizes if threshold > 0 and largest > threshold]
            print(f"--free {free / 1e9:g}: --stream-above {threshold / 1e9:.3g} GB, {len(streamed)} of {len(sizes)} "
                  f"simulations streamed ({sum(streamed) / 1e9:.0f} GB), {(reads - sum(streamed)) / 1e9:.0f} GB written; "
                  f"{fixed / 1e9:.0f} GB besides the reads")


if __name__ == "__main__":
    main()
