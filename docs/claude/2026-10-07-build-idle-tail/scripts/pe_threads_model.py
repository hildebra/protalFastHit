#!/usr/bin/env python3
"""The threads collect_training_data.py gives each paired-end design point (pe_threads), with its own cost estimate
(pe_point_seconds: 20 s per sample and 96 us per 150 bp pair) and with one that also counts the sample's genomes
(G seconds per genome and sample: decompress, ART's start and its reference), and each point's time on those threads
by the second estimate, as the simulator splits threads: one worker per sample (at most the threads), a sample's
genomes on threads / workers threads, a worker's samples one after the other.

Usage: pe_threads_model.py COLLECTOR_DIR [--slots 84] [--genome_seconds 0.15] [--pair_seconds 96e-6] [--test]
COLLECTOR_DIR: the scripts folder holding collect_training_data.py and scenarios.py of 9ecf7ee (git archive), whose
estimate ignored the genomes; later collectors count them already (no files are read or written)."""
import argparse
import importlib
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scripts")
    p.add_argument("--slots", type=int, default=84)
    p.add_argument("--genome_seconds", type=float, default=0.15,
                   help="seconds per genome and sample on one node core (v14: 0.15 with contention; laptop 0.085)")
    p.add_argument("--pair_seconds", type=float, default=96e-6,
                   help="seconds per 150 bp read pair on one core (laptop 96e-6; v14 node ~53e-6)")
    p.add_argument("--test", action="store_true", help="the build's test collection instead of the training one")
    a = p.parse_args()
    sys.path.insert(0, a.scripts)
    ctd = importlib.import_module("collect_training_data")
    scenarios = importlib.import_module("scenarios")

    # As build_gtdb_database.py's collect_command at its defaults (training, or the test set with seed + 1000).
    if a.test:
        design = ["--samples", "4", "--read_pairs", "500,2000,10000,50000,200000,1000000,5000000:2",
                  "--species_per_sample", "10-300", "--long_read_samples", "8",
                  "--long_read_bases", "150000,1000000,5000000,25000000,250000000,3000000000:2", "--seed", "1001",
                  "--scenarios", "gut:3,soil:3,soil_shallow:3,host:3"]
    else:
        design = ["--samples", "12", "--read_pairs",
                  "1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1",
                  "--species_per_sample", "20-200", "--long_read_samples", "36",
                  "--long_read_bases", "300000,1500000,6000000,30000000,150000000,1500000000:4,6000000000:2",
                  "--seed", "1", "--scenarios", "gut:6,soil:6,soil_shallow:6,host:6"]
    opts = ctd.parse_args(["--db", "x", "--genome_table", "x", "-o", "x", "--read_types", "pe,se,pb,ont",
                           "--strains_per_species", "0.3,0.1", *design])
    points, units = ctd.units_of(opts)
    strains = 1.4  # 0.3 and 0.1: a second and a third strain

    def genomes(point):
        species = point["definition"]["species"] if point.get("scenario") else opts.species_per_sample
        lo, hi = scenarios.species_bounds(species)
        return (lo + hi) / 2 * strains

    def pairs_of(point):
        if point.get("community_pairs_of"):
            return [c + h for c, h in zip(ctd.community_pairs_of(point), ctd.host_pairs_of(point))]
        return [float(point["read_pairs"])] * point["samples"]

    def sample_seconds(point, pairs):
        return genomes(point) * a.genome_seconds + pairs * int(point["read_length"]) / 150 * a.pair_seconds

    def true_cost(point):
        return sum(sample_seconds(point, p) for p in pairs_of(point))

    def wall(point, threads):
        """The point's time on `threads` threads, split as write_all_reads splits them (static, samples in order)."""
        samples = [sample_seconds(point, p) for p in pairs_of(point)]
        workers = min(len(samples), threads)
        share = [threads // workers + (1 if w < threads % workers else 0) for w in range(workers)]
        ends = [0.0] * workers
        for s in samples:  # the next sample to the worker free first (the atomic counter)
            w = min(range(workers), key=lambda k: ends[k])
            ends[w] += s / share[w]
        return max(ends)

    pe = [(i, pt) for i, pt in enumerate(points) if pt.get("reads", True) and float(pt["read_pairs"]) > 0]
    long_pending = [(i, u) for i, u in enumerate(u for u in units if ctd.drawn(u))]
    pe_work = sum(ctd.pe_point_seconds(pt, 1) for _, pt in pe)
    long_work = sum(sum(ctd.long_read_seconds({**u, "bases": b}) for b in ctd.bases_of(u)) *
                    (2 if u["setup"]["method"] in ctd.PBSIM_METHODS else 1) for _, u in long_pending)
    pe_slots = max(1, round(a.slots * pe_work / ((pe_work + long_work) or 1)))
    now = ctd.pe_threads(pe, pe_slots)

    # The same split with the genomes counted (pe_point_seconds patched for the call).
    original = ctd.pe_point_seconds
    ctd.pe_point_seconds = lambda pt, threads, *rest: true_cost(pt) / max(1, threads)
    pe_work2 = sum(true_cost(pt) for _, pt in pe)
    pe_slots2 = max(1, round(a.slots * pe_work2 / ((pe_work2 + long_work) or 1)))
    fixed = ctd.pe_threads(pe, pe_slots2)
    ctd.pe_point_seconds = original

    print(f"slots {a.slots}; {a.genome_seconds} s per genome and sample, {a.pair_seconds * 1e6:g} us per pair; long-read work {long_work / 3600:.1f} h "
          f"(the collector's estimate)")
    print(f"pe work: estimate now {pe_work / 3600:.1f} h -> {pe_slots} slots; with genomes {pe_work2 / 3600:.1f} h -> "
          f"{pe_slots2} slots")
    print(f"{'point':32} {'samples':>7} {'genomes':>8} {'est. now h':>10} {'with genomes h':>14} "
          f"{'threads now':>11} {'wall now h':>10} {'threads fixed':>13} {'wall fixed h':>12}")
    for _, pt in sorted(pe, key=lambda x: -true_cost(x[1])):
        n = pt["name"]
        print(f"{n:32} {pt['samples']:7d} {genomes(pt):8.0f} {original(pt, 1) / 3600:10.2f} {true_cost(pt) / 3600:14.2f} "
              f"{now[n]:11d} {wall(pt, now[n]) / 3600:10.2f} {fixed[n]:13d} {wall(pt, fixed[n]) / 3600:12.2f}")


if __name__ == "__main__":
    main()
