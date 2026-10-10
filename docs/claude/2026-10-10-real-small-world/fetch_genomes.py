"""download_gtdb.py's genome step alone, on an inputs folder whose release/ already holds the files (make_inputs.sh:
the small world's genomes only, so the species and strains are those): the strains and representatives to simulate
from, from NCBI, genomes.tsv and simulation_species.txt, as download_gtdb.py writes them. The release step (which
would fetch GTDB's whole archives) and the host genome (no scenario here has host reads) are left out.

usage: fetch_genomes.py -o <inputs folder> [download_gtdb.py's genome options: --species, --per_species, ...]
"""
import datetime
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
import download_gtdb  # noqa: E402

if __name__ == "__main__":  # the download workers import this module again (forkserver)
    opts = download_gtdb.parse_args(sys.argv[1:])
    state = download_gtdb.load_state(opts)
    download_gtdb.get_genomes(opts, state, "226")
    state["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    download_gtdb.save_state(opts, state)
    print(f"Inputs: {opts.out} (build_gtdb_database.py --inputs {opts.out})", flush=True)
