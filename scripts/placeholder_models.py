#!/usr/bin/env python3
"""Write placeholder presence models for read types that have no trained model yet.

A database holds one model per read type (model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml),
and protal stops with an error when a run's read type has none. A placeholder fills the slot
without pretending to be a model: it scores every taxon 0, so protal reports no species (--knob 0
lists every taxon with reads, e.g. to test profiling SAMs of that read type), and protal warns
whenever it loads one. build_gtdb_database.py puts placeholders for se, pb and ont into the database
it builds; replace one with a trained model by

    protal --add_model MODEL --read_type se --db DB

Written to OUT/: model_se.xml, model_PB.xml and model_ONT.xml, or those of --read_types. Store one in
an existing database with the command above (it replaces the model stored for that read type).

  python3 scripts/placeholder_models.py -o models
  python3 scripts/placeholder_models.py -o models --read_types pe,se,pb,ont
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_pmml import MODEL_FILES, write_placeholder  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-o", "--out", default=".", help="folder to write the models to (default: .)")
    p.add_argument("--read_types", default="se,pb,ont",
                   help="comma-separated read types: pe, se, pb, ont (default: se,pb,ont)")
    opts = p.parse_args(argv)
    types = [t.strip() for t in opts.read_types.split(",") if t.strip()]
    unknown = [t for t in types if t not in MODEL_FILES]
    if unknown:
        sys.exit(f"unknown read type(s) {', '.join(unknown)}; protal knows {', '.join(MODEL_FILES)}")
    os.makedirs(opts.out, exist_ok=True)
    for read_type in types:
        path = os.path.join(opts.out, MODEL_FILES[read_type])
        write_placeholder(path, read_type)
        print(f"{path}: placeholder for read type {read_type} "
              f"(store: protal --add_model {path} --read_type {read_type} --db DB)")


if __name__ == "__main__":
    main()
