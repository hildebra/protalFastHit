#!/usr/bin/env python3
"""The read templates of long-read samples, as the collector draws them (collect_training_data.long_read_templates),
for the communities of a build's paired-end test samples: one FASTA per sample, to make reads of the same templates
with pbsim3 and with hifi_reads.py.

usage: make_templates.py SCRIPTS_DIR BUILD_OUT PE_POINT BASES OUT_DIR [SAMPLES [SETUP]]   (SETUP: hifi:15000:3000:3)
  e.g. make_templates.py ~/bprof/gnb/src/scripts ~/opw/b_gn2 rl100_p100000 30000000 ~/hifi/templates 2
"""

import os
import sys


def main():
    scripts, build, point, bases, out = sys.argv[1:6]
    samples = int(sys.argv[6]) if len(sys.argv) > 6 else 2
    setup_text = sys.argv[7] if len(sys.argv) > 7 else "hifi:15000:3000:3"
    sys.path.insert(0, scripts)
    import collect_training_data as collect
    os.makedirs(out, exist_ok=True)
    community_dir = os.path.join(build, "test", "points", point)
    by_sample = {}
    for row in collect.manifest_rows(community_dir):
        by_sample.setdefault(row["sample"], []).append(row)
    setup = collect.parse_long_setup(setup_text)
    for s, (community, genomes) in enumerate(sorted(by_sample.items())[:samples]):
        task = {"seed": 1000 + s, "setup": setup, "bases": int(float(bases)),
                "templates": os.path.join(out, f"{community}.fa"),
                "genomes": [{"genome": g["genome"], "fasta": g["fasta_path"],
                             "weight": float(g["relative_abundance"]) * float(g["genome_length"])} for g in genomes]}
        names, error = collect.long_read_templates(task)
        if error:
            sys.exit(error)
        print(f"{task['templates']}: {len(names)} templates of {len(genomes)} genomes")


if __name__ == "__main__":
    main()
