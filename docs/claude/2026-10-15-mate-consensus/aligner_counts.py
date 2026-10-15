#!/usr/bin/env python3
"""The aligner's counts of mate guidance and of the mates' vote, summed over the samples of each protal log given
(one line per sample in the log: "N fragments had one mate sure ...", "Mate consensus: ...", "Gene neighbours: N
fragments paired across ...").

usage: aligner_counts.py LABEL:protal.log [...]     prints a table, a row per label
Standard library only.
"""
import re
import sys

GUIDE = re.compile(r"^(\d+) fragments had one mate sure of its alignment and the other without a candidate of its taxon: "
                   r"(\d+) aligned from its own anchor of that taxon, (\d+) of (\d+) found on the mate's gene")
VOTE = re.compile(r"^Mate consensus: (\d+) fragments whose best candidate had one mate were voted on \((\d+) with their "
                  r"mates' best alignments on two taxa\): (\d+) went to another taxon than their best candidate's, (\d+) got "
                  r"a MAPQ under \d+; (\d+) of (\d+) mates looked for on a taxon were found \((\d+) from their own anchor, "
                  r"(\d+) on the other mate's gene\); (\d+) mates kept an alignment of another taxon")
ACROSS = re.compile(r"^Gene neighbours: (\d+) fragments paired across two neighbouring genes; (\d+) guided mates found")
COLUMNS = ["guided", "guide_from_anchor", "guide_found", "guide_looked_for", "voted", "vote_disagreed", "vote_moved",
           "vote_mapq_under_4", "vote_found", "vote_looked_for", "vote_from_anchor", "vote_on_gene", "vote_kept",
           "paired_across_genes", "found_on_neighbour"]


def main():
    print("\t".join(["run"] + COLUMNS))
    for spec in sys.argv[1:]:
        label, path = spec.split(":", 1)
        total = dict.fromkeys(COLUMNS, 0)
        with open(path, errors="replace") as fh:
            for line in fh:
                if m := GUIDE.match(line):
                    for key, value in zip(COLUMNS[0:4], (m[1], m[2], m[3], m[4])):
                        total[key] += int(value)
                elif m := VOTE.match(line):
                    for key, value in zip(COLUMNS[4:13], m.groups()):
                        total[key] += int(value)
                elif m := ACROSS.match(line):
                    total["paired_across_genes"] += int(m[1])
                    total["found_on_neighbour"] += int(m[2])
        print("\t".join([label] + [str(total[c]) for c in COLUMNS]))


if __name__ == "__main__":
    main()
