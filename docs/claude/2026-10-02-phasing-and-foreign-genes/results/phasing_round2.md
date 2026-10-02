# Phasing benchmark (phasing_score.py)

Per run: the mixture samples' strains (12 species x 8 samples: 4 at 70:30, 2 at 50:50, one 50:30:20, one 85:15). Shares of the strains; SNPs summed over the strains' matched rows.

| run | strains | recovered | resolved (mean) | wrong calls (mean) | nearest pure sample is the strain's | kept by qcmsa | SNP precision | SNP recall | pure samples given strain rows |
|---|---|---|---|---|---|---|---|---|---|
| ont.nophase | 204 | 0.09 | 0.166 | 0.078 | 0.93 | 0.41 | 0.9127 | 0.2682 | 0 |
| ont.phase | 204 | 0.25 | 0.410 | 0.114 | 0.92 | 0.59 | 0.8623 | 0.4656 | 2 |
| pb.nophase | 204 | 0.09 | 0.126 | 0.064 | 0.94 | 0.41 | 0.9482 | 0.2590 | 0 |
| pb.phase | 204 | 0.45 | 0.550 | 0.086 | 0.97 | 0.68 | 0.9129 | 0.5979 | 0 |

By design and coverage (strains recovered):

| run | design | strain share | coverage < 8x | 8-20x | > 20x |
|---|---|---|---|---|---|
| ont.nophase | 2 strains | 0.85 | 0.50 (4) | 1.00 (3) | 1.00 (5) |
| ont.nophase | 2 strains | 0.7 | 0.06 (16) | 0.42 (12) | 0.00 (20) |
| ont.nophase | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| ont.nophase | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| ont.nophase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| ont.nophase | 3 strains | 0.5 | 0.25 (4) | 0.67 (3) | 0.00 (5) |
| ont.nophase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| ont.nophase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| ont.phase | 2 strains | 0.85 | 0.50 (4) | 1.00 (3) | 1.00 (5) |
| ont.phase | 2 strains | 0.7 | 0.12 (16) | 0.50 (12) | 0.95 (20) |
| ont.phase | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.20 (20) |
| ont.phase | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.15 (20) |
| ont.phase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| ont.phase | 3 strains | 0.5 | 0.25 (4) | 0.67 (3) | 0.80 (5) |
| ont.phase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
| ont.phase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.nophase | 2 strains | 0.85 | 1.00 (4) | 1.00 (3) | 0.80 (5) |
| pb.nophase | 2 strains | 0.7 | 0.25 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.5 | 0.50 (4) | 0.33 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.phase | 2 strains | 0.85 | 1.00 (4) | 0.67 (3) | 1.00 (5) |
| pb.phase | 2 strains | 0.7 | 0.50 (16) | 0.83 (12) | 1.00 (20) |
| pb.phase | 2 strains | 0.5 | 0.00 (16) | 0.17 (12) | 0.40 (20) |
| pb.phase | 2 strains | 0.3 | 0.00 (16) | 0.08 (12) | 0.80 (20) |
| pb.phase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.phase | 3 strains | 0.5 | 0.50 (4) | 1.00 (3) | 1.00 (5) |
| pb.phase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 1.00 (5) |
| pb.phase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
