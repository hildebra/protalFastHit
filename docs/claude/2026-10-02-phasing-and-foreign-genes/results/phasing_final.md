# Phasing benchmark (phasing_score.py)

Per run: the mixture samples' strains (12 species x 8 samples: 4 at 70:30, 2 at 50:50, one 50:30:20, one 85:15). Shares of the strains; SNPs summed over the strains' matched rows.

| run | strains | recovered | resolved (mean) | wrong calls (mean) | nearest pure sample is the strain's | kept by qcmsa | SNP precision | SNP recall | pure samples given strain rows |
|---|---|---|---|---|---|---|---|---|---|
| ont.nophase | 204 | 0.09 | 0.166 | 0.078 | 0.93 | 0.41 | 0.9127 | 0.2682 | 0 |
| ont.phase | 204 | 0.27 | 0.386 | 0.097 | 0.95 | 0.75 | 0.9279 | 0.3984 | 0 |
| pb.nophase | 204 | 0.09 | 0.126 | 0.064 | 0.94 | 0.41 | 0.9482 | 0.2590 | 0 |
| pb.phase | 204 | 0.52 | 0.552 | 0.050 | 0.96 | 0.86 | 0.9734 | 0.5408 | 0 |

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
| ont.phase | 2 strains | 0.85 | 0.25 (4) | 0.00 (3) | 1.00 (5) |
| ont.phase | 2 strains | 0.7 | 0.00 (16) | 0.25 (12) | 0.95 (20) |
| ont.phase | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.20 (20) |
| ont.phase | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.85 (20) |
| ont.phase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
| ont.phase | 3 strains | 0.5 | 0.00 (4) | 0.33 (3) | 0.60 (5) |
| ont.phase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
| ont.phase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
| pb.nophase | 2 strains | 0.85 | 1.00 (4) | 1.00 (3) | 0.80 (5) |
| pb.nophase | 2 strains | 0.7 | 0.25 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.00 (20) |
| pb.nophase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.5 | 0.50 (4) | 0.33 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.3 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.nophase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.00 (5) |
| pb.phase | 2 strains | 0.85 | 0.75 (4) | 0.33 (3) | 1.00 (5) |
| pb.phase | 2 strains | 0.7 | 0.38 (16) | 0.92 (12) | 1.00 (20) |
| pb.phase | 2 strains | 0.5 | 0.00 (16) | 0.33 (12) | 0.40 (20) |
| pb.phase | 2 strains | 0.3 | 0.06 (16) | 0.92 (12) | 1.00 (20) |
| pb.phase | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 1.00 (5) |
| pb.phase | 3 strains | 0.5 | 0.50 (4) | 1.00 (3) | 0.40 (5) |
| pb.phase | 3 strains | 0.3 | 0.00 (4) | 0.33 (3) | 0.60 (5) |
| pb.phase | 3 strains | 0.2 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
