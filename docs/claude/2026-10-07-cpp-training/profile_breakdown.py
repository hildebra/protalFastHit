#!/usr/bin/env python3
"""Where the trainer's time goes, from a cProfile of machine_learning_cmdline.py (python3 -m cProfile -o seq.prof ...):
the callees of train() and of the boosting fit, by cumulative seconds.

    python3 profile_breakdown.py seq.prof
"""

import pstats
import sys


def callees(stats, match, minimum=0.3):
    """The function whose name contains `match`, and its callees as (seconds, ncalls, name), largest first."""
    found = [f for f in stats.stats if match in f"{f[0]}:{f[1]}({f[2]})"]
    if not found:
        return []
    target = found[0]
    rows = []
    for func, (cc, nc, tt, ct, callers) in stats.stats.items():
        if target in callers:
            c_cc, c_nc, c_tt, c_ct = callers[target]
            if c_ct >= minimum:
                rows.append((c_ct, c_nc, f"{func[0].split('site-packages/')[-1]}:{func[1]}({func[2]})"))
    return target, sorted(rows, reverse=True)


def main():
    stats = pstats.Stats(sys.argv[1])
    print(f"total {stats.total_tt:.1f} s")
    for match in ("machine_learning_cmdline.py:1765(train)", "machine_learning_cmdline.py:692(study_evaluation)",
                  "gradient_boosting.py:393(fit)", "grower.py:486(split_next)", "grower.py:402(_initialize_root)"):
        result = callees(stats, match)
        if not result:
            continue
        target, rows = result
        print(f"\n{match}: {stats.stats[target][3]:.1f} s cumulative, {stats.stats[target][2]:.1f} s own")
        for ct, nc, name in rows[:14]:
            print(f"  {ct:8.1f} s  {nc:>7}  {name}")


if __name__ == "__main__":
    main()
