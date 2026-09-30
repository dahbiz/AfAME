#!/usr/bin/env python3
"""Reproduce the frozen sector-recombination benchmark summary.

The expensive 150-run search is represented by the supplied per-seed CSV.
This script validates every row, recomputes the seed-bootstrap intervals, and
writes a machine-readable summary. To regenerate new trajectories, compile the
instrumented C++ source in this folder against the AfAME revision described in
benchmark_protocol.md, then export the same CSV schema.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=here / "recombination_benchmark.csv")
    parser.add_argument("--output", type=Path,
                        default=here / "benchmark_summary.json")
    args = parser.parse_args()

    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"N", "d", "seed", "proposals", "joint_cost", "recombined_cost", "gain"}
    if not rows or not required <= set(rows[0]):
        raise SystemExit(f"missing benchmark columns; need {sorted(required)}")
    if len(rows) != 150:
        raise SystemExit(f"expected 150 held-out rows, found {len(rows)}")

    rng = np.random.default_rng(20260923)
    summary = []
    for n in (10, 11, 12):
        group = [r for r in rows if int(r["N"]) == n]
        if len(group) != 50:
            raise SystemExit(f"N={n}: expected 50 rows, found {len(group)}")
        gains = np.array([int(r["gain"]) for r in group], dtype=int)
        if any(int(r["joint_cost"]) - int(r["recombined_cost"]) != int(r["gain"])
               for r in group):
            raise SystemExit(f"N={n}: gain column is inconsistent")
        if any(int(r["proposals"]) != 2000 for r in group):
            raise SystemExit(f"N={n}: proposal budget is inconsistent")
        bootstrap = rng.choice(gains, size=(50000, 50)).mean(axis=1)
        summary.append({
            "N": n,
            "rows": 50,
            "joint_mean": float(np.mean([int(r["joint_cost"]) for r in group])),
            "recombined_mean": float(np.mean([int(r["recombined_cost"]) for r in group])),
            "mean_gain": float(np.mean(gains)),
            "median_gain": float(np.median(gains)),
            "strict_gains": int(np.count_nonzero(gains > 0)),
            "bootstrap_95_percentile_interval": [float(x) for x in np.quantile(bootstrap, [.025, .975])],
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"source_rows": len(rows), "groups": summary}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
