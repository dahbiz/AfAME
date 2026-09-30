#!/usr/bin/env python3
# AfAME - parallel-tempering search and algebraic certification of
# absolutely maximally entangled (AME) states.
# Copyright (C) 2026 Zakaria Dahbi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
"""Evaluate all 817 qutrit classes using AfAME's independent Python verifier.

The complete cut-rank vector is invariant under partywise Clifford operations
up to permutation of parties. Completeness is inherited from Danielsen's
published classification, not inferred from stochastic searches.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from finite_field_certificate import Field, verify


def read_records(path):
    records = []
    for index, line in enumerate(path.read_text().splitlines(), 1):
        parts = line.split("\t")
        if len(parts) != 6 or parts[0] != "8" or parts[1] not in ("D", "I"):
            raise ValueError(f"Malformed record {index}")
        weights = [int(x) for x in parts[3].strip("{}").split(",")]
        if len(weights) != 9 or sum(weights) != 3**8 or weights[0] != 1:
            raise ValueError(f"Malformed weight enumerator {index}")
        matrix = [[0] * 8 for _ in range(8)]
        for edge in ([] if parts[5] == "-" else parts[5].split(",")):
            match = re.fullmatch(r"(?:(\d+)\*)?(\d+)-(\d+)", edge)
            if not match:
                raise ValueError(f"Malformed edge {edge}")
            w, i, j = int(match[1] or 1), int(match[2]), int(match[3])
            if not (0 <= i < j < 8 and 1 <= w < 3) or matrix[i][j]:
                raise ValueError(f"Invalid/repeated edge {edge}")
            matrix[i][j] = matrix[j][i] = w
        records.append({"index": index, "kind": parts[1], "distance": int(parts[2]),
                        "weights": weights, "automorphisms": int(parts[4]), "matrix": matrix})
    if len(records) != 817 or Counter(r["kind"] for r in records) != {"I": 659, "D": 158}:
        raise ValueError("Classification count mismatch")
    # Uses published automorphism orders as a transcription/completeness checksum.
    observed_mass = sum(Fraction(24**8 * math.factorial(8), r["automorphisms"]) for r in records)
    expected_mass = math.prod(3**i + 1 for i in range(1, 9))
    if observed_mass != expected_mass:
        raise ValueError("Published mass formula does not match extracted records")
    return records, expected_mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    records, mass = read_records(args.database)
    args.output.mkdir(parents=True, exist_ok=True)
    field = Field(3, 1, [])
    cuts = [(k, s, tuple(i for i in range(8) if i not in s), sum(1 << i for i in s))
            for k in range(1, 5) for s in itertools.combinations(range(8), k)]
    all_results, best = [], None
    with (args.output / "python_cut_ranks.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["class", "mask", "size", "rank", "deficit"])
        for record in records:
            hist, cost, linear, failing = Counter(), 0, 0, 0
            for k, s, complement, mask in cuts:
                rank = field.rank([[record["matrix"][i][j] for j in complement] for i in s])
                deficit = k - rank
                hist[(k, deficit)] += 1
                cost += deficit**2
                linear += deficit
                failing += deficit != 0
                writer.writerow([record["index"], mask, k, rank, deficit])
            record["cost"] = cost
            record["histogram"] = {f"{k},{d}": count for (k, d), count in sorted(hist.items())}
            all_results.append({"class": record["index"], "kind": record["kind"],
                                "distance": record["distance"], "cost": cost,
                                "linear_deficit": linear, "failing_cuts": failing})
            if best is None or cost < best["cost"]:
                best = record
    with (args.output / "class_costs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_results[0]))
        writer.writeheader()
        writer.writerows(all_results)
    minima = [r["index"] for r in records if r["cost"] == best["cost"]]
    optimal_records = [r for r in records if r["cost"] == best["cost"]]
    (args.output / "optimal_classes.json").write_text(json.dumps(optimal_records, indent=2) + "\n")
    profiles = Counter(tuple(r["histogram"].get(f"4,{d}", 0) for d in range(5)) for r in optimal_records)
    certificate = {"schema": 1, "N": 8, "d": 3, "ame": False,
                   "verified_cost": best["cost"], "cuts": len(cuts),
                   "failing_cuts": sum(v for k, v in best["histogram"].items() if int(k.split(",")[1]) > 0),
                   "fields": [{"p": 3, "m": 1, "q": 3, "modulus": [], "matrix": best["matrix"]}],
                   "provenance": {"classification_record": best["index"], "method": "exhaustive classification evaluation"}}
    verify(certificate)
    (args.output / "optimal_8_3_certificate.json").write_text(json.dumps(certificate, indent=2) + "\n")
    (args.output / "optimal_8_3_matrix.txt").write_text("\n".join(" ".join(map(str, row)) for row in best["matrix"]) + "\n")
    summary = {"N": 8, "q": 3, "classes": len(records), "indecomposable": 659,
               "decomposable": 158, "cuts_per_class": len(cuts), "cut_evaluations": len(cuts)*len(records),
               "mass_formula": str(mass), "minimum_cost": best["cost"], "minimizing_classes": minima,
               "optimal_balanced_profiles": [{"counts_by_deficit": list(profile), "classes": count}
                                             for profile, count in sorted(profiles.items())],
               "witness_histogram_size_deficit": best["histogram"],
               "minimum_by_distance": {str(d): min(r["cost"] for r in records if r["distance"] == d) for d in range(1, 5)},
               "minimum_by_kind": {kind: min(r["cost"] for r in records if r["kind"] == kind) for kind in ["D", "I"]},
               "cost_histogram": dict(sorted(Counter(r["cost"] for r in records).items())),
               "input_sha256": hashlib.sha256(args.database.read_bytes()).hexdigest()}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "cost_histogram"}, indent=2))


if __name__ == "__main__":
    main()
