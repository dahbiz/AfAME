#!/usr/bin/env python3
"""Compare every independently obtained cut rank and check the paper witness."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from finite_field_certificate import Field, verify


def read_ranks(path):
    rows = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = int(row["class"]), int(row["mask"])
            if key in rows:
                raise ValueError(f"Duplicate cut {key}")
            rows[key] = tuple(int(row[column]) for column in ["size", "rank", "deficit"])
    if len(rows) != 817 * 162:
        raise ValueError("Incomplete cut-rank output")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("composite_certificate", type=Path)
    args = parser.parse_args()
    first = read_ranks(args.results / "python_cut_ranks.csv")
    second = read_ranks(args.results / "stabilizer_cut_ranks.csv")
    if first != second:
        raise ValueError("Independent calculations disagree")
    summary = json.loads((args.results / "summary.json").read_text())
    if summary["minimum_cost"] != 16:
        raise ValueError("Unexpected ternary minimum")
    data = json.loads(args.composite_certificate.read_text())
    if (data["N"], data["d"]) != (8, 6):
        raise ValueError("Expected the eight-party dimension-six paper certificate")
    verified = verify(data)
    import itertools
    sector_costs, sector_histograms = [], []
    for item in data["fields"]:
        field = Field(item["p"], item["m"], item["modulus"])
        matrix = item["matrix"]
        cost, histogram = 0, {}
        for k in range(1, 5):
            for subset in itertools.combinations(range(8), k):
                complement = [j for j in range(8) if j not in subset]
                deficit = k - field.rank([[matrix[i][j] for j in complement] for i in subset])
                cost += deficit**2
                key = f"{k},{deficit}"
                histogram[key] = histogram.get(key, 0) + 1
        sector_costs.append(cost)
        sector_histograms.append({"q": field.q, "histogram_size_deficit": histogram})
    if sector_costs != [28, 16] or verified["verified_cost"] != 44:
        raise ValueError("Paper certificate fails to attain the sector minima")
    report = {"independent_cut_ranks_match": True, "matched_cut_ranks": len(first),
              "paper_certificate": verified, "paper_sector_costs": sector_costs,
              "paper_sector_histograms": sector_histograms,
              "paper_certificate_sha256": hashlib.sha256(args.composite_certificate.read_bytes()).hexdigest(),
              "scope": "C*_3(8)=16 from the complete published classification; combined with the separately verified binary minimum 28, C*_6(8)=44 in the product family. No unrestricted AME nonexistence conclusion."}
    (args.results / "independent_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
