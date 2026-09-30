#!/usr/bin/env python3
"""Reproduce both exact sector minima and independent checks without a search.

Requires Python >= 3.8 and a C++17 compiler. No third-party Python modules.
The bundled classification records are sufficient; downloading the full source
database is optional and can be checked with extract_database.py.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def run(command, log, allowed=(0,)):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True)
    log.write_text(result.stdout + result.stderr)
    if result.returncode not in allowed:
        raise RuntimeError(f"Command failed ({result.returncode}): {command}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New directory; never overwrites a prior run")
    parser.add_argument("--cxx", default="c++")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    database = HERE / "source_records_length8.tsv"
    provenance = json.loads((HERE / "input_provenance.json").read_text())
    if hashlib.sha256(database.read_bytes()).hexdigest() != provenance["extracted_records_sha256"]:
        raise ValueError("Classification input checksum mismatch")
    run([sys.executable, HERE / "analyze_classification.py", database, output], output / "python.log")
    print("Python: all 817 ternary classes evaluated", flush=True)
    with tempfile.TemporaryDirectory(prefix="afame-exact-") as temp:
        build = Path(temp)
        run([args.cxx, "-O3", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
             HERE / "verify_stabilizer_support.cpp", "-o", build / "stabilizer"], output / "stabilizer_build.log")
        run([build / "stabilizer", database, output / "stabilizer_cut_ranks.csv"], output / "stabilizer.log")
        run([sys.executable, HERE / "compare_independent_checks.py", output, HERE / "composite_certificate.json"], output / "comparison.log")
        print("C++ support counting agrees on every one of 132,354 cut ranks", flush=True)
        for method, name in [("elimination", "AfAME_binary8_exact_enumeration.cpp"), ("span", "AfAME_binary8_verify_by_span.cpp")]:
            exe = build / method
            source = HERE / ("binary_eight_party_exhaustive.cpp" if method == "elimination" else "binary_eight_party_span_check.cpp")
            run([args.cxx, "-O3", "-std=c++17", source, "-o", exe], output / f"binary_{method}_build.log")
            result = run([exe, HERE / "binary_graph_atlas_rows.txt"], output / f"binary_{method}.log")
            if "atlas7=1044 extensions=133632" not in result or "min_cost=28" not in result:
                raise ValueError(f"Unexpected binary {method} result")
        run([args.cxx, "-O3", "-std=c++17", "-pthread", HERE / "afame_search_verifier.cpp", "-o", build / "afame"], output / "afame_build.log")
        run([build / "afame", "-N", "8", "-d", "3", "--input", output / "optimal_8_3_matrix.txt", "--verify-only", "--output", output / "native_witness"], output / "afame_witness.log", allowed=(2,))
        native = json.loads((output / "native_witness/certificate.json").read_text())
        if native["verified_cost"] != 16 or native["failing_cuts"] != 16:
            raise ValueError("Native AfAME witness check failed")
        # Negative controls make sure corrupted input/output is rejected.
        corrupt = output / "negative_control.tsv"
        lines = database.read_text().splitlines()
        parts = lines[0].split("\t")
        coefficients = parts[3][1:-1].split(",")
        coefficients[0] = "2"
        parts[3] = "{" + ",".join(coefficients) + "}"
        lines[0] = "\t".join(parts)
        corrupt.write_text("\n".join(lines) + "\n")
        result = subprocess.run([str(build / "stabilizer"), str(corrupt), str(output / "negative_control.csv")], capture_output=True, text=True)
        if result.returncode != 1 or "weight enumerator mismatch" not in result.stderr:
            raise ValueError("Corruption negative control was not rejected")
        (output / "negative_control.log").write_text(result.stderr)
        corrupt.unlink()
        (output / "negative_control.csv").unlink()
    summary = {"binary_minimum": 28, "ternary_minimum": 16, "product_minimum": 44,
               "matched_ternary_cut_ranks": 132354, "native_afame_witness_passed": True,
               "corrupt_input_rejected": True,
               "classification_coverage": "Inherited from the published Danielsen classification and Read-Wilson graph atlas; not regenerated here"}
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
