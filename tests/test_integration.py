#!/usr/bin/env python3
import cmath
import copy
import csv
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("independent", ROOT / "verify_certificate.py")
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
BINARY = Path(sys.argv[1]).resolve()


def rows(path):
    with open(path) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


with tempfile.TemporaryDirectory(prefix="ame-pt-tests-") as temp:
    base = Path(temp)

    def run(name, args, expected=None):
        directory = base / name
        proc = subprocess.run([str(BINARY), *args, "--output", str(directory)],
                              text=True, capture_output=True, timeout=90)
        assert proc.returncode in ([expected] if expected is not None else [0, 2]), proc.stderr + proc.stdout
        cert = json.loads((directory / "certificate.json").read_text())
        result = verifier.verify(cert)
        assert result["ame"] == (proc.returncode == 0)
        # Every factor histogram accounts for all cuts of its size.
        hist = rows(directory / "rank_histogram.tsv")
        for k in range(1, cert["N"] // 2 + 1):
            for factor in range(len(cert["fields"])):
                counts = sum(int(r["count"]) for r in hist if int(r["k"]) == k and int(r["factor"]) == factor)
                assert counts == math.comb(cert["N"], k)
        for row in rows(directory / "entropy_bips.tsv"):
            assert int(row["rank_sum"]) <= int(row["max_rank_sum"])
        timings = cert["timings_seconds"]
        assert math.isclose(sum(timings[k] for k in ("setup", "search", "verification", "data_export")),
                            timings["total_through_data_export"], abs_tol=1e-8)
        return cert, directory

    args = ["-N", "4", "-d", "2", "--replicas", "4", "--steps", "240",
            "--restarts", "2", "--swap-interval", "5", "--log-interval", "7",
            "--keep-going", "--debug-checks", "--seed", "42"]
    first, a = run("threaded-a", args, 2)
    second, b = run("threaded-b", args, 2)
    first.pop("timings_seconds")
    second.pop("timings_seconds")
    assert first == second, "Fixed-seed multithreaded certificate differs"
    for filename in ("convergence.tsv", "swaps.tsv", "restart_stats.tsv"):
        assert (a / filename).read_bytes() == (b / filename).read_bytes(), filename
    swaps = rows(a / "swaps.tsv")
    assert {(int(r["i"]), int(r["j"])) for r in swaps} == {(0, 1), (1, 2), (2, 3)}
    assert any(int(r["accepted"]) for r in swaps)
    stats = rows(a / "restart_stats.tsv")
    assert len(stats) == 8 and sum(int(r["uphill"]) for r in stats) > 0
    assert first["proposals"] == 4 * 240 * 2
    assert len(swaps) == 144
    for r in swaps:
        i, j = int(r["i"]), int(r["j"])
        ti = .25 * (50 / .25) ** (i / 3)
        tj = .25 * (50 / .25) ** (j / 3)
        want = (1 / ti - 1 / tj) * (int(r["energy_i"]) - int(r["energy_j"]))
        assert math.isclose(float(r["log_acceptance"]), want, abs_tol=1e-10)
    for replicas in [1, 2, 3, 5]:
        cert, folder = run("replicas-" + str(replicas),
                          ["-N", "4", "-d", "2", "--replicas", str(replicas),
                           "--steps", "20", "--restarts", "1", "--swap-interval", "2",
                           "--keep-going", "--debug-checks", "--guide", "0"], 2)
        pairs = {(int(r["i"]), int(r["j"])) for r in rows(folder / "swaps.tsv")}
        assert pairs == {(i, i + 1) for i in range(replicas - 1)}

    for dimension in [4, 6, 12, 121, 257, 343, 512, 625]:
        cert, folder = run("dimension-" + str(dimension),
                          ["-N", "2", "-d", str(dimension), "--replicas", "2",
                           "--steps", "20", "--restarts", "1", "--seed", "71",
                           "--debug-checks", "--diagonal"], 0)
        assert len(rows(folder / "restart_stats.tsv")) == 2
        squarefree = all(f["m"] == 1 for f in cert["fields"])
        assert (folder / "crt_matrix.txt").exists() == squarefree
    for dimension in [4, 9, 121]:
        run("extension-cuts-" + str(dimension),
            ["-N", "4", "-d", str(dimension), "--replicas", "3",
             "--steps", "40", "--restarts", "1", "--keep-going",
             "--swap-interval", "4", "--seed", "7", "--debug-checks"])
    zero = base / "zero.txt"
    zero.write_text("0 0\n0 0\n")
    cert, _ = run("noname-input", ["-N", "2", "-d", "2", "--input", str(zero), "--verify-only"], 2)
    assert cert["timings_seconds"]["search"] == 0
    forged = copy.deepcopy(cert)
    forged["ame"] = True
    try:
        verifier.verify(forged)
        raise AssertionError("Forged positive verdict accepted")
    except ValueError:
        pass
    # The GF(4) nonzero element labelled 2 is a full-rank edge.
    bell4 = base / "bell4.txt"
    bell4.write_text("0 2\n2 0\n")
    cert, folder = run("bell4", ["-N", "2", "-d", "4", "--input", str(bell4), "--verify-only"], 0)
    f = verifier.Field(2, 2, cert["fields"][0]["modulus"])
    def trace(x):
        return f.sub(x, f.sub(0, f.power(x, 2)))
    amp = [[cmath.exp(2j * math.pi * trace(f.mul(2, f.mul(x, y))) / 2) / 4
            for y in range(4)] for x in range(4)]
    rho = [[sum(amp[x][y] * amp[z][y].conjugate() for y in range(4)) for z in range(4)] for x in range(4)]
    purity = sum(abs(v)**2 for row in rho for v in row)
    assert math.isclose(purity, .25, abs_tol=1e-12)
    assert math.isclose(float(rows(folder / "entropy_bips.tsv")[0]["entropy_bits"]), 2)
    # Validate failures before expensive work and preserve previous outputs.
    for i, invalid in enumerate([
        ["--log-interval", "0"], ["--log-interval"], ["-N", "64"], ["-d", "1"],
        ["--replicas", "0"], ["--steps", "0"], ["--tmin", "nan"], ["--guide", "1"],
        ["-d", "1031"], ["-N", "24", "--memory-mb", "1"], ["--unknown"],
    ]):
        path = base / ("invalid-" + str(i))
        proc = subprocess.run([str(BINARY), *invalid, "--output", str(path)], capture_output=True, timeout=20)
        assert proc.returncode == 1 and not path.exists(), (invalid, proc.stderr)
    before = (a / "certificate.json").read_bytes()
    proc = subprocess.run([str(BINARY), "--output", str(a)], capture_output=True, timeout=20)
    assert proc.returncode == 1 and before == (a / "certificate.json").read_bytes()
    print("PASS: threaded reproducibility; odd/even ladders; accepted swaps/uphill moves; independent certificates; extension fields; entropy amplitudes; invalid inputs; output preservation")
