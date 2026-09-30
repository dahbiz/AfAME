<p align="center">
  <img src="AfAME_logo.png" alt="AfAME logo" width="480">
</p>

# AfAME

Search for and certify absolutely maximally entangled (AME) states using finite-field cut ranks and parallel tempering.

Developed by **Dr. Zakaria Dahbi**.

[![Build and verify](https://github.com/dahbiz/AfAME/actions/workflows/tests.yml/badge.svg)](https://github.com/dahbiz/AfAME/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

AfAME combines a concurrent C++17 search with exact algebraic checks and a standalone Python certificate verifier. It supports prime fields, extension fields, and tensor products of field sectors. Every reported verdict is decided by exhaustive final verification, not by the search heuristic.

- **Parallel search.** Persistent CPU workers run parallel-tempering replicas with alternating adjacent exchanges, following [Earl and Deem, Eq. (2)](https://arxiv.org/pdf/physics/0508111).
- **Exact certification.** All cuts are enumerated with explicit field arithmetic.
- **Reproducibility.** Seeds, field moduli, search parameters, swap logs, and restart records are exported with every run.
- **Independent checking.** A separate, dependency-free Python implementation re-verifies exported certificates.

## Requirements

A C++17 compiler and CPU threads. Python 3 is needed only for the independent verifier and the integration tests. There are no MPI, PETSc, or OpenMP dependencies. The search program and verifier use no third-party Python packages; the figure and benchmark scripts under `data/` additionally use NumPy and Matplotlib.

## Build and run

    git clone https://github.com/dahbiz/AfAME.git
    cd AfAME
    make CXX=clang++
    ./AfAME -N 5 -d 2 --replicas 4 --steps 5000 --restarts 20 --seed 42

The launcher `run.sh` builds automatically if necessary:

    bash run.sh -N 17 -d 10001 --replicas 4 --steps 5000 --seed 42

Run the executable directly; mpirun would launch independent copies, not a distributed version of this algorithm. Each replica runs on its own persistent thread.

Each run creates a new output directory and prints its path. Choose a location with `--output NEW_DIRECTORY`; existing directories are rejected to preserve previous results. The historical `-outprefix` spelling is accepted as an alias.

Exit codes: **0** for a verified AME state, **2** for a verified non-AME best candidate, **1** for an error. A non-AME result is neither a nonexistence proof nor a global-optimality certificate.

The default temperature ladder spans 0.25 to 50 (`--tmin`, `--tmax`). It is fixed within a restart; restarts initialize all replicas afresh, with no hidden cooling or reheating schedule. Defaults are starting values, not an optimized prescription for all dimensions. Swap-acceptance and uphill-move counts are recorded for tuning.

## Command-line options

| Option | Default | Description |
|---|---|---|
| `-N` | 7 | Number of parties (2–24). |
| `-d` | 2 | Local dimension, factored into prime-power sectors automatically. |
| `--replicas` | 4 | Tempering replicas, one thread each (1–64). |
| `--steps` | 5000 | Proposals per replica per restart. |
| `--restarts` | 20 | Full reinitializations of all replicas. |
| `--seed-tries` | 1 | Total restart opportunities equal `--restarts * --seed-tries`. |
| `--swap-interval` | 25 | Adjacent-pair exchange attempts; even and odd pairs alternate. |
| `--log-interval` | 200 | Convergence-log entries are written on this interval. |
| `--tmin`, `--tmax` | 0.25, 50 | Temperature ladder bounds. |
| `--guide` | 0.85 | Guided-proposal mixture weight; 0 selects ordinary uniform Metropolis proposals. |
| `--memory-mb` | 512 | Estimated cut/cache preallocation check; not a hard process-RSS limit. |
| `--seed` | 1 | Base seed. |
| `--output` | new directory | Run directory, created fresh; existing directories are rejected. Legacy alias: `-outprefix`. |
| `--input` | none | Plain matrix file, used with `--verify-only`. |
| `--verify-only` | off | Verify the input matrix and export a certificate; no search. |
| `--keep-going` | off | Disable early stopping, for exchange diagnostics. |
| `--diagonal` | off | Random local diagonal phases at initialization; fixed during search. |
| `--debug-checks` | off | Extra internal invariant checks. |
| `-print_search_entropy` | off | Print the location of the entropy output files. |
| `-tensor_phases [1]` | | Legacy flag, accepted and ignored; factorization is automatic. |

Bounds enforced by the program: 2 ≤ N ≤ 24, d ≥ 2, field size q ≤ 1024 (16-bit elements), 1 ≤ replicas ≤ 64, 1e-6 ≤ tmin ≤ tmax ≤ 1e9, 0 ≤ guide < 1, and restarts × seed-tries ≤ 1,000,000. Every extension modulus is generated deterministically and tested for irreducibility; every nonzero inverse is checked at setup.

Guided proposals mix uniform edges with weights proportional to the number of failing cuts, and the reverse/forward proposal ratio is included in the Metropolis-Hastings acceptance.

## Certification and outputs

Final exhaustive verification determines every result label. Certification enumerates all nonempty subsets of size at most floor(N/2), including both halves of balanced cuts. This is exponential in N; the program limits N to 24 and checks an estimated cut/cache memory budget before allocation (`--memory-mb`, default 512). The estimate is not a hard process-RSS limit.

Outputs include `certificate.json`, `run_summary.txt`, factor matrices, per-cut entropies, per-size statistics, per-factor rank histograms, edge failure counts, convergence logs, restart records, and every attempted swap. Statistics, rank histograms, and seed records are corrected counts. Timing includes data export but excludes writing the final certificate and summary metadata.

Early stopping is synchronized at dispatch boundaries. `--steps` counts proposals per replica per restart. The logged winning-slot seed is not sufficient by itself to replay an exchanged trajectory; use the base seed with the full replica and run parameters. Repeated runs were deterministic with the same executable and parameters; cross-toolchain floating-point identity is not promised.

## Independent verification

    python3 verify_certificate.py PATH_TO_RUN/certificate.json --progress

The verifier is a separate implementation with its own polynomial arithmetic, trial-division irreducibility checks, and independent elimination. It checks dimensions, symmetry, factorization, cut counts, cost, and the claimed verdict. For extension fields it accepts the explicitly specified valid modulus in its JSON input.

Recheck the manuscript's matrix:

    ./AfAME -N 17 -d 10001 \
      --input examples/paper_17_10001.txt --verify-only --output paper_recheck
    python3 verify_certificate.py paper_recheck/certificate.json --progress

Plain matrix input is supported for one field or for square-free CRT dimensions. The matrix must contain exactly N*N integers; `#` comments are allowed. For a single extension field, labels use this program's generated modulus.

## Included example

The fresh 17-party example can be reproduced with:

    ./AfAME -N 17 -d 10001 --replicas 4 --steps 500 \
      --restarts 1 --seed 42 --swap-interval 5 --log-interval 50

It reached zero cost by the 60-step dispatch boundary in the validation run. Its certificate and independent verification are in `examples/search_17_10001/`. This observation does not establish comparative performance.

## Scientific conventions

For each field GF(p^m), the state phase uses the absolute trace of sum(i ≤ j) P[i,j] x[i] x[j], with character exp(2·π·i·Tr(···)/p). Labels are polynomial-basis coefficients encoded in base p. The exact modulus is included in every certificate and factor-matrix header.

The search runs over the tensor product of the field sectors of the distinct prime-power factors of d. Its objective is the sum over cuts and factors of squared rank deficits; entropy is sum_f rank_f·log(q_f).

**GF(p^m) and Z/(p^m) are distinct models when m > 1.** CRT residue matrices are exported only for square-free d. Their sector character weights can differ from those of the directly defined tensor state; their cut ranks and entropies agree. For nonsquarefree d, use the factor matrices and JSON certificate.

The generated modulus may differ from the original program's table entries. Do not reinterpret an old extension-field matrix's integer labels under a new modulus.

## Documentation

- [VALIDATION.md](VALIDATION.md) — numerical, threading, and sanitizer checks.
- [data/README.md](data/README.md) — reproducibility data and figure scripts.
- [AUTHORS.md](AUTHORS.md) — development credits.
- [CONTRIBUTING.md](CONTRIBUTING.md) — report and change guidelines.
- `examples/` — recorded runs with certificates and logs.

## Citation

Please credit **Dr. Zakaria Dahbi** when referring to AfAME. GitHub can display a software citation using [CITATION.cff](CITATION.cff).

## License

AfAME is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](LICENSE).
