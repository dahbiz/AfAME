# AfAME

**Developed by Dr. Zakaria Dahbi**

Search for and certify absolutely maximally entangled (AME) states using finite-field cut ranks and parallel tempering.

AfAME combines concurrent C++17 search replicas with exact algebraic checks and a standalone Python certificate verifier. It supports prime fields, extension fields, and tensor products of field sectors.

- **Parallel search:** persistent CPU workers and alternating replica exchanges.
- **Exact certification:** full cut enumeration with explicit field arithmetic.
- **Reproducible results:** seeds, field moduli, search parameters, and swap logs.
- **Independent checks:** a separate Python implementation verifies exported certificates.

Start with the commands below. See [validation results](VALIDATION.md), [developer credits](AUTHORS.md), and [citation metadata](CITATION.cff).

## Build and run

Requires a C++17 compiler and CPU threads. Python 3 is needed only for the independent verifier and integration tests. There are no MPI, PETSc, OpenMP, or third-party Python dependencies.

The repository contains portable source. Clone it, then build the executable on your platform or use the automatic launcher.

    git clone https://github.com/dahbiz/AfAME.git
    cd AfAME

From this directory:

    make CXX=clang++
    ./AME7DIA_PT -N 5 -d 2 --replicas 4 --steps 5000 --restarts 20 --seed 42

Or use the launcher, which builds automatically if necessary:

    bash run.sh -N 17 -d 10001 --replicas 4 --steps 5000 --seed 42

The default temperature ladder spans 0.25 to 50. Tune it with --tmin and --tmax. Each replica runs on its own persistent CPU thread. Run the executable directly: mpirun would launch independent copies, not a distributed version of this algorithm.

Each run creates a new directory and prints its path. Choose a location with --output NEW_DIRECTORY; existing directories are rejected to preserve previous results. The historical -outprefix spelling is accepted as an alias for the output directory, not a filename prefix.

Exit codes are **0** for a verified AME state, **2** for a verified non-AME best candidate, and **1** for an error. A non-AME result is neither a nonexistence proof nor a global-optimality certificate.

## What changed

- Field arithmetic uses 16-bit elements, with an enforced maximum field size of 1024. Every extension modulus is generated deterministically and tested for irreducibility; every nonzero inverse is checked. The invalid polynomial table and null-pointer fallback are removed.
- Each accepted move updates the failure counts of every affected edge. Rejected moves restore ranks, costs, and guidance exactly. Immutable cut geometry is shared across replicas.
- Actual parallel tempering replaces the serial, nearly zero-temperature search. Adjacent exchanges alternate (0,1),(2,3),... and (1,2),(3,4),...; a two-replica run attempts its sole pair every round.
- The acceptance exponent is (1/Ti - 1/Tj)*(Ei - Ej). Configurations, caches, and walker identities move together. Temperatures and random generators remain at their slots.
- Guided proposals use a mixture of uniform edges and weights proportional to the number of failing cuts. Their reverse/forward proposal ratio is included in Metropolis-Hastings acceptance. --guide 0 gives ordinary uniform Metropolis proposals. The old biased tournament selection and greedy descent are not retained.
- Final exhaustive verification determines every result label. Statistics, factor rank histograms, seed records, and successful-restart records are corrected. Runtime is divided into setup, search, verification, and data export.

The temperature ladder is fixed within a restart. Restarts initialize all replicas afresh; there is no undocumented cooling or reheating schedule. Defaults are starting values, not an optimized prescription for all dimensions. Swap acceptance and uphill counts are recorded for tuning. The exchange rule follows [Earl and Deem, Eq. (2)](https://arxiv.org/pdf/physics/0508111).

The included fresh 17-party example can be reproduced with:

    ./AME7DIA_PT -N 17 -d 10001 --replicas 4 --steps 500 \
      --restarts 1 --seed 42 --swap-interval 5 --log-interval 50

It reached zero cost by the 60-step dispatch boundary in the validation run. Its certificate and independent verification are in examples/search_17_10001. This observation does not establish comparative performance.

## Scientific conventions

For each field GF(p^m), the state phase uses the absolute trace of sum(i <= j) P[i,j] x[i] x[j], with character exp(2*pi*i*Tr(...)/p).
Labels are polynomial-basis coefficients encoded in base p. The exact modulus is included in every certificate and factor-matrix header.

The program searches a tensor product of field sectors for the distinct prime-power factors of d. Its objective is the sum over cuts and factors of squared rank deficits. Entropy is sum_f rank_f*log(q_f).

**GF(p^m) and Z/(p^m) are distinct models when m > 1.** CRT residue matrices are exported only for square-free d. Their sector character weights can differ from those of the directly defined tensor state; their cut ranks and entropies agree. For nonsquarefree d, use the factor matrices and JSON certificate.

The generated modulus may differ from the original program's valid table entries. Do not reinterpret an old extension-field matrix's integer labels under a new modulus. The independent verifier accepts the explicitly specified valid modulus in its JSON input.

--diagonal randomizes local diagonal phases at initialization. They remain fixed during search because they do not affect any cut rank. No computation is spent optimizing them.

## Verification and outputs

Verify a generated certificate independently:

    python3 verify_certificate.py PATH_TO_RUN/certificate.json --progress

The Python verifier uses separate polynomial arithmetic, trial-division irreducibility checks, and independent elimination. It checks dimensions, symmetry, factorization, cut counts, cost, and the claimed verdict.

Recheck the manuscript's matrix:

    ./AME7DIA_PT -N 17 -d 10001 \
      --input examples/paper_17_10001.txt --verify-only --output paper_recheck
    python3 verify_certificate.py paper_recheck/certificate.json --progress

Plain matrix input is supported for one field or for square-free CRT dimensions. The matrix must contain exactly N*N integers; # comments are allowed. For a single extension field, labels use this program's generated modulus.

Outputs include certificate.json, run_summary.txt, factor matrices, per-cut entropies, per-size statistics, per-factor rank histograms, edge failure counts, convergence logs, restart records, and every attempted swap. The logged winning slot seed is not sufficient by itself to replay an exchanged trajectory: use the base seed and the full replica/run parameters. Repeated runs were deterministic with the same executable and parameters; cross-toolchain floating-point identity is not promised.

--steps counts proposals per replica per restart. The total number of restart opportunities is --restarts * --seed-tries. Early stopping is synchronized at dispatch boundaries. --keep-going disables early stopping for diagnostics; it is useful for observing exchanges in easy cases.

Certification enumerates all nonempty subsets of size at most floor(N/2), including both halves of balanced cuts. This is exponential in N. The program limits N to 24 and checks an estimated cut/cache memory budget before allocation (--memory-mb, default 512). This estimate is not a hard process-RSS limit. Practical limits can be much smaller.

Timing includes data export but excludes writing the final certificate/summary metadata. Verification-only runs have zero search time. Timings in the included examples are single local runs, not performance benchmarks.

## Tests and provenance

    make test CXX=clang++ PYTHON=python3
    make sanitize CXX=clang++ PYTHON=python3

See VALIDATION.md for the executed checks. The original supplied source had SHA-256:
7737b63732425d55a05a5314ea4ea4b685d87ee4f739509651481300ef7f8f50.

This implementation develops Dr. Zakaria Dahbi's AME7DIA research code. Refactoring, parallel-tempering implementation, and tests were produced with Codex assistance; details are recorded in [AUTHORS.md](AUTHORS.md). Software validation does not establish scientific novelty or reproduce historical runtime claims.

## Citation and licensing

Please credit **Dr. Zakaria Dahbi** when referring to AfAME. GitHub can display a software citation using [CITATION.cff](CITATION.cff).

Licensing terms have not yet been specified by the developer. See [COPYRIGHT.md](COPYRIGHT.md).
