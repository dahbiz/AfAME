# Validation record

Executed on the user's Mac, 15 September 2026. These are correctness and integration checks, not a comparison of search speed with other algorithms.

## Numerical and mathematical checks

- Rejected all three reducible polynomials identified in the audit.
- Checked every nonzero inverse for GF(2), GF(73), GF(137), GF(257), GF(4), GF(9), GF(121), GF(343), GF(625), GF(512), and GF(1024).
- Tested distributivity and 2-by-2 determinant/rank agreement on 400 random cases per field.
- Confirmed that the rank of [256] over GF(257) is 1.
- Compared full and incremental ranks, costs, edge failure counts, and total guidance weights after 3,200 Metropolis-Hastings moves, including rejected proposals.
- Checked guided-move detailed balance exactly for every neighboring pair of three-party binary graph states.
- Checked the exchange example T=(1,2), E=(0,10): log acceptance is -5, not +5.
- Confirmed the GF(4) edge labelled 2 has reduced purity 1/4 using independently constructed state amplitudes. Its reported entanglement is two bits.

## Parallel execution and output checks

- Identical seeds and parameters produced identical matrices, convergence logs, swap logs, restart logs, and certificate contents apart from timings in repeated four-thread runs.
- One-, two-, three-, four-, and five-replica cases completed. Every adjacent pair in each ladder was exercised.
- In the included four-replica diagnostic, 32 of 36 exchanges were accepted, with 36 accepted uphill local moves across the replicas. The final best state was correctly reported as non-AME.
- A five-qubit search found a state passing all 15 enumerated cuts.
- A fresh four-replica search at N=17, d=10001, seed=42 found a candidate by the 60-step dispatch boundary; 10 of 17 exchanges were accepted. Exhaustive C++ verification passed all 65,535 cuts. The independently recomputed result is saved with this example. This is a functional test, not a comparative benchmark or a claim of a novel state.
- Successful initialization and successful final restarts appear in restart statistics.
- Composite-dimension rank histograms account for all cuts in every factor. Output no longer presents a rank sum as a single-field rank with an inconsistent maximum.
- Invalid arguments and field sizes fail cleanly. The memory estimate rejects oversized requests before cut allocation. Existing output directories are preserved.
- Thread-pool exceptions propagate to the main thread without hanging the remaining workers.

## Independent verification

The standalone Python verifier uses its own arithmetic and a different irreducibility test. It rechecks exported certificates in prime, extension-field, and tensor-product cases, and rejects a deliberately forged positive verdict. Four-party extension-field runs exercise nontrivial 2-by-2 cuts in GF(4), GF(9), and GF(121).

The manuscript's AME(17,10001) matrix was rechecked in both the C++ program and the independent Python verifier: all 65,535 subsets passed, with both prime sectors checked per subset. Its matrix and the C++ certificate are included under examples.

## Runtime checks

The C++ core and integration suites passed with address and undefined-behavior sanitizers. A separate four-thread, two-restart exchange run completed under ThreadSanitizer with no race report.

Run the checks again using:

    make test CXX=clang++ PYTHON=python3
    make sanitize CXX=clang++ PYTHON=python3

The test suite uses temporary output directories and does not overwrite research results.

## Scope

No historical binary equivalence, cross-machine reproducibility, full parameter-space validation, or comparative performance study is claimed. Exhaustive rank certification remains exponential in party count. Default temperatures and replica counts require tuning for demanding searches. The examples distinguish verification of the manuscript matrix from a fresh search; neither is a discovery-time benchmark.
