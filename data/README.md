# Reproducibility data

This directory contains the machine-readable inputs used for the manuscript
figures, the plotting and benchmark scripts, and the files needed to document
the held-out recombination test.

From the package root, regenerate the four external figures with:

```text
python3 data/plot_paper_figures.py
```

The script reads the CSV and JSON files in this directory and writes
`fig_tradeoffs.pdf`, `fig_recombination.pdf`, `fig_four_six.pdf`, and
`fig_seven_qubit.pdf` (plus PNG previews) in the package root, where the
LaTeX source expects them.

To validate and summarize the frozen 150-run benchmark, run:

```text
python3 data/summarize_recombination_benchmark.py
```

The supplied CSV is the frozen held-out result. The instrumented C++ source and independent verifier are included. Regenerate all 150 search runs, verify every exported certificate, compare integer results with the frozen CSV, and check nine baseline runs with:

```text
python3 data/run_recombination_benchmark.py --output /tmp/afame_benchmark --cxx clang++
```

Use a new output directory. Timings are machine dependent. The `benchmark_certificates/` directory contains the fresh 30 September 2026 reproductions. The protocol file records the original design; its prespecification date is a provenance assertion, not established by rerunning code.

Data files:

- `ternary_class_profiles.json`: all 817 classified eight-qutrit profiles.
- `near_ame_entropies.csv`: cut entropies for the four-quhex and seven-qubit examples.
- `open_frontier_entropies.csv`: open-frontier entropy summaries.
- `recombination_benchmark.csv`: 150 held-out sector-recombination outcomes.
- `figure_validation.json`: validation summary written by the plotting script.

The complete eight-qutrit reproducibility check is in
`eight_qutrit_classification/`. Run it from the package root with:

```text
python3 data/eight_qutrit_classification/run_exact_checks.py \
  --output /tmp/eight_qutrit_check --cxx clang++
```

This evaluates all 817 extracted length-eight records, checks the input hash,
compares finite-field elimination with independent stabilizer-support counting
on all 132,354 cuts, repeats the binary exhaustive checks, verifies the native
AfAME witness, and runs a corrupted-input negative control. The directory
contains the extracted source records, provenance hashes, evaluation programs,
independent C++ verifier, AfAME verification source, and machine-readable
results. The original compressed Danielsen archive is not redistributed; its
URL and SHA-256 hash are recorded in `input_provenance.json`.


## Trade-offs and eleven qutrits

```text
python3 data/tradeoffs/run_checks.py --output /tmp/afame_tradeoffs --cxx clang++
```

This checks both binary enumerations, all 817 ternary profiles, both aligned
product optima, and the published eleven-qutrit matrix against independent
stabilizer-support counts. Fresh checked outputs are in `tradeoffs/results/`.
The eight-party fresh outputs are in `eight_qutrit_classification/rechecked_results/`.
Classification coverage is inherited from Danielsen's published classification;
these programs do not regenerate that classification.

## Requirements and examples

Use Python 3.10 or newer and a C++17 compiler. Only the figure and bootstrap
summary scripts require third-party packages: NumPy and Matplotlib.
Install them in your preferred environment with `python3 -m pip install numpy matplotlib`.
The exact checks use the Python standard library.
`example_certificates/` contains the printed 17-party construction, small-system
examples, and original search baselines. The newer optima are in the check results.
The 17-party file was reconstructed from the manuscript's printed matrices and
verified across all cuts; it is not an original search-run log.

These files are prepared for release under the AfAME repository's `data/`
directory. This package does not establish that the repository upload has happened.
Preserve provenance and third-party attribution when publishing these files.
