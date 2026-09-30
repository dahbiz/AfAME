# Exact eight-party entropy-deficit minima

This research addition evaluates the manuscript's squared cut-rank deficit
objective over the complete published eight-qutrit classification. It establishes
`C*_3(8) = 16`. Together with the supplied binary enumeration, `C*_2(8) = 28`,
sector separability gives `C*_6(8) = 44` in the binary–ternary product family.
The existing manuscript certificate attains this value. This is not a claim about
unrestricted AME(8,6) existence or an unrestricted entanglement optimum.

## Reproduce all checks

From the AfAME repository root, with Python 3.8+ and a C++17 compiler:

```sh
python3 data/eight_qutrit_classification/run_exact_checks.py --output exact_optima_check --cxx clang++
```

Choose a new output directory each time. The command:

1. Checks the supplied input checksum and evaluates all 817 classes by finite-field
   elimination using AfAME's separate Python certificate verifier.
2. Independently parses the same records in C++ and enumerates all 6,561 stabilizer
   labels per class. It counts stabilizers supported in each subsystem without
   using Gaussian elimination. All published weight enumerators and distances must
   match, and all 132,354 cut ranks must agree between implementations.
3. Repeats both existing binary eight-party enumerations (133,632 graph extensions).
4. Verifies the qutrit witness with the unmodified native AfAME executable and checks
   the existing composite certificate with the Python verifier.
5. Confirms that a deliberately corrupted weight enumerator is rejected.

Exit code 2 from native AfAME means an exactly verified non-AME state; the runner
handles that expected outcome. Research enumeration does not change the stochastic
search algorithm or its defaults.

## Classification provenance and scope

The source is Lars Eirik Danielsen's
[Database of Nonbinary Self-Dual Quantum Codes](https://www.codetables.de/larsed/nonbinary/),
specifically `selfdualcodes3.txt.bz2`, retrieved on 29 September 2026.
`data/danielsen_length8.tsv` contains every length-eight record: 659 indecomposable
and 158 decomposable classes. Each line has six tab-separated fields: length,
decomposability, minimum distance, weight enumerator, automorphism group order,
and weighted graph edges. Vertex indices start at zero. `2*0-1` is a weight-two
edge; `0-1` is weight one; `-` is an empty edge set.

The full archive SHA-256 is
`869d056cbd48378bacd954a48db72ab9c8f68140e3791355bc461383ee07d618`.
The extracted-record SHA-256 is
`0644d313a153149893f7760d70473eaea03dceb4b30eac1444348b5b4c81d9dc`.
`data/provenance.json` also records each original source line number. To repeat
the extraction after downloading the approximately 32 MB archive:

```sh
python3 data/eight_qutrit_classification/extract_danielsen_database.py selfdualcodes3.txt.bz2 extracted
```

The mass identity using the published automorphism orders is checked exactly:
`sum(24^8 * 8! / |Aut(C)|) = product(3^i + 1, i=1..8) = 234870301468364800`.
This checks consistency of the input; it does not regenerate the classification
or independently compute the automorphism groups. Completeness relies on the
published classification. Local Clifford transformations preserve all subsystem
spectra, and party permutations only reorder the subsets, so the objective is
constant on each equivalence class. Graph representatives cover every qutrit
stabilizer state up to these operations.

References:

- L. E. Danielsen, *Graph-based classification of self-dual additive codes over
  finite fields*, Advances in Mathematics of Communications 3, 329–348 (2009),
  [doi:10.3934/amc.2009.3.329](https://doi.org/10.3934/amc.2009.3.329).
- L. E. Danielsen, *On the classification of Hermitian self-dual additive codes
  over GF(9)*, IEEE Transactions on Information Theory 58, 5500–5511 (2012),
  [doi:10.1109/TIT.2012.2196255](https://doi.org/10.1109/TIT.2012.2196255).

The classification records are attributed third-party research data; their
inclusion does not assert authorship of that classification or change its rights.

## Findings

All cuts of sizes 1–4 are counted, including both complementary balanced subsets.
The minimum ternary cost is 16. One-based records 741–750 attain it:

| Balanced deficits (zero, one, two) | Number of classes | Total balanced deficit |
|---|---:|---:|
| (54,16,0) | 7 | 16 |
| (60,8,2) | 2 | 12 |
| (66,0,4) | 1 | 8 |

All ten have full-rank smaller cuts. Equal squared cost does not imply equal mean
entropy, worst-cut deficit, or number of maximally mixed reductions. For the paper's
composite witness, the binary and ternary sector costs are 28 and 16; 34 of 70
balanced reductions are maximally mixed. This does not prove 34 is the largest
possible number of such reductions.

`results/` contains the complete outputs, all ten minimizing matrices, provenance
checks, a standard AfAME certificate for record 741, and the independent agreement
report. `results/run_summary.json` records the full reproducibility run. These new
programs are manuscript supplements, not part of the earlier public revision
`6e6ec6088aad` unless separately released.
