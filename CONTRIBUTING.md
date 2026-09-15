# Contributing to AfAME

AfAME is developed by Dr. Zakaria Dahbi.

For a bug report, include the exact command, compiler version, platform, base seed, replica count, and relevant output. Include the certificate when reporting an incorrect AME verdict. Remove private local paths or unrelated research data before posting.

For a proposed change:

1. Describe the mathematical or implementation problem.
2. Keep field conventions and output schemas explicit.
3. Add a regression test that reproduces the defect or validates the new behavior.
4. Run make test and, where supported, make sanitize.
5. Report the validation performed and any change to search trajectories.

Do not claim nonexistence or global optimality from an unsuccessful finite search. Separate field arithmetic from residue-ring arithmetic. Changes to state conventions, certification, or acceptance probabilities need a clear mathematical justification.

Licensing terms are currently unspecified; discuss contribution and reuse terms with the developer before submitting code.
