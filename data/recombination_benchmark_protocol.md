# AfAME sector-recombination benchmark protocol

Frozen 23 September 2026, before the held-out runs below.

## Hypothesis and scope

For square-free `d=6`, the squared cut-rank cost equals the sum of independent GF(2) and GF(3) sector costs. At a fixed sequence of accepted states, storing the best matrix seen for each field and recombining them after the search cannot increase the final cost relative to storing only the best jointly encountered matrix. The test measures how often and by how much this exact postprocessing improves AfAME's output. It makes no claim that replica exchange is faster than another optimizer.

## Exploratory data already inspected

Seeds 200–219, `N=10,11,12`, `d=6`, four replicas, 500 steps each, one restart. Strict cost decreases appeared in 3/20, 9/20, and 12/20 runs respectively. These seeds are excluded from the confirmation set.

## Held-out confirmation design

Use seeds 1000–1049 at each of `N=10,11,12`, giving 150 runs. Fix four replicas, 500 steps per replica, one restart, temperature ladder 0.25–50, guidance 0.85, swap interval 10, and log interval 500. Every run has at most 2000 local proposals. The instrumented implementation records both the best jointly encountered cost and the recombined cost on the same trajectory. Its source derives from local AfAME revision `6e6ec6088aad22ba0e6dde8fa7ec5bbbbcacbc2c` with the published `--no-exchange` benchmark patch; the new tracking uses no random draws. Compare final costs pairwise within each run. The primary setting is `N=12`; `N=10,11` are corroboration settings.

Primary endpoint: mean paired reduction in squared cost at `N=12`. Report a seed-bootstrap 95% interval, the median paired reduction, and the fraction of runs with a strict reduction. Report the same descriptive outcomes for `N=10,11`, with no pooled significance claim. Independently recompute every final certificate using AfAME's Python verifier. Check that a sample of the instrumented code's recorded joint costs equals costs from the uninstrumented baseline at identical seeds and options. The result is a methods effect within the quadratic-phase ansatz; it is not a general AME construction advantage or a comparison against other optimizers.
