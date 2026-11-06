# Registration of the W3 process-state rules and association test (pass 6, 2026-10-26)

registered_sha256: e34bff2c1693bfade1758db70a897f5d9676e7c34500509e91beaa2b97edf270

This file freezes everything that decides a process state or the NH4 association test in `src/bnr_neq/states.py`. The hash is the SHA-256 of `states.registration()`: the three rules, every threshold, window, seed and resample count, the sensor exclusions, the indicator partitions, and the source code of every decision function (`states._REGISTERED_FUNCTIONS`). `tests/test_pass6_fixes.py` fails if the code no longer matches this hash. A change to the rules must update this file in its own commit, with the reason, before any output computed with the new rules is reported.

## Status: revised after first results, not pre-registered

The rules were fixed from sensor resolution and physical scales before the first states were computed (pass 2, `docs/progress.md`). They were not committed before their first results, and two parts changed after results existed:

1. **Pass 2:** the rule that each bootstrap group must span at least 2 blocks (`MIN_BLOCKS_PER_GROUP`). It was added after the first run showed a percentile interval that did not cover its own estimate, with 89 and 66 bins against a 101-bin block.
2. **Pass 3 (X2):** the Theil-Sen slope intervals are widened by √(n/n_eff) for the lag-1 autocorrelation of the window residuals. This changed which bins are classified, after the pass-2 states and associations had been computed.

Both changes were disclosed at the time and were made for validity reasons, not to move a result. Still, the rules now in use are not the ones registered before the data were seen, so every output says `REVISED AFTER FIRST RESULTS (not pre-registered)` (`states_rules.csv`, `states_comparison.csv`, `states.md`).

## Added in pass 6 before it was run on the record: null coverage of the block bootstrap (G2)

`null_response`, `bootstrap_null_coverage`, `null_coverage` and `interval_reading`, and the constants `NULL_COVERAGE_PHIS` (0.5, 0.9, 0.95, 0.98), `NULL_COVERAGE_REPLICATES` (200) and `NULL_COVERAGE_MIN_OK` (0.90).

**Design.** The check keeps the record's eligible bins and group labels. It uses a response with no group effect: the 101-bin least-squares slope of AR(1) noise, so neighbouring values share 100 of their 101 bins, as in the real response. It applies the registered block rule with 999 resamples.

**Decision rule.** An interval is "adequate" only if it covers the true difference 0 in at least 0.90 of the replicates at every φ. Otherwise the comparison is labelled `UNDER-COVERS`, and an interval excluding 0 is not read as evidence of a difference.

**What the check cannot change.** The block rule stays as registered whatever the check finds.
