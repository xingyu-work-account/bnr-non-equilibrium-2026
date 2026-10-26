# Independent re-review prompt for remediation pass 5 (2026-10-25)

You are an independent reviewer who did not implement the change. The previous re-review (`docs/rereview_2026-10-22_pass4.md`) decided GO WITH DOCUMENTED NON-BLOCKING DEBT and left five P2 items (P2-A to P2-E). The implementer then made remediation pass 5, a **report-only** pass, committed as `fccd38a` on branch `remediation/rereview-2026-10-19` (parent `2097d75`, whose code and outputs equal `8e10d0e`). Decide whether that pass closes P2-A to P2-E without changing anything that was accepted.

Treat everything in `docs/progress.md`, section "Remediation pass 5", as claims to verify, not as facts. Use `git diff 2097d75 fccd38a` to see exactly what pass 5 changed.

Reference documents:
- `docs/remediation_prompt_2026-10-22_pass5.md`: work items Z1 to Z5, the report-only rule, tests, acceptance criteria and failure conditions.
- `docs/rereview_2026-10-22_pass4.md`: P2-A to P2-E and the nine judgment calls of §10. Do not reopen those calls.

## Procedure

1. **Run the commands yourself** on a clean checkout of `fccd38a`, and record the full results. A failed command is NO-GO.
   - `uv run ruff check`
   - `rm -rf outputs && uv run python -m bnr_neq run-all`. The implementer reports exit 0 after 1:42:28 under load. Use a worker count different from the implementer's default (all cores).
   - `uv run pytest -q -rs`, after run-all and not at the same time. The implementer reports `458 passed`, none skipped.
   - Diff the regenerated `outputs/tables`, `outputs/reports` and figures against a copy of the committed outputs taken before the run. They must be byte-identical.

2. **The report-only guard, the core check.** Against `8e10d0e`, on your own regenerated outputs:
   - Every figure and 63 of 66 tables must be byte-identical. That includes `pulses.csv`, every `nh4_truth_*.csv`, `model_identifiability.csv`, `model_volume_scenarios.csv` and the `states_*` tables.
   - In `model_fit.csv`, `model_calibration_starts.csv` and `model_volume_scenario_summary.csv`, every pre-existing column must be identical as raw CSV text, and the only new columns may be those of Z1 and Z2.
   - The only other changes allowed are `analysis_report.md`, `model.md`, `progress.md`, tests and fixtures.
   - Check `git diff 2097d75 fccd38a -- src` for any change to a fit, bound, start, budget, selection rule or primary computation.

3. **Close P2-A to P2-E item by item**, with file and line evidence, and check every sub-item of Z1 to Z5.
   - **Z1.** Recompute one ratio and its φ by hand, for example the end window at 1 mL/d: 569.54/293.85 against the median V0/V(t) over the last 0.1 d. Confirm that one function writes the statement in §1 point 9, §4 and `model.md`, and that the branch is decided by the table.
   - **Z2.** Confirm that `cal_at_bounds` is computed from the in-memory fit vectors and not from a refit. Check one flag against a refit of your own, for example `haldane_FA_SCENARIO_T20C` Ks. Confirm that §4 and `model.md` name every flagged candidate. The selected `biomass_decay` calibration is now flagged (C0 on the 300 mg N/L cap). Judge whether the report states what this means for the selected model and its MODEL-BASED quantities, and whether "its identifiability is that of its whole-record refit" is correct and sufficient.
   - **Z3.** Confirm that the condition number was reproduced (4.2e6), and that the correction is dated and placed in the pass-4 section without rewriting the original text.
   - **Z4.** Confirm that §1 point 7 keeps the five required conclusions, that every number traces to `model_fit.csv`, and that the §4 table carries what point 7 dropped.
   - **Z5.** Confirm that the fixtures equal `git show ad51a6a:…` and `git show 8e10d0e:…` extracts, by recomputing at least one, and that a missing fixture fails rather than skips.

4. **Judge test adequacy.**
   - Check `tests/test_pass5_fixes.py` (13 tests) and the new test in `tests/test_pass4_fixes.py` against the prompt's positive, negative, boundary, regression, integration and end-to-end rows.
   - For the changed existing test, `test_pass4_fixes.test_real_report_states_the_skill_range_and_the_limits`, compare it with `git show 2097d75:tests/test_pass4_fixes.py`. Decide whether it was changed because Z4 moved the ranges, or whether it was weakened.

## Rule explicitly on the implementer's judgment calls

1. **Z1 tolerance.** |ratio/φ − 1| ≤ 0.05, fixed when the ratios were already known (disclosed), argued from the refit residual scatter and `MULTIMODAL_LOG_TOL`. The start window is declared uninformative.
2. **Z2 tolerance.** 1e-3 × the bound width in fitted coordinates, applied to every free parameter, including the reference model's `a`.
3. **Z4 word limit.** The test asserts ≤ 160 words. The limit was written after a preview of the point-7 text (144 words, disclosed).
4. **Test corrections during development.** The Z1 positive test changed from "ratio = φ to 1e-9" to amount conservation plus the tolerance. The Z2 positive test changed from C0 at its cap to kd at zero. No source changed for either.
5. **Baseline not regenerated.** The comparison is against the committed `8e10d0e` outputs, which the pass-4 review regenerated byte-identically.

## Decision

Classify each remaining issue as P0, P1 or P2. Any of the following gives NO-GO:
- an unresolved P0 or P1;
- a failed or non-reproducible command;
- any pre-existing value, figure or `processed/` file that changed;
- an unverified acceptance criterion or a met failure condition of `docs/remediation_prompt_2026-10-22_pass5.md`;
- a rejected judgment call rated P0 or P1;
- a weakened test.

Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT listing each P2 item. If the decision is GO or GO WITH DEBT and no P2 item needs code, say whether the remediation branch is ready to merge into `main`. Write the review to `docs/rereview_2026-10-25_pass5.md`, using the structure of `docs/rereview_2026-10-22_pass4.md` and shortened to fit a report-only pass.
