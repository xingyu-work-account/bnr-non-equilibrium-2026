# Independent re-review of remediation pass 5 (2026-10-25)

This review covers commit `fccd38a` on branch `remediation/rereview-2026-10-19` (parent `2097d75`, whose code and outputs equal `8e10d0e`). Pass 5 is report-only: Z1 to Z5 of `docs/remediation_prompt_2026-10-22_pass5.md`, closing P2-A to P2-E of `docs/rereview_2026-10-22_pass4.md`. The claims in `docs/progress.md`, section "Remediation pass 5", were treated as claims to verify. The procedure is `docs/rereview_prompt_2026-10-25_pass5.md`.

Every command was run on clean detached worktrees of `fccd38a`: one for run-all, one for pytest and one for scratch checks. The untracked `src/bnr_non_equilibrium_2026/` stub and the review prompt are not part of the commit and were not present.

## Decision

**GO WITH DOCUMENTED NON-BLOCKING DEBT.** Pass 5 closes P2-A to P2-E of the pass-4 review and keeps its report-only promise:

- Against `8e10d0e`, every figure and 63 of 66 tables are byte-identical.
- In the three tables that changed, every pre-existing column is identical as raw text, and the only additions are the Z1 and Z2 columns.
- The source diff adds code and changes no computation.

All three commands pass on a clean checkout of `fccd38a`, and the outputs regenerate byte-identically at a different worker count. No acceptance criterion fails, no failure condition is met, no test is weakened, and all five judgment calls are accepted.

Two P2 items remain:
- **P2-1** is optional analysis debt. The selected model's calibration has C0 on its 300 mg N/L cap; the MODEL-BASED quantities are unaffected.
- **P2-2** is two wording nits.

Neither needs code before merging.

## 1. Commands (run by the reviewer, 2026-10-25)

| Command | Result |
| --- | --- |
| `uv run ruff check` | `All checks passed!`, exit 0 |
| `rm -rf outputs && uv run python -m bnr_neq run-all` (`BNR_NEQ_WORKERS=6`) | Exit 0, wall 1:38:55. Stage messages: verdicts NH4/EC/ORP FAIL, pH PASS; segments references 11/13; thermo 9/11; dynamics 28/28; model selected `biomass_decay`, 13/15 identifiable; states 12 % classified, idle between 5 and 93 %; pulses 274 computed, 78 not computed. These match `progress.md`. |
| `diff -rq` against a copy of the committed `outputs/` taken before the run | No differences in `tables/` (66 files), `reports/` or `figures/` (PNG and SVG). `processed/` is byte-identical to the implementer's final run in the main checkout. The outputs therefore do not depend on the worker count (6 here, all cores for the implementer). |
| `uv run pytest -q -rs` (`BNR_NEQ_WORKERS=6`, after run-all) | `458 passed in 24764.94s (6:52:44)`, exit 0, none skipped (`-rs` lists none). This matches the implementer's count (444 + 13 + 1). The suite includes the two fresh run-alls with byte-identical CSV and SVG outputs, and the pass-3 regression test, which now reads fixtures. |

Environment: CPython 3.14.8, uv 0.13.1, numpy 2.5.4, pandas 3.0.6, scipy 1.18.1. The implementer used the default worker count (all cores); the reviewer used 6.

## 2. The report-only guard

- **Committed outputs, `8e10d0e` → `fccd38a`.** Exactly five output files changed: `analysis_report.md`, `model.md`, `model_fit.csv`, `model_calibration_starts.csv` and `model_volume_scenario_summary.csv`. Every figure and the other 63 tables are unchanged. That includes `pulses.csv`, every `nh4_truth_*.csv`, `model_identifiability.csv`, `model_volume_scenarios.csv` and the `states_*` tables.
- **Column-level check.** In the three changed tables, every pre-existing column is identical **as raw CSV text**, row for row (9, 33 and 3 rows). The only new columns are `cal_at_bounds` (Z2, in two tables) and `C_start_ratio_to_constant_volume`, `phi_start_median`, `C_end_ratio_to_constant_volume`, `phi_end_median` and `ratio_within_phi_tolerance` (Z1).
- **Source.** `git diff 2097d75 fccd38a -- src` touches only `model.py` and `report.py`, and adds without changing any existing computation:
  - `at_bounds` and `bounds_text`;
  - `phi_scaling` and `volume_scaling_text`;
  - the constants `BOUND_TOL_FRACTION`, `PHI_SCALING_TOL` and `HALDANE_CONDITIONING_NOTE`;
  - the report functions `mechanism_summary_txt`, `mechanism_table` and `alk0_sentence`.

  `cal_at_bounds` reads the fit vectors already in memory (`at_bounds(lay, f.x, f.free)`), not a refit. The one edited line in `volume_scenario_refits` binds the constant-volume C(t) to `c_cv` and computes the same medians. No fit, bound, start, budget, tolerance or selection rule changed.
- **Regenerated outputs.** The reviewer's regenerated tables, reports and figures are byte-identical to the committed `fccd38a` outputs, so the guard verified above on the committed outputs holds for a fresh run as well. `processed/` (not committed) is byte-identical to the implementer's final run, which `progress.md` reports as identical to the `8e10d0e` baseline. The reviewer could not compare `processed/` directly with a `8e10d0e` run, because the pass-4 regeneration was not kept. The tables derived from it (every `nh4_truth_*.csv`, `pulses.csv`, `states_*`) are byte-identical to `8e10d0e`.

## 3. P2-A to P2-E, item by item

| Previous item | Work item | Status | Evidence |
| --- | --- | --- | --- |
| P2-A: scenario C(t) "above the band" unexplained | Z1 | **Resolved** | Hand recomputation, end window: 569.5445/293.8491 = **1.93822** at 1 mL/d. The median of V0/V(t) over the last 0.1 d bins gives **1.91828** at 1 mL/d, 1.31467 at 0.5 and 1.13594 at 0.25; every value equals the table to the last digit. The largest \|ratio/φ − 1\| is 0.020, inside the registered 0.05. One function (`volume_scaling_text`) writes the statement in §1 point 9, §4 and `model.md`, and its branch is decided by `ratio_within_phi_tolerance`. The negative and boundary tests exercise the "differs" branch and the exact-0.05 case. |
| P2-B: bound-constrained calibrations unflagged | Z2 | **Resolved** (P2-1 below) | Calibration flags: `null_no_consumption` a; `biomass_decay` C0; haldane T20C Ks, a; haldane T25C a; combined T20C and T25C Ks, C0, a; no alkalinity candidate. These agree with the reviewer's independent refits in the pass-4 review (haldane T20C: Ks 991.7 and `a` on the bound; combined T25C: Ks, C0 and `a`; alkalinity_FA T25C: none). `bounds_text` in §4 and `model.md` names every flagged candidate and names the reference separately. |
| P2-C: conditioning not stated | Z3 | **Resolved** | Dated correction "Correction (2026-10-22, pass 5, Z3)" added to the pass-4 section of `progress.md`. The diff of `progress.md` deletes no line, so the original text is intact. The implementer's reproduced condition number of 4.224e6 (kinetic block 1.441e4, \|corr(Ks, S_opt)\| 0.8925) matches the reviewer's own 4.223e6 from the pass-4 review. `model.md` carries the note, marked as a quoted diagnostic. |
| P2-D: §1 point 7 too long | Z4 | **Resolved** | 142 words (was 601). It keeps the five required conclusions: the selected model and its skills; "holds" with the two selection skills; one ranking and test range per other-mechanism candidate; Alk0 in one sentence; and the time-confounding caveat. Every number traces to `model_fit.csv`. The ranges, rankings, Alk0 minima and bound flags are in one §4 table (`mechanism_table`), and the test asserts the table text in §4. |
| P2-E: regression test could skip | Z5 | **Resolved** | `_pass3` reads `tests/fixtures/pass3/` and calls `pytest.fail` if a file is missing (`test_a_missing_pass3_fixture_fails_the_regression_test`). The reviewer recomputed every fixture from git: the three pass-3 tables, the selected row of `model_fit.csv` and the 16 MODEL_BASED rows of `pulses.csv` are byte-identical to the `git show ad51a6a:…` extracts. The three pass-4 fixtures are byte-identical to `git show 8e10d0e:…`. |

## 4. Scientific notes

- **The selected model's calibration is bound-constrained.** Z2 makes visible that every `biomass_decay` start inside Δ < 3.84 has C0 on the 300 mg N/L cap, and so do the combined candidates. This was already in the pass-3 diagnosis (C0 = 299.9994), but it was never in the outputs.

  What it means:
  - The **MODEL-BASED quantities are not affected.** They come from the whole-record refit, whose C0 is 185.8 mg N/L, inside the bounds (`model_identifiability.csv`, `optimizer_value_not_an_estimate`), and reported as not identifiable. So `bounds_text`'s "the identifiability of the selected model is that of its whole-record refit" is correct.
  - The **selection skill** (0.958) does come from a C0-capped calibration. C0 lies along the non-identifiable C0-a-b0 valley, so the cap pins a nearly flat direction, and the skill is unlikely to move much. But the sensitivity of the selection to that arbitrary cap has never been tested.
  - The margin to the best other start (0.935) is 0.023. Testing the cap would mean changing a bound, which this pass correctly forbids. **P2-1** (optional analysis debt; not a defect of pass 5).
- **Wording.** `bounds_text` says of every flagged candidate that "their natural-unit kinetic parameters are not estimates". For `biomass_decay` the flagged parameter is C0, not a kinetic one. The next sentence ("the identifiability of the selected model … is that of its whole-record refit") points the reader to the reported v0 and kd, so this is a nit. §1 point 9 says the N amount "is unchanged", which is true within the 0.05 tolerance (≤ 2 % here). "Unchanged within the tolerance" would be exact. Both nits are recorded under P2-2.

## 5. Tests

- **New tests against the prompt's table.** `tests/test_pass5_fixes.py` has 13 tests, and `tests/test_pass4_fixes.py` gains 1.

  | Type | Tests |
  | --- | --- |
  | Positive | Z1 ratio = φ on records generated under 1 mL/d (two laws), with N amount conservation to 1e-9. Z2 flag for a fit whose true kd is 0. |
  | Negative | A ratio outside the tolerance gives the "differs" text naming the scenario and window. An interior fit gives no flag. A missing fixture fails. |
  | Boundary | \|ratio/φ − 1\| = 0.05 exactly is inside and the next float is outside. A parameter at exactly 1e-3 × width is flagged and the next float inside is not. A fixed parameter on its bound is not flagged. |
  | Regression | `test_real_pre_existing_columns_equal_8e10d0e` against the `8e10d0e` fixtures. |
  | Integration | The table's φ and ratios recomputed from the grid; the one Z1 statement in all three places; the flags equal to the chosen start's and named in §4 and `model.md`; point 7 ≤ 160 words, with the required content, and the §4 table. |

- **Changed existing test.** `test_pass4_fixes.test_real_report_states_the_skill_range_and_the_limits` looks for the `alkalinity_limited` "min to max" strings in §4 instead of §1 point 7, because Z4 moved them there. Every other assertion is unchanged: the starts inside, the minima against the table, the chosen start reproducing every skill, the "holds" statement in point 7, the limits text and the status row. **Not weakened.**
- The `_pass3` rewrite (Z5) keeps every assertion of the regression test and removes its only escape (`pytest.skip`). **Strengthened.**
- **No weakened test found.**

## 6. Rulings on the implementer's judgment calls

| # | Judgment call | Ruling | Reason |
| --- | --- | --- | --- |
| 1 | Z1 tolerance 0.05, fixed when the ratios were already known (disclosed) | **Accept** | The prompt itself proposed 0.05. The argument does not use the ratios: it rests on half the refit RMSE and the existing 5 % `MULTIMODAL_LOG_TOL`. The observed maximum of 0.020 would pass any tolerance down to 0.02, so 0.05 was not chosen at the margin. The start window is correctly declared uninformative (φ − 1 ≤ 0.011). |
| 2 | Z2 tolerance 1e-3 × bound width in fitted coordinates, every free parameter | **Accept** | One rule for every parameter. It catches the solver's near-bound stops (Ks 991.7 is 0.0083 log units from 1e3, against a tolerance of 0.016) and flags only the outer 0.1 % of each range. |
| 3 | Z4 test limit ≤ 160 words, set after a 144-word preview (disclosed) | **Accept** | The prompt set "about 150". 160 is a stated margin for candidate-name lengths, and the requirement is the content, which the test also asserts. |
| 4 | Two new tests corrected during development (Z1 exact identity → conservation plus tolerance; Z2 C0 cap → kd at zero) | **Accept** | Both first versions asserted something false about the model, and no source changed: the state starts at day 1, when V is already 1 % below V0, and pushing C0 also drove v0 and kd onto bounds. The corrected tests assert what the prompt requires. |
| 5 | Baseline not regenerated | **Accept** | The pass-4 review regenerated `8e10d0e` byte-identically, and the comparison is with the committed tables. The guard was rechecked on the reviewer's own regenerated outputs (§2). |

## Remaining issues by severity

| Severity | Issue |
| --- | --- |
| P2-1 | The selection skill of `biomass_decay` (0.958, margin 0.023) comes from calibrations with C0 on the 300 mg N/L cap, and the combined candidates are capped too. The MODEL-BASED quantities come from the uncapped whole-record refit and are not affected. The sensitivity of the selection to the cap is untested. It needs a bound change, so it is optional analysis debt for any future model work, not a defect of this pass. |
| P2-2 | Wording nits: `bounds_text` applies "natural-unit kinetic parameters are not estimates" to `biomass_decay`, whose only flag is C0, a state parameter. §1 point 9 says the N amount "is unchanged" where "unchanged within the 0.05 tolerance (≤ 2 %)" is exact. |

## Readiness to merge

**Ready to merge.** The remediation plan of `docs/remediation_prompt_2026-10-19.md` (W1 to W8) and every item raised in the five re-reviews are now implemented and accepted. The remaining items are P2-1, which is optional analysis that would need a bound change, and P2-2, which is wording. Neither blocks a merge.

The scientific NO-GOs that remain in the status table are blocked by measurements, not by code:
- `lab_NH4_N_mg_L`;
- `initial_alkalinity_meq_L`;
- `DO_mg_L`;
- `NO3_N_mg_L` and `NO2_N_mg_L`;
- the vessel masses;
- `kLa_per_d`.

Recommended next steps: open a pull request of `remediation/rereview-2026-10-19` into `main`, with the five review documents as the record. Then write a measurement request that maps each missing quantity to the NO-GO it would unlock.

## What is sound and should be kept

- A report-only pass that proved its own guard: raw-text column equality against stored `8e10d0e` fixtures, in a test.
- Report statements generated from one function each and decided by a registered tolerance on a table, so the text cannot drift from the numbers.
- A bound flag that made a long-known fact visible (the selected calibration's C0 on its cap) instead of hiding it.
- A regression test that can no longer skip, with fixtures traceable to commits by a recorded command.
- Complete disclosure of every development correction, including the one found by the first full run.
