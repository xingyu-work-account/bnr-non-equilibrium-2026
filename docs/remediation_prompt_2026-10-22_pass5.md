# Remediation prompt, pass 5: report clarifications and the last P2 items (report-only)

You are the implementer for `bnr-non-equilibrium-2026`, branch `remediation/rereview-2026-10-19`, starting from commit `8e10d0e`. First read `docs/rereview_2026-10-22_pass4.md`, the latest decision: **GO WITH DOCUMENTED NON-BLOCKING DEBT**, with five P2 items (P2-A to P2-E). Then read `docs/remediation_prompt_2026-10-21_pass4.md`. Its invariants, non-goals and acceptance criteria still apply.

This is a **report-only pass**. It adds explanations, flags and tests. It must not change any fitted value, start, selection, C(t), primary row, MODEL-BASED number, scenario value or figure. Every existing column of every table keeps its values. Only the new columns named below, the report texts and `progress.md` may change. If an item cannot be done without changing a value, stop and report it.

Treat each finding as a claim to reconcile against the code, not as an established fact. Before you edit anything:

- Copy the committed `outputs/` to a scratch location as the baseline.
- Verify every path, line and number cited below. Report any that differ.

Rules for the whole pass:

- Make the smallest coherent change for each item, and do not refactor unrelated code.
- Do not weaken or delete existing tests.
- Never fabricate a measurement, and never tune code to reproduce or avoid a number.
- Fix every tolerance you introduce, with its argument, in `docs/progress.md` before you compute the outputs it affects.
- You are not the approver; an independent reviewer decides.

The nine judgment calls ruled on in §10 of the pass-4 review were accepted. Do not revisit them.

## Work items

| Item | Change | Affected files (verify) |
| --- | --- | --- |
| Z1 (P2-A) | **Explain the volume-scenario C(t) against the band.** Report §1 point 9 shows the scenario C(t) "above the band" at 0.5 and 1 mL/d without explanation. The review found that each scenario C(t) is the constant-volume C(t) × φ = V0/V(t) to within 1-2 %: end ratio 1.138, 1.321, 1.938 against φ = 1.136, 1.315, 1.918. (a) Add to `model_volume_scenario_summary.csv`, for the start and end windows, the ratio of the scenario to the constant-volume median C and the median φ over the same window. (b) Fix a tolerance before computing, for example \|ratio/φ − 1\| ≤ 0.05, with its argument. (c) If every ratio is within the tolerance, state in §1 point 9, §4 and `model.md` that the scenario C(t) is the constant-volume C(t) scaled by V0/V(t): the N amount per initial litre and η are unchanged, concentrations scale by φ, and "above the band" reflects evaporative concentration of a band that is conditional on constant volume. If any ratio is outside the tolerance, say so instead. The text must follow the table, not this prompt. | `model.py` (`volume_scenario_refits`, about lines 1198-1245), `report.py` (`volume_refit_txt`, line 576; §4), `model.md` |
| Z2 (P2-B) | **Flag bound-constrained calibrations.** (a) Add `cal_at_bounds` to `model_fit.csv` and `model_calibration_starts.csv`: the names of the free parameters within a fixed tolerance of a bound in the fitted coordinates (for example 1e-3 × the bound width), or empty. Fix the tolerance in `progress.md` first. Apply it to every candidate and every start, including the observation parameters. (b) Report §4 and `model.md` name every candidate whose calibration is bound-constrained, and state that its natural-unit kinetic parameters are not estimates. The review found Ks at 1e3 for `haldane_FA_SCENARIO_T20C` and `combined_SCENARIO_T25C`, and `a` and C0 on bounds. The selected `biomass_decay` calibration may also be flagged: its C0 sat on the 300 mg N/L cap in the pass-3 diagnosis. Report whatever the flag shows, and do not change a bound to clear a flag. | `model.py` (`split_evaluation`, start rows), `report.py`, `model.md` |
| Z3 (P2-C) | **State what the reparametrization did and did not do.** In `progress.md` (pass 4 section, as a dated correction) and `model.md`: the (α, Ks, S_opt) form removed the near-perfect kinetic correlations (0.99998 → about 0.89), but not the ill-conditioning. The reviewer found a Jacobian condition number of 4.2e6 at the converged alkalinity_FA T25C calibration, and the alkalinity × FA natural parameters (vmax about 2,600 mg N/L/d, Ki about 0.005 mg NH3/L) are not estimates. Reproduce the condition number with a scratch script before quoting it; quote your own value. | `docs/progress.md`, `model.md` text in `model.py` (`write_report`) |
| Z4 (P2-D) | **Shorten report §1 point 7.** It is one paragraph of about 600 words. Move the per-candidate skill ranges and Alk0 profile statuses into one table in §4: candidate, starts inside Δ < 3.84, selection, test and forecast ranges, ranking against the selected model, and Alk0 minima. Keep §1 point 7 to about 150 words with the conclusions only: (1) the selected model and its skills; (2) the selection holds against every indistinguishable start; (3) the ranking against each other-mechanism candidate (undetermined or better or worse, with one range each); (4) the Alk0 evidence in one sentence; (5) the time-confounding caveat. Every number in point 7 must still come from the tables. | `report.py` (§1 point 7 at line 531; §4) |
| Z5 (P2-E) | **Make the pass-3 regression test independent of git history.** `test_real_selection_fit_and_model_based_numbers_equal_pass3` skips when `git show ad51a6a:…` is unavailable, so it would pass silently in a source tarball. Store the reference values it compares under `tests/fixtures/pass3/`: the selected row of `model_fit.csv`, `model_identifiability.csv`, `nh4_truth_summary.csv`, `nh4_truth_vectors.csv`, and the MODEL_BASED rows of `pulses.csv`. Extract them once from `ad51a6a` with `git show`, never by regenerating them, and record the command in `progress.md`. The test compares against the fixtures and fails, never skips, if a fixture is missing. Keep every assertion it has. | `tests/test_pass4_fixes.py`, `tests/fixtures/pass3/` |

## Invariants

All invariants of passes 1 to 4 hold. New for this pass:

- Every pre-existing column of every output table is byte-identical to `8e10d0e`, and so is every figure and `processed/` file.
- A statement in the report about the scenario C(t) or about bound-constrained calibrations follows its table and its registered tolerance.

## Non-goals

- No change to bounds, parametrizations, starts, budgets, tolerances of the fits, the selection rule, the split, or any registered threshold.
- No new candidate, no refit, no new scenario, and no choice of a primary volume scenario.
- No change to the sensor verdicts, the S2 logic, the W2 bound, the energy family or the state classifier.
- No new or assumed measurements.

## Tests

| Type | Required tests |
| --- | --- |
| Positive | Z1: on a synthetic record generated at a known evaporation rate, the ratio equals φ within the registered tolerance. Z2: a synthetic fit pushed onto a bound is flagged with that parameter's name. |
| Negative | Z1: a ratio outside the tolerance produces the "differs" statement, not the "scaled" one. Z2: a fit with every parameter inside its bounds has an empty flag. Z5: a missing fixture makes the regression test fail, not skip. |
| Boundary | Z1: a ratio exactly at the tolerance takes the documented branch. Z2: a parameter exactly at the tolerance distance from a bound is flagged as documented. |
| Regression | Existing tests still pass. `test_real_report_states_the_skill_range_and_the_limits` may look for the ranges in the §4 table instead of §1 point 7 (Z4 moves them), but must keep every assertion on the values. Explain that change in `progress.md`. A new test asserts that every pre-existing column of `model_fit.csv`, `model_calibration_starts.csv` and `model_volume_scenario_summary.csv` equals `8e10d0e`. |
| Integration | On the real run-all outputs: the summary table carries the ratios and φ; §1 point 9 contains the statement that its tolerance decides; `cal_at_bounds` exists and §4 names every flagged candidate; §1 point 7 is at most about 150 words and the §4 table holds the ranges. |
| End-to-end | `uv run python -m bnr_neq run-all` regenerates all outputs from the raw files with fixed seeds; two runs are byte-identical in CSV and SVG (the existing determinism test). |

## Validation commands

```bash
uv run ruff check
```

```bash
uv run pytest -q
```

```bash
rm -rf outputs && uv run python -m bnr_neq run-all
```

Run pytest and run-all one after the other, not at the same time. Then diff `outputs/tables`, `outputs/reports`, the figures and `processed/` against the baseline copy. The only differences allowed are the new columns of Z1 and Z2 and the report texts. Explain each one in `docs/progress.md`, under "Remediation pass 5", by work item. Run the full suite once on the final code, after the last edit. Do not regenerate single stages after the final full run. Commit when all three commands pass, and do not push.

## Acceptance criteria and failure conditions

Accept only if all of the following hold:

- Z1 to Z5 are implemented and tested.
- All three commands pass on the final code.
- The outputs are reproducible.
- The diff against `8e10d0e` contains only the new columns, the report texts and `progress.md`.

The work fails if any of the following occurs:

- any pre-existing value, figure or `processed/` file changes;
- a report statement contradicts its table or its registered tolerance;
- a bound, parametrization or fit setting is changed;
- the regression test can still skip;
- an existing test is weakened.

Report exactly what was changed, the commands run and their full results.
