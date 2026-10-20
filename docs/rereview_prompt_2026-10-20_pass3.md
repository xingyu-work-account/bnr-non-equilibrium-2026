# Independent re-review prompt for remediation pass 3 (2026-10-20)

You are an independent reviewer who did not implement the change. The previous re-review (`docs/rereview_2026-10-19_pass2.md`) gave NO-GO on two P1 defects: N1 (the W5 model comparison rested on unconverged calibrations) and N2 (the idle fraction was a biased lower bound presented as a value). The implementer then made remediation pass 3, committed as `adda251` on branch `remediation/rereview-2026-10-19` (parent `954f87a`). Decide whether that pass resolves the previous decision.

Treat everything in `docs/progress.md`, section "Remediation pass 3", as claims to verify, not as facts. Review the cumulative tree. Use `git diff 954f87a adda251` to see exactly what pass 3 changed.

Reference documents:
- `docs/remediation_prompt_2026-10-19_pass3.md`: work items X1 to X8, invariants, tests, acceptance criteria and failure conditions.
- `docs/rereview_2026-10-19_pass2.md`: the previous decision, N1, N2, the P2 table and the five accepted judgment calls. Do not reopen those calls unless pass 3 changed their basis.
- `docs/remediation_prompt_2026-10-19.md`: its invariants, non-goals and acceptance criteria still apply.

## Procedure

1. **Run the commands yourself** on a clean checkout of `adda251`, and record the full results. A failed command is NO-GO.
   - `uv run ruff check`
   - `uv run pytest -q`. This takes 1.5 to 2.5 h: the suite contains two complete run-alls and an M6 compute, and the model stage now runs 33 calibration starts. The implementer reports 418 passed.
   - `rm -rf outputs && uv run python -m bnr_neq run-all`. This takes about 1 h. `BNR_NEQ_WORKERS` sets the process count; results must not depend on it.
   - Compare the regenerated `outputs/tables`, `outputs/reports` and the SVG figures with a copy of the committed outputs taken before the run, using `diff -rq`. They must be byte-identical.
   - Optional determinism check: rerun the model stage with `BNR_NEQ_WORKERS=1` and confirm that `model_fit.csv` and `model_calibration_starts.csv` are identical.

2. **Close the previous re-review item by item.** For N1, N2 and every P2 row of the "Remaining issues by severity" table in `docs/rereview_2026-10-19_pass2.md`, state resolved, partly resolved or unresolved, with file and line evidence. For each of X1 to X8, check every lettered sub-item of the pass-3 prompt.

3. **N1 / X1, the converged comparison.**
   - Verify the diagnosis claims in `progress.md`: the original stop reason, the evaluation counts, the parameters at bounds, Alk0, and the warm-start runs. Reproduce at least one.
   - Confirm that the convergence policy (`model.CONVERGENCE_POLICY`, `calibration_starts`, `CAL_MAX_NFEV`) is identical for every candidate, and that no candidate gets a budget the others do not. Judge whether "4 starts, kinetic start x 0.5/1/2/4, 8,000 evaluations each" is a fair common rule. The null model has a single start.
   - From `model_calibration_starts.csv`, check that the chosen start of each candidate is the converged start with the lowest RSS. Check that only converged candidates can be selected (`select_converged`), and that unconverged rows carry `not computed: calibration did not converge (status …)` with their skills kept.
   - Judge the result. `alkalinity_limited` converges; its lowest-RSS start has Alk0 = 80 (no limitation inside the calibration window, selection skill −0.91), while other starts converge at Alk0 = 22.6 (selection skill 0.93) at Δ = 3.0. Both alkalinity × FA candidates fail to converge in 4 x 8,000 evaluations. Is "decay selected among converged candidates; alkalinity × FA undecided; the pass-2 claim of better out-of-sample alkalinity skill withdrawn" the correct conclusion? Is selecting the calibration by lowest RSS, rather than by any skill, correct here?
   - Verify the fitted Alk0 reporting and the comparison with the S2 requirement (`alk0_comparison` in the report, `nh4_truth_w2_consistency.csv`).
   - **Disclosed post-hoc rule:** `alk0_profile` withholds the Alk0 interval when another converged start lies below Δ = 3.84 (`MULTIMODAL_LOG_TOL`). It was added after the first full run showed an invalid one-sided limit. Rule whether it is a validity condition or tuning. Recompute Δ for the 22.6 minimum by hand (n_eff from the calibration residuals).
   - Check the performance change in `model.fit`: the simulation is truncated at the window end because the model is causal. Confirm that the residuals inside the window are unchanged.
   - Confirm that the small changes in C(t) (+0.015 to +0.054 %) and in every MODEL-BASED and C(t)-derived quantity are explained by the new biomass_decay start, as `progress.md` claims.

4. **N2 / X2, the idle fraction.**
   - Check the n_eff-scaled interval in `states.rolling_theil_sen` against the prompt's equation: ρ is the lag-1 autocorrelation of the window residuals, capped at 0.999, and the same scaling applies to C(t), pH and ORP. Confirm it is stated in `states_rules.csv` and the docstring.
   - Verify that the method and the 0.5 coverage threshold were recorded in `progress.md` before the outputs they affect. The repository history has a single pass-3 commit, so judge from the log and the code.
   - The implementer reports that on synthetic AR(1) data at φ = 0.9 the scaled interval covers about 0.85, not 0.95, because ρ̂ is biased low on 101 points. They kept the estimator the prompt specifies. Judge whether the test (`test_autocorrelation_scaled_interval_covers_an_ar1_slope`) is an adequate implementation of the prompt's "about 95 %" test, and whether the residual undercoverage is stated correctly in the report.
   - Recompute the idle bounds of one interval from `states_bins.parquet` by hand. Confirm that every idle metric in `pulses.csv` carries LOWER_BOUND, UPPER_BOUND or CONDITIONAL in its name and its coverage in its basis. Confirm that report §1 point 8, the §2 table, Fig. 8c and the status table show bounds and coverage, and that no status row reads GO for the idle fraction.
   - Judge whether the state layer still supports any claim the report makes from it.

5. **X3, the checker.** Try to construct an output that breaks the one-observation-model invariant and is not caught by `model.check_single_observation_model`:
   - an empty or unknown model;
   - a derived metric labelled `electrode_reading`;
   - a reading-based column without SENSITIVITY in a wide table;
   - a series on an NH4-independent row;
   - "not MODEL-BASED" text.
   Confirm that the checker runs before `pulses.csv`, `nh4_truth_windows.csv`, `states_*`, `segments.csv` and `thermo_oxygen.csv` are written. Check that the 28 `NH4-independent` rows of `pulses.csv` really rest on no NH4 quantity.

6. **X4 to X7.**
   - X4: no reading-based column or metric may lack SENSITIVITY in its name. `consumption_uncorrected_*` stays as a measured raw slope; judge whether that is acceptable.
   - X5: the η band lies in [0, 1], and the baseline windows carry `w2_bound_check`.
   - X6: the R8 rule was extended symmetrically instead of being argued away. Is that justified? Pulse 1's rate goes from 7.15 to 4.23 mg N/L/d, and 7.15 is kept as `consumption_weighted_mgN_L_d_R8_CASE_incl_fragments`. Confirm that only window 14 changed status, and that the downstream changes in dynamics, kinetic fits, thermo and the M4 references are as `progress.md` states (reference pass counts 11/13 and 9/11 unchanged).
   - X7: check rounding in §4, the band wording, the QUALIFIED accumulation statement, and `dose_removal_time_d`. Recompute `eta_gap_inversion_misfit_MODEL_BASED` for one pulse and confirm the stated identity (ratio − MODEL-BASED η = gap, within the 0.1 d medians). Check the AE suffix, the DO = 0 against DO = 2 mg/L statement, and the quoted source caveat.

7. **Trace pulse 1 and pulse 2 end to end.** Pulse 1 changed through X6, pulse 2 through X5. Every C(t)-derived row must rest on one `nh4_series` identifier, or be labelled MODEL-BASED, SCENARIO or SENSITIVITY. Every `w2_bound_check` must be consistent with `nh4_truth_w2_consistency.csv`.

8. **Judge test adequacy.** The new tests are in `tests/test_pass3_fixes.py`, `tests/test_states.py` and `tests/test_report.py`. Check them against the prompt's positive, negative, boundary, regression, integration and end-to-end rows. For every existing test changed in pass 3, compare it with `git show 954f87a:<file>` and the rationale in `progress.md`, and decide whether it was changed because a requirement changed or was weakened. These are:
   - `test_truth`: pinned η branch; NH4-independent fact rows;
   - `test_rereview_fixes`: the R8 assertions inverted by X6;
   - `test_states`: the all-artifact interval;
   - `test_pulses` and `test_truth` fixtures;
   - the mechanical renames in `test_pulses`, `test_segments`, `test_thermo_m4`, `test_energy` and `test_report`.

## Rule explicitly on the implementer's judgment calls

1. **X1 start rule.** The kinetic start vector is scaled by 0.5, 1, 2 and 4 for every candidate, rather than an Alk0-only grid or a warm start from biomass_decay. The warm start was rejected on a diagnosis run.
2. **X1 budget.** 8,000 evaluations per start, fixed from the diagnosis (more than twice the 3,620 a converging run needed).
3. **X1 calibration choice.** The calibration is the lowest-RSS converged start, even when another converged start forecasts far better (alkalinity_limited).
4. **X1 post-hoc validity rule.** The Alk0 interval is withheld when the profile is multimodal below Δ = 3.84.
5. **X2 estimator.** ρ̂ is kept as the prompt specifies, without a small-sample bias correction, and the remaining undercoverage is documented.
6. **X2 upper bound.** The upper bound counts unclassified and artifact bins as possibly idle, but treats active-nitrification and oxygen-limited bins as not idle, so it is conditional on the active classification.
7. **X6.** R8 is extended symmetrically to windows governed by a pulse, with the result labelled "R8 CASE" rather than SENSITIVITY, because the checker reserves SENSITIVITY for an alternative observation model.

## Decision

Classify each remaining issue as P0, P1 or P2. Any of the following gives NO-GO:
- an unresolved P0 or P1;
- a failed or non-reproducible command;
- an unverified acceptance criterion of `docs/remediation_prompt_2026-10-19_pass3.md`;
- any of its failure conditions met;
- a rejected judgment call that you rate P0 or P1;
- an unresolved scientific or data-integrity risk;
- a weakened test.

Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT listing each P2 item. Justify the decision with evidence. Write the review to `docs/rereview_2026-10-20_pass3.md`, using the structure of `docs/rereview_2026-10-19_pass2.md`.
