# Independent re-review prompt for remediation pass 4 (2026-10-22)

You are an independent reviewer who did not implement the change. The previous re-review (`docs/rereview_2026-10-20_pass3.md`) decided GO WITH DOCUMENTED NON-BLOCKING DEBT and left six P2 items (P2-1 to P2-6) plus the carried W8 debt. The implementer then made remediation pass 4, committed as `c70cf2d` on branch `remediation/rereview-2026-10-19` (parent `083dde2`, whose code and outputs equal `ad51a6a`). Decide whether that pass closes the previous items without breaking what was accepted.

Treat everything in `docs/progress.md`, section "Remediation pass 4", as claims to verify, not as facts. Review the cumulative tree. Use `git diff 083dde2 c70cf2d` to see exactly what pass 4 changed.

Reference documents:
- `docs/remediation_prompt_2026-10-21_pass4.md`: work items Y1 to Y7, the Step 1 / Step 2 split, invariants, tests, acceptance criteria and failure conditions.
- `docs/rereview_2026-10-20_pass3.md`: the previous decision, P2-1 to P2-6, and the seven judgment calls of §10. Do not reopen those calls unless pass 4 changed their basis.
- `docs/remediation_prompt_2026-10-19_pass3.md` and `docs/remediation_prompt_2026-10-19.md`: their invariants, non-goals and acceptance criteria still apply.

## Procedure

1. **Run the commands yourself** on a clean checkout of `c70cf2d`, and record the full results. A failed command is NO-GO.
   - `uv run ruff check`
   - `uv run pytest -q`. The implementer reports `444 passed in 17257.89s (4:47:37)` on a heavily loaded machine. The suite contains two complete run-alls, and each run-all now includes three volume-scenario refits of every candidate.
   - `rm -rf outputs && uv run python -m bnr_neq run-all`. The implementer reports exit 0 after 4:23:58 under load. Run it in its own worktree, not at the same time as pytest. `BNR_NEQ_WORKERS` sets the process count; results must not depend on it.
   - Compare the regenerated `outputs/tables` (66 files), `outputs/reports` and the figures with a copy of the committed outputs taken before the run, using `diff -rq`. They must be byte-identical.
   - Optional determinism check: run with a different `BNR_NEQ_WORKERS` from the implementer's default and confirm byte identity. This also covers the scenario refits, which run concurrently in threads, each with its own process pool.

2. **Close the previous re-review item by item.** For P2-1 to P2-6 and the carried W8 debt in "Remaining issues by severity" of `docs/rereview_2026-10-20_pass3.md`, state resolved, partly resolved or unresolved, with file and line evidence. For each of Y1 to Y7, check every lettered sub-item of the pass-4 prompt.

3. **The Step 1 guard.** The prompt forbids any Step 1 change to a fitted value, a selection or a MODEL-BASED number.
   - Verify the claim that, against `ad51a6a`, `pulses.csv`, `nh4_truth_*.csv`, `model_identifiability.csv`, `model_checks.csv`, `model_response_fractions.csv`, `volume_constraints.csv`, `processed/` and every figure are byte-identical in the final outputs. Do this on your own regenerated outputs, not the committed ones.
   - Confirm that `null_no_consumption`, `biomass_decay` and `alkalinity_limited` are identical in every pre-existing column of `model_fit.csv` except the deduplicated `alk0_profile_status`.
   - Check that the regression test asserting this is real: it must compare against pass-3 values, not against values recomputed by the code under test.

4. **Y1 / Y2, the skill ranges and the Alk0 evidence.**
   - Recompute Δ and the test and forecast skill of the `alkalinity_limited` x0.5 start by hand or with your own script. The implementer reports 0.812 / 0.921, not the 0.817 / 0.925 the pass-3 review quoted. Confirm the reason: those were the values of the pass-2 fit, cut off at 1,300 evaluations.
   - Confirm that the per-start skills are DIAGNOSTIC only, and never choose a calibration or a model.
   - Judge the registered report rule for "holds" and "undetermined" (`progress.md`, Y1). Does report §1 point 7 now state the ranking correctly for every candidate of the other mechanism, including the two alkalinity × FA candidates that converge in pass 4?
   - Confirm that the Alk0 evidence range uses converged minima only, and that each distinct minimum appears once.

5. **Y3, the limits of the idle bounds.**
   - The new `states_interval_coverage.csv` gives production-interval coverage 0.853 at φ = 0.9 and **0.653** at φ = 0.98. The φ = 0.98 grid point was added in this pass, and its median ρ̂ (0.91) is the one nearest the record's median ρ̂ of C(t) (0.926). This is lower than the "roughly 0.8" in the pass-3 review.
   - Judge whether adding φ = 0.98 and choosing the nearest-ρ̂ grid point is a fixed rule or a choice made after seeing results.
   - Judge whether 0.65 coverage changes the pass-3 conclusion that the bounds stay usable (lower bound conservative, upper bound low by at most the classified non-idle share). Recompute that share from `states_fractions.csv`.
   - Confirm that report §1 point 8, §3 and `states.md` take the numbers from one function (`states.interval_limits_text`), so no number is written twice.

6. **Y4 / Y5, the checker and the tests.** Repeat the pass-3 probes (`docs/rereview_2026-10-20_pass3.md` §6) and the two that were not caught:
   - an NH4-derived metric declared `NH4-independent`;
   - a reading-derived column in a wide table with `M6_log` rows.

   Both must now raise. Try to break the new declarations:
   - a column declared C(t)-derived that is in fact reading-derived;
   - an `o2_supply_capacity_*` name that the kLa pattern wrongly accepts;
   - a wide `NH4-independent` table with an extra column.

   Check that the replacement `at_threshold` case really has a classified fraction of exactly 0.5.

7. **Y6, the Haldane reparametrization.** This is the main scientific change of the pass.
   - Verify the diagnosis: the stopping point, the condition numbers (8.1e5 overall, 2.7e3 for the kinetic block), the three-parameter valley in (log vmax, log Ks, −log Ki), and the correlations −0.99998 / +0.99993 / −0.99984. Reproduce at least the Jacobian condition number.
   - The implementer did not use the prompt's example (vmax, vmax/Ks) but (α = vmax/Ks, Ks, S_opt = √(Ks Ki)). The stated reason is that the example leaves corr(vmax, Ki) at 0.99998. Judge whether the choice was fixed from the diagnosis before any refit under it.
   - Confirm that r = α S/(1 + S/Ks + S²/S_opt²) is identical to the old rate law. The implementer reports agreement to 1e-12 relative; check that the test asserts this.
   - Confirm that the new parametrization was applied to every candidate that shares the Haldane pair (haldane_FA, combined, alkalinity_FA, each at T20C and T25C), and that the policy (four starts × 0.5/1/2/4 on the log kinetic block, 8,000 evaluations, default tolerances, lowest-RSS converged start) is unchanged.
   - **Bounds in the new coordinates.** The pass-3 haldane and combined optima sat on the 1e3 bound of vmax or v0. In the new coordinates that bound no longer caps vmax = α Ks, and the combined optimum now sits on the 1e3 bound of Ks instead. Judge whether "the same bound rule in new coordinates" treats the candidates equally, or whether it changes their effective parameter space. Check whether any reported skill or the selection depends on a bound-constrained optimum, and whether the report says so.
   - The alkalinity × FA candidates now converge (T20C 4 of 4 starts, T25C 3 of 4) with skills 0.792 / 0.642 / 0.802 and 0.793 / 0.644 / 0.803, and are not selected. Confirm that they enter the selection on equal terms, that their one-sided lower Alk0 limits (16.4 and 16.3 meq/L) are valid, and that the model status row may now read GO WITH DOCUMENTED NON-BLOCKING DEBT.
   - The M3 kinetic fits of `segments.py` are also Haldane in FA but were left unchanged. Judge whether that is consistent with "the same parametrization for every candidate that shares it".

8. **Y7, M6 per volume scenario.**
   - Check the scenario state model against the prompt and physics:
     - V(t) = V0 − E t;
     - dC/dt = −r + C E/V;
     - a pulse adds dose × V0/V(t_p);
     - the alkalinity state A = φ (Alk0 − 2n/M_N), with dn/dt = r/φ.
   - **The model also assumes that the biomass activity of every rate law scales with φ(t) = V0/V(t)** (biomass concentrates as water evaporates). Judge whether that assumption is stated, physically defensible, and was fixed before the refits.
   - Confirm that the model at E = 0 equals the constant-volume model (test), and that the constant-volume path itself is untouched.
   - Confirm that every scenario row carries SCENARIO and is never read by a primary row, and that `check_scenario_rows` raises before the tables are written.
   - **The scenario C(t) leaves the primary band.** At 0.5 and 1 mL/d the refitted C(t) ends at 388 and 570 mg N/L, above the constant-volume band (end 294, band 182-357). Those are the scenarios that `volume_constraints.csv` finds consistent with N conservation (≥ 0.314 mL/d, CONDITIONAL on C(t)), while constant volume is INCONSISTENT. Judge whether the primary C(t) band, and every primary number derived from it, now understates its uncertainty. Judge whether the report states this prominently enough (§1 points 2 and 9, the status table), or whether this is a scientific-integrity risk that needs more than a scenario table. The prompt forbade choosing a primary volume scenario in this pass. Do not demand one, but rule on what the report must say.
   - Confirm "the selected model does not change under any scenario" from `model_volume_scenario_summary.csv`.

9. **Trace pulse 1 and pulse 2 end to end** as in pass 3 §7. Every C(t)-derived row must rest on one `nh4_series` identifier, or be labelled MODEL-BASED, SCENARIO or SENSITIVITY. Every `w2_bound_check` must be consistent with `nh4_truth_w2_consistency.csv`. These rows should be byte-identical to pass 3: confirm it.

10. **Judge test adequacy.** The new tests are in `tests/test_pass4_fixes.py` (26 tests). Check them against the prompt's positive, negative, boundary, regression, integration and end-to-end rows. For every existing test changed in pass 4, compare it with `git show 083dde2:<file>` and the rationale in `progress.md`, and decide whether it was changed because a requirement changed or was weakened:
    - `test_states.test_conditional_idle_fraction_needs_the_coverage_but_the_bounds_do_not` (Y5);
    - `test_pass3_fixes.test_every_checked_table_passes_the_hardened_checker` (Y4 declared columns).

    Also judge the development change to the new test `test_scenario_refit_recovers_the_generating_parameters_better_than_constant_volume`: its margin went from "below a third" to "better, below 2 % and a lower RSS" after a first failure. Is this a new test fixed to match the prompt, or a margin tuned to pass?

## Rule explicitly on the implementer's judgment calls

1. **Y1 report rule.** "Holds" if the selected model's lowest selection skill exceeds the highest of every other candidate's starts inside Δ < 3.84. "Undetermined" if the test or forecast ranges overlap; better or worse only if they are disjoint.
2. **Y2 distinct minima.** Converged starts within 5 % in Alk0 (`MULTIMODAL_LOG_TOL`) are one minimum, represented by its lowest-Δ start.
3. **Y3 coverage grid.** φ = 0.98 added; the report quotes φ = 0.9 and the grid point whose median ρ̂ is nearest the record's.
4. **Y4 indirect dependence.** A column whose dependence on C(t) is indirect (a state, a run, a count of eligible bins) is declared C(t)-derived.
5. **Y6 parametrization.** (α, Ks, S_opt) instead of the prompt's (vmax, vmax/Ks), applied to every W5 Haldane candidate, with the M3 Haldane fits unchanged.
6. **Y6 bounds.** The same 1e-4 to 1e3 box in the new coordinates, which uncaps vmax and moves the combined optimum to the Ks bound.
7. **Y7 activity scaling.** The biomass activity scales with V0/V(t), and the dose and alkalinity are concentrated likewise.
8. **Y7 concurrency.** The three scenario evaluations run in threads; determinism is argued from byte-identical runs.
9. **Baseline pytest not rerun.** Replaced by a byte-identical regeneration of the baseline outputs.

## Decision

Classify each remaining issue as P0, P1 or P2. Any of the following gives NO-GO:
- an unresolved P0 or P1;
- a failed or non-reproducible command;
- an unverified acceptance criterion of `docs/remediation_prompt_2026-10-21_pass4.md`;
- any of its failure conditions met;
- a rejected judgment call that you rate P0 or P1;
- an unresolved scientific or data-integrity risk;
- a weakened test.

Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT listing each P2 item. Justify the decision with evidence. Write the review to `docs/rereview_2026-10-22_pass4.md`, using the structure of `docs/rereview_2026-10-20_pass3.md`.
