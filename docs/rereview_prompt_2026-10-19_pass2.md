# Independent re-review prompt for remediation pass 2 (2026-10-19)

You are an independent reviewer who did not implement the change. The previous re-review (`docs/rereview_2026-10-19.md`) gave NO-GO. The implementer then made remediation pass 2, committed as `954f87a` on branch `remediation/rereview-2026-10-19`. Decide whether that pass resolves the previous decision.

Treat everything in `docs/progress.md`, sections "Remediation pass 2" and "SAE verification", as claims to verify, not as facts. Note that `954f87a` is the first commit after `eb74b94`, so it contains the whole project: the original implementation, pass 1 and pass 2. Review the cumulative state of the tree. Use `docs/progress.md` to locate what pass 2 changed.

Reference documents: `docs/review_2026-10-19.md` (findings R1 to R12), `docs/remediation_prompt_2026-10-19.md` (work items W1 to W8, tests, acceptance criteria) and `docs/rereview_2026-10-19.md` (the previous decision, S1 to S3, and the remaining-issues table).

## Procedure

1. **Run the commands yourself** on a clean checkout of `954f87a` and record the full results. A failed command is NO-GO.
   - `uv run ruff check`
   - `uv run pytest -q`. This takes about 60 to 80 min, because it contains two complete run-alls and the M6 fits.
   - `rm -rf outputs && uv run python -m bnr_neq run-all`. This takes about 25 min.
   - Compare the regenerated `outputs/tables` and `outputs/reports` with the committed ones using `diff -rq` against a copy taken before the run. The committed outputs must be reproducible. The SVG figures are also expected to be byte-identical.
2. **Close the previous re-review item by item.** For every row of the "Remaining issues by severity" table in `docs/rereview_2026-10-19.md`, and for R1 to R12 and W1 to W8, state resolved, partly resolved or unresolved, with file and line evidence. Pay particular attention to these:
   - **S1:** the ΔA in the status table must come from the primary (M6_log) rows of `pulses.csv`. No reading-based column may lack SENSITIVITY in its name. Check `thermo_affinity_steps.csv` and `segments.csv`.
   - **W1 provenance:** `check_single_observation_model` must now reject mixed series and MODEL-BASED mislabels, as well as mixed models. Try to construct an output that breaks the one-observation-model invariant and is not caught.
   - **S2:** every C(t)-derived row of a pulse whose checkpoints contradict the W2 bound must carry a `w2_bound_check` value starting with QUALIFIED. Verify this on the committed `pulses.csv` and `nh4_truth_windows.csv`. Recompute one Alk0 requirement in `nh4_truth_w2_consistency.csv` by hand.
   - **S3:** no computed η may lie outside [0, 1]. Check that `removal_to_dose_ratio` and η agree wherever η is computed.
   - **W8 / P2:** check the `top_up_mL` status, the `alkalinity_mmol_L` handling (`alk0_source`), and `volume_constraints.csv`.
3. **Trace pulse 2 end to end.** Its consumption, FA, η, O2 demand, energy metrics and idle fraction must rest on the same `nh4_series` identifier, or be labelled MODEL-BASED or SCENARIO. The `w2_bound_check` value must be consistent with `nh4_truth_w2_consistency.csv`.
4. **Recompute by hand**, checking units throughout:
   - the W2 positive-test values;
   - one aeration-equivalent energy, `aeration_energy_supply_kWh_per_kgN_SCENARIO_SAE1_kLa10_T25C` for pulse 1, from kLa, C*, V, Δt, N removed and SAE;
   - one demand floor (4.57/SAE);
   - one `relaxation_time_d`.
5. **Verify the SAE source independently.** The configuration now uses field (process-water) aeration efficiency of 0.3, 1.0 and 2.6 kg O2/kWh. It is cited to Stenstrom and Rosso (2010, UCLA/UCI "Aeration" notes, p. 56), with Moore (2021, ORNL training series no. 11) as a cross-check. Confirm that the values are process-water rather than clean-water values, and that the range matches the cited tables. Note that the config key keeps the prompt's name "SAE" for a process-water AE, and judge whether that naming is acceptable.
6. **Check the new modules and stages:** `states.py`, `energy.py`, `figures.py`, the W5 candidates in `model.py`, the R8 rule in `segments.py`, and the new `states` stage between `model` and `pulses`. Look for silent fallbacks: a missing input must give a status string, never a default number, and every scenario value must carry its label. Confirm that every new or renamed column is documented and that its readers (report, figures, tests) were updated.
7. **Judge test adequacy.** For each new module, check that positive, negative, boundary, integration and end-to-end tests exercise the production path. Check that tests changed in pass 2 (`test_main`, `test_report`, `test_segments`, `test_truth`) were changed because a requirement changed and were not weakened. Compare each edit against `git show 954f87a` and the rationale in `docs/progress.md`.
8. **Judge scientific correctness:**
   - W3 rules: were they fixed from stated physical and resolution arguments before any association was computed? Is the 15 mV ORP threshold defensible? Is excluding the NH4 slope from the association partition sufficient to avoid circularity? Does the block length (at least the 1.01 d slope window, and at least the first lag with ACF < 1/e) respect autocorrelation? Is the rule "each group must span at least 2 blocks", added after the first run, a validity condition rather than tuning?
   - Coverage: the classifier assigns a state to only 12 % of the record. Judge whether the state layer still supports the claims the report makes from it, including the idle fraction used in W4.
   - W5: biomass_decay is still selected, but `alkalinity_limited` has better test and forecast skill. Is the report's conclusion, that decay cannot be discriminated from alkalinity limitation, correct and adequately prominent?
   - R8: is the electrode-fragment rule (3 × the longest delayed electrode response) principled, and does it remove only artifacts?
   - Report: confirm that no consistency check (coexistence affinity, ΔA, σ ≥ 0) and no association is presented as a finding or a cause, and that every number in the summary (section 1 of `analysis_report.md`) traces to a table.

## Rule explicitly on the implementer's judgment calls

The implementer chose these instead of the literal wording of the prompt or the previous review. Accept or reject each one with reasons:

1. **W1, "r bounded by the alkalinity result of W2".** The inverted rate is flagged QUALIFIED where it contradicts the bound, rather than constrained to obey it. The stated reason is that forcing the bound would alter the data to fit an assumption.
2. **S3.** η is withheld above 1, and the removal-to-dose ratio is reported as a separate metric that may exceed 1 (carry-over TAN removed), rather than treating every ratio above 1 as physically impossible.
3. **Infeasible kLa scenarios.** Where the O2 required exceeds the transfer capacity, utilization and supply-based energy are reported as not computed, rather than as values.
4. **Conditional EC bound.** The flat conductivity record is reported as an evaporation bound conditional on probe validity, and is not used, because the probe failed M2.
5. **Remaining P2 debt.** The M6 fit is not refitted per evaporation scenario.

## Decision

Classify each remaining issue as P0, P1 or P2. Any of the following gives NO-GO:
- an unresolved P0 or P1;
- a failed or non-reproducible command;
- an unverified acceptance criterion of `docs/remediation_prompt_2026-10-19.md`;
- a rejected judgment call that you rate P0 or P1;
- an unresolved scientific or data-integrity risk;
- a weakened test.

Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT listing each P2 item. Justify the decision with evidence. Write the review to `docs/rereview_2026-10-19_pass2.md`, using the structure of `docs/rereview_2026-10-19.md`.
