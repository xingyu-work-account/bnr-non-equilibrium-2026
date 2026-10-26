# Remediation prompt, pass 3: model convergence, idle fraction and provenance hardening

You are the implementer for `bnr-non-equilibrium-2026`, branch `remediation/rereview-2026-10-19`, starting from commit `954f87a`. First read `docs/rereview_2026-10-19_pass2.md`, the latest decision (NO-GO). Then read `docs/remediation_prompt_2026-10-19.md`: its invariants, non-goals and acceptance criteria still apply.

Treat each finding as a claim to reconcile against the code, not as an established fact. Before you edit anything:

- Run the pipeline and the tests on a clean checkout.
- Verify every path, line, column and number cited below. Report any that differ.

Rules for the whole pass:

- Make the smallest coherent change for each item, and do not refactor unrelated code.
- Do not weaken or delete existing tests.
- Never fabricate a measurement, and never tune code to reproduce or avoid a number.
- Fix every threshold, rule or numerical policy you introduce in code, with its physical or statistical argument, before you compute the outputs it affects. Record it in `docs/progress.md` before the run.
- You are not the approver; an independent reviewer decides.

The five judgment calls ruled on in §9 of the review were accepted: QUALIFIED instead of constraining, η withheld above 1, infeasible kLa not computed, the conditional EC bound, and no M6 refit per volume scenario. Do not revisit them.

## Problem and root cause

Pass 2 closed S1, S2 and S3 and implemented W3 to W7. Two P1 defects block acceptance:

- **N1.** The W5 model comparison is invalid. All three alkalinity-limited candidates finish calibration with `cal_converged = False` (`model_fit.csv`; `model.py:384` returns `sol.success`). `biomass_decay` converged and is selected against them by minimum selection RMSE (`model.py:431`). The selected model, and so every MODEL-BASED quantity, rests on an unequal comparison. Neither the report nor `progress.md` mentions this, and the fitted Alk0 is not reported.
- **N2.** The idle fraction is a biased lower bound presented as a value. There are two causes:
  - Only 12 % of the record carries a state, and `idle_fraction_of_stirred_time` divides by all bins (`states.py:277`).
  - The per-window Theil-Sen intervals that decide "C(t) flat" (`states.py:153`, `sign_class` at `states.py:176`) treat 101 autocorrelated bins as independent. The median lag-1 autocorrelation of detrended C(t) inside the slope windows is 0.93, so n_eff ≈ 3.7. The intervals are therefore about 5× too narrow, and only 379 of 2,031 sloped bins are flat.

  So the inactive and oxygen-limited states are under-assigned by construction. The report nevertheless quotes the idle fraction without caveat in §1 point 8, the §2 energy table, Fig. 8c and the W3 status row (GO WITH DEBT).

The root causes are (1) a fitting procedure without a uniform convergence policy, and no convergence gate on the selection, and (2) significance tests that ignore the autocorrelation the pipeline measures elsewhere: `model.effective_n` already does this for the profile likelihood.

## Work items

| Item | Change | Affected files (verify) |
| --- | --- | --- |
| X1 (N1, P1) | **Converge every candidate under one policy, and gate the selection on convergence.** (a) Diagnose why each W5 calibration stops: the `least_squares` status, `nfev` against the default limit, parameters at bounds, and the Alk0 start of 10 meq/L against the S2 requirement of about 20 meq/L. Record this in `progress.md`. (b) Fix one convergence policy for **all** candidates before refitting: the same evaluation budget, and the same multi-start or warm-start rule, for example starts on a fixed Alk0 grid of 5, 10, 20 and 40 meq/L, or warm starts from the converged `biomass_decay` observation parameters. No candidate may get a budget the others do not. (c) A candidate whose calibration still fails gets `selected = False` and a `selection_status` of `not computed: calibration did not converge (status …)`. It is excluded from selection, but its skills stay in the table, labelled. (d) Report `cal_converged`, the fitted Alk0 of each alkalinity candidate, and Alk0's profile interval or `not identifiable` in `model_fit.csv` (or a new `model_parameters_by_candidate.csv`), in report §4 and in §1 point 7. Compare the fitted Alk0 with the S2 requirements in `nh4_truth_w2_consistency.csv` and say whether they agree. (e) If the selected model changes, regenerate every MODEL-BASED quantity, and explain each changed value in `progress.md`. The report's conclusion must follow the converged comparison, whatever it is. | `model.py`, `report.py`, `pulses.py` (MODEL-BASED rows), `docs/progress.md` |
| X2 (N2, P1) | **Make the state slopes autocorrelation-aware, and report the idle fraction as a bounded quantity.** (a) In `rolling_theil_sen`, widen each window's interval for autocorrelation with a method fixed before rerunning. Preferred: scale the half-widths by √(n/n_eff), with n_eff = n (1 − ρ)/(1 + ρ) from the lag-1 autocorrelation ρ of the window's detrended residuals, the same formula as `model.effective_n`. A block-based interval is acceptable if justified. Apply it to the C(t), pH and ORP slopes alike. State the method in `states_rules.csv` and the `states.py` docstring. (b) Recompute the states and comparisons. Do not change the 15 mV threshold, the 1.01 d window, the block length rule or the ≥ 2 blocks rule. (c) Replace the single idle value with a bounded pair: `idle_fraction_lower_bound` = inactive bins / all bins, and `idle_fraction_upper_bound` = (inactive + unclassified + artifact bins) / all bins. Keep `idle_fraction_of_classified_time`, labelled CONDITIONAL on the classified bins being representative, with `classified_fraction` beside it. In `pulses.csv` the metric names must carry LOWER_BOUND, UPPER_BOUND and CONDITIONAL. (d) Fix a minimum coverage, with its reason, before the run, for example classified fraction ≥ 0.5 of the interval. Below it, the conditional idle fraction is `not computed: classified fraction … below …`; the bounds are still reported. (e) Report §1 point 8, the §2 table, Fig. 8c and the W3 status row must show the bounds and the coverage. The W3 status must not read GO for the idle fraction. | `states.py`, `pulses.py`, `report.py`, `figures.py` |
| X3 (P2) | **Harden the one-observation-model checker.** (a) Classify rows positively: every row of a checked table declares an NH4 basis (`M6_log`, a `_SENSITIVITY` model, `electrode_reading`) or an explicit `NH4-independent` marker. An empty model cell raises. (b) Any metric computed from removal, rate or concentration that is labelled `electrode_reading` raises, including `consumption_*` (other than `consumption_uncorrected_*`), `n_removed`, `removal_to_dose_ratio`, `o2_*`, `aeration_energy_*`, `gibbs_*`, `idle_fraction*`, `relaxation_time*` and `alk0_required*`. Prefer an allow-list of the measured raw metrics over a deny-list of derived prefixes. (c) Match the MODEL-BASED label as a token (`MODEL-BASED` or `MODEL_BASED` not preceded by "not"), so "not MODEL-BASED" is no longer a false positive. (d) Run the checker before writing `nh4_truth_windows.csv`, the `states_*` tables that carry C(t) quantities, `segments.csv` and `thermo_oxygen.csv`. Give each of them the provenance columns if it lacks them. | `model.py`, `pulses.py`, `states.py`, `segments.py`, `thermo.py` |
| X4 (P2) | **Put SENSITIVITY in every reading-based column name.** `segments.csv` `consumption_gain_corrected_*` becomes `consumption_SENSITIVITY_linear_gain_*`. In `thermo_oxygen.csv`, `consumption_gain_corrected_mgN_L_d`, `o2_demand_mgO2_L_d_GAIN_CORRECTED` and `kla_required_per_d_SCENARIO_T25C` become `*_SENSITIVITY_linear_gain_*`. That `kla_required_per_d_SCENARIO_T25C` column currently shares the primary metric's name. Update every reader: the report, figures, `pulses.py` and the tests. | `segments.py`, `thermo.py`, `pulses.py`, `report.py`, tests |
| X5 (P2) | **Clean up η and the W2 labels.** (a) An η row must not carry a band edge outside [0, 1]: pulse 2 has `param_band_high` = 1.069. Keep the full band on the `removal_to_dose_ratio` row, and give η the band clipped to [0, 1] with a note, or no band and a pointer to the ratio. (b) Fill `w2_bound_check` for the pre-pulse windows of `nh4_truth_windows.csv` from the checkpoint that bounds them (pre_pulse1), instead of leaving it empty. | `pulses.py` |
| X6 (P2) | **State the scope of the R8 rule.** Sub-segment 14 (8.96 to 9.27 d, 0.31 d) is retained although it starts 0.10 d after the delayed pulse-1 response jump, inside 3τ (M6 τ = 0.26 d). It supplies the highest point in Fig. 4b (18.3 mg N/L/d). Give the physical argument for excluding post-pulse windows from the rule (`segments.py:188-191`), or extend the rule symmetrically, before inspecting the result. Either way, report pulse 1's consumption with and without sub-segment 14 as a labelled sensitivity. | `segments.py`, `pulses.py`, `report.py` |
| X7 (P2) | **Fix the report text and naming.** (a) Round the M6 consistency-check and identifiability tables in report §4. (b) §1 point 2: the C(t) band is open upward over the whole series, not only at the end. (c) §1 point 5: label the accumulation statement as resting on the QUALIFIED C(t). (d) Rename `relaxation_time_d` to a name that says what it is (for example `dose_removal_time_d`, dose / post-pulse rate), or keep it with an explicit note that it is not a fitted relaxation constant. (e) Discuss why the MODEL-BASED η is about half the dose-anchored η for every pulse (0.21 against 0.38, 0.57 against 0.81, 0.035 against 0.20). (f) Aeration efficiency: rename the metric suffix `SAE` to `AE` (field aeration efficiency), or document the mismatch in the column descriptions. State that the supply-based energy divides a DO = 0 capacity by an AE defined at DO = 2 mg/L. Quote the source's caveat that its table is not a design guide (Stenstrom and Rosso 2010, p. 55). | `report.py`, `pulses.py`, `energy.py`, `config/experiment.toml` |
| X8 (P2) | **Close the test gaps listed in review §7.** See Tests. | `tests/` |

## Equations, units and invariants

- **Autocorrelation-scaled interval:** half-width × √(n/n_eff), with n_eff = n (1 − ρ)/(1 + ρ), ρ ∈ [0, 0.999] the lag-1 autocorrelation of the detrended window residuals, and n = 101 bins.
- **Idle bounds:** the lower bound is N_inactive / N_all. The upper bound is (N_inactive + N_unclassified + N_artifact) / N_all. The conditional value is N_inactive / N_classified.
- **Invariants from earlier passes:**
  - one NH4 truth per run;
  - every C(t)-derived row of a W2-contradicted pulse is QUALIFIED;
  - η ∈ [0, 1];
  - every scenario value carries its label;
  - no window spans a pulse or an electrode jump;
  - consistency checks are not findings.
- **New invariants:**
  - only converged candidates can be selected;
  - every row of a checked table declares its NH4 basis;
  - no reading-based column or metric lacks SENSITIVITY in its name;
  - no idle fraction is reported without its coverage.

## Non-goals

- No change to the sensor verdicts.
- No new or assumed measurements.
- No M6 refit per volume scenario (it stays stated debt).
- No change to the S2 qualification logic, the W2 bound, the energy formulas, or the SAE/AE values and source.
- No tuning of the 15 mV threshold, the slope window, the block rule or the fragment threshold.

## Tests

| Type | Required tests |
| --- | --- |
| Positive | On a synthetic model with known parameters, every candidate converges under the common policy and recovers Alk0 within a stated tolerance. On a synthetic AR(1) series with φ = 0.9, the autocorrelation-scaled slope interval covers the true slope in about 95 % of windows (fixed seed); the unscaled interval does not. The idle bounds equal hand counts on a small synthetic state vector. |
| Negative | A candidate with `cal_converged = False` cannot be selected, and its row carries the not-computed selection status. The checker raises on: an empty model on an NH4-derived metric; `aeration_energy_*`, `n_removed_mgN`, `removal_to_dose_ratio` or `consumption_weighted_mgN_L_d` labelled `electrode_reading`; and an unlabelled gain-corrected column. "not MODEL-BASED" in a basis does not raise. |
| Boundary | An interval whose classified fraction is below the coverage threshold yields `not computed` for the conditional idle fraction, while still reporting both bounds. ρ ≥ 0.999 caps n_eff and does not divide by zero. An η band edge above 1 is clipped, or moved to the ratio row. |
| Regression | Existing tests still pass. Any changed reference is explained by X1 to X7. Pin the branch in `test_truth.test_dose_anchored_eta_and_affinity` (assert the branch the synthetic truth must take). Extend `test_reading_based_columns_say_sensitivity` to the gain columns. |
| Integration | On the real run-all outputs: every row of `model_fit.csv` with `selected = True` has `cal_converged = True`; the fitted Alk0 of each alkalinity candidate is reported; the pulse-2 idle metrics carry their bounds and coverage; every checked table passes the hardened checker. |
| End-to-end | `uv run python -m bnr_neq run-all` regenerates all outputs from the raw files with fixed seeds. Two runs are byte-identical in CSV and SVG (the existing determinism test). |

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

Then diff `outputs/tables`, `outputs/reports` and the SVG figures against a copy of the `954f87a` outputs taken before the run. Explain every changed value or column in `docs/progress.md`, under "Remediation pass 3", by its work item. Run the full suite once on the final code, after the last edit. Do not regenerate single stages after the final full run.

## Acceptance criteria and failure conditions

Accept only if all of the following hold:

- N1 and N2 are resolved with evidence in the code and outputs.
- X3 to X8 are implemented and tested.
- All three commands pass on the final code.
- The outputs are reproducible.
- The report states the converged model comparison, the fitted Alk0 against the S2 requirement, and the idle fraction as bounds with coverage.

The work fails if any of the following occurs:

- a non-converged candidate can be selected, or its skills appear without a convergence label;
- the convergence policy differs between candidates;
- the slope-interval method or the coverage threshold is chosen after looking at its effect on the idle fraction;
- any idle fraction appears without its coverage;
- any reading-based column or metric lacks SENSITIVITY in its name;
- an output mixes observation models or series and is not caught;
- an existing test is weakened.

Report exactly what was changed, the commands run and their full results.
