# Remediation prompt, pass 4: report corrections, checker gaps and the remaining model debt

You are the implementer for `bnr-non-equilibrium-2026`, branch `remediation/rereview-2026-10-19`, starting from commit `ad51a6a`. First read `docs/rereview_2026-10-20_pass3.md`, the latest decision: **GO WITH DOCUMENTED NON-BLOCKING DEBT**. Then read `docs/remediation_prompt_2026-10-19.md` and `docs/remediation_prompt_2026-10-19_pass3.md`. Their invariants, non-goals and acceptance criteria still apply.

W1 to W8 of the original plan are implemented and accepted (passes 1 to 3). **Do not rebuild or redesign W3, W4 or W5.** This pass fixes the six P2 items of the pass-3 review and then works off the remaining model debt.

Treat each finding as a claim to reconcile against the code, not as an established fact. Before you edit anything:

- Run the pipeline and the tests on a clean checkout.
- Verify every path, line and number cited below. Report any that differ.

Rules for the whole pass:

- Make the smallest coherent change for each item, and do not refactor unrelated code.
- Do not weaken or delete existing tests.
- Never fabricate a measurement, and never tune code to reproduce or avoid a number.
- Fix every rule or numerical policy you introduce, with its argument, in `docs/progress.md` before you compute the outputs it affects.
- You are not the approver; an independent reviewer decides.

The judgment calls ruled on in §10 of the pass-3 review were accepted. Do not revisit them:

- the common start rule;
- the 8,000-evaluation budget;
- the lowest-RSS calibration;
- the multimodality validity rule for the Alk0 interval;
- ρ̂ without bias correction;
- the upper idle bound treating active bins as not idle;
- the symmetric R8 rule.

## Step 1: fix the pass-3 findings first

Do these before Step 2, and record their output changes separately in `progress.md`. None of them may change a fitted value, a selection or a MODEL-BASED number. If one does, stop and report it.

| Item | Change | Affected files (verify) |
| --- | --- | --- |
| Y1 (P2-1, highest priority) | **Report both alkalinity fits.** `alkalinity_limited` has two converged calibration minima that the calibration window cannot distinguish (Δ = 3.0 < 3.84): Alk0 = 80 (selection / test / forecast skill −0.91 / 0.23 / 0.44) and Alk0 = 22.6 (0.934 / 0.817 / 0.925, better than `biomass_decay` on test and forecast). Report §1 point 7 and §4 quote only the first. (a) For every candidate with more than one converged start inside Δ < 3.84 of its calibration, report the range of its selection, test and forecast skills over those starts. Test skill is not in `model_calibration_starts.csv`, so compute it per start; it is a diagnostic only and is never used to choose a calibration. (b) State in §1 point 7 and §4 that the selection of `biomass_decay` holds against every such start, but that the out-of-sample ranking of decay against alkalinity limitation is undetermined. (c) Correct `progress.md`: the pass-2 skills were not withdrawn "because unconverged". They belong to a minimum that has now converged and cannot be rejected. | `model.py` (per-start diagnostics), `report.py:526-541`, `docs/progress.md` |
| Y2 (P2-2) | **Alk0 evidence from estimates only.** The status-table evidence quotes "fitted Alk0 … 23.6-80.0 meq/L", and the 23.6 is from an unconverged fit that the report itself calls "not an estimate". Build the range from converged candidates only, and name the unconverged ones separately. In `alk0_profile`, list each distinct minimum once: starts within `MULTIMODAL_LOG_TOL` of each other are one minimum (currently "Alk0 22.6 at Δ = 3.01" appears twice). | `report.py:324-327`, `model.py:748-757` |
| Y3 (P2-3) | **State the limits of the idle bounds where they are reported.** (a) The autocorrelation-scaled interval covers about 0.85, not 0.95, at φ = 0.9 on 101 points, roughly 0.8 at the record's median ρ̂ ≈ 0.93, so some truly flat windows are still called rising or falling. (b) Both bounds are conditional on the state classification: the lower bound counts only bins classified inactive, and the upper bound counts active and oxygen-limited bins as not idle. Add these two statements to report §1 point 8, §3 and `states.md`, with the numbers taken from the synthetic check, not hard-coded twice. | `report.py`, `states.py` (`write_report`, line 486) |
| Y4 (P2-4) | **Close the two checker gaps.** (a) `NH4-independent` is a self-declaration and needs an allow-list (`model.NH4_INDEPENDENT_METRICS`): `t_pulse_d`, `segment_end_d`, `median_pH`, `o2_supply_capacity_mgO2_SCENARIO_kLa*_T25C`, `free_energy_per_fully_converted_pulse_J_SCENARIO_T25C` and `stirrer_energy_J`. A metric outside the list that declares `NH4-independent` raises. (b) A wide table whose rows carry the primary model must also declare its NH4-independent columns, and every other column must be on a declared list of C(t)-derived columns. An undeclared column raises, as it already does for reading-based rows. Apply this to `nh4_truth_windows.csv` and the `states_*` tables. | `model.py:87-140`, `pulses.py:894-895`, `states.py:469` |
| Y5 (P2-6) | **Tests.** In `test_conditional_idle_fraction_needs_the_coverage_but_the_bounds_do_not`, replace the misnamed `at_threshold` case (`tests/test_states.py:219`, which is a second below-threshold case) with a true classified fraction of exactly 0.5. Add negative tests for both Y4 gaps: a derived metric declared `NH4-independent`, and an undeclared column in a primary wide table. | `tests/test_states.py`, `tests/test_pass3_fixes.py` or a new `tests/test_pass4_fixes.py` |

## Step 2: the remaining model debt

Start Step 2 only after Step 1 is complete and its outputs are explained.

| Item | Change | Affected files (verify) |
| --- | --- | --- |
| Y6 (P2-5) | **Settle the two alkalinity × FA candidates, or state why they cannot be settled.** All four starts stop at 8,000 evaluations, agreeing to 1e-5 in calibration RMSE near Alk0 = 23.6: a flat, correlated valley. (a) Diagnose which parameters are correlated at the stopping point (Jacobian condition number, parameter correlation), and record it in `progress.md`. (b) If a reparametrization removes the correlation (for example fitting vmax and vmax/Ks instead of vmax and Ks), fix it before refitting and apply the same convergence policy (starts, budget, tolerances). Any parametrization change must also be applied to every other candidate that shares the parameter pair, so the policy stays the same for every candidate, and the changes to those candidates must be explained. (c) If the candidates then converge, they enter the selection on equal terms. If they do not, they stay excluded and labelled, with the diagnosis as the stated reason. Either way the report's §1 point 7 follows the result. | `model.py`, `report.py`, `docs/progress.md` |
| Y7 (W8 debt) | **Refit M6 per volume scenario.** The constant-volume frame is INCONSISTENT with the volume constraints (`volume_constraints.csv`: N conservation needs a net loss ≥ 0.314 mL/d, CONDITIONAL on C(t)). Refit the selected model's observation and kinetic parameters under each SCENARIO evaporation rate already in the config (0.25, 0.5 and 1 mL/d), with the pulse dose and the states diluted or concentrated by V(t). Report, per scenario: the C(t) level against the constant-volume band, the selection skills, and whether the selected model changes. The primary C(t) and every primary row stay constant-volume. The scenario refits are labelled SCENARIO and never feed a primary row. If a scenario changes the selected model, say so in §1 and the status table. Do not choose a primary volume scenario in this pass. | `model.py`, `volume.py`, `report.py` |

## Invariants (unchanged, plus one new)

- One NH4 truth per run; primary rows rest on the one inverted C(t), and MODEL-BASED rows on the selected model state.
- Every C(t)-derived row of a W2-contradicted pulse is QUALIFIED.
- η ∈ [0, 1].
- Every scenario value carries its label.
- No window spans a pulse or an electrode jump.
- Consistency checks are not findings.
- Only converged candidates can be selected.
- Every row of a checked table declares its NH4 basis.
- No reading-based column or metric lacks SENSITIVITY in its name.
- No idle fraction is reported without its coverage.
- **New:** no `NH4-independent` row or column is outside its declared allow-list.

## Non-goals

- No change to W3, W4 or W5 beyond the items above; no new process states, energy metrics or candidate mechanisms.
- No change to the sensor verdicts, the S2 qualification logic, the W2 bound, the energy formulas, the AE values or their source.
- No change to the registered thresholds: 15 mV, the 1.01 d window, the block and ≥ 2 blocks rules, the 0.5 coverage threshold, the fragment rule, the convergence budget and the start factors.
- No new or assumed measurements. Items blocked by `lab_NH4_N_mg_L`, `initial_alkalinity_meq_L`, `DO_mg_L` or the vessel masses stay blocked.

## Tests

| Type | Required tests |
| --- | --- |
| Positive | Y1: on a synthetic record with two minima inside Δ < 3.84, the reported skill range spans both. Y4: every metric on the NH4-independent allow-list passes. Y7: on a synthetic record generated at a known evaporation rate, the scenario refit at that rate recovers the generating parameters better than the constant-volume fit. |
| Negative | Y4: a derived metric declared `NH4-independent` raises, and an undeclared column in a primary wide table raises. Y2: an unconverged candidate's Alk0 never enters the evidence range. Y7: a scenario refit row without SCENARIO in its label raises. |
| Boundary | Y5: a classified fraction of exactly 0.5 computes the conditional value. Y2: two starts within `MULTIMODAL_LOG_TOL` are reported as one minimum. |
| Regression | Existing tests still pass, and every changed reference is explained by Y1 to Y7. Step 1 changes no fitted value, selection or MODEL-BASED number (assert it against the pass-3 outputs). |
| Integration | On the real run-all outputs: §1 point 7 contains the skill range of `alkalinity_limited`; every checked table passes the extended checker; the volume-scenario table exists and its rows are SCENARIO-labelled. |
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

Then diff `outputs/tables`, `outputs/reports` and the SVG figures against a copy of the `ad51a6a` outputs taken before the run. Explain every changed value or column in `docs/progress.md`, under "Remediation pass 4", by its work item, with Step 1 and Step 2 kept apart. Run the full suite once on the final code, after the last edit. Do not regenerate single stages after the final full run.

## Acceptance criteria and failure conditions

Accept only if all of the following hold:

- Y1 to Y5 are implemented and tested.
- Y6 and Y7 are either implemented and tested, or carried as debt with a recorded diagnosis.
- All three commands pass on the final code.
- The outputs are reproducible.
- The report states the skill range of every candidate with indistinguishable minima, the limits of the idle bounds, and the volume-scenario refits as SCENARIO.

The work fails if any of the following occurs:

- a Step 1 item changes a fitted value, a selection or a MODEL-BASED number;
- a non-converged candidate's skills or Alk0 appear as an estimate;
- the convergence policy or a parametrization differs between candidates that share it;
- a scenario refit feeds a primary row, or appears without its SCENARIO label;
- an `NH4-independent` row or column outside the allow-list passes the checker;
- an existing test is weakened.

Report exactly what was changed, the commands run and their full results.
