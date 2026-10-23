# Independent re-review of remediation pass 4 (2026-10-22)

This review covers commit `c70cf2d` on branch `remediation/rereview-2026-10-19` (parent `083dde2`, whose code and outputs equal `ad51a6a`). The claims in `docs/progress.md`, section "Remediation pass 4", were treated as claims to verify. The procedure is `docs/rereview_prompt_2026-10-22_pass4.md`. The reference documents are `docs/remediation_prompt_2026-10-21_pass4.md` (Y1 to Y7), `docs/rereview_2026-10-20_pass3.md` (P2-1 to P2-6, the carried W8 debt, and the judgment calls of §10) and the earlier remediation prompts.

Every command was run on clean detached worktrees of `c70cf2d`: one for run-all, one for pytest and one for scratch scripts. A detached worktree of `083dde2` was used to reproduce the pre-pass-4 diagnosis. The untracked `src/bnr_non_equilibrium_2026/` stub and `docs/rereview_prompt_2026-10-22_pass4.md` are not part of the commit and were not present in the worktrees.

## Decision

**GO WITH DOCUMENTED NON-BLOCKING DEBT.** Pass 4 closes every item the pass-3 review left open, including the W8 refit debt carried since pass 2:

- Step 1 (Y1 to Y5) changed no fitted value, selection or MODEL-BASED number. Every primary table is byte-identical to `ad51a6a`.
- The report now states that the out-of-sample ranking of decay against alkalinity limitation is undetermined.
- Both checker gaps raise.
- The alkalinity × FA candidates converge after a reparametrization that the reviewer reproduced from its diagnosis. They are not selected.
- The selection is unchanged under every volume scenario.

All three commands pass on a clean checkout, and the outputs regenerate byte-identically at a different worker count. No acceptance criterion fails, no failure condition is met, no test is weakened, and all nine judgment calls are accepted.

The question the review prompt raised about Y7, whether the scenario C(t) "above the band" means the primary band understates its uncertainty, resolves without a P1. The scenario C(t) is the constant-volume C(t) times V0/V(t) to within 1-2 %, so amounts and η are unchanged.

Five P2 items remain (last table). The most useful to fix next are P2-A, which explains that ratio in the report, and P2-B, which flags the bound-constrained Haldane calibrations.

## 1. Commands (run by the reviewer, 2026-10-22)

| Command | Result |
| --- | --- |
| `uv run ruff check` | `All checks passed!`, exit 0 |
| `rm -rf outputs && uv run python -m bnr_neq run-all` (`BNR_NEQ_WORKERS=6`) | Exit 0, wall 1:31:48 on a mostly idle machine. Stage messages: verdicts NH4/EC/ORP FAIL, pH PASS; segments references 11/13; thermo 9/11; dynamics 28/28; model selected `biomass_decay`, 13/15 identifiable; states 12 % classified, idle between 5 and 93 %; pulses 274 computed, 78 not computed. These match `progress.md`. |
| `diff -rq` against a copy of the committed `outputs/` taken before the run | No differences in `tables/` (66 files), `reports/` or `figures/` (PNG and SVG). `processed/` is byte-identical to the implementer's final run in the main checkout. |
| Determinism across worker counts | The reviewer used 6 workers and the implementer the default (all cores), and the outputs are byte-identical. This includes the three volume-scenario evaluations, which run concurrently in threads, each with its own process pool. |
| `uv run pytest -q` (`BNR_NEQ_WORKERS=6`) | `444 passed in 20555.92s (5:42:35)`, exit 0, no skips (run with `-rs`). The pass-3 regression test, which reads `ad51a6a` through git, therefore ran rather than skipped. This matches the implementer's count (418 + 26). The run is after run-all, not concurrent with it. A first pytest run was interrupted at about 25 % when the reviewer's session ended, with no failures up to that point. It was rerun in full, and only the complete run is reported. |

Environment: CPython 3.14.8, uv 0.13.1 (warns that `uv-build<0.13.0` does not contain it; harmless), numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2. Scratch scripts (`y6_old.py`, `y6_new.py`, `y4_probe.py`) ran in separate worktrees.

## 2. Previous remaining issues

### 2.1 "Remaining issues by severity" of `docs/rereview_2026-10-20_pass3.md`

| Previous issue | Status | Evidence |
| --- | --- | --- |
| P2-1: the report showed `alkalinity_limited` only at its worse calibration | **Resolved** | `model_calibration_starts.csv` has Δ, `inside_95`, and test and forecast skills per start. These are DIAGNOSTIC columns that are never used to choose a calibration or a model (`model.py:725-790`). `model_fit.csv` has the skill range over the converged starts inside Δ < 3.84. Report §1 point 7: "the out-of-sample ranking of biomass_decay against alkalinity_limited is undetermined (test 0.747 to 0.748 against 0.227 to 0.812; forecast 0.862 to 0.863 against 0.443 to 0.921)". The pass-3 "withdrawn" sentence is corrected in `progress.md`. The reviewer reproduced the x0.5 start: test 0.8124, forecast 0.9214 (§4). |
| P2-2: Alk0 range included a non-estimate; duplicate minimum | **Resolved** | The status row reads "Alk0 of the converged alkalinity candidates 22.6-80.0 meq/L (every converged minimum inside Δ < 3.84)". `alk0_profile_status` lists "Alk0 22.6 meq/L at Δ = 3.01 (2 starts)" once, via `distinct_minima` (`model.py:928`). |
| P2-3: limits of the idle bounds stated only in `progress.md` | **Resolved** | New `states_interval_coverage.csv`. One function, `states.interval_limits_text` (`states.py:244`), writes the limits into report §1 point 8, §3 and `states.md`. §5 has the details. |
| P2-4: two checker gaps | **Resolved** | `NH4_INDEPENDENT_METRICS` (`model.py:75-79`). Wide tables with primary rows declare independent and C(t)-derived columns (`model.py:147-167`; allow-list check at lines 147-151). Both pass-3 probes now raise (§6). |
| P2-5: alkalinity × FA candidates unconverged | **Resolved** | After the reparametrization both converge (T20C 4 of 4 starts, T25C 3 of 4) and are not selected. §7 has the details, including a new P2 on parameter bounds. |
| P2-6: misnamed boundary test; no tests for the checker gaps | **Resolved** | `tests/test_states.py:219-224`: classified fraction exactly 0.5, conditional value 2/5. New negative tests `test_a_derived_metric_declared_nh4_independent_raises` (7 cases) and `test_an_undeclared_column_in_a_primary_wide_table_raises`. |
| Carried W8 debt: M6 not refitted per volume scenario | **Resolved** (P2 wording remnant, §8) | `simulate_volume` (`model.py:413`), `volume_scenario_refits` (`model.py:1198`), `model_volume_scenarios.csv` and `model_volume_scenario_summary.csv`. The selected model is unchanged under 0.25, 0.5 and 1 mL/d. |

### 2.2 Work items Y1 to Y7

| Item | Sub-item | Status | Evidence |
| --- | --- | --- | --- |
| Y1 | (a) skill range over starts inside Δ < 3.84 | Verified | `skill_range`, `start_delta`, `calibration_n_eff`. The chosen start of every candidate reproduces its `model_fit.csv` skills to 1e-12 (test `test_real_report_states_the_skill_range_and_the_limits`). |
| Y1 | (b) "holds" and "undetermined" in §1 point 7 and §4 | Verified | `selection_robustness_txt`, `ranking` (`report.py:643-675`). Decay's lowest selection skill, 0.958, exceeds the highest of every other candidate's starts inside the region, 0.935. The ranking against `alkalinity_limited` is undetermined, and decay is better than both alkalinity × FA candidates over all their starts (disjoint ranges). |
| Y1 | (c) `progress.md` corrected | Verified | `progress.md` line 969. |
| Y2 | converged-only range; one entry per minimum | Verified | §2.1. |
| Y3 | (a), (b) stated where the bounds are reported | Verified | §5. |
| Y4 | (a) allow-list; (b) declared columns | Verified | §6. |
| Y5 | tests | Verified | §2.1, §9. |
| Y6 | (a) diagnosis; (b) reparametrization fixed before the refits; (c) consequence | Verified, with P2-B and P2-C | §7. |
| Y7 | per-scenario refit, SCENARIO-only, selection reported | Verified, with P2-A | §8. |

## 3. Acceptance criteria and failure conditions of the pass-4 prompt

| Criterion / failure condition | Result |
| --- | --- |
| Y1 to Y5 implemented and tested | Met. |
| Y6 and Y7 implemented and tested, or carried as debt with a diagnosis | Met (both implemented). |
| All three commands pass on the final code | Met (§1). |
| Outputs reproducible | Met: byte-identical regeneration at a different worker count. |
| Report states the skill ranges, the limits of the idle bounds, and the scenario refits as SCENARIO | Met. The meaning of "above the band" is not explained (P2-A). |
| A Step 1 item changes a fitted value, a selection or a MODEL-BASED number | Not met. Against `ad51a6a`, `pulses.csv`, every `nh4_truth_*.csv`, `model_identifiability.csv`, `model_checks.csv`, `model_response_fractions.csv`, `volume_constraints.csv` and every figure are byte-identical in the reviewer's regenerated outputs. `null_no_consumption`, `biomass_decay` and `alkalinity_limited` are identical in every pre-existing column of `model_fit.csv` except the deduplicated `alk0_profile_status`. `states_fractions.csv` is identical in every old column. |
| A non-converged candidate's skills or Alk0 appear as an estimate | Not met. Every candidate now converges, and `skill_range` returns `not computed` for an unconverged one. |
| The convergence policy or a parametrization differs between candidates that share it | Not met. All six M6 Haldane candidates use (α, Ks, S_opt) with the same starts, budget and tolerances (test asserts this). The M3 kinetic fits are not M6 candidates (ruling 5). |
| A scenario refit feeds a primary row, or lacks its SCENARIO label | Not met. `check_scenario_rows` runs before the writes. No primary table names a scenario refit (test), and the primary tables are byte-identical to pass 3. |
| An NH4-independent row or column outside the allow-list passes | Not met (§6). |
| An existing test weakened | Not met (§9). |

## 4. Y1 / Y2: skill ranges and Alk0 evidence

- **Reproduced.** The reviewer fitted `alkalinity_limited` from the x0.5 start with the committed code: test skill **0.8124**, forecast skill **0.9214**, against 0.812 and 0.921 in the table. The pass-3 review's 0.817 and 0.925 were the values of the pass-2 fit cut off at 1,300 evaluations, as `progress.md` says. Δ = 3.006 was reproduced in pass 3, and the input is unchanged.
- **Diagnostic only.** The per-start skills enter only `skill_range` and the report text. `best_start` still picks the lowest-RSS converged start, and `select_converged` still selects on the calibrations' own selection RMSE.
- **Report rule (ruling 1).**
  - "Holds" compares the selected model's worst start against every other candidate's best start inside the region. That is conservative.
  - "Undetermined" whenever the test or forecast ranges overlap, and "better" only when both are disjoint. That is symmetric.
  - Applied to the outputs: decay vs `alkalinity_limited` is undetermined. Decay vs alkalinity_FA T20C and T25C is better on both, with disjoint ranges: test 0.747-0.748 against ≤ 0.646, forecast 0.862-0.863 against ≤ 0.804.
  - The concluding sentence, "the data therefore do not discriminate activity decay from substrate-side (alkalinity) limitation", follows from the `alkalinity_limited` result.
- **Readability (P2-D).** §1 point 7 is now one paragraph of about 600 words, in which three candidate ranges and the Alk0 profiles are written inline. The content is correct, but the summary has lost its "one page" character. A three-row table in §4 with a two-sentence §1 point would read faster.

## 5. Y3: the limits of the idle bounds

- **`states_interval_coverage.csv`** (production `rolling_theil_sen`, 600 synthetic windows per φ; values as committed and regenerated):

  | φ | Unscaled | Production (ρ̂) | True ρ | Median ρ̂ |
  | --- | --- | --- | --- | --- |
  | 0.50 | 0.750 | 0.932 | 0.953 | 0.46 |
  | 0.90 | 0.365 | 0.853 | 0.943 | 0.85 |
  | 0.95 | 0.272 | 0.797 | 0.945 | 0.88 |
  | 0.98 | 0.198 | **0.653** | 0.968 | 0.91 |

- **Grid and reference point (ruling 3).** φ = 0.98 was added so that the grid spans the record's median ρ̂ (0.926, known since pass 3), and the rule "quote φ = 0.9 and the grid point whose median ρ̂ is nearest the record's" is recorded under "Step 1 decisions". This is a design rule fixed from a known input, not from the result. The result (0.65) is less favourable than the pass-3 estimate of about 0.8, which argues against tuning. **Accepted.** The pass-3 review's "roughly 0.8" was an extrapolation, and it is superseded by this table.
- **Does 0.65 change the pass-3 conclusion?** No.
  - Undercoverage can only make truly flat windows read as rising (→ `unclassified`, still counted in the upper bound) or falling (→ `active_nitrification`, not counted).
  - So the upper bound is low by at most the share of bins classified active or oxygen-limited: 0.0558 + 0.0134 = **0.069** over the whole record, as the report says (6.9 %). The argument depends on which bins are excluded, not on how often the interval covers.
  - The lower bound stays conservative.
  - The bounds themselves (0.05-0.93) are unchanged.
- **One source.** Report §1 point 8, §3 and `states.md` contain the same sentence from `interval_limits_text`, and the test asserts the phrase in all three.

## 6. Y4: the checker

The pass-3 probes were rerun on the committed tables with the new code (scratch `y4_probe.py`):

| Probe | Result |
| --- | --- |
| Empty model; unknown `M7`; `consumption_weighted_mgN_L_d`, `aeration_energy_supply_…`, `n_removed_mgN` as `electrode_reading`; series on an NH4-independent row; second `C_inv:` series; unlabelled gain column in `segments.csv` | All caught (unchanged from pass 3) |
| **Pass-3 hole 1:** `consumption_weighted_mgN_L_d` declared `NH4-independent`, no series | **Caught:** "NH4-independent metrics outside the allow-list" |
| **Pass-3 hole 2:** reading-derived column in `nh4_truth_windows.csv` (M6 rows); the same in `states_fractions.csv` | **Caught:** "undeclared columns in a wide table with primary rows" |
| Allow-list pattern: `…kLa10_T25C_from_C`, `…kLa-1_T25C`, `…kLa_T25C`, `t_pulse_d_x`, `"median_pH "` | Caught |
| Allow-list pattern: `…kLa1e3_T25C` | Passes, correctly: a numeric kLa, as the rule intends |
| NH4-independent wide table (`states_interval_coverage.csv`) with an extra `C_inv_mgN_L` column | Caught |
| `nh4_truth_windows.csv` checked without declarations | Caught: "must declare its NH4-independent and its C(t)-derived columns" |
| A column declared C(t)-derived but filled with reading values | Passes. This is inherent: the checker vets names and declarations, not the arithmetic behind them. The declarations (`segments.truth_window_columns`, `states.STATE_TABLE_COLUMNS`) are short, explicit and reviewed here. Not a defect. |

- The declared lists match the tables: the committed tables pass, and the integration test also asserts that the `states_*` tables raise without declarations.
- Ruling 4 (indirect dependence declared C(t)-derived) is the safe direction. **Accepted.**

## 7. Y6: the Haldane reparametrization

**Diagnosis reproduced.** The reviewer fitted `alkalinity_FA_SCENARIO_T25C` from start x2 with the 8,000-evaluation budget on the old code (`083dde2`):
- status 0 after 8,000 evaluations, RSS **29.03025**;
- vmax 281.8, Ks 5.287, Ki 0.0452, Alk0 23.60, C0 77.05, no parameter at a bound;
- finite-difference Jacobian: condition number **8.07e5**, kinetic block **2.67e3**, smallest singular value 0.0031 against 0.348 for the next;
- corr(vmax, Ks) = **+0.99993**, corr(vmax, Ki) = **−0.99998**, corr(Ks, Ki) = **−0.99984**.

Every number matches `progress.md`.

**The rate law is unchanged.**
- r = α S/(1 + S/Ks + S²/S_opt²), with α = vmax/Ks and S_opt² = Ks Ki, multiplied through, is vmax S/(Ks + S + S²/Ki).
- `test_haldane_reparametrization_is_the_same_rate_law_for_every_haldane_candidate` simulates every Haldane candidate both ways and asserts agreement to rtol 1e-12. It also asserts the transformed default start and the × 0.5/1/2/4 policy starts.

**Choice (ruling 5).**
- The prompt's example (vmax, vmax/Ks) leaves corr(vmax, Ki) at 0.99998 by the implementer's own Jacobian, so rejecting it was a reason, not a preference.
- The chosen coordinates are the two quantities a Haldane curve in its inhibition branch determines: the initial slope and the substrate of maximum rate. They were registered before any refit under them.
- At the new alkalinity_FA T25C optimum the reviewer finds the largest kinetic correlation |corr(Ks, S_opt)| = **0.892**, which matches the claimed 0.90.
- The M3 kinetic fits in `segments.py` are a separate descriptive analysis, not candidates in the M6 comparison, so leaving them unchanged does not create unequal treatment within the comparison.
- **Accepted.**

**What the reparametrization did not do (P2-C).** It removed the near-perfect correlations but did not improve conditioning:
- At the converged alkalinity_FA T25C optimum (x0.5, 5,094 evaluations, status 2, RMSE 0.110916 as in the table), the Jacobian condition number is **4.2e6**, the kinetic block **1.4e4**, and the smallest singular value 0.0006.
- Convergence is by ftol on a still-flat valley. That is legitimate under the registered policy, but the natural-unit parameters are not estimates: vmax 2,601 mg N/L/d and Ki 0.0047 mg NH3/L, far from any physical range.
- The report quotes only skills and Alk0 for these candidates, which is correct. `progress.md` should not describe the result as resolving the ill-conditioning.

**Bounds in the new coordinates (ruling 6, P2-B).** The reviewer refitted two chosen calibrations under the committed code:
- `haldane_FA_SCENARIO_T20C` (x0.5): Ks = **991.7**, within 0.01 log units of the 1e3 bound; implied vmax = **9.4e4** mg N/L/d and Ki = **1.8e-4** mg NH3/L; observation intercept `a` on its bound.
- `combined_SCENARIO_T25C` (x4): Ks = **1000** (on the bound), and C0 and `a` also on their bounds.

`progress.md` mentions "Ks 300 to 990 for haldane", but neither report says that these calibrations are bound-constrained, or that their natural parameters are outside any physical range.

The rule itself (the same 1e-4 to 1e3 box on every log kinetic parameter) is applied uniformly, and the box was never physically motivated. Its natural-unit image for the Haldane candidates did change: vmax is no longer capped. No selection depends on it: the Haldane and combined candidates score 0.836-0.838 and 0.526-0.530 against decay's 0.958, before and after.

**Accepted, with P2-B:** `model.md` and §4 should flag bound-constrained calibrations (a `cal_at_bounds` column is the obvious form) and say that their natural parameters are not estimates.

**Consequence (Y6 c).**
- The two alkalinity × FA candidates converge and enter the selection on equal terms ("candidate (converged, not selected)").
- Their one-sided lower Alk0 limits, 16.4 and 16.3 meq/L, come from profiles that close on the low side and stay open on the high side. With every start of these candidates within 5 % in Alk0 (23.64-23.71), the multimodality rule does not apply, so a one-sided limit is valid.
- "Compatible, not a test" against the S2 requirement of 20.8 meq/L is correct.
- The model status row may read GO WITH DOCUMENTED NON-BLOCKING DEBT: 9 of 9 converged, and the selection is robust (§4).

## 8. Y7: M6 per volume scenario

**State model.**
- V(t) = V0 − E t; dC/dt = −r + C E/V.
- A pulse adds dose × V0/V(t_p).
- n, the N consumed per initial litre, grows at r/φ.
- The alkalinity state is φ (Alk0 − 2n/M_N).
- The biomass activity is multiplied by φ = V0/V(t).

Each term is the correct concentration balance for conserved amounts in a shrinking volume (`model.py:413-468`).

**Activity scaling (ruling 7).** The scaling assumes that biomass is conserved and concentrates with the liquid, and that the specific activity is unchanged by the concentration. This is the standard assumption, it is the consistent counterpart of concentrating every other conserved amount, and it is stated in the code, `model.md` and the registered decision. **Accepted.** An alternative, activity per vessel rather than per litre, would only rescale v0.

**Tests.**
- At E = 0 the scenario simulation equals the constant-volume one to 1e-12 for every candidate.
- N amount is conserved under evaporation to 1e-6.
- Every pulse step equals dose × φ(t_p).
- A synthetic record generated at 1 mL/d is recovered better by the scenario fit than by the constant-volume fit (max log parameter error 0.0065 against 0.0166, and lower RSS).
- The constant-volume path (`simulate`) is untouched, which is why every primary table is byte-identical.

**Labels.** Every row of the two new tables carries `volume_scenario` with SCENARIO and the role "SCENARIO refit: never feeds a primary row". `check_scenario_rows` raises before the write. No primary table names a refit.

**Selection.** `biomass_decay` is selected under every scenario (selection skill 0.958, 0.957, 0.951), and every candidate converges. Confirmed from `model_volume_scenario_summary.csv`.

**The scenario C(t) "above the band" (P2-A).** The prompt asked the reviewer to rule on whether the primary band now understates its uncertainty. It does not, in the sense that matters. Here is the scenario end level divided by the constant-volume end level (293.85 mg N/L), against the concentration factor φ = V0/V(t) at the end window (t ≈ 47.87 d):

| Scenario | End C(t) (mg N/L) | Ratio to constant volume | φ at the end |
| --- | --- | --- | --- |
| 0.25 mL/d | 334.4 | 1.138 | 1.136 |
| 0.5 mL/d | 388.3 | 1.321 | 1.315 |
| 1 mL/d | 569.5 | 1.938 | 1.918 |

At the start the ratios are 1.013-1.031, against φ of 1.003-1.011.

- So the refitted C(t) is, to within about 1-2 %, the constant-volume C(t) times V0/V(t). The **amount** of N per initial litre that the electrode record implies is essentially the same under every scenario, and only its concentration changes with the liquid volume.
- Amount-based quantities (η, `removal_to_dose_ratio`, `n_removed_mgN`, the O2 requirement) are therefore insensitive to the volume scenario at this order. Concentration-based quantities scale by φ, and W8 already carries them as SCENARIO columns (`consumption_mgN_L_d_SCENARIO…` in `nh4_truth_windows.csv`, and the scenario balances).
- The constant-volume band is a concentration band conditional on constant volume, and §1 point 2 already states that it is open upward. "Above the band" at 0.5 and 1 mL/d is the evaporative concentration factor, not new evidence about the amount of N.
- This is not a P1 integrity risk. But §1 point 9 reports "388 at the end (above the band)" and "570 … (above the band)" without this explanation, and a reader is likely to read it as a contradiction of the primary C(t). **P2-A:** §1 point 9 and `model.md` should state that the scenario C(t) ≈ φ(t) × the constant-volume C(t), so amounts and η are unchanged and concentrations scale by V0/V(t). The summary table could carry the ratio and φ as columns.

## 9. Pulse traces and test adequacy

**Traces.**
- `pulses.csv`, `nh4_truth_windows.csv`, `nh4_truth_w2_consistency.csv` and `nh4_truth_summary.csv` are byte-identical to `ad51a6a` in the reviewer's regenerated outputs.
- The pulse-1 and pulse-2 traces of the pass-3 review (§7 there) therefore hold unchanged: one series `C_inv:c29e26a0fb89`, MODEL-BASED rows on `M6_state:biomass_decay`, and `w2_bound_check` consistent with the checkpoints (pre_pulse1 and pre_pulse2 consistent, pre_pulse3 contradicted).

**New tests (`tests/test_pass4_fixes.py`, 26) against the prompt's table.**

| Type | Tests |
| --- | --- |
| Positive | Skill range spans both minima on a synthetic two-minimum record; every allow-listed metric passes; the scenario refit recovers the generating parameters better. |
| Negative | 7 derived metrics declared NH4-independent raise; undeclared primary wide column raises; unconverged Alk0 never enters the range; scenario row without SCENARIO raises. |
| Boundary | Exactly 0.5 coverage computes; two starts within the tolerance are one minimum. |
| Regression | `test_real_selection_fit_and_model_based_numbers_equal_pass3` compares the selected model's fit, the selection, `model_identifiability.csv`, `nh4_truth_summary/vectors.csv` and every MODEL-BASED value of `pulses.csv` with the files at `ad51a6a`, read with `git show`. It compares against stored pass-3 values, not values recomputed by the code under test ✓. One weakness: it skips instead of failing when git history is unavailable. That is fine in a checkout, but it would pass silently in a source tarball (P2-E). |
| Integration | Skill range and limits in the report; extended checker on the real states tables; scenario tables labelled, with one selection per scenario. |
| End-to-end | The existing determinism test, plus the reviewer's byte-identical run. |

**Changed existing tests, compared with `git show 083dde2:<file>`.**
- `test_states.test_conditional_idle_fraction_needs_the_coverage_but_the_bounds_do_not`: the removed `at_threshold` assertion repeated the first case's 4-of-10 input exactly. Its replacement asserts a true 5-of-10 case (fraction 0.5, conditional value 2/5, status computed). **Strengthened.**
- `test_pass3_fixes.test_every_checked_table_passes_the_hardened_checker`: it now passes the declared column lists that Y4 requires, and the call without them would raise. The assertion is unchanged, against a stricter checker. **Not weakened.**

**The development change to a new test.**
- `test_scenario_refit_recovers_the_generating_parameters_better_than_constant_volume` first asserted a scenario error below a third of the constant-volume error. That margin was the implementer's own, and it failed (0.0065 against 0.0166, a ratio of 0.39).
- It now asserts "better", an absolute error below 2 %, and a lower RSS. The prompt requires "better".
- The test is new in this pass, so no existing requirement was relaxed. The disclosure is complete, and the remaining assertions are still meaningful (both an absolute and a relative criterion). **Accepted.** It would have been cleaner to fix the margin from the prompt before the first run.

**No weakened test found.**

## 10. Rulings on the implementer's judgment calls

| # | Judgment call | Ruling | Reason |
| --- | --- | --- | --- |
| 1 | Y1 report rule ("holds", "undetermined", disjoint for better/worse) | **Accept** | Conservative for "holds", symmetric for ranking (§4). |
| 2 | Y2: starts within 5 % in Alk0 are one minimum | **Accept** | Reuses the registered `MULTIMODAL_LOG_TOL`. It only affects listing, never a value. |
| 3 | Y3: φ = 0.98 added; quote the grid point nearest the record's ρ̂ | **Accept** | Fixed from a known input (ρ̂ = 0.926), and the result is less favourable (§5). |
| 4 | Y4: indirect dependence declared C(t)-derived | **Accept** | The safe direction. |
| 5 | Y6: (α, Ks, S_opt), not the prompt's example; M3 fits unchanged | **Accept, with P2-C** | Justified by the reproduced diagnosis. Conditioning is not improved (§7). |
| 6 | Y6: the same log box in the new coordinates | **Accept, with P2-B** | A uniform rule, and no selection depends on it, but the bound-constrained calibrations must be flagged (§7). |
| 7 | Y7: biomass activity × V0/V(t) | **Accept** | The standard conserved-biomass assumption, stated and registered (§8). |
| 8 | Y7: scenarios in threads | **Accept** | Byte-identical outputs at a different worker count (§1). |
| 9 | Baseline pytest not rerun | **Accept** | The baseline outputs were regenerated byte-identically, and the pass-3 review ran the full suite on the same code (`418 passed`). |

## Remaining issues by severity

| Severity | Issue |
| --- | --- |
| P2-A | Report §1 point 9 and `model.md` show the scenario C(t) "above the band" at 0.5 and 1 mL/d without saying that it is the constant-volume C(t) × V0/V(t) to within 1-2 %: N amounts and η are unchanged, and only concentrations scale. Add the explanation, and the ratio and φ to `model_volume_scenario_summary.csv`. |
| P2-B | The chosen calibrations of `haldane_FA_SCENARIO_T20C` (Ks 991.7; `a` on its bound) and `combined_SCENARIO_T25C` (Ks, C0 and `a` on their bounds) are bound-constrained, with natural parameters far outside physical ranges (vmax 9.4e4 mg N/L/d, Ki 1.8e-4 mg NH3/L). Neither report says so. Flag bound-constrained calibrations in `model_fit.csv`, `model.md` and §4. Selection unaffected. |
| P2-C | The reparametrization removed the near-perfect correlations (0.99998 → 0.89) but not the ill-conditioning (Jacobian condition 4.2e6 at the new alkalinity × FA optimum). The alkalinity × FA natural parameters (vmax 2,601, Ki 0.0047) are not estimates. `progress.md` should say so; the report correctly quotes only skills and Alk0. |
| P2-D | §1 point 7 is a single paragraph of about 600 words with three skill ranges and the Alk0 profiles inline. Move the ranges to a table in §4 and keep §1 to the conclusions. |
| P2-E | `test_real_selection_fit_and_model_based_numbers_equal_pass3` skips instead of failing when git history is unavailable. |

The W8 debt carried since pass 2 is closed. No earlier P2 item remains open.

## What is sound and should be kept

- A Step 1 that demonstrably changed no fitted value, selection or MODEL-BASED number, with a regression test anchored to the stored pass-3 outputs.
- Skill ranges over statistically indistinguishable calibrations, with a symmetric ranking rule. The report now says "undetermined" where the data cannot decide.
- A provenance checker with allow-lists both for metrics and for declared columns, and negative tests for each gap the previous review found.
- An honest synthetic coverage table that reported a worse number (0.65) than the review had estimated, and used it.
- A reparametrization chosen from a reproduced Jacobian diagnosis, proved to be the same rate law, and applied to every candidate that shares it.
- A volume-scenario model whose conservation laws are tested, whose outputs never reach a primary row, and whose results leave the selection unchanged.
