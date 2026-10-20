# Independent re-review of remediation pass 3 (2026-10-20)

This review covers commit `adda251` on branch `remediation/rereview-2026-10-19` (parent `954f87a`). The claims in `docs/progress.md`, section "Remediation pass 3", were treated as claims to verify. The procedure is `docs/rereview_prompt_2026-10-20_pass3.md`, and the reference documents are `docs/remediation_prompt_2026-10-19_pass3.md` (X1 to X8), `docs/rereview_2026-10-19_pass2.md` (N1, N2, the P2 table) and `docs/remediation_prompt_2026-10-19.md`.

Every command was run on clean detached worktrees of `adda251`. One worktree was used for pytest, one for run-all and one for the reviewer's scratch scripts, so no two runs could interfere. The untracked `src/bnr_non_equilibrium_2026/` stub and `docs/rereview_prompt_2026-10-20_pass3.md` are not part of the commit, and the worktrees did not contain them.

## Decision

**GO WITH DOCUMENTED NON-BLOCKING DEBT.** Pass 3 resolves both P1 defects of the previous decision:

- **N1:** every candidate is calibrated under one convergence policy, and only converged candidates can be selected. The selected `biomass_decay` beats every alkalinity calibration, converged or not, so the selection and every MODEL-BASED quantity no longer rest on an unequal comparison.
- **N2:** the slope intervals are autocorrelation-scaled. The idle fraction is reported only as bounds with its coverage (whole record 0.05 to 0.93 at 12 % classified), with "NO-GO for a single value".

All three commands pass on a clean checkout, and the outputs regenerate byte-identically. No acceptance criterion fails, no failure condition is met, no test is weakened, and every judgment call is accepted. Six P2 items remain (last table). The most important is P2-1: the report presents the competing alkalinity mechanism only at its worse-forecasting calibration, although the data cannot tell it apart from one that forecasts better than decay. It should be fixed at the start of the next pass.

## 1. Commands (run by the reviewer, 2026-10-20)

| Command | Result |
| --- | --- |
| `uv run ruff check` | `All checks passed!`, exit 0 |
| `uv run pytest -q` (`BNR_NEQ_WORKERS=5`) | `418 passed in 19633.81s (5:27:13)`, exit 0. This matches the implementer's count (381 from `954f87a` plus 37 new). It includes the two fresh run-all fixtures and the byte-identical CSV/SVG determinism test. Wall clock is long because pytest shared the CPU with run-all and the reviewer's scratch fits. |
| `rm -rf outputs && uv run python -m bnr_neq run-all` (`BNR_NEQ_WORKERS=5`) | Exit 0, wall 2:27:36 (shared CPU). Stage messages: verdicts NH4/EC/ORP FAIL, pH PASS; segments references 11/13; thermo 9/11; dynamics 28/28; model selected `biomass_decay`, 13/15 identifiable; states 12 % classified, idle between 5 and 93 %; pulses 274 computed, 78 not computed. These match `progress.md`. |
| `diff -rq` against a copy of the committed `outputs/` taken before the run | No differences in `tables/` (63 files), `reports/` or `figures/` (SVG and PNG). |
| Determinism across worker counts | Run-all used 5 workers, while the implementer used the default (all cores), and the outputs are byte-identical. A separate `BNR_NEQ_WORKERS=1` model-stage rerun was not made. Its purpose, showing that the results do not depend on the worker count, is covered by this cross-count byte identity. |

Environment: CPython 3.14.8, uv 0.13.1 (warns that `uv-build<0.13.0` does not contain it; harmless), numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2. Scratch scripts (`x1_check.py`, `x3_probe.py`) ran in a third clean worktree.

## 2. Previous remaining issues

### 2.1 "Remaining issues by severity" of `docs/rereview_2026-10-19_pass2.md`

| Previous issue | Status | Evidence |
| --- | --- | --- |
| P1 N1: W5 calibrations unconverged; unequal comparison; Alk0 unreported | **Resolved** (P2 reporting remnants, §4) | One policy for every candidate: `model.py:46-53` (`CAL_START_FACTORS`, `CAL_MAX_NFEV`, `CONVERGENCE_POLICY`), `calibration_starts` (`model.py:465`), `calibrate_all`, `best_start` (`model.py:484`). Selection is gated on convergence: `select_converged` (`model.py:588`) raises `SelectionError` when no candidate converged. `model_fit.csv` carries `cal_converged`, `cal_status`, `cal_nfev`, `cal_start`, `selection_status` and `alk0_*`, and `model_calibration_starts.csv` lists all 33 starts. The two unconverged candidates read `not computed: calibration did not converge (status 0: …); skills not comparable` and keep their skills. The selected model is `biomass_decay` (converged), and that choice is robust to every unresolved point (§4). |
| P1 N2: idle fraction a biased lower bound presented as a value | **Resolved** (P2 reporting remnants, §5) | n_eff-scaled Theil-Sen intervals (`states.py:159-190`), bounds and coverage (`states.py:298-331`), `pulses.csv` metrics `idle_fraction_LOWER_BOUND/UPPER_BOUND/of_classified_time_CONDITIONAL` and `classified_fraction`, each basis citing the coverage. Report §1 point 8, the §2 table, §3, Fig. 8c and a separate status row read "BOUNDS ONLY; NO-GO for a single value". Hand recomputation in §5. |
| P2 W1 checker holes | **Resolved** (two new latent holes, P2, §6) | Every row declares a basis, and an empty or unknown model raises. Raw-reading rows are checked against an allow-list. The MODEL-BASED label is matched as a token. The checker runs before `pulses.csv`, `nh4_truth_windows.csv`, `states_*`, `segments.csv` and `thermo_oxygen.csv` are written. |
| P2 gain-corrected columns without SENSITIVITY | **Resolved** | `segments.csv`: `consumption_SENSITIVITY_linear_gain_*` (`segments.py:281-283`). `thermo_oxygen.csv`: `consumption_/o2_demand_/kla_required_/kla_measured_over_required_SENSITIVITY_linear_gain*`. `pulses.csv`: `consumption_SENSITIVITY_linear_gain_{first_window,weighted}_mgN_L_d`. The checker enforces the names in wide tables (probe in §6). |
| P2 η band above 1; empty baseline `w2_bound_check` | **Resolved** | Pulse 2 `eta_conversion_fraction` band 0.618-1.000, with the basis naming the unclipped 0.618-1.069 that stays on `removal_to_dose_ratio` (`pulses.py:140`). All ten baseline windows of `nh4_truth_windows.csv` read `pre_pulse1: consistent: the removal needs Alk0 ≥ 9.5 meq/L if nitrified` (`pulses.py:769-790`, k = 0 in `qualify` and `check_w2_qualified`). |
| P2 R8 scope (window 14) | **Resolved** by a symmetric extension (judgment call 7) | `segments.py:36, 192-196`. Window 14 is demoted, and pulse 1's rate with it is kept as `consumption_weighted_mgN_L_d_R8_CASE_incl_fragments` = 7.153 (5.51-8.52). Details in §7. |
| P2 R11 (rounding, wording, `relaxation_time_d`) | **Resolved** | §4 tables rounded (`65.5`, `1.71e-13`, `20`). §1 point 2: "open upward over the whole series". §1 point 5: the accumulation rests "on the QUALIFIED C(t)". `dose_removal_time_d` carries "not a fitted relaxation constant". |
| P2 AE naming, DO = 0 against DO = 2 mg/L, source caveat | **Resolved** | Suffix `AE` in every metric. The supply-energy basis and report §2 state the DO inconsistency, and the caveat is quoted in report §2 and the config comment (`config/experiment.toml:37-41`; values unchanged). The quote is taken from the pass-2 review, where it was verified at the source, and is marked as such. |
| P2 W8 M6 not refitted per volume scenario | **Unchanged, stated debt** (accepted in pass 2) | Not in scope. |
| P2 `test_truth` η test branches; no convergence test; no gain-name test | **Resolved** | §9. |
| P2 MODEL-BASED η about half the dose-anchored η, undiscussed | **Resolved** | New `eta_gap_inversion_misfit_MODEL_BASED` and report §1 point 5. The identity was recomputed by hand (§7). |

### 2.2 Work items X1 to X8, sub-item by sub-item

| Item | Sub-item | Status | Evidence |
| --- | --- | --- | --- |
| X1 | (a) diagnosis | Verified | `progress.md` X1 (a): status 0 at scipy's default 100 × n_free, Alk0 22.5-22.7, no parameter at a bound. Reproduced by the reviewer for `alkalinity_limited` (§4). |
| X1 | (b) one policy | Verified | `CONVERGENCE_POLICY` is applied by `calibrate_all` to every spec. The same factors and the same `CAL_MAX_NFEV` are used, the observation start is shared, and only the null model (no kinetic parameter) has a single start. `test_policy_starts_are_the_same_rule_for_every_candidate` and `test_real_calibration_starts_follow_the_policy`. |
| X1 | (c) gate and labels | Verified | `select_converged`. `model_fit.csv`: the unconverged rows carry the not-computed status and keep their skills; `selected` is True only for `biomass_decay` (`cal_converged = True`). |
| X1 | (d) Alk0 reported, against S2 | Verified, with P2 wording remnants | `alk0_fitted_meq_L`, `alk0_ci95_*`, `alk0_profile_status` in `model_fit.csv`; report §4 table `alk0_comparison` (S2 requirement 20.8 at pre_pulse3, 20.9 at the end); §1 point 7. Remnants in §4. |
| X1 | (e) MODEL-BASED quantities regenerated and explained | Verified | Selected model unchanged. The changes come from biomass_decay's chosen start moving from x1 to x0.5 (§4). |
| X2 | (a) n_eff-scaled interval, C(t), pH, ORP; docstring and `states_rules.csv` | Verified | `states.py:159-190`; module docstring; `states_rules.csv` `thresholds`. |
| X2 | (b) states recomputed, registered thresholds unchanged | Verified | 15 mV, 1.01 d, block rule and ≥ 2 blocks unchanged (`git diff`). |
| X2 | (c) bounded pair, CONDITIONAL value, names | Verified | `states.py:316-330`; `pulses.py:283-306`. |
| X2 | (d) coverage threshold | Verified | `MIN_CLASSIFIED_FRACTION = 0.5` (`states.py:59`), with its reason in the code and in `progress.md`. |
| X2 | (e) report, figure, status | Verified | §1 point 8, §2 table, §3, Fig. 8c (inspected), status row "BOUNDS ONLY; NO-GO for a single value". The W3 row no longer quotes an idle fraction. |
| X3 | (a)-(d) | Verified | §6. |
| X4 | renames and readers | Verified | §2.1 row; `test_reading_based_columns_say_sensitivity` extended. |
| X5 | (a), (b) | Verified | §2.1 row; `test_eta_band_edge_above_one_is_clipped`, `test_real_eta_bands_and_baseline_labels`. |
| X6 | rule and labelled case | Verified | §7. |
| X7 | (a)-(f) | Verified | §2.1 rows; (e) recomputed in §7. |
| X8 | tests | Verified, with P2 remarks | §9. |

## 3. Acceptance criteria and failure conditions of the pass-3 prompt

| Criterion / failure condition | Result |
| --- | --- |
| N1 and N2 resolved with evidence in code and outputs | Met (§4, §5). |
| X3 to X8 implemented and tested | Met (§6, §7, §9). |
| All three commands pass on the final code | Met (§1). |
| Outputs reproducible | Met: byte-identical regeneration, also at a different worker count (§1). |
| Report states the converged comparison, the fitted Alk0 against S2, and the idle fraction as bounds with coverage | Met. The converged comparison is stated incompletely (P2-1, §4). |
| A non-converged candidate can be selected, or its skills appear without a convergence label | Not met. `select_converged` excludes it, `test_unconverged_candidate_cannot_be_selected` covers it, and every skill of an unconverged row is labelled `skills not comparable`. |
| Convergence policy differs between candidates | Not met (one rule; §4). |
| Interval method or coverage threshold chosen after looking at its effect | Not met, as far as can be judged (§5). |
| An idle fraction without its coverage | Not met. Every idle row of `pulses.csv` cites `coverage: classified fraction …`, and so do the §1 and §2 tables, Fig. 8c and the status row. |
| A reading-based column or metric without SENSITIVITY in its name | Not met. Every `linear_gain_SENSITIVITY` and `identity_…_SENSITIVITY` row and every reading-based wide column says SENSITIVITY. `consumption_uncorrected_*` is the measured raw slope (§8). |
| An output mixes observation models or series and is not caught | Not met. No committed output does, and the checker passes on every checked table. Two latent holes remain (§6, P2). |
| An existing test weakened | Not met (§9). |

## 4. N1 / X1: the converged comparison

**Diagnosis reproduced.** With the committed code on the regenerated grid (scratch script `x1_check.py`; calibration window 1.0-24.59 d, 2,359 bins, 13 free parameters), the reviewer reproduced three runs:
- `alkalinity_limited` from the default start with scipy's default budget (the pass-2 call) stops with **status 0 after exactly 1,300 evaluations** (100 × 13), at Alk0 = 22.73 meq/L and RSS 29.3253. That is the diagnosis of `progress.md` (Alk0 22.7; RSS 29.325 after 1,300). No parameter is at a bound.
- Under the policy, start x0.5 converges (status 2) after **990** evaluations at Alk0 = 22.637 and calibration RMSE 0.1114908.
- Start x4 converges after **59** evaluations at Alk0 = 79.965, C0 = 33.6 and RMSE 0.1092763.

All three runs match `model_calibration_starts.csv` to every printed digit.

**The policy is identical for every candidate.**
- `calibration_starts` shifts the whole kinetic block by log(f), f ∈ {0.5, 1, 2, 4}. The kinetic parameters are the first `len(spec.kinetic)` entries of the layout and are all on the log scale (`model.py:339-347`).
- The observation start, the tolerances and `CAL_MAX_NFEV = 8000` are shared, and the whole-record refit gets the same budget (`model.py:932`).
- For the alkalinity candidates the starts are Alk0 = 5, 10, 20 and 40 meq/L, which is the grid the prompt proposed.
- Only the null model, which has no kinetic parameter, has one start, and it is only the reference.
- No candidate gets a budget the others do not.

**Selection and labels.**
- In `model_calibration_starts.csv` the chosen start of every candidate is its lowest-RSS converged start:
  - biomass_decay: x0.5 (0.1137772), 5e-9 below x1 and x2;
  - haldane T20C and T25C: x1;
  - combined T20C: x4; combined T25C: x0.5;
  - alkalinity_limited: x4 (0.109276).
- The two alkalinity × FA candidates have no converged start (all four starts at 8,000 evaluations, status 0). They keep the lowest-RSS start, labelled `cal_converged = False`, and are excluded from selection.
- `test_real_calibration_starts_follow_the_policy` asserts this on the real run.

**Judgment of the result.**
- **The selection of `biomass_decay` is robust.** Its selection skill of 0.958 exceeds every value any alkalinity calibration reaches:
  - `alkalinity_limited` at Alk0 = 80: −0.908;
  - `alkalinity_limited` at the converged 22.6 minimum: 0.934;
  - the unconverged alkalinity × FA starts: 0.797-0.801.

  The selected model, and so every MODEL-BASED quantity, therefore does not depend on how the alkalinity calibrations are resolved. This is what N1 required. **N1 is resolved.**
- **Lowest RSS, not lowest skill, is the right calibration rule (judgment call 3, accepted).** The RSS minimizer is the calibration-window estimate. Choosing among starts by selection or forecast skill would fit the calibration to the selection window and make the selection comparison circular.
- **The conclusion is only partly correct as written.** "Decay selected among converged candidates" and "alkalinity × FA undecided" are correct. "The pass-2 claim of better out-of-sample alkalinity skill withdrawn" (`progress.md`) is not accurate:
  - The pass-2 skills (selection 0.934, test 0.817, forecast 0.925) now belong to a **converged** minimum (starts x0.5 and x1, Alk0 = 22.6).
  - That minimum lies at Δ = 3.0 < 3.84 from the chosen calibration, so the calibration window cannot reject it.
  - Within its own 95 % region, single alkalinity limitation therefore ranges from far worse than decay out of sample (test 0.227, forecast 0.443) to better than decay (test 0.817 against 0.748, forecast 0.925 against 0.863).
  - The correct statement is that the out-of-sample ranking of decay and alkalinity limitation is **undetermined**, not that the claim was an artifact of non-convergence.
  - Report §1 point 7 quotes only the Alk0 = 80 skills, and the 22.6 minimum appears only as an Alk0 value, without its skills. A reader is led to conclude that alkalinity limitation forecasts worse. This is **P2-1**. It does not change the selection or any MODEL-BASED number, and the report keeps "time is confounded …" and the QUALIFIED status. The fix is to state the skill range of `alkalinity_limited` over its starts inside Δ < 3.84 in §1 point 7 and §4.
- **Alk0 reporting.**
  - Alk0 = 80.0 is not at a bound (log bounds 1e-4 to 1e3; `model.py:187`). It is a plateau: the limitation term A/(K + A) is near 1 throughout the calibration window, and the local profile stays below 3.84 on the high side.
  - The `alk0_comparison` verdict "compatible, not a test" is correct: both converged minima, 22.6 and 80, exceed the S2 requirement of 20.8, but no interval exists.
  - Remnants, both P2-2:
    - The status-table evidence (`report.py:324`) quotes "fitted Alk0 of the alkalinity candidates 23.6-80.0 meq/L". The 23.6 is the unconverged value that the same report calls "not an estimate".
    - The profile status lists "Alk0 22.6 meq/L at Δ = 3.01" twice (one entry per converged start).
- **Post-hoc multimodality rule (judgment call 4).** The rule withholds the interval when another converged start lies inside Δ < 3.84 (`model.py:731-768`). It is a validity condition, not tuning:
  - A profile warm-started from one minimum cannot describe a 95 % region that contains two minima.
  - The limit it replaced (one-sided lower limit 80.0 = the estimate) was invalid on its face.
  - The rule only withholds an interval: it changes no fit, selection or other output. It is disclosed, and it is tested in both directions (`test_alk0_interval_is_withheld_…`).
  - **Accepted.** The hand recomputation of Δ: n_eff is taken from the residuals of the chosen (x4) calibration: ρ = 0.9384, n_eff = 2359 × (1 − 0.9384)/(1 + 0.9384) = **74.91**. Then Δ = 74.91 × ln(29.32282/28.16955) = 74.91 × 0.040125 = **3.006** < 3.84 ✓. This matches the 3.01 in `model_fit.csv` and the n_eff 74.9 and ρ 0.94 quoted in `progress.md`.
- **Truncated simulation in `fit`.** `simulate` integrates forward (RK4), and `observe` applies a causal `lfilter` with initial state from the first bin. Nothing after the window end enters a residual inside it. The reviewer compared the calibration-window residuals of `alkalinity_limited` at its start vector with full and truncated data: **bit-identical** (`np.array_equal`, 2,359 residuals). The truncated fits also reproduce the committed per-start RSS and evaluation counts exactly (above), so the change did not alter any fit.
- **Changes in C(t).**
  - `nh4_truth_windows.csv` median C rises by a factor of 1.000157 to 1.000498, inside `progress.md`'s whole-series range of 1.00015 to 1.00054.
  - The cause is that biomass_decay's chosen start moved from x1 to x0.5. Its selection skill moved from 0.957776 to 0.957804, the other converged candidates changed in the 6th to 7th digit, and the whole-record refit (C0 is not identifiable, so the refit sits in a flat valley) lands 1e-4 relative away. The new truth identifier is `C_inv:c29e26a0fb89`.
  - `nh4_truth_summary.csv`: start 177.115 → 177.158, end 293.799 → 293.849, band edges in the 5th digit. All downstream 4th-digit changes trace to this.

## 5. N2 / X2: the idle fraction

**The interval matches the prompt's equation.**
- `rolling_theil_sen` multiplies each half-width by `interval_scale(n, ρ)` = √(n/n_eff), with n_eff = n (1 − ρ)/(1 + ρ) and n = 101.
- ρ is the lag-1 autocorrelation of the residuals about the window's own Theil-Sen line (that is, detrended), clipped to [0, 0.999]. It is 0 when the residuals have no variance.
- The same function is used for C(t), pH and ORP. The method is stated in the docstring and in `states_rules.csv`.
- Median ρ on the regenerated bins: C(t) 0.926, pH 0.772, ORP 0.985 (from `states_bins.parquet`, matching `progress.md`).
- Of the 2,031 bins with a C(t) slope, 1,472 are flat, 332 falling and 227 rising, as claimed.

**Pre-registration.** The history has one pass-3 commit, so the order of decisions cannot be read from git. Three things make tuning implausible:
1. The interval method is the one the remediation prompt prescribed, and the code uses it without modification.
2. The synthetic validation that `progress.md` records "before any real-data output" kept the estimator even though the result was unfavourable (0.85 coverage). A tuned choice would have adopted the bias correction.
3. The coverage threshold is not decision-relevant. Every interval has a classified fraction ≤ 0.209, so any threshold between 0.21 and 1 gives the same outputs. Only a threshold below 0.21 would have changed anything, and 0.5 is not one.

**Residual undercoverage (judgment call 5).**
- On AR(1) data at φ = 0.9, the registered ρ̂ (biased low on 101 detrended points) gives coverage of about 0.85. The test asserts:
  - 0.92 to 0.98 for the formula with the true ρ;
  - ≥ 0.80 for the production estimator;
  - ≤ 0.50 for the unscaled interval.
- This is an honest implementation of a prompt whose two requirements (the ρ̂ estimator and "about 95 %") are incompatible at n = 101. The test proves the formula and bounds the shortfall. **Adequate; judgment call 5 accepted.**
- The shortfall, and its consequence that at ρ̂ ≈ 0.93 roughly 1 in 5 truly flat windows can still be called rising or falling, is stated only in `progress.md`. Neither `analysis_report.md` nor `states.md` mentions it. The same holds for the conditionality of the bounds on the classification:
  - the upper bound counts `active_nitrification` and `oxygen_limited` bins as not idle;
  - the lower bound counts only bins classified inactive.

  This is **P2-3**. The bounds stay usable: undercoverage can only move bins from flat to falling or rising, which keeps the lower bound conservative and lowers the upper bound by at most the active share (whole record 0.056).

**Hand recomputation, pulse 2 (10.57-24.59 d, `states_bins.parquet`).**
- 1,402 bins: 901 sensor artifact, 333 unclassified, 121 active, 38 inactive, 9 oxygen-limited.
- Lower bound 38/1402 = **0.02710**, upper bound (38 + 333 + 901)/1402 = 1272/1402 = **0.90728**, coverage 168/1402 = **0.1198**. These match `states_fractions.csv` and `pulses.csv` (0.027104, 0.907275, 0.119829).
- The CONDITIONAL value is `not computed: classified fraction 0.120 below 0.5 …`.

**Names, coverage and report.**
- The three idle metrics carry LOWER_BOUND, UPPER_BOUND and CONDITIONAL, and every basis carries `coverage: classified fraction …` (test `test_real_idle_metrics_carry_bounds_and_coverage`).
- §1 point 8 gives "P1 0.00-1.00 (coverage 0.00) … P4 0.10-0.96 (coverage 0.15)". The §2 column reads "idle fraction: bounds (coverage)".
- Fig. 8c is titled "bounds, not a value" and prints the coverage above each bar.
- The status row reads "BOUNDS ONLY; NO-GO for a single value". No row reads GO for the idle fraction.

**What the state layer supports.** At 12 % coverage, the idle fraction spans 0.05 to 0.93 over the whole record. That is an honest "not determined", not a measure. The layer supports two claims:
1. At the registered thresholds, NH4-free indicators do not separate consumption: pH rising vs not +1.01 mg N/L/d, 95 % block bootstrap −0.82 to +2.61; ORP not computed.
2. A circular, descriptive active-vs-inactive difference (+3.57).

The report makes no stronger claim. The W3 row "GO WITH DOCUMENTED NON-BLOCKING DEBT" refers to the classifier and the association test, not to the idle fraction, and that is acceptable.

## 6. X3: the checker (`model.check_single_observation_model`)

Probed by appending rows or columns to the committed tables (`pulses.csv`, `nh4_truth_windows.csv`, `segments.csv`, `thermo_oxygen.csv`, `states_fractions.csv`) and calling the checker as the writers do. Scratch script; reproducible.

| Probe | Result |
| --- | --- |
| Empty model on `fa_median_…` (pulses), empty model (truth windows) | Caught: "rows without a declared NH4 basis" |
| Unknown model `M7`; trailing-space `M6_log ` | Caught |
| `consumption_weighted_mgN_L_d`, `aeration_energy_supply_…`, `idle_fraction_LOWER_BOUND`, `n_removed_mgN` labelled `electrode_reading` | Caught (allow-list `RAW_MEASURED_METRICS`) |
| Unlabelled gain column in `segments.csv`; `o2_demand_mgO2_L_d` in `segments.csv`; `kla_required_per_d_SCENARIO_T25C` in `thermo_oxygen.csv` | Caught: "reading-based columns without SENSITIVITY" |
| Long-table SENSITIVITY metric without SENSITIVITY in its name | Caught |
| A series on an `NH4-independent` row | Caught |
| Second `C_inv:` series (pulses, `states_fractions`) | Caught |
| An idle row as `electrode_reading` in a wide `states_*` table | Caught (the wide table declares no NH4-independent columns) |
| "not MODEL-BASED" in a primary basis | Passes, correctly (the token rule fixes the false positive) |
| "NOT MODEL-BASED" (upper case) on a primary row | Raises: a false positive in the safe direction |
| **An NH4-derived metric (`consumption_weighted_mgN_L_d`) declared `NH4-independent`, without a series** | **Not caught** |
| **A reading-derived column added to `nh4_truth_windows.csv`, whose rows are `M6_log`** | **Not caught** (column names of wide tables are checked only on reading-based rows) |

- The two holes are latent. The `NH4-independent` marker is a self-declaration without an allow-list, the analogue of the empty-model hole it replaced. Wide tables with primary rows do not vet their columns.
- No committed output exploits either hole, so the failure condition "an output mixes … and is not caught" is not met. This is **P2-4**. The fix is an allow-list of NH4-independent metrics (the seven names below) and the `independent_cols` declaration for primary wide tables too.

**Order.** The checker runs before every write:
- `pulses.compute` (`pulses.py:894-895`) before `pulses.csv` and `nh4_truth_windows.csv`;
- `states.compute` (`states.py:469`) before `states_bins/timeline/fractions/comparison`;
- `segments.slope_table` (`segments.py:276`) before `segments.csv`;
- `thermo.oxygen_table` (`thermo.py:534`) before `thermo_oxygen.csv`.

**The 28 `NH4-independent` rows of `pulses.csv`.** Seven metrics × 4 pulses, each checked against its basis and inputs:
- `t_pulse_d`, `segment_end_d` (config);
- `median_pH` (grid pH);
- `o2_supply_capacity_mgO2_SCENARIO_kLa{1,10,100}_T25C` (kLa × C* × V × Δt; DO = 0, config);
- `free_energy_per_fully_converted_pulse_J_SCENARIO_T25C` (dose × ΔG°′, config and constants).

None rests on an NH4 measurement or C(t) ✓. The stirrer rows are `not computed` and now declare M6_log with the model-state series, consistent with their computed form (MODEL-BASED nitrogen).

## 7. X4 to X7, and the pulse traces

- **X4.**
  - No reading-based column or metric lacks SENSITIVITY: the checker enforces it, and `test_reading_based_columns_say_sensitivity` asserts it on `segments.csv` and `thermo_oxygen.csv`.
  - `consumption_uncorrected_*` stays on the allow-list as the measured raw slope (−slope × the NH4+ → N mass factor, basis "measured (electrode reading, uncorrected)"). A unit conversion of a measured slope is not an observation model, so it is **acceptable**. The word "consumption" is borrowed but qualified by "uncorrected" and the basis.
- **X5.** Every computed η row has its band in [0, 1] (`test_real_eta_bands_and_baseline_labels`). The baseline windows carry the pre_pulse1 label, and `check_w2_qualified` now also guards k = 0.
- **X6 (judgment call 7).**
  - Only window 14 changed its computed/not-computed status. Windows 2, 3, 5, 7 and 9 changed only the status text ("baseline window" → "window").
  - The pulse boundaries (`pulse k`) do not match `ELECTRODE_BOUNDARY_PREFIXES = ("jump (", "candidate (")`, so windows 10, 15, 19 and 21 stay computed.
  - Downstream changes, verified against `954f87a`: `dynamics_windows.csv` 680 → 679 rows; `thermo_co2.csv` 12 → 11; `segments_window_check.csv` 10 → 9; `thermo_oxygen.csv` 11 → 10; references `segments_references.csv` 11/13 and `thermo_references.csv` 9/11 unchanged.
  - Pulse 1 by hand:
    - The duration-weighted rate is window 15 alone: 4.2272 mg N/L/d.
    - O2 demand 4.57 × 4.2272 = 19.318 ✓.
    - Dose-removal time 65.464/4.2272 = 15.486 d, and /2.04 = 7.591 ✓.
    - R8 CASE (0.3102 × 18.29 + 1.1805 × 4.2272)/1.4907 = **7.153** ✓ (old primary value), interval 5.51-8.52.
  - The symmetric extension is justified. A delayed pulse-response jump is a step of the electrode reading, like an offset jump, and a window that starts 0.10 d after it lies inside the electrode's own 3τ (τ = 0.26 d). The governing pulse does not make the electrode settle faster.
  - The cost is a known downward bias of the earliest post-pulse rate. It is carried as the labelled R8 CASE. "R8 CASE" rather than SENSITIVITY is correct, because SENSITIVITY is reserved by the checker for an alternative observation model. **Accepted.**
- **X7.**
  - (a)-(d) as in §2.1.
  - (e) Recomputed from `nh4_truth.parquet` and `model_state.parquet` with 0.1 d pre-medians:
    - Pulse 1: C_inv 110.865 → 151.432, ratio 0.38032, gap 0.17389 (table 0.17389 ✓).
    - Pulse 2: ratio 0.80921, gap **0.23692** (table ✓).
    - Ratio − gap = 0.20643 and 0.57229 equals the median-based model removal per dose exactly. Against the integral-based MODEL-BASED η (0.20512, 0.56860) the difference is 0.0013 and 0.0037, which is the "within 0.004, up to the 0.1 d medians" that `progress.md` states ✓.
  - (f) AE suffix, DO statement and quoted caveat as in §2.1.

**Trace, pulses 1 and 2 (`pulses.csv`).**

| Pulse | Metric | Value | Model / series | `w2_bound_check` |
| --- | --- | --- | --- | --- |
| 1 | `consumption_weighted_mgN_L_d` | 4.227 | M6_log / C_inv:c29e26a0fb89 | pre_pulse1: consistent |
| 1 | `consumption_weighted_mgN_L_d_R8_CASE_incl_fragments` | 7.153 | M6_log / C_inv:c29e26a0fb89 | pre_pulse1: consistent |
| 1 | `eta_conversion_fraction` = `removal_to_dose_ratio` | 0.3803 (band 0.125-0.849) | M6_log / C_inv | pre_pulse1: consistent |
| 1 | `eta_conversion_fraction_MODEL_BASED`; `eta_gap_…_MODEL_BASED` | 0.2051; 0.1739 | M6_log / M6_state:biomass_decay | pre_pulse1: consistent |
| 1 | `aeration_energy_supply_…_AE1_kLa10_T25C` | 6.768 | M6_log / C_inv | pre_pulse1: consistent |
| 1 | `idle_fraction_LOWER/UPPER_BOUND`; `classified_fraction` | 0; 1; 0 | M6_log / C_inv | pre_pulse1: consistent |
| 1 | `o2_supply_capacity_…_kLa10_T25C` | 16.85 | NH4-independent / – | – |
| 2 | `consumption_weighted_mgN_L_d` | 2.599 | M6_log / C_inv | QUALIFIED (pre_pulse3) |
| 2 | `eta_conversion_fraction` | 0.8092 (band 0.618-**1.000**) | M6_log / C_inv | QUALIFIED |
| 2 | `removal_to_dose_ratio` | 0.8092 (band 0.618-1.069) | M6_log / C_inv | QUALIFIED |
| 2 | `n_removed_mgN` → `o2_required_mgO2` | 5.2975 → 24.209 (= 4.57 × 5.2975 ✓) | M6_log / C_inv | QUALIFIED |
| 2 | `eta_…_MODEL_BASED`; `eta_gap_…` | 0.5686; 0.2369 | M6_log / M6_state | QUALIFIED |
| 2 | `idle_fraction_LOWER/UPPER_BOUND`; `classified_fraction` | 0.0271; 0.9073; 0.1198 | M6_log / C_inv | QUALIFIED |
| 2 | `dose_removal_time_d` | 25.19 (= 65.464/2.5988 ✓) | M6_log / C_inv | QUALIFIED |
| 2 | `consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d` | 3.407 | linear_gain_SENSITIVITY / electrode_reading | – |

Every C(t)-derived row rests on the one series `C_inv:c29e26a0fb89` or is MODEL-BASED on `M6_state:biomass_decay`. Every `w2_bound_check` agrees with `nh4_truth_w2_consistency.csv`: pre_pulse1 and pre_pulse2 are consistent (9.53 and 13.1 meq/L), and pre_pulse3 is contradicted (20.8 meq/L), which bounds pulse 2.

## 8. Scientific correctness of the remaining interpretation

- The headline chain of conclusions is unchanged in direction:
  - pH the only validated sensor;
  - W2 bound;
  - S2 contradiction and QUALIFIED;
  - η/ratio split;
  - consistency checks not findings;
  - the closed frame not established.
- Pass 3 removes one overstatement (the idle value) and adds one understatement (alkalinity limitation, P2-1).
- The alkalinity × FA candidates stay unconverged after 4 × 8,000 evaluations. Their four starts agree to 1e-5 in calibration RMSE (Alk0 23.56-23.61), which points to a flat valley rather than a missed optimum. Excluding them and calling them undecided is the conservative reading. A reparametrization (for example vmax/Ks jointly) would be the way to settle them. It is listed as debt (P2-5), not a blocker, because their skills could not make them selected either way (§4).

## 9. Test adequacy

- **New tests against the prompt's table.**
  - Positive:
    - all nine candidates converge under the policy on a noise-free synthetic record, and Alk0 is recovered within 1 % for `alkalinity_limited` and `alkalinity_FA` on their own synthetic records;
    - AR(1) coverage (with the deviation discussed in §5);
    - idle bounds equal hand counts.
  - Negative:
    - an unconverged candidate is not selectable and is labelled;
    - no converged candidate raises;
    - the checker raises on an empty model, on 11 derived metrics as `electrode_reading`, on an unknown model, on a series on an independent row, on an unlabelled gain column, and on a wide reading table without declared columns;
    - "not MODEL-BASED" passes.
  - Boundary:
    - coverage gate below and at 0.5;
    - ρ cap;
    - η band clipping.
  - Regression: pinned η branch; gain names.
  - Integration (real run-all): converged selection, Alk0 reported, every checked table passes, η bands, baseline labels, R8 CASE, renamed metrics, policy starts.
  - End-to-end: `test_run_all_is_deterministic`.
- **Remarks (P2-6).**
  - In `test_conditional_idle_fraction_needs_the_coverage_but_the_bounds_do_not`, the case named `at_threshold` (`tests/test_states.py:219`) is a second below-threshold case (4 of 10, as its own comment says). The "at 0.5" boundary is covered only by the `half` case two lines later. The test is correct but misnamed.
  - No test covers the two checker holes of §6.
- **Changed existing tests, compared with `git show 954f87a:<file>`.**
  - `test_truth.test_dose_anchored_eta_and_affinity`: the `if/elif/else` on its own result is replaced by asserted branches (pulse 2 ratio 1.96 → "removal exceeds the dose"; pulse 1 −0.77 → "negative"). **Strengthened.**
  - `test_truth.test_check_accepts_empty_cells_read_from_csv` and `test_labelled_sensitivity_and_measured_rows_pass`: the fact row declares `NH4-independent` (X3 (a) makes an empty model raise), and the gain metric is renamed (X4). The NaN-label behaviour under test is unchanged. A requirement changed, nothing was weakened.
  - `test_rereview_fixes.test_electrode_fragment_rule` and `test_real_fragments_are_demoted`: the two "not a fragment" assertions for a pulse-response boundary and a pulse-governed window are inverted by the X6 decision. New assertions cover pulse boundaries (not electrode steps), window 14 (fragment) and window 15 (not). The baseline fragments 2, 3, 5 and 7 are still asserted. A requirement changed, nothing was weakened.
  - `test_states.test_all_artifact_interval_yields_no_state`: the conditional value is still NaN and not computed, and the bounds (0, 1) and coverage 0 are now asserted. **Strengthened.**
  - `test_states.test_states_on_the_real_record`: `idle_fraction_LOWER_BOUND` = `idle_fraction_lower_bound`, the renamed form of the same equality.
  - Mechanical renames in `test_pulses`, `test_segments`, `test_thermo_m4`, `test_energy` and the `test_report` fixture: every removed assertion line reappears with the new name and the same tolerance. `test_energy`'s family filter `"SAE|kLa"` → `"AE|kLa"` selects a superset, so the "every scenario value carries SCENARIO" assertion is at least as strict. The `test_pulses` and `test_truth` fixtures gain the model-state column the X7 (e) row reads, and no assertion changed.
  - The one test changed after the first full run (`test_real_idle_metrics_carry_bounds_and_coverage`, new in this pass) now accepts either correct not-computed reason. That fixes the test, not the code.
  - **No weakened test found.**

## 10. Rulings on the implementer's judgment calls

| # | Judgment call | Ruling | Reason |
| --- | --- | --- | --- |
| 1 | X1 start rule: the kinetic start × 0.5/1/2/4 for every candidate | **Accept** | One rule, applied identically, contains the prompt's Alk0 grid (5, 10, 20, 40) and moves no observation parameter. The warm start from biomass_decay was rejected on a recorded diagnosis (it lands at the far end of the C0 valley, RSS 29.360 unconverged), before any refit under the policy. Scaling every kinetic parameter, rather than Alk0 alone, also tests the other candidates' sensitivity to their start, which revealed biomass_decay's x4 local minimum. |
| 2 | X1 budget: 8,000 evaluations per start | **Accept** | More than twice the 3,620 evaluations of the converging diagnostic run, and fixed before the refit. The two candidates that fail stop on a plateau where four starts agree to 1e-5 in RMSE, so a larger budget is unlikely to change their skills, and none of their possible values would change the selection. |
| 3 | X1: calibration = lowest-RSS converged start, even when another converged start forecasts better | **Accept** | The RSS minimizer is the estimate. Choosing by forecast skill would use the selection window to calibrate. The consequence (two indistinguishable minima with opposite out-of-sample behaviour) must be reported, not resolved by a skill-based choice: P2-1. |
| 4 | X1 post-hoc rule: the Alk0 interval withheld when the profile is multimodal below Δ = 3.84 | **Accept** | A validity condition. A single-minimum profile interval is undefined for a 95 % region with two minima, and the limit it replaced equalled the estimate. It only withholds, is disclosed, and is tested both ways. Same class as the pass-2 "≥ 2 blocks" rule. |
| 5 | X2: ρ̂ without small-sample bias correction | **Accept, with P2-3** | The prompt specified the estimator. Correcting it after seeing 0.85 coverage would have been an unregistered change. The shortfall is quantified and tested, but must also be stated in the report and `states.md`. |
| 6 | X2: upper bound treats active and oxygen-limited bins as not idle | **Accept, with P2-3** | A bin with a significant C(t) fall is the one positive sign of activity the layer has. Counting it as possibly idle would push the upper bound to 0.93-1.00 in every interval and make it uninformative. The conditionality must be stated in the report. |
| 7 | X6: R8 extended symmetrically, labelled "R8 CASE" | **Accept** | §7. |

## Remaining issues by severity

| Severity | Issue |
| --- | --- |
| P2-1 | Report §1 point 7 and §4 quote `alkalinity_limited`'s skills only at the chosen Alk0 = 80 calibration (−0.91 / 0.23 / 0.44). The converged Alk0 = 22.6 minimum, inside Δ < 3.84, has 0.934 / 0.817 / 0.925 (better than decay on test and forecast) and is shown only as an Alk0 value. The out-of-sample ranking of decay against alkalinity limitation is undetermined. `progress.md` calls the pass-2 claim "withdrawn" because it "rested on an unconverged fit", but that minimum has now converged. The selection and every MODEL-BASED number are unaffected. |
| P2-2 | The status-table evidence quotes "fitted Alk0 … 23.6-80.0 meq/L", including the unconverged 23.6 that the same report calls "not an estimate" (`report.py:324`). The Alk0 profile status lists the 22.6 minimum twice. |
| P2-3 | The residual undercoverage of the scaled slope interval (about 0.85 at φ = 0.9, roughly 0.8 at ρ̂ ≈ 0.93), and the conditionality of both idle bounds on the state classification, are stated only in `progress.md`, not in `analysis_report.md` or `states.md`. |
| P2-4 | Checker: an NH4-derived metric declared `NH4-independent` (no series) passes, and the columns of wide tables with primary rows are not vetted. No committed output is affected. |
| P2-5 | The alkalinity × FA candidates did not converge under the policy (flat valley). Their mechanism is undecided, which is stated. A reparametrization would be needed to settle it. |
| P2-6 | `test_states.py:219`: the `at_threshold` case is a second below-threshold case. No test covers the P2-4 holes. |
| P2 (carried) | W8: M6 not refitted per volume scenario (stated debt, accepted in pass 2). |

## What is sound and should be kept

- A convergence policy that is one rule for every candidate, a selection that is gated on convergence, and a per-start table that makes multimodality visible instead of hiding it.
- A selected model that is robust to every unresolved calibration, which an independent reader can check from `model_calibration_starts.csv`.
- Autocorrelation-scaled slope intervals with a synthetic test that reports the shortfall of the estimator instead of tuning it away.
- An idle fraction that is reported as bounds with its coverage, with a NO-GO for a single value.
- A provenance checker that classifies positively, uses an allow-list for raw quantities, and runs before every table that carries C(t).
- A symmetric R8 rule, with the excluded value kept as a labelled case.
