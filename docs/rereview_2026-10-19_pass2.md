# Independent re-review of remediation pass 2 (2026-10-19)

This review covers commit `954f87a` on branch `remediation/rereview-2026-10-19`. That commit is the first after `eb74b94`, so it contains the whole project: the original implementation, pass 1 and pass 2. The review checks the cumulative tree. The claims in `docs/progress.md`, sections "Remediation pass 2" and "SAE verification", were treated as claims to verify. The procedure is `docs/rereview_prompt_2026-10-19_pass2.md`.

Every command was run on two clean detached worktrees of `954f87a`, one for pytest and one for run-all, so the two runs could not interfere. The untracked `src/bnr_non_equilibrium_2026/` stub and `docs/rereview_prompt_2026-10-19_pass2.md` are not part of the commit and were not present.

## Decision

**NO-GO.** Pass 2 closes the three P1 defects of the previous decision:

- **S1:** the headline ΔA now comes from the primary rows.
- **S2:** every C(t)-derived row of a contradicted pulse is QUALIFIED.
- **S3:** η is computed only on [0, 1].

Pass 2 also implements W3 to W7. The commands pass, and the committed outputs regenerate byte-identically. Two new P1 defects block acceptance:

1. **N1. The W5 model comparison is invalid as reported.**
   - All three alkalinity-limited candidates end their calibration with `cal_converged = False` in `model_fit.csv`: `least_squares` stopped without success.
   - `biomass_decay` converged and is "selected" against them.
   - Neither the report nor `progress.md` mentions the non-convergence.
   - The selected model, and so every MODEL-BASED quantity, is not established.
2. **N2. The idle fraction used in W4 is presented as a value, but it is a biased lower bound.**
   - Only 12 % of the record carries a state, and the denominator includes the 88 % without one.
   - The per-window Theil-Sen intervals that decide "C(t) flat" assume independent points. The effective sample size is about 3.7 of 101 points, so few bins can pass as flat and the inactive state is under-assigned.
   - The report quotes the idle fraction without a caveat in §1, §2, Fig. 8 and a GO status row.

Both are scientific-integrity risks in the headline, so the decision rule gives NO-GO. The other remaining items are P2 and are listed in the last table.

## 1. Commands (run by the reviewer, 2026-10-19)

| Command | Result |
| --- | --- |
| `uv run ruff check` | `All checks passed!`, exit 0 |
| `uv run pytest -q` | `381 passed in 4209.31s (1:10:09)`, exit 0. This includes the two fresh run-all fixtures and the byte-identical CSV/SVG determinism test. Wall clock under `time` was 11,602 s while sharing the CPU with run-all; pytest's own timer is quoted. |
| `rm -rf outputs && uv run python -m bnr_neq run-all` | Exit 0, 1990 s (33 min; it ran beside pytest). Stage messages: verdicts NH4/EC/ORP FAIL, pH PASS; segments references 11/13; thermo 9/11; dynamics 28/28; model selected `biomass_decay`, 13/15 identifiable; states 12 % classified, idle 3 %; pulses 265 computed, 71 not computed. These match `progress.md`. |
| `diff -rq` against a copy of the committed `outputs/` taken before the run | No differences in `tables/` (62 files), `reports/` or `figures/`. The SVGs and the PNGs are byte-identical. The committed outputs are reproducible, although `progress.md` says they were partly produced by stage re-runs (report, pulses) after the last full run. |

Environment: CPython 3.14.8, uv 0.13.1 (warns that `uv-build<0.13.0` does not contain it; harmless), numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2.

## 2. Previous remaining issues, R1–R12 and W1–W8

### 2.1 "Remaining issues by severity" of `docs/rereview_2026-10-19.md`

| Previous issue | Status | Evidence |
| --- | --- | --- |
| P0 R1: no sensor → state → energy → efficiency chain | **Resolved** (chain exists), with N2 against the idle link | `states.py` (stage between `model` and `pulses`, `__main__.IMPLEMENTED`); `energy.py`; `pulses.py:181-277` (energy family); Figs. 1, 4 and 8 in `figures.py`; report §1–§3. |
| P0 R2 / S1: headline ΔA from the raw reading | **Resolved** | `report.py:228-235` (`primary_delta_a`) filters `pulses.csv` on `nh4_observation_model == M6_log`. The status row reads ΔA = 0.60 to 1.15 kJ/mol, matching pulses.csv (1.151, 0.891, 0.833, 0.595). The SENSITIVITY values 2.36 to 4.17 are quoted beside it. `thermo_affinity_steps.csv` columns are now `delta_A_SENSITIVITY_reading_as_TAN_*` with a `role` column. |
| P1 S2: C(t) contradicts W2, unflagged | **Resolved** | `pulses.py:645-707` (`w2_consistency`, `qualify`, `check_w2_qualified`), called at `pulses.py:806-809`. In the committed `pulses.csv`, 0 M6_log rows of pulses 2–4 lack a QUALIFIED prefix. In `nh4_truth_windows.csv`, every row governed by pulses 2–4 is QUALIFIED. Hand recomputation in §4. |
| P1 S3: η > 1 reported as computed | **Resolved** (one P2 remnant) | `pulses.py:495-510`. η is computed only for 0 ≤ ratio ≤ 1. Pulses 2 and 4 at 1 mL/d show ratios 1.180 and 1.154, with η `not computed: removal exceeds the dose`. Wherever η is computed, `removal_to_dose_ratio` equals it (all 12 dose-anchored cases checked). Remnant: the pulse-2 η band edge is 1.069 (`param_band_high`), P2. |
| P1 R4 (W4), R5 (W3), R7 (W7), R8, R10 (W5); R6 not demoted | **Partly resolved** | See §2.2. R4 and R5 are implemented. R7 and R6 are resolved. R8 is resolved with a scope note (P2). R10 is **not** resolved: N1. |
| P1: aeration-equivalent energy from a stated SAE unverifiable | **Resolved** | `energy.aeration_energy_kWh`; config `sae_field_kgO2_per_kWh = [0.3, 1.0, 2.6]` with `sae_source`. Hand check and source check in §4–§5. |
| P2: `segments.csv` reading-based FA unlabelled | **Resolved for FA, not for the gain columns** | `fa_median_SENSITIVITY_reading_as_TAN_*`. But `segments.csv` `consumption_gain_corrected_mgN_L_d` (and its CI columns) and `thermo_oxygen.csv` `consumption_gain_corrected_mgN_L_d`, `o2_demand_mgO2_L_d_GAIN_CORRECTED` and `kla_required_per_d_SCENARIO_T25C` still lack SENSITIVITY in their names. The last of these has exactly the same name as the primary metric in `pulses.csv` (`thermo.py:502-523`). They are labelled only by a `role` or `gain_corrected_role` column. P2. |
| P2: `top_up_mL` silently 0; duplicate alkalinity keys | **Resolved** | `volume.py:74-84`: the net loss is MEASURED without the top-up, and the split is `not computed: top_up_mL missing`. `alkalinity.py:150-162` (`alk0_source`): `alkalinity_mmol_L` is never used as Alk0, and the status says so in `balance_alkalinity.csv`, `alkalinity.md` and report §3. Tests `test_missing_top_up_is_a_status_not_zero` and `test_older_alkalinity_key_is_never_silently_ignored`. |
| P2: evaporation bound loose; M6 not refitted per volume | **Partly resolved** | `volume_constraints.csv` adds N conservation (≥ 0.314 mL/d over pulse 3, CONDITIONAL on C(t)), liquid present (< 2.09 mL/d) and the EC bound (0.040 mL/d, CONDITIONAL, not used). The M6 refit per scenario is still not done, which is stated debt (P2, accepted in §9). |
| P2: R11 readability; R12 placeholders | **R12 resolved; R11 partly** | `pyproject.toml:4,21,28` hold a real description and `bnr-neq = "bnr_neq.__main__:main"`. The report is narrative-first and the new tables are rounded. The M6 consistency-check and identifiability tables in report §4 are still unrounded (`65.46428371127854`, `4.547473508864641e-13`, `19.999828959312286`). P2. |

### 2.2 Findings R1–R12 and work items W1–W8

| Item | Status | Evidence |
| --- | --- | --- |
| W1 single NH4 model | **Resolved in the outputs; guard incomplete (P2)** | Every primary row of `pulses.csv` uses `M6_log`/`C_inv:58628090efe8` (244 rows). The 16 MODEL-BASED rows use `M6_state:biomass_decay`. SENSITIVITY rows use `electrode_reading`. The 36 rows without a model are not NH4-derived: `t_pulse_d`, `segment_end_d`, `median_pH`, O2 supply capacity, free energy per fully converted pulse, stirrer status. See §6 for the checker probes. |
| W2 alkalinity module | **Resolved** | Unchanged from pass 1. Positive values re-checked by hand in §4. |
| W3 state classifier | **Partly resolved (N2)** | Rules, excluded-feature guard (`states.validate_rules`), artifact mask, block bootstrap and time-ordered table all exist and are tested. The coverage and the flat-slope test are judged in §8. |
| W4 energy family | **Resolved, except the idle link (N2)** | η, the 4.57 and 3.43 O2 cases, supply capacity, utilization, demand floor, supply-based energy, kWh/kg N, ΔG/energy, idle fraction and relaxation time are all in `pulses.csv`. Every SAE- or kLa-dependent metric carries SCENARIO in its name. A negative or zero SAE raises (`energy.py:69-73`, config). |
| W5 alkalinity / IC candidate | **Unresolved (N1)** | Candidates are in `model.py:148-166` and `model_fit.csv`. Calibration did not converge for any of them (`cal_converged = False`; `model.py:384` returns `sol.success`), so the selection at `model.py:431` compares converged with unconverged fits. |
| W6 Figs. 1, 4, 8 | **Resolved** | `figures.py`; inspected visually. Fig. 1 reads the AE range from the configuration. |
| W7 narrative report | **Resolved (R11 remnants P2)** | §1 "What the data say" first, energy §2, sensors §5, estimators in Supplement S1. Coexistence affinity, ΔA and σ ≥ 0 are `CONSISTENCY CHECK (not a finding)`. |
| W8 volume | **Resolved (stated P2 debt)** | See §2.1. |
| R1 | Resolved (N2 qualifies the idle link) | — |
| R2 | Resolved | — |
| R3 | Resolved | The bound is stated in §1 point 3. The contradiction is stated and qualified (§1 point 4, `nh4_truth_w2_consistency.csv`). |
| R4 | Resolved, except N2 | — |
| R5 | Resolved as a test | Relative ORP and pH are used as indicators, with an NH4-free association test. The result is stated as "no indicator separates the consumption reliably". |
| R6 | Resolved | The estimators are moved to Supplement S1, not removed. |
| R7 | Resolved | Relabelled as consistency checks. A trajectory-based time scale, `relaxation_time_d` = dose / post-pulse rate, is added. The name overstates it (it is not a fitted relaxation constant), P2. |
| R8 | Resolved, with a scope note (P2) | `segments.py:28-32, 188-214`. Sub-segments 2, 3, 5, 7 and 9 are demoted, and their slopes are kept in `slope_fragment_diagnostic_mgNH4_L_d`. See §8. |
| R9 | Partly resolved (stated P2 debt) | — |
| R10 | **Unresolved (N1)** | Time confounding is stated. The competing candidate exists, but its fit is unconverged. |
| R11 | Partly resolved (P2) | — |
| R12 | Resolved | — |

## 3. End-to-end trace, pulse 2 (`outputs/tables/pulses.csv`)

| Metric | Value | `nh4_observation_model` | `nh4_series` | `w2_bound_check` | Hand check |
| --- | --- | --- | --- | --- | --- |
| `consumption_weighted_mgN_L_d` | 2.599 | M6_log | C_inv:58628090efe8 | QUALIFIED (pre_pulse3) | — |
| `fa_median_mgNH3_L_SCENARIO_T25C` | 3.783 | M6_log | C_inv:58628090efe8 | QUALIFIED | — |
| `eta_conversion_fraction` = `removal_to_dose_ratio` | 0.809 | M6_log | C_inv:58628090efe8 | QUALIFIED | (C_pre + 65.464 − C_next)/65.464 |
| `o2_demand_mgO2_L_d` | 11.878 | M6_log | C_inv:58628090efe8 | QUALIFIED | 4.57 × 2.5992 = 11.878 ✓ |
| `n_removed_mgN` → `o2_required_mgO2` | 5.298 → 24.21 | M6_log | C_inv:58628090efe8 | QUALIFIED | 0.80930 × 6.5464 = 5.298; × 4.57 = 24.21 ✓ |
| `aeration_energy_supply_kWh_per_kgN_SCENARIO_SAE1_kLa10_T25C` | 21.86 | M6_log | C_inv:58628090efe8 | QUALIFIED | 115.805 mg O2 / 1 / 5.298 mg N = 21.86 ✓ |
| `idle_fraction_of_stirred_time` | 0.010 | M6_log | C_inv:58628090efe8 | QUALIFIED | 14 inactive of 1402 bins (`states_fractions.csv`); a lower bound, see N2 |
| `relaxation_time_d` | 25.19 | M6_log | C_inv:58628090efe8 | QUALIFIED | 65.464 / 2.5992 = 25.19 ✓ |
| `eta_conversion_fraction_MODEL_BASED` | 0.569 | M6_log | M6_state:biomass_decay | QUALIFIED | labelled MODEL-BASED ✓ (rests on the unconverged comparison, N1) |
| `o2_supply_capacity_mgO2_SCENARIO_kLa10_T25C` | 115.8 | (none) | (none) | (none) | 10 × 8.26 × 0.1 × 14.02 = 115.8 ✓; independent of NH4 |

Every C(t)-derived quantity of pulse 2 rests on one series, or is labelled MODEL-BASED. The `w2_bound_check` of pulse 2 cites checkpoint pre_pulse3, the end of its interval. In `nh4_truth_w2_consistency.csv` that checkpoint has `consistent_at_best_fit = False`, while pre_pulse2 is consistent (13.1 meq/L). The label is consistent with the table ✓.

## 4. Hand recomputation

The constants used were taken independently of the code: 25 °C, pK1 = 6.352, pK2 = 10.329, log KH = −1.468, pCO2 = 4.2e-4 atm, M_N = 14.007, Kw = 1e-14.

- **W2 positive tests.**
  - N_nit = M_N × Alk0/2 gives 14.007, 35.017 and 70.035 mg N/L for Alk0 = 2, 5 and 10 meq/L ✓.
  - HCO3⁻ = KH·pCO2·10^(pH−pK1) gives 0.4011 mM at pH 7.8 and 0.6357 mM at pH 8.0 ✓.
  - FA = 17/14 × 200 × 10^7.8/(e^(6344/298.15) + 10^7.8) = 8.492 mg NH3/L ✓.
- **Alk0 requirement, pre_pulse3** (`nh4_truth_w2_consistency.csv`).
  - Removed per initial litre: 177.115 + 2 × 65.464 − 163.876 = 144.168 mg N/L ✓.
  - Alk_min at the smoothed pH 7.5754 is 0.2403 meq/L (HCO3⁻ + 2 CO3²⁻ + OH⁻ − H⁺).
  - Requirement if nitrified: 2 × 144.168/14.007 + 0.2403 = **20.825 meq/L** (table 20.8255 ✓). For a 1 eq/mol N sink: 10.533 (table 10.5330 ✓).
- **Aeration-equivalent energy, pulse 1, kLa 10/d, SAE 1.**
  - O2 capacity: kLa × C* × V × Δt = 10 d⁻¹ × 8.26 mg/L × 0.100 L × 2.04 d = 16.850 mg O2.
  - N removed: 0.380305 × 6.54643 mg = 2.48964 mg N.
  - Energy: (16.850e-6 kg O2 / 1 kg O2 kWh⁻¹) / 2.48964e-6 kg N = **6.768 kWh/kg N** (table 6.768207 ✓).
  - At SAE 0.3 and 2.6: 22.561 and 2.603 ✓.
- **Demand floor.** 4.57/SAE = 15.233, 4.570 and 1.758 kWh/kg N at SAE 0.3, 1 and 2.6 ✓. Units: g O2/g N ÷ kg O2/kWh = kWh/kg N.
- **`relaxation_time_d`.** Pulse 1: 65.4643/7.15314 = 9.152 d, and 9.152/2.04 = 4.486 ✓. Pulse 2: 25.187 d ✓.

## 5. SAE source

Both documents were retrieved and read on the pages cited.

- **Stenstrom M.K. and Rosso D., *Aeration*, UCLA/UCI, 8 Feb 2010, p. 56 ("Efficiency").**
  - Clean-water SAE: high-speed surface 0.9–1.3, low-speed surface 1.5–2.1, coarse bubble 0.6–1.5, turbines/jets 1.2–1.8, fine pore 3.6–4.8 kg O2/kWh.
  - Process-water AE at 2 mg/L DO, low / high SRT: 0.4–0.8; 0.7–1.5; 0.3–0.7 / 0.4–0.9; 0.4–0.6 / 0.6–0.8; 0.7–1.0 / 2.0–2.6.
  - The configured 0.3 and 2.6 are the extremes of the process-water columns ✓. They are process-water, not clean-water, values ✓.
  - p. 55 warns that the table "should not be used as a general guide line and not for design". Using it as a SCENARIO range is compatible with that, but the config comment should say so.
- **Moore L.W., *Biological Wastewater Treatment Training Series #11: Aeration Systems*, ORNL, July 2021, slides 41–42.**
  - Field transfer rates of 1.0–2.3 lb O2/(hp h) at α = 0.84, β = 0.92, DO = 2 mg/L.
  - 1 lb/(hp h) = 0.4536/0.7457 = 0.6083 kg/kWh, so the range is 0.61–1.40 kg O2/kWh ✓. Slide 2 lists Metcalf & Eddy 4th ed. as a reference ✓.
- **Naming.** SAE conventionally means the clean-water standard efficiency (20 °C, zero DO). The values here are field AE. The key `sae_field_kgO2_per_kWh` carries "field", and the config comment and report §2 state the distinction, so this is **acceptable**, but the metric suffixes `SCENARIO_SAE0.3` etc. would read better as `AE`. One conceptual inconsistency should be stated: the supply-based energy divides a vessel capacity computed at DO = 0 by an AE defined at DO = 2 mg/L. P2.

## 6. New modules, stages, fallbacks and schema

- **Stage order.** `__main__.IMPLEMENTED` places `states` between `model` and `pulses`, and `test_every_command_of_section_7_is_implemented` asserts it ✓.
- **Missing inputs.**
  - `stirrer_power_W`, `kLa_per_d`, `states_fractions.csv`, `initial_alkalinity_meq_L` and `top_up_mL` give status strings ✓.
  - The one default number is DO = 0 in the supply capacity (`pulses.py:188`). It is stated in the basis of every row as "the transfer capacity … an upper bound on the supply", which makes it a labelled bound, not a silent fallback ✓.
- **Scenario labels.** Every kLa or SAE metric name carries `SCENARIO_` ✓. The infeasible kLa = 1/d cases (pulses 1 and 2) are `not computed: SCENARIO infeasible` ✓ (judgment call 3).
- **W1 provenance checker** (`model.check_single_observation_model`, `model.py:61-112`). It was probed by appending rows to the committed `pulses.csv` (scratch script; results are reproducible).
  - **Caught:** a second `C_inv:` series; an M6_log row with the reading as series; a second `M6_state:` series; a model-state series without a MODEL-BASED label; a MODEL-BASED label on an inverted series; a trailing-space model name; a derived ΔA from the raw reading.
  - **Not caught:**
    1. A derived row (`fa_median_mgNH3_L_SCENARIO_T25C`, `idle_fraction_of_stirred_time`) with an empty model column, because `is_primary` excludes empty models (`model.py:86`).
    2. W4 and W1 quantities from the raw reading: `aeration_energy_supply_*`, `n_removed_mgN`, `removal_to_dose_ratio` and `consumption_weighted_mgN_L_d`, each labelled `electrode_reading`. `DERIVED_PREFIXES` (`model.py:52`) was not extended when pass 2 added the energy family.
  - **False positive:** a basis containing "not MODEL-BASED".
  - The checker runs only on `pulses.csv` (`pulses.py:578`), not on `nh4_truth_windows.csv`, `states_*` or `segments.csv`.
  - No committed output violates the invariant. The holes are latent, so they are rated **P2** with a concrete fix: raise on an empty model when the metric is NH4-derived, extend the prefixes, and apply the checker to every C(t)-carrying table.
- **`check_w2_qualified`** is applied right after `qualify` on the same frame, so on the production path it guards only against later edits. The pre-pulse window 1 of `nh4_truth_windows.csv` (1.1–2.63 d, computed) has an empty `w2_bound_check`, because qualification is keyed on the governing pulse. pre_pulse1 is consistent, so nothing is mislabelled, but the cell should say so. P2.
- **Documentation.** All 62 tables are described in `report.py` `TABLE_DESCRIPTIONS` ✓. The new columns `w2_bound_check`, `removal_to_dose_ratio`, `slope_fragment_diagnostic_mgNH4_L_d`, the band columns and `role` are documented in `progress.md`, and their readers (report, figures, tests) use the new names ✓.

## 7. Test adequacy

- **New modules.**
  - `test_states.py`: excluded sensor raises, failed pH excludes everything, all-artifact interval gives no state, a block-bootstrap shift is recovered, AR(1) widening, determinism, small-group guard, real-record held spans.
  - `test_energy.py`: hand values, non-positive SAE raises in the function and the config, a missing source raises, real-record recomputation, labels, no utilization above 1.
  - `test_rereview_fixes.py`: series mixing, MODEL-BASED rules, S2 detection, qualification and raise on synthetic and real tables, η on [0, 1], top-up status, N-conservation hand value, EC bound, old alkalinity key, R8 rule and real fragments, W5 specs and stall.
  - Together these cover the positive, negative, boundary, integration and end-to-end cases on the production functions and the `run_all_dir` fixture ✓.
- **Gaps.**
  1. No test checks that the W5 calibrations converge (`test_model_fit_lists_the_alkalinity_candidates` checks only presence), which is why N1 passed.
  2. No test checks the idle fraction against coverage, or that the slope interval respects autocorrelation (N2).
  3. `test_reading_based_columns_say_sensitivity` covers ΔA and FA but not the gain-corrected columns.
  4. The checker tests do not probe empty-model or W4-metric rows.
  5. `test_truth.test_dose_anchored_eta_and_affinity` branches on the computed ratio (`if 0 ≤ η ≤ 1 … elif η > 1 … else`) instead of asserting the branch the synthetic truth must take. It cannot fail by picking the wrong branch. That makes it weaker than a pinned assertion, but it does not relax a requirement.
- **Changed tests.** Pre-pass-2 versions are not in history (`eb74b94` holds only `.gitignore` and `README.md`), so each edit was judged against the requirement and the `progress.md` rationale:
  - `test_main`: stage list and order follow the new `states` stage ✓.
  - `test_report`: consistency-check labels, new W2/S2/W8 rows, 8 figures, and the fixture extended. The S1 test proves that `thermo_affinity_steps.csv` no longer feeds the headline (×10 tampering) ✓.
  - `test_segments::test_scenario_labels`: the renamed FA columns, asserted by exact name ✓.
  - `test_truth`: computed η must lie in [0, 1], and non-computed reasons are restricted to two strings ✓.
  - None of these lowers a threshold or removes an assertion that the old requirement still demands. **No weakened test found** (the branching η test above is a weakness, P2).

## 8. Scientific correctness

- **W3 rules.**
  - The ±15 mV ORP threshold is RT ln10/4F, one decade of p_O2 on a Nernstian 4-electron couple (`states.py:43`). It is a defensible physical scale, and the code says Pt is not strictly Nernstian. It was fixed from physics, not from the outcome, as far as can be judged: the tree has one commit, so the order of decisions cannot be verified from history.
  - Excluding the NH4 slope from the association partition removes the direct circularity ✓. The comparison of full states is labelled DESCRIPTIVE (circular) ✓.
  - **Block length.** The response ACF on the regenerated `states_bins.parquet` is −0.055 at lag 101, so the 101-bin block (≥ the 1.01 d window, ≥ the first lag with ACF < 1/e) captures the autocorrelation of the response ✓.
  - **"≥ 2 blocks per group".** It was added after the first run, which is disclosed. It is a validity condition: with fewer than two blocks the percentile interval cannot represent between-block variability, and it only withholds a result (the ORP comparison). It does not change the direction of any claim. Accepted, but two blocks is a minimal bar.
- **N2. The slope test that defines the states ignores autocorrelation (P1 with the coverage).**
  - `states.py:153` uses scipy `theilslopes` intervals, which treat the 101 bins as independent.
  - On the regenerated bins, the detrended C(t) inside slope windows has a median lag-1 autocorrelation of 0.93, so n_eff ≈ 101 × 0.07/1.93 ≈ 3.7. The intervals are therefore about √(101/3.7) ≈ 5× too narrow.
  - Of 2,031 bins with a slope, only 379 (19 %) are "flat". 838 are "falling" and 814 "rising"; the rising bins go to `unclassified`.
  - "Inactive, over-supplied" and "oxygen-limited" both require a flat C(t), so they are under-assigned by construction. Together with 12 % coverage, `idle_fraction_of_stirred_time` (whole record 0.03; P1 0.00, P2 0.01, P3 0.008, P4 0.08) is a lower bound with an unknown, downward bias.
  - The report states it as a value in §1 point 8, the §2 energy table, Fig. 8c and the W3 status row (GO WITH DEBT, "idle 3 % of stirred time").
  - **Required:** an autocorrelation-aware slope interval (n_eff-scaled or block-based), and the idle fraction reported as a lower bound with its coverage beside it, or not computed below a stated coverage.
- **Coverage.** The state layer supports one claim: at the registered thresholds, relative ORP and pH do not separate NH4 consumption (pH rising vs not: +1.35, 95 % block bootstrap −0.53 to +2.57; ORP: not computed). It does not support an idle fraction as an efficiency metric.
- **W5 (N1, P1).**
  - `alkalinity_limited` has better test skill (0.817 against 0.748) and forecast skill (0.925 against 0.862), but lower selection skill (0.934 against 0.958). The report's conclusion that decay cannot be discriminated from alkalinity limitation is stated prominently (§1 point 7) and is, if anything, conservative.
  - But all three alkalinity candidates have `cal_converged = False`. The selection rule (`model.py:431`, minimum selection RMSE) therefore compares a converged model with fits stopped at the `least_squares` evaluation limit.
  - The "selected" label, and every MODEL-BASED quantity (MODEL-BASED η 0.205, 0.569, 0.070, 0.035), rest on that unequal comparison. The fitted Alk0 of the candidates, which would be a direct cross-check of S2 (requirement ≥ 20.8 meq/L), is not reported.
  - **Required:** converge the W5 fits (more evaluations, warm starts, or a reparametrization), report `cal_converged` and the fitted Alk0 in the report, and treat the selection as provisional until then.
- **R8.**
  - The rule ("baseline window between two electrode steps shorter than 3 × the longest delayed electrode response", 3 × 0.33 d) is physically motivated. Three settling spans of a first-order transient is an analogy here, since the 0.33 d is a delay, not a time constant (M6 τ = 0.26 d gives 3τ = 0.79 d, similar).
  - It removes only pre-pulse windows bounded by offset jumps or candidates, which are artifacts by construction ✓.
  - Its scope excludes post-pulse windows. Sub-segment 14 (8.96–9.27 d, 0.31 d) starts 0.10 d after the delayed pulse-1 response jump, which is within 3τ by the code's own `INVERSION_NOTE`. It is retained as the first window of pulse 1, and it supplies the highest point of Fig. 4b (18.3 mg N/L/d). The asymmetry should be justified, or window 14 shown as a sensitivity. P2.
- **Report.**
  - The coexistence affinity, ΔA and σ ≥ 0 are consistency checks, not findings ✓. The associations are labelled "not a cause" ✓.
  - Every number in §1 traces to a table: verdicts → `sensor_verdicts.csv`; 177/294 and bands → `nh4_truth_summary.csv`; 34/104 and FA 2.4–3.9 → `alkalinity_bound.csv`; 145 and 20.9/11.5/10.5 → `nh4_truth_w2_consistency.csv`; η and relaxation ratios → `pulses.csv`; state fractions and bootstrap → `states_*.csv`; skills → `model_fit.csv`; 1.76–15.23, 6.8, 21.9 and 99.7 → `pulses.csv`; 0.31 and 2.09 → `volume_constraints.csv` ✓.
  - Two wording issues: §1 point 2 says "open upward" only for the end, but `band_high_open = True` applies to the whole series. §1 point 5 ("TAN accumulated and the batch was driven away from its inactive reference") is an interpretation of QUALIFIED numbers, covered only by the blanket statement in point 4. P2.
- **Other.** The MODEL-BASED η is about half the dose-anchored η for every pulse (0.21 against 0.38, 0.57 against 0.81, 0.035 against 0.20). Both rest on the same observation model, and the disagreement is not discussed. P2.

## 9. Rulings on the implementer's judgment calls

| # | Judgment call | Ruling | Reason |
| --- | --- | --- | --- |
| 1 | W1: the inverted rate is flagged QUALIFIED where it contradicts the W2 bound, not constrained to obey it | **Accept** | Forcing the bound would fit the data to an assumption: the scan maximum of 15 meq/L and nitrification as the only sink. The table shows that the contradiction is resolvable within the open C(t) band (≥ 10.7 meq/L) or by a 1 eq/mol N sink. Note that "consistent" at pre_pulse1 and pre_pulse2 still needs 9.5 and 13.1 meq/L, which is high for most wastewaters. The verdict states the requirement, so this is acceptable. |
| 2 | S3: η withheld above 1; removal-to-dose ratio reported separately | **Accept** | A ratio above 1 has a physical reading (carry-over TAN removed), and η as a fraction of one pulse is undefined there. Remnant: the η band edge of 1.069 for pulse 2 should be clipped or moved to the ratio row (P2). |
| 3 | Infeasible kLa scenarios not computed | **Accept** | A utilization above 1 or a supply energy below the demand floor would be physically meaningless. Infeasibility is itself informative and is stated in the status (required/capacity ratio). |
| 4 | Conditional EC bound reported, not used | **Accept** | The probe failed M2 (no response to the pulses), so using its flat record to cap evaporation at 0.04 mL/d would reuse a rejected sensor. Reporting it as CONDITIONAL is correct. |
| 5 | M6 not refitted per evaporation scenario | **Accept as P2 debt** | The scenarios are carried through every balance, η and the rates. The fit itself keeps the constant dose. This changes the C(t) level by at most the concentration factor (≤ 1.92 at 1 mL/d by day 48), which is inside the open band. It must stay listed as debt. |

## Remaining issues by severity

| Severity | Issue |
| --- | --- |
| P1 | N1: the W5 calibrations did not converge (`cal_converged = False` for all three alkalinity candidates). The selection of `biomass_decay` and all MODEL-BASED quantities rest on an unequal comparison, and this is unreported. The fitted Alk0 is not reported. |
| P1 | N2: the idle fraction (W4) is reported as a value, but it is a biased lower bound: 12 % coverage, and slope intervals about 5× too narrow (n_eff ≈ 3.7 of 101) under-assign the flat states. It appears in §1, §2, Fig. 8c and a GO status row. |
| P2 | W1 checker holes: an empty model column passes; W4/W1 metrics (`aeration_energy_*`, `n_removed`, `removal_to_dose_ratio`, `consumption_*`) from the raw reading pass; the checker runs on `pulses.csv` only. |
| P2 | Gain-corrected columns without SENSITIVITY in their names (`segments.csv`, `thermo_oxygen.csv`). `thermo_oxygen.csv` `kla_required_per_d_SCENARIO_T25C` shares the primary metric's name. |
| P2 | η parameter band above 1 (pulse 2, 1.069). The pre-pulse truth window has an empty `w2_bound_check`. |
| P2 | R8 scope: the 0.31 d post-pulse window 14 is retained, though it starts within 3τ of the delayed response. |
| P2 | R11: unrounded M6 consistency-check and identifiability tables in report §4. §1 wording on the open band and on accumulation. `relaxation_time_d` naming. |
| P2 | AE naming (`SAE` suffix for field AE); supply energy mixes a DO = 0 capacity with an AE at DO = 2 mg/L; the source's "not for design" caveat is not quoted. |
| P2 | W8: the M6 fit is not refitted per volume scenario (stated debt). |
| P2 | `test_truth` η test branches on its own result; no convergence test for the W5 fits; no test of the gain-column names. |
| P2 | MODEL-BASED η about half the dose-anchored η, not discussed. |

## What is sound and should be kept

- Byte-reproducible outputs from a clean checkout, including the PNG figures.
- One inverted C(t) with a hashed identifier, and a checker that now also checks series and MODEL-BASED labels.
- The S2 qualification with per-checkpoint verdicts that state their three possible resolutions, and a guard that raises.
- The η/ratio split, the infeasibility statuses, and the explicit `not computed` paths for every missing measurement.
- An energy family whose values reproduce by hand, with a verified, correctly characterized process-water AE source.
- An NH4-free association test with block bootstrap and a validity guard, and a report that calls associations associations.
