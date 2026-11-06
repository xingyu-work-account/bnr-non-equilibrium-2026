# Independent re-review of remediation pass 6 and W6/W7 (2026-11-05)

This review covers the five commits on branch `remediation/pass6-2026-10-26` from `2dd5a2b` (code and outputs equal to the merged pass 5, `0f7abfc`) to `0178cd8`: `f1744fc` (G1/G2), `ba394df` (G3, G4, G6–G8), `d6d9e5d` (G5 rule), `da49e95` (G5 range) and `0178cd8` (W6/W7). The procedure is `docs/rereview_prompt_2026-11-05_pass6.md`. The claims in `docs/progress.md` ("Remediation pass 6", "Remediation pass 6, G5", "W6 and W7") were treated as claims to verify.

Every command was run on a clean detached worktree of `0178cd8`. Scratch checks (mutations, reproductions) ran on separate `git archive` extracts of the relevant commits, never in the worktree used for run-all and pytest. The untracked `src/bnr_non_equilibrium_2026/`, `reference/` and `docs/reports/` are not part of any commit and were not present in the worktree.

## Decision

**GO WITH DEBT at `47dfd18`.** At `0178cd8` the decision was **NO-GO**, for one P1 only; the P1 was fixed in `47dfd18` and the fix is verified in the addendum at the end of this review. The debt is P2-1 to P2-12 below and the two open pass-5 P2 items.

At `0178cd8`, every command passed: ruff clean, run-all exit 0 with the expected stage results, regenerated outputs byte-identical to the committed ones at a different worker count, and 491 tests passed with no skips on code whose hashes did not change. The accepted values are protected by the guards in §2, and G1–G8, W6 and W7 close in substance (§3).

The blocking issue at `0178cd8` was **P1-1**: Fig. 1 drew the idle-fraction State box ("idle 5-93 %") and the W2-bound box with the blue border its own legend defines as "measured", although both values depend on SCENARIO settings (the ORP threshold and the model-inverted C(t); the Alk0, temperature and kLa_CO2 SCENARIOs). That is a SCENARIO-dependent value presented as measured, which `docs/remediation_prompt_2026-10-19.md` names as a failure condition. `47dfd18` draws both boxes as scenario and adds a guard and tests; it changes no table.

## 1. Commands (run by the reviewer, 2026-11-05)

| Command | Result |
| --- | --- |
| `uv sync` in a clean detached worktree of `0178cd8` | `git status --short` empty before and after the run |
| SHA-256 of `src/`, `tests/`, `config/experiment.toml`, `docs/states_registration.md`, `pyproject.toml`, `uv.lock` | 60 files hashed before run-all; compared after pytest (row below) |
| `uv run ruff check` | `All checks passed!`, exit 0 |
| `rm -rf outputs && uv run python -m bnr_neq run-all` (`BNR_NEQ_WORKERS=6`) | Exit 0, wall 2:51:43 (real 10303.44 s). Stage messages: verdicts NH4/EC/ORP FAIL, pH PASS; segments references 11/13; thermo 9/11; dynamics 28/28; model selected `biomass_decay`, consistency checks pass, 13/15 identifiable; states 12 % classified, idle between 5 and 93 %; pulses 274 computed, 78 not computed. These match `progress.md` and pass 5. The only stderr lines are matplotlib's `findfont` weight-600 substitutions, disclosed by the implementer (ruling 6). |
| `diff -rq -x .DS_Store` against a copy of the committed `outputs/` taken before the run | Exit 0, no differences: `tables/` (67 files), `reports/` (10) and `figures/` (116 PNG and SVG, including `no_text/`). `processed/` (9 files, not committed) is byte-identical to the implementer's last run in the main checkout. The outputs do not depend on the worker count (6 here, all 10 cores for the implementer). |
| `uv run pytest -q -rs` (`BNR_NEQ_WORKERS=6`, after run-all) | `491 passed in 30959.44s (8:35:59)`, exit 0, none skipped (`-rs` lists none). This matches the implementer's count. The wall time is long because another full test suite ran in the main checkout at the same time (load average about 150 on 10 cores); that run was not the reviewer's and did not touch the worktree. |
| SHA-256 recomputed after pytest | Identical to the hashes taken before run-all (`diff` exit 0, 60 files). `git status --short` in the worktree is empty. The code that passed is the committed code of `0178cd8`. |

Environment: CPython 3.14.8, uv 0.13.6, numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2, 10 cores. The implementer used the default worker count (all cores); the reviewer used 6. Ruff, run-all and pytest ran one after another, never concurrently.

## 2. Guards: nothing accepted may change

### 2.1 G1–G4, G6–G8 against `2dd5a2b` (committed outputs, `2dd5a2b` → `ba394df`)

Raw-text, column-by-column comparison of every table (`csv` module, no type parsing):

| Table | Rows | New columns | Pre-existing columns changed |
| --- | --- | --- | --- |
| `model_fit.csv` | 9 = 9 | `nh4_observation_model`, `nh4_series`, `assumptions` | none |
| `model_calibration_starts.csv` | 33 = 33 | same three | none |
| `model_identifiability.csv` | 15 = 15 | same three | none |
| `model_checks.csv` | 8 = 8 | same three | none |
| `model_volume_scenarios.csv` | 27 = 27 | same three | none |
| `model_volume_scenario_summary.csv` | 3 = 3 | `nh4_observation_model`, `nh4_series`, `constant_volume_nh4_series`, `assumptions` | none |
| `model_response_fractions.csv` | 4 = 4 | `model_response_fraction_nh4_observation_model`, `model_response_fraction_nh4_series`, `observed_response_fraction_nh4_observation_model`, `assumptions` | none |
| `nh4_truth_vectors.csv` | 86 = 86 | `nh4_observation_model`, `nh4_series`, `assumptions` | none |
| `nh4_truth_w2_consistency.csv` | 5 = 5 | `nh4_observation_model`, `nh4_truth_id`, `assumptions` | none |
| `pulses.csv` | 352 = 352 | `assumptions` | none |
| `states_fractions.csv` | 6 = 6 | `assumptions` | none |
| `states_comparison.csv` | 7 = 7 | `registration_sha256`, `registration_status`, `null_coverage_min`, `interval_reading`, `assumptions` | none |
| `states_rules.csv` | 5 = 5 | `registration_sha256`, `registration_status` | `thresholds` only (5 rows): the old text plus the suffix " (SCENARIO 25 °C: RT ln10/(4F) at 298.15 K = 14.8 mV, rounded; the liquid temperature is not measured)" |

- One new table, `states_bootstrap_null_coverage.csv` (12 rows). No table removed. Every other table byte-identical.
- The new columns match the list in `progress.md`, including `model_checks.csv`, which `progress.md` lists in the "Changed outputs" table but not among the seven G8 tables.
- **Column order.** In 11 tables the new columns are appended after the old ones. In `states_fractions.csv` and `states_comparison.csv` they are inserted **before** the last two pre-existing columns (`nh4_observation_model`, `nh4_series`), which therefore move from positions 19–20 to 20–21 and from 15–16 to 20–21. The relative order of the old columns is kept, and every reader in the repository selects by name, but `progress.md`'s "column order kept" is not literally true for these two tables (P2-4).
- Figures: only `fig4_constraints.{png,svg}` changed (panel-c title), as disclosed.
- **`processed/` cannot be compared with the `2dd5a2b` commit:** `outputs/processed/` is in `.gitignore`, so no baseline exists in git. It was compared with the implementer's last run in the main checkout instead (§1).
- **Source.** An AST comparison of `git diff 2dd5a2b ba394df -- src` (docstrings stripped) shows only additions in `states.py` (`registration`, `registration_sha256`, the null-check functions and constants), `model.py` (assumption constants, `m6_tables`, `check_m6_provenance`, `with_candidate_provenance`), `pulses.py` (`row_assumptions`, `label_assumptions` and constants), and changes in `compute`, `run`, `write_report`, `rules_table`, `comparisons`, the report and fig4. I read each changed function: the changes assign columns, call the null check, render text, or check provenance before the existing `to_csv`. `comparisons` gains two `assign` lines only. No fit, bound, start, selection rule, state rule or primary computation changed. `figures.fig4` changes only the panel-c title string.

### 2.2 G5 against `ba394df` (`ba394df` → `da49e95`)

- Only `pulses.csv` changed among the tables; columns and row count identical.
- **The 244 rows without an AE suffix are byte-identical** (whole row as raw text).
- **The 108 AE rows:** the metric name equals the old name with AE0.3/AE1/AE2.6 mapped to AE0.7/AE1.6/AE3.5, row for row; every column other than `metric`, `value` and `basis` is byte-identical; `basis` differs only in the AE source and DO wording.
- **Recomputed scaling**, all 69 computed AE values (the other 39 are "not computed" in both versions with identical text): AE_old/AE_new for energy, AE_new/AE_old for the Gibbs-to-aeration ratios. Largest relative error **4.2e-16**. Per kind, for example:

  | Kind (rows) | Example | Old | × factor | New |
  | --- | --- | --- | --- | --- |
  | `aeration_energy_demand_Wh` (12) | AE0.3 → AE0.7 | 0.0379268 | 0.3/0.7 | 0.0162544 |
  | `aeration_energy_demand_kWh_per_kgN` (12) | AE1 → AE1.6 | 4.57 | 1/1.6 | 2.85625 |
  | `aeration_energy_supply_kWh_per_kgN` (36) | AE1, kLa 10 → AE1.6 | 6.76797 | 1/1.6 | 4.22998 |
  | `gibbs_to_aeration_energy_demand` (12) | AE0.3, 25 °C → AE0.7 | 0.454082 | 0.7/0.3 | 1.05952 |
  | `gibbs_to_aeration_energy_supply` (36) | AE1, kLa 10 → AE1.6 | 1.02205 | 1.6/1 | 1.63527 |
- Figures: `fig1_linkage` and `fig8_efficiency` only (fig8 not predicted by the registered decision, disclosed). Reports: `pulses.md`, `analysis_report.md`.
- Source: `pulses.py` and `report.py` change only text (`pulse_metrics`'s `sae_basis` string, `write_report`, `assumptions`, `assemble_report`); `config/experiment.toml` changes the range and source text.

### 2.3 W6/W7 against `da49e95` (`da49e95` → `0178cd8`)

- **No table changed** (raw-text comparison of all tables: no difference).
- **Function by function** (AST comparison, docstrings stripped): in `dynamics.py`, `model.py`, `pulses.py`, `segments.py` and `validation.py`, the only removed top-level objects are `figure`, `_style`, `_optional_table` and the colour constants; the only changed function is `run`, where the old `figure(...)` call becomes `figures.<module>_figures(cfg, out_dir)`, placed after every `to_csv`/`to_parquet` of that stage. The import lines gain `figures` (and `validation.py` drops `pyplot`). No computation changed. The two `model.py` docstring edits are "§4" → "§5" cross-references. `report.py` changes `FIGURES`, `summary_lines` and `assemble_report` (order, captions, rounding), and `__main__.py` adds the `figures` command.
- Because the figures are now drawn from the saved tables rather than from in-memory results, a figure can only show what a table holds. That is an improvement for traceability.
- **Old stems.** `git grep -E 'fig1_linkage[."]|fig2_overview|fig3_sensor_validation|fig4_constraints|fig5_dynamics|fig6_model|fig7_pulses|fig8_efficiency' 0178cd8 -- src tests outputs README.md` finds nothing (and no old file remains under `outputs/figures/`). The old stems survive only in `docs/progress.md` and earlier review documents, as history.

## 3. Gap-by-gap closure

| Gap | Status | Evidence |
| --- | --- | --- |
| G1 | **Closed** (P2-1, P2-2) | See 3.1. |
| G2 | **Closed** | See 3.2. |
| G3 | **Closed** | RT ln10/(4F) at 298.15 K = 14.7898 mV → "14.8 mV, rounded", and `ORP_LEVEL_MV = round(...)` = 15. The SCENARIO text is in `states_rules.csv` `thresholds`, in `assumptions` of `states_fractions.csv` and `states_comparison.csv`, in `states.md`, in methods point 16 and on the 16 state-derived rows of `pulses.csv`. |
| G4 | **Closed** (P2-3) | See 3.3. |
| G5 | **Closed** | See 3.5. |
| G6 | **Closed** | `K_NH = 1 mg N/L fixed (ASM1 default, not fitted)` names `biomass_decay`, `alkalinity_limited` and the combined candidates in `model_fit.csv`; `K_ALK = 0.5 meq/L (ASM3)` and the closed balance name the three `alkalinity_*` candidates; every M6 table carries `assumptions` (§2.1). |
| G7 | **Closed** | `model.md`: "these candidates test an alkalinity limitation, not an inorganic-carbon limitation", and the label on every `alkalinity_*` row says the same. |
| G8 | **Closed** | See 3.4. |
| W6 | **Not closed: one P1** (P1-1; P2-5 to P2-7) | See 3.6. |
| W7 | **Closed** (P2-8, P2-9) | See 3.7. |

### 3.1 G1: registration

- **Hash.** I computed `states.registration_sha256()` from `git archive` extracts of `f1744fc`, `ba394df`, `da49e95` and `0178cd8`: all four give `e34bff2c1693bfade1758db70a897f5d9676e7c34500509e91beaa2b97edf270`, equal to `docs/states_registration.md`. So none of the later commits touched a registered object; the "pre-registered" → "registered" docstring edits in hashed functions (`rolling_theil_sen`, `classify`) are already in `f1744fc`.
- **Coverage of the hash input.** `registration()` serialises the three `RULES`, 22 parameters (windows, thresholds, `N_BOOT`, `SEED`, `MIN_BLOCKS_PER_GROUP`, `CI_ALPHA`, `RHO_MAX`, the null-check constants, `S.EXCLUDE_AFTER_BOUNDARY_D`), `EXCLUDED_BY_FAIL`, `INDICATOR_PARTITIONS` and `inspect.getsource` of the 18 functions in `_REGISTERED_FUNCTIONS`. That matches the document. It does **not** cover:
  - `segments.boundaries`, which `artifact_reasons` calls to place the "NH4 step region" exclusions;
  - `states.excluded_features` (verdicts → excluded features) and `validate_rules`;
  - `state_fractions`, which computes the idle bounds, and the eligibility filter and wiring in `compute`.

  A change in any of these can change a state or an idle bound without changing the hash (P2-1). The document's own definition ("the source code of every decision function (`states._REGISTERED_FUNCTIONS`)") is self-consistent, so this is a scope gap, not a false statement.
- **Mutation.** On a scratch extract of `0178cd8` I changed `sign_class` from `lo > 0` to `lo >= 0`: `test_code_matches_the_registered_rules` failed (hash `1a97184e…`). Reverted: it passes, and the file is byte-identical to `git show 0178cd8:src/bnr_neq/states.py`.
- **Timing.** `f1744fc` (2026-10-26 16:02:10 −04:00) contains `null_response`, `bootstrap_null_coverage`, `null_coverage`, `interval_reading`, `NULL_COVERAGE_MIN_OK = 0.90` and the registration document, and only `states.py`, `tests/test_pass6_fixes.py` and the document. `states_bootstrap_null_coverage.csv` first appears in `ba394df`; no earlier commit contains it or any `null_coverage` column.
- **Label.** `REVISED AFTER FIRST RESULTS (not pre-registered)` is in `states_rules.csv`, `states_comparison.csv`, `states.md` and the report (§1 point 6, §4, methods point 16, the fig4c caption, the table index). It is **not** on other outputs that read the rules: `states_fractions.csv`, `states_timeline.csv`, the 16 state-derived rows of `pulses.csv`, `pulses.md`, Fig. 1 and Fig. 8c (P2-2). Those outputs point to the rules only through the ORP SCENARIO label.
- **History.** The two revisions in the registration document match `progress.md` at `2dd5a2b`: line 711 (pass 2, the "≥ 2 blocks" rule added after "a percentile interval not covering its own estimate with 89 vs 66 bins against a 101-bin block") and line 800 (pass 3, X2, the √(n/n_eff) widening).

### 3.2 G2: null coverage

- **Reproduction.** From the processed inputs (identical to my regenerated ones, §1) I rebuilt the eligible bins: 2,031 bins, pH rising 325, not rising 1,706. The real comparison reproduces exactly (diff 1.05987, interval −0.80981 to 2.50545, block 101). With the project's `bootstrap_null_coverage` and φ = (0.5, 0.9), the cells are **0.85 and 0.87**, equal to the table.
- **Independent check.** My own implementation (AR(1) by `lfilter`, the 101-point least-squares slope by convolution, the 1/e ACF block rule with a minimum of 101, my own moving-block bootstrap with 999 resamples, different seeds) gives **0.865 ± 0.024** at φ = 0.9 for 200 replicates. That agrees with 0.87.
- **Is the null fair?** It keeps what drives the coverage: the real eligible bins, the real run structure of the group labels, the 100-of-101 overlap of neighbouring responses, and the registered block rule, which lands on its 101-bin minimum (median 101) as in the real data. It differs in two ways:
  - It joins non-contiguous eligible bins and uses an LS slope instead of Theil-Sen (both stated in `basis`).
  - Its noise is stationary AR(1) with φ ≤ 0.98, while the real C(t) slopes carry trends and pulse responses.

  Both differences make the null **easier** than the real record, so the under-coverage found (0.805–0.90) is, if anything, optimistic for the real data. Concluding UNDER-COVERS from it is therefore sound. Concluding "adequate" would have needed more caution, but that branch did not occur.
- **Reading.** The state-vs-state interval (1.02 to 5.83) appears in `states_comparison.csv`, §4 and `states.md` with "UNDER-COVERS … an interval excluding 0 is not evidence of a difference", and the row is already labelled DESCRIPTIVE (circular). §1 does not quote it. The report nowhere reads it as a difference.

### 3.3 G4: assumption labels in `pulses.csv`

- **Counts** at `0178cd8`: 248 constant volume; 28 constant volume + K_NH; 16 constant volume + ORP SCENARIO; 12 each for evaporation 0.25, 0.5 and 1 mL/d; 24 "none beyond the metric name". This matches `progress.md`.
- **Per category, from the code path and `nh4_series`:**
  - *None (24 rows):* `t_pulse_d`, `segment_end_d`, `median_pH` (no N, no volume); `o2_windows_computed` (a count); `fa_median_SENSITIVITY_reading_as_TAN_*` (series `electrode_reading`). FA from the electrode reading is a concentration read directly, with no N amount, dose or C_inv, so it is free of the volume frame. Correct.
  - *Constant volume (248 rows):* for example `o2_supply_capacity_*` (kLa (C* − C) V0, series empty), `free_energy_per_fully_converted_pulse_J` (dose/V0), `consumption_uncorrected_*` (a concentration slope read as consumption), and every C_inv row. Correct.
  - *Constant volume + K_NH (28 rows):* exactly the seven `M6_state:biomass_decay` metrics × 4 pulses (`eta_conversion_fraction_MODEL_BASED`, `eta_gap_inversion_misfit_MODEL_BASED`, `eta_uncertainty_interval`, `free_energy_released_MODEL_BASED_J_SCENARIO_T25C`, `free_energy_to_stirrer_energy_ratio`, `mean_consumption_rate_MODEL_BASED_mgN_L_d`, `specific_stirring_energy_J_per_gN`). Correct.
  - *Constant volume + ORP SCENARIO (16 rows):* `classified_fraction` and the three idle-fraction rows × 4. Correct.
  - *Evaporation (36 rows):* consumption, η and removal/dose under each SCENARIO; the constant-volume ASSUMPTION is correctly absent.
- **The 24 rows relabelled after the first run** (`fa_median_mgNH3_L_SCENARIO_T20C/T25C`, `delta_A_kJ_mol_SCENARIO_T20C/T25C`, `delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T20C/T25C`) now carry the constant-volume ASSUMPTION. This is correct: C_inv scales by V0/V(t), and ΔA uses ΔC = dose/V0.
- **Incomplete propagation of K_NH (P2-3).** `model.md` says "the inverted C(t) inherits K_NH through the joint fit", and `nh4_truth_vectors.csv` (the source of C_inv) carries the K_NH ASSUMPTION on all 86 rows. But the 200 `pulses.csv` rows and the states tables computed from `C_inv:c29e26a0fb89` do not name it. Their `nh4_series` points to the vectors table, where the label is, so nothing is presented as measured. Still, the label does not follow C(t) downstream, as the model paragraph implies it should.

### 3.4 G8: provenance

- **Columns.** `model_fit`, `model_calibration_starts`, `model_identifiability`, `model_checks`, `model_volume_scenarios`, `model_volume_scenario_summary`, `model_response_fractions` (per column), `nh4_truth_vectors` and `nh4_truth_w2_consistency` all carry the observation-model and series columns (§2.1). Values: `M6_log`; `M6_state:biomass_decay` or the candidate's own state; the scenario C_inv text in the summary; `C_inv:c29e26a0fb89` in the W2-consistency table.
- **Write sites.** `grep -rn to_csv src/bnr_neq` finds one write site per table. The eight M6 tables (`M6_PROVENANCE_TABLES`) are written only in the loop at `model.py:1566-1568`, where `check_m6_provenance(name, df)` precedes each `to_csv`. `nh4_truth_w2_consistency.csv` is written only at `pulses.py:1078`, from `r.consistency`, which `compute` checked with `check_single_observation_model` at `pulses.py:936`; nothing modifies it in between.
- **Deletion test.** On the committed tables (equal to my regenerated ones), deleting any one of `nh4_observation_model`, `nh4_series` or `assumptions` from each M6 table (and each of the four per-column provenance columns of `model_response_fractions.csv`) makes `check_m6_provenance` raise `ObservationModelError`. Two minor exceptions:
  - Deleting `constant_volume_nh4_series` from the summary table does not raise. It is a supplementary column; the primary `nh4_series` is checked.
  - Deleting `nh4_observation_model` from the W2-consistency table raises `KeyError`, not `ObservationModelError`. It still stops the run.

  I did not rerun the model stage with a dropped column (it takes hours), so this check exercises the function the stage calls, not the stage itself.

### 3.5 G5: the AE range

- **Source obtained.** I downloaded Crawford G. and Sandino J. (2010), *Energy Efficiency in Wastewater Treatment in North America*, WERF OWSO4R07e, CH2M HILL Canada, from the DOE Better Buildings Solution Center copy (`EnergyEfficiency_WWTreatment_BMPs.pdf`, 122 pages, SHA-256 `8b7f10c4…bd710def`). Title page: George Crawford, Julian Sandino, CH2M HILL Canada Limited, 2010, © WERF.
- **Table 2-8 "Aeration Efficiency", p. 2-20** (section 2.4.6.1, "Table 2-8 (WEF, 2009)"; "Source: Based on WEF MOP FD-13 Table 5.1 with additional content"). All 10 rows and their field values match `progress.md` exactly. The footnote "1 – Includes blower power requirements" marks fine-pore (2.5–3.5), coarse bubble (1.0–2.0; spelled "Course" in the source) and submerged turbine, sparger (0.7–1.1). The table gives no field DO or alpha.
- **Rule recomputed.** min = 0.7, max = 3.5, √(0.7 × 3.5) = 1.5652 → 1.6 at 2 significant figures. The config `[0.7, 1.6, 3.5]` follows from Table 2-8 under the registered rule.
- **Order.** `d6d9e5d` (2026-10-27 12:57) changes only `progress.md` and contains the rule, the values and the DO wording. `da49e95` (2026-10-28 01:03) contains the code and config. No later commit deletes a line of `progress.md` (`git diff … | grep -c '^-[^-]'` = 0 for each later commit).
- **DO wording.** Registered: "field AE at the plant's operating DO (not stated in Table 2-8; the cross-check sources state 2 mg/L)". Generated: "at the plant's operating DO" in `basis`; "(> 0; not stated in the source table)" in methods; "field DO not stated" in the source string. The cross-check parenthesis is not in the generated text. That is exactly what `progress.md` discloses.
- **Other references.** `reference/1-s2.0-S1944398624019593-main.pdf` is Azis@Aziz N.F., Ramli N.A. and Hamid M.F.A. (2019), Desalination and Water Treatment 156, 38–45, doi:10.5004/dwt.2019.24134. `reference/1-s2.0-S004313541000638X-main.pdf` is Gresch et al. (2011), Water Research. Both are correctly not used.

### 3.6 W6: Figures 1, 4 and 8

| W6 element | Figure | Finding |
| --- | --- | --- |
| Linkage chain sensors → state → rates → O2 demand → energy → efficiency | Fig. 1 | Present, with C(t) between sensors and state and the W2 bound attached. Numbers: "pH valid" (`sensor_verdicts.csv`), "idle 5-93 %" (whole-record row of `states_fractions.csv`, 0.0537 and 0.9307), "4.57 g/g N" (constant), "AE 0.7-3.5" (config), "C(t) fails" (`nh4_truth_w2_consistency.csv`, best fit inconsistent). All equal their sources. **But the State box, which prints "idle 5-93 %", has the blue border that the figure's own legend defines as "measured".** The idle bounds rest on the ORP SCENARIO threshold (labelled so on the same rows of `pulses.csv` by G3) and on the model-inverted C(t). This presents a SCENARIO- and model-dependent value as measured (P1-1). The W2-bound box is also blue, although its curve depends on the scanned Alk0 and on temperature and kLa_CO2 SCENARIOs (part of P1-1). |
| N_nit(Alk0) bound | Fig. 4a | Bound at 25 and 20 °C, the 1 eq/mol "any proton sink" curve, the C(t) removal with its band, and the dosed N. The 20 °C curve lies under the 25 °C curve and is invisible, and the shaded band has no legend entry (it is in the caption) (P2-6). |
| Rate against FA from W1 | Fig. 4b | −dC/dt with 95 % intervals against median FA from C(t), per W1 slope window, coloured by pulse. The temperature SCENARIO is in the caption only (P2-6). |
| State timeline | Fig. 4c | Clear; pulses dotted; "idle" is the inactive-over-supplied state. |
| η per pulse | Fig. 8a | Bars, band and "n.c." at P3, matching `pulses.csv`. The 0.25 mL/d evaporation SCENARIO is in `pulses.csv` but not drawn (only 0.5 and 1) (P2-5). |
| O2 utilization under the kLa SCENARIOs | Fig. 8b | Hatched bars at P1 and P2 for kLa 1/d, where `pulses.csv` says "not computed: SCENARIO infeasible". They are drawn at a fixed height of 5.0 (`figures.py:717`), a placeholder on a data axis that a reader could take for a value (P2-7). "n.c." at P3 is correct. |
| Idle fraction | Fig. 8c | Bounds per pulse equal `pulses.csv`. No conditional marker, because no interval reaches coverage 0.5. The marker and its legend entry correctly appear only when a value exists. |
| Energy under the AE and kLa SCENARIOs | Fig. 8d | Bars at AE 1.6 with AE 0.7–3.5 whiskers, the floor at 4.57/1.6. "n.c." at P3 only; the missing kLa 1/d bars at P1 and P2 carry no "n.c.", although `progress.md` says "d writes n.c. at pulses where no series was computed", which is literally true but not useful there (P2-7). |

- **Readability.** At 25 pt with Arial drawn at weight 700, every label is readable. The fixed subplot box (0.23/0.95/0.18/0.92) leaves a wide empty left margin on Fig. 1, which has no axis labels, but nothing is cut off.
- **Plotting rules.** All 29 figures have one axes, (10, 6), four spines, no title, no grid and no layout engine. The PNGs are 3000 × 1800 at 300 dpi, RGBA with a transparent corner, each with an SVG and a `no_text/` copy (`tests/test_figures.py`, which passed in my run).

### 3.7 W7: report order, labels and precision

- **Order.** §1 What the data say about the process → §2 Energy and efficiency → §3 Sensor validity → §4 Process state and constraints → §5 Reduced model → §6 Status → §7 Facts and assumptions → S1 Stochastic-thermodynamic estimators → S2 Detailed tables.
- **The possible conflict.** W7 lists "a one-page summary …, then energy and efficiency, then sensor validity, with stochastic-thermodynamic results moved to a supplementary section". The final order meets that sequence literally: summary (§1), energy (§2), sensor validity (§3). W7 does not place process state or the model; putting them after sensor validity follows the dependency (states and C(t) use the verdicts), and both are summarised in §1. The move "§5 → §3" in `progress.md` is the move that **made** sensor validity come right after energy; before it, sensor validity came after the model. **Ruling: the order meets W7.**
- **Summary length.** §1 is 1,593 words (points 1–9 plus Fig. 1), about three pages, not one. Points 5, 8 and 9 alone are about 900 words (P2-8).
- **CONSISTENCY CHECK labels.** The coexistence affinity, ΔA and σ ≥ 0 are labelled "CONSISTENCY CHECK (not a finding)" in the §6 status table, the S2 sentence above the ΔA table and the fig7b caption, and σ ≥ 0 sits under "Consistency checks of the selected model (… not findings)" in §5. Two mentions carry no such label (P2-9):
  - ΔA appears as a descriptor in the EXPLORATORY Spearman table of S2 (`consumption_weighted` vs `delta_A_kJ_mol_SCENARIO_T25C`, ρ = 1, n = 4).
  - ΔA and the coexistence affinity appear in the reference-reproduction table (as reproduction checks).

  Neither is presented as a finding.
- **Rounding.** A scan of every number in the report for 5 or more significant figures finds only:
  - the configured constants 0.77649, −2.89855, 0.52087 and 298.15 K;
  - integer row counts (414,757 …) and optimiser evaluation counts (25,249; 26,955);
  - window boundaries in days written to 0.001 d (10.983, 17.875 …), which are time stamps, not estimates.

  No estimated or derived value has more than 4 significant figures. The previously unrounded balance values (for example `13.009538747950954`) now print at 3 significant figures.

## 4. Scientific notes

- **G2 changes nothing and closes the question honestly.** Both pH intervals already included 0. The finding is that they are too narrow (coverage 0.84–0.90 at nominal 0.95), so the W3 association question has no answer from these data. Because the null is easier than the real record (§3.2), the true under-coverage is probably worse. The registered response, keeping the block rule and labelling the intervals, is the right one. Widening the block after the fact would have been a new, unregistered rule.
- **The G5 range is a full-scale equipment range, correctly kept as a SCENARIO.** The central value moved from 1 to 1.6, so every supply energy fell by 1.6 and every Gibbs-to-aeration ratio rose by 1.6. No conclusion depends on the AE value: the ordering across pulses and kLa SCENARIOs is unchanged, and the demand floor is stated as a range. Restricting the rows to blower-inclusive diffused air matches what the energy metric converts (blower energy of a BNR plant). It excludes the mechanical aerators, whose rows are not marked as including blower power. Their field values (0.5–1.4) lie mostly inside the chosen range and extend it downward only to 0.5–0.6, so including them would have moved the low end, not the conclusions.
- **K_NH inheritance.** `model.md`'s statement that C_inv inherits K_NH is correct: the log observation model is fitted jointly with the selected rate law. It implies that every C_inv-derived number rests on K_NH = 1 mg N/L too. At C ≫ K_NH (C(t) ≈ 180–290 mg N/L), C/(K_NH + C) ≈ 0.995–0.997, so the sensitivity is negligible. The label gap (P2-3) is formal, not numerical.
- **Fig. 1 is the only place where the basis of a number is drawn wrongly.** The State box's number is the most SCENARIO-dependent quantity in the chain: the idle bounds run from 5 to 93 %, a width that comes entirely from the unclassified bins. A blue "measured" border is the opposite of what §1 point 8 says ("bounded, not measured").

## 5. Tests

- **New tests.** `tests/test_pass6_fixes.py` has 21 test functions (31 cases with parametrization); `tests/test_figures.py` has 2.

  | Type | Tests | Adequate? |
  | --- | --- | --- |
  | Positive | Hash equals the document; the status text; the null-response structure; coverage counting; `candidate_assumptions`; `m6_tables` keeps values; `row_assumptions` per kind; labels on the committed `pulses.csv`; the G5 rule recomputed from literal Table 2-8 rows; the figure rules. | Yes |
  | Negative | Parameter and rule mutations change the hash; small groups give no replicate; a missing or wrong provenance column raises; a W2 table with two series raises. | Mostly. `progress.md` says "parametrized mutations of a rule, a threshold, **a seed or a decision function** change the hash", but the parameter list is `SLOPE_HALF_BINS`, `ORP_LEVEL_MV`, `N_BOOT`, `MIN_BLOCKS_PER_GROUP` and `NULL_COVERAGE_MIN_OK`: no seed and no decision-function mutation. My manual mutation shows the hash does catch function changes, so the gap is in the test, not the code (P2-10). |
  | Boundary | `interval_reading` at exactly 0.90 and at `nextafter` below. | Yes |
  | Regression | Extended column-equality tests against the pass-3 and `8e10d0e` fixtures. | Yes, see below |
  | Integration | Real M6 tables carry provenance; the W2 table names one series; every pulse row carries assumptions; states registered and read against the null. | Yes |
  | End-to-end | `test_main.test_run_all_is_deterministic` (two fresh run-alls) plus the figure tests on redrawn files. | Yes |

- **The strengthened G4 test** (every row with `nh4_series` starting `C_inv:` or `M6_state:` names its volume frame, independent of the regex, with only `o2_windows_computed` exempt by name) would have caught the 24-row error. It is adequate for the volume frame. It does not check K_NH on C_inv rows (P2-3).
- **Changed existing tests,** each compared with `git show 2dd5a2b:<file>`:
  - `test_pass4_fixes.py`: two section slices renumbered (§4 → §5 for the model, §3 → §4 for states). `test_real_selection_fit_and_model_based_numbers_equal_pass3` now asserts that the extra columns are exactly the three pass-6 columns (none for `nh4_truth_summary.csv`) and compares `got[ref.columns]` with the fixture by `assert_frame_equal`. Every fixture column is still compared exactly, so any change to a pre-existing column still fails. Column **position** is no longer asserted, only the set of extras. **Not weakened** in value terms.
  - `test_pass5_fixes.py`: `PASS6_COLS` appended to the expected new columns; the value comparison of pre-existing columns is unchanged; section slices renumbered. **Not weakened.**
  - `test_report.py`: `len(R.FIGURES) == 8` → `29`, plus the group order and the presence of fig1/fig4/fig8 groups. **Not weakened** (it follows the split).
  - `test_segments.py` and `test_validation.py`: the expected file list follows the new stems, and `test_validation.py` adds `len(figs) == 12`. **Not weakened.**
  - `test_pass6_fixes.py` itself changed in `0178cd8` by one section slice.
- **`test_figures.py`.** The rules on spines, titles, grid, layout, size and rcParams are checked on the figure **objects** captured at `_save`. Size, dpi, mode, transparency and the SVG's existence are checked on the **saved files**, for both the normal and `no_text/` copies. It does not check that `no_text/` files actually contain no text, or anything inside the SVGs, and it checks files redrawn by `make_all` from a copy of the tables, not the files run-all wrote (those are covered by the determinism test). Adequate for W6's plotting rules; noted under P2-10.
- **No weakened test found.**

## 6. Rulings on the implementer's judgment calls

| # | Judgment call | Ruling | Reason |
| --- | --- | --- | --- |
| 1 | G1: freeze after the fact and label "revised after first results" | **Accept** (P2-1, P2-2) | It is the only honest option: a pre-registration cannot be created retroactively. The hash is reproducible, matches at every later commit, and catches a function edit. The label is used consistently wherever the rules are **described**, but not on every output that **uses** them, and the hash scope omits `segments.boundaries`, `excluded_features` and `state_fractions`. Both are P2. |
| 2 | G2 design and threshold (φ ∈ {0.5, 0.9, 0.95, 0.98}, 200 replicates, 999 resamples, ≥ 0.90 at every φ; block rule unchanged; two stage previews before the validation run) | **Accept** | Code and rule are in `f1744fc` before any null output exists. The previews ran on the committed code, and only report rendering changed between them (no `states.py` change: hash identical). 0.90 is a lenient threshold against nominal 0.95, and the result fails it anyway, so the threshold was not chosen at a margin. 200 replicates give SE ≈ 0.02, enough to separate 0.85 from 0.95. |
| 3 | G4 correction after the first full run | **Accept** | The error was found by the implementer's own check, disclosed, and fixed in code and test together. The strengthened test checks the property directly (rows from C_inv or the model state), not through the regex that encoded the error. |
| 4 | G6/G7: no rename, no CO2 term | **Accept** | Renaming would break fixtures and earlier reports for no gain in content. A CO2 term is a new model and a new fit, outside a labelling pass. The scope is stated in the labels and in `model.md`. |
| 5 | G5 source rule; shorter DO wording; Stenstrom caveat dropped; Baeten unused | **Accept** | Verified against the PDF (§3.5). The geometric mean is the natural centre for a quantity that enters as 1/AE. The shorter wording drops a parenthesis, not a fact (the table states no DO), and is disclosed. Stenstrom is no longer the source and is kept only as a cross-check in `sae_source`, where "not for design" still appears. Leaving an unverifiable number unused is correct. |
| 6 | W6: 29 single-axes figures with `no_text/` copies; Arial at 700 | **Accept** | It follows the plotting rules as written. matplotlib's nearest-weight substitution is disclosed. The cost is that W6's "Fig. 4" and "Fig. 8" are now groups of single-axes files, which the report presents under one caption each. That is acceptable. |
| 7 | No pytest rerun after the ruff rename; no code hashes for the W6/W7 run | **Accept for this review; the omission matters beyond it** | My run on `0178cd8` (§1) is the first full test of the committed code, and it passed. The hashes recorded before run-all and after pytest are identical. In general, a change after the last full test, however trivial, must be followed by at least the affected test file. "Code hashes not recorded" removes the evidence that the validated code is the committed code. Future passes should record both. |
| 8 | P2-1 and P2-2 of pass 5 | **Both still open, untouched** | `model_fit.csv` still shows `cal_at_bounds` = C0 for `biomass_decay` (no bound change). `bounds_text` still says "natural-unit kinetic parameters are not estimates" for a C0-only flag, and §1 point 9 still says the N amount "is therefore unchanged" without "within the tolerance". Neither was in scope. |

## Remaining issues by severity

| Severity | Issue | Fix |
| --- | --- | --- |
| **P1-1 (fixed in `47dfd18`, verified)** | Fig. 1 draws the State box ("idle 5-93 %") and the W2-bound box with the blue border that its legend defines as "measured". The idle bounds rest on the ORP SCENARIO threshold and on the model-inverted C(t); the W2 curve rests on the scanned Alk0 and on temperature and kLa_CO2 SCENARIOs. This is a SCENARIO-dependent value presented as measured, a failure condition of `docs/remediation_prompt_2026-10-19.md`. | In `figures.fig1` (`figures.py:240-246`): State → `ORANGE` (scenario; it also depends on the model C(t)); W2 bound → `ORANGE`, or a legend entry saying "measured input, SCENARIO axis". Add a test that every box whose quantity carries a SCENARIO or ASSUMPTION label in `pulses.csv`/`states_*.csv` is not drawn as "measured". Redraw with `python -m bnr_neq figures`; no table changes. |
| P2-1 | The registration hash does not cover `segments.boundaries` (used by `artifact_reasons`), `states.excluded_features`, `validate_rules`, `state_fractions` or the eligibility filter in `compute`. | Add them to `_REGISTERED_FUNCTIONS` (and `segments.boundaries` by module) in a separate commit that updates `docs/states_registration.md`. |
| P2-2 | `REVISED AFTER FIRST RESULTS (not pre-registered)` is missing on `states_fractions.csv`, `states_timeline.csv`, the 16 state-derived `pulses.csv` rows, `pulses.md`, Fig. 1 and Fig. 8c. | Add `registration_status` (or include it in `STATE_ASSUMPTIONS`). |
| P2-3 | K_NH is named on `nh4_truth_vectors.csv` and in `model.md` as inherited by C_inv, but not on the 200 `pulses.csv` rows or the states tables computed from C_inv. Numerically negligible (C ≫ K_NH). | Append the K_NH ASSUMPTION to C_inv-derived rows, or state in `model.md` that it is named only at the source table. |
| P2-4 | In `states_fractions.csv` and `states_comparison.csv`, the new columns are inserted before the last two pre-existing columns; `progress.md` says "column order kept". | Append at the end, or correct the log. |
| P2-5 | Fig. 8a omits the 0.25 mL/d evaporation SCENARIO, which `pulses.csv` carries. | Draw it or name the omission in the caption. |
| P2-6 | Fig. 4a: the 20 °C curve is hidden under the 25 °C curve, and the band has no legend entry. Fig. 4b: the FA temperature SCENARIO is not on the axes. | Offset or annotate; put "SCENARIO 25 °C" in the axis label. |
| P2-7 | Fig. 8b draws infeasible kLa bars at a fixed placeholder height of 5.0 on a data axis. Fig. 8d leaves those pulses empty without "n.c.". | Draw the hatch to the axis top or as a marker; add "n.c." per missing series. |
| P2-8 | §1 is 1,593 words, not the one page W7 asks for. | Shorten points 5, 8 and 9 to their conclusions and move the detail to §2, §4 and §5. |
| P2-9 | ΔA appears without a CONSISTENCY CHECK label in the S2 EXPLORATORY Spearman table and in the reference-reproduction table. | Add the label to the descriptor name or a footnote. |
| P2-10 | Tests: no seed or decision-function mutation in the hash tests (contrary to `progress.md`); `test_figures.py` does not check that `no_text/` files have no text; the W2-consistency deletion raises `KeyError`, not `ObservationModelError`. | Add the two parametrized cases; check that `no_text` SVGs have no `<text>`/glyph elements; check the column before indexing. |
| P2-11 | `check_fig1_bases` skips a box when its source table has no `assumptions` column (`figures.py`, `if "assumptions" not in tab.columns: continue`), so a measured State box passes the production guard if that column is lost (reviewer mutation 3). The tests catch it (`_labelled` must be true for State and W2), and G8 makes the column hard to lose, but the guard itself falls back silently. | Raise when a source table lacks `assumptions`, with a test. |
| P2-12 | The W2-bound → C(t) arrow is still blue, the legend's "measured" colour, although it leaves a scenario box. Arrows are outside the legend, as `progress.md` says, but the colour invites the reading. | Draw it in `INK2` like the other arrows, or add an arrow entry to the caption. |
| P2 (pass 5) | P2-1 (C0 cap of the selected calibration) and P2-2 (two wording nits) of `docs/rereview_2026-10-25_pass5.md`: still open, unchanged. | As in pass 5. |

## Readiness to merge

Ready to merge at `47dfd18` once the working tree is cleaned:

1. ~~P1-1 is fixed~~: done in `47dfd18` and verified (addendum).
2. Delete the untracked uv stub `src/bnr_non_equilibrium_2026/`. For `reference/`, commit the summary note and a WERF note with the PDF's SHA-256 (`8b7f10c4…bd710def`), not the publisher PDFs.

The P2 items can go into a follow-up pass. P2-1 and P2-2 (registration-hash coverage and the missing "revised after first results" label) should come first, since they bear on how far the W3 state results can be trusted.

## What is sound and should be kept

- A frozen, hashed rule set with an honest status line, instead of a retroactive "pre-registered" claim, and a test that fails on any edit to a hashed function.
- A null check committed before it was run, with a decision rule that could only make the reading more cautious, and a conservative conclusion that does not depend on the null's leniency.
- Raw-text column guards and an AST-level "plotting moved, nothing computed" change, which made the W6/W7 refactor easy to verify.
- A primary source for the AE range, a registered selection rule written down before the edit, and a test that recomputes the rule from literal table rows.
- Figures redrawn from the saved tables, so every plotted number traces to a table.
- Complete disclosure of the 24-row labelling error, the two previews, the missing post-rename test run and the font substitution.

## Out of scope

- **`docs/reports/BNR_non_equilibrium_full_report_2026-10-28.{md,html}` (untracked).** Spot-checked against `analysis_report.md` and the tables: the selection skill 0.958/0.748/0.863, the demand floor 1.31–6.53 kWh/kg N, the supply energies 4.2/13.7/62.3 kWh/kg N, the idle bounds 5–93 %, and the state interval +1.02 to +5.83 read as "not evidence of a difference" all match. Two sentences of its "Bottom line" go beyond the analysis:
  - "conversion slowed with time": η rises from P1 0.38 to P2 0.81 before falling.
  - "because aeration capacity was spent on biomass that was idle or limited": the idle fraction is only bounded (5–93 %), so this is not established.

  Fix those before circulating it.
- **Untracked files.**
  - `src/bnr_non_equilibrium_2026/` is the `uv init` stub (`print("Hello …")`) and should be deleted, not committed.
  - `reference/` should not be committed as is: the two journal PDFs are publisher copies, and the WERF source is not in it. For reproducibility, commit `reference/wastewater_aeration_efficiency_literature.md` together with a short note giving the WERF report's citation, where a copy can be obtained (the DOE Better Buildings copy), its SHA-256 (`8b7f10c48995dfaf9a5dd91650cd43519780ff71ebadc317dc24eb15bd710def`) and a transcription of Table 2-8 with its footnote. The G5 test already holds the three rows used.

## Addendum: verification of the P1-1 fix, `47dfd18` (2026-11-06)

`47dfd18` changes `src/bnr_neq/figures.py`, `tests/test_figures.py`, `docs/progress.md` and the four Fig. 1 files (PNG and SVG, with their `no_text/` copies), nothing else (`git diff --stat 73f4079 47dfd18`; `73f4079` only adds the re-review prompt to `0178cd8`).

| Check | Result |
| --- | --- |
| Code read | Each box now carries a named basis (`measured`/`model`/`scenario`) instead of a colour; State and W2 bound are `scenario`; `FIG1_BASIS_COLORS` drives both the borders and the legend, so they cannot disagree. `check_fig1_bases` runs before drawing and raises if a box is `measured` while its source rows (`FIG1_BOX_SOURCES`) carry SCENARIO or ASSUMPTION, or if a basis is unknown. No computation changed. |
| `uv run ruff check` (clean worktree at `47dfd18`) | `All checks passed!` |
| `python -m bnr_neq figures` on a copy of the reviewer's regenerated `0178cd8` outputs | Exit 0, 29 figures. `diff -rq` against the committed `47dfd18` `outputs/`: no differences. Against `0178cd8`: only the four Fig. 1 files. The SVG diff is two lines per file, the stroke of the State and W2-bound rectangles, `#2a78d6` (blue) → `#eb6834` (orange). The redrawn PNG shows Sensors blue, C(t)/Rates/O2 violet, State/W2/Energy/Efficiency orange. |
| `tests/test_figures.py` (`run_all_dir` pointed at the reviewer's regenerated tables, byte-identical to the committed ones, to avoid a third full run-all) | 7 passed (2 existing, 5 new). |
| Mutation 1: State basis back to `measured` | 2 failed, 2 errors (the guard raises inside `fig1`); restored. |
| Mutation 2: borders drawn blue regardless of basis | 1 failed (`test_drawn_fig1_borders_and_legend_follow_the_basis`); restored. |
| Mutation 3: `assumptions` column dropped from `states_fractions.csv`, State set to `measured` | The production guard passes silently (P2-11). The test suite would still fail, because it asserts the State and W2 rows are labelled. |
| Hashes | Identical before and after the checks (`diff` exit 0); the worktree is clean. `figures.py` `d9e91ee3…f8fa8` and `test_figures.py` `8c74332d…139ec` match the hashes `progress.md` records for the implementer's 496-test run. |

No full suite was re-run by the reviewer: the only code change after the reviewer's 491-test run of `0178cd8` is in `figures.py` and `test_figures.py`, both exercised above, and the implementer's full run (496 passed, no skips) is on the same hashes.

P1-1 is closed. Two P2 items are new (P2-11, P2-12).
