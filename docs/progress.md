# Progress log

Moved on 2026-10-18 from `outputs/reports/progress.md` to `docs/progress.md`, so that `rm -rf outputs` no longer deletes it and the generated tree contains only generated files. Earlier entries mention parking the file in the scratchpad during gate runs; that is no longer needed.

Implementer log. Every gate result below comes from commands executed in the session that wrote it. This log does not approve anything; the independent reviewer decides.

## M0. Environment and inspection (2026-10-18)

### Inspection findings

Confirmed facts:

- Tree before changes: `.gitignore`, `.python-version`, `README.md`, `pyproject.toml`, `implementation_prompt.md`, `review_prompt.md`, `Non-Equilibrium Analysis Plan for a BNR System.md`, `data/{Cond,NH4,ORP,pH}.csv`, `src/bnr_non_equilibrium_2026/__init__.py`. No `uv.lock`, no `.venv`, no `main.py`, no tests.
- `.python-version` = `3.14`; `requires-python = ">=3.14"`. Python is newer than 3.11, so `tomllib` is in the standard library and `tomli` is not needed.
- `pyproject.toml` **already has a `[build-system]` table** (`uv_build>=0.12.19,<0.13.0`), and a script entry `bnr-non-equilibrium-2026 = "bnr_non_equilibrium_2026:main"` pointing to the existing placeholder package `src/bnr_non_equilibrium_2026/__init__.py` (a "Hello" `main()`). There is no `main.py`.
- `README.md` is a two-line title/description.
- All four CSV files use CRLF line endings; `NH4.csv` starts with a UTF-8 BOM (`EF BB BF`); the other three have no BOM.
- **Header discrepancy with Section 4.** The actual first lines are:

  | File | Header in Section 4 | Actual header |
  | --- | --- | --- |
  | `Cond.csv` | `"Time, d",Cond. A` | `Time(d),Conductivity(A)` |
  | `NH4.csv` | `Time(d),Concentration` | `Time(d),Concentration(mg/L)` (after BOM) |
  | `ORP.csv` | `Time,ORP (mV)` | `Time(d),ORP(mV)` |
  | `pH.csv` | `Time(d),pH` | `Time(d),pH` (matches) |

- Last row of `Cond.csv` is `49.69529086,` (empty value). Last rows of NH4/ORP repeat a rounded stamp (`47.92`, `53.23`); pH ends at `48.00439815`.

Assumptions (to be verified in M1):

- The renamed headers carry the same quantities and units as Section 4 (time in days; current in A; NH4 in mg NH4+/L as confirmed by the experimenter, even though the file only says `mg/L`; ORP in mV vs bare Ag).

### Environment

- `uv 0.12.19`; interpreter CPython 3.14.7 (`/opt/homebrew/opt/python@3.14`). A pyenv shim warning is printed by `uv python list` and is harmless (uv does not use pyenv).
- Build system: kept the existing `uv_build` backend (the prompt's hatchling block applies only when `[build-system]` is absent). Added `[tool.uv.build-backend] module-name = ["bnr_non_equilibrium_2026", "bnr_neq"]` so both the existing placeholder package and the new `bnr_neq` package are installed in editable mode. The existing package and its script entry are left untouched.

### Commands run

```
uv sync                                   # created .venv and uv.lock (CPython 3.14.7)
uv run python --version                   # Python 3.14.7
# created src/bnr_neq/__init__.py; appended [tool.uv.build-backend] to pyproject.toml
uv sync                                   # exit 0
uv run python -c "import bnr_neq; print(bnr_neq.__file__)"
uv add numpy pandas scipy matplotlib statsmodels ruptures pyarrow
uv add --dev pytest ruff
uv sync                                   # Resolved 30 packages, exit 0
uv run python -c "import bnr_neq; print(bnr_neq.__file__)"   # exit 0
```

Resolved versions: numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2, statsmodels 0.15.0, ruptures 1.1.10, pyarrow 26.0.0; dev: pytest 9.1.1, ruff 0.17.1.

### Gate M0: PASS

- `uv sync`: exit 0.
- Both `uv add` commands: succeeded (packages listed above installed).
- Import check: `.../src/bnr_neq/__init__.py`.

## Proposed approach for M1 and M2 (written before implementation)

### Pre-implementation verification of Section 4 (exploratory scratch script, raw pandas/NumPy)

| Property | Section 4 | Found | Status |
| --- | --- | --- | --- |
| Headers | see M0 table | 3 of 4 differ in spelling | **Discrepancy** (same quantities; map actual headers, still raise if absent) |
| Cond rows / dropped | 71,759 of 71,760; last empty | 71,760 rows; 1 non-numeric (index 71759, `49.69529086,`) | Reproduced |
| Cond span | 0.000693 to 49.6946 | 0.000692521 to 49.69459834 | Reproduced |
| Cond sampling | 60 s uniform | 59.834 to 59.835 s (= 1/1444 d), uniform | Minor: nominal 60 s is actually 59.83 s |
| Cond range | −0.366 to −0.05946 A; 0.6932 to 1.5817 mS/cm | same; 0.693218 to 1.581739 | Reproduced |
| NH4 rows, span, range | 414,785; 0 to 47.92; 9.52 to 75.13 | same | Reproduced |
| NH4 stamps | 0.01 d; 22 to 450 per stamp | 4,793 stamps, spacing 0.01 only; 22 to 450 (interior 42 to 450) | Reproduced |
| NH4 sampling | 20 s, ~3 s after 41.05, ~2 s from 45 | median 43 / 290 / 435 rows per stamp = 20.1 / 2.98 / 1.99 s | Reproduced |
| ORP rows, span, range | 459,884; 0 to 53.23; 77.96 to 318.7 | same | Reproduced |
| ORP stamps | 21 to 87 per stamp | 5,324 stamps; 21 to 87 (interior 86 to 87) | Reproduced |
| pH rows, span, range | 414,757; 0 to 48.0044; 6.47 to 8.02 | same | Reproduced |
| pH sampling | 10 s, two 20 s near 4.398 | two intervals of 20.0000 s at 4.39815 and 4.39838 | Reproduced |
| pH tie fraction | 88% | 0.884 | Reproduced |

Observations that are not in Section 4 and that I will carry as open issues, not act on:

1. **Stamp convention.** The first NH4 stamp holds 22 rows and the first ORP stamp 43, about half of the interior 43 and 86. This suggests the stamps are rounded to the nearest 0.01 d rather than truncated. If so, the prescribed `t = stamp + 0.01 j/n` places every reconstructed time 0.005 d (7.2 min) late, and spreads the first stamp's rows over twice their true span. I will implement the prescribed formula (the prompt governs), state the assumption in the outputs, and record this evidence. 7.2 min is below the resolution of every M2 window.
2. **Possible pH clock offset (plan item D7).** At the multi-probe event near day 2.88, a rolling-median step locator gives day 2.8125 for pH, 2.9095 for ORP (reconstructed) and 2.919 for Cond. A pH clock about 0.1 d early would bias the ORP/pH pulse windows. I will add a small diagnostic table of this event's timing per sensor in M1 and report it, without shifting any clock.

### M1 plan

- `config/experiment.toml`: exactly the Section 5 content. `config.py` loads it with `tomllib`, checks required sections and keys and their types, and exposes `require(cfg, dotted_key)` (raises `MissingInputError` naming the key) and `not_computed(key)` (returns `"not computed: <key> missing"`).
- `constants.py`: R, F, T25, Kohlrausch λ values, AgCl/Ag E°, all with units and sources (only what M1/M2 use, extended in M4).
- `io.py`: one explicit header map per file, using the **actual** headers (deviation D1); the header line is compared exactly and an error naming the file is raised if it is absent. Read with `encoding="utf-8-sig"`, `dtype=str`, `keep_default_na=False`, then `pd.to_numeric(errors="coerce")`; dropped rows are returned with file name and 1-based line number and written to `outputs/tables/dropped_rows.csv`. EC conversion raises outside 0 to 5 mS/cm. `nh4_mgN_L` is added alongside `nh4_mgNH4_L`. Non-decreasing time is checked, else raise.
- `timebase.py`: `check_stamp_spacing` (distinct stamps must differ by exactly 0.01 d, tolerance 1e-9 d, else raise); `reconstruct_stamps` (`t_stamp_d` is kept, `t_d` = stamp + 0.01 j/n; strictly increasing, else raise); `gap_inventory`, where the native interval is a **local** rolling median of 101 intervals, because the NH4 rate changes from 20 s to 2 s (a global median of 3 s would flag every 20 s interval as a gap). Intervals > 2x local median (relative tolerance 1e-4 for 1e-9 d timestamp rounding) are classed `gap`; intervals > 1.5x are classed `missing_sample`, so the two pH 20 s intervals are listed. `common_grid`: 4,792 bins of 0.01 d on [0, 47.92] (left-closed, last bin closed), time-weighted means with weight = forward interval, capped at the local median for an interval classed as a gap; per-sensor sample count and empty mask; no interpolation at all, so empty bins stay NaN.
- Outputs: `outputs/processed/native_{cond,nh4,orp,ph}.parquet`, `outputs/processed/grid_0p01d.parquet`, `outputs/tables/{data_properties,gaps,dropped_rows,clock_check}.csv`, `outputs/reports/ingest.md` (properties vs Section 4, assumptions, dropped rows).
- CLI `python -m bnr_neq` with all nine commands; `ingest` and `validate` are implemented, and the others exit non-zero with "not implemented: milestone Mx".
- Tests (production functions only): row counts, spans and dropped count from the real files; header mismatch raises with the file name (on a temp copy); EC −0.1 A gives 0.810725, a mA value raises; stamp-size ranges; strictly increasing reconstruction; non-monotonic input raises; grid bins with no samples are NaN with count 0, and each non-empty bin mean lies within that bin's sample min and max (no fabricated values); config missing key yields `not computed`.

### M2 plan

- `thermo.py` (M2 subset of Section 6): `pulse_dose` (from config: mmol/L, mg NH4+/L, mg N/L, mg Cl/L), `kohlrausch_delta_ec`, `nernst_slope_V`, `chloride_orp_rise_mV`. Each with a unit test against Section 6.
- `validation.py`: pulse table with the exact windows; NH4 predicted 84.31 mg NH4+/L and response fraction = observed/predicted; EC predicted 0.700 (immediate) and 0.274 (after complete nitrification), both labeled limiting (infinite-dilution) values; ORP predicted rise per pulse per chloride scenario (1, 3, 8, 13.5 mM) at each scenario temperature (20, 25 °C), using the cumulative chloride before each pulse (c0 + (k−1)·Δc), labeled `scenario`; observed pH change. Cumulative balances: EC staircase (EC0 = median day 0.3 to 1.3, +0.700 or +0.274 per pulse), EC medians day 0.3 to 1.3 and 48.5 to 49.6; N balance (start/end medians, 4 × 65.46 added, removed = start + added − end, implied rate = removed / (47.4 − 0.8 d)).
- `qc.py`:
  - NH4 jump detector: on the reconstructed native series, step(t) = median over (t, t+w] − median over [t−w, t), w = 0.05 d (72 min, long compared with the 14.4-min stamp, short compared with consumption: 3 mg/L/d gives 0.15 mg/L per window), threshold |step| ≥ 4 mg NH4+/L (more than 4x the largest sub-segment HF residual, 0.92 mg/L) on a 0.001 d evaluation grid; local maxima merged within 0.1 d. Classification: a positive jump within 0.5 d after a pulse is a `pulse_response` (lag reported); everything else is an `electrode_offset`. 10 to 90% transition time from a 2-min rolling median between the pre and post levels, with the stamp-resolution caveat.
  - Held/damped output: HF residual = x − 30-min centred rolling median; σ_HF = 1.4826 x MAD per 0.1 d window; flagged when σ_HF < 0.01 mg NH4+/L (the reporting resolution); runs ≥ 0.5 d merged into segments; σ_HF also tabulated per reference sub-segment.
  - Late EC pulses: residual from a 6 h centred rolling median of EC; peaks with height and prominence ≥ h, separation ≥ 30 min, after day 41.9, with h = 99.9th percentile of the residual before day 41.9 (noise-floor rule fixed before looking at the count); exponential decay fitted to 2 h after each peak (`scipy.optimize.curve_fit`); median τ reported; a sensitivity table of count vs h is reported, so the threshold dependence is visible and not tuned.
- Outputs: `outputs/tables/sensor_validation.csv` (one row per pulse, units in column names), `outputs/tables/{sensor_balances,nh4_jumps,nh4_held_segments,ec_late_pulses,ec_late_pulse_sensitivity}.csv`, `outputs/figures/fig3_sensor_validation.{png,svg}`, `outputs/reports/sensor_validation.md` with a pass/fail verdict per sensor and the reference comparison.

## M1. Ingestion, time base and quality flags (2026-10-18)

### Implemented

- `config/experiment.toml` (Section 5 content verbatim), `src/bnr_neq/{config,constants,io,timebase,pipeline,report,__main__}.py`, `tests/{conftest,test_config,test_io,test_timebase}.py`.
- `pipeline.py` is an extra module (not in Section 7) that holds the stage orchestration, so the CLI and the tests call the same functions.
- Outputs: `outputs/processed/native_{cond,nh4,orp,ph}.parquet`, `outputs/processed/grid_0p01d.parquet`, `outputs/tables/{data_properties,dropped_rows,gaps,grid_coverage}.csv`, `outputs/reports/ingest.md`.

### Commands run

```
uv run python -m bnr_neq ingest            # first run; found per-interval artifact rows in the gap list
# merged contiguous flagged intervals into runs and labelled partial-stamp artifacts
uv run pytest -q                           # 51 passed in 1.30s
uv run ruff check .                        # 7 findings with ruff 0.17 defaults (RUF007, ISC004 x4, RUF046, I001); fixed in code
rm -rf outputs/processed outputs/tables
uv sync                                    # Resolved 30 packages, Checked 28 packages
uv run ruff check .                        # All checks passed!
uv run pytest -q                           # 51 passed in 1.30s
uv run python -m bnr_neq ingest            # ingest: 1 row(s) dropped
uv run python -m bnr_neq segments          # "not implemented: segments belongs to milestone M3", exit 2
```

### Section 4 properties (from `outputs/tables/data_properties.csv`)

| Property | Section 4 | Obtained | Result |
| --- | --- | --- | --- |
| Cond rows | 71,759 of 71,760 | 71,759 of 71,760; dropped line 71761 (`49.69529086,`) | Reproduced |
| Cond span | 0.000693 to 49.6946 | 0.000692521 to 49.694598 | Reproduced |
| Cond sampling | 60 s uniform | 59.8337 to 59.8346 s | Discrepancy (nominal): actual interval is 59.83 s (1/1444 d); uniform, no gaps |
| Cond range | −0.366 to −0.05946 A; 0.6932 to 1.5817 mS/cm | same; 0.693218 to 1.58174 | Reproduced |
| NH4 rows, span, range | 414,785; 0 to 47.92; 9.52 to 75.13 | same | Reproduced |
| NH4 stamps | 0.01 d; 22 to 450 rows | 4,793 stamps, all 0.01 d apart; 22 to 450 | Reproduced |
| NH4 sampling | 20 s, ~3 s after 41.05, ~2 s from 45 | 20.09 s, 2.98 s, 2.00 s (median rows per stamp 43, 290, 432) | Reproduced |
| ORP rows, span, range | 459,884; 0 to 53.23; 77.96 to 318.7 | same | Reproduced |
| ORP stamps | 21 to 87 rows | 5,324 stamps; 21 to 87 | Reproduced |
| pH rows, span, range | 414,757; 0 to 48.0044; 6.47 to 8.02 | same (t_max 48.004398) | Reproduced |
| pH sampling | 10 s, two 20 s near 4.398 | 9.99994 s median; two 20.0000 s intervals at 4.398148 to 4.398611 | Reproduced |
| Headers | Section 4 spellings | 3 of 4 differ | Discrepancy D1, explained in M0; loader maps the actual headers and raises on any other |
| Common window | 0 to 47.92 | grid of 4,792 bins, no empty bin for any sensor | Reproduced |

Gap inventory: no true gaps in any sensor within the common window. pH: one run of two missing-sample intervals (20 s) at day 4.398. NH4/ORP: only partial-stamp artifacts (first NH4 stamp, first ORP stamp boundary, last ORP stamp at 53.23, outside the grid).

### Gate M1: PASS

- All tests pass (51 passed), including: row counts and spans, dropped count = 1, EC(−0.1 A) = 0.810725, mA input raises, stamp-size ranges 22–450 and 21–87, strictly increasing reconstructed times, no fabricated grid values (empty bins NaN, counts equal native counts, every bin mean inside its samples' min–max).
- Every Section 4 property is reproduced, except the header spellings (D1) and the nominal 60 s Cond interval (actually 59.83 s), both explained above.

## M2. Sensor validation with the known pulses (2026-10-18)

### Implemented

- `src/bnr_neq/thermo.py` (M2 subset of Section 6: pulse dose, Kohlrausch steps, Nernst slope, chloride ORP rise), `src/bnr_neq/qc.py` (rolling-median step detector, jump characterization and classification, held-output detector, late-EC-pulse detector with decay fits, step locator), `src/bnr_neq/validation.py` (pulse table, balances, clock check, reference comparison, verdicts, figure, report); `validate` and `run-all` wired in `__main__.py`.
- Tests: `tests/{test_thermo,test_qc,test_validation}.py` (100 tests in total with M1).
- Outputs: `outputs/tables/sensor_validation.csv` (one row per pulse, units in column names, ORP predictions labelled `SCENARIO_Cl0_<c>mM_T<T>C`), `sensor_validation_windows.csv`, `sensor_balances.csv`, `nh4_jumps.csv`, `nh4_held_segments.csv`, `nh4_hf_residual_subsegments.csv`, `ec_late_pulses.csv`, `ec_late_pulse_sensitivity.csv`, `clock_check.csv`, `sensor_validation_references.csv`; `outputs/figures/fig3_sensor_validation.{png,svg}`; `outputs/reports/sensor_validation.md`.

### Commands run

```
uv run python $SP/probe.py                 # scratch probe of the detectors: found (a) argmax on the flat-topped
                                           #   median-difference step gives arbitrary jump times, (b) curve_fit p0 outside bounds
# fix: jumps located at the 50 % crossing; "sustained" criterion added; p0 amplitude clipped to >= 1e-4
uv run python -m bnr_neq validate
# figure legends moved twice after viewing the PNG (collisions)
uv run pytest -q                           # 1 failed, 98 passed: synthetic EC test with the quantile rule found a noise peak
# test split: detector tested with an explicit threshold; quantile rule tested separately; limitation documented
uv run ruff check --fix --unsafe-fixes .   # 12 findings (ISC004 parentheses, RUF046), behaviour-preserving
# gate run (progress.md copied to the scratchpad and restored, because rm -rf outputs deletes it):
uv sync                                    # Resolved 30 packages, Checked 28 packages
uv run ruff check .                        # All checks passed!
uv run pytest -q                           # 100 passed in 8.29s
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all           # second run
diff -rq <run1 tables> outputs/tables      # all 14 CSV tables identical
```

### Reference comparison (from `outputs/tables/sensor_validation_references.csv`)

| Quantity | Reference | Obtained | Tolerance | Result |
| --- | --- | --- | --- | --- |
| NH4 response fraction, P1 to P4 | 37, 29, 26, 17 % | 37.0, 29.0, 26.2, 17.4 % | ±3 pp | Pass |
| Observed NH4 change, P1 to P4 | +31.1, +24.4, +22.1, +14.7 | +31.23, +24.45, +22.07, +14.70 mg NH4+/L | ±2 | Pass |
| Observed EC change, P1 to P4 | +0.005, −0.062, −0.007, +0.016 | +0.0052, −0.0606, −0.0067, +0.0164 mS/cm | ±0.01 | Pass |
| Observed ORP change, P1 to P4 | +0.4, 0.0, +0.1, +0.2 | +0.400, 0.000, +0.125, +0.210 mV | ±0.5 | Pass |
| Observed pH change, P1 to P4 | −0.01, +0.02, 0.00, +0.04 | −0.01, +0.02, 0.00, +0.04 | ±0.02 | Pass |
| EC median day 0.3–1.3 / 48.5–49.6 | 1.069 / 1.071 | 1.0690 / 1.0710 mS/cm | ±0.005 | Pass |
| NH4 median start / end | 47.1 / 43.1 | 47.02 / 43.13 mg NH4+/L | ±0.5 | Pass |
| Unrelated NH4 jumps | 2.59, 17.84, 37.19, 44.63 | 2.626, 17.875, 37.219, 44.660 d | ±0.02 d | **Not reproduced (explained below)** |
| Delayed responses P1, P2 | 8.85, 10.88 (~7.5 h) | 8.858 (7.88 h), 10.883 (7.52 h) | ±0.02 d | Pass |
| Late EC pulses after 41.9 | 29, median decay ~24 min | 29, median decay 23.1 min | count ±3 | Pass |

Derived (no tolerance given): N added 261.86 mg N/L; removed from the NH4 pool 264.88 mg N/L over 46.6 d; implied 5.68 mg N/L/d (plan: "about 5.6"). Predicted end EC: 3.87 (no nitrification) and 2.17 mS/cm (complete nitrification), limiting values; observed 1.071.

**Explanation of the four jump-time discrepancies (evidence).** In the raw NH4 rows, each of the four unrelated jumps occurs between two consecutive rows inside one stamp: 2.62 (last row 26.73 after 42.89 at stamp start), 17.87 (61.09 → 34.01), 37.21 (67.71 → 38.64), 44.66 (36.76 → 42.90). The detected 50 % crossings (2.626, 17.875, 37.219, 44.660) lie inside those stamps; `tests/test_validation.py::test_offset_jumps_lie_in_the_stamp_where_the_raw_value_changes` checks this independently of the detector. The reference days are 0.029 to 0.036 d earlier, outside every stamp holding the change. They coincide within ±0.014 d with the onset of the detector's above-threshold run (2.600, 17.850, 37.194, 44.636), which is reported as `run_start_d`; the preliminary analysis probably reported the window onset. The stamp-rounding question (open observation 1) shifts the true times by at most −0.005 d and does not explain a 0.03 d difference. Nothing was tuned.

**Additional detections not in the reference list.** The detector also finds four sustained level changes classed `electrode_offset` in the undamped day 2.6 to 8.6 segment (3.214 +5.7, 5.571 −7.4, 6.547 −4.3, 8.553 −5.4 mg NH4+/L), each with a 27 to 51 min 10–90 % transition, against ≤1 sample (0.05–0.34 min) for the four reference offsets. It also finds smaller responses to pulse 1 (8.752, +4.3, lag 5.3 h) and pulse 2 at the pulse (10.576, +7.9, lag 0.14 h), and the pulse 3 (24.665, lag 1.8 h, 92 min transition) and pulse 4 (32.396, lag 0.14 h) responses. Classification uses the a priori rule only.

**Held/damped output.** σ_HF < 0.01 mg NH4+/L for ≥ 0.5 d: 0.2–0.7, 0.8–1.9, 13.0–18.0, 20.7–21.7, 21.8–24.7 d. Per sub-segment σ_HF: 0.80 for 2.59–8.82 (plan D3: 0.66 as RMS-like residual; my RMS 0.92), 0.0 for 10.85–17.84 and 17.84–24.60 (plan: zero), 0.015 for 32.36–37.19 (plan: zero; here 1.5 x the output resolution, so not flagged), 0.0 for 0.3–2.59 (not in plan).

**Clock check (D7).** Step at the multi-probe event: Cond 2.9183, ORP 2.9099, pH 2.8125 d; pH leads Cond by 0.106 d. Not corrected; recorded as an open issue. **Superseded in M3:** the 0.01 d step window latched onto a transient pH excursion at 2.81 to 2.835 d; with a 0.05 d window pH and ORP step within 0.001 d of each other (see the M3 entry). No evidence of a 0.1 d pH clock offset.

### Gate M2: PASS WITH FOUR EXPLAINED DISCREPANCIES

- The table and report exist and reproduce 27 of 31 reference values within tolerance. The four not reproduced (unrelated-jump days) are explained above with raw-data evidence and an independent test. The reviewer decides whether the explanation is accepted.
- The report states the P0 findings: NH4 under-response (17 to 37 %, declining), no EC response (−0.061 to +0.016 mS/cm against +0.700 limiting, no cumulative rise), no ORP chloride response (0.00 to +0.40 mV against 4.0 to 44.6 mV across scenarios). Verdicts: NH4 FAIL, conductivity FAIL, ORP (chloride test) FAIL, pH PASS (no contradiction; quantitative prediction not computed: alkalinity_mmol_L missing).
- Stopped here as instructed; M3 not started.

## M3. Segmentation and consumption kinetics (2026-10-18)

User decision before M3: "Accept the M0–M2 deviations as documented and proceed to M3."

### Implemented

- `src/bnr_neq/segments.py`: boundaries (start day 1, pulses, every M2 jump), sub-segments with a 0.1 d exclusion and a 0.25 d minimum fit span, Theil-Sen slopes on grid bins wholly inside each fit window, gain-corrected consumption (labelled, with an explicit `not computed: no governing pulse (baseline)`), free-ammonia medians per temperature SCENARIO, window check, diagnostic fits on the exact preliminary windows, kinetic fits (zero, first, Monod via the exact Lambert-W solution, Haldane in free ammonia via RK4, per SCENARIO), Jacobian-based identifiability, residual lag-1 autocorrelation and Ljung-Box, PELT change points, Lomb-Scargle, fig2, report. `thermo.py`: `emerson_pka`, `free_ammonia_mgNH3_L` with input guards. CLI `segments` implemented; `run-all` runs ingest, validate, segments.
- Tests: `tests/test_segments.py` (28), speciation tests in `tests/test_thermo.py`, clock-check test in `tests/test_validation.py`; 137 in total.
- Outputs: `outputs/tables/{segments,segments_window_check,segments_reference_windows,segments_references,kinetic_fits,kinetic_parameters,changepoints,changepoint_sensitivity,lomb_scargle_summary}.csv`, `outputs/processed/lomb_scargle_spectra.parquet`, `outputs/figures/fig2_overview.{png,svg}`, `outputs/reports/segments.md`.

### Correction to M2 (code change in `validation.py`)

While checking PELT events I found that the M2 clock check was wrong. With a 0.01 d step window it located the pH "step" on a transient +0.45 excursion at days 2.81 to 2.835, which returns to baseline. The sustained pH rise is at 2.86 to 2.91, coincident with the ORP drop at 2.905 to 2.91 (raw 0.01 d medians inspected). The clock check now uses the jump detector's 0.05 d median window, a 2.75 to 3.00 d search range and ORP as reference (EC shows only a slow drift there); its report text is computed, not hard-coded. Result: pH 2.8900, ORP 2.8893 d (offset +0.0007 d). The M2 statement "pH leads by 0.106 d" is marked superseded above; there is no evidence of a gross pH clock offset, though the check cannot resolve offsets below about half its window. New test: `test_clock_check_ph_and_orp_steps_coincide`.

### Method choices fixed before looking at the outcome

- Every M2 jump (offsets and delayed pulse responses) is a boundary; this follows from the gate ("no slope window spans a pulse or an uncorrected jump").
- PELT: the mean-shift cost staircased the drifting series (hundreds of breaks at any penalty, probed in scratch), so the cost is discontinuous piecewise-linear (`model="linear"` on [t, 1]), BIC penalty 3 σ² ln n, σ from 1-day rolling-median residuals, event if the level discontinuity is at least 5 σ.
- Lomb-Scargle: linear detrend over the common window; the full band (1/T to Nyquist) is dominated by the trend (max 0.21 to 0.63 at periods of 11 to 40 d); the periodicity band requires at least 10 cycles in the record (f ≥ 0.209 c/d). Both bands are reported.

### Commands run

```
uv run python $SP/m3probe*.py              # scratch probes of PELT costs/penalties and LS detrending/bands
uv run python $SP/m3run.py                 # scratch run of slope table, reference windows, kinetics
uv run python -m bnr_neq segments          # 11/13 references
uv run python -m bnr_neq validate          # after the clock-check correction
# fig2 legends moved twice after viewing the PNG
uv run pytest -q                           # 1 failed (test used wrong parameter names), 136 passed; test fixed
uv run ruff check --fix --unsafe-fixes .   # ISC004 parentheses, RUF007 pairwise; behaviour-preserving
# gate run (progress.md parked in the scratchpad and restored):
uv sync                                    # Resolved 30 packages, Checked 28 packages
uv run ruff check .                        # All checks passed!
uv run pytest -q                           # 137 passed in 56.31s
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all           # second run
diff -rq <run1 tables> outputs/tables      # all 23 CSV tables identical
```

### Reference comparison (`outputs/tables/segments_references.csv`, `segments_reference_windows.csv`)

| Quantity | Reference | Production sub-segment | Exact preliminary window | Tolerance | Result |
| --- | --- | --- | --- | --- | --- |
| Slope 0.3–2.59 | −1.89 | −1.671 (1.100–2.626) | −1.914 | ±0.189 | **Not reproduced (explained)** |
| Slope 2.59–8.82 | −3.42 | −1.781 (3.314–5.571) | −3.406 (spans pulse 1 and jumps) | ±0.342 | **Not reproduced (explained)** |
| Slope 10.85–17.84 | −1.95 | −1.939 | −1.951 | ±0.195 | Pass |
| Slope 17.84–24.60 | −0.57 | −0.578 | −0.573 | ±0.057 | Pass |
| Slope 24.60–32.36 | −0.48 | −0.484 | −0.480 | ±0.05 | Pass |
| Slope 32.36–37.19 | −0.31 | −0.296 | −0.304 | ±0.05 | Pass |
| Slope 37.19–44.63 | −0.28 | −0.283 | −0.282 | ±0.05 | Pass |
| Slope 44.63–47.92 | +0.12 | +0.126 | +0.115 | ±0.05 | Pass |
| Gain-corrected after pulse 2 | ~5.2 mg N/L/d | 5.191 | | ±10 % (my reading of "about") | Pass |
| Gain-corrected after pulse 4 | ~1.4 mg N/L/d | 1.319 | | ±10 % | Pass |
| Free ammonia 17.84–32.36, 20/25 °C SCENARIOS | 0.7 to 1.2 mg NH3/L | 0.702 to 1.168 | | ±0.05 (rounding) | Pass |
| Max LS normalized power | ~0.05 | 0.0516 (ORP, 0.346 c/d, periodicity band) | full band 0.628 | ±0.01 | Pass (band-dependent, see above) |

The Theil-Sen method itself reproduces all eight reference slopes on the exact preliminary windows. The two production discrepancies are rule-driven: (1) the prompt's segmentation starts at day 1 (start-up excluded), not 0.3; (2) the preliminary 2.59–8.82 window contains pulse 1 (8.53) and the M2 jumps at 3.214, 5.571, 6.547, 8.553, which the M3 gate forbids. The 2.59–8.82 sub-segments give −21.5 (2.73–3.21, recovery after the 2.626 jump), −1.78, −1.13 and +0.99 mg NH4+/L/d.

Other results: kinetic fits leave residual lag-1 autocorrelation 0.82 to 0.998 and Ljung-Box p < 1e-196 in every segment, so AIC differences are not reliable. Identifiable parameters: zero-order 4/5, first-order 4/5 (segment 1, the rising 8.96–10.57 window, fails), Monod 0/10, Haldane 0/30. PELT events: EC 13, NH4 55, ORP 36, pH 18 (all but 10 not tied to a pulse); they include the plan's ORP steps (28.905 +100 mV, 30.505 −44 mV), the post-pulse-2 ORP decline (10.705, tied), EC steps (9.455, 15.705, 35.955) and the 2.905 pH/ORP event. Near 1 c/d the pH power is 0.0053, i.e. no daily cycle above the background.

### Gate M3: PASS WITH TWO EXPLAINED DISCREPANCIES

- No slope window spans a pulse or a detected jump: 12 of 12 computed windows clean (`segments_window_check.csv`; `test_no_slope_window_spans_a_pulse_or_jump`).
- 11 of 13 references reproduced on the production windows; the other two are reproduced on the exact preliminary windows, and their production values differ for the rule-driven reasons above. The reviewer decides whether the explanation is accepted.
- Stopped here as instructed; M4 not started.

## M4. Thermodynamics and conservation balances (2026-10-18)

User decision before M4: "Accept the M3 discrepancies and the M2 clock-check correction as documented and proceed to M4."

### Implemented

- `src/bnr_neq/thermo.py` (Section 6): `davies_gamma`, `e_ref_agcl_V`, `e_ref_ag2o_V`, `eh_from_orp_V`, `p_o2_atm`, `nitrification_dG_kJ_mol`, `nitrification_affinity_kJ_mol`, `affinity_step_kJ_mol`, `alkalinity_consumed_mgCaCO3_L`, `oxygen_demand_mgO2_L_d`, `required_kla_per_d`, `hco3_air_eq_mmol_L`, `omega_co2`, `de_donder_consistent`, `entropy_production_W_L_K`, `jin_bethke_FT` (plus M2/M3 functions). Guard pattern: `NotComputed` and `*_guarded` forms for every function that needs alkalinity, DO, NO3-N, chloride, the ORP offset or kLa (`eh_from_orp_guarded` returns exactly `not computed: S4`). M4 stage: affinity steps, coexistence SCENARIO with DO sensitivity, chloride SCENARIOS and the oxide alternative, CO2, oxygen, nitrogen/chloride/alkalinity balances, references, report. CLI `thermo` implemented.
- `config/experiment.toml`: added `[scenarios.coexistence]` (NH4 40, NO3 100 mg N/L, DO 0.1 mg/L, pH 7, 25 °C, DO sensitivity list), labelled as SCENARIO. **Deviation:** extends the Section 5 file, so that these scenario values live in the configuration, not the code.
- `thermo.run` now loads its inputs before creating output directories (found by a test).
- Tests: `tests/test_thermo_m4.py` (36), 173 in total.
- Outputs: `outputs/tables/{thermo_affinity_steps,thermo_coexistence,thermo_chloride_scenarios,thermo_chloride_observed,thermo_co2,thermo_oxygen,balance_nitrogen,balance_chloride,balance_alkalinity,thermo_references}.csv`, `outputs/reports/thermo.md`. No figure is required for M4.

### Commands run

```
uv run python -c "...affinity steps from sensor_validation_windows.csv..."   # 4.172, 2.471, 3.138, 2.365 kJ/mol
uv run python -c "...spot checks..."     # coexistence 316.89; Davies 0.90210; HCO3 0.1268/0.4011; AgCl 0.2223 V, 59.16 mV/dec
uv run python -m bnr_neq thermo           # 9/11 references
uv run pytest -q tests/test_thermo_m4.py ...   # 1 failed (run created dirs before checking inputs) -> code fixed
uv run ruff check --fix --unsafe-fixes .  # 7 style findings (ISC004, I001), behaviour-preserving
# gate run (progress.md parked in the scratchpad and restored):
uv sync                                   # Resolved 30 packages, Checked 28 packages
uv run ruff check .                       # All checks passed!
uv run pytest -q                          # 173 passed in 54.45s
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all          # second run
diff -rq <run1 tables> outputs/tables     # all 33 CSV tables identical
```

### Reference comparison (`outputs/tables/thermo_references.csv`)

| Quantity | Reference | Obtained | Tolerance | Result |
| --- | --- | --- | --- | --- |
| ΔA pulses 1–4 (25 °C) | 4.2, 2.5, 3.1, 2.4 kJ/mol | 4.172, 2.471, 3.138, 2.365 | ±0.1 | Pass |
| Coexistence affinity (SCENARIO) | 317 kJ/mol | 316.895 | ±2 | Pass |
| HCO3- air equilibrium, pH 7.3 / 7.8 | 0.127 / 0.401 mmol/L | 0.1268 / 0.4011 | ±0.001 (my choice) | Pass |
| O2 demand min / max (GAIN-CORRECTED) | ~6 / ~33 mg O2/L/d | 5.763 / 23.723 | ±10 % (my reading of "about") | min Pass, **max not reproduced (explained)** |
| Required kLa min / max | ~0.7 / ~4 d⁻¹ | 0.698 / 2.872 | ±10 % | min Pass, **max not reproduced (explained)** |

Unit tests required by M4 (all pass): pKa(25 °C) = 9.2464; Davies γ(z = 1, I = 0.01) = 0.90210; AgCl/Ag E = 0.2223 V at unit activity and 59.16 mV per decade; impossible inputs (pH outside 0–14, negative concentrations, DO ≤ 0, zero C*, T ≤ 0, mismatched shapes) raise. Every missing-input path is tested both ways: missing gives the exact marker (`not computed: S4`, `not computed: chloride_mmol_L missing`, `DO_mg_L`, then `NO3_N_mg_L` when DO is supplied, `alkalinity_mmol_L`, `kLa_per_d`), and a temporary configuration with the keys moved to `[measurements]` gives numbers.

**Explanation of the two upper-end discrepancies (evidence).** The production oxygen range uses the gain-corrected rates of M3 slope windows with a governing pulse and positive consumption (1.26 to 5.19 mg N/L/d). The preliminary upper end (~7.2 mg N/L/d) came from the 2.59–8.82 window, which spans pulse 1 and three NH4 jumps, divided by pulse 1's response fraction. The stage computes this as a labelled `DIAGNOSTIC ONLY` row in `thermo_oxygen.csv`: 7.14 mg N/L/d, i.e. 32.6 mg O2/L/d and kLa 3.95 d⁻¹, which reproduces the references 33 and 4. That window is excluded by the M3 segmentation rule and does not enter the reference range (test `test_preliminary_window_diagnostic_explains_upper_references`).

Other results: in situ affinity `not computed: DO_mg_L missing`; DO from 0.01 to 8 mg/L moves the SCENARIO affinity from 305.5 to 338.6 kJ/mol; affinity steps are 0.7 to 1.2 % of |ΔG°′|. Chloride SCENARIOS predict 4.0 to 44.6 mV per pulse and 22.0 to 76.6 mV cumulative; the oxide-controlled E_ref at pH 7.26 is 0.744 V vs SHE against 0.331 to 0.400 V for chloride control; absolute E_H `not computed: S4`. Balances: nitrogen closure, chloride baseline and final, remaining alkalinity and Ω_CO2 are `not computed` with their keys; full nitrification of the dose would consume 1,870 mg CaCO3/L (37.4 meq/L), 467 mg/L (9.34 meq/L) per pulse.

### Gate M4: PASS WITH TWO EXPLAINED DISCREPANCIES

- 9 of 11 references reproduced; the remaining two (upper ends of oxygen demand and kLa) are reproduced by the computed diagnostic on the preliminary window, which the M3 rule excludes. The reviewer decides whether the explanation is accepted.
- Every missing-input path is tested (guards for alkalinity, DO, NO3-N, chloride, ORP offset and kLa; stage outputs carry the markers).
- Stopped here as instructed; M5 not started.

## M5. Dynamical and information-theoretic signatures (2026-10-18)

User decision before M5: "Accept the M4 discrepancies and the config extension as documented and proceed to M5."

### Implemented

- `src/bnr_neq/estimators.py` (NumPy only; a test checks that it imports neither pandas nor SciPy): moment asymmetry, tie-aware ordinal codes (pairwise-sign encoding: equalities are their own symbols, 13 weak orders for d = 3), ordinal KL irreversibility, tie fraction, tie-aware permutation entropy, directed horizontal visibility graph degrees and KL, IAAFT surrogates and rank tests (one shared surrogate set per window), sample entropy, multiscale entropy, DFA, ring Markov chain simulation and analytic entropy production, transition-level entropy production with moving-block bootstrap basic intervals, Bowker-type symmetry G-test, circulation rate and Battle-type flux field.
- `src/bnr_neq/dynamics.py`: the estimator validation suite runs first, and `compute` raises `EstimatorValidationError` before any plant-data result if a required check fails (tested). Then: windowed analysis (0.25 d, raw ADF/KPSS-stationary windows only, IAAFT tests, complexity), Onsager regression test, stability indicators, probability flux, the electrode-dynamics guard, figure, report. CLI `dynamics` implemented.
- Tests: `tests/test_estimators.py` (20), `tests/test_dynamics.py` (15); 208 in total.
- Outputs: `outputs/tables/{dynamics,dynamics_windows,estimator_validation,dynamics_onsager,dynamics_onsager_lags,dynamics_stability,dynamics_flux,dynamics_electrode}.csv`, `outputs/figures/fig5_dynamics.{png,svg}`, `outputs/reports/dynamics.md`.

### Problems found and fixed before any result was reported

1. **Detrending artifact.** The first version tested windows raw and, failing that, linearly detrended. Detrended pH windows gave a median ordinal KL of 4.8 nats (49/49 "significant") against 0.014 for raw stationary windows, and NH4/EC detrended windows were inflated too. Subtracting a ramp from a quantized staircase makes a sawtooth, and leftover curvature in a relaxation is time-asymmetric. Fix: raw level-stationary windows only.
2. **HVG bug with ties (production code).** When the stack top equalled the new value, the code linked to it and also to the taller element behind it, which the equal value blocks. A brute-force cross-check test on tied integer series found it. Before the fix, the tie sweep showed HVG false positives of 15–20/20 at tie fractions of 0.52 to 0.91, which I first took for a property of HVG. After the fix, HVG keeps 0/20 at every tie level up to 0.949. Two of my new tests had encoded the buggy behaviour (HVG range < 0.5; HVG never computed on pH). They now assert the generic invariant (each index is computed only inside its validated tie range).
3. **Validation did not cover pH's tie regime.** The first sweep stopped at a 0.245 tie fraction, while pH has 0.82 to 0.96. The sweep now runs from continuous data to 0.949. Each index gets a validated tie range (contiguous levels with at most 4/20 false positives), enforced per window with an explicit `not computed: tie fraction ... above validated range`.
4. **Onsager SNR.** The first version measured the tail deviation from the late level. For pulse 2 that is −0.83 pH, the segment's 14-day drift, not a pulse response. The SNR now uses the pulse response measured in M2 (−0.01 to +0.04 pH).
5. Smaller items: excluded-pair counting in `entropy_production_from_counts` now counts unordered pairs; one of my tests had the direction of multiscale entropy for white noise wrong (fixed in the test); one test input was a 5 % false KPSS rejection (seed changed, with a comment).

### Estimator validation (`outputs/tables/estimator_validation.csv`): 25 of 25 required checks pass

| Check | Result |
| --- | --- |
| Ring chain, transition-level EP, 99 % basic bootstrap CI, 200,000 steps | (0.3, 0.3): 0.00005 [−0.00012, 0.00010], truth 0; (0.4, 0.2): 0.1365 [0.1308, 0.1421], truth 0.1386; (0.6, 0.1): 0.8920 [0.8784, 0.9050], truth 0.8959; (0.45, 0.05): 0.8826 [0.8674, 0.9032], truth 0.8789 |
| Detailed-balance G-test at p = q | p = 0.175 |
| Noisy observable quantized to 0.01 (17 % misclassified states) | p = q stays zero (CI includes 0); bias at p ≠ q is −0.100, −0.650, −0.699 per step (−72 % to −80 %): it underestimates |
| IAAFT, reversible AR(1), tie fractions 0 to 0.949 | false positives ≤ 2/20 for all three indices at every level; Wilcoxon on index − median(surrogates) p > 0.01 at the required levels; validated tie range 0.949 for every index |
| Positive control (ring 0.6/0.1 + noise) | power 20/20 for all three indices |
| Flux, 2D OU | reversible: CI includes 0; rotational ω = 1: CI > 0 |
| Complexity sanity | PE(iid) = 0.6985 (= ln 6/ln 13); DFA white 0.510, random walk 1.551; SampEn white 2.159 |

### Results on the data (only after validation)

- Stationary 0.25 d windows: pH 12/187, ORP 2/185, NH4 6/135, EC 8/187. Most windows are non-stationary, so results are window-level statements.
- pH (primary; tie fractions 0.82 to 0.96): of the 10 stationary windows within the validated tie range, ordinal KL is significant against IAAFT in 10, moment asymmetry in 8 (median normalized lag-1 asymmetry strongly negative: abrupt drops, slow recoveries) and HVG in 7. Two stationary pH windows exceed the 0.949 tie range and are not computed.
- Relative ORP, offset-corrected NH4 (provisional) and EC (probe behaviour only): see `dynamics.csv`.
- Onsager regression (pH): pulse-response SNR 1.81, 0.19, 0.00, 1.12, all below 3. **Not interpretable; no conclusion drawn (H4 untested).**
- Stability indicators: computed and labelled EXPLORATORY (overlapping windows, four pulses).
- Probability flux: no significant circulation in any of 10 segment × plane combinations.
- Deconvolved NH4 relaxation: `not computed: electrode_step_response missing`.

### Commands run

```
uv run python $SP/m5probe.py              # scratch: Markov CI coverage, quantized bias, IAAFT false positives, power, timing
uv run python -c "...validate_estimators()..."   # first suite: all required pass
uv run python $SP/m5run.py                # first data run: implausible pH ordinal KL (2 to 8 nats)
uv run python $SP/m5diag.py $SP           # diagnosis: detrending artifact; HVG tie failures (later traced to the bug)
# fixes 1, 3, 4 above
uv run python $SP/m5run2.py $SP
uv run python -m bnr_neq dynamics         # 25/25 required; figure adjusted once after viewing
uv run pytest -q                          # 8 failed / 200 passed: HVG brute-force mismatch (bug 2), MSE test direction, pair counting
# fix 2 and 5; validation rerun: HVG 0/20 false positives at all tie levels
uv run pytest -q                          # 1 failed (KPSS false rejection of the test input); seed changed
uv run ruff check --fix --unsafe-fixes .  # 12 style findings, behaviour-preserving
# gate run (progress.md parked in the scratchpad and restored):
uv sync                                   # Resolved 30 packages, Checked 28 packages
uv run ruff check .                       # All checks passed!
uv run pytest -q                          # 208 passed in 90.90s
rm -rf outputs && uv run python -m bnr_neq run-all   # dynamics: estimator validation 25/25 required checks passed
uv run python -m bnr_neq run-all          # second run
diff -rq <run1 tables> outputs/tables     # all 41 CSV tables identical
```

### Gate M5: PASS

- All required estimator validation tests pass (25/25), and the code refuses plant-data results otherwise (`test_gate_blocks_plant_data`).
- ADF and KPSS run on every window before any stationary estimator; no NH4 relaxation is presented as a process time; conductivity appears only as probe behaviour.
- Stopped here as instructed; M6 not started.

## M6. Reduced model with the electrode observation model (2026-10-18)

User decision before M6: "Accept the M5 fixes and results as documented and proceed to M6."

### Implemented

- `src/bnr_neq/model.py`: state model (true NH4 in mg N/L from day 1, +65.464 mg N/L exactly at each pulse, steps split at the pulse time, RK4 on the 0.01 d grid; rates non-negative by construction). Consumption laws: Haldane in free ammonia (per temperature SCENARIO), first-order decay of biomass activity, both combined, plus a no-consumption reference. Observation model: log R = a + (b0 + b1 (t − 1)) log C + offsets, then a first-order lag τ (exact exponential filter, `scipy.signal.lfilter`). `scipy.optimize.least_squares` in log space. Chronological split, whole-record refit of the selected model, profile likelihood, consistency checks, model-implied response fractions, figure, report. CLI `model` implemented.
- `constants.py`: `K_NH_ASM_MG_N_L = 1.0` (ASM1 default), fixed in the biomass-decay rate law and declared as fixed.
- Tests: `tests/test_model.py` (15: exact pulse dose, mass-balance closure, RK4 against `solve_ivp` for the decay and Haldane laws, lag and offset behaviour, parameter recovery, profiler flags non-identifiable parameters, effective n, split evaluation on a synthetic record, consistency checks, scenario labels, missing inputs); 223 in total.
- Outputs: `outputs/tables/{model_fit,model_identifiability,model_checks,model_response_fractions}.csv`, `outputs/processed/model_state.parquet`, `outputs/figures/fig6_model.{png,svg}`, `outputs/reports/model.md`.

### Deviations and choices (decided before the fits unless stated)

- Offsets at all eight electrode-offset jumps detected in M2 (2.626, 3.214, 5.571, 6.547, 8.553, 17.875, 37.219, 44.660), not four, consistent with the accepted M2 detections and the M3 segmentation.
- In the test window, the two offsets that fall there (unseen in calibration) are refitted as the only nuisance parameters with everything else frozen; the pure forecast (nothing refitted) is reported from 32.39 to the first unseen offset (37.219). No other test information enters any fit.
- Profile likelihood uses Δ = n_eff ln(RSS/RSS_min) with n_eff from the residual lag-1 autocorrelation (ρ = 0.954, n_eff = 110 of 4,692 bins), so autocorrelated residuals do not overstate identifiability. Grid: 3 points per side, warm-started; per-point evaluations capped at 40, and any closing point that did not converge is refitted without the cap (4 such refits); a side counts as closed only by a converged point. Added after a first profile run, to rule out unconverged re-optimizations overstating Δ; the conclusions did not change.

### Out-of-sample skill (`model_fit.csv`; skill = 1 − MSE/MSE(no-consumption), log units)

| Model | Cal RMSE | Sel RMSE | Test RMSE | Forecast RMSE (32.39–37.22) | Sel skill | Test skill | Forecast skill |
| --- | --- | --- | --- | --- | --- | --- | --- |
| no consumption | 0.118 | 0.347 | 0.305 | 0.543 | 0 | 0 | 0 |
| **biomass_decay (selected)** | 0.114 | **0.071** | 0.153 | 0.201 | **0.958** | 0.748 | 0.863 |
| Haldane FA, T 20 / 25 °C SCENARIO | 0.112 | 0.139 / 0.138 | 0.216 | 0.282 / 0.281 | 0.839 / 0.841 | 0.499 / 0.500 | 0.730 / 0.731 |
| combined, T 20 / 25 °C SCENARIO | 0.112 | 0.239 / 0.240 | 0.358 / 0.359 | 0.560 / 0.562 | 0.527 / 0.524 | −0.377 / −0.383 | −0.065 / −0.070 |

The combined model fits calibration best but generalizes worst (overfitting). Whole-record refit of biomass_decay: RMSE 0.0946 log units.

### Identifiability (`model_identifiability.csv`)

Identifiable (approximate 95 % profile intervals): v0 = 20.0 mg N/L/d [5.0, 63.3]; kd = 0.130 d⁻¹ [0.034, 0.195]; b0 = 1.324 [1.181, 1.358]; b1 = −0.0009 d⁻¹ [−0.0017, 0.0012]; τ = 0.262 d [0.167, 0.387]; all eight offsets. **Not identifiable:** C0 (profile not closed on the high side) and a (not closed on the low side). These are reported as "not identifiable", with optimizer values only in a column labelled `optimizer_value_not_an_estimate`. Interpretation, conditional on the model: the absolute NH4 level is undetermined, so NO-GO for absolute concentrations stands. The declining response fractions are reproduced (model 0.35, 0.28, 0.18, 0.19 against observed 0.37, 0.29, 0.26, 0.17) by NH4 accumulation under a log-type electrode with a non-drifting gain, together with decaying activity. The lag cannot reproduce the 7.5 h delayed jumps, and the residuals remain strongly autocorrelated.

### Consistency checks (`model_checks.csv`): all pass

Each of the four simulated pulses adds 65.46428 mg N/L (difference from the dose < 1e-13); mass balance C_end − (C0 + pulses − consumed) = 4.5e-13; minimum rate 0.046 mg N/L/d ≥ 0; σ = r A / T ≥ 0 at every grid time with the coexistence SCENARIO affinity (316.9 kJ/mol); in situ σ `not computed: DO_mg_L missing`.

### Commands run

```
uv run python $SP/m6speed.py              # pure-Python RK4: 9 ms per whole-record simulation
uv run python $SP/m6run.py $SP            # split evaluation of six models (82 s)
uv run python $SP/m6prof.py $SP           # refit + profiles; rerun after making profiles convergence-aware
uv run python -m bnr_neq model            # 182 s; checks pass; 13/15 identifiable
uv run pytest -q tests/test_model.py      # 15 passed
uv run ruff check --fix --unsafe-fixes .  # style fixes; B023 closure fixed by passing arguments explicitly
# gate run (background; progress.md parked in the scratchpad and restored):
uv sync                                   # Resolved 30 packages, Checked 28 packages
uv run ruff check .                       # All checks passed!
uv run pytest -q                          # 223 passed in 112.81s
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all          # second run
diff -rq <run1 tables> outputs/tables     # all 45 CSV tables identical
```

### Gate M6: PASS

- Mass-balance tests pass (unit tests and the real-data checks above).
- Out-of-sample skill is reported for every model, from the chronological split, before the whole-record refit.
- The identifiability table is complete: every one of the 15 parameters has a profile verdict, and the two non-identifiable ones are not presented as estimates.
- Stopped here as instructed; M7 not started.

## M7. Pulse-level efficiency and energy descriptors (2026-10-18)

User decision before M7: "Accept the M6 model, deviations and identifiability results as documented and proceed to M7."

### Implemented

- `src/bnr_neq/pulses.py`: tidy pulse table (one row per pulse and metric with value, unit, basis, inputs and status), exploratory associations, figure, report. CLI `pulses` implemented. Energy metrics are computed when `stirrer_power_W` is measured (stirrer energy = P x segment duration; specific stirring energy per MODEL-BASED g N converted; ratio of MODEL-BASED free energy released to stirrer energy), and are `not computed: stirrer_power_W missing` otherwise. Both paths are tested.
- Tests: `tests/test_pulses.py` (11; the integration tests use the real M2–M4 results with a synthetic model state); 234 in total.
- Outputs: `outputs/tables/pulses.csv` (92 rows: 74 computed with inputs listed, 18 not computed), `outputs/tables/pulses_associations.csv`, `outputs/figures/fig7_pulses.{png,svg}`, `outputs/reports/pulses.md`.

### Pulse table (computed values; units in `pulses.csv`)

| Metric | P1 | P2 | P3 | P4 |
| --- | --- | --- | --- | --- |
| ΔA (kJ/mol, 25 °C) | 4.172 | 2.471 | 3.138 | 2.365 |
| Consumption, uncorrected, duration-weighted (mg N/L/d) | −0.336 | 0.988 | 0.376 | 0.157 |
| Consumption, GAIN-CORRECTED, first window / weighted | −0.906 / −0.906 | 5.191 / 3.407 | 1.435 / 1.435 | 1.319 / 0.901 |
| Median pH (slope windows) | 7.005 | 7.471 | 7.547 | 7.355 |
| Median FA, 20 / 25 °C SCENARIO (mg NH3/L) | 0.200 / 0.287 | 0.501 / 0.717 | 0.702 / 1.003 | 0.356 / 0.510 |
| O2 demand (mg O2/L/d), GAIN-CORRECTED | not computed (non-positive) | 15.57 | 6.56 | 4.12 |
| Required kLa (1/d, 25 °C) | not computed | 1.885 | 0.794 | 0.499 |
| η, MODEL-BASED (to next pulse; P4 to 47.92 d) | 0.205 | 0.569 | 0.070 | 0.035 |
| Free energy released, MODEL-BASED (J; full pulse = 163.0 J) | 33.4 | 92.7 | 11.5 | 5.7 |

Not computed (explicitly): measured/required kLa ratio (`kLa_per_d missing`), specific stirring energy and free-energy/stirrer-energy ratio (`stirrer_power_W missing`), the η uncertainty interval (joint parameter uncertainty is not propagated), and pulse 1's oxygen demand and kLa (non-positive gain-corrected consumption: the electrode rises in its only slope window).

Associations (EXPLORATORY, n = 4, no p-values because the smallest attainable two-sided p for Spearman with n = 4 is 0.083): ρ(consumption, ΔA) = −0.4; ρ(consumption, FA 20/25 °C) = +0.8; ρ(consumption, median pH) = +0.8; ρ(consumption, pulse time) = +0.2. The positive FA sign is the opposite of the inhibition hypothesis; no association is claimed. The MODEL-BASED η of pulse 1 (0.21) disagrees with its measured slope window (rising electrode), so η for pulse 1 rests on the model alone.

### Commands run

```
uv run python -m bnr_neq pulses           # 74 computed, 18 not computed; figure adjusted twice after viewing
# energy path for a measured stirrer power implemented; report notes on pulse 1 and the FA sign added
uv run ruff check --fix --unsafe-fixes .  # RUF007, ISC004 x4; behaviour-preserving
uv run pytest -q tests/test_pulses.py     # 11 passed
# gate run (background; progress.md parked in the scratchpad and restored):
uv sync                                   # Resolved 30 packages, Checked 28 packages
uv run ruff check .                       # All checks passed!
uv run pytest -q                          # 234 passed in 98.09s
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all          # second run
diff -rq <run1 tables> outputs/tables     # all 47 CSV tables identical
```

### Gate M7: PASS

- Every metric in `pulses.csv` is either computed with a non-empty `inputs` entry or carries an explicit `not computed: <reason>` status with no value (`test_every_metric_computed_with_inputs_or_not_computed`).
- MODEL-BASED, GAIN-CORRECTED and SCENARIO labels are carried in the metric names or the basis column; every association is marked EXPLORATORY.
- Stopped here as instructed; M8 not started.

## M8. Report and single-command reproduction (2026-10-18)

User decision before M8: "Accept M7 as documented and proceed to M8, including the Section 9 validation commands and the Section 11 final report."

### Implemented

- `src/bnr_neq/report.py`: `assemble_report` builds `outputs/reports/analysis_report.md` from the generated tables only (experimental facts from the configuration, verified data properties, sensor verdicts, every reference value, the main tables, all five figures, an index of every table with its row count, the stated assumptions, and the status label of each conclusion with its evidence and missing items). It also writes `outputs/tables/all_references.csv`. CLI `report`; `run-all` now runs all eight stages and no command is left unimplemented.
- A hard-coded number in the first draft of the status table ("4 to 45 mV") was replaced by the value computed from `thermo_chloride_scenarios.csv` (a new test caught it).
- Section 7 asks for one test file per module: added `tests/test_main.py` (CLI), `tests/test_pipeline.py` and `tests/test_constants.py`, plus `tests/test_report.py`. `test_constants.py` documents that 14.007/18.039 = 0.776484, so the configured 0.77649 is rounded up in the fifth decimal (relative difference 7.5e-6).

### Status labels (from `analysis_report.md`)

- GO: chemical disequilibrium (coexistence SCENARIO affinity 316.9 kJ/mol, 305.5 to 338.6 kJ/mol for DO 0.01 to 8 mg/L); smallness of the affinity steps (2.36 to 4.17 kJ/mol, 0.7 to 1.2 % of |ΔG°′|).
- GO WITH DOCUMENTED NON-BLOCKING DEBT: pH trajectory analyses (12/187 stationary windows; ordinal KL significant in 10/10 computed windows; Onsager not interpretable); relative comparison of pulse responses under the electrode gain model.
- NO-GO: absolute NH4 concentrations, rates and conversion (lab_NH4_N_mg_L, electrode_step_response, NO3_N_mg_L, NO2_N_mg_L); conductivity as solution conductivity (probe validation, chloride_mmol_L); absolute E_H (orp_reference_offset_V, chloride_mmol_L); oxygen and energy claims (DO_mg_L, kLa_per_d, stirrer_power_W, temperature_C_measured).

### Section 9 validation commands (final run, after the last code change)

`progress.md` was copied to the scratchpad before `rm -rf outputs` and restored afterwards (deviation: Section 9 deletes `outputs/`, where this hand-written log was requested).

```
$ uv sync
Resolved 30 packages in 5ms
Checked 28 packages in 1ms
$ uv run ruff check .
All checks passed!
$ uv run pytest -q
253 passed in 98.78s (0:01:38)
$ rm -rf outputs && uv run python -m bnr_neq run-all
ingest: 1 row(s) dropped; ...
validate: ... verdicts: {'nh4': 'FAIL: under-response, 17 to 37% of the dose', 'ec': 'FAIL: no response to the pulses', 'orp': 'FAIL: no chloride response at the bare Ag', 'ph': 'PASS (no contradiction)'}
segments: ... references 11/13 reproduced
thermo: ... references 9/11 reproduced
dynamics: estimator validation 25/25 required checks passed; ...
model: selected biomass_decay; consistency checks pass; 13/15 parameters identifiable; ...
pulses: 74 metrics computed, 18 not computed; ...
report: wrote .../outputs/reports/analysis_report.md
run-all: ran ingest, validate, segments, thermo, dynamics, model, pulses, report
$ uv run python -m bnr_neq run-all   # second run: identical stage messages
--- compare CSV tables run 1 vs run 2
no differences in 48 CSV tables
--- compare analysis_report.md run 1 vs run 2
analysis_report.md identical
```

(An earlier complete Section 9 run, before the three extra test files, gave the same result with 240 tests.)

### Gate M8: PASS

- `run-all` completes from a clean `outputs/` directory and regenerates every output from the raw files.
- A second run produces identical tables (48 CSV files) and an identical `analysis_report.md`; random procedures use fixed seeds.
- Reference values overall: 47 of 55 reproduced; the 8 others are explained with evidence in the M2, M3 and M4 entries.
- Work complete through M8. Not approved: the independent reviewer decides.

## Review fixes (2026-10-18)

Implementer log for the 14 review findings. Not approved; the reviewer decides.

### What changed, per finding

1. **Matched-null IAAFT gate for pH (dynamics.py).** Every IAAFT test now uses 99 surrogates (`N_SURR = 99`, rank-p resolution 0.01; data windows, tie sweep and positive control). New regime grid: reversible AR(1), φ ∈ {0.9, 0.95, 0.99, 0.995, 0.999} × target tie fraction ∈ {0.6, 0.7, 0.8, 0.85, 0.9, 0.95} (30 regimes), n = 2160, scaled analytically so that the 0.01-quantized series has the target tie fraction (`quantized_ar1`). Only realizations that pass the same ADF/KPSS filter as the data enter a regime (up to 100 of at most 4000 draws). Threshold: false positives ≤ the 95 % binomial quantile at rate 0.05 (`fp_threshold`, 9/100). A stationary pH window is matched to every regime whose stationary realizations span its tie fraction, lag-1 AC and DFA α; an index is reported only if ≥ 1 regime matches and every matched regime is validated, else `not computed: outside validated regime (...)`. Regime rows are stored in `estimator_validation.csv` (role `gate: pH windows`); required checks `matched_null_gate_pH_<index>` re-match every reported window from the stored statistics, and the report re-runs that check. Result: the persistent, highly tied regimes fail for `ordinal_kl` and `asym_norm_lag1` (e.g. R24, φ 0.995/tie 0.95: 46/100 and 51/100 false positives), and the high-tie pH windows (lag-1 AC ≈ 0.997) lie outside every regime (stationary null draws reach AC ≈ 0.99 at most), so only 3 of 12 stationary pH windows (10.28, 12.32, 18.57 d) keep reported indices.
   - Performance: one persistent process pool (`BNR_NEQ_WORKERS`, default `os.cpu_count()`); each task carries its own seed, so results do not depend on the worker count or scheduling (checked: the first version, with a new pool per batch and ADF-first filtering, and the final one gave identical validation tables). `passes_stationarity` evaluates KPSS before ADF (same accept decision as `stationarity`, tested); the suite takes about 190 s on 10 cores.
2. **Every step candidate splits windows (segments.py).** `step_candidates` takes every M2 detector candidate with |step| ≥ 4 mg/L. Located sustained jumps exclude [t10, t90] (left window ends at t10, right resumes at t50 + 0.1 d as before); unlocated or non-sustained candidates exclude their whole above-threshold run. `windows_check` tests every candidate region and every pulse. The 3.69, 4.06, 4.82, 5.02, 5.89 and 9.29 d candidates now bound windows; all 15 computed windows are clean.
3. **Non-positive gain-corrected consumption.** `segments.csv` reports `not computed: non-positive gain-corrected consumption (rising electrode signal)` with NaN rate and interval (wording as `thermo.oxygen_table`); `pulses.py` excludes those windows from both duration-weighted rates and names them in `basis`.
4. **Held/damped NH4.** Flagged, not subtracted: `held_overlap_d`, `held_fraction` and a sensitivity slope with the held bins removed (`slope_excl_held_*`) in `segments.csv`; `held_fraction` is carried to `thermo_oxygen.csv`, every rate row of `pulses.csv` and the reports. Reason: subtracting the 13-18 d and 20.7-24.7 d spans would leave 10.98-13.0 d of the window whose rate is the 5.2 mg N/L/d reference; both numbers are now reported (production 5.19 kept, held-excluded slope −0.91 vs −1.94 mg NH4+/L/d shown alongside).
5. **Profile likelihood.** `model.profile_interval` brackets the first converged Δ ≥ 3.84 on the 3-point grid and finds Δ(θ) = 3.84 with Brent's method on the profile itself (endpoint values cached, tolerance 1e-3 x half-width).
6. **Placeholders removed.** In situ affinity uses the observed NH4-N and pH medians of each slope window; E_H uses the observed ORP medians 1-3 h after each pulse; Ω_CO2 uses each row's pH; the pulse kLa margin divides by the computed required kLa; the M6 in situ σ uses the simulated state and observed pH and, when not computed, has `pass = "not computed"`.
7. **ORP stability** excludes 28.84-28.99 and 30.43-30.58 d; rolling windows never span an exclusion, each piece is detrended separately (`stability_pieces`).
8. **Onsager and flux inputs** pass ADF/KPSS: the late fluctuation window (else `not computed: non-stationary ...`, ratios NaN) and each flux input window in both coordinates (non-stationary windows dropped).
9. **Table-derived reports.** `validation.run` writes `sensor_verdicts.csv`; the analysis report reads verdicts and builds every evidence sentence from tables; the Onsager conclusion comes from `dynamics.onsager_conclusion`; `report.check_required` raises `ReportValidationError` (before writing) if any required estimator check, the matched-null gate re-check or an M6 consistency check fails. `dynamics.write_report` raises likewise.
10. **SCENARIO labels.** Every value computed with 25 °C constants carries `SCENARIO_T25C` or a `SCENARIO 25 °C` basis (ΔA, its fraction of ΔG°′, Kohlrausch predictions, required kLa, oxide E_ref and pH coefficient, air-equilibrium bicarbonate and Ω_CO2, standard free energy per pulse); figure axes and legends too.
11. **Config.** M6 split windows from `pulses.times_d` (`split_windows`); η axis from the configured dose; report assumptions from config and tables.
12. **PELT penalty.** σ_p = σ_e/(1 − φ), the AR(1) long-run scale of the residual about the 1-day rolling median (φ = lag-1 AC, σ_e = robust scale after prewhitening). Sensitivity table: all β for the production scale, and the innovation and former marginal scales at β = 3. The innovation scale alone ("variance after prewhitening") would lower the penalty and add breaks, so it is reported as sensitivity, not used.
13. **Intervals in pulses.csv.** Columns `ci95_low`, `ci95_high` (from `segments.csv`; duration-weighted bounds for weighted rates; propagated to O2 demand and kLa) and `held_fraction`.
14. **SVG determinism.** `report.pyplot()` sets `rcParams["svg.hashsalt"]`; every figure uses it.

Also: a slope window with fewer than 25 whole grid bins is now `not computed` (previously a computed window without bins would give NaN silently; found from a warning in a new synthetic test). The implementation log moved to `docs/progress.md`.

### Tests

Added (all call production functions): candidate window check (synthetic unlocated and edge cases, real 3.69/5.89/9.29 d), boundary split, non-positive consumption in segments and pulses, held span flag and sensitivity, PELT scale properties and sensitivity, quadratic profile toy (asymmetric σ, analytic ends), open-side profile, in situ σ with measurements, split windows from config, guarded thermo paths with measurements (not 1.0/7.0/0.0), kLa margin, matched-null regime table, binomial threshold, tie-fraction calibration, persistent quantized AR(1) at tie 0.94 (new seeds), pH gate and tampered-gate failure, KPSS-first decision equivalence, ORP stability exclusion, non-stationary Onsager window, flux/Onsager stationarity columns, Onsager conclusion, real-data smoke tests of `model.compute`, `pulses.compute`, `dynamics.run`, `report.assemble_report` and the run-all CLI, report raising on a failed required check and on a gate violation, determinism of two run-all trees (CSV and SVG byte-identical). The tautological `saw[::-1][::-1]` assertion was replaced by the time-reversal property: the pattern counts of x[::-1] are those of x's mirrored windows, and since mirroring is an involution the index is invariant under reversal (palindromes score zero).

Changed assertions (data facts changed by the fixes, not weakened): pulse 1 is no longer a rising window (the 9.29 d candidate split it: 8.96-9.27 and 9.39-10.57 d both fall), so the pulse-1 non-positive assertion moved to pulse 4's last window plus synthetic tests; `test_oxygen_not_computed_for_non_positive_rates` selects the windows by their status (their rate is now NaN); `test_data_windows` accepts `outside validated regime` only for pH windows; `test_split_evaluation_on_synthetic_record` passes the split windows as an argument instead of monkeypatching constants; column renames for SCENARIO labels.

### Reference values after the fixes

47 of 55 reproduced (unchanged count). The O2 demand and kLa maxima still fail but for a new reason: the 9.29 d candidate split pulse 1's window, and 8.96-9.27 d now gives 16.96 mg N/L/d gain-corrected (77.5 mg O2/L/d, kLa ≥ 9.4 1/d) instead of a non-positive rate over 8.96-10.57 d. The reference (33 mg O2/L/d) still corresponds only to the preliminary 2.59-8.82 d window, which spans pulse 1 and three jumps. Nothing was tuned.

### Validation commands (exact output; absolute paths shortened to `.../`)

```
$ uv sync
Resolved 30 packages in 4ms
Checked 28 packages in 0.76ms
$ uv run ruff check .
All checks passed!
$ uv run pytest -q
293 passed in 1836.78s (0:30:36)
$ rm -rf outputs && uv run python -m bnr_neq run-all
ingest: 1 row(s) dropped; ...
validate: ... verdicts: {'nh4': 'FAIL: under-response, 17 to 37% of the dose', 'ec': 'FAIL: no response to the pulses', 'orp': 'FAIL: no chloride response at the bare Ag', 'ph': 'PASS (no contradiction)'}
segments: ... references 11/13 reproduced
thermo: ... references 9/11 reproduced
dynamics: estimator validation 28/28 required checks passed; ...
model: selected biomass_decay; consistency checks pass; 13/15 parameters identifiable; ...
pulses: 72 metrics computed, 16 not computed; ...
report: wrote .../outputs/reports/analysis_report.md
run-all: ran ingest, validate, segments, thermo, dynamics, model, pulses, report
$ uv run python -m bnr_neq run-all   # identical stage messages
$ diff -rq <run 1>/{tables,processed,reports,figures} outputs/...
(no output for any of the four trees; exit 0)
```

The test suite now takes about 30 min, dominated by the two real `run-all` executions it contains (determinism test) and `model.compute`.

## Remediation pass 1: W1, W2, W8 (2026-10-19)

Scope: work items W1, W2 and W8 of `docs/remediation_prompt_2026-10-19.md` only; W3 to W7 not started. Not approved; the reviewer decides.

### Verification of the review claims touched by this pass

- R2 reproduced: `segments.py` computed free ammonia from the raw reading (`slope_table`, FA loop) and the linear gain correction (−slope / response fraction); `model.py` `observe` uses log R = a + b(t) log C. Whole-record refit: b0 = 1.324, a = −3.288 and C0 = 185.8 mg N/L (both not identifiable), electrode start median 36.5 mg N/L.
- R3 reproduced: `balance_alkalinity.csv` had only "not computed" for the initial and remaining alkalinity; no carbonate state anywhere. The reviewer's 14, 35, 70 mg N/L (Alk0/2) and FA ≈ 8.5 mg NH3/L at 200 mg N/L, pH 7.8, 25 °C are reproduced by the production functions.
- R9 reproduced: no vessel mass, cover or top-up key in the configuration.
- Not reproduced as worded: the reviewer cites b ≈ 1.32 with profile interval [1.181, 1.358] (from the M6 log); the current `model_identifiability.csv` (after the review fixes) gives [0.619, 1.718]. The point estimate agrees.
- Environment: `uv` recreated `.venv` on first use (the previous interpreter path no longer existed); CPython 3.14.8 instead of 3.14.7.

### W1: one NH4 observation model

- `model.py`: `invert_observation` (C = exp((log R − a − offsets)/b(t)), quasi-static: τ not inverted), `nh4_truth` (C(t) of the whole-record refit plus the envelope of the inversions under every converged profile vector with Δ ≤ 3.84; open sides detected from the non-identifiable profiles), `truth_id` (SHA-256 of the series), `check_single_observation_model` and `ObservationModelError`. `profile_likelihood(..., collect=)` records the profile vectors. Outputs: `processed/nh4_truth.parquet`, `tables/nh4_truth_summary.csv`, `tables/nh4_truth_vectors.csv`.
- `segments.py`: `truth_slope_table` (Theil-Sen on the inverted C(t) over the same M3 windows; FA from C(t); O2 demand 4.57 g O2/g N; parameter band from the profile vectors; consumption per volume scenario r = −dC/dt + C E/V). Gain-corrected columns and reading-as-TAN FA are kept, labelled by the new `gain_corrected_role` and `fa_basis` columns.
- `thermo.py`: `thermo_oxygen.csv` gets a `role` column (SENSITIVITY: linear gain model); values unchanged.
- `pulses.py`: rates, FA, η, ΔA and O2 demand now come from the inverted C(t); every row carries `nh4_observation_model` and `nh4_series`, and `check_single_observation_model` runs before the table is returned. Gain-corrected and reading-as-TAN rows are kept as SENSITIVITY rows. New: dose-anchored η (C before the pulse + dose − C before the next pulse)/dose, η and consumption per volume scenario, `alk0_required_*`. New tables `nh4_truth_windows.csv` and `nh4_truth_balance.csv` (dose-anchored TAN bound C_up = (C(t0) V(t0) + Σ m_N)/V(t), C_low = C_up − N_nit,max(t; Alk0) V0/V(t)).
- Results: C(t) = 177 mg N/L at day 1 (band 75 to 282, open upward) and 294 at day 47.9 (182 to 357, open upward). The removal implied by C(t) needs Alk0 ≈ 20.9 meq/L if nitrified, or 10.5 meq/L if all removal released 1 eq/mol N. From pulse 3 on, C(t) lies below even the Alk0 = 15 meq/L lower bound. η for pulse 3 is negative (−0.18) and is `not computed`.

### W2: proton, alkalinity and carbonate bound

- New `alkalinity.py`: carbonate constants (Plummer and Busenberg 1982; Kw Harned–Owen form), carbonate alkalinity at the observed pH with CO2 at or above air saturation, `nitrification_bound_mgN_L` (raises on negative Alk0), reverse-cumulative-minimum `bound_curve`, a self-consistent kLa_CO2 SCENARIO bound, and the bound table over Alk0 = 1 to 15 meq/L × temperature (20, 25 °C) × CO2 scenario × volume scenario, with the minimum residual TAN and FA at the end of the record. Run in the thermo stage; outputs `alkalinity_bound.csv`, `alkalinity_references.csv`, `reports/alkalinity.md`, and bound rows in `balance_alkalinity.csv`.
- Results (air-equilibrium CO2, constant volume, 25 °C): N_nit,max by the end = 6.0, 13.0, 34.0, 69.0, 104.1 mg N/L for Alk0 = 1, 2, 5, 10, 15 meq/L, against 261.86 mg N/L dosed (at most 39.7 %). Minimum residual TAN at the end is 255.9 to 157.8 mg N/L, and FA is 3.9 to 2.4 mg NH3/L (pH 7.35). The pH record alone requires Alk0 ≥ 0.643 meq/L.

### W8: vessel volume

- `config/experiment.toml`: `vessel_mass_start_g`, `vessel_mass_end_g`, `covered`, `top_up_mL` and `initial_alkalinity_meq_L` added to `[missing]`. New scenario lists `evaporation_mL_per_d = [0.0, 0.25, 0.5, 1.0]` and `kLa_CO2_per_d = [1, 10, 100]`. `config.py` validates them once measured (types, signs, masses together).
- New `volume.py`: a MEASURED trajectory from the masses (salt mass and top-up accounted for), a covered vessel as a constant-frame ASSUMPTION, or else the SCENARIO rates. The upper bound on net loss is V0 / last sensor time. `volume_scenarios.csv`. V(t) is carried through the nitrogen, chloride and alkalinity balances (new rows per scenario; the closed-frame rows are unchanged), the W2 bound, the W1 dose-anchored bound, η and the consumption rates.
- Not carried: the M6 fit itself keeps the constant 65.46 mg N/L dose (refitting per scenario costs about 6 min each); stated as a limitation.

### Tests added

`tests/test_alkalinity.py` (30 cases) and `tests/test_truth.py` (20). `conftest.synthetic_truth` provides a known truth (C = 5 × reading). Positive: 14/35/70, HCO3⁻ 0.401/0.636, FA 8.5 ± 0.1. Negative: mixed models, unlabelled sensitivity rows, raw reading feeding a derived quantity, negative Alk0 (function and config), invalid vessel keys, negative evaporation, a tampered truth id. Boundary: Alk0 at the pH floor gives 0 (function and curve), a window shorter than the minimum span gives a status, a rising C(t) is not a consumption. Integration: on the real run-all outputs, pulse 2's consumption, FA, η and O2 demand share the C(t) identifier, and O2 = 4.57 × the weighted window consumption. Existing tests changed (not weakened): the `real_inputs` fixture gains the synthetic truth; the bad-gain check now targets `o2_demand_SENSITIVITY_linear_gain_mgO2_L_d`, because `o2_demand_mgO2_L_d` no longer depends on the gain; `test_run_and_missing_inputs` also writes the grid (W2 needs pH) and first checks that its absence raises; the run-all report list gains `alkalinity`.

Out of scope this pass (belong to W3, W4): "a state rule referencing an excluded sensor raises", "negative SAE raises", "an all-artifact interval yields no state", and Figs. 1, 4 and 8.

### Commands and results

```
$ uv run ruff check
All checks passed!
$ rm -rf outputs && uv run python -m bnr_neq run-all        # exit 0, 10 min
... segments: references 11/13 reproduced; thermo: references 9/11 reproduced;
dynamics: estimator validation 28/28 required checks passed; model: selected biomass_decay; consistency checks pass;
13/15 parameters identifiable; pulses: 140 metrics computed, 20 not computed
$ uv run pytest -q
1 failed, 341 passed in 1818.16s
```

The failure (`test_one_pulse_derives_everything_from_the_same_inverted_series`) was a bug in `check_single_observation_model`: under pandas 3, `astype(str)` keeps NaN, so an empty `basis` cell read back from CSV broke the string join. Fixed (empty cells count as empty text), regression test `test_check_accepts_empty_cells_read_from_csv` added. The function only validates, so no output changes. Re-run afterwards:

```
$ uv run ruff check
All checks passed!
$ uv run pytest -q tests/test_truth.py tests/test_pulses.py tests/test_model.py   # includes a fresh run-all fixture
57 passed in 816.02s
```

The full suite was not re-run after this one-function fix.

## Remediation pass 2: re-review findings, W3 to W7, R6 to R12 (2026-10-19)

Scope: every open item of `docs/rereview_2026-10-19.md` (S1, S2, S3, W1 a-c, W3, W4, W5, W6, W7, W8 gaps, R6, R8, R10, R11, R12). Implementer log; not approved, the reviewer decides. Thresholds and rules were fixed in code before the corresponding outputs were computed; nothing was tuned to a number.

### Verification of the re-review claims

All reproduced on the code before editing: `report.py` read the GO ΔA from `thermo_affinity_steps.csv` (reading-based, 2.36 to 4.17 kJ/mol) while `pulses.csv` gave 0.60 to 1.15; `segments.csv` FA columns had no SENSITIVITY in their names; `nh4_truth_balance.csv` was contradicted for Alk0 ≤ 15 at pre_pulse3, pre_pulse4 and end; η = 1.180 and 1.154 under 1 mL/d were `computed`; `volume.py` used `top_up_mL = 0.0` under MEASURED; `alkalinity_mmol_L` was ignored by W2; `model_specs` had no alkalinity candidate; no states, energy family, Figs. 1/4/8; `pyproject.toml` placeholders. One extra defect found: the balance row "Initial alkalinity" printed `not computed` even if `alkalinity_mmol_L` were measured (fixed, see W8).

### W1 remaining defects (S1, W1 a-c)

- `thermo.py`: reading-based ΔA columns renamed `delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T*` (and `delta_A_fraction_of_dG0_SENSITIVITY_reading_as_TAN_SCENARIO_T25C`) with a `role` column; the M4 references are compared against them and say so.
- `segments.py`: `fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L_SCENARIO_T*` (label in the column name).
- `report.py`: the ΔA in the status table now comes from the primary `pulses.csv` rows (M6_log); the SENSITIVITY range is quoted beside it (`primary_delta_a`; test `test_headline_affinity_step_comes_from_the_primary_model`, which also proves the M4 table no longer feeds it).
- `model.check_single_observation_model` also checks `nh4_series` (`check_single_series`): primary rows share one C(t) id; only rows labelled MODEL-BASED may (and must) use one `M6_state:` series. The new check immediately caught `eta_uncertainty_interval` (model-state series, no label); fixed by labelling it, and the stirrer ratios, MODEL-BASED / LABORATORY-SPECIFIC.

### S2: C(t) against the W2 bound

`nh4_truth_balance.csv` gains the per-band-vector removal and Alk0 requirement; new `nh4_truth_w2_consistency.csv` gives a verdict per checkpoint; every M6-derived row of `pulses.csv` and `nh4_truth_windows.csv` carries `w2_bound_check`, and `check_w2_qualified` raises before writing if a row of a contradicted pulse is not `QUALIFIED`. Result (constant volume): the best-fit removal needs Alk0 ≥ 9.5, 13.1, 20.8, 19.0 and 20.9 meq/L at pre_pulse1, pre_pulse2, pre_pulse3, pre_pulse4 and the end; contradicted from pre_pulse3, so pulses 2 to 4 are QUALIFIED. The contradiction is resolvable within the data: the band minimum is 10.7 to 11.5 meq/L (the upper part of the open band), and a 1 eq/mol N sink (NH3 stripping) needs ≥ 10.5 meq/L. The report states the three readings (level too low, Alk0 > 15, a non-nitrifying proton-releasing sink).

### S3: η

`removal_to_dose_ratio` (and per volume scenario) reports the removal per dose, which exceeds 1 when carry-over TAN is removed too; `eta_conversion_fraction*` is computed only for 0 ≤ ratio ≤ 1, otherwise `not computed: removal exceeds the dose (ratio …)`; negative stays `not computed`. The MODEL-BASED η has the same guard. Changed values: η under 1 mL/d for pulses 2 (1.180) and 4 (1.154) are now not computed; the ratios are reported. `n_removed_mgN` and the free energy released rest on the ratio (unchanged values where η was computed).

### W3: process states (`states.py`, stage `states` between `model` and `pulses`)

Rules, windows and thresholds fixed before computing (`states_rules.csv`): 1.01 d centred Theil-Sen slopes of C(t), pH and ORP; pH needs one 0.01 quantum per window; ORP relative to its centred 5 d median, ±15 mV = RT ln10/4F (one decade of p_O2 on a Nernstian 4-electron couple). States: active nitrification (C falls, pH not rising), oxygen-limited (C flat, ORP_rel ≤ −15, ORP and pH not rising), inactive and over-supplied (C flat, ORP_rel ≥ +15 or pH rising), sensor artifact (pulse transient, NH4 step region, held NH4, ORP events, empty bins: no process state), unclassified. Features excluded by FAIL verdicts (`nh4_reading`, `ec`, `orp_absolute`) cannot be used: `validate_rules` raises. Because the states use the NH4 slope, NH4 consumption is compared between NH4-free indicator groups (ORP low/high, pH rising/not) with a moving-block bootstrap (block = max(101 bins, first lag with ACF < 1/e), 999 resamples, fixed seed, pooled and within pulse intervals); each group needs ≥ 2 blocks of bins, else `not computed` (added after the first run showed a percentile interval not covering its own estimate with 89 vs 66 bins against a 101-bin block: a validity condition of the method, independent of the direction of the result).

Results: 12 % of the record carries a process state (9 % active, 3 % inactive and over-supplied, 0 % oxygen-limited), 27 % sensor artifact (mostly held NH4 output), 61 % unclassified (814 bins of C(t) rising without an N source, 426 bins of consumption with a rising pH, windows crossing artifacts). pH rising vs not: +1.35 mg N/L/d, 95 % block bootstrap −0.53 to +2.57 (within pulse intervals); ORP low vs high: −0.75, interval not computed (89 and 66 bins). No association is claimed; R5 is answered: at the registered thresholds relative ORP and pH do not separate NH4 consumption in this record.

### W4: energy family (`energy.py`; rows in `pulses.csv`)

Per pulse: `n_removed_mgN`, `o2_required_mgO2` (4.57) and the nitritation case (3.43), `o2_supply_capacity_mgO2_SCENARIO_kLa*` = kLa (C* − DO) V Δt (DO = 0 while missing: a capacity), `o2_utilization_SCENARIO_kLa*`, `aeration_energy_demand_Wh/kWh_per_kgN_SCENARIO_SAE*` (floor), `aeration_energy_supply_kWh_per_kgN_SCENARIO_SAE*_kLa*`, `gibbs_to_aeration_energy_*`, `idle_fraction_of_stirred/classified_time` (W3), `relaxation_time_d` and its ratio to the pulse interval (R7: a non-equilibrium measure from the trajectory). New config SCENARIOS: `kLa_O2_per_d = [1, 10, 100]`, `sae_field_kgO2_per_kWh` (first [1, 2, 4], corrected to [0.3, 1.0, 2.6], see below) with `sae_source` (required; non-positive SAE raises in config and in `energy.aeration_energy_kWh`). A kLa scenario whose capacity is below the requirement is `not computed: SCENARIO infeasible` (kLa 1/d for pulses 1 and 2), so no utilization above 1 and no supply energy below the demand floor is reported. Stirrer terms are LABORATORY-SPECIFIC. SAE range and source were corrected after a verification (see "SAE verification" below): the values are now 0.3, 1.0 and 2.6 kg O2/kWh.

Hand check (pulse 1, kLa 10/d, SAE 2, before the SAE correction): capacity 10 × 8.26 × 0.1 L × 2.04 d = 16.85 mg O2; n_removed 2.490 mg N, O2 required 11.38 mg (utilization 0.675); supply energy 16.85e-6 kg / 2 = 8.43e-6 kWh per 2.490e-6 kg N = 3.38 kWh/kg N; floor 4.57/2 = 2.285 kWh/kg N; ΔG released 62.0 J / 30.3 J = 2.04. `tests/test_energy.py::test_energy_family_on_the_real_record` repeats this computation at the central SAE of the configuration.

### W5: alkalinity-limited candidates (`model.py`)

`alkalinity_limited` (v0, Alk0; r = v0 A/(K_ALK + A) C/(K_NH + C), A = Alk0 − 2 n/M_N from the model's own consumption) and `alkalinity_FA_SCENARIO_T20C/T25C` (× Haldane in FA of the state TAN), K_ALK = 0.5 meq/L fixed (ASM3), same chronological split. Selected: still biomass_decay (selection skill 0.958, test 0.748, forecast 0.862); alkalinity_limited scores 0.934, 0.817 and 0.925: better out of sample, so the data do not discriminate activity decay from alkalinity limitation (stated in the report). C(t) is unchanged (`nh4_truth_balance.csv` values identical).

### W6, W7, R6, R11, R12

- New `figures.py`: Fig. 1 (linkage with each link's basis), Fig. 4 (N_nit(Alk0) against the C(t) removal and its band; rate against FA per window; state timeline), Fig. 8 (η with band and evaporation ratios; utilization per kLa with infeasible marks; idle fraction; kWh/kg N per kLa with the SAE range and the demand floor), drawn in the report stage from the tables. Palette validated with the dataviz validator (CVD ΔE ≥ 9.2; aqua < 3:1 contrast, so every series is labelled and the report carries the tables). Each figure was rendered and inspected; Fig. 1 text overflow and Fig. 8 legend collisions were fixed after the first look.
- `analysis_report.md` restructured: 1 what the data say (nine table-driven points), 2 energy and efficiency (compact metric table and Fig. 8), 3 process state, W2/S2, volume, 4 model, 5 sensor validity, 6 status of the conclusions, 7 facts and assumptions; Supplement S1 holds the stochastic-thermodynamic estimators (R6, demoted, not removed), S2 the detailed tables. Coexistence affinity, ΔA and σ ≥ 0 are `CONSISTENCY CHECK (not a finding)`. Values are rounded (3 significant figures in tables, fixed decimals in text).
- `pyproject.toml`: real description; script `bnr-neq = "bnr_neq.__main__:main"`; `module-name = ["bnr_neq"]`. The template stub `src/bnr_non_equilibrium_2026/` (untracked) is no longer packaged; left on disk for the owner to delete.

### W8 remaining gaps and R8

- Masses measured without `top_up_mL`: the net loss stays MEASURED (it does not need the top-up), the evaporation/top-up split is `not computed: top_up_mL missing` (never 0).
- `alkalinity_mmol_L` (no sampling time) is never taken as Alk0, and `alk0_source` says so in `alkalinity.md`, the report and `balance_alkalinity.csv` (it is a lower bound on Alk0 if sampled after day 0).
- New `volume_constraints.csv`: N conservation on C(t) (CONDITIONAL) needs a net loss ≥ 0.31 mL/d (band 0 to 0.45) over pulse 3; liquid present at the end bounds it below 2.09 mL/d; the flat EC record would cap it at 0.04 mL/d, reported CONDITIONAL and not used (EC failed M2). The constant frame and 0.25 mL/d are inconsistent with C(t); 0.5 and 1 mL/d are consistent. Not done (P2 debt): the M6 fit is still not refitted per volume scenario.
- R8: a baseline window between two electrode steps shorter than 3 × the longest delayed electrode response (0.33 d, M2 lag of the pulse-response jumps) is `not computed: electrode fragment`, slope kept in `slope_fragment_diagnostic_mgNH4_L_d`. Demoted: sub-segments 2, 3, 5, 7, 9 (−21.80, +7.97, −8.49, +5.48, −0.87 mg NH4+/L/d). Window 10 (6.65 to 8.53 d, ends at pulse 1) stays.

### Changed outputs against the previous run (every change explained)

New tables: `nh4_truth_w2_consistency`, `volume_constraints`, `states_rules`, `states_fractions`, `states_timeline`, `states_comparison`. Additive columns: `w2_bound_check` (`pulses`, `nh4_truth_windows`), band columns (`nh4_truth_balance`), `slope_fragment_diagnostic_mgNH4_L_d` (`segments`), `role` (`thermo_affinity_steps`). Renamed (W1, SENSITIVITY in the name): the ΔA columns of `thermo_affinity_steps` and the FA columns of `segments`; readers (pulses, report, segments report, tests) updated. `pulses.csv` 160 → 336 rows; of the earlier rows only the two η > 1 values changed (S3). The R8 demotion removes five computed windows, which explains every other change: `nh4_truth_windows` (rows 2, 3, 5, 7, 9 only), `segments_window_check` 15 → 10, `thermo_oxygen` 16 → 11, `thermo_co2` 17 → 12, the coexistence note (15 → 10 windows), `kinetic_fits`/`kinetic_parameters` (segment 0 has fewer sub-segments), `dynamics`, `dynamics_windows` (687 → 680) and `dynamics_flux` (NH4 windows only inside computed M3 windows; 28/28 required estimator checks still pass). `model_fit` 6 → 9 rows (W5). `balance_alkalinity` +1 row (Alk0 source). Reference quantities relabelled (`all_references`, `segments_references`, `thermo_references`); values and pass counts unchanged (11/13, 9/11).

### Tests

New: `tests/test_states.py` (excluded sensor raises, failed pH excludes every rule, synthetic rule table, all-artifact interval yields no state, slope significance and quantum, windows touching an artifact, block bootstrap recovers a shift, blocks widen the interval of an AR(1) with φ = 0.98 by > 3x, determinism, small-group guard, overlap raises, within-stratum comparison removes a between-stratum trend, real-record checks including held spans never carrying a state and the idle fraction in `pulses.csv`), `tests/test_energy.py` (hand values, negative or zero SAE raises in the function and the config, SAE without source raises, DO above saturation, real-record hand recomputation, every SAE/kLa metric labelled SCENARIO, no utilization above 1, infeasible kLa withholds the supply energy), `tests/test_rereview_fixes.py` (series mixing, MODEL-BASED series rules, SENSITIVITY in column names, S2 detection/qualification/raise on the synthetic and real tables, η only on [0, 1] on the real table, top-up status, N-conservation hand value and ordering check, EC bound, older alkalinity key, R8 rule and the real fragments, W5 specs and stall at Alk0 M_N/2, W5 rows in `model_fit.csv`).

Existing tests changed because the remediation prompt changes the requirement (none weakened): `test_main` (stage list includes `states` and orders model < states < pulses; states report expected), `test_report` (status rows relabelled as consistency checks, new W2/S2/W8 rows, figure list of 8, rounded ORP evidence `4 to 45 mV`, fixture extended with the tables the narrative reads), `test_segments::test_scenario_labels` (renamed FA columns), `test_truth` (η: computed values must lie in [0, 1]; non-computed reasons may be "negative" or "removal exceeds the dose").

### Commands and results

```
$ rm -rf outputs && uv run python -m bnr_neq run-all        # exit 0, 23 min 55 s
validate: verdicts nh4/ec/orp FAIL, ph PASS; segments: references 11/13; thermo: references 9/11;
dynamics: estimator validation 28/28 required checks passed; model: selected biomass_decay; consistency checks pass;
13/15 parameters identifiable; states: 12 % of the record classified; idle 3 % of stirred time;
pulses: 265 metrics computed, 71 not computed; report: wrote outputs/reports/analysis_report.md
$ uv run ruff check
All checks passed!
$ uv run pytest -q
381 passed in 3625.64s (1:00:25)
```

The suite includes two fresh `run-all` executions with byte-identical CSV and SVG outputs (determinism test) on the final code. After the `outputs/` run above, only `report.py` text (mechanism sentence, `n.c.` cells) and an unreachable branch of `volume.n_conservation_min_loss_mL_d` changed; the report stage was re-run (`python -m bnr_neq report`), and the volume change does not alter any value. The suite now takes about 60 min (three more model candidates in every `model.compute` and run-all).

### SAE verification (2026-10-19, after the commands above)

The first range, 1 to 4 kg O2/kWh attributed to Metcalf & Eddy (2014, ch. 5), had not been checked against the source. Verification against two documents that could be read in full:

- Stenstrom M.K. and Rosso D. (2010), *Aeration*, lecture notes, UCLA and UC Irvine, slide "Efficiency" (p. 56). Clean-water SAE: high-speed surface aerators 0.9 to 1.3, low-speed surface 1.5 to 2.1, coarse bubble 0.6 to 1.5, turbines or jets 1.2 to 1.8, fine pore 3.6 to 4.8 kg O2/kWh. Process-water AE at 2 mg/L DO: 0.4 to 0.8, 0.7 to 1.5, 0.3 to 0.9, 0.4 to 0.8 and 0.7 to 2.6 kg O2/kWh respectively.
- Moore L.W. (2021), *Biological Wastewater Treatment Training Series no. 11: Aeration Systems*, Oak Ridge National Laboratory (references Metcalf & Eddy, 4th ed.), slides 41 and 42. Field transfer rates of 1.0 to 2.3 lb O2/(hp h), i.e. 0.61 to 1.40 kg O2/kWh, at α = 0.84, β = 0.92, DO = 2 mg/L.

Finding: the old upper value of 4 was a clean-water fine-pore SAE, not a field value, and the textbook citation was not verified. Correction: `sae_field_kgO2_per_kWh = [0.3, 1.0, 2.6]`, which is the extremes of the process-water columns plus a central value, cited to the two documents above. The config comment and the report state that these are process-water AE values, although the key keeps the prompt's name "SAE". Fig. 1 now reads the range from the configuration instead of a hard-coded string. `tests/test_energy.py` derives the central value from the configuration.

Effect: only SAE-dependent rows of `pulses.csv` and the report change. Metric names now carry SAE0.3, SAE1 and SAE2.6. The demand floor becomes 1.76 to 15.2 kWh/kg N (formerly 1.14 to 4.57). At kLa 10/d and SAE 1, the supply-based energy is 6.8, 21.9, n.c. and 99.7 kWh/kg N for pulses 1 to 4. Hand check at pulse 1, kLa 10/d: capacity 16.85 mg O2, N removed 2.490 mg. Supply-based energy is 22.56, 6.77 and 2.60 kWh/kg N at SAE 0.3, 1.0 and 2.6, and the floor is 4.57/SAE. Both match the table. The pulses and report stages were re-run; the fast energy, report and config tests pass (33 passed). Full validation on the final code, after the SAE correction:

```
$ uv run ruff check
All checks passed!
$ uv run pytest -q
381 passed in 4580.76s (1:16:20)
```

This includes the two fresh `run-all` executions (byte-identical CSV and SVG) and the real-record energy test at the central SAE.

## Remediation pass 3: model convergence, idle fraction and provenance hardening (2026-10-20)

Scope: `docs/remediation_prompt_2026-10-19_pass3.md` (N1, N2 and X3 to X8), from commit `954f87a` on `remediation/rereview-2026-10-19`. Implementer log; not approved, the reviewer decides. The five judgment calls accepted in §9 of `docs/rereview_2026-10-19_pass2.md` are not revisited.

### Verification of the cited claims (before any edit)

Checked on the committed code and outputs of `954f87a` (scratch scripts in the session scratchpad; the baseline run-all and pytest of a clean checkout are reported under "Commands and results").

- **N1.** `model.py:384` returns `bool(sol.success)`; `model.py:431` selects by minimum `sel_rmse_log` among all non-null candidates, converged or not. `model_fit.csv`: `cal_converged = False` for `alkalinity_limited`, `alkalinity_FA_SCENARIO_T20C` and `alkalinity_FA_SCENARIO_T25C`; True for the other six. Skills as quoted (alkalinity_limited 0.934 / 0.817 / 0.925; biomass_decay 0.958 / 0.748 / 0.862). `initial_guess` starts Alk0 at 10 meq/L (`model.py:311`). Neither the report nor this log mentioned the non-convergence. Confirmed.
- **N2.** `states.py:153` (`theilslopes`, independent-point interval), `sign_class` at `states.py:176`, idle fraction over all bins at `states.py:277`: confirmed. On the committed `states_bins.parquet`: 2,031 bins with a C(t) slope, 379 flat, 838 falling, 814 rising; classified fraction 0.118; median lag-1 autocorrelation of the Theil-Sen residuals of C(t) inside the 101-bin windows 0.926, median n_eff = 3.9 (the review's 3.7 uses ρ rounded to 0.93), median √(n/n_eff) = 5.1. Confirmed.
- **X5.** Pulse 2 `eta_conversion_fraction` band 0.618 to 1.069 (`param_band_high` = 1.068849). Every baseline window (sub-segments 1 to 10) of `nh4_truth_windows.csv` has an empty `w2_bound_check`, not only window 1. Confirmed, slightly wider than cited.
- **X6.** Sub-segment 14: fit window 8.958 to 9.269 d (0.31 d), starts 0.100 d after the delayed pulse-1 response jump at 8.858 d (lag 7.88 h), M6 τ = 0.262 d (3τ = 0.785 d), consumption 18.29 mg N/L/d, the highest computed window. Rule at `segments.py:188-191`. Confirmed. Also: the M4 reference "O2 demand max (GAIN-CORRECTED)" (33 ± 3.3) currently fails with 77.5 from window 14; without window 14 the maximum is 23.7, which still fails, so the R8 decision below cannot change a reference pass count.
- **X3, X4, X7.** `DERIVED_PREFIXES` (`model.py:52`) lacks `aeration_energy_`, `n_removed`, `removal_to_dose`, `consumption_`, `idle_`, `relaxation_`, `gibbs_`; `is_primary` excludes empty models (`model.py:86`); `text.str.contains("MODEL-BASED|MODEL_BASED")` (`model.py:102`) matches "not MODEL-BASED"; the checker runs only on `pulses.csv` (`pulses.py:578`). Gain column names in `segments.py:260-262` and `thermo.py:498-523` as cited; `pulses.csv` additionally has the metrics `consumption_gain_corrected_first_window_mgN_L_d` and `consumption_gain_corrected_weighted_mgN_L_d` without SENSITIVITY in the name (renamed under X4 as well), and `thermo_oxygen.csv` has `kla_measured_over_required` computed from the gain-corrected kLa (renamed too). Report §4 renders `model_checks.csv` and `model_identifiability.csv` unrounded because their mixed text/number columns are read back as text. MODEL-BASED against dose-anchored η: 0.205 / 0.380, 0.569 / 0.809, 0.070 / n.c., 0.035 / 0.197. Confirmed.

### Pre-registered decisions (fixed before any output of this pass was computed)

- **X2 (a), slope intervals.** In `states.rolling_theil_sen` every window's Theil-Sen 95 % half-widths (below and above the slope) are multiplied by √(n/n_eff), n_eff = n (1 − ρ)/(1 + ρ) (`model.effective_n`'s formula), ρ = lag-1 autocorrelation of the window's residuals about its own Theil-Sen line, clipped to [0, 0.999] (ρ := 0 when the residuals have no variance, e.g. a pH window without a quantum change). n = 101 bins. Applied to C(t), pH and ORP alike. Argument: for an AR(1) error of coefficient ρ, the variance of a trend estimate is inflated by about (1 + ρ)/(1 − ρ) against independent points; this is the same correction the profile likelihood already uses. Nothing else in the classifier changes: the 15 mV threshold, the 1.01 d window, the pH quantum, the block-length rule and the ≥ 2 blocks rule stay.
- **X2 (c, d), idle fraction and coverage.** Per interval: `idle_fraction_lower_bound` = N_inactive / N_all; `idle_fraction_upper_bound` = (N_inactive + N_unclassified + N_artifact) / N_all (every bin without a positive non-idle state could have been idle); `idle_fraction_of_classified_time` = N_inactive / N_classified, CONDITIONAL on the classified bins being representative, reported with `classified_fraction`. Coverage threshold: `classified_fraction` ≥ 0.5, else the conditional value is `not computed: classified fraction … below 0.5`. Reason: a conditional fraction stands for the interval only if the classified bins are the majority of it; below one half the bins without a state, whose mix is unknown, outweigh those with one, and the conditional value is no closer to the interval's idle fraction than the bounds. The bounds are always reported. Fixed before the states were recomputed.
- **X6, scope of R8.** No physical argument was found that a short window between two electrode steps becomes trustworthy because a pulse governs it: a delayed pulse-response jump is a step of the electrode reading like an offset jump (the pulse itself happened hours earlier), and a window starting 0.10 d after it lies inside 3τ of the electrode's own transient by the code's `INVERSION_NOTE`. The rule is therefore extended symmetrically: a computed window whose two boundaries are both electrode steps (offset jumps, pulse-response jumps or step candidates; a pulse or the record start/end is not an electrode step) and whose fit span is shorter than 3 × the longest delayed electrode response is an electrode fragment, whatever its governing pulse. Chosen from the boundary structure of `segments.csv` only. Pulse 1's consumption is reported with and without the demoted window(s) as a labelled case (the word SENSITIVITY is reserved by the W1 checker for an alternative NH4 observation model, so the label is "R8 CASE").
- **X2, synthetic validation of the interval (before any real-data output of this pass).** The production `rolling_theil_sen` on 600 synthetic windows of 101 points, AR(1) errors, fixed seed: coverage of the true slope by the unscaled / scaled (ρ̂ from the window, as registered) / scaled with the true ρ interval = 0.74 / 0.94 / 0.96 at φ = 0.5, 0.34 / 0.85 / 0.95 at φ = 0.9, 0.24 / 0.77 / 0.95 at φ = 0.95; mean ρ̂ = 0.47, 0.84, 0.88. The √(n/n_eff) formula is right (0.95 with the true ρ); the shortfall comes from the small-sample downward bias of the lag-1 estimate on 101 detrended points. The registered estimator is kept as the remediation prompt specifies it (no bias correction invented after this check; a first-order correction ρ̂ + (2 + 4ρ̂)/n would only reach 0.90 at φ = 0.9). Consequence, stated with the results: at the record's median ρ̂ ≈ 0.93 the scaled interval still under-covers (roughly 0.8 instead of 0.95), so a few truly flat windows can still be called rising or falling; the upper idle bound counts bins classified as active nitrification as not idle and is therefore conditional on that classification.
- **X1 (a), diagnosis of the W5 calibrations** (the calibration of `split_evaluation` replicated with the full `least_squares` result; scratch scripts, committed outputs of `954f87a`). All three alkalinity candidates stop with status 0, "the maximum number of function evaluations is exceeded", at exactly scipy's default limit 100 × n_free (1,300 for `alkalinity_limited`, 1,500 for each `alkalinity_FA_*`); no parameter is at a bound; first-order optimality 0.17 to 0.56 (biomass_decay: 0.0016); fitted Alk0 22.7, 22.6 and 22.5 meq/L, i.e. Alk0 left its start of 10 meq/L and was not held below the S2 requirement (about 20.8 meq/L) by the start. `biomass_decay` converges (ftol) after 81 evaluations with C0 = 299.9994 mg N/L, on the 300 mg N/L cap, which pins the C0-a-b0 direction that the profile likelihood finds non-identifiable; the alkalinity fits sit inside that valley (C0 85 to 115) and move slowly along it. Further runs of `alkalinity_limited`: from the default start, 500 evaluations reach RSS 29.330 against 29.325 after 1,300, with C0 and b0 still drifting; from Alk0 = 20 (other starts unchanged) it converges (ftol) after 3,620 evaluations at RSS 29.3226, Alk0 22.62, C0 119.5, selection RMSE 0.0888; warm starts with the observation parameters of the converged `biomass_decay` calibration (C0 near 300, b0 1.30) and Alk0 = 10, 20, 40 are still unconverged after 4,000 evaluations at RSS 29.360 (C0 ≈ 183, Alk0 22.1), and with Alk0 = 5 in a worse region (RSS 31.74, a = −4.93 near its bound, Alk0 12.9). Root cause: an evaluation budget too small for a slow, ill-conditioned valley, not the Alk0 start; a warm start from the decay model's observation parameters is worse than the default observation start; at least one secondary basin exists.
- **X1 (b), convergence policy (fixed before any refit; identical for every candidate, `model.CONVERGENCE_POLICY`).** (1) Starts: four per candidate, the candidate's default kinetic start vector multiplied by 0.5, 1, 2 and 4 on the natural scale (for Alk0: 5, 10, 20 and 40 meq/L, the grid of the remediation prompt; v0 2.5 to 20 mg N/L/d; kd, Ks, Ki likewise), observation parameters at their default start; the null model has no kinetic parameter, so its four starts coincide and it is fitted once. (2) Budget: at most 8,000 function evaluations per start (more than twice the 3,620 the converging diagnostic run needed; about 6 x scipy's default), the same for every candidate and start, scipy's default tolerances (ftol = xtol = gtol = 1e-8), success = status > 0. (3) Result: the converged start with the lowest calibration RSS (ties: the earlier start); if no start converged, the lowest-RSS start, labelled `cal_converged = False`, `cal_status = 0`, and excluded from selection (`selection_status = not computed: calibration did not converge (status …)`), its skills kept and labelled. (4) The whole-record refit of the selected model uses the same budget. The starts run in parallel processes; results do not depend on the number of workers. The warm start from biomass_decay was rejected on the diagnosis above (it lands at the far end of the valley), before any refit under the policy.

### What changed, per work item

- **X1 (N1), `model.py`.** `fit` records the `least_squares` status and evaluation count (`Fit.status`, `Fit.nfev`, `Fit.start`) and simulates only up to the end of the fitted window (the state model and the lag filter are causal: residuals inside the window are bit-identical, checked for all nine specifications; half the cost of a calibration). `calibration_starts`, `calibrate_all` and `best_start` implement the policy above for every candidate, with the (candidate, start) pairs in a process pool (`BNR_NEQ_WORKERS`; results independent of it). `select_converged` selects by minimum selection RMSE among converged candidates only, labels every row (`selection_status`: selected / candidate (converged, not selected) / reference / `not computed: calibration did not converge (status …)`), and raises `SelectionError` if none converged. `alk0_profile` adds `alk0_fitted_meq_L`, `alk0_ci95_low/high_meq_L` and `alk0_profile_status` (95 % n_eff-scaled profile interval on the calibration window, or `not identifiable: …`) for the alkalinity candidates (`profile_likelihood(params=...)` restricts a profile to named parameters and keeps parameters fixed that the fit kept fixed). New `model_fit.csv` columns: `cal_status`, `cal_nfev`, `cal_start`, `selection_status`, `alk0_*`, `whole_record_refit_converged`. The whole-record refit of the selected model has the policy budget. Report §1 point 7, §4 (convergence columns, an Alk0-against-S2 table `alk0_comparison`) and the status table state the converged comparison.
- **X2 (N2), `states.py`, `pulses.py`, `report.py`, `figures.py`.** `rolling_theil_sen` widens each half-width by `interval_scale(n, lag1_rho(residuals))` and returns ρ (`*_slope_resid_lag1` in `states_bins.parquet`); the method is in the module docstring and `states_rules.csv`. `state_fractions` reports `n_classified`, `n_inactive_oversupplied`, `classified_fraction`, `idle_fraction_lower_bound`, `idle_fraction_upper_bound`, `bounds_status`, and `idle_fraction_of_classified_time` with `status` (`computed (CONDITIONAL …)` or `not computed: classified fraction … below 0.5 …`). `idle_fraction_of_stirred_time` is gone (it equalled the lower bound). `pulses.csv`: `idle_fraction_LOWER_BOUND`, `idle_fraction_UPPER_BOUND`, `idle_fraction_of_classified_time_CONDITIONAL`, each with the coverage in its basis, and `classified_fraction`. Report §1 point 8, the §2 energy table (bounds with coverage), §3, Fig. 1 and Fig. 8c (bounds as a bar, coverage printed, conditional value where computed) and a separate status row "Idle fraction of stirred time (W3 into W4)": BOUNDS ONLY; NO-GO for a single value. The W3 row no longer quotes an idle fraction.
- **X3, `model.py` and the writers.** `check_single_observation_model` classifies positively: every row declares `M6_log`, a `_SENSITIVITY` model, `electrode_reading` or `NH4-independent` (`model.NH4_INDEPENDENT`); an empty or unknown model raises; an NH4-independent row may not carry a series. Raw-reading rows are checked against an allow-list (`RAW_MEASURED_METRICS` for long tables, `RAW_MEASURED_COLUMNS` plus the table's declared NH4-independent columns for wide tables) instead of `DERIVED_PREFIXES` (removed). Metric names of `_SENSITIVITY` rows must contain SENSITIVITY. The MODEL-BASED label is matched as a token not preceded by "not" (`MODEL_BASED_TOKEN`). The checker now runs on `pulses.csv` and `nh4_truth_windows.csv` (after qualification, before writing), on `states_fractions`, `states_timeline`, `states_comparison` and `states_bins` (new `nh4_observation_model`, `nh4_series` columns), on `segments.csv` (new columns; rows = `electrode_reading`) and on `thermo_oxygen.csv` (new columns; rows = `linear_gain_SENSITIVITY`). In `pulses.csv` the time, pH, O2-capacity, standard free-energy and stirrer-energy rows are `NH4-independent`; the not-computed stirrer ratios are MODEL-BASED like their computed form.
- **X4.** `segments.csv`: `consumption_SENSITIVITY_linear_gain_mgN_L_d` and its `_ci95_low/high_` columns (`segments.GAIN_COL/LO/HI`). `thermo_oxygen.csv`: `consumption_SENSITIVITY_linear_gain_mgN_L_d`, `o2_demand_SENSITIVITY_linear_gain_mgO2_L_d`, `kla_required_SENSITIVITY_linear_gain_per_d_SCENARIO_T25C`, `kla_measured_over_required_SENSITIVITY_linear_gain`. `pulses.csv`: `consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d`, `consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d`. The M3/M4 reference labels say "SENSITIVITY: linear gain". Readers updated: thermo, pulses, report, segments report, tests.
- **X5, `pulses.py`.** `eta_band_on_unit_interval`: an η row carries its band clipped to [0, 1] with the unclipped band named in the basis; `removal_to_dose_ratio` keeps the full band. `qualify`/`check_w2_qualified` treat rows without a governing pulse as the interval before pulse 1, labelled from checkpoint pre_pulse1 (and checked if that checkpoint were contradicted).
- **X6, `segments.py`, `pulses.py`.** `electrode_fragment` as registered above (`ELECTRODE_BOUNDARY_PREFIXES = ("jump (", "candidate (")`, no governing-pulse condition). `nh4_truth_windows.csv` gains `consumption_fragment_diagnostic_mgN_L_d` (C(t) slope of a fragment, diagnostic). `pulses.csv` gains `consumption_weighted_mgN_L_d_R8_CASE_incl_fragments` per pulse ("R8 CASE (not the primary rate)"), computed for pulse 1 only.
- **X7.** (a) `report.numeric_cells` turns number-valued text cells of mixed columns into floats, so §4's checks and identifiability tables are rounded. (b) §1 point 2 says the band is open upward over the whole series when `band_high_open`. (c) §1 point 5 says the accumulation rests on the QUALIFIED C(t) of the contradicted pulses. (d) `relaxation_time_d` → `dose_removal_time_d` (and `dose_removal_time_over_pulse_interval`), "not a fitted relaxation constant". (e) new MODEL-BASED row `eta_gap_inversion_misfit_MODEL_BASED` = [(C_inv − C_model) before the pulse − (C_inv − C_model) before the next] / dose, which equals the dose-anchored ratio minus the MODEL-BASED η up to the 0.1 d medians (the model state obeys its own mass balance); §1 point 5 discusses it. (f) metric suffix `SAE` → `AE`; the supply-based basis, report §2 and the config comment state that a DO = 0 capacity is divided by an AE defined at DO = 2 mg/L, and quote the source's caveat that its table is not a design guide ("should not be used as a general guide line and not for design", Stenstrom and Rosso 2010, p. 55; quoted as read by the reviewer, not re-read here). Fig. 8's title notes are now derived from the table (the hard-coded "P3: removal negative" and "P2-P4 QUALIFIED" are gone).

### Tests (X8)

New `tests/test_pass3_fixes.py`: every candidate converges under the policy on a synthetic record and Alk0 is recovered within 1 % (alkalinity_limited and alkalinity_FA on their own synthetic records); the policy starts are the same rule for every candidate (Alk0 5, 10, 20, 40); an unconverged candidate cannot be selected and carries the not-computed status, its skills kept; no converged candidate raises; the checker raises on an empty model (three metrics), on eleven derived metrics labelled `electrode_reading` (including `aeration_energy_*`, `n_removed_mgN`, `removal_to_dose_ratio`, `consumption_weighted_mgN_L_d`), on an unknown model, on a series on an NH4-independent row, on a reading-based metric or wide-table column without SENSITIVITY (an unlabelled gain-corrected column), and on a wide reading table that does not declare its NH4-independent columns; "not MODEL-BASED" does not raise; the η band is clipped; baseline windows are labelled from pre_pulse1 and checked; the R8 fragment diagnostic; and on the real outputs: the selected model converged, the fitted Alk0 is reported, every checked table passes the hardened checker, the η bands lie in [0, 1], the baseline labels, the R8 CASE and the renamed metrics. `tests/test_states.py`: AR(1) φ = 0.9 coverage (true ρ: 0.92 to 0.98; estimated ρ ≥ 0.80; unscaled ≤ 0.50), the ρ cap, hand counts of the bounds, the coverage gate (below, at and on 0.5), and the real idle rows with their coverage. `tests/test_report.py`: the idle row is not GO and the model row is QUALIFIED when a candidate did not converge, with coverage and convergence in the evidence.

Existing tests changed, each because pass 3 changes the requirement it encodes (none weakened; each keeps its original assertion strength):
- `test_truth.test_dose_anchored_eta_and_affinity`: pinned. The synthetic truth (C = 5 x reading) must give a pulse-2 ratio of 1.96 ("removal exceeds the dose" branch) and a pulse-1 ratio of −0.77 (negative branch); both branches are asserted, no `if`.
- `test_truth.test_check_accepts_empty_cells_read_from_csv` and `test_labelled_sensitivity_and_measured_rows_pass`: the fact row declares `NH4-independent` (X3 (a): an empty model now raises) and the gain metric is `consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d` (X4); the NaN-label handling they test is unchanged.
- `test_rereview_fixes.test_reading_based_columns_say_sensitivity`: extended to the gain columns of `segments.csv` and `thermo_oxygen.csv`.
- `test_rereview_fixes.test_electrode_fragment_rule` and `test_real_fragments_are_demoted`: the two assertions that a pulse-response-bounded or pulse-governed window is not a fragment are inverted by the X6 decision; new assertions that a pulse is not an electrode step, that window 14 is a fragment and window 15 is not; the baseline fragments 2, 3, 5, 7 stay demoted.
- `test_states.test_all_artifact_interval_yields_no_state`: the conditional value is NaN and `not computed` as before; the bounds are now reported (0 and 1, X2 (c)). `test_states.test_states_on_the_real_record`: the pulse-table idle row is `idle_fraction_LOWER_BOUND` = `idle_fraction_lower_bound`.
- Mechanical renames only (X4, X7 (f)): `test_pulses`, `test_segments`, `test_thermo_m4` (also the two failing M4 reference labels), `test_energy`, `test_report` fixture; `test_truth` fixture gains a model-state column for the X7 (e) row.

### Added after a first full run (disclosed)

- **Alk0 profile validity (X1 d), added after seeing the first run's `model_fit.csv`.** The first full run under the policy reported for `alkalinity_limited` "not identifiable … one-sided 95 % lower limit 80.0 meq/L", equal to its estimate. `model_calibration_starts.csv` shows why that limit is invalid: two other converged starts (Alk0 starts 5 and 10) sit at Alk0 = 22.6 meq/L with Δ = n_eff ln(RSS/RSS_min) = 3.0 < 3.84 (n_eff = 74.9, ρ = 0.94, recomputed by hand on the calibration residuals), so the 95 % region contains both minima and is not one interval; the warm-started profile follows one minimum only. Rule added (`model.alk0_profile`, `MULTIMODAL_LOG_TOL` = Alk0 more than 5 % apart): if another converged start of the policy lies below Δ = 3.84, no limit is reported and the status names both minima. It only withholds a limit; it changes no selection, no fitted value and no other output. Test: `test_alk0_interval_is_withheld_when_another_minimum_lies_inside_the_threshold`. The "≥ 2 blocks" rule of pass 2 is the same kind of validity condition.
- `report.numeric_cells` and two NaN-safe string loops (`figures.fig8`, `report.summary_lines`): pandas 3 reads text columns as the new `str` dtype, whose `astype(str)` keeps NaN as a float, and whose dtype is not `object`. Found by the first full run (the report stage raised) and by rendering the report; no value changes.
- `tests/test_pulses.py` fixture `_state` gains the model-state column that the X7 (e) row reads (the first subset run failed nine tests with a KeyError on it; no assertion changed).
- The R8 CASE row first had no interval; `test_pulses.test_rates_carry_segment_intervals_and_held_fraction` (every computed consumption rate carries its interval and held fraction) failed on it. It now carries the duration-weighted Theil-Sen bounds and held fraction like the primary rate; `nh4_truth_windows.csv` keeps the fragment's interval in `consumption_fragment_diagnostic_ci95_low/high_mgN_L_d`. The test is unchanged.

### Results of the converged comparison (X1; `model_fit.csv`, `model_calibration_starts.csv`)

| Candidate | Converged (4 starts) | Chosen start | Selection / test / forecast skill | Alk0 (meq/L) |
| --- | --- | --- | --- | --- |
| biomass_decay | yes (391 evaluations in total) | x0.5 | 0.958 / 0.748 / 0.863 (selected) | — |
| haldane_FA T20C / T25C | yes | x1 | 0.839 / 0.499 / 0.730; 0.841 / 0.500 / 0.731 | — |
| combined T20C / T25C | yes | x4 / x0.5 | 0.528 / −0.375 / −0.063; 0.524 / −0.383 / −0.070 | — |
| alkalinity_limited | yes (4,884) | x4 (Alk0 start 40) | −0.908 / 0.227 / 0.443 | 80.0; not identifiable (second minimum at 22.6, Δ = 3.0) |
| alkalinity_FA T20C / T25C | **no**: all four starts at the 8,000 limit (status 0) | — | 0.797 / 0.649 / 0.807; 0.798 / 0.651 / 0.808 (not comparable) | 23.6 (not an estimate) |

- The selected model is unchanged: `biomass_decay`. Its four starts reach the same optimum (calibration RMSE 0.1137772 from x0.5, x1 and x2; x4 stops in a worse local minimum, 0.1249), and the x0.5 start is lower than the default by 5e-9 in RMSE, so the calibration moves in the 7th digit and the whole-record refit starts from a slightly different point.
- `alkalinity_limited` converges in three of four starts. Starts x0.5 and x1 converge at Alk0 = 22.6 (calibration RMSE 0.11149, selection skill 0.934: the pass-2 value, reached after 990 and 3,637 evaluations instead of being cut at 1,300); start x2 at Alk0 = 3.1 (RMSE 0.1252); start x4 at Alk0 = 80.0 with the lowest RSS (RMSE 0.10928), i.e. no alkalinity limitation inside the calibration window, which extrapolates an undiminished activity into the selection window (skill −0.91). By the registered rule (lowest-RSS converged start) the calibration is the x4 one. The 22.6 minimum lies at Δ = 3.0 < 3.84, so Alk0 is not identifiable from the calibration window: the data allow both an alkalinity that runs out near 22.6 meq/L and none that limits at all.
- The two alkalinity x FA candidates did not converge within 4 x 8,000 evaluations from any start (all near Alk0 23.6, calibration RMSE 0.1109, the same slow valley as in the diagnosis) and are excluded from selection, with their skills labelled as not comparable.
- Consequence for the report (§1 point 7, §4, status table): the pass-2 statement that alkalinity limitation is better out of sample than decay (test 0.817 against 0.748) rested on an unconverged fit and is withdrawn. Among converged candidates, decay is selected; single alkalinity limitation does not compete under the registered rule; the alkalinity x FA variants remain undecided (not converged). Fitted Alk0 against the S2 requirement (20.8 meq/L at pre_pulse3, 20.9 at the end): "compatible, not a test" for alkalinity_limited (both minima, 22.6 and 80, are above it, but Alk0 is not identifiable); not comparable for the unconverged candidates. The status row reads "QUALIFIED: unconverged candidates excluded".

### Changed outputs against `954f87a` (every change, by work item)

The outputs of `954f87a` were regenerated byte-identically first (tables, reports, figures; see Commands). Against them:

- **X1, the convergence policy, and the C(t) it feeds.** `model_fit.csv`: new columns (above); the six converged-before candidates differ only in the 6th to 7th digit except biomass_decay's chosen start; the three alkalinity rows as described above. New `model_calibration_starts.csv` (33 rows). Because the selected model's calibration and hence its whole-record refit start moved by ~1e-7, the refit converges to a point that differs in the 4th to 5th digit (v0 19.9998 → 20.0023, b0 1.32437 → 1.32468; `model_identifiability.csv`, `nh4_truth_vectors.csv`, `model_checks.csv` values, `model_response_fractions.csv`), and the inverted C(t) changes by a factor 1.00015 to 1.00054 (median 1.0002); new truth id `C_inv:c29e26a0fb89`. This explains every small change in `nh4_truth_summary`, `nh4_truth_balance`, `nh4_truth_w2_consistency` (Alk0 requirements 20.8255 → 20.8253 etc.; same verdicts), `nh4_truth_windows` (median C, FA), `volume_constraints` (0.31352 → 0.31358 mL/d), `states_*` C(t) quantities, and the 4th-digit changes in `pulses.csv` (e.g. ΔA pulse 1 1.151 → 1.150, η pulse 2 0.8093 → 0.8092, MODEL-BASED η pulse 2 0.5687 → 0.5686, FA pulse 1 1.109 → 1.104, n_removed pulse 4 1.287 → 1.286). `model_identifiability.csv` gains `profile_crossing_low/high`.
- **X2.** States with the autocorrelation-scaled intervals (median lag-1 ρ of the window residuals: C(t) 0.926, pH 0.772, ORP 0.985): of the 2,031 bins with a C(t) slope, 1,472 are flat (was 379), 332 falling (838), 227 rising (814). States: active nitrification 262 bins (was 412), inactive and over-supplied 252 (141), oxygen-limited 63 (0), sensor artifact 1,254 (unchanged), unclassified 2,861 (2,885); classified fraction 0.123 (0.118). The idle fraction, whole record: lower bound 0.054, upper bound 0.931 (was a single 0.030); per pulse interval, lower to upper with coverage: pulse 1 0.00 to 1.00 (0.00), pulse 2 0.027 to 0.907 (0.120), pulse 3 0.069 to 0.860 (0.209), pulse 4 0.103 to 0.956 (0.147); the CONDITIONAL value is not computed anywhere (coverage < 0.5 in every interval). `states_comparison.csv`: pH rising vs not (within pulse intervals) +1.01 mg N/L/d, 95 % block bootstrap −0.82 to +2.61 (was +1.35, −0.53 to +2.57); ORP low vs high unchanged (−0.75, not computed: 89 and 66 bins); the descriptive (circular) active vs inactive comparison is now computed (+3.57, 1.02 to 5.83; 262 and 252 bins), the other two descriptive comparisons not (63 oxygen-limited bins < 2 blocks). `states_timeline.csv` 109 → 113 runs. `states_rules.csv` thresholds text. `pulses.csv`: `idle_fraction_of_stirred_time` and `idle_fraction_of_classified_time` replaced by `idle_fraction_LOWER_BOUND`, `idle_fraction_UPPER_BOUND`, `idle_fraction_of_classified_time_CONDITIONAL` and `classified_fraction` (+8 rows).
- **X3.** New provenance columns in `segments.csv`, `thermo_oxygen.csv`, `states_fractions/timeline/comparison.csv`; `pulses.csv` `nh4_observation_model` is never empty (28 rows `NH4-independent`); no value changes.
- **X4.** Column and metric renames as listed above; values unchanged except where X6 removes window 14 (`thermo_oxygen.csv` 11 → 10 rows).
- **X5.** η pulse 2 band 0.618 to 1.000 (was 1.069; `removal_to_dose_ratio` keeps 1.069). `w2_bound_check` of the ten baseline windows now "pre_pulse1: consistent: the removal needs Alk0 ≥ 9.5 meq/L if nitrified".
- **X6.** Sub-segment 14 becomes an electrode fragment (`segments.csv` status and diagnostic slope; `segments_window_check.csv` 10 → 9; `nh4_truth_windows.csv` status; `thermo_co2.csv` 12 → 11; `thermo_coexistence.csv` note "9 windows"; `kinetic_fits.csv`/`kinetic_parameters.csv` for segment 1; `dynamics.csv` (one window less in pulse 1's NH4 segment), `dynamics_windows.csv` 680 → 679, `dynamics_flux.csv`; 28/28 required estimator checks still pass). Pulse 1: first-window consumption 18.29 → 4.23 mg N/L/d (window 15), duration-weighted 7.15 → 4.23, O2 demand 32.7 → 19.3 mg O2/L/d, required kLa 3.96 → 2.34 1/d, gain-corrected SENSITIVITY 6.97 → 4.34, dose-removal time 9.15 → 15.5 d (7.59 x the pulse interval, was 4.49). The R8 CASE with window 14 included: 7.15 mg N/L/d (95 % 5.51 to 8.52), the old value. M4 references: "O2 demand max" 77.5 → 23.7 and "kLa required max" 9.39 → 2.87 (both still fail against 33 and 4.0, explained by the preliminary 2.59-8.82 window as before; pass count 9/11 unchanged; M3 11/13 unchanged). η, n_removed and the energy family of pulse 1 do not depend on the windows and are unchanged.
- **X7.** `pulses.csv`: `dose_removal_time_d`, `dose_removal_time_over_pulse_interval` (renamed), `eta_gap_inversion_misfit_MODEL_BASED` (+4 rows: 0.174, 0.237, −0.252, 0.161; equal to the ratio minus the MODEL-BASED η within 0.004), AE suffixes. Reports: `analysis_report.md` (§1 points 2, 5, 7, 8; §2 table and text; §3 states table; §4 convergence, Alk0 table, rounding; status table rows), `states.md`, `pulses.md`, `model.md`, `segments.md`/`thermo.md`/`dynamics.md` (renamed columns, window 14). Figures: Fig. 1 (idle bounds), Fig. 2 (window 14 demoted), Fig. 4 (window 14 gone from panel b; states timeline), Fig. 6 (refit in the 4th digit), Fig. 7 (pulse 1 rates), Fig. 8 (idle bounds with coverage, AE labels, data-driven titles). Figs. 3 and 5 are byte-identical.
- `pulses.csv` 336 → 352 rows (+4 classified_fraction, +4 idle (three idle metrics instead of two), +4 R8 CASE, +4 η gap); 274 computed, 78 not computed (was 265, 71).

### Commands and results

Environment: CPython 3.14.8, uv 0.13.1 (warns that `uv_build<0.13.0` does not contain it; harmless). The machine was shared with an unrelated multi-core process for the whole session, so wall times are long.

Baseline, on two clean detached worktrees of `954f87a` (one for run-all, one for pytest), before any edit:

```
$ uv run ruff check
All checks passed!
$ rm -rf outputs && uv run python -m bnr_neq run-all        # exit 0, real 2533.6 s
... model: selected biomass_decay; consistency checks pass; 13/15 parameters identifiable
states: 12 % of the record classified; idle (inactive, over-supplied) 3 % of stirred time
pulses: 265 metrics computed, 71 not computed
$ diff -rq <run-all outputs> <committed outputs of 954f87a>   # tables, reports, figures
(no differences)
$ uv run pytest -q
381 passed in 8549.89s (2:22:29)
```

Final code (after the last edit), in the working tree:

```
$ uv run ruff check
All checks passed!
$ rm -rf outputs && uv run python -m bnr_neq run-all        # exit 0, real 3363.2 s
ingest: 1 row(s) dropped
validate: verdicts nh4/ec/orp FAIL, ph PASS (unchanged)
segments: references 11/13 reproduced
thermo: references 9/11 reproduced
dynamics: estimator validation 28/28 required checks passed
model: selected biomass_decay; consistency checks pass; 13/15 parameters identifiable
states: 12 % of the record classified (coverage); idle (inactive, over-supplied) between 5 and 93 % of stirred time
pulses: 274 metrics computed, 78 not computed
report: wrote outputs/reports/analysis_report.md
```

Two development run-alls preceded it (one before the per-start table and the Alk0 multimodality rule, which failed at the report stage on the pandas-3 NaN issue; one stopped by hand to add the R8 CASE interval); their outputs were discarded and nothing was regenerated stage by stage after the final run.

The first full test run on the final code gave `1 failed, 417 passed in 6606.78s (1:50:06)`: the new `test_states.test_real_idle_metrics_carry_bounds_and_coverage` expected the below-coverage status for every interval under 0.5, but pulse 1 has no classified bin at all and correctly reads `not computed: no bin with a process state in this interval`. The assertion now accepts either not-computed reason (a mistake in the new test, not in the code; no source file changed). The suite was then run again in full:

```
$ uv run ruff check
All checks passed!
$ uv run pytest -q
418 passed in 5000.28s (1:23:20)
```

418 = the 381 tests of `954f87a` (none removed) plus 37 new ones. The suite includes the two fresh `run-all` executions with byte-identical CSV and SVG outputs (`test_main.test_run_all_is_deterministic`) and the real-record integration tests of this pass. Pass 3 is an implementer's log, not an approval; the reviewer decides.
