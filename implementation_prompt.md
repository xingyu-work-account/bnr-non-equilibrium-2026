# Implementation prompt: non-equilibrium analysis of an NH₄Cl-pulsed wastewater batch

## 1. Your role and ground rules

You are the implementer of a scientific data-analysis pipeline. An independent reviewer will assess your work afterwards, so you are never the approver of your own results. Every number you report must come from code you executed in this session.

1. **Inspect before you edit.** Before writing any code:
    - list the repository tree;
    - read `pyproject.toml`, `.python-version`, `uv.lock` (if present), `README.md` and any existing source files;
    - print the header and the first and last five lines of each file in `data/`.

   Report what you found, separating confirmed facts from assumptions.
2. **Verify everything you use.** Check every path, column name, unit and library symbol in the repository or the installed package, never from memory.
3. **Keep changes minimal and the raw data untouched.** Make the smallest coherent change for each milestone, and do no unrelated refactoring. Never modify, move or rewrite anything in `data/`; all outputs go to `outputs/`.
4. **No silent fallbacks.** A missing input must produce either an explicit exception naming the input, or an explicit `not computed: <reason>` entry in the outputs. Never substitute a default value that could be mistaken for a measurement. Every dropped row must be counted and reported.
5. **Test the production path.** Tests must call the same functions the pipeline calls, not re-implementations. Never weaken, skip or delete a test to make it pass.
6. **Execute and report exactly.** Run every validation command and report its exact output, including failures.
7. **Gates.** Each milestone ends with a gate. If a gate check fails, stop, report the failure with evidence, and do not continue to the next milestone.
8. **Git.** Do not commit, push or create branches unless asked.

## 2. Environment: uv only

The project folder was initialized with `uv`, but no virtual environment exists yet. Use `uv` for every Python operation: never call `pip` directly, never activate the environment manually, and run everything through `uv run`.

1. Inspect `pyproject.toml` and `.python-version`, and report the Python version and whether a `[build-system]` table exists.
2. Create the environment with `uv sync`, which creates `.venv` and `uv.lock`, then confirm it with `uv run python --version`.
3. **Package layout.** The analysis code lives in an importable package `src/bnr_neq/`.
    - If `pyproject.toml` has no `[build-system]` (the default for `uv init` applications), add a hatchling build system so that `uv sync` installs the project in editable mode:

      ```toml
      [build-system]
      requires = ["hatchling"]
      build-backend = "hatchling.build"

      [tool.hatch.build.targets.wheel]
      packages = ["src/bnr_neq"]
      ```

    - Then run `uv sync` again, and confirm with `uv run python -c "import bnr_neq; print(bnr_neq.__file__)"`.
    - Leave any existing `main.py` untouched unless it conflicts, and report it either way.
4. **Add libraries only with `uv add`, never by editing the dependency list by hand.**
    - Runtime: `uv add numpy pandas scipy matplotlib statsmodels ruptures pyarrow`.
    - Development: `uv add --dev pytest ruff`.

   Add a further library only when a milestone needs it, and report why.
5. **Run scripts and tests only through uv:**
    - `uv run python -m bnr_neq <command>`
    - `uv run pytest -q`
    - `uv run ruff check .`
6. **Configuration.** Use `tomllib` from the standard library to read `config/experiment.toml` (Section 5). If the Python version is older than 3.11, report this and run `uv add tomli`.

## 3. Experimental facts, confirmed by the experimenter

Treat these as ground truth, and encode them in `config/experiment.toml`, not in code.

1. **Vessel.** The system is a batch of 100 mL of wastewater, in a vessel open to air and mixed by a magnetic stirrer, with no forced aeration.
2. **Additions.** Over about 50 days, the only material added was 25 mg of solid NH₄Cl at days 8.53, 10.57, 24.59 and 32.39, all on the logger clock. Nothing else was added or removed.
3. **Probes.** The conductivity, NH₄⁺, pH and ORP probes all sit in the same liquid, at slightly different positions.
4. **Conductivity.** The raw column is a current in amperes, converted as EC (mS/cm) = −2.89855 × I (A) + 0.520870.
5. **NH₄⁺.** The electrode reports mg NH₄⁺/L. Convert to mg N/L with the factor 14.007/18.039 = 0.77649.
6. **ORP.** ORP is the potential of a Pt working electrode against a bare Ag quasi-reference, in mV. Its offset to the standard hydrogen electrode is unknown.

## 4. Data files and verified properties

A preliminary inspection found the properties below. Re-verify every one of them, and report any difference before proceeding. All files use CRLF line endings, and `NH4.csv` starts with a UTF-8 byte-order mark, so read every file with `encoding="utf-8-sig"`.

| File | Header (exact) | Valid rows | Time span (d) | Native sampling | Timestamp resolution |
| --- | --- | --- | --- | --- | --- |
| `data/Cond.csv` | `"Time, d",Cond. A` | 71,759 of 71,760; the last row has an empty value | 0.000693 to 49.6946 | 60 s, uniform | Full precision |
| `data/NH4.csv` | `Time(d),Concentration` | 414,785 | 0 to 47.92 | 20 s until day 41.05, then about 3 s, then about 2 s from day 45 | Rounded to 0.01 d; 22 to 450 rows share a stamp |
| `data/ORP.csv` | `Time,ORP (mV)` | 459,884 | 0 to 53.23 | 10 s | Rounded to 0.01 d; 21 to 87 rows share a stamp |
| `data/pH.csv` | `Time(d),pH` | 414,757 | 0 to 48.0044 | 10 s, with two 20 s intervals near day 4.398 | Full precision, values quantized to 0.01 |

The raw value ranges are:

- **Cond.csv:** −0.366 to −0.05946 A, which converts to 0.6932 to 1.5817 mS/cm;
- **NH4.csv:** 9.52 to 75.13 mg NH₄⁺/L;
- **ORP.csv:** 77.96 to 318.7 mV;
- **pH.csv:** 6.47 to 8.02.

The four-sensor common window is day 0 to 47.92.

## 5. Configuration file

Create `config/experiment.toml` containing the confirmed facts, and null-valued placeholders for measurements that do not yet exist. TOML has no null, so represent a missing value by omitting the key, and list the missing keys under `[missing]`. Code must read these values from the file and must emit `not computed: <key> missing` when a required value is absent.

```toml
[vessel]
volume_L = 0.100
open_to_air = true
stirring = "magnetic"

[pulses]
salt = "NH4Cl"
mass_mg = 25.0
molar_mass_g_mol = 53.491
times_d = [8.53, 10.57, 24.59, 32.39]

[sensors]
ec_slope_mS_cm_per_A = -2.89855
ec_intercept_mS_cm = 0.520870
nh4_unit = "mg NH4+/L"
nh4_to_N_factor = 0.77649
orp_working = "Pt"
orp_reference = "bare Ag"

[scenarios]
temperature_C = [20.0, 25.0]          # labeled scenarios until temperature is measured
baseline_chloride_mM = [1.0, 3.0, 8.0, 13.5]
pCO2_air_atm = 4.2e-4

[missing]
keys = ["temperature_C_measured", "alkalinity_mmol_L", "DO_mg_L", "kLa_per_d",
        "stirrer_power_W", "chloride_mmol_L", "orp_reference_offset_V",
        "lab_NH4_N_mg_L", "NO3_N_mg_L", "NO2_N_mg_L", "electrode_step_response"]
```

Scenario values must always appear labeled as scenarios in tables and figures, and never as measurements.

## 6. Constants and equations

Put every constant in `src/bnr_neq/constants.py`, with its units and literature source in a comment. Use 25 °C values unless a temperature is supplied, and use R = 8.314462618 J/(mol·K).

1. **Pulse dose.** c = mass / (molar mass × V) = 4.6737 mmol/L per pulse, which corresponds to 84.31 mg NH₄⁺/L, 65.46 mg N/L and 165.7 mg Cl⁻/L.
2. **Kohlrausch additivity at infinite dilution.** The limiting ionic conductivities (S cm²/mol) are NH₄⁺ 73.5, K⁺ 73.5, Na⁺ 50.1, Cl⁻ 76.3, NO₃⁻ 71.4 and HCO₃⁻ 44.5. The conductivity change is ΔEC (mS/cm) = Σ c_i (mmol/L) × λ_i × 10⁻³.
    - Immediately after a pulse, ΔEC = 0.700 mS/cm.
    - After complete nitrification, the net change is c × (λCl⁻ + λNO₃⁻ − 2λHCO₃⁻) = 0.274 mS/cm.
    - State in the outputs that values at the actual ionic strength are lower than these limiting values.
3. **Bare Ag, if chloride-controlled.**
    - Reference potential: E_ref = 0.2223 − (RT ln10 / F) log₁₀(γ·[Cl⁻]) V vs SHE.
    - Predicted rise of ORP (Pt vs Ag) for a chloride addition Δc from a baseline c₀: (RT ln10 / F) log₁₀((c₀ + Δc)/c₀).
    - Alternative, oxide-controlled scenario: E_ref = 1.173 − (RT ln10 / F) pH.
4. **Ammonia speciation.**
    - Emerson et al.: pK_a = 0.09018 + 2729.92 / T (K).
    - Anthonisen free ammonia: FA (mg NH₃/L) = (17/14) × TAN (mg N/L) × 10^pH / (exp(6344/T) + 10^pH).
5. **Activity coefficients (Davies).** log₁₀ γ = −A z² (√I/(1+√I) − 0.3 I), with A = 0.509 at 25 °C.
6. **Nitrification free energy.** The reaction is NH₄⁺ + 2 O₂(g) → NO₃⁻ + 2 H⁺ + H₂O, with ΔG°′ = −348.8 kJ/mol N at pH 7. The sub-steps are nitritation, −274.7 kJ/mol, and nitratation, −74.1 kJ/mol, from the electron-equivalent tables of Rittmann and McCarty.
    - In situ: ΔG = ΔG°′ + RT ln([NO₃⁻] / ([NH₄⁺] p_O₂²)) + 2RT ln([H⁺]/10⁻⁷), where p_O₂ = 0.2095 × DO/C*(T) atm.
    - The affinity is A = −ΔG.
7. **Affinity step at a pulse.** ΔA_k = RT ln((C_k⁻ + ΔC)/C_k⁻), with C in mg N/L and ΔC = 65.46 mg N/L.
8. **Alkalinity and oxygen stoichiometry.** Nitrification consumes 7.14 g CaCO₃ and requires 4.57 g O₂ per g N. The required surface transfer coefficient is k_La ≥ demand / C*, with C* ≈ 8.26 mg/L at 25 °C, to be replaced by a temperature function when temperature exists.
9. **CO₂ at the open surface.** pK₁ = 6.35 and log₁₀ K_H = −1.47 (mol/(L·atm)) at 25 °C.
    - Air-equilibrium bicarbonate: [HCO₃⁻]_eq = K₁ K_H p_CO₂,air / [H⁺].
    - Supersaturation, which requires alkalinity: Ω_CO₂ = [HCO₃⁻][H⁺] / (K₁ K_H p_CO₂,air).
10. **Non-equilibrium criteria.**
    - De Donder relation: A_j r_j ≥ 0.
    - Entropy production rate: σ = Σ r_j A_j / T.
    - Thermodynamic factor (Jin and Bethke): F_T = 1 − exp(−(−ΔG − mΔG_ATP)/(χRT)).

## 7. Package structure

Create these modules, keeping functions small and pure where possible:

- `src/bnr_neq/`:
    - `__init__.py`
    - `__main__.py`: command-line interface with the commands `ingest`, `validate`, `segments`, `thermo`, `dynamics`, `model`, `pulses`, `report` and `run-all`;
    - `config.py`: loads and validates `config/experiment.toml`;
    - `constants.py`
    - `io.py`: file loading and unit conversion;
    - `timebase.py`: stamp reconstruction, gap inventory and the common grid;
    - `qc.py`: flags, the event detection the pipeline uses, and held-output detection;
    - `validation.py`: Phase 0 sensor validation;
    - `segments.py`: Phase 1 segmentation and kinetics;
    - `thermo.py`: Phase 2 thermodynamics and balances;
    - `estimators.py`: Phase 3 irreversibility, entropy and surrogate methods;
    - `dynamics.py`: Phase 3 applied to the data;
    - `model.py`: Phase 4 reduced model;
    - `pulses.py`: Phase 5 pulse-level efficiency and energy descriptors;
    - `report.py`: Markdown report assembly.
- `tests/`: one test file per module.
- `outputs/` (generated, never edited by hand):
    - `processed/`: parquet files;
    - `tables/`: CSV files;
    - `figures/`: PNG and SVG files;
    - `reports/`: Markdown files.

Design the estimator functions in `estimators.py` with plain NumPy array inputs and outputs and no pandas dependency. The same algorithms are meant to be ported later to a Rust sensor-entropy library.

## 8. Milestones

### M0. Environment and inspection

1. **Implement.** Carry out Sections 1 and 2.
2. **Gate.** `uv sync`, the `uv add` commands and the import check all succeed.

### M1. Ingestion, time base and quality flags (Phase 0, part 1)

1. **Loading.**
    - Map the verified headers explicitly:
        - `Cond.csv`: `"Time, d"` → `t_d` and `Cond. A` → `I_A`;
        - `NH4.csv`: `Time(d)` → `t_d` and `Concentration` → `nh4_mgNH4_L`;
        - `ORP.csv`: `Time` → `t_d` and `ORP (mV)` → `orp_mV_vs_Ag`;
        - `pH.csv`: `Time(d)` → `t_d` and `pH` → `pH`.
    - Raise an error, naming the file, if any expected header is absent.
    - Convert with `pd.to_numeric(errors="coerce")`, then count, log and drop the non-numeric rows. Exactly one is expected, the final row of `Cond.csv`.
2. **Unit conversions.**
    - Convert EC as in Section 3. Raise an error if any converted value lies outside 0 to 5 mS/cm, which would indicate a unit error.
    - Add `nh4_mgN_L` as a separate column, keeping the original column alongside.
3. **Monotonicity.** Raise an error if any time decreases.
4. **Stamp reconstruction for NH₄⁺ and ORP.**
    - First verify that the spacing between distinct stamps is 0.01 d.
    - Then set t = stamp + 0.01 × j/n, where j is the zero-based row index within the stamp and n the number of rows sharing it.
    - Record this as a stated assumption (uniform sampling within a stamp) in the outputs.
    - The reconstructed times must be strictly increasing.
5. **Gap inventory.** List every interval longer than twice the median native interval, per sensor.
6. **Common grid.**
    - Build 0.01-day bins over [0, 47.92] with time-weighted means, where each sample is weighted by its own sampling interval, so that NH₄⁺ logging-rate changes cannot bias the means.
    - Store the sample count per bin and a mask for empty bins.
    - Never interpolate across gaps longer than twice the native interval.
7. **Outputs.** Write `outputs/processed/native_<sensor>.parquet` and `outputs/processed/grid_0p01d.parquet`.
8. **Tests.**
    - Row counts and spans as in Section 4.
    - The EC conversion: I = −0.1 A must give 0.810725 mS/cm, and a value in mA must raise an error.
    - Stamp-size ranges: 22 to 450 rows for NH₄⁺, 21 to 87 for ORP.
    - Reconstructed times are strictly increasing.
    - Grid bins contain no fabricated values.
    - The dropped-row count equals 1.
9. **Gate.** All tests pass, and every Section 4 property is reproduced or its discrepancy explained.

### M2. Sensor validation with the known pulses (Phase 0, part 2)

The four pulses are standard additions of known size. Compute predicted and observed responses with these exact windows:

- **NH₄⁺ and EC:** the median from 0.4 to 0.6 d after the pulse, minus the median over the 0.1 d before;
- **ORP and pH:** the median from 1 to 3 h after the pulse, minus the median from 3 to 1 h before.

1. **Pulse table.** For each pulse, record:
    - the predicted and observed NH₄⁺ change and the response fraction;
    - the predicted and observed EC change;
    - the predicted ORP change for each baseline chloride scenario;
    - the observed pH change.
2. **Cumulative balances.**
    - **EC staircase.** Plot the predicted EC staircase under no nitrification and under complete nitrification against the observed record.
    - **Nitrogen.** Compute the nitrogen balance: the start median (day 0.3 to 1.3), the end median (day 46.9 to 47.9), the 262 mg N/L added, and the implied conversion rate.
3. **NH₄⁺ jump detection.**
    - Detect abrupt NH₄⁺ jumps with a robust rolling-median step detector on the reconstructed native series, with documented thresholds.
    - Classify each jump as a pulse response (state the lag from the pulse) or an electrode offset unrelated to any pulse.
    - Measure the 10 to 90% transition time of each jump.
4. **Held-output detection.** Detect held or damped NH₄⁺ output as segments with near-zero high-frequency residual.
5. **Outputs.** Write `outputs/tables/sensor_validation.csv`, `outputs/figures/fig3_sensor_validation.png` and `outputs/reports/sensor_validation.md`. The report states plainly which sensors pass and which fail.
6. **Reference values from the preliminary analysis.** Reproduce each within the stated tolerance, or explain the discrepancy with evidence. Never tune code to hit a number.

| Quantity | Reference | Tolerance |
| --- | --- | --- |
| NH₄⁺ response fraction, pulses 1 to 4 | 37, 29, 26, 17% | ±3 percentage points |
| Observed NH₄⁺ change, pulses 1 to 4 | +31.1, +24.4, +22.1, +14.7 mg NH₄⁺/L | ±2 mg/L |
| Observed EC change, pulses 1 to 4 | +0.005, −0.062, −0.007, +0.016 mS/cm | ±0.01 mS/cm |
| Observed ORP change, pulses 1 to 4 | +0.4, 0.0, +0.1, +0.2 mV | ±0.5 mV |
| Observed pH change, pulses 1 to 4 | −0.01, +0.02, 0.00, +0.04 | ±0.02 |
| EC median, day 0.3 to 1.3 and day 48.5 to 49.6 | 1.069 and 1.071 mS/cm | ±0.005 |
| NH₄⁺ median, start and end | 47.1 and 43.1 mg NH₄⁺/L | ±0.5 |
| NH₄⁺ jumps unrelated to pulses | days 2.59, 17.84, 37.19, 44.63 | ±0.02 d |
| Delayed NH₄⁺ responses to pulses 1 and 2 | days 8.85 and 10.88 (about 7.5 h lag) | ±0.02 d |
| Late EC pulses after day 41.9 | 29 events, median exponential decay about 24 min | count ±3 |

7. **Gate.** The table and report reproduce the references, and the report states the P0 findings: the NH₄⁺ under-response, no EC response, and no ORP chloride response.

### M3. Segmentation and consumption kinetics (Phase 1)

1. **Segmentation.**
    - Segment at the four pulses: day 1 to 8.53, 8.53 to 10.57, 10.57 to 24.59, 24.59 to 32.39, and 32.39 to 47.92.
    - Split into sub-segments at the unrelated electrode jumps, which are treated as offsets.
2. **Slopes.** Compute Theil–Sen slopes (`scipy.stats.theilslopes`) with 95% intervals on each sub-segment, excluding the first 0.1 d after each boundary.
    - Report the slopes in mg NH₄⁺/L/d and in mg N/L/d.
    - Report gain-corrected rates, dividing by the response fraction of the governing pulse, as a separate, clearly labeled column.
3. **Free ammonia.** Compute the median free ammonia per sub-segment for each temperature scenario.
4. **Kinetic models.** Fit zero-order, first-order, Monod and Haldane (free-ammonia-inhibited) models per segment. Compare them by AIC, and report residual autocorrelation and any parameters that are not identifiable.
5. **Change points and periodicity.**
    - Run PELT change-point detection (`ruptures`) on the grid, and tabulate the events not tied to pulses.
    - Compute Lomb–Scargle periodograms of the detrended grid series.
6. **Reference values.**
    - Uncorrected slopes (mg NH₄⁺/L/d, ±10% or ±0.05), for the sub-segments 0.3 to 2.59, 2.59 to 8.82, 10.85 to 17.84, 17.84 to 24.60, 24.60 to 32.36, 32.36 to 37.19, 37.19 to 44.63 and 44.63 to 47.92: −1.89, −3.42, −1.95, −0.57, −0.48, −0.31, −0.28, +0.12.
    - Gain-corrected rate after pulse 2: about 5.2 mg N/L/d; after pulse 4: about 1.4 mg N/L/d.
    - Free ammonia from day 17.84 to 32.36: 0.7 to 1.2 mg NH₃/L across the 20 and 25 °C scenarios.
    - Maximum Lomb–Scargle normalized power: about 0.05.
7. **Outputs.** Write `outputs/tables/segments.csv` and `outputs/figures/fig2_overview.png`, showing the pulses, flags and electrode events.
8. **Gate.** References are reproduced, and no slope window spans a pulse or an uncorrected jump.

### M4. Thermodynamics and conservation balances (Phase 2)

1. **Functions.** Implement the Section 6 functions with explicit units in names or docstrings. Use the guard pattern for missing inputs: a function that needs alkalinity, DO, chloride or the ORP offset must return a `not computed` result when the value is missing.
2. **Affinity tables.**
    - Affinity steps at the pulses. References: 4.2, 2.5, 3.1 and 2.4 kJ/mol at 25 °C, ±0.1.
    - Coexistence affinity for the documented scenario (NH₄⁺ 40 mg N/L, NO₃⁻ 100 mg N/L, DO 0.1 mg/L, pH 7, 25 °C). Reference: 317 ± 2 kJ/mol. Also tabulate its sensitivity to DO from 0.01 to 8 mg/L.
3. **Chloride scenarios.** Predict the bare Ag responses under each chloride scenario. Compute no absolute E_H unless `orp_reference_offset_V` is present; otherwise emit `not computed: S4`.
4. **CO₂ at the open surface.** Compute the air-equilibrium bicarbonate at the observed pH. References: 0.127 mmol/L at pH 7.3 and 0.401 mmol/L at pH 7.8. Compute Ω_CO₂ only if alkalinity is present.
5. **Oxygen.** Compute the oxygen demand from the gain-corrected rates and the required k_La. References: about 6 to 33 mg O₂/L/d, and k_La ≥ 0.7 to 4 d⁻¹.
6. **Balances.** Set out the nitrogen, chloride and alkalinity balances as tables. Every term not computable from the data appears as `not computed: <missing key>`.
7. **Unit tests.**
    - pK_a at 25 °C equals 9.246 ± 0.002.
    - The Davies coefficient for z = 1 at I = 0.01 M equals 0.902 ± 0.002.
    - The Ag/AgCl potential at unit chloride activity equals 0.2223 V, with a slope of 59.16 mV per decade at 25 °C.
    - Physically impossible inputs, such as pH outside 0 to 14 or negative concentrations, raise errors.
8. **Gate.** All references are reproduced, and every missing-input path is tested.

### M5. Dynamical and information-theoretic signatures (Phase 3)

1. **Estimators** (in `estimators.py`, NumPy only):
    - **Moment asymmetry.** ⟨(x(t+τ) − x(t))³⟩ as a function of lag τ.
    - **Ordinal-pattern irreversibility.** The Kullback–Leibler divergence between ordinal-pattern frequencies and those of the time-reversed series, with explicit tie handling. The pH series is quantized to 0.01, and 88% of consecutive pH values are equal, so treat equalities as their own pattern symbols, and report the tie fraction.
    - **Visibility-graph divergence.** The divergence between in-degree and out-degree distributions of the horizontal visibility graph.
    - **Surrogates.** Iterated amplitude-adjusted Fourier transform (IAAFT) surrogates.
    - **Complexity.** Permutation entropy with ties, sample entropy and multiscale entropy, and detrended fluctuation analysis.
    - **Probability flux.** Net probability flux on a two-dimensional binned state space, with block-bootstrap confidence intervals, following Battle et al. (2016).
2. **Estimator validation, mandatory before applying to data.**
    - Simulate a three-state ring Markov chain, with forward probability p, backward probability q and stay probability 1 − p − q.
    - The transition-level estimator must recover the analytic entropy production per step, (p − q) ln(p/q), within its confidence interval for at least three (p, q) pairs, including p = q, where it must give zero.
    - Repeat the test on a version of the chain observed through a noisy observable quantized to 0.01 units, and report the resulting bias.
    - IAAFT surrogates of a reversible process must give indices statistically indistinguishable from zero.
3. **Application to the data.**
    - **pH,** the primary signal: native 10 s series, excluding day 0 to 1.
    - **Relative ORP:** excluding the electrode events at days 28.89 and 30.48.
    - **Offset-corrected NH₄⁺:** excluding held segments.
    - **Conductivity:** excluded from process analyses, and characterized only as probe behavior.

   Run the augmented Dickey–Fuller and KPSS tests (`statsmodels`) before using any stationary estimator.
4. **Onsager regression test.** Compare the normalized late-tail relaxation after each pulse with the normalized autocorrelation of spontaneous fluctuations in the same segment.
    - Report the ratio against lag, with block-bootstrap intervals.
    - Report whether the signal-to-noise ratio is sufficient. If it is not, say so; do not force a conclusion.
5. **Stability indicators.** Compute rolling variance and lag-1 autocorrelation before each pulse.
6. **Electrode dynamics.** NH₄⁺ relaxation times must not be interpreted as process relaxation until the electrode step response exists. Without it, emit `not computed: electrode_step_response` for deconvolved NH₄⁺ relaxation.
7. **Outputs.** Write `outputs/tables/dynamics.csv` and `outputs/figures/fig5_dynamics.png`.
8. **Gate.** All estimator validation tests pass before any result on plant data is reported.

### M6. Reduced model with the electrode observation model (Phase 4)

1. **State model.** Model NH₄⁺ in mg N/L, with impulses of +65.46 mg N/L at the pulse times. Compare three consumption laws: Haldane kinetics in free ammonia, a first-order decay of biomass activity, and both combined.
2. **Observation model.** log R = a + b log C(t) + offsets at the four unrelated jumps, with a first-order sensor lag τ. Allow a slow linear drift in b, to represent the declining response.
3. **Fitting.** Fit with `scipy.optimize.least_squares`, using a chronological split:
    - calibrate on days 1 to 24.59;
    - select the model on days 24.59 to 32.39;
    - test on days 32.39 to 47.92.

   Refit the selected model on the whole record only after reporting out-of-sample skill from the split.
4. **Identifiability.** Assess every parameter by profile likelihood. Report non-identifiable parameters as such, and never present them as estimates.
5. **Consistency checks.** Each simulated pulse adds exactly 65.46 mg N/L, mass balance closes, and σ ≥ 0 wherever affinities are computed.
6. **Outputs.** Write `outputs/tables/model_fit.csv` and `outputs/figures/fig6_model.png`.
7. **Gate.** Mass-balance tests pass, the out-of-sample skill is reported, and the identifiability table is complete.

### M7. Pulse-level efficiency and energy descriptors (Phase 5)

1. **Pulse table.** For each pulse, report:
    - ΔA_k;
    - the uncorrected and gain-corrected consumption rates;
    - median pH and free ammonia per temperature scenario;
    - the oxygen demand and the required k_La;
    - the model-based conversion fraction η_k, labeled as model-based.
2. **Energy metrics.** Specific stirring energy per g N, and the ratio of the free energy released (163 J per fully converted pulse in 100 mL) to the stirrer's electrical energy, require `stirrer_power_W`. Without it, emit `not computed: stirrer_power_W`.
3. **Statistics.** With only four pulses, mark every statistical association as exploratory.
4. **Outputs.** Write `outputs/tables/pulses.csv` and `outputs/figures/fig7_pulses.png`.
5. **Gate.** Every metric is either computed with its inputs listed, or explicitly marked as not computed.

### M8. Report and single-command reproduction

1. **Report.** `uv run python -m bnr_neq run-all` regenerates every output from the raw files, with fixed random seeds, and writes `outputs/reports/analysis_report.md`. The report contains:
    - the experimental facts;
    - the verified data properties;
    - the sensor validation verdicts;
    - every table and figure;
    - every stated assumption;
    - a status label for each conclusion.
2. **Status labels.** These mirror the plan:
    - **GO:** chemical disequilibrium from the coexistence criterion, and the smallness of the affinity steps;
    - **GO WITH DOCUMENTED NON-BLOCKING DEBT:** the pH trajectory analyses and the relative comparison of pulse responses under the stated electrode gain model;
    - **NO-GO:** absolute NH₄⁺ rates and conversion, conductivity as solution conductivity, absolute E_H, and oxygen and energy claims, each with the missing item that blocks it.
3. **Gate.** `run-all` completes from a clean `outputs/` directory, and a second run produces identical tables.

## 9. Validation commands

Run all of these at the end, and after any milestone that changes code. Report their exact output.

```bash
uv sync
uv run ruff check .
uv run pytest -q
rm -rf outputs && uv run python -m bnr_neq run-all
uv run python -m bnr_neq run-all   # second run, for determinism
```

Compare the CSV tables from the two runs, and report any difference.

## 10. Out of scope

- Any conclusion that uses conductivity as solution conductivity.
- Absolute redox potentials.
- Energy numbers without stirrer power.
- A full ASM2d model.
- Any change to `data/`.
- Porting the estimators to Rust.

## 11. Final report to the user

End with a report containing exactly these sections:

1. **Inspection findings.** Confirmed facts, and assumptions with the evidence for each.
2. **Environment.** Python version, every `uv` command run, and the resulting dependency list.
3. **Files created or changed.** Each with a one-line purpose.
4. **Milestone gates.** The result of each gate, with evidence.
5. **Reference-value comparison.** A table of reference, obtained value, tolerance and pass or fail, for every reference value above.
6. **Validation output.** The exact output of every command in Section 9.
7. **Deviations.** Every departure from this prompt, with its reason.
8. **Open issues.** Ranked P0, P1 or P2.
9. **What requires new measurements.** The keys under `[missing]`, and which outputs each one unlocks.

Do not declare the work approved. The independent reviewer decides.
