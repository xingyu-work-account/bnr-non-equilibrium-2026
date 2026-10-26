# Remediation prompt: linking sensor data to process state, energy and efficiency

You are the implementer for `bnr-non-equilibrium-2026`. Read `docs/review_2026-10-19.md` first and treat each finding as a claim to reconcile against the current code, not as established fact. Before editing, inspect the repository, run the existing pipeline and tests, and verify every path, symbol, column and number cited below; report any that differ. Make the smallest coherent change for each item, do not refactor unrelated code, do not weaken or delete existing tests, never fabricate a measurement, and never tune code to reproduce a number. You are not the approver; an independent reviewer decides.

## Problem and root cause

The pipeline audits the four probes thoroughly but never translates their signals into process state, oxygen demand, energy and efficiency, which is the scientific question. The root causes are (1) two incompatible NH4 observation models feeding different outputs, (2) no proton and alkalinity balance although pH is the only validated sensor, (3) an energy metric (stirrer power) that is neither measured nor transferable to BNR aeration, (4) no use of relative ORP and pH as process-state indicators, and (5) report and figure structure that foregrounds sensor verdicts and stochastic estimators over the process narrative.

## Work items

| Item | Change | Affected files (verify) |
| --- | --- | --- |
| W1 | Adopt one NH4 observation model for all downstream quantities. Use the M6 log model, log R = a + b(t) log C + offsets, inverted to C(t) with propagated uncertainty, and derive rates, FA, η and oxygen demand from that single C(t). Report the linear gain correction only as a labelled sensitivity case, and add a dose-anchored TAN bound, C_dose(t) = C(t0) + Σ ΔC − ∫ r dt, with r bounded by the alkalinity result of W2. | `segments.py`, `model.py`, `thermo.py`, `pulses.py` |
| W2 | Add a proton, alkalinity and carbonate module. Compute the upper bound on cumulative nitrification N_nit ≤ (Alk0 − Alk_min(pH, pCO2))/2 as a function of Alk0 over 1 to 15 meq/L, using the pH record, carbonate equilibria and a CO2 transfer term with kLa_CO2 as a scenario range; report N_nit(Alk0) as a curve with the dosed 262 mg N/L for comparison, and the implied minimum residual TAN and FA (Anthonisen, 20 and 25 °C scenarios). | new `alkalinity.py`; `thermo.py`; `config.py` |
| W3 | Add a process-state classifier using pre-registered rules on relative ORP level and slope, pH slope and NH4 slope (from W1) on the 0.01 d grid, with states active nitrification, inactive and over-supplied, oxygen-limited, and sensor artifact, plus a time-ordered state table. Test whether NH4 consumption differs between states with block-bootstrap intervals that respect autocorrelation. | new `states.py`; `pipeline.py` |
| W4 | Replace the energy layer with a metric family, each with its inputs and status: η_k per pulse; oxygen demand (4.57 g O2/g N, with 3.43 g O2/g N for nitritation as a sensitivity); oxygen supply capacity kLa·(C* − C)·V under a kLa scenario until measured; oxygen utilization efficiency = demand/supply; aeration-equivalent energy = demand/SAE with an explicit field SAE scenario range taken from a cited source; kWh per kg N converted; fraction of stirred time spent in the inactive state (W3); and Gibbs energy released over aeration-equivalent energy. Keep stirrer power only as a laboratory-specific term. | `pulses.py`, `thermo.py`, `constants.py` |
| W5 | Add an alkalinity or inorganic-carbon limited candidate (and its combination with FA inhibition using W1 TAN) to the M6 model set, with the chronological split unchanged, so that biomass decay competes against substrate-side limitation. | `model.py` |
| W6 | Produce the missing figures: Fig. 1 conceptual linkage (sensors → state → rates → O2 demand → energy → efficiency), Fig. 4 (N_nit(Alk0) bound, rate against FA from W1, state timeline), Fig. 8 (per-pulse η, oxygen utilization and idle fraction under the SAE and kLa scenarios). | `report.py`, figure functions |
| W7 | Restructure the report narrative-first: a one-page summary of what the data say about the process, then energy and efficiency, then sensor validity, with stochastic-thermodynamic results moved to a supplementary section. Relabel the coexistence affinity, ΔA and σ ≥ 0 as consistency checks, not findings. Round values to meaningful precision. | `report.py` |
| W8 | Add `vessel_mass_start_g`, `vessel_mass_end_g`, `covered`, `top_up_mL` and `initial_alkalinity_meq_L` as optional config keys under `[missing]`, and carry a volume scenario V(t) through every concentration balance, so evaporation is either excluded by measurement or bounded. | `config/experiment.toml`, `config.py`, balances |

## Equations, units and invariants

Alkalinity consumption is 2 eq per mol N nitrified (7.14 g CaCO3/g N); air-equilibrium bicarbonate is [HCO3−] = K1·KH·pCO2/[H+]; FA = (17/14)·TAN·10^pH/(exp(6344/T) + 10^pH) in mg NH3/L with TAN in mg N/L and T in K; oxygen demand 4.57 g O2/g N. Invariants: one NH4 truth per run; N_nit never exceeds the W2 bound; every scenario value carries its scenario label and range; no window spans a pulse or an electrode jump; σ ≥ 0 checks are reported as consistency checks.

## Non-goals

No change to sensor verdicts without new measurements, no removal of existing analyses (demote instead), no full ASM2d build, and no claim of N2O emission.

## Tests

| Type | Required tests |
| --- | --- |
| Positive | W2 bound reproduces 14, 35 and 70 mg N/L for Alk0 of 2, 5 and 10 meq/L with pH floor ignored; air-equilibrium HCO3− is 0.401 mM at pH 7.8 and 0.636 mM at pH 8.0 (25 °C); FA at 200 mg N/L TAN, pH 7.8, 25 °C is 8.5 ± 0.1 mg NH3/L. |
| Negative | Mixing observation models in one output raises; a negative Alk0 or SAE raises; a state rule referencing an excluded sensor raises. |
| Boundary | Alk0 at the pH-floor limit gives N_nit = 0; a window shorter than the minimum span returns a status, not a number; an all-artifact interval yields no state. |
| Regression | Existing tests and reference tables still pass, with any changed reference explained by the W1 model switch. |
| Integration | Rates, FA, η and O2 demand for one pulse are computed from the same inverted C(t) (asserted by provenance column). |
| End-to-end | `uv run python -m bnr_neq run-all` regenerates all outputs including Figs. 1, 4 and 8 from the raw files with fixed seeds. |

## Validation commands

`uv run ruff check`, `uv run pytest -q`, `uv run python -m bnr_neq run-all`, followed by a diff of `outputs/tables` against the previous run with every changed value explained.

## Acceptance criteria and failure conditions

Accept only if every work item is implemented and tested, all commands pass, the report states the W2 bound and the energy metric family with their scenario ranges, and every efficiency number traces to one observation model. The work fails if any output still mixes observation models, if the alkalinity bound is omitted or contradicted without evidence, if a scenario value is presented as measured, or if any existing test is weakened. Report exactly what was changed, the commands run and their full results.
