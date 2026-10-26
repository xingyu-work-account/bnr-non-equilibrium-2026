# Measurement request: what the BNR batch data cannot settle without the lab (2026-10-26)

The analysis pipeline for the NH4Cl-pulsed nitrifying batch is complete and independently reviewed (passes 1 to 5, merged as `0f7abfc`). Every conclusion that is still NO-GO or QUALIFIED is blocked by a missing measurement, not by code. This document lists those measurements in priority order: what to measure, how, when and on which sample, and which conclusion each one unlocks.

The source of truth is the status table in `outputs/reports/analysis_report.md` (§6) and the `[missing]` keys of `config/experiment.toml`.

## How the values enter the analysis

A measured value goes into a `[measurements]` table in `config/experiment.toml`, and its key is removed from `[missing]`. Then `uv run python -m bnr_neq run-all` recomputes everything. For example:

```toml
[measurements]
initial_alkalinity_meq_L = 4.2
vessel_mass_start_g = 152.3
vessel_mass_end_g = 148.1
covered = false
top_up_mL = 0.0
```

Two kinds of item below:
- **Ready:** the pipeline already uses the value as soon as it is entered.
- **Needs code:** the value has a config key, but using it fully needs a small pipeline change. This is noted per item.

**Record the time of every sample in days since the start of the run** (the same clock as the sensor files), and keep the probe readings at that moment.

## Priority list

| # | Measurement (config key) | Unlocks | Ready? |
| --- | --- | --- | --- |
| 1 | Laboratory NH4-N (`lab_NH4_N_mg_L`) | The absolute NH4 level C(t), and with it absolute rates and conversion; the alkalinity contradiction (S2); the electrode gain | **Needs code** to anchor the NH4 model |
| 2 | Alkalinity: initial and current (`initial_alkalinity_meq_L`, `alkalinity_mmol_L`) | The W2 nitrification bound as a value instead of a 1-15 meq/L scan; the test of the S2 requirement (Alk0 ≥ 20.8 meq/L) | Ready |
| 3 | Nitrate and nitrite at the end (`NO3_N_mg_L`, `NO2_N_mg_L`) | The N balance; nitrified vs stripped N (the "1 eq/mol N sink" alternative in S2); the in-situ affinity | Ready |
| 4 | Vessel mass at start and end; cover; water added (`vessel_mass_start_g`, `vessel_mass_end_g`, `covered`, `top_up_mL`) | The closed 100 mL frame (W8 NO-GO); which evaporation scenario applies | Ready |
| 5 | Dissolved oxygen (`DO_mg_L`) | The in-situ affinity (one value); a single idle fraction instead of bounds (a logged series) | Ready for one value (affinity); **needs code** for a series (states, idle fraction) |
| 6 | O2 transfer coefficient of the vessel (`kLa_per_d`) | O2 supply and aeration energy as values instead of SCENARIOs; the NH3 stripping term of the N balance | Ready |
| 7 | Liquid temperature (`temperature_C_measured`) | Replaces the 20/25 °C SCENARIOs in free ammonia, O2 saturation and thermodynamics | **Needs code**: the key is named in the report but no computation reads it yet |
| 8 | NH4 electrode step response (`electrode_step_response`) | Explains the delayed electrode responses; the absolute-concentration NO-GO | Needs code to use it beyond the report |
| 9 | Chloride and ORP reference offset (`chloride_mmol_L`, `orp_reference_offset_V`) | Absolute redox potential E_H (NO-GO) | Ready |
| 10 | Conductivity probe check (no key; a KCl standard addition) | Whether the conductivity record can be used at all (NO-GO) | Needs code |
| 11 | Stirrer power (`stirrer_power_W`) | Laboratory-specific stirring energy. Low priority: not transferable to a full-scale plant | Ready |

## Details

### 1. Laboratory NH4-N (highest priority)

- **Why.** The electrode registers only 17-37 % of each pulse. The model that converts its reading into true NH4 cannot fix the absolute level: the initial concentration and the electrode intercept are not identifiable, and the C(t) band is open upward. One independent concentration at a known time pins that level. Everything absolute depends on it: rates, conversion η, free ammonia, oxygen demand, and the alkalinity contradiction.
- **How.** Take a sample from the vessel and filter it at 0.45 µm. Preserve it with H2SO4 to pH < 2 and keep it at 4 °C; analyse within 7 days. Use a reference method (distillation and titration, or a salicylate or phenate colorimetric method, Standard Methods 4500-NH3), and report mg NH4-N/L.
- **When.** At least one sample with the electrode in the liquid, at a time when the electrode is steady (not within about 1 day of a pulse), and the electrode reading noted at that moment. Two or three samples spread over the run (an early one, a mid-run one and the end) would also test the electrode's gain drift.
- **Pipeline.** The key exists, but today it fills only the "NH4-N start, laboratory" row of the closed-frame balance. Using it to anchor the NH4 model needs a small code change, which would be the next remediation pass once the value exists. Per-sample time stamps would also need a format beyond the current single value.

### 2. Alkalinity

- **Why.** Nothing but NH4Cl was added, so nitrification is capped by the initial alkalinity. The inverted C(t) implies a removal that needs Alk0 ≥ 20.8 meq/L if all of it was nitrified, more than the 1-15 meq/L scan. Either Alk0 was high, or the true NH4 level is higher, or part of the N left by stripping. A measured Alk0 decides between these.
- **How.** Titrate to pH 4.5 with standard acid (Standard Methods 2320 B) and report meq/L.
- **When.** `initial_alkalinity_meq_L` is the starting liquid. If no starting sample was kept, measure the source water or medium prepared the same way, and record that it is a proxy. `alkalinity_mmol_L` is the current liquid, useful at the end of the run, where it bounds the alkalinity actually consumed.

### 3. Nitrate and nitrite at the end

- **Why.** These close the N balance (start + added − end NH4 − NO3 − NO2 = losses). They distinguish nitrification from NH3 stripping through the open surface, which is the alternative sink that S2 leaves open. Nitrite also shows whether nitrification stopped at nitrite.
- **How.** Ion chromatography for both, or colorimetric: sulfanilamide/NED for nitrite, cadmium reduction or a nitrate test kit for nitrate. Filter at 0.45 µm, keep at 4 °C, and analyse within 48 h.
- **When.** At the end of the run, from the same sample as the final NH4 if possible.

### 4. Vessel mass, cover and water added

- **Why.** Conditional on C(t), N conservation needs a net water loss of at least 0.31 mL/d, while liquid present at the end caps it below 2.09 mL/d. A constant 100 mL volume is inconsistent with that. Two weighings settle the evaporation scenario, and with it every concentration balance.
- **How.** Weigh the vessel with its liquid (and stir bar) on a 0.1 g balance.
- **When.** At the start and at the end. Also record whether the vessel was covered (`covered = true/false`) and the total water added during the run (`top_up_mL`, 0 if none). If the start mass was not recorded, an end mass with the empty vessel's mass and the start volume can stand in for it. Note that it is a reconstruction.

### 5. Dissolved oxygen

- **Why.** DO separates "inactive with oxygen to spare" from "oxygen-limited". A logged DO series could turn the idle fraction from bounds (5-93 % of the stirred time) into a value. A single value gives the in-situ thermodynamic driving force.
- **How.** An optical DO probe in the stirred vessel, logging every few minutes, and calibrated in air-saturated water at the run temperature.
- **When.** Ideally throughout a run, which means the next run. A few spot readings during the current run, if any remain, help only partly.
- **Pipeline.** One value is used today, but only for the in-situ affinity. The idle fraction comes from the state classifier, which does not read DO. A logged series, and a code change to make the classifier use it, are both needed before the idle fraction becomes a value.

### 6. O2 transfer coefficient kLa

- **Why.** The aeration-equivalent energy and the O2 utilization are currently SCENARIOs over kLa = 1, 10 and 100 1/d. A measured kLa makes them values. It also bounds NH3 loss through the surface.
- **How.** The dynamic reaeration method in the same vessel, with the same liquid volume and stirring speed. Strip the O2 with N2 (or sulfite with a cobalt catalyst) in clean water or the liquid without biomass. Then log the DO recovery, and fit the slope of ln(C* − C) against time.
- **When.** Any time: it depends on the vessel and the stirring, not on the run. Report it in 1/d.

### 7. Temperature

- **How.** A logged probe in the liquid, or at least the room temperature range over the run.
- **Pipeline.** Needs code. Today the key appears only in the report's text, and every temperature-dependent value uses the 20 and 25 °C SCENARIOs. Reading a measured value, single or logged, is a small change.

### 8. NH4 electrode step response

- **Why.** The electrode responds to pulses with delays of up to 0.33 d and with shifting offsets. A controlled step test shows its true time constant and gain.
- **How.** In a stirred, temperature-controlled solution of known NH4-N and matching ionic background, add a known NH4Cl step and log the electrode until it is steady. Repeat at two concentration levels.
- **Pipeline.** Reported today. Using it in the electrode model needs a code change.

### 9. Chloride and ORP reference offset

- **Why.** The ORP pair uses a bare silver wire as reference, whose potential depends on chloride. Both values are needed to convert ORP into E_H. Low priority: only relative ORP is used in the conclusions.
- **How.** Chloride by ion chromatography or argentometric titration, in mmol/L. For the offset, measure the bare Ag wire against a calibrated Ag/AgCl reference electrode in the same solution, or check the pair in ZoBell solution, in volts.

### 10. Conductivity probe check

- **Why.** The probe showed no response to the pulses, so its record is unusable.
- **How.** Add known amounts of KCl to the probe's solution and compare with a calibrated conductivity meter. A replacement probe may be the practical answer.

### 11. Stirrer power

- **How.** A plug-in wattmeter on the stirrer plate at the run's speed, minus the idle draw, in W. Low priority: the result is specific to this laboratory set-up.

## Already finished batch or next run?

| Can still be done on the finished batch, if the vessel and its liquid were kept | Needs the next run |
| --- | --- |
| End NH4-N (#1, one anchor point), current alkalinity (#2), end NO3/NO2 (#3), end vessel mass (#4), chloride and ORP offset (#9), kLa of the vessel (#6), the electrode step test (#8), the conductivity check (#10), stirrer power (#11) | NH4-N at several times (#1); initial alkalinity, unless the source water is a valid proxy (#2); start vessel mass (#4); a DO log (#5); a temperature log (#7) |

If a next run is planned, recording the starting values (#1, #2, #4) and logging DO and temperature (#5, #7) from the start would resolve almost every remaining NO-GO in one go.
