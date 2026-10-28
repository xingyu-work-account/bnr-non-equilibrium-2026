"""Markdown report helpers."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

SVG_HASHSALT = "bnr_neq"  # fixed salt for SVG element ids, so repeated runs write byte-identical SVG files


def pyplot():
    """matplotlib.pyplot on the non-interactive Agg backend with a fixed SVG hash salt."""
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = SVG_HASHSALT
    import matplotlib.pyplot as plt

    return plt


def fmt(value: object, digits: int = 4) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        return f"{value:.{digits}g}"
    return str(value)


PROVENANCE_COLUMNS = ("nh4_observation_model", "nh4_series", "nh4_truth_id", "assumptions")


def provenance_line(df: pd.DataFrame, series_col: str = "nh4_series") -> str:
    """Pass 6 (G8): one sentence naming a table's NH4 observation model, series and assumptions, for tables whose rows
    all share them (the CSV keeps them per row)."""
    parts = [f"{label}: {', '.join(sorted(set(df[c].astype(str))))}" for label, c in
             (("NH4 observation model", "nh4_observation_model"), ("series", series_col), ("assumptions", "assumptions"))
             if c in df.columns]
    return ("Provenance of every row: " + "; ".join(parts) + ".") if parts else ""


def md_table(df: pd.DataFrame, digits: int = 4) -> str:
    """Render a DataFrame as a GitHub Markdown table (no external dependency)."""
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for row in df.itertuples(index=False):
        cells = [fmt(v.item() if hasattr(v, "item") else v, digits).replace("|", "\\|") for v in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# =============================================================================
# Analysis report (M8): assembled from the generated tables, never from hand-entered numbers
# =============================================================================
TABLE_DESCRIPTIONS = {
    "data_properties.csv": "M1: measured file properties against Section 4",
    "dropped_rows.csv": "M1: every dropped row with file and line",
    "gaps.csv": "M1: gap and missing-sample inventory",
    "grid_coverage.csv": "M1: 0.01 d grid coverage per sensor",
    "sensor_validation.csv": "M2: predicted and observed response to each pulse",
    "sensor_validation_windows.csv": "M2: window medians and sample counts",
    "sensor_balances.csv": "M2: cumulative EC and nitrogen balances",
    "sensor_validation_references.csv": "M2: reference comparison",
    "sensor_verdicts.csv": "M2: sensor verdicts with their a priori criteria",
    "nh4_jumps.csv": "M2: NH4 jump candidates, classification and transition times",
    "nh4_held_segments.csv": "M2: held/damped NH4 output segments",
    "nh4_hf_residual_subsegments.csv": "M2: NH4 high-frequency residual per sub-segment",
    "ec_late_pulses.csv": "M2: late conductivity pulses with decay fits",
    "ec_late_pulse_sensitivity.csv": "M2: late-pulse count against threshold",
    "clock_check.csv": "M2 (corrected in M3): step timing at the day-2.89 multi-probe event",
    "segments.csv": ("M3: sub-segments, Theil-Sen slopes of the reading, held fraction; linear-gain consumption and "
                     "reading-as-TAN free ammonia (SENSITIVITY in the column names since W1/X4); R8 electrode fragments"),
    "segments_window_check.csv": "M3: no slope window contains a pulse or any NH4 step candidate",
    "segments_reference_windows.csv": "M3: slopes on the exact preliminary windows (diagnostic)",
    "segments_references.csv": "M3: reference comparison",
    "kinetic_fits.csv": "M3: zero/first-order, Monod, Haldane fits per segment",
    "kinetic_parameters.csv": "M3: kinetic parameters with Jacobian-based identifiability",
    "changepoints.csv": "M3: PELT breaks and events",
    "changepoint_sensitivity.csv": "M3: PELT break count against penalty",
    "lomb_scargle_summary.csv": "M3: Lomb-Scargle maxima per band",
    "thermo_affinity_steps.csv": "M4: affinity step at each pulse",
    "thermo_coexistence.csv": "M4: coexistence SCENARIO affinity and DO sensitivity",
    "thermo_chloride_scenarios.csv": "M4: bare-Ag chloride SCENARIOS",
    "thermo_chloride_observed.csv": "M4: observed ORP/pH changes and the oxide alternative",
    "thermo_co2.csv": "M4: air-equilibrium bicarbonate",
    "thermo_oxygen.csv": "M4: oxygen demand and required kLa under the linear gain model (SENSITIVITY in the names)",
    "balance_nitrogen.csv": "M4: nitrogen balance",
    "balance_chloride.csv": "M4: chloride balance",
    "balance_alkalinity.csv": "M4: alkalinity balance, with the W2 bound and V(t) scenarios",
    "alkalinity_bound.csv": "W2: upper bound on cumulative nitrification against Alk0 (1-15 meq/L), residual TAN and FA",
    "alkalinity_references.csv": "W2: reference checks of the carbonate and bound functions",
    "volume_scenarios.csv": "W8: vessel volume scenarios V(t), pulse doses and the evaporation upper bound",
    "thermo_references.csv": "M4: reference comparison",
    "estimator_validation.csv": "M5: estimator validation suite, matched-null regimes and gate checks",
    "dynamics.csv": "M5: irreversibility and complexity per series and segment",
    "dynamics_windows.csv": "M5: per-window stationarity, indices and surrogate p-values",
    "dynamics_onsager.csv": "M5: Onsager regression test per pulse",
    "dynamics_onsager_lags.csv": "M5: Onsager ratio against lag",
    "dynamics_stability.csv": "M5: stability indicators before each pulse (exploratory)",
    "dynamics_flux.csv": "M5: probability-flux circulation",
    "dynamics_electrode.csv": "M5: electrode-dynamics guard",
    "model_fit.csv": ("M6/W5: model comparison under one convergence policy: convergence, selection status, "
                      "out-of-sample skill, fitted Alk0 and its profile interval, skill range over the starts inside "
                      "Δ < 3.84, free parameters at a bound (`cal_at_bounds`)"),
    "model_calibration_starts.csv": ("X1/Y1: every start of the convergence policy per candidate, its Δ against the "
                                     "calibration and its selection, test and forecast skill (diagnostic), free "
                                     "parameters at a bound"),
    "model_identifiability.csv": "M6: profile-likelihood identifiability (with one-sided crossings)",
    "model_volume_scenarios.csv": ("Y7: SCENARIO refit of every candidate per evaporation scenario (convergence, skills, "
                                   "selection); never feeds a primary row"),
    "model_volume_scenario_summary.csv": ("Y7: per evaporation SCENARIO, the selected model and the primary model's "
                                          "refit and C(t) level against the constant-volume band, its ratio to the "
                                          "constant-volume C(t) against φ = V0/V(t); never a primary row"),
    "model_checks.csv": "M6: mass-balance and σ consistency checks",
    "model_response_fractions.csv": "M6: model-implied against observed response fractions",
    "nh4_truth_summary.csv": "W1: the single NH4 observation model, its inverted C(t) and band",
    "nh4_truth_vectors.csv": "W1: profile parameter vectors whose inversions form the C(t) band",
    "nh4_truth_windows.csv": "W1: consumption, free ammonia and O2 demand of the inverted C(t) per slope window",
    "nh4_truth_balance.csv": "W1: dose-anchored TAN bound (with the W2 bound) against the inverted C(t)",
    "nh4_truth_w2_consistency.csv": "S2: per checkpoint, the Alk0 the C(t)-implied removal needs, and the verdict",
    "volume_constraints.csv": "W8: data constraints on the net water loss and the evaporation SCENARIOS checked",
    "states_rules.csv": ("W3: the state rules, thresholds and excluded sensor features, with their registration hash "
                         "and status (revised after first results, not pre-registered)"),
    "states_fractions.csv": ("W3: time fraction of each state; idle fraction as bounds with its coverage, and the "
                             "CONDITIONAL value over classified time where coverage ≥ 0.5"),
    "states_timeline.csv": "W3: time-ordered runs of one process state",
    "states_comparison.csv": ("W3: NH4 consumption between indicator groups, block-bootstrap intervals, and their "
                              "reading against the synthetic null (G2)"),
    "states_bootstrap_null_coverage.csv": ("G2: synthetic null check of the block bootstrap with the record's group "
                                           "labels (coverage of the true difference 0 per φ)"),
    "states_interval_coverage.csv": ("Y3: synthetic AR(1) check of the autocorrelation-scaled slope interval (its "
                                     "coverage, i.e. the limit of the idle bounds)"),
    "pulses.csv": ("M7/W4: pulse-level metrics (η, removal, O2, energy family, idle-fraction bounds and coverage) with "
                   "inputs and status"),
    "pulses_associations.csv": "M7: exploratory associations",
    "all_references.csv": "M8: every reference value of the implementation prompt",
}

FIGURES = [
    ("fig1_linkage.png", "Figure 1. From sensor signals to process state, oxygen demand, energy and efficiency (W6)"),
    ("fig2_overview.png", "Figure 2. Common-window overview with pulses, flags and electrode events (M3)"),
    ("fig3_sensor_validation.png", "Figure 3. Sensor validation with the known pulses (M2)"),
    ("fig5_dynamics.png", "Figure 5. Estimator validation, irreversibility, Onsager SNR, probability flux (M5)"),
    ("fig6_model.png", "Figure 6. Reduced model with the electrode observation model (M6)"),
    ("fig7_pulses.png", "Figure 7. Pulse-level descriptors (M7)"),
    ("fig4_constraints.png", "Figure 4. Alkalinity bound, rate against free ammonia, state timeline (W6)"),
    ("fig8_efficiency.png", "Figure 8. Per-pulse η, oxygen utilization, idle fraction and energy under SCENARIOS (W6)"),
]

def assumptions(cfg, t: dict[str, pd.DataFrame]) -> list[str]:
    """Stated assumptions; every number comes from the configuration or a generated table."""
    from . import dynamics, segments, states

    nf = cfg.get("sensors.nh4_to_N_factor")
    temps = ", ".join(f"{x:g}" for x in cfg.get("scenarios.temperature_C"))
    ident = t["model_identifiability.csv"]
    n_off = int(ident["parameter"].str.startswith("offset_log_").sum())
    return [
        ("File headers differ in spelling from implementation_prompt.md Section 4 for Cond.csv, NH4.csv and ORP.csv; "
        "the loader maps the headers present in the files and raises on any other (deviation D1)."),
        ("NH4 and ORP timestamps are rounded to 0.01 d; rows sharing a stamp are uniformly spaced within it "
        "(t = stamp + 0.01 j/n). The half-sized first stamps suggest round-to-nearest stamps, which would place these "
        "times 0.005 d late; not corrected."),
        f"NH4 values are in {cfg.get('sensors.nh4_unit')} (experimenter); converted to mg N/L with {nf:g}.",
        (f"Conductivity is converted with EC = {cfg.get('sensors.ec_slope_mS_cm_per_A'):g} I + "
        f"{cfg.get('sensors.ec_intercept_mS_cm'):g} and is not valid as solution conductivity (S3)."),
        (f"Temperature is not measured: every temperature-dependent value is a SCENARIO value, at the labelled SCENARIO "
        f"temperatures {temps} °C; values computed with 25 °C constants carry `SCENARIO_T25C` or a `SCENARIO 25 °C` "
        "basis; ΔG°′ is a 25 °C value."),
        "Kohlrausch predictions use limiting ionic conductivities at 25 °C (SCENARIO; upper bounds at the actual ionic strength).",
        "Bare-Ag chloride SCENARIOS use chloride activity = concentration (γ = 1) and equal γ before and after a pulse.",
        "C* = 8.26 mg/L (SCENARIO 25 °C) for oxygen; p_O2 = 0.2095 x DO/C*. Ω_CO2 would take [HCO3-] ≈ alkalinity (pH < 8.3).",
        ("Single NH4 observation model (W1): rates, free ammonia, η, ΔA and O2 demand are derived from C(t) inverted "
        "from the M6 log observation model (quasi-static; lag not inverted); a non-positive consumption is not "
        "computed and is excluded from duration-weighted rates. Gain-corrected rates (division by the response "
        "fraction, linear electrode gain model) and free ammonia from the reading taken as TAN are SENSITIVITY cases."),
        ("Alkalinity bound (W2): 2 eq per mol N nitrified, no alkalinity source, carbonate equilibria (Plummer and "
        "Busenberg 1982) with activities = concentrations, dissolved CO2 at or above air saturation, pH from the "
        "0.11 d rolling median; kLa_CO2 and Alk0 are SCENARIO ranges."),
        ("Volume (W8): the closed 100 mL frame is an assumption until the vessel masses are measured; evaporation "
        f"SCENARIOS {', '.join(f'{e:g}' for e in cfg.get('scenarios.evaporation_mL_per_d'))} mL/d are carried "
        "through the balances, the η and the consumption rates (r = −dC/dt + C E/V)."),
        (f"Slope windows exclude {segments.EXCLUDE_AFTER_BOUNDARY_D} d after every pulse and every NH4 step candidate "
        "(located or not, sustained or not) and need ≥ "
        f"{segments.MIN_FIT_SPAN_D} d; held/damped NH4 spans inside a window are flagged (`held_fraction`), with a "
        "sensitivity slope that removes them."),
        ("Stationary estimators run only on raw 0.25 d windows passing ADF and KPSS; an index is computed only within "
        "its validated tie range and, for pH, only inside matched null regimes whose false-positive rate is at most "
        f"{dynamics.FP_ALPHA} + binomial margin ({dynamics.N_SURR} IAAFT surrogates)."),
        (f"The M6 model fixes K_NH = 1.0 mg N/L (ASM1), uses offsets at all {n_off} detected electrode-offset jumps, and "
        "refits the test-window offsets as the only nuisance parameters when scoring the test window; the split "
        "windows come from the configured pulse times."),
        ("MODEL-BASED quantities (η, model rates) are conditional on the selected M6 model; C0 and the electrode "
        "intercept are not identifiable."),
        ("Process states (W3): rules revised after their first results and frozen by SHA-256 from pass 6 "
         f"(`{states.REGISTRATION_DOC}`; not pre-registered) on the 1.01 d Theil-Sen slopes of C(t), pH and ORP and on ORP "
         f"relative to its centred 5 d median (±15 mV, {states.ORP_THRESHOLD_SCENARIO}); the slope intervals are widened by √(n/n_eff) for the lag-1 "
         "autocorrelation of each window's residuals; relative ORP assumes a slowly varying bare-Ag reference; QC-flagged "
         "bins carry no state. The idle fraction is reported as bounds with the classified fraction (coverage)."),
        ("Energy (W4): kLa_O2 and the field (process-water) aeration efficiency AE are SCENARIO ranges (AE source: "
         f"{cfg.get('scenarios.sae_source')}); the O2 supply is a capacity at DO = 0 while DO is not measured, whereas "
         "the AE values are field values at the plant's operating DO (> 0; not stated in the source table); stirrer power "
         "is a laboratory-specific term only."),
        ("Model candidates (W5): alkalinity limitation uses the model's own alkalinity state (2 eq/mol N, nothing added) "
         "with K_ALK = 0.5 meq/L fixed (ASM3)."),
        ("Electrode fragments (R8): a slope window between two electrode steps (offset jumps, pulse-response jumps or "
         "step candidates) shorter than 3 x the longest delayed electrode response is a diagnostic, not a rate, whether "
         "or not a pulse governs it; the pulse rates with the fragments are an R8 CASE in `pulses.csv`."),
    ]


class ReportValidationError(RuntimeError):
    """A required validation check failed; the report is not written."""


def check_required(t: dict[str, pd.DataFrame]) -> None:
    """Raise unless every required estimator check and every M6 consistency check passes, and the pH
    matched-null gate re-checks against the stored tables."""
    from . import dynamics

    v = t["estimator_validation.csv"]
    failed = dynamics.required_failures(v)
    if not v["required"].astype(str).str.lower().eq("true").any():
        raise ReportValidationError("estimator_validation.csv has no required check")
    gate = dynamics.verify_null_gate(t["dynamics_windows.csv"], v)
    failed += [c for c, ok in zip(gate["check"], gate["pass"]) if not ok]
    checks = t["model_checks.csv"]
    failed += [f"model: {c}" for c, ok in zip(checks["check"], checks["pass"].astype(str)) if ok == "False"]
    if failed:
        raise ReportValidationError(f"required validation checks failed: {failed}")


def _read(tables: Path, name: str) -> pd.DataFrame:
    p = tables / name
    if not p.is_file():
        raise FileNotFoundError(f"{p} missing: run `python -m bnr_neq run-all` stages before report")
    return pd.read_csv(p)


def all_references(tables: Path) -> pd.DataFrame:
    m2 = _read(tables, "sensor_validation_references.csv")
    m2 = m2.assign(milestone="M2", quantity=m2["quantity"] + m2["pulse"].map(lambda p: "" if pd.isna(p) else f" (pulse {int(p)})"))
    m3 = _read(tables, "segments_references.csv").assign(milestone="M3")
    m4 = _read(tables, "thermo_references.csv").assign(milestone="M4")
    cols = ["milestone", "quantity", "reference", "obtained", "tolerance", "pass", "note"]
    return pd.concat([m2[cols], m3[cols], m4[cols]], ignore_index=True)


CONSISTENCY = "CONSISTENCY CHECK (not a finding)"


def _pm(t: dict[str, pd.DataFrame], metric: str) -> pd.DataFrame:
    pm = t["pulses.csv"]
    return pm[(pm["metric"] == metric)]


def primary_delta_a(t: dict[str, pd.DataFrame]) -> tuple[pd.Series, pd.Series]:
    """ΔA at 25 °C on the single NH4 observation model (pulses.csv, M6_log), and the reading-as-TAN SENSITIVITY."""
    from .model import OBSERVATION_MODEL

    prim = _pm(t, "delta_A_kJ_mol_SCENARIO_T25C")
    prim = prim[prim["nh4_observation_model"] == OBSERVATION_MODEL]
    sens = _pm(t, "delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T25C")
    return prim["value"].astype(float), sens["value"].astype(float)


def status_labels(t: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Status of each conclusion; every evidence sentence is derived from the generated tables (values rounded)."""
    from . import dynamics
    from .constants import DG0_NITRIFICATION_KJ_MOL_N

    coex = t["thermo_coexistence.csv"]
    scen = coex[coex["basis"] == "SCENARIO"].copy()
    scen["affinity_kJ_mol"] = scen["affinity_kJ_mol"].astype(float)
    doc = float(scen.loc[scen["documented_scenario"].astype(bool), "affinity_kJ_mol"].iloc[0])
    insitu = coex.loc[coex["basis"] == "IN SITU", "affinity_kJ_mol"].astype(str)
    insitu_txt = insitu.iloc[0] if len(insitu) and insitu.iloc[0].startswith("not computed") else "computed"
    dA, dA_sens = primary_delta_a(t)
    dyn = t["dynamics.csv"]
    ph = dyn[dyn["series"] == "pH"]
    ons = t["dynamics_onsager.csv"]
    pv = t["sensor_validation.csv"]
    ident = t["model_identifiability.csv"]
    cl = t["thermo_chloride_scenarios.csv"]
    cl_cols = [c for c in cl.columns if c.startswith("orp_pred_rise_mV_pulse")]
    cl_lo, cl_hi = float(cl[cl_cols].min().min()), float(cl[cl_cols].max().max())
    val = t["estimator_validation.csv"]
    req = val[val["required"].astype(str).str.lower().eq("true")]
    n_req_ok = int(req["pass"].astype(str).str.lower().eq("true").sum())
    win = t["dynamics_windows.csv"]
    phw = win[(win["series"] == "pH") & (win["transform"] == "raw")]
    st = phw["ordinal_kl_status"].fillna("") if "ordinal_kl_status" in phw else pd.Series(dtype=str)
    n_out = int(st.str.startswith(dynamics.OUTSIDE_REGIME).sum())
    sig = pd.to_numeric(ph["ordinal_kl_significant_windows"], errors="coerce").fillna(0).sum()
    snr_ok = int(ons["snr_sufficient"].astype(str).str.lower().eq("true").sum())
    f = pv["nh4_response_fraction"]
    ec_pred = pv["ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C"].iloc[0]
    checks = t["model_checks.csv"]
    sigma = checks[checks["check"].astype(str).str.startswith("σ = r A / T")]
    sigma_txt = (f"min σ = {float(sigma['value'].iloc[0]):.2g} W/(L K)" if len(sigma) else "not in model_checks.csv")
    cons = t["nh4_truth_w2_consistency.csv"]
    bad = cons[~cons["consistent_at_best_fit"].astype(bool)]
    bnd = t["alkalinity_bound.csv"]
    from . import alkalinity

    h = alkalinity.headline(bnd)
    n15 = float(h.loc[h["alk0_meq_L"] == h["alk0_meq_L"].max(), "n_nit_max_end_mgN_L"].iloc[0])
    dosed = float(bnd["dosed_total_mgN_L"].iloc[0])
    fr = t["states_fractions.csv"].set_index("interval").loc["whole record"]
    sc = t["states_comparison.csv"]
    ph_test = sc[sc["comparison"].astype(str).str.startswith("pH rising vs not rising (within")]
    ph_txt = (f"pH-rising minus non-rising bins {float(ph_test['diff'].iloc[0]):+.2f} mg N/L/d "
              f"(95 % block bootstrap {float(ph_test['ci95_low'].iloc[0]):+.2f} to {float(ph_test['ci95_high'].iloc[0]):+.2f})"
              if len(ph_test) and ph_test["status"].iloc[0] == "computed" else "pH test not computed")
    sae = ae_levels(t["pulses.csv"])
    fit = t["model_fit.csv"]
    conv = fit["cal_converged"].astype(str).str.lower().eq("true")
    sel = fit.loc[fit["selected"].astype(str).str.lower().eq("true"), "model"]
    alk = fit[fit["model"].astype(str).str.startswith("alkalinity")]
    vc = t["volume_constraints.csv"]
    e_min = vc.loc[vc["constraint"] == "N conservation over every pulse", "value_mL_per_d"]
    rows = [
        ("Upper bound on cumulative nitrification from the pH record (W2)", "GO",
         (f"With Alk0 from 1 to 15 meq/L at most {n15:.0f} of {dosed:.0f} mg N/L dosed could be nitrified by day 47.9 "
          "(air-equilibrium CO2, constant volume, SCENARIO 25 °C)."),
         "Assumes no alkalinity source; initial_alkalinity_meq_L missing (scan instead of a value)."),
        ("Single NH4 observation model against the W2 bound (S2)", "QUALIFIED" if len(bad) else "GO",
         (f"The removal implied by the inverted C(t) needs Alk0 ≥ {float(cons['alk0_required_if_nitrified_meq_L'].max()):.1f} "
          f"meq/L if nitrified; contradicted at {', '.join(bad['checkpoint'])}" if len(bad) else
          "The removal implied by C(t) is compatible with the scanned Alk0.")
         + " Every C(t)-derived row carries `w2_bound_check`.",
         "Blocked by: lab_NH4_N_mg_L, initial_alkalinity_meq_L, NO3_N_mg_L (resolve level, Alk0, sink)."),
        ("Process states from relative ORP, pH and C(t) (W3)", "GO WITH DOCUMENTED NON-BLOCKING DEBT",
         (f"{100 * fr['classified_fraction']:.0f} % of the record carries a state (autocorrelation-scaled slope "
          f"intervals); {ph_txt}."),
         "Debt: low coverage (held NH4 output, electrode steps); relative ORP assumes a slowly varying reference."),
        ("Idle fraction of stirred time (W3 into W4)", "BOUNDS ONLY; NO-GO for a single value",
         (f"Between {100 * fr['idle_fraction_lower_bound']:.0f} and {100 * fr['idle_fraction_upper_bound']:.0f} % of "
          f"the stirred time (whole record) at a coverage of {100 * fr['classified_fraction']:.0f} % classified; the "
          "CONDITIONAL value over classified time: "
          + (f"{100 * fr['idle_fraction_of_classified_time']:.0f} %" if pd.notna(fr["idle_fraction_of_classified_time"])
             else str(fr["status"]).split(" (")[0]) + "."),
         "Blocked by: coverage (held NH4 output, electrode steps), DO_mg_L."),
        ("Model selection, activity decay against alkalinity limitation (M6, W5)",
         "GO WITH DOCUMENTED NON-BLOCKING DEBT" if conv.all() else "QUALIFIED: unconverged candidates excluded",
         (f"{int(conv.sum())} of {len(fit)} calibrations converged under one policy; selected {', '.join(sel)}"
          + (f"; {alk0_evidence_txt(fit)}; S2 requirement "
             f"{float(cons['alk0_required_if_nitrified_meq_L'].max()):.1f} meq/L" if len(alk) else "") + "."),
         "Debt: four pulses; time is confounded with decay, alkalinity loss, FA and electrode drift; Alk0 not measured."),
        ("Oxygen demand and aeration-equivalent energy (W4)", "SCENARIO FAMILY COMPUTED; NO-GO for absolute values",
         (f"O2 demand 4.57 g O2/g N on the removal of C(t); field AE {min(sae):g}-{max(sae):g} kg O2/kWh and kLa_O2 "
          f"SCENARIOS; demand floor {4.57 / max(sae):.2f}-{4.57 / min(sae):.2f} kWh/kg N."),
         "Blocked by: DO_mg_L, kLa_per_d, temperature_C_measured; QUALIFIED by S2 from pulse 2 on."),
        ("Closed 100 mL frame (W8)", "NO-GO",
         (f"Conditional on C(t), N conservation needs a net water loss ≥ {float(e_min.iloc[0]):.2f} mL/d"
          if len(e_min) else "not computed") + "; evaporation SCENARIOS carried through every balance; M6 refitted per "
         "SCENARIO: " + volume_status_txt(t.get("model_volume_scenario_summary.csv", pd.DataFrame())) + ".",
         "Blocked by: vessel_mass_start_g, vessel_mass_end_g, covered, top_up_mL."),
        ("Chemical disequilibrium of the batch (coexistence affinity)", CONSISTENCY,
         (f"Coexistence SCENARIO affinity {doc:.0f} kJ/mol ({scen['affinity_kJ_mol'].min():.0f} to "
          f"{scen['affinity_kJ_mol'].max():.0f} kJ/mol for DO {scen['DO_mg_L'].min():g} to {scen['DO_mg_L'].max():g} mg/L): "
          "computed from assumed concentrations, positive by construction whenever DO > 0."),
         f"In situ affinity {insitu_txt} (needs DO_mg_L, NO3_N_mg_L)."),
        ("Affinity steps at the pulses", CONSISTENCY,
         (f"ΔA = {dA.min():.2f} to {dA.max():.2f} kJ/mol on the M6-inverted C(t) (SCENARIO 25 °C; "
          f"{100 * dA.min() / abs(DG0_NITRIFICATION_KJ_MOL_N):.2f} to {100 * dA.max() / abs(DG0_NITRIFICATION_KJ_MOL_N):.2f} "
          f"% of |ΔG°′|); SENSITIVITY reading as TAN {dA_sens.min():.2f} to {dA_sens.max():.2f} kJ/mol. A log ratio of "
          "concentrations: small by construction for a dose small against the pool."),
         "None blocking."),
        ("Entropy production σ = r A / T ≥ 0 in the model", CONSISTENCY,
         f"{sigma_txt}: a non-negative modelled rate times a positive SCENARIO affinity, true by construction.",
         "None blocking."),
        ("Relative comparison of pulse responses (electrode)", "GO WITH DOCUMENTED NON-BLOCKING DEBT",
         f"Response fractions {', '.join(f'{100 * x:.0f} %' for x in f)} of the dose.",
         "Debt: gain unverified (lab_NH4_N_mg_L missing); delayed responses unexplained (electrode_step_response)."),
        ("Absolute NH4 concentrations, rates and conversion", "NO-GO",
         (f"Electrode registers {100 * f.min():.0f} to {100 * f.max():.0f} % of the dose; C0 and the intercept are not "
          f"identifiable ({int((~ident['identifiable'].astype(bool)).sum())} parameters not identifiable)."),
         "Blocked by: lab_NH4_N_mg_L, electrode_step_response, NO3_N_mg_L, NO2_N_mg_L (and recalibration)."),
        ("Conductivity as solution conductivity", "NO-GO",
         (f"Observed pulse responses {pv['ec_obs_change_mS_cm'].min():+.3f} to {pv['ec_obs_change_mS_cm'].max():+.3f} mS/cm "
          f"against +{ec_pred:.3f} predicted (limiting, SCENARIO 25 °C)."),
         "Blocked by: probe validation (KCl standard addition, calibrated meter), chloride_mmol_L."),
        ("Absolute redox potential E_H", "NO-GO",
         (f"ORP response to the chloride added {pv['orp_obs_change_mV'].min():+.1f} to {pv['orp_obs_change_mV'].max():+.1f} mV "
          f"against {cl_lo:.0f} to {cl_hi:.0f} mV per pulse under chloride SCENARIOS."),
         "Blocked by: orp_reference_offset_V, chloride_mmol_L."),
        ("pH trajectory estimators (supplement S1)", "GO WITH DOCUMENTED NON-BLOCKING DEBT",
         (f"{int(ph['n_stationary'].sum())} of {int(ph['n_windows'].sum())} pH windows stationary; ordinal KL in "
          f"{int(pd.to_numeric(ph['ordinal_kl_windows'], errors='coerce').fillna(0).sum())} validated windows "
          f"({n_out} outside validated regimes), significant in {int(sig)}; {n_req_ok} of {len(req)} required checks pass."),
         (f"Debt: not connected to rates or energy; Onsager test interpretable for {snr_ok} of {len(ons)} pulses.")),
    ]
    return pd.DataFrame(rows, columns=["conclusion", "status", "evidence", "missing items / debt"])


def ae_levels(pm: pd.DataFrame) -> list[float]:
    """The field aeration-efficiency (AE) SCENARIO levels present in pulses.csv."""
    return sorted({float(m.split("_AE")[1].split("_")[0]) for m in pm["metric"]
                   if str(m).startswith("aeration_energy_demand_kWh_per_kgN_SCENARIO_AE")})


def energy_table(pm: pd.DataFrame) -> pd.DataFrame:
    """Compact per-pulse view of the W4 metric family (constant volume; SCENARIO ranges in the column names)."""
    def get(metric):
        r = pm[pm["metric"] == metric].set_index("pulse")
        return r["value"].where(r["status"] == "computed")

    klas = sorted({float(m.split("_kLa")[1].split("_")[0]) for m in pm["metric"] if str(m).startswith("o2_utilization_SCENARIO_kLa")})
    saes = ae_levels(pm)
    mid = saes[len(saes) // 2]
    out = pd.DataFrame({"pulse": sorted(pm["pulse"].unique())}).set_index("pulse")
    out["η"] = get("eta_conversion_fraction")
    out["removal / dose"] = get("removal_to_dose_ratio")
    out["N removed (mg N)"] = get("n_removed_mgN")
    out["O2 required (mg O2)"] = get("o2_required_mgO2")
    for k in klas:
        out[f"O2 utilization, kLa {k:g}/d"] = get(f"o2_utilization_SCENARIO_kLa{k:g}_T25C")
    out[f"kWh/kg N floor, AE {max(saes):g}-{min(saes):g}"] = [
        f"{a:.2f}-{b:.2f}" if pd.notna(a) else "not computed"
        for a, b in zip(get(f"aeration_energy_demand_kWh_per_kgN_SCENARIO_AE{max(saes):g}"),
                        get(f"aeration_energy_demand_kWh_per_kgN_SCENARIO_AE{min(saes):g}"))]
    for k in klas:
        out[f"kWh/kg N supply, AE {mid:g}, kLa {k:g}/d"] = get(f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{mid:g}_kLa{k:g}_T25C")
    lo_, hi_, cov = get("idle_fraction_LOWER_BOUND"), get("idle_fraction_UPPER_BOUND"), get("classified_fraction")
    out["idle fraction: bounds (coverage)"] = [
        f"{a:.2f}-{b:.2f} ({c:.2f})" if pd.notna(a) and pd.notna(b) else "not computed" for a, b, c in zip(lo_, hi_, cov)]
    out["dose-removal time / pulse interval"] = get("dose_removal_time_over_pulse_interval")
    k10 = klas[len(klas) // 2]
    out[f"ΔG released / aeration energy, AE {mid:g}, kLa {k10:g}/d"] = get(
        f"gibbs_to_aeration_energy_supply_SCENARIO_AE{mid:g}_kLa{k10:g}_T25C")
    out["W2 check"] = pm.groupby("pulse")["w2_bound_check"].agg(
        lambda s: "QUALIFIED" if s.astype(str).str.startswith("QUALIFIED").any() else "consistent")
    return out.astype(object).where(out.notna(), "n.c.").reset_index()  # n.c.: not computed, reason in pulses.csv


def summary_lines(cfg, t: dict[str, pd.DataFrame]) -> list[str]:
    """The one-page process summary; every number from a generated table."""
    from . import alkalinity, model, states

    pm = t["pulses.csv"]
    v = t["sensor_verdicts.csv"].set_index("key")["verdict"].astype(str)
    ts = t["nh4_truth_summary.csv"].set_index("key")["value"]
    cons = t["nh4_truth_w2_consistency.csv"].set_index("checkpoint")
    h = alkalinity.headline(t["alkalinity_bound.csv"])
    dosed = float(t["alkalinity_bound.csv"]["dosed_total_mgN_L"].iloc[0])
    n15 = float(h.loc[h["alk0_meq_L"] == h["alk0_meq_L"].max(), "n_nit_max_end_mgN_L"].iloc[0])
    n5 = float(h.loc[h["alk0_meq_L"] == 5.0, "n_nit_max_end_mgN_L"].iloc[0]) if (h["alk0_meq_L"] == 5.0).any() else np.nan
    fa_lo = float(h["fa_residual_min_end_mgNH3_L"].min())
    fa_hi = float(h["fa_residual_min_end_mgNH3_L"].max())
    end = cons.loc["end"]

    def val(metric):
        r = pm[pm["metric"] == metric].set_index("pulse")
        return {int(k): (float(x) if s == "computed" else None) for k, x, s in zip(r.index, r["value"], r["status"])}

    eta, relax = val("eta_conversion_fraction"), val("dose_removal_time_over_pulse_interval")
    eta_mb, gap = val("eta_conversion_fraction_MODEL_BASED"), val("eta_gap_inversion_misfit_MODEL_BASED")
    ratio = val("removal_to_dose_ratio")
    qualified = sorted({int(k) for k, q in zip(pm["pulse"], pm["w2_bound_check"].fillna("").astype(str)) if q.startswith("QUALIFIED")})
    saes = sorted(float(x) for x in cfg.get("scenarios.sae_field_kgO2_per_kWh"))
    klas = sorted(float(x) for x in cfg.get("scenarios.kLa_O2_per_d"))
    s_mid, k_mid = saes[len(saes) // 2], klas[len(klas) // 2]
    sup = val(f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{s_mid:g}_kLa{k_mid:g}_T25C")
    idle_lo, idle_hi, cov = val("idle_fraction_LOWER_BOUND"), val("idle_fraction_UPPER_BOUND"), val("classified_fraction")
    idle_txt = ", ".join(f"P{k} " + (f"{idle_lo[k]:.2f}-{idle_hi[k]:.2f} (coverage {cov[k]:.2f})"
                                     if idle_lo.get(k) is not None and idle_hi.get(k) is not None else "n.c.")
                         for k in idle_lo)
    fmt_d = lambda dct, f="{:.2f}": ", ".join(f"P{k} " + (f.format(x) if x is not None else "n.c.") for k, x in dct.items())
    fr = t["states_fractions.csv"].set_index("interval")
    whole = fr.loc["whole record"]
    sc = t["states_comparison.csv"].set_index("comparison")
    fit = t["model_fit.csv"]
    sel = fit.loc[fit["selected"].astype(bool)].iloc[0]
    conv = fit["cal_converged"].astype(str).str.lower().eq("true")
    open_up = str(ts.get("band_high_open", "False")) == "True"

    alk_cmp = alk0_comparison(fit, t["nh4_truth_w2_consistency.csv"])
    cal_end = [c for c in alk_cmp.columns if c.startswith("S2 requirement at pre_pulse")]
    cal_end = cal_end[0].removeprefix("S2 requirement at ").removesuffix(" (meq/L)") if cal_end else "end"

    vc = t["volume_constraints.csv"].set_index("constraint")
    e_min = float(vc.loc["N conservation over every pulse", "value_mL_per_d"])
    e_max = float(vc.loc[vc.index.str.startswith("liquid present"), "value_mL_per_d"].iloc[0])

    def cmp_txt(name):
        r = sc.loc[name]
        if r["status"] == "computed":
            under = str(r.get("interval_reading", "")).startswith("UNDER-COVERS")
            return (f"{r['diff']:+.2f} mg N/L/d (95 % block bootstrap {r['ci95_low']:+.2f} to {r['ci95_high']:+.2f}"
                    + (f"; the interval under-covers on a synthetic null, coverage down to {r['null_coverage_min']:.2f}"
                       if under else "") + ")")
        return f"{r['diff']:+.2f} mg N/L/d, interval not computed ({str(r['status']).removeprefix('not computed: ')})"

    lines = [
        ("1. **Only pH is a validated sensor.** pH: " + v["ph"] + f"; NH4 electrode: {v['nh4']}; conductivity: {v['ec']}; "
         f"ORP: {v['orp']} (absolute E_H unusable, relative changes used with a stated assumption)."),
        (f"2. **One NH4 truth.** Every NH4-derived number rests on C(t), the inversion of the M6 log observation model. Its "
         f"level is not identifiable: {float(ts['C_inv_start_median_mgN_L']):.0f} mg N/L at day 1 (band "
         f"{float(ts['C_inv_start_band_low_mgN_L']):.0f}-{float(ts['C_inv_start_band_high_mgN_L']):.0f}) and "
         f"{float(ts['C_inv_end_median_mgN_L']):.0f} at day 47.9 (band {float(ts['C_inv_end_band_low_mgN_L']):.0f}-"
         f"{float(ts['C_inv_end_band_high_mgN_L']):.0f})"
         + (", and the band is open upward over the whole series (the data do not bound the level from above; its "
            "upper edge is only the last evaluated profile vector)" if open_up else "")
         + f", against an electrode reading of {float(ts['reading_end_median_mgN_L']):.0f} mg N/L at the end."),
        (f"3. **Alkalinity caps nitrification.** Nothing but NH4Cl was added, so nitrification (2 eq/mol N) is bounded by "
         f"the initial alkalinity and the pH record: by day 47.9 at most {n5:.0f} mg N/L for Alk0 = 5 meq/L and "
         f"{n15:.0f} mg N/L for 15 meq/L, against {dosed:.0f} mg N/L dosed. Unless Alk0 was high, most of the dose stayed "
         f"as TAN, with free ammonia of at least {fa_lo:.1f}-{fa_hi:.1f} mg NH3/L at the end (SCENARIO 25 °C)."),
        (f"4. **The inverted C(t) and the bound disagree (S2).** The best-fit C(t) implies "
         f"{end['removed_per_initial_litre_mgN_L']:.0f} mg N/L removed by the end, which needs Alk0 ≥ "
         f"{end['alk0_required_if_nitrified_meq_L']:.1f} meq/L if nitrified (scan 1-15). It is consistent only if the "
         f"true level lies in the open upper part of the band (≥ {end['alk0_required_if_nitrified_band_min_meq_L']:.1f} "
         f"meq/L there), if Alk0 exceeds 15 meq/L, or if a sink releasing 1 eq/mol N (NH3 stripping from the open, "
         f"stirred vessel; ≥ {end['alk0_required_any_proton_sink_meq_L']:.1f} meq/L) took part. Every C(t)-derived "
         "number from pulse 2 on is therefore QUALIFIED (`w2_bound_check`)."),
        (f"5. **Conversion per pulse and dose-removal time.** η = {fmt_d(eta)} (fraction of each dose removed before "
         "the next; n.c. = not computed, e.g. a negative removal on C(t) or a removal above the dose). The time to "
         "remove one dose at the post-pulse rate (a ratio of observed quantities, not a fitted relaxation constant), "
         f"relative to the pulse interval, was {fmt_d(relax, '{:.1f}')}"
         + (("; every dose arrived before the previous one could be removed at that rate, so, on the QUALIFIED C(t) "
             f"(pulses {', '.join(str(k) for k in qualified)} contradict the W2 bound, point 4), TAN accumulated and "
             "the batch was driven away from its inactive reference." if qualified else
             "; every dose arrived before the previous one could be removed at that rate, so TAN accumulated and the "
             "batch was driven away from its inactive reference.")
            if relax and all(x is not None and x > 1 for x in relax.values()) else ".")
         + f" The MODEL-BASED η of the selected rate law is {fmt_d(eta_mb, '{:.3f}')}, against a dose-anchored "
         f"removal per dose of {fmt_d(ratio)}. Both use the same observation model, and the model state obeys its own "
         "mass balance, so the difference is the decline over each interval of the misfit between the inverted C(t) "
         f"and the simulated state, per dose: {fmt_d(gap)} (`eta_gap_inversion_misfit_MODEL_BASED`). The rate law "
         "does not follow the fast post-pulse fall of the reading (its residuals are strongly autocorrelated), and "
         "the inversion assigns that misfit to removal; neither value is the better estimate while C0 and the "
         "intercept are not identifiable."),
        (f"6. **Process state.** The registered rules (revised after first results, not pre-registered) give a state to {100 * whole['classified_fraction']:.0f} % of the "
         f"record ({100 * whole['fraction_active_nitrification']:.0f} % active nitrification, "
         f"{100 * whole['fraction_inactive_oversupplied']:.0f} % inactive and over-supplied, "
         f"{100 * whole['fraction_oxygen_limited']:.0f} % oxygen-limited); the rest is sensor artifact "
         f"({100 * whole['fraction_sensor_artifact']:.0f} %, mostly held NH4 output) or unclassified "
         f"({100 * whole['fraction_unclassified']:.0f} %, much of it a C(t) that rises without an N source). Between "
         f"NH4-free indicator groups, NH4 consumption differs by: pH rising vs not {cmp_txt('pH rising vs not rising (within pulse intervals)')}; "
         f"ORP low vs high {cmp_txt('ORP_rel low vs high (within pulse intervals)')}. No indicator separates the "
         "consumption reliably; these are associations, not causes."),
        # Z4: the conclusions only; the per-candidate ranges, rankings and Alk0 minima are in the §4 table
        (f"7. **Mechanism.** {int(conv.sum())} of {len(fit)} calibrations converged under one policy"
         + (f" ({', '.join(fit.loc[~conv, 'model'])} did not and are excluded from selection)" if (~conv).any() else "")
         + f". The split selects **{sel['model']}** (selection skill {float(sel['sel_skill_vs_null']):.3f}, test "
         f"{float(sel['test_skill_vs_null']):.3f}, forecast {float(sel['forecast_skill_vs_null']):.3f}), and "
         + mechanism_summary_txt(fit, sel, other_mechanism(fit))
         + " " + alk0_sentence(fit, alk_cmp, float(cons.loc[cal_end, "alk0_required_if_nitrified_meq_L"]), cal_end)
         + " With four pulses, time is confounded with activity decay, alkalinity loss, FA and electrode drift "
         "(ranges, rankings and Alk0 profiles: §4)."),
        (f"8. **Energy and efficiency.** The transferable energy measure is oxygen: 4.57 g O2 per g N removed, i.e. a "
         f"floor of {4.57 / max(saes):.2f}-{4.57 / min(saes):.2f} kWh/kg N at field aeration efficiencies of "
         f"{max(saes):g}-{min(saes):g} kg O2/kWh (SCENARIO). With aeration sized to kLa_O2 = {k_mid:g} 1/d (SCENARIO) "
         f"and field AE = {s_mid:g}, the energy per kg N removed is {fmt_d(sup, '{:.1f}')} kWh. The idle fraction of "
         "stirred time is bounded, not measured: only part of each interval carries a process state, so it lies "
         f"between the inactive share of all bins and that share plus every bin without a state: {idle_txt} (lower-upper "
         f"bound, classified fraction in brackets; the value over classified time is computed only at a coverage ≥ "
         f"{states.MIN_CLASSIFIED_FRACTION:g}; Fig. 8c; QUALIFIED from pulse 2 on). "
         + states.interval_limits_text(t["states_interval_coverage.csv"], t["states_fractions.csv"])),
        (f"9. **The closed frame is not established.** Conditional on C(t), N conservation over pulse 3 needs a net water "
         f"loss of at least {e_min:.2f} mL/d; liquid was present at the end, so the loss was below {e_max:.2f} mL/d. "
         "Vessel masses would settle it. " + volume_refit_txt(t["model_volume_scenario_summary.csv"])
         + (" " + model.volume_scaling_text(t["model_volume_scenario_summary.csv"])
            if len(t["model_volume_scenario_summary.csv"]) else "")),
    ]
    return lines


def volume_refit_txt(vs: pd.DataFrame) -> str:
    """Y7: the M6 refits per SCENARIO evaporation rate, labelled SCENARIO (they never feed a primary row)."""
    if vs.empty:
        return "M6 per volume scenario: not computed (no SCENARIO evaporation rate configured)."
    changed = vs[vs["selected_model_changes"].astype(str).str.lower() == "true"]
    prim = str(vs["primary_selected"].iloc[0])
    per = "; ".join(
        f"{r['net_loss_mL_per_d']:g} mL/d: selects {r['selected_in_scenario']}, {prim} selection skill "
        f"{float(r['primary_model_sel_skill_vs_null']):.3f} (test {float(r['primary_model_test_skill_vs_null']):.3f}, "
        f"forecast {float(r['primary_model_forecast_skill_vs_null']):.3f}), C(t) {float(r['C_start_median_mgN_L']):.0f} "
        f"mg N/L at day 1 ({r['C_start_vs_constant_volume_band']}) and {float(r['C_end_median_mgN_L']):.0f} at the end "
        f"({r['C_end_vs_constant_volume_band']})" for _, r in vs.iterrows())
    return ("M6 refitted per SCENARIO evaporation rate (every candidate under the same policy and split; SCENARIO, never "
            f"a primary row; `model_volume_scenario_summary.csv`): {per}. "
            + (f"**The selected model changes** under {', '.join(f'{x:g} mL/d' for x in changed['net_loss_mL_per_d'])}."
               if len(changed) else f"The selected model ({prim}) does not change under any scenario.")
            + " The primary C(t) and every primary row stay at constant volume; no primary volume scenario is chosen.")


def volume_status_txt(vs: pd.DataFrame) -> str:
    if vs.empty:
        return "not computed"
    changed = vs[vs["selected_model_changes"].astype(str).str.lower() == "true"]
    return (f"the selected model CHANGES under {', '.join(f'{x:g} mL/d' for x in changed['net_loss_mL_per_d'])} "
            f"(to {', '.join(changed['selected_in_scenario'])})" if len(changed) else
            f"selected model unchanged ({vs['primary_selected'].iloc[0]}) at "
            f"{', '.join(f'{x:g}' for x in vs['net_loss_mL_per_d'])} mL/d")


def other_mechanism(fit: pd.DataFrame) -> list:
    """The converged candidates of the other mechanism (alkalinity limitation against activity decay), best selection
    RMSE first."""
    conv = fit["cal_converged"].astype(str).str.lower().eq("true")
    sel = fit.loc[fit["selected"].astype(bool)].iloc[0]
    alk = fit["model"].astype(str).str.startswith("alkalinity")
    other = fit[conv & (~alk if str(sel["model"]).startswith("alkalinity") else alk)
                & (fit["model"].astype(str) != "null_no_consumption")]
    return [r for _, r in other.sort_values("sel_rmse_log", kind="stable").iterrows()]


def best_rival(fit: pd.DataFrame):
    """The best converged candidate of the other mechanism."""
    other = other_mechanism(fit)
    return other[0] if other else None


def _range(r, col: str) -> tuple[float, float]:
    return float(r[f"{col}_skill_min_inside_95"]), float(r[f"{col}_skill_max_inside_95"])


def _fmt_range(lo: float, hi: float) -> str:
    return f"{lo:.3f}" if f"{lo:.3f}" == f"{hi:.3f}" else f"{lo:.3f} to {hi:.3f}"


def rival_range_txt(rival) -> str:
    """Y1: the skills of a candidate over its converged starts inside Δ < 3.84 of its calibration, if more than one."""
    n = int(float(rival.get("n_starts_inside_95", 0) or 0))
    if n <= 1:
        return ""
    minima = str(rival.get("alk0_converged_minima_inside_95_meq_L", ""))
    alk = f" (Alk0 {minima.replace('; ', ' and ')} meq/L)" if minima and minima != "nan" and ";" in minima else ""
    return (f"; but {n} of its converged starts lie inside Δ < 3.84 of that calibration{alk}, which the calibration "
            f"window cannot tell apart, and over them it scores selection {_fmt_range(*_range(rival, 'sel'))}, test "
            f"{_fmt_range(*_range(rival, 'test'))} and forecast {_fmt_range(*_range(rival, 'forecast'))} (DIAGNOSTIC "
            "ranges; `model_calibration_starts.csv`)")


def ranking(sel, rival, col: str) -> str:
    """'better' / 'worse' (the selected model against the rival) only if the two skill ranges are disjoint."""
    s_lo, s_hi = _range(sel, col)
    r_lo, r_hi = _range(rival, col)
    return "better" if s_lo > r_hi else "worse" if s_hi < r_lo else "undetermined"


def selection_margin(fit: pd.DataFrame, sel) -> tuple[float, float]:
    """The selected model's lowest selection skill over its starts inside Δ < 3.84, and the highest over those of
    every other converged candidate (the selection holds if the first exceeds the second)."""
    conv = fit["cal_converged"].astype(str).str.lower().eq("true")
    others = fit[conv & (fit["model"] != sel["model"]) & (fit["model"] != "null_no_consumption")]
    return _range(sel, "sel")[0], others["sel_skill_max_inside_95"].astype(float).max() if len(others) else np.nan


def selection_robustness_txt(fit: pd.DataFrame, sel, rivals) -> str:
    """Y1 (b): whether the selection holds against every converged start inside Δ < 3.84 of every other candidate,
    and whether the out-of-sample ranking against each converged candidate of the other mechanism is determined."""
    s_min, o_max = selection_margin(fit, sel)
    out = (f"The selection of {sel['model']} {'holds' if s_min > o_max else 'does not hold'} against every converged "
           f"start inside Δ < 3.84 of every other candidate (its lowest selection skill {s_min:.3f} against the highest "
           f"of those starts, {o_max:.3f}).")
    if rivals is None:
        return out
    undetermined = False
    for rival in ([rivals] if isinstance(rivals, pd.Series) else rivals):
        detail = "; ".join(f"{c} {_fmt_range(*_range(sel, c))} against {_fmt_range(*_range(rival, c))}"
                           for c in ("test", "forecast"))
        if (rank := overall_rank(sel, rival)) != "undetermined":
            out += (f" Out of sample, {sel['model']} is {rank} than {rival['model']} over all those starts "
                    f"({detail}).")
        else:
            undetermined = True
            out += f" The out-of-sample ranking of {sel['model']} against {rival['model']} is undetermined ({detail})."
    return out + (" The data therefore do not discriminate activity decay from substrate-side (alkalinity) limitation."
                  if undetermined else "")


def overall_rank(sel, rival) -> str:
    """The selected model against a rival out of sample: 'better' / 'worse' only if the test and the forecast ranges
    are both disjoint in the same direction (``ranking``), else 'undetermined'."""
    ranks = {ranking(sel, rival, c) for c in ("test", "forecast")}
    return ranks.pop() if len(ranks) == 1 and "undetermined" not in ranks else "undetermined"


def mechanism_summary_txt(fit: pd.DataFrame, sel, rivals: list) -> str:
    """Z4: §1 point 7 in short: whether the selection holds, and the ranking against each converged candidate of the
    other mechanism with its test-skill range (the rules of ``selection_robustness_txt``)."""
    s_min, o_max = selection_margin(fit, sel)
    out = (f"the selection {'holds' if s_min > o_max else 'does not hold'} against every converged start inside "
           f"Δ < 3.84 of every other candidate (its lowest selection skill {s_min:.3f} against their highest, "
           f"{o_max:.3f}).")
    by = {k: [r for r in rivals if overall_rank(sel, r) == k] for k in ("better", "worse", "undetermined")}
    rng = lambda rs: " and ".join(f"{r['model']} (test {_fmt_range(*_range(r, 'test'))})" for r in rs)
    parts = [f"{sel['model']} is {k} than {rng(by[k])}" for k in ("better", "worse") if by[k]]
    if by["undetermined"]:
        parts.append(f"its ranking against {rng(by['undetermined'])} is undetermined, so the data do not discriminate "
                     "activity decay from substrate-side (alkalinity) limitation")
    return out + (f" Out of sample, {'; '.join(parts)}." if parts else "")


def alk0_sentence(fit: pd.DataFrame, alk_cmp: pd.DataFrame, required: float, checkpoint: str) -> str:
    """Z4: the Alk0 evidence of the converged alkalinity candidates in one sentence (from ``alk0_evidence_txt`` and
    ``alk0_comparison``)."""
    if alk_cmp.empty:
        return ""
    ok = alk_cmp["cal_converged"].astype(str).str.lower().eq("true")
    ident = int(alk_cmp.loc[ok, "alk0_profile_status"].astype(str).eq("identifiable").sum())
    verdicts = {str(a).split(":")[0] for a in alk_cmp.loc[ok, "agreement"]}
    verdict = (verdicts.pop() if len(verdicts) == 1 else
               "; ".join(f"{m}: {str(a).split(':')[0]}" for m, a in zip(alk_cmp.loc[ok, "model"], alk_cmp.loc[ok, "agreement"])))
    return (f"{alk0_evidence_txt(fit)}, identifiable in {f'{ident} of {int(ok.sum())}' if ident else 'none'}; against the {required:.1f} meq/L "
            f"that the C(t) removal needs if nitrified ({checkpoint}, S2): {verdict}.")


def mechanism_table(fit: pd.DataFrame) -> pd.DataFrame:
    """Z4: per candidate, the starts inside Δ < 3.84, the selection, test and forecast ranges over them, the ranking
    of the selected model against it (test; forecast; ``ranking``), the Alk0 minima and the bound flag."""
    sel = fit.loc[fit["selected"].astype(bool)].iloc[0]
    rows = []
    for _, r in fit.iterrows():
        conv = str(r["cal_converged"]).lower() == "true"
        ok = conv and pd.notna(r["sel_skill_min_inside_95"])
        minima = str(r.get("alk0_converged_minima_inside_95_meq_L", ""))
        rows.append({
            "candidate": r["model"], "converged": r["cal_converged"],
            "starts inside Δ < 3.84": int(float(r["n_starts_inside_95"])) if pd.notna(r["n_starts_inside_95"]) else 0,
            **{f"{c} skill": _fmt_range(*_range(r, k)) if ok else "not computed"
               for c, k in (("selection", "sel"), ("test", "test"), ("forecast", "forecast"))},
            f"{sel['model']} against it (test; forecast)": (
                "selected" if r["model"] == sel["model"] else "reference" if r["model"] == "null_no_consumption"
                else f"{ranking(sel, r, 'test')}; {ranking(sel, r, 'forecast')}" if ok else "not computed"),
            "Alk0 minima inside Δ < 3.84 (meq/L)": minima if minima not in ("", "nan") else "",
            "cal_at_bounds": "" if pd.isna(r.get("cal_at_bounds")) else str(r.get("cal_at_bounds"))})
    return pd.DataFrame(rows)


def alk0_evidence_txt(fit: pd.DataFrame) -> str:
    """Y2: the Alk0 range of the converged alkalinity candidates only (every distinct converged minimum inside
    Δ < 3.84); the unconverged candidates are named separately, without a value."""
    alk = fit[fit["model"].astype(str).str.startswith("alkalinity")]
    ok = alk["cal_converged"].astype(str).str.lower().eq("true")
    vals = []
    for _, r in alk[ok].iterrows():
        minima = str(r.get("alk0_converged_minima_inside_95_meq_L", ""))
        vals += [float(v) for v in minima.split("; ")] if minima and minima != "nan" else [float(r["alk0_fitted_meq_L"])]
    txt = (f"Alk0 of the converged alkalinity candidates {min(vals):.1f}-{max(vals):.1f} meq/L (every converged minimum "
           "inside Δ < 3.84)" if vals else "no alkalinity candidate converged")
    if (~ok).any():
        txt += f"; not converged (no estimate): {', '.join(alk.loc[~ok, 'model'])}"
    return txt


def numeric_cells(df: pd.DataFrame) -> pd.DataFrame:
    """Cells that read back from CSV as text but hold a number (columns mixing numbers and statuses) become floats,
    so that ``md_table`` rounds them (X7)."""
    def conv(v):
        if isinstance(v, str):
            try:
                return float(v)
            except ValueError:
                return v
        return v
    return df.apply(lambda col: col if pd.api.types.is_numeric_dtype(col) or pd.api.types.is_bool_dtype(col)
                    else col.astype(object).map(conv))


def alk0_comparison(fit: pd.DataFrame, cons: pd.DataFrame) -> pd.DataFrame:
    """Fitted Alk0 of each alkalinity candidate against the S2 requirement at the end of its calibration window
    (pre_pulse3 for the split from the configured pulses) and at the end of the record."""
    alk = fit[fit["model"].astype(str).str.startswith("alkalinity")]
    c = cons.set_index("checkpoint")["alk0_required_if_nitrified_meq_L"].astype(float)
    cal_end = [k for k in c.index if k.startswith("pre_pulse")][-2] if len(c) > 2 else c.index[0]
    rows = []
    for _, r in alk.iterrows():
        a0 = float(r["alk0_fitted_meq_L"])
        ok = str(r["cal_converged"]).lower() == "true"
        req = float(c[cal_end])
        low = float(r["alk0_ci95_low_meq_L"]) if pd.notna(r["alk0_ci95_low_meq_L"]) else np.nan
        rows.append({"model": r["model"], "cal_converged": r["cal_converged"], "alk0_fitted_meq_L": a0,
                     "alk0_ci95_low_meq_L": r["alk0_ci95_low_meq_L"], "alk0_ci95_high_meq_L": r["alk0_ci95_high_meq_L"],
                     "alk0_profile_status": r["alk0_profile_status"],
                     "alk0_converged_minima_inside_95_meq_L": r.get("alk0_converged_minima_inside_95_meq_L", np.nan),
                     f"S2 requirement at {cal_end} (meq/L)": req, "S2 requirement at end (meq/L)": float(c["end"]),
                     "agreement": ("not comparable: calibration did not converge (the value is not an estimate)" if not ok
                                   else f"agrees: the 95 % lower profile limit {low:.1f} meq/L is above the requirement"
                                   if np.isfinite(low) and low >= req else
                                   "compatible, not a test: the fitted value is above the requirement, but the 95 % "
                                   "profile reaches below it or is open" if a0 >= req else
                                   "disagrees: the fitted Alk0 is below the requirement")})
    return pd.DataFrame(rows)


def assemble_report(cfg, out_dir) -> Path:
    from . import alkalinity, figures, model, states
    from .dynamics import onsager_conclusion

    out_dir = Path(out_dir)
    tables = out_dir / "tables"
    names = [n for n in TABLE_DESCRIPTIONS if n != "all_references.csv"]
    t = {n: _read(tables, n) for n in names}
    t["dynamics_windows.csv"] = _read(tables, "dynamics_windows.csv")
    check_required(t)  # raises before anything is written
    refs = all_references(tables)
    figures.make_all(cfg, t, out_dir)
    for fig, _ in FIGURES:
        if not (out_dir / "figures" / fig).is_file():
            raise FileNotFoundError(f"{out_dir / 'figures' / fig} missing")
    refs.to_csv(tables / "all_references.csv", index=False)
    labels = status_labels(t)
    pv = t["sensor_validation.csv"]
    verdicts = t["sensor_verdicts.csv"][["sensor", "verdict", "criterion"]]
    vessel, pulses, sensors = cfg.raw["vessel"], cfg.raw["pulses"], cfg.raw["sensors"]
    index = pd.DataFrame([(n, TABLE_DESCRIPTIONS[n], len(t[n]) if n in t else len(refs)) for n in TABLE_DESCRIPTIONS],
                         columns=["table (outputs/tables/)", "content", "rows"])
    dp = t["data_properties.csv"][["sensor", "file", "rows_valid", "rows_dropped", "t_min_d", "t_max_d", "value_min",
                                   "value_max", "interval_median_s"]]
    pm = t["pulses.csv"]
    fig = {f: f"![{c}](../figures/{f})\n\n*{c}*\n" for f, c in FIGURES}
    sc = cfg.raw["scenarios"]
    fit = t["model_fit.csv"]
    lines = [
        "# Non-equilibrium analysis of an NH4Cl-pulsed nitrifying batch",
        "",
        ("Generated by `python -m bnr_neq run-all` (report stage) from the raw files in `data/`. Do not edit by hand. "
         "Implementer's report, not an approval. Labels: SCENARIO = assumed value with its range; MODEL-BASED = from "
         "the selected M6 rate law; SENSITIVITY = an alternative NH4 observation model kept for comparison; QUALIFIED = "
         "derived from a C(t) that contradicts the W2 alkalinity bound."),
        "",
        "## 1. What the data say about the process",
        "",
        *[x + "\n" for x in summary_lines(cfg, t)],
        fig["fig1_linkage.png"],
        "## 2. Energy and efficiency",
        "",
        (f"The metric family per pulse (constant volume; all rows with inputs and status in `pulses.csv`). SCENARIO "
         f"ranges: kLa_O2 = {', '.join(f'{x:g}' for x in sc['kLa_O2_per_d'])} 1/d (not measured); field (process-water) "
         f"aeration efficiency AE = {', '.join(f'{x:g}' for x in sc['sae_field_kgO2_per_kWh'])} kg O2/kWh (field "
         "conditions, not the clean-water standard SAE; fine pore: 5.0-6.5 standard against 2.5-3.5 field; metric "
         "suffix `AE`, configuration key `sae_field_kgO2_per_kWh`), source: "
         f"{sc['sae_source']}. Low and high are the extremes of the diffused-air rows that include blower power, the "
         "central value their geometric mean (pass 6, G5, `docs/progress.md`); they are values of full-scale "
         "equipment, not of this vessel, so these are scenario values only. O2 supply = kLa (C* − DO) V Δt with DO = 0 while DO is not measured, so utilization is a "
         "lower bound and the supply-based energy an upper bound; the supply-based energy also divides a capacity at "
         "DO = 0 by an AE at the plant's operating DO (> 0), an inconsistency that inflates it further. The demand-based energy "
         "(4.57/AE kWh/kg N) is the "
         "floor that any aeration system must spend on the N actually removed; the supply-based value adds the aeration "
         "spent while the biomass was idle. The stirrer power is a laboratory-specific term (not measured, "
         "`stirrer_power_W` missing)."),
        "",
        md_table(energy_table(pm), digits=3),
        "",
        fig["fig8_efficiency.png"],
        "## 3. Process state and the constraints on conversion",
        "",
        "### Process states (W3)",
        "",
        md_table(t["states_rules.csv"][["state", "condition", "features"]], digits=3),
        "",
        ("Idle fraction per interval as bounds with its coverage (`classified_fraction`); the CONDITIONAL value over "
         f"classified time only where the coverage is ≥ {states.MIN_CLASSIFIED_FRACTION:g} (`status`). "
         + states.interval_limits_text(t["states_interval_coverage.csv"], t["states_fractions.csv"])),
        "",
        md_table(t["states_fractions.csv"][["interval", "n_bins", "fraction_active_nitrification",
                                            "fraction_inactive_oversupplied", "fraction_oxygen_limited",
                                            "fraction_sensor_artifact", "fraction_unclassified", "classified_fraction",
                                            "idle_fraction_lower_bound", "idle_fraction_upper_bound",
                                            "idle_fraction_of_classified_time", "status"]], digits=2),
        "",
        (f"Rules SHA-256 `{t['states_rules.csv']['registration_sha256'].iloc[0]}` ({states.REGISTRATION_DOC}): "
         f"{t['states_rules.csv']['registration_status'].iloc[0]}. Each comparison interval is read against a "
         "synthetic null with the record's own group labels (`interval_reading`, `states_bootstrap_null_coverage.csv`): "
         f"an interval that covers the true difference 0 in fewer than {states.NULL_COVERAGE_MIN_OK:g} of the null "
         "replicates at some φ is too narrow, and an interval excluding 0 is then not evidence of a difference."),
        "",
        md_table(t["states_comparison.csv"][["comparison", "role", "diff", "ci95_low", "ci95_high", "n_a", "n_b",
                                              "block_bins", "status", "null_coverage_min", "interval_reading"]],
                 digits=3),
        "",
        md_table(t["states_bootstrap_null_coverage.csv"][["comparison", "phi", "n_computed", "coverage_of_zero",
                                                           "block_bins_median"]], digits=3),
        "",
        md_table(t["states_interval_coverage.csv"][["phi", "n_points", "n_windows", "nominal", "coverage_unscaled",
                                                     "coverage_scaled_estimated_rho", "coverage_scaled_true_rho",
                                                     "mean_rho_hat", "median_rho_hat"]], digits=3),
        "",
        "### Alkalinity bound (W2) and the inverted C(t) (S2)",
        "",
        (f"Air-equilibrium CO2 (conservative), constant volume, SCENARIO 25 °C. Alk0 source: "
         f"{alkalinity.alk0_source(cfg)}. Temperature, CO2-transfer and volume scenarios: `alkalinity_bound.csv`, "
         "`alkalinity.md`."),
        "",
        md_table(alkalinity.headline(t["alkalinity_bound.csv"])[["alk0_meq_L", "consistent_with_pH_record",
                                                                  "n_nit_max_end_mgN_L", "n_nit_max_fraction_of_dose",
                                                                  "n_removed_max_any_proton_sink_end_mgN_L",
                                                                  "tan_residual_min_end_mgN_L",
                                                                  "fa_residual_min_end_mgNH3_L"]], digits=3),
        "",
        provenance_line(t["nh4_truth_w2_consistency.csv"], "nh4_truth_id"),
        "",
        md_table(t["nh4_truth_w2_consistency.csv"].drop(columns=["alk0_scan_max_meq_L", *PROVENANCE_COLUMNS],
                                                        errors="ignore"), digits=3),
        "",
        fig["fig4_constraints.png"],
        "### Water volume (W8)",
        "",
        md_table(t["volume_constraints.csv"], digits=3),
        "",
        "## 4. Reduced model and mechanism (M6, W5)",
        "",
        (f"Convergence policy (one for every candidate): {model.CONVERGENCE_POLICY} A candidate whose calibration did "
         "not converge is excluded from selection (`selection_status`); its skills stay in the table and are not "
         "comparable."),
        "",
        md_table(numeric_cells(fit[["model", "n_params_calibration", "cal_converged", "cal_status", "cal_nfev",
                                    "selection_status", "sel_skill_vs_null", "test_skill_vs_null",
                                    "forecast_skill_vs_null", "selected"]]), digits=3),
        "",
        model.bounds_text(fit),
        "",
        ("Skill range over each candidate's converged starts inside Δ < 3.84 of its calibration (Δ = n_eff ln(RSS/RSS of "
         "the calibration) on the calibration window, n_eff from the calibration's residuals): these starts cannot be "
         "told apart by the calibration window. The ranges are DIAGNOSTIC (the calibration stays the lowest-RSS "
         "converged start; no start's skill is used to choose); every start is in `model_calibration_starts.csv`. "
         + selection_robustness_txt(fit, fit.loc[fit["selected"].astype(bool)].iloc[0], other_mechanism(fit))
         + " Per candidate (the ranking of the selected model by the same rule: better or worse only where the two "
         "ranges are disjoint, else undetermined; Alk0 minima: every distinct converged minimum inside Δ < 3.84):"),
        "",
        md_table(mechanism_table(fit), digits=3),
        "",
        ("Fitted initial alkalinity of the alkalinity candidates (calibration window, conditional on each model) against "
         "the Alk0 that the removal implied by C(t) needs (S2, `nh4_truth_w2_consistency.csv`):"),
        "",
        md_table(alk0_comparison(fit, t["nh4_truth_w2_consistency.csv"]), digits=3),
        "",
        ("M6 per volume SCENARIO (Y7; SCENARIO, never a primary row): every candidate calibrated under the same policy "
         "and split with V(t) = V0 − E t (doses, C, the biomass activity and the alkalinity state concentrated by "
         "V0/V(t)), and the primary model refitted on the whole record; its C(t) against the constant-volume band "
         "(medians over the first and last 0.1 d). " + volume_refit_txt(t["model_volume_scenario_summary.csv"])
         + (" " + model.volume_scaling_text(t["model_volume_scenario_summary.csv"])
            if len(t["model_volume_scenario_summary.csv"]) else "")),
        "",
        provenance_line(t["model_volume_scenario_summary.csv"]),
        "",
        md_table(numeric_cells(t["model_volume_scenario_summary.csv"].drop(
            columns=["role", "constant_volume_nh4_series", *PROVENANCE_COLUMNS], errors="ignore")), digits=3),
        "",
        md_table(numeric_cells(t["model_volume_scenarios.csv"][["volume_scenario", "model", "cal_converged",
                                                                "sel_skill_vs_null", "test_skill_vs_null",
                                                                "forecast_skill_vs_null", "selected_in_scenario",
                                                                "selection_status"]]), digits=3),
        "",
        "Consistency checks of the selected model (true by construction where noted; not findings):",
        "",
        provenance_line(t["model_checks.csv"]),
        "",
        md_table(numeric_cells(t["model_checks.csv"].drop(columns=list(PROVENANCE_COLUMNS), errors="ignore")), digits=3),
        "",
        md_table(numeric_cells(t["model_identifiability.csv"][["parameter", "estimate", "ci95_low", "ci95_high",
                                                                "identifiable"]]), digits=3),
        "",
        fig["fig6_model.png"],
        "## 5. Sensor validity",
        "",
        "From `sensor_verdicts.csv` (criteria fixed before computing the observations):",
        "",
        md_table(verdicts),
        "",
        md_table(pv[["pulse", "nh4_obs_change_mgNH4_L", "nh4_response_fraction", "ec_obs_change_mS_cm",
                     "orp_obs_change_mV", "ph_obs_change"]], digits=3),
        "",
        fig["fig3_sensor_validation.png"],
        fig["fig2_overview.png"],
        "## 6. Status of the conclusions",
        "",
        md_table(labels),
        "",
        "## 7. Experimental facts and stated assumptions",
        "",
        (f"- Batch of {vessel['volume_L'] * 1e3:g} mL, open to air: {vessel['open_to_air']}, {vessel['stirring']} stirring, "
         "no forced aeration."),
        (f"- Additions: {pulses['mass_mg']:g} mg {pulses['salt']} at days {', '.join(f'{x:g}' for x in pulses['times_d'])} "
         "(logger clock); nothing else added or removed."),
        (f"- NH4 electrode in {sensors['nh4_unit']}, factor {sensors['nh4_to_N_factor']} to mg N/L; ORP "
         f"{sensors['orp_working']} vs {sensors['orp_reference']}."),
        f"- Missing measurements: {', '.join(cfg.missing_keys)}.",
        "",
        *[f"{i}. {a}" for i, a in enumerate(assumptions(cfg, t), start=1)],
        "",
        "## Supplement S1. Stochastic-thermodynamic and information-theoretic estimators (M5)",
        "",
        ("Demoted from the main analysis (R6): the pH record is highly tied, few windows are stationary and the NH4 "
         "sub-stamp times are reconstructed, so these indices do not connect to rate, energy or efficiency."),
        "",
        md_table(t["dynamics.csv"][["series", "role", "segment", "n_windows", "n_stationary", "tie_fraction_median",
                                    "ordinal_kl_windows", "ordinal_kl_significant_windows",
                                    "asym_norm_lag1_significant_windows", "hvg_kl_significant_windows"]], digits=3),
        "",
        md_table(t["dynamics_onsager.csv"][["pulse", "pulse_response", "fluctuation_sd", "snr_pulse_response",
                                            "late_window_stationary", "verdict"]], digits=3),
        "",
        onsager_conclusion(t["dynamics_onsager.csv"]),
        "",
        fig["fig5_dynamics.png"],
        "## Supplement S2. Detailed tables",
        "",
        "### Verified data properties",
        "",
        md_table(dp, digits=6),
        "",
        "### Consumption slopes on the electrode reading (M3; gain-corrected and reading-as-TAN columns are SENSITIVITY)",
        "",
        md_table(t["segments.csv"][["sub_segment", "fit_start_d", "fit_end_d", "start_boundary", "slope_mgNH4_L_d",
                                    "held_fraction", "consumption_SENSITIVITY_linear_gain_mgN_L_d", "status"]], digits=3),
        "",
        "### Slope windows of the inverted C(t) (W1)",
        "",
        md_table(t["nh4_truth_windows.csv"][["sub_segment", "governing_pulse", "fit_start_d", "fit_end_d",
                                             "median_C_inv_mgN_L", "consumption_mgN_L_d", "consumption_ci95_low_mgN_L_d",
                                             "consumption_ci95_high_mgN_L_d", "fa_median_mgNH3_L_SCENARIO_T25C",
                                             "o2_demand_mgO2_L_d", "status"]], digits=3),
        "",
        "### Thermodynamics and balances (M4)",
        "",
        md_table(t["thermo_affinity_steps.csv"].drop(columns=["note"]), digits=3),
        "",
        md_table(t["balance_nitrogen.csv"], digits=4),
        "",
        md_table(t["balance_alkalinity.csv"], digits=4),
        "",
        "### Pulse-level descriptors (M7)",
        "",
        fig["fig7_pulses.png"],
        md_table(t["pulses_associations.csv"], digits=3),
        "",
        "### Reference values",
        "",
        (f"{int(refs['pass'].sum())} of {len(refs)} reproduced within tolerance; the others are explained in "
         "`docs/progress.md` and in each milestone report. They are regression checks against the preliminary analysis, "
         "not external validation."),
        "",
        md_table(refs, digits=4),
        "",
        "### Every table",
        "",
        md_table(index),
        "",
        ("Milestone reports in `outputs/reports/`: `ingest.md`, `sensor_validation.md`, `segments.md`, `thermo.md`, "
         "`alkalinity.md`, `dynamics.md`, `model.md`, `states.md`, `pulses.md`."),
    ]
    path = out_dir / "reports" / "analysis_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
