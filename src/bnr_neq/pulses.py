"""Phase 5: pulse-level efficiency and energy descriptors (M7).

The output is a tidy table: one row per (pulse, metric) with its value, unit,
basis, the inputs it was computed from, and a status that is either
``computed`` or an explicit ``not computed: <reason>``.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import alkalinity, energy, model, thermo, volume
from . import segments as S
from .config import Config, MissingInputError
from .constants import (
    DG0_NITRIFICATION_KJ_MOL_N,
    DO_SAT_25C_MG_L,
    KELVIN_OFFSET,
    M_N_G_MOL,
    SECONDS_PER_DAY,
)
from .states import MIN_CLASSIFIED_FRACTION

END_D = 47.92
EXPLORATORY = "EXPLORATORY (n = 4 pulses)"
# X2: idle-fraction metrics carry their nature in the name (states_fractions.csv columns in the comments)
IDLE_LOWER = "idle_fraction_LOWER_BOUND"  # idle_fraction_lower_bound
IDLE_UPPER = "idle_fraction_UPPER_BOUND"  # idle_fraction_upper_bound
IDLE_CONDITIONAL = "idle_fraction_of_classified_time_CONDITIONAL"  # idle_fraction_of_classified_time

INPUTS = {
    "segments": "outputs/tables/segments.csv",
    "steps": "outputs/tables/thermo_affinity_steps.csv",
    "oxygen": "outputs/tables/thermo_oxygen.csv",
    "pulses": "outputs/tables/sensor_validation.csv",
    "model": "outputs/processed/model_state.parquet",
    "model_fit": "outputs/tables/model_fit.csv",
    "grid": "outputs/processed/grid_0p01d.parquet",
    "truth": "outputs/processed/nh4_truth.parquet",
    "truth_vectors": "outputs/tables/nh4_truth_vectors.csv",
    "truth_summary": "outputs/tables/nh4_truth_summary.csv",
    "states_fractions": "outputs/tables/states_fractions.csv",
}
OUTPUT_ONLY = {"truth_windows": "outputs/tables/nh4_truth_windows.csv"}
INPUTS.update(OUTPUT_ONLY)


def load_inputs(out_dir: Path) -> dict:
    out = {}
    for k, rel in INPUTS.items():
        if k in OUTPUT_ONLY:
            continue
        p = Path(out_dir) / rel.removeprefix("outputs/")
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run the earlier stages (`python -m bnr_neq run-all`) first")
        out[k] = pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)
    return out


def governed_windows(segments: pd.DataFrame, k: int) -> pd.DataFrame:
    """Computed M3 slope windows whose governing pulse is k."""
    s = segments[(segments["status"] == "computed") & (segments["governing_pulse"] == k)]
    return s.assign(duration_d=s["fit_end_d"] - s["fit_start_d"])


def positive_windows(w: pd.DataFrame) -> pd.DataFrame:
    """Windows with a positive gain-corrected consumption (M3 status ``computed``); the others are
    ``not computed: non-positive ...`` in segments.csv and are excluded from duration-weighted rates."""
    return w[w["gain_correction_status"] == "computed"]


def weighted(values: pd.Series, weights: pd.Series) -> float:
    return float(np.sum(values * weights) / np.sum(weights))


def weighted_rate(w: pd.DataFrame, col: str, lo: str, hi: str, sign: float = 1.0) -> tuple[float, float, float, float]:
    """Duration-weighted value, interval bounds and held fraction over the windows ``w``.

    The bounds are the duration-weighted M3 interval bounds (window errors treated as fully correlated,
    which is conservative); ``sign = -1`` turns a slope into a consumption and swaps the bounds.
    """
    a, b = (lo, hi) if sign > 0 else (hi, lo)
    d = w["duration_d"]
    return (sign * weighted(w[col], d), sign * weighted(w[a], d), sign * weighted(w[b], d),
            weighted(w["held_fraction"], d))


def model_integral(state: pd.DataFrame, t0: float, t1: float) -> tuple[float, float]:
    """∫ r dt (mg N/L) over [t0, t1) from the 0.01 d model state, and the mean rate (mg N/L/d)."""
    m = (state["t_d"] >= t0) & (state["t_d"] < t1)
    dt = float(np.median(np.diff(state["t_d"])))
    integral = float(state.loc[m, "model_rate_mgN_L_d"].sum() * dt)
    return integral, integral / (t1 - t0)


def free_energy_per_pulse_J(cfg: Config) -> float:
    """|ΔG°′| x moles N added by one pulse in the vessel (standard conditions, pH 7)."""
    dose = thermo.pulse_dose_from_config(cfg)
    mol_n = dose.c_mmol_L / 1e3 * cfg.get("vessel.volume_L")
    return abs(DG0_NITRIFICATION_KJ_MOL_N) * 1e3 * mol_n


def stirrer_energy_J(power_W: float, t0_d: float, t1_d: float) -> float:
    """Electrical energy of the stirrer over [t0, t1], J."""
    if power_W < 0 or t1_d < t0_d:
        raise ValueError("need power >= 0 and t1 >= t0")
    return power_W * (t1_d - t0_d) * SECONDS_PER_DAY


PRE_WINDOW_D = 0.1  # C(t) "just before" a pulse or the end: median over the preceding 0.1 d (M2 pre window)


def band_series(truth: pd.DataFrame, vectors: pd.DataFrame) -> list[np.ndarray]:
    """Inversions of the reading under every profile vector of ``nh4_truth_vectors.csv`` (vector 0 = best)."""
    t, reading = truth["t_d"].to_numpy(), truth["electrode_obs_mgN_L"].to_numpy()
    return [model.invert_from_row(r, t, reading) for _, r in vectors.iterrows()]


def pre_median(t: np.ndarray, c: np.ndarray, t_event: float) -> float:
    m = (t >= t_event - PRE_WINDOW_D - 1e-9) & (t < t_event)
    return float(np.median(c[m]))


def truth_windows(cfg: Config, inputs: dict) -> pd.DataFrame:
    """Per-window consumption, free ammonia and O2 demand of the M6-inverted C(t) (``nh4_truth_windows.csv``)."""
    tw = S.truth_slope_table(cfg, inputs["truth"], inputs["segments"],
                             band_series(inputs["truth"], inputs["truth_vectors"]),
                             volume.volume_scenarios(cfg, END_D))
    if tw["nh4_truth_id"].nunique() != 1 or tw["nh4_observation_model"].nunique() != 1:
        raise model.ObservationModelError("nh4_truth_windows mixes NH4 truth series")
    return tw


def eta_band_on_unit_interval(band: tuple[float, float]) -> tuple[tuple[float, float], str]:
    """X5: η is a fraction of one pulse, so its band is clipped to [0, 1]; the unclipped band stays on the
    removal_to_dose_ratio row. Returns the band and a note for the basis (empty when nothing was clipped)."""
    lo, hi = band
    if not (np.isfinite(lo) and np.isfinite(hi)) or (0.0 <= lo and hi <= 1.0):
        return band, ""
    return ((min(max(lo, 0.0), 1.0), min(max(hi, 0.0), 1.0)),
            f" Param band clipped to [0, 1] (unclipped {lo:.3f} to {hi:.3f} on removal_to_dose_ratio).")


def open_band_note(summary: pd.DataFrame | None) -> str:
    if summary is None:
        return ""
    v = summary.set_index("key")["value"].astype(str)
    sides = [s for s in ("low", "high") if v.get(f"band_{s}_open", "False") == "True"]
    return (" Param band: envelope of the 95 % profile vectors" + (f", open on the {' and '.join(sides)} side(s) "
            "(not bounded by the data there)" if sides else ""))


def pulse_metrics(cfg: Config, inputs: dict) -> pd.DataFrame:
    seg, steps, grid, state = inputs["segments"], inputs["steps"], inputs["grid"], inputs["model"]
    truth = inputs["truth"]
    tw = inputs["truth_windows"] if "truth_windows" in inputs else truth_windows(cfg, inputs)
    bands = band_series(truth, inputs["truth_vectors"])
    band_txt = open_band_note(inputs.get("truth_summary"))
    selected = inputs["model_fit"].loc[inputs["model_fit"]["selected"], "model"].iloc[0]
    dose = thermo.pulse_dose_from_config(cfg)
    nf = cfg.get("sensors.nh4_to_N_factor")
    g0 = free_energy_per_pulse_J(cfg)
    times = [*cfg.pulse_times_d, END_D]
    vols = volume.volume_scenarios(cfg, END_D)
    n_mg = volume.pulse_n_mg(cfg)
    t_tr, c_tr = truth["t_d"].to_numpy(), truth["C_inv_mgN_L"].to_numpy()
    tid = str(truth["nh4_truth_id"].iloc[0])
    OBS, RAW, GAIN, IDENT = model.OBSERVATION_MODEL, model.RAW_READING, model.LINEAR_GAIN, model.IDENTITY
    INDEP = model.NH4_INDEPENDENT
    state_series = f"M6_state:{selected}"
    c_start = float(np.median(c_tr[t_tr < t_tr[0] + PRE_WINDOW_D]))
    rows = []

    def add(k, metric, value, unit, basis, inputs_, status="computed", ci=(np.nan, np.nan), held=np.nan,
            obs="", series="", band=(np.nan, np.nan)):
        rows.append({"pulse": k, "metric": metric, "value": value, "ci95_low": ci[0], "ci95_high": ci[1],
                     "param_band_low": band[0], "param_band_high": band[1], "held_fraction": held, "unit": unit,
                     "basis": basis, "inputs": inputs_, "status": status, "nh4_observation_model": obs,
                     "nh4_series": series})

    def nc(k, metric, unit, basis, inputs_, status, obs, series, held=np.nan):
        add(k, metric, np.nan, unit, basis, inputs_, status, held=held, obs=obs, series=series)

    try:
        stirrer_w = float(cfg.measurement("stirrer_power_W"))
        stirrer_status = "computed"
    except MissingInputError as exc:
        stirrer_w, stirrer_status = None, str(exc)

    # --- W4: oxygen and energy metric family (closure over the loop variables above)
    klas = [float(x) for x in cfg.get("scenarios.kLa_O2_per_d")]
    saes = [float(x) for x in cfg.get("scenarios.sae_field_kgO2_per_kWh")]
    sae_src = str(cfg.get("scenarios.sae_source"))
    v0 = float(cfg.get("vessel.volume_L"))
    try:
        do_mg_L, do_txt = float(cfg.measurement("DO_mg_L")), "measured DO"
    except MissingInputError:
        do_mg_L, do_txt = 0.0, "DO = 0, i.e. the transfer capacity (DO_mg_L missing; an upper bound on the supply)"
    fractions = inputs.get("states_fractions")

    def energy_family(k, tk, tn, n_removed):
        removed_ok = np.isfinite(n_removed)
        no_removal = "not computed: negative dose-anchored removal on the M6-inverted C(t)"
        if removed_ok:
            o2 = float(energy.oxygen_required_mgO2(n_removed))
            add(k, "o2_required_mgO2", o2, "mg O2", "4.57 g O2/g N x N removed (nitrification)",
                "n_removed_mgN; constants.O2_PER_N_G_G", obs=OBS, series=tid)
            add(k, "o2_required_nitritation_case_mgO2", float(energy.oxygen_required_mgO2(n_removed, nitritation=True)),
                "mg O2", "STOICHIOMETRY CASE: nitritation only, 3.43 g O2/g N x N removed",
                "n_removed_mgN; constants.O2_PER_N_NITRITATION_G_G", obs=OBS, series=tid)
        else:
            for metric in ("o2_required_mgO2", "o2_required_nitritation_case_mgO2"):
                nc(k, metric, "mg O2", "", "n_removed_mgN", no_removal, OBS, tid)
        g_rel = float(energy.gibbs_released_J(n_removed)) if removed_ok else np.nan
        for sae in saes:
            sl = f"AE{sae:g}"
            sae_basis = (f"SCENARIO field (process-water) aeration efficiency AE {sae:g} kg O2/kWh at about 2 mg/L DO, "
                         f"not the clean-water SAE ({sae_src})")
            if removed_ok and n_removed > 0:
                e = float(energy.aeration_energy_kWh(o2, sae))
                add(k, f"aeration_energy_demand_Wh_SCENARIO_{sl}", 1e3 * e, "Wh",
                    f"o2_required / AE (vessel scale); {sae_basis}", "o2_required_mgO2; config [scenarios]",
                    obs=OBS, series=tid)
                add(k, f"aeration_energy_demand_kWh_per_kgN_SCENARIO_{sl}", float(energy.specific_energy_kWh_per_kgN(e, n_removed)),
                    "kWh/kg N", f"floor: 4.57 g O2/g N / AE; {sae_basis}", "aeration_energy_demand; n_removed_mgN",
                    obs=OBS, series=tid)
                add(k, f"gibbs_to_aeration_energy_demand_SCENARIO_{sl}_T25C", g_rel / float(energy.kwh_to_J(e)), "-",
                    f"free energy released / demand-based aeration energy; {sae_basis}; ΔG°′ SCENARIO 25 °C",
                    "free_energy_released_J_SCENARIO_T25C; aeration_energy_demand", obs=OBS, series=tid)
            else:
                for metric, unit in ((f"aeration_energy_demand_Wh_SCENARIO_{sl}", "Wh"),
                                     (f"aeration_energy_demand_kWh_per_kgN_SCENARIO_{sl}", "kWh/kg N"),
                                     (f"gibbs_to_aeration_energy_demand_SCENARIO_{sl}_T25C", "-")):
                    nc(k, metric, unit, sae_basis, "n_removed_mgN", no_removal if not removed_ok else
                       "not computed: no N removed", OBS, tid)
        for kla in klas:
            kl = f"kLa{kla:g}"
            cap = float(energy.o2_supply_capacity_mgO2(kla, DO_SAT_25C_MG_L, do_mg_L, v0, tn - tk))
            kla_basis = (f"SCENARIO kLa_O2 = {kla:g} 1/d (not measured); C* = {DO_SAT_25C_MG_L} mg/L (SCENARIO 25 °C); "
                         f"{do_txt}")
            add(k, f"o2_supply_capacity_mgO2_SCENARIO_{kl}_T25C", cap, "mg O2",
                f"kLa (C* − C) V Δt over {tk:g}-{tn:g} d; {kla_basis}", "config [scenarios] kLa_O2_per_d, [vessel]",
                obs=INDEP)
            infeasible = removed_ok and o2 > cap
            infeasible_msg = (f"not computed: SCENARIO infeasible: O2 required exceeds the transfer capacity at kLa "
                              f"{kla:g} 1/d (required/capacity = {o2 / cap:.2f})" if infeasible else "")
            if infeasible:
                nc(k, f"o2_utilization_SCENARIO_{kl}_T25C", "-", kla_basis, "o2_required_mgO2; o2_supply_capacity",
                   infeasible_msg, OBS, tid)
            elif removed_ok:
                add(k, f"o2_utilization_SCENARIO_{kl}_T25C", float(energy.utilization(o2, cap)), "-",
                    f"o2_required / o2_supply_capacity (a lower bound while DO is not measured); {kla_basis}",
                    "o2_required_mgO2; o2_supply_capacity", obs=OBS, series=tid)
            else:
                nc(k, f"o2_utilization_SCENARIO_{kl}_T25C", "-", kla_basis, "o2_required_mgO2", no_removal, OBS, tid)
            for sae in saes:
                sl = f"AE{sae:g}"
                metric_e = f"aeration_energy_supply_kWh_per_kgN_SCENARIO_{sl}_{kl}_T25C"
                metric_g = f"gibbs_to_aeration_energy_supply_SCENARIO_{sl}_{kl}_T25C"
                if infeasible:
                    for metric, unit in ((metric_e, "kWh/kg N"), (metric_g, "-")):
                        nc(k, metric, unit, kla_basis, "o2_supply_capacity; n_removed_mgN", infeasible_msg, OBS, tid)
                elif removed_ok and n_removed > 0:
                    e_sup = float(energy.aeration_energy_kWh(cap, sae))
                    add(k, metric_e, float(energy.specific_energy_kWh_per_kgN(e_sup, n_removed)), "kWh/kg N",
                        f"supply capacity / AE / N removed: includes aeration while idle; {kla_basis}; SCENARIO field AE "
                        f"{sae:g} kg O2/kWh ({sae_src}). Inconsistent by construction: the capacity is at DO = 0 while "
                        "the AE is defined at DO = 2 mg/L, so this is an upper bound", "o2_supply_capacity; n_removed_mgN",
                        obs=OBS, series=tid)
                    add(k, metric_g, g_rel / float(energy.kwh_to_J(e_sup)), "-",
                        f"free energy released / supply-based aeration energy; {kla_basis}; SCENARIO field AE {sae:g}",
                        "free_energy_released_J_SCENARIO_T25C; o2_supply_capacity", obs=OBS, series=tid)
                else:
                    for metric, unit in ((metric_e, "kWh/kg N"), (metric_g, "-")):
                        nc(k, metric, unit, kla_basis, "n_removed_mgN", no_removal if not removed_ok else
                           "not computed: no N removed", OBS, tid)
        # X2: the idle fraction is a bounded quantity; every idle row states its coverage (classified fraction)
        idle_metrics = ((IDLE_LOWER, "idle_fraction_lower_bound", "bounds_status",
                         "LOWER_BOUND: inactive, over-supplied bins / all bins"),
                        (IDLE_UPPER, "idle_fraction_upper_bound", "bounds_status",
                         "UPPER_BOUND: (inactive + unclassified + artifact bins) / all bins"),
                        (IDLE_CONDITIONAL, "idle_fraction_of_classified_time", "status",
                         "CONDITIONAL on the classified bins being representative: inactive bins / classified bins"))
        if fractions is None:
            for metric in (*(m for m, *_ in idle_metrics), "classified_fraction"):
                nc(k, metric, "-", "W3 states", INPUTS["states_fractions"],
                   "not computed: states_fractions.csv missing (run the states stage)", OBS, tid)
        else:
            fr = fractions.set_index("interval").loc[f"pulse {k}"]
            cov = float(fr["classified_fraction"])
            states_txt = f"over {tk:g}-{tn:g} d (states use the slope of the M6-inverted C(t))"
            add(k, "classified_fraction", cov, "-", f"W3 coverage: bins with a process state / all bins {states_txt}",
                INPUTS["states_fractions"], obs=OBS, series=tid)
            for metric, col, status_col, what in idle_metrics:
                basis = (f"W3 idle fraction, {what} {states_txt}; coverage: classified fraction {cov:.3f} "
                         f"(conditional value needs ≥ {MIN_CLASSIFIED_FRACTION:g})")
                if str(fr[status_col]).startswith("computed"):
                    add(k, metric, float(fr[col]), "-", basis, INPUTS["states_fractions"], obs=OBS, series=tid)
                else:
                    nc(k, metric, "-", basis, INPUTS["states_fractions"], str(fr[status_col]), OBS, tid)

    def r8_case(k, tw):
        """X6: pulse k's consumption with the R8 electrode-fragment windows it governs (the primary rate excludes them)."""
        gov = pd.to_numeric(tw["governing_pulse"], errors="coerce") == k
        frag = tw[gov & tw["status"].astype(str).str.startswith(S.FRAGMENT_STATUS)]
        metric = "consumption_weighted_mgN_L_d_R8_CASE_incl_fragments"
        if len(frag) == 0:
            nc(k, metric, "mg N/L/d", "R8 CASE", INPUTS["truth_windows"],
               "not computed: no electrode-fragment window governed by this pulse (the primary rate is the whole case)",
               OBS, tid)
            return
        cols = ["sub_segment", "fit_start_d", "fit_end_d", "held_fraction", "consumption_mgN_L_d",
                "consumption_ci95_low_mgN_L_d", "consumption_ci95_high_mgN_L_d"]
        parts = pd.concat([tw.loc[gov & (tw["status"] == "computed"), cols], frag.assign(
            consumption_mgN_L_d=frag["consumption_fragment_diagnostic_mgN_L_d"],
            consumption_ci95_low_mgN_L_d=frag["consumption_fragment_diagnostic_ci95_low_mgN_L_d"],
            consumption_ci95_high_mgN_L_d=frag["consumption_fragment_diagnostic_ci95_high_mgN_L_d"])[cols]])
        parts = parts[parts["consumption_mgN_L_d"] > 0].assign(duration_d=lambda d: d["fit_end_d"] - d["fit_start_d"])
        frag_txt = ", ".join(f"{int(r.sub_segment)} ({r.fit_start_d:.3f}-{r.fit_end_d:.3f} d, "
                             f"{r.consumption_fragment_diagnostic_mgN_L_d:.2f} mg N/L/d)" for r in frag.itertuples(index=False))
        basis = ("R8 CASE (not the primary rate): duration-weighted consumption of the M6-inverted C(t) over the computed "
                 f"slope windows plus the electrode-fragment window(s) {frag_txt} demoted by R8; the primary "
                 "consumption_weighted_mgN_L_d excludes them; interval = duration-weighted Theil-Sen bounds, as for "
                 "the primary rate")
        src = (f"{INPUTS['truth_windows']} (sub-segments {', '.join(str(int(x)) for x in sorted(parts['sub_segment']))}); "
               f"{truth_src}")
        if len(parts) == 0:
            nc(k, metric, "mg N/L/d", basis, src, "not computed: no positive consumption in any window", OBS, tid)
        else:
            v, lo, hi, hf = weighted_rate(parts, "consumption_mgN_L_d", "consumption_ci95_low_mgN_L_d",
                                          "consumption_ci95_high_mgN_L_d")
            add(k, metric, v, "mg N/L/d", basis, src, ci=(lo, hi), held=hf, obs=OBS, series=tid)

    t_ph, ph_s = alkalinity.smoothed_ph(grid)
    pco2 = float(cfg.get("scenarios.pCO2_air_atm"))
    truth_src = f"{INPUTS['truth']} ({tid})"
    for k, (tk, tn) in enumerate(itertools.pairwise(times), start=1):
        st = steps.set_index("pulse").loc[k]
        add(k, "t_pulse_d", tk, "d", "experimental fact", "config/experiment.toml [pulses] times_d", obs=INDEP)
        add(k, "segment_end_d", tn, "d", "next pulse or end of common window", "config; common window 47.92 d",
            obs=INDEP)
        # --- affinity step: primary from C(t) just before the pulse; reading-as-TAN as SENSITIVITY
        c_pre = pre_median(t_tr, c_tr, tk)
        c_pre_b = [pre_median(t_tr, b, tk) for b in bands]
        for T_K, label in [(thermo.T_REF_K, "T25C"), *[(T + KELVIN_OFFSET, f"T{T:g}C")
                                                       for T in cfg.get("scenarios.temperature_C") if T != 25.0]]:
            v = float(thermo.affinity_step_kJ_mol(c_pre, dose.n_mg_L, T_K))
            bv = [float(thermo.affinity_step_kJ_mol(c, dose.n_mg_L, T_K)) for c in c_pre_b]
            add(k, f"delta_A_kJ_mol_SCENARIO_{label}", v, "kJ/mol",
                f"RT ln((C + ΔC)/C), C = M6-inverted NH4 over the 0.1 d before the pulse; SCENARIO temperature.{band_txt}",
                truth_src, obs=OBS, series=tid, band=(min(bv), max(bv)))
        for col in [c for c in steps.columns if c.startswith(f"{thermo.DELTA_A_READING}_SCENARIO_")]:
            add(k, col, float(st[col]), "kJ/mol",
                "SENSITIVITY: electrode reading taken as TAN (no observation model); SCENARIO temperature",
                INPUTS["steps"], obs=IDENT, series=RAW)

        w = governed_windows(seg, k)
        twk = tw[tw["sub_segment"].isin(w["sub_segment"])].merge(w[["sub_segment", "duration_d"]], on="sub_segment")
        rate_w = np.nan
        if len(w) == 0:
            for metric in ("consumption_uncorrected_first_window_mgN_L_d", "consumption_uncorrected_weighted_mgN_L_d"):
                nc(k, metric, "mg N/L/d", "", INPUTS["segments"], "not computed: no slope window after the pulse", RAW, RAW)
            for metric in ("consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d", "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"):
                nc(k, metric, "mg N/L/d", "SENSITIVITY", INPUTS["segments"], "not computed: no slope window after the pulse",
                   GAIN, RAW)
            for metric in ("consumption_first_window_mgN_L_d", "consumption_weighted_mgN_L_d"):
                nc(k, metric, "mg N/L/d", "", truth_src, "not computed: no slope window after the pulse", OBS, tid)
            r8_case(k, tw)
        else:
            first = w.iloc[0]
            pos = positive_windows(w)
            excluded = w[~w.index.isin(pos.index)]
            excl_note = (f"; excludes sub-segment(s) {', '.join(str(s) for s in excluded['sub_segment'])} "
                         "(non-positive gain-corrected consumption)" if len(excluded) else "")
            src = f"{INPUTS['segments']} (sub-segments {', '.join(str(s) for s in pos['sub_segment'])})"
            first_src = (f"{INPUTS['segments']} (sub-segment {first['sub_segment']}, "
                         f"{first['fit_start_d']:.3f}-{first['fit_end_d']:.3f} d)")
            # measured electrode slopes (no observation model)
            add(k, "consumption_uncorrected_first_window_mgN_L_d", -float(first["slope_mgN_L_d"]), "mg N/L/d",
                "measured (electrode reading, uncorrected); consumption = −slope; Theil-Sen 95 % interval", first_src,
                ci=(-float(first["slope_ci95_high_mgN_L_d"]), -float(first["slope_ci95_low_mgN_L_d"])),
                held=float(first["held_fraction"]), obs=RAW, series=RAW)
            # SENSITIVITY: linear gain model
            gain_basis = "SENSITIVITY: GAIN-CORRECTED (linear electrode gain model, S1)"
            if first["gain_correction_status"] == "computed":
                add(k, "consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d",
                    float(first[S.GAIN_COL]), "mg N/L/d", gain_basis,
                    f"{first_src}; {INPUTS['pulses']}",
                    ci=(float(first[S.GAIN_LO]),
                        float(first[S.GAIN_HI])), held=float(first["held_fraction"]),
                    obs=GAIN, series=RAW)
            else:
                nc(k, "consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d", "mg N/L/d", gain_basis, first_src,
                   str(first["gain_correction_status"]), GAIN, RAW, held=float(first["held_fraction"]))
            gain_w = np.nan
            if len(pos) == 0:
                nc(k, "consumption_uncorrected_weighted_mgN_L_d", "mg N/L/d", "duration-weighted over slope windows",
                   INPUTS["segments"], "not computed: non-positive gain-corrected consumption in every slope window", RAW, RAW)
                nc(k, "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d", "mg N/L/d", gain_basis + "; duration-weighted",
                   INPUTS["segments"], "not computed: non-positive gain-corrected consumption in every slope window",
                   GAIN, RAW)
            else:
                v, lo, hi, hf = weighted_rate(pos, "slope_mgN_L_d", "slope_ci95_low_mgN_L_d", "slope_ci95_high_mgN_L_d",
                                              sign=-1.0)
                add(k, "consumption_uncorrected_weighted_mgN_L_d", v, "mg N/L/d",
                    "measured (electrode reading, uncorrected); duration-weighted over slope windows with positive "
                    "gain-corrected consumption" + excl_note, src, ci=(lo, hi), held=hf, obs=RAW, series=RAW)
                gain_w, glo, ghi, ghf = weighted_rate(pos, S.GAIN_COL, S.GAIN_LO, S.GAIN_HI)
                add(k, "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d", gain_w, "mg N/L/d",
                    gain_basis + "; duration-weighted over slope windows with positive consumption" + excl_note,
                    f"{src}; {INPUTS['pulses']}", ci=(glo, ghi), held=ghf, obs=GAIN, series=RAW)
            # --- primary: the M6-inverted C(t) on the same windows
            f1 = twk.iloc[0]
            f1_src = f"{INPUTS['truth_windows']} (sub-segment {f1['sub_segment']}); {truth_src}"
            m6_basis = f"M6 log observation model, inverted C(t); Theil-Sen 95 % interval.{band_txt}"
            if f1["status"] == "computed":
                add(k, "consumption_first_window_mgN_L_d", float(f1["consumption_mgN_L_d"]), "mg N/L/d", m6_basis, f1_src,
                    ci=(float(f1["consumption_ci95_low_mgN_L_d"]), float(f1["consumption_ci95_high_mgN_L_d"])),
                    held=float(f1["held_fraction"]), obs=OBS, series=tid,
                    band=(float(f1["consumption_param_band_low_mgN_L_d"]), float(f1["consumption_param_band_high_mgN_L_d"])))
            else:
                nc(k, "consumption_first_window_mgN_L_d", "mg N/L/d", m6_basis, f1_src, str(f1["status"]), OBS, tid,
                   held=float(f1["held_fraction"]))
            tpos = twk[twk["status"] == "computed"]
            texcl = twk[twk["status"] != "computed"]
            t_note = (f"; excludes sub-segment(s) {', '.join(str(s) for s in texcl['sub_segment'])} "
                      f"({'; '.join(sorted(set(texcl['status'].str.removeprefix('not computed: '))))})" if len(texcl) else "")
            t_src = (f"{INPUTS['truth_windows']} (sub-segments {', '.join(str(s) for s in tpos['sub_segment'])}); "
                     f"{truth_src}")
            if len(tpos) == 0:
                nc(k, "consumption_weighted_mgN_L_d", "mg N/L/d", m6_basis, t_src,
                   "not computed: non-positive consumption of the M6-inverted NH4 in every slope window", OBS, tid)
            else:
                rate_w, lo, hi, hf = weighted_rate(tpos, "consumption_mgN_L_d", "consumption_ci95_low_mgN_L_d",
                                                   "consumption_ci95_high_mgN_L_d")
                d_ = tpos["duration_d"]
                band_w = (weighted(tpos["consumption_param_band_low_mgN_L_d"], d_),
                          weighted(tpos["consumption_param_band_high_mgN_L_d"], d_))
                add(k, "consumption_weighted_mgN_L_d", rate_w, "mg N/L/d",
                    "M6 log observation model, inverted C(t); duration-weighted over slope windows with positive "
                    f"consumption{t_note}.{band_txt}", t_src, ci=(lo, hi), held=hf, obs=OBS, series=tid, band=band_w)
                for v in vols:
                    if v.is_constant:
                        continue
                    col = f"consumption_mgN_L_d_{volume.slug(v.name)}"
                    shift = weighted(tpos[col] - tpos["consumption_mgN_L_d"], d_)  # the C E/V term
                    add(k, f"consumption_weighted_mgN_L_d_{volume.slug(v.name)}", rate_w + shift,
                        "mg N/L/d", f"{v.basis}: r = −dC/dt + C E/V(t) on the M6-inverted C(t); duration-weighted; "
                        "interval = the constant-volume interval shifted by the C E/V term", t_src, ci=(lo + shift, hi + shift),
                        held=hf, obs=OBS, series=tid)

            r8_case(k, tw)

            mask = np.zeros(len(grid), dtype=bool)
            for r in w.itertuples(index=False):
                mask |= (grid["t_start_d"] >= r.fit_start_d - 1e-9).to_numpy() & (grid["t_end_d"] <= r.fit_end_d + 1e-9).to_numpy()
            gsel = grid[mask]
            add(k, "median_pH", float(gsel["pH"].median()), "pH", "measured", f"{INPUTS['grid']} over the slope windows",
                obs=INDEP)
            tmask = np.zeros(len(truth), dtype=bool)
            for r in w.itertuples(index=False):
                tmask |= truth.index.isin(S.truth_bins(truth, r.fit_start_d, r.fit_end_d).index)
            tsel = truth[tmask]
            for T in cfg.get("scenarios.temperature_C"):
                T_K = T + KELVIN_OFFSET
                fa = thermo.free_ammonia_mgNH3_L(tsel["C_inv_mgN_L"].to_numpy(), tsel["pH"].to_numpy(), T_K)
                fb = [float(np.median(thermo.free_ammonia_mgNH3_L(b[tmask], tsel["pH"].to_numpy(), T_K))) for b in bands]
                add(k, f"fa_median_mgNH3_L_SCENARIO_T{T:g}C", float(np.median(fa)), "mg NH3/L",
                    f"SCENARIO temperature; TAN = M6-inverted C(t) over the slope windows.{band_txt}",
                    f"{truth_src}; config [scenarios]", obs=OBS, series=tid, band=(min(fb), max(fb)))
                fa_r = thermo.free_ammonia_mgNH3_L(gsel["nh4_mgNH4_L"].to_numpy() * nf, gsel["pH"].to_numpy(), T_K)
                add(k, f"fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L_SCENARIO_T{T:g}C", float(np.median(fa_r)), "mg NH3/L",
                    "SENSITIVITY: electrode reading taken as TAN (no observation model); SCENARIO temperature",
                    f"{INPUTS['grid']}; config [scenarios]", obs=IDENT, series=RAW)

            # --- oxygen demand: primary from the inverted C(t); linear gain as SENSITIVITY
            if rate_w > 0:
                d = float(thermo.oxygen_demand_mgO2_L_d(rate_w))
                o2_ci = (float(thermo.oxygen_demand_mgO2_L_d(max(lo, 0.0))), float(thermo.oxygen_demand_mgO2_L_d(max(hi, 0.0))))
                o2_band = tuple(float(thermo.oxygen_demand_mgO2_L_d(max(x, 0.0))) for x in band_w)
                kla = float(thermo.required_kla_per_d(d))
                add(k, "o2_demand_mgO2_L_d", d, "mg O2/L/d", f"4.57 g O2/g N x M6-inverted consumption.{band_txt}",
                    "consumption_weighted_mgN_L_d; constants.O2_PER_N_G_G", ci=o2_ci, held=hf, obs=OBS, series=tid,
                    band=o2_band)
                add(k, "kla_required_per_d_SCENARIO_T25C", kla, "1/d",
                    "lower bound k_La >= demand / C*; C* at SCENARIO 25 °C", "o2_demand_mgO2_L_d; constants.DO_SAT_25C_MG_L",
                    ci=tuple(float(thermo.required_kla_per_d(x)) for x in o2_ci), held=hf, obs=OBS, series=tid,
                    band=tuple(float(thermo.required_kla_per_d(x)) for x in o2_band))
                margin = thermo.kla_margin_guarded(cfg, kla)
                if isinstance(margin, thermo.NotComputed):
                    nc(k, "kla_measured_over_required", "-", "", "kLa_per_d (config [missing])", str(margin), OBS, tid)
                else:
                    add(k, "kla_measured_over_required", float(margin), "-", "measured k_La / required (SCENARIO 25 °C)",
                        "config [measurements] kLa_per_d; kla_required_per_d_SCENARIO_T25C", obs=OBS, series=tid)
            else:
                for metric, unit in (("o2_demand_mgO2_L_d", "mg O2/L/d"), ("kla_required_per_d_SCENARIO_T25C", "1/d"),
                                     ("kla_measured_over_required", "-")):
                    nc(k, metric, unit, "M6 log observation model", "consumption_weighted_mgN_L_d",
                       "not computed: non-positive consumption of the M6-inverted NH4", OBS, tid)
            if gain_w > 0:
                gd = float(thermo.oxygen_demand_mgO2_L_d(gain_w))
                add(k, "o2_demand_SENSITIVITY_linear_gain_mgO2_L_d", gd, "mg O2/L/d",
                    "SENSITIVITY: GAIN-CORRECTED; 4.57 g O2/g N", "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d",
                    ci=(float(thermo.oxygen_demand_mgO2_L_d(max(glo, 0.0))), float(thermo.oxygen_demand_mgO2_L_d(max(ghi, 0.0)))),
                    held=ghf, obs=GAIN, series=RAW)
                add(k, "kla_required_SENSITIVITY_linear_gain_per_d_SCENARIO_T25C", float(thermo.required_kla_per_d(gd)),
                    "1/d", "SENSITIVITY: GAIN-CORRECTED; C* at SCENARIO 25 °C", "o2_demand_SENSITIVITY_linear_gain_mgO2_L_d",
                    ci=(float(thermo.required_kla_per_d(thermo.oxygen_demand_mgO2_L_d(max(glo, 0.0)))),
                        float(thermo.required_kla_per_d(thermo.oxygen_demand_mgO2_L_d(max(ghi, 0.0))))),
                    held=ghf, obs=GAIN, series=RAW)
            else:
                for metric, unit in (("o2_demand_SENSITIVITY_linear_gain_mgO2_L_d", "mg O2/L/d"),
                                     ("kla_required_SENSITIVITY_linear_gain_per_d_SCENARIO_T25C", "1/d")):
                    nc(k, metric, unit, "SENSITIVITY: GAIN-CORRECTED", "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d",
                       "not computed: non-positive gain-corrected consumption (rising electrode signal)", GAIN, RAW)
            add(k, "o2_windows_computed", len(tpos), "count", "slope windows with positive M6-inverted consumption",
                INPUTS["truth_windows"], obs=OBS, series=tid)
        # --- R7: a non-equilibrium measure read from the trajectory itself: the time to remove one pulse at the
        # post-pulse rate, against the time to the next pulse
        if rate_w > 0:
            tau = dose.n_mg_L / rate_w
            add(k, "dose_removal_time_d", tau, "d", "dose / duration-weighted consumption of the M6-inverted C(t): time to "
                "remove one pulse at the post-pulse rate (a ratio of two observed quantities, not a fitted relaxation "
                f"constant).{band_txt}", "consumption_weighted_mgN_L_d; dose", obs=OBS, series=tid)
            add(k, "dose_removal_time_over_pulse_interval", tau / (tn - tk), "-",
                "dose_removal_time_d / time to the next pulse (> 1: the next pulse arrives before the previous one could "
                "be removed at the post-pulse rate)", "dose_removal_time_d; config", obs=OBS, series=tid)
        else:
            for metric, unit in (("dose_removal_time_d", "d"), ("dose_removal_time_over_pulse_interval", "-")):
                nc(k, metric, unit, "M6 log observation model", "consumption_weighted_mgN_L_d",
                   "not computed: no positive consumption of the M6-inverted NH4 after this pulse", OBS, tid)

        # --- conversion (S3): the dose-anchored removal of the inverted C(t) between pulses, as a ratio to the dose;
        # η, a fraction of this pulse, exists only for 0 ≤ ratio ≤ 1 (a ratio > 1 means carry-over TAN was removed too)
        c_next = pre_median(t_tr, c_tr, tn) if tn < END_D else float(np.median(c_tr[t_tr >= END_D - PRE_WINDOW_D]))
        b_next = [pre_median(t_tr, b, tn) if tn < END_D else float(np.median(b[t_tr >= END_D - PRE_WINDOW_D]))
                  for b in bands]
        ratio = (c_pre + dose.n_mg_L - c_next) / dose.n_mg_L
        ratio_b = [(a + dose.n_mg_L - b) / dose.n_mg_L for a, b in zip(c_pre_b, b_next)]
        ratio_basis = (f"(C_before + ΔC − C_before next) / ΔC on the M6-inverted C(t) (0.1 d medians), constant volume."
                       f"{band_txt}")
        eta_src = f"{truth_src} over {tk:g}-{tn:g} d; dose {dose.n_mg_L:.2f} mg N/L"
        n_removed = np.nan
        cases = [(None, ratio, ratio_b, ratio_basis, eta_src)]
        for v in vols:
            if not v.is_constant:
                rv = (c_pre * float(v.volume_L(tk)) + n_mg - c_next * float(v.volume_L(tn))) / n_mg
                cases.append((v, rv, None, (f"{v.basis}: amount balance (C_before V(t_k) + m_N − C_before next V(t_next)) "
                                            "/ m_N on the M6-inverted C(t)"), eta_src))
        for vname, r_, rb, basis_, src_ in cases:
            sfx = "" if vname is None else f"_{volume.slug(vname.name)}"
            band = (min(rb), max(rb)) if rb else (np.nan, np.nan)
            if r_ < 0:
                msg = (f"not computed: negative dose-anchored removal ({r_:.3f}): the M6-inverted C(t) rose by more than "
                       "the dose" + (", which no N source allows (inversion inconsistent over this interval)"
                                     if vname is None else " under this volume scenario"))
                nc(k, f"removal_to_dose_ratio{sfx}", "-", basis_, src_, msg, OBS, tid)
                nc(k, f"eta_conversion_fraction{sfx}", "-", basis_, src_, msg, OBS, tid)
                continue
            add(k, f"removal_to_dose_ratio{sfx}", r_, "-", basis_ + " Removal of TAN between the pulses per dose; > 1 "
                "when carry-over TAN is removed as well.", src_, obs=OBS, series=tid, band=band)
            if r_ <= 1:
                eta_band, eta_note = eta_band_on_unit_interval(band)
                add(k, f"eta_conversion_fraction{sfx}", r_, "-", basis_ + eta_note, src_, obs=OBS, series=tid,
                    band=eta_band)
            else:
                nc(k, f"eta_conversion_fraction{sfx}", "-", basis_, src_,
                   f"not computed: removal exceeds the dose (ratio {r_:.3f}): carry-over TAN was removed as well, so a "
                   "conversion fraction of this pulse is not defined (see removal_to_dose_ratio)", OBS, tid)
            if vname is None:
                n_removed = r_ * n_mg
        if np.isfinite(n_removed):
            add(k, "n_removed_mgN", n_removed, "mg N", "removal_to_dose_ratio x N added by one pulse (constant volume)",
                "removal_to_dose_ratio; config [pulses]", obs=OBS, series=tid)
            add(k, "free_energy_released_J_SCENARIO_T25C", float(energy.gibbs_released_J(n_removed)), "J",
                "N removed x |ΔG°′| (standard conditions, SCENARIO 25 °C)", "n_removed_mgN; constants",
                obs=OBS, series=tid)
        else:
            msg = "not computed: negative dose-anchored removal on the M6-inverted C(t)"
            nc(k, "n_removed_mgN", "mg N", ratio_basis, eta_src, msg, OBS, tid)
            nc(k, "free_energy_released_J_SCENARIO_T25C", "J", "N removed x |ΔG°′|", "n_removed_mgN", msg, OBS, tid)
        energy_family(k, tk, tn, n_removed)
        # --- alkalinity needed for the cumulative removal implied by C(t) (W1 against the W2 bound)
        i_next = int(np.searchsorted(t_ph, tn, side="left")) - 1 if tn < END_D else len(t_ph) - 1
        amin = float(alkalinity.alk_min_meq_L(ph_s[i_next], thermo.T_REF_K, pco2))
        removed = c_start + k * dose.n_mg_L - c_next
        for metric, eq, what in (("alk0_required_if_nitrified_meq_L", alkalinity.EQ_PER_MOL_N_NITRIFICATION, "nitrification"),
                                 ("alk0_required_any_proton_sink_meq_L", alkalinity.EQ_PER_MOL_N_ANY_SINK,
                                  "any proton-releasing N sink")):
            add(k, metric, eq * max(removed, 0.0) / M_N_G_MOL + amin, "meq/L",
                f"cumulative removal C(t0) + {k} doses − C(t) = {removed:.1f} mg N/L by {tn:g} d on the M6-inverted C(t), "
                f"{eq:g} eq/mol N ({what}), + Alk_min at the smoothed pH (air-equilibrium CO2, SCENARIO 25 °C)",
                f"{truth_src}; alkalinity.py", obs=OBS, series=tid)

        integral, mean_r = model_integral(state, tk, tn)
        mb_basis = (f"MODEL-BASED ({selected} rate law, whole-record refit; same M6 observation model; comparison with "
                    "the dose-anchored η; no uncertainty interval)")
        mb_src = f"{INPUTS['model']} ∫r dt over {tk:g}-{tn:g} d / {dose.n_mg_L:.2f} mg N/L"
        if 0 <= integral / dose.n_mg_L <= 1:
            add(k, "eta_conversion_fraction_MODEL_BASED", integral / dose.n_mg_L, "-", mb_basis, mb_src, obs=OBS,
                series=state_series)
        else:
            nc(k, "eta_conversion_fraction_MODEL_BASED", "-", mb_basis, mb_src,
               f"not computed: MODEL-BASED removal/dose {integral / dose.n_mg_L:.3f} outside [0, 1]", OBS, state_series)
        # X7 (e): the model state obeys C_model(next) = C_model(before) + dose − ∫r dt, so the dose-anchored ratio minus
        # the MODEL-BASED η is the decline over the interval of the misfit C_inv − C_model, divided by the dose
        t_st, c_st = state["t_d"].to_numpy(), state["model_state_C_mgN_L_conditional_on_C0"].to_numpy()

        def at_(t_, c_, x):
            return pre_median(t_, c_, x) if x < END_D else float(np.median(c_[t_ >= END_D - PRE_WINDOW_D]))
        mis_pre, mis_next = c_pre - at_(t_st, c_st, tk), c_next - at_(t_st, c_st, tn)
        add(k, "eta_gap_inversion_misfit_MODEL_BASED", (mis_pre - mis_next) / dose.n_mg_L, "-",
            f"MODEL-BASED decomposition: [(C_inv − C_model) before the pulse ({mis_pre:.1f} mg N/L) − (C_inv − C_model) "
            f"before the next ({mis_next:.1f} mg N/L)] / dose. Because the model state obeys its own mass balance, this "
            "equals removal_to_dose_ratio − eta_conversion_fraction_MODEL_BASED up to the 0.1 d medians: the part of the "
            "dose-anchored removal that the rate law does not reproduce",
            f"{truth_src}; {INPUTS['model']}", obs=OBS, series=state_series)
        add(k, "mean_consumption_rate_MODEL_BASED_mgN_L_d", mean_r, "mg N/L/d",
            f"MODEL-BASED ({selected}); no interval (joint parameter uncertainty not propagated)", INPUTS["model"],
            obs=OBS, series=state_series)
        add(k, "free_energy_per_fully_converted_pulse_J_SCENARIO_T25C", g0, "J",
            "SCENARIO 25 °C: standard conditions (ΔG°′ = −348.8 kJ/mol N at pH 7 and 25 °C)",
            "config [pulses], [vessel]; constants", obs=INDEP)
        add(k, "free_energy_released_MODEL_BASED_J_SCENARIO_T25C", integral / dose.n_mg_L * g0, "J",
            "MODEL-BASED η x standard-condition free energy (SCENARIO 25 °C)", "eta_conversion_fraction_MODEL_BASED; free energy per pulse",
            obs=OBS, series=state_series)
        nc(k, "eta_uncertainty_interval", "-", "MODEL-BASED η", "model_identifiability.csv",
           "not computed: joint parameter uncertainty of the MODEL-BASED η is not propagated (the dose-anchored η "
           "carries the profile-vector band)", OBS, state_series)
        if stirrer_w is None:
            for metric, unit in (("specific_stirring_energy_J_per_gN", "J/g N"),
                                 ("free_energy_to_stirrer_energy_ratio", "-")):
                nc(k, metric, unit, "LABORATORY-SPECIFIC; MODEL-BASED", "stirrer_power_W (config [missing])",
                   stirrer_status, OBS, state_series)
        else:
            energy_j = stirrer_energy_J(stirrer_w, tk, tn)
            g_n = integral * cfg.get("vessel.volume_L") / 1e3  # mg N/L x L -> g N
            add(k, "stirrer_energy_J", energy_j, "J", "LABORATORY-SPECIFIC: measured power x segment duration",
                "config [measurements] stirrer_power_W; segment duration", obs=INDEP)
            if g_n > 0:
                add(k, "specific_stirring_energy_J_per_gN", energy_j / g_n, "J/g N",
                    "LABORATORY-SPECIFIC: stirrer energy / MODEL-BASED nitrogen converted", "stirrer_energy_J; eta_conversion_fraction_MODEL_BASED; vessel volume",
                    obs=OBS, series=state_series)
            else:
                nc(k, "specific_stirring_energy_J_per_gN", "J/g N", "LABORATORY-SPECIFIC; MODEL-BASED",
                   "eta_conversion_fraction_MODEL_BASED",
                   "not computed: no nitrogen converted in the model", OBS, state_series)
            add(k, "free_energy_to_stirrer_energy_ratio", integral / dose.n_mg_L * g0 / energy_j, "-",
                "LABORATORY-SPECIFIC: MODEL-BASED free energy released (SCENARIO 25 °C) / stirrer energy",
                "free_energy_released_MODEL_BASED_J_SCENARIO_T25C; stirrer_energy_J", obs=OBS, series=state_series)
    out = pd.DataFrame(rows)
    model.check_single_observation_model(out)
    return out


def truth_balance(cfg: Config, inputs: dict) -> pd.DataFrame:
    """Dose-anchored TAN bound (W1): C_up(t) = (C(t0) V(t0) + Σ m_N)/V(t), the most TAN there can be, and
    C_low(t; Alk0) = C_up(t) − N_nit,max(t; Alk0) V0/V(t), the least TAN there can be if nitrification is the only
    sink (W2 bound, air-equilibrium CO2, SCENARIO 25 °C), against the M6-inverted C(t) at each checkpoint."""
    truth, grid = inputs["truth"], inputs["grid"]
    bands = band_series(truth, inputs["truth_vectors"])
    t_tr, c_tr = truth["t_d"].to_numpy(), truth["C_inv_mgN_L"].to_numpy()
    tid = str(truth["nh4_truth_id"].iloc[0])
    t, ph = alkalinity.smoothed_ph(grid)
    pco2 = float(cfg.get("scenarios.pCO2_air_atm"))
    n_mg, v0 = volume.pulse_n_mg(cfg), float(cfg.get("vessel.volume_L"))
    alk0s = (2.0, 5.0, 10.0, 15.0)
    c0_m = t_tr < t_tr[0] + PRE_WINDOW_D
    rows = []
    for scen in volume.volume_scenarios(cfg, END_D):
        curves = {a: alkalinity.bound_curve(t, ph, a, thermo.T_REF_K, pco2, scen) for a in alk0s}
        c0 = float(np.median(c_tr[c0_m]))
        for name, i in alkalinity.checkpoints(cfg, t):
            ti = float(t[i])
            m = (t_tr > ti - PRE_WINDOW_D + 1e-9) & (t_tr <= ti + 1e-9)
            c_now = float(np.median(c_tr[m]))
            cb = [float(np.median(b[m])) for b in bands]
            c0b = [float(np.median(b[c0_m])) for b in bands]
            added = n_mg * sum(tp <= ti for tp in cfg.pulse_times_d)
            v_t, v_s = float(scen.volume_L(ti)), float(scen.volume_L(t_tr[0]))
            upper = (c0 * v_s + added) / v_t
            removed = (c0 * v_s + added - c_now * v_t) / v0
            amin_scaled = float(alkalinity.alk_min_meq_L(ph[i], thermo.T_REF_K, pco2)) * v_t / v0
            row = {"volume_scenario": scen.name, "volume_basis": scen.basis, "checkpoint": name, "t_d": ti,
                   "C_inv_mgN_L": c_now, "C_inv_band_low_mgN_L": min(cb), "C_inv_band_high_mgN_L": max(cb),
                   "C_dose_upper_mgN_L": upper, "C_inv_within_upper": bool(c_now <= upper + 1e-9),
                   "removed_per_initial_litre_mgN_L": removed,
                   "alk0_required_if_nitrified_meq_L": 2 * max(removed, 0.0) / M_N_G_MOL + amin_scaled,
                   "alk0_required_any_proton_sink_meq_L": max(removed, 0.0) / M_N_G_MOL + amin_scaled}
            # the same requirement for every profile vector of the band (each vector shifts C(t0) and C(t) together)
            rem_b = [(a0 * v_s + added - cn * v_t) / v0 for a0, cn in zip(c0b, cb)]
            req_b = [2 * max(r_, 0.0) / M_N_G_MOL + amin_scaled for r_ in rem_b]
            row.update({"removed_band_low_mgN_L": min(rem_b), "removed_band_high_mgN_L": max(rem_b),
                        "alk0_required_if_nitrified_band_min_meq_L": min(req_b),
                        "alk0_required_if_nitrified_band_max_meq_L": max(req_b),
                        "alk0_required_any_proton_sink_band_min_meq_L":
                            min(max(r_, 0.0) / M_N_G_MOL + amin_scaled for r_ in rem_b)})
            for a in alk0s:
                low = upper - float(curves[a][i]) * v0 / v_t
                row[f"C_dose_lower_alk0_{a:g}_mgN_L"] = low
                row[f"C_inv_within_bounds_alk0_{a:g}"] = bool(low - 1e-9 <= c_now <= upper + 1e-9)
            row.update({"pH_smoothed": float(ph[i]), "C_t0_mgN_L": c0,
                        "basis": "W1 dose-anchored bound; C(t0) from the same inverted C(t); W2 bound at air-equilibrium "
                                 "CO2, SCENARIO 25 °C; nitrification as the only N sink for C_dose_lower",
                        "nh4_observation_model": model.OBSERVATION_MODEL, "nh4_truth_id": tid})
            rows.append(row)
    return pd.DataFrame(rows)


# --- S2: the inverted C(t) against the W2 alkalinity bound -------------------------------------------------------
ALK0_SCAN_MAX = max(alkalinity.ALK0_GRID_MEQ_L)
QUALIFIED = "QUALIFIED"


class BoundContradictionError(ValueError):
    """A quantity derived from a C(t) that contradicts the W2 bound is not qualified as such."""


def w2_consistency(balance: pd.DataFrame) -> pd.DataFrame:
    """Per checkpoint (constant volume, the primary frame): is the cumulative removal implied by the inverted C(t)
    compatible with the W2 bound for an Alk0 inside the scanned 1-15 meq/L? Three readings when it is not: the inverted
    level is too low (C0 and a are not identifiable and the band is open upward), Alk0 exceeds the scan, or a sink
    releasing 1 eq/mol N (NH3 stripping, uptake) or none took part."""
    b = balance[balance["volume_scenario"] == volume.CONSTANT].copy()
    rows = []
    for r in b.itertuples(index=False):
        best, bmin, anys = (r.alk0_required_if_nitrified_meq_L, r.alk0_required_if_nitrified_band_min_meq_L,
                            r.alk0_required_any_proton_sink_meq_L)
        if best <= ALK0_SCAN_MAX:
            verdict, ok = f"consistent: the removal needs Alk0 ≥ {best:.1f} meq/L if nitrified", True
        else:
            where = ("only within the open upper part of the C(t) parameter band "
                     f"(≥ {bmin:.1f} meq/L there)" if bmin <= ALK0_SCAN_MAX else
                     f"at no vector of the C(t) parameter band (≥ {bmin:.1f} meq/L)")
            verdict, ok = (f"{QUALIFIED}: contradicts the W2 bound for every scanned Alk0 (1-{ALK0_SCAN_MAX:g} meq/L) at "
                           f"the best fit (needs ≥ {best:.1f} meq/L if nitrified); consistent {where}, or with a 1 eq/mol "
                           f"N sink (≥ {anys:.1f} meq/L), or if Alk0 > {ALK0_SCAN_MAX:g} meq/L"), False
        rows.append({"checkpoint": r.checkpoint, "t_d": r.t_d, "C_inv_mgN_L": r.C_inv_mgN_L,
                     "removed_per_initial_litre_mgN_L": r.removed_per_initial_litre_mgN_L,
                     "removed_band_low_mgN_L": r.removed_band_low_mgN_L, "removed_band_high_mgN_L": r.removed_band_high_mgN_L,
                     "alk0_required_if_nitrified_meq_L": best, "alk0_required_if_nitrified_band_min_meq_L": bmin,
                     "alk0_required_any_proton_sink_meq_L": anys, "alk0_scan_max_meq_L": ALK0_SCAN_MAX,
                     "consistent_at_best_fit": ok, "verdict": verdict})
    return pd.DataFrame(rows)


def pulse_checkpoints(k: int, n_pulses: int) -> tuple[str, str]:
    """The checkpoints bounding pulse k's interval (its start and end)."""
    return f"pre_pulse{k}", (f"pre_pulse{k + 1}" if k < n_pulses else "end")


def pulse_qualifier(consistency: pd.DataFrame, k: int, n_pulses: int) -> str:
    """Qualifier of pulse k from the checkpoints bounding its interval; k = 0 is the interval before the first pulse
    (baseline windows), bounded by checkpoint pre_pulse1 (X5)."""
    c = consistency.set_index("checkpoint")
    pts = [x for x in (("pre_pulse1",) if k == 0 else pulse_checkpoints(k, n_pulses)) if x in c.index]
    bad = [x for x in pts if not bool(c.loc[x, "consistent_at_best_fit"])]
    if bad:
        return f"{c.loc[bad[0], 'verdict']} (checkpoint {bad[0]}; nh4_truth_w2_consistency.csv)"
    return "; ".join(f"{x}: {c.loc[x, 'verdict']}" for x in pts)


def qualify(df: pd.DataFrame, consistency: pd.DataFrame, n_pulses: int, pulse_col: str, model_col: str) -> pd.DataFrame:
    """Add ``w2_bound_check`` to every row derived from the M6 observation model (the qualifier of its pulse)."""
    out = df.copy()
    out["w2_bound_check"] = ""
    pulse = pd.to_numeric(out[pulse_col], errors="coerce").fillna(0)  # no governing pulse: before pulse 1
    for k in range(n_pulses + 1):
        m = (pulse == k) & (out[model_col] == model.OBSERVATION_MODEL)
        out.loc[m, "w2_bound_check"] = pulse_qualifier(consistency, k, n_pulses)
    return out


def check_w2_qualified(df: pd.DataFrame, consistency: pd.DataFrame, n_pulses: int, pulse_col: str = "pulse",
                       model_col: str = "nh4_observation_model") -> None:
    """Raise unless every M6-derived row of a pulse whose C(t) contradicts the W2 bound carries the QUALIFIED label."""
    if "w2_bound_check" not in df.columns:
        raise BoundContradictionError("no w2_bound_check column")
    pulse = pd.to_numeric(df[pulse_col], errors="coerce").fillna(0)
    for k in range(n_pulses + 1):
        if pulse_qualifier(consistency, k, n_pulses).startswith(QUALIFIED):
            m = (pulse == k) & (df[model_col] == model.OBSERVATION_MODEL)
            if not df.loc[m, "w2_bound_check"].astype(str).str.startswith(QUALIFIED).all():
                raise BoundContradictionError(f"pulse {k}: C(t)-derived rows not qualified against the W2 bound")


def volume_constraints(cfg: Config, inputs: dict) -> pd.DataFrame:
    """W8: data constraints on the net water loss rate and the evaporation SCENARIOS checked against them."""
    truth, grid = inputs["truth"], inputs["grid"]
    bands = band_series(truth, inputs["truth_vectors"])
    t_tr, c_tr = truth["t_d"].to_numpy(), truth["C_inv_mgN_L"].to_numpy()
    n_mg, v0 = volume.pulse_n_mg(cfg), float(cfg.get("vessel.volume_L"))
    times = [*cfg.pulse_times_d, END_D]

    def at(c, t_):
        return pre_median(t_tr, c, t_) if t_ < END_D else float(np.median(c[t_tr >= END_D - PRE_WINDOW_D]))

    rows = []
    per_vec = np.zeros(len(bands))
    best_max = 0.0
    for k, (tk, tn) in enumerate(itertools.pairwise(times), start=1):
        e = volume.n_conservation_min_loss_mL_d(at(c_tr, tk), at(c_tr, tn), tk, tn, n_mg, v0)
        eb = [volume.n_conservation_min_loss_mL_d(at(b, tk), at(b, tn), tk, tn, n_mg, v0) for b in bands]
        per_vec = np.maximum(per_vec, eb)
        best_max = max(best_max, e)
        rows.append({"constraint": f"N conservation over pulse {k} ({tk:g}-{tn:g} d)", "bound": "lower",
                     "value_mL_per_d": e, "band_low_mL_per_d": min(eb), "band_high_mL_per_d": max(eb),
                     "basis": "CONDITIONAL on the M6-inverted C(t): smallest net loss for which the amount balance "
                              "removes N ≥ 0 (no N source)"})
    rows.append({"constraint": "N conservation over every pulse", "bound": "lower", "value_mL_per_d": best_max,
                 "band_low_mL_per_d": float(per_vec.min()), "band_high_mL_per_d": float(per_vec.max()),
                 "basis": "CONDITIONAL on the M6-inverted C(t): max over the pulses"})
    t_last = float(grid["t_end_d"].iloc[-1])
    upper = volume.evaporation_upper_bound_mL_d(v0, t_last)
    rows.append({"constraint": "liquid present at the last sensor reading", "bound": "upper", "value_mL_per_d": upper,
                 "band_low_mL_per_d": np.nan, "band_high_mL_per_d": np.nan,
                 "basis": f"experimental facts: V0 / {t_last:g} d (nothing else was added)"})
    ec = grid["EC_mS_cm"]
    ec0 = float(ec[(grid["t_mid_d"] >= 0.3) & (grid["t_mid_d"] < 1.3)].median())
    ec1 = float(ec[(grid["t_mid_d"] >= END_D - 1.0) & (grid["t_mid_d"] < END_D)].median())
    try:
        ec_max = volume.ec_conditional_max_loss_mL_d(ec0, ec1, v0, END_D - 0.8)
        ec_status = "computed (CONDITIONAL; not used: the probe failed M2)"
    except ValueError as exc:
        ec_max, ec_status = np.nan, f"not computed: {exc}"
    rows.append({"constraint": f"flat conductivity record ({ec0:.4g} to {ec1:.4g} mS/cm)", "bound": "upper",
                 "value_mL_per_d": ec_max, "band_low_mL_per_d": np.nan, "band_high_mL_per_d": np.nan,
                 "basis": "CONDITIONAL on the conductivity probe measuring solution conductivity, which M2 rejected "
                          f"(no response to the pulses); {ec_status}"})
    out = pd.DataFrame(rows)
    for scen in volume.volume_scenarios(cfg, END_D):
        e = 1e3 * scen.net_loss_L_per_d
        checks = [f"{'≥' if e >= best_max - 1e-12 else '<'} N-conservation minimum {best_max:.3g}",
                  f"{'<' if e < upper else '≥'} liquid-present maximum {upper:.3g}"]
        ok = e >= best_max - 1e-12 and e < upper
        out = pd.concat([out, pd.DataFrame([{
            "constraint": f"scenario check: {scen.name}", "bound": "scenario", "value_mL_per_d": e,
            "band_low_mL_per_d": np.nan, "band_high_mL_per_d": np.nan,
            "basis": f"{'consistent with' if ok else 'INCONSISTENT with'} the constraints: net loss {e:g} mL/d "
                     f"{checks[0]} and {checks[1]} mL/d (the N-conservation constraint is CONDITIONAL on C(t))"}])],
            ignore_index=True)
    return out


def wide(metrics: pd.DataFrame, values: str = "value") -> pd.DataFrame:
    w = metrics.pivot(index="pulse", columns="metric", values=values)
    return w.reset_index()


RESPONSE = "consumption_weighted_mgN_L_d"  # W1: the M6-inverted consumption


def associations(metrics: pd.DataFrame) -> pd.DataFrame:
    """Spearman ρ of the M6-inverted consumption against candidate descriptors, n = 4: exploratory only."""
    w = wide(metrics).set_index("pulse")
    y = w[RESPONSE]
    rows = []
    for x in ("delta_A_kJ_mol_SCENARIO_T25C", *[c for c in w.columns if c.startswith("fa_median_mgNH3_L_SCENARIO_")],
              "median_pH", "t_pulse_d"):
        ok = y.notna() & w[x].notna()
        rho = float(spearmanr(w.loc[ok, x], y[ok]).statistic) if ok.sum() >= 3 else np.nan
        rows.append({"response": RESPONSE, "descriptor": x, "n": int(ok.sum()),
                     "spearman_rho": rho, "status": EXPLORATORY,
                     "note": "no p-value: with 4 pulses the smallest attainable two-sided p is 0.083"})
    return pd.DataFrame(rows)


@dataclass
class PulseResult:
    metrics: pd.DataFrame
    associations: pd.DataFrame
    truth_windows: pd.DataFrame | None = None
    truth_balance: pd.DataFrame | None = None
    consistency: pd.DataFrame | None = None
    volume_constraints: pd.DataFrame | None = None


def compute(cfg: Config, inputs: dict) -> PulseResult:
    tw = truth_windows(cfg, inputs)
    m = pulse_metrics(cfg, {**inputs, "truth_windows": tw})
    bal = truth_balance(cfg, inputs)
    cons = w2_consistency(bal)
    n = len(cfg.pulse_times_d)
    m = qualify(m, cons, n, "pulse", "nh4_observation_model")
    tw = qualify(tw, cons, n, "governing_pulse", "nh4_observation_model")
    check_w2_qualified(m, cons, n)
    check_w2_qualified(tw, cons, n, "governing_pulse")
    model.check_single_observation_model(m)  # X3: the qualified tables are checked again right before they are written
    model.check_single_observation_model(tw, label_cols=("status", "w2_bound_check"), series_col="nh4_truth_id")
    return PulseResult(m, associations(m), tw, bal, cons, volume_constraints(cfg, inputs))


INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"


def figure(cfg: Config, r: PulseResult, png: Path, svg: Path) -> None:
    from .report import pyplot

    plt = pyplot()

    w = wide(r.metrics)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), facecolor=SURFACE)
    for ax in axes.flat:
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(True, color=GRID, linewidth=0.6, axis="y")
        ax.set_axisbelow(True)
    x = np.arange(len(w))
    labels = [f"P{int(p)}" for p in w["pulse"]]
    ax = axes[0, 0]
    lo, hi = wide(r.metrics, "ci95_low"), wide(r.metrics, "ci95_high")
    for off, col, color, lab in ((-0.2, "consumption_uncorrected_weighted_mgN_L_d", BLUE, "Electrode reading (measured)"),
                                 (0.2, RESPONSE, ORANGE, "M6-inverted NH4 (single observation model)")):
        v = w[col].to_numpy(dtype=float)
        err = np.vstack([v - lo[col].to_numpy(dtype=float), hi[col].to_numpy(dtype=float) - v])
        ax.bar(x + off, v, 0.38, color=color, label=lab, yerr=err, error_kw={"ecolor": INK2, "capsize": 3, "lw": 0.8})
    ax.scatter(x + 0.2, w["consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"], marker="_", s=160, color=INK, zorder=3,
               label="SENSITIVITY: linear gain correction")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Consumption (mg N/L/d)", fontsize=8.5, color=INK)
    ax.set_title("a  Consumption after each pulse (positive slope windows, 95 % CI)", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7.5, frameon=False)
    ax = axes[0, 1]
    ax.bar(x, w["delta_A_kJ_mol_SCENARIO_T25C"], 0.55, color=AQUA)
    for xi, v in zip(x, w["delta_A_kJ_mol_SCENARIO_T25C"]):
        ax.text(xi, v + 0.05, f"{100 * v / abs(DG0_NITRIFICATION_KJ_MOL_N):.1f} % of |ΔG°′|", ha="center", fontsize=7.5,
                color=INK2)
    ax.set_xticks(x, labels)
    ax.set_ylabel("ΔA at the pulse (kJ/mol, SCENARIO 25 °C)", fontsize=8.5, color=INK)
    ax.set_title("b  Affinity step on the M6-inverted NH4", loc="left", fontsize=10, color=INK)
    ax = axes[1, 0]
    ax.grid(True, color=GRID, linewidth=0.6, axis="x")
    fa_cols = [c for c in w.columns if c.startswith("fa_median_mgNH3_L_SCENARIO_")]
    for col, color in zip(fa_cols, (BLUE, VIOLET, AQUA, ORANGE)):
        lab = col.removeprefix("fa_median_mgNH3_L_SCENARIO_T").removesuffix("C") + " °C SCENARIO"
        ax.scatter(w[col], w[RESPONSE], color=color, s=40, label=lab)
        for xi, yi, lb in zip(w[col], w[RESPONSE], labels):
            ax.annotate(lb, (xi, yi), textcoords="offset points", xytext=(5, 4), fontsize=7.5, color=INK2)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("Median free ammonia from the M6-inverted TAN (mg NH3/L)", fontsize=8.5, color=INK2)
    ax.set_ylabel("M6-inverted consumption (mg N/L/d)", fontsize=8.5, color=INK)
    ax.set_title(f"c  Free ammonia vs consumption, {EXPLORATORY}", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7.5, frameon=False)
    ax = axes[1, 1]
    eta = w["eta_conversion_fraction"].to_numpy(dtype=float) if "eta_conversion_fraction" in w else np.full(len(w), np.nan)
    ax.bar(x - 0.2, np.nan_to_num(eta), 0.38, color=ORANGE, label="Dose-anchored, M6-inverted C(t)")
    for xi, e in zip(x, eta):
        if np.isnan(e):
            ax.text(xi - 0.2, 0.01, "not\ncomputed", ha="center", va="bottom", fontsize=6.5, color=INK2)
    ax.bar(x + 0.2, w["eta_conversion_fraction_MODEL_BASED"], 0.38, color=VIOLET, label="MODEL-BASED (M6 rate law)")
    top = np.nanmax(np.r_[eta, w["eta_conversion_fraction_MODEL_BASED"].to_numpy(dtype=float)])
    ax.set_ylim(0.0, 1.45 * top)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.set_xticks(x, labels)
    ax.set_ylabel(f"η (fraction of the {thermo.pulse_dose_from_config(cfg).n_mg_L:.2f} mg N/L pulse)", fontsize=8.5,
                  color=INK)
    ax.set_title("d  Conversion before the next pulse", loc="left", fontsize=10, color=INK)
    ratio = r.metrics[r.metrics["metric"] == "free_energy_to_stirrer_energy_ratio"]["status"]
    ratio_txt = "computed" if (ratio == "computed").all() else str(ratio.iloc[0]).removeprefix("not computed: ")
    ax.text(0.98, 0.97, f"P{len(w)}: until the end of the record ({END_D:g} d).\nEnergy ratios: "
            + ("computed" if ratio_txt == "computed" else f"not computed ({ratio_txt})"),
            transform=ax.transAxes, va="top", ha="right", fontsize=7.5, color=INK2)
    fig.tight_layout()
    fig.savefig(png, dpi=150, facecolor=SURFACE)
    fig.savefig(svg, facecolor=SURFACE, metadata={"Date": None})
    plt.close(fig)


def write_report(cfg: Config, r: PulseResult, path: Path) -> None:
    from .report import md_table

    m = r.metrics
    computed = m[m["status"] == "computed"]
    nc = m[m["status"] != "computed"]
    w = wide(m)
    rates = m[m["metric"].str.startswith(("consumption_", "o2_demand", "kla_required")) & (m["status"] == "computed")]
    weighted_rows = m[m["metric"] == RESPONSE].set_index("pulse")
    window_lines = []
    for k, row in weighted_rows.iterrows():
        excl = row["basis"].split("; excludes ")[1] if "; excludes " in str(row["basis"]) else ""
        held = "" if pd.isna(row["held_fraction"]) else f"; held/damped fraction of the windows {row['held_fraction']:.2f}"
        window_lines.append(f"- Pulse {k}: {row['status'] if row['status'] != 'computed' else 'computed'} "
                            f"({row['inputs'].split(';')[0]}){'; excludes ' + excl if excl else ''}{held}.")
    stirrer = m[m["metric"] == "specific_stirring_energy_J_per_gN"]["status"]
    rho = r.associations.set_index("descriptor")["spearman_rho"]
    fa_rho = ", ".join(f"{d.removeprefix('fa_median_mgNH3_L_')}: ρ = {v:+.2f}" for d, v in rho.items()
                       if d.startswith("fa_median") and pd.notna(v))
    lines = [
        "# Pulse-level efficiency and energy descriptors (M7)",
        "",
        "Generated by `python -m bnr_neq pulses`. Do not edit by hand. Implementer output, not approved.",
        "",
        (f"`outputs/tables/pulses.csv` has one row per pulse and metric: {len(computed)} computed (inputs listed), "
        f"{len(nc)} explicitly not computed. **Single NH4 observation model (W1):** rates, free ammonia, η, ΔA and "
        "oxygen demand are derived from one C(t), the inversion of the M6 log observation model "
        "(`nh4_observation_model` = M6_log, `nh4_series` = the C(t) identifier, checked before writing). GAIN-CORRECTED "
        "values (linear electrode gain model) and values that take the reading as TAN are labelled SENSITIVITY; the "
        "uncorrected electrode slopes are measured quantities, not NH4 truths. MODEL-BASED values come from the M6 "
        "rate law (same observation model). SCENARIO values use assumed temperatures or evaporation rates. "
        "`param_band_low/high` is the envelope over the 95 % profile vectors; C0 and the intercept a are not "
        "identifiable, so the band is open on the side named in the basis."),
        "",
        ("Rates carry the 95 % intervals of `segments.csv` (`ci95_low`, `ci95_high`; duration-weighted rates use "
        "duration-weighted interval bounds, which treats window errors as fully correlated) and the fraction of their "
        "slope windows covered by held/damped NH4 output (`held_fraction`). A slope window whose M6-inverted "
        "consumption is not positive is not computed in `nh4_truth_windows.csv` and is excluded from the "
        "duration-weighted rates:"),
        "",
        *window_lines,
        "",
        "## Pulse table (computed values)",
        "",
        md_table(w),
        "",
        "### Rates with 95 % intervals and held fraction",
        "",
        md_table(rates[["pulse", "metric", "value", "ci95_low", "ci95_high", "held_fraction", "unit"]]),
        "",
        "## Not computed",
        "",
        md_table(nc[["pulse", "metric", "status", "inputs"]].drop_duplicates(subset=["metric", "status"])),
        "",
        (f"Free energy of one fully converted pulse at standard conditions (SCENARIO 25 °C): "
        f"{free_energy_per_pulse_J(cfg):.1f} J in the {cfg.get('vessel.volume_L') * 1e3:g} mL vessel. The ratio to the "
        "stirrer's electrical energy and the specific stirring energy: "
        + ("computed." if (stirrer == "computed").all() else f"{stirrer.iloc[0]}.")),
        "",
        "## Associations (EXPLORATORY: four pulses)",
        "",
        md_table(r.associations),
        "",
        (f"No association is claimed: with {int(r.associations['n'].max())} pulses every rank statistic rests on a "
        f"handful of points (free ammonia {fa_rho}). A positive ρ would be the opposite of the free-ammonia "
        "inhibition hypothesis."),
        "",
        ("`removal_to_dose_ratio` is (C before the pulse + dose − C before the next pulse)/dose on the M6-inverted C(t) "
         "(S3). It exceeds 1 when carry-over TAN is removed as well; the conversion fraction η of the pulse is reported "
         "only for 0 ≤ ratio ≤ 1 and is otherwise `not computed` with the reason. MODEL-BASED η comes from the rate law "
         "of the M6 selected model and carries no uncertainty interval."),
        "",
        "## Oxygen and energy metric family (W4)",
        "",
        ("Per pulse: N removed, O2 required (4.57 g O2/g N; 3.43 as a nitritation case), O2 transfer capacity "
         "kLa (C* − DO) V Δt per kLa_O2 SCENARIO (DO = 0 while not measured), utilization = required/capacity "
         "(a scenario whose capacity is below the requirement is `not computed: SCENARIO infeasible`), aeration-equivalent "
         "energy = O2/AE per field-AE SCENARIO (demand-based floor and supply-based value, in kWh/kg N), the free "
         "energy released over that energy, the W3 idle fraction as bounds (LOWER_BOUND: inactive bins over all bins; "
         "UPPER_BOUND: inactive, unclassified and artifact bins over all bins) with its coverage (`classified_fraction`) "
         "and, only where the coverage is at least one half, the CONDITIONAL value over classified time, and the "
         "dose-removal time (dose / post-pulse rate, not a fitted relaxation constant) against the pulse interval. "
         "The supply-based energy divides a transfer capacity at DO = 0 by a field AE defined at DO = 2 mg/L. The "
         "stirrer terms are LABORATORY-SPECIFIC."),
        "",
        ("R8 CASE: a slope window between two electrode steps shorter than 3 x the longest delayed electrode response is "
         "an electrode fragment, also after a pulse (pass 3); `consumption_weighted_mgN_L_d_R8_CASE_incl_fragments` "
         "repeats the duration-weighted consumption with the fragment windows of the pulse included, beside the primary "
         "rate without them. `eta_gap_inversion_misfit_MODEL_BASED` decomposes the difference between the "
         "dose-anchored removal per dose and the MODEL-BASED η."),
        "",
        "## The inverted C(t) against the W2 bound (S2)",
        "",
        ("Every row derived from the M6 observation model carries `w2_bound_check`: `QUALIFIED: ...` when a checkpoint "
         "bounding the pulse needs more initial alkalinity than the 1-15 meq/L scan (at the best fit), else the "
         "consistent requirement. The stage refuses to write a C(t)-derived row of a contradicted pulse without it."),
        "",
        md_table(r.consistency, digits=3) if r.consistency is not None else "not computed",
        "",
        "## Water-loss constraints (W8)",
        "",
        md_table(r.volume_constraints, digits=3) if r.volume_constraints is not None else "not computed",
        "",
        "## Dose-anchored TAN bound against the inverted C(t) (W1 with the W2 alkalinity bound)",
        "",
        ("C_up(t) = (C(t0) V(t0) + Σ m_N)/V(t) assumes no removal; C_low(t; Alk0) subtracts the W2 upper bound on "
        "cumulative nitrification (air-equilibrium CO2, SCENARIO 25 °C). `alk0_required_*` is the initial alkalinity "
        "that the removal implied by C(t) would need (2 eq/mol N for nitrification, 1 eq/mol N for any other "
        "proton-releasing sink such as NH3 volatilization)."),
        "",
        md_table(r.truth_balance) if r.truth_balance is not None else "not computed",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir) -> PulseResult:
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    inputs = load_inputs(out_dir)
    d = out_dirs(out_dir)
    r = compute(cfg, inputs)
    r.metrics.to_csv(d["tables"] / "pulses.csv", index=False)
    r.associations.to_csv(d["tables"] / "pulses_associations.csv", index=False)
    r.truth_windows.to_csv(d["tables"] / "nh4_truth_windows.csv", index=False)
    r.truth_balance.to_csv(d["tables"] / "nh4_truth_balance.csv", index=False)
    r.consistency.to_csv(d["tables"] / "nh4_truth_w2_consistency.csv", index=False)
    r.volume_constraints.to_csv(d["tables"] / "volume_constraints.csv", index=False)
    figure(cfg, r, d["figures"] / "fig7_pulses.png", d["figures"] / "fig7_pulses.svg")
    write_report(cfg, r, d["reports"] / "pulses.md")
    return r
