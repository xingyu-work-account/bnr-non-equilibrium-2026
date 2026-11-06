"""Phase 0 sensor validation with the four known NH4Cl pulses (M2)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import figures, qc, thermo
from .config import Config, MissingInputError
from .constants import KELVIN_OFFSET
from .pipeline import OUTPUT_DIR, out_dirs
from .report import md_table

# Windows relative to each pulse, in days: (pre_from, pre_to), (post_from, post_to); [from, to)
NH4_EC_WINDOWS = ((-0.1, 0.0), (0.4, 0.6))
ORP_PH_WINDOWS = ((-3 / 24, -1 / 24), (1 / 24, 3 / 24))

EC_START_WINDOW = (0.3, 1.3)
EC_END_WINDOW = (48.5, 49.6)
NH4_START_WINDOW = (0.3, 1.3)
NH4_END_WINDOW = (46.9, 47.9)
CLOCK_EVENT_WINDOW = (2.75, 3.00)  # multi-probe event near day 2.88 (plan, Phase 1)
CLOCK_STEP_WINDOW_D = qc.JUMP_WINDOW_D  # sustained steps only: a median window ignores excursions < w/2

# A priori verdict criteria (fixed before computing the observations).
NH4_PASS_FRACTION = (0.8, 1.2)  # observed / predicted step
EC_PASS_MIN_FRACTION = 0.5  # of the limiting prediction (allows for activity effects)
ORP_PASS_MIN_FRACTION = 0.5  # of the smallest chloride-scenario prediction
PH_MAX_ABS_CHANGE = 0.05  # "a few hundredths" (plan, Phase 0)

# Preliminary-analysis reference values (implementation_prompt.md, M2 item 6).
REFERENCES = [
    # (quantity, pulse or None, reference, tolerance, unit)
    *[("nh4_response_fraction_pct", k + 1, v, 3.0, "%") for k, v in enumerate([37, 29, 26, 17])],
    *[("nh4_obs_change_mgNH4_L", k + 1, v, 2.0, "mg NH4+/L") for k, v in enumerate([31.1, 24.4, 22.1, 14.7])],
    *[("ec_obs_change_mS_cm", k + 1, v, 0.01, "mS/cm") for k, v in enumerate([0.005, -0.062, -0.007, 0.016])],
    *[("orp_obs_change_mV", k + 1, v, 0.5, "mV") for k, v in enumerate([0.4, 0.0, 0.1, 0.2])],
    *[("ph_obs_change", k + 1, v, 0.02, "pH") for k, v in enumerate([-0.01, 0.02, 0.00, 0.04])],
    ("ec_median_day_0.3_1.3_mS_cm", None, 1.069, 0.005, "mS/cm"),
    ("ec_median_day_48.5_49.6_mS_cm", None, 1.071, 0.005, "mS/cm"),
    ("nh4_median_start_mgNH4_L", None, 47.1, 0.5, "mg NH4+/L"),
    ("nh4_median_end_mgNH4_L", None, 43.1, 0.5, "mg NH4+/L"),
    *[("nh4_offset_jump_day", None, v, 0.02, "d") for v in (2.59, 17.84, 37.19, 44.63)],
    ("nh4_delayed_response_day_pulse1", 1, 8.85, 0.02, "d"),
    ("nh4_delayed_response_day_pulse2", 2, 10.88, 0.02, "d"),
    ("ec_late_pulse_count", None, 29, 3, "count"),
]


@dataclass
class ValidationResult:
    pulse_table: pd.DataFrame
    windows: pd.DataFrame
    balances: pd.DataFrame
    jumps: pd.DataFrame
    held: pd.DataFrame
    hf_subsegments: pd.DataFrame
    ec_pulses: pd.DataFrame
    ec_sensitivity: pd.DataFrame
    clock: pd.DataFrame
    references: pd.DataFrame
    verdicts: dict[str, str] = field(default_factory=dict)


def window_median(t: np.ndarray, x: np.ndarray, lo: float, hi: float, label: str) -> tuple[float, int]:
    """Median of x over [lo, hi); raise if the window holds no sample."""
    m = (t >= lo) & (t < hi)
    n = int(m.sum())
    if n == 0:
        raise ValueError(f"{label}: no samples in window [{lo}, {hi}) d")
    return float(np.median(x[m])), n


def window_change(t, x, t0, windows, label) -> tuple[float, list[dict]]:
    (a0, a1), (b0, b1) = windows
    pre, n_pre = window_median(t, x, t0 + a0, t0 + a1, f"{label} pre")
    post, n_post = window_median(t, x, t0 + b0, t0 + b1, f"{label} post")
    rows = [
        {"window": "pre", "t_from_d": t0 + a0, "t_to_d": t0 + a1, "n": n_pre, "median": pre},
        {"window": "post", "t_from_d": t0 + b0, "t_to_d": t0 + b1, "n": n_post, "median": post},
    ]
    return post - pre, rows


def scenario_label(cl_mM: float, T_C: float) -> str:
    return f"SCENARIO_Cl0_{cl_mM:g}mM_T{T_C:g}C"


def pulse_table(cfg: Config, natives: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    dose = thermo.pulse_dose_from_config(cfg)
    ec_imm = thermo.ec_step_immediate_mS_cm(dose.c_mmol_L)
    ec_nit = thermo.ec_step_nitrified_mS_cm(dose.c_mmol_L)
    n_factor = cfg.get("sensors.nh4_to_N_factor")
    series = {
        "nh4": (natives["nh4"]["t_d"].to_numpy(), natives["nh4"]["nh4_mgNH4_L"].to_numpy(), NH4_EC_WINDOWS),
        "ec": (natives["cond"]["t_d"].to_numpy(), natives["cond"]["EC_mS_cm"].to_numpy(), NH4_EC_WINDOWS),
        "orp": (natives["orp"]["t_d"].to_numpy(), natives["orp"]["orp_mV_vs_Ag"].to_numpy(), ORP_PH_WINDOWS),
        "ph": (natives["ph"]["t_d"].to_numpy(), natives["ph"]["pH"].to_numpy(), ORP_PH_WINDOWS),
    }
    try:
        cfg.measurement("alkalinity_mmol_L")
        ph_pred = "not computed: pH prediction not implemented"
    except MissingInputError as exc:
        ph_pred = str(exc)
    rows, win_rows = [], []
    for k, tp in enumerate(cfg.pulse_times_d, start=1):
        obs = {}
        for name, (t, x, wins) in series.items():
            obs[name], w = window_change(t, x, tp, wins, f"pulse {k} {name}")
            win_rows += [{"pulse": k, "sensor": name, **r} for r in w]
        row = {
            "pulse": k,
            "t_pulse_d": tp,
            "dose_mmol_L": dose.c_mmol_L,
            "nh4_pred_change_mgNH4_L": dose.nh4_mg_L,
            "nh4_obs_change_mgNH4_L": obs["nh4"],
            "nh4_response_fraction": obs["nh4"] / dose.nh4_mg_L,
            "nh4_pred_change_mgN_L": dose.n_mg_L,
            "nh4_obs_change_mgN_L": obs["nh4"] * n_factor,
            "ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C": ec_imm,
            "ec_pred_change_nitrified_mS_cm_limiting_SCENARIO_T25C": ec_nit,
            "ec_obs_change_mS_cm": obs["ec"],
            "ec_response_fraction_of_immediate": obs["ec"] / ec_imm,
            "orp_obs_change_mV": obs["orp"],
        }
        for T_C in cfg.get("scenarios.temperature_C"):
            for c0 in cfg.get("scenarios.baseline_chloride_mM"):
                c_before = c0 + (k - 1) * dose.c_mmol_L  # chloride is conservative
                row[f"orp_pred_change_mV_{scenario_label(c0, T_C)}"] = thermo.chloride_orp_rise_mV(
                    c_before, dose.c_mmol_L, T_C + KELVIN_OFFSET
                )
        row["ph_obs_change"] = obs["ph"]
        row["ph_pred_change"] = ph_pred
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(win_rows)


def balances(cfg: Config, natives: dict[str, pd.DataFrame]) -> pd.DataFrame:
    dose = thermo.pulse_dose_from_config(cfg)
    n_pulses = len(cfg.pulse_times_d)
    n_factor = cfg.get("sensors.nh4_to_N_factor")
    tc, ec = natives["cond"]["t_d"].to_numpy(), natives["cond"]["EC_mS_cm"].to_numpy()
    tn, nh4 = natives["nh4"]["t_d"].to_numpy(), natives["nh4"]["nh4_mgNH4_L"].to_numpy()
    ec0, _ = window_median(tc, ec, *EC_START_WINDOW, "EC start")
    ec1, _ = window_median(tc, ec, *EC_END_WINDOW, "EC end")
    n0, _ = window_median(tn, nh4, *NH4_START_WINDOW, "NH4 start")
    n1, _ = window_median(tn, nh4, *NH4_END_WINDOW, "NH4 end")
    added = n_pulses * dose.n_mg_L
    removed = n0 * n_factor + added - n1 * n_factor
    span = np.mean(NH4_END_WINDOW) - np.mean(NH4_START_WINDOW)
    imm, nit = thermo.ec_step_immediate_mS_cm(dose.c_mmol_L), thermo.ec_step_nitrified_mS_cm(dose.c_mmol_L)
    rows = [
        ("ec_median_day_0.3_1.3_mS_cm", ec0, "mS/cm", "measured", "median EC, day 0.3 to 1.3"),
        ("ec_median_day_48.5_49.6_mS_cm", ec1, "mS/cm", "measured", "median EC, day 48.5 to 49.6 (Cond only, outside common window)"),
        ("ec_pred_end_no_nitrification_mS_cm_limiting_SCENARIO_T25C", ec0 + n_pulses * imm, "mS/cm", "predicted",
         f"EC start + {n_pulses} x {imm:.3f} (limiting values, SCENARIO 25 °C)"),
        ("ec_pred_end_complete_nitrification_mS_cm_limiting_SCENARIO_T25C", ec0 + n_pulses * nit, "mS/cm", "predicted",
         f"EC start + {n_pulses} x {nit:.3f} (limiting values, SCENARIO 25 °C)"),
        ("nh4_median_start_mgNH4_L", n0, "mg NH4+/L", "measured", "median NH4, day 0.3 to 1.3"),
        ("nh4_median_end_mgNH4_L", n1, "mg NH4+/L", "measured", "median NH4, day 46.9 to 47.9"),
        ("nh4_median_start_mgN_L", n0 * n_factor, "mg N/L", "measured", "start median x N factor"),
        ("nh4_median_end_mgN_L", n1 * n_factor, "mg N/L", "measured", "end median x N factor"),
        ("n_added_mgN_L", added, "mg N/L", "known dose", f"{n_pulses} x {dose.n_mg_L:.2f}"),
        ("n_removed_from_nh4_pool_mgN_L", removed, "mg N/L", "derived", "start + added - end (electrode values, uncorrected)"),
        ("balance_interval_d", span, "d", "derived", "midpoint of end window - midpoint of start window"),
        ("implied_n_removal_rate_mgN_L_d", removed / span, "mg N/L/d", "derived", "removed / interval"),
        ("chloride_added_mmol_L", n_pulses * dose.c_mmol_L, "mmol/L", "known dose", "conservative; staircase of equal steps"),
    ]
    return pd.DataFrame(rows, columns=["quantity", "value", "unit", "basis", "definition"])


def hf_subsegment_table(t: np.ndarray, resid: np.ndarray) -> pd.DataFrame:
    """Robust HF residual per preliminary-analysis sub-segment (plan item D3)."""
    bounds = [(0.3, 2.59), (2.59, 8.82), (8.82, 10.85), (10.85, 17.84), (17.84, 24.60),
              (24.60, 32.36), (32.36, 37.19), (37.19, 44.63), (44.63, 47.92)]
    rows = []
    for a, b in bounds:
        m = (t >= a) & (t < b)
        rows.append({"from_d": a, "to_d": b, "sigma_hf_mgNH4_L": qc.robust_sigma(resid[m]),
                     "rms_hf_mgNH4_L": float(np.sqrt(np.nanmean(resid[m] ** 2)))})
    return pd.DataFrame(rows)


def clock_check(natives: dict[str, pd.DataFrame]) -> pd.DataFrame:
    cols = {"cond": "EC_mS_cm", "orp": "orp_mV_vs_Ag", "ph": "pH"}
    rows = []
    for key, col in cols.items():
        t, x = natives[key]["t_d"].to_numpy(), natives[key][col].to_numpy()
        ts, step = qc.step_time_in_window(t, x, *CLOCK_EVENT_WINDOW, w_d=CLOCK_STEP_WINDOW_D)
        rows.append({"sensor": key, "variable": col, "t_step_d": ts, "step": step})
    df = pd.DataFrame(rows)
    # ORP has the sharpest step at this event; EC shows only a slow drift there
    df["offset_vs_orp_d"] = df["t_step_d"] - df.loc[df["sensor"] == "orp", "t_step_d"].iloc[0]
    return df


def reference_comparison(pt: pd.DataFrame, bal: pd.DataFrame, jumps: pd.DataFrame,
                         ec_pulses: pd.DataFrame) -> pd.DataFrame:
    bal_v = dict(zip(bal["quantity"], bal["value"]))
    jumps_ok = jumps[jumps["jump"]]
    offsets = jumps_ok[jumps_ok["class"] == "electrode_offset"]
    rows = []
    for qty, pulse, ref, tol, unit in REFERENCES:
        note = ""
        if qty == "nh4_response_fraction_pct":
            val = 100 * pt.loc[pt["pulse"] == pulse, "nh4_response_fraction"].iloc[0]
        elif qty in ("nh4_obs_change_mgNH4_L", "ec_obs_change_mS_cm", "orp_obs_change_mV", "ph_obs_change"):
            val = pt.loc[pt["pulse"] == pulse, qty].iloc[0]
        elif qty in bal_v:
            val = bal_v[qty]
        elif qty == "nh4_offset_jump_day":
            k = int(np.argmin(np.abs(offsets["t_jump_d"] - ref)))
            val = offsets["t_jump_d"].iloc[k]
            note = f"nearest detected offset; detector run onset {offsets['run_start_d'].iloc[k]:.3f} d"
        elif qty.startswith("nh4_delayed_response_day"):
            resp = jumps_ok[(jumps_ok["class"] == "pulse_response") & (jumps_ok["pulse"] == pulse)]
            k = int(np.argmax(resp["level_change"].to_numpy()))
            val = resp["t_jump_d"].iloc[k]
            note = f"largest response to pulse {pulse}; lag {resp['lag_h'].iloc[k]:.2f} h"
        elif qty == "ec_late_pulse_count":
            val = float(len(ec_pulses))
            note = f"median decay {ec_pulses['tau_min'].median():.1f} min (reference about 24 min)"
        else:
            raise KeyError(qty)
        rows.append({"quantity": qty, "pulse": pulse, "reference": ref, "obtained": val, "tolerance": tol,
                     "unit": unit, "abs_diff": abs(val - ref), "pass": abs(val - ref) <= tol + 1e-9, "note": note})
    out = pd.DataFrame(rows)
    out["pulse"] = pd.array(out["pulse"], dtype="Int64")
    return out


def verdicts(pt: pd.DataFrame) -> dict[str, str]:
    lo, hi = NH4_PASS_FRACTION
    f = pt["nh4_response_fraction"]
    orp_pred_cols = [c for c in pt.columns if c.startswith("orp_pred_change_mV_")]
    orp_min_pred = pt[orp_pred_cols].min(axis=1)
    return {
        "nh4": "PASS" if ((f >= lo) & (f <= hi)).all() else
               f"FAIL: under-response, {100 * f.min():.0f} to {100 * f.max():.0f}% of the dose",
        "ec": "PASS" if (pt["ec_response_fraction_of_immediate"] >= EC_PASS_MIN_FRACTION).all() else
              "FAIL: no response to the pulses",
        "orp": "PASS" if (pt["orp_obs_change_mV"] >= ORP_PASS_MIN_FRACTION * orp_min_pred).all() else
               "FAIL: no chloride response at the bare Ag",
        "ph": "PASS (no contradiction)" if (pt["ph_obs_change"].abs() <= PH_MAX_ABS_CHANGE).all() else
              "FAIL: change larger than expected",
    }


def verdict_table(v: dict[str, str], pt: pd.DataFrame) -> pd.DataFrame:
    """The verdicts with their a priori criteria, as written to ``sensor_verdicts.csv``."""
    return pd.DataFrame([
        ("nh4", "NH4+ electrode", v["nh4"],
         f"observed/predicted step within {NH4_PASS_FRACTION[0]} to {NH4_PASS_FRACTION[1]} for every pulse"),
        ("ec", "Conductivity", v["ec"], f"observed step ≥ {EC_PASS_MIN_FRACTION} x limiting Kohlrausch prediction"),
        ("orp", "ORP (Pt vs bare Ag)", v["orp"],
         f"observed rise ≥ {ORP_PASS_MIN_FRACTION} x smallest chloride-scenario prediction"),
        ("ph", "pH", v["ph"], f"abs(ΔpH) ≤ {PH_MAX_ABS_CHANGE}; quantitative prediction {pt['ph_pred_change'].iloc[0]}"),
    ], columns=["key", "sensor", "verdict", "criterion"])


# --- Stage --------------------------------------------------------------------
def compute(cfg: Config, natives: dict[str, pd.DataFrame]) -> ValidationResult:
    pt, windows = pulse_table(cfg, natives)
    bal = balances(cfg, natives)
    tn, xn = natives["nh4"]["t_d"].to_numpy(), natives["nh4"]["nh4_mgNH4_L"].to_numpy()
    jumps = qc.jump_table(tn, xn, cfg.pulse_times_d)
    resid = qc.hf_residual(tn, xn)
    held = qc.held_segments(qc.hf_blocks(tn, resid, 0.0, 47.92))
    tc, ec = natives["cond"]["t_d"].to_numpy(), natives["cond"]["EC_mS_cm"].to_numpy()
    rc = qc.ec_residual(tc, ec)
    thr = qc.ec_pulse_threshold(tc, rc)
    ec_pulses = qc.detect_ec_pulses(tc, ec, rc, thr)
    sens = qc.ec_pulse_sensitivity(tc, rc, [thr, *qc.EC_SENSITIVITY_THRESHOLDS])
    sens["rule"] = ["noise floor (used)"] + ["sensitivity"] * len(qc.EC_SENSITIVITY_THRESHOLDS)
    refs = reference_comparison(pt, bal, jumps, ec_pulses)
    return ValidationResult(pt, windows, bal, jumps, held, hf_subsegment_table(tn, resid), ec_pulses, sens,
                            clock_check(natives), refs, verdicts(pt))


def load_natives(out_dir: Path) -> dict[str, pd.DataFrame]:
    natives = {}
    for key in ("cond", "nh4", "orp", "ph"):
        p = Path(out_dir) / "processed" / f"native_{key}.parquet"
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run `python -m bnr_neq ingest` first")
        natives[key] = pd.read_parquet(p)
    return natives


def write_report(cfg: Config, r: ValidationResult, path: Path) -> None:
    pt, refs = r.pulse_table, r.references
    dose = thermo.pulse_dose_from_config(cfg)
    scen_cols = [c for c in pt.columns if c.startswith("orp_pred_change_mV_")]
    pulse_view = pt[["pulse", "t_pulse_d", "nh4_pred_change_mgNH4_L", "nh4_obs_change_mgNH4_L",
                     "nh4_response_fraction", "ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C",
                     "ec_pred_change_nitrified_mS_cm_limiting_SCENARIO_T25C", "ec_obs_change_mS_cm",
                     "orp_obs_change_mV", "ph_obs_change"]].copy()
    pulse_view["orp_pred_change_mV_SCENARIO_range"] = [
        f"{a:.1f} to {b:.1f}" for a, b in zip(pt[scen_cols].min(axis=1), pt[scen_cols].max(axis=1))]
    orp_scen = pt[["pulse", *scen_cols]].copy()
    orp_scen.columns = ["pulse", *[c.replace("orp_pred_change_mV_SCENARIO_", "") for c in scen_cols]]
    j = r.jumps
    jump_view = j[["run_start_d", "t_jump_d", "level_change", "transition_10_90_min", "sustained", "class",
                   "pulse", "lag_h", "jump"]]
    v = r.verdicts
    n_fail = int((~refs["pass"]).sum())
    sens_txt = ", ".join(f"{h:.3f}: {n}" for h, n in zip(r.ec_sensitivity["threshold_mS_cm"],
                                                          r.ec_sensitivity["n_pulses"]))
    lines = [
        "# Sensor validation with the known NH4Cl pulses (M2)",
        "",
        "Generated by `python -m bnr_neq validate`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Verdicts",
        "",
        "Criteria were fixed before computing the observations (`outputs/tables/sensor_verdicts.csv`).",
        "",
        md_table(verdict_table(v, pt).drop(columns=["key"])),
        "",
        "**P0 findings.**",
        "",
        (f"1. **NH4+ under-response (S1).** The electrode registered {100 * pt['nh4_response_fraction'].min():.0f} to "
        f"{100 * pt['nh4_response_fraction'].max():.0f}% of the {dose.nh4_mg_L:.2f} mg NH4+/L added at each pulse, "
        "declining from pulse 1 to pulse 4. Absolute NH4+ levels, rates and conversion are NO-GO until recalibration "
        "and laboratory NH4-N."),
        (f"2. **No conductivity response (S3).** Each pulse should raise EC by {pt['ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C'].iloc[0]:.3f} "
        f"mS/cm (limiting value); observed changes are {pt['ec_obs_change_mS_cm'].min():+.3f} to "
        f"{pt['ec_obs_change_mS_cm'].max():+.3f} mS/cm, and the record ends where it started "
        f"({r.balances.set_index('quantity').loc['ec_median_day_0.3_1.3_mS_cm', 'value']:.3f} vs "
        f"{r.balances.set_index('quantity').loc['ec_median_day_48.5_49.6_mS_cm', 'value']:.3f} mS/cm) instead of rising to "
        f"at least {r.balances.set_index('quantity').loc['ec_pred_end_complete_nitrification_mS_cm_limiting_SCENARIO_T25C', 'value']:.2f} mS/cm. "
        "The record does not represent bulk solution conductivity."),
        (f"3. **No ORP chloride response (S4).** A chloride-controlled Ag surface predicts a rise of "
        f"{pt[scen_cols].min().min():.1f} to {pt[scen_cols].max().max():.1f} mV per pulse across the chloride and temperature "
        f"scenarios; observed changes are {pt['orp_obs_change_mV'].min():+.2f} to {pt['orp_obs_change_mV'].max():+.2f} mV. "
        "Absolute E_H is NO-GO; ORP remains usable as a relative signal only."),
        "",
        f"{thermo.KOHLRAUSCH_NOTE}",
        "",
        "## Pulse table",
        "",
        ("Windows: NH4+ and EC, median 0.4 to 0.6 d after the pulse minus median over the 0.1 d before; ORP and pH, "
        "median 1 to 3 h after minus median 3 to 1 h before. Native (reconstructed) series. Sample counts per window "
        "are in `outputs/tables/sensor_validation_windows.csv`."),
        "",
        md_table(pulse_view),
        "",
        ("ORP predictions under each SCENARIO (baseline chloride Cl0 and temperature T are scenarios, not measurements). "
        "The chloride before pulse k is Cl0 + (k − 1) x "
        f"{dose.c_mmol_L:.4f} mmol/L. Values in mV:"),
        "",
        md_table(orp_scen, digits=3),
        "",
        "## Cumulative balances",
        "",
        md_table(r.balances, digits=5),
        "",
        "## NH4+ jump detection",
        "",
        (f"Detector: step(t) = median over (t, t + {qc.JUMP_WINDOW_D} d] − median over (t − {qc.JUMP_WINDOW_D} d, t] on the "
        f"reconstructed native series; candidates where abs(step) ≥ {qc.JUMP_THRESHOLD_MG_L} mg NH4+/L, merged within "
        f"{qc.JUMP_MERGE_D} d. Each candidate is located at the 50 % crossing (`t_jump_d`) of a {qc.TRANSITION_SMOOTH} rolling "
        f"median between pre and post levels taken {qc.LEVEL_OFFSET_D} to {qc.LEVEL_OFFSET_D + qc.LEVEL_WINDOW_D} d away. "
        f"A candidate is a jump if sustained (abs(post − pre) ≥ {qc.JUMP_THRESHOLD_MG_L}) and located. A rising jump "
        f"{qc.PULSE_LAG_MIN_D} to {qc.PULSE_LAG_MAX_D} d after a pulse is a pulse response; every other jump is an electrode "
        "offset. Transition times shorter than the 14.4-min stamp rely on the uniform-within-stamp assumption; a "
        "transition of one sample interval (about 0.33 min at 20 s) means the jump falls between two consecutive rows."),
        "",
        md_table(jump_view),
        "",
        "## Held or damped NH4+ output",
        "",
        (f"HF residual = x − centred {qc.HF_MEDIAN_WINDOW} rolling median; σ_HF = 1.4826 x MAD per {qc.HF_BLOCK_D} d block. "
        f"Held/damped where σ_HF < {qc.HELD_SIGMA_MG_L} mg NH4+/L (the output resolution) for at least {qc.HELD_MIN_D} d."),
        "",
        md_table(r.held) if len(r.held) else "No held segment found.",
        "",
        "σ_HF and RMS of the HF residual per preliminary-analysis sub-segment:",
        "",
        md_table(r.hf_subsegments),
        "",
        "## Late conductivity pulses (after day 41.9)",
        "",
        (f"Residual from a centred {qc.EC_BASELINE_WINDOW} rolling median; peaks with height and prominence ≥ threshold, "
        f"≥ {qc.EC_MIN_SEPARATION_MIN:g} min apart; threshold = {qc.EC_THRESHOLD_QUANTILE:.1%} quantile of the residual before "
        f"day {qc.EC_LATE_START_D} = {r.ec_sensitivity['threshold_mS_cm'].iloc[0]:.4f} mS/cm (rule fixed before counting). "
        f"Found {len(r.ec_pulses)} pulses, median amplitude {r.ec_pulses['amplitude_mS_cm'].median():.3f} mS/cm, median "
        f"exponential decay {r.ec_pulses['tau_min'].median():.1f} min ({int(r.ec_pulses['fit_ok'].sum())} of "
        f"{len(r.ec_pulses)} fits converged). Count vs threshold (mS/cm): {sens_txt}. Treated as probe events. Limitation: a quantile threshold admits about 0.1 % of pure-noise samples by construction; the smallest detected amplitude is reported in `outputs/tables/ec_late_pulses.csv` so its margin over the threshold can be checked."),
        "",
        "## Clock check at the multi-probe event near day 2.88 (plan item D7)",
        "",
        md_table(r.clock),
        "",
        (f"Sustained steps located with a {CLOCK_STEP_WINDOW_D} d median window in {CLOCK_EVENT_WINDOW[0]} to "
         f"{CLOCK_EVENT_WINDOW[1]} d. Largest offset between sensors: "
         f"pH vs ORP {r.clock.set_index('sensor').loc['ph', 'offset_vs_orp_d']:+.4f} d; EC has no distinct step at "
         "this event (its largest median difference falls on a slow drift). The step locator cannot resolve offsets below about "
         "half its window, and the sensors respond to the event with different shapes (pH also shows a transient "
         "excursion at about day 2.81 to 2.835), so this check does not establish clock alignment; it only bounds a "
         "gross offset. Not corrected."),
        "",
        "## Reference comparison",
        "",
        (f"{len(refs) - n_fail} of {len(refs)} reference values reproduced within tolerance; {n_fail} not reproduced "
        "(explained in `docs/progress.md`)."),
        "",
        md_table(refs),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir: Path | str = OUTPUT_DIR) -> ValidationResult:
    d = out_dirs(out_dir)
    natives = load_natives(Path(out_dir))
    r = compute(cfg, natives)
    r.pulse_table.to_csv(d["tables"] / "sensor_validation.csv", index=False)
    r.windows.to_csv(d["tables"] / "sensor_validation_windows.csv", index=False)
    r.balances.to_csv(d["tables"] / "sensor_balances.csv", index=False)
    r.jumps.to_csv(d["tables"] / "nh4_jumps.csv", index=False)
    r.held.to_csv(d["tables"] / "nh4_held_segments.csv", index=False)
    r.hf_subsegments.to_csv(d["tables"] / "nh4_hf_residual_subsegments.csv", index=False)
    r.ec_pulses.to_csv(d["tables"] / "ec_late_pulses.csv", index=False)
    r.ec_sensitivity.to_csv(d["tables"] / "ec_late_pulse_sensitivity.csv", index=False)
    r.clock.to_csv(d["tables"] / "clock_check.csv", index=False)
    r.references.to_csv(d["tables"] / "sensor_validation_references.csv", index=False)
    verdict_table(r.verdicts, r.pulse_table).to_csv(d["tables"] / "sensor_verdicts.csv", index=False)
    figures.validation_figures(cfg, out_dir)
    write_report(cfg, r, d["reports"] / "sensor_validation.md")
    return r
