"""Phase 1: segmentation at the pulses, consumption slopes, kinetics, change points, periodicity (M3)."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import ruptures as rpt
from scipy.optimize import least_squares
from scipy.signal import lombscargle
from scipy.special import wrightomega
from scipy.stats import theilslopes
from statsmodels.stats.diagnostic import acorr_ljungbox

from . import qc, thermo, volume
from .config import Config
from .constants import ANTHONISEN_EXP_K, KELVIN_OFFSET, NH3_PER_N

SEGMENT_START_D = 1.0  # start-up transient (day 0 to 1) excluded, plan item D10
COMMON_END_D = 47.92
EXCLUDE_AFTER_BOUNDARY_D = 0.1
MIN_FIT_SPAN_D = 0.25  # a slope needs at least 25 grid bins after the exclusion
MIN_FIT_BINS = 25  # grid bins (0.01 d) lying wholly inside the fit window
THEIL_SEN_ALPHA = 0.95
# R8: a baseline window bounded on both sides by electrode steps (offset jumps or step candidates, not pulses or pulse
# responses) and shorter than FRAGMENT_DELAYS x the longest delayed electrode response to a pulse (M2 lag of the
# pulse-response jumps) cannot separate electrode relaxation from consumption; its slope is kept only as a diagnostic.
# 3 delays: the usual settling span of a first-order transient (95 %). Fixed before inspecting any fragment slope.
FRAGMENT_DELAYS = 3.0
ELECTRODE_BOUNDARY_PREFIXES = ("jump (electrode_offset)", "candidate (")

# Preliminary-analysis sub-segments (implementation_prompt.md, M3 item 6), used diagnostically.
REFERENCE_WINDOWS = [(0.3, 2.59), (2.59, 8.82), (10.85, 17.84), (17.84, 24.60),
                     (24.60, 32.36), (32.36, 37.19), (37.19, 44.63), (44.63, 47.92)]
REFERENCE_SLOPES_MGNH4_L_D = [-1.89, -3.42, -1.95, -0.57, -0.48, -0.31, -0.28, 0.12]

# Kinetic fits
LOG_BOUNDS = (np.log(1e-6), np.log(1e4))
REL_SE_MAX = 1.0  # relative standard error above which a parameter is not identifiable
LJUNG_BOX_LAG = 10

# Change points (PELT)
PELT_MODEL = "linear"  # discontinuous piecewise-linear regression on [t, 1]
PELT_MIN_SIZE = 25  # bins (0.25 d)
PELT_JUMP = 5
PELT_NOISE_WINDOW = 101  # bins (1 d) rolling median; residual about it = noise
PELT_BETA = 3.0  # BIC: 3 new parameters (slope, intercept, location) per break
PELT_SENSITIVITY_BETAS = (3.0, 10.0, 30.0, 100.0)
# Penalty variance: the noise is strongly autocorrelated, so the iid BIC penalty uses the AR(1) long-run
# variance sigma_e^2 / (1 - phi)^2 of the residual (sigma_e: robust scale after AR(1) prewhitening).
PELT_PENALTY_SCALES = ("long_run_ar1", "innovation_ar1", "marginal")  # first = production, others = sensitivity
PELT_PENALTY_SCALE = PELT_PENALTY_SCALES[0]
EVENT_SIGMA = 5.0  # a break is an event if its level discontinuity is >= 5 sigma (marginal robust scale)
PULSE_TIE_WINDOW_D = (-0.02, 0.5)  # same rule as the NH4 pulse-response classification

# Lomb-Scargle
LS_MIN_CYCLES = 10  # periodicity band: at least 10 cycles in the record
LS_OVERSAMPLING = 5
LS_NYQUIST_C_D = 50.0  # 0.01 d grid

GRID_SERIES = {"pH": "pH", "orp_mV_vs_Ag": "ORP", "EC_mS_cm": "EC", "nh4_mgNH4_L": "NH4"}


# --- Inputs -------------------------------------------------------------------
def load_inputs(out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Grid, every NH4 step-detector candidate (M2), the pulse table (M2) and the held/damped NH4 segments (M2)."""
    paths = {
        "grid": Path(out_dir) / "processed" / "grid_0p01d.parquet",
        "jumps": Path(out_dir) / "tables" / "nh4_jumps.csv",
        "pulses": Path(out_dir) / "tables" / "sensor_validation.csv",
        "held": Path(out_dir) / "tables" / "nh4_held_segments.csv",
    }
    for name, p in paths.items():
        if not p.is_file():
            stage = "ingest" if name == "grid" else "validate"
            raise FileNotFoundError(f"{p} missing: run `python -m bnr_neq {stage}` first")
    grid = pd.read_parquet(paths["grid"])
    return grid, pd.read_csv(paths["jumps"]), pd.read_csv(paths["pulses"]), pd.read_csv(paths["held"])


def step_candidates(jumps: pd.DataFrame, threshold: float = qc.JUMP_THRESHOLD_MG_L) -> pd.DataFrame:
    """Every detector candidate with |step| >= threshold, located or not, sustained or not.

    Each candidate excludes a region from the slope windows: [t10, t90] for a located, sustained jump
    (``jump``), and the whole above-threshold run [run_start, run_end] for an unlocated or non-sustained
    candidate, whose transition time is not known. ``t_event_d`` is the 10-90 % midpoint where the
    transition was located, else the time of the rolling-step peak.
    """
    c = jumps[jumps["step_peak"].abs() >= threshold].copy()
    located = c["t_jump_d"].notna()
    c["t_event_d"] = np.where(located, c["t_jump_d"], c["t_step_peak_d"])
    jump = c["jump"].astype(bool)
    run_lo = c["run_start_d"] if "run_start_d" in c else c["t_event_d"]
    run_hi = c["run_end_d"] if "run_end_d" in c else c["t_event_d"]
    t10 = c["t10_d"].fillna(c["t_event_d"]) if "t10_d" in c else c["t_event_d"]
    t90 = c["t90_d"].fillna(c["t_event_d"]) if "t90_d" in c else c["t_event_d"]
    c["region_start_d"] = np.where(jump, t10, np.minimum(run_lo, c["t_event_d"]))
    c["region_end_d"] = np.where(jump, t90, np.maximum(run_hi, c["t_event_d"]))
    # The left window ends where the region starts; the right one resumes at the 10-90 % midpoint of a
    # located jump (as before; the 0.1 d exclusion covers t90) or at the end of an unlocated run.
    c["t_resume_d"] = np.where(jump, c["t_event_d"], c["region_end_d"])
    c["candidate_kind"] = np.where(
        jump, "jump (" + c["class"].astype(str) + ")",
        "candidate (" + c["class"].astype(str) + ", " + np.where(c["sustained"].astype(bool), "sustained",
                                                                  "non-sustained") + ")")
    return c.sort_values("t_event_d", kind="stable").reset_index(drop=True)


# --- Segmentation -------------------------------------------------------------
def boundaries(pulse_times_d: list[float], jumps: pd.DataFrame,
               start: float = SEGMENT_START_D, end: float = COMMON_END_D) -> pd.DataFrame:
    """Every pulse and every NH4 step-detector candidate with |step| >= threshold is a boundary:
    located jumps (offsets, delayed responses) and unlocated or non-sustained candidates alike."""
    cand = step_candidates(jumps)
    rows = [{"t_d": start, "t_resume_d": start, "kind": "start"}, {"t_d": end, "t_resume_d": end, "kind": "end"}]
    rows += [{"t_d": t, "t_resume_d": t, "kind": f"pulse {k}"} for k, t in enumerate(pulse_times_d, start=1)]
    rows += [{"t_d": a, "t_resume_d": b, "kind": k}
             for a, b, k in zip(cand["region_start_d"], cand["t_resume_d"], cand["candidate_kind"]) if start < b and a < end]
    return pd.DataFrame(rows).sort_values("t_d", kind="stable").reset_index(drop=True)


def sub_segments(bounds: pd.DataFrame, pulse_times_d: list[float]) -> pd.DataFrame:
    """Sub-segment between consecutive boundaries: from the previous boundary's resume time to the next
    boundary's time (equal for pulses and located jumps)."""
    pulses = np.asarray(pulse_times_d, dtype=float)
    rows = []
    prev = bounds[["t_resume_d", "kind"]].to_numpy()[:-1]
    nxt = bounds[["t_d", "kind"]].to_numpy()[1:]
    for (a, ka), (b, kb) in zip(prev, nxt):
        gov = int(np.searchsorted(pulses, a, side="right"))  # pulses at or before a
        fit_start, fit_end = a + EXCLUDE_AFTER_BOUNDARY_D, b
        span = fit_end - fit_start
        rows.append({
            "segment": gov,
            "governing_pulse": gov if gov > 0 else pd.NA,
            "t_start_d": a, "t_end_d": b, "start_boundary": ka, "end_boundary": kb,
            "fit_start_d": fit_start, "fit_end_d": fit_end,
            "status": "computed" if span >= MIN_FIT_SPAN_D else
                      f"not computed: fit span {max(span, 0):.3f} d < {MIN_FIT_SPAN_D} d",
        })
    df = pd.DataFrame(rows)
    df["governing_pulse"] = pd.array(df["governing_pulse"], dtype="Int64")
    df.insert(0, "sub_segment", np.arange(1, len(df) + 1))
    return df


def bins_in(grid: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    """Grid bins lying wholly inside [lo, hi]."""
    return grid[(grid["t_start_d"] >= lo - 1e-9) & (grid["t_end_d"] <= hi + 1e-9)]


def theil_sen(t: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    res = theilslopes(y, t, alpha=THEIL_SEN_ALPHA)
    return float(res.slope), float(res.low_slope), float(res.high_slope)


def scenario_T(cfg: Config) -> list[float]:
    return [float(T) for T in cfg.get("scenarios.temperature_C")]


def held_overlap_d(lo: float, hi: float, held: pd.DataFrame) -> float:
    """Total time (d) of [lo, hi] covered by held/damped NH4 segments."""
    return float(sum(max(0.0, min(hi, b) - max(lo, a)) for a, b in zip(held["t_start_d"], held["t_end_d"])))


def bins_outside_held(b: pd.DataFrame, held: pd.DataFrame) -> pd.DataFrame:
    """Grid bins that do not overlap any held/damped NH4 segment."""
    keep = np.ones(len(b), dtype=bool)
    for a, e in zip(held["t_start_d"], held["t_end_d"]):
        keep &= ~((b["t_end_d"].to_numpy() > a + 1e-9) & (b["t_start_d"].to_numpy() < e - 1e-9))
    return b[keep]


NON_POSITIVE = "not computed: non-positive gain-corrected consumption (rising electrode signal)"


def electrode_delay_d(jumps: pd.DataFrame) -> float:
    """Longest delay (d) between a pulse and its pulse-response jump (M2 ``lag_h``); NaN if none is located."""
    if "lag_h" not in jumps.columns:
        return float("nan")
    lag = pd.to_numeric(jumps.loc[jumps["class"] == "pulse_response", "lag_h"], errors="coerce").dropna()
    return float(lag.max() / 24.0) if len(lag) else float("nan")


def electrode_fragment(start_boundary: str, end_boundary: str, governing_pulse, span_d: float, delay_d: float) -> bool:
    """R8: True for a baseline window between two electrode steps that is shorter than FRAGMENT_DELAYS delays."""
    both = all(str(b).startswith(ELECTRODE_BOUNDARY_PREFIXES) for b in (start_boundary, end_boundary))
    return bool(pd.isna(governing_pulse) and both and np.isfinite(delay_d) and span_d < FRAGMENT_DELAYS * delay_d)


def slope_table(cfg: Config, grid: pd.DataFrame, jumps: pd.DataFrame, pulse_table: pd.DataFrame,
                held: pd.DataFrame) -> pd.DataFrame:
    n_factor = cfg.get("sensors.nh4_to_N_factor")
    segs = sub_segments(boundaries(cfg.pulse_times_d, jumps), cfg.pulse_times_d)
    fractions = dict(zip(pulse_table["pulse"], pulse_table["nh4_response_fraction"]))
    delay = electrode_delay_d(jumps)
    rows = []
    for s in segs.itertuples(index=False):
        row = s._asdict()
        b = bins_in(grid, s.fit_start_d, s.fit_end_d)
        row["n_bins"] = len(b)
        status = s.status
        if status == "computed" and len(b) < MIN_FIT_BINS:
            status = row["status"] = f"not computed: {len(b)} grid bins in the fit window (< {MIN_FIT_BINS})"
        row["slope_fragment_diagnostic_mgNH4_L_d"] = np.nan
        if status == "computed" and electrode_fragment(s.start_boundary, s.end_boundary, s.governing_pulse,
                                                       s.fit_end_d - s.fit_start_d, delay):
            row["slope_fragment_diagnostic_mgNH4_L_d"] = theil_sen(b["t_mid_d"].to_numpy(), b["nh4_mgNH4_L"].to_numpy())[0]
            status = row["status"] = (f"not computed: electrode fragment (baseline window between two electrode steps, "
                                      f"fit span {s.fit_end_d - s.fit_start_d:.2f} d < {FRAGMENT_DELAYS:g} x the "
                                      f"longest delayed electrode response {delay:.2f} d; slope kept as a diagnostic)")
        span = s.fit_end_d - s.fit_start_d
        row["held_overlap_d"] = held_overlap_d(s.fit_start_d, s.fit_end_d, held) if span > 0 else 0.0
        row["held_fraction"] = row["held_overlap_d"] / span if span > 0 else np.nan
        row["median_pH"] = float(b["pH"].median()) if len(b) else np.nan
        row["median_nh4_mgN_L"] = float(b["nh4_mgN_L"].median()) if len(b) else np.nan
        for T in scenario_T(cfg):
            col = f"{FA_READING}_SCENARIO_T{T:g}C"
            row[col] = (float(np.median(thermo.free_ammonia_mgNH3_L(
                b["nh4_mgN_L"].to_numpy(), b["pH"].to_numpy(), T + KELVIN_OFFSET))) if len(b) else np.nan)
        if status == "computed":
            slope, lo, hi = theil_sen(b["t_mid_d"].to_numpy(), b["nh4_mgNH4_L"].to_numpy())
        else:
            slope = lo = hi = np.nan
        row.update(slope_mgNH4_L_d=slope, slope_ci95_low_mgNH4_L_d=lo, slope_ci95_high_mgNH4_L_d=hi,
                   slope_mgN_L_d=slope * n_factor, slope_ci95_low_mgN_L_d=lo * n_factor,
                   slope_ci95_high_mgN_L_d=hi * n_factor)
        # Sensitivity: the same window with the held/damped NH4 spans removed
        bh = bins_outside_held(b, held)
        if status != "computed":
            row["slope_excl_held_status"] = status
        elif row["held_overlap_d"] == 0:
            row["slope_excl_held_status"] = "computed (no held span in the window)"
        elif len(bh) < MIN_FIT_SPAN_D / 0.01:
            row["slope_excl_held_status"] = (f"not computed: {len(bh)} bins outside held spans "
                                             f"(< {MIN_FIT_SPAN_D} d)")
        else:
            row["slope_excl_held_status"] = "computed"
        if row["slope_excl_held_status"].startswith("computed"):
            sh, sh_lo, sh_hi = theil_sen(bh["t_mid_d"].to_numpy(), bh["nh4_mgNH4_L"].to_numpy())
        else:
            sh = sh_lo = sh_hi = np.nan
        row.update(slope_excl_held_mgNH4_L_d=sh, slope_excl_held_ci95_low_mgNH4_L_d=sh_lo,
                   slope_excl_held_ci95_high_mgNH4_L_d=sh_hi)
        if pd.isna(s.governing_pulse):
            f, gain_status = np.nan, "not computed: no governing pulse (baseline)"
        elif status != "computed":
            f, gain_status = fractions[int(s.governing_pulse)], status
        elif not -slope * n_factor / fractions[int(s.governing_pulse)] > 0:
            f, gain_status = fractions[int(s.governing_pulse)], NON_POSITIVE
        else:
            f, gain_status = fractions[int(s.governing_pulse)], "computed"
        row["response_fraction_governing_pulse"] = f
        # GAIN-CORRECTED (electrode gain model): consumption = -slope / response fraction; only a positive
        # consumption is reported (as thermo.oxygen_table does), a rising signal is not a consumption rate
        ok = gain_status == "computed"
        row["consumption_gain_corrected_mgN_L_d"] = -slope * n_factor / f if ok else np.nan
        row["consumption_gain_corrected_ci95_low_mgN_L_d"] = -hi * n_factor / f if ok else np.nan
        row["consumption_gain_corrected_ci95_high_mgN_L_d"] = -lo * n_factor / f if ok else np.nan
        row["gain_correction_status"] = gain_status
        rows.append(row)
    out = pd.DataFrame(rows)
    # W1 provenance: these electrode-scale columns are not the NH4 truth (see truth_slope_table)
    out["gain_corrected_role"] = GAIN_ROLE
    out["fa_basis"] = FA_READING_ROLE
    return out


GAIN_ROLE = ("SENSITIVITY: linear electrode gain model (consumption = −slope / response fraction); the single NH4 "
             "observation model is the M6 log model (nh4_truth_windows.csv)")
FA_READING = "fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L"  # W1: the SENSITIVITY label is part of the column name
FA_READING_ROLE = ("SENSITIVITY: electrode reading taken as TAN (no observation model); the primary free ammonia uses "
                   "the M6-inverted TAN (nh4_truth_windows.csv)")


# --- W1: slopes of the single NH4 truth on the same windows ----------------------------
NON_POSITIVE_TRUTH = "not computed: non-positive consumption of the M6-inverted NH4 (C(t) rises in the window)"


def truth_bins(truth: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    """Rows of the inverted series (0.01 d bin mid times) whose grid bin lies wholly inside [lo, hi]."""
    t = truth["t_d"]
    return truth[(t - 0.005 >= lo - 1e-9) & (t + 0.005 <= hi + 1e-9)]


def median_pairwise_slope(t: np.ndarray, y: np.ndarray) -> float:
    """Theil-Sen point estimate (median of pairwise slopes), without the interval."""
    i, j = np.triu_indices(len(t), k=1)
    return float(np.median((y[j] - y[i]) / (t[j] - t[i])))


def truth_slope_table(cfg: Config, truth: pd.DataFrame, table: pd.DataFrame,
                      band_series: list[np.ndarray] | None = None, volumes=None) -> pd.DataFrame:
    """Consumption, free ammonia and oxygen demand of the M6-inverted C(t) on the M3 slope windows (W1).

    ``truth`` is ``nh4_truth.parquet``; ``band_series`` are the inversions under the profile vectors (their slope
    range is the parameter band); ``volumes`` adds the consumption under each non-constant volume scenario,
    r = −dC/dt + C E/V(t) (evaporation concentrates C at the rate C E/V).
    """
    from .model import OBSERVATION_MODEL, truth_id

    c_all = truth["C_inv_mgN_L"].to_numpy()
    tid = str(truth["nh4_truth_id"].iloc[0])
    if tid != truth_id(c_all):
        raise ValueError("nh4_truth_id does not match the inverted series it labels")
    temps = scenario_T(cfg)
    rows = []
    for s in table.itertuples(index=False):
        b = truth_bins(truth, s.fit_start_d, s.fit_end_d)
        row = {"sub_segment": s.sub_segment, "segment": s.segment, "governing_pulse": s.governing_pulse,
               "fit_start_d": s.fit_start_d, "fit_end_d": s.fit_end_d, "n_bins": len(b), "held_fraction": s.held_fraction}
        status = s.status if s.status != "computed" or len(b) >= MIN_FIT_BINS else \
            f"not computed: {len(b)} grid bins in the fit window (< {MIN_FIT_BINS})"
        c = b["C_inv_mgN_L"].to_numpy()
        row["median_C_inv_mgN_L"] = float(np.median(c)) if len(b) else np.nan
        row["median_pH"] = float(b["pH"].median()) if len(b) else np.nan
        for T in temps:
            row[f"fa_median_mgNH3_L_SCENARIO_T{T:g}C"] = (float(np.median(thermo.free_ammonia_mgNH3_L(
                c, b["pH"].to_numpy(), T + KELVIN_OFFSET))) if len(b) else np.nan)
        cons = lo = hi = band_lo = band_hi = np.nan
        if status == "computed":
            tt = b["t_d"].to_numpy()
            slope, s_lo, s_hi = theil_sen(tt, c)
            cons, lo, hi = -slope, -s_hi, -s_lo
            if band_series:
                idx = truth.index.get_indexer(b.index)
                cb = [-median_pairwise_slope(tt, bs[idx]) for bs in band_series]
                band_lo, band_hi = float(min(cb)), float(max(cb))
            if not cons > 0:
                status = NON_POSITIVE_TRUTH
        ok = status == "computed"
        row.update({"consumption_mgN_L_d": cons if ok else np.nan, "consumption_ci95_low_mgN_L_d": lo if ok else np.nan,
                    "consumption_ci95_high_mgN_L_d": hi if ok else np.nan,
                    "consumption_param_band_low_mgN_L_d": band_lo if ok else np.nan,
                    "consumption_param_band_high_mgN_L_d": band_hi if ok else np.nan, "status": status})
        row["o2_demand_mgO2_L_d"] = float(thermo.oxygen_demand_mgO2_L_d(cons)) if ok else np.nan
        row["o2_demand_ci95_low_mgO2_L_d"] = float(thermo.oxygen_demand_mgO2_L_d(max(lo, 0.0))) if ok else np.nan
        row["o2_demand_ci95_high_mgO2_L_d"] = float(thermo.oxygen_demand_mgO2_L_d(max(hi, 0.0))) if ok else np.nan
        row["kla_required_per_d_SCENARIO_T25C"] = float(thermo.required_kla_per_d(row["o2_demand_mgO2_L_d"])) if ok \
            else np.nan
        for v in volumes or []:
            if v.is_constant:
                continue
            col = f"consumption_mgN_L_d_{volume.slug(v.name)}"
            tm = 0.5 * (s.fit_start_d + s.fit_end_d)
            row[col] = (cons + row["median_C_inv_mgN_L"] * v.net_loss_L_per_d / float(v.volume_L(tm))) if ok else np.nan
        row["nh4_observation_model"] = OBSERVATION_MODEL
        row["nh4_truth_id"] = tid
        rows.append(row)
    out = pd.DataFrame(rows)
    out["governing_pulse"] = pd.array(out["governing_pulse"], dtype="Int64")
    return out


def windows_check(table: pd.DataFrame, pulse_times_d: list[float], jumps: pd.DataFrame) -> pd.DataFrame:
    """For each computed slope window: does it contain a pulse or any step-detector candidate with
    |step| >= threshold (located or unlocated, sustained or not)? A window is clean only if neither
    any candidate's excluded region (see :func:`step_candidates`) nor a pulse intersects its interior."""
    cand = step_candidates(jumps)
    pulses = np.asarray(pulse_times_d, dtype=float)
    lo = np.r_[pulses, cand["region_start_d"].to_numpy()]
    hi = np.r_[pulses, cand["region_end_d"].to_numpy()]
    labels = [f"pulse {k}" for k in range(1, len(pulses) + 1)] + cand["candidate_kind"].tolist()
    rows = []
    for r in table[table["status"] == "computed"].itertuples(index=False):
        inside = (hi > r.fit_start_d + 1e-9) & (lo < r.fit_end_d - 1e-9)
        rows.append({"sub_segment": r.sub_segment, "fit_start_d": r.fit_start_d, "fit_end_d": r.fit_end_d,
                     "events_inside": ";".join(f"{a:.3f}-{b:.3f} ({lab})" for a, b, lab, m in zip(lo, hi, labels, inside)
                                               if m),
                     "held_fraction": r.held_fraction, "clean": not inside.any()})
    return pd.DataFrame(rows)


def reference_window_table(cfg: Config, grid: pd.DataFrame, jumps: pd.DataFrame,
                           table: pd.DataFrame) -> pd.DataFrame:
    """Diagnostic: Theil-Sen on the exact preliminary windows (first 0.1 d excluded),
    next to the production sub-segment with the largest overlap."""
    pulses = np.asarray(cfg.pulse_times_d, dtype=float)
    jt = step_candidates(jumps)["t_event_d"].to_numpy()
    comp = table[table["status"] == "computed"]
    rows = []
    for (a, b), ref in zip(REFERENCE_WINDOWS, REFERENCE_SLOPES_MGNH4_L_D):
        lo, hi = a + EXCLUDE_AFTER_BOUNDARY_D, b
        bb = bins_in(grid, lo, hi)
        slope, ci_lo, ci_hi = theil_sen(bb["t_mid_d"].to_numpy(), bb["nh4_mgNH4_L"].to_numpy())
        overlap = (np.minimum(comp["fit_end_d"], hi) - np.maximum(comp["fit_start_d"], lo)).clip(lower=0)
        k = int(np.argmax(overlap.to_numpy()))
        prod = comp.iloc[k]
        tol = max(0.1 * abs(ref), 0.05)
        rows.append({
            "ref_from_d": a, "ref_to_d": b, "reference_slope_mgNH4_L_d": ref, "tolerance": tol,
            "exact_window_slope_mgNH4_L_d": slope, "exact_window_ci95": f"[{ci_lo:.3f}, {ci_hi:.3f}]",
            "exact_window_pass": abs(slope - ref) <= tol + 1e-9,
            "exact_window_spans_pulse": bool(((pulses > lo) & (pulses < hi)).any()),
            "exact_window_spans_jump": bool(((jt > lo) & (jt < hi)).any()),
            "production_sub_segment": int(prod["sub_segment"]),
            "production_fit_window_d": f"{prod['fit_start_d']:.3f} to {prod['fit_end_d']:.3f}",
            "production_slope_mgNH4_L_d": prod["slope_mgNH4_L_d"],
            "production_pass": abs(prod["slope_mgNH4_L_d"] - ref) <= tol + 1e-9,
        })
    return pd.DataFrame(rows)


# --- Kinetic models -----------------------------------------------------------
@dataclass(frozen=True)
class KineticModel:
    name: str
    params: tuple[str, ...]  # kinetic parameters (fitted on log scale)
    T_C: float | None = None  # temperature scenario for the FA-based model


def model_set(cfg: Config) -> list[KineticModel]:
    models = [KineticModel("zero_order", ("k0_mgN_L_d",)),
              KineticModel("first_order", ("k1_per_d",)),
              KineticModel("monod", ("vmax_mgN_L_d", "Ks_mgN_L"))]
    models += [KineticModel(f"haldane_FA_SCENARIO_T{T:g}C", ("vmax_mgN_L_d", "Ks_mgNH3_L", "Ki_mgNH3_L"), T)
               for T in scenario_T(cfg)]
    return models


def _rk4_haldane(c0: float, t: np.ndarray, ph: np.ndarray, vmax: float, ks: float, ki: float, T_K: float) -> np.ndarray:
    """dC/dt = -vmax S/(Ks + S + S^2/Ki), S = FA(C, pH(t), T); RK4 on the bin times."""
    def rate(c, p):
        s = NH3_PER_N * max(c, 0.0) * 10 ** p / (np.exp(ANTHONISEN_EXP_K / T_K) + 10 ** p)
        return -vmax * s / (ks + s + s * s / ki)
    out = np.empty(len(t))
    out[0] = c = c0
    for i in range(1, len(t)):
        h = t[i] - t[i - 1]
        pm = 0.5 * (ph[i - 1] + ph[i])
        k1 = rate(c, ph[i - 1])
        k2 = rate(c + 0.5 * h * k1, pm)
        k3 = rate(c + 0.5 * h * k2, pm)
        k4 = rate(c + h * k3, ph[i])
        c = c + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        out[i] = c
    return out


def simulate(model: KineticModel, kin: np.ndarray, c0: float, t: np.ndarray, ph: np.ndarray) -> np.ndarray:
    """NH4 (mg N/L, electrode scale) from c0 at t[0]."""
    tau = t - t[0]
    if model.name == "zero_order":
        return c0 - kin[0] * tau
    if model.name == "first_order":
        return c0 * np.exp(-kin[0] * tau)
    if model.name == "monod":
        vmax, ks = kin
        # exact solution: C + Ks ln C = C0 + Ks ln C0 - vmax t  ->  C = Ks W(exp(z))
        z = np.log(c0 / ks) + (c0 - vmax * tau) / ks
        return ks * np.real(wrightomega(z))
    return _rk4_haldane(c0, t, ph, *kin, model.T_C + KELVIN_OFFSET)


def fit_segment(model: KineticModel, pieces: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
                slope_guess: float) -> dict:
    """Shared kinetic parameters, one free initial level per sub-segment (absorbs electrode offsets)."""
    c_med = float(np.median(np.concatenate([p[1] for p in pieces])))
    rate = max(abs(slope_guess), 0.05)
    guesses = {"k0_mgN_L_d": rate, "k1_per_d": rate / c_med, "vmax_mgN_L_d": 2 * rate,
               "Ks_mgN_L": c_med, "Ks_mgNH3_L": 0.5, "Ki_mgNH3_L": 1.0}
    x0 = np.r_[np.log([guesses[p] for p in model.params]), [p[1][0] for p in pieces]]
    nk = len(model.params)
    lo = np.r_[np.full(nk, LOG_BOUNDS[0]), np.full(len(pieces), 1e-3)]
    hi = np.r_[np.full(nk, LOG_BOUNDS[1]), np.full(len(pieces), 500.0)]
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    def resid(x):
        kin = np.exp(x[:nk])
        return np.concatenate([simulate(model, kin, x[nk + j], t, ph) - c for j, (t, c, ph) in enumerate(pieces)])

    sol = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", method="trf")
    r = sol.fun
    n, p = len(r), len(sol.x)
    rss = float(r @ r)
    sigma2 = rss / max(n - p, 1)
    jac = sol.jac
    jtj = jac.T @ jac
    cond = float(np.linalg.cond(jtj))
    cov = sigma2 * np.linalg.pinv(jtj)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    params = []
    for i, name in enumerate(model.params):
        at_bound = bool(min(sol.x[i] - LOG_BOUNDS[0], LOG_BOUNDS[1] - sol.x[i]) < 1e-3)
        rel_se = float(se[i])  # SE of the log parameter = relative SE
        reasons = [m for m, bad in (("at bound", at_bound), (f"relative SE {rel_se:.2g} > {REL_SE_MAX}",
                                    not rel_se <= REL_SE_MAX), ("ill-conditioned", cond > 1e12)) if bad]
        params.append({"parameter": name, "estimate": float(np.exp(sol.x[i])), "relative_se": rel_se,
                       "identifiable": not reasons, "reason": "; ".join(reasons)})
    lb = acorr_ljungbox(r, lags=[LJUNG_BOX_LAG], return_df=True)
    return {
        "model": model.name, "n": n, "n_params": p, "rss": rss, "rmse_mgN_L": float(np.sqrt(rss / n)),
        "aic": n * np.log(rss / n) + 2 * p,
        "resid_lag1_autocorr": float(np.corrcoef(r[:-1], r[1:])[0, 1]),
        "ljung_box_p_lag10": float(lb["lb_pvalue"].iloc[0]),
        "converged": bool(sol.success), "params": params,
    }


def kinetics(cfg: Config, grid: pd.DataFrame, table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fits, params = [], []
    comp = table[table["status"] == "computed"]
    for seg, g in comp.groupby("segment"):
        pieces = []
        for r in g.itertuples(index=False):
            b = bins_in(grid, r.fit_start_d, r.fit_end_d)
            pieces.append((b["t_mid_d"].to_numpy(), b["nh4_mgN_L"].to_numpy(), b["pH"].to_numpy()))
        slope_guess = float(np.median(g["slope_mgN_L_d"]))
        for m in model_set(cfg):
            f = fit_segment(m, pieces, slope_guess)
            ps = f.pop("params")
            fits.append({"segment": seg, "sub_segments": ";".join(str(s) for s in g["sub_segment"]), **f,
                         "non_identifiable": ";".join(p["parameter"] for p in ps if not p["identifiable"])})
            params += [{"segment": seg, "model": m.name, **p} for p in ps]
    fits = pd.DataFrame(fits)
    fits["delta_aic"] = fits["aic"] - fits.groupby("segment")["aic"].transform("min")
    return fits, pd.DataFrame(params)


# --- Change points ------------------------------------------------------------
def _robust_scale(r: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(r - np.median(r))))


def noise_residual(y: np.ndarray) -> np.ndarray:
    return y - pd.Series(y).rolling(PELT_NOISE_WINDOW, center=True, min_periods=1).median().to_numpy()


def noise_sigma(y: np.ndarray) -> float:
    """Marginal robust scale of the residual about the rolling median (event criterion)."""
    return _robust_scale(noise_residual(y))


def penalty_scales(y: np.ndarray) -> dict[str, float]:
    """Residual scales for the PELT penalty: AR(1) prewhitening of the residual about the rolling median.

    phi = lag-1 autocorrelation of the residual; e = r_t − phi r_{t−1}; sigma_e = robust scale of e.
    ``long_run_ar1`` = sigma_e / (1 − phi), the scale of the AR(1) long-run variance, which is what a
    piecewise-linear RSS reduction is driven by when the noise is autocorrelated.
    """
    r = noise_residual(y)
    phi = float(np.corrcoef(r[:-1], r[1:])[0, 1])
    phi = min(max(phi, 0.0), 0.999)
    sigma_e = _robust_scale(r[1:] - phi * r[:-1])
    return {"phi": phi, "marginal": _robust_scale(r), "innovation_ar1": sigma_e, "long_run_ar1": sigma_e / (1 - phi)}


def pelt_breaks(t: np.ndarray, y: np.ndarray, beta: float, scale: str = PELT_PENALTY_SCALE) -> list[int]:
    sigma = penalty_scales(y)[scale]
    X = np.column_stack([y, t, np.ones_like(t)])
    algo = rpt.Pelt(model=PELT_MODEL, min_size=PELT_MIN_SIZE, jump=PELT_JUMP).fit(X)
    return algo.predict(pen=beta * sigma**2 * np.log(len(y)))[:-1]


def level_discontinuity(t: np.ndarray, y: np.ndarray, bks: list[int]) -> np.ndarray:
    """Right minus left linear fit, evaluated at each break."""
    edges = [0, *bks, len(y)]
    fits = [np.polyfit(t[a:b], y[a:b], 1) for a, b in itertools.pairwise(edges)]
    return np.array([np.polyval(fits[i + 1], t[k]) - np.polyval(fits[i], t[k]) for i, k in enumerate(bks)])


def tied_pulse(tc: float, pulse_times_d: list[float]) -> int | None:
    lags = tc - np.asarray(pulse_times_d, dtype=float)
    ok = (lags >= PULSE_TIE_WINDOW_D[0]) & (lags <= PULSE_TIE_WINDOW_D[1])
    return int(np.flatnonzero(ok)[-1]) + 1 if ok.any() else None


def change_points(grid: pd.DataFrame, pulse_times_d: list[float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    t = grid["t_mid_d"].to_numpy()
    rows, sens = [], []
    for col, label in GRID_SERIES.items():
        y = grid[col].to_numpy()
        sigma = noise_sigma(y)
        scales = penalty_scales(y)
        bks = pelt_breaks(t, y, PELT_BETA)
        # sensitivity: every beta with the production scale, and the other scales at the production beta
        grid_sens = [(PELT_PENALTY_SCALE, b) for b in PELT_SENSITIVITY_BETAS]
        grid_sens += [(sc, PELT_BETA) for sc in PELT_PENALTY_SCALES if sc != PELT_PENALTY_SCALE]
        for scale, beta in grid_sens:
            production = scale == PELT_PENALTY_SCALE and beta == PELT_BETA
            n = len(bks) if production else len(pelt_breaks(t, y, beta, scale))
            sens.append({"series": label, "penalty_scale": scale, "production": production, "beta": beta,
                         "residual_lag1_phi": scales["phi"], "sigma_penalty": scales[scale], "n_breaks": n})
        jumps = level_discontinuity(t, y, bks)
        for k, d in zip(bks, jumps):
            p = tied_pulse(t[k], pulse_times_d)
            rows.append({"series": label, "variable": col, "t_break_d": t[k], "level_discontinuity": d,
                         "noise_sigma": sigma, "discontinuity_sigma": d / sigma,
                         "event": abs(d) >= EVENT_SIGMA * sigma, "tied_pulse": p})
    cols = ["series", "variable", "t_break_d", "level_discontinuity", "noise_sigma", "discontinuity_sigma", "event",
            "tied_pulse"]
    cp = pd.DataFrame(rows, columns=cols)
    cp["event"] = cp["event"].astype(bool)
    cp["tied_pulse"] = pd.array(cp["tied_pulse"], dtype="Int64")
    ev = cp[cp["event"]]
    cp["co_occurring_series"] = [
        ";".join(sorted(set(ev.loc[(ev["series"] != r.series) & ((ev["t_break_d"] - r.t_break_d).abs() <= 0.05),
                                   "series"]))) if r.event else ""
        for r in cp.itertuples(index=False)]
    return cp, pd.DataFrame(sens)


# --- Periodicity --------------------------------------------------------------
def periodograms(grid: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    t = grid["t_mid_d"].to_numpy()
    span = t[-1] - t[0] + 0.01
    f = np.arange(1 / span, LS_NYQUIST_C_D, 1 / (LS_OVERSAMPLING * span))
    band = f >= LS_MIN_CYCLES / span
    near1 = (f >= 0.9) & (f <= 1.1)
    spectra, summary = {"f_cycles_per_d": f}, []
    for col, label in GRID_SERIES.items():
        y = grid[col].to_numpy()
        yd = y - np.polyval(np.polyfit(t, y, 1), t)  # linear detrend over the common window
        p = lombscargle(t, yd - yd.mean(), 2 * np.pi * f, normalize=True)
        spectra[label] = p
        for name, m in (("full (1/T to Nyquist)", np.ones_like(band)), (f"periodicity (>= {LS_MIN_CYCLES} cycles)", band),
                        ("near 1 c/d (0.9 to 1.1)", near1)):
            k = int(np.argmax(np.where(m, p, -np.inf)))
            summary.append({"series": label, "band": name, "f_min_c_d": float(f[m].min()),
                            "max_normalized_power": float(p[k]), "f_at_max_c_d": float(f[k]),
                            "period_at_max_d": float(1 / f[k])})
    return pd.DataFrame(summary), pd.DataFrame(spectra)


# --- References ---------------------------------------------------------------
def reference_comparison(table: pd.DataFrame, refwin: pd.DataFrame, ls: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in refwin.itertuples(index=False):
        rows.append({"quantity": f"slope_mgNH4_L_d {r.ref_from_d:g}-{r.ref_to_d:g} (production sub-segment)",
                     "reference": r.reference_slope_mgNH4_L_d, "obtained": r.production_slope_mgNH4_L_d,
                     "tolerance": r.tolerance, "note": f"fit window {r.production_fit_window_d} d; exact window gives "
                     f"{r.exact_window_slope_mgNH4_L_d:.3f}"})
    comp = table[table["status"] == "computed"]
    first_after = {k: comp[comp["governing_pulse"] == k].iloc[0] for k in (2, 4)}
    for k, ref in ((2, 5.2), (4, 1.4)):
        r = first_after[k]
        rows.append({"quantity": f"consumption_gain_corrected_mgN_L_d after pulse {k}", "reference": ref,
                     "obtained": r["consumption_gain_corrected_mgN_L_d"], "tolerance": 0.1 * ref,
                     "note": f"sub-segment {r['sub_segment']} ({r['fit_start_d']:.3f} to {r['fit_end_d']:.3f} d); "
                             "tolerance: 10 % (reference says 'about')"})
    fa_cols = [c for c in table.columns if c.startswith(f"{FA_READING}_SCENARIO_")]
    fa = comp[(comp["t_start_d"] >= 17.8) & (comp["t_end_d"] <= 32.4)][fa_cols].to_numpy()
    for name, ref, val in (("min", 0.7, fa.min()), ("max", 1.2, fa.max())):
        rows.append({"quantity": f"free ammonia {name}, day 17.84-32.36, 20/25 C SCENARIOS (mg NH3/L; SENSITIVITY: "
                                 "reading as TAN, as in the preliminary analysis)",
                     "reference": ref, "obtained": float(val), "tolerance": 0.05,
                     "note": "range of sub-segment medians; tolerance: rounding of the reference to 0.1"})
    band = ls[ls["band"].str.startswith("periodicity")]
    k = band["max_normalized_power"].idxmax()
    full = ls[ls["band"].str.startswith("full")]["max_normalized_power"].max()
    rows.append({"quantity": "max Lomb-Scargle normalized power (periodicity band)", "reference": 0.05,
                 "obtained": float(band.loc[k, "max_normalized_power"]), "tolerance": 0.01,
                 "note": f"{band.loc[k, 'series']} at {band.loc[k, 'f_at_max_c_d']:.3f} c/d; full band max {full:.3f}; "
                         "tolerance: 0.01 (reference says 'about')"})
    out = pd.DataFrame(rows)
    out["abs_diff"] = (out["obtained"] - out["reference"]).abs()
    out["pass"] = out["abs_diff"] <= out["tolerance"] + 1e-9
    return out[["quantity", "reference", "obtained", "tolerance", "abs_diff", "pass", "note"]]


# --- Figure -------------------------------------------------------------------
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"


def figure(cfg: Config, grid: pd.DataFrame, table: pd.DataFrame, jumps: pd.DataFrame, cp: pd.DataFrame,
           held: pd.DataFrame | None, ec_pulses: pd.DataFrame | None, png: Path, svg: Path) -> None:
    from .report import pyplot

    plt = pyplot()

    t = grid["t_mid_d"].to_numpy()
    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True, facecolor=SURFACE)
    panels = [("nh4_mgNH4_L", "NH4+ (mg NH4+/L)", "NH4"), ("pH", "pH", "pH"),
              ("orp_mV_vs_Ag", "ORP (mV vs bare Ag, relative)", "ORP"), ("EC_mS_cm", "EC (mS/cm, unvalidated)", "EC")]
    for ax, (col, lab, series) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.axvspan(0, SEGMENT_START_D, color=GRID, alpha=0.6, lw=0)
        ax.plot(t, grid[col], color=BLUE, lw=0.8)
        for tp in cfg.pulse_times_d:
            ax.axvline(tp, color=INK2, lw=0.8, ls=":")
        ev = cp[(cp["series"] == series) & cp["event"] & cp["tied_pulse"].isna()]
        y = np.interp(ev["t_break_d"], t, grid[col])
        ax.scatter(ev["t_break_d"], y, marker="D", s=22, color=ORANGE, edgecolor=SURFACE, linewidth=1, zorder=3,
                   label="Change-point event not tied to a pulse")
        ax.set_ylabel(lab, fontsize=8.5, color=INK)
    ax = axes[0]
    if held is not None:
        for h in held.itertuples(index=False):
            ax.axvspan(h.t_start_d, h.t_end_d, ymin=0.95, ymax=1.0, color=INK2, alpha=0.5, lw=0)
    comp = table[table["status"] == "computed"]
    for r in comp.itertuples(index=False):
        b = bins_in(grid, r.fit_start_d, r.fit_end_d)
        mid_t, mid_y = b["t_mid_d"].median(), b["nh4_mgNH4_L"].median()
        tt = np.array([r.fit_start_d, r.fit_end_d])
        ax.plot(tt, mid_y + r.slope_mgNH4_L_d * (tt - mid_t), color=AQUA, lw=2.2,
                label="Theil-Sen slope window" if r.sub_segment == comp["sub_segment"].iloc[0] else None)
    located = jumps[jumps["jump"].astype(bool)]
    for cls, color, marker, lab in (("pulse_response", AQUA, "^", "NH4 jump: delayed pulse response"),
                                    ("electrode_offset", VIOLET, "v", "NH4 jump: electrode offset")):
        s = located[located["class"] == cls]
        ax.scatter(s["t_jump_d"], s["post_level"], marker=marker, s=40, color=color, edgecolor=SURFACE,
                   linewidth=1.2, zorder=4, label=lab)
    other = step_candidates(jumps)
    other = other[~other["jump"].astype(bool)]
    ax.scatter(other["t_event_d"], other["post_level"], marker="x", s=30, color=ORANGE, zorder=4,
               label="NH4 step candidate (unlocated or non-sustained): window boundary")
    ax.plot([], [], color=INK2, alpha=0.5, lw=5, label="Held/damped output (bar at top)")
    ax.legend(fontsize=7.5, frameon=False, loc="lower center", ncol=2)
    if ec_pulses is not None and len(ec_pulses):
        axes[3].scatter(ec_pulses["t_peak_d"], np.interp(ec_pulses["t_peak_d"], t, grid["EC_mS_cm"]) + 0.05,
                        marker="|", s=40, color=VIOLET, label="Late EC pulse (probe event)")
    axes[3].legend(fontsize=7.5, frameon=False, loc="upper center")
    axes[3].set_xlabel("Time (d, logger clock). Dotted: NH4Cl pulses. Grey band: start-up (day 0 to 1).",
                       fontsize=8, color=INK2)
    axes[0].set_title("Common-window overview (0.01 d grid) with pulses, flags and electrode events",
                      loc="left", fontsize=10, color=INK)
    axes[3].set_xlim(0, COMMON_END_D)
    fig.tight_layout()
    fig.savefig(png, dpi=150, facecolor=SURFACE)
    fig.savefig(svg, facecolor=SURFACE, metadata={"Date": None})
    plt.close(fig)


# --- Stage --------------------------------------------------------------------
@dataclass
class SegmentResult:
    table: pd.DataFrame
    windows: pd.DataFrame
    reference_windows: pd.DataFrame
    fits: pd.DataFrame
    params: pd.DataFrame
    changepoints: pd.DataFrame
    pelt_sensitivity: pd.DataFrame
    ls_summary: pd.DataFrame
    ls_spectra: pd.DataFrame
    references: pd.DataFrame


def compute(cfg: Config, grid: pd.DataFrame, jumps: pd.DataFrame, pulse_table: pd.DataFrame,
            held: pd.DataFrame) -> SegmentResult:
    """``jumps`` is the full M2 candidate table (``nh4_jumps.csv``), not only the located jumps."""
    table = slope_table(cfg, grid, jumps, pulse_table, held)
    windows = windows_check(table, cfg.pulse_times_d, jumps)
    refwin = reference_window_table(cfg, grid, jumps, table)
    fits, params = kinetics(cfg, grid, table)
    cp, sens = change_points(grid, cfg.pulse_times_d)
    ls, spectra = periodograms(grid)
    return SegmentResult(table, windows, refwin, fits, params, cp, sens, ls, spectra,
                         reference_comparison(table, refwin, ls))


def _optional_table(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.is_file() else None


def write_report(cfg: Config, r: SegmentResult, path: Path) -> None:
    from .report import md_table

    t = r.table
    view = t[["sub_segment", "segment", "t_start_d", "t_end_d", "start_boundary", "fit_start_d", "fit_end_d",
              "n_bins", "slope_mgNH4_L_d", "slope_ci95_low_mgNH4_L_d", "slope_ci95_high_mgNH4_L_d",
              "slope_mgN_L_d", "held_fraction", "status"]]
    held_view = t.loc[t["held_overlap_d"] > 0, ["sub_segment", "fit_start_d", "fit_end_d", "held_fraction",
                                                  "slope_mgNH4_L_d", "slope_excl_held_mgNH4_L_d",
                                                  "slope_excl_held_ci95_low_mgNH4_L_d",
                                                  "slope_excl_held_ci95_high_mgNH4_L_d", "slope_excl_held_status"]]
    n_nonpos = int((t["gain_correction_status"] == NON_POSITIVE).sum())
    sens = r.pelt_sensitivity
    prod_scale = sens[sens["penalty_scale"] == PELT_PENALTY_SCALE].drop_duplicates("series")
    sens_view = sens.pivot_table(index=["series", "penalty_scale"], columns="beta", values="n_breaks").reset_index()
    sens_view.columns = [c if isinstance(c, str) else f"breaks_beta{c:g}" for c in sens_view.columns]
    gain = t[["sub_segment", "governing_pulse", "response_fraction_governing_pulse",
              "consumption_gain_corrected_mgN_L_d", "consumption_gain_corrected_ci95_low_mgN_L_d",
              "consumption_gain_corrected_ci95_high_mgN_L_d", "gain_correction_status"]]
    fa_cols = [c for c in t.columns if c.startswith(f"{FA_READING}_")]
    fa = t.loc[t["status"] == "computed", ["sub_segment", "t_start_d", "t_end_d", "median_pH", *fa_cols]]
    best = r.fits.loc[r.fits.groupby("segment")["aic"].idxmin(), ["segment", "model"]]
    events = r.changepoints[r.changepoints["event"]]
    free = events[events["tied_pulse"].isna()][["series", "t_break_d", "level_discontinuity",
                                                 "discontinuity_sigma", "co_occurring_series"]]
    tied = events[events["tied_pulse"].notna()][["series", "t_break_d", "level_discontinuity", "tied_pulse"]]
    refs = r.references
    n_fail = int((~refs["pass"]).sum())
    lines = [
        "# Segmentation and consumption kinetics (M3)",
        "",
        "Generated by `python -m bnr_neq segments`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Segmentation rule",
        "",
        (f"Segments start at day {SEGMENT_START_D} (start-up excluded) and at each pulse "
        f"({', '.join(f'{p:g}' for p in cfg.pulse_times_d)}), and end at {COMMON_END_D}. Every NH4 step-detector "
        f"candidate of M2 with |step| ≥ {qc.JUMP_THRESHOLD_MG_L:g} mg NH4+/L splits a segment: located, sustained jumps "
        "(electrode offsets and delayed pulse responses) at their 10-90 % midpoint, and unlocated or non-sustained "
        "candidates at their rolling-step peak, so no slope window spans a pulse or any detected step. The first "
        f"{EXCLUDE_AFTER_BOUNDARY_D} d after each "
        f"boundary is excluded; a sub-segment needs at least {MIN_FIT_SPAN_D} d of fit window. Slopes are Theil-Sen "
        "(`scipy.stats.theilslopes`) on the 0.01 d grid bins lying wholly inside the fit window."),
        "",
        ("**Caveat.** The Theil-Sen 95 % intervals assume independent residuals; the NH4 residuals are strongly "
        "autocorrelated (lag-1 0.82 to 0.998 in the kinetic fits), so these intervals are too narrow."),
        "",
        md_table(view),
        "",
        (f"Window check: {int(r.windows['clean'].sum())} of {len(r.windows)} computed windows contain no pulse and no "
        "step-detector candidate (event time or step peak), located or not, sustained or not."),
        "",
        "## Held/damped NH4 output inside slope windows",
        "",
        ("`held_fraction` is the fraction of each fit window covered by the held/damped NH4 segments of M2 "
        "(`nh4_held_segments.csv`). The production slope keeps those spans (flagged); the sensitivity slope below "
        "removes them. Downstream tables (`pulses.csv`, `thermo_oxygen.csv` via `segments.csv`) carry the flag."),
        "",
        md_table(held_view) if len(held_view) else "No slope window overlaps a held segment.",
        "",
        "## SENSITIVITY: GAIN-CORRECTED consumption (linear electrode gain model)",
        "",
        ("Since W1 the single NH4 observation model is the M6 log model: rates, free ammonia, conversion and oxygen "
        "demand are derived from its inverted C(t) in `nh4_truth_windows.csv` (pulses stage). The linear gain "
        "correction below is kept only as a labelled sensitivity case. "
        "Consumption = −slope (mg N/L/d) / response fraction of the governing pulse (M2). This is conditional on the "
        "electrode-gain interpretation of the under-response (S1) and is not an absolute rate (NO-GO until "
        "laboratory NH4-N). The interval ignores the uncertainty of the response fraction. A window whose gain-corrected "
        f"consumption is not positive is `{NON_POSITIVE}` ({n_nonpos} window(s)); its electrode slope is still listed "
        "above."),
        "",
        md_table(gain),
        "",
        "## SENSITIVITY: free ammonia per sub-segment from the electrode reading taken as TAN (temperature SCENARIOS)",
        "",
        ("The primary free ammonia uses the M6-inverted TAN (`nh4_truth_windows.csv`); these columns keep the reading "
        "as TAN for comparison with the preliminary analysis (reference values below)."),
        "",
        md_table(fa),
        "",
        "## Kinetic models",
        "",
        ("Per segment: shared kinetic parameters, one free initial level per sub-segment (which absorbs the electrode "
        "offsets), fitted to NH4 in mg N/L on the electrode scale with `scipy.optimize.least_squares` (log-scale "
        "kinetic parameters, bounds 1e-6 to 1e4). Monod uses the exact Lambert-W solution; Haldane uses free ammonia "
        "(observed pH, temperature SCENARIO) as substrate and inhibitor, integrated with RK4 on the grid. A parameter "
        f"is not identifiable if at a bound, if its relative SE (from the Jacobian) exceeds {REL_SE_MAX}, or if JᵀJ is "
        "ill-conditioned."),
        "",
        ("**AIC comparisons are not reliable here**: every model leaves strongly autocorrelated residuals "
        "(Ljung-Box p ≈ 0), which violates the independence that AIC assumes. Lowest-AIC models are listed for "
        "completeness only."),
        "",
        md_table(r.fits.drop(columns=["rss"])),
        "",
        "Lowest AIC per segment (not a model selection): "
        + ", ".join(f"segment {s}: {m}" for s, m in zip(best["segment"], best["model"])) + ".",
        "",
        md_table(r.params),
        "",
        "## Change points (PELT)",
        "",
        (f"`ruptures.Pelt(model='{PELT_MODEL}', min_size={PELT_MIN_SIZE}, jump={PELT_JUMP})` on each grid series with "
        f"[t, 1] as covariates (discontinuous piecewise-linear, so drift does not create breaks). Penalty "
        f"{PELT_BETA:g} σ_p² ln n (BIC, three parameters per break). The residual r about a {PELT_NOISE_WINDOW}-bin "
        "rolling median is strongly autocorrelated, so σ_p is autocorrelation-aware: r is prewhitened with its lag-1 "
        "autocorrelation φ (e = r_t − φ r_{t−1}), σ_e is the robust scale of e, and σ_p = σ_e/(1 − φ) is the AR(1) "
        "long-run scale. A break is an event if its level discontinuity is at least "
        f"{EVENT_SIGMA:g} σ (σ = marginal robust scale of r). Tied to a pulse: {PULSE_TIE_WINDOW_D[0]} to "
        f"+{PULSE_TIE_WINDOW_D[1]} d after it."),
        "",
        md_table(prod_scale[["series", "residual_lag1_phi", "sigma_penalty"]]),
        "",
        ("Sensitivity: breaks per penalty factor β for the production scale (long-run AR(1)), and at the production "
        "β for the prewhitened innovation scale σ_e and for the marginal scale used before this correction (a smaller "
        "scale gives a smaller penalty and more breaks):"),
        "",
        md_table(sens_view),
        "",
        f"Events not tied to a pulse ({len(free)}):",
        "",
        md_table(free) if len(free) else "None.",
        "",
        f"Events tied to a pulse ({len(tied)}):",
        "",
        md_table(tied) if len(tied) else "None.",
        "",
        "## Periodicity (Lomb-Scargle)",
        "",
        ("Each grid series is linearly detrended over the common window; normalized power from "
        "`scipy.signal.lombscargle(normalize=True)`; frequencies from 1/T to the Nyquist frequency (50 c/d), "
        f"oversampling {LS_OVERSAMPLING}. The periodicity band requires at least {LS_MIN_CYCLES} cycles in the record "
        f"(f ≥ {LS_MIN_CYCLES}/T); below it the power reflects the non-stationary trend, not a period."),
        "",
        md_table(r.ls_summary),
        "",
        "## Reference comparison",
        "",
        (f"{len(refs) - n_fail} of {len(refs)} reproduced within tolerance. The exact preliminary windows are "
        "evaluated diagnostically below (production windows follow the segmentation rule above)."),
        "",
        md_table(refs),
        "",
        md_table(r.reference_windows),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir: Path | str) -> SegmentResult:
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    d = out_dirs(out_dir)
    grid, jumps, pt, held = load_inputs(out_dir)
    r = compute(cfg, grid, jumps, pt, held)
    r.table.to_csv(d["tables"] / "segments.csv", index=False)
    r.windows.to_csv(d["tables"] / "segments_window_check.csv", index=False)
    r.reference_windows.to_csv(d["tables"] / "segments_reference_windows.csv", index=False)
    r.fits.to_csv(d["tables"] / "kinetic_fits.csv", index=False)
    r.params.to_csv(d["tables"] / "kinetic_parameters.csv", index=False)
    r.changepoints.to_csv(d["tables"] / "changepoints.csv", index=False)
    r.pelt_sensitivity.to_csv(d["tables"] / "changepoint_sensitivity.csv", index=False)
    r.ls_summary.to_csv(d["tables"] / "lomb_scargle_summary.csv", index=False)
    r.ls_spectra.to_parquet(d["processed"] / "lomb_scargle_spectra.parquet", index=False)
    r.references.to_csv(d["tables"] / "segments_references.csv", index=False)
    figure(cfg, grid, r.table, jumps, r.changepoints, held, _optional_table(d["tables"] / "ec_late_pulses.csv"),
           d["figures"] / "fig2_overview.png", d["figures"] / "fig2_overview.svg")
    write_report(cfg, r, d["reports"] / "segments.md")
    return r
