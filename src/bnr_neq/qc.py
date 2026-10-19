"""Quality flags and the event detection used by the pipeline.

All detectors work on native (or reconstructed native) series with time in
days. Thresholds are module constants, documented here and repeated in the
reports.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from .constants import SECONDS_PER_DAY

# --- NH4 jump detector --------------------------------------------------------
JUMP_WINDOW_D = 0.05  # one-sided median window (72 min): >> 14.4-min stamp, << consumption drift
JUMP_THRESHOLD_MG_L = 4.0  # |step| in mg NH4+/L; > 4x the largest sub-segment HF residual (0.92)
JUMP_MERGE_D = 0.1  # above-threshold samples closer than this belong to one jump
PULSE_LAG_MIN_D = -0.02  # two stamps: tolerance for timestamp resolution before the pulse
PULSE_LAG_MAX_D = 0.5  # a positive jump within 0.5 d after a pulse is a pulse response
LEVEL_OFFSET_D = 0.05  # pre/post levels are taken 0.05 to 0.15 d away from the jump
LEVEL_WINDOW_D = 0.10
TRANSITION_SMOOTH = "2min"  # rolling median before locating the 10 and 90 % crossings

# --- Held / damped NH4 output -------------------------------------------------
HF_MEDIAN_WINDOW = "30min"  # centred rolling median; residual = x - median
HF_BLOCK_D = 0.1
HELD_SIGMA_MG_L = 0.01  # reporting resolution of the NH4 output
HELD_MIN_D = 0.5

# --- Late conductivity pulses -------------------------------------------------
EC_BASELINE_WINDOW = "6h"
EC_LATE_START_D = 41.9
EC_THRESHOLD_QUANTILE = 0.999  # of the residual before EC_LATE_START_D
EC_MIN_SEPARATION_MIN = 30.0
EC_FIT_WINDOW_D = 2.0 / 24.0
EC_SENSITIVITY_THRESHOLDS = (0.02, 0.03, 0.04, 0.05, 0.06)


def _series(t: np.ndarray, x: np.ndarray) -> pd.Series:
    return pd.Series(np.asarray(x, dtype=float), index=pd.to_datetime(np.asarray(t, dtype=float), unit="D"))


def rolling_step(t: np.ndarray, x: np.ndarray, w_d: float) -> np.ndarray:
    """step_i = median(x over (t_i, t_i + w]) − median(x over (t_i − w, t_i]).

    NaN within ``w_d`` of either end of the record.
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    win = pd.Timedelta(days=w_d)
    trail = _series(t, x).rolling(win, closed="right").median().to_numpy()
    rev = pd.Series(x[::-1], index=pd.to_datetime(t[-1] - t[::-1], unit="D"))
    lead_incl = rev.rolling(win, closed="right").median().to_numpy()[::-1]  # [t_i, t_i + w)
    lead = np.r_[lead_incl[1:], np.nan]  # starts at the next sample: (t_i, t_i + w]
    step = lead - trail
    step[(t - t[0] < w_d) | (t[-1] - t < w_d)] = np.nan
    return step


def detect_steps(
    t: np.ndarray,
    x: np.ndarray,
    w_d: float = JUMP_WINDOW_D,
    threshold: float = JUMP_THRESHOLD_MG_L,
    merge_d: float = JUMP_MERGE_D,
) -> pd.DataFrame:
    """Candidate level changes: one row per run of |step| >= threshold.

    The median-difference step is flat-topped (full height whenever more than
    half of each window lies on its own side), so the peak only marks the run;
    the jump time is located afterwards by :func:`characterize_jump`.
    """
    t = np.asarray(t, dtype=float)
    step = rolling_step(t, x, w_d)
    idx = np.flatnonzero(np.abs(np.nan_to_num(step)) >= threshold)
    rows = []
    if len(idx):
        for run in np.split(idx, np.flatnonzero(np.diff(t[idx]) > merge_d) + 1):
            k = run[np.argmax(np.abs(step[run]))]
            rows.append({"run_start_d": t[run[0]], "run_end_d": t[run[-1]], "t_step_peak_d": t[k], "step_peak": step[k]})
    return pd.DataFrame(rows, columns=["run_start_d", "run_end_d", "t_step_peak_d", "step_peak"])


def characterize_jump(t: np.ndarray, x: np.ndarray, t_ref: float) -> dict[str, float]:
    """Levels and 10/50/90 % crossing times of a jump located near ``t_ref``.

    Pre and post levels are medians 0.05 to 0.15 d before and after ``t_ref``.
    Crossings are taken on a 2-min rolling median within ±0.05 d: t90 is the
    first time the normalized signal reaches 0.9, t10 the last time at or below
    0.1 before t90, and t50 the first time at or above 0.5 after t10. NaN where
    a crossing does not exist.
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    pre_m = (t >= t_ref - LEVEL_OFFSET_D - LEVEL_WINDOW_D) & (t < t_ref - LEVEL_OFFSET_D)
    post_m = (t > t_ref + LEVEL_OFFSET_D) & (t <= t_ref + LEVEL_OFFSET_D + LEVEL_WINDOW_D)
    pre, post = float(np.median(x[pre_m])), float(np.median(x[post_m]))
    out = {"pre_level": pre, "post_level": post, "level_change": post - pre,
           "t10_d": np.nan, "t50_d": np.nan, "t90_d": np.nan, "transition_10_90_min": np.nan}
    win = (t >= t_ref - LEVEL_OFFSET_D) & (t <= t_ref + LEVEL_OFFSET_D)
    tw = t[win]
    smooth = _series(t, x).rolling(TRANSITION_SMOOTH, center=True).median().to_numpy()[win]
    y = (smooth - pre) / (post - pre)
    above90 = np.flatnonzero(y >= 0.9)
    if len(above90) == 0:
        return out
    i90 = above90[0]
    below10 = np.flatnonzero(y[: i90 + 1] <= 0.1)
    if len(below10) == 0:
        return out
    i10 = below10[-1]
    i50 = i10 + int(np.flatnonzero(y[i10 : i90 + 1] >= 0.5)[0])
    out.update(t10_d=tw[i10], t50_d=tw[i50], t90_d=tw[i90], transition_10_90_min=(tw[i90] - tw[i10]) * 1440.0)
    return out


def classify_jump(t_jump: float, level_change: float, pulse_times_d: list[float],
                  lag_min_d: float = PULSE_LAG_MIN_D, lag_max_d: float = PULSE_LAG_MAX_D) -> tuple[str, int | None, float]:
    """``pulse_response`` if the jump rises within the lag window after a pulse, else ``electrode_offset``."""
    lags = t_jump - np.asarray(pulse_times_d, dtype=float)
    ok = (lags >= lag_min_d) & (lags <= lag_max_d)
    if level_change > 0 and ok.any():
        k = int(np.flatnonzero(ok)[-1])
        return "pulse_response", k + 1, float(lags[k] * 24.0)
    return "electrode_offset", None, np.nan


def jump_table(t: np.ndarray, x: np.ndarray, pulse_times_d: list[float],
               threshold: float = JUMP_THRESHOLD_MG_L) -> pd.DataFrame:
    """All candidates with their characterization.

    ``sustained`` requires |post − pre| >= threshold with the levels taken 0.05 to
    0.15 d away, so a transient excursion is not a jump. ``t_jump_d`` is the
    50 % crossing; candidates without one keep NaN and are not classified.
    """
    cand = detect_steps(t, x, threshold=threshold)
    rows = []
    for c in cand.itertuples(index=False):
        ch = characterize_jump(t, x, c.t_step_peak_d)
        row = {**c._asdict(), **ch, "t_jump_d": ch["t50_d"], "sustained": abs(ch["level_change"]) >= threshold}
        if np.isfinite(ch["t50_d"]):
            row["class"], row["pulse"], row["lag_h"] = classify_jump(ch["t50_d"], ch["level_change"], pulse_times_d)
        else:
            row["class"], row["pulse"], row["lag_h"] = "unlocated", None, np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    if len(df):
        df["pulse"] = pd.array(df["pulse"], dtype="Int64")
        df["jump"] = df["sustained"] & (df["class"] != "unlocated")
    return df


def hf_residual(t: np.ndarray, x: np.ndarray) -> np.ndarray:
    """High-frequency residual: x minus its centred 30-min rolling median."""
    s = _series(t, x)
    return (s - s.rolling(HF_MEDIAN_WINDOW, center=True).median()).to_numpy()


def robust_sigma(r: np.ndarray) -> float:
    r = r[np.isfinite(r)]
    return float(1.4826 * np.median(np.abs(r - np.median(r)))) if len(r) else np.nan


def hf_blocks(t: np.ndarray, resid: np.ndarray, start: float, end: float, block_d: float = HF_BLOCK_D) -> pd.DataFrame:
    """Robust HF residual scale in consecutive blocks."""
    edges = np.arange(start, end + 1e-12, block_d)
    idx = np.searchsorted(edges, t, side="right") - 1
    rows = []
    for b in range(len(edges) - 1):
        r = resid[idx == b]
        rows.append({"block_start_d": edges[b], "block_end_d": edges[b + 1], "n": len(r), "sigma_hf": robust_sigma(r)})
    return pd.DataFrame(rows)


def held_segments(blocks: pd.DataFrame, threshold: float = HELD_SIGMA_MG_L, min_d: float = HELD_MIN_D) -> pd.DataFrame:
    """Runs of consecutive blocks with sigma_hf below ``threshold`` lasting at least ``min_d``."""
    low = (blocks["sigma_hf"] < threshold).to_numpy()
    rows = []
    idx = np.flatnonzero(low)
    if len(idx):
        for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            start, end = blocks["block_start_d"].iloc[run[0]], blocks["block_end_d"].iloc[run[-1]]
            if end - start >= min_d - 1e-9:
                rows.append({"t_start_d": start, "t_end_d": end, "duration_d": end - start,
                             "sigma_hf_max": blocks["sigma_hf"].iloc[run].max()})
    return pd.DataFrame(rows, columns=["t_start_d", "t_end_d", "duration_d", "sigma_hf_max"])


def ec_residual(t: np.ndarray, ec: np.ndarray) -> np.ndarray:
    s = _series(t, ec)
    return (s - s.rolling(EC_BASELINE_WINDOW, center=True).median()).to_numpy()


def ec_pulse_threshold(t: np.ndarray, resid: np.ndarray, start_d: float = EC_LATE_START_D) -> float:
    """Noise-floor rule: the 99.9th percentile of the residual before ``start_d``."""
    return float(np.nanquantile(resid[t < start_d], EC_THRESHOLD_QUANTILE))


def _decay(tt, base, amp, tau):
    return base + amp * np.exp(-tt / tau)


def detect_ec_pulses(t: np.ndarray, ec: np.ndarray, resid: np.ndarray, threshold: float,
                     start_d: float = EC_LATE_START_D) -> pd.DataFrame:
    """Late EC pulses: residual peaks (height and prominence >= threshold) after ``start_d``,
    each with an exponential decay fitted over the following 2 h (or until the next peak)."""
    dt_min = float(np.median(np.diff(t))) * 1440.0
    r = np.where((t >= start_d) & np.isfinite(resid), resid, 0.0)
    peaks, _ = find_peaks(r, height=threshold, prominence=threshold,
                          distance=max(1, round(EC_MIN_SEPARATION_MIN / dt_min)))
    rows = []
    for n, p in enumerate(peaks):
        t_end = t[p] + EC_FIT_WINDOW_D
        if n + 1 < len(peaks):
            t_end = min(t_end, t[peaks[n + 1]])
        m = (t >= t[p]) & (t < t_end)
        tt = (t[m] - t[p]) * 1440.0
        y = ec[m]
        fit = {"tau_min": np.nan, "fit_amp_mS_cm": np.nan, "fit_base_mS_cm": np.nan, "fit_ok": False}
        if len(y) >= 5:
            try:
                popt, _ = curve_fit(_decay, tt, y, p0=(y[-1], max(y[0] - y[-1], 1e-4), 20.0),
                                    bounds=([-np.inf, 0.0, 0.5], [np.inf, np.inf, 600.0]), maxfev=10000)
                fit = {"tau_min": popt[2], "fit_amp_mS_cm": popt[1], "fit_base_mS_cm": popt[0], "fit_ok": True}
            except RuntimeError:
                pass  # recorded as fit_ok = False with NaN parameters
        rows.append({"t_peak_d": t[p], "amplitude_mS_cm": resid[p], "n_fit_points": len(y), **fit})
    return pd.DataFrame(rows, columns=["t_peak_d", "amplitude_mS_cm", "n_fit_points", "tau_min",
                                       "fit_amp_mS_cm", "fit_base_mS_cm", "fit_ok"])


def ec_pulse_sensitivity(t: np.ndarray, resid: np.ndarray, thresholds, start_d: float = EC_LATE_START_D) -> pd.DataFrame:
    dt_min = float(np.median(np.diff(t))) * 1440.0
    r = np.where((t >= start_d) & np.isfinite(resid), resid, 0.0)
    rows = []
    for h in thresholds:
        peaks, _ = find_peaks(r, height=h, prominence=h, distance=max(1, round(EC_MIN_SEPARATION_MIN / dt_min)))
        rows.append({"threshold_mS_cm": h, "n_pulses": len(peaks)})
    return pd.DataFrame(rows)


def step_time_in_window(t: np.ndarray, x: np.ndarray, lo: float, hi: float, w_d: float = 0.01) -> tuple[float, float]:
    """Time and size of the largest rolling-median step between ``lo`` and ``hi``."""
    m = (t >= lo - 2 * w_d) & (t <= hi + 2 * w_d)
    tt, step = t[m], rolling_step(t[m], x[m], w_d)
    inside = (tt >= lo) & (tt <= hi) & np.isfinite(step)
    k = np.flatnonzero(inside)[np.argmax(np.abs(step[inside]))]
    return float(tt[k]), float(step[k])


def native_interval_s(t: np.ndarray) -> float:
    return float(np.median(np.diff(t))) * SECONDS_PER_DAY
