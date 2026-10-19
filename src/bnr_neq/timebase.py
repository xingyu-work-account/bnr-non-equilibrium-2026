"""Stamp reconstruction, gap inventory and the common 0.01-day grid."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .constants import SECONDS_PER_DAY

STAMP_D = 0.01  # d, timestamp resolution of NH4.csv and ORP.csv
STAMP_TOL_D = 1e-9  # d, tolerance for comparing parsed decimal stamps

GRID_START_D = 0.0
GRID_END_D = 47.92  # four-sensor common window
GRID_STEP_D = 0.01

LOCAL_WINDOW = 101  # intervals in the rolling median that defines the native interval
GAP_FACTOR = 2.0  # interval > GAP_FACTOR x native interval is a gap
MISSING_FACTOR = 1.5  # interval > MISSING_FACTOR x native interval is listed as a missing sample
INTERVAL_RTOL = 1e-4  # relative tolerance for timestamp rounding (1e-9 d on 10 s is 1e-5)

STAMP_ASSUMPTION = (
    "NH4 and ORP timestamps are rounded to 0.01 d; rows sharing a stamp are assumed to be "
    "uniformly spaced within the stamp: t = stamp + 0.01 x j/n (j = zero-based row index "
    "within the stamp, n = rows sharing it)."
)


class StampError(ValueError):
    """Stamp spacing or reconstructed times violate the stated assumptions."""


def check_stamp_spacing(t_stamp: np.ndarray, label: str) -> None:
    """Distinct stamps must be exactly 0.01 d apart."""
    u = np.unique(t_stamp)
    d = np.diff(u)
    bad = np.abs(d - STAMP_D) > STAMP_TOL_D
    if bad.any():
        i = int(np.argmax(bad))
        raise StampError(
            f"{label}: distinct stamps {u[i]} and {u[i + 1]} differ by {d[i]} d, not {STAMP_D} d "
            f"({int(bad.sum())} such pairs)"
        )


def stamp_groups(t_stamp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Start index and size of each run of equal stamps (input must be non-decreasing)."""
    starts = np.r_[0, np.flatnonzero(np.diff(t_stamp) > 0) + 1]
    sizes = np.diff(np.r_[starts, len(t_stamp)])
    return starts, sizes


def reconstruct_stamps(t_stamp: np.ndarray, label: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return reconstructed times, within-stamp index j and stamp size n."""
    t_stamp = np.asarray(t_stamp, dtype=float)
    check_stamp_spacing(t_stamp, label)
    starts, sizes = stamp_groups(t_stamp)
    n = np.repeat(sizes, sizes)
    j = np.arange(len(t_stamp)) - np.repeat(starts, sizes)
    t = t_stamp + STAMP_D * j / n
    if not (np.diff(t) > 0).all():
        i = int(np.argmax(np.diff(t) <= 0))
        raise StampError(f"{label}: reconstructed times not strictly increasing at row {i}")
    return t, j, n


def local_interval_d(t: np.ndarray, window: int = LOCAL_WINDOW) -> np.ndarray:
    """Native sampling interval around each interval: centred rolling median of ``diff(t)``."""
    dt = np.diff(np.asarray(t, dtype=float))
    return pd.Series(dt).rolling(window, center=True, min_periods=1).median().to_numpy()


def classify_intervals(t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (dt, local, is_gap, is_missing) per interval, in days."""
    dt = np.diff(np.asarray(t, dtype=float))
    local = local_interval_d(t)
    is_gap = dt > GAP_FACTOR * local * (1 + INTERVAL_RTOL)
    is_missing = dt > MISSING_FACTOR * local * (1 + INTERVAL_RTOL)
    return dt, local, is_gap, is_missing


def gap_inventory(t: np.ndarray, sensor: str, t_stamp: np.ndarray | None = None) -> pd.DataFrame:
    """Runs of intervals longer than 1.5x the local native interval; ``gap`` if any exceeds 2x.

    Contiguous flagged intervals are merged into one row. For stamped sensors,
    ``t_stamp`` marks runs lying wholly inside the first or last stamp, where the
    uniform-spacing assumption stretches a partial stamp over the full 0.01 d
    (``partial_stamp_artifact``).
    """
    t = np.asarray(t, dtype=float)
    dt, local, is_gap, is_missing = classify_intervals(t)
    idx = np.flatnonzero(is_missing)
    columns = ["sensor", "t_start_d", "t_end_d", "n_intervals", "interval_max_s",
               "native_interval_s", "ratio_max", "class", "partial_stamp_artifact"]
    if len(idx) == 0:
        return pd.DataFrame(columns=columns)
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    if t_stamp is not None:
        first_hi = t_stamp[0] + STAMP_D
        last_lo = t_stamp[-1]
    rows = []
    for r in runs:
        start, end = t[r[0]], t[r[-1] + 1]
        artifact = False
        if t_stamp is not None:
            artifact = bool(end <= first_hi + STAMP_TOL_D or start >= last_lo - STAMP_TOL_D)
        rows.append(
            {
                "sensor": sensor,
                "t_start_d": start,
                "t_end_d": end,
                "n_intervals": len(r),
                "interval_max_s": dt[r].max() * SECONDS_PER_DAY,
                "native_interval_s": float(np.median(local[r])) * SECONDS_PER_DAY,
                "ratio_max": float((dt[r] / local[r]).max()),
                "class": "gap" if is_gap[r].any() else "missing_sample",
                "partial_stamp_artifact": artifact,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def sample_weights_d(t: np.ndarray) -> np.ndarray:
    """Weight of each sample = its own (forward) sampling interval, in days.

    An interval classed as a gap is replaced by the local native interval, so a
    sample before a gap does not stand in for the missing period. The last
    sample takes the last local native interval.
    """
    dt, local, is_gap, _ = classify_intervals(t)
    w = np.where(is_gap, local, dt)
    return np.r_[w, local[-1]]


def grid_edges(start: float = GRID_START_D, end: float = GRID_END_D, step: float = GRID_STEP_D) -> np.ndarray:
    n = round((end - start) / step)
    return np.round(start + step * np.arange(n + 1), 10)


def bin_index(t: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Bin of each time (left-closed bins, last bin closed); -1 outside the grid."""
    nbins = len(edges) - 1
    idx = np.searchsorted(edges, t, side="right") - 1
    idx[t == edges[-1]] = nbins - 1
    idx[(idx < 0) | (idx >= nbins)] = -1
    return idx


def bin_time_weighted(
    t: np.ndarray, x: np.ndarray, w: np.ndarray, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Time-weighted bin means and sample counts; empty bins are NaN, never filled."""
    nbins = len(edges) - 1
    idx = bin_index(np.asarray(t, dtype=float), edges)
    ok = idx >= 0
    count = np.bincount(idx[ok], minlength=nbins)
    sw = np.bincount(idx[ok], weights=w[ok], minlength=nbins)
    swx = np.bincount(idx[ok], weights=w[ok] * x[ok], minlength=nbins)
    mean = np.full(nbins, np.nan)
    np.divide(swx, sw, out=mean, where=count > 0)
    return mean, count


def common_grid(natives: dict[str, pd.DataFrame], variables: dict[str, list[str]]) -> pd.DataFrame:
    """Build the 0.01-day grid over [0, 47.92].

    ``natives`` maps sensor to its native frame (column ``t_d``), ``variables``
    maps sensor to the columns to average. For each sensor the grid holds the
    sample count ``n_<sensor>``, the empty-bin mask ``empty_<sensor>`` and
    ``gap_<sensor>`` (a gap interval starts in the bin).
    """
    edges = grid_edges()
    grid = pd.DataFrame({"t_start_d": edges[:-1], "t_end_d": edges[1:]})
    grid["t_mid_d"] = 0.5 * (grid["t_start_d"] + grid["t_end_d"])
    for sensor, df in natives.items():
        t = df["t_d"].to_numpy()
        w = sample_weights_d(t)
        count = None
        for col in variables[sensor]:
            mean, count = bin_time_weighted(t, df[col].to_numpy(dtype=float), w, edges)
            grid[col] = mean
        grid[f"n_{sensor}"] = count
        grid[f"empty_{sensor}"] = count == 0
        _, _, is_gap, _ = classify_intervals(t)
        gap_bins = bin_index(t[:-1][is_gap], edges)
        flag = np.zeros(len(grid), dtype=bool)
        flag[gap_bins[gap_bins >= 0]] = True
        grid[f"gap_{sensor}"] = flag
    return grid
