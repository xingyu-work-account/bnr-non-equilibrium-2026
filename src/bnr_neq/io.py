"""Raw file loading and unit conversion.

Header maps are explicit. The headers below are the ones present in ``data/``
(verified in M0); they differ in spelling from implementation_prompt.md,
Section 4, for three of the four files (deviation D1, see
``docs/progress.md``). Any other header raises an error naming the
file.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PROJECT_ROOT, Config

DATA_DIR = PROJECT_ROOT / "data"

EC_VALID_RANGE_MS_CM = (0.0, 5.0)


@dataclass(frozen=True)
class SensorSpec:
    key: str
    filename: str
    header: tuple[str, str]  # exact raw column names
    columns: tuple[str, str]  # mapped names
    stamped: bool  # timestamps rounded to 0.01 d, reconstructed in timebase


SENSORS: dict[str, SensorSpec] = {
    "cond": SensorSpec("cond", "Cond.csv", ("Time(d)", "Conductivity(A)"), ("t_d", "I_A"), False),
    "nh4": SensorSpec(
        "nh4", "NH4.csv", ("Time(d)", "Concentration(mg/L)"), ("t_d", "nh4_mgNH4_L"), True
    ),
    "orp": SensorSpec("orp", "ORP.csv", ("Time(d)", "ORP(mV)"), ("t_d", "orp_mV_vs_Ag"), True),
    "ph": SensorSpec("ph", "pH.csv", ("Time(d)", "pH"), ("t_d", "pH"), False),
}

# Headers as written in implementation_prompt.md, Section 4, kept for the
# discrepancy report only. They are not accepted by the loader.
SECTION4_HEADERS: dict[str, tuple[str, str]] = {
    "cond": ("Time, d", "Cond. A"),
    "nh4": ("Time(d)", "Concentration"),
    "orp": ("Time", "ORP (mV)"),
    "ph": ("Time(d)", "pH"),
}


class HeaderError(ValueError):
    """A file does not carry the expected header."""


class MonotonicityError(ValueError):
    """Time decreases within a file."""


@dataclass
class RawLoad:
    spec: SensorSpec
    path: Path
    n_rows_read: int
    data: pd.DataFrame  # valid numeric rows, mapped column names
    dropped: pd.DataFrame  # file, line, raw values of every dropped row


def read_raw(spec: SensorSpec, data_dir: Path | str = DATA_DIR) -> RawLoad:
    """Read one file, check its header, coerce to numbers and drop non-numeric rows."""
    path = Path(data_dir) / spec.filename
    if not path.is_file():
        raise FileNotFoundError(f"{spec.filename}: input file not found at {path}")
    raw = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    absent = [h for h in spec.header if h not in raw.columns]
    if absent or len(raw.columns) != len(spec.header):
        raise HeaderError(
            f"{spec.filename}: expected header {list(spec.header)}, "
            f"found {list(raw.columns)} (absent: {absent})"
        )
    numeric = pd.DataFrame(
        {new: pd.to_numeric(raw[old], errors="coerce") for old, new in zip(spec.header, spec.columns)}
    )
    valid = np.isfinite(numeric.to_numpy(dtype=float)).all(axis=1)
    dropped = raw.loc[~valid, list(spec.header)].copy()
    dropped.columns = ["raw_time", "raw_value"]
    dropped.insert(0, "line", dropped.index.to_numpy() + 2)  # 1-based, header is line 1
    dropped.insert(0, "file", spec.filename)
    data = numeric.loc[valid].reset_index(drop=True)
    check_non_decreasing(data["t_d"].to_numpy(), spec.filename)
    return RawLoad(spec, path, len(raw), data, dropped.reset_index(drop=True))


def check_non_decreasing(t: np.ndarray, label: str) -> None:
    d = np.diff(t)
    if (d < 0).any():
        i = int(np.argmax(d < 0))
        raise MonotonicityError(f"{label}: time decreases after valid row {i} ({t[i]} -> {t[i + 1]})")


def convert_ec(current_A: np.ndarray, slope_mS_cm_per_A: float, intercept_mS_cm: float) -> np.ndarray:
    """EC (mS/cm) = slope x I (A) + intercept; raise if any value is outside 0 to 5 mS/cm.

    A value outside that range indicates a unit error, for example a current in mA.
    """
    ec = slope_mS_cm_per_A * np.asarray(current_A, dtype=float) + intercept_mS_cm
    lo, hi = EC_VALID_RANGE_MS_CM
    bad = (ec < lo) | (ec > hi) | ~np.isfinite(ec)
    if bad.any():
        raise ValueError(
            f"EC conversion gives {int(bad.sum())} value(s) outside {lo} to {hi} mS/cm "
            f"(first: I = {np.asarray(current_A)[bad][0]} -> {ec[bad][0]} mS/cm); "
            "is the current in A?"
        )
    return ec


def nh4_to_n(nh4_mgNH4_L: np.ndarray, factor: float) -> np.ndarray:
    """mg NH4+/L to mg N/L."""
    return np.asarray(nh4_mgNH4_L, dtype=float) * factor


def load_sensor(key: str, cfg: Config, data_dir: Path | str = DATA_DIR) -> RawLoad:
    """Read one sensor and apply its unit conversion from the configuration."""
    load = read_raw(SENSORS[key], data_dir)
    df = load.data
    if key == "cond":
        df["EC_mS_cm"] = convert_ec(
            df["I_A"].to_numpy(),
            cfg.get("sensors.ec_slope_mS_cm_per_A"),
            cfg.get("sensors.ec_intercept_mS_cm"),
        )
    elif key == "nh4":
        df["nh4_mgN_L"] = nh4_to_n(df["nh4_mgNH4_L"].to_numpy(), cfg.get("sensors.nh4_to_N_factor"))
    return load
