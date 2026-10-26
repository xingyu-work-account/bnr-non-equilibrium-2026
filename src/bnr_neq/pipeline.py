"""Pipeline stages: each reads raw files or earlier outputs and writes to ``outputs/``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import io, timebase
from .config import PROJECT_ROOT, Config
from .constants import SECONDS_PER_DAY
from .report import md_table

OUTPUT_DIR = PROJECT_ROOT / "outputs"

GRID_VARIABLES = {
    "cond": ["I_A", "EC_mS_cm"],
    "nh4": ["nh4_mgNH4_L", "nh4_mgN_L"],
    "orp": ["orp_mV_vs_Ag"],
    "ph": ["pH"],
}


def out_dirs(out_dir: Path | str) -> dict[str, Path]:
    root = Path(out_dir)
    dirs = {k: root / k for k in ("processed", "tables", "figures", "reports")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


@dataclass
class Ingested:
    loads: dict[str, io.RawLoad]
    natives: dict[str, pd.DataFrame]
    dropped: pd.DataFrame
    gaps: pd.DataFrame
    grid: pd.DataFrame
    properties: pd.DataFrame


def build_native(load: io.RawLoad) -> pd.DataFrame:
    """Native series with reconstructed time for stamped sensors."""
    df = load.data.copy()
    if load.spec.stamped:
        t, j, n = timebase.reconstruct_stamps(df["t_d"].to_numpy(), load.spec.filename)
        df = df.rename(columns={"t_d": "t_stamp_d"})
        df.insert(1, "t_d", t)
        df.insert(2, "stamp_j", j)
        df.insert(3, "stamp_n", n)
    return df


def data_properties(loads: dict[str, io.RawLoad], natives: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Measured file properties, for comparison with implementation_prompt.md, Section 4."""
    rows = []
    for key, load in loads.items():
        df = natives[key]
        value_col = load.spec.columns[1]
        t_raw = df["t_stamp_d"] if load.spec.stamped else df["t_d"]
        dt_s = np.diff(df["t_d"].to_numpy()) * SECONDS_PER_DAY
        row = {
            "sensor": key,
            "file": load.spec.filename,
            "header_found": ",".join(load.spec.header),
            "header_section4": ",".join(io.SECTION4_HEADERS[key]),
            "header_matches_section4": load.spec.header == io.SECTION4_HEADERS[key],
            "rows_read": load.n_rows_read,
            "rows_valid": len(df),
            "rows_dropped": len(load.dropped),
            "t_min_d": float(t_raw.min()),
            "t_max_d": float(t_raw.max()),
            "value_column": value_col,
            "value_min": float(df[value_col].min()),
            "value_max": float(df[value_col].max()),
            "interval_median_s": float(np.median(dt_s)),
            "interval_min_s": float(dt_s.min()),
            "interval_max_s": float(dt_s.max()),
            "stamp_rows_min": np.nan,
            "stamp_rows_max": np.nan,
            "stamp_rows_interior_min": np.nan,
            "n_stamps": np.nan,
            "EC_min_mS_cm": np.nan,
            "EC_max_mS_cm": np.nan,
            "tie_fraction": float(np.mean(np.diff(df[value_col].to_numpy()) == 0)),
        }
        if load.spec.stamped:
            _, sizes = timebase.stamp_groups(df["t_stamp_d"].to_numpy())
            row.update(
                stamp_rows_min=int(sizes.min()),
                stamp_rows_max=int(sizes.max()),
                stamp_rows_interior_min=int(sizes[1:-1].min()),
                n_stamps=len(sizes),
            )
        if key == "cond":
            row.update(EC_min_mS_cm=float(df["EC_mS_cm"].min()), EC_max_mS_cm=float(df["EC_mS_cm"].max()))
        rows.append(row)
    return pd.DataFrame(rows)


def nh4_rate_regimes(nh4: pd.DataFrame) -> pd.DataFrame:
    """Median rows per stamp and implied interval in the three NH4 logging regimes."""
    sizes = nh4.groupby("t_stamp_d").size()
    rows = []
    for lo, hi in [(0.0, 41.05), (41.05, 45.0), (45.0, 47.92)]:
        s = sizes[(sizes.index > lo + 0.01) & (sizes.index < hi - 0.01)]
        med = float(np.median(s))
        rows.append(
            {"from_d": lo, "to_d": hi, "median_rows_per_stamp": med,
             "implied_interval_s": timebase.STAMP_D * SECONDS_PER_DAY / med}
        )
    return pd.DataFrame(rows)


def ingest(cfg: Config, data_dir: Path | str = io.DATA_DIR) -> Ingested:
    loads = {key: io.load_sensor(key, cfg, data_dir) for key in io.SENSORS}
    natives = {key: build_native(load) for key, load in loads.items()}
    dropped = pd.concat([load.dropped for load in loads.values()], ignore_index=True)
    gaps = pd.concat(
        [
            timebase.gap_inventory(
                df["t_d"].to_numpy(), key, df["t_stamp_d"].to_numpy() if "t_stamp_d" in df else None
            )
            for key, df in natives.items()
        ],
        ignore_index=True,
    )
    grid = timebase.common_grid(natives, GRID_VARIABLES)
    properties = data_properties(loads, natives)
    return Ingested(loads, natives, dropped, gaps, grid, properties)


def grid_coverage(ing: Ingested) -> pd.DataFrame:
    edges = timebase.grid_edges()
    rows = []
    for key, df in ing.natives.items():
        t = df["t_d"].to_numpy()
        idx = timebase.bin_index(t, edges)
        rows.append(
            {
                "sensor": key,
                "samples_in_grid": int((idx >= 0).sum()),
                "samples_outside_grid": int((idx < 0).sum()),
                "empty_bins": int(ing.grid[f"empty_{key}"].sum()),
                "bins_with_gap": int(ing.grid[f"gap_{key}"].sum()),
                "min_samples_per_nonempty_bin": int(ing.grid.loc[~ing.grid[f"empty_{key}"], f"n_{key}"].min()),
            }
        )
    return pd.DataFrame(rows)


def write_ingest(ing: Ingested, out_dir: Path | str = OUTPUT_DIR) -> dict[str, Path]:
    d = out_dirs(out_dir)
    for key, df in ing.natives.items():
        df.to_parquet(d["processed"] / f"native_{key}.parquet", index=False)
    ing.grid.to_parquet(d["processed"] / "grid_0p01d.parquet", index=False)
    ing.dropped.to_csv(d["tables"] / "dropped_rows.csv", index=False)
    ing.gaps.to_csv(d["tables"] / "gaps.csv", index=False)
    ing.properties.to_csv(d["tables"] / "data_properties.csv", index=False)
    coverage = grid_coverage(ing)
    coverage.to_csv(d["tables"] / "grid_coverage.csv", index=False)
    regimes = nh4_rate_regimes(ing.natives["nh4"])

    props = ing.properties
    first_stamp = {
        key: int(ing.natives[key]["stamp_n"].iloc[0]) for key in ("nh4", "orp")
    }
    lines = [
        "# Ingestion and time base (M1)",
        "",
        "Generated by `python -m bnr_neq ingest`. Do not edit by hand.",
        "",
        "## Dropped rows",
        "",
        f"Total dropped: {len(ing.dropped)}. Each is a row whose time or value is not a finite number.",
        "",
        md_table(ing.dropped) if len(ing.dropped) else "None.",
        "",
        "## File properties",
        "",
        md_table(props[["sensor", "file", "header_found", "header_section4", "header_matches_section4",
                        "rows_read", "rows_valid", "rows_dropped", "t_min_d", "t_max_d"]], digits=8),
        "",
        md_table(props[["sensor", "value_column", "value_min", "value_max", "EC_min_mS_cm", "EC_max_mS_cm",
                        "interval_median_s", "interval_min_s", "interval_max_s", "tie_fraction"]], digits=6),
        "",
        md_table(props.loc[props["sensor"].isin(["nh4", "orp"]),
                           ["sensor", "n_stamps", "stamp_rows_min", "stamp_rows_max", "stamp_rows_interior_min"]]),
        "",
        "NH4 logging regimes (median rows per 0.01 d stamp, interior stamps only):",
        "",
        md_table(regimes),
        "",
        "Interval statistics for NH4 and ORP are computed on reconstructed times.",
        "",
        "## Gap inventory",
        "",
        ("The native interval is the centred rolling median of 101 intervals, because the NH4 logging "
        f"rate changes from 20 s to 2 s. Intervals longer than {timebase.GAP_FACTOR}x the native interval "
        f"are classed `gap`; intervals longer than {timebase.MISSING_FACTOR}x are listed as `missing_sample`. "
        f"A relative tolerance of {timebase.INTERVAL_RTOL} absorbs timestamp rounding. Contiguous flagged "
        "intervals are merged into one row. `partial_stamp_artifact` marks runs inside the first or last "
        "NH4/ORP stamp: a partial stamp (for example 22 rows at 20 s, i.e. 0.005 d of data) is stretched over "
        "the full 0.01 d by the prescribed reconstruction, so its intervals look about twice the native interval. "
        "These are consequences of the stated assumption, not missing data."),
        "",
        md_table(ing.gaps, digits=8) if len(ing.gaps) else "No interval exceeds 1.5x the native interval.",
        "",
        "## Common grid",
        "",
        ("Bins of 0.01 d over [0, 47.92] (left-closed, last bin closed). Each bin holds the time-weighted "
        "mean of its native samples, each weighted by its own forward sampling interval (an interval "
        "classed as a gap is replaced by the native interval). Empty bins are NaN with count 0; nothing "
        "is interpolated."),
        "",
        md_table(coverage),
        "",
        "## Stated assumptions",
        "",
        f"1. {timebase.STAMP_ASSUMPTION}",
        ("2. The conductivity record is converted with the configured equation and is flagged as "
        "unvalidated as solution conductivity until M2 and the probe checks decide otherwise."),
        "3. NH4 values are in mg NH4+/L as confirmed by the experimenter (the file header says only `mg/L`).",
        "",
        "## Open observations",
        "",
        (f"- The first NH4 stamp holds {first_stamp['nh4']} rows and the first ORP stamp {first_stamp['orp']}, "
        "about half the interior stamp sizes. This is consistent with stamps rounded to the nearest 0.01 d "
        "(stamp 0 can only collect samples from [0, 0.005) d). If so, the prescribed reconstruction places "
        "every NH4 and ORP time 0.005 d (7.2 min) late. Not corrected; the prescribed formula is used."),
        "- Headers differ from Section 4 for three files (see table). The loader maps the headers present in the files.",
    ]
    report = d["reports"] / "ingest.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return d
