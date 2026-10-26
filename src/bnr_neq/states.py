"""Process states from relative ORP, pH and the inverted NH4 (W3).

Every 0.01 d grid bin from day 1 is assigned one of

- ``active_nitrification``: the inverted NH4 C(t) (W1) falls and pH does not rise (acid production possible);
- ``inactive_oversupplied``: NH4 is available but not consumed while the redox level is high or pH rises (O2 present,
  CO2 stripped without acid production): stirring or aeration spent without conversion;
- ``oxygen_limited``: NH4 is available but not consumed while the redox level is low and neither ORP nor pH rises;
- ``sensor_artifact``: the bin lies in a pulse transient, an NH4 step region, held/damped NH4 output, an ORP electrode
  event or an empty bin; such bins carry **no process state**;
- ``unclassified``: a slope window crosses an artifact or event, or the indicators conflict or stay inside thresholds.

The rules, windows and thresholds below were fixed from sensor resolution and physical scales before any state was
computed or compared (``RULES``, ``rules_table``), but two parts were revised after the first results (the ≥ 2 blocks
rule and the autocorrelation-widened intervals, ``REGISTRATION_STATUS``); they are frozen by ``registration_sha256``
from pass 6 and are not pre-registered. A slope is rising, falling or
flat by its Theil-Sen 95 % interval over the centred 1.01 d window, widened for autocorrelation (pass 3, X2): each
half-width is multiplied by √(n/n_eff), n_eff = n (1 − ρ)/(1 + ρ), ρ the lag-1 autocorrelation of the window's
residuals about its Theil-Sen line (capped at 0.999), the formula of ``model.effective_n``; the same for C(t), pH and
ORP. Only a minority of the record carries a state, so the idle fraction is reported as bounds (inactive over all bins;
inactive, unclassified and artifact over all bins) and, CONDITIONAL on the classified bins being representative, over
the classified bins only when they cover at least half of the interval. The classification uses the NH4
slope, so comparing NH4 consumption between full states is circular. The association test therefore partitions the
bins by ORP and pH alone (``INDICATOR_PARTITIONS``) and compares the NH4 consumption between those groups with
moving-block bootstrap intervals whose block length covers the slope window and the autocorrelation of the response.
Associations are not causes: time, pulse number and electrode state are confounded with every indicator.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

from . import model
from . import segments as S
from .config import Config
from .constants import F_C_PER_MOL, R_J_PER_MOL_K, T_REF_K

T0_D = 1.0  # start-up excluded, as in M3 and M6
END_D = 47.92
BIN_D = 0.01
SLOPE_HALF_BINS = 50  # centred 1.01 d window: one pH quantum (0.01) per day is resolvable over it
PH_QUANTUM = 0.01  # pH logger resolution (tie fractions 0.82-0.96 in M1)
# Relative ORP threshold: one decade of p_O2 moves a Nernstian 4-electron O2/H2O couple by RT ln10/(4F) ≈ 14.8 mV at
# 25 °C (SCENARIO temperature). Rounded to 15 mV. Pt in water is not strictly Nernstian; this only fixes the scale.
ORP_LEVEL_MV = round(1e3 * R_J_PER_MOL_K * T_REF_K * math.log(10) / (4 * F_C_PER_MOL))
ORP_THRESHOLD_SCENARIO = "SCENARIO 25 °C: RT ln10/(4F) at 298.15 K = 14.8 mV, rounded; the liquid temperature is not measured"
ORP_REF_HALF_BINS = 250  # centred 5 d reference median: relative level, insensitive to slow reference drift
ORP_REF_MIN_BINS = 100
ORP_EVENTS_D = (28.89, 30.48)  # dynamics.ORP_ELECTRODE_EVENTS_D (implementation_prompt.md M5 item 3)
ORP_EVENT_MARGIN_D = (0.05, 0.1)
CI_ALPHA = 0.95
RHO_MAX = 0.999  # cap of the lag-1 autocorrelation used to widen the slope intervals (n_eff stays > 0)
# X2 coverage: the idle fraction over classified time stands for the interval only if the classified bins are the
# majority of it; below one half, the unknown mix of the bins without a state outweighs them (fixed before the run)
MIN_CLASSIFIED_FRACTION = 0.5
N_BOOT = 999
MIN_BLOCKS_PER_GROUP = 2  # a block bootstrap needs each group to span at least two blocks' worth of bins
SEED = 20261019

STATES = ("active_nitrification", "inactive_oversupplied", "oxygen_limited")
ARTIFACT, UNCLASSIFIED = "sensor_artifact", "unclassified"
# Y3: synthetic check of the slope interval (AR(1) errors, 101 points), fixed before the run; the grid spans the
# record's median ρ̂ of C(t) (0.93 > the mean ρ̂ 0.88 at φ = 0.95), and φ = 0.9 is the reference of the pass-3 check
COVERAGE_PHIS = (0.5, 0.9, 0.95, 0.98)
COVERAGE_REFERENCE_PHI = 0.9
COVERAGE_WINDOWS = 600
COVERAGE_SLOPE, COVERAGE_NOISE, COVERAGE_BURN_IN = -2.0, 0.05, 300


class StateRuleError(ValueError):
    """A state rule uses a sensor feature that the sensor validation excluded."""


# --- Registered rules (revised after first results; see REGISTRATION_STATUS) -----------------------------------
# Features: nh4_inverted = slope of the W1 C(t) (the observation model absorbs the electrode under-response);
# ph = pH slope (pH PASSED M2); orp_relative = ORP minus its centred 5 d median and the ORP slope (the M2 failure
# concerns the bare-Ag reference potential: absolute E_H is excluded, relative changes over days are kept as an
# ASSUMPTION of a slowly varying reference). Raw NH4 readings, absolute ORP and conductivity are not used.
@dataclass(frozen=True)
class Rule:
    state: str
    condition: str
    features: tuple[str, ...]


RULES = (
    Rule("active_nitrification", "NH4 slope 95 % interval < 0 and pH not rising", ("nh4_inverted", "ph")),
    Rule("oxygen_limited", f"NH4 slope interval spans 0, ORP_rel ≤ −{ORP_LEVEL_MV} mV, ORP not rising, pH not rising",
         ("nh4_inverted", "orp_relative", "ph")),
    Rule("inactive_oversupplied", f"NH4 slope interval spans 0 and (ORP_rel ≥ +{ORP_LEVEL_MV} mV or pH rising), "
         f"ORP_rel not ≤ −{ORP_LEVEL_MV} mV", ("nh4_inverted", "orp_relative", "ph")),
)
# A FAIL verdict of a sensor excludes these features (key = sensor key of sensor_verdicts.csv).
EXCLUDED_BY_FAIL = {"nh4": ("nh4_reading",), "ec": ("ec",), "orp": ("orp_absolute",), "ph": ("ph",)}
# NH4-free partitions for the association test (the NH4 consumption is the response).
INDICATOR_PARTITIONS = {
    "ORP_rel low vs high": ("orp_low", "orp_high"),
    "pH rising vs not rising": ("ph_rising", "ph_not_rising"),
}


# --- Registration (pass 6, G1) ----------------------------------------------------------------------------------
# The rules were fixed before the first states were computed, but they landed in git together with their first
# results, and two parts changed after results existed: the ≥ 2 blocks rule (pass 2, after the first comparison) and
# the autocorrelation-widened slope intervals (pass 3, X2). They are therefore not pre-registered. From pass 6 every
# parameter and decision function of the classifier and of the association test is frozen by one SHA-256, recorded in
# REGISTRATION_DOC and checked by a test; any change to them changes the hash and must be re-registered there.
REGISTRATION_DOC = "docs/states_registration.md"
REGISTRATION_STATUS = ("REVISED AFTER FIRST RESULTS (not pre-registered): the >= 2 blocks rule was added after the "
                       "first comparison (pass 2) and the slope intervals were widened for autocorrelation after the "
                       "first states (pass 3, X2); rules frozen by SHA-256 from pass 6 (" + REGISTRATION_DOC + ")")
_REGISTERED_FUNCTIONS = ("artifact_reasons", "lag1_rho", "interval_scale", "rolling_theil_sen", "relative_orp",
                         "window_median", "sign_class", "features", "classify", "acf_block_length",
                         "block_bootstrap_diff", "indicator_groups", "comparisons",
                         # G2: the null check of the bootstrap, frozen before it was run on the record
                         "ar1", "null_response", "bootstrap_null_coverage", "null_coverage", "interval_reading")


def registration() -> dict:
    """Canonical form of everything that decides a state or the association test (parameters and function sources)."""
    import inspect
    import sys

    mod = sys.modules[__name__]
    return {
        "rules": [[r.state, r.condition, list(r.features)] for r in RULES],
        "parameters": {"T0_D": T0_D, "END_D": END_D, "BIN_D": BIN_D, "SLOPE_HALF_BINS": SLOPE_HALF_BINS,
                       "PH_QUANTUM": PH_QUANTUM, "ORP_LEVEL_MV": ORP_LEVEL_MV, "ORP_REF_HALF_BINS": ORP_REF_HALF_BINS,
                       "ORP_REF_MIN_BINS": ORP_REF_MIN_BINS, "ORP_EVENTS_D": list(ORP_EVENTS_D),
                       "ORP_EVENT_MARGIN_D": list(ORP_EVENT_MARGIN_D), "CI_ALPHA": CI_ALPHA, "RHO_MAX": RHO_MAX,
                       "N_BOOT": N_BOOT, "MIN_BLOCKS_PER_GROUP": MIN_BLOCKS_PER_GROUP, "SEED": SEED,
                       "EXCLUDE_AFTER_BOUNDARY_D": S.EXCLUDE_AFTER_BOUNDARY_D,
                       "NULL_COVERAGE_PHIS": list(NULL_COVERAGE_PHIS),
                       "NULL_COVERAGE_REPLICATES": NULL_COVERAGE_REPLICATES,
                       "NULL_COVERAGE_MIN_OK": NULL_COVERAGE_MIN_OK, "COVERAGE_BURN_IN": COVERAGE_BURN_IN},
        "excluded_by_fail": {k: list(v) for k, v in EXCLUDED_BY_FAIL.items()},
        "indicator_partitions": {k: list(v) for k, v in INDICATOR_PARTITIONS.items()},
        "functions": {name: inspect.getsource(getattr(mod, name)) for name in _REGISTERED_FUNCTIONS},
    }


def registration_sha256() -> str:
    import hashlib
    import json

    return hashlib.sha256(json.dumps(registration(), sort_keys=True, ensure_ascii=True).encode("utf-8")).hexdigest()


def excluded_features(verdicts: pd.DataFrame) -> set[str]:
    out = set()
    for key, verdict in zip(verdicts["key"], verdicts["verdict"].astype(str)):
        if verdict.startswith("FAIL"):
            out.update(EXCLUDED_BY_FAIL.get(key, ()))
    return out


def validate_rules(rules, excluded: set[str]) -> None:
    """Raise :class:`StateRuleError` if a rule references an excluded sensor feature."""
    bad = [(r.state, f) for r in rules for f in r.features if f in excluded]
    if bad:
        raise StateRuleError(f"state rules use excluded sensor features: {bad}")


def rules_table(excluded: set[str]) -> pd.DataFrame:
    rows = [{"state": r.state, "condition": r.condition, "features": ", ".join(r.features),
             "order": i} for i, r in enumerate(RULES, start=1)]
    rows.append({"state": ARTIFACT, "condition": "pulse transient (0.1 d), NH4 step region (to resume + 0.1 d), held/damped "
                 "NH4, ORP electrode event, empty NH4/pH/ORP bin, before day 1", "features": "QC flags", "order": 0})
    rows.append({"state": UNCLASSIFIED, "condition": "slope window crosses an artifact, NH4 rising (no N source), "
                 "consumption with a rising pH, or no indicator beyond its threshold", "features": "-", "order": 99})
    df = pd.DataFrame(rows).sort_values("order", kind="stable")
    df["thresholds"] = (f"slope window {(2 * SLOPE_HALF_BINS + 1) * BIN_D:.2f} d (Theil-Sen {100 * CI_ALPHA:g} % interval, "
                        "half-widths x √(n/n_eff), n_eff = n (1 − ρ)/(1 + ρ), ρ = lag-1 autocorrelation of the window "
                        f"residuals, capped at {RHO_MAX}; C(t), pH and ORP alike); "
                        f"pH rising/falling also needs |slope| x window ≥ {PH_QUANTUM} (one quantum); ORP_rel level "
                        f"±{ORP_LEVEL_MV} mV against the centred {2 * ORP_REF_HALF_BINS * BIN_D:g} d median "
                        f"({ORP_THRESHOLD_SCENARIO})")
    df["excluded_features"] = ", ".join(sorted(excluded)) or "none"
    df["registration_sha256"] = registration_sha256()
    df["registration_status"] = REGISTRATION_STATUS
    return df.drop(columns=["order"])


# --- Artifact mask ----------------------------------------------------------------------------------------------
def artifact_reasons(t_mid: np.ndarray, grid: pd.DataFrame, cfg: Config, jumps: pd.DataFrame,
                     held: pd.DataFrame) -> np.ndarray:
    """First applicable artifact reason per bin ('' for a clean bin)."""
    reason = np.full(len(t_mid), "", dtype=object)

    def mark(mask, why):
        reason[(reason == "") & mask] = why

    mark(t_mid < T0_D, "start-up (before day 1)")
    for tp in cfg.pulse_times_d:
        mark((t_mid >= tp) & (t_mid < tp + S.EXCLUDE_AFTER_BOUNDARY_D), "pulse transient")
    if len(jumps):
        bnd = S.boundaries(cfg.pulse_times_d, jumps, start=0.0, end=END_D)
        cand = bnd[~bnd["kind"].astype(str).str.match(r"^(start|end|pulse \d+)$")]
        for a, b in zip(cand["t_d"], cand["t_resume_d"]):
            mark((t_mid >= a - BIN_D) & (t_mid < b + S.EXCLUDE_AFTER_BOUNDARY_D), "NH4 step region")
    for a, b in zip(held["t_start_d"], held["t_end_d"]):
        mark((t_mid >= a) & (t_mid < b), "NH4 held/damped output")
    for te in ORP_EVENTS_D:
        mark((t_mid >= te - ORP_EVENT_MARGIN_D[0]) & (t_mid < te + ORP_EVENT_MARGIN_D[1]), "ORP electrode event")
    for key in ("nh4", "ph", "orp"):
        mark(grid[f"empty_{key}"].to_numpy(dtype=bool), f"empty {key} bin")
    return reason


# --- Features ---------------------------------------------------------------------------------------------------
def lag1_rho(resid: np.ndarray) -> float:
    """Lag-1 autocorrelation of residuals clipped to [0, RHO_MAX] (0 when the residuals have no variance)."""
    e = np.asarray(resid, dtype=float)
    if len(e) < 3 or np.std(e[:-1]) == 0 or np.std(e[1:]) == 0:
        return 0.0
    return min(max(float(np.corrcoef(e[:-1], e[1:])[0, 1]), 0.0), RHO_MAX)


def interval_scale(n: int, rho: float) -> float:
    """√(n/n_eff) with n_eff = n (1 − ρ)/(1 + ρ) (``model.effective_n``); ρ is capped at RHO_MAX, so n_eff > 0."""
    rho = min(max(rho, 0.0), RHO_MAX)
    return math.sqrt(n / (n * (1 - rho) / (1 + rho)))


def rolling_theil_sen(t: np.ndarray, y: np.ndarray, clean: np.ndarray, half: int = SLOPE_HALF_BINS):
    """Theil-Sen slope and autocorrelation-scaled 95 % interval over the centred window of every bin whose window is
    wholly clean (rows: slope, low, high, lag-1 ρ of the window residuals).

    scipy's interval treats the 2 half + 1 points as independent; each half-width (below and above the slope) is
    multiplied by √(n/n_eff), n_eff = n (1 − ρ)/(1 + ρ), ρ = lag-1 autocorrelation of the residuals about the window's
    own Theil-Sen line (added in pass 3, X2, after the first states had been computed)."""
    n = len(t)
    out = np.full((4, n), np.nan)
    bad = np.concatenate([[0], np.cumsum(~clean | ~np.isfinite(y))])
    for i in range(half, n - half):
        if bad[i + half + 1] - bad[i - half] == 0:
            tt, yy = t[i - half:i + half + 1], y[i - half:i + half + 1]
            r = theilslopes(yy, tt, alpha=CI_ALPHA)
            rho = lag1_rho(yy - (r.intercept + r.slope * tt))
            k = interval_scale(len(tt), rho)
            out[:, i] = (r.slope, r.slope - k * (r.slope - r.low_slope), r.slope + k * (r.high_slope - r.slope), rho)
    return out


def ar1(n: int, phi: float, seed: int) -> np.ndarray:
    """AR(1) series x_i = φ x_{i−1} + e_i with unit normal innovations (fixed seed)."""
    e = np.random.default_rng(seed).normal(size=n)
    x = np.empty(n)
    x[0] = e[0]
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def interval_coverage(phis=COVERAGE_PHIS, n_windows: int = COVERAGE_WINDOWS, seed: int = SEED) -> pd.DataFrame:
    """Y3: coverage of the true slope by the production interval at the centre of 2 half + 1 point windows with AR(1)
    errors: unscaled (independent points), production (ρ̂ from the window) and with the true ρ = φ. Synthetic data
    only (NH4-independent); it quantifies the residual undercoverage of the ρ̂ estimator on 101 points."""
    rng = np.random.default_rng(seed)
    n = 2 * SLOPE_HALF_BINS + 1
    t = BIN_D * np.arange(n)
    c = SLOPE_HALF_BINS
    rows = []
    for phi in phis:
        hit = {"unscaled": 0, "estimated": 0, "true": 0}
        rhos = []
        for _ in range(n_windows):
            y = COVERAGE_SLOPE * t + COVERAGE_NOISE * ar1(n + COVERAGE_BURN_IN, phi, int(rng.integers(1 << 31)))[
                COVERAGE_BURN_IN:]
            out = rolling_theil_sen(t, y, np.ones(n, dtype=bool))
            r = theilslopes(y, t, alpha=CI_ALPHA)
            k = interval_scale(n, phi)
            hit["estimated"] += out[1, c] <= COVERAGE_SLOPE <= out[2, c]
            hit["unscaled"] += r.low_slope <= COVERAGE_SLOPE <= r.high_slope
            hit["true"] += (r.slope - k * (r.slope - r.low_slope) <= COVERAGE_SLOPE
                            <= r.slope + k * (r.high_slope - r.slope))
            rhos.append(out[3, c])
        rows.append({"phi": phi, "n_points": n, "n_windows": n_windows, "nominal": CI_ALPHA,
                     "coverage_unscaled": hit["unscaled"] / n_windows,
                     "coverage_scaled_estimated_rho": hit["estimated"] / n_windows,
                     "coverage_scaled_true_rho": hit["true"] / n_windows,
                     "mean_rho_hat": float(np.mean(rhos)), "median_rho_hat": float(np.median(rhos)),
                     "reference": phi == COVERAGE_REFERENCE_PHI})
    out = pd.DataFrame(rows)
    out["basis"] = ("SYNTHETIC: production rolling_theil_sen on AR(1) errors (fixed seed); coverage_scaled_estimated_rho "
                    "is the interval the states use")
    return out


def interval_limits_text(coverage: pd.DataFrame, fractions: pd.DataFrame, interval: str = "whole record") -> str:
    """Y3: the two limits of the idle bounds, with every number from the synthetic check and the state fractions (one
    source for the report and ``states.md``)."""
    cov = coverage.assign(phi=coverage["phi"].astype(float))
    ref = cov.loc[cov["reference"].astype(str).str.lower() == "true"].iloc[0]
    fr = fractions.set_index("interval").loc[interval]
    rho = float(fr["median_nh4_slope_resid_lag1"])
    near = cov.loc[(cov["median_rho_hat"].astype(float) - rho).abs().idxmin()]
    non_idle = float(fr["fraction_active_nitrification"]) + float(fr["fraction_oxygen_limited"])
    return (f"Limits of the idle bounds: (a) the autocorrelation-scaled slope interval under-covers, because ρ̂ from "
            f"{int(ref['n_points'])} detrended points is biased low (synthetic AR(1) check, "
            f"`states_interval_coverage.csv`): it covers the true slope in {float(ref['coverage_scaled_estimated_rho']):.2f} "
            f"of the windows at φ = {float(ref['phi']):g} instead of {float(ref['nominal']):.2f}, and in "
            f"{float(near['coverage_scaled_estimated_rho']):.2f} at φ = {float(near['phi']):g}, whose median ρ̂ "
            f"({float(near['median_rho_hat']):.2f}) is nearest the record's median ρ̂ of C(t) ({rho:.2f}); so some truly "
            "flat windows are still called rising or falling, which moves bins from a state of no change to "
            "unclassified or active. (b) Both bounds are conditional on the state classification: the lower bound "
            "counts only bins classified inactive, and the upper bound counts the active and oxygen-limited bins as "
            f"not idle ({100 * non_idle:.1f} % of the {interval}), so a truly idle bin misclassified as active lowers "
            "the upper bound by at most that share, while the lower bound stays conservative.")


def relative_orp(t: np.ndarray, orp: np.ndarray, clean: np.ndarray) -> np.ndarray:
    """ORP minus its centred 5 d median, computed within the pieces between ORP electrode events."""
    rel = np.full(len(t), np.nan)
    cuts = [-np.inf, *ORP_EVENTS_D, np.inf]
    for a, b in itertools.pairwise(cuts):
        piece = (t >= a) & (t < b)
        s = pd.Series(np.where(clean & piece, orp, np.nan))[piece]
        ref = s.rolling(2 * ORP_REF_HALF_BINS + 1, center=True, min_periods=ORP_REF_MIN_BINS).median()
        rel[piece] = (s - ref).to_numpy()
    return rel


def window_median(x: np.ndarray, valid: np.ndarray, half: int = SLOPE_HALF_BINS) -> np.ndarray:
    """Median of x over the centred slope window where ``valid`` (a finite slope) holds; NaN elsewhere."""
    m = pd.Series(x).rolling(2 * half + 1, center=True, min_periods=2 * half + 1).median().to_numpy()
    return np.where(valid, m, np.nan)


def sign_class(slope, lo, hi, min_abs: float = 0.0) -> np.ndarray:
    """+1 rising, −1 falling, 0 flat, NaN where no slope; significance = interval excludes 0 and |slope| ≥ min_abs."""
    out = np.where(np.isfinite(slope), 0.0, np.nan)
    out[(lo > 0) & (np.abs(slope) >= min_abs)] = 1.0
    out[(hi < 0) & (np.abs(slope) >= min_abs)] = -1.0
    return out


def features(grid: pd.DataFrame, truth: pd.DataFrame, reason: np.ndarray) -> pd.DataFrame:
    g = grid.reset_index(drop=True)
    t = g["t_mid_d"].to_numpy()
    c = pd.Series(truth["C_inv_mgN_L"].to_numpy(), index=np.round(truth["t_d"].to_numpy(), 6))
    c_inv = c.reindex(np.round(t, 6)).to_numpy()
    clean = reason == ""
    nh4 = rolling_theil_sen(t, c_inv, clean)
    ph = rolling_theil_sen(t, g["pH"].to_numpy(), clean)
    orp = rolling_theil_sen(t, g["orp_mV_vs_Ag"].to_numpy(), clean)
    window_d = (2 * SLOPE_HALF_BINS) * BIN_D
    rel = relative_orp(t, g["orp_mV_vs_Ag"].to_numpy(), clean)
    valid = np.isfinite(nh4[0]) & np.isfinite(ph[0]) & np.isfinite(orp[0])
    return pd.DataFrame({
        "t_mid_d": t, "artifact_reason": reason, "C_inv_mgN_L": c_inv, "pH": g["pH"].to_numpy(),
        "orp_mV_vs_Ag": g["orp_mV_vs_Ag"].to_numpy(),
        "nh4_slope_mgN_L_d": nh4[0], "nh4_slope_ci95_low": nh4[1], "nh4_slope_ci95_high": nh4[2],
        "ph_slope_per_d": ph[0], "ph_slope_ci95_low": ph[1], "ph_slope_ci95_high": ph[2],
        "orp_slope_mV_d": orp[0], "orp_slope_ci95_low": orp[1], "orp_slope_ci95_high": orp[2],
        "nh4_slope_resid_lag1": nh4[3], "ph_slope_resid_lag1": ph[3], "orp_slope_resid_lag1": orp[3],
        "orp_rel_mV": rel, "orp_rel_window_median_mV": window_median(rel, valid),
        "nh4_sign": sign_class(nh4[0], nh4[1], nh4[2]),
        "ph_sign": sign_class(ph[0], ph[1], ph[2], PH_QUANTUM / window_d),
        "orp_sign": sign_class(orp[0], orp[1], orp[2]),
    })


# --- Classification --------------------------------------------------------------------------------------------
def classify(f: pd.DataFrame) -> pd.DataFrame:
    """Apply the registered rules, in order, to the feature table (adds ``state`` and ``state_reason``)."""
    n = len(f)
    state = np.full(n, UNCLASSIFIED, dtype=object)
    why = np.full(n, "", dtype=object)
    art = f["artifact_reason"].astype(str).to_numpy() != ""
    state[art], why[art] = ARTIFACT, f["artifact_reason"].to_numpy()[art]
    nh, ph, os_ = f["nh4_sign"].to_numpy(), f["ph_sign"].to_numpy(), f["orp_sign"].to_numpy()
    lvl = f["orp_rel_window_median_mV"].to_numpy()
    have = ~art & np.isfinite(nh) & np.isfinite(ph) & np.isfinite(os_) & np.isfinite(lvl)
    why[~art & ~have] = "slope window crosses an artifact or event"
    low, high = lvl <= -ORP_LEVEL_MV, lvl >= ORP_LEVEL_MV
    rules = [
        (have & (nh > 0), UNCLASSIFIED, "C(t) rising without an N source (drift or evaporation)"),
        (have & (nh < 0) & (ph > 0), UNCLASSIFIED, "consumption with a rising pH (conflicting indicators)"),
        (have & (nh < 0) & (ph <= 0), "active_nitrification", RULES[0].condition),
        (have & (nh == 0) & low & (os_ <= 0) & (ph <= 0), "oxygen_limited", RULES[1].condition),
        (have & (nh == 0) & ~low & (high | (ph > 0)), "inactive_oversupplied", RULES[2].condition),
        (have & (nh == 0), UNCLASSIFIED, "no indicator beyond its threshold, or low ORP with a rising ORP or pH"),
    ]
    done = ~have
    for mask, st, reason in rules:
        m = mask & ~done
        state[m], why[m] = st, reason
        done |= m
    return f.assign(state=state, state_reason=why)


def timeline(bins: pd.DataFrame) -> pd.DataFrame:
    """Time-ordered runs of one state (start, end, duration, medians)."""
    s = bins["state"].to_numpy()
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    ends = np.r_[starts[1:], len(s)]
    rows = []
    for a, b in zip(starts, ends):
        part = bins.iloc[a:b]
        rows.append({"state": s[a], "t_start_d": float(part["t_mid_d"].iloc[0] - BIN_D / 2),
                     "t_end_d": float(part["t_mid_d"].iloc[-1] + BIN_D / 2), "duration_d": (b - a) * BIN_D,
                     "n_bins": b - a, "median_C_inv_mgN_L": float(part["C_inv_mgN_L"].median()),
                     "median_consumption_mgN_L_d": float((-part["nh4_slope_mgN_L_d"]).median()),
                     "median_pH": float(part["pH"].median()),
                     "median_orp_rel_mV": float(part["orp_rel_window_median_mV"].median()),
                     "reason": "; ".join(sorted(set(part["state_reason"].astype(str)) - {""}))[:300]})
    return pd.DataFrame(rows)


def interval_bounds(cfg: Config) -> list[tuple[str, float, float]]:
    p = cfg.pulse_times_d
    return [("before pulse 1", T0_D, p[0]), *[(f"pulse {k}", a, b) for k, (a, b)
                                                   in enumerate(zip(p, [*p[1:], END_D]), start=1)], ("whole record", T0_D, END_D)]


def state_fractions(bins: pd.DataFrame, lo: float, hi: float, min_classified: float = MIN_CLASSIFIED_FRACTION) -> dict:
    """Time fractions of each state in [lo, hi) and the idle fraction as a bounded quantity (pass 3, X2):

    - ``idle_fraction_lower_bound`` = inactive bins / all bins;
    - ``idle_fraction_upper_bound`` = (inactive + unclassified + artifact bins) / all bins (every bin without a positive
      non-idle state could have been idle);
    - ``idle_fraction_of_classified_time`` = inactive / classified bins, CONDITIONAL on the classified bins being
      representative; ``not computed`` (in ``status``) unless ``classified_fraction`` ≥ ``min_classified``.

    The bounds are reported whenever the interval has bins; ``status`` is the status of the conditional value."""
    m = (bins["t_mid_d"] >= lo) & (bins["t_mid_d"] < hi)
    st = bins.loc[m, "state"]
    total = len(st)
    out = {"t_from_d": lo, "t_to_d": hi, "n_bins": total}
    for s in (*STATES, ARTIFACT, UNCLASSIFIED):
        out[f"fraction_{s}"] = float((st == s).sum() / total) if total else np.nan
    classified = int(st.isin(STATES).sum())
    inactive = int((st == "inactive_oversupplied").sum())
    out["n_classified"], out["n_inactive_oversupplied"] = classified, inactive
    out["classified_fraction"] = classified / total if total else np.nan
    out["idle_fraction_lower_bound"] = inactive / total if total else np.nan
    out["idle_fraction_upper_bound"] = (inactive + int(st.isin((UNCLASSIFIED, ARTIFACT)).sum())) / total if total else np.nan
    out["bounds_status"] = "computed" if total else "not computed: no bin in this interval"
    if classified == 0:
        out["idle_fraction_of_classified_time"] = np.nan
        out["status"] = "not computed: no bin with a process state in this interval"
    elif out["classified_fraction"] < min_classified:
        out["idle_fraction_of_classified_time"] = np.nan
        out["status"] = (f"not computed: classified fraction {out['classified_fraction']:.3f} below {min_classified:g} "
                         "(the classified bins are not representative of the interval; see the bounds)")
    else:
        out["idle_fraction_of_classified_time"] = inactive / classified
        out["status"] = "computed (CONDITIONAL on the classified bins being representative)"
    if "nh4_slope_resid_lag1" in bins.columns:  # Y3: compared with the synthetic coverage check
        out["median_nh4_slope_resid_lag1"] = float(bins.loc[m, "nh4_slope_resid_lag1"].median()) if total else np.nan
    return out


# --- Block bootstrap -------------------------------------------------------------------------------------------
def acf_block_length(y: np.ndarray, minimum: int) -> int:
    """Block length: the first lag at which the autocorrelation of y drops below 1/e, at least ``minimum``."""
    y = np.asarray(y, dtype=float) - np.mean(y)
    n = len(y)
    den = float(y @ y)
    for lag in range(1, n // 2):
        if den == 0 or float(y[:-lag] @ y[lag:]) / den < 1 / math.e:
            return max(lag, minimum)
    return max(n // 2, minimum)


def block_bootstrap_diff(y: np.ndarray, group_a: np.ndarray, group_b: np.ndarray, block: int,
                         n_boot: int = N_BOOT, seed: int = SEED, strata: np.ndarray | None = None) -> dict:
    """mean(y | A) − mean(y | B) with a moving-block bootstrap percentile interval.

    The series is time-ordered; blocks of ``block`` consecutive bins keep their group labels, so the autocorrelation
    inside a block is preserved. With ``strata`` (e.g. the pulse interval), blocks are drawn within each stratum and
    the statistic is the bin-weighted mean of the within-stratum differences (strata lacking a group are skipped).
    """
    y = np.asarray(y, dtype=float)
    a, b = np.asarray(group_a, dtype=bool), np.asarray(group_b, dtype=bool)
    if np.any(a & b):
        raise ValueError("the two groups overlap")
    keep = a | b
    y, a, b = y[keep], a[keep], b[keep]
    st = np.zeros(len(y), dtype=int) if strata is None else np.asarray(strata)[keep]
    labels = [s for s in pd.unique(st) if a[st == s].any() and b[st == s].any()]

    def stat(idx_by_stratum):
        num = den = 0.0
        for idx in idx_by_stratum.values():
            ya, yb = y[idx][a[idx]], y[idx][b[idx]]
            if len(ya) == 0 or len(yb) == 0:
                continue
            w = len(idx)
            num += w * (ya.mean() - yb.mean())
            den += w
        return num / den if den else np.nan

    base = {s: np.flatnonzero(st == s) for s in labels}
    est = stat(base)
    out = {"diff": est, "ci95_low": np.nan, "ci95_high": np.nan, "n_a": int(a.sum()), "n_b": int(b.sum()),
           "block_bins": block, "n_strata": len(labels), "n_boot_valid": 0}
    if not labels:
        out["status"] = "not computed: a group is empty in every stratum"
        return out
    if min(out["n_a"], out["n_b"]) < MIN_BLOCKS_PER_GROUP * block:
        out["status"] = (f"not computed: a group has fewer than {MIN_BLOCKS_PER_GROUP} blocks of {block} bins "
                         f"({out['n_a']} and {out['n_b']} bins); the point difference is descriptive only")
        return out
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        res = {}
        for s, idx in base.items():
            n = len(idx)
            L = min(block, n)
            k = math.ceil(n / L)
            starts = rng.integers(0, n - L + 1, size=k)
            res[s] = idx[(starts[:, None] + np.arange(L)).ravel()[:n]]
        v = stat(res)
        if np.isfinite(v):
            draws.append(v)
    out["n_boot_valid"] = len(draws)
    if len(draws) >= 0.9 * n_boot:
        out["ci95_low"], out["ci95_high"] = (float(x) for x in np.percentile(draws, [2.5, 97.5]))
        out["status"] = "computed"
    else:
        out["status"] = f"not computed: only {len(draws)} of {n_boot} resamples contain both groups"
    return out


def indicator_groups(bins: pd.DataFrame) -> dict[str, np.ndarray]:
    lvl = bins["orp_rel_window_median_mV"].to_numpy()
    ph = bins["ph_sign"].to_numpy()
    return {"orp_low": lvl <= -ORP_LEVEL_MV, "orp_high": lvl >= ORP_LEVEL_MV,
            "ph_rising": ph > 0, "ph_not_rising": ph <= 0}


def comparisons(bins: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """NH4 consumption (−slope of C(t)) between NH4-free indicator groups, pooled and within pulse intervals; and,
    descriptively (circular), between the full states."""
    ok = (bins["artifact_reason"] == "") & bins["nh4_slope_mgN_L_d"].notna() & bins["orp_rel_window_median_mV"].notna() \
        & bins["ph_sign"].notna()
    b = bins[ok].reset_index(drop=True)
    y = -b["nh4_slope_mgN_L_d"].to_numpy()
    if len(y) < 2 * (2 * SLOPE_HALF_BINS + 1):
        return pd.DataFrame([{"comparison": "all", "status": f"not computed: {len(y)} eligible bins"}])
    block = acf_block_length(y, 2 * SLOPE_HALF_BINS + 1)
    strata = np.searchsorted(np.asarray(cfg.pulse_times_d), b["t_mid_d"].to_numpy(), side="right")
    g = indicator_groups(b)
    rows = []
    for name, (ga, gb) in INDICATOR_PARTITIONS.items():
        for scope, stv in (("pooled", None), ("within pulse intervals", strata)):
            r = block_bootstrap_diff(y, g[ga], g[gb], block, strata=stv)
            rows.append({"comparison": f"{name} ({scope})", "group_a": ga, "group_b": gb, "role": "ASSOCIATION TEST "
                         "(NH4-free partition; not a cause: time and electrode state are confounded)", **r})
    st = b["state"].to_numpy()
    for sa, sb in (("active_nitrification", "inactive_oversupplied"), ("active_nitrification", "oxygen_limited"),
                   ("oxygen_limited", "inactive_oversupplied")):
        r = block_bootstrap_diff(y, st == sa, st == sb, block)
        rows.append({"comparison": f"{sa} vs {sb} (pooled)", "group_a": sa, "group_b": sb,
                     "role": "DESCRIPTIVE (circular: the states are defined with the NH4 slope)", **r})
    out = pd.DataFrame(rows)
    out["response"] = "NH4 consumption −dC/dt of the M6-inverted C(t), 1.01 d Theil-Sen (mg N/L/d)"
    out["registration_sha256"] = registration_sha256()
    out["registration_status"] = REGISTRATION_STATUS
    return out


# --- Null coverage of the block bootstrap (pass 6, G2) ---------------------------------------------------------
# Fixed before the check was run: the record's own eligible bins and group labels, a response with no group effect
# built like the real one (the 2 half + 1 bin least-squares slope of AR(1) noise, so neighbouring values share all but
# one bin), the registered block rule and N_BOOT. An interval is read as adequate only if it covers the true
# difference 0 in at least NULL_COVERAGE_MIN_OK of the replicates at every φ; the block rule is not changed by it.
NULL_COVERAGE_PHIS = COVERAGE_PHIS
NULL_COVERAGE_REPLICATES = 200
NULL_COVERAGE_MIN_OK = 0.90
NULL_COVERAGE_COLUMNS = ("comparison", "phi", "n_replicates", "n_computed", "coverage_of_zero", "nominal",
                         "block_bins_median", "basis")


def null_response(n: int, phi: float, seed: int, half: int = SLOPE_HALF_BINS) -> np.ndarray:
    """n values of the centred (2 half + 1)-point least-squares slope of AR(1) noise: no group effect, overlap as in the
    real response."""
    w = np.arange(-half, half + 1, dtype=float)
    w /= w @ w
    x = ar1(n + 2 * half + COVERAGE_BURN_IN, phi, seed)[COVERAGE_BURN_IN:]
    return np.convolve(x, w[::-1], mode="valid")


def bootstrap_null_coverage(group_a: np.ndarray, group_b: np.ndarray, strata: np.ndarray | None = None,
                            phis=NULL_COVERAGE_PHIS, n_rep: int = NULL_COVERAGE_REPLICATES, n_boot: int = N_BOOT,
                            seed: int = SEED, min_block: int = 2 * SLOPE_HALF_BINS + 1) -> list[dict]:
    """Coverage of the true difference 0 by ``block_bootstrap_diff`` (registered block rule) on synthetic responses
    with the given group labels; one row per φ."""
    rng = np.random.default_rng(seed)
    n = len(group_a)
    rows = []
    for phi in phis:
        hits = done = 0
        blocks = []
        for _ in range(n_rep):
            y = null_response(n, phi, int(rng.integers(1 << 31)))
            block = acf_block_length(y, min_block)
            r = block_bootstrap_diff(y, group_a, group_b, block, n_boot=n_boot, seed=int(rng.integers(1 << 31)),
                                     strata=strata)
            if r["status"] != "computed":
                continue
            done += 1
            hits += r["ci95_low"] <= 0.0 <= r["ci95_high"]
            blocks.append(block)
        rows.append({"phi": phi, "n_replicates": n_rep, "n_computed": done,
                     "coverage_of_zero": hits / done if done else np.nan, "nominal": CI_ALPHA,
                     "block_bins_median": float(np.median(blocks)) if blocks else np.nan})
    return rows


def null_coverage(bins: pd.DataFrame, comp: pd.DataFrame, cfg: Config, **kw) -> pd.DataFrame:
    """G2: the null coverage of every computed comparison, with the eligible bins and groups of ``comparisons``."""
    ok = (bins["artifact_reason"] == "") & bins["nh4_slope_mgN_L_d"].notna() & bins["orp_rel_window_median_mV"].notna() \
        & bins["ph_sign"].notna()
    b = bins[ok].reset_index(drop=True)
    strata = np.searchsorted(np.asarray(cfg.pulse_times_d), b["t_mid_d"].to_numpy(), side="right")
    g = indicator_groups(b)
    st = b["state"].to_numpy()
    rows = []
    for r in comp[comp["status"] == "computed"].itertuples(index=False):
        ga, gb = (g[r.group_a], g[r.group_b]) if r.group_a in g else (st == r.group_a, st == r.group_b)
        stv = strata if "within pulse intervals" in r.comparison else None
        for row in bootstrap_null_coverage(ga, gb, stv, **kw):
            rows.append({"comparison": r.comparison, **row})
    out = pd.DataFrame(rows, columns=list(NULL_COVERAGE_COLUMNS[:-1]))
    out["basis"] = (f"SYNTHETIC null: the record's eligible bins and group labels, response = {2 * SLOPE_HALF_BINS + 1}-bin "
                    "least-squares slope of AR(1) noise (no group effect), registered block rule and resamples, fixed "
                    "seed; eligible bins treated as contiguous")
    return out


def interval_reading(comp: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """Adds ``null_coverage_min`` and ``interval_reading`` (decided by NULL_COVERAGE_MIN_OK) to the comparison table."""
    worst = cov.groupby("comparison")["coverage_of_zero"].min() if len(cov) else pd.Series(dtype=float)
    out = comp.copy()
    out["null_coverage_min"] = out["comparison"].map(worst)

    def reading(row):
        if row["status"] != "computed":
            return "no interval"
        c = row["null_coverage_min"]
        if not np.isfinite(c):
            return "not checked: no computed null replicate"
        if c >= NULL_COVERAGE_MIN_OK:
            return f"adequate on the synthetic null (coverage of 0 ≥ {NULL_COVERAGE_MIN_OK:g} at every φ)"
        return (f"UNDER-COVERS on the synthetic null (coverage of 0 down to {c:.2f} < {NULL_COVERAGE_MIN_OK:g}): the "
                "interval is too narrow; an interval excluding 0 is not evidence of a difference")

    out["interval_reading"] = out.apply(reading, axis=1)
    return out


# --- Stage -----------------------------------------------------------------------------------------------------
# Y4: the NH4-independent and the C(t)-derived columns of each checked table (the provenance columns aside); a column
# that depends on C(t) only through a state, a run or the bins eligible for a slope is declared C(t)-derived
STATE_TABLE_COLUMNS = {
    "bins": (("t_mid_d", "artifact_reason", "pH", "orp_mV_vs_Ag", "ph_slope_per_d", "ph_slope_ci95_low",
              "ph_slope_ci95_high", "orp_slope_mV_d", "orp_slope_ci95_low", "orp_slope_ci95_high", "ph_slope_resid_lag1",
              "orp_slope_resid_lag1", "orp_rel_mV", "ph_sign", "orp_sign"),
             ("C_inv_mgN_L", "nh4_slope_mgN_L_d", "nh4_slope_ci95_low", "nh4_slope_ci95_high", "nh4_slope_resid_lag1",
              "orp_rel_window_median_mV", "nh4_sign", "state", "state_reason")),
    "timeline": ((), ("state", "t_start_d", "t_end_d", "duration_d", "n_bins", "median_C_inv_mgN_L",
                      "median_consumption_mgN_L_d", "median_pH", "median_orp_rel_mV", "reason")),
    "fractions": (("interval", "t_from_d", "t_to_d", "n_bins", "fraction_sensor_artifact"),
                  ("fraction_active_nitrification", "fraction_inactive_oversupplied", "fraction_oxygen_limited",
                   "fraction_unclassified", "n_classified", "n_inactive_oversupplied", "classified_fraction",
                   "idle_fraction_lower_bound", "idle_fraction_upper_bound", "bounds_status",
                   "idle_fraction_of_classified_time", "status", "median_nh4_slope_resid_lag1")),
    "comparisons": (("comparison", "group_a", "group_b", "role", "response", "registration_sha256",
                     "registration_status", "null_coverage_min", "interval_reading"),
                    ("diff", "ci95_low", "ci95_high", "n_a", "n_b", "block_bins", "n_strata", "n_boot_valid", "status")),
}
COVERAGE_COLUMNS = ("phi", "n_points", "n_windows", "nominal", "coverage_unscaled", "coverage_scaled_estimated_rho",
                    "coverage_scaled_true_rho", "mean_rho_hat", "median_rho_hat", "reference", "basis")


@dataclass
class StateResult:
    bins: pd.DataFrame
    timeline: pd.DataFrame
    fractions: pd.DataFrame
    comparisons: pd.DataFrame
    rules: pd.DataFrame
    excluded: set = field(default_factory=set)
    coverage: pd.DataFrame | None = None  # Y3: synthetic check of the slope interval
    null_coverage: pd.DataFrame | None = None  # G2: synthetic null check of the block bootstrap


def compute(cfg: Config, inputs: dict) -> StateResult:
    excluded = excluded_features(inputs["verdicts"])
    validate_rules(RULES, excluded)
    grid = inputs["grid"]
    g = grid[(grid["t_start_d"] >= -1e-9) & (grid["t_end_d"] <= END_D + 1e-9)].reset_index(drop=True)
    reason = artifact_reasons(g["t_mid_d"].to_numpy(), g, cfg, inputs["jumps"], inputs["held"])
    f = features(g, inputs["truth"], reason)
    bins = classify(f)
    bins = bins[bins["t_mid_d"] >= T0_D].reset_index(drop=True)
    fr = pd.DataFrame([{"interval": name, **state_fractions(bins, lo, hi)} for name, lo, hi in interval_bounds(cfg)])
    # X3: every table carrying C(t) quantities (slopes, consumption, states built on them) declares its NH4 basis
    prov = {"nh4_observation_model": str(inputs["truth"]["nh4_observation_model"].iloc[0]),
            "nh4_series": str(inputs["truth"]["nh4_truth_id"].iloc[0])}
    comp = comparisons(bins, cfg)
    null = null_coverage(bins, comp, cfg)
    comp = interval_reading(comp, null)
    out = [x.assign(**prov) for x in (bins, timeline(bins), fr, comp)]
    for x, name in zip(out, ("bins", "timeline", "fractions", "comparisons")):
        indep, derived = STATE_TABLE_COLUMNS[name]
        model.check_single_observation_model(x, label_cols=(), independent_cols=indep, derived_cols=derived)
    cov = interval_coverage().assign(nh4_observation_model=model.NH4_INDEPENDENT, nh4_series="")
    model.check_single_observation_model(cov, label_cols=(), independent_cols=COVERAGE_COLUMNS)
    null = null.assign(nh4_observation_model=model.NH4_INDEPENDENT, nh4_series="")
    model.check_single_observation_model(null, label_cols=(), independent_cols=NULL_COVERAGE_COLUMNS)
    return StateResult(*out, rules_table(excluded), excluded, cov, null)


def load_inputs(out_dir: Path) -> dict:
    files = {"grid": "processed/grid_0p01d.parquet", "truth": "processed/nh4_truth.parquet",
             "jumps": "tables/nh4_jumps.csv", "held": "tables/nh4_held_segments.csv",
             "verdicts": "tables/sensor_verdicts.csv"}
    out = {}
    for k, rel in files.items():
        p = Path(out_dir) / rel
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run the earlier stages (`python -m bnr_neq run-all`) first")
        out[k] = pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)
    return out


def write_report(cfg: Config, r: StateResult, path: Path) -> None:
    from .report import md_table

    whole = r.fractions[r.fractions["interval"] == "whole record"].iloc[0]
    lines = [
        "# Process states from relative ORP, pH and the inverted NH4 (W3)",
        "",
        "Generated by `python -m bnr_neq states`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Registered rules (revised after first results; not pre-registered)",
        "",
        (f"Fixed from sensor resolution and physical scales before any state was computed (see `states.py`), then "
         f"revised after the first results: {REGISTRATION_STATUS}. SHA-256 of the rules, parameters and decision "
         f"functions: `{registration_sha256()}`. Rules are "
         "applied in order; a bin flagged by QC carries no process state. Excluded sensor features (from the M2 "
         f"verdicts): {', '.join(sorted(r.excluded)) or 'none'}; a rule that used one would raise."),
        "",
        md_table(r.rules),
        "",
        "## Time fractions per interval",
        "",
        (f"Over the whole record (day {T0_D:g} to {END_D:g}) {100 * whole['classified_fraction']:.0f} % of the time has a "
         f"process state (coverage). The idle fraction (inactive, over-supplied) is bounded: at least "
         f"{100 * whole['idle_fraction_lower_bound']:.0f} % of the time (inactive bins over all bins) and at most "
         f"{100 * whole['idle_fraction_upper_bound']:.0f} % (inactive, unclassified and artifact bins over all bins). "
         f"Over the classified time only, CONDITIONAL on those bins being representative, it is "
         + (f"{100 * whole['idle_fraction_of_classified_time']:.0f} %" if pd.notna(whole["idle_fraction_of_classified_time"])
            else f"{whole['status']}")
         + f" (computed only where the classified fraction is ≥ {MIN_CLASSIFIED_FRACTION:g}). Slope intervals are "
         "widened for the autocorrelation inside each window (`*_slope_resid_lag1` in `states_bins.parquet`)."),
        "",
        interval_limits_text(r.coverage, r.fractions),
        "",
        md_table(r.fractions, digits=3),
        "",
        "Synthetic check of the slope interval (`states_interval_coverage.csv`):",
        "",
        md_table(r.coverage.drop(columns=["nh4_observation_model", "nh4_series"]), digits=3),
        "",
        "## NH4 consumption between indicator groups (block bootstrap)",
        "",
        ("The full states use the NH4 slope, so their NH4 comparison is circular and only descriptive. The association "
         "test uses ORP and pH alone. Moving-block bootstrap (fixed seed, "
         f"{N_BOOT} resamples): the block length is the larger of the slope window and the first lag where the "
         "autocorrelation of the response drops below 1/e, so overlapping windows and autocorrelation are kept inside "
         "blocks. 'Within pulse intervals' compares bins of the same pulse interval only, which removes the "
         "between-interval time trend but not the trend inside an interval. These are associations, not causes."),
        "",
        md_table(r.comparisons, digits=3),
        "",
        (f"Synthetic null check of these intervals (pass 6, G2; `states_bootstrap_null_coverage.csv`): the record's "
         f"eligible bins and group labels with a response that has no group effect ({NULL_COVERAGE_REPLICATES} "
         f"replicates per φ; {2 * SLOPE_HALF_BINS + 1}-bin least-squares slope of AR(1) noise), the registered block "
         f"rule and {N_BOOT} resamples. An interval is read as adequate only if it covers 0 in at least "
         f"{NULL_COVERAGE_MIN_OK:g} of the replicates at every φ (`interval_reading`); the block rule is not changed by "
         "this check."),
        "",
        md_table(r.null_coverage.drop(columns=["nh4_observation_model", "nh4_series", "basis"]), digits=3)
        if r.null_coverage is not None and len(r.null_coverage) else "No computed comparison to check.",
        "",
        "## State timeline (time-ordered runs)",
        "",
        md_table(r.timeline[r.timeline["duration_d"] >= 0.1], digits=3),
        "",
        "Runs shorter than 0.1 d are omitted here; `states_timeline.csv` lists every run.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir) -> StateResult:
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    inputs = load_inputs(out_dir)
    d = out_dirs(out_dir)
    r = compute(cfg, inputs)
    r.bins.to_parquet(d["processed"] / "states_bins.parquet", index=False)
    r.timeline.to_csv(d["tables"] / "states_timeline.csv", index=False)
    r.fractions.to_csv(d["tables"] / "states_fractions.csv", index=False)
    r.comparisons.to_csv(d["tables"] / "states_comparison.csv", index=False)
    r.rules.to_csv(d["tables"] / "states_rules.csv", index=False)
    r.coverage.to_csv(d["tables"] / "states_interval_coverage.csv", index=False)
    r.null_coverage.to_csv(d["tables"] / "states_bootstrap_null_coverage.csv", index=False)
    write_report(cfg, r, d["reports"] / "states.md")
    return r
