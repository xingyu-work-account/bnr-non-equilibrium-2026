"""Phase 4: reduced NH4 model with the electrode observation model (M6).

State: true NH4 in mg N/L, C(t), with impulses of +dose at the pulse times and a
non-negative consumption rate r. Observation: log R = a + b(t) log C + offsets at
the detected electrode jumps, b(t) = b0 + b1 (t − t0), followed by a first-order
sensor lag τ. Fitted with ``scipy.optimize.least_squares`` on the 0.01 d grid,
residuals in log space, with a chronological calibration / selection / test split.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.signal import lfilter

from . import thermo
from .config import Config
from .constants import (
    ANTHONISEN_EXP_K,
    EQ_ALK_PER_MOL_N_NITRIFIED,
    K_ALK_ASM3_MEQ_L,
    K_NH_ASM_MG_N_L,
    KELVIN_OFFSET,
    M_N_G_MOL,
    NH3_PER_N,
    SECONDS_PER_DAY,
)

T0_D = 1.0
END_D = 47.92
CHI2_95_1DOF = 3.841
PROFILE_POINTS_PER_SIDE = 3
PROFILE_MAX_NFEV = 40  # per profile grid point; warm-started from the neighbouring point
PROFILE_XTOL_FRACTION = 1e-3  # root-finding tolerance as a fraction of the profile half-width
PROFILE_MAXITER = 40
MULTIMODAL_LOG_TOL = math.log(1.05)  # another minimum: Alk0 more than 5 % away from the fitted value
# X1: one convergence policy for every candidate's calibration (fixed from the pass-3 diagnosis before any refit)
CAL_START_FACTORS = (0.5, 1.0, 2.0, 4.0)  # x the default kinetic start (Alk0: 5, 10, 20, 40 meq/L)
CAL_MAX_NFEV = 8000  # per start, every candidate; > 2 x the 3,620 evaluations a converging W5 diagnostic run needed
CONVERGENCE_POLICY = (f"{len(CAL_START_FACTORS)} starts per candidate, the default kinetic start multiplied by "
                      f"{', '.join(f'{f:g}' for f in CAL_START_FACTORS)} (Alk0 5, 10, 20, 40 meq/L), observation "
                      f"parameters at their default start; at most {CAL_MAX_NFEV} function evaluations per start; scipy "
                      "least_squares (trf) default tolerances, converged = status > 0; the calibration is the converged "
                      "start with the lowest RSS; the whole-record refit of the selected model has the same budget.")
N_WORKERS = int(os.environ.get("BNR_NEQ_WORKERS", os.cpu_count() or 1))  # results do not depend on this

# --- W1: the single NH4 observation model ---------------------------------------------------
# Every downstream NH4 quantity (rates, free ammonia, conversion fraction, oxygen demand, affinity step) is derived
# from C(t) obtained by inverting this observation model; the linear gain correction survives only as a labelled
# SENSITIVITY case and the raw electrode reading only as a measured slope.
OBSERVATION_MODEL = "M6_log"
RAW_READING = "electrode_reading"  # measured slope of the reading; never an NH4 truth, never feeds a derived quantity
SENSITIVITY_SUFFIX = "_SENSITIVITY"
LINEAR_GAIN = "linear_gain" + SENSITIVITY_SUFFIX
IDENTITY = "identity_reading_as_TAN" + SENSITIVITY_SUFFIX
NH4_INDEPENDENT = "NH4-independent"  # explicit marker of a row that rests on no NH4 quantity (time, pH, capacity, ...)
# X3: allow-list of the measured quantities of the raw reading; anything else on an electrode_reading row is derived
# from removal, rate or concentration of the reading and raises (a deny-list of derived prefixes missed new metrics)
RAW_MEASURED_METRICS = ("consumption_uncorrected_first_window_mgN_L_d", "consumption_uncorrected_weighted_mgN_L_d")
RAW_MEASURED_COLUMNS = ("slope_mgNH4_L_d", "slope_ci95_low_mgNH4_L_d", "slope_ci95_high_mgNH4_L_d", "slope_mgN_L_d",
                        "slope_ci95_low_mgN_L_d", "slope_ci95_high_mgN_L_d", "slope_excl_held_mgNH4_L_d",
                        "slope_excl_held_ci95_low_mgNH4_L_d", "slope_excl_held_ci95_high_mgNH4_L_d",
                        "slope_fragment_diagnostic_mgNH4_L_d", "median_nh4_mgN_L", "response_fraction_governing_pulse")
# Y4: allow-list of the metrics that may declare NH4-independent (time, pH, transfer capacity, standard free energy,
# measured stirrer energy); ``*`` stands for a numeric scenario value only. Any other NH4-independent metric raises.
NH4_INDEPENDENT_METRICS = ("t_pulse_d", "segment_end_d", "median_pH", "o2_supply_capacity_mgO2_SCENARIO_kLa*_T25C",
                           "free_energy_per_fully_converted_pulse_J_SCENARIO_T25C", "stirrer_energy_J")
_NH4_INDEPENDENT_RE = "^(" + "|".join(re.escape(m).replace(r"\*", r"[0-9][0-9.e+-]*") for m in NH4_INDEPENDENT_METRICS) + ")$"
MODEL_BASED_TOKEN = r"(?<![Nn]ot )(?<![A-Za-z0-9])MODEL[-_]BASED(?![A-Za-z0-9])"  # "not MODEL-BASED" is not the label
INVERSION_NOTE =("quasi-static inversion C = exp((log R − a − offsets)/b(t)) of the 0.01 d grid reading; the sensor "
                  "lag τ is not inverted (slopes after ~3τ are unaffected, window starts within 3τ of a boundary carry "
                  "the lag transient)")


class ObservationModelError(ValueError):
    """An output combines NH4 quantities that rest on different observation models."""


def _as_text(s: pd.Series) -> pd.Series:
    """Cells as text; empty cells (NaN when read back from CSV) count as empty text."""
    return s.map(lambda v: "" if pd.isna(v) else str(v))


def check_single_observation_model(df: pd.DataFrame, model_col: str = "nh4_observation_model",
                                   label_cols: tuple[str, ...] = ("metric", "basis"),
                                   series_col: str = "nh4_series",
                                   independent_cols: tuple[str, ...] | None = None,
                                   derived_cols: tuple[str, ...] | None = None) -> str | None:
    """Raise :class:`ObservationModelError` unless every row of ``df`` declares its NH4 basis and every primary row
    uses one observation model and one series. Returns the primary model.

    - Every row declares a basis (X3): ``OBSERVATION_MODEL``, a model ending in ``_SENSITIVITY``, ``RAW_READING`` or
      the explicit marker ``NH4_INDEPENDENT``. An empty or unknown model raises.
    - Rows of a ``_SENSITIVITY`` model say SENSITIVITY in their label columns and, in a long table, in the metric name;
      no other row may say SENSITIVITY.
    - A ``RAW_READING`` row may carry only a measured quantity of the reading: in a long table (``metric`` column) a
      metric of ``RAW_MEASURED_METRICS``. In a wide table (one column per quantity), every column of a reading-based
      row is a measured column (``RAW_MEASURED_COLUMNS``), says SENSITIVITY in its name, or is declared NH4-independent
      in ``independent_cols``.
    - ``NH4_INDEPENDENT`` rows carry no NH4 series. In a long table their metric is on ``NH4_INDEPENDENT_METRICS``
      (Y4); a wide table with such rows declares ``independent_cols`` and carries no other column.
    - A wide table with primary rows declares its NH4-independent columns (``independent_cols``) and its C(t)-derived
      columns (``derived_cols``); an undeclared column raises (Y4).
    - Primary rows share one series (:func:`check_single_series`).
    """
    m = _as_text(df[model_col])
    cols = [c for c in label_cols if c in df.columns]
    text = (df[cols].apply(_as_text).agg(" ".join, axis=1) if cols else pd.Series("", index=df.index))
    metric = _as_text(df["metric"]) if "metric" in df.columns else None
    sens = m.str.endswith(SENSITIVITY_SUFFIX)
    raw, indep = m == RAW_READING, m == NH4_INDEPENDENT

    def where(mask):
        return sorted(set(metric[mask])) if metric is not None else f"{int(mask.sum())} row(s)"

    if (m == "").any():
        raise ObservationModelError(f"rows without a declared NH4 basis (empty {model_col}): {where(m == '')}")
    if (sens & ~text.str.contains("SENSITIVITY")).any():
        raise ObservationModelError(f"sensitivity rows without a SENSITIVITY label: {sorted(set(m[sens]))}")
    if metric is not None and (sens & ~metric.str.contains("SENSITIVITY")).any():
        raise ObservationModelError(f"reading-based metrics without SENSITIVITY in their names: {where(sens & ~metric.str.contains('SENSITIVITY'))}")
    if (~sens & text.str.contains("SENSITIVITY")).any():
        raise ObservationModelError("a row labelled SENSITIVITY carries a primary observation model")
    if metric is not None:
        raw_derived = raw & ~metric.isin(RAW_MEASURED_METRICS)
        if raw_derived.any():
            raise ObservationModelError(f"derived quantities from the raw electrode reading: {where(raw_derived)}")
    elif (raw | sens).any():
        if independent_cols is None:
            raise ObservationModelError("a wide table with reading-based rows must declare its NH4-independent columns")
        skip = {model_col, series_col, *cols, *independent_cols}
        bad = [c for c in df.columns if c not in skip and c not in RAW_MEASURED_COLUMNS and "SENSITIVITY" not in c]
        if bad:
            raise ObservationModelError(f"reading-based columns without SENSITIVITY in their names: {bad}")
    if series_col in df.columns and (indep & (_as_text(df[series_col]) != "")).any():
        raise ObservationModelError(f"{NH4_INDEPENDENT} rows carry an NH4 series: {where(indep & (_as_text(df[series_col]) != ''))}")
    is_primary = ~sens & ~raw & ~indep
    if metric is not None:
        outside = indep & ~metric.str.match(_NH4_INDEPENDENT_RE)
        if outside.any():
            raise ObservationModelError(f"{NH4_INDEPENDENT} metrics outside the allow-list NH4_INDEPENDENT_METRICS: "
                                        f"{where(outside)}")
    else:
        prov = {model_col, series_col, *cols}
        if indep.any():
            if independent_cols is None:
                raise ObservationModelError(f"a wide table with {NH4_INDEPENDENT} rows must declare its NH4-independent "
                                            "columns")
            bad = [c for c in df.columns if c not in prov and c not in independent_cols]
            if bad:
                raise ObservationModelError(f"columns of {NH4_INDEPENDENT} rows not declared NH4-independent: {bad}")
        if is_primary.any():
            if independent_cols is None or derived_cols is None:
                raise ObservationModelError("a wide table with primary rows must declare its NH4-independent and its "
                                            "C(t)-derived columns")
            bad = [c for c in df.columns if c not in prov and c not in independent_cols and c not in derived_cols]
            if bad:
                raise ObservationModelError(f"undeclared columns in a wide table with primary rows: {bad}")
    primary = sorted(set(m[is_primary]))
    if len(primary) > 1:
        raise ObservationModelError(f"one output mixes NH4 observation models: {primary}")
    if primary and primary[0] != OBSERVATION_MODEL:
        raise ObservationModelError(f"undeclared NH4 basis {primary[0]!r}: expected {OBSERVATION_MODEL}, a model ending "
                                    f"in {SENSITIVITY_SUFFIX}, {RAW_READING} or {NH4_INDEPENDENT}")
    if series_col in df.columns and is_primary.any():
        check_single_series(df[is_primary], text[is_primary], series_col)
    return primary[0] if primary else None


MODEL_STATE_PREFIX = "M6_state:"  # series of the simulated model state (MODEL-BASED quantities)


def check_single_series(rows: pd.DataFrame, text: pd.Series, series_col: str = "nh4_series") -> None:
    """Primary rows must share one inverted C(t) series; MODEL-BASED rows (and only they) may instead use one simulated
    model-state series of the same observation model. Raises :class:`ObservationModelError` otherwise."""
    s = _as_text(rows[series_col])
    model_based = text.str.contains(MODEL_BASED_TOKEN, regex=True)
    if (s == "").any():
        raise ObservationModelError("primary NH4 rows without a series identifier")
    state = s.str.startswith(MODEL_STATE_PREFIX)
    if (state & ~model_based).any():
        raise ObservationModelError(f"model-state series on rows not labelled MODEL-BASED: {sorted(set(s[state & ~model_based]))}")
    if (model_based & ~state).any():
        raise ObservationModelError("rows labelled MODEL-BASED carry an inverted-series identifier")
    for what, sel in (("inverted C(t)", ~state), ("model-state", state)):
        ids = sorted(set(s[sel]))
        if len(ids) > 1:
            raise ObservationModelError(f"one output mixes {what} series: {ids}")


def truth_id(c: np.ndarray) -> str:
    """Identifier of one inverted C(t) series (provenance of every quantity derived from it)."""
    return "C_inv:" + hashlib.sha256(np.ascontiguousarray(c, dtype=np.float64).tobytes()).hexdigest()[:12]


@dataclass(frozen=True)
class SplitWindows:
    cal: tuple[float, float]
    sel: tuple[float, float]
    test: tuple[float, float]


def split_windows(cfg: Config) -> SplitWindows:
    """Chronological split from the configured pulse times (M6 item 3): calibrate from day T0 to the
    second-to-last pulse, select between the last two pulses, test from the last pulse to the end."""
    p = cfg.pulse_times_d
    if len(p) < 2:
        raise ValueError("the chronological split needs at least two pulses in config [pulses] times_d")
    return SplitWindows((T0_D, p[-2]), (p[-2], p[-1]), (p[-1], END_D))
LOG_KIN_BOUNDS = (math.log(1e-4), math.log(1e3))


@dataclass(frozen=True)
class ModelSpec:
    name: str
    kinetic: tuple[str, ...]  # fitted on log scale
    fa: bool  # rate depends on free ammonia (needs a temperature SCENARIO)
    decay: bool  # first-order decay of biomass activity
    T_C: float | None = None
    alk: bool = False  # W5: Monod limitation by the alkalinity state Alk0 − 2 n_consumed / M_N (substrate side)


# Y6: the Haldane law r = vmax S/(Ks + S + S²/Ki) in the inhibition branch is a three-parameter valley in (vmax, Ks, Ki)
# (|corr| 0.99998 at the alkalinity x FA stopping point, docs/progress.md); it is fitted as the initial slope
# α = vmax/Ks, Ks and the substrate of maximum rate S_opt = √(Ks Ki): r = α S/(1 + S/Ks + S²/S_opt²)
HALDANE_SHAPE = ("Ks_mgNH3_L", "Sopt_mgNH3_L")
HALDANE_ALPHA = ("alpha_mgN_L_d_per_mgNH3_L", "alpha0_mgN_L_d_per_mgNH3_L")
# Z3: diagnostic of the reparametrized optimum, reproduced in pass 5 by a scratch script (docs/progress.md); not computed
# by the pipeline
HALDANE_CONDITIONING_NOTE = ("finite-difference Jacobian at the converged alkalinity_FA_SCENARIO_T25C calibration, start "
                             "x0.5: condition number 4.2e6, |corr(Ks, S_opt)| 0.89")


def haldane_natural(kin: dict[str, float]) -> dict[str, float]:
    """(vmax or v0, Ks, Ki) of a Haldane law fitted as (α, Ks, S_opt); other parameters unchanged."""
    out = dict(kin)
    for a, v in zip(HALDANE_ALPHA, ("vmax_mgN_L_d", "v0_mgN_L_d")):
        if a in kin:
            out[v] = kin[a] * kin["Ks_mgNH3_L"]
    if "Sopt_mgNH3_L" in kin:
        out["Ki_mgNH3_L"] = kin["Sopt_mgNH3_L"] ** 2 / kin["Ks_mgNH3_L"]
    return out


def model_specs(cfg: Config) -> list[ModelSpec]:
    specs = [ModelSpec("null_no_consumption", (), False, False),
             ModelSpec("biomass_decay", ("v0_mgN_L_d", "kd_per_d"), False, True)]
    # Y6: every Haldane law is fitted as (α = vmax/Ks, Ks, S_opt = √(Ks Ki)), the same function (HALDANE_ALPHA, HALDANE_SHAPE)
    for T in cfg.get("scenarios.temperature_C"):
        specs.append(ModelSpec(f"haldane_FA_SCENARIO_T{T:g}C", ("alpha_mgN_L_d_per_mgNH3_L", *HALDANE_SHAPE), True, False, T))
        specs.append(ModelSpec(f"combined_SCENARIO_T{T:g}C", ("alpha0_mgN_L_d_per_mgNH3_L", "kd_per_d", *HALDANE_SHAPE),
                               True, True, T))
    # W5: substrate-side limitation by alkalinity / inorganic carbon, alone and with FA inhibition (TAN = model state,
    # i.e. the same W1 observation model), competing with biomass decay under the unchanged chronological split
    specs.append(ModelSpec("alkalinity_limited", ("v0_mgN_L_d", "Alk0_meq_L"), False, False, alk=True))
    for T in cfg.get("scenarios.temperature_C"):
        specs.append(ModelSpec(f"alkalinity_FA_SCENARIO_T{T:g}C", ("alpha_mgN_L_d_per_mgNH3_L", *HALDANE_SHAPE,
                                                                  "Alk0_meq_L"), True, False, T, alk=True))
    return specs


def alkalinity_state_meq_L(alk0_meq_L: float, consumed_mgN_L: float) -> float:
    """Alkalinity left after nitrifying ``consumed`` mg N/L (2 eq/mol N, nothing added), never below zero."""
    return max(0.0, alk0_meq_L - EQ_ALK_PER_MOL_N_NITRIFIED * consumed_mgN_L / M_N_G_MOL)


# --- Data ------------------------------------------------------------------------------
@dataclass
class ModelData:
    t: np.ndarray  # bin mid times, d
    log_r: np.ndarray  # log of electrode NH4, mg N/L
    ph: np.ndarray
    pulses: list[float]
    dose: float
    offsets: list[float]  # electrode offset jump times, d
    offset_guess: list[float]  # log level change at each offset (M2)
    volume: object = None  # Y7: a SCENARIO ``volume.VolumeScenario``; None = the constant-volume (primary) model


def build_data(cfg: Config, grid: pd.DataFrame, jumps: pd.DataFrame) -> ModelData:
    g = grid[(grid["t_start_d"] >= T0_D - 1e-9) & (grid["t_end_d"] <= END_D + 1e-9)]
    if g["empty_nh4"].any() or g["empty_ph"].any():
        raise ValueError("model grid has empty NH4 or pH bins")
    off = jumps[jumps["jump"] & (jumps["class"] == "electrode_offset")].sort_values("t_jump_d")
    return ModelData(
        t=g["t_mid_d"].to_numpy(), log_r=np.log(g["nh4_mgN_L"].to_numpy()), ph=g["pH"].to_numpy(),
        pulses=cfg.pulse_times_d, dose=thermo.pulse_dose_from_config(cfg).n_mg_L,
        offsets=off["t_jump_d"].tolist(),
        offset_guess=np.log(off["post_level"].to_numpy() / off["pre_level"].to_numpy()).tolist())


# --- State model -----------------------------------------------------------------------------
def fa_factor(ph: np.ndarray, T_C: float) -> np.ndarray:
    """Free ammonia per unit TAN, mg NH3/L per mg N/L (Anthonisen)."""
    h = 10.0 ** ph
    return NH3_PER_N * h / (np.exp(ANTHONISEN_EXP_K / (T_C + KELVIN_OFFSET)) + h)


def make_rate(spec: ModelSpec, kin: dict[str, float]):
    """r(t, C, f, n, φ) >= 0 in mg N/L/d; f = free-ammonia factor at t; n = N consumed since T0 (mg N/L), which sets the
    alkalinity state of the W5 candidates (ignored by the others); φ = V0/V(t) under a SCENARIO volume (Y7; 1 at
    constant volume, where the arithmetic is unchanged)."""
    if spec.name == "null_no_consumption":
        return lambda t, c, f, n=0.0, phi=1.0: 0.0
    alpha = kin.get(HALDANE_ALPHA[0], kin.get(HALDANE_ALPHA[1]))
    if alpha is not None:  # Y6: the (α, Ks, S_opt) form
        return _haldane_alpha_rate(spec, alpha, kin)
    v = kin.get("v0_mgN_L_d", kin.get("vmax_mgN_L_d"))
    kd = kin.get("kd_per_d", 0.0)
    ks, ki = kin.get("Ks_mgNH3_L"), kin.get("Ki_mgNH3_L")
    alk0 = kin.get("Alk0_meq_L")

    def rate(t, c, f, n=0.0, phi=1.0):
        c = max(0.0, c)
        act = v * math.exp(-kd * (t - T0_D)) if spec.decay else v
        if phi != 1.0:  # Y7: a SCENARIO volume concentrates the biomass by φ = V0/V(t)
            act *= phi
        if spec.alk:
            a = alkalinity_state_meq_L(alk0, n)
            if phi != 1.0:
                a *= phi
            act *= a / (K_ALK_ASM3_MEQ_L + a)
        if spec.fa:
            s = f * c
            return act * s / (ks + s + s * s / ki)
        return act * c / (K_NH_ASM_MG_N_L + c)

    return rate


def _haldane_alpha_rate(spec: ModelSpec, alpha: float, kin: dict[str, float]):
    """r = act α S/(1 + S/Ks + S²/S_opt²), S = FA; act = exp(−kd (t − 1)) for decay, the alkalinity factor for the
    W5 variant, φ under a SCENARIO volume (Y6, Y7)."""
    kd = kin.get("kd_per_d", 0.0)
    ks, so2 = kin["Ks_mgNH3_L"], kin["Sopt_mgNH3_L"] ** 2
    alk0 = kin.get("Alk0_meq_L")

    def rate(t, c, f, n=0.0, phi=1.0):
        c = max(0.0, c)
        act = alpha * math.exp(-kd * (t - T0_D)) if spec.decay else alpha
        if phi != 1.0:
            act *= phi
        if spec.alk:
            a = alkalinity_state_meq_L(alk0, n)
            if phi != 1.0:
                a *= phi
            act *= a / (K_ALK_ASM3_MEQ_L + a)
        s = f * c
        return act * s / (1.0 + s / ks + s * s / so2)

    return rate


@dataclass
class Simulation:
    c: np.ndarray
    rate: np.ndarray
    consumed: float
    pulse_jumps: list[tuple[float, float, float]] = field(default_factory=list)  # (t, before, after)


def simulate(spec: ModelSpec, kin: dict[str, float], c0: float, t: np.ndarray, ph: np.ndarray,
             pulses: list[float], dose: float) -> Simulation:
    """RK4 on the bin times; a pulse inside a step splits the step and adds exactly ``dose``."""
    rate = make_rate(spec, kin)
    f = fa_factor(ph, spec.T_C).tolist() if spec.fa else [0.0] * len(t)
    tl = t.tolist()
    n = len(tl)
    c = c0
    cs = [0.0] * n
    rs = [0.0] * n
    cs[0] = c
    rs[0] = rate(tl[0], c, f[0])
    consumed = 0.0
    jumps = []
    pk = [p for p in pulses if tl[0] < p <= tl[-1]]
    ip = 0
    for i in range(1, n):
        ta, tb, fa, fb = tl[i - 1], tl[i], f[i - 1], f[i]
        span = tb - ta
        stops = [ta]
        while ip < len(pk) and pk[ip] <= tb:
            stops.append(pk[ip])
            ip += 1
        stops.append(tb)
        for j in range(len(stops) - 1):
            s0, s1 = stops[j], stops[j + 1]
            h = s1 - s0
            if h > 0:
                f0 = fa + (fb - fa) * (s0 - ta) / span
                f1 = fa + (fb - fa) * (s1 - ta) / span
                fm = 0.5 * (f0 + f1)
                k1 = rate(s0, c, f0, consumed)
                k2 = rate(s0 + h / 2, c - h / 2 * k1, fm, consumed + h / 2 * k1)
                k3 = rate(s0 + h / 2, c - h / 2 * k2, fm, consumed + h / 2 * k2)
                k4 = rate(s1, c - h * k3, f1, consumed + h * k3)
                d = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                c -= d
                consumed += d
            if j < len(stops) - 2:  # a pulse at s1
                jumps.append((s1, c, c + dose))
                c += dose
        cs[i] = c
        rs[i] = rate(tb, c, fb, consumed)
    return Simulation(np.array(cs), np.array(rs), consumed, jumps)


def simulate_volume(spec: ModelSpec, kin: dict[str, float], c0: float, t: np.ndarray, ph: np.ndarray,
                    pulses: list[float], dose: float, vol) -> Simulation:
    """Y7: :func:`simulate` under a SCENARIO volume V(t) = V0 − E (t − t0) (``volume.VolumeScenario``). Every conserved
    amount is concentrated by φ(t) = V0/V(t): dC/dt = −r + C E/V(t), with the biomass activity and the alkalinity state
    of the rate law multiplied by φ, and n (N consumed per initial litre) growing at r/φ; a pulse adds its N mass into
    V(t_p), i.e. ``dose`` (at V0) x φ(t_p). RK4 on the bin times with the pulse splitting of :func:`simulate`."""
    rate = make_rate(spec, kin)
    f = fa_factor(ph, spec.T_C).tolist() if spec.fa else [0.0] * len(t)
    v0, e, t0 = float(vol.v0_L), float(vol.net_loss_L_per_d), float(vol.t0_d)

    def deriv(s, c, fs, n):
        v = v0 - e * (s - t0)
        if v <= 0:
            raise ValueError(f"volume scenario {vol.name!r} empties the vessel before t = {s:g} d")
        phi = v0 / v
        r = rate(s, c, fs, n, phi)
        return r - max(c, 0.0) * e / v, r / phi, r  # (−dC/dt, dn/dt, volumetric rate)

    tl = t.tolist()
    n_bins = len(tl)
    c, consumed = c0, 0.0
    cs, rs = [0.0] * n_bins, [0.0] * n_bins
    cs[0], rs[0] = c, deriv(tl[0], c, f[0], 0.0)[2]
    jumps = []
    pk = [p for p in pulses if tl[0] < p <= tl[-1]]
    ip = 0
    for i in range(1, n_bins):
        ta, tb, fa, fb = tl[i - 1], tl[i], f[i - 1], f[i]
        span = tb - ta
        stops = [ta]
        while ip < len(pk) and pk[ip] <= tb:
            stops.append(pk[ip])
            ip += 1
        stops.append(tb)
        for j in range(len(stops) - 1):
            s0, s1 = stops[j], stops[j + 1]
            h = s1 - s0
            if h > 0:
                f0 = fa + (fb - fa) * (s0 - ta) / span
                f1 = fa + (fb - fa) * (s1 - ta) / span
                fm = 0.5 * (f0 + f1)
                k1c, k1n, _ = deriv(s0, c, f0, consumed)
                k2c, k2n, _ = deriv(s0 + h / 2, c - h / 2 * k1c, fm, consumed + h / 2 * k1n)
                k3c, k3n, _ = deriv(s0 + h / 2, c - h / 2 * k2c, fm, consumed + h / 2 * k2n)
                k4c, k4n, _ = deriv(s1, c - h * k3c, f1, consumed + h * k3n)
                c -= h * (k1c + 2 * k2c + 2 * k3c + k4c) / 6
                consumed += h * (k1n + 2 * k2n + 2 * k3n + k4n) / 6
            if j < len(stops) - 2:  # a pulse at s1: its N mass into V(s1)
                step = dose * v0 / (v0 - e * (s1 - t0))
                jumps.append((s1, c, c + step))
                c += step
        cs[i] = c
        rs[i] = deriv(tb, c, fb, consumed)[2]
    return Simulation(np.array(cs), np.array(rs), consumed, jumps)


# --- Parameters ----------------------------------------------------------------------------
@dataclass
class Layout:
    names: list[str]
    lo: np.ndarray
    hi: np.ndarray
    log: list[bool]


def layout(spec: ModelSpec, data: ModelData) -> Layout:
    names, lo, hi, log = [], [], [], []
    for k in spec.kinetic:
        names.append(k)
        lo.append(LOG_KIN_BOUNDS[0])
        hi.append(LOG_KIN_BOUNDS[1])
        log.append(True)
    for name, a, b, lg in (("C0_mgN_L", math.log(1.0), math.log(300.0), True), ("a", -5.0, 5.0, False),
                           ("b0", 0.02, 3.0, False), ("b1_per_d", -0.1, 0.1, False),
                           ("tau_d", math.log(1e-3), math.log(3.0), True)):
        names.append(name)
        lo.append(a)
        hi.append(b)
        log.append(lg)
    for t in data.offsets:
        names.append(f"offset_log_{t:.3f}d")
        lo.append(-3.0)
        hi.append(3.0)
        log.append(False)
    return Layout(names, np.array(lo), np.array(hi), log)


def initial_guess(spec: ModelSpec, data: ModelData, lay: Layout) -> np.ndarray:
    g = {"v0_mgN_L_d": 5.0, "vmax_mgN_L_d": 5.0, "kd_per_d": 0.05, "Ks_mgNH3_L": 0.5, "Ki_mgNH3_L": 1.0,
         "Alk0_meq_L": 10.0}
    # Y6: the (α, Ks, S_opt) start is the (vmax, Ks, Ki) default start transformed (α = 10, S_opt = 0.707 mg NH3/L)
    g[HALDANE_ALPHA[0]] = g[HALDANE_ALPHA[1]] = g["vmax_mgN_L_d"] / g["Ks_mgNH3_L"]
    g["Sopt_mgNH3_L"] = math.sqrt(g["Ks_mgNH3_L"] * g["Ki_mgNH3_L"])
    x = [math.log(g[k]) for k in spec.kinetic]
    x += [float(data.log_r[0]), 0.0, 1.0, 0.0, math.log(0.05)]
    x += list(data.offset_guess)
    return np.clip(np.array(x), lay.lo + 1e-9, lay.hi - 1e-9)


def unpack(spec: ModelSpec, data: ModelData, lay: Layout, x: np.ndarray) -> dict:
    v = {n: (math.exp(xi) if lg else xi) for n, xi, lg in zip(lay.names, x, lay.log)}
    v["kin"] = {k: v[k] for k in spec.kinetic}
    v["offset_values"] = [v[f"offset_log_{t:.3f}d"] for t in data.offsets]
    return v


def observe(c: np.ndarray, t: np.ndarray, p: dict, offsets: list[float]) -> np.ndarray:
    """log of the lagged electrode reading (mg N/L)."""
    u = p["a"] + (p["b0"] + p["b1_per_d"] * (t - T0_D)) * np.log(np.maximum(c, 1e-3))
    for tj, dj in zip(offsets, p["offset_values"]):
        u = u + dj * (t >= tj)
    r = np.exp(u)
    alpha = math.exp(-(t[1] - t[0]) / p["tau_d"])
    lagged, _ = lfilter([1 - alpha], [1, -alpha], r, zi=[alpha * r[0]])
    return np.log(lagged)


def invert_observation(log_r: np.ndarray, t: np.ndarray, p: dict, offsets: list[float]) -> np.ndarray:
    """True NH4 C(t) (mg N/L) from the electrode reading under the M6 observation model (lag not inverted):
    C = exp((log R − a − Σ offsets)/b(t))."""
    b = p["b0"] + p["b1_per_d"] * (np.asarray(t, dtype=float) - T0_D)
    if np.any(b <= 0):
        raise ValueError("the observation-model gain b(t) must stay positive to be inverted")
    u = np.asarray(log_r, dtype=float) - p["a"]
    for tj, dj in zip(offsets, p["offset_values"]):
        u = u - dj * (t >= tj)
    return np.exp(u / b)


def predict(spec: ModelSpec, data: ModelData, lay: Layout, x: np.ndarray) -> tuple[np.ndarray, Simulation]:
    p = unpack(spec, data, lay, x)
    if data.volume is None:
        sim = simulate(spec, p["kin"], p["C0_mgN_L"], data.t, data.ph, data.pulses, data.dose)
    else:
        sim = simulate_volume(spec, p["kin"], p["C0_mgN_L"], data.t, data.ph, data.pulses, data.dose, data.volume)
    return observe(sim.c, data.t, p, data.offsets), sim


# --- Fitting ----------------------------------------------------------------------------------
@dataclass
class Fit:
    spec: ModelSpec
    x: np.ndarray
    free: np.ndarray
    window: tuple[float, float]
    rss: float
    n: int
    converged: bool
    status: int = 1  # least_squares status (0: evaluation budget exhausted; > 0: a tolerance was met)
    nfev: int = 0  # function evaluations (all starts, for a calibration under the convergence policy)
    start: str = ""  # the start of the convergence policy that gave this fit


def fit(spec: ModelSpec, data: ModelData, lay: Layout, window: tuple[float, float], x0: np.ndarray,
        free: np.ndarray, max_nfev: int | None = None) -> Fit:
    """Least squares in log space over bins in ``window``; parameters not in ``free`` stay at x0.

    The state model and the lag filter are causal, so the record after the window cannot change a residual inside it:
    the simulation stops at the window end (bit-identical residuals, half the cost for the calibration window)."""
    data = truncated(data, window[1])
    mask = (data.t >= window[0]) & (data.t < window[1])
    x0 = np.clip(x0.copy(), lay.lo + 1e-9, lay.hi - 1e-9)

    def resid(z):
        x = x0.copy()
        x[free] = z
        return (predict(spec, data, lay, x)[0] - data.log_r)[mask]

    if free.sum() == 0:
        r = resid(np.array([]))
        return Fit(spec, x0, free, window, float(r @ r), int(mask.sum()), True)
    sol = least_squares(resid, x0[free], bounds=(lay.lo[free], lay.hi[free]), x_scale="jac", method="trf",
                        max_nfev=max_nfev)
    x = x0.copy()
    x[free] = sol.x
    return Fit(spec, x, free, window, float(sol.fun @ sol.fun), int(mask.sum()), bool(sol.success), int(sol.status),
               int(sol.nfev))


def truncated(data: ModelData, t_end: float) -> ModelData:
    """The bins before ``t_end`` (the same object when none is cut)."""
    keep = data.t < t_end
    if keep.all():
        return data
    return ModelData(t=data.t[keep], log_r=data.log_r[keep], ph=data.ph[keep], pulses=data.pulses, dose=data.dose,
                     offsets=data.offsets, offset_guess=data.offset_guess, volume=data.volume)


def window_skill(spec, data, lay, x, window) -> dict:
    mask = (data.t >= window[0]) & (data.t < window[1])
    lr, _ = predict(spec, data, lay, x)
    e = (lr - data.log_r)[mask]
    lin = (np.exp(lr) - np.exp(data.log_r))[mask]
    return {"rmse_log": float(np.sqrt(np.mean(e**2))), "rmse_mgN_L": float(np.sqrt(np.mean(lin**2))),
            "n": int(mask.sum())}


def calibration_starts(spec: ModelSpec, lay: Layout, x0: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """X1: the starts of the convergence policy (the kinetic parameters are on the log scale, so a factor is a shift)."""
    if not spec.kinetic:
        return [("default (no kinetic parameter)", x0)]
    out = []
    for f in CAL_START_FACTORS:
        x = x0.copy()
        x[:len(spec.kinetic)] += math.log(f)
        out.append((f"kinetic start x{f:g}", np.clip(x, lay.lo + 1e-9, lay.hi - 1e-9)))
    return out


BOUND_TOL_FRACTION = 1e-3  # Z2: a free parameter within 1e-3 x its bound width of a bound (fitted coordinates), inclusive


def at_bounds(lay: Layout, x: np.ndarray, free: np.ndarray) -> str:
    """Z2: the free parameters of a fit at a bound (``BOUND_TOL_FRACTION``; docs/progress.md, pass 5), or ''."""
    tol = BOUND_TOL_FRACTION * (lay.hi - lay.lo)
    return "; ".join(n for n, xi, lo, hi, t, f in zip(lay.names, x, lay.lo, lay.hi, tol, free)
                     if f and (xi - lo <= t or hi - xi <= t))


def bounds_text(fit_df: pd.DataFrame) -> str:
    """Z2: every candidate whose calibration is bound-constrained (`cal_at_bounds`), for §4 and ``model.md``."""
    flags = fit_df["cal_at_bounds"].fillna("").astype(str)
    rule = (f"`cal_at_bounds`: a free parameter within {BOUND_TOL_FRACTION:g} x its bound width of a bound, in the "
            "fitted coordinates")
    if not (flags != "").any():
        return f"No calibration is bound-constrained ({rule})."
    cand = fit_df["model"] != NULL_MODEL
    named = "; ".join(f"{m} ({f.replace('; ', ', ')})" for m, f in zip(fit_df["model"][cand], flags[cand]) if f)
    ref = flags[~cand].iloc[0] if (~cand).any() else ""
    sel = fit_df.loc[fit_df["selected"].astype(bool), "model"].iloc[0] if "selected" in fit_df else None
    return ((f"Bound-constrained calibrations ({rule}): {named}. A bound, not the data, fixes at least one parameter "
             "of these calibrations, so their natural-unit kinetic parameters are not estimates; their skills are "
             "what the comparison uses, and no bound was changed to clear a flag." if named else
             f"No candidate calibration is bound-constrained ({rule}).")
            + (f" The reference {NULL_MODEL} (no kinetic parameter) has {ref.replace('; ', ', ')} at a bound."
               if ref else "")
            + (f" The identifiability of the selected model ({sel}) is that of its whole-record refit (profile "
               "likelihood), not of this calibration." if sel is not None and flags[fit_df["model"] == sel].ne("").any()
               else ""))


def _calibration_job(args) -> Fit:
    spec, data, lay, window, x0, free, start = args
    f = fit(spec, data, lay, window, x0, free, max_nfev=CAL_MAX_NFEV)
    f.start = start
    return f


def best_start(fits: list[Fit]) -> Fit:
    """The converged start with the lowest RSS (ties: the earlier start); else the lowest-RSS start, unconverged.
    ``nfev`` is the total over the starts."""
    conv = [f for f in fits if f.converged]
    best = min(conv or fits, key=lambda f: f.rss)  # min keeps the first of equal values
    return replace(best, nfev=sum(f.nfev for f in fits))


def _map(fn, tasks: list) -> list:
    if N_WORKERS > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(tasks))) as pool:
            return list(pool.map(fn, tasks))
    return [fn(t) for t in tasks]


def calibrate_all(specs: list[ModelSpec], data: ModelData, window: tuple[float, float],
                  setup: dict[str, tuple[Layout, np.ndarray, np.ndarray]]) -> tuple[dict[str, Fit], dict[str, list[Fit]]]:
    """Calibrate every candidate under the convergence policy; every (candidate, start) pair runs in parallel.
    Returns the calibration of each candidate and every start's fit."""
    tasks, owner = [], []
    for spec in specs:
        lay, x0, free = setup[spec.name]
        for label, xs in calibration_starts(spec, lay, x0):
            tasks.append((spec, data, lay, window, xs, free, label))
            owner.append(spec.name)
    done = _map(_calibration_job, tasks)
    starts = {spec.name: [f for o, f in zip(owner, done) if o == spec.name] for spec in specs}
    return {name: best_start(fs) for name, fs in starts.items()}, starts


def split_evaluation(cfg: Config, data: ModelData, windows: SplitWindows | None = None) -> tuple[pd.DataFrame, dict]:
    """Calibrate on CAL, evaluate frozen on SEL, then TEST (unseen offsets as the only nuisance)."""
    w = split_windows(cfg) if windows is None else windows
    CAL, SEL, TEST = w.cal, w.sel, w.test
    rows, fits, setup = [], {}, {}
    specs = model_specs(cfg)
    for spec in specs:
        lay = layout(spec, data)
        x0 = initial_guess(spec, data, lay)
        free = np.ones(len(x0), dtype=bool)
        for i, t in enumerate(data.offsets):
            if t >= CAL[1]:
                free[len(lay.names) - len(data.offsets) + i] = False
                x0[len(lay.names) - len(data.offsets) + i] = 0.0
        setup[spec.name] = (lay, x0, free)
    cals, all_starts = calibrate_all(specs, data, CAL, setup)
    # Y1: out-of-sample skills of every start (DIAGNOSTIC; never used to choose a calibration or a model)
    oos = iter(_map(_start_skill_job, [(spec, data, setup[spec.name][0], f.x, TEST) for spec in specs
                                       for f in all_starts[spec.name]]))
    start_rows = []
    for spec in specs:
        lay = setup[spec.name][0]
        cal = cals[spec.name]
        n_eff = calibration_n_eff(spec, data, lay, cal)
        for f in all_starts[spec.name]:
            d = start_delta(f, cal, n_eff)
            start_rows.append({"model": spec.name, "start": f.start, "converged": f.converged, "status": f.status,
                               "nfev": f.nfev, "cal_rmse_log": math.sqrt(f.rss / f.n),
                               "alk0_meq_L": unpack(spec, data, lay, f.x)["Alk0_meq_L"] if spec.alk else np.nan,
                               "sel_rmse_log_diagnostic": window_skill(spec, data, lay, f.x, SEL)["rmse_log"],
                               "chosen": f.start == cals[spec.name].start, "delta_vs_calibration": d,
                               "inside_95": bool(cal.converged and f.converged and d < CHI2_95_1DOF),
                               **next(oos), "cal_at_bounds": at_bounds(lay, f.x, f.free)})
    bounds = []
    for spec in specs:
        lay, x0, free = setup[spec.name]
        cal = cals[spec.name]
        sel = window_skill(spec, data, lay, cal.x, SEL)
        test, forecast, unseen, fc_window = out_of_sample(spec, data, lay, cal.x, TEST)
        rows.append({"model": spec.name, "n_params_calibration": int(free.sum()),
                     "cal_rmse_log": math.sqrt(cal.rss / cal.n), "cal_converged": cal.converged,
                     "cal_status": cal.status, "cal_nfev": cal.nfev, "cal_start": cal.start,
                     "sel_rmse_log": sel["rmse_log"], "sel_rmse_mgN_L": sel["rmse_mgN_L"],
                     "test_rmse_log": test["rmse_log"], "test_rmse_mgN_L": test["rmse_mgN_L"],
                     "test_nuisance": f"offsets at {', '.join(f'{t:.3f}' for t in unseen)} d refitted on the test window; "
                                      "all other parameters frozen",
                     "forecast_window_d": f"{fc_window[0]:.2f} to {fc_window[1]:.3f}",
                     "forecast_rmse_log": forecast["rmse_log"], "forecast_rmse_mgN_L": forecast["rmse_mgN_L"]})
        fits[spec.name] = (spec, lay, cal)
        bounds.append(at_bounds(lay, cal.x, cal.free))
    df = pd.DataFrame(rows)
    null = df.set_index("model").loc["null_no_consumption"]
    for col in ("sel", "test", "forecast"):
        df[f"{col}_skill_vs_null"] = 1 - df[f"{col}_rmse_log"] ** 2 / null[f"{col}_rmse_log"] ** 2
    df = select_converged(df)
    alk = pd.DataFrame(_map(_alk0_job, [(spec, data, lay, cal, all_starts[spec.name]) for spec, lay, cal in fits.values()]))
    st = pd.DataFrame(start_rows)
    nul = df.set_index("model").loc[NULL_MODEL]
    for col in ("sel", "test", "forecast"):
        st[f"{col}_skill_vs_null_diagnostic"] = 1 - st[f"{col}_rmse_log_diagnostic"] ** 2 / float(nul[f"{col}_rmse_log"]) ** 2
    st["role"] = ("DIAGNOSTIC: every start of the convergence policy; the calibration is the converged start with the "
                  "lowest RSS (chosen); inside_95 = converged and Δ = n_eff ln(RSS/RSS_calibration) < 3.84; the "
                  "selection, test and forecast skills of the starts are reported as ranges and never used to choose")
    out = pd.concat([df, alk, pd.DataFrame([skill_range(st[st["model"] == m]) for m in df["model"]])], axis=1)
    out["cal_at_bounds"] = bounds  # Z2 (rows in the order of ``specs``, as ``df``)
    out.attrs["starts"] = st
    return out, fits


def out_of_sample(spec, data, lay, x, test_window) -> tuple[dict, dict, list[float], tuple[float, float]]:
    """Test skill (offsets unseen in calibration refitted on the test window, all else frozen at ``x``) and pure
    forecast skill (nothing refitted, up to the first unseen offset)."""
    test_off = np.array([n.startswith("offset_log_") and float(n[11:-1]) >= test_window[0] for n in lay.names])
    test_fit = fit(spec, data, lay, test_window, x, test_off)
    test = window_skill(spec, data, lay, test_fit.x, test_window)
    unseen = [t for t in data.offsets if t >= test_window[0]]
    fc_window = (test_window[0], min(unseen) if unseen else test_window[1])
    return test, window_skill(spec, data, lay, x, fc_window), unseen, fc_window


def _start_skill_job(args) -> dict:
    test, forecast, _, _ = out_of_sample(*args)
    return {"test_rmse_log_diagnostic": test["rmse_log"], "forecast_rmse_log_diagnostic": forecast["rmse_log"]}


def calibration_n_eff(spec, data, lay, cal: Fit) -> float:
    """n_eff of the residuals of a calibration on its window (the n_eff of its profile likelihood)."""
    mask = (data.t >= cal.window[0]) & (data.t < cal.window[1])
    return effective_n((predict(spec, data, lay, cal.x)[0] - data.log_r)[mask])[0]


def start_delta(f: Fit, cal: Fit, n_eff: float) -> float:
    """Δ = n_eff ln(RSS_start/RSS_calibration) of one start against the chosen calibration."""
    return n_eff * math.log(f.rss / cal.rss) if f.rss > 0 and cal.rss > 0 else 0.0


SKILL_RANGE_COLS = [f"{c}_skill_{e}_inside_95" for c in ("sel", "test", "forecast") for e in ("min", "max")]


def skill_range(st: pd.DataFrame) -> dict:
    """Y1: range of the selection, test and forecast skills over a candidate's converged starts inside Δ < 3.84 of its
    calibration (DIAGNOSTIC; the calibration stays the lowest-RSS converged start)."""
    inside = st[st["inside_95"].astype(bool)]
    out = {"n_starts_inside_95": len(inside), **dict.fromkeys(SKILL_RANGE_COLS, np.nan)}
    if not st["converged"].astype(bool).any():
        out["skill_range_status"] = "not computed: calibration did not converge (no start is an estimate)"
        return out
    for c in ("sel", "test", "forecast"):
        v = inside[f"{c}_skill_vs_null_diagnostic"]
        out[f"{c}_skill_min_inside_95"], out[f"{c}_skill_max_inside_95"] = float(v.min()), float(v.max())
    out["skill_range_status"] = (f"{len(inside)} converged starts inside Δ < 3.84 of the calibration: DIAGNOSTIC range, "
                                 "never used to choose" if len(inside) > 1 else
                                 "one converged start inside Δ < 3.84 (the calibration): no range")
    return out


def _alk0_job(args) -> dict:
    return alk0_profile(*args)


NULL_MODEL = "null_no_consumption"
LSQ_STATUS = {-1: "improper input", 0: "maximum number of function evaluations exceeded", 1: "gtol satisfied",
              2: "ftol satisfied", 3: "xtol satisfied", 4: "ftol and xtol satisfied"}


class SelectionError(RuntimeError):
    """No candidate model has a converged calibration, so none can be selected."""


def select_converged(df: pd.DataFrame) -> pd.DataFrame:
    """X1: select by minimum selection RMSE among the candidates whose calibration converged; an unconverged candidate
    keeps its skills in the table, labelled by ``selection_status``, and cannot be selected."""
    out = df.copy()
    conv = out["cal_converged"].astype(bool)
    cand = out[(out["model"] != NULL_MODEL) & conv]
    if cand.empty:
        raise SelectionError("no candidate calibration converged under the convergence policy")
    selected = cand.loc[cand["sel_rmse_log"].idxmin(), "model"]
    out["selected"] = out["model"] == selected
    out["selection_status"] = [
        "reference (not a candidate)" if m == NULL_MODEL else
        f"not computed: calibration did not converge (status {st}: {LSQ_STATUS.get(int(st), '?')}); skills not "
        "comparable" if not c else "selected" if m == selected else "candidate (converged, not selected)"
        for m, c, st in zip(out["model"], conv, out["cal_status"])]
    return out


# --- Identifiability -------------------------------------------------------------------------
def effective_n(resid: np.ndarray) -> tuple[float, float]:
    rho = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])
    rho = min(max(rho, 0.0), 0.999)
    return len(resid) * (1 - rho) / (1 + rho), rho


def profile_interval(delta, x_best: float, half: float, lo: float, hi: float,
                     n_points: int = PROFILE_POINTS_PER_SIDE, target: float = CHI2_95_1DOF,
                     xtol: float | None = None) -> dict:
    """Profile-likelihood interval of one parameter by root finding on the profile.

    ``delta(v)`` returns (Δ(v), converged) for the parameter fixed at v and the others re-optimized. On
    each side a grid of ``n_points`` steps up to ``half`` brackets the first converged Δ ≥ ``target``; the
    crossing Δ(v) = target is then found with Brent's method inside that bracket (endpoint values are
    reused, never re-evaluated). A side is closed only by a converged Δ ≥ target within the bounds.
    """
    xtol = PROFILE_XTOL_FRACTION * half if xtol is None else xtol
    out = {}
    for side, sign in (("low", -1), ("high", 1)):
        cache = {x_best: 0.0}
        pts = []
        prev = x_best
        closed, crossing = False, np.nan
        for k in range(1, n_points + 1):
            val = x_best + sign * half * k / n_points
            at_bound = val <= lo or val >= hi
            val = float(np.clip(val, lo + 1e-9, hi - 1e-9))
            d, ok = delta(val)
            cache[val] = d
            pts.append((val, d, ok))
            if d >= target and ok:
                def g(z, cache=cache):
                    if z not in cache:
                        cache[z] = delta(z)[0]
                    return cache[z] - target
                a, b = sorted((prev, val))
                crossing = float(brentq(g, a, b, xtol=xtol, maxiter=PROFILE_MAXITER))
                closed = True
                break
            if d >= target or at_bound:
                break
            prev = val
        out[side] = {"closed": closed, "crossing": crossing, "points": pts}
    return out


def profile_likelihood(spec, data, lay, best: Fit, collect: list | None = None,
                       params: tuple[str, ...] | None = None) -> pd.DataFrame:
    """Profile Δ(θ_i) = n_eff ln(RSS(θ_i)/RSS_min); identifiable if Δ ≥ 3.84 is reached on both sides within
    the bounds; the 95 % interval ends are the roots of Δ(θ_i) = 3.84 (:func:`profile_interval`). ``params`` restricts
    the profiles to those parameters; parameters that ``best`` kept fixed stay fixed.

    ``collect``, if given, receives (parameter, value, Δ, converged, x) for every profile evaluation (W1 uses the
    converged vectors with Δ ≤ 3.84 to propagate parameter uncertainty into the inverted C(t))."""
    mask = (data.t >= best.window[0]) & (data.t < best.window[1])
    r_best = (predict(spec, data, lay, best.x)[0] - data.log_r)[mask]
    n_eff, rho = effective_n(r_best)
    # approximate standard errors to size the grid (n_eff-scaled Jacobian)
    eps = 1e-5
    jac = np.empty((mask.sum(), len(best.x)))
    for i in range(len(best.x)):
        xp = best.x.copy()
        xp[i] += eps
        jac[:, i] = ((predict(spec, data, lay, xp)[0] - data.log_r)[mask] - r_best) / eps
    s2 = best.rss / max(n_eff - len(best.x), 1.0)
    se = np.sqrt(np.clip(np.diag(np.linalg.pinv(jac.T @ jac)) * s2, 0, None))
    rows = []
    n_refits = [0]
    n_evals = [0]
    for i, name in enumerate(lay.names):
        if params is not None and name not in params:
            continue
        half = max(4 * se[i], 1.0 if lay.log[i] else 0.2 if name != "b1_per_d" else 0.01)
        free = best.free.copy()
        free[i] = False
        solved = [(best.x[i], best.x.copy())]  # (value, solution) for warm starts

        def delta(val, i=i, free=free, solved=solved):
            x0 = min(solved, key=lambda vs: abs(vs[0] - val))[1].copy()
            x0[i] = val
            f = fit(spec, data, lay, best.window, x0, free, max_nfev=PROFILE_MAX_NFEV)
            d = n_eff * math.log(f.rss / best.rss) if f.rss > 0 else 0.0
            if not f.converged and d >= CHI2_95_1DOF - 1.0:
                # an unconverged re-optimization overstates Δ near the crossing: refit without an evaluation cap
                f = fit(spec, data, lay, best.window, f.x, free)
                d = n_eff * math.log(f.rss / best.rss) if f.rss > 0 else 0.0
                n_refits[0] += 1
            n_evals[0] += 1
            solved.append((val, f.x))
            if collect is not None:
                collect.append((lay.names[i], float(val), float(d), bool(f.converged), f.x.copy()))
            return d, f.converged

        prof = profile_interval(delta, float(best.x[i]), half, float(lay.lo[i]), float(lay.hi[i]))
        closed = {s: prof[s]["closed"] for s in prof}
        all_pts = [d for side in prof.values() for _, d, _ in side["points"]]
        min_delta = min(all_pts, default=0.0)
        to_nat = (lambda z: math.exp(z)) if lay.log[i] else (lambda z: z)
        identifiable = closed["low"] and closed["high"]
        rows.append({
            "parameter": name, "scale": "log" if lay.log[i] else "linear",
            "estimate": to_nat(best.x[i]) if identifiable else "not identifiable",
            "ci95_low": to_nat(prof["low"]["crossing"]) if identifiable else np.nan,
            "ci95_high": to_nat(prof["high"]["crossing"]) if identifiable else np.nan,
            "identifiable": identifiable,
            "reason": "" if identifiable else
                      "profile stays below 3.84 on the " + " and ".join(s for s in ("low", "high") if not closed[s])
                      + " side(s) within the bounds or grid",
            "optimizer_value_not_an_estimate": np.nan if identifiable else to_nat(best.x[i]),
            "profile_max_delta_low": max((d for _, d, _ in prof["low"]["points"]), default=np.nan),
            "profile_max_delta_high": max((d for _, d, _ in prof["high"]["points"]), default=np.nan),
            "profile_min_delta": min_delta,
            # one-sided 95 % limits: the crossing of a side that closed, also when the other side stays open
            "profile_crossing_low": to_nat(prof["low"]["crossing"]) if closed["low"] else np.nan,
            "profile_crossing_high": to_nat(prof["high"]["crossing"]) if closed["high"] else np.nan,
            "better_optimum_found": bool(min_delta < -0.5),
            "n_eff": n_eff, "resid_lag1": rho,
        })
    out = pd.DataFrame(rows)
    out.attrs["uncapped_refits"] = n_refits[0]
    out.attrs["profile_evaluations"] = n_evals[0]
    return out


def distinct_minima(points: list[tuple[float, float]], tol: float = MULTIMODAL_LOG_TOL) -> list[tuple[float, float, int]]:
    """Y2: (Alk0, Δ) of converged starts grouped into distinct minima: in order of increasing Δ, a start within ``tol``
    (log scale) of an existing minimum joins it, else it starts a new one. Returns (Alk0, Δ, number of starts) of each
    minimum's lowest-Δ start."""
    out: list[list] = []
    for a, d in sorted(points, key=lambda p: p[1]):
        for m in out:
            if abs(math.log(a / m[0])) <= tol:
                m[2] += 1
                break
        else:
            out.append([a, d, 1])
    return [tuple(m) for m in out]


def alk0_profile(spec: ModelSpec, data: ModelData, lay: Layout, cal: Fit, others: list[Fit] | None = None) -> dict:
    """X1 (d): fitted Alk0 of an alkalinity candidate (its calibration) and the 95 % profile interval of Alk0 on the
    calibration window, or ``not identifiable``; empty for the other candidates.

    The profile is warm-started from neighbouring points, so it follows one minimum. Any other converged start of the
    policy (``others``) whose Alk0 differs and whose Δ = n_eff ln(RSS/RSS_min) is below 3.84 lies inside the 95 %
    region as well; the region is then not one interval and no limit is reported (a validity condition). Starts within
    ``MULTIMODAL_LOG_TOL`` of each other are one minimum (:func:`distinct_minima`), listed once with its number of
    starts; ``alk0_converged_minima_inside_95_meq_L`` lists every distinct converged minimum inside Δ < 3.84."""
    out = {"alk0_fitted_meq_L": np.nan, "alk0_ci95_low_meq_L": np.nan, "alk0_ci95_high_meq_L": np.nan,
           "alk0_profile_status": "", "alk0_converged_minima_inside_95_meq_L": ""}
    if not spec.alk:
        return out
    out["alk0_fitted_meq_L"] = unpack(spec, data, lay, cal.x)["Alk0_meq_L"]
    if not cal.converged:
        out["alk0_profile_status"] = "not computed: calibration did not converge (the value is not an estimate)"
        return out
    r = profile_likelihood(spec, data, lay, cal, params=("Alk0_meq_L",)).iloc[0]
    a_best = out["alk0_fitted_meq_L"]
    region = []  # Y2: every converged start inside Δ < 3.84 (``others`` contains the calibration's own start)
    for f in others or []:
        a_f = unpack(spec, data, lay, f.x)["Alk0_meq_L"]
        d = float(r["n_eff"]) * math.log(f.rss / cal.rss) if f.rss > 0 else 0.0
        if f.converged and d < CHI2_95_1DOF:
            region.append((a_f, d))
    if not any(abs(math.log(a / a_best)) <= MULTIMODAL_LOG_TOL for a, _ in region):
        region.append((a_best, 0.0))
    minima = distinct_minima(region)
    out["alk0_converged_minima_inside_95_meq_L"] = "; ".join(f"{a:.1f}" for a, _, _ in sorted(minima))
    inside = [(a, d, k) for a, d, k in minima if abs(math.log(a / a_best)) > MULTIMODAL_LOG_TOL]
    if inside:
        mins = "; ".join(f"Alk0 {a:.1f} meq/L at Δ = {d:.2f} ({k} start{'s' if k > 1 else ''})" for a, d, k in sorted(inside))
        out["alk0_profile_status"] = (f"not identifiable: the profile has more than one minimum below Δ = 3.84 (fitted "
                                      f"{a_best:.1f} meq/L; other converged minima: {mins}); no interval"
                                      + ("" if r["identifiable"] else f"; local profile: {r['reason']}"))
        return out
    if r["identifiable"]:
        out.update(alk0_ci95_low_meq_L=r["ci95_low"], alk0_ci95_high_meq_L=r["ci95_high"], alk0_profile_status="identifiable")
        return out
    out.update(alk0_ci95_low_meq_L=r["profile_crossing_low"], alk0_ci95_high_meq_L=r["profile_crossing_high"])
    one = [f"one-sided 95 % {side} limit {v:.1f} meq/L" for side, v in (("lower", r["profile_crossing_low"]),
                                                                    ("upper", r["profile_crossing_high"])) if np.isfinite(v)]
    out["alk0_profile_status"] = f"not identifiable: {r['reason']}" + (f"; {', '.join(one)}" if one else "")
    return out


# --- Consistency checks ---------------------------------------------------------------------------
def consistency_checks(cfg: Config, spec, data, lay, best: Fit, coexistence_affinity_kJ_mol: float) -> pd.DataFrame:
    _, sim = predict(spec, data, lay, best.x)
    p = unpack(spec, data, lay, best.x)
    rows = []
    for t, before, after in sim.pulse_jumps:
        rows.append({"check": "pulse adds exactly the dose", "where": f"{t:g} d", "value": after - before,
                     "expected": data.dose, "pass": abs(after - before - data.dose) < 1e-9})
    closure = sim.c[-1] - (p["C0_mgN_L"] + data.dose * len(sim.pulse_jumps) - sim.consumed)
    rows.append({"check": "mass balance closes: C_end − (C0 + pulses − consumed)", "where": "whole record",
                 "value": closure, "expected": 0.0, "pass": abs(closure) < 1e-9})
    rows.append({"check": "consumption rate non-negative", "where": "every grid time", "value": float(sim.rate.min()),
                 "expected": ">= 0", "pass": bool(sim.rate.min() >= 0)})
    T_K = cfg.get("scenarios.coexistence.temperature_C") + KELVIN_OFFSET
    sigma = sim.rate / M_N_G_MOL / 1e3 / SECONDS_PER_DAY * coexistence_affinity_kJ_mol * 1e3 / T_K  # W/(L K)
    rows.append({"check": "σ = r A / T ≥ 0 with the coexistence SCENARIO affinity", "where": "every grid time",
                 "value": float(sigma.min()), "expected": ">= 0", "pass": bool(sigma.min() >= 0)})
    # In situ affinity from the simulated state and the observed pH (guarded on DO and NO3-N); 25 °C SCENARIO
    a_insitu = thermo.nitrification_affinity_guarded(cfg, np.maximum(sim.c, 1e-9), data.ph)
    if isinstance(a_insitu, thermo.NotComputed):
        rows.append({"check": "σ in situ", "where": "every grid time (model state C, observed pH)", "value": str(a_insitu),
                     "expected": ">= 0", "pass": "not computed"})
    else:
        sig = sim.rate / M_N_G_MOL / 1e3 / SECONDS_PER_DAY * np.asarray(a_insitu) * 1e3 / thermo.T_REF_K
        rows.append({"check": "σ in situ (SCENARIO 25 °C)", "where": "every grid time (model state C, observed pH)",
                     "value": float(sig.min()), "expected": ">= 0", "pass": bool(sig.min() >= 0)})
    return pd.DataFrame(rows)


def response_fractions(spec, data, lay, best: Fit, pulse_table: pd.DataFrame, n_factor: float) -> pd.DataFrame:
    """Model-implied electrode response to each pulse with the M2 window, against the observed fraction."""
    lr, _ = predict(spec, data, lay, best.x)
    r = np.exp(lr) / n_factor  # back to mg NH4+/L
    rows = []
    for k, tp in enumerate(data.pulses, start=1):
        pre = np.median(r[(data.t >= tp - 0.1) & (data.t < tp)])
        post = np.median(r[(data.t >= tp + 0.4) & (data.t < tp + 0.6)])
        obs = float(pulse_table.set_index("pulse").loc[k, "nh4_response_fraction"])
        rows.append({"pulse": k, "model_response_fraction": (post - pre) / (data.dose / n_factor),
                     "observed_response_fraction": obs})
    return pd.DataFrame(rows)


# --- W1: inverted NH4 truth with propagated parameter uncertainty ----------------------------------
def vector_row(spec, data, lay, x: np.ndarray) -> dict:
    """Natural-scale parameter values of one parameter vector (log parameters exponentiated)."""
    p = unpack(spec, data, lay, x)
    return {n: p[n] for n in lay.names}


def params_from_row(row) -> tuple[dict, list[float]]:
    """Observation-model parameters and offset times from a row of ``nh4_truth_vectors.csv``."""
    names = [k for k in row.index if str(k).startswith("offset_log_")]
    offsets = [float(str(k)[len("offset_log_"):-1]) for k in names]
    p = {"a": float(row["a"]), "b0": float(row["b0"]), "b1_per_d": float(row["b1_per_d"]),
         "offset_values": [float(row[k]) for k in names]}
    return p, offsets


def invert_from_row(row, t: np.ndarray, reading_mgN_L: np.ndarray) -> np.ndarray:
    p, offsets = params_from_row(row)
    return invert_observation(np.log(reading_mgN_L), t, p, offsets)


@dataclass
class Truth:
    series: pd.DataFrame  # t_d, reading, pH, C_inv, band, provenance
    vectors: pd.DataFrame  # parameter vectors whose inversions form the band (vector 0 = best)
    summary: pd.DataFrame  # key, value


def nh4_truth(spec, data, lay, best: Fit, ident: pd.DataFrame, points: list) -> Truth:
    """C(t) from the whole-record refit, and its band over every converged profile vector with Δ ≤ 3.84.

    The band is an envelope of 95 % profile vectors (a lower bound on the joint uncertainty, not a joint region).
    A side is open when a non-identifiable parameter's open profile side moves C(t) that way: the data then do not
    bound C(t) on that side and the band edge there is only the last evaluated vector.
    """
    c = invert_observation(data.log_r, data.t, unpack(spec, data, lay, best.x), data.offsets)
    kept = [(name, val, d, x) for name, val, d, ok, x in points if ok and d <= CHI2_95_1DOF]
    vecs = [("best (whole-record refit)", np.nan, 0.0, best.x)] + kept
    cs = np.array([invert_observation(data.log_r, data.t, unpack(spec, data, lay, x), data.offsets)
                   for *_, x in vecs])
    lo, hi = cs.min(axis=0), cs.max(axis=0)
    open_hi = open_lo = False
    reasons = []
    for r in ident[~ident["identifiable"].astype(bool)].itertuples(index=False):
        i = lay.names.index(r.parameter)
        for side, sign, dmax in (("low", -1, r.profile_max_delta_low), ("high", 1, r.profile_max_delta_high)):
            if not dmax < CHI2_95_1DOF:
                continue
            pts = [(val, x) for name, val, d, ok, x in points if name == r.parameter and sign * (val - best.x[i]) > 0]
            if not pts:
                continue
            _, x_far = max(pts, key=lambda vx: abs(vx[0] - best.x[i]))
            c_far = invert_observation(data.log_r, data.t, unpack(spec, data, lay, x_far), data.offsets)
            up = float(np.median(c_far / c)) > 1.0
            open_hi |= up
            open_lo |= not up
            reasons.append(f"parameter {r.parameter}: profile open on the {side} side moves C(t) {'up' if up else 'down'}")
    tid = truth_id(c)
    series = pd.DataFrame({
        "t_d": data.t, "electrode_obs_mgN_L": np.exp(data.log_r), "pH": data.ph, "C_inv_mgN_L": c,
        "C_inv_band_low_mgN_L": lo, "C_inv_band_high_mgN_L": hi,
        "gain_b": unpack(spec, data, lay, best.x)["b0"] + unpack(spec, data, lay, best.x)["b1_per_d"] * (data.t - T0_D),
        "nh4_observation_model": OBSERVATION_MODEL, "nh4_truth_id": tid})
    vectors = pd.DataFrame([{"vector": k, "profiled_parameter": name.split(" ")[0] if k else "best",
                             "profiled_value_internal": val, "delta": d, **vector_row(spec, data, lay, x)}
                            for k, (name, val, d, x) in enumerate(vecs)])
    t0, t1 = data.t < T0_D + 0.1, data.t >= END_D - 0.1
    summary = pd.DataFrame([
        ("nh4_truth_id", tid), ("nh4_observation_model", OBSERVATION_MODEL), ("model", spec.name),
        ("inversion", INVERSION_NOTE), ("n_vectors_in_band", len(vecs)),
        ("band_basis", "envelope of the whole-record refit and every converged profile vector with Δ ≤ 3.84"),
        ("band_low_open", open_lo), ("band_high_open", open_hi), ("open_reason", "; ".join(reasons) or "none"),
        ("C_inv_start_median_mgN_L", float(np.median(c[t0]))),
        ("C_inv_start_band_low_mgN_L", float(np.median(lo[t0]))), ("C_inv_start_band_high_mgN_L", float(np.median(hi[t0]))),
        ("C_inv_end_median_mgN_L", float(np.median(c[t1]))),
        ("C_inv_end_band_low_mgN_L", float(np.median(lo[t1]))), ("C_inv_end_band_high_mgN_L", float(np.median(hi[t1]))),
        ("reading_start_median_mgN_L", float(np.median(np.exp(data.log_r[t0])))),
        ("reading_end_median_mgN_L", float(np.median(np.exp(data.log_r[t1])))),
    ], columns=["key", "value"])
    return Truth(series, vectors, summary)


# --- Stage ---------------------------------------------------------------------------------------
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"


@dataclass
class ModelResult:
    comparison: pd.DataFrame
    identifiability: pd.DataFrame
    checks: pd.DataFrame
    responses: pd.DataFrame
    state: pd.DataFrame
    selected: str
    windows: SplitWindows | None = None
    offsets: list[float] = field(default_factory=list)
    truth: Truth | None = None
    starts: pd.DataFrame | None = None  # every start of the convergence policy (X1)
    volume_candidates: pd.DataFrame | None = None  # Y7: SCENARIO refits per volume scenario (never a primary row)
    volume_summary: pd.DataFrame | None = None


def load_inputs(out_dir: Path) -> dict:
    files = {"grid": "processed/grid_0p01d.parquet", "jumps": "tables/nh4_jumps.csv",
             "pulses": "tables/sensor_validation.csv", "coexistence": "tables/thermo_coexistence.csv"}
    out = {}
    for k, rel in files.items():
        p = Path(out_dir) / rel
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run the earlier stages (`python -m bnr_neq run-all`) first")
        out[k] = pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)
    return out


def compute(cfg: Config, inputs: dict) -> ModelResult:
    data = build_data(cfg, inputs["grid"], inputs["jumps"])
    comparison, fits = split_evaluation(cfg, data)
    starts = comparison.attrs.get("starts")
    selected = comparison.loc[comparison["selected"], "model"].iloc[0]
    spec, lay, cal = fits[selected]
    full = fit(spec, data, lay, (T0_D, END_D), cal.x, np.ones(len(cal.x), dtype=bool), max_nfev=CAL_MAX_NFEV)
    comparison["whole_record_refit_rmse_log"] = np.where(comparison["selected"], math.sqrt(full.rss / full.n), np.nan)
    comparison["whole_record_refit_converged"] = np.where(comparison["selected"], full.converged, pd.NA)
    points: list = []
    ident = profile_likelihood(spec, data, lay, full, collect=points)
    truth = nh4_truth(spec, data, lay, full, ident, points)
    coex = inputs["coexistence"]
    a_scen = float(coex.loc[coex["documented_scenario"].astype(bool), "affinity_kJ_mol"].iloc[0])
    checks = consistency_checks(cfg, spec, data, lay, full, a_scen)
    resp = response_fractions(spec, data, lay, full, inputs["pulses"], cfg.get("sensors.nh4_to_N_factor"))
    lr_full, sim = predict(spec, data, lay, full.x)
    lr_cal, _ = predict(spec, data, lay, cal.x)
    p = unpack(spec, data, lay, full.x)
    state = pd.DataFrame({
        "t_d": data.t, "electrode_obs_mgN_L": np.exp(data.log_r),
        "model_electrode_wholefit_mgN_L": np.exp(lr_full), "model_electrode_calibration_fit_mgN_L": np.exp(lr_cal),
        "model_state_C_mgN_L_conditional_on_C0": sim.c, "model_rate_mgN_L_d": sim.rate,
        "gain_b": p["b0"] + p["b1_per_d"] * (data.t - T0_D)})
    vol_cands, vol_summary = volume_scenario_refits(cfg, data, selected, truth)
    return ModelResult(comparison, ident, checks, resp, state, selected, split_windows(cfg), list(data.offsets), truth,
                       starts, vol_cands, vol_summary)


# --- Y7: the model per SCENARIO volume (never a primary row) ------------------------------------------------------
SCENARIO_ROLE = "SCENARIO refit: never feeds a primary row"


class ScenarioLabelError(ValueError):
    """A volume-scenario refit row does not say SCENARIO."""


def check_scenario_rows(df: pd.DataFrame, col: str = "volume_scenario") -> None:
    """Y7: every row of a volume-scenario refit table names its SCENARIO and carries the SCENARIO role."""
    lab = _as_text(df[col]) if col in df.columns else pd.Series("", index=df.index)
    role = _as_text(df["role"]) if "role" in df.columns else pd.Series("", index=df.index)
    if (~lab.str.contains("SCENARIO") | (role != SCENARIO_ROLE)).any():
        raise ScenarioLabelError(f"volume-scenario refit rows without a SCENARIO label or role: {int((~lab.str.contains('SCENARIO') | (role != SCENARIO_ROLE)).sum())}")


def level_against_band(c: np.ndarray, lo: np.ndarray, hi: np.ndarray, sel: np.ndarray) -> tuple[float, str]:
    """Median of C over ``sel`` and where it lies against the medians of the band edges."""
    m, a, b = float(np.median(c[sel])), float(np.median(lo[sel])), float(np.median(hi[sel]))
    return m, "inside the band" if a <= m <= b else "above the band" if m > b else "below the band"


def volume_scenario_refits(cfg: Config, data: ModelData, selected: str, truth: Truth) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Y7: for every non-constant SCENARIO volume, the whole comparison under the unchanged policy and split, and the
    whole-record refit of the primary selected model; its C(t) level against the constant-volume band."""
    from . import volume

    scens = [v for v in volume.volume_scenarios(cfg, END_D) if not v.is_constant and "SCENARIO" in v.name]
    cand_rows, summary = [], []
    t0, t1 = data.t < T0_D + 0.1, data.t >= END_D - 0.1
    lo, hi = truth.series["C_inv_band_low_mgN_L"].to_numpy(), truth.series["C_inv_band_high_mgN_L"].to_numpy()

    def evaluate(v):
        d = replace(data, volume=v)
        comp, fits = split_evaluation(cfg, d)
        spec, lay, cal = fits[selected]
        return d, comp, fits, fit(spec, d, lay, (T0_D, END_D), cal.x, np.ones(len(cal.x), dtype=bool), max_nfev=CAL_MAX_NFEV)

    # the scenarios are independent and each fit is deterministic: they run concurrently (each with its own pool)
    with ThreadPoolExecutor(max_workers=max(1, len(scens))) as ex:
        done = list(ex.map(evaluate, scens))
    for v, (d, comp, fits, full) in zip(scens, done):
        sel_v = comp.loc[comp["selected"], "model"].iloc[0]
        for r in comp.itertuples(index=False):
            cand_rows.append({"volume_scenario": v.name, "net_loss_mL_per_d": 1e3 * v.net_loss_L_per_d, "model": r.model,
                              "cal_converged": r.cal_converged, "cal_status": r.cal_status, "cal_start": r.cal_start,
                              "cal_rmse_log": r.cal_rmse_log, "sel_skill_vs_null": r.sel_skill_vs_null,
                              "test_skill_vs_null": r.test_skill_vs_null, "forecast_skill_vs_null": r.forecast_skill_vs_null,
                              "selected_in_scenario": r.selected, "selection_status": r.selection_status,
                              "alk0_fitted_meq_L": r.alk0_fitted_meq_L, "role": SCENARIO_ROLE})
        spec, lay, _ = fits[selected]
        c = invert_observation(d.log_r, d.t, unpack(spec, d, lay, full.x), d.offsets)
        (c_start, where_start), (c_end, where_end) = level_against_band(c, lo, hi, t0), level_against_band(c, lo, hi, t1)
        row = comp.set_index("model").loc[selected]
        c_cv = truth.series["C_inv_mgN_L"].to_numpy()
        summary.append({
            "volume_scenario": v.name, "net_loss_mL_per_d": 1e3 * v.net_loss_L_per_d,
            "selected_in_scenario": sel_v, "primary_selected": selected, "selected_model_changes": sel_v != selected,
            "primary_model_sel_skill_vs_null": row["sel_skill_vs_null"], "primary_model_test_skill_vs_null": row["test_skill_vs_null"],
            "primary_model_forecast_skill_vs_null": row["forecast_skill_vs_null"],
            "primary_model_cal_converged": row["cal_converged"],
            "whole_record_refit_rmse_log": math.sqrt(full.rss / full.n), "whole_record_refit_converged": full.converged,
            "C_start_median_mgN_L": c_start, "C_start_vs_constant_volume_band": where_start,
            "C_end_median_mgN_L": c_end, "C_end_vs_constant_volume_band": where_end,
            "constant_volume_C_start_median_mgN_L": float(np.median(c_cv[t0])),
            "constant_volume_C_end_median_mgN_L": float(np.median(c_cv[t1])),
            **phi_scaling(c, c_cv, d.t, v, t0, t1), "role": SCENARIO_ROLE})
    cands = pd.DataFrame(cand_rows, columns=VOLUME_CANDIDATE_COLS)  # fixed columns, also without a scenario
    summ = pd.DataFrame(summary, columns=VOLUME_SUMMARY_COLS)
    for df in (cands, summ):
        check_scenario_rows(df)
    return cands, summ


VOLUME_CANDIDATE_COLS = ("volume_scenario", "net_loss_mL_per_d", "model", "cal_converged", "cal_status", "cal_start",
                         "cal_rmse_log", "sel_skill_vs_null", "test_skill_vs_null", "forecast_skill_vs_null",
                         "selected_in_scenario", "selection_status", "alk0_fitted_meq_L", "role")
VOLUME_SUMMARY_COLS = ("volume_scenario", "net_loss_mL_per_d", "selected_in_scenario", "primary_selected",
                       "selected_model_changes", "primary_model_sel_skill_vs_null", "primary_model_test_skill_vs_null",
                       "primary_model_forecast_skill_vs_null", "primary_model_cal_converged",
                       "whole_record_refit_rmse_log", "whole_record_refit_converged", "C_start_median_mgN_L",
                       "C_start_vs_constant_volume_band", "C_end_median_mgN_L", "C_end_vs_constant_volume_band",
                       "constant_volume_C_start_median_mgN_L", "constant_volume_C_end_median_mgN_L",
                       "C_start_ratio_to_constant_volume", "phi_start_median", "C_end_ratio_to_constant_volume",
                       "phi_end_median", "ratio_within_phi_tolerance", "role")
PHI_SCALING_TOL = 0.05  # Z1: |ratio − φ|/φ, inclusive (docs/progress.md, pass 5)


def phi_deviation(ratio: float, phi: float) -> float:
    """Z1: |ratio − φ|/φ = |ratio/φ − 1| of a scenario-to-constant-volume C ratio against φ = V0/V(t)."""
    return abs(ratio - phi) / phi


def within_phi_tolerance(ratio: float, phi: float) -> bool:
    return phi_deviation(ratio, phi) <= PHI_SCALING_TOL


def phi_scaling(c: np.ndarray, c_cv: np.ndarray, t: np.ndarray, vol, start: np.ndarray, end: np.ndarray) -> dict:
    """Z1: per window, the ratio of the scenario to the constant-volume median C and the median φ = V0/V(t) over the
    same bins; whether every ratio is within ``PHI_SCALING_TOL`` of its φ."""
    out, ok = {}, True
    for w, sel in (("start", start), ("end", end)):
        ratio, phi = float(np.median(c[sel])) / float(np.median(c_cv[sel])), float(np.median(vol.factor(t[sel])))
        out[f"C_{w}_ratio_to_constant_volume"], out[f"phi_{w}_median"] = ratio, phi
        ok = ok and within_phi_tolerance(ratio, phi)
    out["ratio_within_phi_tolerance"] = ok
    return out


def volume_scaling_text(vs: pd.DataFrame) -> str:
    """Z1: whether every scenario C(t) is the constant-volume C(t) scaled by φ = V0/V(t) (start and end windows of
    `model_volume_scenario_summary.csv`, ``PHI_SCALING_TOL``); one text for §1 point 9, §4 and ``model.md``."""
    if vs.empty:
        return ""
    pairs = [(r["net_loss_mL_per_d"], w, float(r[f"C_{w}_ratio_to_constant_volume"]), float(r[f"phi_{w}_median"]))
             for _, r in vs.iterrows() for w in ("start", "end")]
    out = [(e, w, q, p) for e, w, q, p in pairs if not within_phi_tolerance(q, p)]
    tol = f"|ratio/φ − 1| ≤ {PHI_SCALING_TOL:g}"
    if out:
        return ("The scenario C(t) is not the constant-volume C(t) scaled by φ = V0/V(t) within the registered "
                f"tolerance ({tol}): " + "; ".join(f"{e:g} mL/d, {w}: ratio {q:.3f} against φ {p:.3f}"
                                                   for e, w, q, p in out)
                + ". No scaling statement is made, and the scenario and constant-volume N amounts differ there.")
    end = "; ".join(f"{e:g} mL/d {q:.3f} against {p:.3f}" for e, w, q, p in pairs if w == "end")
    above = [f"{r['net_loss_mL_per_d']:g}" for _, r in vs.iterrows()
             if "above the band" in (str(r["C_start_vs_constant_volume_band"]), str(r["C_end_vs_constant_volume_band"]))]
    start_max = max(abs(p - 1) for e, w, q, p in pairs if w == "start")
    return ("Every scenario C(t) is the constant-volume C(t) scaled by φ = V0/V(t): in the start and end windows the "
            "ratio of the scenario to the constant-volume median C equals the median φ of the same window within "
            f"{tol} (end: {end}). The N amount per initial litre, and with it η, is therefore unchanged, and only "
            "the concentrations scale by φ"
            + (f"; \"above the band\" at {' and '.join(above)} mL/d is the evaporative concentration of a band that is "
               "conditional on constant volume, not more N" if above else "")
            + (f". At the start φ − 1 ≤ {start_max:.3f}, below the tolerance, so the start window cannot tell scaled "
               "from unchanged." if start_max < PHI_SCALING_TOL else "."))


def figure(cfg: Config, r: ModelResult, png: Path, svg: Path) -> None:
    from .report import pyplot

    plt = pyplot()

    def style(ax):
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)

    s = r.state
    fig = plt.figure(figsize=(11, 9), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], hspace=0.35, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    style(ax)
    sw = r.windows or split_windows(cfg)
    for (a, b), lab in ((sw.cal, "calibration"), (sw.sel, "selection"), (sw.test, "test")):
        ax.axvspan(a, b, color=GRID if lab != "calibration" else SURFACE, alpha=0.6 if lab == "selection" else 0.9, lw=0)
        ax.text(0.5 * (a + b), 0.97, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8,
                color=INK2)
    ax.plot(s["t_d"], s["electrode_obs_mgN_L"], color=BLUE, lw=0.8, label="Electrode (observed)")
    ax.plot(s["t_d"], s["model_electrode_calibration_fit_mgN_L"], color=ORANGE, lw=1.4,
            label=f"Calibrated on days {sw.cal[0]:g}-{sw.cal[1]:g}, frozen afterwards (test offsets not refitted here)")
    ax.plot(s["t_d"], s["model_electrode_wholefit_mgN_L"], color=AQUA, lw=1.2, ls="--", label="Whole-record refit")
    for tp in cfg.pulse_times_d:
        ax.axvline(tp, color=INK2, lw=0.8, ls=":")
    ax.set_ylabel("NH4+ electrode (mg N/L)", fontsize=8.5, color=INK)
    ax.set_xlabel("Time (d)", fontsize=8.5, color=INK2)
    ax.set_title(f"a  Selected model ({r.selected}) with the electrode observation model", loc="left", fontsize=10,
                 color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[1, 0])
    style(ax)
    ax.plot(s["t_d"], s["model_rate_mgN_L_d"], color=VIOLET, lw=1.5)
    for tp in cfg.pulse_times_d:
        ax.axvline(tp, color=INK2, lw=0.8, ls=":")
    ax.set_ylabel("Model consumption rate (mg N/L/d)", fontsize=8.5, color=INK)
    ax.set_xlabel("Time (d)", fontsize=8.5, color=INK2)
    ax.set_title("b  Consumption rate of the whole-record refit", loc="left", fontsize=10, color=INK)

    ax = fig.add_subplot(gs[1, 1])
    style(ax)
    idf = r.identifiability
    y = np.arange(len(idf))
    lo = idf["profile_max_delta_low"].clip(upper=40)
    hi = idf["profile_max_delta_high"].clip(upper=40)
    ax.scatter(lo, y - 0.15, marker="<", color=BLUE, s=28, label="Low side")
    ax.scatter(hi, y + 0.15, marker=">", color=ORANGE, s=28, label="High side")
    ax.axvline(CHI2_95_1DOF, color=INK2, lw=1, ls="--")
    ax.set_yticks(y, [n + ("" if ok else "  (not identifiable)") for n, ok in zip(idf["parameter"], idf["identifiable"])],
                  fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Max profile Δ reached (n_eff-scaled; capped at 40)", fontsize=8.5, color=INK2)
    ax.set_title("c  Profile likelihood: identifiability (3.84 line)", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.savefig(png, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    fig.savefig(svg, facecolor=SURFACE, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def write_report(cfg: Config, r: ModelResult, path: Path) -> None:
    from .report import md_table

    c = r.comparison
    sel = c[c["selected"]].iloc[0]
    idf = r.identifiability
    sw = r.windows or split_windows(cfg)
    CAL, SEL, TEST = sw.cal, sw.sel, sw.test
    n_test_off = sum(t >= TEST[0] for t in r.offsets)
    lines = [
        "# Reduced model with the electrode observation model (M6)",
        "",
        "Generated by `python -m bnr_neq model`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Model",
        "",
        (f"State: true NH4 C(t) in mg N/L from day {T0_D}, +{thermo.pulse_dose_from_config(cfg).n_mg_L:.4f} mg N/L "
        "exactly at each pulse, dC/dt = −r ≥ 0. Rate laws: Haldane in free ammonia, r = vmax S/(Ks + S + S²/Ki) with "
        "S = FA(C, pH(t), temperature SCENARIO), fitted as α = vmax/Ks, Ks and S_opt = √(Ks Ki) (r = α S/(1 + S/Ks + "
        "S²/S_opt²), the same law; pass 4, Y6); first-order decay of biomass activity, r = v0 exp(−kd (t − 1)) "
        f"C/(K_NH + C) with K_NH = {K_NH_ASM_MG_N_L} mg N/L fixed (ASM1 default, not fitted); combined, r = v0 "
        "exp(−kd (t − 1)) x Haldane(S); W5 substrate-side candidates: alkalinity-limited, r = v0 A/(K_ALK + A) "
        f"C/(K_NH + C), and alkalinity x Haldane(S), with the alkalinity state A = Alk0 − 2 n/M_N (n = N consumed "
        f"since day 1, 2 eq/mol N, nothing added), Alk0 fitted and K_ALK = {K_ALK_ASM3_MEQ_L} meq/L fixed (ASM3 "
        "default); and a no-consumption reference. Observation: log R = a + (b0 + b1 (t − 1)) "
        f"log C + offsets at every electrode-offset jump detected in M2 ({len(r.offsets)}; the M2 detections are used), "
        "then a first-order lag τ (exact exponential filter). Residuals in log R on the 0.01 d grid."),
        "",
        (f"What the (α, Ks, S_opt) form did and did not do (pass 5, Z3; {HALDANE_CONDITIONING_NOTE}): it removed the "
         "near-perfect correlations of (vmax, Ks, Ki) (largest kinetic |corr| 0.99998 → 0.89 at the alkalinity x FA "
         "optimum), but not the ill-conditioning. The natural-unit Haldane parameters of the alkalinity x FA "
         "candidates (vmax about 2,600 mg N/L/d, Ki about 0.005 mg NH3/L at that calibration) are therefore not "
         "estimates; only their skills and Alk0 are reported."),
        "",
        "## Chronological split and out-of-sample skill",
        "",
        (f"Calibrate on days {CAL[0]}-{CAL[1]} (offsets after day {CAL[1]} not fitted), select on {SEL[0]}-{SEL[1]} with all "
        f"parameters frozen, test on {TEST[0]}-{TEST[1]} (split from config [pulses] times_d). In the test window, the "
        f"{n_test_off} offset(s) detected there (unseen in "
        "calibration) are refitted as the only nuisance parameters; the pure forecast (nothing refitted) is reported up to "
        "the first unseen offset. Skill = 1 − MSE/MSE(no-consumption reference), in log units."),
        "",
        (f"Convergence policy, identical for every candidate (pass 3, X1): {CONVERGENCE_POLICY} Only a candidate whose "
         "calibration converged can be selected (`selection_status`); `alk0_*` is the fitted initial alkalinity of the "
         "alkalinity candidates with its 95 % profile interval on the calibration window (n_eff-scaled, as below)."),
        "",
        md_table(c.drop(columns=["test_nuisance"])),
        "",
        bounds_text(c),
        "",
        (f"Selected on the selection window: **{r.selected}** (selection skill {sel['sel_skill_vs_null']:.3f}, test skill "
        f"{sel['test_skill_vs_null']:.3f}, pure-forecast skill {sel['forecast_skill_vs_null']:.3f}). It was then refitted "
        f"on the whole record (RMSE {sel['whole_record_refit_rmse_log']:.4f} log units)."),
        "",
        "## Identifiability of the whole-record refit (profile likelihood)",
        "",
        ("Δ(θ) = n_eff ln(RSS(θ)/RSS_min), with n_eff = n (1 − ρ)/(1 + ρ) from the residual lag-1 autocorrelation ρ, so "
        f"that autocorrelated residuals do not overstate the information. {PROFILE_POINTS_PER_SIDE} points per side, "
        "warm-started; a side is closed only by a converged re-optimization with Δ ≥ 3.84 "
        f"({idf.attrs.get('uncapped_refits', 0)} points were refitted without an evaluation cap). The interval ends are "
        "the roots of Δ(θ) = 3.84 found with Brent's method on the profile itself inside the bracketing grid step "
        f"(tolerance {PROFILE_XTOL_FRACTION:g} x the grid half-width; {idf.attrs.get('profile_evaluations', 0)} profile "
        "evaluations), not linear interpolations. **Parameters marked not identifiable are not "
        "estimates**; their optimizer values are shown only for traceability."),
        "",
        md_table(idf),
        "",
        ("Interpretation, conditional on the model: the absolute NH4 level (C0) and the electrode intercept a are not "
        "identifiable, so the model does not determine absolute NH4 (NO-GO for absolute concentrations stands). The "
        "dose anchors the response, and the declining response fractions are explained by accumulation of NH4 under a "
        "log-type electrode with a non-drifting gain b, together with decaying activity, not by a drifting gain. The "
        "lag τ cannot reproduce the 7.5 h delayed jumps (a pure delay), and the residuals stay strongly autocorrelated: "
        "the model is misspecified at short time scales."),
        "",
        "## Consistency checks",
        "",
        md_table(r.checks),
        "",
        "## Model-implied electrode response to each pulse (M2 window) vs observed",
        "",
        md_table(r.responses),
        "",
        "## Single NH4 observation model: inverted C(t) (W1)",
        "",
        (f"Every downstream NH4 quantity is derived from C(t) obtained by inverting this observation model with the "
        f"whole-record refit ({INVERSION_NOTE}). Its uncertainty band is the envelope of the inversions under every "
        "converged profile vector with Δ ≤ 3.84 (`nh4_truth_vectors.csv`). Because C0 and the intercept a are not "
        "identifiable, the band is open on the side(s) listed below: the data do not bound the absolute level there, "
        "and the band edge is only the last evaluated vector. `outputs/processed/nh4_truth.parquet` holds the series."),
        "",
        md_table(r.truth.summary) if r.truth is not None else "not computed",
        "",
        "## The model per SCENARIO volume (Y7; never a primary row)",
        "",
        ("Every candidate calibrated and selected under the same policy and split with V(t) = V0 − E t (pulse doses, C, "
         "the biomass activity and the alkalinity state concentrated by V0/V(t)); the primary selected model refitted "
         "on the whole record, and its C(t) compared with the constant-volume band. SCENARIO values only: the primary "
         "C(t) and every primary row stay at constant volume. Every candidate row: `model_volume_scenarios.csv`."
         + (" " + volume_scaling_text(r.volume_summary) if r.volume_summary is not None and len(r.volume_summary)
            else "")),
        "",
        md_table(r.volume_summary.drop(columns=["role"])) if r.volume_summary is not None else "not computed",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir) -> ModelResult:
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    inputs = load_inputs(out_dir)
    d = out_dirs(out_dir)
    r = compute(cfg, inputs)
    r.comparison.to_csv(d["tables"] / "model_fit.csv", index=False)
    r.starts.to_csv(d["tables"] / "model_calibration_starts.csv", index=False)
    r.identifiability.to_csv(d["tables"] / "model_identifiability.csv", index=False)
    r.checks.to_csv(d["tables"] / "model_checks.csv", index=False)
    r.responses.to_csv(d["tables"] / "model_response_fractions.csv", index=False)
    r.state.to_parquet(d["processed"] / "model_state.parquet", index=False)
    r.truth.series.to_parquet(d["processed"] / "nh4_truth.parquet", index=False)
    r.truth.vectors.to_csv(d["tables"] / "nh4_truth_vectors.csv", index=False)
    r.truth.summary.to_csv(d["tables"] / "nh4_truth_summary.csv", index=False)
    r.volume_candidates.to_csv(d["tables"] / "model_volume_scenarios.csv", index=False)
    r.volume_summary.to_csv(d["tables"] / "model_volume_scenario_summary.csv", index=False)
    figure(cfg, r, d["figures"] / "fig6_model.png", d["figures"] / "fig6_model.svg")
    write_report(cfg, r, d["reports"] / "model.md")
    return r
