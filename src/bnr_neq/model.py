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
from dataclasses import dataclass, field
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

# --- W1: the single NH4 observation model ---------------------------------------------------
# Every downstream NH4 quantity (rates, free ammonia, conversion fraction, oxygen demand, affinity step) is derived
# from C(t) obtained by inverting this observation model; the linear gain correction survives only as a labelled
# SENSITIVITY case and the raw electrode reading only as a measured slope.
OBSERVATION_MODEL = "M6_log"
RAW_READING = "electrode_reading"  # measured slope of the reading; never an NH4 truth, never feeds a derived quantity
SENSITIVITY_SUFFIX = "_SENSITIVITY"
LINEAR_GAIN = "linear_gain" + SENSITIVITY_SUFFIX
IDENTITY = "identity_reading_as_TAN" + SENSITIVITY_SUFFIX
DERIVED_PREFIXES = ("fa_", "o2_", "kla_", "eta_", "delta_A", "free_energy", "alk0_required")
INVERSION_NOTE = ("quasi-static inversion C = exp((log R − a − offsets)/b(t)) of the 0.01 d grid reading; the sensor "
                  "lag τ is not inverted (slopes after ~3τ are unaffected, window starts within 3τ of a boundary carry "
                  "the lag transient)")


class ObservationModelError(ValueError):
    """An output combines NH4 quantities that rest on different observation models."""


def check_single_observation_model(df: pd.DataFrame, model_col: str = "nh4_observation_model",
                                   label_cols: tuple[str, ...] = ("metric", "basis"),
                                   series_col: str = "nh4_series") -> str | None:
    """Raise :class:`ObservationModelError` unless every primary NH4 row of ``df`` uses one observation model.

    Rows whose model ends in ``_SENSITIVITY`` must say SENSITIVITY in their label columns (and only they may);
    the raw electrode reading may appear only as a measured, non-derived quantity. When ``series_col`` is present,
    the primary rows must also share one series (:func:`check_single_series`). Returns the primary model.
    """
    m = df[model_col].fillna("").astype(str)
    cols = [c for c in label_cols if c in df.columns]
    # empty cells (NaN when read back from CSV) count as empty text
    text = (df[cols].apply(lambda col: col.map(lambda v: "" if pd.isna(v) else str(v))).agg(" ".join, axis=1)
            if cols else pd.Series("", index=df.index))
    sens = m.str.endswith(SENSITIVITY_SUFFIX)
    if (sens & ~text.str.contains("SENSITIVITY")).any():
        raise ObservationModelError(f"sensitivity rows without a SENSITIVITY label: {sorted(set(m[sens]))}")
    if (~sens & (m != "") & text.str.contains("SENSITIVITY")).any():
        raise ObservationModelError("a row labelled SENSITIVITY carries a primary observation model")
    if "metric" in df.columns:
        raw_derived = (m == RAW_READING) & df["metric"].astype(str).str.startswith(DERIVED_PREFIXES)
        if raw_derived.any():
            raise ObservationModelError(f"derived quantities from the raw electrode reading: "
                                        f"{sorted(set(df.loc[raw_derived, 'metric']))}")
    is_primary = (m != "") & (m != RAW_READING) & ~sens
    primary = sorted(set(m[is_primary]))
    if len(primary) > 1:
        raise ObservationModelError(f"one output mixes NH4 observation models: {primary}")
    if series_col in df.columns and is_primary.any():
        check_single_series(df[is_primary], text[is_primary], series_col)
    return primary[0] if primary else None


MODEL_STATE_PREFIX = "M6_state:"  # series of the simulated model state (MODEL-BASED quantities)


def check_single_series(rows: pd.DataFrame, text: pd.Series, series_col: str = "nh4_series") -> None:
    """Primary rows must share one inverted C(t) series; MODEL-BASED rows (and only they) may instead use one simulated
    model-state series of the same observation model. Raises :class:`ObservationModelError` otherwise."""
    s = rows[series_col].map(lambda v: "" if pd.isna(v) else str(v))
    model_based = text.str.contains("MODEL-BASED|MODEL_BASED", regex=True)
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


def model_specs(cfg: Config) -> list[ModelSpec]:
    specs = [ModelSpec("null_no_consumption", (), False, False),
             ModelSpec("biomass_decay", ("v0_mgN_L_d", "kd_per_d"), False, True)]
    for T in cfg.get("scenarios.temperature_C"):
        specs.append(ModelSpec(f"haldane_FA_SCENARIO_T{T:g}C", ("vmax_mgN_L_d", "Ks_mgNH3_L", "Ki_mgNH3_L"), True, False, T))
        specs.append(ModelSpec(f"combined_SCENARIO_T{T:g}C", ("v0_mgN_L_d", "kd_per_d", "Ks_mgNH3_L", "Ki_mgNH3_L"),
                               True, True, T))
    # W5: substrate-side limitation by alkalinity / inorganic carbon, alone and with FA inhibition (TAN = model state,
    # i.e. the same W1 observation model), competing with biomass decay under the unchanged chronological split
    specs.append(ModelSpec("alkalinity_limited", ("v0_mgN_L_d", "Alk0_meq_L"), False, False, alk=True))
    for T in cfg.get("scenarios.temperature_C"):
        specs.append(ModelSpec(f"alkalinity_FA_SCENARIO_T{T:g}C", ("vmax_mgN_L_d", "Ks_mgNH3_L", "Ki_mgNH3_L", "Alk0_meq_L"),
                               True, False, T, alk=True))
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
    """r(t, C, f, n) >= 0 in mg N/L/d; f = free-ammonia factor at t; n = N consumed since T0 (mg N/L), which sets the
    alkalinity state of the W5 candidates (ignored by the others)."""
    if spec.name == "null_no_consumption":
        return lambda t, c, f, n=0.0: 0.0
    v = kin.get("v0_mgN_L_d", kin.get("vmax_mgN_L_d"))
    kd = kin.get("kd_per_d", 0.0)
    ks, ki = kin.get("Ks_mgNH3_L"), kin.get("Ki_mgNH3_L")
    alk0 = kin.get("Alk0_meq_L")

    def rate(t, c, f, n=0.0):
        c = max(0.0, c)
        act = v * math.exp(-kd * (t - T0_D)) if spec.decay else v
        if spec.alk:
            a = alkalinity_state_meq_L(alk0, n)
            act *= a / (K_ALK_ASM3_MEQ_L + a)
        if spec.fa:
            s = f * c
            return act * s / (ks + s + s * s / ki)
        return act * c / (K_NH_ASM_MG_N_L + c)

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
    sim = simulate(spec, p["kin"], p["C0_mgN_L"], data.t, data.ph, data.pulses, data.dose)
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


def fit(spec: ModelSpec, data: ModelData, lay: Layout, window: tuple[float, float], x0: np.ndarray,
        free: np.ndarray, max_nfev: int | None = None) -> Fit:
    """Least squares in log space over bins in ``window``; parameters not in ``free`` stay at x0."""
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
    return Fit(spec, x, free, window, float(sol.fun @ sol.fun), int(mask.sum()), bool(sol.success))


def window_skill(spec, data, lay, x, window) -> dict:
    mask = (data.t >= window[0]) & (data.t < window[1])
    lr, _ = predict(spec, data, lay, x)
    e = (lr - data.log_r)[mask]
    lin = (np.exp(lr) - np.exp(data.log_r))[mask]
    return {"rmse_log": float(np.sqrt(np.mean(e**2))), "rmse_mgN_L": float(np.sqrt(np.mean(lin**2))),
            "n": int(mask.sum())}


def split_evaluation(cfg: Config, data: ModelData, windows: SplitWindows | None = None) -> tuple[pd.DataFrame, dict]:
    """Calibrate on CAL, evaluate frozen on SEL, then TEST (unseen offsets as the only nuisance)."""
    w = split_windows(cfg) if windows is None else windows
    CAL, SEL, TEST = w.cal, w.sel, w.test
    rows, fits = [], {}
    for spec in model_specs(cfg):
        lay = layout(spec, data)
        x0 = initial_guess(spec, data, lay)
        free = np.ones(len(x0), dtype=bool)
        for i, t in enumerate(data.offsets):
            if t >= CAL[1]:
                free[len(lay.names) - len(data.offsets) + i] = False
                x0[len(lay.names) - len(data.offsets) + i] = 0.0
        cal = fit(spec, data, lay, CAL, x0, free)
        sel = window_skill(spec, data, lay, cal.x, SEL)
        test_off = np.array([n.startswith("offset_log_") and float(n[11:-1]) >= TEST[0] for n in lay.names])
        test_fit = fit(spec, data, lay, TEST, cal.x, test_off)
        test = window_skill(spec, data, lay, test_fit.x, TEST)
        unseen = [t for t in data.offsets if t >= TEST[0]]
        fc_window = (TEST[0], min(unseen) if unseen else TEST[1])
        forecast = window_skill(spec, data, lay, cal.x, fc_window)
        rows.append({"model": spec.name, "n_params_calibration": int(free.sum()),
                     "cal_rmse_log": math.sqrt(cal.rss / cal.n), "cal_converged": cal.converged,
                     "sel_rmse_log": sel["rmse_log"], "sel_rmse_mgN_L": sel["rmse_mgN_L"],
                     "test_rmse_log": test["rmse_log"], "test_rmse_mgN_L": test["rmse_mgN_L"],
                     "test_nuisance": f"offsets at {', '.join(f'{t:.3f}' for t in unseen)} d refitted on the test window; "
                                      "all other parameters frozen",
                     "forecast_window_d": f"{fc_window[0]:.2f} to {fc_window[1]:.3f}",
                     "forecast_rmse_log": forecast["rmse_log"], "forecast_rmse_mgN_L": forecast["rmse_mgN_L"]})
        fits[spec.name] = (spec, lay, cal)
    df = pd.DataFrame(rows)
    null = df.set_index("model").loc["null_no_consumption"]
    for col in ("sel", "test", "forecast"):
        df[f"{col}_skill_vs_null"] = 1 - df[f"{col}_rmse_log"] ** 2 / null[f"{col}_rmse_log"] ** 2
    cand = df[df["model"] != "null_no_consumption"]
    selected = cand.loc[cand["sel_rmse_log"].idxmin(), "model"]
    df["selected"] = df["model"] == selected
    return df, fits


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


def profile_likelihood(spec, data, lay, best: Fit, collect: list | None = None) -> pd.DataFrame:
    """Profile Δ(θ_i) = n_eff ln(RSS(θ_i)/RSS_min); identifiable if Δ ≥ 3.84 is reached on both sides within
    the bounds; the 95 % interval ends are the roots of Δ(θ_i) = 3.84 (:func:`profile_interval`).

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
        half = max(4 * se[i], 1.0 if lay.log[i] else 0.2 if name != "b1_per_d" else 0.01)
        free = np.ones(len(best.x), dtype=bool)
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
            "better_optimum_found": bool(min_delta < -0.5),
            "n_eff": n_eff, "resid_lag1": rho,
        })
    out = pd.DataFrame(rows)
    out.attrs["uncapped_refits"] = n_refits[0]
    out.attrs["profile_evaluations"] = n_evals[0]
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
    selected = comparison.loc[comparison["selected"], "model"].iloc[0]
    spec, lay, cal = fits[selected]
    full = fit(spec, data, lay, (T0_D, END_D), cal.x, np.ones(len(cal.x), dtype=bool))
    comparison["whole_record_refit_rmse_log"] = np.where(comparison["selected"], math.sqrt(full.rss / full.n), np.nan)
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
    return ModelResult(comparison, ident, checks, resp, state, selected, split_windows(cfg), list(data.offsets), truth)


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
        "S = FA(C, pH(t), temperature SCENARIO); first-order decay of biomass activity, r = v0 exp(−kd (t − 1)) "
        f"C/(K_NH + C) with K_NH = {K_NH_ASM_MG_N_L} mg N/L fixed (ASM1 default, not fitted); combined, r = v0 "
        "exp(−kd (t − 1)) x Haldane(S); W5 substrate-side candidates: alkalinity-limited, r = v0 A/(K_ALK + A) "
        f"C/(K_NH + C), and alkalinity x Haldane(S), with the alkalinity state A = Alk0 − 2 n/M_N (n = N consumed "
        f"since day 1, 2 eq/mol N, nothing added), Alk0 fitted and K_ALK = {K_ALK_ASM3_MEQ_L} meq/L fixed (ASM3 "
        "default); and a no-consumption reference. Observation: log R = a + (b0 + b1 (t − 1)) "
        f"log C + offsets at every electrode-offset jump detected in M2 ({len(r.offsets)}; the M2 detections are used), "
        "then a first-order lag τ (exact exponential filter). Residuals in log R on the 0.01 d grid."),
        "",
        "## Chronological split and out-of-sample skill",
        "",
        (f"Calibrate on days {CAL[0]}-{CAL[1]} (offsets after day {CAL[1]} not fitted), select on {SEL[0]}-{SEL[1]} with all "
        f"parameters frozen, test on {TEST[0]}-{TEST[1]} (split from config [pulses] times_d). In the test window, the "
        f"{n_test_off} offset(s) detected there (unseen in "
        "calibration) are refitted as the only nuisance parameters; the pure forecast (nothing refitted) is reported up to "
        "the first unseen offset. Skill = 1 − MSE/MSE(no-consumption reference), in log units."),
        "",
        md_table(c.drop(columns=["test_nuisance"])),
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
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir) -> ModelResult:
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    inputs = load_inputs(out_dir)
    d = out_dirs(out_dir)
    r = compute(cfg, inputs)
    r.comparison.to_csv(d["tables"] / "model_fit.csv", index=False)
    r.identifiability.to_csv(d["tables"] / "model_identifiability.csv", index=False)
    r.checks.to_csv(d["tables"] / "model_checks.csv", index=False)
    r.responses.to_csv(d["tables"] / "model_response_fractions.csv", index=False)
    r.state.to_parquet(d["processed"] / "model_state.parquet", index=False)
    r.truth.series.to_parquet(d["processed"] / "nh4_truth.parquet", index=False)
    r.truth.vectors.to_csv(d["tables"] / "nh4_truth_vectors.csv", index=False)
    r.truth.summary.to_csv(d["tables"] / "nh4_truth_summary.csv", index=False)
    figure(cfg, r, d["figures"] / "fig6_model.png", d["figures"] / "fig6_model.svg")
    write_report(cfg, r, d["reports"] / "model.md")
    return r
