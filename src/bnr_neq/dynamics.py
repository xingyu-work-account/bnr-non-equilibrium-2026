"""Phase 3 applied: estimator validation, then irreversibility, complexity, Onsager regression,
stability indicators and probability flux on the plant data (M5).

The estimator validation suite runs first; if any required check fails, the stage
writes the validation table and stops without computing anything on plant data.

pH windows are additionally gated by a matched null: reversible quantized AR(1)
regimes spanning the tie fraction, lag-1 autocorrelation and DFA exponent of the
analysed pH windows are tested with the same IAAFT procedure; an index is reported
for a pH window only if every regime it falls in keeps the false-positive rate at
or below 0.05 plus a binomial margin.
"""

from __future__ import annotations

import itertools
import math
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.signal import lfilter
from scipy.stats import binom, binomtest, chi2, kendalltau, norm, wilcoxon
from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.stattools import adfuller, kpss

from . import estimators as E
from . import timebase
from .config import Config, MissingInputError

BASE_SEED = 20261018
LEVEL = 0.95

# --- Validation suite parameters ------------------------------------------------------
MARKOV_PAIRS = ((0.3, 0.3), (0.4, 0.2), (0.6, 0.1), (0.45, 0.05))
MARKOV_STEPS = 200_000
MARKOV_BLOCK = 100
MARKOV_BOOT = 300
MARKOV_LEVEL = 0.99  # declared before running: fewer chance misses when the suite is rerun with other seeds
QUANT_LEVEL_SPACING = 0.02  # state s observed at 0.02 s + N(0, 0.007), quantized to 0.01
QUANT_NOISE_SD = 0.007
QUANT_STEP = 0.01
SURR_REALIZATIONS = 20
SURR_LENGTH = 2160  # one pH window
N_SURR = 99  # IAAFT surrogates per test, in the validation suite and on the data (rank p resolution 0.01)
AR_PHI = 0.95
# Tie sweep: AR(1) scaled to this sd, then quantized to 0.01 (None = continuous). Tie fractions span 0 to
# beyond the pH record (0.82 to 0.96). The first N_REQUIRED_TIE_LEVELS levels are required checks; the
# others define each index's validated tie range, which the data analysis enforces.
TIE_SWEEP_SD = (None, 0.05, 0.02, 0.01, 0.006, 0.004, 0.003)
N_REQUIRED_TIE_LEVELS = 2
MAX_FALSE_POSITIVES = 4  # of 20 at alpha 0.05: P(X >= 5 | 0.05, 20) = 0.016
MIN_POWER = 16  # of 20 for the positive control
OU_STEPS = 20_000
OU_DT = 0.01
OU_OMEGA = 1.0
# Matched null for the pH windows: reversible AR(1) at these phi, scaled so that the series quantized to 0.01
# has the target tie fraction. The grid spans the analysed pH windows (tie fraction 0.6 to 0.96, lag-1
# autocorrelation 0.82 to 0.998, DFA alpha 0.9 to 1.55). Only realizations that pass the same ADF/KPSS filter as
# the data windows enter a regime (the data indices are computed only on such windows).
NULL_PHI = (0.9, 0.95, 0.99, 0.995, 0.999)
NULL_TIE = (0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
NULL_TARGET_STATIONARY = 100  # stationary realizations per regime
NULL_MAX_REALIZATIONS = 4000  # realizations drawn per regime at most (high phi rarely passes ADF/KPSS)
NULL_BATCH = 200
NULL_MIN_STATIONARY = 30  # fewer stationary realizations: the regime is not validated
FP_ALPHA = 0.05
FP_MARGIN_LEVEL = 0.95  # threshold: the 95 % binomial quantile of the false-positive count at rate 0.05
N_WORKERS = int(os.environ.get("BNR_NEQ_WORKERS", os.cpu_count() or 1))  # results do not depend on this

# --- Data analysis parameters ------------------------------------------------------------
START_D = 1.0  # start-up excluded (plan item D10)
END_D = 47.92
WINDOW_D = 0.25
MIN_WINDOW_FILL = 0.9
ORDINAL_D = 3
ASYM_LAGS = (1, 3, 10, 30, 100)  # samples
SAMPEN_M = 2
SAMPEN_R_SD = 0.2
MSE_SCALES = (1, 2, 4)
ADF_ALPHA = 0.05
KPSS_ALPHA = 0.05
NH4_BIN_S = 20.0  # NH4 resampled to 20 s time-weighted bins (logging rate changes from 20 s to 2 s)
ORP_ELECTRODE_EVENTS_D = (28.89, 30.48)  # implementation_prompt.md M5 item 3
EVENT_MARGIN_D = (0.05, 0.1)  # excluded before / after each event
ONSAGER_TAIL_START_D = 0.25
ONSAGER_LATE_FRACTION = 0.4
ONSAGER_MAX_LAG_D = 1.0
ONSAGER_MIN_C = 0.2
ONSAGER_MIN_SNR = 3.0
ONSAGER_BLOCK = 10
ONSAGER_BOOT = 500
STABILITY_SPAN_D = 2.0
STABILITY_WINDOW = 50  # bins (0.5 d)
FLUX_BLOCK = 25
FLUX_BOOT = 999
FLUX_MIN_INCREMENTS = 100

INDEX_FNS = {
    "asym_norm_lag1": lambda z: E.moment_asymmetry(z, [1])[1][0],
    "ordinal_kl": lambda z: E.ordinal_irreversibility(z, ORDINAL_D),
    "hvg_kl": E.hvg_irreversibility,
}
TWO_SIDED = frozenset({"asym_norm_lag1"})


def _rng(*keys: int) -> np.random.Generator:
    return np.random.default_rng([BASE_SEED, *keys])


# =============================================================================
# Estimator validation
# =============================================================================
def _ar1(n: int, phi: float, rng: np.random.Generator) -> np.ndarray:
    """Stationary AR(1) with unit innovations: x_0 ~ N(0, 1/(1 − phi²)), x_i = phi x_{i−1} + e_i."""
    e = rng.normal(size=n)
    x0 = e[0] / np.sqrt(1 - phi**2)
    if n == 1:
        return np.array([x0])
    rest, _ = lfilter([1.0], [1.0, -phi], e[1:], zi=[phi * x0])
    return np.r_[x0, rest]


def _tie_probability(sd_increment: float, q: float = QUANT_STEP) -> float:
    """P(two consecutive values quantize equally) for a Gaussian increment of this sd and a uniform grid phase:
    E[(1 − |Δ|/q)+] = (2Φ(a) − 1) − (2/(a√(2π)))(1 − exp(−a²/2)), a = q/sd."""
    a = q / sd_increment
    return float((2 * norm.cdf(a) - 1) - (2 / (a * math.sqrt(2 * math.pi))) * (1 - math.exp(-a * a / 2)))


def quantized_ar1(n: int, phi: float, tie_target: float, rng: np.random.Generator) -> np.ndarray:
    """Reversible AR(1) scaled so that its 0.01-quantized version has about ``tie_target`` ties."""
    sd_inc = brentq(lambda s: _tie_probability(s) - tie_target, 1e-6, 10.0)
    z = _ar1(n, phi, rng) * sd_inc / math.sqrt(2 / (1 + phi))  # increment sd of unit-innovation AR(1): √(2/(1+φ))
    return np.round(z / QUANT_STEP) * QUANT_STEP


def dfa_scales(n: int) -> np.ndarray:
    return np.unique(np.logspace(1, np.log10(n // 4), 10).astype(int))


def window_stats(z: np.ndarray) -> dict[str, float]:
    """Tie fraction, lag-1 autocorrelation and DFA exponent: the coordinates used to match a null regime."""
    z = np.asarray(z, dtype=float)
    if np.std(z) == 0:
        return {"tie_fraction": E.tie_fraction(z), "ac1": np.nan, "dfa_alpha": np.nan}
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = E.dfa(z, dfa_scales(len(z)))[0]
    return {"tie_fraction": E.tie_fraction(z), "ac1": float(np.corrcoef(z[:-1], z[1:])[0, 1]),
            "dfa_alpha": float(alpha) if np.isfinite(alpha) else np.nan}


def _ou_2d(n: int, dt: float, omega: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """dX = (−X + ω J X) dt + dW, J the 90° rotation; detailed balance holds only for ω = 0."""
    x = np.zeros((n, 2))
    a = np.array([[-1.0, -omega], [omega, -1.0]])
    noise = rng.normal(size=(n, 2)) * np.sqrt(dt)
    for i in range(1, n):
        x[i] = x[i - 1] + a @ x[i - 1] * dt + noise[i]
    return x[:, 0], x[:, 1]


def null_regimes() -> list[tuple[str, float, float]]:
    return [(f"R{i:02d}", phi, tie) for i, (phi, tie) in enumerate(itertools.product(NULL_PHI, NULL_TIE), start=1)]


def _null_realization(rid: int, phi: float, tie: float, m: int) -> np.ndarray:
    return quantized_ar1(SURR_LENGTH, phi, tie, _rng(40, rid, m))


def _null_stationarity_task(task: tuple[int, float, float, int]) -> tuple[int, int, bool, dict]:
    rid, phi, tie, m = task
    z = _null_realization(rid, phi, tie, m)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            ok = passes_stationarity(z)
        except (np.linalg.LinAlgError, ValueError):
            ok = False
        return rid, m, ok, window_stats(z) if ok else {}


def _null_surrogate_task(task: tuple[int, float, float, int]) -> tuple[int, int, dict[str, float]]:
    rid, phi, tie, m = task
    z = _null_realization(rid, phi, tie, m)
    out = E.surrogate_indices(z, INDEX_FNS, _rng(41, rid, m), N_SURR, TWO_SIDED)
    return rid, m, {name: pv for name, (_, _, pv) in out.items()}


def _map(fn, tasks: list, pool: ProcessPoolExecutor | None) -> list:
    """Deterministic map: every task carries its own seed, results come back in task order."""
    if pool is None or len(tasks) < 2:
        return [fn(t) for t in tasks]
    return list(pool.map(fn, tasks, chunksize=max(1, len(tasks) // (8 * N_WORKERS))))


def _pool() -> ProcessPoolExecutor | None:
    return ProcessPoolExecutor(max_workers=N_WORKERS) if N_WORKERS > 1 else None


def fp_threshold(n: int) -> int:
    """Largest false-positive count consistent with rate 0.05: the 95 % binomial quantile (0.05 + margin)."""
    return int(binom.ppf(FP_MARGIN_LEVEL, n, FP_ALPHA))


def matched_null_table(pool: ProcessPoolExecutor | None = None) -> pd.DataFrame:
    """False-positive rate of each index under each reversible quantized AR(1) regime (stationary realizations)."""
    regimes = null_regimes()
    stats = {rid: [] for rid, _, _ in regimes}
    drawn = dict.fromkeys(stats, 0)
    open_ = {i: r for i, r in enumerate(regimes, start=1)}
    while open_:
        tasks = [(i, phi, tie, m) for i, (rid, phi, tie) in open_.items()
                 for m in range(drawn[rid], min(drawn[rid] + NULL_BATCH, NULL_MAX_REALIZATIONS))]
        for i, m, ok, st in _map(_null_stationarity_task, tasks, pool):
            if ok:
                stats[regimes[i - 1][0]].append((m, st))
        for i, (rid, _, _) in list(open_.items()):
            drawn[rid] = min(drawn[rid] + NULL_BATCH, NULL_MAX_REALIZATIONS)
            if len(stats[rid]) >= NULL_TARGET_STATIONARY or drawn[rid] >= NULL_MAX_REALIZATIONS:
                del open_[i]
    chosen = {rid: sorted(v, key=lambda ms: ms[0])[:NULL_TARGET_STATIONARY] for rid, v in stats.items()}
    tasks = [(i, phi, tie, m) for i, (rid, phi, tie) in enumerate(regimes, start=1) for m, _ in chosen[rid]]
    pvals = {}
    for i, m, pv in _map(_null_surrogate_task, tasks, pool):
        pvals[(i, m)] = pv
    rows = []
    for i, (rid, phi, tie) in enumerate(regimes, start=1):
        sel = chosen[rid]
        n = len(sel)
        st = pd.DataFrame([x for _, x in sel], columns=["tie_fraction", "ac1", "dfa_alpha"])
        k_max = fp_threshold(n) if n else 0
        for name in INDEX_FNS:
            fp = int(sum(pvals[(i, m)][name] <= FP_ALPHA for m, _ in sel))
            enough = n >= NULL_MIN_STATIONARY
            ok = enough and fp <= k_max
            rows.append({
                "regime": rid, "phi": phi, "tie_target": tie, "index": name, "n_realizations": drawn[rid],
                "n_stationary": n, "tie_min": st["tie_fraction"].min(), "tie_max": st["tie_fraction"].max(),
                "ac1_min": st["ac1"].min(), "ac1_max": st["ac1"].max(), "dfa_min": st["dfa_alpha"].min(),
                "dfa_max": st["dfa_alpha"].max(), "false_positives": fp, "fp_rate": fp / n if n else np.nan,
                "fp_threshold_count": k_max, "fp_threshold_rate": k_max / n if n else np.nan, "validated": ok,
                "reason": "" if ok else (f"{n} stationary realizations < {NULL_MIN_STATIONARY}" if not enough else
                                         f"false-positive rate {fp / n:.3f} > {k_max / n:.3f}"),
            })
    return pd.DataFrame(rows)


def match_regimes(stats: dict[str, float], null: pd.DataFrame) -> list[str]:
    """Regimes whose stationary realizations span the window's tie fraction, lag-1 AC and DFA exponent."""
    reg = null.drop_duplicates("regime")
    if not all(np.isfinite(stats[k]) for k in ("tie_fraction", "ac1", "dfa_alpha")):
        return []
    inside = ((reg["tie_min"] <= stats["tie_fraction"]) & (stats["tie_fraction"] <= reg["tie_max"])
              & (reg["ac1_min"] <= stats["ac1"]) & (stats["ac1"] <= reg["ac1_max"])
              & (reg["dfa_min"] <= stats["dfa_alpha"]) & (stats["dfa_alpha"] <= reg["dfa_max"]))
    return reg.loc[inside.fillna(False).astype(bool), "regime"].tolist()


def null_gate(name: str, matched: list[str], null: pd.DataFrame) -> tuple[bool, str]:
    """An index is reported only if the window lies in at least one regime and every matched regime is validated."""
    if not matched:
        return False, "no matched null regime"
    g = null[(null["index"] == name) & null["regime"].isin(matched)]
    bad = g[~g["validated"].astype(bool)]
    if len(bad):
        return False, "; ".join(f"{r.regime} (phi {r.phi:g}, tie {r.tie_target:g}): {r.reason}"
                                for r in bad.itertuples(index=False))
    return True, ", ".join(matched)


def null_rows(null: pd.DataFrame) -> pd.DataFrame:
    """Matched-null regimes as validation-table rows (role: gate for the pH windows)."""
    out = null.copy()
    out.insert(0, "check", "iaaft_matched_null_false_positives_" + out["index"])
    out["parameters"] = [f"{r.regime}: quantized AR(1) phi={r.phi:g}, target tie {r.tie_target:g}, n={SURR_LENGTH}; "
                         f"{r.n_stationary} stationary of {r.n_realizations} realizations"
                         for r in out.itertuples(index=False)]
    out["value"] = out["false_positives"]
    out["expected"] = [f"<= {k} of {n} (rate <= {FP_ALPHA} + binomial margin)"
                       for k, n in zip(out["fp_threshold_count"], out["n_stationary"])]
    out["required"] = False
    out["pass"] = out["validated"].astype(bool)
    out["role"] = "gate: pH windows"
    out["tie_fraction"] = np.nan
    out["note"] = [f"{N_SURR} IAAFT surrogates, rank test at {FP_ALPHA}; envelope tie {r.tie_min:.3f}-{r.tie_max:.3f}, "
                   f"lag-1 AC {r.ac1_min:.4f}-{r.ac1_max:.4f}, DFA {r.dfa_min:.3f}-{r.dfa_max:.3f}"
                   + (f"; NOT VALIDATED: {r.reason}" if r.reason else "") for r in out.itertuples(index=False)]
    return out


def null_from_validation(validation: pd.DataFrame) -> pd.DataFrame:
    cols = ["regime", "phi", "tie_target", "index", "n_realizations", "n_stationary", "tie_min", "tie_max", "ac1_min",
            "ac1_max", "dfa_min", "dfa_max", "false_positives", "fp_rate", "fp_threshold_count", "fp_threshold_rate",
            "validated", "reason"]
    g = validation[validation["check"].str.startswith("iaaft_matched_null_false_positives_")]
    if not len(g):
        raise EstimatorValidationError("estimator validation has no matched-null regime rows")
    out = g[cols].reset_index(drop=True)
    out["validated"] = out["validated"].astype(str).str.lower().eq("true")
    out["reason"] = out["reason"].fillna("").astype(str)
    return out


def _tie_sweep_task(task: tuple[int, float | None, int]) -> tuple[float, dict]:
    level, sd, m = task
    rng = _rng(2, level, m)
    z = _ar1(SURR_LENGTH, AR_PHI, rng)
    if sd is not None:
        z = np.round(z * sd / z.std() / QUANT_STEP) * QUANT_STEP
    out = E.surrogate_indices(z, INDEX_FNS, rng, N_SURR, TWO_SIDED)
    return E.tie_fraction(z), {name: (pv, obs - np.median(vals)) for name, (obs, vals, pv) in out.items()}


def _power_task(m: int) -> dict[str, float]:
    rng = _rng(3, m)
    z = E.simulate_ring(0.6, 0.1, SURR_LENGTH, rng) + rng.normal(0, 0.3, SURR_LENGTH)
    return {name: pv for name, (_, _, pv) in E.surrogate_indices(z, INDEX_FNS, rng, N_SURR, TWO_SIDED).items()}


def validate_estimators() -> pd.DataFrame:
    """The validation suite; the surrogate-heavy parts run on a process pool (results do not depend on it)."""
    pool = _pool()
    try:
        return _validate_estimators(pool)
    finally:
        if pool is not None:
            pool.shutdown()


def _validate_estimators(pool: ProcessPoolExecutor | None) -> pd.DataFrame:
    rows = []

    def add(check, params, value, expected, lo=np.nan, hi=np.nan, required=True, passed=True, note=""):
        rows.append({"check": check, "parameters": params, "value": value, "expected": expected,
                     "ci_low": lo, "ci_high": hi, "required": required, "pass": bool(passed), "note": note})

    for k, (p, q) in enumerate(MARKOV_PAIRS):
        rng = _rng(1, k)
        s = E.simulate_ring(p, q, MARKOV_STEPS, rng)
        true = E.ring_entropy_production(p, q)
        r = E.markov_ep_with_ci(s, 3, MARKOV_BLOCK, MARKOV_BOOT, MARKOV_LEVEL, rng)
        add("markov_ring_entropy_production", f"p={p}, q={q}", r["estimate"], true, r["ci_low"], r["ci_high"],
            passed=r["ci_low"] <= true <= r["ci_high"],
            note=f"{MARKOV_LEVEL:.0%} moving-block basic bootstrap; {MARKOV_STEPS} steps")
        if p == q:
            pval = float(chi2.sf(r["g_stat"], r["g_dof"]))
            add("markov_ring_detailed_balance_g_test", f"p={p}, q={q}", pval, "> 0.01", passed=pval > 0.01,
                note="Bowker-type symmetry G-test of transition counts (conservative)")
        y = np.round((QUANT_LEVEL_SPACING * s + rng.normal(0, QUANT_NOISE_SD, len(s))) / QUANT_STEP) * QUANT_STEP
        sh = np.clip(np.round(y / QUANT_LEVEL_SPACING), 0, 2).astype(np.int64)
        rq = E.markov_ep_with_ci(sh, 3, MARKOV_BLOCK, MARKOV_BOOT, MARKOV_LEVEL, rng)
        params = (f"p={p}, q={q}; levels {QUANT_LEVEL_SPACING} s + N(0, {QUANT_NOISE_SD}), quantized {QUANT_STEP}; "
                  f"misclassified {np.mean(sh != s):.3f}; tie fraction {E.tie_fraction(y):.3f}")
        if p == q:
            add("quantized_observable_zero_at_p_eq_q", params, rq["estimate"], 0.0, rq["ci_low"], rq["ci_high"],
                passed=rq["ci_low"] <= 0.0 <= rq["ci_high"], note="quantization and noise must not create irreversibility")
        else:
            add("quantized_observable_bias", params, rq["estimate"] - true, "reported", rq["ci_low"] - true,
                rq["ci_high"] - true, required=False,
                note=f"relative bias {(rq['estimate'] - true) / true:+.1%} (misclassification makes transitions look reversible)")
        add("quantized_observable_ordinal_kl", params, E.ordinal_irreversibility(y, ORDINAL_D), "reported",
            required=False, note="tie-aware ordinal index on the quantized observable")

    sweep = _map(_tie_sweep_task, [(level, sd, m) for level, sd in enumerate(TIE_SWEEP_SD)
                                   for m in range(SURR_REALIZATIONS)], pool)
    for level, sd in enumerate(TIE_SWEEP_SD):
        res = {name: {"p": [], "delta": []} for name in INDEX_FNS}
        ties = []
        for tie, out in sweep[level * SURR_REALIZATIONS: (level + 1) * SURR_REALIZATIONS]:
            ties.append(tie)
            for name, (pv, delta) in out.items():
                res[name]["p"].append(pv)
                res[name]["delta"].append(delta)
        required = level < N_REQUIRED_TIE_LEVELS
        variant = "continuous" if sd is None else f"quantized {QUANT_STEP} at sd {sd}"
        for name, r in res.items():
            fp = int(np.sum(np.array(r["p"]) <= 0.05))
            params = f"AR(1) phi={AR_PHI}, n={SURR_LENGTH}, {variant}"
            add(f"iaaft_reversible_false_positives_{name}", params, fp, f"<= {MAX_FALSE_POSITIVES} of {SURR_REALIZATIONS}",
                required=required, passed=fp <= MAX_FALSE_POSITIVES,
                note=f"tie fraction {np.mean(ties):.3f}; {N_SURR} IAAFT surrogates, rank test at 0.05"
                     + ("" if required else "; defines the validated tie range"))
            rows[-1]["tie_fraction"] = float(np.mean(ties))
            rows[-1]["index"] = name
            if required:
                d = np.array(r["delta"])
                wp = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
                add(f"iaaft_reversible_index_minus_surrogates_{name}", params, float(np.median(d)),
                    "0 (Wilcoxon p > 0.01)", passed=wp > 0.01,
                    note=f"tie fraction {np.mean(ties):.3f}; Wilcoxon signed-rank p = {wp:.3g} on index − median(surrogates)")

    power = dict.fromkeys(INDEX_FNS, 0)
    for pvs in _map(_power_task, list(range(SURR_REALIZATIONS)), pool):
        for name, pv in pvs.items():
            power[name] += pv <= 0.05
    for name, k in power.items():
        add(f"positive_control_power_{name}", "ring p=0.6, q=0.1 + N(0, 0.3)", k,
            f">= {MIN_POWER} of {SURR_REALIZATIONS}" if name == "ordinal_kl" else "reported",
            required=name == "ordinal_kl", passed=k >= MIN_POWER if name == "ordinal_kl" else True)

    for omega, label in ((0.0, "reversible"), (OU_OMEGA, "rotational")):
        rng = _rng(4, int(omega * 10))
        x, y = _ou_2d(OU_STEPS, OU_DT, omega, rng)
        valid = np.ones(OU_STEPS - 1, dtype=bool)
        r = E.circulation_with_ci(x, y, valid, 200, 500, LEVEL, rng)
        expected = "CI includes 0" if omega == 0 else "CI > 0"
        ok = r["ci_low"] <= 0 <= r["ci_high"] if omega == 0 else r["ci_low"] > 0
        add(f"flux_circulation_{label}_OU", f"2D OU, omega={omega}, {OU_STEPS} steps", r["estimate"], expected,
            r["ci_low"], r["ci_high"], passed=ok)

    rng = _rng(5)
    pe = E.permutation_entropy_ties(rng.random(20_000))
    target = np.log(6) / np.log(13)
    add("permutation_entropy_iid", "uniform iid, d=3", pe, round(target, 4), passed=abs(pe - target) < 0.01,
        note="iid continuous data use 6 of the 13 tie-aware patterns: ln 6 / ln 13")
    sc = np.unique(np.logspace(1, np.log10(500), 12).astype(int))
    a_w = E.dfa(rng.normal(size=8000), sc)[0]
    a_r = E.dfa(np.cumsum(rng.normal(size=8000)), sc)[0]
    add("dfa_white_noise", "n=8000", a_w, 0.5, passed=abs(a_w - 0.5) < 0.1)
    add("dfa_random_walk", "n=8000", a_r, 1.5, passed=abs(a_r - 1.5) < 0.1)
    se = E.sample_entropy(rng.normal(size=2000), SAMPEN_M, 0.2)
    add("sample_entropy_white_noise", "m=2, r=0.2 sd, n=2000", se, 2.2, passed=abs(se - 2.2) < 0.2,
        note="expected about 2.2 for Gaussian white noise")
    v = pd.DataFrame(rows)
    v["role"] = np.where(v["required"], "required", "reported")
    v.loc[v["check"].str.startswith("iaaft_reversible_false_positives_") & ~v["required"], "role"] = "gate: tie range"
    return pd.concat([v, null_rows(matched_null_table(pool))], ignore_index=True)


def validated_tie_range(validation: pd.DataFrame) -> dict[str, float]:
    """Largest tie fraction up to which an index passes the false-positive check at every lower level."""
    out = {}
    fp = validation[validation["check"].str.startswith("iaaft_reversible_false_positives_")]
    for name in INDEX_FNS:
        g = fp[fp["index"] == name].sort_values("tie_fraction")
        limit = -1.0
        for tie, passed in zip(g["tie_fraction"], g["pass"]):
            if not passed:
                break
            limit = float(tie)
        out[name] = limit
    return out


# =============================================================================
# Data windows
# =============================================================================
def pulse_segments(cfg: Config) -> list[tuple[int, float, float]]:
    b = [START_D, *cfg.pulse_times_d, END_D]
    return [(k, a, c) for k, (a, c) in enumerate(itertools.pairwise(b))]


def subtract_intervals(intervals: list[tuple[float, float]], cuts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    out = []
    for a, b in intervals:
        pieces = [(a, b)]
        for c0, c1 in cuts:
            nxt = []
            for p0, p1 in pieces:
                if c1 <= p0 or c0 >= p1:
                    nxt.append((p0, p1))
                    continue
                if c0 > p0:
                    nxt.append((p0, c0))
                if c1 < p1:
                    nxt.append((c1, p1))
            pieces = nxt
        out += pieces
    return out


def allowed_intervals(cfg: Config, series: str, held: pd.DataFrame, seg_table: pd.DataFrame) -> list[tuple[int, float, float]]:
    """(segment, start, end) intervals each series may be analysed in."""
    out = []
    for k, a, b in pulse_segments(cfg):
        if series in ("pH", "EC"):
            pieces = [(a, b)]
        elif series == "ORP":
            cuts = [(t - EVENT_MARGIN_D[0], t + EVENT_MARGIN_D[1]) for t in ORP_ELECTRODE_EVENTS_D]
            pieces = subtract_intervals([(a, b)], cuts)
        elif series == "NH4":
            comp = seg_table[(seg_table["status"] == "computed") & (seg_table["segment"] == k)]
            base = [(r.fit_start_d, r.fit_end_d) for r in comp.itertuples(index=False)]
            pieces = subtract_intervals(base, list(zip(held["t_start_d"], held["t_end_d"])))
        else:
            raise KeyError(series)
        out += [(k, p0, p1) for p0, p1 in pieces]
    return out


def nh4_uniform(nh4: pd.DataFrame) -> pd.DataFrame:
    """NH4 on 20 s time-weighted bins (empty bins dropped, never filled)."""
    step = NH4_BIN_S / 86400.0
    edges = np.round(np.arange(0.0, END_D + step, step), 9)
    t = nh4["t_d"].to_numpy()
    mean, count = timebase.bin_time_weighted(t, nh4["nh4_mgNH4_L"].to_numpy(), timebase.sample_weights_d(t), edges)
    ok = count > 0
    return pd.DataFrame({"t_d": 0.5 * (edges[:-1] + edges[1:])[ok], "value": mean[ok]})


def stationarity(x: np.ndarray) -> dict:
    """ADF (unit root) and KPSS (level stationarity) on the raw window.

    No detrending: subtracting a trend from a quantized series turns plateaus into a
    sawtooth, and leftover curvature in a relaxation is itself time-asymmetric; both
    create spurious irreversibility (found in M5, see docs/progress.md).
    """
    if np.std(x) == 0:
        return {"transform": "not computed: constant window"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", InterpolationWarning)
        adf_p = float(adfuller(x, autolag="AIC", result_object=False)[1])
        kpss_p = float(kpss(x, regression="c", nlags="auto", result_object=False)[1])
    out = {"adf_p": adf_p, "kpss_p": kpss_p}
    if adf_p < ADF_ALPHA and kpss_p > KPSS_ALPHA:
        out["transform"] = "raw"
        out["series"] = x
    else:
        out["transform"] = "not computed: non-stationary (ADF p >= 0.05 or KPSS p <= 0.05)"
    return out


OUTSIDE_REGIME = "not computed: outside validated regime"


def passes_stationarity(x: np.ndarray) -> bool:
    """The accept decision of :func:`stationarity` (ADF p < 0.05 and KPSS p > 0.05), evaluated KPSS first:
    KPSS costs ~0.1 ms against ~14 ms for the autolag ADF, and most persistent null draws fail it."""
    if np.std(x) == 0:
        return False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", InterpolationWarning)
        if not float(kpss(x, regression="c", nlags="auto", result_object=False)[1]) > KPSS_ALPHA:
            return False
        return float(adfuller(x, autolag="AIC", result_object=False)[1]) < ADF_ALPHA


def analyse_window(x: np.ndarray, rng: np.random.Generator, tie_range: dict[str, float],
                   null: pd.DataFrame | None = None) -> dict:
    """Stationarity, then the irreversibility indices against IAAFT surrogates and complexity measures.

    With ``null`` (the matched-null regime table), an index is reported only if the window's tie fraction,
    lag-1 autocorrelation and DFA exponent fall in at least one null regime and every regime they fall in
    is validated for that index; otherwise ``not computed: outside validated regime (...)``.
    """
    st = stationarity(x)
    row = {k: v for k, v in st.items() if k != "series"}
    row["tie_fraction"] = E.tie_fraction(x)
    if "series" not in st:
        return row
    z = st["series"]
    ws = window_stats(z)
    row["ac1"] = ws["ac1"]
    matched = match_regimes(ws, null) if null is not None else []
    if null is not None:
        row["null_regimes"] = ";".join(matched) if matched else "none"
    for name, (obs, vals, pv) in E.surrogate_indices(z, INDEX_FNS, rng, N_SURR, TWO_SIDED).items():
        if row["tie_fraction"] > tie_range[name]:
            row[f"{name}_status"] = (f"not computed: tie fraction {row['tie_fraction']:.3f} above validated "
                                     f"range ({tie_range[name]:.3f})")
            continue
        if null is not None:
            ok, why = null_gate(name, matched, null)
            if not ok:
                row[f"{name}_status"] = f"{OUTSIDE_REGIME} ({why})"
                continue
            row[f"{name}_null_regimes"] = why
        row[f"{name}_status"] = "computed"
        row[name] = obs
        row[f"{name}_surr_median"] = float(np.median(vals))
        row[f"{name}_p"] = pv
    _raw, norm = E.moment_asymmetry(z, list(ASYM_LAGS))
    for lag, v in zip(ASYM_LAGS, norm):
        row[f"asym_norm_lag{lag}"] = v
    row["perm_entropy_ties"] = E.permutation_entropy_ties(z, ORDINAL_D)
    r = SAMPEN_R_SD * np.std(z)
    mse = E.multiscale_entropy(z, MSE_SCALES, SAMPEN_M, r) if r > 0 else [np.nan] * len(MSE_SCALES)
    for s, v in zip(MSE_SCALES, mse):
        row[f"sampen_scale{s}"] = v
    row["dfa_alpha"] = ws["dfa_alpha"]
    return row


def window_table(series: str, role: str, t: np.ndarray, x: np.ndarray, intervals, series_id: int,
                 tie_range: dict[str, float], null: pd.DataFrame | None = None) -> pd.DataFrame:
    dt = float(np.median(np.diff(t)))
    expected = WINDOW_D / dt
    rows = []
    w_id = 0
    for seg, a, b in intervals:
        n_win = int(np.floor((b - a) / WINDOW_D + 1e-9))
        for i in range(n_win):
            w0, w1 = a + i * WINDOW_D, a + (i + 1) * WINDOW_D
            m = (t >= w0) & (t < w1)
            row = {"series": series, "role": role, "segment": seg, "window": w_id, "t_start_d": w0, "t_end_d": w1,
                   "n": int(m.sum())}
            if m.sum() < MIN_WINDOW_FILL * expected:
                row["transform"] = f"not computed: window {m.sum() / expected:.0%} filled"
            else:
                row.update(analyse_window(x[m], _rng(10, series_id, w_id), tie_range, null))
            rows.append(row)
            w_id += 1
    return pd.DataFrame(rows)


def verify_null_gate(windows: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    """Required checks, recomputed from the tables: every stationary pH window was matched against the null,
    and every pH index reported lies inside validated regimes only (re-matched from the stored statistics)."""
    null = null_from_validation(validation)
    ph = windows[(windows["series"] == "pH") & (windows["transform"] == "raw")]
    rows = []
    for name in INDEX_FNS:
        status = ph.get(f"{name}_status", pd.Series(index=ph.index, dtype=object)).fillna("")
        reported = ph[status == "computed"]
        violations = []
        for r in reported.itertuples(index=False):
            matched = match_regimes({"tie_fraction": r.tie_fraction, "ac1": r.ac1, "dfa_alpha": r.dfa_alpha}, null)
            ok, _ = null_gate(name, matched, null)
            if not ok:
                violations.append(f"{r.t_start_d:.2f}")
        gated = status.eq("computed") | status.str.startswith(OUTSIDE_REGIME) | status.str.startswith(
            "not computed: tie fraction")
        unmatched = int((~gated).sum())
        rows.append({"check": f"matched_null_gate_pH_{name}", "parameters": f"{len(ph)} stationary pH windows",
                     "value": len(violations) + unmatched, "expected": 0, "ci_low": np.nan, "ci_high": np.nan,
                     "required": True, "pass": not violations and unmatched == 0, "role": "required",
                     "note": (f"{len(reported)} reported inside validated regimes, "
                              f"{int(status.str.startswith(OUTSIDE_REGIME).sum())} outside validated regime"
                              + (f"; violations at {', '.join(violations)} d" if violations else "")
                              + (f"; {unmatched} window(s) not gated" if unmatched else ""))})
    return pd.DataFrame(rows)


def summarize(windows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (series, role, seg), g in windows.groupby(["series", "role", "segment"], sort=False):
        ok = g[g["transform"] == "raw"]
        row = {"series": series, "role": role, "segment": seg, "n_windows": len(g), "n_stationary": len(ok),
               "stationary_fraction": len(ok) / len(g) if len(g) else np.nan,
               "tie_fraction_median": float(g["tie_fraction"].median()) if "tie_fraction" in g else np.nan}
        for name in INDEX_FNS:
            done = ok[ok[f"{name}_status"] == "computed"] if f"{name}_status" in ok else ok.iloc[0:0]
            row[f"{name}_windows"] = len(done)
            if len(done):
                k = int((done[f"{name}_p"] <= 0.05).sum())
                row[f"{name}_median"] = float(done[name].median())
                row[f"{name}_minus_surr_median"] = float((done[name] - done[f"{name}_surr_median"]).median())
                row[f"{name}_significant_windows"] = k
                row[f"{name}_binomial_p"] = float(binomtest(k, len(done), 0.05, alternative="greater").pvalue)
            elif len(ok) and f"{name}_status" in ok and ok[f"{name}_status"].str.startswith(OUTSIDE_REGIME).any():
                row[f"{name}_median"] = f"{OUTSIDE_REGIME} or tie fraction above validated range in every stationary window"
            elif len(ok):
                row[f"{name}_median"] = "not computed: tie fraction above validated range in every stationary window"
            else:
                row[f"{name}_median"] = "not computed: no stationary window"
        for col in ("perm_entropy_ties", "sampen_scale1", "sampen_scale2", "sampen_scale4", "dfa_alpha"):
            row[f"{col}_median"] = float(ok[col].median()) if len(ok) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Onsager regression, stability, flux (0.01 d grid)
# =============================================================================
def _acf(f: np.ndarray, max_lag: int) -> np.ndarray:
    f = f - f.mean()
    v = np.dot(f, f) / len(f)
    return np.array([np.dot(f[: len(f) - k], f[k:]) / len(f) / v for k in range(max_lag + 1)])


def onsager(cfg: Config, grid: pd.DataFrame, pulse_table: pd.DataFrame, col: str = "pH") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Late-tail relaxation R(τ) vs normalized autocorrelation C(τ) of late fluctuations, per pulse.

    The signal-to-noise ratio uses the pulse response measured in M2 (1 to 3 h after minus 3 to 1 h
    before), because the deviation of the tail from the late level also contains the slow drift of
    the segment, which cannot be attributed to the pulse.
    """
    response = dict(zip(pulse_table["pulse"], pulse_table["ph_obs_change"]))
    t = grid["t_mid_d"].to_numpy()
    x = grid[col].to_numpy()
    lag_rows, summary = [], []
    bounds = [*cfg.pulse_times_d, END_D]
    for k, (t0, t1) in enumerate(itertools.pairwise(bounds), start=1):
        seg_len = t1 - t0
        late_start = t1 - max(ONSAGER_LATE_FRACTION * seg_len, 0.5)
        late = (t >= late_start) & (t < t1)
        tl = t[late]
        f = x[late] - np.polyval(np.polyfit(tl, x[late], 1), tl)
        x_late = float(np.median(x[late]))
        max_lag = round(min(ONSAGER_MAX_LAG_D, (t1 - late_start) / 3) / 0.01)
        tail = (t >= t0 + ONSAGER_TAIL_START_D) & (t < t0 + ONSAGER_TAIL_START_D + (max_lag + 1) * 0.01 - 1e-9)
        dx = x[tail][: max_lag + 1] - x_late
        snr = abs(response[k]) / np.std(f)
        drift_ratio = abs(dx[0]) / abs(response[k]) if response[k] != 0 else np.inf
        st = stationarity(f)
        stationary = "series" in st
        c = _acf(f, max_lag)
        r = dx / dx[0]
        rng = _rng(20, k)
        smooth = pd.Series(dx).rolling(5, center=True, min_periods=1).median().to_numpy()
        e = dx - smooth
        boot = np.full((ONSAGER_BOOT, max_lag + 1), np.nan)
        for bi in range(ONSAGER_BOOT):
            cb = _acf(f[E.block_bootstrap_rows(len(f), ONSAGER_BLOCK, rng)], max_lag)
            db = smooth + e[E.block_bootstrap_rows(len(e), 5, rng)]
            with np.errstate(divide="ignore", invalid="ignore"):
                boot[bi] = (db / db[0]) / cb
        for lag in range(max_lag + 1):
            valid = c[lag] >= ONSAGER_MIN_C and stationary
            lo, hi = np.nanquantile(boot[:, lag], [0.025, 0.975]) if valid else (np.nan, np.nan)
            lag_rows.append({"pulse": k, "lag_d": lag * 0.01, "relaxation_R": r[lag], "acf_C": c[lag],
                             "ratio_R_over_C": r[lag] / c[lag] if valid else np.nan, "ratio_ci_low": lo,
                             "ratio_ci_high": hi})
        rows_k = pd.DataFrame(lag_rows)
        rows_k = rows_k[(rows_k["pulse"] == k) & rows_k["ratio_R_over_C"].notna()]
        sufficient = bool(snr >= ONSAGER_MIN_SNR)
        tests = "" if "adf_p" not in st else f"; ADF p {st['adf_p']:.3g}, KPSS p {st['kpss_p']:.3g}"
        nonstat = f"not computed: non-stationary late fluctuation window ({st['transform'].removeprefix('not computed: ')}{tests})"
        if not sufficient:
            verdict = f"not interpretable: pulse-response SNR {snr:.2f} < {ONSAGER_MIN_SNR}"
            if not stationary:
                verdict += f"; {nonstat}"
        elif not stationary:
            verdict = nonstat
        elif len(rows_k) < 3:
            verdict = "not interpretable: fewer than 3 lags with C >= 0.2"
        else:
            excl = ((rows_k["ratio_ci_low"] > 1) | (rows_k["ratio_ci_high"] < 1)).mean()
            verdict = f"regression violated at {excl:.0%} of lags" if excl > 0 else "consistent with regression"
        summary.append({"pulse": k, "series": col, "tail_start_d": t0 + ONSAGER_TAIL_START_D, "late_start_d": late_start,
                        "late_end_d": t1, "pulse_response": response[k], "tail_deviation_x0": dx[0],
                        "tail_deviation_over_pulse_response": drift_ratio, "fluctuation_sd": float(np.std(f)),
                        "snr_pulse_response": snr, "snr_sufficient": sufficient, "max_lag_d": max_lag * 0.01,
                        "late_adf_p": st.get("adf_p", np.nan), "late_kpss_p": st.get("kpss_p", np.nan),
                        "late_window_stationary": stationary, "verdict": verdict})
    return pd.DataFrame(lag_rows), pd.DataFrame(summary)


def orp_event_cuts() -> list[tuple[float, float]]:
    """Intervals excluded around the ORP electrode events (M5 item 3)."""
    return [(t - EVENT_MARGIN_D[0], t + EVENT_MARGIN_D[1]) for t in ORP_ELECTRODE_EVENTS_D]


def stability_pieces(series: str, a: float, b: float) -> list[tuple[float, float]]:
    """Contiguous pieces of [a, b) a stability series is computed on: ORP excludes its electrode events."""
    return subtract_intervals([(a, b)], orp_event_cuts()) if series == "ORP" else [(a, b)]


def stability(cfg: Config, grid: pd.DataFrame) -> pd.DataFrame:
    """Rolling variance and lag-1 autocorrelation before each pulse; rolling windows never span an excluded
    interval, and each piece is linearly detrended on its own."""
    t = grid["t_mid_d"].to_numpy()
    rows = []
    prev = START_D
    for k, tp in enumerate(cfg.pulse_times_d, start=1):
        a = max(prev, tp - STABILITY_SPAN_D)
        for col, label in (("pH", "pH"), ("orp_mV_vs_Ag", "ORP")):
            pieces = stability_pieces(label, a, tp)
            excluded = [(c0, c1) for c0, c1 in (orp_event_cuts() if label == "ORP" else []) if c1 > a and c0 < tp]
            var, ac1, tc = [], [], []
            for p0, p1 in pieces:
                m = (t >= p0) & (t < p1)
                if m.sum() < STABILITY_WINDOW + 2:
                    continue
                tt, x = t[m], grid[col].to_numpy()[m]
                z = x - np.polyval(np.polyfit(tt, x, 1), tt)
                for i in range(len(z) - STABILITY_WINDOW + 1):
                    w = z[i: i + STABILITY_WINDOW]
                    var.append(np.var(w))
                    ac1.append(np.corrcoef(w[:-1], w[1:])[0, 1] if np.std(w) > 0 else np.nan)
                    tc.append(tt[i + STABILITY_WINDOW - 1])
            row = {"pulse": k, "series": label, "from_d": a, "to_d": tp, "n_pieces": len(pieces),
                   "excluded_d": "; ".join(f"{c0:.2f}-{c1:.2f}" for c0, c1 in excluded),
                   "rolling_end_min_d": min(tc, default=np.nan), "rolling_end_max_d": max(tc, default=np.nan),
                   "n_rolling": len(var)}
            if len(var) >= 3:
                kv = kendalltau(tc, var)
                ka = kendalltau(tc, ac1, nan_policy="omit")
                row.update(variance_kendall_tau=kv.statistic, variance_p=kv.pvalue, ac1_kendall_tau=ka.statistic,
                           ac1_p=ka.pvalue)
            else:
                row.update(variance_kendall_tau=f"not computed: {len(var)} rolling windows", variance_p=np.nan,
                           ac1_kendall_tau=np.nan, ac1_p=np.nan)
            row["note"] = "EXPLORATORY: rolling windows overlap, so p-values are anti-conservative"
            rows.append(row)
        prev = tp
    return pd.DataFrame(rows)


def flux(cfg: Config, grid: pd.DataFrame, intervals: dict) -> pd.DataFrame:
    t = grid["t_mid_d"].to_numpy()
    rows = []
    for plane, ycol, series in (("pH-ORP", "orp_mV_vs_Ag", "ORP"), ("pH-NH4", "nh4_mgNH4_L", "NH4")):
        for k, a, b in pulse_segments(cfg):
            pieces = [(p0, p1) for s, p0, p1 in intervals[series] if s == k]
            xs, ys, valid = [], [], []
            n_pieces = n_stationary = 0
            for p0, p1 in pieces:
                m = (grid["t_start_d"].to_numpy() >= p0 - 1e-9) & (grid["t_end_d"].to_numpy() <= p1 + 1e-9)
                if m.sum() < 3:
                    continue
                n_pieces += 1
                tt = t[m]
                dx, dy = (grid[col].to_numpy()[m] - np.polyval(np.polyfit(tt, grid[col].to_numpy()[m], 1), tt)
                          for col in ("pH", ycol))
                # ADF/KPSS on each input window (both coordinates, as passed to the estimator)
                if m.sum() < FLUX_MIN_INCREMENTS or not all("series" in stationarity(v) for v in (dx, dy)):
                    continue
                n_stationary += 1
                xs.append(dx)
                ys.append(dy)
                valid.append(np.r_[np.ones(m.sum() - 1, dtype=bool), False])
            row = {"plane": plane, "segment": k, "from_d": a, "to_d": b, "n_input_windows": n_pieces,
                   "n_stationary_input_windows": n_stationary}
            if not xs or sum(v.sum() for v in valid) < FLUX_MIN_INCREMENTS:
                row["circulation"] = (f"not computed: fewer than {FLUX_MIN_INCREMENTS} increments in stationary input "
                                      f"windows ({n_stationary} of {n_pieces} windows pass ADF/KPSS)")
                rows.append(row)
                continue
            x = np.concatenate(xs)
            y = np.concatenate(ys)
            x, y = (x - x.mean()) / x.std(), (y - y.mean()) / y.std()
            vmask = np.concatenate(valid)[:-1]
            r = E.circulation_with_ci(x, y, vmask, FLUX_BLOCK, FLUX_BOOT, LEVEL, _rng(30, k, len(plane)))
            row.update(circulation=r["estimate"], ci_low=r["ci_low"], ci_high=r["ci_high"],
                       n_increments=r["n_increments"], significant=bool(r["ci_low"] > 0 or r["ci_high"] < 0))
            rows.append(row)
    return pd.DataFrame(rows)


def electrode_dynamics(cfg: Config) -> pd.DataFrame:
    try:
        cfg.measurement("electrode_step_response")
        status = "not computed: deconvolution not implemented"
    except MissingInputError as exc:
        status = str(exc)
    return pd.DataFrame([{"quantity": "NH4 relaxation times, deconvolved (process time scales)", "status": status,
                          "note": "NH4 relaxation is electrode-limited (S5) until the step response is measured"}])


# =============================================================================
# Stage
# =============================================================================
@dataclass
class DynamicsResult:
    validation: pd.DataFrame
    windows: pd.DataFrame | None = None
    summary: pd.DataFrame | None = None
    onsager_lags: pd.DataFrame | None = None
    onsager_summary: pd.DataFrame | None = None
    stability: pd.DataFrame | None = None
    flux: pd.DataFrame | None = None
    electrode: pd.DataFrame | None = None


class EstimatorValidationError(RuntimeError):
    """A required estimator validation check failed; no plant-data result is computed."""


def load_inputs(out_dir: Path) -> dict:
    files = {"ph": "processed/native_ph.parquet", "orp": "processed/native_orp.parquet",
             "nh4": "processed/native_nh4.parquet", "cond": "processed/native_cond.parquet",
             "grid": "processed/grid_0p01d.parquet", "held": "tables/nh4_held_segments.csv",
             "segments": "tables/segments.csv", "pulses": "tables/sensor_validation.csv"}
    out = {}
    for k, rel in files.items():
        p = Path(out_dir) / rel
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run `python -m bnr_neq run-all` stages before dynamics")
        out[k] = pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)
    return out


def required_failures(validation: pd.DataFrame) -> list[str]:
    req = validation["required"].astype(str).str.lower().eq("true")
    ok = validation["pass"].astype(str).str.lower().eq("true")
    return validation.loc[req & ~ok, "check"].tolist()


def compute(cfg: Config, inputs: dict, validation: pd.DataFrame | None = None) -> DynamicsResult:
    val = validate_estimators() if validation is None else validation
    val = val[~val["check"].str.startswith("matched_null_gate_")].reset_index(drop=True)  # recomputed below
    failed = required_failures(val)
    if failed:
        raise EstimatorValidationError(f"estimator validation failed: {failed}")
    tie_range = validated_tie_range(val)
    null = null_from_validation(val)
    intervals = {s: allowed_intervals(cfg, s, inputs["held"], inputs["segments"]) for s in ("pH", "ORP", "NH4", "EC")}
    nh4u = nh4_uniform(inputs["nh4"])
    data = {
        "pH": ("process", inputs["ph"]["t_d"].to_numpy(), inputs["ph"]["pH"].to_numpy()),
        "ORP": ("process (relative)", inputs["orp"]["t_d"].to_numpy(), inputs["orp"]["orp_mV_vs_Ag"].to_numpy()),
        "NH4": ("process (offset-corrected, provisional)", nh4u["t_d"].to_numpy(), nh4u["value"].to_numpy()),
        "EC": ("probe behaviour only (S3)", inputs["cond"]["t_d"].to_numpy(), inputs["cond"]["EC_mS_cm"].to_numpy()),
    }
    windows = pd.concat([window_table(s, role, t, x, intervals[s], i, tie_range, null if s == "pH" else None)
                         for i, (s, (role, t, x)) in enumerate(data.items())], ignore_index=True)
    val = pd.concat([val, verify_null_gate(windows, val)], ignore_index=True)
    failed = required_failures(val)
    if failed:
        raise EstimatorValidationError(f"matched-null gate failed: {failed}")
    lags, osum = onsager(cfg, inputs["grid"], inputs["pulses"])
    return DynamicsResult(val, windows, summarize(windows), lags, osum, stability(cfg, inputs["grid"]),
                          flux(cfg, inputs["grid"], intervals), electrode_dynamics(cfg))


# =============================================================================
# Figure, report, stage
# =============================================================================
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def figure(cfg: Config, r: DynamicsResult, png: Path, svg: Path) -> None:
    from .report import pyplot

    plt = pyplot()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), facecolor=SURFACE)
    # (a) Markov validation
    ax = axes[0, 0]
    _style(ax)
    v = r.validation
    clean = v[v["check"] == "markov_ring_entropy_production"]
    true = clean["expected"].astype(float).to_numpy()
    ax.errorbar(true, clean["value"], yerr=[clean["value"] - clean["ci_low"], clean["ci_high"] - clean["value"]],
                fmt="o", color=BLUE, ms=6, capsize=3, label="Transition estimator (99 % CI)")
    q = v[v["check"] == "quantized_observable_bias"]
    q_true = true[1:]
    ax.scatter(q_true, q_true + q["value"].astype(float), marker="s", color=ORANGE, s=36, zorder=3,
               label="Noisy observable quantized to 0.01")
    ax.plot([0, 1], [0, 1], color=INK2, lw=0.8, ls="--")
    ax.set_xlabel("Analytic (p − q) ln(p/q) per step", fontsize=8.5, color=INK2)
    ax.set_ylabel("Estimated entropy production per step", fontsize=8.5, color=INK2)
    ax.set_title("a  Estimator validation: three-state ring chain", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    # (b) pH ordinal KL per window
    ax = axes[0, 1]
    _style(ax)
    w = r.windows
    for series, color in (("pH", BLUE), ("ORP", AQUA), ("NH4", VIOLET)):
        g = w[(w["series"] == series) & (w.get("ordinal_kl_status", pd.Series(dtype=str)) == "computed")]
        mid = 0.5 * (g["t_start_d"] + g["t_end_d"])
        sig = g["ordinal_kl_p"] <= 0.05
        ax.scatter(mid[sig], g.loc[sig, "ordinal_kl"], color=color, s=30, label=f"{series}: significant vs IAAFT")
        ax.scatter(mid[~sig], g.loc[~sig, "ordinal_kl"], facecolor="none", edgecolor=color, s=30,
                   label=f"{series}: not significant")
    for tp in cfg.pulse_times_d:
        ax.axvline(tp, color=INK2, lw=0.8, ls=":")
    ax.set_xlim(0, END_D)
    ax.set_xlabel("Time (d); stationary 0.25 d windows only", fontsize=8.5, color=INK2)
    ax.set_ylabel("Ordinal-pattern KL (tie-aware, nats)", fontsize=8.5, color=INK2)
    ax.set_title("b  Time irreversibility in stationary windows", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="upper right", ncol=1)
    # (c) Onsager
    ax = axes[1, 0]
    _style(ax)
    s = r.onsager_summary
    ax.bar(s["pulse"].astype(str), s["snr_pulse_response"], color=BLUE, width=0.55)
    ax.axhline(ONSAGER_MIN_SNR, color=ORANGE, lw=1.5, ls="--", label=f"Required SNR ({ONSAGER_MIN_SNR:g})")
    ax.set_xlabel("Pulse", fontsize=8.5, color=INK2)
    ax.set_ylabel("|pH pulse response| / fluctuation SD", fontsize=8.5, color=INK2)
    ax.set_title("c  Onsager regression test: signal-to-noise", loc="left", fontsize=10, color=INK)
    ax.text(0.02, 0.86, "Below the line: regression test not interpretable", transform=ax.transAxes, fontsize=7.5,
            color=INK2, va="top")
    for xi, (val, resp) in enumerate(zip(s["snr_pulse_response"], s["pulse_response"])):
        ax.text(xi, val + 0.05, f"ΔpH {resp:+.2f}", ha="center", fontsize=7.5, color=INK2)
    ax.legend(fontsize=7.5, frameon=False, loc="center right")
    # (d) Flux
    ax = axes[1, 1]
    _style(ax)
    f = r.flux[pd.to_numeric(r.flux["circulation"], errors="coerce").notna()].copy()
    f["circulation"] = f["circulation"].astype(float)
    for i, (plane, color) in enumerate((("pH-ORP", AQUA), ("pH-NH4", VIOLET))):
        g = f[f["plane"] == plane]
        xpos = g["segment"].to_numpy() + (i - 0.5) * 0.25
        ax.errorbar(xpos, g["circulation"], yerr=[g["circulation"] - g["ci_low"], g["ci_high"] - g["circulation"]],
                    fmt="o", color=color, capsize=3, label=plane)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(5), ["base", "P1", "P2", "P3", "P4"])
    ax.set_xlabel("Segment", fontsize=8.5, color=INK2)
    ax.set_ylabel("Circulation rate (standardized units/step, 95 % CI)", fontsize=8.5, color=INK2)
    ax.set_title("d  Probability flux (Battle et al.)", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(png, dpi=150, facecolor=SURFACE)
    fig.savefig(svg, facecolor=SURFACE, metadata={"Date": None})
    plt.close(fig)


def onsager_conclusion(osum: pd.DataFrame) -> str:
    """The Onsager conclusion, derived from the per-pulse verdicts."""
    v = osum["verdict"].astype(str)
    usable = osum[~v.str.startswith(("not interpretable", "not computed"))]
    snr_bad = int((~osum["snr_sufficient"].astype(bool)).sum())
    nonstat = int((~osum["late_window_stationary"].astype(bool)).sum()) if "late_window_stationary" in osum else 0
    if len(usable) == 0:
        return (f"**No conclusion about Onsager regression (H4) is drawn**: the pulse-response SNR is insufficient for "
                f"{snr_bad} of {len(osum)} pulses (best {osum['snr_pulse_response'].max():.2f} < {ONSAGER_MIN_SNR:g}) and "
                f"the late fluctuation window is non-stationary for {nonstat} of {len(osum)}.")
    return ("Interpretable pulses: " + "; ".join(f"pulse {int(r.pulse)}: {r.verdict}" for r in usable.itertuples(index=False))
            + f". Not interpretable or not computed: {len(osum) - len(usable)} of {len(osum)}.")


def write_report(cfg: Config, r: DynamicsResult, path: Path) -> None:
    from .report import md_table

    v = r.validation
    failed = required_failures(v)
    if failed:
        raise EstimatorValidationError(f"required estimator validation checks failed: {failed}")
    req = v[v["required"].astype(bool)]
    tie_range = validated_tie_range(v)
    null = null_from_validation(v)
    core = v[~v["check"].str.startswith("iaaft_matched_null_false_positives_")]
    core_cols = ["check", "parameters", "value", "expected", "ci_low", "ci_high", "required", "role", "pass", "note"]
    regimes = null.pivot_table(index=["regime", "phi", "tie_target", "n_realizations", "n_stationary", "tie_min",
                                      "tie_max", "ac1_min", "ac1_max", "dfa_min", "dfa_max", "fp_threshold_count"],
                               columns="index", values="false_positives").reset_index()
    regimes.columns = [c if c not in INDEX_FNS else f"fp_{c}" for c in regimes.columns]
    valid = null.pivot_table(index="regime", columns="index", values="validated", aggfunc="first")
    regimes["validated_indices"] = [", ".join(n for n in INDEX_FNS if bool(valid.loc[rid, n])) or "none"
                                    for rid in regimes["regime"]]
    w = r.windows
    ph = w[(w["series"] == "pH") & (w["transform"] == "raw")]
    gate_lines = []
    for name in INDEX_FNS:
        st = ph[f"{name}_status"].fillna("") if f"{name}_status" in ph else pd.Series(dtype=str)
        gate_lines.append(f"- {name}: {int((st == 'computed').sum())} reported, "
                          f"{int(st.str.startswith(OUTSIDE_REGIME).sum())} `{OUTSIDE_REGIME}`, "
                          f"{int(st.str.startswith('not computed: tie fraction').sum())} above the tie range.")
    ph_view = ph[["t_start_d", "tie_fraction", "ac1", "dfa_alpha", "null_regimes",
                  *[c for name in INDEX_FNS for c in (f"{name}_status", f"{name}_p") if c in ph]]]
    s = r.summary
    view_cols = ["series", "role", "segment", "n_windows", "n_stationary", "tie_fraction_median"]
    for name in INDEX_FNS:
        view_cols += [f"{name}_windows", f"{name}_median", f"{name}_significant_windows", f"{name}_binomial_p"]
    view = s[[c for c in view_cols if c in s.columns]]
    cx = s[["series", "segment", "perm_entropy_ties_median", "sampen_scale1_median", "sampen_scale2_median",
            "sampen_scale4_median", "dfa_alpha_median"]]
    sweep_max = float(v.loc[v["check"].str.startswith("iaaft_reversible_false_positives_"), "tie_fraction"].max())
    flux_ok = pd.to_numeric(r.flux["circulation"], errors="coerce").notna()
    lines = [
        "# Dynamical and information-theoretic signatures (M5)",
        "",
        "Generated by `python -m bnr_neq dynamics`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Estimator validation (run first; the stage stops if a required check fails)",
        "",
        f"{int(req['pass'].astype(bool).sum())} of {len(req)} required checks pass; {len(v) - len(req)} further rows "
        "are reported or gate the data. Validated tie ranges (largest tie fraction up to which reversible quantized "
        f"AR(1) surrogate tests keep ≤ {MAX_FALSE_POSITIVES}/{SURR_REALIZATIONS} false positives at every lower level, "
        f"{N_SURR} surrogates): " + ", ".join(f"{k} ≤ {val:.3f}" for k, val in tie_range.items())
        + f". An index is not computed on a window whose tie fraction exceeds its range (the sweep reaches {sweep_max:.3f}).",
        "",
        md_table(core[[c for c in core_cols if c in core.columns]]),
        "",
        "Permutation entropy uses tie-aware patterns (13 symbols for d = 3); iid continuous data reach ln 6 / ln 13 = 0.699.",
        "",
        "### Matched null for the pH windows",
        "",
        (f"Reversible AR(1) with φ ∈ {{{', '.join(f'{x:g}' for x in NULL_PHI)}}}, scaled so that the series quantized to "
        f"{QUANT_STEP} has a target tie fraction in {{{', '.join(f'{x:g}' for x in NULL_TIE)}}}, n = {SURR_LENGTH} (one pH "
        f"window). Realizations are drawn until {NULL_TARGET_STATIONARY} pass the same ADF/KPSS filter as the data windows "
        f"(at most {NULL_MAX_REALIZATIONS}); each is tested with {N_SURR} IAAFT surrogates exactly as a data window. A regime "
        f"is validated for an index if it has ≥ {NULL_MIN_STATIONARY} stationary realizations and its false-positive count "
        f"at α = {FP_ALPHA} is at most the {FP_MARGIN_LEVEL:.0%} binomial quantile (rate ≤ {FP_ALPHA} + binomial margin). A pH "
        "window is matched to every regime whose stationary realizations span its tie fraction, lag-1 autocorrelation and "
        "DFA α; an index is reported only if at least one regime matches and every matched regime is validated, "
        f"otherwise `{OUTSIDE_REGIME}`. The required `matched_null_gate_pH_*` checks re-match every reported window."),
        "",
        md_table(regimes, digits=4),
        "",
        "Stationary pH windows (gating):",
        "",
        *gate_lines,
        "",
        md_table(ph_view),
        "",
        "## Application to the data",
        "",
        (f"Native series in {WINDOW_D} d windows inside allowed intervals: pH 10 s from day {START_D} (primary signal); "
        f"relative ORP 10 s, excluding {EVENT_MARGIN_D[0]} d before to {EVENT_MARGIN_D[1]} d after the electrode events at "
        f"{', '.join(f'{t:g}' for t in ORP_ELECTRODE_EVENTS_D)} d; NH4 on {NH4_BIN_S:g} s time-weighted bins, only inside the M3 "
        "offset-free slope windows and outside held/damped segments (provisional: electrode under-response, S1); "
        "conductivity as probe behaviour only (S3), never as solution conductivity. Windows never span a pulse."),
        "",
        (f"Stationarity: ADF p < {ADF_ALPHA} and KPSS p > {KPSS_ALPHA} on the raw window; nothing is detrended (detrending a "
        "quantized series creates a sawtooth and spurious irreversibility; an earlier version of this stage did so and "
        "was corrected before reporting). Each stationary window is tested against "
        f"{N_SURR} IAAFT surrogates (rank p, one-sided for the KL indices, two-sided for the asymmetry); per segment, "
        "the number of significant windows is compared with the 5 % expected under reversibility (binomial). The matched "
        "null gates pH only; ORP, NH4 and EC windows are gated by the tie range."),
        "",
        md_table(view),
        "",
        "Complexity (medians over stationary windows):",
        "",
        md_table(cx),
        "",
        "## Onsager regression test (pH, 0.01 d grid)",
        "",
        (f"R(τ): deviation of the tail ({ONSAGER_TAIL_START_D} d after the pulse) from the late level, normalized at τ = 0; "
        f"C(τ): normalized autocorrelation of the detrended late {ONSAGER_LATE_FRACTION:.0%} of the segment, which must pass "
        "ADF/KPSS (otherwise not computed). The SNR is the "
        "pulse response measured in M2 divided by the fluctuation SD; the tail deviation also contains the segment's "
        f"slow drift, which cannot be attributed to the pulse. Required SNR ≥ {ONSAGER_MIN_SNR:g}. Ratios with "
        f"moving-block bootstrap intervals are in `dynamics_onsager_lags.csv` (only where C ≥ {ONSAGER_MIN_C} and the late "
        "window is stationary)."),
        "",
        md_table(r.onsager_summary),
        "",
        onsager_conclusion(r.onsager_summary),
        "",
        "## Stability indicators before each pulse (EXPLORATORY)",
        "",
        md_table(r.stability.drop(columns=["note"])),
        "",
        (f"Kendall τ of the rolling ({STABILITY_WINDOW * 0.01:g} d) variance and lag-1 autocorrelation against time over up "
        f"to {STABILITY_SPAN_D:g} d before each pulse, on the linearly detrended grid series. ORP excludes "
        + ", ".join(f"{a:.2f}-{b:.2f} d" for a, b in orp_event_cuts())
        + " (electrode events); rolling windows never span an exclusion and each piece is detrended separately. "
        f"Overlapping windows make the p-values anti-conservative; with {len(cfg.pulse_times_d)} pulses these are "
        "exploratory."),
        "",
        "## Probability flux (0.01 d grid)",
        "",
        md_table(r.flux),
        "",
        ("Circulation rate ½⟨(x − x̄)Δy − (y − ȳ)Δx⟩ in standardized, linearly detrended coordinates, within allowed "
        "intervals; each input window (both coordinates) must pass ADF/KPSS, non-stationary windows are dropped; "
        f"moving-block bootstrap ({FLUX_BLOCK} bins, {FLUX_BOOT} resamples) basic 95 % intervals. Computed for "
        f"{int(flux_ok.sum())} of {len(r.flux)} plane-segment combinations."),
        "",
        "## Electrode dynamics",
        "",
        md_table(r.electrode),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(cfg: Config, out_dir, validation: pd.DataFrame | None = None) -> DynamicsResult:
    """``validation`` may be passed to reuse a suite computed by :func:`validate_estimators` (tests)."""
    from .pipeline import out_dirs

    out_dir = Path(out_dir)
    inputs = load_inputs(out_dir)
    d = out_dirs(out_dir)
    val = validate_estimators() if validation is None else validation
    val.to_csv(d["tables"] / "estimator_validation.csv", index=False)
    r = compute(cfg, inputs, val)  # raises EstimatorValidationError before any plant-data result
    r.validation.to_csv(d["tables"] / "estimator_validation.csv", index=False)  # with the matched-null gate checks
    r.summary.to_csv(d["tables"] / "dynamics.csv", index=False)
    r.windows.to_csv(d["tables"] / "dynamics_windows.csv", index=False)
    r.onsager_lags.to_csv(d["tables"] / "dynamics_onsager_lags.csv", index=False)
    r.onsager_summary.to_csv(d["tables"] / "dynamics_onsager.csv", index=False)
    r.stability.to_csv(d["tables"] / "dynamics_stability.csv", index=False)
    r.flux.to_csv(d["tables"] / "dynamics_flux.csv", index=False)
    r.electrode.to_csv(d["tables"] / "dynamics_electrode.csv", index=False)
    figure(cfg, r, d["figures"] / "fig5_dynamics.png", d["figures"] / "fig5_dynamics.svg")
    write_report(cfg, r, d["reports"] / "dynamics.md")
    return r
