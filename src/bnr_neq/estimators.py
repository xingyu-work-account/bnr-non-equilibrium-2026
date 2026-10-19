"""Irreversibility, complexity and surrogate estimators (Phase 3).

Plain NumPy arrays in and out, no pandas, no SciPy: these functions are meant
to be ported to a Rust sensor-entropy library. Every random function takes an
explicit ``numpy.random.Generator``.
"""

from __future__ import annotations

import math

import numpy as np

# Number of weak orderings (Fubini numbers) of d elements: patterns with ties.
FUBINI = (1, 1, 3, 13, 75, 541, 4683)


# --- Moment asymmetry -----------------------------------------------------------
def moment_asymmetry(x: np.ndarray, lags: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """<(x(t+τ) − x(t))³> for each lag τ (in samples), raw and normalized by <(·)²>^{3/2}."""
    x = np.asarray(x, dtype=float)
    raw = np.full(len(lags), np.nan)
    norm = np.full(len(lags), np.nan)
    for i, tau in enumerate(lags):
        if tau < 1 or tau >= len(x):
            raise ValueError(f"lag {tau} outside 1 to {len(x) - 1}")
        d = x[tau:] - x[:-tau]
        raw[i] = np.mean(d**3)
        m2 = np.mean(d**2)
        norm[i] = raw[i] / m2**1.5 if m2 > 0 else np.nan
    return raw, norm


# --- Ordinal patterns with ties ---------------------------------------------------
def ordinal_codes(x: np.ndarray, d: int = 3, lag: int = 1) -> np.ndarray:
    """Code of the weak-order (tie-aware) pattern of each embedding vector.

    A pattern is the tuple of signs sign(x_k − x_j) over all position pairs
    j < k (−1, 0, +1), so equalities are their own symbols. Code = Σ (sign+1) 3^i.
    """
    x = np.asarray(x, dtype=float)
    n = len(x) - (d - 1) * lag
    if n < 1:
        raise ValueError("series too short for the embedding")
    emb = np.stack([x[k * lag: k * lag + n] for k in range(d)], axis=1)
    code = np.zeros(n, dtype=np.int64)
    i = 0
    for j in range(d):
        for k in range(j + 1, d):
            code += (np.sign(emb[:, k] - emb[:, j]).astype(np.int64) + 1) * 3**i
            i += 1
    return code


def tie_fraction(x: np.ndarray) -> float:
    """Fraction of consecutive samples that are equal."""
    x = np.asarray(x, dtype=float)
    return float(np.mean(x[1:] == x[:-1]))


def _kl(cf: np.ndarray, cr: np.ndarray, pseudocount: float) -> float:
    support = (cf + cr) > 0
    p = cf[support] + pseudocount
    q = cr[support] + pseudocount
    p = p / p.sum()
    q = q / q.sum()
    return float(np.sum(p * np.log(p / q)))


def ordinal_irreversibility(x: np.ndarray, d: int = 3, lag: int = 1, pseudocount: float = 0.5) -> float:
    """KL divergence between ordinal-pattern frequencies of x and of x reversed in time.

    Pseudocount (Jeffreys, 0.5) on the observed support avoids infinite terms.
    """
    nsym = 3 ** (d * (d - 1) // 2)
    cf = np.bincount(ordinal_codes(x, d, lag), minlength=nsym).astype(float)
    cr = np.bincount(ordinal_codes(np.asarray(x)[::-1], d, lag), minlength=nsym).astype(float)
    return _kl(cf, cr, pseudocount)


def permutation_entropy_ties(x: np.ndarray, d: int = 3, lag: int = 1) -> float:
    """Shannon entropy of tie-aware ordinal patterns, normalized by ln(number of weak orders)."""
    c = np.bincount(ordinal_codes(x, d, lag)).astype(float)
    p = c[c > 0] / c.sum()
    return float(-np.sum(p * np.log(p)) / math.log(FUBINI[d]))


# --- Horizontal visibility graph ------------------------------------------------------
def hvg_degrees(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """In- and out-degrees of the directed horizontal visibility graph (Luque et al.; Lacasa et al. 2012).

    i -> j (i < j) if x_k < min(x_i, x_j) for every i < k < j; an equal intermediate value blocks.
    Out-degree counts links to the future, in-degree links from the past.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    k_in = np.zeros(n, dtype=np.int64)
    k_out = np.zeros(n, dtype=np.int64)
    stack: list[int] = []
    for i in range(n):
        xi = x[i]
        while stack and x[stack[-1]] < xi:
            j = stack.pop()
            k_out[j] += 1
            k_in[i] += 1
        if stack and x[stack[-1]] == xi:
            # an equal value is visible but blocks every taller value behind it
            j = stack.pop()
            k_out[j] += 1
            k_in[i] += 1
        elif stack:
            k_out[stack[-1]] += 1
            k_in[i] += 1
        stack.append(i)
    return k_in, k_out


def hvg_irreversibility(x: np.ndarray, pseudocount: float = 0.5) -> float:
    """KL divergence between the out-degree and in-degree distributions of the HVG."""
    k_in, k_out = hvg_degrees(x)
    m = int(max(k_in.max(), k_out.max())) + 1
    return _kl(np.bincount(k_out, minlength=m).astype(float), np.bincount(k_in, minlength=m).astype(float),
               pseudocount)


# --- Surrogates -----------------------------------------------------------------------
def iaaft(x: np.ndarray, rng: np.random.Generator, max_iter: int = 100) -> np.ndarray:
    """Iterated amplitude-adjusted Fourier transform surrogate (Schreiber and Schmitz 2000).

    Keeps the exact value distribution and approximately the power spectrum; it is
    time-reversible by construction. Iterates until the rank order stops changing.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    sorted_x = np.sort(x)
    amp = np.abs(np.fft.rfft(x))
    s = rng.permutation(x)
    ranks = np.argsort(np.argsort(s, kind="stable"), kind="stable")
    for _ in range(max_iter):
        phase = np.angle(np.fft.rfft(s))
        s = np.fft.irfft(amp * np.exp(1j * phase), n)
        new_ranks = np.argsort(np.argsort(s, kind="stable"), kind="stable")
        s = sorted_x[new_ranks]
        if np.array_equal(new_ranks, ranks):
            break
        ranks = new_ranks
    return s


def surrogate_test(x: np.ndarray, index_fn, rng: np.random.Generator, n_surr: int = 19,
                   two_sided: bool = False) -> tuple[float, np.ndarray, float]:
    """Index of x, indices of IAAFT surrogates, and the rank p-value (1 + #{surr ≥ obs})/(1 + n)."""
    obs = float(index_fn(x))
    surr = np.array([index_fn(iaaft(x, rng)) for _ in range(n_surr)])
    if two_sided:
        exceed = np.sum(np.abs(surr) >= abs(obs))
    else:
        exceed = np.sum(surr >= obs)
    return obs, surr, float((1 + exceed) / (1 + n_surr))


def surrogate_indices(x: np.ndarray, index_fns: dict, rng: np.random.Generator, n_surr: int = 19,
                      two_sided: frozenset = frozenset()) -> dict:
    """Evaluate several indices on x and on one shared set of IAAFT surrogates.

    Returns {name: (obs, surrogate_values, rank_p)}; names in ``two_sided`` use |index|.
    """
    surr = [iaaft(x, rng) for _ in range(n_surr)]
    out = {}
    for name, fn in index_fns.items():
        obs = float(fn(x))
        vals = np.array([fn(s) for s in surr])
        exceed = np.sum(np.abs(vals) >= abs(obs)) if name in two_sided else np.sum(vals >= obs)
        out[name] = (obs, vals, float((1 + exceed) / (1 + n_surr)))
    return out


# --- Complexity -------------------------------------------------------------------------
def _match_count(templates: np.ndarray, r: float, chunk: int = 512) -> int:
    """Number of pairs i < j with Chebyshev distance <= r."""
    n = len(templates)
    total = 0
    for a in range(0, n, chunk):
        block = templates[a: a + chunk]
        dist = np.max(np.abs(block[:, None, :] - templates[None, :, :]), axis=2)
        rows = np.arange(a, a + len(block))[:, None]
        total += int(np.sum((dist <= r) & (np.arange(n)[None, :] > rows)))
    return total


def sample_entropy(x: np.ndarray, m: int, r: float) -> float:
    """SampEn(m, r) = −ln(A/B) (Richman and Moorman 2000), Chebyshev distance, r absolute.

    NaN when no template pairs match (undefined).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if r <= 0 or n <= m + 1:
        raise ValueError("need r > 0 and len(x) > m + 1")
    tm = np.stack([x[k: n - m + k] for k in range(m)], axis=1)
    tm1 = np.stack([x[k: n - m + k] for k in range(m + 1)], axis=1)
    b = _match_count(tm, r)
    a = _match_count(tm1, r)
    if a == 0 or b == 0:
        return float("nan")
    return float(-math.log(a / b))


def coarse_grain(x: np.ndarray, scale: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = len(x) // scale
    return x[: n * scale].reshape(n, scale).mean(axis=1)


def multiscale_entropy(x: np.ndarray, scales, m: int, r: float) -> np.ndarray:
    """SampEn of coarse-grained series (Costa et al. 2002), r fixed from the original series."""
    return np.array([sample_entropy(coarse_grain(x, s), m, r) for s in scales])


def dfa(x: np.ndarray, scales: np.ndarray, order: int = 1) -> tuple[float, np.ndarray]:
    """Detrended fluctuation analysis: exponent α and fluctuation F(s) for each scale (samples)."""
    x = np.asarray(x, dtype=float)
    y = np.cumsum(x - x.mean())
    f = np.empty(len(scales))
    for i, s in enumerate(scales):
        nseg = len(y) // s
        if nseg < 2:
            raise ValueError(f"scale {s} too large for {len(y)} samples")
        segs = y[: nseg * s].reshape(nseg, s)
        t = np.arange(s, dtype=float)
        v = np.vander(t, order + 1)
        coef, *_ = np.linalg.lstsq(v, segs.T, rcond=None)
        resid = segs.T - v @ coef
        f[i] = math.sqrt(np.mean(resid**2))
    alpha = float(np.polyfit(np.log(scales), np.log(f), 1)[0])
    return alpha, f


# --- Markov chains --------------------------------------------------------------------------
def simulate_ring(p: float, q: float, n: int, rng: np.random.Generator, n_states: int = 3) -> np.ndarray:
    """Ring Markov chain: +1 with probability p, −1 with q, stay with 1 − p − q."""
    if p < 0 or q < 0 or p + q > 1:
        raise ValueError("need p, q >= 0 and p + q <= 1")
    u = rng.random(n - 1)
    step = np.where(u < p, 1, np.where(u < p + q, -1, 0))
    return (rng.integers(n_states) + np.concatenate([[0], np.cumsum(step)])) % n_states


def ring_entropy_production(p: float, q: float) -> float:
    """Analytic entropy production per step of the ring chain: (p − q) ln(p/q)."""
    if p == q:
        return 0.0
    return float((p - q) * math.log(p / q))


def transition_counts(pairs: np.ndarray, n_states: int) -> np.ndarray:
    """Counts N_ij from an (n, 2) array of consecutive state pairs."""
    return np.bincount(pairs[:, 0] * n_states + pairs[:, 1], minlength=n_states * n_states).reshape(n_states, n_states)


def entropy_production_from_counts(c: np.ndarray) -> tuple[float, int]:
    """Plug-in entropy production per step Σ_{i≠j} (N_ij/N) ln(N_ij/N_ji).

    Unordered pairs observed in one direction only are excluded and counted.
    """
    n = c.sum()
    total, excluded = 0.0, 0
    k = len(c)
    for i in range(k):
        for j in range(i + 1, k):
            a, b = c[i, j], c[j, i]
            if a == 0 and b == 0:
                continue
            if a == 0 or b == 0:
                excluded += 1  # unordered pair observed in one direction only
                continue
            total += (a - b) / n * math.log(a / b)
    return float(total), excluded


def symmetry_g_test(c: np.ndarray) -> tuple[float, int]:
    """Bowker-type G statistic for N_ij = N_ji (detailed balance of transition counts) and its dof.

    Conservative for a stationary chain, where flux conservation removes degrees of freedom.
    """
    g, dof = 0.0, 0
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            s = c[i, j] + c[j, i]
            if s == 0:
                continue
            dof += 1
            for a in (c[i, j], c[j, i]):
                if a > 0:
                    g += 2 * a * math.log(2 * a / s)
    return float(g), dof


def block_bootstrap_rows(n: int, block: int, rng: np.random.Generator) -> np.ndarray:
    """Row indices of one moving-block bootstrap resample of length n."""
    if block < 1 or block > n:
        raise ValueError("block must be between 1 and n")
    starts = rng.integers(0, n - block + 1, size=-(-n // block))
    return (starts[:, None] + np.arange(block)[None, :]).ravel()[:n]


def basic_bootstrap_ci(estimate: float, boot: np.ndarray, level: float) -> tuple[float, float]:
    """Basic (reverse-percentile) interval [2θ − q_hi, 2θ − q_lo]; corrects the plug-in bias."""
    a = (1 - level) / 2
    lo, hi = np.quantile(boot, [a, 1 - a])
    return float(2 * estimate - hi), float(2 * estimate - lo)


def markov_ep_with_ci(states: np.ndarray, n_states: int, block: int, n_boot: int, level: float,
                      rng: np.random.Generator) -> dict:
    """Transition-level entropy production per step with a moving-block bootstrap basic interval."""
    s = np.asarray(states, dtype=np.int64)
    pairs = np.stack([s[:-1], s[1:]], axis=1)
    c = transition_counts(pairs, n_states)
    est, excluded = entropy_production_from_counts(c)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        boot[b] = entropy_production_from_counts(transition_counts(pairs[block_bootstrap_rows(len(pairs), block, rng)],
                                                                   n_states))[0]
    lo, hi = basic_bootstrap_ci(est, boot, level)
    g, dof = symmetry_g_test(c)
    return {"estimate": est, "ci_low": lo, "ci_high": hi, "excluded_pairs": excluded, "g_stat": g, "g_dof": dof}


# --- Probability flux (Battle et al. 2016) ----------------------------------------------------
def circulation_rate(x: np.ndarray, y: np.ndarray, valid: np.ndarray | None = None) -> float:
    """Mean signed area rate ½ <(x_t − x̄) Δy_t − (y_t − ȳ) Δx_t> over valid increments.

    Zero in expectation under detailed balance; its sign gives the direction of the
    net probability circulation. ``valid[t]`` marks increments t -> t+1 to use.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx, dy = np.diff(x), np.diff(y)
    if valid is None:
        valid = np.ones(len(dx), dtype=bool)
    xc, yc = x[:-1] - x.mean(), y[:-1] - y.mean()
    a = 0.5 * (xc * dy - yc * dx)
    return float(a[valid].mean())


def circulation_with_ci(x: np.ndarray, y: np.ndarray, valid: np.ndarray, block: int, n_boot: int,
                        level: float, rng: np.random.Generator) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx, dy = np.diff(x), np.diff(y)
    a = 0.5 * ((x[:-1] - x.mean()) * dy - (y[:-1] - y.mean()) * dx)
    a = a[valid]
    est = float(a.mean())
    boot = np.array([a[block_bootstrap_rows(len(a), block, rng)].mean() for _ in range(n_boot)])
    lo, hi = basic_bootstrap_ci(est, boot, level)
    return {"estimate": est, "ci_low": lo, "ci_high": hi, "n_increments": len(a)}


def flux_field(x: np.ndarray, y: np.ndarray, edges_x: np.ndarray, edges_y: np.ndarray,
               valid: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Battle et al. flux estimate per bin: ½ (mean displacement leaving + entering the bin), per step."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx, dy = np.diff(x), np.diff(y)
    if valid is None:
        valid = np.ones(len(dx), dtype=bool)
    nx, ny = len(edges_x) - 1, len(edges_y) - 1
    ix = np.clip(np.searchsorted(edges_x, x, side="right") - 1, 0, nx - 1)
    iy = np.clip(np.searchsorted(edges_y, y, side="right") - 1, 0, ny - 1)
    jx = np.zeros((nx, ny))
    jy = np.zeros((nx, ny))
    count = np.zeros((nx, ny))
    for src in (slice(None, -1), slice(1, None)):
        bx, by = ix[src][valid], iy[src][valid]
        np.add.at(jx, (bx, by), dx[valid] / 2)
        np.add.at(jy, (bx, by), dy[valid] / 2)
        np.add.at(count, (bx, by), 0.5)
    n = valid.sum()
    return jx / n, jy / n, count
