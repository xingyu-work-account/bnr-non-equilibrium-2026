import itertools
import math

import numpy as np
import pytest

from bnr_neq import estimators as E

RNG = np.random.default_rng(99)


def brute_hvg(x):
    n = len(x)
    k_in, k_out = np.zeros(n, int), np.zeros(n, int)
    for i, j in itertools.combinations(range(n), 2):
        if all(x[k] < min(x[i], x[j]) for k in range(i + 1, j)):
            k_out[i] += 1
            k_in[j] += 1
    return k_in, k_out


def test_hvg_small_examples():
    k_in, k_out = E.hvg_degrees(np.array([1.0, 3.0, 2.0, 4.0]))
    assert k_out.tolist() == [1, 2, 1, 0] and k_in.tolist() == [0, 1, 1, 2]
    k_in, k_out = E.hvg_degrees(np.array([2.0, 1.0, 2.0, 2.0]))  # equal intermediate value blocks
    assert k_out.tolist() == [2, 1, 1, 0] and k_in.tolist() == [0, 1, 2, 1]


@pytest.mark.parametrize("seed", range(5))
def test_hvg_matches_brute_force_with_ties(seed):
    x = np.random.default_rng(seed).integers(0, 4, size=60).astype(float)
    k_in, k_out = E.hvg_degrees(x)
    b_in, b_out = brute_hvg(x)
    assert k_in.tolist() == b_in.tolist() and k_out.tolist() == b_out.tolist()


def test_ordinal_codes_distinguish_ties_and_count_weak_orders():
    codes = E.ordinal_codes(RNG.integers(0, 3, size=5000).astype(float), d=3)
    assert len(np.unique(codes)) == 13  # Fubini number for d = 3
    assert E.ordinal_codes(np.array([1.0, 1.0, 2.0]))[0] != E.ordinal_codes(np.array([1.0, 2.0, 2.0]))[0]


def test_ordinal_irreversibility_properties():
    assert E.ordinal_irreversibility(np.ones(100)) == 0.0
    saw = np.tile(np.r_[np.linspace(0, 1, 9), 0.0], 50)  # slow rise, abrupt fall
    assert E.ordinal_irreversibility(saw) > 0.5
    assert E.ordinal_irreversibility(saw[::-1]) > 0.5  # abrupt rise, slow fall: just as irreversible


def test_ordinal_irreversibility_time_reversal_property():
    """Reversing time maps each embedding window to its mirror image, so the forward pattern counts of x[::-1]
    are the counts of x's mirrored windows; the index compares exactly these two distributions. Reversal swaps
    the arguments of the divergence, and since mirroring is an involution on the symbols, D(p||πp) = D(πp||p):
    the index of x[::-1] equals that of x (a non-trivial invariance), and a palindrome scores exactly zero."""
    x = RNG.integers(0, 4, size=3000).astype(float) + np.tile(np.r_[np.linspace(0, 2, 7), 0.0], 375)
    d = 3
    mirrored = np.array([E.ordinal_codes(x[i: i + d][::-1], d)[0] for i in range(len(x) - d + 1)])[::-1]
    np.testing.assert_array_equal(E.ordinal_codes(x[::-1], d), mirrored)
    nsym = 3 ** (d * (d - 1) // 2)
    cf = np.bincount(E.ordinal_codes(x, d), minlength=nsym).astype(float)
    cr = np.bincount(E.ordinal_codes(x[::-1], d), minlength=nsym).astype(float)
    assert E.ordinal_irreversibility(x) == pytest.approx(E._kl(cf, cr, 0.5))
    assert E.ordinal_irreversibility(x[::-1]) == pytest.approx(E._kl(cr, cf, 0.5))
    assert E.ordinal_irreversibility(x) > 0.01
    assert E.ordinal_irreversibility(x[::-1]) == pytest.approx(E.ordinal_irreversibility(x), rel=1e-12)
    perm = np.random.default_rng(0).permutation(nsym)  # a non-involutive relabelling breaks the equality
    assert E._kl(cf, cf[perm], 0.5) != pytest.approx(E._kl(cf[perm], cf, 0.5))
    pal = np.r_[x[:1500], x[:1500][::-1]]
    assert E.ordinal_irreversibility(pal) == 0.0


def test_tie_fraction_and_permutation_entropy():
    assert E.tie_fraction(np.array([1.0, 1.0, 2.0, 2.0, 3.0])) == 0.5
    assert E.permutation_entropy_ties(RNG.random(20000)) == pytest.approx(math.log(6) / math.log(13), abs=0.01)
    assert E.permutation_entropy_ties(np.ones(50)) == 0.0


def test_moment_asymmetry():
    saw = np.tile(np.r_[np.linspace(0, 1, 9), 0.0], 50)
    raw, norm = E.moment_asymmetry(saw, [1, 2])
    assert (raw < 0).all() and (norm < 0).all()  # abrupt falls dominate the third moment
    with pytest.raises(ValueError):
        E.moment_asymmetry(saw, [0])


def test_iaaft_preserves_values_and_spectrum():
    x = np.cumsum(RNG.normal(size=1024)) * 0.1
    s = E.iaaft(x, np.random.default_rng(1))
    np.testing.assert_array_equal(np.sort(s), np.sort(x))
    a, b = np.abs(np.fft.rfft(x))[1:], np.abs(np.fft.rfft(s))[1:]
    assert np.corrcoef(np.log(a), np.log(b))[0, 1] > 0.9


def test_surrogate_indices_p_values():
    x = RNG.normal(size=500)
    out = E.surrogate_indices(x, {"o": E.ordinal_irreversibility}, np.random.default_rng(2), n_surr=19)
    _obs, vals, p = out["o"]
    assert len(vals) == 19 and 1 / 20 <= p <= 1.0


def test_sample_entropy():
    assert E.sample_entropy(np.tile([1.0, 2.0, 3.0], 200), 2, 0.1) == pytest.approx(0.0, abs=1e-9)
    assert E.sample_entropy(RNG.normal(size=2000), 2, 0.2) == pytest.approx(2.2, abs=0.2)
    with pytest.raises(ValueError):
        E.sample_entropy(np.ones(10), 2, 0.0)


def test_multiscale_entropy_and_dfa():
    x = RNG.normal(size=4000)
    mse = E.multiscale_entropy(x, [1, 2], 2, 0.2)
    assert mse[1] < mse[0]  # white noise: coarse-graining lowers the variance at fixed r (Costa et al.)
    sc = np.unique(np.logspace(1, np.log10(500), 10).astype(int))
    assert E.dfa(x, sc)[0] == pytest.approx(0.5, abs=0.1)
    with pytest.raises(ValueError):
        E.dfa(x[:30], np.array([20]))


def test_ring_chain_and_entropy_production():
    s = E.simulate_ring(0.5, 0.2, 50_000, np.random.default_rng(3))
    steps = (np.diff(s) % 3)
    assert np.mean(steps == 1) == pytest.approx(0.5, abs=0.01) and np.mean(steps == 2) == pytest.approx(0.2, abs=0.01)
    assert E.ring_entropy_production(0.3, 0.3) == 0.0
    assert E.ring_entropy_production(0.6, 0.1) == pytest.approx(0.5 * math.log(6))
    with pytest.raises(ValueError):
        E.simulate_ring(0.7, 0.5, 10, RNG)


def test_entropy_production_from_counts():
    sym = np.array([[5, 3, 2], [3, 5, 4], [2, 4, 5]])
    assert E.entropy_production_from_counts(sym) == (0.0, 0)
    one_way = np.array([[0, 4, 0], [0, 0, 4], [4, 0, 0]])
    est, excluded = E.entropy_production_from_counts(one_way)
    assert est == 0.0 and excluded == 3
    g, dof = E.symmetry_g_test(sym)
    assert g == pytest.approx(0.0) and dof == 3


def test_markov_ep_ci_covers_truth():
    rng = np.random.default_rng(4)
    r = E.markov_ep_with_ci(E.simulate_ring(0.4, 0.2, 100_000, rng), 3, 100, 200, 0.99, rng)
    assert r["ci_low"] <= E.ring_entropy_production(0.4, 0.2) <= r["ci_high"]


def test_bootstrap_helpers():
    idx = E.block_bootstrap_rows(100, 10, np.random.default_rng(5))
    assert len(idx) == 100 and idx.min() >= 0 and idx.max() < 100
    with pytest.raises(ValueError):
        E.block_bootstrap_rows(10, 20, RNG)
    lo, hi = E.basic_bootstrap_ci(1.0, np.array([1.0, 1.2, 1.4]), 0.5)
    assert lo < hi and lo == pytest.approx(2 - np.quantile([1.0, 1.2, 1.4], 0.75))


def test_circulation_and_flux_field_sign():
    th = np.linspace(0, 20 * np.pi, 4000)
    x, y = np.cos(th), np.sin(th)  # counter-clockwise
    assert E.circulation_rate(x, y) > 0
    assert E.circulation_rate(x[::-1], y[::-1]) < 0
    edges = np.linspace(-1.1, 1.1, 5)
    jx, jy, count = E.flux_field(x, y, edges, edges)
    assert jx.shape == (4, 4) and count.sum() == pytest.approx(len(x) - 1)
    assert jy[3, 1:3].mean() > 0 and jy[0, 1:3].mean() < 0  # up on the right, down on the left


def test_estimators_module_has_no_pandas_or_scipy():
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(E))
    names = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))
             for a in n.names} | {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    assert "pandas" not in names and "scipy" not in names
