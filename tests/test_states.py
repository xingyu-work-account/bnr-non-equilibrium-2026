"""W3: process states from relative ORP, pH and the inverted NH4, and the block-bootstrap comparison."""

import numpy as np
import pandas as pd
import pytest

from bnr_neq import states as ST

VERDICTS = pd.DataFrame({"key": ["nh4", "ec", "orp", "ph"],
                         "verdict": ["FAIL: under-response", "FAIL: no response", "FAIL: no chloride response",
                                     "PASS (no contradiction)"]})


# --- negative: a rule that references an excluded sensor raises ---------------------------------------------
def test_rule_referencing_an_excluded_sensor_raises():
    excluded = ST.excluded_features(VERDICTS)
    assert {"ec", "orp_absolute", "nh4_reading"} <= excluded and "ph" not in excluded
    bad = (*ST.RULES, ST.Rule("active_nitrification", "conductivity falls", ("ec",)))
    with pytest.raises(ST.StateRuleError, match="ec"):
        ST.validate_rules(bad, excluded)
    with pytest.raises(ST.StateRuleError, match="orp_absolute"):
        ST.validate_rules([ST.Rule("oxygen_limited", "E_H < 0", ("orp_absolute",))], excluded)
    ST.validate_rules(ST.RULES, excluded)  # the registered rules use none of them


def test_a_failed_ph_sensor_excludes_every_registered_rule():
    v = VERDICTS.assign(verdict=lambda d: np.where(d["key"] == "ph", "FAIL: drift", d["verdict"]))
    with pytest.raises(ST.StateRuleError, match="ph"):
        ST.validate_rules(ST.RULES, ST.excluded_features(v))


# --- rules on synthetic features -------------------------------------------------------------------------------
def _features(nh4, ph, orp, level, artifact=""):
    n = len(nh4)
    return pd.DataFrame({"t_mid_d": 1.005 + 0.01 * np.arange(n), "artifact_reason": [artifact] * n,
                         "nh4_sign": nh4, "ph_sign": ph, "orp_sign": orp, "orp_rel_window_median_mV": level,
                         "C_inv_mgN_L": 100.0, "pH": 7.5, "nh4_slope_mgN_L_d": -1.0})


def test_registered_rules_on_synthetic_bins():
    f = _features(nh4=[-1, -1, 0, 0, 0, 1, 0, np.nan], ph=[0, 1, 0, 1, 0, 0, 0, 0], orp=[0, 0, -1, 0, 0, 0, 0, 0],
                  level=[0, 0, -20, 0, 20, 0, 0, 0])
    st = ST.classify(f)["state"].tolist()
    assert st == ["active_nitrification", "unclassified", "oxygen_limited", "inactive_oversupplied",
                  "inactive_oversupplied", "unclassified", "unclassified", "unclassified"]


def test_orp_level_threshold_is_the_nernstian_decade():
    assert ST.ORP_LEVEL_MV == 15  # RT ln10 / 4F at 25 °C = 14.8 mV, rounded


# --- boundary: an all-artifact interval yields no state ------------------------------------------------------
def test_all_artifact_interval_yields_no_state():
    f = _features(nh4=[-1] * 50, ph=[0] * 50, orp=[0] * 50, level=[0] * 50, artifact="NH4 held/damped output")
    bins = ST.classify(f)
    assert (bins["state"] == ST.ARTIFACT).all() and not bins["state"].isin(ST.STATES).any()
    fr = ST.state_fractions(bins, 0.0, 10.0)
    # pass 3 (X2): the idle fraction over classified time is not computed; the bounds are still reported and are the
    # uninformative [0, 1] of an interval without any process state
    assert fr["status"].startswith("not computed") and np.isnan(fr["idle_fraction_of_classified_time"])
    assert (fr["idle_fraction_lower_bound"], fr["idle_fraction_upper_bound"], fr["classified_fraction"]) == (0.0, 1.0, 0.0)
    assert fr["fraction_sensor_artifact"] == 1.0


def test_sign_needs_the_interval_to_exclude_zero_and_the_quantum():
    s = ST.sign_class(np.array([0.02, 0.02, -0.5, 0.3]), np.array([0.01, -0.01, -0.6, 0.1]),
                      np.array([0.03, 0.05, -0.4, 0.5]), min_abs=0.05)
    assert s.tolist() == [0.0, 0.0, -1.0, 1.0]


def test_rolling_slope_skips_windows_that_touch_an_artifact():
    t = 0.01 * np.arange(300)
    clean = np.ones(300, dtype=bool)
    clean[150] = False
    out = ST.rolling_theil_sen(t, 2.0 * t, clean, half=10)
    assert np.isnan(out[0, 140:161]).all() and np.allclose(out[0, 20:139], 2.0) and np.isnan(out[0, :10]).all()


# --- block bootstrap -------------------------------------------------------------------------------------------
def _ar1(n, phi, seed):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0]
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def test_block_bootstrap_recovers_a_shift():
    y = _ar1(2000, 0.5, 1)
    a = np.zeros(2000, dtype=bool)
    a[::2] = True
    r = ST.block_bootstrap_diff(y + 3.0 * a, a, ~a, block=20, n_boot=299)
    assert r["status"] == "computed" and r["ci95_low"] < 3.0 < r["ci95_high"] and r["ci95_low"] > 2.0


def test_blocks_widen_the_interval_of_an_autocorrelated_series():
    y = _ar1(3000, 0.98, 2)
    a = np.zeros(3000, dtype=bool)
    a[:1500] = True
    iid = ST.block_bootstrap_diff(y, a, ~a, block=1, n_boot=299)
    blk = ST.block_bootstrap_diff(y, a, ~a, block=ST.acf_block_length(y, 1), n_boot=299)
    assert ST.acf_block_length(y, 1) > 20
    assert (blk["ci95_high"] - blk["ci95_low"]) > 3 * (iid["ci95_high"] - iid["ci95_low"])


def test_block_bootstrap_is_deterministic_and_guards_small_groups():
    y = _ar1(500, 0.5, 3)
    a = np.zeros(500, dtype=bool)
    a[:100] = True
    r1 = ST.block_bootstrap_diff(y, a, ~a, block=10, n_boot=199)
    r2 = ST.block_bootstrap_diff(y, a, ~a, block=10, n_boot=199)
    assert r1 == r2
    small = ST.block_bootstrap_diff(y, a, ~a, block=60, n_boot=199)
    assert small["status"].startswith("not computed: a group has fewer than 2 blocks") and np.isnan(small["ci95_low"])
    with pytest.raises(ValueError, match="overlap"):
        ST.block_bootstrap_diff(y, a, a, block=10)


def test_within_stratum_comparison_removes_a_between_stratum_trend():
    rng = np.random.default_rng(4)
    strata = np.repeat([0, 1], 1000)
    a = np.tile([True, False], 1000)
    a[1000:] = rng.random(1000) < 0.9  # stratum 1 is mostly group A and sits 10 higher
    y = 10.0 * strata + rng.normal(size=2000)
    pooled = ST.block_bootstrap_diff(y, a, ~a, block=5, n_boot=199)
    within = ST.block_bootstrap_diff(y, a, ~a, block=5, n_boot=199, strata=strata)
    assert pooled["diff"] > 3 and abs(within["diff"]) < 0.3 and within["n_strata"] == 2


# --- real data (production stage on the run-all outputs) ---------------------------------------------------------
def test_states_on_the_real_record(cfg, run_all_dir):
    bins = pd.read_parquet(run_all_dir / "processed" / "states_bins.parquet")
    assert set(bins["state"]) <= {*ST.STATES, ST.ARTIFACT, ST.UNCLASSIFIED}
    art = bins["state"] == ST.ARTIFACT
    assert (bins.loc[art, "artifact_reason"] != "").all() and (bins.loc[~art, "artifact_reason"] == "").all()
    held = pd.read_csv(run_all_dir / "tables" / "nh4_held_segments.csv")
    for a, b in zip(held["t_start_d"], held["t_end_d"]):
        inside = (bins["t_mid_d"] >= a) & (bins["t_mid_d"] < b)
        assert not bins.loc[inside, "state"].isin(ST.STATES).any()  # held NH4 output never carries a process state
    fr = pd.read_csv(run_all_dir / "tables" / "states_fractions.csv").set_index("interval")
    cols = [f"fraction_{s}" for s in (*ST.STATES, ST.ARTIFACT, ST.UNCLASSIFIED)]
    assert np.allclose(fr[cols].sum(axis=1), 1.0)
    comp = pd.read_csv(run_all_dir / "tables" / "states_comparison.csv")
    assert comp["role"].str.startswith(("ASSOCIATION TEST (NH4-free", "DESCRIPTIVE (circular")).all()
    assert (comp["block_bins"] >= 2 * ST.SLOPE_HALF_BINS + 1).all()
    # the idle fraction in the pulse table is the W3 one (pass 3: as bounds with the coverage beside them)
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    idle = pm[pm["metric"] == "idle_fraction_LOWER_BOUND"].set_index("pulse")["value"]
    for k in (1, 2, 3, 4):
        assert idle[k] == pytest.approx(fr.loc[f"pulse {k}", "idle_fraction_lower_bound"])


# --- pass 3, X2: autocorrelation-aware slope intervals and the bounded idle fraction ----------------------------------
def _window_coverage(phi, n_windows, seed, rho=None):
    """Coverage of the true slope by the production interval at the centre of 101-point windows with AR(1) errors."""
    rng = np.random.default_rng(seed)
    n, slope = 2 * ST.SLOPE_HALF_BINS + 1, -2.0
    t = 0.01 * np.arange(n)
    scaled = unscaled = known = 0
    for _ in range(n_windows):
        y = slope * t + 0.05 * _ar1(n + 300, phi, int(rng.integers(1 << 31)))[300:]
        out = ST.rolling_theil_sen(t, y, np.ones(n, dtype=bool))
        c = ST.SLOPE_HALF_BINS
        scaled += out[1, c] <= slope <= out[2, c]
        r = ST.theilslopes(y, t, alpha=ST.CI_ALPHA)
        unscaled += r.low_slope <= slope <= r.high_slope
        k = ST.interval_scale(n, phi)
        known += r.slope - k * (r.slope - r.low_slope) <= slope <= r.slope + k * (r.high_slope - r.slope)
    return scaled / n_windows, unscaled / n_windows, known / n_windows


def test_autocorrelation_scaled_interval_covers_an_ar1_slope():
    scaled, unscaled, known = _window_coverage(0.9, 400, 20261020)
    # the √(n/n_eff) formula holds its 95 % with the true ρ; with ρ estimated from 101 detrended points (biased low,
    # mean ≈ 0.84 for φ = 0.9) the coverage is ≈ 0.85; the independent-point interval covers ≈ 0.35 (docs/progress.md)
    assert 0.92 <= known <= 0.98
    assert scaled >= 0.80 and unscaled <= 0.50 and scaled - unscaled >= 0.35


def test_rho_cap_keeps_n_eff_positive():
    n = 2 * ST.SLOPE_HALF_BINS + 1
    k_cap = ST.interval_scale(n, 1.0)  # ρ = 1 would give n_eff = 0
    assert np.isfinite(k_cap) and k_cap == pytest.approx(ST.interval_scale(n, ST.RHO_MAX))
    assert k_cap == pytest.approx(np.sqrt((1 + ST.RHO_MAX) / (1 - ST.RHO_MAX)))
    assert ST.interval_scale(n, -0.3) == 1.0 and ST.lag1_rho(np.ones(n)) == 0.0  # no variance: no widening
    assert ST.lag1_rho(np.arange(n, dtype=float)) == ST.RHO_MAX  # lag-1 autocorrelation 1, capped
    t = 0.01 * np.arange(n)
    out = ST.rolling_theil_sen(t, 2.0 * t + np.repeat([0.0, 1.0], [n // 2, n - n // 2]), np.ones(n, dtype=bool))
    assert np.isfinite(out[:, ST.SLOPE_HALF_BINS]).all() and out[3, ST.SLOPE_HALF_BINS] <= ST.RHO_MAX


def _state_bins(states):
    return pd.DataFrame({"t_mid_d": 1.005 + 0.01 * np.arange(len(states)), "state": states})


def test_idle_bounds_equal_hand_counts():
    st = (["inactive_oversupplied"] * 3 + ["active_nitrification"] * 4 + ["oxygen_limited"] * 1 + [ST.UNCLASSIFIED] * 1
          + [ST.ARTIFACT] * 1)  # 10 bins, 8 classified
    fr = ST.state_fractions(_state_bins(st), 0.0, 10.0)
    assert fr["n_bins"] == 10 and fr["n_classified"] == 8 and fr["n_inactive_oversupplied"] == 3
    assert fr["idle_fraction_lower_bound"] == pytest.approx(3 / 10)
    assert fr["idle_fraction_upper_bound"] == pytest.approx((3 + 1 + 1) / 10)
    assert fr["classified_fraction"] == pytest.approx(0.8)
    assert fr["idle_fraction_of_classified_time"] == pytest.approx(3 / 8)
    assert fr["status"].startswith("computed (CONDITIONAL")


def test_conditional_idle_fraction_needs_the_coverage_but_the_bounds_do_not():
    st = ["inactive_oversupplied"] * 2 + ["active_nitrification"] * 2 + [ST.UNCLASSIFIED] * 6  # coverage 0.4 < 0.5
    fr = ST.state_fractions(_state_bins(st), 0.0, 10.0)
    assert ST.MIN_CLASSIFIED_FRACTION == 0.5
    assert np.isnan(fr["idle_fraction_of_classified_time"])
    assert fr["status"] == ("not computed: classified fraction 0.400 below 0.5 (the classified bins are not "
                            "representative of the interval; see the bounds)")
    assert (fr["idle_fraction_lower_bound"], fr["idle_fraction_upper_bound"]) == pytest.approx((0.2, 0.8))
    assert fr["bounds_status"] == "computed"
    # pass 4, Y5: a classified fraction of exactly 0.5 (5 of 10) computes the conditional value
    at_threshold = ST.state_fractions(_state_bins(["inactive_oversupplied"] * 2 + ["active_nitrification"] * 3
                                                  + [ST.UNCLASSIFIED] * 5), 0.0, 10.0)
    assert at_threshold["classified_fraction"] == 0.5
    assert at_threshold["idle_fraction_of_classified_time"] == pytest.approx(2 / 5)
    assert at_threshold["status"].startswith("computed (CONDITIONAL")
    half = ST.state_fractions(_state_bins(["inactive_oversupplied"] * 5 + [ST.ARTIFACT] * 5), 0.0, 10.0)
    assert half["classified_fraction"] == 0.5 and half["idle_fraction_of_classified_time"] == 1.0  # ≥ 0.5 computes


def test_real_idle_metrics_carry_bounds_and_coverage(run_all_dir):
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    fr = pd.read_csv(run_all_dir / "tables" / "states_fractions.csv").set_index("interval")
    idle = pm[pm["metric"].str.startswith("idle_fraction")]
    assert set(idle["metric"]) == {"idle_fraction_LOWER_BOUND", "idle_fraction_UPPER_BOUND",
                                   "idle_fraction_of_classified_time_CONDITIONAL"}
    assert idle["basis"].str.contains("coverage: classified fraction").all()  # no idle fraction without its coverage
    m = pm.set_index(["pulse", "metric"])
    for k in (1, 2, 3, 4):
        cov = fr.loc[f"pulse {k}", "classified_fraction"]
        assert m.loc[(k, "classified_fraction"), "value"] == pytest.approx(cov)
        lo, hi = m.loc[(k, "idle_fraction_LOWER_BOUND"), "value"], m.loc[(k, "idle_fraction_UPPER_BOUND"), "value"]
        assert 0 <= lo <= hi <= 1
        cond = m.loc[(k, "idle_fraction_of_classified_time_CONDITIONAL")]
        if cov < ST.MIN_CLASSIFIED_FRACTION:  # below the coverage, or (coverage 0) no process state at all
            assert cond["status"].startswith(("not computed: classified fraction",
                                              "not computed: no bin with a process state")) and np.isnan(cond["value"])
        else:
            assert cond["status"] == "computed"
    bins = pd.read_parquet(run_all_dir / "processed" / "states_bins.parquet")
    have = bins["nh4_slope_mgN_L_d"].notna()
    assert (bins.loc[have, "nh4_slope_resid_lag1"].between(0, ST.RHO_MAX)).all()
