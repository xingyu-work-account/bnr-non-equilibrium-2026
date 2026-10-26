import numpy as np
import pandas as pd
import pytest

from bnr_neq import qc

RNG = np.random.default_rng(20261018)


def _series(days=2.0, dt_s=20.0, noise=0.05):
    t = np.arange(0, days, dt_s / 86400)
    return t, 40.0 + noise * RNG.standard_normal(len(t))


def test_rolling_step_nan_near_edges():
    t, x = _series()
    s = qc.rolling_step(t, x, 0.05)
    assert np.isnan(s[t < 0.05]).all() and np.isnan(s[t > t[-1] - 0.05]).all()
    assert np.isfinite(s[(t > 0.06) & (t < t[-1] - 0.06)]).all()


def test_instant_step_located_and_classified_as_offset():
    t, x = _series()
    x[t >= 1.0] -= 20.0
    j = qc.jump_table(t, x, pulse_times_d=[0.3])
    j = j[j["jump"]]
    assert len(j) == 1
    row = j.iloc[0]
    assert row["t_jump_d"] == pytest.approx(1.0, abs=25 / 86400)  # within about one sample
    assert row["level_change"] == pytest.approx(-20.0, abs=0.1)
    assert row["transition_10_90_min"] <= 20 / 60 + 1e-6
    assert row["class"] == "electrode_offset"


def test_delayed_rise_after_pulse_is_pulse_response():
    t, x = _series()
    x[t >= 1.3] += 30.0
    j = qc.jump_table(t, x, pulse_times_d=[1.0])
    j = j[j["jump"]]
    assert j["class"].tolist() == ["pulse_response"]
    assert j["pulse"].iloc[0] == 1
    assert j["lag_h"].iloc[0] == pytest.approx(7.2, abs=0.05)


def test_rise_outside_lag_window_is_offset():
    t, x = _series(days=3.0)
    x[t >= 1.8] += 10.0  # 0.8 d after the pulse
    j = qc.jump_table(t, x, pulse_times_d=[1.0])
    assert j.loc[j["jump"], "class"].tolist() == ["electrode_offset"]


def test_transient_excursion_is_not_a_sustained_jump():
    t, x = _series()
    x[(t >= 1.0) & (t < 1.06)] += 10.0  # 86-min excursion that returns
    j = qc.jump_table(t, x, pulse_times_d=[])
    assert len(j) >= 1 and not j["jump"].any()


def test_slow_ramp_transition_time():
    t, x = _series(noise=0.01)
    ramp = (t >= 1.0) & (t < 1.0 + 1 / 24)  # 60-min linear rise
    x[ramp] += 10.0 * (t[ramp] - 1.0) * 24
    x[t >= 1.0 + 1 / 24] += 10.0
    ch = qc.characterize_jump(t, x, 1.0 + 0.5 / 24)
    assert ch["transition_10_90_min"] == pytest.approx(48.0, abs=3.0)  # 80 % of 60 min
    assert ch["t50_d"] == pytest.approx(1.0 + 0.5 / 24, abs=2 / 1440)


def test_held_segments_found_only_where_output_is_smooth():
    t = np.arange(0, 3.0, 20 / 86400)
    x = 40.0 + 0.5 * RNG.standard_normal(len(t))
    held = (t >= 1.0) & (t < 2.0)
    x[held] = np.round(40.0 - 0.5 * (t[held] - 1.0), 2)  # smooth, quantized to 0.01
    blocks = qc.hf_blocks(t, qc.hf_residual(t, x), 0.0, 3.0)
    seg = qc.held_segments(blocks)
    assert len(seg) == 1
    assert seg["t_start_d"].iloc[0] == pytest.approx(1.0, abs=0.11)
    assert seg["t_end_d"].iloc[0] == pytest.approx(2.0, abs=0.11)


def test_ec_pulses_recovered_with_decay_time():
    t = np.arange(0, 3.0, 60 / 86400)
    ec = 1.0 + 0.002 * RNG.standard_normal(len(t))
    peaks = [2.2, 2.4, 2.6, 2.8]
    for tp in peaks:
        m = t >= tp
        ec[m] += 0.09 * np.exp(-(t[m] - tp) * 1440 / 24.0)
    resid = qc.ec_residual(t, ec)
    found = qc.detect_ec_pulses(t, ec, resid, 0.03, start_d=2.0)
    assert len(found) == 4
    np.testing.assert_allclose(found["t_peak_d"], peaks, atol=2 / 1440)
    assert found["fit_ok"].all()
    assert found["tau_min"].median() == pytest.approx(24.0, rel=0.15)


def test_ec_threshold_rule_is_the_pre_period_quantile():
    t = np.arange(0, 3.0, 60 / 86400)
    resid = RNG.standard_normal(len(t))
    resid[t >= 2.0] += 100.0  # the post period must not influence the threshold
    expected = np.quantile(resid[t < 2.0], qc.EC_THRESHOLD_QUANTILE)
    assert qc.ec_pulse_threshold(t, resid, start_d=2.0) == pytest.approx(expected)


def test_ec_sensitivity_counts_decrease_with_threshold():
    t = np.arange(0, 3.0, 60 / 86400)
    ec = 1.0 + 0.002 * RNG.standard_normal(len(t))
    for tp, a in [(2.2, 0.03), (2.5, 0.09)]:
        m = t >= tp
        ec[m] += a * np.exp(-(t[m] - tp) * 1440 / 24.0)
    resid = qc.ec_residual(t, ec)
    sens = qc.ec_pulse_sensitivity(t, resid, [0.02, 0.05], start_d=2.0)
    assert sens["n_pulses"].tolist() == [2, 1]


def test_robust_sigma_of_constant_is_zero():
    assert qc.robust_sigma(np.zeros(10)) == 0.0
    assert np.isnan(qc.robust_sigma(np.array([np.nan])))


def test_step_time_in_window():
    t, x = _series()
    x[t >= 1.2] += 1.0
    ts, step = qc.step_time_in_window(t, x, 1.1, 1.3)
    assert step == pytest.approx(1.0, abs=0.05)
    assert abs(ts - 1.2) <= 0.006  # median-difference plateau is +-w/2 wide


def test_jump_table_columns_are_stable():
    t, x = _series()
    assert isinstance(qc.jump_table(t, x, []), pd.DataFrame)
