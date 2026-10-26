import numpy as np
import pandas as pd
import pytest
from conftest import synthetic_truth

from bnr_neq import pulses as P
from bnr_neq import segments as S
from bnr_neq import thermo
from bnr_neq.config import CONFIG_PATH, load_config


def test_free_energy_per_pulse(cfg):
    assert P.free_energy_per_pulse_J(cfg) == pytest.approx(163.0, abs=0.1)


def test_model_integral_constant_rate():
    t = np.round(np.arange(1.005, 5.0, 0.01), 10)
    state = pd.DataFrame({"t_d": t, "model_rate_mgN_L_d": 2.0})
    integral, mean = P.model_integral(state, 2.0, 3.0)
    assert integral == pytest.approx(2.0, abs=1e-9) and mean == pytest.approx(2.0)


def test_weighted_and_governed_windows():
    seg = pd.DataFrame({"status": ["computed", "computed", "not computed: x"], "governing_pulse": [1.0, 1.0, 1.0],
                        "fit_start_d": [0.0, 2.0, 5.0], "fit_end_d": [1.0, 5.0, 5.1], "sub_segment": [1, 2, 3]})
    w = P.governed_windows(seg, 1)
    assert w["sub_segment"].tolist() == [1, 2] and w["duration_d"].tolist() == [1.0, 3.0]
    assert P.weighted(pd.Series([1.0, 3.0]), w["duration_d"]) == pytest.approx(2.5)


def test_stirrer_energy():
    assert P.stirrer_energy_J(0.5, 1.0, 2.0) == pytest.approx(0.5 * 86400)
    with pytest.raises(ValueError):
        P.stirrer_energy_J(-1.0, 1.0, 2.0)


def _state(rate=1.0):
    t = np.round(np.arange(1.005, 47.92, 0.01), 10)
    return pd.DataFrame({"t_d": t, "model_rate_mgN_L_d": rate, "model_state_C_mgN_L_conditional_on_C0": 100.0})


@pytest.fixture(scope="module")
def real_inputs(cfg, ingested, validated, segmented):
    th = thermo.compute_stage(cfg, {"windows": validated.windows, "pulses": validated.pulse_table,
                                    "balances": validated.balances, "segments": segmented.table,
                                    "reference_windows": segmented.reference_windows})
    truth, vectors, summary = synthetic_truth(ingested.grid)  # W1: a known inverted NH4 truth (5 x reading)
    return {"segments": segmented.table, "steps": th["steps"], "oxygen": th["oxygen"], "pulses": validated.pulse_table,
            "model": _state(1.0), "model_fit": pd.DataFrame({"model": ["biomass_decay"], "selected": [True]}),
            "grid": ingested.grid, "truth": truth, "truth_vectors": vectors, "truth_summary": summary}


def test_every_metric_computed_with_inputs_or_not_computed(cfg, real_inputs):
    m = P.pulse_metrics(cfg, real_inputs)
    assert set(m["pulse"]) == {1, 2, 3, 4}
    ok = m["status"] == "computed"
    assert (ok | m["status"].str.startswith("not computed")).all()
    assert (m.loc[ok, "inputs"].str.len() > 0).all() and m.loc[ok, "value"].notna().all()
    assert m.loc[~ok, "value"].isna().all()


def test_required_metrics_present(cfg, real_inputs):
    m = P.pulse_metrics(cfg, real_inputs)
    for metric in ("delta_A_kJ_mol_SCENARIO_T25C", "consumption_uncorrected_weighted_mgN_L_d",
                   "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d", "median_pH", "fa_median_mgNH3_L_SCENARIO_T20C",
                   "fa_median_mgNH3_L_SCENARIO_T25C", "o2_demand_mgO2_L_d", "kla_required_per_d_SCENARIO_T25C",
                   "eta_conversion_fraction_MODEL_BASED", "specific_stirring_energy_J_per_gN",
                   "free_energy_to_stirrer_energy_ratio"):
        assert (m["metric"] == metric).sum() == 4, metric


def test_labels_and_missing_markers(cfg, real_inputs):
    m = P.pulse_metrics(cfg, real_inputs).set_index(["pulse", "metric"])
    assert "MODEL-BASED" in m.loc[(2, "eta_conversion_fraction_MODEL_BASED"), "basis"]
    assert "GAIN-CORRECTED" in m.loc[(2, "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"), "basis"]
    assert "SCENARIO" in m.loc[(2, "fa_median_mgNH3_L_SCENARIO_T25C"), "basis"]
    assert m.loc[(3, "specific_stirring_energy_J_per_gN"), "status"] == "not computed: stirrer_power_W missing"
    assert m.loc[(3, "kla_measured_over_required"), "status"] == "not computed: kLa_per_d missing"
    # pulse 4's last slope window (44.76-47.92 d) has a rising electrode signal: excluded from the weighted rates
    assert "excludes sub-segment(s)" in m.loc[(4, "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"), "basis"]


def test_eta_from_model_state(cfg, real_inputs):
    m = P.pulse_metrics(cfg, real_inputs).set_index(["pulse", "metric"])["value"]
    dose = thermo.pulse_dose_from_config(cfg).n_mg_L
    assert m[(2, "eta_conversion_fraction_MODEL_BASED")] == pytest.approx((24.59 - 10.57) / dose, rel=1e-3)
    assert m[(2, "free_energy_released_MODEL_BASED_J_SCENARIO_T25C")] == pytest.approx(
        m[(2, "eta_conversion_fraction_MODEL_BASED")] * 163.018, rel=1e-3)


def test_energy_metrics_with_measured_stirrer_power(tmp_path, real_inputs):
    text = CONFIG_PATH.read_text().replace('"stirrer_power_W", ', "") + "\n[measurements]\nstirrer_power_W = 0.5\n"
    p = tmp_path / "experiment.toml"
    p.write_text(text)
    c = load_config(p)
    m = P.pulse_metrics(c, real_inputs).set_index(["pulse", "metric"])
    energy = 0.5 * (24.59 - 10.57) * 86400
    assert m.loc[(2, "stirrer_energy_J"), "value"] == pytest.approx(energy)
    g_n = (24.59 - 10.57) * 1.0 * 0.1 / 1e3
    assert m.loc[(2, "specific_stirring_energy_J_per_gN"), "value"] == pytest.approx(energy / g_n, rel=1e-3)
    assert m.loc[(2, "free_energy_to_stirrer_energy_ratio"), "status"] == "computed"


def test_associations_are_exploratory(cfg, real_inputs):
    a = P.associations(P.pulse_metrics(cfg, real_inputs))
    assert (a["status"] == P.EXPLORATORY).all() and "p_value" not in a.columns
    assert (a["n"] <= 4).all()


def test_load_inputs_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        P.load_inputs(tmp_path)


def test_negative_consumption_excluded_and_not_computed(cfg, real_inputs):
    seg = real_inputs["segments"]
    m = P.pulse_metrics(cfg, real_inputs).set_index(["pulse", "metric"])
    # pulse 4: the non-positive window is excluded from the duration-weighted rate
    w = P.governed_windows(seg, 4)
    pos = w[w["gain_correction_status"] == "computed"]
    assert len(pos) < len(w) and (w.loc[~w.index.isin(pos.index), "gain_correction_status"] == S.NON_POSITIVE).all()
    expected = float((pos["consumption_SENSITIVITY_linear_gain_mgN_L_d"] * pos["duration_d"]).sum() / pos["duration_d"].sum())
    assert m.loc[(4, "consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"), "value"] == pytest.approx(expected)
    # every window of pulse 3 made non-positive: weighted rates and O2 are not computed
    bad = seg.copy()
    k3 = (bad["governing_pulse"] == 3) & (bad["status"] == "computed")
    bad.loc[k3, ["slope_mgN_L_d", "slope_mgNH4_L_d"]] = 0.5
    bad.loc[k3, "consumption_SENSITIVITY_linear_gain_mgN_L_d"] = np.nan
    bad.loc[k3, "gain_correction_status"] = S.NON_POSITIVE
    m3 = P.pulse_metrics(cfg, {**real_inputs, "segments": bad}).set_index(["pulse", "metric"])
    for metric in ("consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d", "consumption_uncorrected_weighted_mgN_L_d",
                   "consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d", "o2_demand_SENSITIVITY_linear_gain_mgO2_L_d"):
        assert m3.loc[(3, metric), "status"].startswith("not computed: non-positive"), metric
        assert np.isnan(m3.loc[(3, metric), "value"])


def test_rates_carry_segment_intervals_and_held_fraction(cfg, real_inputs):
    seg = real_inputs["segments"].set_index("sub_segment")
    m = P.pulse_metrics(cfg, real_inputs).set_index(["pulse", "metric"])
    first = P.governed_windows(real_inputs["segments"], 2).iloc[0]
    r = m.loc[(2, "consumption_SENSITIVITY_linear_gain_first_window_mgN_L_d")]
    s = seg.loc[first["sub_segment"]]
    assert (r["ci95_low"], r["ci95_high"]) == (s["consumption_SENSITIVITY_linear_gain_ci95_low_mgN_L_d"],
                                               s["consumption_SENSITIVITY_linear_gain_ci95_high_mgN_L_d"])
    assert r["held_fraction"] == pytest.approx(s["held_fraction"]) and r["held_fraction"] > 0.5
    u = m.loc[(2, "consumption_uncorrected_first_window_mgN_L_d")]
    assert u["ci95_low"] == pytest.approx(-s["slope_ci95_high_mgN_L_d"]) and u["ci95_high"] == pytest.approx(-s["slope_ci95_low_mgN_L_d"])
    rates = m[m.index.get_level_values("metric").str.startswith(("consumption_", "o2_demand", "kla_required"))]
    done = rates[rates["status"] == "computed"]
    assert len(done) > 0 and done["ci95_low"].notna().all() and done["ci95_high"].notna().all()
    assert (done["ci95_low"] <= done["value"]).all() and (done["value"] <= done["ci95_high"]).all()
    assert done["held_fraction"].notna().all()


def test_kla_margin_uses_required_kla(measured_cfg, real_inputs):
    c = measured_cfg(kLa_per_d=4.0)
    m = P.pulse_metrics(c, real_inputs).set_index(["pulse", "metric"])
    req = m.loc[(2, "kla_required_per_d_SCENARIO_T25C"), "value"]
    assert m.loc[(2, "kla_measured_over_required"), "value"] == pytest.approx(4.0 / req)
    assert m.loc[(2, "kla_measured_over_required"), "value"] != pytest.approx(4.0)  # not divided by a 1.0 placeholder


def test_pulses_compute_on_real_run_all_outputs(cfg, run_all_dir):
    r = P.compute(cfg, P.load_inputs(run_all_dir))
    assert set(r.metrics["pulse"]) == {1, 2, 3, 4}
    ok = r.metrics["status"] == "computed"
    assert (ok | r.metrics["status"].str.startswith("not computed")).all()
    assert r.metrics.to_csv(index=False) == (run_all_dir / "tables" / "pulses.csv").read_text()
