"""W1: one NH4 observation model (the M6 log model) for every downstream quantity."""

import math

import numpy as np
import pandas as pd
import pytest
from conftest import synthetic_truth

from bnr_neq import model as M
from bnr_neq import pulses as P
from bnr_neq import segments as S
from bnr_neq import thermo, volume

T = np.round(np.arange(1.005, 7.0, 0.01), 10)


# --- inversion ----------------------------------------------------------------------------
def test_inversion_recovers_the_state_without_lag():
    c = 50.0 + 30.0 * np.sin(T)
    p = {"a": -1.2, "b0": 1.3, "b1_per_d": -0.002, "tau_d": 1e-9, "offset_values": [0.3, -0.2]}
    offsets = [2.5, 4.0]
    log_r = M.observe(c, T, p, offsets)
    np.testing.assert_allclose(M.invert_observation(log_r, T, p, offsets), c, rtol=1e-9)


def test_inversion_rejects_a_non_positive_gain():
    p = {"a": 0.0, "b0": 0.01, "b1_per_d": -0.01, "offset_values": []}
    with pytest.raises(ValueError, match="gain"):
        M.invert_observation(np.zeros(len(T)), T, p, [])


def test_inverted_series_from_profile_vectors():
    truth, vectors, _ = synthetic_truth(pd.DataFrame({
        "t_start_d": T - 0.005, "t_end_d": T + 0.005, "t_mid_d": T, "nh4_mgN_L": 10.0 + T, "pH": 7.3}))
    series = P.band_series(truth, vectors)
    np.testing.assert_allclose(series[0], truth["C_inv_mgN_L"])
    np.testing.assert_allclose(series[1], 4.0 * truth["electrode_obs_mgN_L"])


def test_truth_id_identifies_the_series():
    c = np.linspace(1.0, 2.0, 50)
    assert M.truth_id(c) == M.truth_id(c.copy()) and M.truth_id(c) != M.truth_id(c * (1 + 1e-12))


# --- negative: mixing observation models raises ----------------------------------------------
def _rows(*specs):
    return pd.DataFrame([{"metric": m, "basis": b, "nh4_observation_model": o} for m, b, o in specs])


def test_mixing_observation_models_in_one_output_raises():
    with pytest.raises(M.ObservationModelError, match="mixes"):
        M.check_single_observation_model(_rows(("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL),
                                               ("o2_demand_mgO2_L_d", "gain", "linear_gain")))


def test_unlabelled_sensitivity_row_raises():
    with pytest.raises(M.ObservationModelError, match="SENSITIVITY"):
        M.check_single_observation_model(_rows(("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL),
                                               ("consumption_gain_corrected_weighted_mgN_L_d", "GAIN-CORRECTED",
                                                M.LINEAR_GAIN)))


def test_derived_quantity_from_the_raw_reading_raises():
    with pytest.raises(M.ObservationModelError, match="raw electrode reading"):
        M.check_single_observation_model(_rows(("fa_median_mgNH3_L_SCENARIO_T25C", "reading", M.RAW_READING)))


def test_check_accepts_empty_cells_read_from_csv(tmp_path):
    # pass 3 (X3): the model cell must declare a basis; empty label cells (NaN read back from CSV) are still accepted
    df = _rows(("t_pulse_d", "", M.NH4_INDEPENDENT), ("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL))
    df.to_csv(tmp_path / "p.csv", index=False)
    back = pd.read_csv(tmp_path / "p.csv")
    assert back["basis"].isna().any()
    assert M.check_single_observation_model(back) == M.OBSERVATION_MODEL


def test_labelled_sensitivity_and_measured_rows_pass():
    # pass 3 (X3, X4): the gain metric says SENSITIVITY in its name, the fact row declares itself NH4-independent
    df = _rows(("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL),
               ("consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d", "SENSITIVITY: GAIN-CORRECTED", M.LINEAR_GAIN),
               ("consumption_uncorrected_weighted_mgN_L_d", "measured", M.RAW_READING),
               ("t_pulse_d", "fact", M.NH4_INDEPENDENT))
    assert M.check_single_observation_model(df) == M.OBSERVATION_MODEL


# --- boundary: a short window returns a status ------------------------------------------------
def test_short_window_returns_a_status_not_a_number(cfg):
    grid = pd.DataFrame({"t_start_d": T - 0.005, "t_end_d": T + 0.005, "t_mid_d": T, "nh4_mgN_L": 40.0 - T,
                         "pH": 7.3})
    truth, vectors, _ = synthetic_truth(grid)
    table = pd.DataFrame({"sub_segment": [1, 2], "segment": [1, 1], "governing_pulse": [1, 1],
                          "fit_start_d": [2.0, 4.0], "fit_end_d": [2.2, 6.0], "held_fraction": [0.0, 0.0],
                          "status": ["computed", "computed"]})
    tw = S.truth_slope_table(cfg, truth, table, P.band_series(truth, vectors))
    short, ok = tw.iloc[0], tw.iloc[1]
    assert short["status"].startswith("not computed: ") and "(< 25)" in short["status"]
    assert np.isnan(short["consumption_mgN_L_d"]) and np.isnan(short["o2_demand_mgO2_L_d"])
    assert ok["status"] == "computed" and ok["consumption_mgN_L_d"] == pytest.approx(5.0)  # C = 5 x (40 − t)
    assert (ok["consumption_param_band_low_mgN_L_d"], ok["consumption_param_band_high_mgN_L_d"]) == pytest.approx((4.0, 6.0))
    assert ok["o2_demand_mgO2_L_d"] == pytest.approx(4.57 * 5.0)


def test_rising_truth_is_not_a_consumption(cfg):
    grid = pd.DataFrame({"t_start_d": T - 0.005, "t_end_d": T + 0.005, "t_mid_d": T, "nh4_mgN_L": 10.0 + T,
                         "pH": 7.3})
    truth, _, _ = synthetic_truth(grid)
    table = pd.DataFrame({"sub_segment": [1], "segment": [1], "governing_pulse": [1], "fit_start_d": [2.0],
                          "fit_end_d": [5.0], "held_fraction": [0.0], "status": ["computed"]})
    tw = S.truth_slope_table(cfg, truth, table)
    assert tw["status"].iloc[0] == S.NON_POSITIVE_TRUTH and np.isnan(tw["consumption_mgN_L_d"].iloc[0])


def test_tampered_truth_id_raises(cfg):
    grid = pd.DataFrame({"t_start_d": T - 0.005, "t_end_d": T + 0.005, "t_mid_d": T, "nh4_mgN_L": 40.0 - T,
                         "pH": 7.3})
    truth, _, _ = synthetic_truth(grid)
    truth.loc[5, "C_inv_mgN_L"] *= 1.01
    table = pd.DataFrame({"sub_segment": [1], "segment": [1], "governing_pulse": [1], "fit_start_d": [2.0],
                          "fit_end_d": [5.0], "held_fraction": [0.0], "status": ["computed"]})
    with pytest.raises(ValueError, match="does not match"):
        S.truth_slope_table(cfg, truth, table)


def test_evaporation_adds_the_concentration_term(cfg):
    grid = pd.DataFrame({"t_start_d": T - 0.005, "t_end_d": T + 0.005, "t_mid_d": T, "nh4_mgN_L": 40.0 - T,
                         "pH": 7.3})
    truth, _, _ = synthetic_truth(grid)
    table = pd.DataFrame({"sub_segment": [1], "segment": [1], "governing_pulse": [1], "fit_start_d": [2.0],
                          "fit_end_d": [5.0], "held_fraction": [0.0], "status": ["computed"]})
    vols = volume.volume_scenarios(cfg, 47.92)
    tw = S.truth_slope_table(cfg, truth, table, volumes=vols).iloc[0]
    e1 = next(v for v in vols if v.net_loss_L_per_d == 0.001)
    col = f"consumption_mgN_L_d_{volume.slug(e1.name)}"
    expected = 5.0 + tw["median_C_inv_mgN_L"] * 0.001 / float(e1.volume_L(3.5))
    assert tw[col] == pytest.approx(expected)


# --- stage: the pulse table rests on one inverted C(t) --------------------------------------
@pytest.fixture(scope="module")
def synth_inputs(cfg, ingested, validated, segmented):
    th = thermo.compute_stage(cfg, {"windows": validated.windows, "pulses": validated.pulse_table,
                                    "balances": validated.balances, "segments": segmented.table,
                                    "reference_windows": segmented.reference_windows})
    truth, vectors, summary = synthetic_truth(ingested.grid)
    t = np.round(np.arange(1.005, 47.92, 0.01), 10)
    return {"segments": segmented.table, "steps": th["steps"], "oxygen": th["oxygen"], "pulses": validated.pulse_table,
            "model": pd.DataFrame({"t_d": t, "model_rate_mgN_L_d": 1.0, "model_state_C_mgN_L_conditional_on_C0": 100.0}),
            "model_fit": pd.DataFrame({"model": ["biomass_decay"], "selected": [True]}), "grid": ingested.grid,
            "truth": truth, "truth_vectors": vectors, "truth_summary": summary}


def test_primary_rates_scale_with_the_truth(cfg, synth_inputs):
    m = P.pulse_metrics(cfg, synth_inputs).set_index(["pulse", "metric"])
    # C = 5 x reading: Theil-Sen is scale-equivariant, so the primary first-window rate is 5x the measured one
    for k in (2, 3):
        assert m.loc[(k, "consumption_first_window_mgN_L_d"), "value"] == pytest.approx(
            5.0 * m.loc[(k, "consumption_uncorrected_first_window_mgN_L_d"), "value"], rel=1e-9)
        assert m.loc[(k, "consumption_first_window_mgN_L_d"), "param_band_low"] == pytest.approx(
            4.0 * m.loc[(k, "consumption_uncorrected_first_window_mgN_L_d"), "value"], rel=1e-9)


def test_dose_anchored_eta_and_affinity(cfg, synth_inputs):
    m = P.pulse_metrics(cfg, synth_inputs).set_index(["pulse", "metric"])
    truth = synth_inputs["truth"]
    t, c = truth["t_d"].to_numpy(), truth["C_inv_mgN_L"].to_numpy()
    dose = thermo.pulse_dose_from_config(cfg).n_mg_L
    pre1, pre2, pre3 = P.pre_median(t, c, 8.53), P.pre_median(t, c, 10.57), P.pre_median(t, c, 24.59)
    eta = (pre2 + dose - pre3) / dose
    row = m.loc[(2, "eta_conversion_fraction")]
    ratio = m.loc[(2, "removal_to_dose_ratio")]
    # pinned (pass 3, X8): C = 5 x reading falls by more than the dose over pulse 2 and rises by more than it over pulse 1
    assert eta == pytest.approx(1.96, abs=0.01)
    assert row["status"].startswith("not computed: removal exceeds the dose") and np.isnan(row["value"])
    assert ratio["status"] == "computed" and ratio["value"] == pytest.approx(eta)
    eta1 = (pre1 + dose - pre2) / dose
    assert eta1 == pytest.approx(-0.77, abs=0.01)
    assert m.loc[(1, "eta_conversion_fraction"), "status"].startswith("not computed: negative dose-anchored removal")
    assert m.loc[(1, "removal_to_dose_ratio"), "status"].startswith("not computed: negative dose-anchored removal")
    assert m.loc[(2, "delta_A_kJ_mol_SCENARIO_T25C"), "value"] == pytest.approx(
        float(thermo.affinity_step_kJ_mol(pre2, dose)))
    assert m.loc[(2, "delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T25C"), "nh4_observation_model"] == M.IDENTITY


def test_pulse_metrics_carry_provenance(cfg, synth_inputs):
    m = P.pulse_metrics(cfg, synth_inputs)
    tid = synth_inputs["truth"]["nh4_truth_id"].iloc[0]
    primary = m[m["nh4_observation_model"] == M.OBSERVATION_MODEL]
    derived = primary[primary["metric"].isin(["consumption_weighted_mgN_L_d", "fa_median_mgNH3_L_SCENARIO_T25C",
                                              "o2_demand_mgO2_L_d", "delta_A_kJ_mol_SCENARIO_T25C"])]
    assert len(derived) == 16 and (derived["nh4_series"] == tid).all()
    sens = m[m["nh4_observation_model"].str.endswith(M.SENSITIVITY_SUFFIX)]
    assert len(sens) > 0 and (sens["metric"] + sens["basis"]).str.contains("SENSITIVITY").all()
    assert M.check_single_observation_model(m) == M.OBSERVATION_MODEL
    # a dose-anchored conversion is reported only when positive, under every volume scenario
    eta = m[m["metric"].str.startswith("eta_conversion_fraction") & (m["metric"] != "eta_conversion_fraction_MODEL_BASED")]
    computed = eta.loc[eta["status"] == "computed", "value"]
    assert len(eta) == 16 and ((computed >= 0) & (computed <= 1)).all()  # S3: a fraction of the pulse
    assert eta.loc[eta["status"] != "computed", "status"].str.startswith(
        ("not computed: negative", "not computed: removal exceeds the dose")).all()


def test_pulse_metrics_refuse_a_mixed_output(cfg, synth_inputs, monkeypatch):
    monkeypatch.setattr(M, "LINEAR_GAIN", "linear_gain")  # the gain rows lose their SENSITIVITY suffix
    with pytest.raises(M.ObservationModelError):
        P.pulse_metrics(cfg, synth_inputs)


def test_dose_anchored_bound_table(cfg, synth_inputs):
    b = P.truth_balance(cfg, synth_inputs)
    const = b[b["volume_scenario"] == volume.CONSTANT].set_index("checkpoint")
    assert list(const.index) == ["pre_pulse1", "pre_pulse2", "pre_pulse3", "pre_pulse4", "end"]
    dose = thermo.pulse_dose_from_config(cfg).n_mg_L
    c0 = float(const["C_t0_mgN_L"].iloc[0])
    assert const.loc["end", "C_dose_upper_mgN_L"] == pytest.approx(c0 + 4 * dose)
    assert const.loc["pre_pulse1", "C_dose_upper_mgN_L"] == pytest.approx(c0)
    assert (const["C_dose_lower_alk0_2_mgN_L"] >= const["C_dose_lower_alk0_15_mgN_L"]).all()
    assert (const["C_dose_lower_alk0_15_mgN_L"] <= const["C_dose_upper_mgN_L"]).all()
    removed = max(float(const.loc["end", "removed_per_initial_litre_mgN_L"]), 0.0)
    assert (const.loc["end", "alk0_required_if_nitrified_meq_L"] - const.loc["end", "alk0_required_any_proton_sink_meq_L"]
            == pytest.approx(removed / 14.007))


# --- integration on the real run-all outputs ------------------------------------------------------
def test_one_pulse_derives_everything_from_the_same_inverted_series(run_all_dir):
    truth = pd.read_parquet(run_all_dir / "processed" / "nh4_truth.parquet")
    tid = M.truth_id(truth["C_inv_mgN_L"].to_numpy())
    assert (truth["nh4_truth_id"] == tid).all() and (truth["nh4_observation_model"] == M.OBSERVATION_MODEL).all()
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    p2 = pm[(pm["pulse"] == 2) & pm["metric"].isin(["consumption_weighted_mgN_L_d", "fa_median_mgNH3_L_SCENARIO_T25C",
                                                     "eta_conversion_fraction", "o2_demand_mgO2_L_d"])]
    assert len(p2) == 4
    assert (p2["nh4_series"] == tid).all() and (p2["nh4_observation_model"] == M.OBSERVATION_MODEL).all()
    tw = pd.read_csv(run_all_dir / "tables" / "nh4_truth_windows.csv")
    assert (tw["nh4_truth_id"] == tid).all()
    # the O2 demand is 4.57 x the duration-weighted consumption of the same windows
    w = tw[(tw["governing_pulse"] == 2) & (tw["status"] == "computed")]
    d = w["fit_end_d"] - w["fit_start_d"]
    rate = float((w["consumption_mgN_L_d"] * d).sum() / d.sum())
    v = p2.set_index("metric")["value"]
    assert v["consumption_weighted_mgN_L_d"] == pytest.approx(rate)
    assert v["o2_demand_mgO2_L_d"] == pytest.approx(4.57 * rate)
    assert M.check_single_observation_model(pm) == M.OBSERVATION_MODEL
    summary = pd.read_csv(run_all_dir / "tables" / "nh4_truth_summary.csv").set_index("key")["value"]
    assert summary["nh4_truth_id"] == tid


def test_real_band_contains_the_best_series(run_all_dir):
    truth = pd.read_parquet(run_all_dir / "processed" / "nh4_truth.parquet")
    assert (truth["C_inv_band_low_mgN_L"] <= truth["C_inv_mgN_L"] + 1e-9).all()
    assert (truth["C_inv_mgN_L"] <= truth["C_inv_band_high_mgN_L"] + 1e-9).all()
    v = pd.read_csv(run_all_dir / "tables" / "nh4_truth_vectors.csv")
    assert v["vector"].iloc[0] == 0 and (v["delta"] <= M.CHI2_95_1DOF).all()
    assert math.isfinite(float(v["a"].iloc[0]))
