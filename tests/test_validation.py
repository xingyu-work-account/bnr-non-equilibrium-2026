import numpy as np
import pandas as pd
import pytest

from bnr_neq import validation


def test_window_median_empty_window_raises():
    with pytest.raises(ValueError, match="no samples"):
        validation.window_median(np.array([0.0, 1.0]), np.array([1.0, 2.0]), 0.2, 0.8, "x")


def test_window_is_left_closed():
    t = np.array([0.0, 0.1, 0.2])
    assert validation.window_median(t, np.array([1.0, 2.0, 3.0]), 0.1, 0.2, "x") == (2.0, 1)


def test_pulse_table_shape_and_labels(validated, cfg):
    pt = validated.pulse_table
    assert pt["pulse"].tolist() == [1, 2, 3, 4]
    assert pt["t_pulse_d"].tolist() == cfg.pulse_times_d
    scen = [c for c in pt.columns if c.startswith("orp_pred_change_mV_")]
    assert len(scen) == 8 and all("SCENARIO_" in c for c in scen)
    assert (pt["ph_pred_change"] == "not computed: alkalinity_mmol_L missing").all()
    np.testing.assert_allclose(pt["nh4_pred_change_mgNH4_L"], 84.31, atol=5e-3)
    np.testing.assert_allclose(pt["ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C"], 0.700, atol=5e-4)


def test_orp_scenario_predictions_shrink_with_cumulative_chloride(validated):
    col = "orp_pred_change_mV_SCENARIO_Cl0_1mM_T25C"
    assert validated.pulse_table[col].iloc[0] == pytest.approx(44.6, abs=0.05)
    assert validated.pulse_table[col].is_monotonic_decreasing


def test_windows_have_samples(validated):
    assert (validated.windows["n"] > 0).all()
    assert len(validated.windows) == 4 * 4 * 2


REFERENCE_QUANTITIES_EXPECTED_TO_PASS = [
    "nh4_response_fraction_pct", "nh4_obs_change_mgNH4_L", "ec_obs_change_mS_cm", "orp_obs_change_mV",
    "ph_obs_change", "ec_median_day_0.3_1.3_mS_cm", "ec_median_day_48.5_49.6_mS_cm",
    "nh4_median_start_mgNH4_L", "nh4_median_end_mgNH4_L", "nh4_delayed_response_day_pulse1",
    "nh4_delayed_response_day_pulse2", "ec_late_pulse_count",
]


@pytest.mark.parametrize("quantity", REFERENCE_QUANTITIES_EXPECTED_TO_PASS)
def test_reference_values_reproduced(validated, quantity):
    r = validated.references[validated.references["quantity"] == quantity]
    assert len(r) >= 1
    bad = r[(r["obtained"] - r["reference"]).abs() > r["tolerance"] + 1e-9]
    assert bad.empty, bad.to_string()


def test_offset_jumps_lie_in_the_stamp_where_the_raw_value_changes(validated, ingested):
    """Independent of the detector: find, in the raw stamped rows, the consecutive pair with the
    largest absolute change near each reference day, and require the detected jump in that stamp."""
    nh4 = ingested.natives["nh4"]
    offsets = validated.jumps[validated.jumps["jump"] & (validated.jumps["class"] == "electrode_offset")]
    for ref in (2.59, 17.84, 37.19, 44.63):
        near = nh4[(nh4["t_stamp_d"] > ref - 0.06) & (nh4["t_stamp_d"] < ref + 0.06)]
        k = int(np.argmax(np.abs(np.diff(near["nh4_mgNH4_L"].to_numpy()))))
        stamp = near["t_stamp_d"].iloc[k + 1]
        detected = offsets["t_jump_d"].to_numpy()
        tj = detected[np.argmin(np.abs(detected - ref))]
        assert stamp <= tj < stamp + 0.01, (ref, stamp, tj)


def test_ec_late_pulse_threshold_rule(validated):
    sens = validated.ec_sensitivity
    assert sens["rule"].iloc[0] == "noise floor (used)"
    assert sens["n_pulses"].iloc[0] == len(validated.ec_pulses)


def test_verdicts(validated):
    v = validated.verdicts
    assert v["nh4"].startswith("FAIL") and v["ec"].startswith("FAIL") and v["orp"].startswith("FAIL")
    assert v["ph"].startswith("PASS")


def test_balances_n(validated):
    b = validated.balances.set_index("quantity")["value"]
    assert b["n_added_mgN_L"] == pytest.approx(4 * 65.46, abs=0.02)
    assert b["n_removed_from_nh4_pool_mgN_L"] == pytest.approx(
        b["nh4_median_start_mgN_L"] + b["n_added_mgN_L"] - b["nh4_median_end_mgN_L"])


def test_run_writes_outputs(cfg, ingested, tmp_path):
    from bnr_neq import pipeline

    pipeline.write_ingest(ingested, tmp_path)
    validation.run(cfg, tmp_path)
    for rel in ("tables/sensor_validation.csv", "figures/fig3_sensor_validation.png",
                "figures/fig3_sensor_validation.svg", "reports/sensor_validation.md"):
        assert (tmp_path / rel).is_file()
    pt = pd.read_csv(tmp_path / "tables/sensor_validation.csv")
    assert len(pt) == 4


def test_run_without_ingest_raises(cfg, tmp_path):
    with pytest.raises(FileNotFoundError, match="ingest"):
        validation.run(cfg, tmp_path)


def test_clock_check_ph_and_orp_steps_coincide(validated):
    c = validated.clock.set_index("sensor")
    assert abs(c.loc["ph", "offset_vs_orp_d"]) < 0.02
