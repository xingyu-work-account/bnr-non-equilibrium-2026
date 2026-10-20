"""Section 6 functions (M4), guard paths and the thermo stage."""

import numpy as np
import pytest

from bnr_neq import thermo as th
from bnr_neq.config import CONFIG_PATH, load_config
from bnr_neq.constants import LOG10_KH_CO2_25C, PK1_CARBONATE_25C

RT25 = 8.314462618 * 298.15 / 1e3


def cfg_with(tmp_path, **measured):
    """Configuration with the given keys moved from [missing] to [measurements]."""
    text = CONFIG_PATH.read_text()
    for k in measured:
        assert f'"{k}"' in text
        text = text.replace(f'"{k}", ', "").replace(f', "{k}"', "")
    text += "\n[measurements]\n" + "".join(f"{k} = {v}\n" for k, v in measured.items())
    p = tmp_path / "experiment.toml"
    p.write_text(text)
    return load_config(p)


# --- unit tests required by M4 ----------------------------------------------------
def test_davies_gamma_z1_I001():
    assert th.davies_gamma(1, 0.01) == pytest.approx(0.902, abs=0.002)
    assert th.davies_gamma(1, 0.0) == 1.0


def test_agcl_reference_potential_and_slope():
    assert th.e_ref_agcl_V(1.0) == pytest.approx(0.2223, abs=1e-12)
    slope_mV = 1e3 * (th.e_ref_agcl_V(0.01) - th.e_ref_agcl_V(0.1))
    assert slope_mV == pytest.approx(59.16, abs=0.01)


def test_oxide_reference():
    assert th.e_ref_ag2o_V(7.3) == pytest.approx(1.173 - 0.0591593 * 7.3, abs=1e-6)


def test_pka_25C():
    assert th.emerson_pka(298.15) == pytest.approx(9.246, abs=0.002)


@pytest.mark.parametrize(
    "call",
    [
        lambda: th.davies_gamma(1, -0.01),
        lambda: th.e_ref_agcl_V(0.0),
        lambda: th.e_ref_ag2o_V(15.0),
        lambda: th.nitrification_affinity_kJ_mol(40, 100, 0.1, 14.2),
        lambda: th.nitrification_affinity_kJ_mol(-40, 100, 0.1, 7.0),
        lambda: th.nitrification_affinity_kJ_mol(40, -1, 0.1, 7.0),
        lambda: th.nitrification_affinity_kJ_mol(40, 100, 0.0, 7.0),
        lambda: th.affinity_step_kJ_mol(0.0, 65.46),
        lambda: th.affinity_step_kJ_mol(10.0, -1.0),
        lambda: th.hco3_air_eq_mmol_L(-1.0, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, 4.2e-4),
        lambda: th.omega_co2(-0.5, 7.0, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, 4.2e-4),
        lambda: th.oxygen_demand_mgO2_L_d(-1.0),
        lambda: th.required_kla_per_d(5.0, 0.0),
        lambda: th.alkalinity_consumed_mgCaCO3_L(-1.0),
        lambda: th.entropy_production_W_L_K([1.0], [1.0], 0.0),
        lambda: th.entropy_production_W_L_K([1.0, 2.0], [1.0], 298.15),
        lambda: th.jin_bethke_FT(-1e5, 1.0, 5e4, 0.0, 298.15),
    ],
)
def test_impossible_inputs_raise(call):
    with pytest.raises(ValueError):
        call()


# --- equations ----------------------------------------------------------------------
def test_coexistence_affinity_reference():
    assert th.nitrification_affinity_kJ_mol(40, 100, 0.1, 7.0) == pytest.approx(317.0, abs=2.0)


def test_affinity_dependencies():
    a = th.nitrification_affinity_kJ_mol
    assert a(40, 100, 0.01, 7.0) == pytest.approx(a(40, 100, 0.1, 7.0) - 2 * RT25 * np.log(10), abs=1e-9)
    assert a(40, 100, 0.1, 8.0) == pytest.approx(a(40, 100, 0.1, 7.0) + 2 * RT25 * np.log(10), abs=1e-9)
    assert a(4000, 100, 0.1, 7.0) == pytest.approx(a(40, 100, 0.1, 7.0) + RT25 * np.log(100), abs=1e-9)
    assert th.nitrification_dG_kJ_mol(40, 100, 0.1, 7.0) == -a(40, 100, 0.1, 7.0)


def test_affinity_step():
    assert th.affinity_step_kJ_mol(15.1, 65.46) == pytest.approx(4.2, abs=0.05)
    assert th.affinity_step_kJ_mol(10.0, 0.0) == 0.0


def test_hco3_air_equilibrium_references():
    f = th.hco3_air_eq_mmol_L
    assert f(7.3, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, 4.2e-4) == pytest.approx(0.127, abs=0.001)
    assert f(7.8, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, 4.2e-4) == pytest.approx(0.401, abs=0.001)
    assert th.omega_co2(1.27, 7.3, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, 4.2e-4) == pytest.approx(10.0, rel=0.01)


def test_stoichiometry():
    assert th.oxygen_demand_mgO2_L_d(1.0) == pytest.approx(4.57)
    assert th.required_kla_per_d(8.26) == pytest.approx(1.0)
    assert th.alkalinity_consumed_mgCaCO3_L(65.46) / th.CACO3_MG_PER_MEQ == pytest.approx(9.35, abs=0.02)


def test_non_equilibrium_criteria():
    assert th.de_donder_consistent([1.0, -1.0, 2.0], [0.5, -0.1, -0.1]).tolist() == [True, True, False]
    assert th.entropy_production_W_L_K([2e-9, 1e-9], [3e5, 1e5], 300.0) == pytest.approx((6e-4 + 1e-4) / 300)
    assert th.jin_bethke_FT(-1e5, 2.0, 5e4, 1.0, 298.15) == pytest.approx(0.0, abs=1e-12)
    assert th.jin_bethke_FT(-3e5, 1.0, 5e4, 1.0, 298.15) == pytest.approx(1.0, abs=1e-12)


# --- every missing-input path -----------------------------------------------------------
def test_guards_report_missing_inputs(cfg):
    assert str(th.eh_from_orp_guarded(cfg, 150.0)) == "not computed: S4"
    assert str(th.e_ref_agcl_guarded(cfg, 0.9)) == "not computed: chloride_mmol_L missing"
    assert str(th.nitrification_affinity_guarded(cfg, 30.0, 7.3)) == "not computed: DO_mg_L missing"
    assert str(th.omega_co2_guarded(cfg, 7.3, PK1_CARBONATE_25C, LOG10_KH_CO2_25C)) == \
        "not computed: alkalinity_mmol_L missing"
    assert str(th.kla_margin_guarded(cfg, 1.0)) == "not computed: kLa_per_d missing"


def test_affinity_guard_second_key(tmp_path):
    c = cfg_with(tmp_path, DO_mg_L=2.0)
    assert str(th.nitrification_affinity_guarded(c, 30.0, 7.3)) == "not computed: NO3_N_mg_L missing"


def test_guards_compute_when_measured(tmp_path):
    c = cfg_with(tmp_path, DO_mg_L=0.1, NO3_N_mg_L=100.0, chloride_mmol_L=1000.0, orp_reference_offset_V=0.2,
                 alkalinity_mmol_L=1.27, kLa_per_d=4.0)
    assert th.nitrification_affinity_guarded(c, 40.0, 7.0) == pytest.approx(316.9, abs=0.1)
    assert th.e_ref_agcl_guarded(c, 1.0) == pytest.approx(0.2223)
    assert th.eh_from_orp_guarded(c, 150.0) == pytest.approx(0.35)
    assert th.omega_co2_guarded(c, 7.3, PK1_CARBONATE_25C, LOG10_KH_CO2_25C) == pytest.approx(10.0, rel=0.01)
    assert th.kla_margin_guarded(c, 2.0) == pytest.approx(2.0)


# --- stage on the real data ----------------------------------------------------------
@pytest.fixture(scope="module")
def stage(cfg, validated, segmented):
    inputs = {"windows": validated.windows, "pulses": validated.pulse_table, "balances": validated.balances,
              "segments": segmented.table, "reference_windows": segmented.reference_windows}
    return th.compute_stage(cfg, inputs)


def test_stage_references(stage):
    r = stage["references"].set_index("quantity")
    failing = sorted(r.index[~r["pass"]])
    assert failing == ["O2 demand max mg O2/L/d (SENSITIVITY: linear gain, as in the preliminary analysis)",
                       "kLa required max 1/d (SENSITIVITY: linear gain, as in the preliminary analysis)"]


def test_preliminary_window_diagnostic_explains_upper_references(stage):
    diag = stage["oxygen"][stage["oxygen"]["status"].str.startswith("DIAGNOSTIC")]
    assert len(diag) == 1
    assert diag["o2_demand_SENSITIVITY_linear_gain_mgO2_L_d"].iloc[0] == pytest.approx(33.0, rel=0.1)
    assert diag["kla_required_SENSITIVITY_linear_gain_per_d_SCENARIO_T25C"].iloc[0] == pytest.approx(4.0, rel=0.1)


def test_stage_missing_markers(stage):
    text = " ".join(str(v) for key in ("nitrogen", "chloride_balance", "alkalinity")
                    for v in stage[key]["value"])
    for key in ("lab_NH4_N_mg_L", "NO3_N_mg_L", "NO2_N_mg_L", "kLa_per_d", "chloride_mmol_L", "alkalinity_mmol_L"):
        assert f"not computed: {key} missing" in text
    assert (stage["chloride_observed"]["eh_absolute"] == "not computed: S4").all()
    assert (stage["co2"]["omega_co2_SCENARIO_T25C"] == "not computed: alkalinity_mmol_L missing").all()
    insitu = stage["coexistence"][stage["coexistence"]["basis"] == "IN SITU"]
    assert insitu["affinity_kJ_mol"].iloc[0] == "not computed: DO_mg_L missing"


def test_stage_scenarios_labelled(stage):
    assert stage["chloride_scenarios"]["scenario"].str.startswith("SCENARIO_").all()
    assert (stage["coexistence"].loc[stage["coexistence"]["basis"] != "IN SITU", "basis"] == "SCENARIO").all()
    assert any("SCENARIO_T20C" in c for c in stage["steps"].columns)


def test_oxygen_not_computed_for_non_positive_rates(stage, segmented):
    # non-positive gain-corrected consumption is "not computed" in segments.csv and stays so in thermo_oxygen.csv
    o = stage["oxygen"]
    nonpos = segmented.table.loc[segmented.table["gain_correction_status"].str.startswith("not computed: non-positive"),
                                 "sub_segment"]
    neg = o[o["sub_segment"].isin(nonpos)]
    assert len(neg) >= 1 and neg["o2_demand_SENSITIVITY_linear_gain_mgO2_L_d"].isna().all()
    assert neg["status"].str.startswith("not computed: non-positive").all()
    assert neg["consumption_SENSITIVITY_linear_gain_mgN_L_d"].isna().all()


def test_run_and_missing_inputs(cfg, ingested, validated, segmented, tmp_path):
    with pytest.raises(FileNotFoundError, match="validate"):
        th.run(cfg, tmp_path)
    t = tmp_path / "tables"
    t.mkdir()
    validated.windows.to_csv(t / "sensor_validation_windows.csv", index=False)
    validated.pulse_table.to_csv(t / "sensor_validation.csv", index=False)
    validated.balances.to_csv(t / "sensor_balances.csv", index=False)
    segmented.table.to_csv(t / "segments.csv", index=False)
    segmented.reference_windows.to_csv(t / "segments_reference_windows.csv", index=False)
    with pytest.raises(FileNotFoundError, match="ingest"):  # W2 needs the pH grid
        th.run(cfg, tmp_path)
    (tmp_path / "processed").mkdir()
    ingested.grid.to_parquet(tmp_path / "processed" / "grid_0p01d.parquet", index=False)
    th.run(cfg, tmp_path)
    for name in ("thermo_affinity_steps", "thermo_coexistence", "balance_nitrogen", "thermo_references",
                 "alkalinity_bound", "alkalinity_references", "volume_scenarios"):
        assert (t / f"{name}.csv").is_file()
    assert (tmp_path / "reports" / "thermo.md").is_file() and (tmp_path / "reports" / "alkalinity.md").is_file()


def test_guarded_stage_paths_use_observed_data(measured_cfg, validated, segmented):
    """With the measurements present, the guarded stage paths are evaluated on the observed series,
    never on the former placeholders (NH4 1.0 mg N/L, pH 7.0, ORP 0.0 mV)."""
    c = measured_cfg(DO_mg_L=0.1, NO3_N_mg_L=100.0, alkalinity_mmol_L=1.27, orp_reference_offset_V=0.2)
    inputs = {"windows": validated.windows, "pulses": validated.pulse_table, "balances": validated.balances,
              "segments": segmented.table, "reference_windows": segmented.reference_windows}
    st = th.compute_stage(c, inputs)
    comp = segmented.table[segmented.table["status"] == "computed"]
    insitu = st["coexistence"][st["coexistence"]["basis"] == "IN SITU"]
    assert len(insitu) == len(comp)
    expected = th.nitrification_affinity_kJ_mol(comp["median_nh4_mgN_L"].to_numpy(), 100.0, 0.1, comp["median_pH"].to_numpy())
    np.testing.assert_allclose(insitu["affinity_kJ_mol"].astype(float), expected)
    np.testing.assert_allclose(insitu["nh4_mgN_L"], comp["median_nh4_mgN_L"])
    placeholder = float(th.nitrification_affinity_kJ_mol(1.0, 100.0, 0.1, 7.0))
    assert not np.isclose(insitu["affinity_kJ_mol"].astype(float), placeholder).any()
    obs = st["chloride_observed"]
    post = validated.windows[(validated.windows["sensor"] == "orp") & (validated.windows["window"] == "post")]
    np.testing.assert_allclose(obs["eh_absolute"].astype(float), post.set_index("pulse").loc[obs["pulse"], "median"] / 1e3 + 0.2)
    assert not np.isclose(obs["eh_absolute"].astype(float), 0.2).any()  # 0.0 mV placeholder would give the offset
    co2 = st["co2"]
    pco2 = c.get("scenarios.pCO2_air_atm")
    np.testing.assert_allclose(co2["omega_co2_SCENARIO_T25C"].astype(float),
                               th.omega_co2(1.27, co2["pH"].to_numpy(), PK1_CARBONATE_25C, LOG10_KH_CO2_25C, pco2))
    at7 = float(th.omega_co2(1.27, 7.0, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, pco2))
    assert not np.isclose(co2.loc[co2["basis"] != "reference pH", "omega_co2_SCENARIO_T25C"].astype(float), at7).any()
