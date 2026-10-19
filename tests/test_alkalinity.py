"""W2 (proton, alkalinity and carbonate bound) and W8 (vessel volume scenarios)."""

import math

import numpy as np
import pytest

from bnr_neq import alkalinity as A
from bnr_neq import thermo, volume
from bnr_neq.config import CONFIG_PATH, ConfigError, load_config
from bnr_neq.constants import KELVIN_OFFSET, LOG10_KH_CO2_25C, PK1_CARBONATE_25C

T25 = 25.0 + KELVIN_OFFSET
PCO2 = 4.2e-4


def _cfg(tmp_path, measurements: str):
    """Config with the given TOML lines as [measurements] (keys removed from [missing])."""
    text = CONFIG_PATH.read_text()
    for line in measurements.strip().splitlines():
        k = line.split("=")[0].strip()
        text = text.replace(f'"{k}", ', "").replace(f', "{k}"', "")
    p = tmp_path / "experiment.toml"
    p.write_text(text + "\n[measurements]\n" + measurements.strip() + "\n")
    return load_config(p)


# --- positive -----------------------------------------------------------------------------
@pytest.mark.parametrize(("alk0", "expected"), [(2.0, 14.0), (5.0, 35.0), (10.0, 70.0)])
def test_bound_without_ph_floor(alk0, expected):
    assert A.nitrification_bound_mgN_L(alk0) == pytest.approx(expected, abs=0.05)


@pytest.mark.parametrize(("ph", "expected"), [(7.8, 0.401), (8.0, 0.636)])
def test_air_equilibrium_bicarbonate(ph, expected):
    assert A.hco3_air_eq_mmol_L(ph, T25, PCO2) == pytest.approx(expected, abs=0.001)
    # same value as the M4 constants at 25 °C
    assert thermo.hco3_air_eq_mmol_L(ph, PK1_CARBONATE_25C, LOG10_KH_CO2_25C, PCO2) == pytest.approx(expected, abs=0.001)


def test_free_ammonia_reference():
    assert thermo.free_ammonia_mgNH3_L(200.0, 7.8, T25) == pytest.approx(8.5, abs=0.1)


def test_carbonate_constants_temperature_dependence():
    k25, k20 = A.carbonate_constants(T25), A.carbonate_constants(20.0 + KELVIN_OFFSET)
    assert k25["pK1"] == pytest.approx(6.35, abs=0.01) and k25["log10_KH"] == pytest.approx(-1.47, abs=0.01)
    assert k25["pK2"] == pytest.approx(10.33, abs=0.01) and k25["pKw"] == pytest.approx(14.0, abs=0.01)
    assert k20["log10_KH"] > k25["log10_KH"]  # CO2 more soluble when colder
    assert k20["pKw"] == pytest.approx(14.17, abs=0.01)


def test_alkalinity_is_bicarbonate_at_neutral_ph():
    hco3 = A.hco3_air_eq_mmol_L(7.3, T25, PCO2)
    alk = A.alk_min_meq_L(7.3, T25, PCO2)
    assert alk == pytest.approx(hco3, rel=0.01)  # CO3, OH, H negligible at pH 7.3


def test_co2_excess_tightens_the_bound():
    air = float(A.nitrification_bound_mgN_L(10.0, A.alk_min_meq_L(7.8, T25, PCO2)))
    vals = [A.bound_with_co2(10.0, 7.8, T25, PCO2, 1.0, 40.0, k) for k in (1.0, 10.0, 100.0)]
    assert vals[0] < vals[1] < vals[2] < air
    assert vals[2] == pytest.approx(air, rel=0.02)


def test_bound_table_on_real_grid(cfg, ingested):
    b = A.bound_table(cfg, ingested.grid)
    assert set(b["alk0_meq_L"]) == set(A.ALK0_GRID_MEQ_L) and min(A.ALK0_GRID_MEQ_L) == 1 and max(A.ALK0_GRID_MEQ_L) == 15
    assert b["temperature_basis"].str.startswith("SCENARIO").all()
    assert set(b["volume_scenario"]) == {volume.scenario_label(e) for e in cfg.get("scenarios.evaporation_mL_per_d")}
    # invariants: the bound never exceeds the pH-floor-free value and never decreases in time
    cols = [c for c in b.columns if c.startswith("n_nit_max_") and c.endswith("_mgN_L")]
    v = b[cols].to_numpy()
    assert (np.diff(v, axis=1) >= -1e-9).all()
    no_floor = A.nitrification_bound_mgN_L(b["alk0_meq_L"].to_numpy(), 0.0, 1.0)
    assert (b.loc[b["volume_scenario"] == volume.CONSTANT, "n_nit_max_end_mgN_L"]
            <= no_floor[b["volume_scenario"] == volume.CONSTANT] + 1e-9).all()
    h = A.headline(b)
    assert len(h) == 15 and (h["n_nit_max_end_mgN_L"] < A.nitrification_bound_mgN_L(h["alk0_meq_L"].to_numpy())).all()
    # the dose is far beyond the bound at every Alk0 in the scan
    assert (b["n_nit_max_fraction_of_dose"] < 0.41).all()
    assert (b["tan_residual_min_end_mgN_L"] > 150).all()


def test_reference_checks_pass(cfg):
    r = A.reference_checks(cfg)
    assert r["pass"].all(), r[~r["pass"]]


# --- negative -----------------------------------------------------------------------------
def test_negative_initial_alkalinity_raises():
    with pytest.raises(ValueError, match="initial alkalinity"):
        A.nitrification_bound_mgN_L(-1.0)
    with pytest.raises(ValueError):
        A.nitrification_bound_mgN_L(5.0, alk_min_meq_L=-0.1)


def test_negative_measured_initial_alkalinity_raises(tmp_path):
    with pytest.raises(ConfigError, match="initial_alkalinity_meq_L"):
        _cfg(tmp_path, "initial_alkalinity_meq_L = -1.0")


# --- boundary -----------------------------------------------------------------------------
def test_alk0_at_the_ph_floor_gives_zero():
    floor = float(A.alk_min_meq_L(8.0, T25, PCO2))
    assert A.nitrification_bound_mgN_L(floor, floor) == 0.0
    assert A.nitrification_bound_mgN_L(0.5 * floor, floor) == 0.0  # below the floor: zero, never negative


def test_bound_curve_is_zero_before_the_ph_maximum_when_alk0_is_the_floor():
    t = np.arange(1.005, 10.0, 0.01)
    ph = np.where(t < 5.0, 7.2 + 0.1 * (t - 1.0), 7.3)  # peaks at 7.6 just before day 5
    scen = volume.VolumeScenario("c", "ASSUMPTION", 0.1, 0.0)
    floor = float(np.max(A.alk_min_meq_L(ph, T25, PCO2)))
    curve = A.bound_curve(t, ph, floor, T25, PCO2, scen)
    assert np.all(curve[t < 5.0] == 0.0) and curve[-1] > 0.0
    assert np.all(np.diff(curve) >= 0.0)


# --- W8: volume -----------------------------------------------------------------------------
def test_optional_keys_are_listed_as_missing(cfg):
    for key in ("vessel_mass_start_g", "vessel_mass_end_g", "covered", "top_up_mL", "initial_alkalinity_meq_L"):
        assert key in cfg.missing_keys and not cfg.has_measurement(key)


def test_default_scenarios_are_constant_plus_evaporation(cfg):
    s = volume.volume_scenarios(cfg, 47.92)
    assert s[0].is_constant and s[0].name == volume.CONSTANT
    assert [1e3 * x.net_loss_L_per_d for x in s] == cfg.get("scenarios.evaporation_mL_per_d")
    assert all(x.basis == "SCENARIO" for x in s[1:])


def test_measured_vessel_masses_give_one_measured_trajectory(tmp_path):
    c = _cfg(tmp_path, "vessel_mass_start_g = 150.0\nvessel_mass_end_g = 140.1\ntop_up_mL = 5.0")
    (s,) = volume.volume_scenarios(c, 50.0)
    assert s.basis == "MEASURED"
    # net loss 9.9 g plus 0.1 g of dosed salt inside the end mass = 10 mL over 50 d
    assert 1e3 * s.net_loss_L_per_d == pytest.approx(10.0 / 50.0)
    assert s.volume_L(50.0) == pytest.approx(0.090)


def test_covered_vessel_keeps_the_constant_frame(tmp_path):
    (s,) = volume.volume_scenarios(_cfg(tmp_path, "covered = true"), 47.92)
    assert s.is_constant and s.basis == "ASSUMPTION" and "covered" in s.name


@pytest.mark.parametrize("bad", ["vessel_mass_start_g = 150.0", "covered = 1", "top_up_mL = -2.0",
                                 "vessel_mass_start_g = -1.0\nvessel_mass_end_g = 140.0"])
def test_invalid_vessel_measurements_raise(tmp_path, bad):
    with pytest.raises(ConfigError):
        _cfg(tmp_path, bad)


def test_negative_evaporation_scenario_raises(tmp_path):
    text = CONFIG_PATH.read_text().replace("evaporation_mL_per_d = [0.0, 0.25, 0.5, 1.0]",
                                           "evaporation_mL_per_d = [0.0, -0.5]")
    p = tmp_path / "e.toml"
    p.write_text(text)
    with pytest.raises(ConfigError, match="evaporation_mL_per_d"):
        load_config(p)


def test_volume_that_empties_the_vessel_raises():
    s = volume.VolumeScenario("x", "SCENARIO", 0.1, 0.003)
    with pytest.raises(ValueError, match="empties"):
        s.volume_L(40.0)


def test_scenario_table_doses_and_bound(cfg):
    v = volume.scenario_table(cfg, 47.92, 53.23).set_index("volume_scenario")
    const = v.loc[volume.CONSTANT]
    assert const["dose_pulse1_mgN_L"] == pytest.approx(thermo.pulse_dose_from_config(cfg).n_mg_L)
    one = v.loc[volume.scenario_label(1.0)]
    assert one["V_at_pulse4_mL"] == pytest.approx(100.0 - 32.39)
    assert one["dose_pulse4_mgN_L"] == pytest.approx(const["dose_pulse4_mgN_L"] * 100.0 / (100.0 - 32.39))
    assert one["net_loss_upper_bound_mL_per_d"] == pytest.approx(100.0 / 53.23) and bool(one["within_bound"])


def test_balances_carry_the_volume_scenarios(cfg, validated):
    vols = volume.volume_scenarios(cfg, 47.92)
    nitrogen, chloride, alk = thermo.balance_tables(cfg, validated.balances, vols)
    assert (nitrogen.loc[nitrogen["volume_scenario"] == volume.CONSTANT, "term"].str.len() > 0).all()
    s1 = volume.scenario_label(1.0)
    cl = chloride[(chloride["volume_scenario"] == s1) & chloride["term"].str.contains("baseline 3 mM")]
    base = chloride[(chloride["volume_scenario"] == volume.CONSTANT) & chloride["term"].str.contains("baseline 3 mM")]
    f_end = 0.1 / (0.1 - 0.001 * thermo.NH4_BALANCE_END_D)
    assert float(cl["value"].iloc[0]) == pytest.approx(float(base["value"].iloc[0]) * f_end)
    added = nitrogen[(nitrogen["volume_scenario"] == s1) & nitrogen["term"].str.contains("concentration at the end")]
    assert float(added["value"].iloc[0]) == pytest.approx(4 * thermo.pulse_dose_from_config(cfg).n_mg_L * f_end)
    assert set(alk["volume_scenario"]) == {v.name for v in vols}


def test_alkalinity_balance_states_the_bound(cfg, validated, ingested):
    r = A.compute(cfg, ingested.grid)
    _, _, alk = thermo.balance_tables(cfg, validated.balances, volume.volume_scenarios(cfg, 47.92), r)
    rows = alk[alk["term"].str.startswith("Upper bound on cumulative nitrification")
               & (alk["volume_scenario"] == volume.CONSTANT)].set_index("term")["value"].astype(float)
    assert rows["Upper bound on cumulative nitrification by the end, Alk0 = 10 meq/L"] < 70.04
    h = A.headline(r.bound).set_index("alk0_meq_L")["n_nit_max_end_mgN_L"]
    assert rows["Upper bound on cumulative nitrification by the end, Alk0 = 2 meq/L"] == pytest.approx(h[2.0])
    assert "not computed: initial_alkalinity_meq_L missing" in " ".join(alk["value"].astype(str))


def test_measured_initial_alkalinity_is_evaluated(tmp_path, ingested):
    c = _cfg(tmp_path, "initial_alkalinity_meq_L = 4.2")
    b = A.bound_table(c, ingested.grid)
    m = b[b["alk0_basis"] == "MEASURED"]
    assert len(m) > 0 and (m["alk0_meq_L"] == 4.2).all()
    scan4 = b[(b["alk0_meq_L"] == 4.0) & (b["alk0_basis"] != "MEASURED")]["n_nit_max_end_mgN_L"].to_numpy()
    assert (m["n_nit_max_end_mgN_L"].to_numpy() > scan4).all()


def test_math_of_the_any_sink_bound():
    assert A.nitrification_bound_mgN_L(3.0, 0.0, 1.0, A.EQ_PER_MOL_N_ANY_SINK) == pytest.approx(
        2 * A.nitrification_bound_mgN_L(3.0))
    assert math.isclose(float(A.co2_excess_mM(14.007, 10.0, 2.0)), 2 * 1.0 / 10.0 / 2.0)
