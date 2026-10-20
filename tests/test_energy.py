"""W4: oxygen and aeration-equivalent energy metric family."""

import numpy as np
import pandas as pd
import pytest

from bnr_neq import energy as E
from bnr_neq.config import CONFIG_PATH, ConfigError, load_config
from bnr_neq.constants import DO_SAT_25C_MG_L


def test_hand_values():
    assert E.oxygen_required_mgO2(10.0) == pytest.approx(45.7)
    assert E.oxygen_required_mgO2(10.0, nitritation=True) == pytest.approx(34.3)
    # 10 1/d x 8.26 mg/L x 0.1 L x 2 d = 16.52 mg O2
    assert E.o2_supply_capacity_mgO2(10.0, 8.26, 0.0, 0.1, 2.0) == pytest.approx(16.52)
    assert E.o2_supply_capacity_mgO2(10.0, 8.26, 2.0, 0.1, 2.0) == pytest.approx(10 * 6.26 * 0.2)
    assert E.utilization(8.26, 16.52) == pytest.approx(0.5)
    # 4.57 kg O2 per kg N at 2 kg O2/kWh = 2.285 kWh per kg N
    e = E.aeration_energy_kWh(4.57e6, 2.0)
    assert e == pytest.approx(2.285) and E.specific_energy_kWh_per_kgN(e, 1e6) == pytest.approx(2.285)
    assert E.gibbs_released_J(14.007) == pytest.approx(348.8)  # 1 mmol N x 348.8 kJ/mol
    assert E.kwh_to_J(1.0) == 3.6e6


@pytest.mark.parametrize("sae", [0.0, -1.0, np.nan])
def test_negative_or_zero_sae_raises(sae):
    with pytest.raises(ValueError, match="SAE"):
        E.aeration_energy_kWh(1.0, sae)


def test_negative_sae_in_the_config_raises(tmp_path):
    lines = CONFIG_PATH.read_text().splitlines()
    assert sum(x.startswith("sae_field_kgO2_per_kWh =") for x in lines) == 1
    text = "\n".join("sae_field_kgO2_per_kWh = [-1.0, 2.0]" if x.startswith("sae_field_kgO2_per_kWh =") else x for x in lines)
    p = tmp_path / "c.toml"
    p.write_text(text)
    with pytest.raises(ConfigError, match="sae_field_kgO2_per_kWh"):
        load_config(p)


def test_sae_scenario_needs_a_source(tmp_path):
    lines = [x for x in CONFIG_PATH.read_text().splitlines() if not x.startswith("sae_source")]
    p = tmp_path / "c.toml"
    p.write_text("\n".join(lines))
    with pytest.raises(ConfigError, match="sae_source"):
        load_config(p)


def test_do_above_saturation_and_bad_inputs_raise():
    with pytest.raises(ValueError, match="saturation"):
        E.o2_supply_capacity_mgO2(10.0, 8.26, 9.0, 0.1, 1.0)
    with pytest.raises(ValueError):
        E.o2_supply_capacity_mgO2(-1.0, 8.26, 0.0, 0.1, 1.0)
    with pytest.raises(ValueError):
        E.oxygen_required_mgO2(-1.0)


def test_energy_family_on_the_real_record(cfg, run_all_dir):
    """Hand recomputation of one aeration-equivalent energy from its stated SCENARIOS, and the feasibility guards."""
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    m = pm.set_index(["pulse", "metric"])
    t1, t2 = cfg.pulse_times_d[0], cfg.pulse_times_d[1]
    n = m.loc[(1, "n_removed_mgN"), "value"]
    cap = 10.0 * DO_SAT_25C_MG_L * cfg.get("vessel.volume_L") * (t2 - t1)
    assert m.loc[(1, "o2_supply_capacity_mgO2_SCENARIO_kLa10_T25C"), "value"] == pytest.approx(cap)
    sae = sorted(cfg.get("scenarios.sae_field_kgO2_per_kWh"))[1]  # the central SCENARIO value
    kwh_per_kg = cap / 1e6 / sae / (n / 1e6)
    assert m.loc[(1, f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{sae:g}_kLa10_T25C"), "value"] == pytest.approx(kwh_per_kg)
    assert m.loc[(1, f"aeration_energy_demand_kWh_per_kgN_SCENARIO_AE{sae:g}"), "value"] == pytest.approx(4.57 / sae)
    assert m.loc[(1, "o2_required_mgO2"), "value"] == pytest.approx(4.57 * n)
    # every scenario value carries its label; no computed utilization exceeds 1 (an infeasible kLa gets a status)
    fam = pm[pm["metric"].str.contains("AE|kLa")]
    assert fam["metric"].str.contains("SCENARIO").all()
    u = pm[pm["metric"].str.startswith("o2_utilization")]
    assert (u.loc[u["status"] == "computed", "value"] <= 1).all()
    infeasible = u[u["status"].str.startswith("not computed: SCENARIO infeasible")]
    for k, metric in zip(infeasible["pulse"], infeasible["metric"]):
        kl = metric.split("SCENARIO_")[1].removesuffix("_T25C")
        sup = pm[(pm["pulse"] == k) & pm["metric"].str.startswith("aeration_energy_supply") & pm["metric"].str.contains(f"_{kl}_")]
        assert (sup["status"] != "computed").all()
    # the stirrer stays a laboratory-specific term
    assert m.loc[(1, "specific_stirring_energy_J_per_gN"), "status"] == "not computed: stirrer_power_W missing"
