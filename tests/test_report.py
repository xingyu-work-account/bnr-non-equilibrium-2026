import pandas as pd
import pytest

from bnr_neq import report as R


def test_md_table_formats_and_escapes():
    df = pd.DataFrame({"a": [1.23456789, float("nan")], "b": ["x|y", True]})
    out = R.md_table(df, digits=3).splitlines()
    assert out[0] == "| a | b |" and out[1] == "| --- | --- |"
    assert out[2] == "| 1.23 | x\\|y |" and out[3] == "| NaN | True |"


def test_all_references_combines_milestones(tmp_path):
    pd.DataFrame({"quantity": ["q"], "pulse": [1], "reference": [1.0], "obtained": [1.1], "tolerance": [0.2],
                  "pass": [True], "note": [""]}).to_csv(tmp_path / "sensor_validation_references.csv", index=False)
    for name in ("segments_references.csv", "thermo_references.csv"):
        pd.DataFrame({"quantity": ["r"], "reference": [2.0], "obtained": [3.0], "tolerance": [0.5], "pass": [False],
                      "note": ["n"]}).to_csv(tmp_path / name, index=False)
    refs = R.all_references(tmp_path)
    assert refs["milestone"].tolist() == ["M2", "M3", "M4"]
    assert refs["quantity"].iloc[0] == "q (pulse 1)" and refs["pass"].tolist() == [True, False, False]


def _bound():
    from bnr_neq import alkalinity, volume

    rows = [{"alk0_meq_L": a, "alk0_basis": "SCAN 1-15 meq/L", "consistent_with_pH_record": True, "temperature_C": 25.0,
             "co2_scenario": alkalinity.AIR_EQ, "volume_scenario": volume.CONSTANT, "n_nit_max_end_mgN_L": 7.0 * a - 1,
             "n_nit_max_fraction_of_dose": 0.1, "n_removed_max_any_proton_sink_end_mgN_L": 14.0 * a,
             "tan_residual_min_end_mgN_L": 200.0, "fa_residual_min_end_mgNH3_L": 3.0,
             "tan_residual_min_any_proton_sink_end_mgN_L": 100.0, "fa_residual_min_any_proton_sink_end_mgNH3_L": 1.5,
             "dosed_total_mgN_L": 261.86} for a in (1.0, 5.0, 15.0)]
    return pd.DataFrame(rows)


def _tables():
    return {
        "thermo_coexistence.csv": pd.DataFrame({"basis": ["SCENARIO", "SCENARIO", "IN SITU"],
                                                "affinity_kJ_mol": ["305.5", "316.9", "not computed: DO_mg_L missing"],
                                                "DO_mg_L": [0.01, 0.1, float("nan")],
                                                "documented_scenario": [False, True, False]}),
        "thermo_affinity_steps.csv": pd.DataFrame({"delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T25C": [4.17, 2.36]}),
        "dynamics.csv": pd.DataFrame({"series": ["pH"], "n_stationary": [12], "n_windows": [187],
                                      "ordinal_kl_significant_windows": [10], "ordinal_kl_windows": [10]}),
        "dynamics_onsager.csv": pd.DataFrame({"snr_pulse_response": [1.8, 0.2], "snr_sufficient": [False, False]}),
        "pulses.csv": pd.DataFrame({
            "pulse": [1, 2, 1, 2, 1, 2],
            "metric": ["delta_A_kJ_mol_SCENARIO_T25C"] * 2 + ["delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T25C"] * 2
            + ["aeration_energy_demand_kWh_per_kgN_SCENARIO_SAE1", "aeration_energy_demand_kWh_per_kgN_SCENARIO_SAE4"],
            "value": [0.89, 0.61, 4.17, 2.36, 4.57, 1.1425], "status": ["computed"] * 6,
            "nh4_observation_model": ["M6_log", "M6_log", "identity_reading_as_TAN_SENSITIVITY",
                                      "identity_reading_as_TAN_SENSITIVITY", "M6_log", "M6_log"]}),
        "model_checks.csv": pd.DataFrame({"check": ["σ = r A / T ≥ 0 with the coexistence SCENARIO affinity"],
                                          "value": [0.0], "pass": [True]}),
        "nh4_truth_w2_consistency.csv": pd.DataFrame({"checkpoint": ["pre_pulse2", "end"], "consistent_at_best_fit": [True, False],
                                                      "alk0_required_if_nitrified_meq_L": [13.1, 20.9]}),
        "alkalinity_bound.csv": _bound(),
        "states_fractions.csv": pd.DataFrame({"interval": ["whole record"], "classified_fraction": [0.12],
                                              "idle_fraction_of_stirred_time": [0.03]}),
        "states_comparison.csv": pd.DataFrame({"comparison": ["pH rising vs not rising (within pulse intervals)"],
                                               "diff": [1.3], "ci95_low": [-0.5], "ci95_high": [2.6],
                                               "status": ["computed"]}),
        "volume_constraints.csv": pd.DataFrame({"constraint": ["N conservation over every pulse"],
                                                "value_mL_per_d": [0.31]}),
        "estimator_validation.csv": pd.DataFrame({"check": ["a", "b", "c"], "required": [True, True, False],
                                                  "pass": [True, True, False]}),
        "dynamics_windows.csv": pd.DataFrame({"series": ["pH", "pH"], "transform": ["raw", "raw"],
                                              "ordinal_kl_status": ["computed",
                                                                    "not computed: outside validated regime (x)"]}),
        "sensor_validation.csv": pd.DataFrame({"nh4_response_fraction": [0.37, 0.17], "ec_obs_change_mS_cm": [0.005, -0.06],
                                               "ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C": [0.7, 0.7],
                                               "orp_obs_change_mV": [0.4, 0.0]}),
        "model_identifiability.csv": pd.DataFrame({"identifiable": [True, False, False]}),
        "thermo_chloride_scenarios.csv": pd.DataFrame({"orp_pred_rise_mV_pulse1": [44.6, 7.6],
                                                       "orp_pred_rise_mV_pulse4": [7.0, 4.0]}),
    }


def test_status_labels_follow_the_prompt():
    labels = R.status_labels(_tables()).set_index("conclusion")["status"]
    # W7 / R7: the coexistence affinity, ΔA and σ ≥ 0 are consistency checks, not findings
    for c in ("Chemical disequilibrium of the batch (coexistence affinity)", "Affinity steps at the pulses",
              "Entropy production σ = r A / T ≥ 0 in the model"):
        assert labels[c] == R.CONSISTENCY
    assert "GO" not in {labels[c] for c in labels.index if "affinity" in c.lower() or "σ" in c}
    assert labels["pH trajectory estimators (supplement S1)"] == "GO WITH DOCUMENTED NON-BLOCKING DEBT"
    assert labels["Upper bound on cumulative nitrification from the pH record (W2)"] == "GO"
    assert labels["Single NH4 observation model against the W2 bound (S2)"] == "QUALIFIED"
    for c in ("Absolute NH4 concentrations, rates and conversion", "Conductivity as solution conductivity",
              "Absolute redox potential E_H", "Closed 100 mL frame (W8)"):
        assert labels[c] == "NO-GO"


def test_headline_affinity_step_comes_from_the_primary_model():
    """S1: the ΔA in the status table is the M6-inverted one (pulses.csv); the reading-as-TAN value is SENSITIVITY."""
    ev = R.status_labels(_tables()).set_index("conclusion").loc["Affinity steps at the pulses", "evidence"]
    assert ev.startswith("ΔA = 0.61 to 0.89 kJ/mol on the M6-inverted C(t)")
    assert "SENSITIVITY reading as TAN 2.36 to 4.17 kJ/mol" in ev
    t = _tables()
    t["thermo_affinity_steps.csv"] = t["thermo_affinity_steps.csv"] * 10  # the M4 table no longer feeds the headline
    assert R.status_labels(t).set_index("conclusion").loc["Affinity steps at the pulses", "evidence"] == ev


def test_no_go_rows_name_their_missing_items():
    lab = R.status_labels(_tables())
    no_go = lab[lab["status"] == "NO-GO"]
    assert no_go["missing items / debt"].str.startswith("Blocked by:").all()
    assert "orp_reference_offset_V" in no_go.set_index("conclusion").loc["Absolute redox potential E_H", "missing items / debt"]
    assert "4 to 45 mV" in no_go.set_index("conclusion").loc["Absolute redox potential E_H", "evidence"]  # R11 rounding


def test_assemble_report_requires_every_table(cfg, tmp_path):
    (tmp_path / "tables").mkdir()
    with pytest.raises(FileNotFoundError, match="data_properties.csv"):
        R.assemble_report(cfg, tmp_path)


def test_every_figure_and_table_indexed():
    assert len(R.FIGURES) == 8 and all(f.endswith(".png") for f, _ in R.FIGURES)
    assert {"fig1_linkage.png", "fig4_constraints.png", "fig8_efficiency.png"} <= {f for f, _ in R.FIGURES}  # W6
    assert "pulses.csv" in R.TABLE_DESCRIPTIONS and "model_fit.csv" in R.TABLE_DESCRIPTIONS


def test_evidence_is_derived_from_tables():
    t = _tables()
    lab = R.status_labels(t).set_index("conclusion")
    assert "17 to 37 %" in lab.loc["Absolute NH4 concentrations, rates and conversion", "evidence"]
    assert "2 of 2 required checks pass" in lab.loc["pH trajectory estimators (supplement S1)", "evidence"]
    assert "1 outside validated regimes" in lab.loc["pH trajectory estimators (supplement S1)", "evidence"]
    t["sensor_validation.csv"]["nh4_response_fraction"] = [0.5, 0.2]
    lab = R.status_labels(t).set_index("conclusion")
    assert "20 to 50 %" in lab.loc["Absolute NH4 concentrations, rates and conversion", "evidence"]
    assert "demand floor 1.14-4.57 kWh/kg N" in lab.loc["Oxygen demand and aeration-equivalent energy (W4)", "evidence"]
    assert "needs Alk0 ≥ 20.9 meq/L if nitrified; contradicted at end" in lab.loc[
        "Single NH4 observation model against the W2 bound (S2)", "evidence"]


def _copy(run_all_dir, tmp_path):
    import shutil

    out = tmp_path / "out"
    shutil.copytree(run_all_dir, out)
    return out


def test_assemble_report_on_real_outputs(cfg, run_all_dir, tmp_path):
    out = _copy(run_all_dir, tmp_path)
    path = R.assemble_report(cfg, out)
    text = path.read_text()
    v = pd.read_csv(out / "tables" / "sensor_verdicts.csv")
    for verdict in v["verdict"]:
        assert verdict in text
    assert "docs/progress.md" in text and "outputs/reports/progress.md" not in text
    # verdicts are read from the table, not written by hand
    v.loc[v["key"] == "ph", "verdict"] = "PASS (edited in the table)"
    v.to_csv(out / "tables" / "sensor_verdicts.csv", index=False)
    assert "PASS (edited in the table)" in R.assemble_report(cfg, out).read_text()


def test_report_raises_if_a_required_check_fails(cfg, run_all_dir, tmp_path):
    out = _copy(run_all_dir, tmp_path)
    (out / "reports" / "analysis_report.md").unlink()
    v = pd.read_csv(out / "tables" / "estimator_validation.csv")
    i = v.index[v["required"]][0]
    v.loc[i, "pass"] = False
    v.to_csv(out / "tables" / "estimator_validation.csv", index=False)
    with pytest.raises(R.ReportValidationError, match=v.loc[i, "check"]):
        R.assemble_report(cfg, out)
    assert not (out / "reports" / "analysis_report.md").exists()


def test_report_raises_on_a_matched_null_gate_violation(cfg, run_all_dir, tmp_path):
    out = _copy(run_all_dir, tmp_path)
    w = pd.read_csv(out / "tables" / "dynamics_windows.csv")
    i = w.index[(w["series"] == "pH") & (w["transform"] == "raw")][0]
    w.loc[i, "ac1"] = -0.5
    w["ordinal_kl_status"] = w["ordinal_kl_status"].astype(object)
    w.loc[i, "ordinal_kl_status"] = "computed"
    w.to_csv(out / "tables" / "dynamics_windows.csv", index=False)
    with pytest.raises(R.ReportValidationError, match="matched_null_gate_pH_ordinal_kl"):
        R.assemble_report(cfg, out)


def test_assumptions_read_the_configuration(cfg, run_all_dir):
    t = {"model_identifiability.csv": pd.read_csv(run_all_dir / "tables" / "model_identifiability.csv")}
    text = " ".join(R.assumptions(cfg, t))
    assert f"{cfg.get('sensors.nh4_to_N_factor'):g}" in text
    n_off = int(t["model_identifiability.csv"]["parameter"].str.startswith("offset_log_").sum())
    assert f"all {n_off} detected electrode-offset jumps" in text
