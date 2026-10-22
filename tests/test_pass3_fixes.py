"""Remediation pass 3: converged model comparison (X1), hardened provenance checker (X3), names (X4), η band and W2
labels (X5), R8 scope (X6) and the real-record integration checks."""

import numpy as np
import pandas as pd
import pytest

from bnr_neq import model as M
from bnr_neq import pulses as P
from bnr_neq import segments as S
from bnr_neq import thermo


def _rows(*specs):
    return pd.DataFrame([{"metric": m, "basis": b, "nh4_observation_model": o, "nh4_series": s} for m, b, o, s in specs])


PRIMARY = ("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL, "C_inv:aaa")


# --- X3 (a): every row declares its NH4 basis --------------------------------------------------------------------------
def test_empty_model_on_an_nh4_derived_metric_raises():
    for metric in ("fa_median_mgNH3_L_SCENARIO_T25C", "idle_fraction_LOWER_BOUND", "t_pulse_d"):
        with pytest.raises(M.ObservationModelError, match="without a declared NH4 basis"):
            M.check_single_observation_model(_rows(PRIMARY, (metric, "", "", "")))


def test_unknown_model_and_series_on_an_independent_row_raise():
    with pytest.raises(M.ObservationModelError, match="undeclared NH4 basis"):
        M.check_single_observation_model(_rows(("o2_demand_mgO2_L_d", "M7", "M7_linear", "C_inv:aaa")))
    with pytest.raises(M.ObservationModelError, match="NH4-independent rows carry an NH4 series"):
        M.check_single_observation_model(_rows(PRIMARY, ("t_pulse_d", "fact", M.NH4_INDEPENDENT, "C_inv:aaa")))
    ok = _rows(PRIMARY, ("t_pulse_d", "fact", M.NH4_INDEPENDENT, ""), ("median_pH", "measured", M.NH4_INDEPENDENT, ""))
    assert M.check_single_observation_model(ok) == M.OBSERVATION_MODEL


# --- X3 (b): an allow-list of measured raw metrics -----------------------------------------------------------------------
@pytest.mark.parametrize("metric", ["aeration_energy_supply_kWh_per_kgN_SCENARIO_AE1_kLa10_T25C",
                                    "aeration_energy_demand_kWh_per_kgN_SCENARIO_AE1", "n_removed_mgN",
                                    "removal_to_dose_ratio", "consumption_weighted_mgN_L_d", "consumption_first_window_mgN_L_d",
                                    "o2_required_mgO2", "gibbs_to_aeration_energy_demand_SCENARIO_AE1_T25C",
                                    "idle_fraction_UPPER_BOUND", "dose_removal_time_d", "alk0_required_if_nitrified_meq_L"])
def test_derived_metric_labelled_as_the_raw_reading_raises(metric):
    with pytest.raises(M.ObservationModelError, match="raw electrode reading"):
        M.check_single_observation_model(_rows(PRIMARY, (metric, "reading", M.RAW_READING, M.RAW_READING)))


def test_measured_raw_metrics_pass():
    rows = _rows(PRIMARY, *[(m, "measured", M.RAW_READING, M.RAW_READING) for m in M.RAW_MEASURED_METRICS])
    assert M.check_single_observation_model(rows) == M.OBSERVATION_MODEL


def test_reading_based_metric_names_say_sensitivity():
    with pytest.raises(M.ObservationModelError, match="without SENSITIVITY in their names"):
        M.check_single_observation_model(_rows(PRIMARY, ("consumption_gain_corrected_weighted_mgN_L_d",
                                                         "SENSITIVITY: GAIN-CORRECTED", M.LINEAR_GAIN, M.RAW_READING)))


# --- X3 (b), wide tables: an unlabelled gain-corrected column raises ---------------------------------------------------
def _wide(gain_col):
    return pd.DataFrame({"sub_segment": [1, 2], "status": ["computed"] * 2, "slope_mgNH4_L_d": [-1.0, -2.0],
                         gain_col: [3.0, 4.0], "nh4_observation_model": M.RAW_READING, "nh4_series": M.RAW_READING})


def test_unlabelled_gain_corrected_column_raises():
    with pytest.raises(M.ObservationModelError, match="reading-based columns without SENSITIVITY"):
        M.check_single_observation_model(_wide("consumption_gain_corrected_mgN_L_d"), label_cols=(),
                                         independent_cols=("sub_segment", "status"))
    with pytest.raises(M.ObservationModelError, match="declare its NH4-independent columns"):
        M.check_single_observation_model(_wide(S.GAIN_COL), label_cols=())
    assert M.check_single_observation_model(_wide(S.GAIN_COL), label_cols=(),
                                            independent_cols=("sub_segment", "status")) is None


# --- X3 (c): the MODEL-BASED label is a token ------------------------------------------------------------------------
def test_not_model_based_in_a_basis_does_not_raise():
    rows = _rows(("consumption_weighted_mgN_L_d", "M6 inverted C(t), not MODEL-BASED", M.OBSERVATION_MODEL, "C_inv:aaa"),
                 ("eta_conversion_fraction_MODEL_BASED", "MODEL-BASED", M.OBSERVATION_MODEL, "M6_state:x"))
    assert M.check_single_observation_model(rows) == M.OBSERVATION_MODEL
    with pytest.raises(M.ObservationModelError, match="MODEL-BASED carry an inverted"):  # the token still counts
        M.check_single_observation_model(_rows(("x_MODEL_BASED", "", M.OBSERVATION_MODEL, "C_inv:aaa")))


# --- X5 ------------------------------------------------------------------------------------------------------------------
def test_eta_band_edge_above_one_is_clipped():
    band, note = P.eta_band_on_unit_interval((0.618, 1.069))
    assert band == (0.618, 1.0) and "clipped to [0, 1]" in note and "1.069" in note
    assert P.eta_band_on_unit_interval((-0.1, 0.5))[0] == (0.0, 0.5)
    assert P.eta_band_on_unit_interval((0.2, 0.9)) == ((0.2, 0.9), "")


def test_baseline_windows_are_labelled_from_pre_pulse1():
    cons = pd.DataFrame({"checkpoint": ["pre_pulse1", "pre_pulse2", "end"], "consistent_at_best_fit": [True, False, False],
                         "verdict": ["consistent: needs 9.5", "QUALIFIED: contradicts", "QUALIFIED: contradicts"]})
    tw = pd.DataFrame({"governing_pulse": pd.array([pd.NA, pd.NA, 1, 2], dtype="Int64"),
                       "nh4_observation_model": M.OBSERVATION_MODEL})
    q = P.qualify(tw, cons, 2, "governing_pulse", "nh4_observation_model")
    assert (q["w2_bound_check"].iloc[:2] == "pre_pulse1: consistent: needs 9.5").all()
    assert q["w2_bound_check"].iloc[2].startswith("QUALIFIED") and q["w2_bound_check"].iloc[3].startswith("QUALIFIED")
    bad = cons.assign(consistent_at_best_fit=False, verdict="QUALIFIED: contradicts")
    with pytest.raises(P.BoundContradictionError):  # a contradicted baseline must be qualified too
        P.check_w2_qualified(q, bad, 2, "governing_pulse")


# --- X6 ------------------------------------------------------------------------------------------------------------------
def test_r8_case_on_the_synthetic_truth(cfg, segmented, ingested):
    from conftest import synthetic_truth

    truth, vectors, _ = synthetic_truth(ingested.grid)
    tw = S.truth_slope_table(cfg, truth, segmented.table, P.band_series(truth, vectors))
    frag = tw[tw["status"].str.startswith(S.FRAGMENT_STATUS)]
    assert frag["consumption_mgN_L_d"].isna().all() and frag["consumption_fragment_diagnostic_mgN_L_d"].notna().all()
    assert tw.loc[tw["status"] == "computed", "consumption_fragment_diagnostic_mgN_L_d"].isna().all()


# --- real record (run-all outputs) ---------------------------------------------------------------------------------------
def test_every_checked_table_passes_the_hardened_checker(run_all_dir):
    tab = run_all_dir / "tables"
    pm = pd.read_csv(tab / "pulses.csv")
    assert M.check_single_observation_model(pm) == M.OBSERVATION_MODEL
    # pass 4 (Y4): a wide table with primary rows declares its NH4-independent and C(t)-derived columns
    from bnr_neq import states as ST
    from bnr_neq import volume
    from bnr_neq.config import load_config

    cfg = load_config()
    indep, derived = S.truth_window_columns(cfg, volume.volume_scenarios(cfg, P.END_D))
    assert M.check_single_observation_model(pd.read_csv(tab / "nh4_truth_windows.csv"),
                                            label_cols=("status", "w2_bound_check"), series_col="nh4_truth_id",
                                            independent_cols=indep, derived_cols=derived) == M.OBSERVATION_MODEL
    for name, key in (("states_fractions.csv", "fractions"), ("states_timeline.csv", "timeline"),
                      ("states_comparison.csv", "comparisons")):
        indep, derived = ST.STATE_TABLE_COLUMNS[key]
        assert M.check_single_observation_model(pd.read_csv(tab / name), label_cols=(), independent_cols=indep,
                                                derived_cols=derived) == M.OBSERVATION_MODEL
    assert M.check_single_observation_model(pd.read_csv(tab / "segments.csv"), label_cols=(),
                                            independent_cols=S.SEGMENTS_NH4_INDEPENDENT) is None
    assert M.check_single_observation_model(pd.read_csv(tab / "thermo_oxygen.csv"), label_cols=("role",),
                                            independent_cols=thermo.OXYGEN_NH4_INDEPENDENT) is None
    assert (pm["nh4_observation_model"].fillna("") != "").all()


def test_real_eta_bands_and_baseline_labels(run_all_dir):
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    eta = pm[pm["metric"].str.startswith("eta_conversion_fraction") & (pm["status"] == "computed")]
    assert ((eta["param_band_low"].fillna(0) >= 0) & (eta["param_band_high"].fillna(0) <= 1)).all()
    tw = pd.read_csv(run_all_dir / "tables" / "nh4_truth_windows.csv")
    base = tw[tw["governing_pulse"].isna()]
    assert len(base) and base["w2_bound_check"].fillna("").str.len().gt(0).all()
    assert base["w2_bound_check"].str.startswith(("pre_pulse1: ", "QUALIFIED")).all()


def test_real_r8_case_and_renamed_metrics(run_all_dir):
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    seg = pd.read_csv(run_all_dir / "tables" / "segments.csv").set_index("sub_segment")
    assert seg.loc[14, "status"].startswith(S.FRAGMENT_STATUS)
    case = pm[pm["metric"] == "consumption_weighted_mgN_L_d_R8_CASE_incl_fragments"].set_index("pulse")
    assert case.loc[1, "status"] == "computed" and "14 (" in case.loc[1, "basis"]
    assert case.loc[1, "basis"].startswith("R8 CASE (not the primary rate)")
    prim = pm[(pm["metric"] == "consumption_weighted_mgN_L_d") & (pm["pulse"] == 1)].iloc[0]
    assert "14" not in prim["inputs"].split("(sub-segments ")[1].split(")")[0].split(", ")
    assert not pm["metric"].str.contains("relaxation_time|_SAE|gain_corrected").any()
    assert {"dose_removal_time_d", "eta_gap_inversion_misfit_MODEL_BASED"} <= set(pm["metric"])


# --- X1: only converged candidates can be selected ----------------------------------------------------------------------
def _comparison(converged, sel_rmse):
    names = ["null_no_consumption", "biomass_decay", "alkalinity_limited"]
    return pd.DataFrame({"model": names, "cal_converged": converged, "cal_status": [2 if c else 0 for c in converged],
                         "sel_rmse_log": sel_rmse, "sel_skill_vs_null": [0.0, 0.9, 0.95]})


def test_unconverged_candidate_cannot_be_selected():
    out = M.select_converged(_comparison([True, True, False], [0.3, 0.07, 0.05])).set_index("model")
    assert out.loc["biomass_decay", "selected"] and not out.loc["alkalinity_limited", "selected"]
    assert out.loc["alkalinity_limited", "selection_status"] == (
        "not computed: calibration did not converge (status 0: maximum number of function evaluations exceeded); "
        "skills not comparable")
    assert out.loc["alkalinity_limited", "sel_skill_vs_null"] == 0.95  # its skills stay in the table, labelled
    assert out.loc["null_no_consumption", "selection_status"] == "reference (not a candidate)"
    assert out.loc["biomass_decay", "selection_status"] == "selected"
    # once it converges, the better selection skill wins
    out = M.select_converged(_comparison([True, True, True], [0.3, 0.07, 0.05])).set_index("model")
    assert out.loc["alkalinity_limited", "selected"] and out["selected"].sum() == 1
    assert out.loc["biomass_decay", "selection_status"] == "candidate (converged, not selected)"


def test_no_converged_candidate_raises():
    with pytest.raises(M.SelectionError):
        M.select_converged(_comparison([True, False, False], [0.3, 0.07, 0.05]))


def test_real_selected_model_converged_and_alk0_reported(run_all_dir):
    fit = pd.read_csv(run_all_dir / "tables" / "model_fit.csv")
    assert fit["selected"].sum() == 1 and fit.loc[fit["selected"], "cal_converged"].all()
    assert (fit.loc[~fit["cal_converged"].astype(bool), "selection_status"].str.startswith(
        "not computed: calibration did not converge")).all()
    alk = fit[fit["model"].str.startswith("alkalinity")]
    assert len(alk) == 3 and alk["alk0_fitted_meq_L"].notna().all()
    assert alk["alk0_profile_status"].str.startswith(("identifiable", "not identifiable", "not computed")).all()
    ok = alk["alk0_profile_status"] == "identifiable"
    assert (alk.loc[ok, "alk0_ci95_low_meq_L"] <= alk.loc[ok, "alk0_fitted_meq_L"]).all()
    assert (alk.loc[ok, "alk0_fitted_meq_L"] <= alk.loc[ok, "alk0_ci95_high_meq_L"]).all()
    assert fit["cal_start"].notna().all() and (fit["cal_nfev"] > 0).all()


# --- X1 positive: on a synthetic record of known parameters every candidate converges and Alk0 is recovered -------------
T_SYN = np.round(np.arange(1.01, 16.0, 0.02), 10)
PH_SYN = 7.6 - 0.02 * (T_SYN - 1)
OBS_SYN = {"a": 0.1, "b0": 0.8, "b1_per_d": 0.0, "tau_d": 0.05, "offset_values": []}
ALK0_TRUE = 15.0  # meq/L: the alkalinity runs out (105 mg N/L consumed) before the last of three 60 mg N/L pulses


def _synthetic(spec, kin):
    sim = M.simulate(spec, kin, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0)
    assert sim.consumed == pytest.approx(ALK0_TRUE / 2 * 14.007, rel=1e-3)  # the bound is reached: Alk0 identifiable
    return M.ModelData(t=T_SYN, log_r=M.observe(sim.c, T_SYN, OBS_SYN, []), ph=PH_SYN, pulses=[4.0, 8.0, 12.0],
                       dose=60.0, offsets=[], offset_guess=[])


def _calibrate(specs, data):
    setup = {}
    for s in specs:
        lay = M.layout(s, data)
        x0 = M.initial_guess(s, data, lay)
        setup[s.name] = (lay, x0, np.ones(len(x0), dtype=bool))
    return M.calibrate_all(specs, data, (1.0, 16.0), setup)[0], setup


def test_every_candidate_converges_under_the_policy_and_alk0_is_recovered(cfg):
    specs = {s.name: s for s in M.model_specs(cfg)}
    data = _synthetic(specs["alkalinity_limited"], {"v0_mgN_L_d": 12.0, "Alk0_meq_L": ALK0_TRUE})
    fits, setup = _calibrate(list(specs.values()), data)
    for name, f in fits.items():
        assert f.converged and f.status > 0 and f.nfev <= len(M.CAL_START_FACTORS) * M.CAL_MAX_NFEV, name
    lay = setup["alkalinity_limited"][0]
    alk0 = M.unpack(specs["alkalinity_limited"], data, lay, fits["alkalinity_limited"].x)["Alk0_meq_L"]
    assert alk0 == pytest.approx(ALK0_TRUE, rel=0.01)  # stated tolerance: 1 %
    # the FA variant recovers Alk0 on a record generated by its own rate law
    fa = specs["alkalinity_FA_SCENARIO_T25C"]
    data_fa = _synthetic(fa, {"vmax_mgN_L_d": 15.0, "Ks_mgNH3_L": 0.5, "Ki_mgNH3_L": 8.0, "Alk0_meq_L": ALK0_TRUE})
    f, setup = _calibrate([fa], data_fa)
    assert f[fa.name].converged
    assert M.unpack(fa, data_fa, setup[fa.name][0], f[fa.name].x)["Alk0_meq_L"] == pytest.approx(ALK0_TRUE, rel=0.01)


def test_policy_starts_are_the_same_rule_for_every_candidate(cfg):
    data = _synthetic(M.ModelSpec("alkalinity_limited", ("v0_mgN_L_d", "Alk0_meq_L"), False, False, alk=True),
                      {"v0_mgN_L_d": 12.0, "Alk0_meq_L": ALK0_TRUE})
    for spec in M.model_specs(cfg):
        lay = M.layout(spec, data)
        x0 = M.initial_guess(spec, data, lay)
        starts = M.calibration_starts(spec, lay, x0)
        assert len(starts) == (len(M.CAL_START_FACTORS) if spec.kinetic else 1)
        for _, xs in starts:  # only the kinetic block moves; the observation start is shared
            np.testing.assert_array_equal(xs[len(spec.kinetic):], x0[len(spec.kinetic):])
        if spec.alk:
            k = lay.names.index("Alk0_meq_L")
            assert [round(float(np.exp(xs[k])), 6) for _, xs in starts] == [5.0, 10.0, 20.0, 40.0]


def test_real_calibration_starts_follow_the_policy(run_all_dir):
    st = pd.read_csv(run_all_dir / "tables" / "model_calibration_starts.csv")
    fit = pd.read_csv(run_all_dir / "tables" / "model_fit.csv").set_index("model")
    for model, g in st.groupby("model"):
        assert len(g) == (1 if model == M.NULL_MODEL else len(M.CAL_START_FACTORS))
        assert g["chosen"].sum() == 1 and (g["nfev"] <= M.CAL_MAX_NFEV).all()
        chosen = g[g["chosen"]].iloc[0]
        pool = g[g["converged"]] if g["converged"].any() else g
        assert chosen["cal_rmse_log"] == pytest.approx(pool["cal_rmse_log"].min())  # the lowest-RSS converged start
        assert bool(chosen["converged"]) == bool(fit.loc[model, "cal_converged"])
        assert fit.loc[model, "cal_nfev"] == g["nfev"].sum()


def test_alk0_interval_is_withheld_when_another_minimum_lies_inside_the_threshold(cfg):
    from dataclasses import replace

    spec = {s.name: s for s in M.model_specs(cfg)}["alkalinity_limited"]
    rng = np.random.default_rng(3)
    data = _synthetic(spec, {"v0_mgN_L_d": 12.0, "Alk0_meq_L": ALK0_TRUE})
    data = M.ModelData(t=data.t, log_r=data.log_r + 0.02 * rng.normal(size=len(data.t)), ph=data.ph, pulses=data.pulses,
                       dose=data.dose, offsets=[], offset_guess=[])
    fits, setup = _calibrate([spec], data)
    cal, lay = fits[spec.name], setup[spec.name][0]
    alone = M.alk0_profile(spec, data, lay, cal, [cal])
    assert alone["alk0_profile_status"] == "identifiable"
    assert alone["alk0_ci95_low_meq_L"] < ALK0_TRUE < alone["alk0_ci95_high_meq_L"]
    k = lay.names.index("Alk0_meq_L")
    far = cal.x.copy()
    far[k] += np.log(5.0)  # a second converged minimum at 5 x Alk0 ...
    inside = replace(cal, x=far, rss=cal.rss * 1.001)  # ... with Δ far below 3.84
    r = M.alk0_profile(spec, data, lay, cal, [cal, inside])
    assert r["alk0_profile_status"].startswith("not identifiable: the profile has more than one minimum below Δ = 3.84")
    assert np.isnan(r["alk0_ci95_low_meq_L"]) and np.isnan(r["alk0_ci95_high_meq_L"])
    outside = replace(cal, x=far, rss=cal.rss * 10.0)  # Δ far above 3.84: the interval stands
    assert M.alk0_profile(spec, data, lay, cal, [cal, outside])["alk0_profile_status"] == "identifiable"
