"""Remediation pass 4: skill range over indistinguishable starts (Y1), Alk0 evidence and distinct minima (Y2), the
limits of the idle bounds (Y3), the NH4-independent allow-list and declared wide-table columns (Y4), and the real-record
integration checks."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bnr_neq import model as M
from bnr_neq import report as R
from bnr_neq import states as ST

ROOT = Path(__file__).resolve().parents[1]
PASS3_COMMIT = "ad51a6a"  # outputs of pass 3, accepted by docs/rereview_2026-10-20_pass3.md
# Z5 (pass 5): the reference values of PASS3_COMMIT, extracted once with `git show` (command in docs/progress.md)
PASS3_FIXTURES = ROOT / "tests" / "fixtures" / "pass3"

T_SYN = np.round(np.arange(1.01, 16.0, 0.02), 10)
PH_SYN = 7.6 - 0.02 * (T_SYN - 1)
OBS_SYN = {"a": 0.1, "b0": 0.8, "b1_per_d": 0.0, "tau_d": 0.05, "offset_values": []}


def _synthetic_data(spec, kin, noise=0.02, seed=3):
    sim = M.simulate(spec, kin, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0)
    lr = M.observe(sim.c, T_SYN, OBS_SYN, []) + noise * np.random.default_rng(seed).normal(size=len(T_SYN))
    return M.ModelData(t=T_SYN, log_r=lr, ph=PH_SYN, pulses=[4.0, 8.0, 12.0], dose=60.0, offsets=[], offset_guess=[])


# --- Y1: the skill range spans every converged start inside Δ < 3.84 -------------------------------------------------
def _starts(rows):
    return pd.DataFrame([{"model": "m", "converged": c, "inside_95": i, "sel_skill_vs_null_diagnostic": s,
                          "test_skill_vs_null_diagnostic": t, "forecast_skill_vs_null_diagnostic": f}
                         for c, i, s, t, f in rows])


def test_skill_range_spans_both_minima_inside_the_region():
    st = _starts([(True, True, -0.9, 0.23, 0.44), (True, True, 0.93, 0.82, 0.93),  # two minima inside
                  (True, False, -1.2, -0.5, 0.0), (False, False, 0.99, 0.99, 0.99)])  # outside; unconverged
    r = M.skill_range(st)
    assert r["n_starts_inside_95"] == 2 and r["skill_range_status"].startswith("2 converged starts inside")
    assert (r["sel_skill_min_inside_95"], r["sel_skill_max_inside_95"]) == (-0.9, 0.93)
    assert (r["test_skill_min_inside_95"], r["test_skill_max_inside_95"]) == (0.23, 0.82)
    assert (r["forecast_skill_min_inside_95"], r["forecast_skill_max_inside_95"]) == (0.44, 0.93)
    one = M.skill_range(st.iloc[[0, 2]])
    assert one["n_starts_inside_95"] == 1 and one["skill_range_status"].startswith("one converged start")
    none = M.skill_range(st.iloc[[3]].assign(converged=False))
    assert none["skill_range_status"].startswith("not computed") and np.isnan(none["sel_skill_min_inside_95"])


def test_two_minima_on_a_synthetic_record_give_a_range_that_spans_both(cfg):
    """On a synthetic record, a second calibration point inside Δ < 3.84 is inside the region and its out-of-sample
    skills enter the range; one far above Δ = 3.84 does not."""
    spec = {s.name: s for s in M.model_specs(cfg)}["alkalinity_limited"]
    data = _synthetic_data(spec, {"v0_mgN_L_d": 12.0, "Alk0_meq_L": 15.0})
    lay = M.layout(spec, data)
    x0 = M.initial_guess(spec, data, lay)
    cal = M.fit(spec, data, lay, (1.0, 8.0), x0, np.ones(len(x0), dtype=bool), max_nfev=M.CAL_MAX_NFEV)
    assert cal.converged
    n_eff = M.calibration_n_eff(spec, data, lay, cal)
    k = lay.names.index("Alk0_meq_L")
    other_x = cal.x.copy()
    other_x[k] += np.log(3.0)
    inside = replace(cal, x=other_x, rss=cal.rss * np.exp(2.0 / n_eff))  # Δ = 2 < 3.84
    outside = replace(cal, x=other_x, rss=cal.rss * np.exp(10.0 / n_eff))  # Δ = 10
    assert M.start_delta(inside, cal, n_eff) == pytest.approx(2.0)
    rows = []
    for f in (cal, inside, outside):
        d = M.start_delta(f, cal, n_eff)
        rows.append({"model": spec.name, "converged": True, "inside_95": d < M.CHI2_95_1DOF,
                     "sel_rmse_log_diagnostic": M.window_skill(spec, data, lay, f.x, (8.0, 12.0))["rmse_log"],
                     **M._start_skill_job((spec, data, lay, f.x, (12.0, 16.0)))})
    st = pd.DataFrame(rows)
    for c in ("sel", "test", "forecast"):  # skills against an arbitrary reference RMSE of 1
        st[f"{c}_skill_vs_null_diagnostic"] = 1 - st[f"{c}_rmse_log_diagnostic"] ** 2
    r = M.skill_range(st)
    assert r["n_starts_inside_95"] == 2
    for c in ("sel", "test", "forecast"):
        v = st[f"{c}_skill_vs_null_diagnostic"].iloc[:2]
        assert r[f"{c}_skill_min_inside_95"] == v.min() and r[f"{c}_skill_max_inside_95"] == v.max()
        assert v.iloc[0] != v.iloc[1]  # the two points forecast differently, and the range spans both


def _fit_rows(alk_rival_test=(0.227, 0.817), alk_converged=True):
    return pd.DataFrame([
        {"model": "null_no_consumption", "cal_converged": True, "selected": False, "sel_rmse_log": 0.3,
         "sel_skill_min_inside_95": 0.0, "sel_skill_max_inside_95": 0.0},
        {"model": "biomass_decay", "cal_converged": True, "selected": True, "sel_rmse_log": 0.07,
         "sel_skill_min_inside_95": 0.9577, "sel_skill_max_inside_95": 0.9578, "test_skill_min_inside_95": 0.748,
         "test_skill_max_inside_95": 0.748, "forecast_skill_min_inside_95": 0.862, "forecast_skill_max_inside_95": 0.863,
         "sel_skill_vs_null": 0.958},
        {"model": "alkalinity_limited", "cal_converged": alk_converged, "selected": False, "sel_rmse_log": 0.48,
         "sel_skill_min_inside_95": -0.908, "sel_skill_max_inside_95": 0.935, "test_skill_min_inside_95": alk_rival_test[0],
         "test_skill_max_inside_95": alk_rival_test[1], "forecast_skill_min_inside_95": 0.443,
         "forecast_skill_max_inside_95": 0.925, "n_starts_inside_95": 3, "alk0_fitted_meq_L": 80.0,
         "alk0_converged_minima_inside_95_meq_L": "22.6; 80.0"},
        {"model": "alkalinity_FA_SCENARIO_T25C", "cal_converged": False, "selected": False, "sel_rmse_log": 0.15,
         "sel_skill_min_inside_95": np.nan, "sel_skill_max_inside_95": np.nan, "alk0_fitted_meq_L": 23.6,
         "alk0_converged_minima_inside_95_meq_L": np.nan},
    ])


def test_selection_holds_but_the_out_of_sample_ranking_is_undetermined():
    fit = _fit_rows()
    sel, rival = fit.iloc[1], R.best_rival(fit)
    assert rival["model"] == "alkalinity_limited"
    txt = R.selection_robustness_txt(fit, sel, rival)
    assert "biomass_decay holds against every converged start" in txt and "0.958 against the highest" in txt
    assert "out-of-sample ranking of biomass_decay against alkalinity_limited is undetermined" in txt
    assert "test 0.748 against 0.227 to 0.817" in txt
    assert "selection -0.908 to 0.935, test 0.227 to 0.817 and forecast 0.443 to 0.925" in R.rival_range_txt(rival)
    assert "Alk0 22.6 and 80.0 meq/L" in R.rival_range_txt(rival)
    worse = _fit_rows().assign(sel_skill_max_inside_95=[0.0, 0.9578, 0.97, np.nan])
    assert "does not hold" in R.selection_robustness_txt(worse, worse.iloc[1], R.best_rival(worse))


def test_every_converged_candidate_of_the_other_mechanism_is_ranked():
    fit = _fit_rows()
    fa = fit.index[fit["model"] == "alkalinity_FA_SCENARIO_T25C"][0]
    fit.loc[fa, ["cal_converged", "sel_rmse_log", "sel_skill_min_inside_95", "sel_skill_max_inside_95",
                 "test_skill_min_inside_95", "test_skill_max_inside_95", "forecast_skill_min_inside_95",
                 "forecast_skill_max_inside_95", "n_starts_inside_95"]] = [True, 0.15, 0.79, 0.793, 0.641, 0.644, 0.801,
                                                                          0.803, 3]
    others = R.other_mechanism(fit)
    assert [r["model"] for r in others] == ["alkalinity_FA_SCENARIO_T25C", "alkalinity_limited"]  # best first
    txt = R.selection_robustness_txt(fit, fit.iloc[1], others)
    assert "biomass_decay is better than alkalinity_FA_SCENARIO_T25C over all those starts" in txt
    assert "ranking of biomass_decay against alkalinity_limited is undetermined" in txt
    assert txt.endswith("do not discriminate activity decay from substrate-side (alkalinity) limitation.")


def test_disjoint_ranges_give_a_determined_ranking():
    fit = _fit_rows(alk_rival_test=(0.1, 0.3)).assign(forecast_skill_max_inside_95=[np.nan, 0.863, 0.5, np.nan])
    txt = R.selection_robustness_txt(fit, fit.iloc[1], R.best_rival(fit))
    assert "biomass_decay is better than alkalinity_limited over all those starts" in txt


# --- Y2: Alk0 evidence from converged candidates only; one entry per distinct minimum --------------------------------
def test_unconverged_alk0_never_enters_the_evidence_range():
    txt = R.alk0_evidence_txt(_fit_rows())
    assert txt.startswith("Alk0 of the converged alkalinity candidates 22.6-80.0 meq/L")
    assert "23.6" not in txt and "not converged (no estimate): alkalinity_FA_SCENARIO_T25C" in txt
    none = R.alk0_evidence_txt(_fit_rows(alk_converged=False))
    assert none.startswith("no alkalinity candidate converged") and "80.0" not in none and "23.6" not in none


def test_two_starts_within_the_tolerance_are_one_minimum(cfg):
    assert M.distinct_minima([(22.637, 3.02), (22.614, 3.01), (80.0, 0.0)]) == [(80.0, 0.0, 1), (22.614, 3.01, 2)]
    a = 20.0
    edge = a * np.exp(M.MULTIMODAL_LOG_TOL) * 0.999  # just inside the 5 % tolerance
    assert len(M.distinct_minima([(a, 0.0), (edge, 1.0)])) == 1
    assert len(M.distinct_minima([(a, 0.0), (a * np.exp(M.MULTIMODAL_LOG_TOL) * 1.01, 1.0)])) == 2
    # in the profile status: two starts near the same second minimum are listed once
    spec = {s.name: s for s in M.model_specs(cfg)}["alkalinity_limited"]
    data = _synthetic_data(spec, {"v0_mgN_L_d": 12.0, "Alk0_meq_L": 15.0})
    lay = M.layout(spec, data)
    x0 = M.initial_guess(spec, data, lay)
    cal = M.fit(spec, data, lay, (1.0, 16.0), x0, np.ones(len(x0), dtype=bool), max_nfev=M.CAL_MAX_NFEV)
    k = lay.names.index("Alk0_meq_L")
    far1, far2 = cal.x.copy(), cal.x.copy()
    far1[k] += np.log(5.0)
    far2[k] += np.log(5.0 * 1.01)
    twins = [replace(cal, x=far1, rss=cal.rss * 1.001), replace(cal, x=far2, rss=cal.rss * 1.0012)]
    r = M.alk0_profile(spec, data, lay, cal, [cal, *twins])
    assert r["alk0_profile_status"].count("meq/L at Δ") == 1 and "(2 starts)" in r["alk0_profile_status"]
    assert len(r["alk0_converged_minima_inside_95_meq_L"].split("; ")) == 2


# --- Y3: the limits of the idle bounds come from the synthetic check --------------------------------------------------
def test_interval_coverage_check_and_the_limits_text():
    cov = ST.interval_coverage(phis=(0.0, 0.9), n_windows=150, seed=7)
    assert list(cov["phi"]) == [0.0, 0.9] and cov["reference"].tolist() == [False, True]
    r = cov.set_index("phi")
    assert r.loc[0.0, "coverage_scaled_true_rho"] == pytest.approx(0.95, abs=0.05)
    assert r.loc[0.9, "coverage_scaled_true_rho"] == pytest.approx(0.95, abs=0.05)
    assert r.loc[0.9, "coverage_unscaled"] < r.loc[0.9, "coverage_scaled_estimated_rho"] < 0.95  # residual shortfall
    assert r.loc[0.9, "median_rho_hat"] < 0.9  # ρ̂ is biased low on 101 points
    fake = pd.DataFrame({"phi": [0.9, 0.98], "n_points": 101, "n_windows": 600, "nominal": 0.95,
                         "coverage_scaled_estimated_rho": [0.851, 0.733], "median_rho_hat": [0.85, 0.93],
                         "reference": [True, False]})
    fr = pd.DataFrame({"interval": ["whole record"], "median_nh4_slope_resid_lag1": [0.926],
                       "fraction_active_nitrification": [0.056], "fraction_oxygen_limited": [0.013]})
    txt = ST.interval_limits_text(fake, fr)
    assert "0.85 of the windows at φ = 0.9 instead of 0.95" in txt and "0.73 at φ = 0.98" in txt
    assert "record's median ρ̂ of C(t) (0.93)" in txt and "(6.9 % of the whole record)" in txt
    assert "lower bound counts only bins classified inactive" in txt and "active and oxygen-limited bins" in txt


# --- Y4: the allow-list of NH4-independent metrics and the declared columns of primary wide tables --------------------
def _rows(*specs):
    return pd.DataFrame([{"metric": m, "basis": b, "nh4_observation_model": o, "nh4_series": s} for m, b, o, s in specs])


PRIMARY = ("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL, "C_inv:aaa")
ALLOWED = ("t_pulse_d", "segment_end_d", "median_pH", "o2_supply_capacity_mgO2_SCENARIO_kLa1_T25C",
           "o2_supply_capacity_mgO2_SCENARIO_kLa10_T25C", "o2_supply_capacity_mgO2_SCENARIO_kLa100_T25C",
           "o2_supply_capacity_mgO2_SCENARIO_kLa0.5_T25C", "free_energy_per_fully_converted_pulse_J_SCENARIO_T25C",
           "stirrer_energy_J")


def test_every_allow_listed_nh4_independent_metric_passes():
    rows = _rows(PRIMARY, *[(m, "fact", M.NH4_INDEPENDENT, "") for m in ALLOWED])
    assert M.check_single_observation_model(rows) == M.OBSERVATION_MODEL


@pytest.mark.parametrize("metric", ["consumption_weighted_mgN_L_d", "eta_conversion_fraction", "n_removed_mgN",
                                    "idle_fraction_LOWER_BOUND", "o2_supply_capacity_mgO2_SCENARIO_kLa_derived_T25C",
                                    "median_pH_of_C", "aeration_energy_supply_kWh_per_kgN_SCENARIO_AE1_kLa10_T25C"])
def test_a_derived_metric_declared_nh4_independent_raises(metric):
    with pytest.raises(M.ObservationModelError, match="outside the allow-list"):
        M.check_single_observation_model(_rows(PRIMARY, (metric, "", M.NH4_INDEPENDENT, "")))


def _primary_wide(extra=None):
    df = pd.DataFrame({"sub_segment": [1, 2], "fit_start_d": [1.0, 2.0], "consumption_mgN_L_d": [3.0, 4.0],
                       "status": ["computed"] * 2, "nh4_observation_model": M.OBSERVATION_MODEL,
                       "nh4_truth_id": "C_inv:aaa"})
    return df.assign(**{extra: [1.0, 2.0]}) if extra else df


def test_an_undeclared_column_in_a_primary_wide_table_raises():
    kw = {"label_cols": ("status",), "series_col": "nh4_truth_id", "independent_cols": ("sub_segment", "fit_start_d"),
          "derived_cols": ("consumption_mgN_L_d",)}
    assert M.check_single_observation_model(_primary_wide(), **kw) == M.OBSERVATION_MODEL
    with pytest.raises(M.ObservationModelError, match="undeclared columns in a wide table with primary rows"):
        M.check_single_observation_model(_primary_wide("consumption_SENSITIVITY_gain_mgN_L_d_x"), **kw)
    with pytest.raises(M.ObservationModelError, match="undeclared columns"):
        M.check_single_observation_model(_primary_wide("reading_slope_mgN_L_d"), **kw)
    with pytest.raises(M.ObservationModelError, match="must declare its NH4-independent and its C\\(t\\)-derived"):
        M.check_single_observation_model(_primary_wide(), label_cols=("status",), series_col="nh4_truth_id",
                                         independent_cols=("sub_segment", "fit_start_d"))
    indep = pd.DataFrame({"phi": [0.9], "coverage": [0.85], "nh4_observation_model": M.NH4_INDEPENDENT,
                          "nh4_series": ""})
    assert M.check_single_observation_model(indep, label_cols=(), independent_cols=("phi", "coverage")) is None
    with pytest.raises(M.ObservationModelError, match="not declared NH4-independent"):
        M.check_single_observation_model(indep.assign(C_inv_mgN_L=1.0), label_cols=(),
                                         independent_cols=("phi", "coverage"))


# --- real record (run-all outputs) ---------------------------------------------------------------------------------------
def test_real_report_states_the_skill_range_and_the_limits(run_all_dir):
    fit = pd.read_csv(run_all_dir / "tables" / "model_fit.csv").set_index("model")
    st = pd.read_csv(run_all_dir / "tables" / "model_calibration_starts.csv")
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    point7 = next(line for line in rep.splitlines() if line.startswith("7. **Mechanism.**"))
    alk = fit.loc["alkalinity_limited"]
    inside = st[(st["model"] == "alkalinity_limited") & st["inside_95"]]
    assert int(alk["n_starts_inside_95"]) == len(inside)
    sec4 = rep.split("## 4.")[1].split("## 5.")[0]
    if len(inside) > 1:  # the range of the indistinguishable starts is in the §4 table (pass 5, Z4; was §1 point 7)
        for c in ("sel", "test", "forecast"):
            v = inside[f"{c}_skill_vs_null_diagnostic"]
            assert alk[f"{c}_skill_min_inside_95"] == pytest.approx(v.min())
            assert f"{v.min():.3f} to {v.max():.3f}" in sec4
    assert "holds against every converged start" in point7 or "does not hold against" in point7
    for name, g in st.groupby("model"):  # the chosen start reproduces the candidate's own skills
        ch = g[g["chosen"]].iloc[0]
        for c in ("sel", "test", "forecast"):
            assert ch[f"{c}_skill_vs_null_diagnostic"] == pytest.approx(fit.loc[name, f"{c}_skill_vs_null"], abs=1e-12)
    assert not st.loc[~st["converged"], "inside_95"].any()
    point8 = next(line for line in rep.splitlines() if line.startswith("8. **Energy and efficiency.**"))
    states_md = (run_all_dir / "reports" / "states.md").read_text()
    cov = pd.read_csv(run_all_dir / "tables" / "states_interval_coverage.csv")
    ref = cov[cov["reference"]].iloc[0]
    phrase = f"{ref['coverage_scaled_estimated_rho']:.2f} of the windows at φ = {ref['phi']:g}"
    for text in (point8, rep.split("## 3.")[1].split("## 4.")[0], states_md):
        assert phrase in text and "conditional on the state classification" in text
    row = next(line for line in rep.splitlines() if line.startswith("| Model selection, activity decay"))
    for m in fit.index[~fit["cal_converged"].astype(bool)]:
        assert f"{fit.loc[m, 'alk0_fitted_meq_L']:.1f}" not in row.split("Alk0 of the converged")[1].split(";")[0]


def test_real_states_tables_pass_the_extended_checker(run_all_dir):
    tab = run_all_dir / "tables"
    for name, key in (("states_fractions.csv", "fractions"), ("states_timeline.csv", "timeline"),
                      ("states_comparison.csv", "comparisons")):
        df = pd.read_csv(tab / name)
        indep, derived = ST.STATE_TABLE_COLUMNS[key]
        assert M.check_single_observation_model(df, label_cols=(), independent_cols=indep, derived_cols=derived)
        with pytest.raises(M.ObservationModelError):  # no declaration: raises
            M.check_single_observation_model(df, label_cols=())
    bins = pd.read_parquet(run_all_dir / "processed" / "states_bins.parquet")
    indep, derived = ST.STATE_TABLE_COLUMNS["bins"]
    assert M.check_single_observation_model(bins, label_cols=(), independent_cols=indep, derived_cols=derived)
    cov = pd.read_csv(tab / "states_interval_coverage.csv")
    assert M.check_single_observation_model(cov, label_cols=(), independent_cols=ST.COVERAGE_COLUMNS) is None
    pm = pd.read_csv(tab / "pulses.csv")
    ind = pm[pm["nh4_observation_model"] == M.NH4_INDEPENDENT]
    assert len(ind) and ind["metric"].str.match(M._NH4_INDEPENDENT_RE).all()


def _pass3(name: str, fixtures: Path = PASS3_FIXTURES) -> pd.DataFrame:
    """A stored pass-3 reference table; a missing fixture fails the test (never skips)."""
    path = fixtures / name
    if not path.is_file():
        pytest.fail(f"pass-3 fixture {path} missing (extracted from {PASS3_COMMIT}; docs/progress.md, pass 5)")
    return pd.read_csv(path)


def test_a_missing_pass3_fixture_fails_the_regression_test(tmp_path):
    with pytest.raises(pytest.fail.Exception, match="fixture .* missing"):
        _pass3("model_fit_selected.csv", fixtures=tmp_path)
    for name in ("model_fit_selected.csv", "model_identifiability.csv", "nh4_truth_summary.csv",
                 "nh4_truth_vectors.csv", "pulses_model_based.csv"):
        assert not _pass3(name).empty


def test_real_selection_fit_and_model_based_numbers_equal_pass3(run_all_dir):
    """Regression (Step 1): no fitted value of the selected model, no selection and no MODEL-BASED number changes."""
    tab = run_all_dir / "tables"
    new, old = pd.read_csv(tab / "model_fit.csv"), _pass3("model_fit_selected.csv")
    assert new.loc[new["selected"], "model"].tolist() == old.loc[old["selected"], "model"].tolist()
    sel = new.loc[new["selected"], "model"].iloc[0]
    cols = ["cal_rmse_log", "cal_start", "sel_skill_vs_null", "test_skill_vs_null", "forecast_skill_vs_null",
            "whole_record_refit_rmse_log"]
    pd.testing.assert_frame_equal(new.set_index("model").loc[[sel], cols], old.set_index("model").loc[[sel], cols])
    for name in ("model_identifiability.csv", "nh4_truth_summary.csv", "nh4_truth_vectors.csv"):
        # pass 6 (G8) appends provenance columns only: every column of the fixture is still compared exactly
        got, ref = pd.read_csv(tab / name), _pass3(name)
        extra = [c for c in got.columns if c not in ref.columns]
        assert extra == ([] if name == "nh4_truth_summary.csv" else ["nh4_observation_model", "nh4_series", "assumptions"])
        pd.testing.assert_frame_equal(got[list(ref.columns)], ref)
    pm, pm_old = pd.read_csv(tab / "pulses.csv"), _pass3("pulses_model_based.csv")
    key = ["pulse", "metric"]
    mb = pm[pm["metric"].str.contains("MODEL_BASED")].set_index(key).sort_index()
    mb_old = pm_old[pm_old["metric"].str.contains("MODEL_BASED")].set_index(key).sort_index()
    pd.testing.assert_series_equal(mb["value"], mb_old["value"])


# --- Y6: the Haldane laws are fitted as (α, Ks, S_opt), the same function for every candidate that has one ---------------
def test_haldane_reparametrization_is_the_same_rate_law_for_every_haldane_candidate(cfg):
    specs = M.model_specs(cfg)
    hald = [s for s in specs if s.fa]
    assert {s.name.split("_SCENARIO")[0] for s in hald} == {"haldane_FA", "combined", "alkalinity_FA"}
    for s in hald:  # one parametrization for every candidate that shares the Haldane pair
        assert s.kinetic[0] in M.HALDANE_ALPHA and set(M.HALDANE_SHAPE) <= set(s.kinetic)
        assert "vmax_mgN_L_d" not in s.kinetic and "Ki_mgNH3_L" not in s.kinetic
    old = {"vmax_mgN_L_d": 282.0, "v0_mgN_L_d": 282.0, "kd_per_d": 0.03, "Ks_mgNH3_L": 5.3, "Ki_mgNH3_L": 0.045,
           "Alk0_meq_L": 23.6}
    for s in hald:
        new = {k: old[k] for k in s.kinetic if k in old}
        new[s.kinetic[0]] = old["vmax_mgN_L_d"] / old["Ks_mgNH3_L"]
        new["Sopt_mgNH3_L"] = np.sqrt(old["Ks_mgNH3_L"] * old["Ki_mgNH3_L"])
        nat = M.haldane_natural(new)
        assert nat["Ki_mgNH3_L"] == pytest.approx(0.045) and nat.get("vmax_mgN_L_d", nat.get("v0_mgN_L_d")) == pytest.approx(282)
        ref = {k: old[k] for k in ("vmax_mgN_L_d" if not s.decay else "v0_mgN_L_d", "kd_per_d", "Ks_mgNH3_L",
                                   "Ki_mgNH3_L", "Alk0_meq_L") if (k != "kd_per_d" or s.decay) and (k != "Alk0_meq_L" or s.alk)}
        a = M.simulate(s, ref, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0).c
        b = M.simulate(s, new, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0).c
        np.testing.assert_allclose(b, a, rtol=1e-12)
    data = _synthetic_data(specs[0], {})
    for s in hald:  # the policy starts: the transformed default start x 0.5, 1, 2, 4, observation start shared
        lay = M.layout(s, data)
        x0 = M.initial_guess(s, data, lay)
        k = lay.names.index("Sopt_mgNH3_L")
        assert np.exp(x0[k]) == pytest.approx(np.sqrt(0.5)) and np.exp(x0[0]) == pytest.approx(10.0)
        assert [round(float(np.exp(xs[k])), 6) for _, xs in M.calibration_starts(s, lay, x0)] == [
            round(f * np.sqrt(0.5), 6) for f in M.CAL_START_FACTORS]


# --- Y7: the model per SCENARIO volume --------------------------------------------------------------------------------
def test_volume_simulation_at_zero_loss_is_the_constant_volume_model(cfg):
    from bnr_neq import volume

    for s in M.model_specs(cfg):
        kin = {"v0_mgN_L_d": 8.0, "kd_per_d": 0.05, "alpha_mgN_L_d_per_mgNH3_L": 10.0, "alpha0_mgN_L_d_per_mgNH3_L": 10.0,
               "Ks_mgNH3_L": 0.5, "Sopt_mgNH3_L": 0.7, "Alk0_meq_L": 15.0}
        kin = {k: kin[k] for k in s.kinetic}
        a = M.simulate(s, kin, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0)
        b = M.simulate_volume(s, kin, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0,
                              volume.VolumeScenario("SCENARIO 0", "SCENARIO", 0.1, 0.0))
        np.testing.assert_allclose(b.c, a.c, rtol=1e-12)
        assert b.consumed == pytest.approx(a.consumed, rel=1e-12)


def test_volume_simulation_conserves_the_nitrogen_amount(cfg):
    from bnr_neq import volume

    spec = {s.name: s for s in M.model_specs(cfg)}["biomass_decay"]
    v = volume.VolumeScenario("SCENARIO evaporation 1 mL/d", "SCENARIO", 0.1, 1e-3)
    sim = M.simulate_volume(spec, {"v0_mgN_L_d": 8.0, "kd_per_d": 0.05}, 80.0, T_SYN, PH_SYN, [4.0, 8.0, 12.0], 60.0, v)
    amount_end = sim.c[-1] * v.volume_L(T_SYN[-1])  # mg N in the vessel at the end
    amount_in = 80.0 * v.volume_L(T_SYN[0]) + 3 * 60.0 * v.v0_L  # initial + three pulse masses (dose at V0)
    assert amount_end == pytest.approx(amount_in - sim.consumed * v.v0_L, rel=1e-6)  # consumed is per initial litre
    for t, before, after in sim.pulse_jumps:
        assert after - before == pytest.approx(60.0 * float(v.factor(t)))


def test_scenario_refit_recovers_the_generating_parameters_better_than_constant_volume(cfg):
    from bnr_neq import volume

    spec = {s.name: s for s in M.model_specs(cfg)}["biomass_decay"]
    v = volume.VolumeScenario("SCENARIO evaporation 1 mL/d", "SCENARIO", 0.1, 1e-3)
    true = {"v0_mgN_L_d": 8.0, "kd_per_d": 0.05}
    t = np.round(np.arange(1.01, 30.0, 0.02), 10)
    ph = 7.6 - 0.01 * (t - 1)
    pulses = [4.0, 10.0, 18.0, 25.0]
    sim = M.simulate_volume(spec, true, 80.0, t, ph, pulses, 60.0, v)
    lr = M.observe(sim.c, t, OBS_SYN, []) + 0.01 * np.random.default_rng(5).normal(size=len(t))
    errors, rss = {}, {}
    for vol in (None, v):
        data = M.ModelData(t=t, log_r=lr, ph=ph, pulses=pulses, dose=60.0, offsets=[], offset_guess=[], volume=vol)
        lay = M.layout(spec, data)
        x0 = M.initial_guess(spec, data, lay)
        f = M.fit(spec, data, lay, (1.0, 30.0), x0, np.ones(len(x0), dtype=bool), max_nfev=M.CAL_MAX_NFEV)
        p = M.unpack(spec, data, lay, f.x)["kin"]
        errors[vol is None] = max(abs(np.log(p[k] / true[k])) for k in true)
        rss[vol is None] = f.rss
    assert errors[False] < 0.02 and errors[False] < errors[True]  # the scenario at the generating rate recovers better
    assert rss[False] < rss[True]


def test_a_scenario_refit_row_without_scenario_raises():
    ok = pd.DataFrame({"volume_scenario": ["SCENARIO evaporation 0.5 mL/d"], "role": [M.SCENARIO_ROLE]})
    M.check_scenario_rows(ok)
    for bad in (ok.assign(volume_scenario="evaporation 0.5 mL/d"), ok.assign(role="primary"),
                ok.drop(columns=["volume_scenario"])):
        with pytest.raises(M.ScenarioLabelError):
            M.check_scenario_rows(bad)


def test_real_volume_scenario_tables_are_scenario_labelled(run_all_dir, cfg):
    from bnr_neq import volume

    tab = run_all_dir / "tables"
    cands, summ = pd.read_csv(tab / "model_volume_scenarios.csv"), pd.read_csv(tab / "model_volume_scenario_summary.csv")
    rates = sorted(1e3 * v.net_loss_L_per_d for v in volume.volume_scenarios(cfg, M.END_D) if not v.is_constant)
    assert sorted(summ["net_loss_mL_per_d"]) == pytest.approx(rates) and len(rates) == 3
    for df in (cands, summ):
        M.check_scenario_rows(df)
        assert df["volume_scenario"].str.contains("SCENARIO").all()
    fit = pd.read_csv(tab / "model_fit.csv")
    assert set(cands["model"]) == set(fit["model"]) and len(cands) == len(rates) * len(fit)
    assert (cands.groupby("volume_scenario")["selected_in_scenario"].sum() == 1).all()
    sel = fit.loc[fit["selected"], "model"].iloc[0]
    assert (summ["primary_selected"] == sel).all()
    for _, r in summ.iterrows():
        chosen = cands[(cands["volume_scenario"] == r["volume_scenario"]) & cands["selected_in_scenario"]]["model"].iloc[0]
        assert r["selected_in_scenario"] == chosen and bool(r["selected_model_changes"]) == (chosen != sel)
        assert r["C_start_vs_constant_volume_band"] in ("inside the band", "above the band", "below the band")
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    assert "M6 refitted per SCENARIO evaporation rate" in rep
    # no primary table names a volume-scenario refit
    for name in ("pulses.csv", "nh4_truth_windows.csv", "nh4_truth_summary.csv"):
        assert "SCENARIO refit" not in (tab / name).read_text()
