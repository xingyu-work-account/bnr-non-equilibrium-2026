"""Remediation pass 5 (report-only): the volume-scenario C(t) against φ = V0/V(t) (Z1), bound-constrained calibrations
(Z2), the shortened §1 point 7 and its §4 table (Z4), and the regression of every pre-existing column against `8e10d0e`."""

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bnr_neq import model as M
from bnr_neq import report as R
from bnr_neq import volume

PASS4 = Path(__file__).resolve().parent / "fixtures" / "pass4"  # tables of 8e10d0e (git show; docs/progress.md, pass 5)
Z1_COLS = ["C_start_ratio_to_constant_volume", "phi_start_median", "C_end_ratio_to_constant_volume", "phi_end_median",
           "ratio_within_phi_tolerance"]
NEW_COLS = {"model_fit.csv": ["cal_at_bounds"], "model_calibration_starts.csv": ["cal_at_bounds"],
            "model_volume_scenario_summary.csv": Z1_COLS}

T = np.round(np.arange(1.0, 20.0, 0.01), 10)
PH = 7.6 - 0.02 * (T - 1)
PULSES, DOSE = [4.0, 8.0, 12.0], 60.0


def _summary(rows):
    """A volume-scenario summary with the columns the Z1 text reads."""
    return pd.DataFrame([{"net_loss_mL_per_d": e, "C_start_vs_constant_volume_band": "inside the band",
                          "C_end_vs_constant_volume_band": band, "C_start_ratio_to_constant_volume": rs,
                          "phi_start_median": ps, "C_end_ratio_to_constant_volume": re, "phi_end_median": pe}
                         for e, band, rs, ps, re, pe in rows])


# --- Z1: the scenario C(t) is the constant-volume C(t) x φ within the registered tolerance ----------------------------
@pytest.mark.parametrize("name,kin", [("null_no_consumption", {}),
                                      ("biomass_decay", {"v0_mgN_L_d": 5.0, "kd_per_d": 0.05})])
def test_a_record_generated_at_a_known_evaporation_rate_has_ratio_phi(cfg, name, kin):
    """Positive: C(t) generated under E = 1 mL/d against the same law at constant volume. With conserved amounts the
    ratio is φ (exactly for no consumption; to the Monod factor C/(K_NH + C) with consumption)."""
    spec = {s.name: s for s in M.model_specs(cfg)}[name]
    vol = volume.VolumeScenario("SCENARIO evaporation 1 mL/d", "SCENARIO", 0.1, 1e-3)
    c_cv = M.simulate(spec, kin, 80.0, T, PH, PULSES, DOSE).c
    c = M.simulate_volume(spec, kin, 80.0, T, PH, PULSES, DOSE, vol).c
    start, end = T < 1.1, T >= T[-1] - 0.1
    z = M.phi_scaling(c, c_cv, T, vol, start, end)
    assert z["phi_end_median"] == pytest.approx(float(np.median(0.1 / (0.1 - 1e-3 * T[end]))), rel=1e-12)
    assert z["phi_end_median"] > 1.2  # the end window tests the scaling (φ − 1 far above the tolerance)
    for w in ("start", "end"):
        assert M.within_phi_tolerance(z[f"C_{w}_ratio_to_constant_volume"], z[f"phi_{w}_median"])
    assert z["ratio_within_phi_tolerance"]
    if name == "null_no_consumption":  # the amount is conserved exactly; C0 is the concentration at T[0], in V(T[0])
        amount = c[end] * vol.volume_L(T[end])
        assert amount == pytest.approx(80.0 * vol.volume_L(T[0]) + len(PULSES) * DOSE * vol.v0_L, rel=1e-9)
    txt = M.volume_scaling_text(_summary([(1.0, "above the band", z["C_start_ratio_to_constant_volume"],
                                           z["phi_start_median"], z["C_end_ratio_to_constant_volume"],
                                           z["phi_end_median"])]))
    assert txt.startswith("Every scenario C(t) is the constant-volume C(t) scaled by φ = V0/V(t)")
    assert "\"above the band\" at 1 mL/d is the evaporative concentration" in txt


def test_a_ratio_outside_the_tolerance_says_it_differs(cfg):
    """Negative: an unscaled C(t) (the constant-volume one taken as the scenario's) is not within the tolerance, and
    the text says so instead of the scaling statement."""
    vol = volume.VolumeScenario("SCENARIO evaporation 1 mL/d", "SCENARIO", 0.1, 1e-3)
    spec = {s.name: s for s in M.model_specs(cfg)}["null_no_consumption"]
    c_cv = M.simulate(spec, {}, 80.0, T, PH, PULSES, DOSE).c
    z = M.phi_scaling(c_cv, c_cv, T, vol, T < 1.1, T >= T[-1] - 0.1)
    assert not z["ratio_within_phi_tolerance"]
    txt = M.volume_scaling_text(_summary([(0.5, "inside the band", 1.01, 1.005, 1.10, 1.31),
                                          (1.0, "above the band", 1.02, 1.01, 1.94, 1.92)]))
    assert txt.startswith("The scenario C(t) is not the constant-volume C(t) scaled by φ")
    assert "0.5 mL/d, end: ratio 1.100 against φ 1.310" in txt and "1 mL/d" not in txt.split("tolerance")[1]
    assert "Every scenario C(t)" not in txt and "η" not in txt


def test_a_ratio_exactly_at_the_tolerance_takes_the_scaled_branch():
    """Boundary: |ratio − φ|/φ = 0.05 exactly (21/20 is exact in binary arithmetic) is inside (inclusive)."""
    assert M.phi_deviation(21.0, 20.0) == M.PHI_SCALING_TOL == 0.05
    assert M.within_phi_tolerance(21.0, 20.0) and M.within_phi_tolerance(19.0, 20.0)
    assert not M.within_phi_tolerance(math.nextafter(21.0, 22.0), 20.0)
    txt = M.volume_scaling_text(_summary([(1.0, "above the band", 19.0, 20.0, 21.0, 20.0)]))
    assert txt.startswith("Every scenario C(t) is the constant-volume C(t) scaled")
    assert "start window cannot tell" not in txt  # φ − 1 at the start is far above the tolerance here
    assert M.volume_scaling_text(pd.DataFrame()) == ""


# --- Z2: bound-constrained calibrations are flagged ------------------------------------------------------------------------
def _fit_kd(cfg, kd_true):
    """A synthetic biomass-decay record with the observation parameters fixed at their true values; v0, kd and C0 are
    free. A true decay constant of zero pushes the fitted kd onto its lower bound (1e-4 /d); v0 and C0 stay inside."""
    spec = {s.name: s for s in M.model_specs(cfg)}["biomass_decay"]
    kin = {"v0_mgN_L_d": 12.0, "kd_per_d": kd_true}
    obs = {"a": 0.1, "b0": 0.8, "b1_per_d": 0.0, "tau_d": 0.05, "offset_values": []}
    sim = M.simulate(spec, kin, 80.0, T, PH, PULSES, DOSE)
    lr = M.observe(sim.c, T, obs, []) + 0.01 * np.random.default_rng(5).normal(size=len(T))
    data = M.ModelData(t=T, log_r=lr, ph=PH, pulses=PULSES, dose=DOSE, offsets=[], offset_guess=[])
    lay = M.layout(spec, data)
    x0 = M.initial_guess(spec, data, lay)
    x0[lay.names.index("a")], x0[lay.names.index("b0")] = obs["a"], obs["b0"]
    x0[lay.names.index("b1_per_d")], x0[lay.names.index("tau_d")] = 0.0, math.log(obs["tau_d"])
    free = np.array([n in ("v0_mgN_L_d", "kd_per_d", "C0_mgN_L") for n in lay.names])
    f = M.fit(spec, data, lay, (T[0], T[-1] + 0.01), x0, free, max_nfev=2000)
    return lay, f


def test_a_fit_pushed_onto_a_bound_is_flagged_with_that_parameter(cfg):
    lay, f = _fit_kd(cfg, 0.0)
    i = lay.names.index("kd_per_d")
    assert f.converged and f.x[i] - lay.lo[i] <= M.BOUND_TOL_FRACTION * (lay.hi[i] - lay.lo[i])
    assert M.at_bounds(lay, f.x, f.free) == "kd_per_d"


def test_a_fit_inside_its_bounds_has_an_empty_flag(cfg):
    lay, f = _fit_kd(cfg, 0.05)
    assert f.converged and M.at_bounds(lay, f.x, f.free) == ""


def test_a_parameter_exactly_at_the_tolerance_distance_is_flagged():
    """Boundary: bound width 1000, so the tolerance is exactly 1.0; inclusive on both sides; a fixed parameter on
    its bound is not flagged (only free parameters are)."""
    lay = M.Layout(["p", "q"], np.array([0.0, 0.0]), np.array([1000.0, 1000.0]), [False, False])
    assert M.BOUND_TOL_FRACTION * 1000.0 == 1.0
    free = np.array([True, False])
    assert M.at_bounds(lay, np.array([1.0, 0.0]), free) == "p"
    assert M.at_bounds(lay, np.array([999.0, 0.0]), free) == "p"
    assert M.at_bounds(lay, np.array([math.nextafter(1.0, 2.0), 0.0]), free) == ""
    assert M.at_bounds(lay, np.array([math.nextafter(999.0, 0.0), 1000.0]), free) == ""
    assert M.at_bounds(lay, np.array([0.0, 0.0]), np.array([True, True])) == "p; q"


def test_bounds_text_names_every_flagged_calibration():
    fit = pd.DataFrame({"model": ["null_no_consumption", "biomass_decay", "haldane_FA_SCENARIO_T20C"],
                        "selected": [False, True, False], "cal_at_bounds": [np.nan, "C0_mgN_L", "Ks_mgNH3_L; a"]})
    txt = M.bounds_text(fit)
    assert "biomass_decay (C0_mgN_L); haldane_FA_SCENARIO_T20C (Ks_mgNH3_L, a)" in txt
    assert "natural-unit kinetic parameters are not estimates" in txt and "null_no_consumption" not in txt
    assert "selected model (biomass_decay) is that of its whole-record refit" in txt
    assert M.bounds_text(fit.assign(cal_at_bounds="")).startswith("No calibration is bound-constrained")
    # the reference has no kinetic parameter: named in its own clause, never among the candidates
    ref = M.bounds_text(fit.assign(cal_at_bounds=["a", "", ""]))
    assert ref.startswith("No candidate calibration is bound-constrained") and "not estimates" not in ref
    assert "The reference null_no_consumption (no kinetic parameter) has a at a bound." in ref
    both = M.bounds_text(fit.assign(cal_at_bounds=["a", "C0_mgN_L", ""]))
    assert "): biomass_decay (C0_mgN_L)." in both and both.count("null_no_consumption") == 1


# --- Z4: §1 point 7 in short, with the ranking rule of pass 4 ----------------------------------------------------------
def _fit_table():
    base = {"cal_converged": True, "n_starts_inside_95": 3, "alk0_converged_minima_inside_95_meq_L": np.nan,
            "cal_at_bounds": ""}
    rows = [("null_no_consumption", 0.0, (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), False),
            ("biomass_decay", 0.11, (0.958, 0.958), (0.747, 0.748), (0.862, 0.863), True),
            ("alkalinity_limited", 0.13, (-0.908, 0.935), (0.227, 0.812), (0.443, 0.921), False),
            ("alkalinity_FA_SCENARIO_T25C", 0.12, (0.790, 0.793), (0.641, 0.644), (0.801, 0.803), False)]
    return pd.DataFrame([{**base, "model": m, "sel_rmse_log": r, "selected": s,
                          "sel_skill_min_inside_95": a[0], "sel_skill_max_inside_95": a[1],
                          "test_skill_min_inside_95": b[0], "test_skill_max_inside_95": b[1],
                          "forecast_skill_min_inside_95": c[0], "forecast_skill_max_inside_95": c[1]}
                         for m, r, a, b, c, s in rows])


def test_point7_summary_states_the_ranking_with_one_range_each():
    fit = _fit_table()
    sel = fit.iloc[1]
    txt = R.mechanism_summary_txt(fit, sel, R.other_mechanism(fit))
    assert txt.startswith("the selection holds against every converged start")
    assert "better than alkalinity_FA_SCENARIO_T25C (test 0.641 to 0.644)" in txt
    assert "ranking against alkalinity_limited (test 0.227 to 0.812) is undetermined" in txt
    assert R.overall_rank(sel, fit.iloc[3]) == "better" and R.overall_rank(sel, fit.iloc[2]) == "undetermined"
    tab = R.mechanism_table(fit).set_index("candidate")
    assert tab.loc["alkalinity_limited", "test skill"] == "0.227 to 0.812"
    assert tab.loc["alkalinity_limited", "biomass_decay against it (test; forecast)"] == "undetermined; undetermined"
    assert tab.loc["alkalinity_FA_SCENARIO_T25C", "biomass_decay against it (test; forecast)"] == "better; better"
    assert tab.loc["biomass_decay", "biomass_decay against it (test; forecast)"] == "selected"


# --- real record (run-all outputs) ---------------------------------------------------------------------------------------
def test_real_pre_existing_columns_equal_8e10d0e(run_all_dir):
    """Regression: every pre-existing column of the three tables that gain columns in this pass equals `8e10d0e`; the
    only new columns are those of Z1 and Z2."""
    for name, new_cols in NEW_COLS.items():
        ref = PASS4 / name
        if not ref.is_file():
            pytest.fail(f"fixture {ref} missing (extracted with git show 8e10d0e; docs/progress.md, pass 5)")
        old, new = pd.read_csv(ref), pd.read_csv(run_all_dir / "tables" / name)
        assert [c for c in new.columns if c not in old.columns] == new_cols
        pd.testing.assert_frame_equal(new[list(old.columns)], old)


def test_real_volume_scaling_columns_and_statement(run_all_dir, cfg):
    tab = run_all_dir / "tables"
    summ = pd.read_csv(tab / "model_volume_scenario_summary.csv")
    grid = pd.read_parquet(run_all_dir / "processed" / "grid_0p01d.parquet")
    data = M.build_data(cfg, grid, pd.read_csv(tab / "nh4_jumps.csv"))
    start, end = data.t < M.T0_D + 0.1, data.t >= M.END_D - 0.1
    vols = {v.name: v for v in volume.volume_scenarios(cfg, M.END_D) if not v.is_constant}
    for _, r in summ.iterrows():
        v = vols[r["volume_scenario"]]
        assert r["phi_start_median"] == pytest.approx(float(np.median(v.factor(data.t[start]))), rel=1e-12)
        assert r["phi_end_median"] == pytest.approx(float(np.median(v.factor(data.t[end]))), rel=1e-12)
        for w in ("start", "end"):
            assert r[f"C_{w}_ratio_to_constant_volume"] == pytest.approx(
                r[f"C_{w}_median_mgN_L"] / r[f"constant_volume_C_{w}_median_mgN_L"], rel=1e-12)
        assert bool(r["ratio_within_phi_tolerance"]) == all(
            M.within_phi_tolerance(r[f"C_{w}_ratio_to_constant_volume"], r[f"phi_{w}_median"]) for w in ("start", "end"))
    txt = M.volume_scaling_text(summ)
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    point9 = next(line for line in rep.splitlines() if line.startswith("9. **The closed frame"))
    sec4 = rep.split("## 4.")[1].split("## 5.")[0]
    model_md = (run_all_dir / "reports" / "model.md").read_text()
    for text in (point9, sec4, model_md):  # one statement, decided by the table and the tolerance
        assert txt in text
    scaled = bool(summ["ratio_within_phi_tolerance"].all())
    assert txt.startswith("Every scenario C(t) is the constant-volume C(t) scaled") == scaled
    assert txt.startswith("The scenario C(t) is not the constant-volume C(t) scaled") == (not scaled)


def test_real_bound_flags_are_named_in_the_report(run_all_dir):
    tab = run_all_dir / "tables"
    fit = pd.read_csv(tab / "model_fit.csv")
    st = pd.read_csv(tab / "model_calibration_starts.csv")
    assert "cal_at_bounds" in fit.columns and "cal_at_bounds" in st.columns
    for _, r in fit.iterrows():  # the calibration's flag is its chosen start's flag
        ch = st[(st["model"] == r["model"]) & st["chosen"]].iloc[0]
        assert str(ch["cal_at_bounds"]) == str(r["cal_at_bounds"])
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    sec4 = rep.split("## 4.")[1].split("## 5.")[0]
    model_md = (run_all_dir / "reports" / "model.md").read_text()
    txt = M.bounds_text(fit)
    assert txt in sec4 and txt in model_md
    for m in fit.loc[fit["cal_at_bounds"].fillna("").astype(str) != "", "model"]:
        assert f"{m} (" in txt
    assert "Jacobian condition number" in model_md or "condition number 4.2e6" in model_md


def test_real_point7_is_short_and_the_ranges_are_in_the_section4_table(run_all_dir):
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    point7 = next(line for line in rep.splitlines() if line.startswith("7. **Mechanism.**"))
    assert len(point7.split()) <= 160  # "about 150 words" (registered in docs/progress.md, pass 5)
    fit = pd.read_csv(run_all_dir / "tables" / "model_fit.csv")
    sel = fit.loc[fit["selected"]].iloc[0]
    assert f"**{sel['model']}** (selection skill {sel['sel_skill_vs_null']:.3f}" in point7
    assert "holds against every converged start" in point7 or "does not hold against" in point7
    for r in R.other_mechanism(fit):
        assert f"{r['model']} (test {R._fmt_range(*R._range(r, 'test'))})" in point7
    assert "time is confounded" in point7
    sec4 = rep.split("## 4.")[1].split("## 5.")[0]
    assert R.md_table(R.mechanism_table(fit), digits=3) in sec4
