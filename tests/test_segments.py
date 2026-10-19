import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp

from bnr_neq import segments as S
from bnr_neq.constants import ANTHONISEN_EXP_K, NH3_PER_N

RNG = np.random.default_rng(7)


# --- segmentation -------------------------------------------------------------
def _jumps(times, classes):
    """Located, sustained jumps in the M2 candidate-table layout."""
    return pd.DataFrame({"t_jump_d": times, "t_step_peak_d": times, "step_peak": 10.0, "class": classes,
                         "sustained": True, "jump": True, "post_level": 30.0})


def test_boundaries_and_governing_pulse():
    b = S.boundaries([8.53, 10.57], _jumps([2.6, 8.86], ["electrode_offset", "pulse_response"]), 1.0, 12.0)
    assert b["t_d"].tolist() == [1.0, 2.6, 8.53, 8.86, 10.57, 12.0]
    segs = S.sub_segments(b, [8.53, 10.57])
    assert segs["segment"].tolist() == [0, 0, 1, 1, 2]
    assert segs["governing_pulse"].isna().tolist() == [True, True, False, False, False]
    assert (segs["fit_start_d"] - segs["t_start_d"]).round(9).eq(0.1).all()
    assert segs.loc[2, "status"].startswith("not computed: fit span")  # 8.53-8.86: 0.23 d after exclusion


def test_theil_sen_recovers_slope():
    t = np.linspace(0, 5, 501)
    y = 40 - 0.5 * t + 0.05 * RNG.standard_normal(len(t))
    slope, lo, hi = S.theil_sen(t, y)
    assert lo < slope < hi and slope == pytest.approx(-0.5, abs=0.01)


def test_bins_in_excludes_straddling_bins():
    g = pd.DataFrame({"t_start_d": [0.0, 0.01, 0.02], "t_end_d": [0.01, 0.02, 0.03]})
    assert S.bins_in(g, 0.005, 0.03)["t_start_d"].tolist() == [0.01, 0.02]


# --- production results on the real data (gate) -----------------------------
def test_no_slope_window_spans_a_pulse_or_jump(segmented, cfg, validated):
    w = segmented.windows
    assert len(w) > 0 and w["clean"].all()
    events = np.r_[cfg.pulse_times_d, validated.jumps.loc[validated.jumps["jump"], "t_jump_d"]]
    for r in segmented.table[segmented.table["status"] == "computed"].itertuples():
        assert not ((events > r.fit_start_d) & (events < r.fit_end_d)).any()
        assert r.fit_start_d - r.t_start_d == pytest.approx(0.1)


def test_no_slope_window_contains_any_step_candidate(segmented, validated):
    """Every M2 candidate with |step| >= threshold, including unlocated and non-sustained ones (e.g. near
    3.69, 5.89 and 9.29 d), is outside every computed slope window, run and step peak included."""
    j = validated.jumps
    cand = j[j["step_peak"].abs() >= S.qc.JUMP_THRESHOLD_MG_L]
    assert (~cand["jump"].astype(bool)).sum() >= 3
    for t in (3.69, 5.89, 9.29):
        assert ((cand["t_step_peak_d"] - t).abs() < 0.03).any()
    comp = segmented.table[segmented.table["status"] == "computed"]
    for c in cand[~cand["jump"].astype(bool)].itertuples():
        for r in comp.itertuples():
            assert not (c.run_end_d > r.fit_start_d and c.run_start_d < r.fit_end_d), (c.t_step_peak_d, r.sub_segment)


def _candidates(rows):
    cols = ["run_start_d", "run_end_d", "t_step_peak_d", "step_peak", "post_level", "t10_d", "t90_d", "t_jump_d",
            "sustained", "class", "jump"]
    return pd.DataFrame(rows, columns=cols)


def test_window_check_catches_unlocated_step_inside_a_window():
    table = pd.DataFrame({"sub_segment": [1], "status": ["computed"], "fit_start_d": [1.1], "fit_end_d": [3.0],
                          "held_fraction": [0.0]})
    unlocated = _candidates([(1.98, 2.02, 2.0, 6.0, 30.0, np.nan, np.nan, np.nan, True, "unlocated", False)])
    w = S.windows_check(table, [8.53], unlocated)
    assert not w["clean"].iloc[0] and "unlocated" in w["events_inside"].iloc[0]
    # a non-sustained candidate whose run overlaps the window edge is caught too
    edge = _candidates([(2.95, 3.05, 3.02, -5.0, 30.0, np.nan, np.nan, np.nan, False, "unlocated", False)])
    assert not S.windows_check(table, [8.53], edge)["clean"].iloc[0]
    # below the detector threshold: not a candidate
    small = _candidates([(1.98, 2.02, 2.0, 1.0, 30.0, np.nan, np.nan, np.nan, True, "unlocated", False)])
    assert S.windows_check(table, [8.53], small)["clean"].iloc[0]


def test_boundaries_split_around_unlocated_candidate():
    unlocated = _candidates([(3.69, 3.715, 3.705, 4.7, 37.9, np.nan, np.nan, np.nan, True, "unlocated", False)])
    b = S.boundaries([8.53], unlocated, 1.0, 9.0)
    row = b[b["kind"].str.startswith("candidate")].iloc[0]
    assert row["t_d"] == pytest.approx(3.69) and row["t_resume_d"] == pytest.approx(3.715)
    segs = S.sub_segments(b, [8.53])
    assert segs["t_end_d"].iloc[0] == pytest.approx(3.69) and segs["t_start_d"].iloc[1] == pytest.approx(3.715)


def _synthetic_grid(t0, t1, nh4_fn):
    t = np.round(np.arange(t0, t1, 0.01), 10) + 0.005
    y = nh4_fn(t)
    return pd.DataFrame({"t_start_d": t - 0.005, "t_end_d": t + 0.005, "t_mid_d": t, "pH": 7.3,
                         "nh4_mgNH4_L": y, "nh4_mgN_L": y * 0.77649})


NO_JUMPS = _candidates([])
NO_HELD = pd.DataFrame({"t_start_d": [], "t_end_d": []})


def test_negative_consumption_not_computed(cfg):
    p1 = cfg.pulse_times_d[0]
    grid = _synthetic_grid(p1 - 1.0, p1 + 1.5, lambda t: 40.0 + np.where(t > p1, 2.0 * (t - p1), 0.0))  # rising
    pt = pd.DataFrame({"pulse": [1, 2, 3, 4], "nh4_response_fraction": [0.37, 0.29, 0.26, 0.17]})
    tab = S.slope_table(cfg, grid, NO_JUMPS, pt, NO_HELD)
    r = tab[(tab["governing_pulse"] == 1) & (tab["status"] == "computed")].iloc[0]
    assert r["slope_mgNH4_L_d"] > 0
    assert r["gain_correction_status"] == S.NON_POSITIVE and r["gain_correction_status"].startswith("not computed: non-positive")
    assert np.isnan(r["consumption_gain_corrected_mgN_L_d"]) and np.isnan(r["consumption_gain_corrected_ci95_low_mgN_L_d"])
    falling = S.slope_table(cfg, _synthetic_grid(p1 - 1.0, p1 + 1.5, lambda t: 40.0 - 2.0 * (t - p1)), NO_JUMPS, pt, NO_HELD)
    f = falling[(falling["governing_pulse"] == 1) & (falling["status"] == "computed")].iloc[0]
    assert f["gain_correction_status"] == "computed"
    assert f["consumption_gain_corrected_mgN_L_d"] == pytest.approx(2.0 * 0.77649 / 0.37, rel=1e-6)


def test_held_span_flagged_and_excluded_in_sensitivity(cfg):
    p1 = cfg.pulse_times_d[0]
    # falling at 2 mg/L/d, but the output is held flat over [p1 + 0.5, p1 + 1.0]
    def nh4(t):
        y = 40.0 - 2.0 * (t - p1)
        held = (t >= p1 + 0.5) & (t < p1 + 1.0)
        return np.where(held, 40.0 - 2.0 * 0.5, np.where(t >= p1 + 1.0, y + 1.0, y))
    grid = _synthetic_grid(p1, p1 + 1.6, nh4)
    held = pd.DataFrame({"t_start_d": [p1 + 0.5], "t_end_d": [p1 + 1.0]})
    pt = pd.DataFrame({"pulse": [1, 2, 3, 4], "nh4_response_fraction": [0.37, 0.29, 0.26, 0.17]})
    tab = S.slope_table(cfg, grid, NO_JUMPS, pt, held)
    r = tab[(tab["governing_pulse"] == 1) & (tab["status"] == "computed")].iloc[0]
    span = r["fit_end_d"] - r["fit_start_d"]
    assert r["held_overlap_d"] == pytest.approx(0.5) and r["held_fraction"] == pytest.approx(0.5 / span)
    assert r["slope_excl_held_status"] == "computed"
    b = S.bins_outside_held(S.bins_in(grid, r["fit_start_d"], r["fit_end_d"]), held)
    assert not ((b["t_end_d"] > p1 + 0.5 + 1e-9) & (b["t_start_d"] < p1 + 1.0 - 1e-9)).any()
    assert r["slope_excl_held_mgNH4_L_d"] == pytest.approx(S.theil_sen(b["t_mid_d"].to_numpy(), b["nh4_mgNH4_L"].to_numpy())[0])


def test_real_held_fraction_reported(segmented, validated):
    t = segmented.table
    comp = t[t["status"] == "computed"]
    assert comp["held_fraction"].between(0, 1).all()
    for r in comp.itertuples():
        assert r.held_overlap_d == pytest.approx(S.held_overlap_d(r.fit_start_d, r.fit_end_d, validated.held))
    assert (comp["held_fraction"] > 0.5).sum() >= 2  # 10.98-17.88 and 17.98-24.59 overlap 13-18 and 20.7-24.7 d


def test_pelt_penalty_is_autocorrelation_aware():
    rng = np.random.default_rng(3)
    e = rng.normal(size=4000)
    ar = np.zeros(4000)
    for i in range(1, 4000):
        ar[i] = 0.9 * ar[i - 1] + e[i]
    sc = S.penalty_scales(ar)
    assert sc["phi"] > 0.5
    assert sc["long_run_ar1"] == pytest.approx(sc["innovation_ar1"] / (1 - sc["phi"]))
    assert sc["long_run_ar1"] > sc["marginal"] > sc["innovation_ar1"]
    white = S.penalty_scales(rng.normal(size=4000))
    assert white["long_run_ar1"] == pytest.approx(white["marginal"], rel=0.25)


def test_pelt_sensitivity_reported(segmented):
    sens = segmented.pelt_sensitivity
    assert set(sens["penalty_scale"]) == set(S.PELT_PENALTY_SCALES)
    prod = sens[sens["production"]]
    assert len(prod) == len(S.GRID_SERIES) and (prod["penalty_scale"] == "long_run_ar1").all()
    for series, g in sens[sens["beta"] == S.PELT_BETA].groupby("series"):
        n = g.set_index("penalty_scale")["n_breaks"]
        assert n["innovation_ar1"] >= n["marginal"] >= n["long_run_ar1"], series


def test_baseline_gain_correction_not_computed(segmented):
    base = segmented.table[segmented.table["segment"] == 0]
    assert (base["gain_correction_status"] == "not computed: no governing pulse (baseline)").all()
    assert base["consumption_gain_corrected_mgN_L_d"].isna().all()


def test_scenario_labels(segmented):
    fa = [c for c in segmented.table.columns if c.startswith("fa_median")]
    # W1: reading-as-TAN free ammonia is a SENSITIVITY case, and the label is part of the column name
    assert fa == ["fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L_SCENARIO_T20C",
                  "fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L_SCENARIO_T25C"]
    assert segmented.fits["model"].str.contains("SCENARIO").sum() == 2 * segmented.fits["segment"].nunique()


def test_exact_reference_windows_reproduce_all_slopes(segmented):
    assert segmented.reference_windows["exact_window_pass"].all()


@pytest.mark.parametrize("ref_from", [10.85, 17.84, 24.60, 32.36, 37.19, 44.63])
def test_production_slopes_reproduce_references(segmented, ref_from):
    r = segmented.reference_windows.set_index("ref_from_d").loc[ref_from]
    assert r["production_pass"]


def test_unreproducible_reference_windows_are_the_explained_ones(segmented):
    rw = segmented.reference_windows
    failing = rw.loc[~rw["production_pass"], "ref_from_d"].tolist()
    assert failing == [0.3, 2.59]
    assert rw.set_index("ref_from_d").loc[2.59, "exact_window_spans_pulse"]


@pytest.mark.parametrize("label", ["after pulse 2", "after pulse 4", "free ammonia min", "free ammonia max",
                                   "Lomb-Scargle"])
def test_other_references(segmented, label):
    r = segmented.references[segmented.references["quantity"].str.contains(label)]
    assert len(r) == 1 and r["pass"].iloc[0], r.to_string()


# --- kinetic solvers ------------------------------------------------------------
def test_monod_exact_solution_matches_ode():
    t = np.linspace(0, 5, 101)
    m = S.KineticModel("monod", ("vmax_mgN_L_d", "Ks_mgN_L"))
    exact = S.simulate(m, np.array([3.0, 8.0]), 40.0, t, np.full_like(t, 7.0))
    ode = solve_ivp(lambda _, c: -3.0 * c / (8.0 + c), (0, 5), [40.0], t_eval=t, rtol=1e-10, atol=1e-10).y[0]
    np.testing.assert_allclose(exact, ode, rtol=1e-6)


def test_haldane_rk4_matches_ode():
    t = np.linspace(0, 3, 301)
    ph = 7.2 + 0.3 * np.sin(t)
    T_K = 298.15
    m = S.KineticModel("haldane", ("vmax", "Ks", "Ki"), 25.0)
    rk4 = S.simulate(m, np.array([4.0, 0.3, 1.0]), 50.0, t, ph)

    def rhs(tt, c):
        p = 7.2 + 0.3 * np.sin(tt)
        s = NH3_PER_N * c[0] * 10**p / (np.exp(ANTHONISEN_EXP_K / T_K) + 10**p)
        return [-4.0 * s / (0.3 + s + s * s / 1.0)]

    ode = solve_ivp(rhs, (0, 3), [50.0], t_eval=t, rtol=1e-10, atol=1e-10).y[0]
    np.testing.assert_allclose(rk4, ode, rtol=1e-5)


def _pieces(model, kin, levels, noise=0.02):
    out = []
    for k, c0 in enumerate(levels):
        t = 2.0 * k + np.arange(0, 1.5, 0.01)
        ph = np.full_like(t, 7.3)
        c = S.simulate(model, np.asarray(kin), c0, t, ph) + noise * RNG.standard_normal(len(t))
        out.append((t, c, ph))
    return out


def test_fit_recovers_zero_and_first_order():
    for name, par, kin in (("zero_order", ("k0_mgN_L_d",), [2.0]), ("first_order", ("k1_per_d",), [0.05])):
        m = S.KineticModel(name, par)
        f = S.fit_segment(m, _pieces(m, kin, [40.0, 25.0]), -1.0)
        assert f["converged"] and f["params"][0]["estimate"] == pytest.approx(kin[0], rel=0.03)
        assert f["params"][0]["identifiable"]


def test_monod_half_saturation_not_identifiable_when_saturated():
    m = S.KineticModel("monod", ("vmax_mgN_L_d", "Ks_mgN_L"))
    f = S.fit_segment(m, _pieces(m, [2.0, 0.01], [40.0, 30.0]), -2.0)
    ks = next(p for p in f["params"] if p["parameter"] == "Ks_mgN_L")
    assert not ks["identifiable"] and ks["reason"]


# --- change points and periodicity ---------------------------------------------
def _grid(n=2000, f=None):
    t = np.arange(n) * 0.01 + 0.005
    base = {c: 0.01 * RNG.standard_normal(n) for c in S.GRID_SERIES}
    g = pd.DataFrame({"t_start_d": t - 0.005, "t_end_d": t + 0.005, "t_mid_d": t, **base})
    return g, t


def test_pelt_finds_step_but_not_drift():
    g, t = _grid()
    g["pH"] += 0.2 * t + np.where(t >= 12.0, 0.5, 0.0)  # drift everywhere, one step
    cp, _ = S.change_points(g, pulse_times_d=[5.0])
    ev = cp[(cp["series"] == "pH") & cp["event"]]
    assert len(ev) == 1 and ev["t_break_d"].iloc[0] == pytest.approx(12.0, abs=0.06)
    assert pd.isna(ev["tied_pulse"].iloc[0])


def test_tied_pulse_rule():
    assert S.tied_pulse(8.6, [8.53]) == 1
    assert S.tied_pulse(8.52, [8.53]) == 1  # within the -0.02 d stamp tolerance
    assert S.tied_pulse(9.2, [8.53]) is None


def test_periodogram_finds_daily_cycle():
    g, t = _grid()
    g["pH"] += 0.2 * np.sin(2 * np.pi * t)
    summary, spectra = S.periodograms(g)
    row = summary[(summary["series"] == "pH") & summary["band"].str.startswith("periodicity")].iloc[0]
    assert row["f_at_max_c_d"] == pytest.approx(1.0, abs=0.01)
    assert row["max_normalized_power"] > 0.9
    assert len(spectra) > 0


def test_load_inputs_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="ingest"):
        S.load_inputs(tmp_path)


def test_run_writes_outputs(cfg, ingested, tmp_path):
    from bnr_neq import pipeline, validation

    pipeline.write_ingest(ingested, tmp_path)
    validation.run(cfg, tmp_path)
    S.run(cfg, tmp_path)
    for rel in ("tables/segments.csv", "figures/fig2_overview.png", "figures/fig2_overview.svg",
                "reports/segments.md", "tables/kinetic_fits.csv", "tables/changepoints.csv"):
        assert (tmp_path / rel).is_file()
