import numpy as np
import pandas as pd
import pytest

from bnr_neq import dynamics as D
from bnr_neq import estimators as E


def test_validation_required_checks_pass(estimator_validation):
    v = estimator_validation
    assert v["required"].sum() >= 20
    assert v.loc[v["required"], "pass"].all(), v.loc[v["required"] & ~v["pass"], "check"].tolist()
    markov = v[v["check"] == "markov_ring_entropy_production"]
    assert len(markov) >= 3 and (markov["expected"] == 0.0).any()


def test_quantized_bias_reported(estimator_validation):
    q = estimator_validation[estimator_validation["check"] == "quantized_observable_bias"]
    assert len(q) == 3 and (q["value"] < 0).all()  # misclassification biases the estimate low


def test_validated_tie_ranges(estimator_validation):
    r = D.validated_tie_range(estimator_validation)
    assert all(r[name] >= 0.9 for name in ("ordinal_kl", "asym_norm_lag1", "hvg_kl"))


def test_gate_blocks_plant_data(cfg, dynamics_inputs, estimator_validation):
    bad = estimator_validation.copy()
    bad.loc[bad.index[0], "pass"] = False
    with pytest.raises(D.EstimatorValidationError):
        D.compute(cfg, dynamics_inputs, bad)


def test_stationarity_raw_only():
    rng = np.random.default_rng(0)
    assert D.stationarity(rng.normal(size=2000))["transform"] == "raw"
    rw = D.stationarity(np.cumsum(rng.normal(size=2000)))
    assert rw["transform"].startswith("not computed: non-stationary")
    trend = D.stationarity(np.arange(2000) * 0.01 + rng.normal(size=2000))
    assert trend["transform"].startswith("not computed")  # no detrending path
    assert D.stationarity(np.ones(100))["transform"] == "not computed: constant window"


def test_tie_range_enforced_per_window():
    # seed 2: a draw that KPSS accepts as stationary (seed 1 is one of the 5 % false KPSS rejections)
    rng = np.random.default_rng(2)
    z = np.round(rng.normal(size=2160) * 0.004 / 0.01) * 0.01
    row = D.analyse_window(z, rng, {"asym_norm_lag1": 1.0, "ordinal_kl": 1.0, "hvg_kl": 0.25})
    assert row["transform"] == "raw" and row["tie_fraction"] > 0.25
    assert row["hvg_kl_status"].startswith("not computed: tie fraction")
    assert "hvg_kl" not in row and row["ordinal_kl_status"] == "computed"


def test_subtract_intervals():
    assert D.subtract_intervals([(0, 10)], [(2, 3), (5, 6)]) == [(0, 2), (3, 5), (6, 10)]
    assert D.subtract_intervals([(0, 1)], [(-1, 2)]) == []


def test_allowed_intervals(cfg, dynamics_inputs):
    orp = D.allowed_intervals(cfg, "ORP", dynamics_inputs["held"], dynamics_inputs["segments"])
    for t in D.ORP_ELECTRODE_EVENTS_D:
        assert not any(a <= t <= b for _, a, b in orp)
    nh4 = D.allowed_intervals(cfg, "NH4", dynamics_inputs["held"], dynamics_inputs["segments"])
    held = dynamics_inputs["held"]
    for _, a, b in nh4:
        assert not any((a < h1) and (b > h0) for h0, h1 in zip(held["t_start_d"], held["t_end_d"]))
    for series in ("pH", "ORP", "NH4", "EC"):
        for _, a, b in D.allowed_intervals(cfg, series, held, dynamics_inputs["segments"]):
            assert a >= D.START_D and not any(a < p < b for p in cfg.pulse_times_d)


def test_data_windows(cfg, dynamics_result):
    w = dynamics_result.windows
    assert set(w["transform"].dropna().unique()) <= {"raw"} | set(w["transform"][w["transform"].str.startswith("not computed")])
    for r in w.itertuples(index=False):
        assert not any(r.t_start_d < p < r.t_end_d for p in cfg.pulse_times_d)
    ranges = D.validated_tie_range(dynamics_result.validation)
    raw = w[w["transform"] == "raw"]
    assert len(raw[raw["series"] == "pH"]) > 0
    for name, limit in ranges.items():
        done = raw[raw[f"{name}_status"] == "computed"]
        skipped = raw[raw[f"{name}_status"] != "computed"]
        tie_skipped = skipped[skipped[f"{name}_status"].str.startswith("not computed: tie fraction")]
        assert (done["tie_fraction"] <= limit).all() and (tie_skipped["tie_fraction"] > limit).all()
        # every other skip is the pH matched-null gate, with the window below the tie limit
        other = skipped.drop(tie_skipped.index)
        assert other[f"{name}_status"].str.startswith(D.OUTSIDE_REGIME).all()
        assert (other["series"] == "pH").all() and (other["tie_fraction"] <= limit).all()
    assert (w.loc[w["series"] == "EC", "role"] == "probe behaviour only (S3)").all()


def test_onsager_uses_pulse_response(dynamics_result, validated):
    s = dynamics_result.onsager_summary
    np.testing.assert_allclose(s["pulse_response"], validated.pulse_table["ph_obs_change"])
    np.testing.assert_allclose(s["snr_pulse_response"], s["pulse_response"].abs() / s["fluctuation_sd"])
    assert s.loc[~s["snr_sufficient"], "verdict"].str.startswith("not interpretable").all()


def test_flux_and_stability_tables(dynamics_result):
    f = dynamics_result.flux
    assert set(f["plane"]) == {"pH-ORP", "pH-NH4"} and len(f) == 10
    assert dynamics_result.stability["note"].str.startswith("EXPLORATORY").all()


def test_electrode_marker(dynamics_result):
    assert dynamics_result.electrode["status"].iloc[0] == "not computed: electrode_step_response missing"


def test_nh4_uniform_bins_not_filled():
    nh4 = pd.DataFrame({"t_d": np.r_[np.arange(0, 0.01, 20 / 86400), 0.02 + np.arange(0, 0.01, 2 / 86400)],
                        "nh4_mgNH4_L": 1.0})
    u = D.nh4_uniform(nh4)
    assert u["value"].eq(1.0).all()
    assert not ((u["t_d"] > 0.0102) & (u["t_d"] < 0.0198)).any()  # the gap stays empty


def test_load_inputs_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        D.load_inputs(tmp_path)


def test_ar1_and_ou_helpers():
    rng = np.random.default_rng(3)
    x = D._ar1(20000, 0.9, rng)
    assert np.corrcoef(x[:-1], x[1:])[0, 1] == pytest.approx(0.9, abs=0.02)
    a, b = D._ou_2d(20000, 0.01, 1.0, rng)
    assert E.circulation_rate(a, b) > 0


# --- matched null for the pH windows (review item 1) --------------------------------------
def test_matched_null_regimes(estimator_validation):
    null = D.null_from_validation(estimator_validation)
    assert D.N_SURR >= 99
    assert len(null) == len(D.NULL_PHI) * len(D.NULL_TIE) * len(D.INDEX_FNS)
    assert set(null["phi"]) == set(D.NULL_PHI) and {0.95, 0.99, 0.995, 0.999} <= set(null["phi"])
    assert {0.8, 0.85, 0.9, 0.95} <= set(null["tie_target"])
    assert (null["n_stationary"] <= D.NULL_TARGET_STATIONARY).all()
    for r in null.itertuples(index=False):
        k = D.fp_threshold(r.n_stationary) if r.n_stationary else 0
        assert r.fp_threshold_count == k
        assert r.validated == (r.n_stationary >= D.NULL_MIN_STATIONARY and r.false_positives <= k)
    # the null reaches the persistent, highly tied pH regime and shows why the gate is needed there
    hi = null[(null["phi"] >= 0.99) & (null["tie_target"] == 0.95) & (null["index"] == "ordinal_kl")]
    assert (hi["tie_max"] >= 0.94).all() and (~hi["validated"]).any()
    rows = estimator_validation[estimator_validation["check"].str.startswith("iaaft_matched_null_false_positives_")]
    assert (~rows["required"].astype(bool)).all() and (rows["role"] == "gate: pH windows").all()


def test_fp_threshold_is_005_plus_binomial_margin():
    from scipy.stats import binom

    assert D.fp_threshold(100) == 9 and D.fp_threshold(60) == 6
    for n in (30, 57, 100):
        k = D.fp_threshold(n)
        assert binom.sf(k, n, 0.05) <= 0.05 < binom.sf(k - 1, n, 0.05)


def test_quantized_ar1_hits_target_tie_fraction():
    rng = np.random.default_rng(11)
    for phi, tie in ((0.95, 0.8), (0.99, 0.94), (0.995, 0.9)):
        ties = [E.tie_fraction(D.quantized_ar1(2160, phi, tie, rng)) for _ in range(20)]
        assert np.mean(ties) == pytest.approx(tie, abs=0.02)


def test_persistent_quantized_ar1_at_tie_094_is_gated(estimator_validation):
    """Reversible, persistent, quantized AR(1) at a tie fraction of 0.94 (new seeds, not the null's): among the
    stationary windows, a reported index keeps <= 5 % false positives (+ binomial margin); every other window
    is marked outside the validated regime."""
    null = D.null_from_validation(estimator_validation)
    ranges = D.validated_tie_range(estimator_validation)
    rng = np.random.default_rng(940)
    rows = []
    for phi in (0.99, 0.995):
        n_stat = 0
        for _ in range(600):
            z = D.quantized_ar1(D.SURR_LENGTH, phi, 0.94, rng)
            if "series" not in D.stationarity(z):
                continue
            rows.append(D.analyse_window(z, rng, ranges, null))
            n_stat += 1
            if n_stat >= 20:
                break
    assert len(rows) >= 25
    w = pd.DataFrame(rows)
    assert w["tie_fraction"].between(0.9, 0.97).all()
    for name in D.INDEX_FNS:
        st = w[f"{name}_status"]
        reported = w[st == "computed"]
        assert (st == "computed").sum() + st.str.startswith(D.OUTSIDE_REGIME).sum() + st.str.startswith(
            "not computed: tie fraction").sum() == len(w)
        if len(reported):
            fp = int((reported[f"{name}_p"] <= 0.05).sum())
            assert fp <= D.fp_threshold(len(reported)), (name, fp, len(reported))


def test_ph_windows_gated_by_matched_null(dynamics_result):
    w = dynamics_result.windows
    ph = w[(w["series"] == "pH") & (w["transform"] == "raw")]
    assert len(ph) > 0 and ph["null_regimes"].notna().all()
    null = D.null_from_validation(dynamics_result.validation)
    for name in D.INDEX_FNS:
        for r in ph.itertuples(index=False):
            status = getattr(r, f"{name}_status")
            ok, _ = D.null_gate(name, D.match_regimes({"tie_fraction": r.tie_fraction, "ac1": r.ac1,
                                                       "dfa_alpha": r.dfa_alpha}, null), null)
            if status == "computed":
                assert ok
            elif not status.startswith("not computed: tie fraction"):
                assert not ok and status.startswith(D.OUTSIDE_REGIME)
    gate = dynamics_result.validation[dynamics_result.validation["check"].str.startswith("matched_null_gate_pH_")]
    assert len(gate) == len(D.INDEX_FNS) and gate["required"].all() and gate["pass"].all()


def test_gate_check_fails_on_a_reported_window_outside_regime(dynamics_result):
    w = dynamics_result.windows.copy()
    i = w[(w["series"] == "pH") & (w["transform"] == "raw")].index[0]
    w.loc[i, ["ac1", "ordinal_kl_status"]] = [-0.5, "computed"]  # anti-persistent: inside no null regime
    gate = D.verify_null_gate(w, dynamics_result.validation).set_index("check")
    assert not gate.loc["matched_null_gate_pH_ordinal_kl", "pass"]
    assert gate.loc["matched_null_gate_pH_ordinal_kl", "value"] >= 1


# --- stability, Onsager and flux inputs (review items 7 and 8) ----------------------------
def test_orp_stability_window_excludes_30_48(cfg, ingested):
    pieces = D.stability_pieces("ORP", 30.39, 32.39)
    assert len(pieces) == 2 and not any(a <= 30.48 <= b for a, b in pieces)
    assert D.stability_pieces("pH", 30.39, 32.39) == [(30.39, 32.39)]
    st = D.stability(cfg, ingested.grid).set_index(["pulse", "series"])
    k = cfg.pulse_times_d.index(32.39) + 1
    row = st.loc[(k, "ORP")]
    assert "30.43-30.58" in row["excluded_d"] and row["n_pieces"] == 2
    assert st.loc[(k, "pH"), "excluded_d"] == ""
    # 28.89 lies before pulse 4's 2-day window; every other pulse window is unaffected
    assert (st.xs("ORP", level="series").drop(index=k)["excluded_d"] == "").all()


def test_onsager_not_computed_on_non_stationary_late_window(cfg):
    rng = np.random.default_rng(5)
    t = np.round(np.arange(1.0, 47.92, 0.01), 10) + 0.005
    grid = pd.DataFrame({"t_mid_d": t, "pH": 7.3 + 0.02 * np.cumsum(rng.normal(size=len(t)))})
    pt = pd.DataFrame({"pulse": [1, 2, 3, 4], "ph_obs_change": [5.0, 5.0, 5.0, 5.0]})  # SNR sufficient
    lags, osum = D.onsager(cfg, grid, pt)
    bad = osum[~osum["late_window_stationary"]]
    assert len(bad) >= 3
    assert bad["verdict"].str.startswith("not computed: non-stationary").all()
    assert lags[lags["pulse"].isin(bad["pulse"])]["ratio_R_over_C"].isna().all()


def test_onsager_and_flux_report_stationarity(dynamics_result):
    o = dynamics_result.onsager_summary
    assert {"late_adf_p", "late_kpss_p", "late_window_stationary"} <= set(o.columns)
    assert o.loc[~o["late_window_stationary"].astype(bool), "verdict"].str.contains("non-stationary").all()
    f = dynamics_result.flux
    assert (f["n_stationary_input_windows"] <= f["n_input_windows"]).all()
    done = pd.to_numeric(f["circulation"], errors="coerce").notna()
    assert (f.loc[done, "n_stationary_input_windows"] >= 1).all()
    assert f.loc[~done, "circulation"].str.startswith("not computed").all()


def test_onsager_conclusion_derived_from_table():
    osum = pd.DataFrame({"pulse": [1, 2], "snr_pulse_response": [0.5, 4.0], "snr_sufficient": [False, True],
                         "late_window_stationary": [True, True],
                         "verdict": ["not interpretable: pulse-response SNR 0.50 < 3.0", "consistent with regression"]})
    assert D.onsager_conclusion(osum).startswith("Interpretable pulses: pulse 2: consistent with regression")
    osum["verdict"] = ["not interpretable: x", "not computed: non-stationary late fluctuation window"]
    assert "No conclusion" in D.onsager_conclusion(osum)


def test_dynamics_run_on_real_outputs(cfg, run_all_dir, estimator_validation, tmp_path):
    import shutil

    out = tmp_path / "out"
    shutil.copytree(run_all_dir, out)
    (out / "tables" / "dynamics.csv").unlink()
    r = D.run(cfg, out, validation=estimator_validation)
    assert (out / "tables" / "dynamics.csv").is_file() and (out / "reports" / "dynamics.md").is_file()
    assert not D.required_failures(r.validation)
    # same seeds, same validation suite: identical to the run-all table
    assert (out / "tables" / "dynamics.csv").read_bytes() == (run_all_dir / "tables" / "dynamics.csv").read_bytes()
    assert (out / "tables" / "estimator_validation.csv").read_bytes() == (
        run_all_dir / "tables" / "estimator_validation.csv").read_bytes()


def test_passes_stationarity_matches_stationarity_decision():
    rng = np.random.default_rng(21)
    for phi, tie in ((0.9, 0.6), (0.99, 0.9), (0.999, 0.95), (0.995, 0.8)):
        for _ in range(15):
            z = D.quantized_ar1(D.SURR_LENGTH, phi, tie, rng)
            assert D.passes_stationarity(z) == ("series" in D.stationarity(z))
    assert not D.passes_stationarity(np.ones(500))
