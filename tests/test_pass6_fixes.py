"""Pass 6 (2026-10-26): G1 registration hash, G2 null coverage of the block bootstrap, G3/G4/G6/G7 SCENARIO and
ASSUMPTION labels, G8 provenance columns of the M6 and W2-consistency tables."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from bnr_neq import states as ST

ROOT = Path(__file__).resolve().parents[1]


# --- G1: registration -------------------------------------------------------------------------------------------
def _registered_hash() -> str:
    text = (ROOT / ST.REGISTRATION_DOC).read_text(encoding="utf-8")
    m = re.search(r"^registered_sha256: ([0-9a-f]{64})$", text, flags=re.MULTILINE)
    if m is None:
        pytest.fail(f"{ST.REGISTRATION_DOC} has no registered_sha256 line")
    return m.group(1)


def test_code_matches_the_registered_rules():
    assert ST.registration_sha256() == _registered_hash()


@pytest.mark.parametrize("name, value", [("SLOPE_HALF_BINS", 49), ("ORP_LEVEL_MV", 14), ("N_BOOT", 998),
                                         ("MIN_BLOCKS_PER_GROUP", 1), ("NULL_COVERAGE_MIN_OK", 0.85)])
def test_any_change_to_a_registered_parameter_changes_the_hash(monkeypatch, name, value):
    before = ST.registration_sha256()
    monkeypatch.setattr(ST, name, value)
    assert ST.registration_sha256() != before


def test_a_changed_rule_changes_the_hash(monkeypatch):
    before = ST.registration_sha256()
    monkeypatch.setattr(ST, "RULES", (*ST.RULES[:2], ST.Rule("inactive_oversupplied", "changed", ("ph",))))
    assert ST.registration_sha256() != before


def test_registration_status_says_not_pre_registered():
    assert ST.REGISTRATION_STATUS.startswith("REVISED AFTER FIRST RESULTS (not pre-registered)")
    rules = ST.rules_table(set())
    assert (rules["registration_sha256"] == ST.registration_sha256()).all()
    assert (rules["registration_status"] == ST.REGISTRATION_STATUS).all()


# --- G2: null coverage of the block bootstrap -------------------------------------------------------------------
def test_null_response_has_the_overlap_dependence_and_no_group_effect():
    y = ST.null_response(3000, 0.0, seed=1)
    assert len(y) == 3000
    assert np.corrcoef(y[:-1], y[1:])[0, 1] > 0.95  # neighbouring values share 100 of 101 points
    assert abs(y.mean()) < 4 * y.std() / np.sqrt(3000 / ST.acf_block_length(y, 1))


def test_null_coverage_counts_the_replicates_whose_interval_covers_zero():
    n = 2000
    a = (np.arange(n) // 250) % 2 == 0  # long runs, as the real labels
    rows = ST.bootstrap_null_coverage(a, ~a, phis=(0.0,), n_rep=40, n_boot=199, seed=3)
    (r,) = rows
    assert r["n_replicates"] == 40 and r["n_computed"] == 40
    assert 0.0 <= r["coverage_of_zero"] <= 1.0
    assert r["block_bins_median"] >= 2 * ST.SLOPE_HALF_BINS + 1
    # determinism (fixed seeds)
    assert ST.bootstrap_null_coverage(a, ~a, phis=(0.0,), n_rep=40, n_boot=199, seed=3) == rows


def test_small_groups_give_no_computed_replicate():
    a = np.zeros(1000, dtype=bool)
    a[:50] = True  # fewer than 2 blocks of 101 bins
    (r,) = ST.bootstrap_null_coverage(a, ~a, phis=(0.5,), n_rep=5, n_boot=99)
    assert r["n_computed"] == 0 and np.isnan(r["coverage_of_zero"])


def _comp(status="computed"):
    import pandas as pd

    return pd.DataFrame([{"comparison": "x (pooled)", "status": status}])


def _cov(values):
    import pandas as pd

    return pd.DataFrame({"comparison": "x (pooled)", "coverage_of_zero": values})


@pytest.mark.parametrize("values, start", [([0.97, 0.95], "adequate"), ([0.95, ST.NULL_COVERAGE_MIN_OK], "adequate"),
                                           ([0.95, np.nextafter(ST.NULL_COVERAGE_MIN_OK, 0)], "UNDER-COVERS"),
                                           ([0.6, 0.95], "UNDER-COVERS")])
def test_interval_reading_is_decided_by_the_worst_phi(values, start):
    out = ST.interval_reading(_comp(), _cov(values))
    assert out["null_coverage_min"].iloc[0] == min(values)
    assert out["interval_reading"].iloc[0].startswith(start)


def test_a_comparison_without_an_interval_is_not_read():
    out = ST.interval_reading(_comp("not computed: small"), _cov([0.5]))
    assert out["interval_reading"].iloc[0] == "no interval"


# --- G6, G7: fixed constants and scope of the M6 candidates --------------------------------------------------------
from bnr_neq import model as M
from bnr_neq import pulses as P

FIXTURES = ROOT / "tests" / "fixtures" / "pass4"


def _specs(cfg):
    return {s.name: s for s in M.model_specs(cfg)}


def test_candidate_assumptions_name_every_fixed_constant(cfg):
    sp = _specs(cfg)
    assert M.candidate_assumptions(sp["null_no_consumption"]) == M.CONSTANT_VOLUME_ASSUMPTION
    bd = M.candidate_assumptions(sp["biomass_decay"])
    assert M.K_NH_ASSUMPTION in bd and "K_ALK" not in bd and bd.endswith(M.CONSTANT_VOLUME_ASSUMPTION)
    alk = M.candidate_assumptions(sp["alkalinity_limited"])
    assert M.K_NH_ASSUMPTION in alk and M.ALK_ASSUMPTION in alk
    assert "not an inorganic-carbon limitation" in alk and "no CO2 exchange" in alk
    for T in cfg.get("scenarios.temperature_C"):
        h = M.candidate_assumptions(sp[f"haldane_FA_SCENARIO_T{T:g}C"])
        assert f"SCENARIO {T:g} °C" in h and "K_NH" not in h  # the Haldane law has no K_NH
        assert M.ALK_ASSUMPTION in M.candidate_assumptions(sp[f"alkalinity_FA_SCENARIO_T{T:g}C"])


def test_a_volume_scenario_replaces_the_constant_volume_assumption(cfg):
    a = M.candidate_assumptions(_specs(cfg)["biomass_decay"], "SCENARIO evaporation 1 mL/d")
    assert a.endswith("SCENARIO evaporation 1 mL/d") and M.CONSTANT_VOLUME_ASSUMPTION not in a


# --- G8: provenance of the M6 tables --------------------------------------------------------------------------------
def _fake_result(cfg):
    """A ModelResult from committed M6 tables (no fit): enough for the provenance layer, which changes no value."""
    import pandas as pd

    t = ROOT / "outputs" / "tables"
    read = lambda n: pd.read_csv(t / n).drop(columns=["nh4_observation_model", "nh4_series", "assumptions",
                                                      "constant_volume_nh4_series",
                                                      "model_response_fraction_nh4_observation_model",
                                                      "model_response_fraction_nh4_series",
                                                      "observed_response_fraction_nh4_observation_model"],
                                             errors="ignore")
    fit = read("model_fit.csv")
    sel = fit.loc[fit["selected"].astype(bool), "model"].iloc[0]
    truth = M.Truth(pd.DataFrame(), read("nh4_truth_vectors.csv"), pd.read_csv(t / "nh4_truth_summary.csv"))
    return M.ModelResult(fit, read("model_identifiability.csv"), read("model_checks.csv"),
                         read("model_response_fractions.csv"), pd.DataFrame(), sel, truth=truth,
                         starts=read("model_calibration_starts.csv"), volume_candidates=read("model_volume_scenarios.csv"),
                         volume_summary=read("model_volume_scenario_summary.csv"))


def test_every_m6_table_gets_provenance_and_keeps_its_values(cfg):
    r = _fake_result(cfg)
    tabs = M.m6_tables(cfg, r)
    assert set(tabs) == set(M.M6_PROVENANCE_TABLES)
    sp = _specs(cfg)
    originals = {"model_fit.csv": r.comparison, "model_calibration_starts.csv": r.starts,
                 "model_identifiability.csv": r.identifiability, "model_checks.csv": r.checks,
                 "model_response_fractions.csv": r.responses, "nh4_truth_vectors.csv": r.truth.vectors,
                 "model_volume_scenarios.csv": r.volume_candidates, "model_volume_scenario_summary.csv": r.volume_summary}
    for name, df in tabs.items():
        M.check_m6_provenance(name, df)
        import pandas as pd

        pd.testing.assert_frame_equal(df[list(originals[name].columns)], originals[name])  # nothing else changes
    fit = tabs["model_fit.csv"]
    assert (fit["nh4_series"] == M.MODEL_STATE_PREFIX + fit["model"]).all()
    assert all(a == M.candidate_assumptions(sp[m]) for m, a in zip(fit["model"], fit["assumptions"]))
    vc = tabs["model_volume_scenarios.csv"]
    assert all(s == M.model_state_series(m, v) for s, m, v in zip(vc["nh4_series"], vc["model"], vc["volume_scenario"]))
    assert all(a.endswith(v) for a, v in zip(vc["assumptions"], vc["volume_scenario"]))
    for name in ("model_identifiability.csv", "model_checks.csv", "nh4_truth_vectors.csv"):
        assert (tabs[name]["nh4_series"] == M.model_state_series(r.selected)).all()
    resp = tabs["model_response_fractions.csv"]
    assert (resp["observed_response_fraction_nh4_observation_model"] == M.RAW_READING).all()
    summ = tabs["model_volume_scenario_summary.csv"]
    tid = r.truth.summary.set_index("key")["value"]["nh4_truth_id"]
    assert (summ["constant_volume_nh4_series"] == tid).all()


@pytest.mark.parametrize("name, col, value", [("model_fit.csv", "nh4_observation_model", M.RAW_READING),
                                              ("model_fit.csv", "nh4_series", "C_inv:abc"),
                                              ("model_fit.csv", "assumptions", ""),
                                              ("model_response_fractions.csv",
                                               "observed_response_fraction_nh4_observation_model", M.OBSERVATION_MODEL)])
def test_a_wrong_or_missing_provenance_raises(cfg, name, col, value):
    df = M.m6_tables(cfg, _fake_result(cfg))[name].copy()
    df.loc[df.index[0], col] = value
    with pytest.raises(M.ObservationModelError):
        M.check_m6_provenance(name, df)
    with pytest.raises(M.ObservationModelError):
        M.check_m6_provenance(name, df.drop(columns=[col]))


# --- G3, G4: SCENARIO and ASSUMPTION labels of pulses.csv ----------------------------------------------------------
def test_row_assumptions_per_kind_of_metric(cfg):
    sel = _specs(cfg)["biomass_decay"]
    assert P.row_assumptions("t_pulse_d", "", sel) == "none beyond the metric name"
    assert P.row_assumptions("fa_median_SENSITIVITY_reading_as_TAN_mgNH3_L_SCENARIO_T25C", "", sel) == (
        "none beyond the metric name")
    for metric in ("fa_median_mgNH3_L_SCENARIO_T25C", "delta_A_kJ_mol_SCENARIO_T25C",
                   "delta_A_SENSITIVITY_reading_as_TAN_kJ_mol_SCENARIO_T25C"):  # C(t) or dose/V0: volume-dependent
        assert P.row_assumptions(metric, "", sel) == M.CONSTANT_VOLUME_ASSUMPTION
    assert P.row_assumptions("eta_conversion_fraction", "", sel) == M.CONSTANT_VOLUME_ASSUMPTION
    evap = P.row_assumptions("eta_conversion_fraction_SCENARIO_evaporation_0.5_mL_per_d", "", sel)
    assert evap.startswith("SCENARIO: evaporation 0.5 mL/d") and M.CONSTANT_VOLUME_ASSUMPTION not in evap
    idle = P.row_assumptions(P.IDLE_LOWER, "", sel)
    assert ST.ORP_THRESHOLD_SCENARIO in idle and M.CONSTANT_VOLUME_ASSUMPTION in idle
    mb = P.row_assumptions("eta_conversion_fraction_MODEL_BASED", "MODEL-BASED η", sel)
    assert M.K_NH_ASSUMPTION in mb and "selected model biomass_decay" in mb
    assert M.K_NH_ASSUMPTION not in P.row_assumptions("o2_required_mgO2", "not MODEL-BASED; C_inv:x", sel)


def test_every_committed_pulse_row_gets_its_labels(cfg):
    import pandas as pd

    pm = pd.read_csv(ROOT / "outputs" / "tables" / "pulses.csv").drop(columns=["assumptions"], errors="ignore")
    out = P.label_assumptions(pm, _specs(cfg)["biomass_decay"])
    assert (out["assumptions"].fillna("") != "").all()
    volume_dependent = ~out["metric"].str.match(P._VOLUME_FREE_RE) & ~out["metric"].str.contains("evaporation")
    assert out.loc[volume_dependent, "assumptions"].str.contains(M.CONSTANT_VOLUME_ASSUMPTION, regex=False).all()
    # independent of the regex: every row computed from the inverted C(t) or the model state names its volume frame
    from_c = (out["nh4_series"].fillna("").str.startswith(("C_inv:", M.MODEL_STATE_PREFIX))
              & (out["metric"] != "o2_windows_computed"))  # a count of windows, not a value of C
    volume_label = out["assumptions"].str.contains("constant volume V0|SCENARIO: evaporation")
    assert from_c.sum() > 0 and volume_label[from_c].all()
    assert out.loc[out["metric"].isin(P.STATE_METRICS), "assumptions"].str.contains("SCENARIO 25 °C").all()


# --- G8: provenance of the W2-consistency table ---------------------------------------------------------------------
def test_w2_consistency_with_two_series_raises():
    import pandas as pd

    cons = pd.read_csv(ROOT / "outputs" / "tables" / "nh4_truth_w2_consistency.csv")
    cons = cons.assign(nh4_observation_model=M.OBSERVATION_MODEL, nh4_truth_id="C_inv:a",
                       assumptions=P.W2_CONSISTENCY_ASSUMPTIONS)
    kw = {"label_cols": ("verdict",), "series_col": "nh4_truth_id", "independent_cols": P.W2_CONSISTENCY_INDEPENDENT,
          "derived_cols": P.W2_CONSISTENCY_DERIVED}
    assert M.check_single_observation_model(cons, **kw) == M.OBSERVATION_MODEL
    cons.loc[0, "nh4_truth_id"] = "C_inv:b"
    with pytest.raises(M.ObservationModelError):
        M.check_single_observation_model(cons, **kw)


# --- real record (run-all outputs) ----------------------------------------------------------------------------------
def test_real_m6_tables_carry_provenance_and_assumptions(run_all_dir, cfg):
    import pandas as pd

    tab = run_all_dir / "tables"
    sp = _specs(cfg)
    for name in M.M6_PROVENANCE_TABLES:
        M.check_m6_provenance(name, pd.read_csv(tab / name))
    fit = pd.read_csv(tab / "model_fit.csv")
    assert all(a == M.candidate_assumptions(sp[m]) for m, a in zip(fit["model"], fit["assumptions"]))
    for m in ("alkalinity_limited", *(f"alkalinity_FA_SCENARIO_T{T:g}C" for T in cfg.get("scenarios.temperature_C"))):
        assert "not an inorganic-carbon limitation" in fit.set_index("model").loc[m, "assumptions"]
    tid = pd.read_csv(tab / "nh4_truth_summary.csv").set_index("key")["value"]["nh4_truth_id"]
    assert (pd.read_csv(tab / "model_volume_scenario_summary.csv")["constant_volume_nh4_series"] == tid).all()


def test_real_w2_consistency_names_the_one_inverted_series(run_all_dir):
    import pandas as pd

    tab = run_all_dir / "tables"
    cons = pd.read_csv(tab / "nh4_truth_w2_consistency.csv")
    assert M.check_single_observation_model(cons, label_cols=("verdict",), series_col="nh4_truth_id",
                                            independent_cols=P.W2_CONSISTENCY_INDEPENDENT,
                                            derived_cols=P.W2_CONSISTENCY_DERIVED) == M.OBSERVATION_MODEL
    tid = pd.read_csv(tab / "nh4_truth_summary.csv").set_index("key")["value"]["nh4_truth_id"]
    assert (cons["nh4_truth_id"] == tid).all()
    assert (cons["assumptions"] == P.W2_CONSISTENCY_ASSUMPTIONS).all()


def test_real_pulse_rows_carry_their_assumptions(run_all_dir, cfg):
    import pandas as pd

    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv")
    sel = pd.read_csv(run_all_dir / "tables" / "model_fit.csv").query("selected")["model"].iloc[0]
    expected = P.label_assumptions(pm.drop(columns=["assumptions"]), _specs(cfg)[sel])["assumptions"]
    assert (pm["assumptions"] == expected).all() and (pm["assumptions"].fillna("") != "").all()
    assert pm.loc[pm["metric"].isin(P.STATE_METRICS), "assumptions"].str.contains("SCENARIO 25 °C").all()


def test_real_states_are_registered_and_every_interval_is_read_against_the_null(run_all_dir):
    import pandas as pd

    tab = run_all_dir / "tables"
    h = _registered_hash()
    rules, comp = pd.read_csv(tab / "states_rules.csv"), pd.read_csv(tab / "states_comparison.csv")
    null = pd.read_csv(tab / "states_bootstrap_null_coverage.csv")
    for df in (rules, comp):
        assert (df["registration_sha256"] == h).all() and (df["registration_status"] == ST.REGISTRATION_STATUS).all()
    computed = comp[comp["status"] == "computed"]
    assert set(null["comparison"]) == set(computed["comparison"])
    assert (null.groupby("comparison").size() == len(ST.NULL_COVERAGE_PHIS)).all()
    assert (null["n_replicates"] == ST.NULL_COVERAGE_REPLICATES).all()
    worst = null.groupby("comparison")["coverage_of_zero"].min()
    for r in computed.itertuples(index=False):
        assert r.null_coverage_min == pytest.approx(worst[r.comparison])
        assert r.interval_reading.startswith("adequate" if worst[r.comparison] >= ST.NULL_COVERAGE_MIN_OK
                                             else "UNDER-COVERS")
    assert (comp.loc[comp["status"] != "computed", "interval_reading"] == "no interval").all()
    rep = (run_all_dir / "reports" / "analysis_report.md").read_text()
    states_md = (run_all_dir / "reports" / "states.md").read_text()
    for text in (rep.split("## 3.")[1].split("## 4.")[0], states_md):
        assert ST.REGISTRATION_STATUS in text and h in text
    assert "Pre-registered rules" not in rep and "pre-registered rules" not in rep
