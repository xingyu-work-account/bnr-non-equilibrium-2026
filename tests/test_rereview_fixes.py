"""Fixes for the 2026-10-19 re-review: S1-S3, W1 series provenance, W5, W8 and R8."""

import numpy as np
import pandas as pd
import pytest
from conftest import config_with

from bnr_neq import alkalinity as A
from bnr_neq import model as M
from bnr_neq import pulses as P
from bnr_neq import segments as S
from bnr_neq import thermo, volume


# --- W1: one series, checked as well as one model ---------------------------------------------------------------
def _rows(*specs):
    return pd.DataFrame([{"metric": m, "basis": b, "nh4_observation_model": o, "nh4_series": s} for m, b, o, s in specs])


def test_two_inverted_series_in_one_output_raise():
    with pytest.raises(M.ObservationModelError, match="inverted C"):
        M.check_single_observation_model(_rows(("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL, "C_inv:aaa"),
                                               ("o2_demand_mgO2_L_d", "M6", M.OBSERVATION_MODEL, "C_inv:bbb")))


def test_model_state_series_must_be_labelled_model_based():
    with pytest.raises(M.ObservationModelError, match="not labelled MODEL-BASED"):
        M.check_single_observation_model(_rows(("eta_uncertainty_interval", "", M.OBSERVATION_MODEL, "M6_state:x")))
    with pytest.raises(M.ObservationModelError, match="MODEL-BASED carry an inverted"):
        M.check_single_observation_model(_rows(("eta_conversion_fraction_MODEL_BASED", "MODEL-BASED", M.OBSERVATION_MODEL,
                                                "C_inv:aaa")))
    ok = _rows(("consumption_weighted_mgN_L_d", "M6", M.OBSERVATION_MODEL, "C_inv:aaa"),
               ("eta_conversion_fraction_MODEL_BASED", "MODEL-BASED", M.OBSERVATION_MODEL, "M6_state:x"))
    assert M.check_single_observation_model(ok) == M.OBSERVATION_MODEL


# --- S1 / W1 (b): the reading-based columns are labelled in their names ----------------------------------------------
def test_reading_based_columns_say_sensitivity(cfg, segmented, validated):
    steps = thermo.affinity_step_table(cfg, validated.windows, validated.pulse_table)
    assert all("SENSITIVITY" in c for c in steps.columns if c.startswith("delta_A"))
    assert (steps["role"].str.startswith("SENSITIVITY")).all()
    assert not any(c.startswith("delta_A_kJ_mol") for c in steps.columns)
    fa = [c for c in segmented.table.columns if c.startswith("fa_median")]
    assert fa and all("SENSITIVITY" in c for c in fa)


# --- S2: the inverted C(t) against the W2 bound ---------------------------------------------------------------------
def _balance(best_end, band_min_end):
    rows = []
    for name, t, best, bmin in (("pre_pulse1", 8.5, 5.0, 2.0), ("pre_pulse2", 10.5, 9.0, 4.0), ("end", 47.9, best_end,
                                                                                               band_min_end)):
        rows.append({"volume_scenario": volume.CONSTANT, "checkpoint": name, "t_d": t, "C_inv_mgN_L": 100.0,
                     "removed_per_initial_litre_mgN_L": 7 * best, "removed_band_low_mgN_L": 1.0,
                     "removed_band_high_mgN_L": 2.0, "alk0_required_if_nitrified_meq_L": best,
                     "alk0_required_if_nitrified_band_min_meq_L": bmin, "alk0_required_any_proton_sink_meq_L": best / 2})
    return pd.DataFrame(rows)


def test_contradiction_is_detected_and_qualifies_the_pulse():
    cons = P.w2_consistency(_balance(20.9, 11.5))
    assert cons["consistent_at_best_fit"].tolist() == [True, True, False]
    assert cons["verdict"].iloc[-1].startswith("QUALIFIED") and "open upper part" in cons["verdict"].iloc[-1]
    assert P.pulse_qualifier(cons, 1, 2).startswith("pre_pulse1: consistent")
    assert P.pulse_qualifier(cons, 2, 2).startswith("QUALIFIED")
    m = pd.DataFrame({"pulse": [1, 2, 2, 2], "metric": ["a", "b", "c", "t_pulse_d"],
                      "nh4_observation_model": [M.OBSERVATION_MODEL] * 3 + [""]})
    q = P.qualify(m, cons, 2, "pulse", "nh4_observation_model")
    assert q["w2_bound_check"].iloc[1].startswith("QUALIFIED") and q["w2_bound_check"].iloc[3] == ""
    P.check_w2_qualified(q, cons, 2)
    q.loc[2, "w2_bound_check"] = "consistent"
    with pytest.raises(P.BoundContradictionError, match="pulse 2"):
        P.check_w2_qualified(q, cons, 2)


def test_contradiction_at_every_band_vector_says_so():
    cons = P.w2_consistency(_balance(30.0, 25.0))
    assert "at no vector of the C(t) parameter band" in cons["verdict"].iloc[-1]


def test_real_outputs_qualify_every_contradicted_quantity(cfg, run_all_dir):
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv", keep_default_na=False, na_values=[""])
    cons = pd.read_csv(run_all_dir / "tables" / "nh4_truth_w2_consistency.csv")
    pm["w2_bound_check"] = pm["w2_bound_check"].fillna("")
    P.check_w2_qualified(pm, cons, len(cfg.pulse_times_d))
    tw = pd.read_csv(run_all_dir / "tables" / "nh4_truth_windows.csv")
    tw["w2_bound_check"] = tw["w2_bound_check"].fillna("")
    P.check_w2_qualified(tw, cons, len(cfg.pulse_times_d), "governing_pulse")
    # the verdict follows from the numbers: best-fit requirement against the 15 meq/L scan maximum
    assert (cons["consistent_at_best_fit"] == (cons["alk0_required_if_nitrified_meq_L"] <= max(A.ALK0_GRID_MEQ_L))).all()
    assert (cons["alk0_required_if_nitrified_band_min_meq_L"] <= cons["alk0_required_if_nitrified_meq_L"] + 1e-9).all()


# --- S3: η is a fraction of the pulse ---------------------------------------------------------------------------------
def test_eta_is_computed_only_between_zero_and_one(run_all_dir):
    pm = pd.read_csv(run_all_dir / "tables" / "pulses.csv").set_index(["pulse", "metric"])
    for (k, metric), row in pm.iterrows():
        if not metric.startswith("removal_to_dose_ratio"):
            continue
        eta = pm.loc[(k, metric.replace("removal_to_dose_ratio", "eta_conversion_fraction"))]
        if row["status"] != "computed":
            assert eta["status"].startswith("not computed: negative")
        elif row["value"] > 1:
            assert eta["status"].startswith("not computed: removal exceeds the dose") and np.isnan(eta["value"])
        else:
            assert eta["status"] == "computed" and eta["value"] == pytest.approx(row["value"])
    eta = pm[pm.index.get_level_values("metric").str.startswith("eta_conversion_fraction")]
    assert (eta.loc[eta["status"] == "computed", "value"].between(0, 1)).all()


# --- W8 ---------------------------------------------------------------------------------------------------------------
def test_missing_top_up_is_a_status_not_zero(tmp_path):
    c = config_with(tmp_path / "v.toml", vessel_mass_start_g=150.0, vessel_mass_end_g=140.1)
    (s,) = volume.volume_scenarios(c, 50.0)
    assert s.basis == "MEASURED" and "not computed: top_up_mL missing" in s.name and "top-up 0" not in s.name
    assert 1e3 * s.net_loss_L_per_d == pytest.approx(10.0 / 50.0)  # the net loss itself is measured


def test_n_conservation_minimum_loss_hand_value():
    # C rises 160 -> 240 mg N/L over 24.6-32.4 d in 0.1 L with 6.546 mg N dosed: deficit 8 - 6.546 = 1.454 mg N
    e = volume.n_conservation_min_loss_mL_d(160.0, 240.0, 24.6, 32.4, 6.546, 0.1)
    assert e == pytest.approx(1e3 * 1.454 / (240 * 32.4 - 160 * 24.6))
    v = 0.1 - e / 1e3 * np.array([24.6, 32.4])
    assert 160 * v[0] + 6.546 - 240 * v[1] == pytest.approx(0.0, abs=1e-9)
    assert volume.n_conservation_min_loss_mL_d(200.0, 150.0, 10.0, 20.0, 6.5, 0.1) == 0.0
    # with t_k < t_next a feasible loss below V0 / t_next always exists (the constant-volume deficit is at most V0 dC)
    assert volume.n_conservation_min_loss_mL_d(100.0, 900.0, 10.0, 11.0, 0.0, 0.1) < 1e3 * 0.1 / 11.0
    with pytest.raises(ValueError, match="t_next"):
        volume.n_conservation_min_loss_mL_d(100.0, 900.0, 11.0, 10.0, 0.0, 0.1)


def test_ec_conditional_bound():
    assert volume.ec_conditional_max_loss_mL_d(1.0, 1.1, 0.1, 50.0) == pytest.approx(100 * (1 - 1 / 1.1) / 50)
    assert volume.ec_conditional_max_loss_mL_d(1.1, 1.0, 0.1, 50.0) == 0.0
    with pytest.raises(ValueError):
        volume.ec_conditional_max_loss_mL_d(0.0, 1.0, 0.1, 50.0)


def test_real_volume_constraints(run_all_dir):
    vc = pd.read_csv(run_all_dir / "tables" / "volume_constraints.csv")
    assert vc.loc[vc["bound"] == "lower", "basis"].str.startswith("CONDITIONAL").all()
    ec = vc[vc["constraint"].str.startswith("flat conductivity")]
    assert "not used" in ec["basis"].iloc[0]
    assert (vc["bound"] == "scenario").sum() == 4


def test_older_alkalinity_key_is_never_silently_ignored(tmp_path):
    c = config_with(tmp_path / "a.toml", alkalinity_mmol_L=4.0)
    assert A.alk0_source(c).startswith("not computed: initial_alkalinity_meq_L missing; alkalinity_mmol_L = 4")
    assert all(b == "SCAN 1-15 meq/L" for _, b in A.alk0_values(c))
    c2 = config_with(tmp_path / "b.toml", initial_alkalinity_meq_L=4.2)
    assert A.alk0_source(c2).startswith("MEASURED initial_alkalinity_meq_L = 4.2")


# --- R8 -------------------------------------------------------------------------------------------------------------
def test_electrode_fragment_rule():
    j, c = "jump (electrode_offset)", "candidate (unlocated, sustained)"
    assert S.electrode_fragment(j, c, pd.NA, 0.5, 0.33)
    assert not S.electrode_fragment(j, c, pd.NA, 1.2, 0.33)  # longer than 3 delays
    assert not S.electrode_fragment("jump (pulse_response)", c, pd.NA, 0.3, 0.33)  # a pulse response bounds it
    assert not S.electrode_fragment(j, j, 4, 0.3, 0.33)  # governed by a pulse
    assert not S.electrode_fragment(j, j, pd.NA, 0.3, float("nan"))


def test_real_fragments_are_demoted(segmented):
    t = segmented.table.set_index("sub_segment")
    frag = t[t["status"].str.contains("electrode fragment")]
    assert set(frag.index) >= {2, 3, 5, 7}  # the -21.8, +7.97, -8.49, +5.48 mg NH4+/L/d windows of the review
    assert frag["slope_mgNH4_L_d"].isna().all() and frag["slope_fragment_diagnostic_mgNH4_L_d"].notna().all()
    assert frag["governing_pulse"].isna().all()
    assert t.loc[10, "status"] == "computed"  # 6.65-8.53 d ends at pulse 1: not a fragment


# --- W5 -------------------------------------------------------------------------------------------------------------
def test_alkalinity_candidates_are_in_the_model_set(cfg):
    specs = {s.name: s for s in M.model_specs(cfg)}
    assert specs["alkalinity_limited"].alk and "Alk0_meq_L" in specs["alkalinity_limited"].kinetic
    for T in cfg.get("scenarios.temperature_C"):
        s = specs[f"alkalinity_FA_SCENARIO_T{T:g}C"]
        assert s.alk and s.fa and s.T_C == T
    assert "biomass_decay" in specs  # still competing


def test_alkalinity_limitation_stalls_at_the_bound():
    spec = M.ModelSpec("alkalinity_limited", ("v0_mgN_L_d", "Alk0_meq_L"), False, False, alk=True)
    t = np.round(np.arange(1.005, 30.0, 0.01), 10)
    sim = M.simulate(spec, {"v0_mgN_L_d": 5.0, "Alk0_meq_L": 4.0}, 100.0, t, np.full(len(t), 7.5), [10.0], 60.0)
    assert sim.consumed == pytest.approx(4.0 / 2 * 14.007, rel=1e-3)
    assert sim.rate[-1] < 1e-6 and np.all(sim.rate >= 0)
    assert M.alkalinity_state_meq_L(1.0, 100.0) == 0.0


def test_model_fit_lists_the_alkalinity_candidates(run_all_dir):
    fit = pd.read_csv(run_all_dir / "tables" / "model_fit.csv")
    assert {"alkalinity_limited", "alkalinity_FA_SCENARIO_T20C", "alkalinity_FA_SCENARIO_T25C"} <= set(fit["model"])
    assert fit["selected"].sum() == 1
