import math

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp

from bnr_neq import model as M
from bnr_neq import thermo
from bnr_neq.constants import K_NH_ASM_MG_N_L

DOSE = 65.46428371127853
T = np.round(np.arange(1.005, 7.0, 0.01), 10)


def data(t=T, pulses=(2.53,), offsets=(), log_r=None, ph=7.3):
    return M.ModelData(t=t, log_r=np.zeros(len(t)) if log_r is None else log_r, ph=np.full(len(t), ph),
                       pulses=list(pulses), dose=DOSE, offsets=list(offsets), offset_guess=[0.0] * len(offsets))


NULL = M.ModelSpec("null_no_consumption", (), False, False)
DECAY = M.ModelSpec("biomass_decay", ("v0_mgN_L_d", "kd_per_d"), False, True)
HALD = M.ModelSpec("haldane_FA_SCENARIO_T25C", ("vmax_mgN_L_d", "Ks_mgNH3_L", "Ki_mgNH3_L"), True, False, 25.0)


def test_pulse_adds_exactly_the_dose_and_mass_balance_closes():
    d = data(pulses=(2.53, 4.57))
    sim = M.simulate(DECAY, {"v0_mgN_L_d": 5.0, "kd_per_d": 0.1}, 30.0, d.t, d.ph, d.pulses, d.dose)
    assert [round(t, 6) for t, _, _ in sim.pulse_jumps] == [2.53, 4.57]
    for _, before, after in sim.pulse_jumps:
        assert after - before == pytest.approx(DOSE, abs=1e-12)
    assert sim.c[-1] == pytest.approx(30.0 + 2 * DOSE - sim.consumed, abs=1e-9)
    assert (sim.rate >= 0).all()


def test_null_model_only_adds_pulses():
    d = data()
    sim = M.simulate(NULL, {}, 10.0, d.t, d.ph, d.pulses, d.dose)
    assert sim.consumed == 0.0 and sim.c[-1] == pytest.approx(10.0 + DOSE)


def test_decay_model_matches_ode():
    d = data(pulses=())
    sim = M.simulate(DECAY, {"v0_mgN_L_d": 4.0, "kd_per_d": 0.2}, 25.0, d.t, d.ph, [], d.dose)
    ode = solve_ivp(lambda t, c: [-4.0 * math.exp(-0.2 * (t - 1.0)) * c[0] / (K_NH_ASM_MG_N_L + c[0])],
                    (d.t[0], d.t[-1]), [25.0], t_eval=d.t, rtol=1e-10, atol=1e-10).y[0]
    np.testing.assert_allclose(sim.c, ode, rtol=1e-7)


def test_haldane_model_matches_ode():
    d = data(pulses=())
    kin = {"vmax_mgN_L_d": 6.0, "Ks_mgNH3_L": 0.3, "Ki_mgNH3_L": 0.8}
    sim = M.simulate(HALD, kin, 40.0, d.t, d.ph, [], d.dose)
    f = float(M.fa_factor(np.array([7.3]), 25.0)[0])

    def rhs(t, c):
        s = f * c[0]
        return [-6.0 * s / (0.3 + s + s * s / 0.8)]

    ode = solve_ivp(rhs, (d.t[0], d.t[-1]), [40.0], t_eval=d.t, rtol=1e-10, atol=1e-10).y[0]
    np.testing.assert_allclose(sim.c, ode, rtol=1e-7)


def test_rate_never_negative_when_state_depleted():
    d = data(pulses=())
    sim = M.simulate(DECAY, {"v0_mgN_L_d": 50.0, "kd_per_d": 0.0}, 1.0, d.t, d.ph, [], d.dose)
    assert (sim.rate >= 0).all()


def test_observation_lag_and_offsets():
    t = np.round(np.arange(0.005, 3.0, 0.01), 10)
    c = np.where(t >= 1.0, 20.0, 10.0)
    p = {"a": 0.0, "b0": 1.0, "b1_per_d": 0.0, "tau_d": 0.1, "offset_values": [0.5]}
    lr = M.observe(c, t + (M.T0_D - 0.005), p, [2.5 + M.T0_D - 0.005])
    r = np.exp(lr)
    i = np.searchsorted(t, 1.1)  # one time constant after the step
    assert (r[i] - 10) / 10 == pytest.approx(1 - math.exp(-1), abs=0.06)
    assert r[-1] == pytest.approx(20 * math.exp(0.5), rel=0.02)  # offset applied after its time


def _synthetic(spec, x_true, d, noise, seed=0):
    lay = M.layout(spec, d)
    lr, _ = M.predict(spec, d, lay, x_true)
    d.log_r = lr + noise * np.random.default_rng(seed).normal(size=len(lr))
    return lay


def test_fit_recovers_known_parameters():
    d = data(pulses=(2.53, 4.57))
    lay = M.layout(DECAY, d)
    x_true = np.array([math.log(6.0), math.log(0.3), math.log(30.0), 0.1, 0.9, 0.0, math.log(0.05)])
    _synthetic(DECAY, x_true, d, 0.002)
    free = np.ones(len(x_true), dtype=bool)
    free[[2, 3]] = False  # fix C0 and a (not jointly identifiable with b0)
    x0 = x_true.copy()
    x0[[0, 1, 4, 6]] += [0.3, -0.3, 0.1, 0.5]
    f = M.fit(DECAY, d, lay, (1.0, 7.0), x0, free)
    assert f.converged
    np.testing.assert_allclose(f.x[[0, 1, 4]], x_true[[0, 1, 4]], atol=0.05)


def test_profile_flags_non_identifiable_parameters():
    d = data(pulses=(2.53,))
    lay = M.layout(NULL, d)
    x_true = np.array([math.log(20.0), 0.2, 0.8, 0.0, math.log(0.1)])
    _synthetic(NULL, x_true, d, 0.005)
    best = M.fit(NULL, d, lay, (1.0, 7.0), x_true, np.ones(len(x_true), dtype=bool))
    prof = M.profile_likelihood(NULL, d, lay, best).set_index("parameter")
    assert prof.loc["tau_d", "identifiable"]
    level = prof.loc[["C0_mgN_L", "a", "b0"]]
    assert (~level["identifiable"]).any()  # one jump ratio and one level cannot fix three parameters
    assert (level.loc[~level["identifiable"], "estimate"] == "not identifiable").all()
    assert level.loc[~level["identifiable"], "ci95_low"].isna().all()


def test_effective_n():
    n_eff, _rho = M.effective_n(np.random.default_rng(1).normal(size=1000))
    assert n_eff == pytest.approx(1000, rel=0.1)
    ar = np.cumsum(np.random.default_rng(2).normal(size=1000)) * 0.01
    assert M.effective_n(ar)[0] < 20


def test_split_evaluation_on_synthetic_record(cfg):
    d = data(pulses=(2.53, 4.57), offsets=(3.2, 5.6))
    lay = M.layout(DECAY, d)
    x_true = np.array([math.log(6.0), math.log(0.2), math.log(30.0), 0.0, 0.9, 0.0, math.log(0.05), -0.2, 0.15])
    _synthetic(DECAY, x_true, d, 0.003)
    df, fits = M.split_evaluation(cfg, d, M.SplitWindows((1.0, 4.0), (4.0, 5.0), (5.0, 7.0)))
    assert df["selected"].sum() == 1 and not df.loc[df["selected"], "model"].iloc[0].startswith("null")
    assert "5.600" in df["test_nuisance"].iloc[0] and "3.200" not in df["test_nuisance"].iloc[0]
    sel = df[df["selected"]].iloc[0]
    assert sel["sel_skill_vs_null"] > 0.5 and sel["forecast_window_d"] == "5.00 to 5.600"
    _spec, lay, cal = fits["biomass_decay"]
    off_test = lay.names.index("offset_log_5.600d")
    assert cal.x[off_test] == 0.0 and not cal.free[off_test]  # unseen offset never calibrated


def test_consistency_checks_pass(cfg):
    d = data(pulses=(2.53, 4.57))
    lay = M.layout(DECAY, d)
    x = np.array([math.log(6.0), math.log(0.2), math.log(30.0), 0.0, 0.9, 0.0, math.log(0.05)])
    best = M.Fit(DECAY, x, np.ones(len(x), dtype=bool), (1.0, 7.0), 1.0, len(d.t), True)
    checks = M.consistency_checks(cfg, DECAY, d, lay, best, 316.9)
    computed = checks[checks["pass"] != "not computed"]
    assert len(computed) == len(checks) - 1 and computed["pass"].astype(bool).all()
    assert checks.loc[checks["check"] == "σ in situ", "pass"].iloc[0] == "not computed"
    assert (checks["check"] == "pulse adds exactly the dose").sum() == 2
    assert checks.loc[checks["check"] == "σ in situ", "value"].iloc[0] == "not computed: DO_mg_L missing"


def test_model_specs_label_scenarios(cfg):
    names = [s.name for s in M.model_specs(cfg)]
    assert "null_no_consumption" in names and "biomass_decay" in names
    assert all("SCENARIO" in n for n in names if n.startswith(("haldane", "combined")))


def test_load_inputs_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        M.load_inputs(tmp_path)


def test_build_data_uses_detected_offsets(cfg, ingested, validated):
    d = M.build_data(cfg, ingested.grid, validated.jumps)
    off = validated.jumps[validated.jumps["jump"] & (validated.jumps["class"] == "electrode_offset")]
    assert d.offsets == sorted(off["t_jump_d"].tolist())
    assert d.t[0] > M.T0_D and d.t[-1] < M.END_D and d.dose == pytest.approx(DOSE)
    assert isinstance(d.log_r, np.ndarray) and np.isfinite(d.log_r).all()


def test_response_fraction_table(cfg):
    d = data(t=np.round(np.arange(1.005, 12.0, 0.01), 10), pulses=(2.53, 8.53))
    lay = M.layout(NULL, d)
    x = np.array([math.log(20.0), 0.0, 1.0, 0.0, math.log(0.01)])
    best = M.Fit(NULL, x, np.ones(len(x), dtype=bool), (1.0, 12.0), 1.0, len(d.t), True)
    pt = pd.DataFrame({"pulse": [1, 2], "nh4_response_fraction": [0.3, 0.2]})
    r = M.response_fractions(NULL, d, lay, best, pt, 0.77649)
    assert r["model_response_fraction"].iloc[0] == pytest.approx(1.0, abs=0.02)  # b = 1, no consumption: full response


@pytest.mark.parametrize("mu, sigma_low, sigma_high", [(0.0, 1.0, 1.0), (2.5, 0.3, 0.7), (-1.0, 2.0, 0.5)])
def test_profile_interval_matches_analytic_quadratic(mu, sigma_low, sigma_high):
    """Δ(θ) = ((θ − μ)/σ)², with σ different on each side: the 95 % ends are μ ∓ √3.841 σ exactly
    (linear interpolation between the 3-point grid would miss them by about 0.1 σ)."""
    def delta(v):
        s = sigma_low if v < mu else sigma_high
        return ((v - mu) / s) ** 2, True

    half = 4 * max(sigma_low, sigma_high)
    out = M.profile_interval(delta, mu, half, -1e6, 1e6)
    z = math.sqrt(M.CHI2_95_1DOF)
    assert out["low"]["closed"] and out["high"]["closed"]
    assert out["low"]["crossing"] == pytest.approx(mu - z * sigma_low, abs=2e-3 * half)
    assert out["high"]["crossing"] == pytest.approx(mu + z * sigma_high, abs=2e-3 * half)


def test_profile_interval_open_side_and_bounds():
    out = M.profile_interval(lambda v: (min(v, 0.0) ** 2, True), 0.0, 4.0, -10.0, 10.0)  # flat above the optimum
    assert out["low"]["closed"] and not out["high"]["closed"] and np.isnan(out["high"]["crossing"])
    unconverged = M.profile_interval(lambda v: (v * v, False), 0.0, 4.0, -10.0, 10.0)
    assert not unconverged["low"]["closed"] and not unconverged["high"]["closed"]


def test_consistency_check_in_situ_sigma_uses_model_state(measured_cfg):
    c = measured_cfg(DO_mg_L=0.1, NO3_N_mg_L=100.0)
    d = data(pulses=(2.53, 4.57), ph=7.6)
    lay = M.layout(DECAY, d)
    x = np.array([math.log(6.0), math.log(0.2), math.log(30.0), 0.0, 0.9, 0.0, math.log(0.05)])
    best = M.Fit(DECAY, x, np.ones(len(x), dtype=bool), (1.0, 7.0), 1.0, len(d.t), True)
    checks = M.consistency_checks(c, DECAY, d, lay, best, 316.9).set_index("check")
    row = checks.loc["σ in situ (SCENARIO 25 °C)"]
    assert isinstance(row["pass"], (bool, np.bool_)) and bool(row["pass"])  # a bool, not "not computed"
    _, sim = M.predict(DECAY, d, lay, x)
    a = thermo.nitrification_affinity_kJ_mol(sim.c, 100.0, 0.1, d.ph)
    sig = sim.rate / 14.007 / 1e3 / 86400 * a * 1e3 / 298.15
    assert row["value"] == pytest.approx(float(sig.min()))
    a_placeholder = float(thermo.nitrification_affinity_kJ_mol(1.0, 100.0, 0.1, 7.0))
    assert row["value"] != pytest.approx(float((sim.rate / 14.007 / 1e3 / 86400 * a_placeholder * 1e3 / 298.15).min()))


def test_split_windows_from_config(cfg):
    w = M.split_windows(cfg)
    p = cfg.pulse_times_d
    assert w.cal == (M.T0_D, p[2]) and w.sel == (p[2], p[3]) and w.test == (p[3], M.END_D)


def test_model_compute_on_real_data(cfg, run_all_dir):
    r = M.compute(cfg, M.load_inputs(run_all_dir))
    assert r.selected in set(r.comparison["model"]) and r.comparison["selected"].sum() == 1
    done = r.checks[r.checks["pass"] != "not computed"]
    assert done["pass"].astype(bool).all()
    assert (r.checks.loc[r.checks["check"] == "σ in situ", "pass"] == "not computed").all()
    ident = r.identifiability
    ok = ident[ident["identifiable"]]
    assert len(ok) > 0 and (ok["ci95_low"].astype(float) < ok["ci95_high"].astype(float)).all()
    on_disk = pd.read_csv(run_all_dir / "tables" / "model_identifiability.csv")
    np.testing.assert_allclose(on_disk["ci95_low"], ident["ci95_low"].astype(float), equal_nan=True)
