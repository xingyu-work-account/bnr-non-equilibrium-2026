import numpy as np
import pytest

from bnr_neq import thermo


def test_pulse_dose(cfg):
    d = thermo.pulse_dose_from_config(cfg)
    assert d.c_mmol_L == pytest.approx(4.6737, abs=5e-5)
    assert d.nh4_mg_L == pytest.approx(84.31, abs=5e-3)
    assert d.n_mg_L == pytest.approx(65.46, abs=5e-3)
    assert d.cl_mg_L == pytest.approx(165.7, abs=5e-2)


def test_kohlrausch_predictions():
    c = 25.0 / (53.491 * 0.100)
    assert thermo.ec_step_immediate_mS_cm(c) == pytest.approx(0.700, abs=5e-4)
    assert thermo.ec_step_nitrified_mS_cm(c) == pytest.approx(0.274, abs=5e-4)


def test_kohlrausch_unknown_ion_raises():
    with pytest.raises(KeyError):
        thermo.kohlrausch_delta_ec_mS_cm({"SO4--": 1.0})


def test_nernst_slope_25C():
    assert 1e3 * thermo.nernst_slope_V(298.15) == pytest.approx(59.16, abs=0.005)


@pytest.mark.parametrize("c0, expected_mV", [(1.0, 44.6), (13.5, 7.6)])
def test_chloride_orp_rise(c0, expected_mV):
    assert thermo.chloride_orp_rise_mV(c0, 4.6737, 298.15) == pytest.approx(expected_mV, abs=0.05)


def test_chloride_orp_rise_cumulative_totals():
    # four pulses together: 22 to 77 mV (plan, finding S4)
    assert thermo.chloride_orp_rise_mV(13.5, 4 * 4.6737, 298.15) == pytest.approx(22.3, abs=0.1)
    assert thermo.chloride_orp_rise_mV(1.0, 4 * 4.6737, 298.15) == pytest.approx(76.6, abs=0.1)


@pytest.mark.parametrize(
    "call",
    [
        lambda: thermo.pulse_dose(-25.0, 53.491, 0.1),
        lambda: thermo.pulse_dose(25.0, 53.491, 0.0),
        lambda: thermo.chloride_orp_rise_mV(0.0, 4.67, 298.15),
        lambda: thermo.chloride_orp_rise_mV(1.0, -4.67, 298.15),
        lambda: thermo.nernst_slope_V(-1.0),
        lambda: thermo.ec_step_immediate_mS_cm(-1.0),
    ],
)
def test_impossible_inputs_raise(call):
    with pytest.raises(ValueError):
        call()


def test_emerson_pka_25C():
    assert thermo.emerson_pka(298.15) == pytest.approx(9.246, abs=0.002)


def test_free_ammonia_matches_pka_speciation():
    # Anthonisen's exp(6344/T) is Ka^-1 in its own fit; check it against Emerson within 2 %
    T, pH, tan = 298.15, 7.5, 30.0
    frac = 1.0 / (1.0 + 10 ** (thermo.emerson_pka(T) - pH))
    assert thermo.free_ammonia_mgNH3_L(tan, pH, T) == pytest.approx(17 / 14 * tan * frac, rel=0.02)


def test_free_ammonia_vectorized_and_monotonic():
    fa = thermo.free_ammonia_mgNH3_L(np.array([10.0, 10.0]), np.array([7.0, 8.0]), 298.15)
    assert fa[1] > fa[0] > 0


@pytest.mark.parametrize(
    "call",
    [
        lambda: thermo.free_ammonia_mgNH3_L(-1.0, 7.0, 298.15),
        lambda: thermo.free_ammonia_mgNH3_L(10.0, 14.5, 298.15),
        lambda: thermo.free_ammonia_mgNH3_L(10.0, -0.1, 298.15),
        lambda: thermo.free_ammonia_mgNH3_L(10.0, float("nan"), 298.15),
        lambda: thermo.emerson_pka(0.0),
    ],
)
def test_impossible_speciation_inputs_raise(call):
    with pytest.raises(ValueError):
        call()
