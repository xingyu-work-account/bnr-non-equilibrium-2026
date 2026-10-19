import math

import pytest

from bnr_neq import constants as C


def test_section6_values():
    assert C.R_J_PER_MOL_K == 8.314462618
    assert C.LAMBDA0_S_CM2_MOL == {"NH4+": 73.5, "K+": 73.5, "Na+": 50.1, "Cl-": 76.3, "NO3-": 71.4, "HCO3-": 44.5}
    assert C.E0_AGCL_AG_V == 0.2223 and C.E0_AG2O_AG_V == 1.173
    assert (C.EMERSON_A, C.EMERSON_B_K, C.ANTHONISEN_EXP_K) == (0.09018, 2729.92, 6344.0)
    assert C.DAVIES_A_25C == 0.509
    assert C.DG0_NITRIFICATION_KJ_MOL_N == -348.8
    assert C.DG0_NITRITATION_KJ_MOL_N + C.DG0_NITRATATION_KJ_MOL_N == pytest.approx(C.DG0_NITRIFICATION_KJ_MOL_N)
    assert (C.ALK_CACO3_PER_N_G_G, C.O2_PER_N_G_G, C.DO_SAT_25C_MG_L) == (7.14, 4.57, 8.26)
    assert (C.PK1_CARBONATE_25C, C.LOG10_KH_CO2_25C) == (6.35, -1.47)


def test_molar_masses_consistent():
    assert C.M_NH4_G_MOL == pytest.approx(C.M_N_G_MOL + 4 * 1.008, abs=1e-9)
    # 14.007/18.039 = 0.776484; the configured factor 0.77649 (implementation_prompt.md Section 3) is this value
    # rounded up in the fifth decimal (relative difference 7.5e-6), documented as a minor discrepancy.
    assert C.M_N_G_MOL / C.M_NH4_G_MOL == pytest.approx(0.77649, abs=1e-5)


def test_nernst_slope_from_constants():
    slope = C.R_J_PER_MOL_K * C.T_REF_K * math.log(10) / C.F_C_PER_MOL
    assert slope * 1e3 == pytest.approx(59.16, abs=0.005)
