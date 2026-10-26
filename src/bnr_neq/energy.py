"""Oxygen and energy metric family (W4).

The laboratory vessel was stirred, not aerated, and its stirrer power is not measured, so the energy terms that
transfer to a BNR plant are built on oxygen:

- oxygen required by the nitrogen actually removed, 4.57 g O2/g N (nitrification; 3.43 g O2/g N for nitritation as a
  stoichiometry case);
- oxygen supply capacity kLa (C* − C) V Δt under a kLa SCENARIO (C = measured DO, else 0: the capacity);
- oxygen utilization = required / supply capacity;
- aeration-equivalent energy = O2 / SAE under a field aeration-efficiency SCENARIO with a cited source, converted to
  kWh per kg N removed; the demand-based value is the floor, the supply-based value includes idle aeration;
- Gibbs energy released by the removal (|ΔG°′| at pH 7, SCENARIO 25 °C) over the aeration-equivalent energy.

Units: mg, L, d, kWh; 1 kg O2 = 1e6 mg; 1 kWh = 3.6e6 J.
"""

from __future__ import annotations

import numpy as np

from .constants import (
    DG0_NITRIFICATION_KJ_MOL_N,
    J_PER_KWH,
    M_N_G_MOL,
    O2_PER_N_G_G,
    O2_PER_N_NITRITATION_G_G,
)

MG_PER_KG = 1e6


def _positive(name: str, v) -> None:
    a = np.asarray(v, dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a <= 0):
        raise ValueError(f"{name} must be positive and finite, got {v!r}")


def _non_negative(name: str, v) -> None:
    a = np.asarray(v, dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError(f"{name} must be non-negative and finite, got {v!r}")


def oxygen_required_mgO2(n_removed_mgN, nitritation: bool = False):
    """O2 needed to oxidize the removed N: 4.57 g O2/g N (nitrification) or 3.43 g O2/g N (nitritation)."""
    _non_negative("n_removed_mgN", n_removed_mgN)
    return (O2_PER_N_NITRITATION_G_G if nitritation else O2_PER_N_G_G) * np.asarray(n_removed_mgN, dtype=float)


def o2_supply_capacity_mgO2(kla_per_d, c_star_mg_L: float, do_mg_L: float, volume_L: float, duration_d: float):
    """kLa (C* − C) V Δt, mg O2 transferable over Δt."""
    _positive("kla_per_d", kla_per_d)
    _positive("c_star_mg_L", c_star_mg_L)
    _non_negative("do_mg_L", do_mg_L)
    _positive("volume_L", volume_L)
    _positive("duration_d", duration_d)
    if do_mg_L > c_star_mg_L:
        raise ValueError(f"DO {do_mg_L} mg/L above saturation {c_star_mg_L} mg/L: no transfer into the liquid")
    return np.asarray(kla_per_d, dtype=float) * (c_star_mg_L - do_mg_L) * volume_L * duration_d


def utilization(required_mgO2, supply_mgO2):
    """Oxygen utilization efficiency = required / supplied (≤ 1 is feasible; > 1 means the supply cannot carry it)."""
    _non_negative("required_mgO2", required_mgO2)
    _positive("supply_mgO2", supply_mgO2)
    return np.asarray(required_mgO2, dtype=float) / np.asarray(supply_mgO2, dtype=float)


def aeration_energy_kWh(o2_mg, sae_kgO2_per_kWh):
    """Aeration-equivalent energy of transferring ``o2_mg`` at a field aeration efficiency SAE (kg O2/kWh)."""
    _positive("SAE", sae_kgO2_per_kWh)
    _non_negative("o2_mg", o2_mg)
    return np.asarray(o2_mg, dtype=float) / MG_PER_KG / np.asarray(sae_kgO2_per_kWh, dtype=float)


def specific_energy_kWh_per_kgN(energy_kWh, n_removed_mgN):
    """kWh per kg N removed."""
    _non_negative("energy_kWh", energy_kWh)
    _positive("n_removed_mgN", n_removed_mgN)
    return np.asarray(energy_kWh, dtype=float) / (np.asarray(n_removed_mgN, dtype=float) / MG_PER_KG)


def gibbs_released_J(n_removed_mgN):
    """|ΔG°′| (nitrification, pH 7, 25 °C) x moles of N removed, J."""
    _non_negative("n_removed_mgN", n_removed_mgN)
    return abs(DG0_NITRIFICATION_KJ_MOL_N) * 1e3 * np.asarray(n_removed_mgN, dtype=float) / 1e3 / M_N_G_MOL


def kwh_to_J(kwh):
    return np.asarray(kwh, dtype=float) * J_PER_KWH
