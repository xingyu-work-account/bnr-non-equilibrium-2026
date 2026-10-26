"""Liquid volume V(t) of the open, stirred vessel (W8).

Every earlier milestone assumed a closed 100 mL frame. An open beaker stirred for about 50 days can lose
millilitres per day by evaporation, which concentrates every dissolved species. This module turns the
optional vessel measurements into a volume trajectory, or, while they are missing, into labelled
evaporation SCENARIOS, so that each concentration balance can be evaluated per scenario:

- ``vessel_mass_start_g`` and ``vessel_mass_end_g``: MEASURED net loss (evaporation minus top-up), spread linearly
  in time; the split into evaporation and top-up needs ``top_up_mL`` and is otherwise ``not computed`` (never 0);
- ``covered = true`` without masses: evaporation ASSUMED negligible (cover measured, loss not);
- otherwise: the SCENARIO rates of ``[scenarios] evaporation_mL_per_d``.

Water density is taken as 1 g/mL. Amounts (mg) are conserved by evaporation; concentrations are
amount / V(t).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import thermo
from .config import Config, not_computed
from .constants import M_N_G_MOL

CONSTANT = "constant volume (closed-frame assumption)"
WATER_G_PER_ML = 1.0


@dataclass(frozen=True)
class VolumeScenario:
    name: str
    basis: str  # MEASURED, ASSUMPTION or SCENARIO
    v0_L: float
    net_loss_L_per_d: float  # evaporation minus top-up, spread uniformly in time
    t0_d: float = 0.0

    def volume_L(self, t_d) -> np.ndarray:
        v = self.v0_L - self.net_loss_L_per_d * (np.asarray(t_d, dtype=float) - self.t0_d)
        if np.any(v <= 0):
            raise ValueError(f"volume scenario {self.name!r} empties the vessel before t = {np.max(t_d):g} d")
        return v

    def factor(self, t_d) -> np.ndarray:
        """Concentration factor V0 / V(t) of a conserved amount."""
        return self.v0_L / self.volume_L(t_d)

    @property
    def is_constant(self) -> bool:
        return self.net_loss_L_per_d == 0.0


def slug(name: str) -> str:
    """Scenario name as a column / metric suffix."""
    return name.replace("/", "_per_").replace(" ", "_")


def scenario_label(e_mL_d: float) -> str:
    return CONSTANT if e_mL_d == 0 else f"SCENARIO evaporation {e_mL_d:g} mL/d"


def evaporation_upper_bound_mL_d(v0_L: float, t_last_d: float) -> float:
    """Net loss rate that would empty the vessel by the last sensor reading: liquid was still present then
    (the probes recorded), so the net loss rate is below this value."""
    return 1e3 * v0_L / t_last_d


def volume_scenarios(cfg: Config, t_end_d: float) -> list[VolumeScenario]:
    """Volume trajectories over [0, t_end_d]: one MEASURED trajectory when the vessel masses exist, else the
    constant frame plus the configured evaporation SCENARIOS (a covered vessel keeps only the constant frame)."""
    v0 = float(cfg.get("vessel.volume_L"))
    if cfg.has_measurement("vessel_mass_start_g"):
        m0, m1 = float(cfg.measurement("vessel_mass_start_g")), float(cfg.measurement("vessel_mass_end_g"))
        salt_g = len(cfg.pulse_times_d) * float(cfg.get("pulses.mass_mg")) / 1e3  # dosed solid is inside m1
        net_L = (m0 - m1 + salt_g) / WATER_G_PER_ML / 1e3  # top-up water is inside m1 already
        if cfg.has_measurement("top_up_mL"):
            top_up_L = float(cfg.measurement("top_up_mL")) / 1e3
            split = f"evaporation {1e3 * (net_L + top_up_L):g} mL, top-up {1e3 * top_up_L:g} mL"
        else:  # the net loss is measured either way; only its split needs the top-up record
            split = f"evaporation/top-up split {not_computed('top_up_mL')}"
        name = f"MEASURED net loss {1e3 * net_L:g} mL ({split}) over {t_end_d:g} d, spread linearly (ASSUMPTION)"
        return [VolumeScenario(name, "MEASURED", v0, net_L / t_end_d)]
    if cfg.has_measurement("covered") and cfg.measurement("covered"):
        return [VolumeScenario(CONSTANT + "; vessel covered (measured), loss not measured", "ASSUMPTION", v0, 0.0)]
    out = []
    for e in cfg.get("scenarios.evaporation_mL_per_d"):
        out.append(VolumeScenario(scenario_label(float(e)), "ASSUMPTION" if e == 0 else "SCENARIO", v0, float(e) / 1e3))
    return out


def pulse_n_mg(cfg: Config) -> float:
    """Nitrogen mass added by one pulse, mg N (independent of the volume)."""
    return float(cfg.get("pulses.mass_mg")) / float(cfg.get("pulses.molar_mass_g_mol")) * M_N_G_MOL


def dose_mgN_L(cfg: Config, scen: VolumeScenario, t_d: float) -> float:
    """Concentration step of one pulse added at t_d into V(t_d), mg N/L."""
    return pulse_n_mg(cfg) / float(scen.volume_L(t_d))


def scenario_table(cfg: Config, t_end_d: float, t_last_d: float) -> pd.DataFrame:
    """One row per volume scenario: V at each pulse and at the end, the concentration factor and the pulse dose
    in mg N/L at the moment of each pulse. ``t_last_d`` is the last sensor reading (liquid still present)."""
    bound = evaporation_upper_bound_mL_d(float(cfg.get("vessel.volume_L")), t_last_d)
    rows = []
    for s in volume_scenarios(cfg, t_end_d):
        row = {"volume_scenario": s.name, "basis": s.basis, "net_loss_mL_per_d": 1e3 * s.net_loss_L_per_d,
               "net_loss_upper_bound_mL_per_d": bound,
               "within_bound": bool(1e3 * s.net_loss_L_per_d < bound)}
        for k, tp in enumerate(cfg.pulse_times_d, start=1):
            row[f"V_at_pulse{k}_mL"] = 1e3 * float(s.volume_L(tp))
            row[f"dose_pulse{k}_mgN_L"] = dose_mgN_L(cfg, s, tp)
        row["V_end_mL"] = 1e3 * float(s.volume_L(t_end_d))
        row["concentration_factor_end"] = float(s.factor(t_end_d))
        row["t_end_d"] = t_end_d
        rows.append(row)
    df = pd.DataFrame(rows)
    df["note"] = (f"upper bound = V0 / {t_last_d:g} d (liquid present at the last sensor reading; the stated facts say "
                  "nothing else was added); constant-volume dose = "
                  f"{thermo.pulse_dose_from_config(cfg).n_mg_L:.4f} mg N/L")
    return df


# --- Data constraints on the net loss rate (W8) ---------------------------------------------------------------
def n_conservation_min_loss_mL_d(c_before: float, c_next: float, t_k: float, t_next: float, n_mg: float,
                                 v0_L: float) -> float:
    """Smallest uniform net loss rate E (mL/d) for which the amount balance of one pulse interval removes N ≥ 0:
    C_before V(t_k) + m_N − C_next V(t_next) ≥ 0 with V(t) = V0 − E t (no N source). Returns 0 when the constant
    frame already conserves N. For t_k < t_next and m_N ≥ 0 the root lies below V0/t_next, so a loss rate always
    exists; whether it is plausible is judged against the other constraints."""
    if not t_next > t_k:
        raise ValueError(f"need t_next > t_k, got {t_k!r} and {t_next!r}")
    deficit = v0_L * (c_next - c_before) - n_mg  # mg N that the constant frame would have to create
    if deficit <= 0:
        return 0.0
    slope = c_next * t_next - c_before * t_k  # d(removal)/dE, mg N per (L/d)
    return 1e3 * deficit / slope


def ec_conditional_max_loss_mL_d(ec_start: float, ec_end: float, v0_L: float, t_d: float) -> float:
    """Net loss rate at which evaporation alone would raise a conserved conductivity from ``ec_start`` to ``ec_end``
    (EC ∝ V0/V, all ions conserved, the NH4Cl additions ignored, which only tightens it). Valid only if the probe
    measures solution conductivity, which M2 rejected; reported as CONDITIONAL, never used as a constraint."""
    if not (ec_start > 0 and ec_end > 0 and t_d > 0):
        raise ValueError("need positive conductivities and duration")
    return max(0.0, 1e3 * v0_L * (1.0 - ec_start / ec_end) / t_d)
