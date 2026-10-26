"""Proton, alkalinity and carbonate balance: an upper bound on cumulative nitrification (W2).

Nitrification, NH4+ + 2 O2 -> NO3- + 2 H+ + H2O, releases 2 eq of acidity per mol N (7.14 g CaCO3/g N).
Nothing but NH4Cl was added to the batch, so with no alkalinity-producing process (denitrification,
ammonification of organic N, mineral dissolution: assumptions, not measured), the alkalinity at time t is

    Alk(t) V(t) = Alk0 V0 − 2 n_nit(t),

per litre of initial volume. The pH electrode (the only validated sensor) bounds Alk(t) from below:
dissolved CO2 is produced in the liquid (the released H+ converts HCO3- to CO2) and can only leave
through the surface, so [CO2] ≥ K_H p_CO2,air, and at the observed pH

    Alk(t) ≥ Alk_min(pH(t)) = [HCO3-] + 2 [CO3 2-] + [OH-] − [H+] at [CO2] = K_H p_CO2 (+ excess).

Hence N_nit(t) ≤ M_N (Alk0 − Alk_min(t) V(t)/V0) / 2, and because n_nit never decreases, the bound at t is
the minimum of that expression over every later time (reverse cumulative minimum). Air-equilibrium CO2
gives the conservative (loosest) bound; a finite CO2 transfer coefficient kLa_CO2 (SCENARIO range) keeps
CO2 above saturation and tightens it. Any proton-releasing N sink other than nitrification (NH3
volatilization NH4+ -> NH3(g) + H+, NH4+ uptake) releases 1 eq/mol N; that looser bound is reported too.

Carbonate constants: Plummer and Busenberg (1982), Geochim. Cosmochim. Acta 46, 1011-1040; Kw: Harned and
Owen form pKw = 4470.99/T − 6.0875 + 0.01706 T. Activities are taken equal to concentrations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from . import thermo, volume
from .config import Config
from .constants import KELVIN_OFFSET, M_N_G_MOL

ALK0_GRID_MEQ_L = tuple(float(x) for x in range(1, 16))  # W2: 1 to 15 meq/L
EQ_PER_MOL_N_NITRIFICATION = 2.0
EQ_PER_MOL_N_ANY_SINK = 1.0
PH_SMOOTH_BINS = 11  # centred rolling median of the 0.01 d pH grid (0.11 d), so a single bin cannot set the bound
T0_D = 1.0  # start of the nitrification record (start-up excluded, as in M3 and M6)
AIR_EQ = "air-equilibrium CO2 (conservative)"


# --- Carbonate equilibria ----------------------------------------------------------------
def carbonate_constants(T_K: float) -> dict[str, float]:
    """pK1, pK2, log10 K_H (mol/(L atm)) and pKw at temperature T (K)."""
    if not T_K > 0:
        raise ValueError(f"T_K must be positive, got {T_K!r}")
    lg = math.log10(T_K)
    log_kh = 108.3865 + 0.01985076 * T_K - 6919.53 / T_K - 40.45154 * lg + 669365.0 / T_K**2
    log_k1 = -356.3094 - 0.06091964 * T_K + 21834.37 / T_K + 126.8339 * lg - 1684915.0 / T_K**2
    log_k2 = -107.8871 - 0.03252849 * T_K + 5151.79 / T_K + 38.92561 * lg - 563713.9 / T_K**2
    pkw = 4470.99 / T_K - 6.0875 + 0.01706 * T_K
    return {"pK1": -log_k1, "pK2": -log_k2, "log10_KH": log_kh, "pKw": pkw}


def co2_saturation_mM(T_K: float, pco2_atm: float) -> float:
    thermo._require_positive("pco2_atm", pco2_atm)
    return 1e3 * 10.0 ** carbonate_constants(T_K)["log10_KH"] * pco2_atm


def hco3_air_eq_mmol_L(pH, T_K: float, pco2_atm: float):
    """[HCO3-] = K1 K_H p_CO2 / [H+], mmol/L."""
    thermo._require_ph(pH)
    k = carbonate_constants(T_K)
    return co2_saturation_mM(T_K, pco2_atm) * 10.0 ** (np.asarray(pH, dtype=float) - k["pK1"])


def carbonate_alkalinity_meq_L(pH, co2_aq_mM, T_K: float):
    """Alk = [HCO3-] + 2 [CO3 2-] + [OH-] − [H+] (meq/L) at the given pH and dissolved CO2 (mmol/L)."""
    thermo._require_ph(pH)
    thermo._require_non_negative("co2_aq_mM", co2_aq_mM)
    k = carbonate_constants(T_K)
    ph = np.asarray(pH, dtype=float)
    hco3 = np.asarray(co2_aq_mM, dtype=float) * 10.0 ** (ph - k["pK1"])
    co3 = hco3 * 10.0 ** (ph - k["pK2"])
    return hco3 + 2.0 * co3 + 1e3 * 10.0 ** (ph - k["pKw"]) - 1e3 * 10.0 ** (-ph)


def alk_min_meq_L(pH, T_K: float, pco2_atm: float, co2_excess_mM: float = 0.0):
    """Lowest alkalinity compatible with the observed pH, for dissolved CO2 at air saturation plus an excess."""
    return carbonate_alkalinity_meq_L(pH, co2_saturation_mM(T_K, pco2_atm) + co2_excess_mM, T_K)


# --- Bound -------------------------------------------------------------------------------
def nitrification_bound_mgN_L(alk0_meq_L, alk_min_meq_L=0.0, volume_factor=1.0,
                              eq_per_mol_N: float = EQ_PER_MOL_N_NITRIFICATION):
    """N_nit ≤ M_N (Alk0 − Alk_min / f) / eq, mg N per litre of initial volume; f = V0/V(t); never negative.

    ``alk_min_meq_L = 0`` ignores the pH floor (then 2, 5, 10 meq/L give 14, 35, 70 mg N/L).
    """
    a0 = np.asarray(alk0_meq_L, dtype=float)
    if not np.all(np.isfinite(a0)) or np.any(a0 < 0):
        raise ValueError(f"initial alkalinity must be non-negative and finite, got {alk0_meq_L!r}")
    thermo._require_non_negative("alk_min_meq_L", alk_min_meq_L)
    thermo._require_positive("volume_factor", volume_factor)
    thermo._require_positive("eq_per_mol_N", eq_per_mol_N)
    return M_N_G_MOL * np.clip(a0 - np.asarray(alk_min_meq_L, dtype=float) / volume_factor, 0.0, None) / eq_per_mol_N


def bound_curve(t_d: np.ndarray, ph: np.ndarray, alk0_meq_L: float, T_K: float, pco2_atm: float,
                scen: volume.VolumeScenario, eq_per_mol_N: float = EQ_PER_MOL_N_NITRIFICATION) -> np.ndarray:
    """Upper bound on cumulative nitrification at each time (air-equilibrium CO2): the reverse cumulative
    minimum of the pointwise bound, because n_nit(t) ≤ n_nit(s) for every s ≥ t."""
    point = nitrification_bound_mgN_L(alk0_meq_L, alk_min_meq_L(ph, T_K, pco2_atm), scen.factor(t_d), eq_per_mol_N)
    return np.minimum.accumulate(point[::-1])[::-1]


def co2_excess_mM(n_nit_mgN_L: float, duration_d: float, kla_co2_per_d: float) -> float:
    """Pseudo-steady CO2 excess: production 2 n_nit / (M_N Δt) (one CO2 per released H+, below pH 8.3)
    divided by kLa_CO2. Uses the mean rate over Δt, the smallest excess for a given n_nit (conservative)."""
    thermo._require_positive("kla_co2_per_d", kla_co2_per_d)
    thermo._require_positive("duration_d", duration_d)
    return EQ_PER_MOL_N_NITRIFICATION * n_nit_mgN_L / M_N_G_MOL / duration_d / kla_co2_per_d


def bound_with_co2(alk0_meq_L: float, ph: float, T_K: float, pco2_atm: float, factor: float,
                   duration_d: float, kla_co2_per_d: float) -> float:
    """Self-consistent bound N = M_N (Alk0 − Alk_min(pH, CO2_sat + excess(N))/f)/2, solved for N."""
    hi = float(nitrification_bound_mgN_L(alk0_meq_L, alk_min_meq_L(ph, T_K, pco2_atm), factor))
    if hi == 0.0:
        return 0.0

    def g(n):
        amin = alk_min_meq_L(ph, T_K, pco2_atm, co2_excess_mM(n, duration_d, kla_co2_per_d))
        return n - float(nitrification_bound_mgN_L(alk0_meq_L, amin, factor))

    return float(brentq(g, 0.0, hi, xtol=1e-9))


def smoothed_ph(grid: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    g = grid[grid["t_start_d"] >= T0_D - 1e-9]
    ph = g["pH"].rolling(PH_SMOOTH_BINS, center=True, min_periods=1).median().to_numpy()
    return g["t_mid_d"].to_numpy(), ph


def checkpoints(cfg: Config, t: np.ndarray) -> list[tuple[str, int]]:
    """Grid index just before each pulse and at the end of the record."""
    out = [(f"pre_pulse{k}", int(np.searchsorted(t, tp, side="left")) - 1) for k, tp in enumerate(cfg.pulse_times_d, 1)]
    return [*out, ("end", len(t) - 1)]


def scenario_temperatures(cfg: Config) -> list[float]:
    return [float(T) for T in cfg.get("scenarios.temperature_C")]


def alk0_source(cfg: Config) -> str:
    """Which configuration value, if any, the W2 bound takes as Alk0, stated so that no measurement is silently ignored.

    ``alkalinity_mmol_L`` (the older key, used for Ω_CO2 and the pH prediction) carries no sampling time. Alkalinity
    only falls during nitrification, so a later sample is a lower bound on Alk0 and would make the upper bound on
    N_nit too tight; it is therefore never taken as Alk0, and the status says so."""
    if cfg.has_measurement("initial_alkalinity_meq_L"):
        return f"MEASURED initial_alkalinity_meq_L = {float(cfg.measurement('initial_alkalinity_meq_L')):g} meq/L (plus the 1-15 meq/L scan)"
    if cfg.has_measurement("alkalinity_mmol_L"):
        return (f"not computed: initial_alkalinity_meq_L missing; alkalinity_mmol_L = "
                f"{float(cfg.measurement('alkalinity_mmol_L')):g} is measured but has no sampling time (a sample after "
                "day 0 is only a lower bound on Alk0), so it is not used as Alk0; only the 1-15 meq/L scan is reported")
    return "not computed: initial_alkalinity_meq_L missing; only the 1-15 meq/L scan is reported"


def alk0_values(cfg: Config) -> list[tuple[float, str]]:
    vals = [(a, "SCAN 1-15 meq/L") for a in ALK0_GRID_MEQ_L]
    if cfg.has_measurement("initial_alkalinity_meq_L"):
        vals.append((float(cfg.measurement("initial_alkalinity_meq_L")), "MEASURED"))
    return vals


def bound_table(cfg: Config, grid: pd.DataFrame) -> pd.DataFrame:
    """N_nit upper bound against Alk0 (1-15 meq/L and any measured value), per temperature SCENARIO, CO2 scenario
    and volume scenario, with the implied minimum residual TAN and free ammonia at the end of the record.

    Residual TAN takes the initial TAN ≥ 0 (conservative: the electrode cannot fix it) and no N source:
    TAN_end ≥ (dosed − N_removed,max) V0/V_end.
    """
    t, ph = smoothed_ph(grid)
    pco2 = float(cfg.get("scenarios.pCO2_air_atm"))
    marks = checkpoints(cfg, t)
    i_end = marks[-1][1]
    dosed = len(cfg.pulse_times_d) * volume.pulse_n_mg(cfg) / float(cfg.get("vessel.volume_L"))  # per V0 litre
    klas = [float(k) for k in cfg.get("scenarios.kLa_CO2_per_d")]
    rows = []
    for scen in volume.volume_scenarios(cfg, float(t[-1])):
        f = scen.factor(t)
        for T in scenario_temperatures(cfg):
            T_K = T + KELVIN_OFFSET
            amin = alk_min_meq_L(ph, T_K, pco2)
            alk0_floor = float(np.max(amin / f))
            for alk0, basis in alk0_values(cfg):
                air = bound_curve(t, ph, alk0, T_K, pco2, scen)
                any_sink = bound_curve(t, ph, alk0, T_K, pco2, scen, EQ_PER_MOL_N_ANY_SINK)
                for co2 in [AIR_EQ, *klas]:
                    if co2 == AIR_EQ:
                        at = {name: float(air[i]) for name, i in marks}
                        label = AIR_EQ
                    else:
                        # CO2 SCENARIO at the checkpoints only, capped by the air-equilibrium curve, reverse minimum
                        vals = [min(float(air[i]), bound_with_co2(alk0, float(ph[i]), T_K, pco2, float(f[i]),
                                                                  float(t[i] - T0_D), co2)) for _, i in marks]
                        vals = np.minimum.accumulate(np.array(vals)[::-1])[::-1]
                        at = {name: float(v) for (name, _), v in zip(marks, vals)}
                        label = f"SCENARIO kLa_CO2 = {co2:g} 1/d (pseudo-steady, mean rate)"
                    n_end = at["end"]
                    row = {"alk0_meq_L": alk0, "alk0_basis": basis, "temperature_C": T,
                           "temperature_basis": f"SCENARIO {T:g} °C", "co2_scenario": label,
                           "volume_scenario": scen.name, "volume_basis": scen.basis,
                           "alk0_min_for_pH_record_meq_L": alk0_floor,
                           "consistent_with_pH_record": bool(alk0 >= alk0_floor)}
                    row.update({f"n_nit_max_{name}_mgN_L": v for name, v in at.items()})
                    n_any = float(any_sink[i_end]) if co2 == AIR_EQ else np.nan
                    tan_min = max(dosed - n_end, 0.0) * float(f[i_end])
                    tan_any = max(dosed - n_any, 0.0) * float(f[i_end]) if co2 == AIR_EQ else np.nan
                    row.update({
                        "dosed_total_mgN_L": dosed, "n_nit_max_fraction_of_dose": n_end / dosed,
                        "n_removed_max_any_proton_sink_end_mgN_L": n_any,
                        "pH_end_smoothed": float(ph[i_end]), "alk_min_end_meq_L": float(amin[i_end]),
                        "tan_residual_min_end_mgN_L": tan_min,
                        "tan_residual_min_any_proton_sink_end_mgN_L": tan_any,
                        "fa_residual_min_end_mgNH3_L": float(thermo.free_ammonia_mgNH3_L(tan_min, ph[i_end], T_K)),
                        "fa_residual_min_any_proton_sink_end_mgNH3_L":
                            float(thermo.free_ammonia_mgNH3_L(tan_any, ph[i_end], T_K)) if co2 == AIR_EQ else np.nan,
                    })
                    rows.append(row)
    return pd.DataFrame(rows)


def reference_checks(cfg: Config) -> pd.DataFrame:
    """Reference values of the remediation prompt (W2 positive tests), evaluated by the production functions."""
    T25 = 25.0 + KELVIN_OFFSET
    pco2 = float(cfg.get("scenarios.pCO2_air_atm"))
    rows = [(f"N_nit bound, Alk0 {a:g} meq/L, pH floor ignored (mg N/L)", ref,
             float(nitrification_bound_mgN_L(a)), 0.1) for a, ref in ((2.0, 14.0), (5.0, 35.0), (10.0, 70.0))]
    rows += [(f"air-equilibrium HCO3- at pH {p:g}, 25 °C (mmol/L)", ref, float(hco3_air_eq_mmol_L(p, T25, pco2)), 0.001)
             for p, ref in ((7.8, 0.401), (8.0, 0.636))]
    rows.append(("FA at 200 mg N/L TAN, pH 7.8, 25 °C (mg NH3/L)", 8.5,
                 float(thermo.free_ammonia_mgNH3_L(200.0, 7.8, T25)), 0.1))
    k = carbonate_constants(T25)
    rows += [("pK1 at 25 °C (Plummer-Busenberg vs constants.py)", 6.35, k["pK1"], 0.01),
             ("log10 K_H at 25 °C (Plummer-Busenberg vs constants.py)", -1.47, k["log10_KH"], 0.01)]
    df = pd.DataFrame(rows, columns=["quantity", "reference", "obtained", "tolerance"])
    df["abs_diff"] = (df["obtained"] - df["reference"]).abs()
    df["pass"] = df["abs_diff"] <= df["tolerance"] + 1e-9
    return df


@dataclass
class AlkalinityResult:
    bound: pd.DataFrame
    references: pd.DataFrame
    volumes: pd.DataFrame


def compute(cfg: Config, grid: pd.DataFrame, t_last_d: float | None = None) -> AlkalinityResult:
    t_end = float(grid["t_end_d"].iloc[-1])
    return AlkalinityResult(bound_table(cfg, grid), reference_checks(cfg),
                            volume.scenario_table(cfg, t_end, t_end if t_last_d is None else t_last_d))


def headline(bound: pd.DataFrame, T: float = 25.0) -> pd.DataFrame:
    """Air-equilibrium bound at constant volume, one temperature: the compact view used in the reports."""
    b = bound[(bound["co2_scenario"] == AIR_EQ) & (bound["volume_scenario"].str.startswith(volume.CONSTANT))
              & (bound["temperature_C"] == T)]
    cols = ["alk0_meq_L", "alk0_basis", "consistent_with_pH_record",
            *[c for c in b.columns if c.startswith("n_nit_max_") and c.endswith("_mgN_L")],
            "n_nit_max_fraction_of_dose", "n_removed_max_any_proton_sink_end_mgN_L", "tan_residual_min_end_mgN_L",
            "fa_residual_min_end_mgNH3_L", "tan_residual_min_any_proton_sink_end_mgN_L",
            "fa_residual_min_any_proton_sink_end_mgNH3_L"]
    return b[cols]


def write_report(cfg: Config, r: AlkalinityResult, path: Path) -> None:
    from .report import md_table

    b = r.bound
    h25, h20 = headline(b, 25.0), headline(b, 20.0)
    floor = b["alk0_min_for_pH_record_meq_L"]
    co2 = b[(b["volume_scenario"].str.startswith(volume.CONSTANT)) & (b["temperature_C"] == 25.0)].pivot_table(
        index="alk0_meq_L", columns="co2_scenario", values="n_nit_max_end_mgN_L").reset_index()
    vol = b[(b["co2_scenario"] == AIR_EQ) & (b["temperature_C"] == 25.0)].pivot_table(
        index="alk0_meq_L", columns="volume_scenario", values="n_nit_max_end_mgN_L").reset_index()
    dosed = float(b["dosed_total_mgN_L"].iloc[0])
    lines = [
        "# Proton, alkalinity and carbonate balance (W2)",
        "",
        "Generated by `python -m bnr_neq thermo`. Do not edit by hand. Implementer output, not approved.",
        "",
        "## Bound",
        "",
        ("N_nit(t) ≤ M_N (Alk0 − Alk_min(pH(s)) V(s)/V0)/2 for every s ≥ t (2 eq per mol N nitrified, nothing added "
        "but NH4Cl). Alk_min is the carbonate alkalinity at the observed pH (0.11 d rolling median of the 0.01 d grid) "
        "with dissolved CO2 at air saturation (conservative) or above it under a kLa_CO2 SCENARIO. Assumptions: no "
        "alkalinity source (denitrification, ammonification of organic N, mineral dissolution), activities equal "
        "concentrations, pH electrode accurate (the only sensor that passed M2). A proton-releasing N sink other than "
        "nitrification (NH3 volatilization, NH4+ uptake: 1 eq/mol N) would allow up to twice the removal; that "
        "looser bound is listed as `any_proton_sink`."),
        "",
        f"Alk0 source: {alk0_source(cfg)}.",
        "",
        (f"Dosed in total: {dosed:.2f} mg N/L (constant volume). The pH record requires Alk0 ≥ "
        f"{floor.min():.3f} to {floor.max():.3f} meq/L (over the temperature and volume scenarios) at air-equilibrium "
        "CO2: a lower Alk0 cannot hold the highest observed pH."),
        "",
        "### Air-equilibrium CO2, constant volume, SCENARIO 25 °C",
        "",
        md_table(h25),
        "",
        "### Air-equilibrium CO2, constant volume, SCENARIO 20 °C",
        "",
        md_table(h20),
        "",
        "Residual TAN assumes an initial TAN ≥ 0 and no N source; FA at the smoothed end-of-record pH (Anthonisen).",
        "",
        "## CO2 transfer SCENARIOS (end of record, constant volume, 25 °C)",
        "",
        ("Pseudo-steady CO2 excess = 2 N_nit / (M_N Δt kLa_CO2), with the mean rate over Δt (from day "
        f"{T0_D:g}); the bound is solved self-consistently at the checkpoints and capped by the air-equilibrium curve."),
        "",
        md_table(co2),
        "",
        "## Volume scenarios (end of record, air-equilibrium CO2, 25 °C)",
        "",
        md_table(vol),
        "",
        md_table(r.volumes),
        "",
        "## Reference checks",
        "",
        md_table(r.references),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write(cfg: Config, r: AlkalinityResult, d: dict[str, Path]) -> None:
    r.bound.to_csv(d["tables"] / "alkalinity_bound.csv", index=False)
    r.references.to_csv(d["tables"] / "alkalinity_references.csv", index=False)
    r.volumes.to_csv(d["tables"] / "volume_scenarios.csv", index=False)
    write_report(cfg, r, d["reports"] / "alkalinity.md")
