"""Thermodynamics and balances (Phase 2).

Section 6 of implementation_prompt.md. Units are in names or docstrings.
Functions that need a measurement listed under ``[missing]`` (alkalinity, DO,
chloride, ORP reference offset, kLa) have a ``*_guarded`` form that takes the
configuration and returns :class:`NotComputed` when the measurement is absent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import Config, not_computed
from .constants import (
    ALK_CACO3_PER_N_G_G,
    ANTHONISEN_EXP_K,
    DAVIES_A_25C,
    DG0_NITRIFICATION_KJ_MOL_N,
    DO_SAT_25C_MG_L,
    E0_AG2O_AG_V,
    E0_AGCL_AG_V,
    EMERSON_A,
    EMERSON_B_K,
    F_C_PER_MOL,
    KELVIN_OFFSET,
    KOHLRAUSCH_MMOL_L_TO_MS_CM,
    LAMBDA0_S_CM2_MOL,
    LOG10_KH_CO2_25C,
    M_CL_G_MOL,
    M_N_G_MOL,
    M_NH4_G_MOL,
    NH3_PER_N,
    O2_PER_N_G_G,
    P_O2_AIR_ATM,
    PK1_CARBONATE_25C,
    R_J_PER_MOL_K,
    T_REF_K,
)

CACO3_MG_PER_MEQ = 50.04  # mg CaCO3 per meq of alkalinity (100.09 / 2)

KOHLRAUSCH_NOTE = (
    "Kohlrausch predictions use limiting (infinite-dilution) ionic conductivities; "
    "at the actual ionic strength the conductivity changes are lower than these values."
)


def _require_positive(name: str, value: float) -> None:
    if not np.all(np.isfinite(value)) or np.any(np.asarray(value) <= 0):
        raise ValueError(f"{name} must be positive and finite, got {value!r}")


def _require_non_negative(name: str, value: float) -> None:
    if not np.all(np.isfinite(value)) or np.any(np.asarray(value) < 0):
        raise ValueError(f"{name} must be non-negative and finite, got {value!r}")


@dataclass(frozen=True)
class PulseDose:
    c_mmol_L: float  # mmol NH4Cl per litre
    nh4_mg_L: float  # mg NH4+/L
    n_mg_L: float  # mg N/L
    cl_mg_L: float  # mg Cl-/L


def pulse_dose(mass_mg: float, molar_mass_g_mol: float, volume_L: float) -> PulseDose:
    """c = mass / (molar mass x V); mg / (g/mol x L) = mmol/L."""
    for name, v in (("mass_mg", mass_mg), ("molar_mass_g_mol", molar_mass_g_mol), ("volume_L", volume_L)):
        _require_positive(name, v)
    c = mass_mg / (molar_mass_g_mol * volume_L)
    return PulseDose(c, c * M_NH4_G_MOL, c * M_N_G_MOL, c * M_CL_G_MOL)


def pulse_dose_from_config(cfg: Config) -> PulseDose:
    return pulse_dose(cfg.get("pulses.mass_mg"), cfg.get("pulses.molar_mass_g_mol"), cfg.get("vessel.volume_L"))


def kohlrausch_delta_ec_mS_cm(changes_mmol_L: dict[str, float]) -> float:
    """ΔEC (mS/cm) = Σ Δc_i (mmol/L) x λ_i (S cm²/mol) x 1e-3, at infinite dilution.

    ``changes_mmol_L`` maps an ion in ``LAMBDA0_S_CM2_MOL`` to its signed change.
    """
    unknown = sorted(set(changes_mmol_L) - set(LAMBDA0_S_CM2_MOL))
    if unknown:
        raise KeyError(f"no limiting conductivity for {unknown}")
    return sum(dc * LAMBDA0_S_CM2_MOL[ion] for ion, dc in changes_mmol_L.items()) * KOHLRAUSCH_MMOL_L_TO_MS_CM


def ec_step_immediate_mS_cm(c_mmol_L: float) -> float:
    """Conductivity step right after adding c mmol/L NH4Cl."""
    _require_positive("c_mmol_L", c_mmol_L)
    return kohlrausch_delta_ec_mS_cm({"NH4+": c_mmol_L, "Cl-": c_mmol_L})


def ec_step_nitrified_mS_cm(c_mmol_L: float) -> float:
    """Net change after complete nitrification: c (λCl + λNO3 − 2 λHCO3)."""
    _require_positive("c_mmol_L", c_mmol_L)
    return kohlrausch_delta_ec_mS_cm({"Cl-": c_mmol_L, "NO3-": c_mmol_L, "HCO3-": -2 * c_mmol_L})


def nernst_slope_V(T_K: float) -> float:
    """RT ln10 / F in V per decade."""
    _require_positive("T_K", T_K)
    return R_J_PER_MOL_K * T_K * np.log(10.0) / F_C_PER_MOL


def chloride_orp_rise_mV(c0_mM: float, dc_mM: float, T_K: float) -> float:
    """Predicted rise of ORP (Pt vs Ag) when chloride rises from c0 to c0 + Δc.

    If the bare Ag is a chloride-controlled AgCl/Ag electrode, E_ref falls by
    (RT ln10/F) log10((c0 + Δc)/c0), so Pt-vs-Ag rises by the same amount.
    Activity coefficients are taken as equal before and after.
    """
    _require_positive("c0_mM", c0_mM)
    _require_non_negative("dc_mM", dc_mM)
    return 1e3 * nernst_slope_V(T_K) * np.log10((c0_mM + dc_mM) / c0_mM)


def _require_ph(pH) -> None:
    a = np.asarray(pH, dtype=float)
    if not np.all(np.isfinite(a)) or np.any((a < 0) | (a > 14)):
        raise ValueError("pH must be finite and within 0 to 14")


def emerson_pka(T_K: float) -> float:
    """pKa of NH4+/NH3 (Emerson et al. 1975): 0.09018 + 2729.92 / T."""
    _require_positive("T_K", T_K)
    return EMERSON_A + EMERSON_B_K / np.asarray(T_K, dtype=float)


def free_ammonia_mgNH3_L(tan_mgN_L, pH, T_K: float):
    """Free ammonia (Anthonisen et al. 1976), mg NH3/L.

    FA = (17/14) x TAN (mg N/L) x 10^pH / (exp(6344/T) + 10^pH).
    """
    _require_non_negative("tan_mgN_L", tan_mgN_L)
    _require_ph(pH)
    _require_positive("T_K", T_K)
    tan = np.asarray(tan_mgN_L, dtype=float)
    h = 10.0 ** np.asarray(pH, dtype=float)
    return NH3_PER_N * tan * h / (np.exp(ANTHONISEN_EXP_K / T_K) + h)



# --- Guard pattern ------------------------------------------------------------
@dataclass(frozen=True)
class NotComputed:
    """Explicit result for a quantity whose input is missing."""

    key: str
    message: str = ""

    def __str__(self) -> str:
        return self.message or not_computed(self.key)


S4_NOT_COMPUTED = NotComputed("orp_reference_offset_V", "not computed: S4")


def _measured(cfg: Config, *keys: str) -> dict | NotComputed:
    for k in keys:
        if not cfg.has_measurement(k):
            return NotComputed(k)
    return {k: cfg.measurement(k) for k in keys}


# --- Activity coefficients ------------------------------------------------------
def davies_gamma(z: int, ionic_strength_M: float, A: float = DAVIES_A_25C) -> float:
    """Davies: log10 γ = −A z² (√I/(1+√I) − 0.3 I); A = 0.509 at 25 °C."""
    _require_non_negative("ionic_strength_M", ionic_strength_M)
    sq = np.sqrt(ionic_strength_M)
    return 10.0 ** (-A * z * z * (sq / (1 + sq) - 0.3 * ionic_strength_M))


# --- Bare Ag reference ------------------------------------------------------------
def e_ref_agcl_V(cl_activity_M: float, T_K: float = T_REF_K) -> float:
    """Chloride-controlled AgCl/Ag: E_ref = 0.2223 − (RT ln10/F) log10(a_Cl), V vs SHE."""
    _require_positive("cl_activity_M", cl_activity_M)
    return E0_AGCL_AG_V - nernst_slope_V(T_K) * np.log10(cl_activity_M)


def e_ref_ag2o_V(pH, T_K: float = T_REF_K):
    """Oxide-controlled Ag2O/Ag alternative: E_ref = 1.173 − (RT ln10/F) pH, V vs SHE."""
    _require_ph(pH)
    return E0_AG2O_AG_V - nernst_slope_V(T_K) * np.asarray(pH, dtype=float)


def eh_from_orp_V(orp_mV_vs_Ag, offset_V: float):
    """E_H (V vs SHE) = ORP (V vs bare Ag) + reference offset (V)."""
    return np.asarray(orp_mV_vs_Ag, dtype=float) / 1e3 + offset_V


def eh_from_orp_guarded(cfg: Config, orp_mV_vs_Ag):
    m = _measured(cfg, "orp_reference_offset_V")
    if isinstance(m, NotComputed):
        return S4_NOT_COMPUTED
    return eh_from_orp_V(orp_mV_vs_Ag, m["orp_reference_offset_V"])


def e_ref_agcl_guarded(cfg: Config, gamma_cl: float, T_K: float = T_REF_K):
    """E_ref from the measured chloride (mmol/L) and an activity coefficient."""
    m = _measured(cfg, "chloride_mmol_L")
    if isinstance(m, NotComputed):
        return m
    return e_ref_agcl_V(gamma_cl * m["chloride_mmol_L"] / 1e3, T_K)


# --- Nitrification energetics ----------------------------------------------------
def p_o2_atm(do_mg_L, c_star_mg_L: float = DO_SAT_25C_MG_L):
    """p_O2 = 0.2095 x DO / C*(T), atm. C* = 8.26 mg/L at 25 °C until temperature is measured."""
    _require_positive("DO_mg_L", do_mg_L)
    _require_positive("c_star_mg_L", c_star_mg_L)
    return P_O2_AIR_ATM * np.asarray(do_mg_L, dtype=float) / c_star_mg_L


def nitrification_dG_kJ_mol(nh4_mgN_L, no3_mgN_L, do_mg_L, pH, T_K: float = T_REF_K,
                            c_star_mg_L: float = DO_SAT_25C_MG_L):
    """In situ ΔG of NH4+ + 2 O2(g) -> NO3- + 2 H+ + H2O, kJ/mol N.

    ΔG = ΔG°′ + RT ln([NO3-]/([NH4+] p_O2²)) + 2RT ln([H+]/1e-7). NO3- and NH4+
    are both in mg N/L, so their ratio is a molar ratio. ΔG°′ (−348.8 kJ/mol N)
    is the 25 °C value; at another T only the RT terms change.
    """
    _require_positive("nh4_mgN_L", nh4_mgN_L)
    _require_positive("no3_mgN_L", no3_mgN_L)
    _require_ph(pH)
    _require_positive("T_K", T_K)
    rt = R_J_PER_MOL_K * T_K / 1e3
    po2 = p_o2_atm(do_mg_L, c_star_mg_L)
    h = 10.0 ** -np.asarray(pH, dtype=float)
    return (DG0_NITRIFICATION_KJ_MOL_N
            + rt * np.log(np.asarray(no3_mgN_L, dtype=float) / (np.asarray(nh4_mgN_L, dtype=float) * po2**2))
            + 2 * rt * np.log(h / 1e-7))


def nitrification_affinity_kJ_mol(nh4_mgN_L, no3_mgN_L, do_mg_L, pH, T_K: float = T_REF_K,
                                  c_star_mg_L: float = DO_SAT_25C_MG_L):
    """A = −ΔG, kJ/mol N."""
    return -nitrification_dG_kJ_mol(nh4_mgN_L, no3_mgN_L, do_mg_L, pH, T_K, c_star_mg_L)


def nitrification_affinity_guarded(cfg: Config, nh4_mgN_L, pH, T_K: float = T_REF_K):
    """In situ affinity: needs measured DO and NO3-N."""
    m = _measured(cfg, "DO_mg_L", "NO3_N_mg_L")
    if isinstance(m, NotComputed):
        return m
    return nitrification_affinity_kJ_mol(nh4_mgN_L, m["NO3_N_mg_L"], m["DO_mg_L"], pH, T_K)


def affinity_step_kJ_mol(c_before_mgN_L, dc_mgN_L, T_K: float = T_REF_K):
    """ΔA_k = RT ln((C_k⁻ + ΔC)/C_k⁻), kJ/mol."""
    _require_positive("c_before_mgN_L", c_before_mgN_L)
    _require_non_negative("dc_mgN_L", dc_mgN_L)
    _require_positive("T_K", T_K)
    c = np.asarray(c_before_mgN_L, dtype=float)
    return R_J_PER_MOL_K * T_K * np.log((c + dc_mgN_L) / c) / 1e3


# --- Stoichiometry ------------------------------------------------------------------
def alkalinity_consumed_mgCaCO3_L(n_nitrified_mgN_L):
    """7.14 g CaCO3 per g N nitrified."""
    _require_non_negative("n_nitrified_mgN_L", n_nitrified_mgN_L)
    return ALK_CACO3_PER_N_G_G * np.asarray(n_nitrified_mgN_L, dtype=float)


def oxygen_demand_mgO2_L_d(rate_mgN_L_d):
    """4.57 g O2 per g N nitrified."""
    _require_non_negative("rate_mgN_L_d", rate_mgN_L_d)
    return O2_PER_N_G_G * np.asarray(rate_mgN_L_d, dtype=float)


def required_kla_per_d(demand_mgO2_L_d, c_star_mg_L: float = DO_SAT_25C_MG_L):
    """Lower bound k_La >= demand / C* (supply at zero DO), 1/d."""
    _require_non_negative("demand_mgO2_L_d", demand_mgO2_L_d)
    _require_positive("c_star_mg_L", c_star_mg_L)
    return np.asarray(demand_mgO2_L_d, dtype=float) / c_star_mg_L


def kla_margin_guarded(cfg: Config, required_kla):
    """Measured k_La divided by the required k_La."""
    m = _measured(cfg, "kLa_per_d")
    if isinstance(m, NotComputed):
        return m
    return m["kLa_per_d"] / np.asarray(required_kla, dtype=float)


# --- CO2 at the open surface ---------------------------------------------------------
def hco3_air_eq_mmol_L(pH, pK1: float, log10_KH: float, pco2_atm: float):
    """[HCO3-]_eq = K1 K_H p_CO2 / [H+], mmol/L."""
    _require_ph(pH)
    _require_positive("pco2_atm", pco2_atm)
    h = 10.0 ** -np.asarray(pH, dtype=float)
    return 1e3 * 10.0 ** -pK1 * 10.0 ** log10_KH * pco2_atm / h


def omega_co2(hco3_mmol_L, pH, pK1: float, log10_KH: float, pco2_atm: float):
    """Ω_CO2 = [HCO3-][H+] / (K1 K_H p_CO2,air)."""
    _require_non_negative("hco3_mmol_L", hco3_mmol_L)
    return np.asarray(hco3_mmol_L, dtype=float) / hco3_air_eq_mmol_L(pH, pK1, log10_KH, pco2_atm)


def omega_co2_guarded(cfg: Config, pH, pK1: float, log10_KH: float):
    """Needs alkalinity; [HCO3-] ≈ alkalinity below pH 8.3 (carbonate alkalinity assumption)."""
    m = _measured(cfg, "alkalinity_mmol_L")
    if isinstance(m, NotComputed):
        return m
    return omega_co2(m["alkalinity_mmol_L"], pH, pK1, log10_KH, cfg.get("scenarios.pCO2_air_atm"))


# --- Non-equilibrium criteria ---------------------------------------------------------
def de_donder_consistent(affinity, rate) -> np.ndarray:
    """A_j r_j >= 0 (affinity and net rate share sign)."""
    return np.asarray(affinity, dtype=float) * np.asarray(rate, dtype=float) >= 0


def entropy_production_W_L_K(rates_mol_L_s, affinities_J_mol, T_K: float):
    """σ = Σ r_j A_j / T, W/(L K), for rates in mol/(L s) and affinities in J/mol."""
    _require_positive("T_K", T_K)
    r = np.atleast_1d(np.asarray(rates_mol_L_s, dtype=float))
    a = np.atleast_1d(np.asarray(affinities_J_mol, dtype=float))
    if r.shape != a.shape:
        raise ValueError("rates and affinities must have the same shape")
    return float(np.sum(r * a) / T_K)


def jin_bethke_FT(dG_J_mol: float, m_atp: float, dG_atp_J_mol: float, chi: float, T_K: float) -> float:
    """F_T = 1 − exp(−(−ΔG − m ΔG_ATP)/(χ R T)). Negative when the driving force is below m ΔG_ATP."""
    _require_positive("chi", chi)
    _require_positive("T_K", T_K)
    _require_non_negative("m_atp", m_atp)
    return float(1.0 - np.exp(-(-dG_J_mol - m_atp * dG_atp_J_mol) / (chi * R_J_PER_MOL_K * T_K)))


# ==============================================================================
# Stage (M4): tables and balances from the M1 to M3 outputs
# ==============================================================================
def _nc(x) -> object:
    """Value for a table cell: the number, or the not-computed marker string."""
    return str(x) if isinstance(x, NotComputed) else x


def load_stage_inputs(out_dir) -> dict:
    from pathlib import Path

    import pandas as pd

    files = {
        "windows": ("tables/sensor_validation_windows.csv", "validate"),
        "pulses": ("tables/sensor_validation.csv", "validate"),
        "balances": ("tables/sensor_balances.csv", "validate"),
        "segments": ("tables/segments.csv", "segments"),
        "reference_windows": ("tables/segments_reference_windows.csv", "segments"),
        "grid": ("processed/grid_0p01d.parquet", "ingest"),
    }
    out = {}
    for name, (rel, stage) in files.items():
        p = Path(out_dir) / rel
        if not p.is_file():
            raise FileNotFoundError(f"{p} missing: run `python -m bnr_neq {stage}` first")
        out[name] = pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)
    return out


# W1: the affinity step from the raw reading is a SENSITIVITY case (reading taken as TAN); the primary ΔA is computed on
# the M6-inverted C(t) in pulses.csv (`delta_A_kJ_mol_SCENARIO_T*`). The label sits in the column name itself.
DELTA_A_READING = "delta_A_SENSITIVITY_reading_as_TAN_kJ_mol"
DELTA_A_READING_ROLE = ("SENSITIVITY: electrode reading taken as TAN (no observation model); the primary ΔA uses the "
                        "M6-inverted C(t) (pulses.csv, delta_A_kJ_mol_SCENARIO_T*)")


def affinity_step_table(cfg: Config, windows, pulses):
    import pandas as pd

    dose = pulse_dose_from_config(cfg)
    nf = cfg.get("sensors.nh4_to_N_factor")
    pre = windows[(windows["sensor"] == "nh4") & (windows["window"] == "pre")].set_index("pulse")
    rows = []
    for k in pulses["pulse"]:
        c = float(pre.loc[k, "median"]) * nf
        # 25 °C is not measured: the 25 °C value is a SCENARIO like the others (computed at T_REF_K)
        row = {"pulse": int(k), "t_pulse_d": float(pulses.set_index("pulse").loc[k, "t_pulse_d"]),
               "baseline_before_mgNH4_L": float(pre.loc[k, "median"]), "baseline_before_mgN_L": c,
               "dose_mgN_L": dose.n_mg_L,
               f"{DELTA_A_READING}_SCENARIO_T25C": float(affinity_step_kJ_mol(c, dose.n_mg_L, T_REF_K))}
        for T in cfg.get("scenarios.temperature_C"):
            row[f"{DELTA_A_READING}_SCENARIO_T{T:g}C"] = float(affinity_step_kJ_mol(c, dose.n_mg_L, T + KELVIN_OFFSET))
        row["delta_A_fraction_of_dG0_SENSITIVITY_reading_as_TAN_SCENARIO_T25C"] = (
            row[f"{DELTA_A_READING}_SCENARIO_T25C"] / abs(DG0_NITRIFICATION_KJ_MOL_N))
        row["role"] = DELTA_A_READING_ROLE
        row["note"] = ("baseline = median NH4 over the 0.1 d before the pulse (electrode, uncorrected); "
                       f"a 4x gain error changes ΔA by at most RT ln 4 = {R_J_PER_MOL_K * T_REF_K * np.log(4) / 1e3:.1f} kJ/mol")
        rows.append(row)
    return pd.DataFrame(rows)


def coexistence_table(cfg: Config, segments_tab):
    """Documented coexistence SCENARIO with its DO sensitivity, and the in situ affinity guarded on DO and NO3-N.

    The in situ affinity is evaluated with the observed electrode NH4-N and pH medians of each M3 slope window
    (never a placeholder), at the 25 °C SCENARIO.
    """
    import pandas as pd

    sc = cfg.get("scenarios.coexistence")
    T_K = sc["temperature_C"] + KELVIN_OFFSET
    rows = []
    for do in sc["DO_sensitivity_mg_L"]:
        rows.append({"basis": "SCENARIO", "nh4_mgN_L": sc["nh4_mgN_L"], "no3_mgN_L": sc["no3_mgN_L"],
                     "DO_mg_L": do, "pH": sc["pH"], "temperature_C": sc["temperature_C"],
                     "p_O2_atm": float(p_o2_atm(do)),
                     "affinity_kJ_mol": float(nitrification_affinity_kJ_mol(sc["nh4_mgN_L"], sc["no3_mgN_L"], do,
                                                                            sc["pH"], T_K)),
                     "documented_scenario": bool(np.isclose(do, sc["DO_mg_L"])), "note": ""})
    comp = segments_tab[segments_tab["status"] == "computed"]
    nh4_obs = comp["median_nh4_mgN_L"].to_numpy(dtype=float)
    ph_obs = comp["median_pH"].to_numpy(dtype=float)
    a = nitrification_affinity_guarded(cfg, nh4_obs, ph_obs, T_REF_K)
    if isinstance(a, NotComputed):
        rows.append({"basis": "IN SITU", "nh4_mgN_L": np.nan,
                     "no3_mgN_L": np.nan, "DO_mg_L": np.nan, "pH": np.nan, "temperature_C": np.nan, "p_O2_atm": np.nan,
                     "affinity_kJ_mol": str(a), "documented_scenario": False,
                     "note": f"inputs: observed NH4-N and pH medians of {len(comp)} slope windows; needs DO_mg_L, NO3_N_mg_L"})
    else:
        m = cfg.raw["measurements"]
        rows += [{"basis": "IN SITU", "nh4_mgN_L": c, "no3_mgN_L": m["NO3_N_mg_L"],
                  "DO_mg_L": m["DO_mg_L"], "pH": ph, "temperature_C": T_REF_K - KELVIN_OFFSET,
                  "p_O2_atm": float(p_o2_atm(m["DO_mg_L"])), "affinity_kJ_mol": float(v), "documented_scenario": False,
                  "note": f"sub-segment {int(sub)}: observed NH4-N and pH medians; measured DO, NO3-N; SCENARIO 25 °C"}
                 for sub, c, ph, v in zip(comp["sub_segment"], nh4_obs, ph_obs, np.atleast_1d(a))]
    return pd.DataFrame(rows)


def chloride_table(cfg: Config, pulses, segments_tab, windows):
    """Chloride SCENARIOS, observed changes, and absolute E_H guarded on the reference offset.

    E_H is evaluated on the observed ORP medians 1 to 3 h after each pulse (M2 ``post`` window), never on a
    placeholder value.
    """
    import pandas as pd

    dose = pulse_dose_from_config(cfg)
    n = len(cfg.pulse_times_d)
    rows = []
    for T in cfg.get("scenarios.temperature_C"):
        T_K = T + KELVIN_OFFSET
        for c0 in cfg.get("scenarios.baseline_chloride_mM"):
            row = {"scenario": f"SCENARIO_Cl0_{c0:g}mM_T{T:g}C", "baseline_chloride_mM": c0, "temperature_C": T}
            for k in range(1, n + 1):
                row[f"orp_pred_rise_mV_pulse{k}"] = chloride_orp_rise_mV(c0 + (k - 1) * dose.c_mmol_L, dose.c_mmol_L, T_K)
            row["orp_pred_rise_mV_cumulative"] = chloride_orp_rise_mV(c0, n * dose.c_mmol_L, T_K)
            row["e_ref_V_vs_SHE_start_gamma1"] = e_ref_agcl_V(c0 / 1e3, T_K)
            row["e_ref_V_vs_SHE_end_gamma1"] = e_ref_agcl_V((c0 + n * dose.c_mmol_L) / 1e3, T_K)
            rows.append(row)
    scen = pd.DataFrame(rows)
    obs = pulses[["pulse", "orp_obs_change_mV", "ph_obs_change"]].copy()
    obs["oxide_scenario_pred_orp_change_mV_SCENARIO_T25C"] = 1e3 * nernst_slope_V(T_REF_K) * obs["ph_obs_change"]
    # gamma_Cl = 1 is the stated activity assumption of the chloride SCENARIOS (report, assumptions)
    obs["e_ref_measured_chloride_V"] = _nc(e_ref_agcl_guarded(cfg, 1.0))
    post = windows[(windows["sensor"] == "orp") & (windows["window"] == "post")].set_index("pulse")["median"]
    obs["orp_post_median_mV_vs_Ag"] = obs["pulse"].map(post).astype(float)
    eh = eh_from_orp_guarded(cfg, obs["orp_post_median_mV_vs_Ag"].to_numpy())
    obs["eh_absolute"] = str(eh) if isinstance(eh, NotComputed) else eh
    med_ph = float(np.nanmedian(segments_tab.loc[segments_tab["status"] == "computed", "median_pH"]))
    oxide = {"median_pH_of_slope_windows": med_ph,
             "e_ref_ag2o_V_vs_SHE_SCENARIO_T25C": float(e_ref_ag2o_V(med_ph, T_REF_K))}
    return scen, obs, oxide


def co2_table(cfg: Config, segments_tab):
    import pandas as pd

    pk1, lkh, pco2 = PK1_CARBONATE_25C, LOG10_KH_CO2_25C, cfg.get("scenarios.pCO2_air_atm")
    rows = [{"basis": "reference pH", "sub_segment": pd.NA, "pH": ph,
             "hco3_air_eq_mmol_L_SCENARIO_T25C": float(hco3_air_eq_mmol_L(ph, pk1, lkh, pco2))} for ph in (7.3, 7.8)]
    comp = segments_tab[segments_tab["status"] == "computed"]
    rows += [{"basis": "median pH of slope window", "sub_segment": int(r.sub_segment), "pH": r.median_pH,
              "hco3_air_eq_mmol_L_SCENARIO_T25C": float(hco3_air_eq_mmol_L(r.median_pH, pk1, lkh, pco2))}
             for r in comp.itertuples(index=False)]
    df = pd.DataFrame(rows)
    df["temperature"] = "SCENARIO 25 °C (pK1, K_H at 25 °C)"
    # Ω_CO2 at each row's own pH (observed for the slope windows), guarded on alkalinity
    df["omega_co2_SCENARIO_T25C"] = [_nc(omega_co2_guarded(cfg, ph, pk1, lkh)) for ph in df["pH"]]
    return df


OXYGEN_SENSITIVITY_ROLE = ("SENSITIVITY: linear electrode gain model (W1); the primary O2 demand comes from the M6 "
                           "log observation model in nh4_truth_windows.csv")
# X4: every reading-based column of thermo_oxygen.csv says SENSITIVITY in its name
OXY_CONS = "consumption_SENSITIVITY_linear_gain_mgN_L_d"
OXY_DEMAND = "o2_demand_SENSITIVITY_linear_gain_mgO2_L_d"
OXY_KLA = "kla_required_SENSITIVITY_linear_gain_per_d_SCENARIO_T25C"
OXY_MARGIN = "kla_measured_over_required_SENSITIVITY_linear_gain"
OXYGEN_NH4_INDEPENDENT = ("sub_segment", "fit_start_d", "fit_end_d", "governing_pulse", "held_fraction", "status", "role")


def oxygen_table(cfg: Config, segments_tab, reference_windows=None, pulses=None):
    """O2 demand from the GAIN-CORRECTED consumption: a labelled SENSITIVITY case since W1."""
    import pandas as pd

    from .model import LINEAR_GAIN, RAW_READING, check_single_observation_model
    from .segments import GAIN_COL

    rows = []
    for r in segments_tab[segments_tab["status"] == "computed"].itertuples(index=False):
        rate = getattr(r, GAIN_COL)
        row = {"sub_segment": int(r.sub_segment), "fit_start_d": r.fit_start_d, "fit_end_d": r.fit_end_d,
               "governing_pulse": r.governing_pulse, OXY_CONS: rate, "held_fraction": r.held_fraction,
               OXY_DEMAND: np.nan, OXY_KLA: np.nan, "status": ""}
        if pd.isna(rate):
            row["status"] = r.gain_correction_status
        elif rate <= 0:
            row["status"] = "not computed: non-positive gain-corrected consumption (rising electrode signal)"
        else:
            d = float(oxygen_demand_mgO2_L_d(rate))
            row.update({OXY_DEMAND: d, OXY_KLA: float(required_kla_per_d(d)), "status": "computed"})
        row[OXY_MARGIN] = (_nc(kla_margin_guarded(cfg, row[OXY_KLA])) if row["status"] == "computed" else np.nan)
        row["role"] = OXYGEN_SENSITIVITY_ROLE
        rows.append(row)
    if reference_windows is not None and pulses is not None:
        # Diagnostic only: the preliminary 2.59-8.82 window (spans pulse 1) divided by pulse 1's response fraction.
        w = reference_windows.set_index("ref_from_d").loc[2.59]
        f1 = float(pulses.set_index("pulse").loc[1, "nh4_response_fraction"])
        rate = -float(w["exact_window_slope_mgNH4_L_d"]) * cfg.get("sensors.nh4_to_N_factor") / f1
        d = float(oxygen_demand_mgO2_L_d(rate))
        rows.append({"sub_segment": pd.NA, "fit_start_d": 2.69, "fit_end_d": 8.82, "governing_pulse": 1,
                     OXY_CONS: rate, "held_fraction": np.nan, OXY_DEMAND: d, OXY_KLA: float(required_kla_per_d(d)),
                     "status": "DIAGNOSTIC ONLY: preliminary window 2.59-8.82 spans pulse 1 and three NH4 jumps; "
                               "excluded by the M3 segmentation rule",
                     OXY_MARGIN: np.nan, "role": OXYGEN_SENSITIVITY_ROLE})
    out = pd.DataFrame(rows)
    out["nh4_observation_model"], out["nh4_series"] = LINEAR_GAIN, RAW_READING  # X3: declared basis, checked
    check_single_observation_model(out, label_cols=("role",), independent_cols=OXYGEN_NH4_INDEPENDENT)
    return out


NH4_BALANCE_START_D, NH4_BALANCE_END_D = 0.8, 47.4  # midpoints of the M2 start and end windows


def balance_tables(cfg: Config, sensor_balances, volumes=None, alk=None):
    """Nitrogen, chloride and alkalinity balances. ``volumes`` (list of :class:`volume.VolumeScenario`) carries
    V(t) through each concentration balance (W8): the closed-frame rows are labelled with the constant-volume
    scenario, and every other scenario adds amount-basis rows. ``alk`` (``alkalinity.AlkalinityResult``) adds
    the W2 bound rows to the alkalinity balance."""
    import pandas as pd

    from . import alkalinity, volume

    nitrogen, chloride, alk_tab = _closed_frame_balances(cfg, sensor_balances)
    for df in (nitrogen, chloride, alk_tab):
        df["volume_scenario"] = volume.CONSTANT
    b = sensor_balances.set_index("quantity")["value"]
    n_mg = volume.pulse_n_mg(cfg)
    n = len(cfg.pulse_times_d)
    v0 = float(cfg.get("vessel.volume_L"))
    dose = pulse_dose_from_config(cfg)
    extra_n, extra_cl, extra_alk = [], [], []
    for s in volumes or []:
        if s.is_constant:
            continue
        fs, fe = float(s.factor(NH4_BALANCE_START_D)), float(s.factor(NH4_BALANCE_END_D))
        removed = b["nh4_median_start_mgN_L"] / fs + n * n_mg / v0 - b["nh4_median_end_mgN_L"] / fe
        extra_n += [
            ("N added (4 pulses), sum of the concentration steps at each pulse",
             sum(volume.dose_mgN_L(cfg, s, tp) for tp in cfg.pulse_times_d), "mg N/L", s.basis, s.name),
            ("N added (4 pulses), concentration at the end if nothing were removed", n * n_mg / v0 * fe, "mg N/L",
             s.basis, s.name),
            ("Removed from the NH4 pool, electrode (uncorrected), amount basis per initial litre", removed, "mg N/L",
             "derived: (C_start V_start + added − C_end V_end) / V0; NO-GO as an absolute value (S1, S2)", s.name)]
        extra_cl += [(f"Final chloride, SCENARIO baseline {c0:g} mM", (c0 + n * dose.c_mmol_L) * fe, "mmol/L",
                      "SCENARIO; amount / V_end", s.name) for c0 in cfg.get("scenarios.baseline_chloride_mM")]
        extra_alk.append(("Consumed if all added N were nitrified, concentration at the end",
                          float(alkalinity_consumed_mgCaCO3_L(n * dose.n_mg_L)) / CACO3_MG_PER_MEQ * fe, "meq/L",
                          "amount per initial litre x V0/V_end", s.name))
    if alk is not None:
        bnd = alk.bound
        base = bnd[(bnd["co2_scenario"] == alkalinity.AIR_EQ) & (bnd["temperature_C"] == T_REF_K - KELVIN_OFFSET)
                   & (bnd["alk0_basis"] != "MEASURED")]
        for name, g in base.groupby("volume_scenario", sort=False):
            extra_alk.append(("Minimum initial alkalinity that can hold the observed pH (air-equilibrium CO2)",
                              float(g["alk0_min_for_pH_record_meq_L"].iloc[0]), "meq/L",
                              "W2: max over t of Alk_min(pH(t)) V(t)/V0; SCENARIO 25 °C", name))
            for a0 in (2.0, 5.0, 10.0, 15.0):
                r = g[g["alk0_meq_L"] == a0].iloc[0]
                extra_alk.append((f"Upper bound on cumulative nitrification by the end, Alk0 = {a0:g} meq/L",
                                  float(r["n_nit_max_end_mgN_L"]), "mg N/L",
                                  "W2 bound, air-equilibrium CO2 (conservative), SCENARIO 25 °C; alkalinity_bound.csv",
                                  name))
        for name, g in bnd[bnd["alk0_basis"] == "MEASURED"].groupby("volume_scenario", sort=False):
            r = g[(g["co2_scenario"] == alkalinity.AIR_EQ) & (g["temperature_C"] == T_REF_K - KELVIN_OFFSET)].iloc[0]
            extra_alk.append(("Upper bound on cumulative nitrification by the end, measured Alk0",
                              float(r["n_nit_max_end_mgN_L"]), "mg N/L", "W2 bound with initial_alkalinity_meq_L", name))
    cols = ["term", "value", "unit", "basis", "volume_scenario"]
    nitrogen = pd.concat([nitrogen, pd.DataFrame(extra_n, columns=cols)], ignore_index=True)
    chloride = pd.concat([chloride, pd.DataFrame(extra_cl, columns=cols)], ignore_index=True)
    alk_tab = pd.concat([alk_tab, pd.DataFrame(extra_alk, columns=cols)], ignore_index=True)
    return nitrogen, chloride, alk_tab


def _closed_frame_balances(cfg: Config, sensor_balances):
    import pandas as pd

    b = sensor_balances.set_index("quantity")["value"]
    dose = pulse_dose_from_config(cfg)
    n = len(cfg.pulse_times_d)
    added_n = n * dose.n_mg_L
    nc = lambda key: str(NotComputed(key))
    nitrogen = pd.DataFrame([
        ("NH4-N start, electrode (uncorrected), day 0.3-1.3", b["nh4_median_start_mgN_L"], "mg N/L", "measured (electrode)"),
        ("NH4-N start, laboratory", nc("lab_NH4_N_mg_L"), "mg N/L", ""),
        ("N added (4 pulses)", added_n, "mg N/L", "known dose"),
        ("NH4-N end, electrode (uncorrected), day 46.9-47.9", b["nh4_median_end_mgN_L"], "mg N/L", "measured (electrode)"),
        ("NO3-N end", nc("NO3_N_mg_L"), "mg N/L", ""),
        ("NO2-N end", nc("NO2_N_mg_L"), "mg N/L", ""),
        ("NH3 volatilized through the surface", nc("kLa_per_d"), "mg N/L", "needs a measured surface transfer coefficient"),
        ("N assimilated into biomass", "not computed: no biomass measurement (MLVSS is not among the [missing] keys)",
         "mg N/L", ""),
        ("N2 / N2O loss (closure term)", nc("NO3_N_mg_L"), "mg N/L", "requires NO3-N and NO2-N"),
        ("Removed from the NH4 pool, electrode (uncorrected)", b["n_removed_from_nh4_pool_mgN_L"], "mg N/L",
         "derived: start + added − end; NO-GO as an absolute value (S1, S2)"),
        ("Balance closure", nc("NO3_N_mg_L"), "mg N/L", ""),
    ], columns=["term", "value", "unit", "basis"])
    cl_rows = [
        ("Baseline chloride, measured", nc("chloride_mmol_L"), "mmol/L", ""),
        ("Chloride added (4 pulses)", n * dose.c_mmol_L, "mmol/L", "known dose; conservative"),
        ("Chloride added (4 pulses)", n * dose.cl_mg_L, "mg/L", "known dose"),
        ("Final chloride, measured basis", nc("chloride_mmol_L"), "mmol/L", ""),
    ]
    cl_rows += [(f"Final chloride, SCENARIO baseline {c0:g} mM", c0 + n * dose.c_mmol_L, "mmol/L", "SCENARIO")
                for c0 in cfg.get("scenarios.baseline_chloride_mM")]
    cl_rows.append(("Check against conductivity", ("not computed: S3 (conductivity record not valid as solution "
                    "conductivity)"), "", ""))
    chloride = pd.DataFrame(cl_rows, columns=["term", "value", "unit", "basis"])
    alk_all = float(alkalinity_consumed_mgCaCO3_L(added_n))
    alk_pulse = float(alkalinity_consumed_mgCaCO3_L(dose.n_mg_L))
    alk_removed = float(alkalinity_consumed_mgCaCO3_L(b["n_removed_from_nh4_pool_mgN_L"]))
    alk0 = (cfg.measurement("initial_alkalinity_meq_L") if cfg.has_measurement("initial_alkalinity_meq_L")
            else nc("initial_alkalinity_meq_L"))
    from .alkalinity import alk0_source

    alk_any = cfg.measurement("alkalinity_mmol_L") if cfg.has_measurement("alkalinity_mmol_L") else nc("alkalinity_mmol_L")
    alkalinity = pd.DataFrame([
        ("Alkalinity, sampling time not recorded (alkalinity_mmol_L)", alk_any, "mmol/L",
         "measured" if not isinstance(alk_any, str) else ""),
        ("Initial alkalinity (W2 key)", alk0, "meq/L", "measured" if not isinstance(alk0, str) else ""),
        ("Alk0 used by the W2 bound", alk0_source(cfg), "", "W2: alkalinity_bound.csv"),
        ("Consumed per fully nitrified pulse", alk_pulse, "mg CaCO3/L", "7.14 g CaCO3 per g N"),
        ("Consumed per fully nitrified pulse", alk_pulse / CACO3_MG_PER_MEQ, "meq/L", ""),
        ("Consumed if all added N were nitrified", alk_all, "mg CaCO3/L", "upper bound from the dose"),
        ("Consumed if all added N were nitrified", alk_all / CACO3_MG_PER_MEQ, "meq/L", ""),
        ("Consumed if the electrode-implied removal were nitrified", alk_removed, "mg CaCO3/L",
         "conditional on the uncorrected electrode balance"),
        ("Remaining alkalinity", nc("alkalinity_mmol_L"), "mmol/L", ""),
        ("CO2 supersaturation Ω_CO2", nc("alkalinity_mmol_L"), "-", ""),
    ], columns=["term", "value", "unit", "basis"])
    return nitrogen, chloride, alkalinity


M4_REFERENCES = [
    *[(f"delta_A_kJ_mol pulse {k} (SENSITIVITY: reading as TAN, as in the preliminary analysis)", v, 0.1)
      for k, v in enumerate([4.2, 2.5, 3.1, 2.4], start=1)],
    ("coexistence affinity kJ/mol (documented SCENARIO)", 317.0, 2.0),
    ("HCO3- air-equilibrium mmol/L at pH 7.3", 0.127, 0.001),
    ("HCO3- air-equilibrium mmol/L at pH 7.8", 0.401, 0.001),
    ("O2 demand min mg O2/L/d (SENSITIVITY: linear gain, as in the preliminary analysis)", 6.0, 0.6),
    ("O2 demand max mg O2/L/d (SENSITIVITY: linear gain, as in the preliminary analysis)", 33.0, 3.3),
    ("kLa required min 1/d (SENSITIVITY: linear gain, as in the preliminary analysis)", 0.7, 0.07),
    ("kLa required max 1/d (SENSITIVITY: linear gain, as in the preliminary analysis)", 4.0, 0.4),
]


def reference_comparison(steps, coex, co2, oxy):
    import pandas as pd

    ok = oxy[oxy["status"] == "computed"]
    vals = [*steps[f"{DELTA_A_READING}_SCENARIO_T25C"].tolist(),
            float(coex.loc[coex["documented_scenario"], "affinity_kJ_mol"].iloc[0]),
            *co2.loc[co2["basis"] == "reference pH", "hco3_air_eq_mmol_L_SCENARIO_T25C"].tolist(),
            ok[OXY_DEMAND].min(), ok[OXY_DEMAND].max(), ok[OXY_KLA].min(), ok[OXY_KLA].max()]
    df = pd.DataFrame(M4_REFERENCES, columns=["quantity", "reference", "tolerance"])
    df["obtained"] = vals
    df["abs_diff"] = (df["obtained"] - df["reference"]).abs()
    df["pass"] = df["abs_diff"] <= df["tolerance"] + 1e-9
    df["note"] = ""
    df.loc[df["quantity"].str.startswith(("O2", "kLa")), "note"] = (
        "tolerance 10 % (reference says 'about'); range over slope windows with positive gain-corrected consumption")
    return df[["quantity", "reference", "obtained", "tolerance", "abs_diff", "pass", "note"]]


def write_report(cfg: Config, res: dict, path) -> None:
    from .report import md_table

    refs = res["references"]
    n_fail = int((~refs["pass"]).sum())
    lines = [
        "# Thermodynamics and conservation balances (M4)",
        "",
        "Generated by `python -m bnr_neq thermo`. Do not edit by hand. Implementer output, not approved.",
        "",
        ("Constants: 25 °C values (`constants.py`) unless another temperature SCENARIO is named; R = 8.314462618 J/(mol K). "
        "Temperature is not measured (`temperature_C_measured` missing), so every temperature-dependent value is a "
        "SCENARIO value: columns computed with 25 °C constants carry `SCENARIO_T25C` (or a `SCENARIO 25 °C` basis), "
        "the others name their SCENARIO temperature. ΔG°′ is a 25 °C value; at 20 °C only the RT terms change."),
        "",
        "## Affinity step at each pulse",
        "",
        md_table(res["steps"].drop(columns=["note", "role"])),
        "",
        f"Note: {res['steps']['note'].iloc[0]}. Role: {DELTA_A_READING_ROLE}.",
        "",
        "## Coexistence affinity (documented SCENARIO) and DO sensitivity",
        "",
        ("A = −ΔG for NH4+ + 2 O2(g) → NO3- + 2 H+ + H2O, with p_O2 = 0.2095 x DO / 8.26 mg/L (C* at 25 °C). "
        "The in situ value needs measured DO and NO3-N."),
        "",
        md_table(res["coexistence"]),
        "",
        "## Bare Ag reference: chloride SCENARIOS",
        "",
        ("Predicted ORP (Pt vs Ag) rise per pulse if the Ag surface were chloride-controlled, with the chloride before "
        "pulse k equal to Cl0 + (k − 1) x dose; E_ref with activity = concentration (γ = 1). Absolute E_H is "
        "`not computed: S4` because `orp_reference_offset_V` is missing."),
        "",
        md_table(res["chloride_scenarios"], digits=4),
        "",
        (f"Observed changes and the oxide-controlled alternative (E_ref = {E0_AG2O_AG_V} − (RT ln10/F) pH: Pt-vs-Ag would "
        f"rise by {1e3 * nernst_slope_V(T_REF_K):.2f} mV per pH unit at 25 °C if the Pt potential did not itself depend on pH):"),
        "",
        md_table(res["chloride_observed"]),
        "",
        (f"Oxide-controlled E_ref at the median pH of the slope windows ({res['oxide']['median_pH_of_slope_windows']:.2f}): "
        f"{res['oxide']['e_ref_ag2o_V_vs_SHE_SCENARIO_T25C']:.3f} V vs SHE, against "
        f"{res['chloride_scenarios']['e_ref_V_vs_SHE_start_gamma1'].min():.3f} to "
        f"{res['chloride_scenarios']['e_ref_V_vs_SHE_start_gamma1'].max():.3f} V for the chloride SCENARIOS: the surface state alone moves "
        "E_ref by hundreds of mV, so E_H is undetermined."),
        "",
        "## CO2 at the open surface",
        "",
        ("[HCO3-]_eq = K1 K_H p_CO2,air / [H+] with pK1 = 6.35, log10 K_H = −1.47 (SCENARIO 25 °C), p_CO2 = 4.2e-4 atm. "
        "Ω_CO2 is guarded on alkalinity and evaluated at each row's pH."),
        "",
        md_table(res["co2"]),
        "",
        "## SENSITIVITY: oxygen demand under the linear electrode gain model (GAIN-CORRECTED, S1)",
        "",
        ("Demand = 4.57 x gain-corrected consumption; required k_La ≥ demand / C*, C* = 8.26 mg/L (SCENARIO 25 °C, to be "
        "replaced by a temperature function). Since W1 the single NH4 observation model is the M6 log model; this table "
        "is a labelled sensitivity case and the primary O2 demand is in `nh4_truth_windows.csv` (pulses stage). Oxygen "
        "claims are NO-GO until DO and k_La are measured."),
        "",
        md_table(res["oxygen"]),
        "",
        "## Balances",
        "",
        "### Nitrogen",
        "",
        md_table(res["nitrogen"], digits=5),
        "",
        "### Chloride",
        "",
        md_table(res["chloride_balance"], digits=5),
        "",
        "### Alkalinity",
        "",
        ("Rows labelled with a volume scenario other than the constant frame carry V(t) (W8); the W2 bound is detailed in "
        "`alkalinity.md` and `alkalinity_bound.csv`."),
        "",
        md_table(res["alkalinity"], digits=5),
        "",
        "## Reference comparison",
        "",
        f"{len(refs) - n_fail} of {len(refs)} reproduced within tolerance.",
        "",
        md_table(refs),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def compute_stage(cfg: Config, inputs: dict) -> dict:
    """M4 tables. With ``inputs["grid"]`` (the stage always passes it) the W2 alkalinity bound and the W8 volume
    scenarios are computed and carried into the balances."""
    from . import alkalinity, volume

    alk = vols = None
    if "grid" in inputs:
        alk = alkalinity.compute(cfg, inputs["grid"])
        vols = volume.volume_scenarios(cfg, float(inputs["grid"]["t_end_d"].iloc[-1]))
    steps = affinity_step_table(cfg, inputs["windows"], inputs["pulses"])
    coex = coexistence_table(cfg, inputs["segments"])
    scen, obs, oxide = chloride_table(cfg, inputs["pulses"], inputs["segments"], inputs["windows"])
    co2 = co2_table(cfg, inputs["segments"])
    oxy = oxygen_table(cfg, inputs["segments"], inputs["reference_windows"], inputs["pulses"])
    nitrogen, chloride, alk_tab = balance_tables(cfg, inputs["balances"], vols, alk)
    return {"steps": steps, "coexistence": coex, "chloride_scenarios": scen, "chloride_observed": obs,
            "oxide": oxide, "co2": co2, "oxygen": oxy, "nitrogen": nitrogen, "chloride_balance": chloride,
            "alkalinity": alk_tab, "references": reference_comparison(steps, coex, co2, oxy), "alk": alk}


def run(cfg: Config, out_dir) -> dict:
    from pathlib import Path

    from .pipeline import out_dirs

    inputs = load_stage_inputs(Path(out_dir))  # fail before creating any output directory
    d = out_dirs(out_dir)
    res = compute_stage(cfg, inputs)
    names = {"steps": "thermo_affinity_steps", "coexistence": "thermo_coexistence",
             "chloride_scenarios": "thermo_chloride_scenarios", "chloride_observed": "thermo_chloride_observed",
             "co2": "thermo_co2", "oxygen": "thermo_oxygen", "nitrogen": "balance_nitrogen",
             "chloride_balance": "balance_chloride", "alkalinity": "balance_alkalinity",
             "references": "thermo_references"}
    for key, name in names.items():
        res[key].to_csv(d["tables"] / f"{name}.csv", index=False)
    from . import alkalinity

    alkalinity.write(cfg, res["alk"], d)
    write_report(cfg, res, d["reports"] / "thermo.md")
    return res
