"""Physical and chemical constants, with units and sources.

Values are at 25 °C unless a temperature function is given. Experimental facts
(dose, volume, sensor calibrations) are not constants; they live in
``config/experiment.toml``.
"""

# --- Fundamental (CODATA 2018, exact SI values) -------------------------------
R_J_PER_MOL_K = 8.314462618  # J/(mol K), molar gas constant (CODATA 2018)
F_C_PER_MOL = 96485.33212  # C/mol, Faraday constant (CODATA 2018)
T_REF_K = 298.15  # K, 25 °C reference temperature
KELVIN_OFFSET = 273.15  # K, 0 °C
SECONDS_PER_DAY = 86400.0  # s/d

# --- Molar masses (IUPAC 2013 standard atomic weights, conventional) ---------
M_N_G_MOL = 14.007  # g/mol, nitrogen
M_NH4_G_MOL = 18.039  # g/mol, NH4+ (14.007 + 4 x 1.008)
M_CL_G_MOL = 35.453  # g/mol, chlorine

# --- Limiting ionic molar conductivities at 25 °C, S cm^2/mol ----------------
# CRC Handbook of Chemistry and Physics, "Equivalent conductivity of
# electrolytes in aqueous solution"; values as specified in
# implementation_prompt.md, Section 6.2. Infinite-dilution values: at the
# actual ionic strength the conductivity changes are lower.
LAMBDA0_S_CM2_MOL = {
    "NH4+": 73.5,
    "K+": 73.5,
    "Na+": 50.1,
    "Cl-": 76.3,
    "NO3-": 71.4,
    "HCO3-": 44.5,
}
# mmol/L x S cm^2/mol x 1e-3 = mS/cm (1 mmol/L = 1e-6 mol/cm^3; S/cm -> mS/cm x 1e3)
KOHLRAUSCH_MMOL_L_TO_MS_CM = 1e-3

# --- Silver reference electrode ----------------------------------------------
# AgCl/Ag standard potential at 25 °C, V vs SHE (Bard and Faulkner,
# Electrochemical Methods, Appendix C).
E0_AGCL_AG_V = 0.2223
# Ag2O/Ag standard potential at 25 °C, V vs SHE, for the oxide-controlled
# alternative E_ref = 1.173 - (RT ln10/F) pH (Bard and Faulkner; Section 6.3).
E0_AG2O_AG_V = 1.173

# --- Ammonia speciation -------------------------------------------------------
# Emerson et al. (1975): pKa = 0.09018 + 2729.92 / T(K).
EMERSON_A = 0.09018
EMERSON_B_K = 2729.92  # K
# Anthonisen et al. (1976): FA = (17/14) TAN 10^pH / (exp(6344/T) + 10^pH).
ANTHONISEN_EXP_K = 6344.0  # K
NH3_PER_N = 17.0 / 14.0  # g NH3 per g N (Anthonisen's rounded ratio)

# --- Activity coefficients ----------------------------------------------------
DAVIES_A_25C = 0.509  # (L/mol)^0.5, Davies constant at 25 °C (Stumm and Morgan)

# --- Nitrification energetics (Rittmann and McCarty, electron-equivalent tables)
DG0_NITRIFICATION_KJ_MOL_N = -348.8  # kJ/mol N, NH4+ + 2 O2 -> NO3- + 2 H+ + H2O, pH 7
DG0_NITRITATION_KJ_MOL_N = -274.7  # kJ/mol N, NH4+ -> NO2-, pH 7
DG0_NITRATATION_KJ_MOL_N = -74.1  # kJ/mol N, NO2- -> NO3-, pH 7
P_O2_AIR_ATM = 0.2095  # atm, O2 partial pressure of dry air at 1 atm

# --- Stoichiometry (Metcalf and Eddy, Wastewater Engineering) ----------------
ALK_CACO3_PER_N_G_G = 7.14  # g CaCO3 consumed per g N nitrified
O2_PER_N_G_G = 4.57  # g O2 required per g N nitrified
O2_PER_N_NITRITATION_G_G = 3.43  # g O2 per g N oxidized to nitrite only (W4 sensitivity)
J_PER_KWH = 3.6e6
DO_SAT_25C_MG_L = 8.26  # mg/L, air-saturated DO at 25 °C, 1 atm, fresh water

# --- Carbonate system at 25 °C (Stumm and Morgan) ----------------------------
PK1_CARBONATE_25C = 6.35  # -log10 K1 for CO2(aq) + H2O = HCO3- + H+
LOG10_KH_CO2_25C = -1.47  # log10 K_H, mol/(L atm)

# --- Reduced model (M6) -----------------------------------------------------------------
# Ammonium half-saturation of autotrophs, ASM1/ASM3 default K_NH = 1.0 g N/m3 (Henze et al. 2000).
# Fixed (not fitted) in the biomass-decay rate law, because the electrode cannot identify it (M3).
K_NH_ASM_MG_N_L = 1.0
# Alkalinity half-saturation of autotrophs, ASM3 K_A,ALK = 0.5 mol HCO3-/m3 (Gujer et al. 1999, Water Sci. Technol.
# 39(1), 183-193), taken as 0.5 meq/L of carbonate alkalinity (pH < 8.3). Fixed (not fitted) in the W5 candidates.
K_ALK_ASM3_MEQ_L = 0.5
EQ_ALK_PER_MOL_N_NITRIFIED = 2.0  # 2 eq of alkalinity consumed per mol N nitrified (7.14 g CaCO3/g N)
