# Aeration Efficiency in Wastewater Treatment: Literature Summary

**Target unit:** kg O₂/kWh  
**Prepared:** October 27, 2026  
**Focus:** Oxygen transfer per unit of aeration energy, with emphasis on biological nutrient removal (BNR) applications.

## Definition and interpretation

Aeration efficiency is the mass of oxygen transferred to the liquid per unit of electrical energy consumed:

$$
\mathrm{AE\;(kg\ O_2/kWh)}=\frac{\mathrm{oxygen\ transfer\ rate\;(kg\ O_2/h)}}{\mathrm{aeration\ power\;(kW)}}.
$$

**Standard aeration efficiency (SAE)** normally refers to clean-water reference test conditions (commonly normalized to 20 °C, zero dissolved oxygen, and standard pressure). **Field or process aeration efficiency (AE)** reflects oxygen transfer under actual wastewater operating conditions. These are **not directly interchangeable**: wastewater characteristics, diffuser fouling, tank geometry, dissolved oxygen, air delivery, and blower efficiency can alter field performance. Comparisons should also establish whether reported energy includes blowers, mixers, and other auxiliaries.

## 1. Engineering reference values

### 1.1 Standard (clean-water) efficiency

Source: Ontario Ministry of the Environment, *Water and Energy Conservation Guidance Manual for Sewage Works*, Table 5.1, adapted from WEF/ASCE (1998). [Source](https://www.ontario.ca/document/water-and-energy-conservation-guidance-manual-sewage-works/energy-conservation-measures)

| Aeration technology | Standard aeration efficiency (kg O₂/kWh) |
|:--|--:|
| Fine-bubble diffusers | 1.9–6.6 |
| Fine-bubble jet aerators | 2.2–3.5 |
| Mechanical surface aerators | 1.1–2.5 |
| Coarse-bubble diffusers | 1.3–1.9 |

### 1.2 Standard versus field efficiency

Source: Crawford and Sandino (2010), *Energy Efficiency in Wastewater Treatment in North America: A Compendium of Best Practices and Case Studies of Novel Approaches*, WERF, Table 2-8, drawing on WEF (2009). [Report PDF](https://betterbuildingssolutioncenter.energy.gov/sites/default/files/tools/EnergyEfficiency_WWTreatment_BMPs.pdf)

| Aeration technology | Standard (kg O₂/kWh) | Field (kg O₂/kWh) |
|:--|--:|--:|
| Fine-pore diffusion | 5.0–6.5 | 2.5–3.5 |
| Coarse-bubble diffusion | 2.5–3.5 | 1.0–2.0 |
| Low-speed centrifugal surface aeration | 1.2–3.0 | 0.7–1.4 |
| High-speed axial surface aeration | 1.2–2.2 | 0.7–1.2 |
| Submerged turbine with sparger | 1.2–2.0 | 0.7–1.1 |

*Note:* Differences between the Ontario and WERF standard ranges reflect differences in the reference compilations and equipment categories; they should not be interpreted as conflicting measurements of a single system.

## 2. Research studies with numerical values

| Reference | System and basis | Reported efficiency (kg O₂/kWh) | Interpretation |
|:--|:--|--:|:--|
| **Baeten et al. (2021)** | Full-scale aerobic granular sludge reactor; off-gas-based estimate under process conditions | **1.53** | Average actual-process AE; highly relevant full-scale wastewater comparison |
| **Côté et al. (2015)** | Pilot-scale membrane-aerated biofilm reactor (MABR), municipal primary effluent | **3–5** | Aeration efficiency calculated using observed oxygen transfer and *modeled blower energy*; not directly metered full-scale electricity |
| **Adel et al. (2019)** | Laboratory mechanical surface aerator with nine-blade configuration | **2.60** | Maximum *standard* efficiency under optimized test conditions |
| **Thakre et al. (2009)** | Laboratory oxidation ditch, curved-blade rotor | **2.269** | Maximum in clean-water reaeration tests |
| **Roy et al. (2021)** | Laboratory propeller-diffused aeration system | **1.711** | Experimentally validated optimum *standard* efficiency |
| **Mahmud et al. (2024)** | Laboratory confined-tube/Venturi aeration, three parallel injectors | **0.542** | Maximum *standard* efficiency in the tested configuration |

### Publication links and methodological notes

1. **Baeten, J. E., et al. (2021).** “Potential of off-gas analyses for sequentially operated reactors demonstrated on full-scale aerobic granular sludge technology.” *Science of the Total Environment*, **787**, 147651. DOI: [10.1016/j.scitotenv.2021.147651](https://doi.org/10.1016/j.scitotenv.2021.147651). The authors report **1.53 kg O₂/kWh on average under process conditions**, derived from full-scale off-gas analysis.

2. **Côté, P., et al. (2015).** “A New Membrane-Aerated Biofilm Reactor for Low Energy Wastewater Treatment: Pilot Results.” *Proceedings of the Water Environment Federation / WEFTEC*. DOI: [10.2175/193864715819540883](https://doi.org/10.2175/193864715819540883). [Technical paper PDF](https://www.mcilvainecompany.com/Decision_Tree/subscriber/Tree/DescriptionTextLinks/Hybrid_membrane.pdf). The authors describe **3–4 kg O₂/kWh after biofilm establishment**, with values **up to 5 kg O₂/kWh** under changed test conditions; electrical demand was calculated for a modeled blower, while **greater than 6 kg O₂/kWh** was a projected optimization, not a demonstrated result.

3. **Adel, M., et al. (2019).** “A comparative study of impeller aerators configurations.” *Alexandria Engineering Journal*, **58**(4), 1431–1438. DOI: [10.1016/j.aej.2019.11.014](https://doi.org/10.1016/j.aej.2019.11.014). Maximum SAE: **2.60 kg O₂/kWh** for a nine-blade aerator at 500 rpm.

4. **Thakre, S. B., Bhuyar, L. B., & Deshmukh, S. J. (2009).** “Oxidation ditch process using curved blade rotor as aerator.” *International Journal of Environmental Science and Technology*, **6**(1), 113–122. DOI: [10.1007/BF03326065](https://doi.org/10.1007/BF03326065). [Article PDF](https://applications.emro.who.int/imemrf/ijes_2009_6_1_113.pdf). Laboratory maximum: **2.269 kg O₂/kWh** at 48 rpm, 5.5 cm immersion depth, and a 47° blade-tip angle.

5. **Roy, S. M., et al. (2021).** “Prediction of standard aeration efficiency of a propeller diffused aeration system using response surface methodology and an artificial neural network.” *Water Supply*, **21**(8), 4534–4547. DOI: [10.2166/ws.2021.199](https://doi.org/10.2166/ws.2021.199). Experimentally validated optimum SAE: **1.711 kg O₂/kWh**, compared with a model prediction of 2.020 kg O₂/kWh.

6. **Mahmud, R., Carpenter, J., & MacPhee, D. W. (2024).** “Performance Assessment of Confined Tube Aerators in Parallel Configuration.” *Green Energy and Environmental Technology*. DOI: [10.5772/geet.20240041](https://doi.org/10.5772/geet.20240041). Maximum SAE: **0.542 kg O₂/kWh** using three parallel injectors; this was a laboratory configuration rather than a full-scale operating wastewater plant.

## 3. Practical interpretation for BNR systems

For an initial **field-condition comparison**, the most useful sources are:

- **1.53 kg O₂/kWh:** a published **full-scale process-condition** result for an aerobic granular sludge reactor (Baeten et al., 2021).
- **2.5–3.5 kg O₂/kWh:** WERF/WEF **field-condition benchmark for fine-pore diffusion**; use only if that technology is comparable to the system under study.
- **1.0–2.0 kg O₂/kWh:** WERF/WEF **field-condition benchmark for coarse-bubble diffusion**.

An illustrative **1.5–3.5 kg O₂/kWh** interval spans the full-scale Baeten result and the WERF fine-pore field benchmark. **It is not a universal measured range for all BNR reactors**, nor a direct measurement of any specific experimental tank. In particular, a continuously or intermittently aerated, oxygen-limited pilot with unusual airflows may have substantially different efficiency.

**Important:** Low dissolved oxygen (DO), even near zero, does **not** show that oxygen transfer is negligible: simultaneous microbial oxygen consumption can keep DO low during substantial oxygen transfer. Air flow alone also cannot establish AE because it describes oxygen *supplied* in air, not oxygen actually *transferred* into wastewater. A defensible process AE estimate needs either off-gas measurements, an appropriate oxygen mass balance, or calibrated mass-transfer measurements, plus aeration-system electricity use over a matched interval.

## Recommended citation strategy

For a manuscript or experimental report, cite **Baeten et al. (2021)** for a quantitative full-scale process observation and **Crawford and Sandino (2010), Table 2-8** for technology-specific field benchmarks. Use the **Ontario guidance manual, Table 5.1** for general clean-water standard reference ranges. Keep the laboratory optimization studies in a separate comparison group and label them as *standard*, not *field*, efficiency.
