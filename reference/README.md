# Reference material for the field aeration-efficiency range (pass 6, G5)

This folder documents the sources behind `scenarios.sae_field_kgO2_per_kWh = [0.7, 1.6, 3.5]` in `config/experiment.toml`. Publisher PDFs are not committed (`reference/*.pdf` is in `.gitignore`). Each is identified below by citation and SHA-256 so that a copy can be checked.

## 1. Primary source used: Crawford and Sandino (2010), Table 2-8

- **Citation.** Crawford, G. and Sandino, J. (2010). *Energy Efficiency in Wastewater Treatment in North America: A Compendium of Best Practices and Case Studies of Novel Approaches.* Water Environment Research Foundation (WERF), report OWSO4R07e. CH2M HILL Canada Limited. © WERF.
- **Where to obtain it.** The U.S. DOE Better Buildings Solution Center copy: <https://betterbuildingssolutioncenter.energy.gov/sites/default/files/tools/EnergyEfficiency_WWTreatment_BMPs.pdf> (file `EnergyEfficiency_WWTreatment_BMPs.pdf`, 122 pages).
- **SHA-256 of that copy:** `8b7f10c48995dfaf9a5dd91650cd43519780ff71ebadc317dc24eb15bd710def`. Recorded by the independent reviewer in `docs/rereview_2026-11-05_pass6.md`.
- **Location in the report.** Table 2-8 "Aeration Efficiency", p. 2-20, section 2.4.6.1. It is cited in the text as "Table 2-8 (WEF, 2009)", with the source line "Based on WEF MOP FD-13 Table 5.1 with additional content".

### Transcription of Table 2-8, field column (kg O2/kWh)

| Aerator type | Field |
| --- | --- |
| fine-pore¹ | 2.5–3.5 |
| coarse bubble¹ | 1.0–2.0 |
| surface centrifugal, low speed | 0.7–1.4 |
| surface centrifugal, draft tube | 0.7–1.3 |
| surface axial, high speed | 0.7–1.2 |
| downdraft open turbine | 0.6–1.2 |
| downdraft closed turbine | 0.7–1.3 |
| submerged turbine, sparger¹ | 0.7–1.1 |
| submerged impeller | 0.7–1.1 |
| surface brush and blade | 0.5–1.1 |

¹ Footnote in the source: "1 – Includes blower power requirements". The source spells "coarse" as "Course".

- **Field conditions.** The table states no operating DO or alpha for the field values.
- **Standard column.** The table also has a standard (clean-water) column, which is not used and not transcribed here. Only the fine-pore value, 5.0–6.5, is quoted in the report.
- **Verification.** The ten field rows were transcribed in `docs/progress.md` ("Remediation pass 6, G5") and checked row by row against the PDF by the reviewer.

### How the range follows from it (registered rule, G5)

- **Rows used:** the three rows marked "Includes blower power requirements", that is, fine-pore, coarse bubble, and submerged turbine with sparger.
- **Range:** low = 0.7, high = 3.5, central = √(0.7 × 3.5) = 1.565, rounded to 1.6.
- **Test:** `tests/test_pass6_fixes.py::test_field_ae_range_is_the_registered_rule_on_table_2_8` holds those three rows as literals and recomputes the rule.

## 2. Secondary summary: `wastewater_aeration_efficiency_literature.md`

A literature summary dated 2026-10-27, supplied by the user. It is context, not a source; only Table 2-8 was taken from it, after checking against the primary PDF. Known limits, from the G5 verification in `docs/progress.md`:

- Its "Table 2-8" reproduces 5 of the 10 rows (§1 above has all ten field values).
- Baeten et al. (2021), 1.53 kg O2/kWh: **not verified.** The full text is closed access and the abstract gives no kg O2/kWh value. Not used.
- Côté (2015), Adel (2019), Thakre (2009), Roy (2021), Mahmud (2024) and the Ontario Table 5.1 are laboratory, pilot or clean-water (standard) values. Not used as field values.
- Its §3 "illustrative 1.5–3.5 kg O2/kWh" interval is **not** the project's range.

## 3. Publisher PDFs read during G5 and not used (not committed)

| Local file | Citation | SHA-256 | Outcome |
| --- | --- | --- | --- |
| `1-s2.0-S004313541000638X-main.pdf` | Gresch, M., Armbruster, M., Braun, D. and Gujer, W. (2011). *Water Research* 45, 810–818. doi:10.1016/j.watres.2010.09.009 | `1657ae1c9837c151a589fc55e53cf9026d47239d63b7a2f824082b9aa6a87796` | No aeration efficiency in kg O2/kWh. Only a power density of 10–30 W/m³, quoted from Metcalf & Eddy (2003). Not used. |
| `1-s2.0-S1944398624019593-main.pdf` | Azis@Aziz, Ramli and Hamid (2019). *Desalination and Water Treatment* 156, 38–45. doi:10.5004/dwt.2019.24134 | `98373d42a9e5b4f3822065cac2e5b7cf60f8cfacba29ae4db89ad2500703331c` | Blower power and kWh only, no oxygen transferred. Not used. |
