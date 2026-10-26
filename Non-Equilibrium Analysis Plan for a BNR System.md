# Non-Equilibrium Analysis Plan for a BNR System

Oct 18, 2026 · @Xingyu Wang

## Scope, verified data status and hypotheses

This plan establishes how far a 100 mL batch of wastewater departs from thermodynamic equilibrium when it is driven by four known NH₄Cl pulses, and how that departure relates to nitrogen conversion and energy use. It also uses the known pulses to validate the four sensors, because the first comparison of predicted and recorded responses shows that three of the four sensors did not respond as the dose requires.

**Experimental system.**

1. **Batch closed to solutes.** The system is a 100 mL batch of wastewater in a vessel open to air and mixed by a magnetic stirrer, without forced aeration. Over the 50 days, the only material added was 25 mg of NH₄Cl at days 8.53, 10.57, 24.59 and 32.39 (times on the logger clock), and nothing else was added or removed. Oxygen enters, and CO₂ and NH₃ can leave, only through the stirred air–water surface.
2. **Known pulses.** Each pulse adds 4.67 mmol/L of NH₄Cl, which is 84.3 mg NH₄⁺/L (65.5 mg N/L) and 166 mg Cl⁻/L. If fully nitrified, each pulse releases 163 J of Gibbs energy in the 100 mL (1.63 kJ/L at standard conditions). The four pulses add 262 mg N/L and 18.7 mmol/L of chloride in total.
3. **Chloride as a conservative tracer.** Chloride is neither produced nor consumed, so it must rise as a staircase of four equal steps and stay there. This gives a built-in test of the conductivity and ORP sensors.
4. **Sensors.** Conductivity, NH₄⁺, pH and ORP probes sit in the same liquid at slightly different positions. ORP is measured with a Pt working electrode against a bare Ag quasi-reference. Conductivity is EC (mS/cm) = −2.89855 × I (A) + 0.520870.
5. **Common window.** The four records overlap from day 0 to 47.92.

**Verified file properties.**

| File | Valid rows | Span (d) | Native sampling | Timestamp resolution | Value range | Open points |
| --- | --- | --- | --- | --- | --- | --- |
| Cond.csv | 71,759 of 71,760 (last row empty) | 0.0007 to 49.70 | 60 s, uniform, no gaps | Full precision | 0.693 to 1.582 mS/cm after conversion | No response to the NH₄Cl pulses (Phase 0, S2) |
| NH4.csv | 414,785 | 0 to 47.92 | 20 s to day 41.05, then about 3 s, then about 2 s from day 45 | Rounded to 0.01 d (14.4 min); 43 to 450 rows per stamp | 9.5 to 75.1 mg NH₄⁺/L (7.4 to 58.3 mg N/L) | Responds at 17 to 37% of the dose, with delays of up to 7.7 h; four abrupt jumps unrelated to any addition (S1, S2) |
| ORP.csv | 459,884 | 0 to 53.23 | 10 s | Rounded to 0.01 d; 86 to 87 rows per stamp | +78 to +319 mV vs bare Ag | No response to the chloride added with each pulse (S3) |
| pH.csv | 414,757 | 0 to 48.00 | 10 s; two samples missing at day 4.40 | Full precision | 6.47 to 8.02 | Quantized to 0.01 pH |

All time columns are monotonic. In the NH₄⁺ and ORP files, rows sharing a stamp are in chronological order, so their true times are reconstructed by spacing them evenly within each 14.4-minute stamp, under the stated assumption of uniform sampling.

**Hypotheses.**

1. **H0, sensor validity.** Each sensor reproduces the response that the known pulses require. This hypothesis is currently rejected for conductivity and ORP and only partly supported for NH₄⁺, so it is the gate for every quantitative conclusion.
2. **H1, chemical disequilibrium.** The batch remains far from equilibrium throughout, because NH₄⁺ coexists with oxygen supplied through the open, stirred surface, and because it is consumed at a nonzero net rate. Both criteria are independent of the ORP reference.
3. **H2, broken detailed balance.** After sensor artifacts are excluded, the trajectories are statistically time-irreversible, with irreversibility concentrated in the relaxation after each pulse.
4. **H3, conversion is kinetically, not thermodynamically, limited.** Each pulse raises the nitrification affinity by only 2.4 to 4.2 kJ/mol against about 349 kJ/mol N. The rate response therefore reflects substrate kinetics, inhibition and biomass state, not driving force.
5. **H4, Onsager regression fails.** The relaxation after the pulses does not follow the regression of spontaneous fluctuations, which is the signature of a system far from equilibrium.
6. **H5, conservation closes.** Once the sensors are validated, nitrogen, chloride, charge and alkalinity balances close across the four pulses.
7. **H6, efficiency and energy.** Non-equilibrium descriptors computed for each relaxation segment predict the conversion of each pulse and, once surface oxygen transfer and stirrer power are measured, the oxygen and energy intensity.

## Phase 0: sensor validation with the known NH₄Cl pulses

The four NH₄Cl pulses are standard additions of known size into a known, stirred volume, so they test every sensor directly. Only pH responded as predicted. The NH₄⁺ electrode registered 17 to 37% of the dose, its sensitivity declined from pulse to pulse, and for pulses 1 and 2 most of its response arrived about 7.5 h late. Magnetic stirring mixes 100 mL within minutes, so that delay is electrode or transmitter behavior. Conductivity showed neither the step at each pulse nor the cumulative rise that the added chloride requires, and ORP showed no chloride response at the bare Ag reference.

These results are P0 findings for any quantitative use of the conductivity record, the NH₄⁺ magnitudes and absolute ORP. They supersede the earlier interpretation of the NH₄⁺ steps as additions and removals, and they withdraw the chloride-based redox-gap result.

**Predicted and observed responses to each pulse.** NH₄⁺ and conductivity changes are taken as the median 0.4 to 0.6 d after the pulse minus the median 0.1 d before. ORP and pH changes are taken 1 to 3 h after minus 3 to 1 h before. NH₄⁺ is reported by the instrument in mg NH₄⁺/L. The ORP prediction assumes a chloride-controlled Ag surface and a baseline chloride of 1 to 13.5 mM.

| Pulse (day) | NH₄⁺ predicted / observed (mg NH₄⁺/L) | NH₄⁺ response | EC predicted / observed (mS/cm) | ORP predicted / observed (mV) | pH observed |
| --- | --- | --- | --- | --- | --- |
| 8.53 | +84.3 / +31.1, almost all in a jump 7.7 h later | 37% | +0.65 / +0.005 | +7.6 to +45 / +0.4 | −0.01 |
| 10.57 | +84.3 / +24.4 (+8.1 at the pulse, +15.9 at 7.4 h) | 29% | +0.65 / −0.062, during an ongoing decline | +7.6 to +45 / 0.0 | +0.02 |
| 24.59 | +84.3 / +22.1, rising over about 2 h | 26% | +0.65 / −0.007 | +7.6 to +45 / +0.1 | 0.00 |
| 32.39 | +84.3 / +14.7 | 17% | +0.65 / +0.016 | +7.6 to +45 / +0.2 | +0.04 |

The predicted pH change is at most a few hundredths of a unit, because only about 1% of the added NH₄⁺ deprotonates at pH 7.3. The observed pH changes agree with this prediction.

**Cumulative balances.** The cumulative balances reach the same verdict as the individual pulses:

- **Conductivity.** Starting from 1.07 mS/cm, the end-of-record conductivity should be 2.16 mS/cm even if all added NH₄⁺ were nitrified (chloride and nitrate replacing two bicarbonates), or 3.65 mS/cm with no nitrification. Evaporation from the open vessel would raise it further. The record ends at 1.07 mS/cm.
- **Nitrogen.** With 36.6 mg N/L at the start (47.1 mg NH₄⁺/L), 262 mg N/L added and 33.5 mg N/L at the end, 265 mg N/L must have left the NH₄⁺ pool, about 5.6 mg N/L/d. The NH₄⁺ record itself contains only 39 mg N/L of gradual decline and 53 mg N/L of abrupt drops.

**Findings.**

| ID | Severity | Finding | Required action |
| --- | --- | --- | --- |
| S1 | P0 for NH₄⁺ magnitudes | **NH₄⁺ electrode under-response.** In its reporting unit of mg NH₄⁺/L, the electrode registered 37, 29, 26 and 17% of the four identical pulses of 84.3 mg NH₄⁺/L. The steady decline points to a falling electrode slope, as in membrane aging, and to calibration error, rather than to mixing. Any gain error scales all consumption rates by the same factor, but leaves their ratios and the logarithmic affinity changes nearly unaffected. | Review the calibration record. Recalibrate the electrode in standards matched to the sample ionic strength. Measure NH₄⁺-N in the laboratory on stored or new samples, and perform a multi-point standard addition (Gran plot) in the vessel. |
| S2 | P0 for NH₄⁺ levels | **Abrupt NH₄⁺ jumps.** Four jumps are unrelated to any addition, at days 2.59 (−12), 17.84 (−27), 37.19 (−29) and 44.63 (+6 mg NH₄⁺/L). Two more, at days 8.85 (+33) and 10.88 (+16), are delayed responses to pulses 1 and 2 (S5). Nothing was removed, and the three decreases are more than 200 times faster than the consumption measured in the same liquid, so they are not concentration changes. At a Nernstian slope, the 29 mg NH₄⁺/L drop at day 37.19 corresponds to a potential shift of only about 14 mV, typical of a reference junction or membrane shift. | Treat each of the four unrelated jumps as an electrode offset, estimated jointly with the kinetics (Phase 4). Check the maintenance and handling record for these times. |
| S3 | P0 for conductivity | **Conductivity does not register the pulses.** Each pulse should raise conductivity by 0.65 mS/cm, and the chloride should keep it elevated. The observed changes are between −0.06 and +0.02 mS/cm, and the record shows no cumulative rise. No mixing delay can explain a missing cumulative rise over 48 days in 100 mL, so the record does not represent bulk solution conductivity. Its other features (steps at days 9.41, 15.67 and 35.93, the 29 late pulses, and the overall range of 0.69 to 1.58 mS/cm) therefore cannot be attributed to ionic composition until the probe is validated. | Check that the cell is fully immersed and free of trapped air. Measure the vessel with a calibrated handheld meter. Add a known KCl standard to test the response, and verify the current-to-conductivity equation against standards spanning 0.5 to 4 mS/cm. |
| S4 | P0 for absolute E\_H | **No chloride response at the bare Ag reference.** If the Ag surface were chloride-controlled, each pulse would raise ORP by 7.6 to 45 mV, and the four pulses together by 22 to 77 mV. Observed changes are 0.0 to 0.4 mV within 3 h. Either the Ag surface is not chloride-controlled, or the chloride did not reach it within hours. An oxide-controlled surface (Ag₂O/Ag, about 0.74 V vs SHE at pH 7.3) or a sulfide-controlled one would give an E\_ref hundreds of mV from the chloride case, so E\_ref, and with it absolute E\_H, is undetermined. The earlier redox-gap table and the chloride inference about the late conductivity pulses are withdrawn. | Spike a known NaCl amount under stirring and record ORP. Check the Pt and bare Ag pair in ZoBell's solution. Inspect the Ag surface, and preferably replace the bare Ag with a double-junction Ag/AgCl reference. |
| S5 | P1 | **Delayed electrode response.** The pulse times are confirmed on the logger clock. The electrode registered pulse 1 almost entirely 7.7 h late, and pulse 2 partly at once (+8.1) and mostly 7.4 h late (+15.9 mg NH₄⁺/L), while pulse 3 took about 2 h to register and pulse 4 registered at once. With magnetic stirring in 100 mL, mixing takes minutes, so these delays reflect a slow or intermittently updating electrode or transmitter, not the liquid. | Measure the electrode's step response in a stirred standard, and check the transmitter's averaging, hold and update settings. |
| D1 | Closed | **Conductivity unit.** EC (mS/cm) = −2.89855 × I (A) + 0.520870 gives 0.693 to 1.582 mS/cm. | Validate the equation against standards (S3). |
| D3 | P1 | **Held NH₄⁺ output.** The high-frequency residual is 0.66 mg NH₄⁺/L from day 2.59 to 8.82 and zero from day 10.85 to 24.60 and from day 32.36 to 37.19, which indicates held or strongly damped output. | Confirm the output filter settings; exclude held segments from fluctuation estimators. |
| D4 | P2 | **NH₄⁺ logging rate.** The logging rate rises from one sample per 20 s to about one per 3 s at day 41.05, and to about one per 2 s near day 45. | Use time-weighted bin means. |
| D7 | P1 | **Time base.** The four-sensor overlap is day 0 to 47.92, and pH alignment with the other loggers is not yet verified. | Confirm logger start times. |
| D8 | P1 | **Late conductivity pulses.** From day 42, 29 conductivity pulses occur (median amplitude 0.090 mS/cm, 24-minute decay), with no response in any other probe. Since nothing was added, they cannot be ionic additions, and they are treated as conductivity-probe events, for example bubbles in the cell. | Inspect the probe; exclude the pulses from process analyses. |
| D9 | P2 | **pH quantization and dips.** pH is quantized to 0.01, two samples are missing, and 111 dips exceed 0.08 below the 1-hour median. | Use tie-aware estimators; investigate the dips. |
| D10 | P2 | **Start-up transient.** pH and ORP are unsettled during day 0 to 1. | Exclude day 0 to 1 from steady-state statistics. |

**Dataset construction.**

1. **Reconstruct times.** Reconstruct NH₄⁺ and ORP times within stamps, and convert conductivity with the stated equation, flagged as unvalidated (S3).
2. **Build a common grid.** Build a 0.01-day grid of time-weighted means over day 0 to 47.92, with a gap mask and per-value QC flags. Keep the native series for pH dynamics.
3. **Build the event table.** Build an event table of the four known pulses (time, 4.67 mmol/L NH₄Cl) and of the sensor artifacts (S2 jumps, D8 pulses), with offsets carried as parameters.
4. **Segment the record.** Segment the record at the four pulses, so that every relaxation segment starts at a known input, and correct the S2 offsets within segments.
5. **Generate outputs.** Produce the flagged dataset, the data dictionary, the event table and a validation report with the predicted and observed responses above.

## Theoretical framework and computable quantities

Non-equilibrium is defined operationally by nonzero thermodynamic forces and by positive entropy production. Each reaction j and each gas–liquid transfer contributes an affinity, A\_j, that is computed from in situ activities. The analysis therefore rests on four connected layers:

1. **Equilibrium chemistry** supplies the reference state.
2. **Bioenergetics** relates affinities to microbial rates and yields.
3. **Irreversible thermodynamics** converts forces and fluxes into dissipation.
4. **Stochastic thermodynamics** bounds dissipation directly from the statistics of measured trajectories.

**Affinity and entropy production.**

```latex
A_j = -\Delta_r G_j = -\left(\Delta_r G_j^{\circ\prime} + RT \ln Q_j\right), \qquad \sigma = \frac{1}{T}\sum_j r_j A_j + \frac{1}{T}\sum_k J_k\, RT \ln\frac{C_k^{*}}{C_k} \;\ge\; 0
```

Here r\_j is the volumetric rate of reaction j. J\_k = k\_La,k (C\*\_k − C\_k) is the gas transfer flux of species k (O₂ absorption, CO₂ and N₂O stripping) driven by its departure from saturation C\*\_k. Equilibrium corresponds to all A\_j = 0 and all C\_k = C\*\_k, so every measurable departure from these conditions is a quantitative signature of non-equilibrium.

**Disequilibrium criteria independent of ORP.** Two criteria establish non-equilibrium without any electrode potential:

- **Net rate criterion.** The De Donder relation requires each reaction's affinity and net rate to share sign, A\_j r\_j ≥ 0, with A\_j = 0 only where r\_j = 0. A sustained nonzero net consumption of NH₄⁺ therefore proves a positive nitrification affinity, provided laboratory analysis confirms that the decline is real and not electrode drift (S1, S2).
- **Coexistence criterion.** NH₄⁺ coexisting with dissolved O₂ is a state of very large affinity. With NH₄⁺ at 40 mg N/L, NO₃⁻ at 100 mg N/L and DO at only 0.1 mg/L, the in situ affinity of complete nitrification at pH 7 and 25 °C is still about 317 kJ/mol. Even a hundredfold error in any concentration changes this by at most RT ln 100 = 11 kJ/mol. This criterion requires only that oxygen is present, which here is ensured by the open, stirred air–water surface.

**Redox couples and apparent electron activity.** For a reduction half-reaction with n electrons and m protons, the Nernst potential and the apparent electron activity derived from ORP (converted to the hydrogen scale, E\_H) are:

```latex
E_{H,i} = E_i^{\circ} - \frac{2.303RT}{n_i F}\log\frac{\prod a_{\mathrm{red}}^{\nu}}{\prod a_{\mathrm{ox}}^{\nu}} - \frac{2.303RT\, m_i}{n_i F}\,\mathrm{pH}, \qquad p\varepsilon_{\mathrm{app}} = \frac{F\,E_H}{2.303RT}
```

The Pt electrode potential measured against the bare Ag quasi-reference converts to the hydrogen scale as E\_H = E\_meas + E\_ref(t). In chloride-containing water, the Ag surface acquires an AgCl film and behaves as an AgCl/Ag electrode of the second kind, whose potential follows the chloride activity of the liquid. Sulfide would instead form Ag₂S, whose standard potential is about −0.69 V, and shift E\_ref strongly negative. An oxidized surface would follow the Ag₂O/Ag couple, about 0.74 V vs SHE at pH 7.3, and a clean surface would take an ill-defined mixed potential. The absent chloride response to the NH₄Cl pulses (S4) shows that the chloride case cannot currently be assumed. Every E\_H therefore requires a chloride record, or an in situ calibration, and a statement of the assumed surface state:

```latex
E_{\mathrm{ref}} = E^{\circ}_{\mathrm{AgCl/Ag}} - \frac{RT}{F}\ln a_{\mathrm{Cl^-}}, \qquad E^{\circ}_{\mathrm{AgCl/Ag}} = 0.2223\ \mathrm{V}\ (25\,^{\circ}\mathrm{C})
```

Once the conductivity probe is validated, Kohlrausch additivity bounds the chloride concentration from above, since chloride cannot carry more current than the whole solution. A chloride standard addition measures it directly: a known spike Δc changes the ORP signal by ΔE = 0.05916 log\[(c₀ + Δc)/c₀\], which can be solved for the baseline c₀. The same measurement verifies the Nernstian slope of the Ag surface.

The relevant formal potentials at pH 7 and 25 °C, derived from the electron-equivalent free energies tabulated by Rittmann and McCarty, are:

| Couple (reduction) | n (e⁻) | E°′ at pH 7 (V vs SHE) |
| --- | --- | --- |
| O₂ / H₂O | 4 | +0.816 |
| NO₃⁻ / N₂ | 5 | +0.748 |
| NO₃⁻ / NO₂⁻ | 2 | +0.432 |
| NO₃⁻ / NH₄⁺ | 8 | +0.364 |
| NO₂⁻ / NH₄⁺ | 6 | +0.341 |

The corresponding standard free energies are −274.7 kJ/mol for nitritation (NH₄⁺ to NO₂⁻ with O₂) and −74.1 kJ/mol for nitratation (NO₂⁻ to NO₃⁻), giving −348.8 kJ/mol N for complete nitrification. These values must be corrected to in situ temperature and activities before use.

**Speciation from the available sensors.** Conductivity, pH and NH₄⁺ suffice to compute four quantities:

- ionic strength, from Russell's approximation I ≈ 1.6 × 10⁻⁵ × EC (µS/cm);
- activity coefficients, from the Davies equation;
- the ammonium dissociation constant, pK\_a = 0.09018 + 2729.92/T (Emerson et al.);
- free ammonia (Anthonisen et al.), which with free nitrous acid governs the selective inhibition of nitrite oxidizers. Free nitrous acid itself requires a nitrite measurement.

```latex
\log\gamma_i = -A z_i^2\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3I\right), \qquad \mathrm{FA} = \frac{17}{14}\,\frac{\mathrm{TAN}\cdot 10^{\mathrm{pH}}}{\exp(6344/T) + 10^{\mathrm{pH}}}
```

**Proton and ionic balances.** The proton production rate follows from r\_H = −β (dpH/dt), where the buffer capacity β requires alkalinity or a carbonate system model. Kohlrausch's additivity, κ ≈ Σ λ\_i |z\_i| c\_i, predicts the conductivity signature of nitrification:

```latex
\mathrm{NH_4^+} + 2\,\mathrm{HCO_3^-} + 2\,\mathrm{O_2} \rightarrow \mathrm{NO_3^-} + 2\,\mathrm{CO_2} + 3\,\mathrm{H_2O}, \qquad \Delta\Lambda = \lambda_{\mathrm{NO_3^-}} - \lambda_{\mathrm{NH_4^+}} - 2\lambda_{\mathrm{HCO_3^-}} \approx -91\ \mathrm{S\,cm^2\,mol^{-1}}
```

At infinite dilution this corresponds to a decrease of about 6.5 µS/cm per mg/L of NH₄⁺-N nitrified. Departures from this ratio expose concurrent ionic processes.

**Kinetic and thermodynamic coupling.** The rate of a microbial reaction approaching equilibrium is reduced by the thermodynamic factor of Jin and Bethke:

```latex
r = r_{\mathrm{kin}}\,F_T, \qquad F_T = 1 - \exp\!\left(-\frac{-\Delta_r G - m\,\Delta G_{\mathrm{ATP}}}{\chi RT}\right)
```

Here m is the number of ATP synthesized per reaction and χ is the average stoichiometric number. This factor links the kinetic ASM framework to the thermodynamic driving forces calculated above.

**Stochastic thermodynamics.** For a stationary process observed through a coarse-grained variable x, the Kullback–Leibler divergence between the probabilities of forward and time-reversed trajectory segments gives a lower bound on the entropy production rate (Kawai, Parrondo and Van den Broeck; Roldán and Parrondo):

```latex
\dot S \;\ge\; k_B \lim_{n\to\infty}\frac{1}{n\,\tau}\, D_{KL}\!\left[P(x_1,\dots,x_n)\,\Vert\, P(x_n,\dots,x_1)\right]
```

The bound is rigorous but loose for a macroscopic reactor, because four sensors observe only a tiny fraction of the dissipative degrees of freedom. It is therefore used as a relative, phase-resolved index of irreversibility, never as an estimate of the physical entropy production.

**Pulsed closed batch.** In the 100 mL batch, NH₄⁺-N concentration C is driven only by the four pulses of known size, ΔC = 65.5 mg N/L at times t\_k, and by the net biological consumption rate r. Chloride C\_Cl is conservative and therefore rises as a staircase:

```latex
\frac{dC}{dt} = \sum_{k=1}^{4} \Delta C\,\delta(t-t_k) \;-\; r(C,\mathrm{pH},\mathrm{FA},T,X), \qquad C_{\mathrm{Cl}}(t) = C_{\mathrm{Cl},0} + \sum_{t_k<t} 4.67\ \mathrm{mmol\,L^{-1}}
```

Here X is the active nitrifier biomass. Without further input, a batch closed to solutes relaxes toward equilibrium, and only gas exchange through the stirred open surface keeps it driven. Oxygen enters at J = k\_La(C\* − C), with k\_La set by the stirring, while CO₂, NH₃, N₂ and N₂O can leave. Each pulse injects a known Gibbs energy of 163 J into the 100 mL, available for nitrification if oxygen is supplied. The pulses are therefore ideal for measuring how a driven biological system dissipates a known free-energy input.

Each pulse raises the nitrification affinity by ΔA = RT ln\[(C + ΔC)/C\]. Using the electrode baselines before each pulse, converted to mg N/L, this gives 2.4 to 4.2 kJ/mol at 25 °C, about 1% of the 349 kJ/mol standard driving force. Because the dependence is logarithmic, even a fourfold calibration error in C moves ΔA by at most RT ln 4 = 3.4 kJ/mol. Any change in consumption rate after a pulse is therefore governed by substrate kinetics, free ammonia inhibition and biomass state, for instance through a Haldane rate law, rather than by thermodynamic driving force. This is hypothesis H3.

**Onsager regression and fluctuation–dissipation.** Close to equilibrium, the spontaneous fluctuations of a variable regress on average exactly like a small imposed perturbation. Equivalently, the normalized autocorrelation of spontaneous fluctuations equals the normalized relaxation function after a small perturbation:

```latex
\frac{\langle \delta x(t+\tau)\,\delta x(t)\rangle}{\langle \delta x^2\rangle} \;=\; \frac{\Delta x(\tau)}{\Delta x(0)} \quad \text{(near equilibrium)}
```

The pulses supply the right-hand side of this equality, while the quasi-steady stretches late in each segment supply the left-hand side. Each pulse is large relative to the background (ΔC/C between about 1.6 and 4.3 on the electrode baselines), so the early relaxation is expected to be nonlinear. The test is therefore applied to the late tail of each relaxation, where the remaining deviation is small enough for linear response. There, a systematic, significant violation of the equality is a model-free signature of a far-from-equilibrium state. The test requires approximate stationarity, correction for sensor noise, and exclusion of the S2 electrode jumps.

**What is computable now.** From the four sensors alone, the following can be computed:

- ionic strength and activity coefficients;
- free ammonia;
- apparent pε, together with its offset from the formal potentials;
- proton and ionic stoichiometric ratios;
- all dynamical and information-theoretic signatures.

The complete affinities of nitritation, nitratation, denitrification and oxygen transfer, together with σ in physical units, additionally require temperature, DO, NO₂⁻, NO₃⁻ and alkalinity. These can be supplied either by measurement or by state estimates from the models in Phase 4.

## Phase 1: pulse segmentation and relaxation characterization

The record divides into a baseline period and four relaxation segments, each starting at a known NH₄Cl pulse: day 1 to 8.53, 8.53 to 10.57, 10.57 to 24.59, 24.59 to 32.39, and 32.39 to 47.92. Within these segments, the four unrelated S2 jumps are treated as electrode offsets, and the delayed responses to pulses 1 and 2 as electrode dynamics (S5), not as process events.

Lomb–Scargle analysis shows no dominant period (maximum normalized power 0.053), apart from a weak component near 1.0 d in pH. ORP never falls below +78 mV against the bare Ag, so no persistent anoxic phase is visible, although the absolute redox level is undetermined (S4).

**Consumption between electrode events.** Net NH₄⁺ trends are computed with Theil–Sen slopes on sub-segments bounded by pulses and S2 jumps, excluding the first 0.1 d after each boundary. They are converted from the electrode's mg NH₄⁺/L to mg N/L (factor 0.7765). If the under-response (S1) is a gain error, true rates are larger by the inverse of the response fraction, about 2.7 to 5.9 times. Free ammonia is computed with the Anthonisen expression at two assumed temperatures, pending the temperature record.

| Sub-segment (d) | Relaxation after | Net slope (mg N/L/d, uncorrected for gain) | Median pH | Median FA, 20 °C / 25 °C (mg NH₃/L) |
| --- | --- | --- | --- | --- |
| 0.3 to 2.59 | Baseline | −1.47 | 6.92 | 0.14 / 0.20 |
| 2.59 to 8.82 | Baseline, then pulse 1 at 8.53 | −2.66 | 7.15 | 0.15 / 0.21 |
| 8.82 to 10.85 | Pulse 1, contains pulse 2 | Not interpretable | 7.01 | 0.20 / 0.29 |
| 10.85 to 17.84 | Pulse 2 | −1.51 | 7.18 | 0.36 / 0.53 |
| 17.84 to 24.60 | Pulse 2 | −0.44 | 7.81 | 0.82 / 1.16 |
| 24.60 to 32.36 | Pulse 3 | −0.37 | 7.55 | 0.71 / 1.02 |
| 32.36 to 37.19 | Pulse 4 | −0.24 | 7.28 | 0.47 / 0.68 |
| 37.19 to 44.63 | Pulse 4 | −0.22 | 7.37 | 0.33 / 0.47 |
| 44.63 to 47.92 | Pulse 4 | +0.09 | 7.37 | 0.37 / 0.54 |

Two features stand out:

1. **Slowing consumption.** Consumption slows from 1.5 to 2.7 mg N/L/d before day 18 to 0.2 to 0.4 mg N/L/d afterward. Part of this apparent decline comes from the falling electrode sensitivity (S1). After correction by the pulse response fractions, the rate after pulse 2 (1.51/0.29 ≈ 5.2 mg N/L/d) still exceeds that after pulse 4 (0.24/0.17 ≈ 1.4 mg N/L/d) by about fourfold, although both pulses were identical.
2. **Rising pH and free ammonia.** From day 17.84 to 32.36, pH rises to 7.5 to 7.8 and free ammonia reaches 0.7 to 1.2 mg NH₃/L.

Time, pH, free ammonia, biomass state and electrode gain are confounded in this record, and separating them is the central modeling task of Phase 4. Slopes in segments with held output (D3) carry wider uncertainty.

**Kinetic signature of each relaxation.** Each segment is fitted with zero-order, first-order, Monod and Haldane (free-ammonia-inhibited) models, with the S2 offsets and an electrode gain as nuisance parameters. The models are compared by AIC, and residual autocorrelation is checked.

**Other events.** None of these events coincides with an addition, so each is either an internal process change or a sensor or handling event:

| Day | Signature | Provisional interpretation |
| --- | --- | --- |
| 2.88 | pH +0.5, ORP −29 mV, conductivity +0.036 mS/cm together | Handling or maintenance event affecting several probes, to be checked against the lab record |
| 10.69 | ORP decays from 195 to about 100 mV over about 0.5 d, starting 0.12 d after pulse 2 | Only ORP feature that follows a pulse; candidate oxygen depletion as nitrification of the pulse raises oxygen demand |
| About 24.5 and 32.5 | pH drops of about 0.15, about 0.1 d before pulse 3 and 0.14 d after pulse 4 | Candidate acidification by nitrification, or a pH-logger clock offset (D7) |
| 28.89, 30.48 | ORP steps of +132 and −48 mV | Pt or Ag electrode event (S4) |
| 9.41, 15.67, 35.93 | Conductivity steps of +0.15, −0.18 and +0.19 mS/cm | Not interpretable until the conductivity probe is validated (S3) |

**Segment-level features.** For each relaxation segment, the following features form the matrix used in Phases 3 to 5:

- the mean, variance and drift of pH and relative ORP;
- the integrated autocorrelation time;
- the spectral exponent;
- the kinetic parameters above;
- the conversion fraction of the pulse;
- phase-averaged daily profiles, which separate any diurnal component from relaxation.

**Validation.** The electrode's step response is measured to deconvolve the delayed responses (S5), every artifact is checked against the handling record, and gain-corrected rates are compared with laboratory NH₄⁺-N.

## Phase 2: stoichiometric closure and thermodynamic disequilibrium

Disequilibrium is established here without relying on ORP. NH₄⁺ coexists with an oxidant at an affinity of about 300 kJ/mol or more, and the pulses change that affinity by only about 1%. The absolute redox gap derived from ORP is withdrawn until the bare Ag reference is calibrated (S4), but relative ORP dynamics remain usable.

**Affinity of the batch.** The vessel is open to air and stirred, so oxygen is supplied continuously through the surface. As long as dissolved O₂ stays above zero, the affinity of complete nitrification remains above about 300 kJ/mol for any plausible NH₄⁺, NO₃⁻ and DO levels (theoretical framework). Lowering DO from 0.1 to 0.01 mg/L reduces it by only about 11 kJ/mol, so the batch is far from equilibrium throughout. A single DO measurement confirms the premise, and the net rate criterion adds independent support once laboratory NH₄⁺-N confirms that the electrode trends are real.

**Affinity step at each pulse.** The step is computed as follows, using the electrode baseline before each pulse and the known dose of 65.5 mg N/L:

```latex
\Delta A_k = RT\ln\frac{C(t_k^{-})+\Delta C}{C(t_k^{-})}, \qquad \Delta C = 65.5\ \mathrm{mg\,N\,L^{-1}}
```

| Pulse (day) | Baseline before pulse (mg N/L, from electrode mg NH₄⁺/L) | ΔA at 25 °C (kJ/mol) |
| --- | --- | --- |
| 8.53 | 15.1 | 4.2 |
| 10.57 | 38.3 | 2.5 |
| 24.59 | 25.7 | 3.1 |
| 32.39 | 41.0 | 2.4 |

All four steps are about 1% of the 349 kJ/mol standard driving force, and they remain small under any plausible correction of the baseline. Identical pulses that produce different rate responses (Phase 1) therefore probe kinetics and biomass state, not thermodynamic driving force (H3).

**Redox analysis once the reference is calibrated.** Until then, ORP is analyzed only as a relative signal, for example the 95 mV decline that began 0.12 d after pulse 2. After calibration, the analysis follows Lindberg and Runnells:

1. **Convert ORP.** Convert ORP to E\_H with the calibrated E\_ref(t), propagating its uncertainty (10 mV corresponds to 0.17 pε units at 25 °C).
2. **Compute couple potentials.** Compute Nernst potentials for the nitrogen and oxygen couples in each relaxation segment, with NO₃⁻, NO₂⁻ and O₂ from samples or Phase 4 estimates.
3. **Form disequilibrium indices.** Compute D\_i = E\_H − E\_H,i and the spread D\_spread, and track their response to each pulse.
4. **Compute equilibrium-implied compositions.** Report the equilibrium-implied NO₃⁻/NH₄⁺ ratio with its uncertainty.

**Speciation.** Ionic strength cannot be taken from the unvalidated conductivity record (S3). It is estimated instead from a major-ion analysis or a validated conductivity measurement, raised by 4.67 mmol/L at each pulse. Over the plausible range of 0.01 to 0.04 M, Davies coefficients for monovalent ions fall from about 0.90 to 0.83. Free ammonia reaches 0.7 to 1.2 mg NH₃/L between days 17.84 and 32.36 at 20 to 25 °C, within the range Anthonisen et al. report for inhibition of nitrite oxidizers. It scales with the electrode gain (S1) and with temperature, both of which must be fixed.

**Conservation balances of the closed batch.** Three balances constrain the process independently of the electrode problems, and each becomes quantitative once the missing measurements exist:

- **Chloride and conductivity.** Chloride must rise by 4.67 mmol/L at each pulse and stay there. This predicts a conductivity staircase of 0.65 mS/cm per pulse, falling to a net 0.27 mS/cm per pulse after complete nitrification. Its absence is the S3 finding; once the probe is validated, the observed staircase gives the extent of nitrification directly.
- **Nitrogen.** By the end of the record, 266 mg N/L must have left the NH₄⁺ pool. Even after gain correction, the electrode trends account for only about half of this. Closure requires NO₃⁻, NO₂⁻ and total nitrogen, and should consider ammonia volatilization at pH 7.5 to 7.8, assimilation into biomass, and gaseous losses as N₂ or N₂O.
- **Protons and alkalinity.** Complete nitrification of one pulse consumes 9.35 mmol/L of alkalinity, 468 mg/L as CaCO₃, so all four pulses would consume 1,870 mg/L as CaCO₃. Without base addition, pH stayed between 6.5 and 8.0. Either the initial alkalinity was very high, or only part of the added nitrogen was nitrified, or denitrification returned part of the alkalinity. A single initial alkalinity measurement, combined with the pH record and a carbonate model, therefore bounds the cumulative nitrification independently of the NH₄⁺ electrode.

**Reaction affinities and dissipation.** In situ affinities for nitritation, nitratation, denitrification and oxygen reduction, together with σ\_chem(t) = Σ r\_j A\_j / T, require temperature, DO, NO₂⁻ and NO₃⁻ from measurement or Phase 4 estimates. Because consumption rates are small, chemical dissipation per unit volume is small. Most of the imposed disequilibrium may therefore be dissipated through aeration and gas transfer rather than through nitrogen conversion, which links directly to energy efficiency in Phase 5.

**Gas–liquid disequilibrium at the open surface.** The open, stirred surface makes gas exchange the only sustained driving force, so its disequilibria are central descriptors:

- **Oxygen deficit.** The deficit (C\* − C) sets the oxygen supply k\_La(C\* − C). C\* is about 8.3 mg/L at 25 °C, so the gain-corrected nitrification demand of about 6 to 33 mg O₂/L/d requires k\_La of at least 0.7 to 4 d⁻¹ even at zero DO. Surface reaeration alone may therefore limit nitrification whenever biomass respiration also draws on the oxygen supply. DO and k\_La, the latter measured by a reaeration test, decide this.
- **CO₂ supersaturation.** At pH 7.3, water in CO₂ equilibrium with air (about 420 µatm) can hold only about 0.13 mmol/L of bicarbonate alkalinity, or 0.40 mmol/L at pH 7.8. Any realistic alkalinity therefore implies CO₂ supersaturation sustained by respiration, and its degree, Ω\_CO₂ = \[HCO₃⁻\]\[H⁺\]/(K₁K\_H pCO₂,air), follows from the pH record once alkalinity is known.
- **Ammonia volatilization.** Free ammonia of 0.7 to 1.2 mg/L against an essentially ammonia-free atmosphere drives slow loss of NH₃ through the stirred surface. This loss belongs in the nitrogen balance.

Every index carries a Monte Carlo interval that propagates the pH, NH₄⁺ gain, alkalinity and temperature uncertainties.

## Phase 3: dynamical and information-theoretic non-equilibrium signatures

Phase 3 tests non-equilibrium from trajectory statistics, without assumed activities. Irreversibility estimators cannot by themselves distinguish process dynamics from instrument dynamics, such as electrode drift, jumps and slow responses. For that reason, Phase 3 relies first on pH, the one sensor that passed the pulse test, and treats NH₄⁺, ORP and conductivity results as provisional until their sensors are validated. Stationary estimators are applied within relaxation segments, and response estimators across the four pulses. The analysis variables, and the segments excluded from each, are:

- **pH:** native 10 s series, excluding day 0 to 1;
- **Conductivity:** excluded from process analyses until validated (S3), and analyzed only to characterize the probe's own behavior;
- **ORP and NH₄⁺:** 0.01-day grid, offset-corrected for the S2 jumps, excluding held NH₄⁺ segments (D3) and the ORP electrode events (days 28.89 and 30.48), and with ORP used only as a relative signal (S4).

**Time irreversibility.** Three estimators are computed within each stationary segment, where stationarity is verified by augmented Dickey–Fuller and KPSS tests:

- **Moment asymmetry.** The asymmetry ⟨(x(t+τ) − x(t))³⟩ as a function of lag τ.
- **Ordinal-pattern divergence.** The Kullback–Leibler divergence between ordinal-pattern frequencies and those of their time-reversed counterparts.
- **Visibility-graph divergence.** The in-degree and out-degree divergence of the horizontal visibility graph (Lacasa et al.).

Because 88% of consecutive pH readings are identical, ordinal patterns must be computed with explicit tie handling, as in modified permutation entropy, or on bin-averaged series. Otherwise ties bias every pattern statistic.

Significance is tested against iterated amplitude-adjusted Fourier transform surrogates, which are time-reversible by construction.

**Relaxation after the four pulses.** The four identical NH₄Cl pulses are the strongest material in the record, because they apply the same known input at four different system states. Three analyses use them:

1. **Response functions.** For each pulse, extract the response of NH₄⁺ (gain- and offset-corrected), pH and relative ORP. Fit single or multiple relaxation times, and test whether the relaxation is exponential, as in linear response, or substrate-limited and inhibited, as expected far from equilibrium.
2. **Onsager regression test.** Compare the normalized late-tail response of each pulse with the normalized autocorrelation of spontaneous fluctuations in the same segment. Report their ratio as a function of lag, with block-bootstrap intervals (H4).
3. **State dependence.** Because the input is identical, differences between the four responses measure how the system's state changed, through biomass activity, free ammonia, pH and accumulated products. A system near equilibrium with fixed parameters would respond identically to identical inputs; systematic divergence traces the evolution of the driven state.

Magnetic stirring mixes the 100 mL vessel within minutes, so the NH₄⁺ response delays of up to 7.7 h belong to the electrode (S5). The electrode's step response must therefore be measured and deconvolved before NH₄⁺ relaxation times are interpreted. pH has no such delay and carries the primary relaxation analysis.

**Probability currents in state space.** Following Battle et al., net probability fluxes are estimated in the pH–ORP and pH–NH₄⁺ planes on the 0.01-day grid. This is done within relaxation segments, with block-bootstrap confidence intervals. A circulating flux with nonzero curl is model-free evidence of broken detailed balance.

**Stability.** Rolling variance and lag-1 autocorrelation are evaluated before each NH₄⁺ addition and each classified process event, to test for critical slowing down (Scheffer et al.). Integrated autocorrelation times are compared across relaxation segments.

**Complexity and directed coupling.** The following measures are computed per relaxation segment:

- permutation entropy (tie-aware), sample entropy and multiscale entropy;
- detrended fluctuation exponents;
- recurrence quantification;
- transfer entropy and Granger causality, conditioned on the remaining variables.

These Shannon-type measures describe the organization of the driven dynamics, not thermodynamic entropy. They can be implemented and unit-tested in your Rust sensor entropy library.

Transfer entropy is computed only between validated signals, and its asymmetry is interpreted with care, since a sensor with slow dynamics can appear to lag or lead artificially.

**Estimator validation.** Each estimator must return values indistinguishable from zero on reversible surrogates. It must also recover the analytic entropy production of a synthetic driven Markov chain, including a version quantized to 0.01 units like the pH record, before it is applied to plant data.

## Phase 4: mechanistic and data-driven modeling

Modeling serves three purposes:

- reconstructing unmeasured states (DO, NO₂⁻, NO₃⁻, alkalinity) needed for complete affinities;
- providing a thermodynamically consistent description of dissipation;
- separating the effects of driving forces from those of loading and temperature.

Models are calibrated chronologically, with every split falling at a pulse so that each relaxation segment belongs wholly to one set. A suitable split calibrates on days 1 to 24.59 (baseline and pulses 1 and 2), selects on days 24.59 to 32.39 (pulse 3) and tests on days 32.39 to 47.92 (pulse 4). Only the selected model is refitted to the whole record, and out-of-sample skill is reported from the split.

**Reduced gray-box model.** Ordinary differential equations for the closed 100 mL batch are written for six states, forced only by the four known NH₄Cl pulses and by gas exchange, following the mass balance in the theoretical framework:

- NH₄⁺ (observed through the calibrated electrode model);
- NO₂⁻, NO₃⁻ and DO (latent unless measured);
- alkalinity, linked to pH through a charge balance with the carbonate system;
- chloride, which is conservative and fixed by the pulses, with conductivity computed from the ionic inventory by Kohlrausch additivity once the probe is validated.

The kinetics treat nitrification in two steps, using Monod terms with oxygen switching functions and with free ammonia and free nitrous acid inhibition. Because ORP is a mixed potential, its observation model is empirical, for example a regression on log DO and pH calibrated on aerated and anoxic periods, and it is never forced to obey a single Nernst equation.

Because the dose and volume are known, each pulse calibrates the NH₄⁺ electrode in situ. The observation model, log R = a + b log C plus an offset at each S2 jump, has its gain b and offsets estimated jointly with the kinetics from the four pulses, and is anchored by laboratory NH₄⁺-N where available. A slow drift in b is allowed, to represent the declining response from pulse 1 to pulse 4. The relaxation between pulses identifies the kinetics, with shared parameters across segments. Time, pH, free ammonia, biomass state and electrode gain are confounded in this record (Phase 1), so three nested model variants are compared:

- one with a slowly drifting biomass activity;
- one with free ammonia inhibition;
- one with both.

The comparison determines whether the roughly fourfold, gain-corrected fall in consumption rate between pulse 2 and pulse 4 reflects inhibition, loss of nitrifier activity, or both. Free nitrous acid inhibition is added once NO₂⁻ data are available. The chloride and alkalinity states are updated at each pulse by the known dose, and the pH observation model is constrained by the alkalinity balance of Phase 2.

Structural identifiability is checked before fitting, and practical identifiability afterward through profile likelihoods. Parameters that the four sensors cannot identify are fixed at literature values and declared as such.

**Full activated sludge model.** ASM2d or ASM3 is coupled to a physicochemical module for pH, ion pairing and gas exchange (Flores-Alsina et al.). It is implemented in open-source QSDsan or in a commercial simulator, and calibrated following the IWA good modeling practice guidelines (Rieger et al.). This model provides plant-scale mass balances, oxygen demand and the influence of sludge retention time, which the reduced model cannot resolve.

**Thermodynamically consistent extension.** Each microbial rate expression is multiplied by the Jin–Bethke factor F\_T computed from the model's own state, and σ(t) is evaluated from the simulated rates and affinities. Two consistency tests follow. First, σ(t) ≥ 0 must hold at every time step, which is a mandatory check of model and parameter consistency. Second, the predicted Gibbs energy dissipation per C-mol biomass is compared with the Heijnen correlation to test whether calibrated yields are thermodynamically plausible.

**State estimation.** An extended or unscented Kalman filter, or moving horizon estimation, assimilates the four sensors into the reduced model to estimate the latent states continuously. An observability analysis determines in advance which states are recoverable from the available measurements. Grab samples are withheld from the filter and used only for validation.

**Data-driven models.** Four complementary approaches extract structure that the mechanistic models may miss:

1. **Dynamic mode decomposition.** DMD, a finite-dimensional Koopman approximation, identifies dominant modes and their decay rates, which are compared with the relaxation constants from Phase 3.
2. **Sparse identification.** SINDy identifies governing equations from data; their structure is compared with the mechanistic kinetics.
3. **Gradient-boosted and temporal convolutional models.** These predict the consumption rate and relaxation time after each pulse from segment features; with only four pulses they serve for exploration, and they become predictive only with replicate experiments. SHAP attribution is used to check that the learned dependencies are physically interpretable.
4. **Hybrid model.** A machine-learning term learns only the residual of the mechanistic model, preserving mass balance by construction.

**Model selection.** Candidate models are compared by AIC and BIC, out-of-sample forecast skill on the held-out days, and residual whiteness. A model that improves fit while violating mass balance or σ ≥ 0 is rejected regardless of its score.

## Phase 5: linking non-equilibrium descriptors to efficiency and energy

The four identical NH₄Cl pulses turn the efficiency question into a controlled experiment. Each pulse injects the same nitrogen load (65.5 mg N/L) and the same free-energy input (163 J in 100 mL), so differences in how each pulse is converted measure changes in the system's capacity, not in its input. In this vessel the only external energy input is the magnetic stirrer, which also sets surface oxygen transfer. The energy question therefore becomes how much stirring power, and how much oxygen transfer, each gram of nitrogen converted requires. These are the quantities that carry over to aeration energy in a full-scale BNR plant.

**Conversion of each pulse.** For pulse k, three quantities are computed from the relaxation segment that follows it, using the known dose rather than the electrode step:

```latex
\eta_k = \frac{\int_{t_k}^{t_{k+1}} r\,dt}{\Delta C}, \qquad \bar r_k = \frac{1}{t_{k+1}-t_k}\int_{t_k}^{t_{k+1}} r\,dt, \qquad \Phi_k = \frac{\bar r_k\,\bar A_k}{T}, \qquad \Delta C = 65.5\ \mathrm{mg\,N\,L^{-1}}
```

Here η\_k is the fraction of the pulse converted before the next pulse, r̄\_k is the mean volumetric consumption rate, and Φ\_k is the mean chemical dissipation, with Ā\_k the mean in situ affinity.

The data already show a clear contrast between two identical pulses. After gain correction, consumption following pulse 2 (about 5.2 mg N/L/d) exceeds that following pulse 4 (about 1.4 mg N/L/d) by about fourfold, although the affinity steps are almost equal (2.5 and 2.4 kJ/mol). The same input converted at very different rates under the same driving force is the central efficiency finding to be explained. The finding is conditional on the electrode gain model (S1) and must be confirmed with laboratory NH₄⁺-N.

**Metrics and their data status.**

| Metric | Definition | Status |
| --- | --- | --- |
| Net NH₄⁺ consumption rate | Slope within each relaxation segment (mg N/L/d) | 0.2 to 2.7 uncorrected for gain; about 1 to 7 after gain correction; laboratory confirmation needed |
| Conversion fraction η\_k | Share of each 65.5 mg N/L pulse consumed before the next pulse | Computable after the electrode gain and offsets are fixed (Phase 4) |
| Affinity step ΔA\_k | RT ln\[(C + 65.5)/C\] at each pulse | 2.4 to 4.2 kJ/mol (Phase 2) |
| Relative ORP response | ORP change in the relaxation after each pulse | Computable now in relative terms; absolute after reference calibration (S4) |
| Total nitrogen conversion | Balance of NH₄⁺, NO₂⁻, NO₃⁻, gaseous and assimilated N against the 262 mg N/L added | Requires NO₃⁻, NO₂⁻ and total N |
| Alkalinity-based nitrification extent | Alkalinity consumed, from the initial alkalinity and the pH record | Requires one initial alkalinity measurement |
| Specific oxygen and energy demand | O₂ transferred through the surface and electrical energy of the stirrer per g N converted | Requires k\_La from a reaeration test, DO, and stirrer power |
| Thermodynamic ratio | Free energy released by nitrification (163 J per fully converted pulse) over the electrical energy of the stirrer during the same relaxation | Requires stirrer power and run time |

**Relationships to test.**

| Descriptor | Expected association | Mechanism |
| --- | --- | --- |
| Free ammonia in each relaxation segment | Higher free ammonia, lower consumption rate; FA of 0.7 to 1.2 mg/L coincides with low rates between days 17.84 and 32.36, but rates stayed low after pulse 4 when FA fell to 0.3 to 0.7 mg/L, so inhibition alone may not explain the decline | Free ammonia inhibition of nitrifiers |
| Affinity step ΔA\_k | Little or no effect on rate | Conversion is kinetically limited (H3) |
| Cumulative chloride and ionic strength | Possible slowing as chloride reaches 18.7 mmol/L above background | Salt stress on nitrifiers |
| Pulse order and elapsed time | Declining conversion with successive pulses | Biomass decay, product accumulation, alkalinity depletion |
| Departure from Onsager regression | Larger departure in segments with inhibited or nonlinear kinetics | Far-from-equilibrium relaxation |
| Relaxation time and autocorrelation before each pulse | Lengthening precedes falling conversion | Loss of resilience in the driven state |

**Statistical design.** The experiment provides four replicated inputs and five segments, so the design has four elements:

1. **Pulse-level analysis.** Compare the kinetic and non-equilibrium descriptors of the four pulses, with pulse order and background state as factors and bootstrap confidence intervals throughout.
2. **Windowed mixed models.** Divide each segment into windows of about 6 hours and fit mixed-effects models with segment random effects and autocorrelation-robust errors.
3. **Nested comparison for H6.** Compare models with and without the non-equilibrium descriptors out of sample, using the Phase 4 split.
4. **Causal graph.** Use an explicit graph (pulse to NH₄⁺ to free ammonia and affinity, aeration to DO to ORP, both to rate, and rate to oxygen and energy demand) to choose adjustment sets, so that mediators are not conditioned on by mistake.

With only four pulses, these analyses are hypothesis-generating. A replicate experiment with more pulses, varied doses and validated sensors is the way to confirm them.

**Mechanistic attribution and scenarios.** Once DO, k\_La and stirrer power are measured, the calibrated model compares the oxygen the surface can transfer with the stoichiometric requirement of 4.57 g O₂ per g N nitrified. The gain-corrected consumption of about 1 to 7 mg N/L/d needs about 6 to 33 mg O₂/L/d, or 0.6 to 3.3 mg O₂/d in the 100 mL vessel. If transfer at the operating stir speed falls short, oxygen limits conversion; if it far exceeds the demand, the extra stirring power is dissipated without nitrogen benefit. Counterfactual simulations then quantify oxygen and energy use against residual NH₄⁺ for three strategies:

- stirring speed matched to the oxygen demand after each pulse;
- pulse size and timing matched to the measured nitrification capacity;
- pH control to keep free ammonia below inhibitory levels.

The specific oxygen demand and the dissipation per gram of nitrogen are the quantities that transfer to a full-scale BNR plant, where oxygen is supplied by blowers and aeration energy dominates.

**Nitrous oxide risk.** With DO and NO₂⁻ available, periods combining high free ammonia, residual NH₄⁺ and low or fluctuating DO after each addition are mapped as N₂O risk windows. Without an N₂O sensor these remain risk proxies and are never reported as emissions.

## Uncertainty, validation and acceptance criteria

The current decisions are as follows:

- **GO.** Chemical disequilibrium of the batch, which follows from the coexistence of NH₄⁺ with oxygen at an affinity of about 300 kJ/mol or more, given the open, stirred vessel; one DO measurement confirms that dissolved O₂ is above zero. Also the smallness of the affinity steps at the pulses, which is insensitive to electrode gain.
- **GO WITH DOCUMENTED NON-BLOCKING DEBT.** The pH-based trajectory analyses, and relative comparisons of consumption between pulses under a stated electrode gain model.
- **NO-GO** for four uses until the open items are closed:
  - absolute NH₄⁺ concentrations, rates and conversion fractions, until the electrode is recalibrated and confirmed by laboratory analysis (S1, S2);
  - any use of the conductivity record as solution conductivity (S3);
  - absolute E\_H and couple potentials (S4);
  - oxygen and energy claims, until DO, k\_La and stirrer power are measured.

Every quantitative claim must carry an uncertainty interval and be traceable to a validated computation. The implementer builds and tests; an independent reviewer inspects the complete cumulative diff and the executed results before any figure enters a manuscript.

**Required tests.**

| Component | Positive and boundary tests | Negative tests |
| --- | --- | --- |
| File parsing | UTF-8 BOM and CRLF are handled; row counts match the verified table; the empty final row of Cond.csv is reported, not dropped silently | Malformed rows raise an error that names the file and row |
| Conductivity conversion | EC = −2.89855 × I (A) + 0.520870 reproduces 0.693 and 1.582 mS/cm at the extreme raw values; a KCl standard addition in the vessel produces the Kohlrausch-predicted step within 10% | A raw value in milliamperes, or a converted value outside 0 to 5 mS/cm, raises an error |
| Time reconstruction | Rows sharing a 0.01-day stamp are spaced evenly and stay monotonic; stamp counts match the verified ranges (43 to 450 for NH₄⁺, 86 to 87 for ORP) | Duplicate times after reconstruction raise an error |
| Segmentation | The four pulses are placed at days 8.53, 10.57, 24.59 and 32.39 on the logger clock, and the four unrelated S2 jumps are flagged as electrode offsets | No relaxation rate, derivative or estimator window spans a pulse or an uncorrected electrode jump |
| Speciation and Nernst | pK\_a at 25 °C reproduces 9.25; formal potentials reproduce the table; the bare Ag (AgCl film) conversion reproduces 0.2223 V at unit chloride activity and a 59.2 mV shift per decade of chloride at 25 °C; the predicted chloride staircase rises by 4.67 mmol/L per pulse; Davies coefficients match published values | Physically impossible inputs are rejected |
| Irreversibility estimators | Analytic entropy production of a driven Markov chain is recovered, including at 0.01-unit quantization | Reversible surrogates give values indistinguishable from zero |
| Models | Mass and charge balances close; σ ≥ 0 at every step; each simulated pulse adds exactly 65.5 mg N/L and 4.67 mmol/L chloride; the fitted electrode gain reproduces the observed 22 to 47% responses | Parameter sets violating either are rejected |

**Acceptance criteria.**

1. **Data integrity.** Each sensor reproduces a known standard addition in the vessel: NH₄Cl for the NH₄⁺ electrode, KCl for conductivity, and NaCl for the Ag reference. Laboratory NH₄⁺-N, chloride and alkalinity confirm the calibrated records (S1 to S4).
2. **Time base.** The logger start times are confirmed, and multivariate analyses are restricted to the verified common window.
3. **Event classification.** The electrode's step response explains its delayed responses to pulses 1 to 3 (S5), and every other change point is attributed to a sensor, handling or process event.
4. **Disequilibrium.** Disequilibrium indices exceed their Monte Carlo intervals wherever non-equilibrium is claimed.
5. **Validation.** Latent-state estimates are validated against withheld grab samples before they enter any affinity calculation.
6. **Generalization.** Reported relationships survive chronological out-of-sample tests and false discovery rate control.
7. **Reproducibility.** One command regenerates all results from the raw files, in a pinned environment with fixed seeds.

**Severity and decision rules.**

- **P0** findings invalidate conclusions, for example unconfirmed units, rates computed across a perturbation, misaligned clocks, or a second-law violation.
- **P1** findings undermine specific claims, for example an unvalidated estimator, missing uncertainty, or leakage between training and test periods.
- **P2** findings concern documentation.

Any unresolved P0 or P1 finding, failed test or unverified criterion results in NO-GO. The implementer is never the sole approver, and decisions are GO, GO WITH DOCUMENTED NON-BLOCKING DEBT, or NO-GO, based on evidence.

**Scientific failure conditions.** The non-equilibrium claim is restricted if indices fall within uncertainty or cannot be distinguished from surrogates. Conclusions that depend on a sensor are withdrawn if that sensor fails its standard-addition test after recalibration. The pulse-to-pulse rate contrast is withdrawn if laboratory NH₄⁺-N does not reproduce it. H3 is rejected if the consumption rate scales with the affinity step rather than with substrate, free ammonia and biomass state. The efficiency and energy conclusions are reported as associations only if H6 fails out of sample.

## Records and measurements needed

The highest-value next step is a short validation run in the same vessel with the same probes, in which known standard additions test each sensor directly. Together with a few laboratory analyses and a description of the vessel conditions, it would convert this dataset from qualitative to quantitative and make a replicate experiment trustworthy.

| Priority | Item | What it resolves or enables |
| --- | --- | --- |
| 1 | Standard additions in the stirred vessel: KCl, recorded by conductivity and checked with a calibrated handheld meter; NaCl, recorded by ORP; and NH₄Cl at two or three levels, recorded by the NH₄⁺ electrode and checked in the laboratory | Tests S1, S3 and S4 directly; gives the NH₄⁺ electrode gain (Gran plot), the conductivity response and the chloride response of the bare Ag |
| 1 | NH₄⁺ electrode step response in a stirred standard, with the transmitter's averaging, hold and update settings and its calibration history | Explains the 7.5 h delays and the jumps (S1, S2, S5) |
| 1 | DO, and k\_La at the operating stir speed by a reaeration test (deoxygenate, then record the DO recovery) | Confirms the oxygen premise of the affinity criterion; oxygen supply and its possible limitation of nitrification |
| 1 | Stirrer power, temperature, and vessel weight at start and end | Energy metrics, equilibrium constants and free ammonia, and the evaporation loss of the open vessel |
| 2 | Laboratory NH₄⁺-N, NO₃⁻, NO₂⁻, total N, chloride and alkalinity, on any stored samples and at the end of the current batch | Nitrogen, chloride and alkalinity balances; confirms the consumption trends |
| 2 | Initial alkalinity, or that of the source wastewater | Bounds cumulative nitrification from the pH record and gives the CO₂ supersaturation |
| 2 | Ag surface inspection, a ZoBell's check, and preferably replacement of the bare Ag with a double-junction Ag/AgCl reference | Resolves S4 and gives absolute E\_H in future runs |
| 3 | Intensive sampling every 15 to 30 min after a future pulse | Ground truth for the relaxation analysis and the Onsager regression test |
| 3 | MLSS and MLVSS, and NH₃ and N₂O in the headspace | Biomass-specific rates and gaseous nitrogen losses |

## Deliverables and manuscript structure

The work proceeds as a gated sequence. Phases 2 and 3 run in parallel on the Phase 1 perturbation segmentation, and Phase 4 feeds estimated latent states back into Phase 2 so that complete affinities can be computed.

&#91;embedded content: analysis workflow · six phases, one feedback loop, two review gates\]

The thermodynamic and the trajectory-based evidence are deliberately independent. Agreement between them is therefore the strongest available demonstration of non-equilibrium.

**Deliverables.**

1. **Verified dataset.** A verified, flagged dataset with its data dictionary, gap and event inventory, and QC report.
2. **Analysis code.** A version-controlled codebase with the tests specified above, a pinned environment and a single-command pipeline.
3. **Feature matrix.** The perturbation-level and segment-level matrix of features and non-equilibrium descriptors, with uncertainty intervals.
4. **Calibrated models.** Calibrated reduced and full models, with their identifiability analysis and validated latent-state estimates.
5. **Operational recommendations.** These cover:
   - the aeration control strategy;
   - DO set points;
   - alkalinity dosing triggers;
   - nitrite shunt feasibility;
   - N₂O risk windows.
6. **Manuscript.** A manuscript built on the figure plan below.

**Proposed figures.**

| Figure | Content | Phase |
| --- | --- | --- |
| 1 | Conceptual framework: the closed 100 mL batch, the four NH₄Cl pulses as free-energy inputs, gas exchange, and the part of the system each sensor observes | Theory |
| 2 | Common-window overview (day 0 to 47.92) with the four pulses, QC flags and electrode artifacts | 0 and 1 |
| 3 | Sensor validation: predicted and observed responses of each sensor to each pulse, and the predicted chloride and conductivity staircase against the record | 0 |
| 4 | Relaxation after each pulse, gain-corrected consumption rate against free ammonia and affinity step, and the alkalinity-bounded extent of nitrification | 1 and 2 |
| 5 | Onsager regression test, pulse-to-pulse comparison of responses, irreversibility indices against surrogates, and probability flux maps | 3 |
| 6 | Model fit with the electrode observation model, latent states against laboratory samples, and σ(t) with F\_T | 4 |
| 7 | Conversion fraction and dissipation for each pulse, and partial effects of non-equilibrium descriptors on consumption rate and oxygen demand | 5 |
| 8 | Counterfactual aeration and pH control strategies: oxygen and energy saving against residual NH₄⁺ | 5 |

## References

The references below are cited from established literature without online verification in this session. Volume and page details should be confirmed before they enter a manuscript.

- Anthonisen, A. C., et al. (1976). Inhibition of nitrification by ammonia and nitrous acid. *Journal of the Water Pollution Control Federation*.
- Bandt, C., and Pompe, B. (2002). Permutation entropy: a natural complexity measure for time series. *Physical Review Letters*.
- Battle, C., et al. (2016). Broken detailed balance at mesoscopic scales in active biological systems. *Science*.
- Emerson, K., et al. (1975). Aqueous ammonia equilibrium calculations: effect of pH and temperature. *Journal of the Fisheries Research Board of Canada*.
- Flores-Alsina, X., et al. (2015). A plant-wide aqueous phase chemistry module describing pH variations and ion speciation/pairing in wastewater treatment process models. *Water Research*.
- Heijnen, J. J., and van Dijken, J. P. (1992). In search of a thermodynamic description of biomass yields for the chemotrophic growth of microorganisms. *Biotechnology and Bioengineering*.
- Henze, M., et al. (2000). *Activated Sludge Models ASM1, ASM2, ASM2d and ASM3*. IWA Publishing.
- Jin, Q., and Bethke, C. M. (2005). Predicting the rate of microbial respiration in geochemical environments. *Geochimica et Cosmochimica Acta*.
- Kawai, R., Parrondo, J. M. R., and Van den Broeck, C. (2007). Dissipation: the phase-space perspective. *Physical Review Letters*.
- Kleerebezem, R., and van Loosdrecht, M. C. M. (2010). A generalized method for thermodynamic state analysis of environmental systems. *Critical Reviews in Environmental Science and Technology*.
- Lacasa, L., et al. (2012). Time series irreversibility: a visibility graph approach. *European Physical Journal B*.
- Lindberg, R. D., and Runnells, D. D. (1984). Ground water redox reactions: an analysis of equilibrium state applied to Eh measurements and geochemical modeling. *Science*.
- Rieger, L., et al. (2012). *Guidelines for Using Activated Sludge Models*. IWA Scientific and Technical Report No. 22.
- Rittmann, B. E., and McCarty, P. L. *Environmental Biotechnology: Principles and Applications*. McGraw-Hill.
- Roldán, É., and Parrondo, J. M. R. (2010). Estimating dissipation from single stationary trajectories. *Physical Review Letters*.
- Russell, L. L. (1976). Chemical aspects of groundwater recharge with wastewaters. PhD thesis, University of California, Berkeley.
- Scheffer, M., et al. (2009). Early-warning signals for critical transitions. *Nature*.
- Schreiber, T. (2000). Measuring information transfer. *Physical Review Letters*.
- Schreiber, T., and Schmitz, A. (2000). Surrogate time series. *Physica D*.
- Stumm, W., and Morgan, J. J. *Aquatic Chemistry*. Wiley.
- Onsager, L. (1931). Reciprocal relations in irreversible processes. II. *Physical Review*.
- Kubo, R. (1966). The fluctuation–dissipation theorem. *Reports on Progress in Physics*.
- Bard, A. J., and Faulkner, L. R. *Electrochemical Methods: Fundamentals and Applications*. Wiley.
- Haldane, J. B. S. (1930). *Enzymes*. Longmans, Green and Co.
