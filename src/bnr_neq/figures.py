"""Figures 1, 4 and 8 (W6), drawn from the generated tables only.

Fig. 1: the linkage sensors -> process state -> rates -> O2 demand -> energy -> efficiency, with each link's basis.
Fig. 4: the W2 bound N_nit(Alk0) against the removal implied by C(t); rate against free ammonia (W1); state timeline (W3).
Fig. 8: per-pulse conversion, oxygen utilization, idle-fraction bounds with coverage and aeration energy under the
field-AE and kLa SCENARIOS.

Palette: the project's fixed categorical order (validated: CVD ΔE ≥ 9 for adjacent pairs; aqua is below 3:1 contrast on
the surface, so every series is also named in a legend or direct label, and the report carries the tables).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
MUTED, FAINT = "#9a988f", "#d6d4cd"
STATE_COLORS = {"active_nitrification": AQUA, "inactive_oversupplied": ORANGE, "oxygen_limited": VIOLET,
                "sensor_artifact": MUTED, "unclassified": FAINT}
FIGURES_W6 = ("fig1_linkage", "fig4_constraints", "fig8_efficiency")


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _save(fig, d: Path, name: str, plt) -> None:
    fig.savefig(d / f"{name}.png", dpi=150, facecolor=SURFACE, bbox_inches="tight")
    fig.savefig(d / f"{name}.svg", facecolor=SURFACE, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def _metric(pm: pd.DataFrame, metric: str, col: str = "value") -> pd.Series:
    return pm[pm["metric"] == metric].set_index("pulse")[col]


# --- Fig. 1 ----------------------------------------------------------------------------------------------------
def fig1(t: dict[str, pd.DataFrame], d: Path, saes: list[float]) -> None:
    from .report import pyplot

    plt = pyplot()
    v = t["sensor_verdicts.csv"].set_index("key")["verdict"].astype(str)
    fr = t["states_fractions.csv"].set_index("interval").loc["whole record"]
    cons = t["nh4_truth_w2_consistency.csv"]
    first_bad = cons.loc[~cons["consistent_at_best_fit"].astype(bool), "checkpoint"]
    boxes = [
        ("Sensors", (f"pH: {v['ph'].split(' (')[0]}\nNH4: {v['nh4'].split(':')[0]} (under-response)\n"
                     f"ORP: {v['orp'].split(':')[0]} (absolute);\nrelative ORP used\nEC: {v['ec'].split(':')[0]}"), BLUE),
        ("NH4 truth C(t)", "M6 log observation\nmodel, inverted\n(level not identifiable)", VIOLET),
        ("Process state", (f"rules on ORP_rel, pH\nand C(t) slopes:\n{100 * fr['classified_fraction']:.0f} % classified;\n"
                           f"idle {100 * fr['idle_fraction_lower_bound']:.0f}-{100 * fr['idle_fraction_upper_bound']:.0f} % "
                           "of\nstirred time (bounds)"), BLUE),
        ("Rates", "−dC/dt per slope window,\nremoval per pulse, η", VIOLET),
        ("O2 demand", "4.57 g O2/g N\n(3.43 nitritation case)", VIOLET),
        ("Energy", f"O2 / field AE\n(SCENARIO {min(saes):g}-{max(saes):g}\nkg O2/kWh);\nsupply kLa (C*−C) V\n(SCENARIO kLa)",
         ORANGE),
        ("Efficiency", "η, O2 utilization,\nidle-fraction bounds,\nkWh/kg N,\nΔG released / energy", ORANGE),
    ]
    fig, ax = plt.subplots(figsize=(15, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 7)
    ax.set_ylim(-1.35, 1.3)
    ax.axis("off")
    for i, (title, body, color) in enumerate(boxes):
        x = i + 0.5
        ax.add_patch(plt.Rectangle((x - 0.44, -0.55), 0.88, 1.35, facecolor=SURFACE, edgecolor=color, lw=2))
        ax.text(x, 0.66, title, ha="center", va="center", fontsize=10, color=INK, weight="bold")
        ax.text(x, 0.08, body, ha="center", va="center", fontsize=7.2, color=INK2, linespacing=1.3)
        if i < len(boxes) - 1:
            ax.annotate("", xy=(x + 0.58, 0.1), xytext=(x + 0.42, 0.1),
                        arrowprops={"arrowstyle": "-|>", "color": INK2, "lw": 1.2})
    # the pH-based W2 constraint acts on the rates
    ax.add_patch(plt.Rectangle((0.06, -1.3), 3.5, 0.5, facecolor=SURFACE, edgecolor=BLUE, lw=1.2, ls="--"))
    bad = (f"C(t) contradicts it from {first_bad.iloc[0]} on (QUALIFIED)" if len(first_bad)
           else "C(t) consistent with it")
    ax.text(1.81, -1.05, f"W2 alkalinity bound from the pH record: N_nit ≤ M_N (Alk0 − Alk_min(pH))/2\n{bad}",
            ha="center", va="center", fontsize=7.6, color=INK2)
    ax.annotate("", xy=(3.5, -0.56), xytext=(3.2, -0.8), arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 1.1})
    for j, (lab, c) in enumerate((("measured / validated sensor input", BLUE), ("model-based (one observation model)", VIOLET),
                                  ("SCENARIO-dependent", ORANGE))):
        ax.plot([4.1], [-0.85 - 0.2 * j], marker="s", ms=9, color=c, ls="none")
        ax.text(4.2, -0.85 - 0.2 * j, lab, va="center", fontsize=7.6, color=INK2)
    ax.set_title("Figure 1. From sensor signals to process state, oxygen demand, energy and efficiency (box border = "
                 "basis)", loc="left", fontsize=10.5, color=INK)
    _save(fig, d, "fig1_linkage", plt)


# --- Fig. 4 ----------------------------------------------------------------------------------------------------
def fig4(t: dict[str, pd.DataFrame], states_bins: pd.DataFrame, pulse_times: list[float], d: Path) -> None:
    from . import alkalinity, volume
    from .report import pyplot

    plt = pyplot()
    fig = plt.figure(figsize=(12, 8.6), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.75], hspace=0.38, wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    _style(ax)
    b = t["alkalinity_bound.csv"]
    base = b[(b["volume_scenario"] == volume.CONSTANT) & (b["alk0_basis"] != "MEASURED")]
    for T, color, ls in ((25.0, BLUE, "-"), (20.0, BLUE, ":")):
        g = base[(base["co2_scenario"] == alkalinity.AIR_EQ) & (base["temperature_C"] == T)]
        ax.plot(g["alk0_meq_L"], g["n_nit_max_end_mgN_L"], color=color, lw=2, ls=ls,
                label=f"W2 bound by day 47.9, air-equilibrium CO2, SCENARIO {T:g} °C")
    any_s = base[(base["co2_scenario"] == alkalinity.AIR_EQ) & (base["temperature_C"] == 25.0)]
    ax.plot(any_s["alk0_meq_L"], any_s["n_removed_max_any_proton_sink_end_mgN_L"], color=VIOLET, lw=1.5, ls="--",
            label="bound if every removed N released 1 eq/mol (e.g. NH3 stripping)")
    dosed = float(b["dosed_total_mgN_L"].iloc[0])
    ax.axhline(dosed, color=INK2, lw=1, ls="-.")
    ax.text(14.9, dosed + 4, f"dosed {dosed:.0f} mg N/L", fontsize=7.5, color=INK2, ha="right")
    cons = t["nh4_truth_w2_consistency.csv"].set_index("checkpoint").loc["end"]
    ax.axhspan(cons["removed_band_low_mgN_L"], cons["removed_band_high_mgN_L"], color=ORANGE, alpha=0.12, lw=0)
    ax.axhline(cons["removed_per_initial_litre_mgN_L"], color=ORANGE, lw=2)
    ax.text(1.1, cons["removed_per_initial_litre_mgN_L"] + 4,
            f"removal implied by C(t), best fit {cons['removed_per_initial_litre_mgN_L']:.0f} mg N/L "
            "(band shaded, open upward in C)", fontsize=7.5, color=INK2)
    ax.set_xlabel("Initial alkalinity Alk0 (meq/L, scanned; not measured)", fontsize=8.5, color=INK2)
    ax.set_ylabel("Cumulative N removed by day 47.9 (mg N/L)", fontsize=8.5, color=INK)
    ax.set_title("a  Alkalinity bound against the inverted C(t) (constant volume)", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="lower right")

    ax = fig.add_subplot(gs[0, 1])
    _style(ax)
    tw = t["nh4_truth_windows.csv"]
    tw = tw[tw["status"] == "computed"]
    fa_col = "fa_median_mgNH3_L_SCENARIO_T25C"
    for k, color in zip((1, 2, 3, 4), (BLUE, ORANGE, AQUA, VIOLET)):
        g = tw[tw["governing_pulse"] == k]
        if not len(g):
            continue
        err = np.vstack([g["consumption_mgN_L_d"] - g["consumption_ci95_low_mgN_L_d"],
                         g["consumption_ci95_high_mgN_L_d"] - g["consumption_mgN_L_d"]])
        ax.errorbar(g[fa_col], g["consumption_mgN_L_d"], yerr=err, fmt="o", ms=8, color=color, ecolor=color,
                    elinewidth=1, capsize=2, mec=SURFACE, mew=1.5, label=f"after pulse {k}")
        for x, y in zip(g[fa_col], g["consumption_mgN_L_d"]):
            ax.annotate(f"P{k}", (x, y), textcoords="offset points", xytext=(6, 4), fontsize=7.5, color=INK2)
    ax.set_xlabel("Median free ammonia from C(t) (mg NH3/L, SCENARIO 25 °C)", fontsize=8.5, color=INK2)
    ax.set_ylabel("Consumption −dC/dt (mg N/L/d, Theil-Sen 95 %)", fontsize=8.5, color=INK)
    ax.set_title("b  Rate against free ammonia, slope windows (W1; pulses 2-4 QUALIFIED)", loc="left", fontsize=10,
                 color=INK)
    ax.legend(fontsize=7.5, frameon=False)

    ax = fig.add_subplot(gs[1, :])
    _style(ax)
    ax.grid(False)
    tl = states_bins
    for i, (state, color) in enumerate(STATE_COLORS.items()):
        m = (tl["state"] == state).to_numpy()
        ax.broken_barh([(x - 0.005, 0.01) for x in tl.loc[m, "t_mid_d"]], (len(STATE_COLORS) - 1 - i - 0.4, 0.8),
                       facecolors=color, linewidth=0)
    for tp in pulse_times:
        ax.axvline(tp, color=INK2, lw=0.8, ls=":")
    ax.set_yticks(range(len(STATE_COLORS)), [s.replace("_", " ") for s in reversed(STATE_COLORS)], fontsize=8)
    ax.set_xlabel("Time (d); dotted lines: NH4Cl pulses", fontsize=8.5, color=INK2)
    ax.set_title("c  Process-state timeline (registered rules, revised after first results; 0.01 d bins)", loc="left", fontsize=10, color=INK)
    fig.suptitle("Figure 4. Constraints on the conversion: alkalinity, free ammonia and process state", x=0.06,
                 ha="left", fontsize=11, color=INK)
    _save(fig, d, "fig4_constraints", plt)


# --- Fig. 8 ----------------------------------------------------------------------------------------------------
def fig8(t: dict[str, pd.DataFrame], d: Path) -> None:
    from .report import pyplot

    plt = pyplot()
    pm = t["pulses.csv"]
    pulses = sorted(pm["pulse"].unique())
    x = np.arange(len(pulses))
    labels = [f"P{int(p)}" for p in pulses]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), facecolor=SURFACE)
    for ax in axes.flat:
        _style(ax)

    ax = axes[0, 0]
    eta = _metric(pm, "eta_conversion_fraction").reindex(pulses)
    lo = _metric(pm, "eta_conversion_fraction", "param_band_low").reindex(pulses)
    hi = _metric(pm, "eta_conversion_fraction", "param_band_high").reindex(pulses)
    ax.bar(x, np.nan_to_num(eta.to_numpy(dtype=float)), 0.5, color=ORANGE, label="η, constant volume (best fit)")
    ax.vlines(x, lo, hi, color=INK, lw=1.2, label="profile-vector band")
    for metric, marker, lab in (("removal_to_dose_ratio_SCENARIO_evaporation_0.5_mL_per_d", "^", "removal/dose, 0.5 mL/d"),
                                ("removal_to_dose_ratio_SCENARIO_evaporation_1_mL_per_d", "v", "removal/dose, 1 mL/d")):
        ax.plot(x + 0.33, _metric(pm, metric).reindex(pulses), marker, ms=8, color=VIOLET, ls="none",
                label=f"SCENARIO evaporation: {lab}")
    for xi, e in zip(x, eta):
        if np.isnan(e):
            ax.text(xi, 0.03, "not\ncomputed", ha="center", fontsize=7, color=INK2)
    ax.axhline(1.0, color=INK2, lw=0.8, ls="--")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Fraction of the pulse removed before the next one", fontsize=8.5, color=INK)
    ax.set_ylim(0, 1.55)
    qual = sorted({int(k) for k, q in zip(pm["pulse"], pm["w2_bound_check"].fillna("").astype(str)) if q.startswith("QUALIFIED")})
    ax.set_title("a  Conversion η per pulse (W1 C(t)" + (f"; {', '.join(f'P{k}' for k in qual)} QUALIFIED by W2)"
                                                        if qual else ")"), loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)

    ax = axes[0, 1]
    klas = sorted({float(m.split("_kLa")[1].split("_")[0]) for m in pm["metric"] if m.startswith("o2_utilization_SCENARIO_kLa")})
    for j, (kla, color) in enumerate(zip(klas, (BLUE, AQUA, VIOLET))):
        u = _metric(pm, f"o2_utilization_SCENARIO_kLa{kla:g}_T25C").reindex(pulses).to_numpy(dtype=float)
        ax.bar(x + (j - 1) * 0.26, np.nan_to_num(u), 0.24, color=color, label=f"kLa {kla:g} 1/d (SCENARIO)")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 5)
    ax.axhline(1.0, color=INK2, lw=0.8, ls="--")
    for j, kla in enumerate(klas):
        u = _metric(pm, f"o2_utilization_SCENARIO_kLa{kla:g}_T25C", "status").reindex(pulses).astype(str)
        for xi, st in zip(x, u):
            if st.startswith("not computed: SCENARIO infeasible"):
                ax.text(xi + (j - 1) * 0.26, 1.15, "infeasible", rotation=90, ha="center", va="bottom", fontsize=6.5,
                        color=INK2)
    ax.set_xticks(x, labels)
    ax.set_ylabel("O2 required / O2 transfer capacity (log)", fontsize=8.5, color=INK)
    ax.set_title("b  Oxygen utilization (DO not measured: capacity at DO = 0)", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="lower left")

    ax = axes[1, 0]
    lo_i = _metric(pm, "idle_fraction_LOWER_BOUND").reindex(pulses).to_numpy(dtype=float)
    hi_i = _metric(pm, "idle_fraction_UPPER_BOUND").reindex(pulses).to_numpy(dtype=float)
    cov = _metric(pm, "classified_fraction").reindex(pulses).to_numpy(dtype=float)
    cond = _metric(pm, "idle_fraction_of_classified_time_CONDITIONAL").reindex(pulses).to_numpy(dtype=float)
    ax.bar(x, np.nan_to_num(hi_i - lo_i), 0.42, bottom=np.nan_to_num(lo_i), color=ORANGE, alpha=0.35,
           label="idle fraction: between the lower and upper bound")
    ax.hlines(lo_i, x - 0.21, x + 0.21, color=ORANGE, lw=2, label="lower bound (inactive bins / all bins)")
    ax.plot(x, cond, "D", ms=6, color=BLUE, ls="none", label="CONDITIONAL, over classified time (coverage ≥ 0.5)")
    for xi, c, h in zip(x, cov, hi_i):
        ax.text(xi, min(np.nan_to_num(h) + 0.03, 1.02), f"coverage\n{c:.2f}", ha="center", va="bottom", fontsize=7,
                color=INK2)
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.6)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_ylabel("Idle fraction of stirred time (W3 states)", fontsize=8.5, color=INK)
    ax.set_title("c  Time stirred without conversion: bounds, not a value", loc="left", fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=2)

    ax = axes[1, 1]
    saes = sorted({float(m.split("_AE")[1].split("_")[0]) for m in pm["metric"]
                   if m.startswith("aeration_energy_supply_kWh_per_kgN_SCENARIO_AE")})
    mid = saes[len(saes) // 2]
    for j, (kla, color) in enumerate(zip(klas, (BLUE, AQUA, VIOLET))):
        vals = {s: _metric(pm, f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{s:g}_kLa{kla:g}_T25C").reindex(pulses)
                .to_numpy(dtype=float) for s in saes}
        xx = x + (j - 1) * 0.26
        ax.bar(xx, np.nan_to_num(vals[mid]), 0.24, color=color, label=f"supply-based, kLa {kla:g} 1/d, AE {mid:g}")
        ax.vlines(xx, vals[max(saes)], vals[min(saes)], color=INK, lw=1)
    floor = _metric(pm, f"aeration_energy_demand_kWh_per_kgN_SCENARIO_AE{mid:g}").dropna()
    if len(floor):
        ax.axhline(float(floor.iloc[0]), color=ORANGE, lw=1.5, ls="--",
                   label=f"demand-based floor 4.57/AE, AE {mid:g}")
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Aeration-equivalent energy (kWh/kg N removed, log)", fontsize=8.5, color=INK)
    ax.set_title(f"d  Energy per N removed (whiskers: field AE {min(saes):g}-{max(saes):g} kg O2/kWh)", loc="left",
                 fontsize=10, color=INK)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    eta_st = _metric(pm, "eta_conversion_fraction", "status").reindex(pulses).fillna("").astype(str)
    nc_txt = "; ".join(f"P{int(p)}: η not computed, {'removal negative' if 'negative' in st else 'removal above the dose'}"
                       for p, st in eta_st.items() if st != "computed")
    fig.suptitle("Figure 8. Per-pulse efficiency under the field-AE and kLa SCENARIOS" + (f" ({nc_txt})" if nc_txt else ""),
                 x=0.06, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    _save(fig, d, "fig8_efficiency", plt)


def make_all(cfg, t: dict[str, pd.DataFrame], out_dir: Path) -> None:
    d = Path(out_dir) / "figures"
    states_bins = pd.read_parquet(Path(out_dir) / "processed" / "states_bins.parquet")
    fig1(t, d, [float(x) for x in cfg.get("scenarios.sae_field_kgO2_per_kWh")])
    fig4(t, states_bins, cfg.pulse_times_d, d)
    fig8(t, d)
