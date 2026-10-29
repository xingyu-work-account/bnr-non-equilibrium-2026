"""Every figure of the pipeline, drawn from the files under ``outputs/`` only (no recomputation).

One canvas and one axes per file, no title and no tight layout; what a figure shows is in its file name, and the
report carries the full caption (:data:`FIGURE_GROUPS`). Every figure uses :data:`PARAMS`, a 10 x 6 inch canvas,
all four spines and no grid, and is saved at 300 dpi on a transparent background as PNG and SVG; a copy without any
text label goes to ``figures/no_text/`` under the same name.

``python -m bnr_neq figures`` redraws all of them from the saved tables; each stage also draws its own after writing
its tables, and the report stage draws figures 1, 4 and 8.

Palette: the project's fixed categorical order (validated: CVD ΔE ≥ 9 for adjacent pairs; aqua is below 3:1 contrast on
a light background, so every series is also named in a legend, and the report carries the tables).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d6d4cd"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
MUTED, FAINT = "#9a988f", "#d6d4cd"
STATE_COLORS = {"active_nitrification": AQUA, "inactive_oversupplied": ORANGE, "oxygen_limited": VIOLET,
                "sensor_artifact": MUTED, "unclassified": FAINT}
STATE_LABELS = {"active_nitrification": "active", "inactive_oversupplied": "idle", "oxygen_limited": "O2-lim.",
                "sensor_artifact": "artifact", "unclassified": "none"}

FONT_SIZE = 25
PARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "font.size": FONT_SIZE,
    "font.weight": 600,
    "axes.labelsize": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "axes.labelweight": 600,
    "axes.titleweight": 600,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "axes.linewidth": 3,
    "figure.dpi": 300,
    "figure.subplot.left": 0.23,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.92,
    "lines.linewidth": 3,
    "savefig.dpi": 300,
}
FIGSIZE = (10, 6)
NO_TEXT_DIR = "no_text"  # the same figures without any text label, same file names

# (group key used by the report, group caption, [(file stem, panel caption), ...]); the stems are the file names
FIGURE_GROUPS: list[tuple[str, str, list[tuple[str, str]]]] = [
    ("fig1", "Figure 1. From sensor signals to process state, oxygen demand, energy and efficiency (W6)", [
        ("fig1_linkage_sensors_state_O2_energy_efficiency",
         ("Sensors (verdicts, §3) → NH4 truth C(t) (the M6 log observation model, inverted) → process state (rules on "
         "relative ORP, pH and C(t) slopes; whole-record idle-fraction bounds) → rates (−dC/dt per slope window, "
         "removal and η per pulse) → O2 demand (4.57 g O2/g N) → aeration-equivalent energy (field AE SCENARIO range, "
         "kg O2/kWh) → efficiency (kWh/kg N, O2 utilization, idle bounds). Dashed: the W2 alkalinity bound from the pH "
         "record, which C(t) contradicts or holds. Box border = basis: blue measured or validated sensor input, violet "
         "model-based (one observation model), orange SCENARIO-dependent")),
    ]),
    ("fig2", "Figure 2. Common-window overview (0.01 d grid) with pulses, flags and electrode events (M3)", [
        ("fig2a_overview_NH4_slope_windows_jumps_held",
         ("NH4+ electrode: Theil-Sen slope windows, located jumps (pulse response, electrode offset), step candidates "
         "and held/damped output (bar at top); dotted: NH4Cl pulses; grey band: start-up (day 0 to 1)")),
        ("fig2b_overview_pH", "pH with change-point events not tied to a pulse"),
        ("fig2c_overview_ORP_relative", "ORP (mV vs bare Ag, relative) with change-point events not tied to a pulse"),
        ("fig2d_overview_EC_unvalidated", "EC (unvalidated) with change-point events and late EC pulses (probe events)"),
    ]),
    ("fig3", "Figure 3. Sensor validation with the known pulses (M2)", [
        ("fig3a_validation_EC_predicted_staircase_vs_record",
         ("Conductivity: predicted staircase (limiting λ°, SCENARIO 25 °C) for no and for complete nitrification "
         "against the record")),
        ("fig3b_validation_NH4_jumps_held_segments",
         "NH4+ electrode: detected jumps and held/damped segments (σ_HF < 0.01, shaded)"),
        ("fig3c_validation_NH4_change_per_pulse", "NH4+ change per pulse: predicted from the dose vs observed"),
        ("fig3d_validation_EC_change_per_pulse",
         "EC change per pulse: predicted (limiting λ°, SCENARIO 25 °C) vs observed"),
        ("fig3e_validation_ORP_change_per_pulse",
         "ORP change per pulse: predicted range over the Cl- and temperature SCENARIOS vs observed"),
        ("fig3f_validation_pH_change_per_pulse", ("pH change per pulse, observed (prediction not computed: alkalinity "
         "not measured)")),
    ]),
    ("fig5", "Figure 5. Estimator validation, irreversibility, Onsager SNR, probability flux (M5)", [
        ("fig5a_estimator_validation_ring_chain",
         ("Estimator validation on a three-state ring chain: transition estimator (99 % CI) and a noisy observable "
         "quantized to 0.01, against the analytic (p − q) ln(p/q)")),
        ("fig5b_ordinal_KL_irreversibility_stationary_windows",
         ("Time irreversibility (tie-aware ordinal-pattern KL) in stationary 0.25 d windows; filled: significant "
         "against IAAFT surrogates")),
        ("fig5c_Onsager_SNR_per_pulse",
         "Onsager regression test: |pH pulse response| / fluctuation SD; below the line the test is not interpretable"),
        ("fig5d_probability_flux_circulation", "Probability-flux circulation rate per segment (Battle et al.; 95 % CI)"),
    ]),
    ("fig6", "Figure 6. Reduced model with the electrode observation model (M6)", [
        ("fig6a_model_selected_vs_electrode_cal_sel_test",
         ("Selected model with the electrode observation model: calibrated on the calibration window and frozen "
         "afterwards, and the whole-record refit; shaded: selection and test windows")),
        ("fig6b_model_consumption_rate_whole_record_refit", "Consumption rate of the whole-record refit"),
        ("fig6c_profile_likelihood_identifiability",
         ("Profile likelihood: maximum Δ reached on each side (n_eff-scaled, capped at 40) against the 3.84 line; "
         "hollow: not identifiable; o1-o8: electrode offsets in time order")),
    ]),
    ("fig7", "Figure 7. Pulse-level descriptors (M7)", [
        ("fig7a_consumption_after_each_pulse_95CI",
         ("Consumption after each pulse (positive slope windows, 95 % CI): electrode reading and M6-inverted NH4; "
         "bar: SENSITIVITY linear gain correction")),
        ("fig7b_affinity_step_M6_inverted_NH4_25C",
         ("CONSISTENCY CHECK (not a finding): affinity step ΔA at each pulse on the M6-inverted NH4 (SCENARIO 25 °C), "
          "with its share of |ΔG°′|")),
        ("fig7c_free_ammonia_vs_consumption_EXPLORATORY",
         ("Free ammonia from the M6-inverted TAN against the M6-inverted consumption, per temperature SCENARIO "
         "(EXPLORATORY, n = 4 pulses)")),
        ("fig7d_conversion_eta_before_next_pulse",
         ("Conversion η before the next pulse: dose-anchored on the M6-inverted C(t) and MODEL-BASED (M6 rate law); "
         "the last pulse runs to the end of the record")),
    ]),
    ("fig4", "Figure 4. Alkalinity bound, rate against free ammonia, state timeline (W6)", [
        ("fig4a_alkalinity_bound_vs_Ct_removal_day47p9",
         ("W2 alkalinity bound by day 47.9 (air-equilibrium CO2, SCENARIO 25 and 20 °C; dashed: every removed N "
         "releasing 1 eq/mol) against the removal implied by C(t) (best fit, band shaded); dash-dot: dosed N")),
        ("fig4b_rate_vs_free_ammonia_W1_slope_windows",
         ("Consumption −dC/dt (Theil-Sen, 95 %) against median free ammonia from C(t) (SCENARIO 25 °C), per slope "
         "window (W1)")),
        ("fig4c_process_state_timeline_W3",
         "Process-state timeline (registered rules, revised after first results; 0.01 d bins); dotted: NH4Cl pulses"),
    ]),
    ("fig8", "Figure 8. Per-pulse η, oxygen utilization, idle fraction and energy under SCENARIOS (W6)", [
        ("fig8a_conversion_eta_per_pulse_evaporation_SCENARIO",
         ("Conversion η per pulse (constant volume, best fit) with the profile-vector band, and removal/dose under "
         "the evaporation SCENARIOS")),
        ("fig8b_O2_utilization_kLa_SCENARIO_DO0",
         ("O2 required / O2 transfer capacity per kLa SCENARIO (DO not measured: capacity at DO = 0); hatched: "
         "infeasible")),
        ("fig8c_idle_fraction_bounds_coverage",
         ("Idle fraction of stirred time: between the lower and upper bound, with the CONDITIONAL value over "
         "classified time (coverage ≥ 0.5)")),
        ("fig8d_aeration_energy_per_N_removed_AE_kLa_SCENARIO",
         ("Aeration-equivalent energy per N removed at the central field AE per kLa SCENARIO (whiskers: field AE "
         "range); dashed: the demand-based floor 4.57/AE")),
    ]),
]
FIGURES = [(f"{stem}.png", f"{caption}: {panel}") for _, caption, panels in FIGURE_GROUPS for stem, panel in panels]


# --- Canvas ----------------------------------------------------------------------------------------------------
def _plt():
    from .report import pyplot

    plt = pyplot()
    plt.rcParams.update(PARAMS)
    return plt


def _new(plt):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for side in ("top", "bottom", "left", "right"):
        ax.spines[side].set_visible(True)
    ax.grid(False)
    return fig, ax


def _save(fig, d: Path, stem: str, plt) -> None:
    """Save the figure, then the same figure with every text label removed to ``d / NO_TEXT_DIR`` under the same
    name: the canvas, axes, data, ticks, legend handles and arrows stay where they are."""
    for out, strip in ((d, False), (d / NO_TEXT_DIR, True)):
        if strip:
            _strip_text(fig)
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{stem}.png", dpi=300, transparent=True)
        fig.savefig(out / f"{stem}.svg", dpi=300, transparent=True, metadata={"Date": None})
    plt.close(fig)


def _strip_text(fig) -> None:
    """Hide every text label (axis labels, tick labels, offset text, legend text, annotations, box text). Hidden
    legend text keeps its space, so the legend handles do not move; an annotation keeps its arrow."""
    from matplotlib.text import Annotation, Text

    for ax in fig.axes:
        ax.tick_params(which="both", labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.xaxis.offsetText.set_visible(False)
        ax.yaxis.offsetText.set_visible(False)
    for t in fig.findobj(Text):
        if isinstance(t, Annotation) and t.arrow_patch is not None:
            t.set_text("")
        else:
            t.set_visible(False)


def _pulse_lines(ax, times) -> None:
    for tp in times:
        ax.axvline(tp, color=INK2, lw=1.5, ls=":")


def _legend(ax, **kw) -> None:
    ax.legend(frameon=False, handlelength=1.2, borderaxespad=0.2, **kw)


def _headroom(ax, factor: float) -> None:
    """Extend the y-axis above the data so the legend sits clear of it."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, lo + factor * (hi - lo))


def _metric(pm: pd.DataFrame, metric: str, col: str = "value") -> pd.Series:
    return pm[pm["metric"] == metric].set_index("pulse")[col]


def _read(out_dir: Path, rel: str) -> pd.DataFrame:
    p = Path(out_dir) / rel
    if not p.is_file():
        raise FileNotFoundError(f"{p} missing: run the stage that writes it first")
    return pd.read_parquet(p) if rel.endswith(".parquet") else pd.read_csv(p)


# --- Fig. 1 (W6) -----------------------------------------------------------------------------------------------
def fig1(t: dict[str, pd.DataFrame], d: Path, saes: list[float]) -> None:
    """The linkage diagram on a 3 x 3 grid, every label at the rcParams size: a title and one short line per box; the
    report caption and §1 carry the detail."""
    from matplotlib.patches import Patch

    plt = _plt()
    v = t["sensor_verdicts.csv"].set_index("key")["verdict"].astype(str)
    fr = t["states_fractions.csv"].set_index("interval").loc["whole record"]
    cons = t["nh4_truth_w2_consistency.csv"]
    contradicted = (~cons["consistent_at_best_fit"].astype(bool)).any()
    valid = [k for k in ("ph", "nh4", "orp", "ec") if v[k].startswith("PASS")]
    sensors_line = (" ".join({"ph": "pH", "nh4": "NH4", "orp": "ORP", "ec": "EC"}[k] for k in valid) + " valid"
                    if valid else "none valid")
    # (cell, title, line, border); the chain runs Sensors → C(t) → State ↓ Rates → O2 ↓ Energy → Efficiency
    boxes = {
        "sensors": ((0, 2), "Sensors", sensors_line, BLUE),
        "ct": ((1, 2), "C(t)", "M6 inverse", VIOLET),
        "state": ((2, 2), "State", (f"idle {100 * fr['idle_fraction_lower_bound']:.0f}-"
                                    f"{100 * fr['idle_fraction_upper_bound']:.0f} %"), BLUE),
        "rates": ((2, 1), "Rates", "−dC/dt, η", VIOLET),
        "o2": ((1, 1), "O2 demand", "4.57 g/g N", VIOLET),
        "energy": ((1, 0), "Energy", f"AE {min(saes):g}-{max(saes):g}", ORANGE),
        "eff": ((2, 0), "Efficiency", "kWh/kg N", ORANGE),
        "w2": ((0, 1), "W2 bound", "C(t) fails" if contradicted else "C(t) holds", BLUE),
    }
    fig, ax = _new(plt)
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 3)
    ax.set_xticks([])
    ax.set_yticks([])
    for key, ((i, j), title, line, color) in boxes.items():
        ax.add_patch(plt.Rectangle((i + BOX_X, j + BOX_Y), 1 - 2 * BOX_X, 1 - 2 * BOX_Y, facecolor="none",
                                   edgecolor=color, lw=3,
                                   ls="--" if key == "w2" else "-"))
        ax.text(i + 0.5, j + 0.68, title, ha="center", va="center", color=INK)
        ax.text(i + 0.5, j + 0.30, line, ha="center", va="center", color=INK2)
    for a, b in (("sensors", "ct"), ("ct", "state"), ("state", "rates"), ("rates", "o2"), ("o2", "energy"),
                 ("energy", "eff"), ("sensors", "w2")):
        _arrow(ax, boxes[a][0], boxes[b][0], BLUE if b == "w2" else INK2)
    ax.annotate("", xy=(1.12, 2.12), xytext=(0.88, 1.88), arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 3})
    ax.legend(handles=[Patch(facecolor="none", edgecolor=c, lw=3, label=lab)
                       for c, lab in ((BLUE, "measured"), (VIOLET, "model"), (ORANGE, "scenario"))],
              loc="lower left", frameon=False, handlelength=0.8, handletextpad=0.4, labelspacing=0.2,
              borderaxespad=0.2)
    _save(fig, d, "fig1_linkage_sensors_state_O2_energy_efficiency", plt)


BOX_X, BOX_Y = 0.08, 0.1  # fig1 box inset from its grid cell; the gaps between boxes hold the arrows


def _arrow(ax, a: tuple[int, int], b: tuple[int, int], color: str) -> None:
    """An arrow between two adjacent grid cells, from the edge of box a to the facing edge of box b."""
    (i, j), (i2, j2) = a, b
    if j == j2:
        x0, x1 = (i + 1 - BOX_X, i2 + BOX_X) if i2 > i else (i + BOX_X, i2 + 1 - BOX_X)
        start, end = (x0, j + 0.5), (x1, j + 0.5)
    else:
        start, end = (i + 0.5, j + BOX_Y), (i + 0.5, j2 + 1 - BOX_Y)
    ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "-|>", "color": color, "lw": 3})


# --- Fig. 2 (M3) -----------------------------------------------------------------------------------------------
def segments_figures(cfg, out_dir) -> None:
    from .segments import COMMON_END_D, SEGMENT_START_D, bins_in, step_candidates

    plt = _plt()
    out_dir = Path(out_dir)
    d = out_dir / "figures"
    grid = _read(out_dir, "processed/grid_0p01d.parquet")
    table = _read(out_dir, "tables/segments.csv")
    jumps = _read(out_dir, "tables/nh4_jumps.csv")
    cp = _read(out_dir, "tables/changepoints.csv")
    held = _read(out_dir, "tables/nh4_held_segments.csv")
    ec_path = out_dir / "tables" / "ec_late_pulses.csv"
    ec_pulses = pd.read_csv(ec_path) if ec_path.is_file() else None
    t = grid["t_mid_d"].to_numpy()
    panels = [("fig2a_overview_NH4_slope_windows_jumps_held", "nh4_mgNH4_L", "NH4+ (mg/L)", "NH4"),
              ("fig2b_overview_pH", "pH", "pH", "pH"),
              ("fig2c_overview_ORP_relative", "orp_mV_vs_Ag", "ORP (mV vs Ag)", "ORP"),
              ("fig2d_overview_EC_unvalidated", "EC_mS_cm", "EC (mS/cm)", "EC")]
    for stem, col, lab, series in panels:
        fig, ax = _new(plt)
        ax.axvspan(0, SEGMENT_START_D, color=GRID, alpha=0.6, lw=0)
        ax.plot(t, grid[col], color=BLUE, lw=1.5)
        _pulse_lines(ax, cfg.pulse_times_d)
        ev = cp[(cp["series"] == series) & cp["event"].astype(bool) & cp["tied_pulse"].isna()]
        ax.scatter(ev["t_break_d"], np.interp(ev["t_break_d"], t, grid[col]), marker="D", s=90, color=ORANGE,
                   edgecolor="white", linewidth=1, zorder=3, label="event")
        if series == "NH4":
            for h in held.itertuples(index=False):
                ax.axvspan(h.t_start_d, h.t_end_d, ymin=0.95, ymax=1.0, color=INK2, alpha=0.5, lw=0)
            comp = table[table["status"] == "computed"]
            for k, r in enumerate(comp.itertuples(index=False)):
                b = bins_in(grid, r.fit_start_d, r.fit_end_d)
                mid_t, mid_y = b["t_mid_d"].median(), b["nh4_mgNH4_L"].median()
                tt = np.array([r.fit_start_d, r.fit_end_d])
                ax.plot(tt, mid_y + r.slope_mgNH4_L_d * (tt - mid_t), color=AQUA, lw=4,
                        label="slope" if k == 0 else None)
            located = jumps[jumps["jump"].astype(bool)]
            for cls, color, marker, lab_j in (("pulse_response", AQUA, "^", "jump: pulse"),
                                              ("electrode_offset", VIOLET, "v", "jump: offset")):
                s = located[located["class"] == cls]
                ax.scatter(s["t_jump_d"], s["post_level"], marker=marker, s=150, color=color, edgecolor="white",
                           linewidth=1.2, zorder=4, label=lab_j)
            other = step_candidates(jumps)
            other = other[~other["jump"].astype(bool)]
            ax.scatter(other["t_event_d"], other["post_level"], marker="x", s=110, color=ORANGE, linewidth=2.5,
                       zorder=4, label="step cand.")
            lo, hi = np.nanmin(grid[col]), np.nanmax(grid[col])
            ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.75 * (hi - lo))
            _legend(ax, loc="upper center", ncol=3, columnspacing=0.8, handletextpad=0.3)
        else:
            if series == "EC" and ec_pulses is not None and len(ec_pulses):
                ax.scatter(ec_pulses["t_peak_d"], np.interp(ec_pulses["t_peak_d"], t, grid["EC_mS_cm"]) + 0.05,
                           marker="|", s=300, linewidth=3, color=VIOLET, label="late EC pulse")
            _legend(ax, loc="best")
        ax.set_xlim(0, COMMON_END_D)
        ax.set_xlabel("Time (d)")
        ax.set_ylabel(lab)
        _save(fig, d, stem, plt)


# --- Fig. 3 (M2) -----------------------------------------------------------------------------------------------
def validation_figures(cfg, out_dir) -> None:
    from .validation import load_natives

    plt = _plt()
    out_dir = Path(out_dir)
    d = out_dir / "figures"
    natives = load_natives(out_dir)
    pt = _read(out_dir, "tables/sensor_validation.csv")
    bal = _read(out_dir, "tables/sensor_balances.csv")
    jumps = _read(out_dir, "tables/nh4_jumps.csv")
    held = _read(out_dir, "tables/nh4_held_segments.csv")
    pulses = cfg.pulse_times_d
    bal_v = dict(zip(bal["quantity"], bal["value"]))

    fig, ax = _new(plt)
    c = natives["cond"]
    ax.plot(c["t_d"], c["EC_mS_cm"], color=BLUE, lw=1.5, label="observed")
    t_grid = np.linspace(0, c["t_d"].max(), 4000)
    n_before = np.searchsorted(pulses, t_grid, side="right")
    imm = pt["ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C"].iloc[0]
    nit = pt["ec_pred_change_nitrified_mS_cm_limiting_SCENARIO_T25C"].iloc[0]
    ec0 = float(bal_v["ec_median_day_0.3_1.3_mS_cm"])
    ax.plot(t_grid, ec0 + n_before * imm, color=ORANGE, label="pred., no nitrification")
    ax.plot(t_grid, ec0 + n_before * nit, color=AQUA, label="pred., full nitrification")
    _pulse_lines(ax, pulses)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("EC (mS/cm)")
    _headroom(ax, 1.45)
    _legend(ax, loc="upper left")
    _save(fig, d, "fig3a_validation_EC_predicted_staircase_vs_record", plt)

    fig, ax = _new(plt)
    n = natives["nh4"]
    for _, h in held.iterrows():
        ax.axvspan(h["t_start_d"], h["t_end_d"], color=GRID, alpha=0.9, lw=0)
    ax.plot(n["t_d"], n["nh4_mgNH4_L"], color=BLUE, lw=1.5, label="observed")
    j = jumps[jumps["jump"].astype(bool)]
    for cls, color, marker, lab in (("pulse_response", AQUA, "^", "jump: pulse"),
                                    ("electrode_offset", VIOLET, "v", "jump: offset")):
        s = j[j["class"] == cls]
        ax.scatter(s["t_jump_d"], s["post_level"], color=color, marker=marker, s=150, zorder=3, edgecolor="white",
                   linewidth=1.5, label=lab)
    _pulse_lines(ax, pulses)
    ax.plot([], [], color=GRID, lw=10, label="held")
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("NH4+ (mg/L)")
    _headroom(ax, 1.6)
    _legend(ax, loc="upper right", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig3b_validation_NH4_jumps_held_segments", plt)

    k = pt["pulse"].to_numpy()
    x = np.arange(len(k))
    wbar = 0.38
    panels = [
        ("fig3c_validation_NH4_change_per_pulse", "ΔNH4+ (mg/L)", pt["nh4_pred_change_mgNH4_L"],
         pt["nh4_obs_change_mgNH4_L"], None),
        ("fig3d_validation_EC_change_per_pulse", "ΔEC (mS/cm)",
         pt["ec_pred_change_immediate_mS_cm_limiting_SCENARIO_T25C"], pt["ec_obs_change_mS_cm"], None),
        ("fig3e_validation_ORP_change_per_pulse", "ΔORP (mV)", None, pt["orp_obs_change_mV"],
         pt[[col for col in pt.columns if col.startswith("orp_pred_change_mV_")]]),
        ("fig3f_validation_pH_change_per_pulse", "ΔpH", None, pt["ph_obs_change"], None),
    ]
    for stem, lab, pred, obs, scen in panels:
        fig, ax = _new(plt)
        if pred is not None:
            ax.bar(x - wbar / 2 - 0.01, pred, wbar, color=ORANGE, label="predicted")
        if scen is not None:
            lo_, hi_ = scen.min(axis=1), scen.max(axis=1)
            ax.bar(x - wbar / 2 - 0.01, hi_ - lo_, wbar, bottom=lo_, color=ORANGE, alpha=0.45, label="pred. range")
        ax.bar(x + wbar / 2 + 0.01, obs, wbar, color=BLUE, label="observed")
        ax.axhline(0, color=INK2, lw=1.5)
        ax.set_xticks(x, [f"P{v}" for v in k])
        ax.set_xlabel("Pulse")
        ax.set_ylabel(lab)
        _headroom(ax, 1.4)
        _legend(ax, loc="upper right", ncol=2, columnspacing=0.8)
        _save(fig, d, stem, plt)


# --- Fig. 5 (M5) -----------------------------------------------------------------------------------------------
def dynamics_figures(cfg, out_dir) -> None:
    from .dynamics import END_D, ONSAGER_MIN_SNR

    plt = _plt()
    out_dir = Path(out_dir)
    d = out_dir / "figures"
    v = _read(out_dir, "tables/estimator_validation.csv")
    w = _read(out_dir, "tables/dynamics_windows.csv")
    s = _read(out_dir, "tables/dynamics_onsager.csv")
    flux = _read(out_dir, "tables/dynamics_flux.csv")

    fig, ax = _new(plt)
    clean = v[v["check"] == "markov_ring_entropy_production"]
    true = clean["expected"].astype(float).to_numpy()
    val = clean["value"].astype(float)
    ax.errorbar(true, val, yerr=[val - clean["ci_low"].astype(float), clean["ci_high"].astype(float) - val],
                fmt="o", color=BLUE, ms=10, capsize=5, elinewidth=2.5, label="estimator")
    q = v[v["check"] == "quantized_observable_bias"]
    q_true = true[1:]
    ax.scatter(q_true, q_true + q["value"].astype(float).to_numpy(), marker="s", color=ORANGE, s=110, zorder=3,
               label="quantized 0.01")
    ax.plot([0, 1], [0, 1], color=INK2, lw=2, ls="--")
    ax.set_xlabel("Analytic (p − q) ln(p/q)")
    ax.set_ylabel("Estimated EP/step")
    _legend(ax, loc="upper left")
    _save(fig, d, "fig5a_estimator_validation_ring_chain", plt)

    fig, ax = _new(plt)
    status = w["ordinal_kl_status"] if "ordinal_kl_status" in w else pd.Series("", index=w.index)
    for series, color in (("pH", BLUE), ("ORP", AQUA), ("NH4", VIOLET)):
        g = w[(w["series"] == series) & (status == "computed")]
        mid = 0.5 * (g["t_start_d"] + g["t_end_d"])
        sig = g["ordinal_kl_p"] <= 0.05
        ax.scatter(mid[sig], g.loc[sig, "ordinal_kl"], color=color, s=90, label=series)
        ax.scatter(mid[~sig], g.loc[~sig, "ordinal_kl"], facecolor="none", edgecolor=color, s=90, linewidth=2)
    _pulse_lines(ax, cfg.pulse_times_d)
    ax.set_xlim(0, END_D)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("Ordinal KL (nats)")
    _legend(ax, loc="upper right")
    _save(fig, d, "fig5b_ordinal_KL_irreversibility_stationary_windows", plt)

    fig, ax = _new(plt)
    ax.bar(s["pulse"].astype(str), s["snr_pulse_response"], color=BLUE, width=0.55)
    ax.axhline(ONSAGER_MIN_SNR, color=ORANGE, ls="--", label=f"required ({ONSAGER_MIN_SNR:g})")
    ax.set_xlabel("Pulse")
    ax.set_ylabel("Pulse SNR")
    ax.set_ylim(0, 1.35 * max(ONSAGER_MIN_SNR, float(s["snr_pulse_response"].max())))
    _legend(ax, loc="upper right")
    _save(fig, d, "fig5c_Onsager_SNR_per_pulse", plt)

    fig, ax = _new(plt)
    f = flux[pd.to_numeric(flux["circulation"], errors="coerce").notna()].copy()
    f["circulation"] = f["circulation"].astype(float)
    for i, (plane, color) in enumerate((("pH-ORP", AQUA), ("pH-NH4", VIOLET))):
        g = f[f["plane"] == plane]
        xpos = g["segment"].to_numpy() + (i - 0.5) * 0.25
        ax.errorbar(xpos, g["circulation"], yerr=[g["circulation"] - g["ci_low"], g["ci_high"] - g["circulation"]],
                    fmt="o", color=color, ms=10, capsize=5, elinewidth=2.5, label=plane)
    ax.axhline(0, color=INK2, lw=1.5)
    ax.set_xticks(range(5), ["base", "P1", "P2", "P3", "P4"])
    ax.set_xlim(-0.5, 4.5)
    ax.set_xlabel("Segment")
    ax.set_ylabel("Circulation (1/step)")
    _legend(ax, loc="best")
    _save(fig, d, "fig5d_probability_flux_circulation", plt)


# --- Fig. 6 (M6) -----------------------------------------------------------------------------------------------
def model_figures(cfg, out_dir) -> None:
    from .model import CHI2_95_1DOF, split_windows

    plt = _plt()
    out_dir = Path(out_dir)
    d = out_dir / "figures"
    s = _read(out_dir, "processed/model_state.parquet")
    idf = _read(out_dir, "tables/model_identifiability.csv")
    sw = split_windows(cfg)

    fig, ax = _new(plt)
    for (a, b), lab in ((sw.sel, "selection"), (sw.test, "test")):
        ax.axvspan(a, b, color=GRID, alpha=0.6 if lab == "selection" else 0.9, lw=0)
    ax.plot(s["t_d"], s["electrode_obs_mgN_L"], color=BLUE, lw=1.5, label="electrode")
    ax.plot(s["t_d"], s["model_electrode_calibration_fit_mgN_L"], color=ORANGE, label="calibrated, frozen")
    ax.plot(s["t_d"], s["model_electrode_wholefit_mgN_L"], color=AQUA, ls="--", label="whole-record refit")
    _pulse_lines(ax, cfg.pulse_times_d)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("NH4+ (mg N/L)")
    _headroom(ax, 1.5)
    _legend(ax, loc="upper left")
    _save(fig, d, "fig6a_model_selected_vs_electrode_cal_sel_test", plt)

    fig, ax = _new(plt)
    ax.plot(s["t_d"], s["model_rate_mgN_L_d"], color=VIOLET)
    _pulse_lines(ax, cfg.pulse_times_d)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("Rate (mg N/L/d)")
    _save(fig, d, "fig6b_model_consumption_rate_whole_record_refit", plt)

    fig, ax = _new(plt)
    x = np.arange(len(idf))
    ok = idf["identifiable"].astype(str).str.lower().eq("true").to_numpy()
    lo = idf["profile_max_delta_low"].astype(float).clip(upper=40).to_numpy()
    hi = idf["profile_max_delta_high"].astype(float).clip(upper=40).to_numpy()
    for vals, off, marker, color, lab in ((lo, -0.15, "v", BLUE, "low side"), (hi, 0.15, "^", ORANGE, "high side")):
        ax.scatter(x[ok] + off, vals[ok], marker=marker, color=color, s=120, label=lab)
        ax.scatter(x[~ok] + off, vals[~ok], marker=marker, facecolor="none", edgecolor=color, s=120, linewidth=2)
    ax.axhline(CHI2_95_1DOF, color=INK2, lw=2, ls="--")
    ax.set_xticks(x, _short_params(idf["parameter"]), rotation=90)
    ax.set_ylim(0, 48)
    ax.set_xlabel("Parameter")
    ax.set_ylabel("Max profile Δ")
    _legend(ax, loc="upper right", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig6c_profile_likelihood_identifiability", plt)


def _short_params(names: pd.Series) -> list[str]:
    """Tick labels that fit at 25 pt: the electrode offsets become o1, o2, ... in time order."""
    short = {"v0_mgN_L_d": "v0", "kd_per_d": "kd", "C0_mgN_L": "C0", "b1_per_d": "b1", "tau_d": "τ"}
    out, k = [], 0
    for n in names:
        if str(n).startswith("offset_"):
            k += 1
            out.append(f"o{k}")
        else:
            out.append(short.get(n, str(n)))
    return out


# --- Fig. 7 (M7) -----------------------------------------------------------------------------------------------
def pulses_figures(cfg, out_dir) -> None:
    from . import thermo
    from .constants import DG0_NITRIFICATION_KJ_MOL_N
    from .pulses import RESPONSE, wide

    plt = _plt()
    out_dir = Path(out_dir)
    d = out_dir / "figures"
    metrics = _read(out_dir, "tables/pulses.csv")
    w = wide(metrics)
    x = np.arange(len(w))
    labels = [f"P{int(p)}" for p in w["pulse"]]

    fig, ax = _new(plt)
    lo, hi = wide(metrics, "ci95_low"), wide(metrics, "ci95_high")
    for off, col, color, lab in ((-0.2, "consumption_uncorrected_weighted_mgN_L_d", BLUE, "electrode"),
                                 (0.2, RESPONSE, ORANGE, "M6-inverted")):
        v = w[col].to_numpy(dtype=float)
        err = np.vstack([v - lo[col].to_numpy(dtype=float), hi[col].to_numpy(dtype=float) - v])
        ax.bar(x + off, v, 0.38, color=color, label=lab, yerr=err,
               error_kw={"ecolor": INK2, "capsize": 5, "lw": 2})
    ax.scatter(x + 0.2, w["consumption_SENSITIVITY_linear_gain_weighted_mgN_L_d"], marker="_", s=600, linewidth=4,
               color=INK, zorder=3, label="linear gain")
    ax.axhline(0, color=INK2, lw=1.5)
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel("Rate (mg N/L/d)")
    _legend(ax, loc="best")
    _save(fig, d, "fig7a_consumption_after_each_pulse_95CI", plt)

    fig, ax = _new(plt)
    da = w["delta_A_kJ_mol_SCENARIO_T25C"].to_numpy(dtype=float)
    ax.bar(x, da, 0.55, color=AQUA)
    for xi, v in zip(x, da):
        ax.text(xi, v, f"{100 * v / abs(DG0_NITRIFICATION_KJ_MOL_N):.1f} %", ha="center", va="bottom")
    ax.set_ylim(0, 1.25 * np.nanmax(da))
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel("ΔA (kJ/mol)")
    _save(fig, d, "fig7b_affinity_step_M6_inverted_NH4_25C", plt)

    fig, ax = _new(plt)
    fa_cols = [c for c in w.columns if c.startswith("fa_median_mgNH3_L_SCENARIO_")]
    for col, color in zip(fa_cols, (BLUE, VIOLET, AQUA, ORANGE)):
        lab = col.removeprefix("fa_median_mgNH3_L_SCENARIO_T").removesuffix("C") + " °C"
        ax.scatter(w[col], w[RESPONSE], color=color, s=110, label=lab)
    ax.axhline(0, color=INK2, lw=1.5)
    ax.set_xlabel("Free ammonia (mg NH3/L)")
    ax.set_ylabel("Rate (mg N/L/d)")
    _legend(ax, loc="best")
    _save(fig, d, "fig7c_free_ammonia_vs_consumption_EXPLORATORY", plt)

    fig, ax = _new(plt)
    eta = w["eta_conversion_fraction"].to_numpy(dtype=float) if "eta_conversion_fraction" in w else np.full(len(w), np.nan)
    mb = w["eta_conversion_fraction_MODEL_BASED"].to_numpy(dtype=float)
    ax.bar(x - 0.2, np.nan_to_num(eta), 0.38, color=ORANGE, label="dose-anchored")
    for xi, e in zip(x, eta):
        if np.isnan(e):
            ax.text(xi - 0.2, 0.01, "n.c.", ha="center", va="bottom", color=INK2)
    ax.bar(x + 0.2, mb, 0.38, color=VIOLET, label="MODEL-BASED")
    ax.set_ylim(0.0, 1.6 * np.nanmax(np.r_[eta, mb]))
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel(f"η (of {thermo.pulse_dose_from_config(cfg).n_mg_L:.0f} mg N/L)")
    _legend(ax, loc="upper left")
    _save(fig, d, "fig7d_conversion_eta_before_next_pulse", plt)


# --- Fig. 4 (W6) -----------------------------------------------------------------------------------------------
def fig4(t: dict[str, pd.DataFrame], states_bins: pd.DataFrame, pulse_times: list[float], d: Path) -> None:
    from . import alkalinity, volume

    plt = _plt()
    fig, ax = _new(plt)
    b = t["alkalinity_bound.csv"]
    base = b[(b["volume_scenario"] == volume.CONSTANT) & (b["alk0_basis"] != "MEASURED")]
    for T, ls in ((25.0, "-"), (20.0, ":")):
        g = base[(base["co2_scenario"] == alkalinity.AIR_EQ) & (base["temperature_C"] == T)]
        ax.plot(g["alk0_meq_L"], g["n_nit_max_end_mgN_L"], color=BLUE, ls=ls, label=f"W2 bound, {T:g} °C")
    any_s = base[(base["co2_scenario"] == alkalinity.AIR_EQ) & (base["temperature_C"] == 25.0)]
    ax.plot(any_s["alk0_meq_L"], any_s["n_removed_max_any_proton_sink_end_mgN_L"], color=VIOLET, ls="--",
            label="any proton sink")
    dosed = float(b["dosed_total_mgN_L"].iloc[0])
    ax.axhline(dosed, color=INK2, lw=2, ls="-.", label="dosed")
    cons = t["nh4_truth_w2_consistency.csv"].set_index("checkpoint").loc["end"]
    ax.axhspan(cons["removed_band_low_mgN_L"], cons["removed_band_high_mgN_L"], color=ORANGE, alpha=0.15, lw=0)
    ax.axhline(cons["removed_per_initial_litre_mgN_L"], color=ORANGE, label="C(t) removal")
    ax.set_xlabel("Alk0 (meq/L, scanned)")
    ax.set_ylabel("N removed (mg N/L)")
    _headroom(ax, 1.6)
    _legend(ax, loc="upper center", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig4a_alkalinity_bound_vs_Ct_removal_day47p9", plt)

    fig, ax = _new(plt)
    tw = t["nh4_truth_windows.csv"]
    tw = tw[tw["status"] == "computed"]
    fa_col = "fa_median_mgNH3_L_SCENARIO_T25C"
    for k, color in zip((1, 2, 3, 4), (BLUE, ORANGE, AQUA, VIOLET)):
        g = tw[tw["governing_pulse"] == k]
        if not len(g):
            continue
        err = np.vstack([g["consumption_mgN_L_d"] - g["consumption_ci95_low_mgN_L_d"],
                         g["consumption_ci95_high_mgN_L_d"] - g["consumption_mgN_L_d"]])
        ax.errorbar(g[fa_col], g["consumption_mgN_L_d"], yerr=err, fmt="o", ms=12, color=color, ecolor=color,
                    elinewidth=2, capsize=4, mec="white", mew=1.5, label=f"P{k}")
    ax.set_xlabel("Free ammonia (mg NH3/L)")
    ax.set_ylabel("−dC/dt (mg N/L/d)")
    _legend(ax, loc="best", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig4b_rate_vs_free_ammonia_W1_slope_windows", plt)

    fig, ax = _new(plt)
    tl = states_bins
    for i, (state, color) in enumerate(STATE_COLORS.items()):
        m = (tl["state"] == state).to_numpy()
        ax.broken_barh([(x - 0.005, 0.01) for x in tl.loc[m, "t_mid_d"]], (len(STATE_COLORS) - 1 - i - 0.4, 0.8),
                       facecolors=color, linewidth=0)
    _pulse_lines(ax, pulse_times)
    ax.set_yticks(range(len(STATE_COLORS)), [STATE_LABELS[s] for s in reversed(STATE_COLORS)])
    ax.set_ylim(-0.6, len(STATE_COLORS) - 0.4)
    ax.set_xlabel("Time (d)")
    ax.set_ylabel("State")
    _save(fig, d, "fig4c_process_state_timeline_W3", plt)


# --- Fig. 8 (W6) -----------------------------------------------------------------------------------------------
def fig8(t: dict[str, pd.DataFrame], d: Path) -> None:
    plt = _plt()
    pm = t["pulses.csv"]
    pulses = sorted(pm["pulse"].unique())
    x = np.arange(len(pulses))
    labels = [f"P{int(p)}" for p in pulses]

    fig, ax = _new(plt)
    eta = _metric(pm, "eta_conversion_fraction").reindex(pulses)
    lo = _metric(pm, "eta_conversion_fraction", "param_band_low").reindex(pulses)
    hi = _metric(pm, "eta_conversion_fraction", "param_band_high").reindex(pulses)
    ax.bar(x, np.nan_to_num(eta.to_numpy(dtype=float)), 0.5, color=ORANGE, label="η")
    ax.vlines(x, lo, hi, color=INK, lw=3, label="band")
    for metric, marker, lab in (("removal_to_dose_ratio_SCENARIO_evaporation_0.5_mL_per_d", "^", "0.5 mL/d"),
                                ("removal_to_dose_ratio_SCENARIO_evaporation_1_mL_per_d", "v", "1 mL/d")):
        ax.plot(x + 0.33, _metric(pm, metric).reindex(pulses), marker, ms=13, color=VIOLET, ls="none", label=lab)
    for xi, e in zip(x, eta):
        if np.isnan(e):
            ax.text(xi, 0.03, "n.c.", ha="center", color=INK2)
    ax.axhline(1.0, color=INK2, lw=2, ls="--")
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel("Removed / dose")
    ax.set_ylim(0, 2.0)
    _legend(ax, loc="upper left", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig8a_conversion_eta_per_pulse_evaporation_SCENARIO", plt)

    fig, ax = _new(plt)
    klas = sorted({float(m.split("_kLa")[1].split("_")[0]) for m in pm["metric"]
                   if m.startswith("o2_utilization_SCENARIO_kLa")})
    for j, (kla, color) in enumerate(zip(klas, (BLUE, AQUA, VIOLET))):
        u = _metric(pm, f"o2_utilization_SCENARIO_kLa{kla:g}_T25C").reindex(pulses).to_numpy(dtype=float)
        st = _metric(pm, f"o2_utilization_SCENARIO_kLa{kla:g}_T25C", "status").reindex(pulses).astype(str)
        xx = x + (j - 1) * 0.26
        ax.bar(xx, np.nan_to_num(u), 0.24, color=color, label=f"kLa {kla:g}")
        inf = st.str.startswith("not computed: SCENARIO infeasible").to_numpy()
        ax.bar(xx[inf], np.full(inf.sum(), 5.0), 0.24, facecolor="none", edgecolor=color, hatch="//", lw=2)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 50)
    _not_computed(ax, x, [_metric(pm, f"o2_utilization_SCENARIO_kLa{k:g}_T25C").reindex(pulses) for k in klas])
    ax.axhline(1.0, color=INK2, lw=2, ls="--")
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel("O2 utilization")
    _legend(ax, loc="upper center", ncol=3, columnspacing=0.8)
    _save(fig, d, "fig8b_O2_utilization_kLa_SCENARIO_DO0", plt)

    fig, ax = _new(plt)
    lo_i = _metric(pm, "idle_fraction_LOWER_BOUND").reindex(pulses).to_numpy(dtype=float)
    hi_i = _metric(pm, "idle_fraction_UPPER_BOUND").reindex(pulses).to_numpy(dtype=float)
    cond = _metric(pm, "idle_fraction_of_classified_time_CONDITIONAL").reindex(pulses).to_numpy(dtype=float)
    ax.bar(x, np.nan_to_num(hi_i - lo_i), 0.42, bottom=np.nan_to_num(lo_i), color=ORANGE, alpha=0.35, label="bounds")
    ax.hlines(lo_i, x - 0.21, x + 0.21, color=ORANGE, lw=4, label="lower")
    if np.isfinite(cond).any():  # only at a coverage >= 0.5; no legend entry without a marker
        ax.plot(x, cond, "D", ms=12, color=BLUE, ls="none", label="conditional")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.5)
    ax.set_yticks(np.arange(0, 1.01, 0.25))
    ax.set_xlabel("Pulse")
    ax.set_ylabel("Idle fraction")
    _legend(ax, loc="upper center", ncol=3, columnspacing=0.8)
    _save(fig, d, "fig8c_idle_fraction_bounds_coverage", plt)

    fig, ax = _new(plt)
    saes = sorted({float(m.split("_AE")[1].split("_")[0]) for m in pm["metric"]
                   if m.startswith("aeration_energy_supply_kWh_per_kgN_SCENARIO_AE")})
    mid = saes[len(saes) // 2]
    for j, (kla, color) in enumerate(zip(klas, (BLUE, AQUA, VIOLET))):
        vals = {s: _metric(pm, f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{s:g}_kLa{kla:g}_T25C").reindex(pulses)
                .to_numpy(dtype=float) for s in saes}
        xx = x + (j - 1) * 0.26
        ax.bar(xx, np.nan_to_num(vals[mid]), 0.24, color=color, label=f"kLa {kla:g}")
        ax.vlines(xx, vals[max(saes)], vals[min(saes)], color=INK, lw=2.5)
    floor = _metric(pm, f"aeration_energy_demand_kWh_per_kgN_SCENARIO_AE{mid:g}").dropna()
    if len(floor):
        ax.axhline(float(floor.iloc[0]), color=ORANGE, lw=3, ls="--", label="floor")
    ax.set_yscale("log")
    ax.set_ylim(top=ax.get_ylim()[1] * 30)
    _not_computed(ax, x, [_metric(pm, f"aeration_energy_supply_kWh_per_kgN_SCENARIO_AE{mid:g}_kLa{k:g}_T25C")
                          .reindex(pulses) for k in klas], y=0.15)  # clear of the floor line
    ax.set_xticks(x, labels)
    ax.set_xlabel("Pulse")
    ax.set_ylabel("kWh/kg N")
    _legend(ax, loc="upper center", ncol=2, columnspacing=0.8)
    _save(fig, d, "fig8d_aeration_energy_per_N_removed_AE_kLa_SCENARIO", plt)


def _not_computed(ax, x, series: list[pd.Series], y: float = 0.03) -> None:
    """Write "n.c." at a pulse where no series has a value (reason in ``pulses.csv``), at height ``y`` (axes
    fraction)."""
    empty = np.all([s.isna().to_numpy() for s in series], axis=0)
    for xi in np.asarray(x)[empty]:
        ax.text(xi, y, "n.c.", ha="center", va="bottom", color=INK2, transform=ax.get_xaxis_transform())


def make_w6(cfg, t: dict[str, pd.DataFrame], out_dir: Path) -> None:
    """Figures 1, 4 and 8, drawn by the report stage from its tables."""
    d = Path(out_dir) / "figures"
    states_bins = pd.read_parquet(Path(out_dir) / "processed" / "states_bins.parquet")
    fig1(t, d, [float(x) for x in cfg.get("scenarios.sae_field_kgO2_per_kWh")])
    fig4(t, states_bins, cfg.pulse_times_d, d)
    fig8(t, d)


def make_all(cfg, out_dir) -> list[Path]:
    """Every figure from the saved outputs; no stage is rerun. Returns the PNG paths written."""
    out_dir = Path(out_dir)
    validation_figures(cfg, out_dir)
    segments_figures(cfg, out_dir)
    dynamics_figures(cfg, out_dir)
    model_figures(cfg, out_dir)
    pulses_figures(cfg, out_dir)
    names = ("sensor_verdicts.csv", "states_fractions.csv", "nh4_truth_w2_consistency.csv", "alkalinity_bound.csv",
             "nh4_truth_windows.csv", "pulses.csv")
    make_w6(cfg, {n: _read(out_dir, f"tables/{n}") for n in names}, out_dir)
    return [out_dir / "figures" / png for png, _ in FIGURES]
