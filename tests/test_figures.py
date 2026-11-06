"""The plotting rules hold for every figure: one canvas and one axes, a 10 x 6 inch canvas, all four spines, no title,
no tight layout, no grid, the project rcParams, and 300 dpi PNG/SVG on a transparent background."""

import shutil

import numpy as np
import pytest
from PIL import Image

from bnr_neq import figures as F


@pytest.fixture(scope="module")
def drawn(cfg, run_all_dir, tmp_path_factory):
    """Redraw every figure from a copy of the run-all tables, keeping each figure object as it is saved."""
    out = tmp_path_factory.mktemp("figures")
    for sub in ("tables", "processed"):
        shutil.copytree(run_all_dir / sub, out / sub)
    seen, save = {}, F._save

    def capture(fig, d, stem, plt):
        seen[stem] = {"n_axes": len(fig.axes), "size": tuple(fig.get_size_inches()),
                      "layout": fig.get_layout_engine(), "suptitle": fig._suptitle,
                      "spines": [s.get_visible() for s in fig.axes[0].spines.values()],
                      "titles": [fig.axes[0].get_title(loc) for loc in ("left", "center", "right")],
                      "grid": any(g.get_visible() for g in fig.axes[0].get_xgridlines() + fig.axes[0].get_ygridlines()),
                      "rc": {k: plt.rcParams[k] for k in F.PARAMS}}
        save(fig, d, stem, plt)

    F._save = capture
    try:
        F.make_all(cfg, out)
    finally:
        F._save = save
    return out, seen


def test_every_figure_follows_the_plotting_rules(drawn):
    _, seen = drawn
    assert set(seen) == {png[:-4] for png, _ in F.FIGURES}
    for stem, s in seen.items():
        assert s["n_axes"] == 1, stem
        assert s["size"] == F.FIGSIZE == (10, 6), stem
        assert s["layout"] is None and s["suptitle"] is None, stem
        assert len(s["spines"]) == 4 and all(s["spines"]), stem
        assert s["titles"] == ["", "", ""], stem
        assert not s["grid"], stem
        # matplotlib stores the font family names as lists and the sizes as floats
        assert s["rc"] == {k: [v] if k.startswith("font.") and isinstance(v, str) else v
                           for k, v in F.PARAMS.items()}, stem


def test_saved_files_are_300_dpi_and_transparent(drawn):
    out, seen = drawn
    for stem in seen:
        for d in (out / "figures", out / "figures" / F.NO_TEXT_DIR):
            assert (d / f"{stem}.svg").is_file(), stem
            im = Image.open(d / f"{stem}.png")
            assert im.size == (3000, 1800), stem  # 10 x 6 in at 300 dpi
            assert round(im.info["dpi"][0]) == 300, stem
            assert im.mode == "RGBA" and np.asarray(im)[0, 0, 3] == 0, stem  # transparent background


# --- P1-1 (pass 6): no Fig. 1 box whose number carries a SCENARIO or ASSUMPTION label is drawn as "measured" ----
FIG1_TABLES = ("sensor_verdicts.csv", "states_fractions.csv", "nh4_truth_w2_consistency.csv")


@pytest.fixture(scope="module")
def fig1_inputs(cfg, run_all_dir):
    t = {n: F._read(run_all_dir, f"tables/{n}") for n in FIG1_TABLES}
    return t, [float(x) for x in cfg.get("scenarios.sae_field_kgO2_per_kWh")]


def _labelled(t, key) -> bool:
    """True when the rows a box's number is read from carry a SCENARIO or ASSUMPTION label (read independently of
    figures.check_fig1_bases)."""
    name, row = F.FIG1_BOX_SOURCES[key]
    tab = t[name]
    if "assumptions" not in tab.columns:
        return False
    if row is not None:
        tab = tab[tab[row[0]].astype(str) == row[1]]
        assert len(tab) == 1, (key, row)
    return any(("SCENARIO" in s or "ASSUMPTION" in s) for s in tab["assumptions"].fillna("").astype(str))


def test_fig1_labelled_boxes_are_not_drawn_as_measured(fig1_inputs):
    t, saes = fig1_inputs
    boxes = F.fig1_boxes(t, saes)
    assert set(F.FIG1_BOX_SOURCES) <= set(boxes)
    # the two boxes the pass-6 re-review found drawn as measured rest on labelled rows in the run-all tables
    assert _labelled(t, "state") and _labelled(t, "w2")
    for key in F.FIG1_BOX_SOURCES:
        if _labelled(t, key):
            assert boxes[key][3] != "measured", key
    assert boxes["state"][3] == "scenario" and boxes["w2"][3] == "scenario"
    # energy and efficiency rest on the field AE SCENARIO range; only the sensor box is measured
    assert {k for k, b in boxes.items() if b[3] == "measured"} == {"sensors"}
    F.check_fig1_bases(t, boxes)  # the production guard accepts the real boxes


@pytest.mark.parametrize("key", ["state", "w2"])
def test_fig1_refuses_a_labelled_box_drawn_as_measured(fig1_inputs, key):
    t, saes = fig1_inputs
    boxes = F.fig1_boxes(t, saes)
    cell, title, line, _ = boxes[key]
    boxes[key] = (cell, title, line, "measured")
    with pytest.raises(ValueError, match="drawn as measured"):
        F.check_fig1_bases(t, boxes)


def test_fig1_refuses_an_unknown_basis(fig1_inputs):
    t, saes = fig1_inputs
    boxes = F.fig1_boxes(t, saes)
    cell, title, line, _ = boxes["ct"]
    boxes["ct"] = (cell, title, line, "validated")
    with pytest.raises(ValueError, match="unknown basis"):
        F.check_fig1_bases(t, boxes)


def test_drawn_fig1_borders_and_legend_follow_the_basis(fig1_inputs, tmp_path, monkeypatch):
    t, saes = fig1_inputs
    boxes = F.fig1_boxes(t, saes)
    seen = {}

    def capture(fig, d, stem, plt):
        ax = fig.axes[0]
        seen["edges"] = [p.get_edgecolor() for p in ax.patches]
        leg = ax.get_legend()
        seen["legend"] = {txt.get_text(): h.get_edgecolor()
                          for txt, h in zip(leg.get_texts(), leg.legend_handles, strict=True)}
        plt.close(fig)

    monkeypatch.setattr(F, "_save", capture)
    F.fig1(t, tmp_path, saes)
    from matplotlib.colors import to_rgba

    assert seen["edges"] == [to_rgba(F.FIG1_BASIS_COLORS[b[3]]) for b in boxes.values()]
    assert seen["legend"] == {lab: to_rgba(c) for lab, c in F.FIG1_BASIS_COLORS.items()}
    measured = to_rgba(F.FIG1_BASIS_COLORS["measured"])
    for key, edge in zip(boxes, seen["edges"], strict=True):
        if key in F.FIG1_BOX_SOURCES and _labelled(t, key):
            assert edge != measured, key
