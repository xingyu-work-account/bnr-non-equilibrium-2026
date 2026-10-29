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
