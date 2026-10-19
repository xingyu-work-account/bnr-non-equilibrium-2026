import numpy as np
import pandas as pd
import pytest

from bnr_neq import timebase


@pytest.mark.parametrize("key, lo, hi", [("nh4", 22, 450), ("orp", 21, 87)])
def test_stamp_size_ranges(ingested, key, lo, hi):
    _, sizes = timebase.stamp_groups(ingested.natives[key]["t_stamp_d"].to_numpy())
    assert sizes.min() == lo and sizes.max() == hi


@pytest.mark.parametrize("key", ["nh4", "orp"])
def test_reconstructed_times_strictly_increasing(ingested, key):
    df = ingested.natives[key]
    assert (np.diff(df["t_d"].to_numpy()) > 0).all()
    # each reconstructed time stays inside its own stamp
    off = df["t_d"] - df["t_stamp_d"]
    assert (off >= 0).all() and (off < timebase.STAMP_D).all()


@pytest.mark.parametrize("key", ["cond", "ph"])
def test_full_precision_times_strictly_increasing(ingested, key):
    assert (np.diff(ingested.natives[key]["t_d"].to_numpy()) > 0).all()


def test_reconstruction_formula():
    t, j, n = timebase.reconstruct_stamps(np.array([0.0, 0.0, 0.01, 0.01, 0.01, 0.01]), "x")
    np.testing.assert_allclose(t, [0.0, 0.005, 0.01, 0.0125, 0.015, 0.0175])
    assert j.tolist() == [0, 1, 0, 1, 2, 3] and n.tolist() == [2, 2, 4, 4, 4, 4]


def test_stamp_spacing_violation_raises():
    with pytest.raises(timebase.StampError, match="differ by"):
        timebase.reconstruct_stamps(np.array([0.0, 0.0, 0.02, 0.02]), "NH4.csv")


def test_gap_inventory_synthetic():
    t = np.arange(0, 2000) * 10.0 / 86400
    t = np.r_[t[:1000], t[1001:]]  # one missing sample -> 20 s, exactly 2x: missing_sample
    t = np.r_[t[:1500], t[1503:]]  # three missing -> 40 s: gap
    inv = timebase.gap_inventory(t, "s")
    assert inv["class"].tolist() == ["missing_sample", "gap"]
    assert inv["interval_max_s"].round(6).tolist() == [20.0, 40.0]


def test_rate_change_is_not_a_gap():
    t = np.r_[np.arange(0, 500) * 20.0, 500 * 20.0 + np.arange(0, 5000) * 2.0] / 86400
    assert timebase.gap_inventory(t, "s").empty


def test_ph_missing_samples_listed(ingested):
    ph = ingested.gaps[ingested.gaps["sensor"] == "ph"]
    assert len(ph) == 1 and ph["n_intervals"].iloc[0] == 2
    assert ph["t_start_d"].iloc[0] == pytest.approx(4.398, abs=1e-3)
    assert (ph["class"] == "missing_sample").all()


def test_only_partial_stamp_artifacts_for_stamped_sensors(ingested):
    g = ingested.gaps[ingested.gaps["sensor"].isin(["nh4", "orp"])]
    assert g["partial_stamp_artifact"].all()


def test_time_weighting():
    # one bin: one sample representing 20 s, then ten representing 2 s each
    t = np.r_[0.0, 20.0 + np.arange(10) * 2.0, 40.0] / 86400
    x = np.r_[1.0, np.full(10, 3.0), 99.0]
    w = np.r_[20.0, np.full(10, 2.0), 2.0] / 86400
    mean, count = timebase.bin_time_weighted(t, x, w, np.array([0.0, 40.0 / 86400 - 1e-12]))
    assert count[0] == 11
    assert mean[0] == pytest.approx((20 * 1 + 20 * 3) / 40)  # an unweighted mean would be 2.82


def test_empty_bins_are_nan_not_filled():
    edges = np.array([0.0, 0.01, 0.02, 0.03])
    t = np.array([0.001, 0.002, 0.021, 0.025])
    x = np.array([1.0, 2.0, 5.0, 6.0])
    mean, count = timebase.bin_time_weighted(t, x, np.ones(4), edges)
    assert count.tolist() == [2, 0, 2]
    assert np.isnan(mean[1]) and mean[0] == 1.5 and mean[2] == 5.5


def test_grid_has_no_fabricated_values(ingested):
    grid = ingested.grid
    edges = timebase.grid_edges()
    assert len(grid) == 4792 and grid["t_end_d"].iloc[-1] == 47.92
    for sensor, cols in [("cond", ["EC_mS_cm"]), ("nh4", ["nh4_mgNH4_L"]), ("orp", ["orp_mV_vs_Ag"]),
                         ("ph", ["pH"])]:
        df = ingested.natives[sensor]
        idx = timebase.bin_index(df["t_d"].to_numpy(), edges)
        inside = idx >= 0
        for col in cols:
            stats = pd.DataFrame({"b": idx[inside], "x": df[col].to_numpy()[inside]}).groupby("b")["x"].agg(
                ["min", "max", "size"]
            )
            stats = stats.reindex(range(len(grid)))
            empty = grid[f"empty_{sensor}"].to_numpy()
            assert np.isnan(grid[col].to_numpy()[empty]).all()
            assert (grid[f"n_{sensor}"].to_numpy()[~empty] == stats["size"].to_numpy()[~empty]).all()
            m = grid[col].to_numpy()[~empty]
            assert (m >= stats["min"].to_numpy()[~empty] - 1e-9).all()
            assert (m <= stats["max"].to_numpy()[~empty] + 1e-9).all()
