import numpy as np
import pytest

from bnr_neq import io, pipeline


def test_build_native_reconstructs_only_stamped_sensors(ingested):
    assert {"t_stamp_d", "stamp_j", "stamp_n"} <= set(ingested.natives["nh4"].columns)
    assert {"t_stamp_d", "stamp_j", "stamp_n"} <= set(ingested.natives["orp"].columns)
    assert "t_stamp_d" not in ingested.natives["ph"].columns and "t_stamp_d" not in ingested.natives["cond"].columns


def test_data_properties_flag_header_discrepancies(ingested):
    p = ingested.properties.set_index("sensor")
    assert p.loc["ph", "header_matches_section4"] and not p.loc["cond", "header_matches_section4"]
    assert p.loc["cond", "rows_dropped"] == 1 and p["rows_dropped"].sum() == 1
    assert p.loc["ph", "tie_fraction"] == pytest.approx(0.884, abs=0.001)


def test_nh4_rate_regimes(ingested):
    r = pipeline.nh4_rate_regimes(ingested.natives["nh4"])
    np.testing.assert_allclose(r["implied_interval_s"], [20.1, 3.0, 2.0], atol=0.15)


def test_grid_variables_cover_every_sensor():
    assert set(pipeline.GRID_VARIABLES) == set(io.SENSORS)


def test_grid_coverage(ingested):
    c = pipeline.grid_coverage(ingested).set_index("sensor")
    assert (c["empty_bins"] == 0).all()
    assert c.loc["orp", "samples_outside_grid"] > 0  # ORP runs to day 53.23
