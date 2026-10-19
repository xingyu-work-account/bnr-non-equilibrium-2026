import shutil

import numpy as np
import pytest

from bnr_neq import io


@pytest.mark.parametrize(
    "key, rows_read, rows_valid, t_min, t_max",
    [
        ("cond", 71760, 71759, 0.000693, 49.6946),
        ("nh4", 414785, 414785, 0.0, 47.92),
        ("orp", 459884, 459884, 0.0, 53.23),
        ("ph", 414757, 414757, 0.0, 48.0044),
    ],
)
def test_row_counts_and_spans(ingested, key, rows_read, rows_valid, t_min, t_max):
    load = ingested.loads[key]
    assert load.n_rows_read == rows_read
    assert len(load.data) == rows_valid
    t = load.data["t_d"]
    assert t.min() == pytest.approx(t_min, abs=5e-7)
    assert t.max() == pytest.approx(t_max, abs=5e-5)


def test_dropped_row_count_is_one(ingested):
    assert len(ingested.dropped) == 1
    row = ingested.dropped.iloc[0]
    assert row["file"] == "Cond.csv"
    assert row["line"] == 71761  # last line of the file
    assert row["raw_time"] == "49.69529086" and row["raw_value"] == ""


@pytest.mark.parametrize(
    "key, lo, hi",
    [("cond", -0.366, -0.05946), ("nh4", 9.52, 75.13), ("orp", 77.96, 318.7), ("ph", 6.47, 8.02)],
)
def test_raw_value_ranges(ingested, key, lo, hi):
    v = ingested.loads[key].data[io.SENSORS[key].columns[1]]
    assert v.min() == pytest.approx(lo) and v.max() == pytest.approx(hi)


def test_ec_conversion_reference_value(cfg):
    ec = io.convert_ec(np.array([-0.1]), cfg.get("sensors.ec_slope_mS_cm_per_A"),
                       cfg.get("sensors.ec_intercept_mS_cm"))
    assert ec[0] == pytest.approx(0.810725, abs=1e-9)


def test_ec_conversion_range_on_real_data(ingested):
    ec = ingested.natives["cond"]["EC_mS_cm"]
    assert ec.min() == pytest.approx(0.6932, abs=5e-5)
    assert ec.max() == pytest.approx(1.5817, abs=5e-5)


def test_ec_conversion_in_milliamperes_raises(cfg):
    with pytest.raises(ValueError, match="outside 0.0 to 5.0 mS/cm"):
        io.convert_ec(np.array([-144.6]), cfg.get("sensors.ec_slope_mS_cm_per_A"),
                      cfg.get("sensors.ec_intercept_mS_cm"))


def test_nh4_n_column_kept_alongside(ingested, cfg):
    df = ingested.natives["nh4"]
    np.testing.assert_allclose(df["nh4_mgN_L"], df["nh4_mgNH4_L"] * 0.77649)
    assert cfg.get("sensors.nh4_to_N_factor") == 0.77649


def _copy_data(tmp_path, *names):
    for n in names:
        shutil.copy(io.DATA_DIR / n, tmp_path / n)
    return tmp_path


def test_wrong_header_raises_naming_file(tmp_path):
    p = tmp_path / "ORP.csv"
    p.write_bytes(b"Time,ORP (mV)\r\n0,123.6\r\n0,123.9\r\n")
    with pytest.raises(io.HeaderError, match="ORP.csv"):
        io.read_raw(io.SENSORS["orp"], tmp_path)


def test_missing_file_raises_naming_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="pH.csv"):
        io.read_raw(io.SENSORS["ph"], tmp_path)


def test_bom_and_crlf_are_handled(tmp_path):
    p = tmp_path / "NH4.csv"
    p.write_bytes(b"\xef\xbb\xbfTime(d),Concentration(mg/L)\r\n0,44.9\r\n0.01,45.0\r\n")
    load = io.read_raw(io.SENSORS["nh4"], tmp_path)
    assert list(load.data.columns) == ["t_d", "nh4_mgNH4_L"]
    assert load.data["nh4_mgNH4_L"].tolist() == [44.9, 45.0]


def test_non_numeric_rows_are_counted_with_line_numbers(tmp_path):
    p = tmp_path / "pH.csv"
    p.write_bytes(b"Time(d),pH\r\n0,7.4\r\n0.1,bad\r\n0.2,7.5\r\n0.3,\r\n")
    load = io.read_raw(io.SENSORS["ph"], tmp_path)
    assert len(load.data) == 2
    assert load.dropped["line"].tolist() == [3, 5]
    assert load.dropped["file"].unique().tolist() == ["pH.csv"]


def test_decreasing_time_raises(tmp_path):
    p = tmp_path / "pH.csv"
    p.write_bytes(b"Time(d),pH\r\n0,7.4\r\n0.2,7.5\r\n0.1,7.5\r\n")
    with pytest.raises(io.MonotonicityError, match="pH.csv"):
        io.read_raw(io.SENSORS["ph"], tmp_path)
