import pytest

from bnr_neq.config import (
    CONFIG_PATH,
    ConfigError,
    MissingInputError,
    load_config,
    not_computed,
)


def test_confirmed_facts_are_encoded(cfg):
    assert cfg.get("vessel.volume_L") == 0.100
    assert cfg.get("vessel.open_to_air") is True
    assert cfg.get("vessel.stirring") == "magnetic"
    assert cfg.get("pulses.salt") == "NH4Cl"
    assert cfg.get("pulses.mass_mg") == 25.0
    assert cfg.pulse_times_d == [8.53, 10.57, 24.59, 32.39]
    assert cfg.get("sensors.nh4_unit") == "mg NH4+/L"
    assert cfg.get("sensors.nh4_to_N_factor") == pytest.approx(14.007 / 18.039, abs=1e-5)
    assert cfg.get("sensors.orp_working") == "Pt"
    assert cfg.get("sensors.orp_reference") == "bare Ag"


@pytest.mark.parametrize(
    "key",
    ["temperature_C_measured", "alkalinity_mmol_L", "DO_mg_L", "kLa_per_d", "stirrer_power_W",
     "chloride_mmol_L", "orp_reference_offset_V", "lab_NH4_N_mg_L", "NO3_N_mg_L", "NO2_N_mg_L",
     "electrode_step_response"],
)
def test_missing_measurement_is_not_computed(cfg, key):
    assert key in cfg.missing_keys
    assert not cfg.has_measurement(key)
    with pytest.raises(MissingInputError) as exc:
        cfg.measurement(key)
    assert str(exc.value) == f"not computed: {key} missing" == not_computed(key)


def test_absent_config_key_raises(cfg):
    with pytest.raises(MissingInputError):
        cfg.get("vessel.temperature_C")


def _write(tmp_path, text):
    p = tmp_path / "experiment.toml"
    p.write_text(text)
    return p


def test_missing_section_raises(tmp_path):
    text = CONFIG_PATH.read_text().replace("[vessel]", "[vessel_renamed]")
    with pytest.raises(ConfigError, match=r"\[vessel\]"):
        load_config(_write(tmp_path, text))


def test_wrong_type_raises(tmp_path):
    text = CONFIG_PATH.read_text().replace("volume_L = 0.100", "volume_L = true")
    with pytest.raises(ConfigError, match="volume_L"):
        load_config(_write(tmp_path, text))


def test_non_increasing_pulse_times_raise(tmp_path):
    text = CONFIG_PATH.read_text().replace("[8.53, 10.57, 24.59, 32.39]", "[8.53, 24.59, 10.57, 32.39]")
    with pytest.raises(ConfigError, match="strictly increasing"):
        load_config(_write(tmp_path, text))


def test_key_both_missing_and_measured_raises(tmp_path):
    text = CONFIG_PATH.read_text() + "\n[measurements]\nDO_mg_L = 2.0\n"
    with pytest.raises(ConfigError, match="DO_mg_L"):
        load_config(_write(tmp_path, text))


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "nope.toml")
