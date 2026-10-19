"""Load and validate ``config/experiment.toml``.

The configuration holds the experimenter-confirmed facts. Values that have not
been measured are absent from the file and listed under ``[missing] keys``.
Measured values, once they exist, belong in a ``[measurements]`` table under the
same key names. Code asks for them with :meth:`Config.measurement`, which never
substitutes a default.
"""

from __future__ import annotations

import itertools
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "experiment.toml"

# Required keys and their accepted Python types.
_SCHEMA: dict[str, dict[str, tuple[type, ...]]] = {
    "vessel": {"volume_L": (float,), "open_to_air": (bool,), "stirring": (str,)},
    "pulses": {
        "salt": (str,),
        "mass_mg": (float,),
        "molar_mass_g_mol": (float,),
        "times_d": (list,),
    },
    "sensors": {
        "ec_slope_mS_cm_per_A": (float,),
        "ec_intercept_mS_cm": (float,),
        "nh4_unit": (str,),
        "nh4_to_N_factor": (float,),
        "orp_working": (str,),
        "orp_reference": (str,),
    },
    "scenarios": {
        "temperature_C": (list,),
        "baseline_chloride_mM": (list,),
        "pCO2_air_atm": (float,),
    },
    "missing": {"keys": (list,)},
}


class ConfigError(ValueError):
    """The configuration file is malformed or internally inconsistent."""


class MissingInputError(KeyError):
    """A required input is absent from the configuration."""

    def __init__(self, key: str):
        super().__init__(key)
        self.key = key

    def __str__(self) -> str:
        return not_computed(self.key)


def not_computed(key: str) -> str:
    """The exact marker written to outputs when ``key`` is absent."""
    return f"not computed: {key} missing"


@dataclass(frozen=True)
class Config:
    raw: dict[str, Any]
    path: Path
    missing_keys: tuple[str, ...] = field(default=())

    def get(self, dotted: str) -> Any:
        """Return ``section.key``; raise :class:`MissingInputError` if absent."""
        node: Any = self.raw
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                raise MissingInputError(dotted)
            node = node[part]
        return node

    def has_measurement(self, key: str) -> bool:
        return key in self.raw.get("measurements", {})

    def measurement(self, key: str) -> Any:
        """Return a measured value from ``[measurements]``.

        Raises :class:`MissingInputError` (whose ``str`` is the
        ``not computed: <key> missing`` marker) when it has not been measured.
        """
        if not self.has_measurement(key):
            raise MissingInputError(key)
        return self.raw["measurements"][key]

    @property
    def pulse_times_d(self) -> list[float]:
        return [float(t) for t in self.get("pulses.times_d")]


def _check_types(raw: dict[str, Any], path: Path) -> None:
    for section, keys in _SCHEMA.items():
        if section not in raw:
            raise ConfigError(f"{path}: missing section [{section}]")
        for key, types in keys.items():
            if key not in raw[section]:
                raise ConfigError(f"{path}: missing key {section}.{key}")
            value = raw[section][key]
            # bool is a subclass of int; never accept it for a numeric field.
            if not isinstance(value, types) or (bool not in types and isinstance(value, bool)):
                raise ConfigError(
                    f"{path}: {section}.{key} must be {types[0].__name__}, got {value!r}"
                )


# Optional measurements (W8, W2): accepted type and the bound they must satisfy once measured.
OPTIONAL_MEASUREMENTS: dict[str, tuple[tuple[type, ...], str]] = {
    "vessel_mass_start_g": ((int, float), "positive"),
    "vessel_mass_end_g": ((int, float), "positive"),
    "covered": ((bool,), "bool"),
    "top_up_mL": ((int, float), "non-negative"),
    "initial_alkalinity_meq_L": ((int, float), "non-negative"),
}


def _check_optional_measurements(measured: dict[str, Any], path: Path) -> None:
    for key, (types, rule) in OPTIONAL_MEASUREMENTS.items():
        if key not in measured:
            continue
        v = measured[key]
        if rule == "bool":
            if not isinstance(v, bool):
                raise ConfigError(f"{path}: measurements.{key} must be true or false, got {v!r}")
            continue
        if isinstance(v, bool) or not isinstance(v, types):
            raise ConfigError(f"{path}: measurements.{key} must be a number, got {v!r}")
        if (rule == "positive" and not v > 0) or (rule == "non-negative" and not v >= 0):
            raise ConfigError(f"{path}: measurements.{key} must be {rule}, got {v!r}")
    if ("vessel_mass_start_g" in measured) != ("vessel_mass_end_g" in measured):
        raise ConfigError(f"{path}: vessel_mass_start_g and vessel_mass_end_g must be measured together")


def _check_values(raw: dict[str, Any], path: Path) -> None:
    if raw["vessel"]["volume_L"] <= 0:
        raise ConfigError(f"{path}: vessel.volume_L must be positive")
    if raw["pulses"]["mass_mg"] <= 0 or raw["pulses"]["molar_mass_g_mol"] <= 0:
        raise ConfigError(f"{path}: pulse mass and molar mass must be positive")
    times = raw["pulses"]["times_d"]
    if not times or any(not isinstance(t, (int, float)) for t in times):
        raise ConfigError(f"{path}: pulses.times_d must be a non-empty list of numbers")
    if any(b <= a for a, b in itertools.pairwise(times)):
        raise ConfigError(f"{path}: pulses.times_d must be strictly increasing")
    for key in ("temperature_C", "baseline_chloride_mM"):
        values = raw["scenarios"][key]
        if not values or any(not isinstance(v, (int, float)) for v in values):
            raise ConfigError(f"{path}: scenarios.{key} must be a non-empty list of numbers")
    if any(c <= 0 for c in raw["scenarios"]["baseline_chloride_mM"]):
        raise ConfigError(f"{path}: scenarios.baseline_chloride_mM must be positive")
    for key in ("evaporation_mL_per_d", "kLa_CO2_per_d"):
        if key in raw["scenarios"]:
            values = raw["scenarios"][key]
            if (not isinstance(values, list) or not values
                    or any(isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0 for v in values)):
                raise ConfigError(f"{path}: scenarios.{key} must be a non-empty list of non-negative numbers")
    if any(v <= 0 for v in raw["scenarios"].get("kLa_CO2_per_d", [1.0])):
        raise ConfigError(f"{path}: scenarios.kLa_CO2_per_d must be positive")
    for key in ("kLa_O2_per_d", "sae_field_kgO2_per_kWh"):  # W4 SCENARIO ranges: strictly positive
        if key in raw["scenarios"]:
            values = raw["scenarios"][key]
            if (not isinstance(values, list) or not values
                    or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not v > 0 for v in values)):
                raise ConfigError(f"{path}: scenarios.{key} must be a non-empty list of positive numbers")
    if "sae_field_kgO2_per_kWh" in raw["scenarios"] and not str(raw["scenarios"].get("sae_source", "")).strip():
        raise ConfigError(f"{path}: scenarios.sae_field_kgO2_per_kWh needs a cited scenarios.sae_source")
    _check_optional_measurements(raw.get("measurements", {}), path)
    missing = raw["missing"]["keys"]
    measured = set(raw.get("measurements", {}))
    both = sorted(measured.intersection(missing))
    if both:
        raise ConfigError(f"{path}: keys listed as missing but also measured: {both}")


def load_config(path: Path | str = CONFIG_PATH) -> Config:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"configuration file not found: {path}")
    with path.open("rb") as fh:
        raw = tomllib.load(fh)
    _check_types(raw, path)
    _check_values(raw, path)
    return Config(raw=raw, path=path, missing_keys=tuple(raw["missing"]["keys"]))
