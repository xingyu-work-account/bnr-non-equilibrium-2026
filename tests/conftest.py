import pytest

from bnr_neq import pipeline
from bnr_neq.config import CONFIG_PATH, load_config


def config_with(path, **measured):
    """Configuration with the given keys moved from [missing] to [measurements], written to ``path``."""
    text = CONFIG_PATH.read_text()
    for k in measured:
        assert f'"{k}"' in text
        text = text.replace(f'"{k}", ', "").replace(f', "{k}"', "")
    text += "\n[measurements]\n" + "".join(f"{k} = {v}\n" for k, v in measured.items())
    path.write_text(text)
    return load_config(path)


@pytest.fixture
def measured_cfg(tmp_path):
    """Factory: ``measured_cfg(DO_mg_L=0.1, ...)`` returns a Config with those measurements present."""
    return lambda **measured: config_with(tmp_path / "experiment_measured.toml", **measured)


@pytest.fixture(scope="session")
def cfg():
    return load_config()


@pytest.fixture(scope="session")
def ingested(cfg):
    """The production ingest of the real raw files (read-only)."""
    return pipeline.ingest(cfg)


@pytest.fixture(scope="session")
def validated(cfg, ingested):
    """The production M2 computation on the ingested natives (no files written)."""
    from bnr_neq import validation

    return validation.compute(cfg, ingested.natives)


@pytest.fixture(scope="session")
def segmented(cfg, ingested, validated):
    """The production M3 computation on the M2 results (no files written)."""
    from bnr_neq import segments

    return segments.compute(cfg, ingested.grid, validated.jumps, validated.pulse_table, validated.held)


@pytest.fixture(scope="session")
def estimator_validation():
    from bnr_neq import dynamics

    return dynamics.validate_estimators()


@pytest.fixture(scope="session")
def dynamics_inputs(ingested, validated, segmented):
    n = ingested.natives
    return {"ph": n["ph"], "orp": n["orp"], "nh4": n["nh4"], "cond": n["cond"], "grid": ingested.grid,
            "held": validated.held, "segments": segmented.table, "pulses": validated.pulse_table}


@pytest.fixture(scope="session")
def dynamics_result(cfg, dynamics_inputs, estimator_validation):
    from bnr_neq import dynamics

    return dynamics.compute(cfg, dynamics_inputs, estimator_validation)


@pytest.fixture(scope="session")
def run_all_dir(tmp_path_factory):
    """Output tree of the production ``run-all`` CLI on the real data (first run)."""
    from bnr_neq.__main__ import main

    out = tmp_path_factory.mktemp("run_all_1")
    assert main(["run-all", "--out-dir", str(out)]) == 0
    return out


def synthetic_truth(grid, scales=(5.0, 4.0, 6.0)):
    """A known NH4 truth for stage tests: C = scale x reading, written in the W1 file formats.

    Vector k has a = −ln(scale_k), b0 = 1, b1 = 0 and no offsets, so the M6 inversion of the reading gives exactly
    scale_k x reading; vector 0 is the 'best' series."""
    import numpy as np
    import pandas as pd

    from bnr_neq import model as M

    g = grid[(grid["t_start_d"] >= M.T0_D - 1e-9) & (grid["t_end_d"] <= M.END_D + 1e-9)]
    reading = g["nh4_mgN_L"].to_numpy()
    c = scales[0] * reading
    truth = pd.DataFrame({"t_d": g["t_mid_d"].to_numpy(), "electrode_obs_mgN_L": reading, "pH": g["pH"].to_numpy(),
                          "C_inv_mgN_L": c, "C_inv_band_low_mgN_L": min(scales) * reading,
                          "C_inv_band_high_mgN_L": max(scales) * reading, "gain_b": 1.0,
                          "nh4_observation_model": M.OBSERVATION_MODEL, "nh4_truth_id": M.truth_id(c)})
    vectors = pd.DataFrame([{"vector": k, "profiled_parameter": "best" if k == 0 else "a",
                             "profiled_value_internal": np.nan, "delta": 0.0, "a": -np.log(s), "b0": 1.0, "b1_per_d": 0.0}
                            for k, s in enumerate(scales)])
    summary = pd.DataFrame([("nh4_truth_id", M.truth_id(c)), ("band_low_open", False), ("band_high_open", True)],
                           columns=["key", "value"])
    return truth, vectors, summary
