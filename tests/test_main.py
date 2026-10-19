import pytest

from bnr_neq.__main__ import IMPLEMENTED, NOT_IMPLEMENTED, main


def test_every_command_of_section_7_is_implemented():
    assert set(IMPLEMENTED) == {"ingest", "validate", "segments", "thermo", "dynamics", "model", "states", "pulses",
                                "report"}
    assert IMPLEMENTED.index("model") < IMPLEMENTED.index("states") < IMPLEMENTED.index("pulses")  # W3 feeds W4
    assert NOT_IMPLEMENTED == {}


def test_unknown_command_exits_with_error():
    with pytest.raises(SystemExit) as exc:
        main(["nonsense"])
    assert exc.value.code == 2


def test_missing_config_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="configuration file not found"):
        main(["ingest", "--config", str(tmp_path / "none.toml"), "--out-dir", str(tmp_path / "out")])


def test_stage_without_inputs_names_the_missing_stage(tmp_path):
    with pytest.raises(FileNotFoundError, match="ingest"):
        main(["validate", "--out-dir", str(tmp_path)])


def test_ingest_command_writes_outputs(tmp_path, capsys):
    assert main(["ingest", "--out-dir", str(tmp_path)]) == 0
    assert "1 row(s) dropped" in capsys.readouterr().out
    assert (tmp_path / "processed" / "grid_0p01d.parquet").is_file()
    assert (tmp_path / "reports" / "ingest.md").is_file()


def test_run_all_cli_on_tmp_out_dir(run_all_dir):
    from bnr_neq.report import FIGURES, TABLE_DESCRIPTIONS

    for name in TABLE_DESCRIPTIONS:
        assert (run_all_dir / "tables" / name).is_file(), name
    for fig, _ in FIGURES:
        assert (run_all_dir / "figures" / fig).is_file() and (run_all_dir / "figures" / fig.replace(".png", ".svg")).is_file()
    for rep in ("ingest", "sensor_validation", "segments", "thermo", "alkalinity", "dynamics", "model", "states",
                "pulses", "analysis_report"):
        assert (run_all_dir / "reports" / f"{rep}.md").is_file(), rep
    assert not (run_all_dir / "reports" / "progress.md").exists()  # the implementation log lives in docs/


def test_run_all_is_deterministic(run_all_dir, tmp_path_factory):
    """A second run-all into another directory writes byte-identical CSV tables (and SVG figures)."""
    out = tmp_path_factory.mktemp("run_all_2")
    assert main(["run-all", "--out-dir", str(out)]) == 0
    first = sorted(p.relative_to(run_all_dir) for p in (run_all_dir / "tables").glob("*.csv"))
    second = sorted(p.relative_to(out) for p in (out / "tables").glob("*.csv"))
    assert first == second and len(first) > 40
    differ = [str(p) for p in first if (run_all_dir / p).read_bytes() != (out / p).read_bytes()]
    assert differ == []
    svg = [p.relative_to(run_all_dir) for p in (run_all_dir / "figures").glob("*.svg")]
    assert [str(p) for p in svg if (run_all_dir / p).read_bytes() != (out / p).read_bytes()] == []
