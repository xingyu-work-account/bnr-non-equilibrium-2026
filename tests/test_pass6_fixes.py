"""Pass 6 (2026-10-26): G1 registration hash, G2 null coverage of the block bootstrap, G3/G4/G6/G7 SCENARIO and
ASSUMPTION labels, G8 provenance columns of the M6 and W2-consistency tables."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from bnr_neq import states as ST

ROOT = Path(__file__).resolve().parents[1]


# --- G1: registration -------------------------------------------------------------------------------------------
def _registered_hash() -> str:
    text = (ROOT / ST.REGISTRATION_DOC).read_text(encoding="utf-8")
    m = re.search(r"^registered_sha256: ([0-9a-f]{64})$", text, flags=re.MULTILINE)
    if m is None:
        pytest.fail(f"{ST.REGISTRATION_DOC} has no registered_sha256 line")
    return m.group(1)


def test_code_matches_the_registered_rules():
    assert ST.registration_sha256() == _registered_hash()


@pytest.mark.parametrize("name, value", [("SLOPE_HALF_BINS", 49), ("ORP_LEVEL_MV", 14), ("N_BOOT", 998),
                                         ("MIN_BLOCKS_PER_GROUP", 1), ("NULL_COVERAGE_MIN_OK", 0.85)])
def test_any_change_to_a_registered_parameter_changes_the_hash(monkeypatch, name, value):
    before = ST.registration_sha256()
    monkeypatch.setattr(ST, name, value)
    assert ST.registration_sha256() != before


def test_a_changed_rule_changes_the_hash(monkeypatch):
    before = ST.registration_sha256()
    monkeypatch.setattr(ST, "RULES", (*ST.RULES[:2], ST.Rule("inactive_oversupplied", "changed", ("ph",))))
    assert ST.registration_sha256() != before


def test_registration_status_says_not_pre_registered():
    assert ST.REGISTRATION_STATUS.startswith("REVISED AFTER FIRST RESULTS (not pre-registered)")
    rules = ST.rules_table(set())
    assert (rules["registration_sha256"] == ST.registration_sha256()).all()
    assert (rules["registration_status"] == ST.REGISTRATION_STATUS).all()


# --- G2: null coverage of the block bootstrap -------------------------------------------------------------------
def test_null_response_has_the_overlap_dependence_and_no_group_effect():
    y = ST.null_response(3000, 0.0, seed=1)
    assert len(y) == 3000
    assert np.corrcoef(y[:-1], y[1:])[0, 1] > 0.95  # neighbouring values share 100 of 101 points
    assert abs(y.mean()) < 4 * y.std() / np.sqrt(3000 / ST.acf_block_length(y, 1))


def test_null_coverage_counts_the_replicates_whose_interval_covers_zero():
    n = 2000
    a = (np.arange(n) // 250) % 2 == 0  # long runs, as the real labels
    rows = ST.bootstrap_null_coverage(a, ~a, phis=(0.0,), n_rep=40, n_boot=199, seed=3)
    (r,) = rows
    assert r["n_replicates"] == 40 and r["n_computed"] == 40
    assert 0.0 <= r["coverage_of_zero"] <= 1.0
    assert r["block_bins_median"] >= 2 * ST.SLOPE_HALF_BINS + 1
    # determinism (fixed seeds)
    assert ST.bootstrap_null_coverage(a, ~a, phis=(0.0,), n_rep=40, n_boot=199, seed=3) == rows


def test_small_groups_give_no_computed_replicate():
    a = np.zeros(1000, dtype=bool)
    a[:50] = True  # fewer than 2 blocks of 101 bins
    (r,) = ST.bootstrap_null_coverage(a, ~a, phis=(0.5,), n_rep=5, n_boot=99)
    assert r["n_computed"] == 0 and np.isnan(r["coverage_of_zero"])


def _comp(status="computed"):
    import pandas as pd

    return pd.DataFrame([{"comparison": "x (pooled)", "status": status}])


def _cov(values):
    import pandas as pd

    return pd.DataFrame({"comparison": "x (pooled)", "coverage_of_zero": values})


@pytest.mark.parametrize("values, start", [([0.97, 0.95], "adequate"), ([0.95, ST.NULL_COVERAGE_MIN_OK], "adequate"),
                                           ([0.95, np.nextafter(ST.NULL_COVERAGE_MIN_OK, 0)], "UNDER-COVERS"),
                                           ([0.6, 0.95], "UNDER-COVERS")])
def test_interval_reading_is_decided_by_the_worst_phi(values, start):
    out = ST.interval_reading(_comp(), _cov(values))
    assert out["null_coverage_min"].iloc[0] == min(values)
    assert out["interval_reading"].iloc[0].startswith(start)


def test_a_comparison_without_an_interval_is_not_read():
    out = ST.interval_reading(_comp("not computed: small"), _cov([0.5]))
    assert out["interval_reading"].iloc[0] == "no interval"
