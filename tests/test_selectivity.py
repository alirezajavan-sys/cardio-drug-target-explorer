import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cardiotarget.selectivity import (  # noqa: E402
    TRIALS, independent, load_snapshot, paired_selectivity, pooled_summary, selectivity_table,
)


@pytest.fixture(scope="module")
def snap():
    return load_snapshot()


def test_snapshot_integrity(snap):
    assert set(snap["drug"]) == set(TRIALS)
    assert set(snap["target"]) == {"SGLT1", "SGLT2"}
    assert snap["activity_id"].is_unique
    # pChEMBL reported by ChEMBL must equal -log10(IC50 in M): catches unit/transcription errors
    assert np.allclose(snap["pchembl_value"], snap["pIC50"], atol=0.01)


def test_independent_removes_duplicates(snap):
    ind = independent(snap)
    assert (ind["potential_duplicate"] == 0).all()
    assert not ind.duplicated(subset=["record_id", "target", "ic50_nM"]).any()
    assert len(ind) < len(snap)


def test_selectivity_maths():
    toy = pd.DataFrame({
        "activity_id": [1, 2, 3, 4], "record_id": [1, 1, 2, 2], "drug": ["a"] * 4,
        "target": ["SGLT2", "SGLT1", "SGLT2", "SGLT1"], "ic50_nM": [1.0, 1000.0, 10.0, 1000.0],
        "pchembl_value": [9, 6, 8, 6], "document_chembl_id": ["d1", "d1", "d2", "d2"], "potential_duplicate": 0,
    })
    toy["pIC50"] = 9 - np.log10(toy["ic50_nM"])
    paired = paired_selectivity(toy).set_index("document_chembl_id")["fold"]
    assert paired["d1"] == pytest.approx(1000)
    assert paired["d2"] == pytest.approx(100)
    pooled = pooled_summary(toy)
    assert pooled.loc["a", "selectivity_fold_pooled"] == pytest.approx(10 ** (8.5 - 6))


def test_sotagliflozin_is_the_dual_inhibitor(snap):
    tab = selectivity_table(snap)
    # Sotagliflozin is a dual SGLT1/2 inhibitor: by far the lowest selectivity in the panel
    assert tab["selectivity_fold_paired"].idxmin() == "sotagliflozin"
    assert tab.loc["sotagliflozin", "selectivity_fold_paired"] < 50
    others = tab.drop("sotagliflozin")["selectivity_fold_paired"]
    assert (others > 100).all()
