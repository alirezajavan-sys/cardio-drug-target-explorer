"""SGLT2-vs-SGLT1 potency and selectivity of SGLT2 inhibitors from ChEMBL bioactivity data.

Selectivity is expressed as a fold ratio:  IC50(SGLT1) / IC50(SGLT2) = 10 ** (pIC50_SGLT2 - pIC50_SGLT1).

Two estimators are provided because IC50 values vary between laboratories and assay formats:

* **pooled**: median pIC50 per target over all independent measurements, then the ratio;
* **paired (within-document)**: ratio computed inside each publication that measured *both* transporters under the
  same conditions, then the median across publications. This cancels between-lab offsets and is the preferred estimate.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA_PATH = Path(__file__).resolve().parents[2] / "data" / "chembl_sglt_ic50_snapshot.csv"

TARGETS = {"SGLT2": "CHEMBL3884", "SGLT1": "CHEMBL4979"}

# Outcome trials in which each drug was tested (links this analysis to the sglt2-inhibitor-meta-analysis repo).
TRIALS = {
    "empagliflozin": "EMPA-REG OUTCOME, EMPEROR-Reduced, EMPEROR-Preserved, EMPA-KIDNEY",
    "dapagliflozin": "DECLARE-TIMI 58, DAPA-HF, DAPA-CKD, DELIVER",
    "canagliflozin": "CANVAS, CREDENCE",
    "ertugliflozin": "VERTIS CV",
    "sotagliflozin": "SCORED, SOLOIST-WHF",
}


def load_snapshot(path: Path | str = DATA_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"activity_id", "record_id", "drug", "target", "ic50_nM", "pchembl_value", "document_chembl_id",
                "potential_duplicate"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if not set(df["target"]).issubset(TARGETS):
        raise ValueError("Unexpected target labels")
    df["pIC50"] = 9 - np.log10(df["ic50_nM"])
    return df


def independent(df: pd.DataFrame) -> pd.DataFrame:
    """Drop ChEMBL-flagged citations of earlier data and exact re-depositions of the same record."""
    out = df[df["potential_duplicate"] == 0]
    return out.drop_duplicates(subset=["record_id", "target", "ic50_nM"]).copy()


def pooled_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Median (IQR) pIC50 per drug and target, and the pooled selectivity fold."""
    g = df.groupby(["drug", "target"])["pIC50"]
    s = g.agg(n="size", median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75)).reset_index()
    wide = s.pivot(index="drug", columns="target")
    out = pd.DataFrame({
        "n_SGLT2": wide[("n", "SGLT2")].astype(int),
        "pIC50_SGLT2": wide[("median", "SGLT2")],
        "n_SGLT1": wide[("n", "SGLT1")].astype(int),
        "pIC50_SGLT1": wide[("median", "SGLT1")],
    })
    out["IC50_SGLT2_nM"] = 10 ** (9 - out["pIC50_SGLT2"])
    out["IC50_SGLT1_nM"] = 10 ** (9 - out["pIC50_SGLT1"])
    out["selectivity_fold_pooled"] = 10 ** (out["pIC50_SGLT2"] - out["pIC50_SGLT1"])
    return out


def paired_selectivity(df: pd.DataFrame) -> pd.DataFrame:
    """Selectivity fold within each publication that reports both SGLT1 and SGLT2 for the same drug."""
    m = df.groupby(["drug", "document_chembl_id", "target"])["pIC50"].mean().unstack("target").dropna()
    m["fold"] = 10 ** (m["SGLT2"] - m["SGLT1"])
    return m.reset_index()


def selectivity_table(df: pd.DataFrame) -> pd.DataFrame:
    ind = independent(df)
    pooled = pooled_summary(ind)
    paired = paired_selectivity(ind)
    pstats = paired.groupby("drug")["fold"].agg(n_documents="size", selectivity_fold_paired="median",
                                                 fold_min="min", fold_max="max")
    out = pooled.join(pstats)
    out["outcome_trials"] = out.index.map(TRIALS)
    return out.sort_values("selectivity_fold_paired", ascending=False)
