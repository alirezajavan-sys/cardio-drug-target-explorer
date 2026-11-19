"""Refresh the ChEMBL snapshot from the public ChEMBL REST API.

    python scripts/fetch_chembl.py            # writes data/chembl_sglt_ic50_snapshot.csv

Pulls every human SGLT1 / SGLT2 IC50 measurement ("=" relation, nM) for the five SGLT2 inhibitors tested in
cardiovascular and kidney outcome trials. The schema matches the committed snapshot, so `run_analysis.py` works
unchanged on refreshed data. Requires `requests` and internet access.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import requests

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
OUT = Path(__file__).resolve().parents[1] / "data" / "chembl_sglt_ic50_snapshot.csv"

DRUGS = {
    "dapagliflozin": "CHEMBL429910",
    "empagliflozin": "CHEMBL2107830",
    "canagliflozin": "CHEMBL2048484",
    "ertugliflozin": "CHEMBL1770248",
    "sotagliflozin": "CHEMBL3039507",
}
TARGETS = {"SGLT2": "CHEMBL3884", "SGLT1": "CHEMBL4979"}
HEADER = ["activity_id", "record_id", "drug", "molecule_chembl_id", "target", "target_chembl_id", "ic50_nM",
          "pchembl_value", "assay_chembl_id", "document_chembl_id", "document_year", "potential_duplicate"]


def fetch(molecule: str, target: str) -> list[dict]:
    params = {"molecule_chembl_id": molecule, "target_chembl_id": target, "standard_type": "IC50",
              "standard_relation": "=", "standard_units": "nM", "limit": 1000, "format": "json"}
    r = requests.get(API, params=params, timeout=60)
    r.raise_for_status()
    return r.json()["activities"]


def main() -> int:
    rows = []
    for drug, mol in DRUGS.items():
        for label, tgt in TARGETS.items():
            acts = fetch(mol, tgt)
            for a in acts:
                if a.get("standard_value") is None:
                    continue
                rows.append([a["activity_id"], a["record_id"], drug, mol, label, tgt, float(a["standard_value"]),
                             a.get("pchembl_value"), a["assay_chembl_id"], a["document_chembl_id"],
                             a.get("document_year"), a.get("potential_duplicate", 0)])
            print(f"{drug:14s} {label}: {len(acts)} measurements", file=sys.stderr)
            time.sleep(0.5)  # be polite to the public API
    rows.sort(key=lambda r: (r[2], r[4], r[0]))
    with OUT.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {OUT}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
