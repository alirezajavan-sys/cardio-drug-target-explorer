"""Reproduce the tables and figures in the README:  python run_analysis.py"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from cardiotarget.selectivity import independent, load_snapshot, paired_selectivity, selectivity_table  # noqa: E402

FIG, RES = ROOT / "figures", ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
COL = {"SGLT2": "#1f6feb", "SGLT1": "#cf222e"}


def main() -> None:
    raw = load_snapshot()
    ind = independent(raw)
    tab = selectivity_table(raw)
    order = list(tab.index)  # most to least selective

    out = tab.copy()
    out.to_csv(RES / "selectivity.csv", float_format="%.3g")
    pretty = pd.DataFrame({
        "SGLT2 IC50, nM (n)": [f"{a:.2g} ({n})" for a, n in zip(tab["IC50_SGLT2_nM"], tab["n_SGLT2"])],
        "SGLT1 IC50, nM (n)": [f"{a:,.0f} ({n})" for a, n in zip(tab["IC50_SGLT1_nM"], tab["n_SGLT1"])],
        "Selectivity, paired (range)": [f"{f:,.0f}× ({lo:,.0f}–{hi:,.0f})" for f, lo, hi in
                                        zip(tab["selectivity_fold_paired"], tab["fold_min"], tab["fold_max"])],
        "Papers with both": tab["n_documents"].astype(int),
        "Selectivity, pooled": [f"{f:,.0f}×" for f in tab["selectivity_fold_pooled"]],
        "Outcome trials": tab["outcome_trials"],
    }, index=tab.index)
    md = ["# Results", "",
          f"Snapshot: {len(raw)} ChEMBL IC50 records → {len(ind)} independent measurements after removing "
          f"citations/re-depositions.", "",
          "IC50 values are medians of independent measurements; selectivity = IC50(SGLT1) / IC50(SGLT2).", "",
          pretty.to_markdown(), ""]
    (RES / "results.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))

    # ---- Figure 1: potency at each transporter -------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    rng = np.random.default_rng(1)
    for i, drug in enumerate(order):
        for tgt, dx in (("SGLT2", -0.15), ("SGLT1", 0.15)):
            v = ind.loc[(ind.drug == drug) & (ind.target == tgt), "pIC50"]
            ax.scatter(np.full(len(v), i + dx) + rng.uniform(-0.05, 0.05, len(v)), v, s=14, alpha=0.6, color=COL[tgt],
                       label=tgt if i == 0 else None)
            ax.hlines(v.median(), i + dx - 0.12, i + dx + 0.12, color=COL[tgt], lw=2.5)
    ax.set_xticks(range(len(order)), [d.capitalize() for d in order])
    ax.set_ylabel("pIC50  (−log10 M)")
    sec = ax.secondary_yaxis("right", functions=(lambda p: p, lambda p: p))
    ticks = [5, 6, 7, 8, 9]
    sec.set_yticks(ticks, ["10 µM", "1 µM", "100 nM", "10 nM", "1 nM"])
    ax.set_title("Potency at SGLT2 vs SGLT1 (each dot = one independent measurement; bar = median)")
    ax.title.set_fontsize(9.5)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "potency_by_transporter.png")
    plt.close(fig)

    # ---- Figure 2: selectivity ------------------------------------------------------------
    paired = paired_selectivity(ind)
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    y = np.arange(len(order))[::-1]
    for yi, drug in zip(y, order):
        f = paired.loc[paired.drug == drug, "fold"]
        ax.scatter(f, np.full(len(f), yi) + rng.uniform(-0.12, 0.12, len(f)), s=16, color="#8c959f", zorder=2)
        med = tab.loc[drug, "selectivity_fold_paired"]
        ax.barh(yi, med, color="#cf222e" if drug == "sotagliflozin" else "#1f6feb", alpha=0.85, height=0.55, zorder=1)
        ax.text(med * 1.15, yi + 0.3, f"{med:,.0f}×", va="center", fontsize=9)
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}×"))
    ax.set_yticks(y, [d.capitalize() for d in order])
    ax.set_xlim(5, 10000)
    ax.set_xlabel("SGLT2 selectivity over SGLT1 (IC50 ratio, log scale)")
    ax.set_title("Within-paper selectivity (bars = median; dots = individual papers)", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "selectivity.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
