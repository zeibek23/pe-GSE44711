"""
06_volcano.py  --  GSE44711: volcano plot, EOPET vs Control (8 vs 7 primary analysis)

Input : results/04_stats.csv          primary (GSM1089242 excluded)
        results/06_robust_set.csv     probes at padj < 0.05 in BOTH 8v7 and 8v8 runs
Output: figures/06_volcano.png
        results/06_labelled_genes.csv

Colouring
  grey   : not significant
  orange : padj < 0.05 in the primary run only
  red    : padj < 0.05 in both runs (robust set)
Vertical lines at |log2FC| = 1, horizontal at padj = 0.05.

Run from project root:  python scripts/06_volcano.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths
PROJECT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT / "results"
FIG = PROJECT / "figures"
FIG.mkdir(exist_ok=True)

PADJ, LFC = 0.05, 1.0

# ---------------------------------------------------------------- load + seatbelts
res = pd.read_csv(RESULTS / "04_stats.csv", index_col=0)
assert {"log2FC", "padj", "symbol"} <= set(res.columns), res.columns
assert res["padj"].between(0, 1).all()
robust_path = RESULTS / "06_robust_set.csv"
robust = set(pd.read_csv(robust_path, index_col=0).index) if robust_path.exists() else set()
print("probes:", len(res), "| robust set:", len(robust))

res["neg_log10_padj"] = -np.log10(res["padj"])
res["sig"] = res["padj"] < PADJ
res["robust"] = res.index.isin(robust)
res["big"] = res["log2FC"].abs() > LFC

# ---------------------------------------------------------------- what to label
# canonical EOPE genes, H2S enzymes, and the stromal set from the top-20 list
CANONICAL = ["HTRA4", "LEP", "FSTL3", "PAPPA2", "INHA", "ENG", "PDPN", "COL17A1"]
H2S = ["CBS", "MPST"]
STROMAL = ["HAND2", "OSR1", "ANGPTL2", "CNN2", "IRX3"]
named = res[res["symbol"].isin(CANONICAL + H2S + STROMAL)]

# plus the ten largest |log2FC| among significant probes not already named
top_effect = (res[res["sig"] & ~res.index.isin(named.index)]
              .assign(abs_lfc=lambda d: d["log2FC"].abs())
              .sort_values("abs_lfc", ascending=False)
              .head(10))
labelled = pd.concat([named, top_effect]).drop_duplicates()
# one label per gene: keep the probe with the smallest padj
labelled = labelled.sort_values("padj").drop_duplicates("symbol")
print("labelled genes:", len(labelled))

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(8.5, 7))
ns = res[~res["sig"]]
primary_only = res[res["sig"] & ~res["robust"]]
both = res[res["sig"] & res["robust"]]

ax.scatter(ns["log2FC"], ns["neg_log10_padj"], s=6, c="#c8c8c8", alpha=0.6, linewidths=0, label=f"n.s. ({len(ns)})")
ax.scatter(primary_only["log2FC"], primary_only["neg_log10_padj"], s=10, c="#f39c12", alpha=0.8,
           linewidths=0, label=f"padj<{PADJ} primary only ({len(primary_only)})")
ax.scatter(both["log2FC"], both["neg_log10_padj"], s=14, c="#c0392b", alpha=0.9,
           linewidths=0, label=f"padj<{PADJ} in both runs ({len(both)})")

ax.axhline(-np.log10(PADJ), color="k", lw=0.7, ls="--")
ax.axvline(-LFC, color="k", lw=0.7, ls="--")
ax.axvline(LFC, color="k", lw=0.7, ls="--")

for pid, row in labelled.iterrows():
    color = "#1a5276" if row["symbol"] in H2S else ("#6c3483" if row["symbol"] in STROMAL else "k")
    ax.annotate(row["symbol"], (row["log2FC"], row["neg_log10_padj"]),
                fontsize=8, color=color, fontweight="bold" if row["symbol"] in H2S else "normal",
                xytext=(5, 3), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.4))

ax.set_xlabel("log2 fold change (EOPET / Control)")
ax.set_ylabel("-log10 BH-adjusted p (Welch)")
ax.set_title("GSE44711  early-onset PE vs preterm control, 8 vs 7\n"
             "19,601 detected probes; blue = H2S enzymes, purple = stromal set")
ax.legend(loc="upper left", fontsize=8, frameon=False)
xlim = max(abs(res["log2FC"].min()), abs(res["log2FC"].max())) * 1.05
ax.set_xlim(-xlim, xlim)
fig.tight_layout()
fig.savefig(FIG / "06_volcano.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- table
cols = ["symbol", "log2FC", "pval", "padj", "robust", "big"]
out = labelled[cols].sort_values("log2FC", ascending=False)
print(out.round(3).to_string())
out.to_csv(RESULTS / "06_labelled_genes.csv")
print("\nsaved:", FIG / "06_volcano.png", "|", RESULTS / "06_labelled_genes.csv")