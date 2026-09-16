"""
07_gene_panel.py  --  GSE44711: per-sample expression panel for the genes that matter

Rows
  1. H2S enzymes : CTH (detection p-values, floored out on all 3 probes), CBS, MPST
  2. canonical EOPE : HTRA4, LEP, FSTL3, PAPPA2, INHA, ENG
  3. stromal set (down in EOPE) : PDPN, HAND2, OSR1, CNN2, ANGPTL2

Markers
  red = EOPET, blue = Control
  hollow circle = GSM1089242 (technical outlier, excluded from the primary analysis)
  triangle      = GSM1089235 / GSM1089236 (EOPET with LEP at floor; ..36 also least PE-like on PCA)
Panel titles carry log2FC and BH padj from the 8 vs 7 primary run (results/04_stats.csv).
One probe per gene: the probe with the smallest padj.

Run from project root:  python scripts/07_gene_panel.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths
PROJECT = Path(__file__).resolve().parent.parent
PROC = PROJECT / "data" / "processed"
RESULTS = PROJECT / "results"
FIG = PROJECT / "figures"
FIG.mkdir(exist_ok=True)

EXCLUDED = ["GSM1089242"]
FLOOR_LEP = ["GSM1089235", "GSM1089236"]
DET_P = 0.05

# ---------------------------------------------------------------- load + seatbelts
expr = pd.read_csv(PROC / "expr_log2.csv", index_col=0)
meta = pd.read_csv(PROC / "meta.csv", index_col=0)
annot = pd.read_csv(PROC / "annotation.csv", index_col=0)
detp = pd.read_csv(PROC / "detection_pvals.csv", index_col=0)
stats = pd.read_csv(RESULTS / "04_stats.csv", index_col=0)
assert all(meta.index == expr.columns) and all(meta.index == detp.columns)
assert set(stats.index) <= set(expr.index)

is_pe = meta["diagnosis"] == "EOPET"
x_pos = np.where(is_pe, 1.0, 0.0)
rng = np.random.default_rng(0)
jitter = rng.uniform(-0.12, 0.12, size=len(meta))

ROWS = [
    ("H2S enzymes", ["CTH", "CBS", "MPST"]),
    ("canonical early-onset PE genes", ["HTRA4", "LEP", "FSTL3", "PAPPA2", "INHA", "ENG"]),
    ("stromal / mesenchymal set", ["PDPN", "HAND2", "OSR1", "CNN2", "ANGPTL2"]),
]


def best_probe(gene):
    """Probe with the smallest padj among tested probes for this gene."""
    hits = stats[stats["symbol"] == gene]
    assert len(hits) > 0, f"{gene}: no tested probe"
    return hits["padj"].idxmin()


def draw_points(ax, y):
    for i, g in enumerate(meta.index):
        color = "#c0392b" if is_pe[g] else "#2c7fb8"
        marker = "^" if g in FLOOR_LEP else "o"
        hollow = g in EXCLUDED
        ax.scatter(x_pos[i] + jitter[i], y[g], s=42, marker=marker,
                   facecolors="none" if hollow else color, edgecolors=color, linewidths=1.2, zorder=3)


def draw_means(ax, y):
    kept = ~meta.index.isin(EXCLUDED)
    for grp, xc in [("Control", 0.0), ("EOPET", 1.0)]:
        m = (meta["diagnosis"] == grp).values & kept
        ax.hlines(y[m].mean(), xc - 0.25, xc + 0.25, color="k", lw=1.5, zorder=4)


# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(3, 6, hspace=0.55, wspace=0.45)

for r, (row_title, genes) in enumerate(ROWS):
    span = 2 if r == 0 else 1
    for c, gene in enumerate(genes):
        ax = fig.add_subplot(gs[r, c * span:(c + 1) * span])

        if gene == "CTH":
            # floored out: show detection p-values per sample for the three probes
            probes = [p for p in annot.index[annot["symbol"] == "CTH"] if p in detp.index]
            for k, pid in enumerate(probes):
                y = detp.loc[pid]
                for i, g in enumerate(meta.index):
                    color = "#c0392b" if is_pe[g] else "#2c7fb8"
                    ax.scatter(x_pos[i] + jitter[i] + (k - 1) * 0.0, y[g], s=18, color=color,
                               alpha=0.7, marker="^" if g in FLOOR_LEP else "o", linewidths=0)
            ax.axhline(DET_P, color="k", lw=0.8, ls="--")
            n_det = (detp.loc[probes] < DET_P).sum(axis=1)
            ax.set_title("CTH  --  below detection on 3 probes\n"
                         "detected in " + ", ".join(f"{int(n)}/16" for n in n_det) + " samples",
                         fontsize=9)
            ax.set_ylabel("Illumina detection p-value", fontsize=8)
            ax.set_ylim(0, 1.02)
        else:
            pid = best_probe(gene)
            y = expr.loc[pid]
            draw_points(ax, y)
            draw_means(ax, y)
            s = stats.loc[pid]
            ax.set_title(f"{gene}  ({pid})\nlog2FC {s['log2FC']:+.2f}   padj {s['padj']:.3f}",
                         fontsize=9)
            ax.set_ylabel("log2 QN expression", fontsize=8)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Control", "EOPET"], fontsize=8)
        ax.set_xlim(-0.5, 1.5)
        ax.tick_params(axis="y", labelsize=8)
        ax.axhline(6.5, color="grey", lw=0.5, ls=":", zorder=1) if gene != "CTH" else None
        if c == 0:
            ax.annotate(row_title, (0, 1.32), xycoords="axes fraction", fontsize=11,
                        fontweight="bold", ha="left")

fig.suptitle("GSE44711  Illumina HT-12 v4, chorionic villi: 8 early-onset PE vs 8 GA-matched preterm controls\n"
             "red = EOPET, blue = Control; hollow = GSM1089242 (excluded); triangle = EOPET with LEP at floor; "
             "dotted = detection floor (~6.5); bars = group mean (8 vs 7)",
             fontsize=10, y=0.995)
fig.savefig(FIG / "07_gene_panel.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- table of what was plotted
rows = []
for _, genes in ROWS:
    for gene in genes:
        if gene == "CTH":
            continue
        pid = best_probe(gene)
        s = stats.loc[pid]
        rows.append({"gene": gene, "probe": pid, "log2FC": round(s["log2FC"], 3),
                     "pval": s["pval"], "padj": round(s["padj"], 4)})
table = pd.DataFrame(rows).set_index("gene")
print(table.to_string())
table.to_csv(RESULTS / "07_gene_panel.csv")
print("\nsaved:", FIG / "07_gene_panel.png", "|", RESULTS / "07_gene_panel.csv")