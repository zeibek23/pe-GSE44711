from pathlib import Path


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

WATCH = ["GSM1089240", "GSM1089242", "GSM1089235", "GSM1089236"]


expression = pd.read_csv("data/processed/expr_log2.csv", index_col=0)
meta = pd.read_csv("data/processed/meta.csv", index_col=0)
detp = pd.read_csv("data/processed/detection_pvals.csv", index_col=0)
assert all(meta.index == expression.columns)
assert all(meta.index == detp.columns)
assert expression.min().min() > 0 
print("matrix", expression.shape)

is_pe = meta["diagnosis"] == "EOPET"
colors = np.where(is_pe, "#c0392b", "#2c7fb8") # Preeclampsia RED, Control BLUE
short = [g.replace("GSM10892", "..") for g in meta.index]

#__________________SAMPLE-SAMPLE PEARSON CORRELATION__________________
corr = expression.corr(method="pearson")
corr = corr.mask(np.eye(len(corr), dtype=bool))   # self-correlation -> NaN

mean_corr= corr.mean(axis=1)
print ("sample-sample Pearson correlation (excluding self)")
print("overall min:", round(np.nanmin(corr.values), 4), "median: ", round(np.nanmedian(corr.values), 4))
z = (mean_corr - mean_corr.mean()) / mean_corr.std()
flag_corr = z < -2
print("per-sample mean correlation:")
for g in meta.index:
    tag = "  <-- FLAG (z < -2)" if flag_corr[g] else ("  <-- watch" if g in WATCH else "")
    print(f"    {g}  {meta.loc[g, 'diagnosis']:8s} {mean_corr[g]:.4f}  z={z[g]:+.2f}{tag}")

fig, ax = plt.subplots(figsize=(7.5, 6.5))
im = ax.imshow(corr.values, cmap="viridis", vmin=np.nanmin(corr.values), vmax=1)
ax.set_xticks(range(16))
ax.set_yticks(range(16))
ax.set_xticklabels(short, rotation=90, fontsize=8)
ax.set_yticklabels(short, fontsize=8)
for t, c in zip(ax.get_xticklabels(), colors): t.set_color(c)
for t, c in zip(ax.get_yticklabels(), colors): t.set_color(c)
plt.colorbar(im, ax=ax, label="Pearson r")
ax.set_title("GSE44711  sample correlation, 19,601 probes (red = EOPET, blue = Control)")
fig.tight_layout()
fig.savefig("results/05_correlation_heatmap.png", dpi=150)
plt.close(fig)


#__________________PCA by SVD (probes centred, samples as points)__________________
X = expression.sub(expression.mean(axis=1), axis=0).T.values          # samples x probes, probe-centred
U, S, Vt = np.linalg.svd(X, full_matrices=False)
scores = U * S                                             # samples x components
var_expl = S**2 / (S**2).sum()
pcs = pd.DataFrame(scores[:, :4], index=meta.index, columns=["PC1", "PC2", "PC3", "PC4"])
print("\nvariance explained: " + "  ".join(f"PC{i+1} {v:.1%}" for i, v in enumerate(var_expl[:4])))
 
# what do the PCs track?
det_rate = (detp < 0.05).mean()
print("\nPC associations:")
for pc in ["PC1", "PC2", "PC3"]:
    t_diag = stats.ttest_ind(pcs.loc[is_pe.values, pc], pcs.loc[~is_pe.values, pc], equal_var=False)
    r_ga = stats.pearsonr(pcs[pc], meta["gest_age"])
    r_det = stats.pearsonr(pcs[pc], det_rate)
    print(f"  {pc}: EOPET vs Control p={t_diag.pvalue:.3g} | "
          f"gest_age r={r_ga[0]:+.2f} (p={r_ga[1]:.2g}) | "
          f"detection rate r={r_det[0]:+.2f} (p={r_det[1]:.2g})")

print("\nPC2 vs detection rate within each group:")
for grp in ["EOPET", "Control"]:
    m = (meta["diagnosis"] == grp).values
    r = stats.pearsonr(pcs.loc[m, "PC2"], det_rate[m])
    print(f"  {grp:8s} r={r[0]:+.2f}  p={r[1]:.2g}")
 
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(pcs["PC1"], pcs["PC2"], c=colors, s=70, edgecolor="k", linewidth=0.5)
for g, s in zip(meta.index, short):
    ax.annotate(s, (pcs.loc[g, "PC1"], pcs.loc[g, "PC2"]), fontsize=8,
                xytext=(4, 4), textcoords="offset points",
                fontweight="bold" if g in WATCH else "normal")
ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
ax.set_xlabel(f"PC1 ({var_expl[0]:.1%})"); ax.set_ylabel(f"PC2 ({var_expl[1]:.1%})")
ax.set_title("GSE44711  PCA, 19,601 probes (red = EOPET, blue = Control; bold = watch list)")
fig.tight_layout()
fig.savefig("results/05_pca.png", dpi=150)
plt.close(fig)

#__________________per-sample QC table__________________
qc = pd.DataFrame({
    "diagnosis": meta["diagnosis"],
    "gest_age": meta["gest_age"],
    "detection_rate": det_rate.round(3),
    "mean_corr": mean_corr.round(4),
    "corr_z": z.round(2),
    "PC1": pcs["PC1"].round(2),
    "PC2": pcs["PC2"].round(2),
    "watch": meta.index.isin(WATCH),
    "flag_low_corr": flag_corr,
})
print("\nper-sample QC:")
print(qc.to_string())
qc.to_csv("results/05_sample_qc.csv")
 
print("\nsaved: 05_correlation_heatmap.png, 05_pca.png, 05_sample_qc.csv")