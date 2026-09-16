import pandas as pd
from pathlib import Path
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

expression = pd.read_csv("data/processed/expr_log2.csv", index_col=0)
meta = pd.read_csv("data/processed/meta.csv", index_col=0)
annotation = pd.read_csv("data/processed/annotation.csv", index_col=0)
detp = pd.read_csv("data/processed/detection_pvals.csv", index_col=0)
assert all(meta.index == detp.columns), "detection table misaligned"

print(expression.shape)
print(expression.head(10))
print(type(meta))
print(len(meta))
print(all(meta.index == expression.columns))
#---------------------GROUP---------------------
is_pe = meta["diagnosis"] == "EOPET"
group_pe = expression.loc[:, is_pe.values]
group_control = expression.loc[:, (~is_pe).values]
assert group_pe.shape[1] == 8 and group_control.shape[1] == 8, (group_pe.shape, group_control.shape)
#---------------------TESTS---------------------
tstat, pval = stats.ttest_ind(group_pe, group_control, axis=1, equal_var=False)
assert len(pval) == len(expression) and not np.isnan(pval).any()
 
mwu = stats.mannwhitneyu(group_pe, group_control, axis=1, alternative="two-sided")
 
results = pd.DataFrame(
    {
        "mean_PE": group_pe.mean(axis=1),
        "mean_Control": group_control.mean(axis=1),
        "t": tstat,
        "pval": pval,
        "pval_mwu": mwu.pvalue,
    },
    index=expression.index,
)
results["log2FC"] = results["mean_PE"] - results["mean_Control"]
results["padj"] = multipletests(results["pval"], method="fdr_bh")[1]
results["padj_mwu"] = multipletests(results["pval_mwu"], method="fdr_bh")[1]
results["neg_log10_padj"] = -np.log10(results["padj"])

#---------------------ANNOTATE---------------------
results = results.join(annotation[["symbol", "probe_start"]], how="left")
assert len(results) == len(expression)
frac_annotated = results["symbol"].notna().mean()
print(f"Probes with a symbol: {frac_annotated:.1%}")
assert frac_annotated > 0.90 



#---------------------POSITIVE CONTROL---------------------
lep = results.loc["ILMN_2207504"]
assert lep["symbol"] == "LEP", lep["symbol"]
assert lep["log2FC"] > 1 and lep["pval"] < 0.01, lep[["log2FC", "pval"]]
htra4 = results.loc["ILMN_2099277"]
assert htra4["symbol"] == "HTRA4", htra4["symbol"]
assert htra4["log2FC"] > 1 and htra4["padj"] < 0.05, htra4[["log2FC", "padj"]]
print("controls OK | LEP   log2FC %+.2f  p=%.2g  padj=%.3f  padj_mwu=%.3f"
      % (lep["log2FC"], lep["pval"], lep["padj"], lep["padj_mwu"]))
print("            | HTRA4 log2FC %+.2f  p=%.2g  padj=%.3f  padj_mwu=%.3f"
      % (htra4["log2FC"], htra4["pval"], htra4["padj"], htra4["padj_mwu"]))
 
# FLT1 single probe is isoform-blind and was floored out; it must not be here.
assert "ILMN_1752307" not in results.index, "FLT1 probe should have been floored out"


#---------------------SUMMARY---------------------
sig = results["padj"] < 0.05
sig_mwu = results["padj_mwu"] < 0.05
print("\nprobes tested:", len(results))
print("Welch  padj < 0.05:", int(sig.sum()),
      "| up:", int((sig & (results["log2FC"] > 0)).sum()),
      "| down:", int((sig & (results["log2FC"] < 0)).sum()),
      "| padj < 0.10:", int((results["padj"] < 0.10).sum()))
print("MWU    padj < 0.05:", int(sig_mwu.sum()),
      "| in both:", int((sig & sig_mwu).sum()))
print("Welch padj < 0.05 & |log2FC| > 1:", int((sig & (results["log2FC"].abs() > 1)).sum()))
 
results = results.sort_values("padj")
cols = ["symbol", "mean_PE", "mean_Control", "log2FC", "pval", "padj", "padj_mwu"]
print("\ntop 20 by Welch padj:")
print(results[cols].head(20).round(3).to_string())
 
#---------------------HYDROGEN SULFIDE ENZYMES & SGLT2---------------------
panel = ["CTH", "CBS", "MPST", "SLC5A2"]
print("\nH2S enzymes + SGLT2 -- probes that passed the detection floor:")
tested = results[results["symbol"].isin(panel)]
print(tested[cols].round(3).to_string() if len(tested) else "  (none)")
 
print("\nH2S enzymes + SGLT2 -- probes floored out (reported by detection, not fold change):")
for gene in panel:
    for pid in annotation.index[annotation["symbol"] == gene]:
        if pid in results.index or pid not in detp.index:
            continue
        n_pe = int((detp.loc[pid, is_pe.values] < 0.05).sum())
        n_ct = int((detp.loc[pid, (~is_pe).values] < 0.05).sum())
        print(f"  {gene:7s} {pid}  detected in {n_pe}/8 PE, {n_ct}/8 Control  -> not tested")

 
results.to_csv("results/04_stats.csv")
print("results saved!")