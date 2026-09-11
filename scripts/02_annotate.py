import pandas as pd

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 250)

annot = pd.read_csv("data/GPL10558.annot", sep="\t", skiprows=28, skipfooter=1, engine="python")
results = pd.read_csv("results/01_log2fc.csv", index_col=0)

print ("Loaded annotation Data:")
print(annot.shape)
print(annot.columns.to_list())
print(annot.head(10))

print("Loaded Results data for Log2FC")
print(results.shape)
print(results.head(10))

columns_to_keep = ["ID", "Gene title", "Gene symbol", "Gene ID"]
annot = annot[columns_to_keep]
print("After selecting important columnds:")
print(annot.shape)
print(annot.head(10))

print(annot["Gene symbol"].isna().sum())
print(annot["ID"].dtype, results.index.dtype)
print(annot["ID"].duplicated().sum())


merged = results.merge(annot, left_index=True, right_on="ID", how="left")
print(merged.shape)
print(merged.head(10))
print(merged.columns.to_list())
print(merged["Gene symbol"].isna().sum())

top20 = merged.sort_values("abs_log2FC", ascending=False).head(20)
print(top20[["ID", "Gene symbol", "log2FC", "Gene title"]])

merged.to_csv("results/02_annotated_log2fc.csv", index=False)
print(merged[merged["Gene symbol"] == "SLC5A2"] [["ID", "log2FC", "mean_PE", "mean_nonPE"]])
print(merged[merged["Gene symbol"] == "CTH"] [["ID", "log2FC", "mean_PE", "mean_nonPE"]])
print(merged[merged["Gene symbol"] == "CBS"] [["ID", "log2FC", "mean_PE", "mean_nonPE"]])
print(merged[merged["Gene symbol"] == "MPST"] [["ID", "log2FC", "mean_PE", "mean_nonPE"]])
print(merged[merged["Gene symbol"] == "FLT1"] [["ID", "log2FC", "mean_PE", "mean_nonPE"]])

overall_mean = (merged["mean_PE"] + merged["mean_nonPE"]) / 2
print(overall_mean.describe())
print(overall_mean.shape)
print(overall_mean.head(5))


genes = ["SLC5A2", "LEP", "FLT1", "SLC5A1", "SLC2A1", "CTH", "CBS", "MPST"]
print(merged[merged["Gene symbol"].isin(genes)][["Gene symbol", "mean_PE", "mean_nonPE", "log2FC"]])