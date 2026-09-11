import pandas as pd

pe_pathway = "data/raw/GSE44711_series_matrix.txt"

pe_data = pd.read_csv(pe_pathway, sep="\t", skiprows=69, skipfooter=1, engine="python" )

print("Loading new Preeclampsia data: ")
print(pe_data.head(10))
print(pe_data.shape)
print(pe_data.columns.to_list())

pe_data = pe_data.set_index("ID_REF")
print(pe_data.dtypes.value_counts())
print(pe_data.min().min(), pe_data.max().max())

with open("data/raw/GSE44711_series_matrix.txt") as file:
    meta_lines = [line.rstrip("\n") for line in file][:69]
print(len(meta_lines))

gsm_line = None
for line in meta_lines:
    if line.startswith("!Sample_geo_accession"):
        gsm_line = line

diag_line = None
for line in meta_lines:
    if line.startswith("!Sample_characteristics_ch1") and ('"condition:'):
        diag_line=line
print(diag_line[0:100])


matches = [m for m in meta_lines
            if m.startswith("!Sample_characteristics_ch1") and ('"condition:' in m) ]
print(len(matches))
diag_line = matches[0]
print(diag_line[0:150])


gsm_ids = [p.strip('"') for p in gsm_line.split("\t")[1:]]
diagnoses = [p.strip('"').split(": ", 1)[1] for p in diag_line.split("\t")[1:]]
print(len(gsm_ids), len(diagnoses))

diagnosis = pd.Series(diagnoses, index=gsm_ids, name="diagnosis")
print(diagnosis.value_counts())
print(len(diagnosis))
print(pe_data.shape[1])
print(all(diagnosis.index == pe_data.columns))

is_pe = diagnosis == "EOPET"
print(is_pe.sum())
print(~is_pe.sum())

pe_samples = pe_data.loc[:, is_pe.values]
non_pe_samples = pe_data.loc[:, (~is_pe.values)]
print(pe_samples.shape)
print(non_pe_samples.shape)

mean_pe = pe_samples.mean(axis=1)
mean_non_pe = non_pe_samples.mean(axis=1)
print(mean_pe.shape)

log2fc = mean_pe - mean_non_pe
results = pd.DataFrame({
    "mean_PE": mean_pe,
    "mean_nonPE": mean_non_pe,
    "log2FC": log2fc
})

results["abs_log2FC"] = results["log2FC"].abs()
top = results.sort_values("abs_log2FC", ascending=False).head(25)
print(top)

import os
os.makedirs("results", exist_ok=True)
results.to_csv("results/01_log2fc.csv")
print("Saved the following:", results.shape)

pe_data.to_csv("results/01_expression.csv")
diagnosis.to_csv("results/01_diagnosis.csv")
print("Saved the pe_data and diagnosis as CSV")