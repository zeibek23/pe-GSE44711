"""
02_detection_floor.py  --  GSE44711: detection floor, log2, quantile normalisation

Starts from the NON-normalised Illumina export (raw AVG_Signal + Detection Pval per sample),
not from the series matrix: the series matrix is background-subtracted and quantile-normalised
with ~25% of cells negative, which makes log2 unsafe and hides per-sample QC.

Steps
  1. load raw signal + detection p-values, split the paired columns
  2. map the authors' sample labels (PM12, PL-113, ...) to GSM IDs via !Sample_description
     -- TRAP: series matrix writes 'PL113', the export writes 'PL-113'; hyphens are stripped
  3. detection call per probe per sample; per-sample detection rate (real QC, unlike script 01)
  4. probe filter: keep probes detected in >= MIN_DETECTED samples
  5. log2 -> quantile normalise -> apply filter
  6. positive / negative controls on the final matrix
  7. save data/processed/expr_log2.csv, detection_pvals.csv, probe_filter.csv

Run from project root:  python scripts/02_detection_floor.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- DECISIONS (yours to defend)
DET_P = 0.05          # detection p-value threshold. Illumina convention is 0.05; 0.01 is stricter.
MIN_DETECTED = 8      # keep a probe if detected in >= this many of 16 samples.
                      # 8 = smallest group size: a gene expressed in all of one group but none
                      # of the other survives; a gene detected in 3 random samples does not.
# ---------------------------------------------------------------- paths
PROJECT = Path(__file__).resolve().parent.parent
NONNORM = PROJECT / "data" / "raw" / "GSE44711_non-normalized_data.txt"
SERIES = PROJECT / "data" / "raw" / "GSE44711_series_matrix.txt"
OUT = PROJECT / "data" / "processed"
for p in (NONNORM, SERIES, OUT / "meta.csv", OUT / "annotation.csv", OUT / "expr_linear.csv"):
    assert p.exists(), f"missing: {p}"

meta = pd.read_csv(OUT / "meta.csv", index_col=0)
annot = pd.read_csv(OUT / "annotation.csv", index_col=0)
assert len(meta) == 16

# ---------------------------------------------------------------- 1. raw export
raw = pd.read_csv(NONNORM, sep="\t", index_col=0)
raw.index.name = "ID_REF"
print("raw export shape:", raw.shape)

signal = raw.filter(like=".AVG_Signal").copy()
detp = raw.filter(like=".Detection Pval").copy()
signal.columns = signal.columns.str.replace(".AVG_Signal", "", regex=False)
detp.columns = detp.columns.str.replace(".Detection Pval", "", regex=False)

assert list(signal.columns) == list(detp.columns), "signal/p-value columns do not pair up"
assert signal.shape[1] == 16, signal.shape
assert (signal > 0).all().all(), "non-positive raw signal found; log2 would be unsafe"
assert ((detp >= 0) & (detp <= 1)).all().all(), "detection p-values outside [0, 1]"
print("raw signal min:", round(signal.min().min(), 2), "| max:", round(signal.max().max(), 2))

# same probe set as the deposited matrix?
series_ids = pd.read_csv(OUT / "expr_linear.csv", index_col=0, usecols=[0]).index
only_raw = raw.index.difference(series_ids)
only_series = series_ids.difference(raw.index)
print("probes only in raw export:", len(only_raw), "| only in series matrix:", len(only_series))
assert len(only_series) == 0, "series matrix has probes the raw export lacks"

# ---------------------------------------------------------------- 2. author label -> GSM
with open(SERIES) as f:
    header = [l.rstrip("\n") for l in f if l.startswith("!Sample_")]


def find_line(prefix, contains=""):
    hits = [l for l in header if l.startswith(prefix) and contains in l]
    assert len(hits) == 1, f"{prefix} {contains!r}: found {len(hits)}"
    return hits[0]


def cells(line):
    return [c.strip('"') for c in line.split("\t")[1:]]


gsm_ids = cells(find_line("!Sample_geo_accession"))
labels = [c.split(": ", 1)[1] for c in cells(find_line("!Sample_description"))]


def norm(s):
    """'PL-113' and 'PL113' are the same sample."""
    return s.replace("-", "").strip()


label_to_gsm = {norm(l): g for l, g in zip(labels, gsm_ids)}
unmapped = [c for c in signal.columns if norm(c) not in label_to_gsm]
assert not unmapped, f"unmapped sample labels: {unmapped}"

new_cols = [label_to_gsm[norm(c)] for c in signal.columns]
assert len(set(new_cols)) == 16, "two labels mapped to the same GSM"
signal.columns = new_cols
detp.columns = new_cols

# order columns exactly as meta, then seatbelt
signal = signal[meta.index]
detp = detp[meta.index]
assert all(signal.columns == meta.index) and all(detp.columns == meta.index)
print("label -> GSM mapping:")
print(pd.Series(label_to_gsm).to_string())

# ---------------------------------------------------------------- 3. detection calls
detected = detp < DET_P
print(f"\nper-sample detection rate (p < {DET_P}):")
print(detected.mean().round(3).to_string())

n_det = detected.sum(axis=1)
print("\nprobes by number of samples detected in:")
print(n_det.value_counts().sort_index().to_string())

# ---------------------------------------------------------------- 4. probe filter
keep = n_det >= MIN_DETECTED
print(f"\nprobes kept (detected in >= {MIN_DETECTED}/16): {int(keep.sum())} of {len(keep)}")

# ---------------------------------------------------------------- 5. log2 + quantile normalisation
log2 = np.log2(signal)


def quantile_normalize(df):
    """Classic quantile normalisation: every column gets the mean sorted distribution."""
    sorted_mean = np.sort(df.values, axis=0).mean(axis=1)
    ranks = df.rank(method="average")
    positions = np.arange(1, len(sorted_mean) + 1)
    return ranks.apply(lambda r: np.interp(r, positions, sorted_mean))


qn = quantile_normalize(log2)
assert qn.shape == log2.shape
print("\nlog2 raw    median per sample:", log2.median().round(2).to_dict())
print("log2 QN     median per sample:", qn.median().round(2).to_dict())

expr_log2 = qn.loc[keep]
print("\nfinal matrix:", expr_log2.shape, "| min:", round(expr_log2.min().min(), 2),
      "| max:", round(expr_log2.max().max(), 2))

# ---------------------------------------------------------------- 6. controls
is_pe = meta["diagnosis"] == "EOPET"


def probes_for(gene):
    ids = annot.index[annot["symbol"] == gene]
    return [i for i in ids if i in raw.index]


print("\ncontrols  (n_detected/16, kept, mean log2 PE, mean log2 Ctrl, diff):")
for gene in ["LEP", "PAPPA2", "INHA", "FSTL3", "HTRA4", "ENG", "FLT1", "CTH", "CBS", "MPST", "SLC5A2"]:
    for pid in probes_for(gene):
        row = qn.loc[pid]
        d = row[is_pe].mean() - row[~is_pe].mean()
        print(f"{gene:7s} {pid}  {int(n_det[pid]):2d}/16  {'KEPT ' if keep[pid] else 'DROP '} "
              f"{row[is_pe].mean():5.2f} {row[~is_pe].mean():5.2f}  {d:+5.2f}")

# positive control must survive the filter; FLT1 is NOT a valid control on this platform
lep = probes_for("LEP")
assert any(keep[p] for p in lep), "LEP filtered out -- floor is wrong"

# CTH: the definitive statement, per probe per sample
cth = probes_for("CTH")
print("\nCTH detection p-values (rows = probes, cols = samples):")
print(detp.loc[cth].round(3).to_string())
print("CTH samples detected per probe:", n_det[cth].to_dict())

# ---------------------------------------------------------------- 7. save
expr_log2.to_csv(OUT / "expr_log2.csv")
detp.to_csv(OUT / "detection_pvals.csv")
pd.DataFrame({"n_detected": n_det, "kept": keep}).to_csv(OUT / "probe_filter.csv")
print("\nsaved:", sorted(p.name for p in OUT.iterdir()))