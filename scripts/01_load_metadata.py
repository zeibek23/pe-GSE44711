"""
01_load_metadata.py  --  GSE44711 (Illumina HT-12 v4; 8 EOPET vs 8 GA-matched preterm controls)

Loads the series matrix, parses sample metadata, checks alignment, and saves the
LINEAR (untransformed) expression matrix plus metadata to data/processed/.

Deliberately does NOT transform or compute fold changes: the deposited matrix is
linear, quantile-normalised, background-subtracted and contains negative values.
The floor decision (using detection p-values) is script 02.

Run from project root:  python scripts/01_load_metadata.py
"""
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------- paths
PROJECT = Path(__file__).resolve().parent.parent
RAW = PROJECT / "data" / "raw" / "GSE44711_series_matrix.txt"
OUT = PROJECT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)
assert RAW.exists(), f"missing input: {RAW}"

# ---------------------------------------------------------------- header block
with open(RAW) as f:
    lines = [line.rstrip("\n") for line in f]

table_begin = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
header = lines[:table_begin]
print("header lines:", len(header), "| table header is file line:", table_begin + 2)


def find_line(prefix, contains=""):
    """Return the ONE header line starting with `prefix` and containing `contains`."""
    hits = [l for l in header if l.startswith(prefix) and contains in l]
    assert len(hits) == 1, f"{prefix} {contains!r}: expected 1 line, found {len(hits)}"
    return hits[0]


def cells(line):
    """Tab-split a header line, drop the tag, strip quotes."""
    return [c.strip('"') for c in line.split("\t")[1:]]


def values(line):
    """For 'key: value' cells, keep the value."""
    return [c.split(": ", 1)[1] for c in cells(line)]


gsm_ids = cells(find_line("!Sample_geo_accession"))
diagnosis = values(find_line("!Sample_characteristics_ch1", '"condition:'))
gest_age = values(find_line("!Sample_characteristics_ch1", '"gestational age'))

meta = pd.DataFrame(
    {"diagnosis": diagnosis, "gest_age": gest_age},
    index=pd.Index(gsm_ids, name="gsm"),
)
meta["gest_age"] = meta["gest_age"].astype(float)  # everything from GEO is a string
assert meta["gest_age"].dtype == float
assert meta["diagnosis"].isin(["EOPET", "Control"]).all(), meta["diagnosis"].unique()
assert len(meta) == 16, len(meta)

print("\nmetadata:")
print(meta)
print("\ngestational age by group:")
print(meta.groupby("diagnosis")["gest_age"].agg(["count", "mean", "min", "max"]).round(2))

# ---------------------------------------------------------------- expression table
expr = pd.read_csv(
    RAW,
    sep="\t",
    skiprows=table_begin + 1,   # header row of the table
    skipfooter=1,               # drops '!series_matrix_table_end'
    engine="python",            # required by skipfooter
    index_col=0,
)
expr.index.name = "ID_REF"

print("\nexpression shape:", expr.shape)
assert expr.shape == (47227, 16), expr.shape
assert (expr.dtypes == float).all(), expr.dtypes.value_counts()
assert all(meta.index == expr.columns), "metadata/matrix misaligned"

# ---------------------------------------------------------------- scale & floor reconnaissance (report only)
print("min:", expr.min().min(), "| max:", expr.max().max(), "  -> linear scale")
print("negative cells:", int((expr < 0).sum().sum()))
print("zero cells:    ", int((expr == 0).sum().sum()))
print("\nfraction of cells <= 0 per sample (detection-rate proxy):")
print((expr <= 0).mean().round(3).to_string())

is_pe = meta["diagnosis"] == "EOPET"
print("\nPE samples:", int(is_pe.sum()), "| control samples:", int((~is_pe).sum()))

# ---------------------------------------------------------------- save
meta.to_csv(OUT / "meta.csv")
expr.to_csv(OUT / "expr_linear.csv")
print("\nsaved:", sorted(p.name for p in OUT.iterdir()))