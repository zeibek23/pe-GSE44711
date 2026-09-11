"""
03_annotate.py  --  probe -> symbol table for GPL10558 (Illumina HT-12 v4)

Uses Illumina's own platform table (GEO 'Download full table', GPL10558-50081.txt),
NOT the GEO .annot file: the .annot re-annotates via GenBank and leaves ~16.8k probes
without a symbol; Illumina's table also carries Probe_Start, which is what showed the
single FLT1 probe (ILMN_1752307) sits at nt 5279 in the 3' UTR of full-length NM_002019
and is therefore blind to the sFlt-1 isoforms (sFlt1-i13 / sFlt1-e15a).

Output: data/processed/annotation.csv   (index = ILMN probe ID)
This script is independent of the expression values; merge on ID_REF in the stats script.

Run from project root:  python scripts/03_annotate.py
"""
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------- paths
PROJECT = Path(__file__).resolve().parent.parent
GPL = PROJECT / "data" / "raw" / "GPL10558-50081.txt"
OUT = PROJECT / "data" / "processed"
assert GPL.exists(), f"missing input: {GPL}"
assert (OUT / "expr_linear.csv").exists(), "run 01_load_metadata.py first"

# ---------------------------------------------------------------- load Illumina table
# The file starts with '#' comment lines; the header is the first line beginning 'ID<tab>'.
with open(GPL) as f:
    header_row = next(i for i, l in enumerate(f) if l.startswith("ID\t"))
print("header row (0-based):", header_row)

annot = pd.read_csv(GPL, sep="\t", skiprows=header_row, dtype=str, low_memory=False)
print("raw table shape:", annot.shape)

# keep only real probe rows (guards against any trailing non-probe lines)
annot = annot[annot["ID"].str.startswith("ILMN_", na=False)]

keep = {
    "ID": "ID",
    "Symbol": "symbol",
    "RefSeq_ID": "refseq",
    "Probe_Start": "probe_start",
    "Definition": "definition",
}
annot = annot[list(keep)].rename(columns=keep).set_index("ID")
assert annot.index.is_unique, "duplicate probe IDs in platform table"
print("probes:", len(annot), "| without symbol:", int(annot["symbol"].isna().sum()))

# ---------------------------------------------------------------- platform sanity checks
# FLT1: exactly one probe on HT-12 v4, in the 3' UTR of the full-length transcript.
flt1 = annot.index[annot["symbol"] == "FLT1"].tolist()
assert flt1 == ["ILMN_1752307"], flt1
print("FLT1 probes:", flt1, "| probe_start:", annot.loc["ILMN_1752307", "probe_start"])

# LEP must exist (it is the positive control on this platform, FLT1 is not).
lep = annot.index[annot["symbol"] == "LEP"].tolist()
assert len(lep) >= 1, "no LEP probe found"
print("LEP probes:", lep)

for gene in ["CTH", "CBS", "MPST", "PAPPA2", "INHA", "FSTL3", "HTRA4", "ENG", "SLC5A2"]:
    ids = annot.index[annot["symbol"] == gene].tolist()
    print(f"{gene:7s} {len(ids)} probe(s): {ids}")

# ---------------------------------------------------------------- coverage vs expression matrix
expr_ids = pd.read_csv(OUT / "expr_linear.csv", index_col=0, usecols=[0]).index
missing = expr_ids.difference(annot.index)
print("expression probes with no annotation row:", len(missing))
assert len(missing) == 0, list(missing[:10])

# ---------------------------------------------------------------- save
annot.to_csv(OUT / "annotation.csv")
print("saved:", OUT / "annotation.csv")