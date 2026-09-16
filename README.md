# pe-GSE44711 — early-onset preeclampsia vs gestational-age-matched preterm controls

Third independent placental cohort in the CSE/H₂S preeclampsia re-analysis series
(after GSE75010, Affymetrix, and GSE25906, Illumina HumanWG-6 v2).

**Dataset.** GSE44711, Blair et al. 2013 (Robinson lab, UBC). Illumina HumanHT-12 v4
(GPL10558). Chorionic villi from 8 early-onset preeclampsia (EOPET) and 8 preterm
controls matched for gestational age and fetal sex. Companion 450K methylation
dataset from the same placentas: GSE44667.

**Question.** Does the PE transcriptomic signature reproduce in a third cohort, and is
CTH (cystathionine γ-lyase, the H₂S-producing enzyme) detectable in villous tissue?

## Pipeline

| script | does | writes |
|---|---|---|
| `01_load_metadata.py` | parse series matrix; diagnosis + gestational age; alignment seatbelt; **no transformation** | `data/processed/meta.csv`, `expr_linear.csv` |
| `02_detection_floor.py` | non-normalised export → detection calls (p < 0.05), probe filter (≥ 8/16 detected), log2, quantile normalisation; author-label → GSM mapping (hyphen trap) | `expr_log2.csv`, `detection_pvals.csv`, `probe_filter.csv` |
| `03_annotate.py` | Illumina's own GPL10558 table (3,270 unannotated probes vs 16,841 in GEO `.annot`); FLT1 single-probe assert | `annotation.csv` |
| `04_statistics.py` | Welch t + Mann-Whitney per probe, BH; `--all16` for the sensitivity run; robust set = padj < 0.05 in both | `results/04_stats.csv`, `04_stats_all16.csv`, `06_robust_set.csv` |
| `05_qc.py` | sample correlation, PCA, PC associations with diagnosis / GA / detection rate | `figures/05_*.png`, `results/05_sample_qc.csv` |
| `06_volcano.py` | volcano, 8 vs 7 primary, robust set marked | `figures/06_volcano.png` |
| `07_gene_panel.py` | per-sample panels: H₂S enzymes, canonical EOPE genes, stromal set | `figures/07_gene_panel.png` |

Run from the project root, in order. `data/` and large results are git-ignored.

## Decisions (and why)

- **Start from the non-normalised export, not the series matrix.** The deposited
  matrix is background-subtracted and quantile-normalised with 24.6% of cells negative
  (185,714 of 755,632) — `log2` would produce silent NaNs, and per-sample QC is erased
  by normalisation (every sample has 24.5–24.7% cells ≤ 0). The raw export has
  positive signal plus Illumina detection p-values.
- **Detection floor.** Probe kept if detection p < 0.05 in ≥ 8 of 16 samples (8 =
  smallest group). 19,601 of 47,227 probes kept. Sensitivity: at p < 0.01, [fill] probes kept.
  Floor on the log2-QN scale ≈ 6.5.
- **Positive control is LEP/HTRA4, not FLT1.** HT-12 v4 carries a single FLT1 probe
  (ILMN_1752307) at nt 5279 of NM_002019 — in the 3′ UTR of full-length FLT1, ~3 kb
  downstream of where the soluble isoforms sFlt1-i13 / sFlt1-e15a terminate. The
  platform cannot measure sFlt-1; the probe is at the floor (4/16 detected) and was
  filtered. Any re-analysis reporting FLT1 as unchanged on this platform has misread it.
- **GSM1089242 excluded from the primary analysis.** Pre-specified rule: mean
  sample-sample correlation z < −2 (observed z = −3.32, r = 0.892 vs 0.93–0.96 for all
  others). Two further independent criteria agree: lowest detection rate on the chip
  (0.380) and PC1 extreme (65 vs ≤ 42). Primary = 8 vs 7; full 8 vs 8 kept as sensitivity.
- **GSM1089240 retained.** Second-lowest detection rate (0.374), PC2 extreme, but
  correlation z = −0.56 does not meet the pre-specified rule; excluding it would be post hoc.
- **Probe-level analysis.** Probes of the same gene disagree in magnitude (LEP +3.4 vs
  +1.1; HTRA4 +3.1 vs +0.5); never averaged into genes.

## Results

**QC.** PC1 (35.7%) is not diagnosis (p = 0.5) — it is GSM1089242 and two other
lower-quality arrays. PC2 (20.5%) separates EOPET from Control (p = 1.5 × 10⁻⁵) with
every PE below zero and every control above; within the EOPET group PC2 is
independent of detection rate (r = −0.16), within Controls it is not (r = −0.75,
driven by ..40 and ..42).

**Differential expression** (Welch, BH):

| | 8 vs 7 (primary) | 8 vs 8 (all) |
|---|---|---|
| padj < 0.05 | 778 (352 up / 426 down) | 335 (179 / 156) |
| padj < 0.05 and \|log2FC\| > 1 | 98 | 26 |
| robust (padj < 0.05 in both) | [fill] | — |
| HTRA4 ILMN_2099277 | +3.23, padj 0.024 | +3.13, padj 0.041 |
| LEP ILMN_2207504 | +3.38, padj 0.056 | +3.37, padj 0.073 |
| CBS ILMN_1804735 | −0.29, padj 0.031 | −0.25, padj 0.085 |

One excluded array more than doubled the discovery set: at n = 8 conclusions are
sensitive to single samples, so paired numbers are reported throughout.
Mann-Whitney saturates at this n (minimum exact two-sided p = 2/12,870; all BH-adjusted
values 0.087–0.112) and is not usable for discovery; kept as a robustness column only.
The most significant Welch probes have |log2FC| < 1 (low-variance artefact); the
interpretable set is the 98 / 26 large-effect probes, and a limma moderated-t
cross-check is planned.

**Canonical EOPE signature reproduces.** HTRA4, LEP, FSTL3, PAPPA2, INHA, ENG, PRG2,
ABP1, EBI3, PSG4, ADAM12, PAPPA all up. A coordinated stromal/mesenchymal set is down
(PDPN −1.75, ANGPTL2 −1.33, IRX3 −1.09, CNN2 −1.01, OSR1 −0.99, HAND2 −0.72).

**H₂S enzymes.**
- **CTH: not expressed in villous tissue.** Three independent probes detected in
  1/8, 1/8, 0/8 EOPET and 0/8, 2/8, 1/8 Control (p < 0.05); scattered calls consistent
  with chance (48 tests). Third cohort, third platform, same result — consistent with
  CTH localising to stem villous artery smooth muscle rather than villous parenchyma.
- **CBS:** expressed just above floor, 16–18% lower in EOPET (padj 0.031 / 0.085).
- **MPST:** two tested probes flat; a third probe at the detection boundary detected in
  4/8 EOPET vs 1/8 Control (not tested).
- **SLC5A2 (SGLT2):** 0/8 and 0/8, 1/8 — absent, as expected for a kidney transporter.

**Sample-level observations.** Two EOPET placentas have LEP at floor: GSM1089235
(32.4 wk) and GSM1089236 (37.3 wk delivery). On unsupervised PCA, ..36 is the least
PE-like PE sample (PC2 = +1.0, at the group boundary). Probable non-canonical
molecular subclass (cf. Leavey et al. 2016); both retained.

## Limitations

- n = 8 vs 8 (7 after QC); discovery set sensitive to single samples.
- Fetal sex and matched-pair ID are not in the GEO deposit (paper states matched);
  labour status not recorded — preterm controls were delivered early for other
  reasons, EOPET often by indicated caesarean.
- Illumina detection-rate differs by group (EOPET ~0.45, Control ~0.42); PC2 partly
  tracks array quality within controls.
- Bulk villous tissue: a smooth-muscle-restricted transcript would be diluted even
  where present; the stromal-signature decrease in EOPET would dilute it further.

## Open items

- limma (`lmFit` → `eBayes`) cross-check on `expr_log2.csv`.
- FLT1 probe positions on GPL6102 (GSE25906) to close the three-platform sFlt-1 note.
- GSE44667 (450K, same placentas): CTH promoter methylation vs non-detection.