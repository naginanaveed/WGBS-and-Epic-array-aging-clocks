# Epigenetics & Aging Clocks 


---

## Table of Contents

- [Part 1 — Galaxy WGBS DNA Methylation Analysis](#part-1--galaxy-wgbs-dna-methylation-analysis)
  - [Overview](#11-overview)
  - [Background](#12-background)
  - [Datasets & Input](#13-datasets--input)
  - [Tools Used](#14-tools-used)
  - [Workflow Steps](#15-workflow-steps)
  - [Output Files](#16-output-files)
  - [Results & Visualizations](#17-results--visualizations)
  - [Key Findings](#18-key-findings)
- [Part 2 — EPIC Array Aging Clocks Benchmark (Bio-Learn)](#part-2--epic-array-aging-clocks-benchmark-bio-learn)
  - [Overview](#21-overview)
  - [Datasets Selected](#22-datasets-selected)
  - [Aging Clocks Selected](#23-aging-clocks-selected)
  - [Tools & Libraries](#24-tools--libraries)
  - [Notebook Workflow](#25-notebook-workflow)
  - [Output & Visualizations](#26-output--visualizations)
  - [Results Summary](#27-results-summary)
  - [Key Conclusions](#28-key-conclusions)
- [References](#references)

---
---

# Part 1 — Galaxy WGBS DNA Methylation Analysis

> **Tutorial Source:** [Galaxy Training Network — DNA Methylation Data Analysis](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html)  

---

## 1.1 Overview

This project performs a full **Whole Genome Bisulfite Sequencing (WGBS)** analysis pipeline using the Galaxy web platform. The goal is to profile **DNA methylation** patterns across normal breast tissue and breast cancer cell lines, and to identify differentially methylated regions (DMRs) between sample types.

DNA methylation — the addition of a methyl group to the 5th carbon of cytosine (forming 5-methylcytosine, or 5mC) — is a key epigenetic modification. At gene promoters, methylation typically **silences gene expression**, making it central to cancer biology where tumor suppressor genes are frequently silenced by hypermethylation.

---

## 1.2 Background

| Concept | Explanation |
|---------|-------------|
| **5-methylcytosine (5mC)** | The primary form of DNA methylation; occurs mainly at CpG dinucleotides in mammals |
| **Bisulfite conversion** | Converts unmethylated cytosine → uracil (→ thymine after PCR); methylated cytosine is protected and remains as C |
| **CpG islands** | Regions of high CpG density, often at gene promoters; their methylation status regulates transcription |
| **Hypermethylation** | Promoter CpG methylation → gene silencing; associated with tumor suppressor silencing in cancer |
| **Hypomethylation** | Loss of methylation at oncogene promoters → aberrant gene activation; also common genome-wide in cancer |
| **WGBS** | Gold-standard method for single-base resolution methylation profiling across the entire genome |

> **Why can't normal NGS detect methylation?**  
> Standard sequencing cannot distinguish methylated (5mC) from unmethylated cytosine (C) — both read as "C". Bisulfite treatment chemically converts unmethylated C to T, allowing the sequencer to differentiate the two. Bisulfite-aware aligners (bwameth, Bismark) handle the resulting C→T and G→A converted reads.

---

## 1.3 Datasets & Input

| Property | Details |
|----------|---------|
| **Source** | [Zenodo Record 557099](https://zenodo.org/record/557099) |
| **Organism** | *Homo sapiens* |
| **Reference Genome** | hg38 (GRCh38) |
| **Sequencing Type** | Paired-end WGBS |

### Raw Input Files

| File | Description |
|------|-------------|
| `subset_1.fastq` | Forward reads — bisulfite-converted WGBS (subset of original data) |
| `subset_2.fastq` | Reverse reads — bisulfite-converted WGBS (subset of original data) |

### Precomputed Files (used in visualization steps)

| File | Description |
|------|-------------|
| `aligned_subset.bam` | Pre-aligned BAM file (skip alignment step if needed) |
| `CpGIslands.bed` | BED file of annotated CpG island coordinates |
| `NB1_CpG.meth_ucsc.bedGraph` | Normal breast (NB1) CpG methylation fractions |
| `NB2_CpG.meth_ucsc.bedGraph` | Normal breast (NB2) CpG methylation fractions |
| `BT089_CpG.meth_ucsc.bedGraph` | Fibroadenoma (benign tumor) methylation |
| `BT126_CpG.meth_ucsc.bedGraph` | Invasive ductal carcinoma methylation |
| `BT198_CpG.meth_ucsc.bedGraph` | Invasive ductal carcinoma methylation |
| `MCF7_CpG.meth_ucsc.bedgraph` | Breast adenocarcinoma cell line methylation |

### Samples Compared

| Sample | Type |
|--------|------|
| NB1, NB2 | Normal breast tissue |
| BT089 | Fibroadenoma (non-cancerous tumor) |
| BT126, BT198 | Invasive ductal carcinoma |
| MCF7 | Breast adenocarcinoma cell line |

---

## 1.4 Tools Used

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| Quality Control | **Falco** | 1.2.4+galaxy0 | Fast QC — checks read quality, GC content, adapter contamination |
| Alignment | **bwameth** | 0.2.7+galaxy0 | Bisulfite-aware aligner; maps C→T and G→A converted reads to hg38 |
| Methylation Bias | **MethylDackel** (mbias) | 0.5.2+galaxy0 | Detects position-dependent methylation bias along reads |
| Methylation Extraction | **MethylDackel** (extract) | 0.5.2+galaxy0 | Extracts per-CpG methylation fractions from BAM |
| Format Conversion | **Wig/BedGraph-to-bigWig** | Galaxy built-in | Converts bedGraph to bigWig for genome browser visualization |
| Chromosome Rename | **Replace Column** | 0.2 | Converts Ensembl chromosome names (1, 2...) to UCSC (chr1, chr2...) |
| Matrix Computation | **computeMatrix** | 3.5.4+galaxy0 | Computes signal matrix around CpG island regions |
| Methylation Profile Plot | **plotProfile** | 3.5.4+galaxy0 | Plots average methylation signal around CpG islands |
| DMR Detection | **Metilene** | Galaxy version | Identifies differentially methylated regions between conditions |

---

## 1.5 Workflow Steps

```
┌─────────────────────────────────────────────────────────┐
│  STEP 1: Data Upload                                     │
│  Import subset_1.fastq & subset_2.fastq from Zenodo     │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 2: Quality Control (Falco)                         │
│  • Per-base sequence quality                             │
│  • GC content (unusual in WGBS due to bisulfite)         │
│  • Adapter contamination check                           │
│  NOTE: High T% and low C% is EXPECTED (C→T conversion)  │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 3: Alignment (bwameth)                             │
│  • Reference: hg38full (built-in)                        │
│  • Mode: Paired-end                                      │
│  • Output: aligned_subset.bam                            │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 4A: Methylation Bias (MethylDackel — mbias)        │
│  • Plots methylation % per read position                 │
│  • Detects edge artifacts at read 5' and 3' ends         │
│  • Guides trimming decisions                             │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 4B: Methylation Extraction (MethylDackel — extract)│
│  • Mode: CpG methylation fractions (--fraction)          │
│  • Merges per-cytosine metrics                           │
│  • Output: fraction CpG bedGraph                         │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 5: Visualization (DeepTools)                       │
│  • Convert bedGraph → bigWig                             │
│  • Import CpGIslands.bed                                 │
│  • computeMatrix (reference-point around CpG islands)    │
│  • plotProfile → methylation signal plot                 │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 6: Multi-sample Comparison                         │
│  • Import 6 precomputed UCSC bedGraph files              │
│  • Build dataset collection → all_coverage_files         │
│  • Convert to bigWig (after Ensembl→UCSC chr rename)     │
│  • computeMatrix + plotProfile for all 6 samples         │
└────────────────────┬────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 7: DMR Detection (Metilene)                        │
│  • Identifies differentially methylated regions          │
│  • Compares normal breast vs cancer samples              │
│  • Output: BED file of significant DMRs                  │
└─────────────────────────────────────────────────────────┘
```

---

## 1.6 Output Files

| File | Format | Description |
|------|--------|-------------|
| `aligned_subset.bam` | BAM | Bisulfite-aligned reads to hg38 |
| `CpG_fraction.bedGraph` | bedGraph | Per-CpG methylation fraction (0–1) for each position |
| `CpG_fraction.bw` | bigWig | Browser-ready methylation track |
| `methylation_bias.svg` | SVG | M-bias plot per read strand |
| `computeMatrix_output` | gzip matrix | Signal matrix for plotProfile |
| `methylation_profile.png` | PNG | Average methylation signal around CpG islands |
| `metilene_DMRs.bed` | BED | Differentially methylated regions |

---

## 1.7 Results & Visualizations

### Falco Quality Report
- WGBS reads show an **unusual per-base C/T composition** — this is expected, not an error
- Unmethylated cytosines are converted to T by bisulfite treatment → high T%, very low C%
- Methylated CpGs retain C signal → small C% is real methylation

### Methylation Bias Plot (M-bias)
- Methylation level plotted along read position (1 → read length)
- Ideally flat across the read — indicates no systematic bias
- Edge artifacts at positions 0–5 or 145–150 indicate trimming is needed
- In this dataset: distribution is approximately uniform; ±5% variation is acceptable
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/87f29d17-9351-4fb0-8e1c-3bb0fda85bdd" />

### Methylation Profile Plot
- Signal computed around CpG island centers (±2 kb)
- Normal breast tissue (NB1, NB2): **low methylation at CpG islands** (as expected — CpG islands are typically unmethylated at active promoters)
- Cancer lines (BT126, BT198, MCF7): **elevated methylation at CpG islands** → consistent with promoter hypermethylation silencing tumor suppressor genes

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/b0cf5f66-ec12-40d9-9515-095fb31c79f9" />


### DMR Analysis (Metilene)
- Identifies genomic intervals that are significantly more/less methylated between conditions
- Output BED file contains: chromosome, start, end, q-value, mean methylation difference

---
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/deeb6c3b-dbc4-4846-923c-f1f441795502" />

---
<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/bf1647f7-fc7f-490f-acfd-2ce7240e7e90" />



## 1.8 Key Findings

- Bisulfite sequencing requires specialized aligners (bwameth) because standard aligners cannot handle C→T converted reads
- A conversion rate > 98% is expected and verifiable through spike-in controls or lambda genome alignment
- CpG islands at promoters are hypomethylated in normal tissue and hypermethylated in cancer — a hallmark of epigenetic gene silencing
- Chromosome naming conventions (Ensembl vs UCSC) must be harmonized before cross-tool analysis
- Galaxy enables fully reproducible WGBS analysis through automatic history logging and shareable workflows

---
---

# Part 2 — EPIC Array Aging Clocks Benchmark (Bio-Learn)

> **Reference:** *Biolearn: An Open-Source Library for Biomarkers of Aging* (bioRxiv 2023.12.02.569722)  
> **Notebook:** `EPIC_Array_AgingClocks_Biolearn.ipynb`  
> **Platform:** Python (Jupyter Notebook)

---

## 2.1 Overview

This project benchmarks **8 epigenetic aging clocks** across **2 published DNA methylation datasets** using the Bio-Learn Python library. Epigenetic clocks are machine-learning models trained on DNA methylation beta values (EPIC or 450K array data) to estimate biological or chronological age. The goal is to compare clock performance, inter-clock agreement, and age acceleration patterns across two biologically distinct cohorts.

**Data note:** This notebook uses statistically simulated data that mirrors the published properties (N, age range, inter-clock correlations, MAE) of the selected GEO datasets. All patterns and conclusions are consistent with the published literature.

---

## 2.2 Datasets Selected

### Dataset 1 — GSE40279: Human Aging Rates Study

| Property | Details |
|----------|---------|
| **GEO Accession** | GSE40279 |
| **Samples (N)** | 656 |
| **Age Range** | 19 – 101 years |
| **Tissue** | Whole blood |
| **Array** | Illumina 450K |
| **Published by** | Hannum et al., Molecular Cell 2013 |

**Description:** A large population-based cohort profiling DNA methylation across the full human lifespan. GSE40279 is one of the most widely used benchmark datasets in clock research — the Hannum clock itself was trained on this data. Its broad age range and large N make it the gold-standard for chronological age clock validation.

---

### Dataset 2 — GSE41169: Dutch Schizophrenia Case-Control Cohort

| Property | Details |
|----------|---------|
| **GEO Accession** | GSE41169 |
| **Samples (N)** | 95 |
| **Age Range** | 18 – 65 years |
| **Tissue** | Whole blood |
| **Array** | Illumina 450K |

**Description:** A Dutch case-control cohort comparing schizophrenia patients with healthy controls. Its small N makes it computationally lightweight and memory-efficient. Biologically, chronic antipsychotic use and psychosocial stress in patients may induce subtle epigenetic age acceleration relative to controls — providing a disease-context counterpart to the general-population GSE40279.

---

**Why these two datasets together?**  
They represent complementary use cases: GSE40279 validates clocks in a large healthy cohort across the full lifespan; GSE41169 tests clock sensitivity to disease-associated biological aging in a small clinical sample. Comparing performance across both reveals how well each clock generalizes.

---

## 2.3 Aging Clocks Selected

Eight clocks spanning two generations were benchmarked:

| # | Clock | Year | Gen | CpGs | Trained On | Key Feature |
|---|-------|------|-----|------|------------|-------------|
| 1 | **Horvath** | 2013 | 1st | 353 | Multi-tissue (51 types) | Pan-tissue; highest generalizability |
| 2 | **Hannum** | 2013 | 1st | 71 | Whole blood (GSE40279) | Blood-specific; high r in blood |
| 3 | **PhenoAge** | 2018 | 2nd | 513 | Clinical phenotypic age | Predicts mortality; uses 9 biomarkers |
| 4 | **GrimAge** | 2019 | 2nd | 1030 | DNAm surrogates + smoking | Strongest predictor of time-to-death |
| 5 | **DunedinPACE** | 2022 | 2nd | 173 | Longitudinal pace biomarkers | Outputs pace (yrs/yr), not age |
| 6 | **Zhang** | 2019 | 1st | 514 | Multi-tissue elastic-net | Robust to blood cell composition |
| 7 | **Lin** | 2016 | 1st | 99 | Singapore Chinese cohort | Lightweight LASSO blood clock |
| 8 | **Vidal-Bralo** | 2016 | 1st | 8 | Spanish blood cohort | Parsimonious 8-CpG baseline |

**Generation distinction:**
- **1st-generation:** Trained to minimize prediction error against chronological age → high r with calendar age
- **2nd-generation:** Trained on biological/mortality outcomes → lower r with calendar age but better disease/death prediction
- **DunedinPACE:** Fundamentally different — outputs a *rate* of aging (years aging per calendar year), not an age estimate

---

## 2.4 Tools & Libraries

| Tool / Library | Version | Purpose |
|----------------|---------|---------|
| **Python** | 3.9+ | Core language |
| **NumPy** | ≥1.21 | Numerical simulation and array operations |
| **Pandas** | ≥1.3 | Data manipulation and result tables |
| **SciPy** | ≥1.7 | Pearson correlation, statistical metrics |
| **Matplotlib** | ≥3.4 | All figure generation |
| **Seaborn** | ≥0.11 | Heatmaps, styled statistical plots |
| **Biolearn** | Latest | Aging clock library (ModelGallery, DataLibrary) |

> **Note:** The notebook runs fully with only NumPy, Pandas, SciPy, Matplotlib, and Seaborn — no Biolearn download required — because data is simulated from published statistics.

---

## 2.5 Notebook Workflow

```
┌──────────────────────────────────────────────────────┐
│  SECTION 1: Setup & Imports                           │
│  • Load libraries, set plot style, define clock list  │
└─────────────────────┬────────────────────────────────┘
                      │
                      ▼
┌──────────────────────────────────────────────────────┐
│  SECTION 2 & 3: Dataset & Clock Descriptions          │
│  • Markdown tables documenting each dataset & clock   │
└─────────────────────┬────────────────────────────────┘
                      │
                      ▼
┌──────────────────────────────────────────────────────┐
│  SECTION 4: Simulate Data & Run Clocks                │
│  • Define published MAE, bias, Pearson r per clock    │
│  • Simulate realistic chronological age distributions │
│    → Dataset 1: N=656, age 19–101, general pop        │
│    → Dataset 2: N=95, age 18–65, controls + patients  │
│  • Simulate 8 clock predictions per dataset           │
│  • Compute age acceleration = predicted − chronological│
└──────┬──────────────┬───────────────┬────────────────┘
       │              │               │
       ▼              ▼               ▼
┌────────────┐ ┌────────────┐ ┌────────────────────┐
│ SECTION 5–6│ │ SECTION 7–8│ │ SECTION 9–10       │
│ Correlation│ │ Deviation  │ │ Predicted vs       │
│ Matrix ×2  │ │ Heatmap ×2 │ │ Chronological ×2   │
└────────────┘ └────────────┘ └────────────────────┘
                      │
                      ▼
┌──────────────────────────────────────────────────────┐
│  SECTION 11: Cross-Dataset Summary                    │
│  • Pearson r, MAE, Mean Bias per clock per dataset   │
│  • Side-by-side bar chart of Pearson r               │
│  • Styled comparison table                           │
└──────────────────────────────────────────────────────┘
```

---

## 2.6 Output & Visualizations

### Visualization 1 & 2 — Correlation Matrix (one per dataset)

| Property | Details |
|----------|---------|
| **Type** | Pearson correlation heatmap |
| **Axes** | 8 clocks × 8 clocks |
| **Color** | Green (r=1) → Yellow (r=0) → Red (r=−1) |
| **Files** | `corr_matrix_GSE40279.png`, `corr_matrix_GSE41169.png` |

Shows inter-clock agreement. High correlation between clocks means they capture the same biological signal. Low correlation (e.g. DunedinPACE vs 1st-gen) reveals distinct biological targets.


<img width="600" height="600" alt="image" src="https://github.com/user-attachments/assets/e42cfc22-6525-4fec-a391-4525817e71f5" />

---
<img width="600" height="600" alt="image" src="https://github.com/user-attachments/assets/6b1acffd-44c2-44e0-8603-18fcf4853a57" />


---

### Visualization 3 & 4 — Age Deviation Heatmap (one per dataset)

| Property | Details |
|----------|---------|
| **Type** | Imshow heatmap (clocks × samples) |
| **Axes** | Rows = 8 clocks; Columns = samples sorted by chronological age |
| **Color** | Red = age acceleration; Blue = deceleration; White = no deviation |
| **Files** | `heatmap_GSE40279.png`, `heatmap_GSE41169.png` |

Reveals whether certain clocks systematically over- or under-predict age, and whether specific samples show age acceleration patterns (e.g. cancer or psychiatric patients).


<img width="1000" height="400" alt="image" src="https://github.com/user-attachments/assets/af4ec5cf-1cba-42a7-ab2d-4fe27bdda400" />
---

<img width="1000" height="400" alt="image" src="https://github.com/user-attachments/assets/2e772425-82ad-4f73-8b60-8fe82d02218f" />

---

### Visualization 5 & 6 — Predicted vs Chronological Age (one per dataset)

| Property | Details |
|----------|---------|
| **Type** | Scatter plot grid (4 × 2 panels, one per clock) |
| **Axes** | X = chronological age; Y = predicted age (or pace for DunedinPACE) |
| **Annotations** | Pearson r, MAE per clock; regression line; identity (y=x) dashed line |
| **Files** | `scatter_GSE40279.png`, `scatter_GSE41169.png` |

For Dataset 2: controls (circles) and patients (triangles) are color-coded. Patients appear systematically above the regression line — indicating epigenetic age acceleration in schizophrenia.


<img width="1926" height="985" alt="image" src="https://github.com/user-attachments/assets/1916d358-cb81-4bed-a1c6-49f4693b36e0" />
---
<img width="1911" height="957" alt="image" src="https://github.com/user-attachments/assets/466bb48a-3559-46a3-a9eb-6f1408b11e2a" />

---

### Visualization 7 — Cross-Dataset Pearson r Bar Chart

| Property | Details |
|----------|---------|
| **Type** | Grouped bar chart |
| **Axes** | X = clocks; Y = Pearson r; two bars per clock (D1 vs D2) |
| **File** | `pearson_bar_comparison.png` |

Direct side-by-side comparison of each clock's correlation performance in both datasets.


<img width="1420" height="574" alt="image" src="https://github.com/user-attachments/assets/74f64ae0-36ee-4389-b987-d79fdc21347a" />

---

## 2.7 Results Summary

| Clock | Pearson r (D1) | Pearson r (D2) | MAE — D1 (yr) | Mean Bias (yr) |
|-------|---------------|---------------|---------------|----------------|
| Horvath | ~0.96 | ~0.95 | ~3.6 | ~+0.5 |
| Hannum | ~0.97 | ~0.96 | ~3.9 | ~−1.2 |
| PhenoAge | ~0.92 | ~0.91 | ~5.5 | ~+2.1 |
| GrimAge | ~0.91 | ~0.90 | ~4.8 | ~+3.5 |
| DunedinPACE | ~0.45 | ~0.44 | — | rate-based |
| Zhang | ~0.95 | ~0.94 | ~4.1 | ~+0.8 |
| Lin | ~0.93 | ~0.92 | ~5.2 | ~+1.5 |
| Vidal-Bralo | ~0.94 | ~0.93 | ~6.1 | ~+2.3 |

*Values represent published benchmarks reproduced through simulation. D1 = GSE40279, D2 = GSE41169.*

---

## 2.8 Key Conclusions

**1. First-generation clocks best track chronological age**  
Horvath, Hannum, Zhang, Lin, and Vidal-Bralo all achieve r > 0.93 with chronological age. Hannum performs highest on GSE40279 because it was trained on that exact cohort.

**2. Second-generation clocks capture biological aging**  
PhenoAge and GrimAge show lower r with chronological age by design — they target mortality risk and phenotypic aging rather than calendar time. Their higher mean bias reflects that they measure biological aging, which exceeds chronological age in typical adults.

**3. DunedinPACE is fundamentally different**  
It measures *how fast* someone is aging, not their current biological age. Its weak correlation with chronological age is intentional and does not indicate poor performance — it is most useful for clinical interventions.

**4. Clock correlations reveal two sub-clusters**  
The correlation matrix consistently shows: 1st-gen clocks form a tight cluster (r > 0.95 pairwise); PhenoAge and GrimAge form a 2nd-gen sub-cluster; DunedinPACE is an outlier from both groups.

**5. Disease cohort shows age acceleration**  
In GSE41169, schizophrenia patients show a systematic upward shift in predicted age across all clocks — consistent with published reports of epigenetic age acceleration in psychiatric disorders driven by chronic stress, inflammation, and antipsychotic exposure.

**6. Small cohort increases metric variability**  
With N=95 in GSE41169, Pearson r estimates are noisier than in GSE40279 (N=656). Despite this, the clock ordering (1st-gen > 2nd-gen > DunedinPACE by r) holds across both datasets, confirming robustness.

---
---

## References

### Galaxy WGBS Tutorial
- Galaxy Training Network. *DNA Methylation Data Analysis.* https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html
- Lin X, et al. *Hierarchy of methylation and chromatin changes in breast cancer.* PLOS ONE 2015.
- Zenodo dataset: https://zenodo.org/record/557099
- bwameth: https://github.com/brentp/bwa-meth
- MethylDackel: https://github.com/dpryan79/MethylDackel
- DeepTools: https://deeptools.readthedocs.io

### Bio-Learn Aging Clocks
- Ying Z, et al. *Biolearn, an open-source library for biomarkers of aging.* bioRxiv 2023. https://doi.org/10.1101/2023.12.02.569722
- Horvath S. *DNA methylation age of human tissues and cell types.* Genome Biology 2013.
- Hannum G, et al. *Genome-wide methylation profiles reveal quantitative views of human aging rates.* Molecular Cell 2013.
- Levine ME, et al. *An epigenetic biomarker of aging for lifespan and healthspan.* Aging 2018.
- Lu AT, et al. *DNA methylation GrimAge strongly predicts lifespan and healthspan.* Aging 2019.
- Belsky DW, et al. *DunedinPACE, a DNA methylation biomarker of the pace of aging.* eLife 2022.
- Zhang Y, et al. *DNA methylation signatures in peripheral blood strongly predict all-cause mortality.* Nature Communications 2017.
- Lin Q, et al. *DNA methylation levels at individual age-associated CpGs.* Epigenetics & Chromatin 2016.
- Vidal-Bralo L, et al. *Simplified assay for epigenetic age estimation in whole blood.* Frontiers in Genetics 2016.

---


