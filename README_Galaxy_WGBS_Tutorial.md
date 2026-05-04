# README: Galaxy WGBS Methylation-Seq Tutorial

> **Tutorial:** [Methylation-Seq Analysis — Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html)  
> **Background Slides:** [Introduction to DNA Methylation](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/introduction-dna-methylation/slides-plain.html)

---

## Overview

This hands-on Galaxy tutorial walks through a complete **WGBS (Whole Genome Bisulfite Sequencing)** analysis pipeline using the Galaxy web platform — no command-line required. The workflow reproduces the breast cancer methylome analysis from Bock et al. (PLOS ONE, 2015), comparing methylation profiles between normal and cancerous breast tissue.

**Galaxy** is an open, web-based bioinformatics platform that makes tools accessible without coding, while ensuring reproducibility through automatic metadata capture.

---

## Background: DNA Methylation Basics

| Concept | Explanation |
|---------|-------------|
| **5-methylcytosine (5mC)** | Primary form of DNA methylation in mammals; occurs mainly at CpG dinucleotides |
| **Bisulfite conversion** | Converts unmethylated cytosine → uracil (→ thymine); methylated cytosine remains unchanged |
| **CpG islands** | Regions with high CpG density; often at gene promoters |
| **Hypermethylation** | Silences gene expression — common at tumor suppressor genes in cancer |
| **Hypomethylation** | Can activate oncogenes — also common in cancer |

---

## Input Data

| File | Description |
|------|-------------|
| `subset_1.fastq` | Forward reads — bisulfite-converted WGBS |
| `subset_2.fastq` | Reverse reads — bisulfite-converted WGBS |
| **Source** | Zenodo record 557099 |
| **Organism** | *Homo sapiens* |
| **Reference** | hg38 (chr7 subset for tutorial efficiency) |

```
# Galaxy Data Upload — paste these URLs into Galaxy Upload tool:
https://zenodo.org/record/557099/files/subset_1.fastq
https://zenodo.org/record/557099/files/subset_2.fastq
```

---

## Tools Used (Galaxy Tool IDs)

| Step | Tool | Purpose |
|------|------|---------|
| 1 | **FastQC** | Raw read quality assessment |
| 2 | **Trim Galore** | Adapter removal + quality trimming (bisulfite-aware) |
| 3 | **FastQC** (post-trim) | Verify trimming quality |
| 4 | **Bismark** (genome prep + align) | Align bisulfite reads to reference; handles C→T and G→A converted reads |
| 5 | **Samtools flagstat** | Alignment statistics |
| 6 | **Bismark deduplicate** | Remove PCR duplicates |
| 7 | **Bismark methylation extractor** | Call methylation at each CpG position |
| 8 | **MethylDackel** | Additional methylation bias assessment |
| 9 | **bismark2bedGraph** | Convert methylation calls to bedGraph/bigWig |
| 10 | **DeepTools bamCompare** | Methylation visualization tracks |

---

## Step-by-Step Workflow

### Step 1 — Quality Control
- Upload FASTQ files to Galaxy history
- Run **FastQC** on both reads
- Check: sequence quality, GC content, adapter contamination, per-base quality

### Step 2 — Trimming
- Run **Trim Galore** (paired-end mode, bisulfite-specific settings)
- Removes Illumina adapters and low-quality bases (Q < 20)
- Re-run FastQC to confirm improvement

### Step 3 — Alignment
- Run **Bismark** with hg38 reference (chr7 for tutorial)
- Bismark builds a C→T and G→A converted reference for bisulfite-aware alignment
- Output: `.bam` file with methylation tags

### Step 4 — Deduplication & Sorting
- **Bismark deduplicate** removes PCR clonal reads
- **Samtools sort + index** for downstream use

### Step 5 — Methylation Extraction
- **Bismark methylation extractor** produces CpG-context files
- Output: position, strand, count methylated, count unmethylated
- Generates M-bias plot to check for end-of-read artifacts

### Step 6 — Visualization
- Convert to **bedGraph / bigWig** format
- Load in **UCSC Genome Browser** or Galaxy's built-in display
- Compare methylation tracks between samples

---

## Key Output Files

| File | Format | Content |
|------|--------|---------|
| `*.bam` | BAM | Aligned reads |
| `CpG_OT_*.txt` | Text | Original top-strand methylation calls |
| `*.bedGraph` | bedGraph | Per-CpG methylation percentage |
| `*.bw` | bigWig | Genome browser-ready methylation track |
| M-bias plot | PNG | Read-position methylation bias (QC) |
| Summary report | HTML | Bismark alignment + methylation statistics |

---

## Visualization Examples

- **M-bias plot**: Checks if methylation level is uniform across read position; if edges are biased, clip those bases
- **Genome browser track**: Per-CpG methylation % displayed as continuous signal over chr7
- **Coverage histogram**: Distribution of read depth across CpG sites

---

## Expected Results

| Metric | Expected Value |
|--------|----------------|
| Bisulfite conversion rate | > 98% |
| Mapping efficiency | 50–80% (lower than WGS due to bisulfite complexity) |
| Unique reads after deduplication | Depends on library complexity |
| CpG coverage (at 1×) | Variable; tutorial uses small subset |

---

## Galaxy-Specific Tips

- Use **Galaxy Histories** to track every step — automatically logged for reproducibility
- Use **Shared Workflows** from GTN to import the complete workflow in one click
- Intermediate files can be hidden to keep histories clean
- Use **Collections** for batch processing multiple samples

---

## Related Resources

| Resource | Link |
|----------|------|
| Galaxy Training Network | https://training.galaxyproject.org |
| Bismark documentation | https://www.bioinformatics.babraham.ac.uk/projects/bismark/ |
| WGBS review | PMC8963483 |
| Zenodo dataset | https://zenodo.org/record/557099 |

---
