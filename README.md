# Population Structure Analysis — 1000 Genomes Project (Phase 3)

![PLINK](https://img.shields.io/badge/PLINK-1.9-blue)
![ADMIXTURE](https://img.shields.io/badge/ADMIXTURE-1.3-green)
![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r)
![License](https://img.shields.io/badge/License-MIT-yellow)

A reproducible population genomics pipeline analysing human genetic diversity across 26 populations using **1000 Genomes Phase 3** data. Covers quality control, principal component analysis (PCA), ancestry inference (ADMIXTURE), and runs-of-homozygosity (ROH) — the same methods used in ancient DNA and modern population genetics research.

---

## 🎯 Objectives

1. Characterise global human population structure using PCA
2. Infer ancestry proportions across superpopulations using ADMIXTURE (K = 2–7)
3. Quantify inbreeding and recent relatedness signals using ROH analysis
4. Produce publication-quality visualisations in R

---

## 🗂 Repository Structure

```
population-genomics-1000genomes/
├── README.md
├── data/               # Raw data not tracked (see download script)
│   └── .gitkeep
├── scripts/
│   ├── 01_download_data.sh       # Download 1000G Phase 3 data
│   ├── 02_qc_filtering.sh        # QC and variant filtering with PLINK
│   ├── 03_pca_analysis.R         # PCA computation and visualisation
│   ├── 04_admixture.sh           # ADMIXTURE ancestry inference
│   └── 05_roh_analysis.R         # ROH analysis and visualisation
├── results/
│   └── figures/        # Output plots (PCA, ADMIXTURE bar chart, ROH)
└── docs/
    └── methods.md      # Detailed methods description
```

---

## 🔬 Methods Overview

### 1. Data
- **Source**: 1000 Genomes Project Phase 3 (2,504 individuals, 26 populations, 5 superpopulations: AFR, AMR, EAS, EUR, SAS)
- **Format**: VCF → converted to PLINK binary format (BED/BIM/FAM)

### 2. Quality Control (Script 02)
| Filter | Threshold |
|---|---|
| Minor allele frequency (MAF) | > 0.05 |
| Genotype missingness per SNP | < 0.05 |
| Genotype missingness per individual | < 0.10 |
| Hardy-Weinberg equilibrium | p > 1×10⁻⁶ |
| LD pruning (window/step/r²) | 50 / 10 / 0.2 |

### 3. PCA (Script 03)
- Computed using PLINK `--pca` on LD-pruned SNP set
- Visualised in R using `ggplot2` — coloured by superpopulation and population

### 4. ADMIXTURE (Script 04)
- Run at K = 2 to 7 with cross-validation to identify optimal K
- Results plotted as stacked bar charts per individual, ordered by population

### 5. ROH Analysis (Script 05)
- Called using PLINK `--homozyg` with minimum ROH length 500 kb
- Distributions of total ROH length (FROH) plotted by superpopulation
- Connects to IBD and kinship inference methods used in ancient DNA research

---

## ⚙️ Requirements

### Software
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [ADMIXTURE 1.3](https://dalexander.github.io/admixture/)
- R ≥ 4.0 with packages: `ggplot2`, `dplyr`, `tidyr`, `RColorBrewer`, `cowplot`

### Install R packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "cowplot"))
```

---

## 🚀 How to Reproduce

```bash
# Step 1 — Download data
bash scripts/01_download_data.sh

# Step 2 — QC and filtering
bash scripts/02_qc_filtering.sh

# Step 3 — PCA
Rscript scripts/03_pca_analysis.R

# Step 4 — ADMIXTURE
bash scripts/04_admixture.sh

# Step 5 — ROH analysis
Rscript scripts/05_roh_analysis.R
```

---

## 📊 Results

> Figures will be added upon pipeline completion.

| Analysis | Output |
|---|---|
| PCA | PC1 vs PC2 coloured by superpopulation |
| ADMIXTURE | Stacked bar chart at optimal K |
| ROH | FROH distribution by superpopulation |

---

## 🔗 Connection to Ancient DNA Research

The ROH and IBD methods implemented here are foundational to **ancient DNA pedigree reconstruction** — the same statistical framework underlies tools like `hapROH`, `ancIBD`, and `KIN` used in archaeogenomics to infer kinship and social structure in prehistoric populations.

---

## 👤 Author

**Niharika Joshi**  
MSc Biological Sciences, Seoul National University  
[LinkedIn](https://linkedin.com/in/niharika-joshi-003a48194) · [Email](mailto:joshiniharika243@gmail.com)

---

## 📜 License
MIT License — free to use and adapt with attribution.