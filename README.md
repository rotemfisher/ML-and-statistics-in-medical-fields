# 🧬 From Gene Expression to Clinical Outcomes
### A Comparative Analysis of Statistical and Machine Learning Models in Pan-Cancer Survival Prediction

> **Course Project** — Machine Learning & Statistics in Medical Applications  
> **Authors:** Rotem Fisher & Orian Aziz  
> **Dataset:** TCGA Pan-Cancer Atlas  

---

## 📋 Table of Contents

- [Overview](#overview)
- [Research Questions](#research-questions)
- [Results Summary](#results-summary)
- [Repository Structure](#repository-structure)
- [Dataset](#dataset)
- [Pipeline](#pipeline)
- [Installation & Setup](#installation--setup)
- [Reproducing the Analysis](#reproducing-the-analysis)
- [Key Findings](#key-findings)
- [Limitations](#limitations)

---

## Overview

Traditional cancer prognosis relies primarily on clinical staging systems; however, patients with identical stages often exhibit drastically different survival outcomes. This project investigates whether **transcriptome-wide RNA-Seq gene expression profiles** can improve survival prediction beyond conventional clinical variables in a pan-cancer cohort.

We design a structured machine learning pipeline for survival analysis, comparing:
- **Cox Proportional Hazards (CoxPH)** — classical statistical baseline with L2 regularization
- **Random Survival Forest (RSF)** — non-linear, tree-based survival model with Optuna hyperparameter tuning and permutation-based feature selection

The full analysis is grounded in survival analysis theory (concordance index, Kaplan–Meier curves, Log-rank tests), dimensionality reduction (PCA), and rigorous model evaluation.

---

## Research Questions

| # | Question |
|---|----------|
| **Primary** | Can genomic features improve survival prediction beyond traditional clinical variables in Pan-Cancer patients? |
| **Secondary** | Does a non-linear ML model (RSF) outperform a classical statistical model (Cox PH)? |
| **Hypothesis** | High-dimensional genomic expression patterns encode prognostic information not captured by clinical variables alone. |

---

## Results Summary

| Model | C-Index |
|---|---|
| CoxPH — Clinical only | 0.6667 |
| CoxPH — Clinical + Genomic (PCA) | 0.6633 |
| Naïve RSF (default hyperparameters) | 0.7246 |
| Optimized RSF (Optuna tuning) | 0.7427 |
| **Refined RSF — Top-100 Features** | **0.7827** |

**Δ Improvement (Clinical CoxPH → Refined RSF): +0.116**

Risk stratification using refined RSF scores produced a pronounced Kaplan–Meier separation between high- and low-risk patient groups, confirmed by a highly significant Log-rank test (χ² = 400.88, **p = 3.55 × 10⁻⁸⁹**).

---

## Repository Structure

```
ML-and-statistics-in-medical-fields/
│
├── data/                                  # Data directory (see Dataset section)
│   ├── EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena   # RNA-Seq gene expression
│   ├── Survival_SupplementalTable_S1_20171025_xena_sp          # Clinical & survival metadata
│   └── processed_pancan_pca.csv                                # Preprocessed output (generated)
│
├── PreProcess.ipynb                       # Data loading, merging, PCA, feature engineering
├── EDAGen.ipynb                           # Exploratory analysis of gene expression
├── EDAMetaData.ipynb                      # Exploratory analysis of clinical metadata
├── CoxHazard.ipynb                        # Cox PH baseline models (clinical-only & genomic)
├── RSF.ipynb                              # Random Survival Forest + Optuna tuning + evaluation
│
├── permutation_importance.csv             # Feature importance scores from optimized RSF
│
├── From_Gene_Expression_to_Clinical_Outcomes.pdf    # Full written report (paper format)
├── _Final_ML_and_statistics_in_medical_fields_slides.pdf   # Presentation slides
│
└── README.md
```

---

## Dataset

This project uses two publicly available files from the **TCGA Pan-Cancer Atlas** via the [UCSC Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)):

| File | Description | Source |
|------|-------------|--------|
| `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena` | Batch-corrected whole-transcriptome RNA-Seq gene expression (~20,530 genes × ~11,000 patients) | UCSC Xena |
| `Survival_SupplementalTable_S1_20171025_xena_sp` | Patient-level clinical annotations and overall survival endpoints | UCSC Xena |

### Downloading the Data

1. Visit [UCSC Xena PANCAN Hub](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN))
2. Download the two files above and place them in the `data/` directory
3. Run `PreProcess.ipynb` to generate `processed_pancan_pca.csv`

> **Note:** The raw data files are large (~1 GB for gene expression). The processed file `processed_pancan_pca.csv` is generated locally and is not tracked in this repository.

---

## Pipeline

```
Raw Data (RNA-Seq + Clinical)
         │
         ▼
┌─────────────────────┐
│   PreProcess.ipynb  │  Merge datasets, filter, impute, z-score normalize,
│                     │  PCA (95% variance) → 3,390 components, label encode
└────────┬────────────┘
         │  processed_pancan_pca.csv (10,952 × 3,399)
         ▼
┌─────────────────────┐
│  EDAGen.ipynb       │  Distribution of gene expression, variance structure,
│  EDAMetaData.ipynb  │  survival time, censoring rate, cancer type breakdown
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  CoxHazard.ipynb    │  Baseline CoxPH (clinical-only) → C-index: 0.6667
│                     │  Extended CoxPH (clinical + PCA) → C-index: 0.6633
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  RSF.ipynb          │  Naïve RSF → C-index: 0.7246
│                     │  Optuna-tuned RSF → C-index: 0.7427
│                     │  Permutation importance → Top-100 features
│                     │  Refined RSF → C-index: 0.7827
│                     │  Kaplan–Meier + Log-rank stratification
└─────────────────────┘
```

---

## Installation & Setup

### Prerequisites

- Python 3.9+
- Jupyter Notebook or JupyterLab

### Install Dependencies

```bash
pip install pandas numpy scikit-learn scikit-survival matplotlib seaborn optuna lifelines
```

Or using a requirements file (if provided):

```bash
pip install -r requirements.txt
```

### Clone the Repository

```bash
git clone https://github.com/rotemfisher/ML-and-statistics-in-medical-fields.git
cd ML-and-statistics-in-medical-fields
```

---

## Reproducing the Analysis

Run the notebooks **in the following order**:

### Step 1 — Preprocessing
```bash
jupyter notebook PreProcess.ipynb
```
Outputs: `data/processed_pancan_pca.csv`

### Step 2 — Exploratory Data Analysis
```bash
jupyter notebook EDAGen.ipynb
jupyter notebook EDAMetaData.ipynb
```

### Step 3 — Baseline Cox Models
```bash
jupyter notebook CoxHazard.ipynb
```
Outputs: C-index scores for clinical-only and genomic CoxPH models.

### Step 4 — Random Survival Forest
```bash
jupyter notebook RSF.ipynb
```
Outputs: C-index scores, permutation importance plot, Kaplan–Meier curves, Log-rank test results.

> ⚠️ **Runtime note:** Optuna hyperparameter tuning in `RSF.ipynb` may take 30–90 minutes depending on hardware. The number of trials can be reduced in the notebook by modifying `n_trials`.

---

## Key Findings

1. **Clinical variables alone** provide moderate predictive performance (C-index ≈ 0.667), consistent with known limitations of stage-based prognosis.

2. **Linear genomic modeling (CoxPH + PCA)** does not improve over the clinical baseline — transcriptomic information is not linearly separable under proportional hazards assumptions.

3. **Non-linear survival modeling (RSF)** substantially improves performance (+0.076 C-index over CoxPH), demonstrating that molecular-clinical interactions require non-parametric modeling.

4. **Feature refinement** (Top-100 permutation-importance features) further improves generalization (+0.04 C-index), confirming that survival signals are concentrated in a compact feature subset.

5. **Risk stratification** using refined RSF scores produces biologically meaningful patient separation, with a Log-rank p-value of 3.55 × 10⁻⁸⁹.

---

## Limitations

- **Pan-cancer heterogeneity:** The model does not explicitly account for cancer-type–specific baseline hazard functions.
- **PCA interpretability:** Transcriptomic features are represented as principal components, limiting gene-level biological interpretation.
- **Zero imputation:** Missing RNA-Seq values were imputed with zero, which may not fully represent measurement noise.
- **Single dataset evaluation:** No external validation cohort was used; performance estimates may reflect dataset-specific characteristics.

---

## Citation

If you use this work, please cite:

```
Fisher, R. & Aziz, O. (2025). From Gene Expression to Clinical Outcomes: 
A Comparative Analysis of Statistical and Machine Learning Models in 
Pan-Cancer Survival Prediction. Course Project, ML & Statistics in Medical Applications.
```

---

<div align="center">
  <sub>Built with ❤️ using Python, scikit-survival, and the TCGA Pan-Cancer Atlas</sub>
</div>
