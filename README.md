# Microbiome GAMLSS Analysis
This repository contains an R-based workflow for performing univariable beta zero-inflated regression (GAMLSS-BEZI) on gut microbiome composition data at different taxonomic levels: Phylum, Class, Family, and Genus.

---

## 📁 Repository Structure

```
microbiome-gamlss/
├── README.md
├── .gitignore
├── scripts/
│   └── microbiome_gamlss_analysis.R
├── data/
│   ├── phyla&covariates.RData
│   ├── class&covariates.RData
│   ├── family&covariates.RData
│   └── genus&covariates.RData
```

---

## 📦 Data Requirements

Each `.RData` file must contain a `data.frame` with **taxa proportions (0–1)** and **covariates** as columns:

* Example for `class&covariates.RData`: must contain an object named `class` with columns like:

```
$ Clostridia %          : num  0.822 0.638 ...
$ Bacteroidia %         : num  0.0923 0.2799 ...
...
$ Female gender         : num 0 1 1 ...
$ Age, BMI, Antibiotics, ...
```

### Required covariates:

```
"Female gender", "Age", "BMI", "Antibiotics",
"Anti-inflammatories", "Proton pump inhibitors", "Probiotics/Prebiotics",
"Disorders", "Hours of sleep", "Tobacco",
"Pathologies", "Physical activity", "Allergies"
```

---

## ▶️ Running the Script

### Prerequisites

Install required R packages if not already installed:

```r
install.packages(c("gamlss", "broom", "dplyr", "tidyr", "ggplot2", "gtsummary"))
```

### Run

```r
# Set working directory to the repository root
setwd("path/to/microbiome-gamlss")

# Source the analysis script
source("scripts/microbiome_gamlss_analysis.R")
```

---

## 📊 Output

The script produces heatmaps with coefficient estimates for each taxon-covariate pair, highlighting significant results (FDR-adjusted p < 0.05).

---


## 👩‍🔬 Author

Maria De Martino
PhD student in Biostatistics

---

For issues, suggestions, or improvements, feel free to open a GitHub issue or pull reque

