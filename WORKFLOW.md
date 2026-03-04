# Workflow

This document describes the order in which scripts are to be executed to reproduce the results of the paper.

---

## 0. Environment setup

Install the Python dependencies:

```bash
pip install -r scripts/requirements.txt
```

Install the R dependencies (cmdstanr, tidyverse, ape, posterior, loo, BBMV, ggdist, geiger):

```r
install.packages(c("tidyverse", "ape", "posterior", "loo", "ggdist", "geiger", "cmdstanr"))
# BBMV is installed from GitHub:
devtools::install_github("fcboucher/BBMV")
```

The file `consonant_inventories.yml` specifies a conda environment with all Python dependencies.

---

## 1. Data preparation

**Input:** Raw PHOIBLE data (download separately from https://phoible.org/; place `phoible.csv` in `phoible/data/`).

```bash
python scripts/process-inventories.py
```
Selects one inventory per language and computes markedness-weighted segment counts. Writes `processed/inventories.csv` and `processed/segments.csv`.

```bash
python scripts/generate-series.py
```
Generates phonological charts for each language and computes obstruent series counts. Writes `charts/charts_full/` and `data/series_counts.tsv`.

---

## 2. Tree preparation

```r
Rscript R/replace_subtrees.r
```
Takes the global language tree (`data/tree_pruned.nex`) and replaces its family-level clades with dated family-specific trees from `data/families/`. Writes `data/tree_families_replaced.nex`.

---

## 3. Global OU model fitting

```r
Rscript R/OU_stan.r
```
Fits the Bayesian Ornstein–Uhlenbeck model (`stan/OU.stan` and `stan/OU_logtransformed.stan`) to the global tree for all four variables (plosive, fricative, affricate, series). Writes fitted model objects to `data/family_models/`.

Model comparison between the log-transformed and identity-transformed versions is performed here; the preferred transformation for each variable is selected based on LOOIC.

---

## 4. Sensitivity analysis

```bash
bash R/run_sensitivity_analysis.sh
```
Calls `R/OU_sensitivity_analysis.r` once per variable and model specification (12 specifications × 4 variables = 48 runs). Each run uses the Stan template `stan/OU_template.stan` with priors filled in from the model bundle. Results are written to `sensitivity_analysis/{variable}/{model}/`.

```r
Rscript R/evaluate_sensitivity_analysis.r
```
Reads the per-model outputs and produces the summary table `sensitivity_analysis/sensitivity_analysis_results.csv`.

---

## 5. Simulation study

```r
Rscript R/simulated_data.R
```
Simulates trait evolution on the global tree under Brownian Motion (10 datasets) and the OU model (120 datasets across a grid of σ and t½ values). Writes simulated datasets to `data/simulated/`.

```bash
bash R/fit_simulated_data.sh
```
Calls `R/OU_stan_simulations.r` for each simulated dataset, fitting the Bayesian OU model and writing summary CSVs to `data/simulations_fitted/`.

```r
Rscript R/evaluate_simulations.r
```
Reads the fitted summaries and produces parameter-recovery statistics.

---

## 6. BBMV landscape analysis

```r
Rscript R/BBMV_consonants_forest.r
```
Fits the Fokker–Planck–Kolmogorov (FPK/BBMV) model to the series variable on the global tree to estimate the fitness landscape (equilibrium distribution).

---

## 7. Family-level reconstructions

Run each family script independently. Each reads the family-specific tree from `data/families/<family>/phylogeny/` and fits the OU model using `R/bayesian_fit.r`.

```r
Rscript R/families/Indo-European_OU.r
Rscript R/families/Austronesian_OU.r
Rscript R/families/Pama-Nyungan_OU.r
Rscript R/families/Sino-Tibetan_OU.r
Rscript R/families/Tupi-Guarani_OU.r
Rscript R/families/Turkic_OU.r
Rscript R/families/Uralic_OU.r
```

---

## 8. Grambank comparison

```r
Rscript R/grambank_rates.r
```
Fits discrete-character ARD models to binarized Grambank features and to the binarized obstruent series variables. Produces characteristic-time estimates for comparison across typological features. Requires Grambank data (download separately from https://grambank.clld.org/).

---

## 9. Figure generation

```r
Rscript R/code_for_figures.r
```
Reads fitted model outputs and produces all publication figures. Writes PDFs and SVGs to `figures/`.
