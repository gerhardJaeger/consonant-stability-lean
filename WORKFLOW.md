# Workflow

This document describes the order in which scripts are to be executed to reproduce the results of the paper.

---

## 0. Environment setup

Install the R dependencies (cmdstanr, tidyverse, ape, posterior, loo, ggdist, geiger, treeio, phangorn):

```r
install.packages(c("tidyverse", "ape", "posterior", "loo", "ggdist", "geiger", "cmdstanr", "treeio", "phangorn"))
```

The file `consonant_inventories.yml` specifies a mamba environment with all Python dependencies:

```bash
mamba env create -f consonant_inventories.yml
mamba activate consonant_inventories
```

---

## 1. Data preparation

**Input:** Raw PHOIBLE data (`phoible/data/phoible.csv`), included in this repository. Retrieved from https://github.com/phoible/dev (commit `7030ae02`, 2023-04-16) on 2024-04-19. Glottolog language metadata (`data/languages.csv`) is also included.

```bash
python scripts/process-inventories.py -l data/languages.csv
```
Selects one inventory per language and computes markedness-weighted segment counts. Writes `processed/inventories.csv` and `processed/segments.csv`.

```bash
python scripts/generate-series.py \
  -i processed/inventories.csv \
  -s processed/segments.csv \
  -t data/transformation_rules.tsv \
  -o data/
```
Generates phonological charts for each language and computes obstruent series counts. Writes `data/charts_full/` and `data/series_counts.tsv`.

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
Fits the Bayesian Ornstein–Uhlenbeck model (`stan/OU.stan` and `stan/OU_logtransformed.stan`) to the global tree for all four variables (plosive, fricative, affricate, series). Selects the preferred transformation per variable via LOOIC. Also:
- Exports per-node posterior ASR samples to `data/asr/asr_{variable}.csv`
- Generates `figures/halflive.svg` (posterior distributions of t½)
- Generates `figures/landscapes.svg` (inferred equilibrium distributions)

---

## 4. Ancestral state reconstruction (world and family trees)

```r
Rscript R/familywise_asr.r
```
Reads the ASR sample CSVs from `data/asr/`, computes per-node posterior summaries (mean, 2.5th–97.5th percentile), and writes annotated NEXUS files for the world tree and for each of the 8 family subtrees to `data/asr/`.

---

## 5. Sensitivity analysis

```bash
bash R/run_sensitivity_analysis.sh
```
Calls `R/OU_sensitivity_analysis.r` once per variable and model specification (12 specifications × 4 variables = 48 runs). Each run uses the Stan template `stan/OU_template.stan` with priors filled in from the model bundle. Results are written to `sensitivity_analysis/{variable}/{model}/`.

```r
Rscript R/evaluate_sensitivity_analysis.r
```
Reads the per-model outputs and produces the summary table `sensitivity_analysis_results.csv`.

---

## 6. Simulation study

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
Reads the fitted summaries and produces the parameter-recovery figure `figures/simulation_recovery.svg`.

---

## 7. Figure generation

```r
Rscript R/code_for_figures.r
```
Reads `data/data_pruned.csv` and produces the distribution and geographic map figures. Writes PDFs and SVGs to `figures/`.

```r
Rscript R/plot_sensitivity_analysis.r
```
Reads `sensitivity_analysis_results.csv` and writes `figures/sensitivity_analysis.pdf/.svg`.

```r
Rscript R/plot_simulations.r
```
Reads `data/simulations_fitted/` (OU datasets) and writes `figures/simulation_recovery.pdf/.svg`.

```r
Rscript R/plot_equilibrium_recovery.r
```
Reads `data/simulations_fitted/` (OU and BBM datasets) and writes `figures/equilibrium_recovery.pdf/.svg`.
