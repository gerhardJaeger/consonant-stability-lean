# Consonant Inventory Stability: Replication Repository

This repository contains the code and data to reproduce the analyses and figures for the paper *[title omitted for review]*.

## Overview

We model the diachronic stability of consonant inventories across language families using Bayesian Ornstein–Uhlenbeck (OU) processes fitted to language phylogenies via Stan. The primary output is the **half-life** t½ = ln(2)/λ, which quantifies how quickly phonological traits revert toward an equilibrium. Four traits are analyzed: plosive series count, fricative series count, affricate series count, and combined obstruent series count.

## Requirements

A [mamba](https://mamba.readthedocs.io/) environment with all Python and R dependencies is specified in `consonant_inventories.yml`:

```bash
mamba env create -f consonant_inventories.yml
mamba activate consonant_inventories
```

Key R packages: `tidyverse`, `ape`, `cmdstanr`, `posterior`, `loo`, `ggdist`, `treeio`, `phangorn`, `patchwork`, `OUwie`.

CmdStan must be installed separately; see the [cmdstanr documentation](https://mc-stan.org/cmdstanr/).

## Reproduction

See [`WORKFLOW.md`](WORKFLOW.md) for the full step-by-step reproduction guide. In brief:

| Step | Script | Output |
|------|--------|--------|
| 1. Data preparation | `scripts/process-inventories.py`, `scripts/generate-series.py` | `processed/`, `data/series_counts.tsv` |
| 2. Tree preparation | `R/replace_subtrees.r` | `data/tree_families_replaced.nex` |
| 3. Global OU model | `R/OU_stan.r` | `data/asr/`, `figures/halflive.svg`, `figures/landscapes.svg` |
| 4. Ancestral state reconstruction | `R/familywise_asr.r` | `data/asr/*.nexus` |
| 5. Sensitivity analysis | `R/run_sensitivity_analysis.sh`, `R/evaluate_sensitivity_analysis.r` | `sensitivity_analysis_results.csv` |
| 6. Simulation study | `R/simulated_data.R`, `R/fit_simulated_data.sh`, `R/evaluate_simulations.r` | `data/simulations_fitted/` |
| 7. Figures | `R/code_for_figures.r`, `R/plot_sensitivity_analysis.r`, `R/plot_simulations.r`, `R/plot_equilibrium_recovery.r` | `figures/` |

## Data

- **PHOIBLE** (`phoible/data/phoible.csv`): phonological inventory database, commit `7030ae02` (2023-04-16), retrieved from [https://github.com/phoible/dev](https://github.com/phoible/dev).
- **Glottolog** (`data/languages.csv`): language metadata from [https://glottolog.org](https://glottolog.org).
- **Family phylogenies** (`data/families/`): dated trees for 8 language families (Indo-European, Austronesian, Bantu, Pama-Nyungan, Sino-Tibetan, Tupi-Guarani, Turkic, Uralic).

## Stan models

- `stan/OU.stan` — identity-scale OU model
- `stan/OU_logtransformed.stan` — log-scale OU model
- `stan/OU_template.stan` — parameterized template used by the sensitivity analysis

## Repository structure

```
.
├── data/               # Input data and generated outputs
│   ├── families/       # Family-specific phylogenies
│   ├── asr/            # Ancestral state reconstruction outputs
│   ├── simulated/      # Simulated datasets
│   └── simulations_fitted/  # Fitted summaries for simulation study
├── figures/            # All figures (SVG and PDF)
├── phoible/            # Bundled PHOIBLE data
├── processed/          # Intermediate processed data
├── R/                  # R analysis scripts
├── scripts/            # Python data preparation scripts
├── sensitivity_analysis/  # Per-model sensitivity analysis outputs
├── stan/               # Stan model files
└── WORKFLOW.md         # Step-by-step reproduction guide
```
