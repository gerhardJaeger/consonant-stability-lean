library(BBMV)
library(geiger)
library(tidyverse)
library(coda)
library(parallel)
library(bayesplot)

# Read the data and tree
d <- read.csv("../data/data_pruned.csv")
tree <- read.nexus("../data/tree_pruned.nex")
d <- d[match(tree$tip.label, d$glottocode), ]

# Define a function to prepare data and run MCMC
run_analysis <- function(variable_name, seeds, tree, d, params_to_keep) {
  trait <- d[[variable_name]]
  names(trait) <- d$glottocode

  run_mcmc <- function(seed) {
    set.seed(seed)
    tryCatch({
      MH_MCMC_FPK(
        tree,
        trait = trait,
        bounds = range(trait),
        Nsteps = 100000,
        record_every = 10,
        Npts = 50,
        pars_init = c(0, -4, -4, 0, 1),
        prob_update = c(0.2, 0.25, 0.25, 0.25, 0.05),
        verbose = TRUE,
        plot = FALSE,
        save_to = paste0(variable_name, "_MCMC_FPK_", seed, ".Rdata"),
        save_every = 10,
        type_priors = c(rep("Normal", 4), "Uniform"),
        shape_priors = list(
          c(0, 10),
          c(0, 10),
          c(0, 10),
          c(0, 10),
          NA
        ),
        proposal_type = "Uniform",
        proposal_sensitivity = c(0.1, 0.1, 0.1, 0.1, 1),
        prior.only = FALSE
      )
    }, error = function(e) {
      message("Error in run_mcmc with seed ", seed, " for ", variable_name, ": ", e$message)
      NULL
    })
  }

  # Run MCMC for all seeds in parallel
  results <- mclapply(seeds, run_mcmc, mc.cores = 2)

  # Subset, filter, and process results
  results <- lapply(results, function(chain) {
    if (is.null(chain)) return(NULL)
    chain[, params_to_keep, drop = FALSE]
  })
  results <- Filter(Negate(is.null), results)

  # Discard the first half of each chain
  results <- lapply(results, function(chain) {
    chain[(nrow(chain) %/% 2 + 1):nrow(chain), ]
  })

  results
}

# List of variables to analyze
variables <- c(
  "series_markedness_fullness",
  "plosive_markedness_fullness",
  "fricative_markedness_fullness",
  "affricate_markedness_fullness"
)

# Parameters to retain
params_to_keep <- c("sigsq", "a", "b", "c", "root", "lnprior", "lnlik")

# Run analyses in parallel
start_time <- Sys.time()
results_all <- mclapply(
  variables,
  function(var) run_analysis(var, seeds = c(123, 456), tree = tree, d = d, params_to_keep = params_to_keep),
  mc.cores = length(variables)
)
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)

# Save results for each variable
names(results_all) <- variables
saveRDS(results_all, file = "all_mcmc_results.rds")


# Function to compute Gelman diagnostics
compute_gelman_diagnostics <- function(chains) {
  # Combine chains into an mcmc.list
  mcmc_list <- lapply(chains, function(chain) as.mcmc(chain))
  mcmc_list <- as.mcmc.list(mcmc_list)
  
  # Compute Gelman diagnostics
  gelman_results <- gelman.diag(mcmc_list, autoburnin = FALSE)
  return(gelman_results)
}

# Apply the function to each element of results_all
gelman_diagnostics <- lapply(results_all, function(variable_results) {
  if (length(variable_results) > 1) {  # Ensure there are multiple chains
    compute_gelman_diagnostics(variable_results)
  } else {
    warning("Not enough chains for Gelman diagnostics for this variable.")
    NULL
  }
})

# Inspect Gelman diagnostics
gelman_diagnostics



