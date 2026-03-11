suppressMessages({
    library(tidyverse)
    library(ape)
    library(cmdstanr)
    library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript OU_stan_simulations.r <data_csv>")
}
data_csv <- args[1]

message("Using data CSV: ", data_csv)

d <- suppressMessages(read_csv(data_csv))
colnames(d) <- c("glottocode", "series_markedness_fullness")


tree <- read.nexus("data/tree_pruned.nex")
tree <- reorder(tree, "postorder")

# Align to tree order
d <- d %>% filter(glottocode %in% tree$tip.label)
d <- d[match(tree$tip.label, d$glottocode), ]


stan_data <- list(
  Ndata        = nrow(d),
  Nnodes       = length(unique(as.vector(tree$edge))),
  Nedges       = nrow(tree$edge),
  edges        = tree$edge,
  edge_lengths = tree$edge.length,
  root_node    = length(tree$tip.label) + 1,
  y            = d$series_markedness_fullness
)
storage.mode(stan_data$edges) <- "integer"

# Floor zero-length edges (also guarded in Stan)
eps <- 1e-9
stan_data$edge_lengths[stan_data$edge_lengths == 0] <- eps



mod <- cmdstan_model("stan/OU.stan")
fit <- mod$sample(
  data = stan_data,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  chains          = 2,
  parallel_chains = 2,
  save_warmup     = FALSE,
  seed            = 42
)



# Summaries (omit the big arrays; here only globals of interest)
fit$summary(variables = c("mu", "sigma", "t_half", "tau"))


summary <- fit$summary(variables = c("mu", "sigma", "t_half", "tau")) %>% as_tibble


outdir <- file.path("data", "simulations_fitted")

if (file.exists(outdir) && !dir.exists(outdir)) {
    stop(sprintf("Path '%s' exists and is not a directory", outdir))
}

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message("Created directory: ", outdir)
} else {
    message("Directory already exists: ", outdir)
}

write_csv(summary,
          file.path(outdir,
                    paste0(tools::file_path_sans_ext(basename(data_csv)),
                           ".summary.csv")))