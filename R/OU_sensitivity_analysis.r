library(tidyverse)
library(ape)
library(cmdstanr)
library(posterior)
library(ggdist)
library(loo)
library(stringr)
library(bayesplot)

# library(future)
# library(furrr)

# plan(multisession)


traits <- c("plosive", "fricative", "affricate", "series")

args <- commandArgs(trailingOnly = TRUE)
trait_no <- if (length(args) > 0) as.integer(args[1]) else 1
trait <- traits[trait_no]

model_no <- if (length(args) > 1) as.integer(args[2]) else 1


d <- read_csv("../data/data_pruned.csv", show_col_types = FALSE)
tree <- read.nexus("tree_families_replaced.nex")
tree <- reorder(tree, "postorder")

# Align to tree order
d <- d %>% filter(glottocode %in% tree$tip.label)
d <- d[match(tree$tip.label, d$glottocode), ]

# Build data for Stan

stan_data <- list(
  Ndata        = nrow(d),
  Nnodes       = length(unique(as.vector(tree$edge))),
  Nedges       = nrow(tree$edge),
  edges        = tree$edge,
  edge_lengths = tree$edge.length,
  root_node     = length(tree$tip.label) + 1
)
storage.mode(stan_data$edges) <- "integer"

# Floor zero-length edges (also guarded in Stan)
eps <- 1e-9
stan_data$edge_lengths[stan_data$edge_lengths == 0] <- eps



if (!dir.exists("sensitivity_analysis")) {
    dir.create("sensitivity_analysis")
}



# __MU_PRIOR__
# __MU_PRIOR_GENERATOR__
# __V_PRIOR__
# __V_PRIOR_GENERATOR__
# __TAU_PRIOR__
# __TAU_PRIOR_GENERATOR__
# __T_HALF_PRIOR__
# __T_HALF_PRIOR_GENERATOR__
# __TRANSFORMATION__
# __Y_PRIOR_GENERATOR__
# __Y_REP_GENERATOR__
# __Y_LIKELIHOOD__

baseline <- list(
    name = "baseline",
    mu_prior = "normal(0, 10)",
    mu_prior_generator = "normal_rng(0, 10)",
    V_prior = "lognormal(0, 1)",
    V_prior_generator = "lognormal_rng(0, 1)",
    tau_prior = "lognormal(-1, 0.5)",
    tau_prior_generator = "lognormal_rng(-1, 0.5)",
    t_half_prior = "lognormal(2, 0.5)",
    t_half_prior_generator = "lognormal_rng(2, 0.5)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)

# exponential on V

model_01 <- list(
    name = "model_01",
    mu_prior = "normal(1, 0.5)",
    mu_prior_generator = "normal_rng(1, 0.5)",
    V_prior = "exponential(1)",
    V_prior_generator = "exponential_rng(1)",
    tau_prior = "lognormal(1, 0.5)",
    tau_prior_generator = "lognormal_rng(1, 0.5)",
    t_half_prior = "lognormal(2, 0.5)",
    t_half_prior_generator = "lognormal_rng(2, 0.5)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)


# half-normal on tau, wide 

model_02 <- list(
    name = "model_02",
    mu_prior = "normal(1, 0.5)",
    mu_prior_generator = "normal_rng(1, 0.5)",
    V_prior = "lognormal(0, 1)",
    V_prior_generator = "lognormal_rng(0, 1)",
    tau_prior = "normal(0, 1) T[0,]",
    tau_prior_generator = "abs(normal_rng(0, 1))",
    t_half_prior = "lognormal(2, 0.5)",
    t_half_prior_generator = "lognormal_rng(2, 0.5)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)

# half-normal on tau, narrow 

model_03 <- list(
    name = "model_03",
    mu_prior = "normal(1, 0.5)",
    mu_prior_generator = "normal_rng(1, 0.5)",
    V_prior = "lognormal(0, 1)",
    V_prior_generator = "lognormal_rng(0, 1)",
    tau_prior = "normal(0, 0.15) T[0,]",
    tau_prior_generator = "abs(normal_rng(0, 0.15))",
    t_half_prior = "lognormal(2, 0.5)",
    t_half_prior_generator = "lognormal_rng(2, 0.5)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)



# half-t on tau

model_04 <- list(
    name = "model_04",
    mu_prior = "normal(1, 0.5)",
    mu_prior_generator = "normal_rng(1, 0.5)",
    V_prior = "lognormal(0, 1)",
    V_prior_generator = "lognormal_rng(0, 1)",
    tau_prior = "student_t(3, 0, 0.2) T[0,]",
    tau_prior_generator = "abs(student_t_rng(3, 0, 0.2))",
    t_half_prior = "lognormal(2, 0.5)",
    t_half_prior_generator = "lognormal_rng(2, 0.5)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)


# uniform on t_half

model_05 <- list(
    name = "model_05",
    mu_prior = "normal(1, 0.5)",
    mu_prior_generator = "normal_rng(1, 0.5)",
    V_prior = "lognormal(0, 1)",
    V_prior_generator = "lognormal_rng(0, 1)",
    tau_prior = "exponential(10)",
    tau_prior_generator = "exponential_rng(10)",
    t_half_prior = "uniform(0.1, 100)",
    t_half_prior_generator = "uniform_rng(0.1, 100)",
    transformation = "log(y[i] + 0.001)",
    y_prior_generator = "exp(normal_rng(z_prior[i], tau_prior))-0.001",
    y_rep_generator = "exp(normal_rng(z[i], tau))-0.001",
    y_likelihood = "lognormal_lpdf(y[i] + eps | z[i], tau)"
)

as_identity <- function(m, newname) {
  m$name              <- newname
  m$transformation    <- "y[i]"
  m$y_prior_generator <- "normal_rng(z_prior[i], tau_prior)"
  m$y_rep_generator   <- "normal_rng(z[i], tau)"
  m$mu_prior          <- "normal(2, 1)"
  m$mu_prior_generator <- "normal_rng(2, 1)"
  m$y_likelihood     <- "normal_lpdf(y_transformed[i] | z[i], tau)"
  m
}

model_06 <- as_identity(baseline, "model_06")
model_07 <- as_identity(model_01, "model_07")
model_08 <- as_identity(model_02, "model_08")
model_09 <- as_identity(model_03, "model_09")
model_10 <- as_identity(model_04, "model_10")
model_11 <- as_identity(model_05, "model_11")

models <- list(
    baseline = baseline,
    model_01 = model_01,
    model_02 = model_02,
    model_03 = model_03,
    model_04 = model_04,
    model_05 = model_05,
    model_06 = model_06,
    model_07 = model_07,
    model_08 = model_08,
    model_09 = model_09,
    model_10 = model_10,
    model_11 = model_11
)



stan_data$y <- d %>% pull(str_c(trait, "_markedness_fullness"))

model_text <- readLines("sensitivity_analysis/OU_template.stan")

bundle_to_placeholders <- function(bundle) {
  c(
    "__MU_PRIOR__"               = bundle$mu_prior,
    "__MU_PRIOR_GENERATOR__"     = bundle$mu_prior_generator,
    "__V_PRIOR__"                = bundle$V_prior,
    "__V_PRIOR_GENERATOR__"      = bundle$V_prior_generator,
    "__TAU_PRIOR__"              = bundle$tau_prior,
    "__TAU_PRIOR_GENERATOR__"    = bundle$tau_prior_generator,
    "__T_HALF_PRIOR__"           = bundle$t_half_prior,
    "__T_HALF_PRIOR_GENERATOR__" = bundle$t_half_prior_generator,
    "__TRANSFORMATION__"         = bundle$transformation,
    "__Y_REP_GENERATOR__"        = bundle$y_rep_generator,
    "__Y_PRIOR_GENERATOR__"      = bundle$y_prior_generator,
    "__Y_LIKELIHOOD__"           = bundle$y_likelihood
  )
}


run_model <- function(
    model_name,
    stan_data_in = stan_data,
    model_text_in = model_text,
    trait_in = trait
) {
    stan_code <- str_replace_all(
        paste(model_text_in, collapse = "\n"),
        bundle_to_placeholders(model_name)
    )
    subdir <- file.path("sensitivity_analysis", trait, model_name$name)
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

    stan_file <- str_c(subdir, "/modelcode.stan")
    writeLines(stan_code, stan_file)
    mod <- cmdstan_model(stan_file)
    fit <- mod$sample(
        data = stan_data_in,
        iter_warmup     = 2000,
        iter_sampling   = 2000,
        chains          = 4,
        parallel_chains = 4,
        seed            = 123,
    )

    fit$save_object(str_c(subdir, "/fit.rds"))
    jsonlite::write_json(model_name, str_c(subdir, "/bundle.json"), pretty = TRUE, auto_unbox = TRUE)

    summ <- fit$summary()
    readr::write_csv(summ, file.path(subdir, "summary.csv"))

    log_lik <- fit$draws("log_lik", format = "draws_matrix")

    loo_res <- loo::loo(log_lik)

    # 1) Overall LOO estimates table
    loo_est <- as.data.frame(loo_res$estimates)
    loo_est <- tibble::rownames_to_column(loo_est, var = "metric")
    readr::write_csv(loo_est, file.path(subdir, "loo_estimates.csv"))

    # 2) Pointwise elpd/p_loo etc. (+ Pareto k)
    pw <- as.data.frame(loo_res$pointwise)
    if (!is.null(loo_res$diagnostics$pareto_k)) {
        pw$pareto_k <- loo_res$diagnostics$pareto_k
    }
    readr::write_csv(pw, file.path(subdir, "loo_pointwise.csv"))

    # 3) Keep the full object for later comparisons
    saveRDS(loo_res, file.path(subdir, "loo.rds"))

    # Optional quick text summary
    bad_k_07 <- mean(loo_res$diagnostics$pareto_k > 0.7)
    bad_k_10 <- mean(loo_res$diagnostics$pareto_k > 1.0)
    writeLines(
        sprintf("bad_k>0.7=%.3f, bad_k>1=%.3f", bad_k_07, bad_k_10),
        file.path(subdir, "loo_summary.txt")
    )
    invisible(subdir)
}




run_model(models[[model_no]], stan_data, model_text, trait)

