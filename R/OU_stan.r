## Global OU model fitting

#| vscode: {languageId: r}
library(tidyverse)
library(ape)
library(cmdstanr)
library(posterior)
library(ggdist)
library(patchwork)
library(loo)
options(repr.plot.width = 20, repr.plot.height = 12.4)
theme_set(theme_minimal(base_size = 24))

#| vscode: {languageId: r}
d <- read_csv("data/data_pruned.csv", show_col_types = FALSE)
tree <- read.nexus("data/tree_families_replaced.nex")
tree <- reorder(tree, "postorder")

# Align to tree order
d <- d %>% filter(glottocode %in% tree$tip.label)
d <- d[match(tree$tip.label, d$glottocode), ]

#| vscode: {languageId: r}
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

#| vscode: {languageId: r}

## Plosives

stan_data$y <- d$plosive_markedness_fullness

# Compile & sample
mod <- cmdstan_model("stan/OU.stan")
fit_plosive <- mod$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_plosive$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}

fit_plosive$diagnostic_summary()

#| vscode: {languageId: r}

fit_plosive$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_plosive$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_plosive  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_plosive)

#| vscode: {languageId: r}

## plot
draws_plosive <- posterior::as_draws_df(fit_plosive$draws("t_half"))

ggplot(draws_plosive, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (median and 95% HDI indicated)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Plosives: log-transformed

#| vscode: {languageId: r}

##  log-transformed

mod_logtransformed <- cmdstan_model("stan/OU_logtransformed.stan")
fit_logplosive <- mod_logtransformed$sample(
  data = stan_data,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_logplosive$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}



fit_logplosive$diagnostic_summary()

#| vscode: {languageId: r}

fit_logplosive$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_logplosive$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_logplosive  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_logplosive)

loo_compare(loo_plosive, loo_logplosive)

#| vscode: {languageId: r}


## Plot

draws_logplosive <- posterior::as_draws_df(fit_logplosive$draws("t_half"))

ggplot(draws_logplosive, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (median and 95% HDI indicated)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Fricatives: untransformed

#| vscode: {languageId: r}


################
## Fricatives 
###############

stan_data$y <- d$fricative_markedness_fullness

# Compile & sample
mod <- cmdstan_model("stan/OU.stan")
fit_fricative <- mod$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_fricative$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}


fit_fricative$diagnostic_summary()

#| vscode: {languageId: r}

fit_fricative$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_fricative$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_fricative  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_fricative)

#| vscode: {languageId: r}
## plot
draws_fricative <- posterior::as_draws_df(fit_fricative$draws("t_half"))


ggplot(draws_fricative, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (95% HDI shaded)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Fricatives: log-transformed

#| vscode: {languageId: r}

##  log-transformed

mod_logtransformed <- cmdstan_model("stan/OU_logtransformed.stan")
fit_logfricative <- mod_logtransformed$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

fit_logfricative$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}


fit_logfricative$diagnostic_summary()

#| vscode: {languageId: r}

fit_logfricative$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_logfricative$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_logfricative  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_logfricative)

loo_compare(loo_fricative, loo_logfricative)

#| vscode: {languageId: r}


## Plot

draws_logfricative <- posterior::as_draws_df(fit_logfricative$draws("t_half"))

ggplot(draws_logfricative, aes(x = t_half, y = "")) +
  stat_halfeye(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (95% HDI shaded)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Affricates: untransformed

#| vscode: {languageId: r}

################
## Affricates
###############

stan_data$y <- d$affricate_markedness_fullness

# Compile & sample
mod <- cmdstan_model("stan/OU.stan")
fit_affricate <- mod$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_affricate$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}


fit_affricate$diagnostic_summary()

#| vscode: {languageId: r}

fit_affricate$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_affricate$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_affricate  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_affricate)

#| vscode: {languageId: r}

## plot
draws_affricate <- posterior::as_draws_df(fit_affricate$draws("t_half"))

ggplot(draws_affricate, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (95% HDI shaded)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Affricates: log-transformed

#| vscode: {languageId: r}


##  log-transformed

mod_logtransformed <- cmdstan_model("stan/OU_logtransformed.stan")
fit_logaffricate <- mod_logtransformed$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_logaffricate$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}


fit_logaffricate$diagnostic_summary()

#| vscode: {languageId: r}

fit_logaffricate$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_logaffricate$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_logaffricate  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_logaffricate)

loo_compare(loo_affricate, loo_logaffricate)

#| vscode: {languageId: r}

## Plot

draws_logaffricate <- posterior::as_draws_df(fit_logaffricate$draws("t_half"))

ggplot(draws_logaffricate, aes(x = t_half, y = "")) +
  stat_halfeye(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (95% HDI shaded)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Obstruent series: untransformed

#| vscode: {languageId: r}

################
## Series
###############

stan_data$y <- d$series_markedness_fullness

# Compile & sample
mod <- cmdstan_model("stan/OU.stan")
fit_series <- mod$sample(
  data = stan_data,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

# Summaries (omit the big arrays; here only globals of interest)
fit_series$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}


fit_series$diagnostic_summary()

#| vscode: {languageId: r}

fit_series$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_series$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_series  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_series)

#| vscode: {languageId: r}

## plot
draws_series <- posterior::as_draws_df(fit_series$draws("t_half"))

ggplot(draws_series, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (95% HDI shaded)", x = "t_half", y = NULL) +
  theme_minimal()

## ## Obstruent series: log-transformed

#| vscode: {languageId: r}



##  log-transformed

mod_logtransformed <- cmdstan_model("stan/OU_logtransformed.stan")
fit_logseries <- mod_logtransformed$sample(
  data = stan_data,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  chains          = 4,
  parallel_chains = 4,
  seed            = 123,
  adapt_delta     = 0.95,
  max_treedepth   = 13
)

fit_logseries$summary(variables = c("mu", "sigma", "lambda", "tau", "t_half"))

#| vscode: {languageId: r}

fit_logseries$diagnostic_summary()

#| vscode: {languageId: r}

fit_logseries$cmdstan_diagnose()

#| vscode: {languageId: r}

## LOO

ll_draws <- fit_logseries$draws("log_lik")
ll_array <- posterior::as_draws_array(ll_draws)

loo_logseries  <- loo::loo(ll_array, cores = parallel::detectCores())
print(loo_logseries)

loo_compare(loo_series, loo_logseries)

#| vscode: {languageId: r}

## Plot

draws_logseries <- posterior::as_draws_df(fit_logseries$draws("t_half"))

ggplot(draws_plosive, aes(x = t_half, y = "")) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  labs(title = "t_half posterior (median and 95% HDI indicated)", x = "t_half", y = NULL) +
  theme_minimal()

#| vscode: {languageId: r}
draws_combined <- bind_rows(
    draws_logplosive %>% mutate(variable = "plosive"),
    draws_fricative %>% mutate(variable = "fricative"),
    draws_affricate %>% mutate(variable = "affricate"),
    draws_logseries %>% mutate(variable = "series")
)
draws_combined <- draws_combined %>%
  mutate(variable = factor(variable,
                           levels = c("plosive", "fricative", "affricate", "series")))

#| vscode: {languageId: r}
custom_colors <- c(plosive = "#4C72B0", fricative = "#55A868", affricate = "#E69F00", series = "#CC79A7")
p <- ggplot(draws_combined, aes(x = t_half, y = variable, fill = variable)) +
  stat_slabinterval(.width = 0.95, point_interval = median_hdi) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "phylogenetic half-life", x = "kya", y = NULL) +
  theme_gray(base_size = 24) +
  theme(legend.position = "none")
p
ggsave("figures/halflive.svg", plot = p, width = 18, height = 12, device = "svg")

#| vscode: {languageId: r}
draws_combined %>%
  group_by(variable) %>%
  summarize(
    hdi = list(median_hdi(t_half)),
    .groups = "drop"
  ) %>%
  unnest(hdi)

## ASR export

backtransform <- function(x) exp(x) + 0.001

dir.create("data/asr", recursive = TRUE, showWarnings = FALSE)

asr_plosive <- posterior::as_draws_df(fit_logplosive$draws("z")) %>%
  select(-.chain, -.iteration, -.draw) %>%
  t() %>%
  backtransform() %>%
  as.data.frame()

write_csv(asr_plosive, "data/asr/asr_plosive.csv", col_names = TRUE)

asr_fricative <- posterior::as_draws_df(fit_fricative$draws("z")) %>%
  select(-.chain, -.iteration, -.draw) %>%
  t() %>%
  as.data.frame()

write_csv(asr_fricative, "data/asr/asr_fricative.csv", col_names = TRUE)

asr_affricate <- posterior::as_draws_df(fit_logaffricate$draws("z")) %>%
  select(-.chain, -.iteration, -.draw) %>%
  t() %>%
  as.data.frame()

write_csv(asr_affricate, "data/asr/asr_affricate.csv", col_names = TRUE)

asr_series <- posterior::as_draws_df(fit_logseries$draws("z")) %>%
  select(-.chain, -.iteration, -.draw) %>%
  t() %>%
  backtransform() %>%
  as.data.frame()

write_csv(asr_series, "data/asr/asr_series.csv", col_names = TRUE)

## Equilibrium distributions (landscapes)

# Compute density matrix across posterior draws for transformed models (Y = exp(X) - 0.001)
compute_transformed_density_matrix <- function(fit, y_grid, n_draws = 200) {
  draws_eq <- posterior::as_draws_df(fit$draws(c("mu", "V")))
  if (!"V" %in% names(draws_eq)) stop("Parameter 'V' not found in fit draws for transformed plot")
  mu_draws <- draws_eq$mu
  V_draws  <- draws_eq$V
  S <- length(mu_draws)
  sel <- sample(seq_len(S), size = min(n_draws, S))
  x_grid <- log(y_grid + 0.001)
  dens_mat <- sapply(sel, function(i) {
    dnorm(x_grid, mean = mu_draws[i], sd = sqrt(V_draws[i])) * (1 / (y_grid + 0.001))
  })
  list(mat = dens_mat, sel = sel)
}

# Compute density matrix across posterior draws for untransformed models (X ~ Normal(mu, V))
compute_untransformed_density_matrix <- function(fit, x_grid, n_draws = 200) {
  all_draws <- posterior::as_draws_df(fit$draws())
  if ("V" %in% names(all_draws)) {
    draws_eq <- posterior::as_draws_df(fit$draws(c("mu", "V")))
    mu_draws <- draws_eq$mu
    V_draws  <- draws_eq$V
  } else if (all(c("sigma", "lambda") %in% names(all_draws))) {
    draws_eq <- posterior::as_draws_df(fit$draws(c("mu", "sigma", "lambda")))
    mu_draws <- draws_eq$mu
    V_draws  <- (draws_eq$sigma^2) / (2 * draws_eq$lambda)
  } else if (all(c("sigma", "t_half") %in% names(all_draws))) {
    draws_eq <- posterior::as_draws_df(fit$draws(c("mu", "sigma", "t_half")))
    mu_draws <- draws_eq$mu
    lambda_draws <- log(2) / draws_eq$t_half
    V_draws <- (draws_eq$sigma^2) / (2 * lambda_draws)
  } else {
    stop("Required parameters not found in fit draws for untransformed density.")
  }
  S <- length(mu_draws)
  sel <- sample(seq_len(S), size = min(n_draws, S))
  dens_mat <- sapply(sel, function(i) dnorm(x_grid, mean = mu_draws[i], sd = sqrt(V_draws[i])))
  list(mat = dens_mat, sel = sel)
}

# Helper to build a plot from a quantile summary dataframe and apply shared limits
build_plot_from_df <- function(df, obs_vals, fill_col, title = "", xlim = NULL, ylim = NULL) {
  obs_df <- tibble(obs = obs_vals[!is.na(obs_vals)])
  p <- ggplot(df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = fill_col, alpha = 0.12) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = fill_col, alpha = 0.28) +
    geom_line(aes(y = median), color = "black", linewidth = 1) +
    geom_line(aes(y = mean), color = fill_col, linewidth = 0.8, linetype = "dashed") +
    geom_rug(data = obs_df, aes(x = jitter(obs, amount = 0.1)), inherit.aes = FALSE,
             sides = "b", color = "black", alpha = 0.25, linewidth = 0.6) +
    labs(title = title, x = NULL, y = NULL) +
    theme_gray(base_size = 24) +
    theme(panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"),
          panel.background = element_rect(fill = "grey90", colour = NA),
          axis.title.y = element_blank())
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  p
}

# Observations
obs_plosive   <- d$plosive_markedness_fullness
obs_fricative <- d$fricative_markedness_fullness
obs_affricate <- d$affricate_markedness_fullness
obs_series    <- d$series_markedness_fullness

# Shared top-row x-range
top_x_min <- min(c(min(obs_plosive, na.rm = TRUE),
                   min(obs_fricative, na.rm = TRUE),
                   min(obs_affricate, na.rm = TRUE)), na.rm = TRUE) - 0.5
top_x_max <- max(c(max(obs_plosive, na.rm = TRUE),
                   max(obs_fricative, na.rm = TRUE),
                   max(obs_affricate, na.rm = TRUE)), na.rm = TRUE) + 0.5
top_x_min <- max(top_x_min, -0.0009)  # ensure > -0.001 for transformed models

# Grids
y_grid_plosive    <- seq(top_x_min, top_x_max, length.out = 400)
x_grid_fricative  <- y_grid_plosive
y_grid_affricate  <- seq(top_x_min, top_x_max, length.out = 400)
series_min        <- max(min(obs_series, na.rm = TRUE) - 0.5, -0.0009)
series_max        <- max(obs_series, na.rm = TRUE) + 0.5
x_grid_series     <- seq(series_min, series_max, length.out = 600)

# Density matrices and quantiles
dm_plosive   <- compute_transformed_density_matrix(fit_logplosive,  y_grid = y_grid_plosive,   n_draws = 300)
dm_fricative <- compute_untransformed_density_matrix(fit_fricative, x_grid = x_grid_fricative, n_draws = 300)
dm_affricate <- compute_transformed_density_matrix(fit_logaffricate, y_grid = y_grid_affricate, n_draws = 300)
dm_series    <- compute_transformed_density_matrix(fit_logseries,   y_grid = x_grid_series,    n_draws = 300)

qs_plosive    <- apply(dm_plosive$mat,    1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
qs_fricative  <- apply(dm_fricative$mat,  1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
qs_affricate  <- apply(dm_affricate$mat,  1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
qs_series     <- apply(dm_series$mat,     1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

mean_plosive   <- rowMeans(dm_plosive$mat)
mean_fricative <- rowMeans(dm_fricative$mat)
mean_affricate <- rowMeans(dm_affricate$mat)
mean_series    <- rowMeans(dm_series$mat)

df_plosive <- tibble(
  x = y_grid_plosive,
  lower95 = qs_plosive[1,], lower50 = qs_plosive[2,],
  median  = qs_plosive[3,], upper50 = qs_plosive[4,],
  upper95 = qs_plosive[5,], mean = mean_plosive
)

df_fricative <- tibble(
  x = x_grid_fricative,
  lower95 = qs_fricative[1,], lower50 = qs_fricative[2,],
  median  = qs_fricative[3,], upper50 = qs_fricative[4,],
  upper95 = qs_fricative[5,], mean = mean_fricative
)

df_affricate <- tibble(
  x = y_grid_affricate,
  lower95 = qs_affricate[1,], lower50 = qs_affricate[2,],
  median  = qs_affricate[3,], upper50 = qs_affricate[4,],
  upper95 = qs_affricate[5,], mean = mean_affricate
)

df_series <- tibble(
  x = x_grid_series,
  lower95 = qs_series[1,], lower50 = qs_series[2,],
  median  = qs_series[3,], upper50 = qs_series[4,],
  upper95 = qs_series[5,], mean = mean_series
)

# Global y-max across all panels (use upper95) + small margin
y_max_global <- max(
  max(df_plosive$upper95,   na.rm = TRUE),
  max(df_fricative$upper95, na.rm = TRUE),
  max(df_affricate$upper95, na.rm = TRUE),
  max(df_series$upper95,    na.rm = TRUE),
  na.rm = TRUE
) * 1.05

# Build plots with shared limits
p_plosive <- build_plot_from_df(
  df_plosive, obs_plosive, fill_col = "#4C72B0",
  title = "plosive", xlim = c(top_x_min, top_x_max), ylim = c(0, 1)
)

p_fricative <- build_plot_from_df(
  df_fricative, obs_fricative, fill_col = "#55A868",
  title = "fricative", xlim = c(top_x_min, top_x_max), ylim = c(0, 1)
)

p_affricate <- build_plot_from_df(
  df_affricate, obs_affricate, fill_col = "#E69F00",
  title = "affricate", xlim = c(top_x_min, top_x_max), ylim = c(0, 1)
)

p_series <- build_plot_from_df(
  df_series, obs_series, fill_col = "#CC79A7",
  title = "series", xlim = c(series_min, series_max), ylim = c(0, 1)
)

combined_plot <- (p_plosive | p_fricative | p_affricate) / p_series +
  plot_annotation(theme = theme(text = element_text(size = 18)))

combined_plot

ggsave("figures/landscapes.svg", combined_plot, width = 20, height = 12.4)
