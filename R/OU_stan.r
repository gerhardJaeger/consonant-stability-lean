## Global OU model fitting

#| vscode: {languageId: r}
library(tidyverse)
library(ape)
library(cmdstanr)
library(posterior)
library(ggdist)
library(loo)
options(repr.plot.width = 20, repr.plot.height = 12.4)
theme_set(theme_minimal(base_size = 24))

#| vscode: {languageId: r}
d <- read_csv("../data/data_pruned.csv", show_col_types = FALSE)
tree <- read.nexus("tree_families_replaced.nex")
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
mod <- cmdstan_model("OU.stan")
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

mod_logtransformed <- cmdstan_model("OU_logtransformed.stan")
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
mod <- cmdstan_model("OU.stan")
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

mod_logtransformed <- cmdstan_model("OU_logtransformed.stan")
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
mod <- cmdstan_model("OU.stan")
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

mod_logtransformed <- cmdstan_model("OU_logtransformed.stan")
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
mod <- cmdstan_model("OU.stan")
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

mod_logtransformed <- cmdstan_model("OU_logtransformed.stan")
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
ggsave("../figures/halflive.svg", plot = p, width = 18, height = 12, device = "svg")

#| vscode: {languageId: r}
draws_combined %>%
  group_by(variable) %>%
  summarize(
    hdi = list(median_hdi(t_half)),
    .groups = "drop"
  ) %>%
  unnest(hdi)

