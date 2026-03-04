library(tidyverse)
library(ggdist)


results <- tibble(
  variable = character(),
  trait = character(),
  model_no = numeric(),
  statistic = character(),
  value = numeric()
)



traits <- c('plosive', 'fricative', 'affricate', 'series')

vars <- c("V", "t_half", "mu", "tau")


for (trait in traits) {
    print(trait)
    models <- list.files(str_c('sensitivity_analysis/', trait, '/'), full.names=TRUE)

    for (model_no in seq_along(models)) {
        print(model_no)
        fit <- readRDS(str_c(models[model_no], '/fit.rds'))
        fit_summary <- fit$summary(variables = vars)
        draws <- fit$draws(vars, format = "draws_df")
        draws_long <- draws %>%
        pivot_longer(
            cols = -c(.chain, .iteration, .draw),  # keep metadata columns
            names_to = "variable",
            values_to = "value"
        )
        local_results <- draws_long %>% 
            group_by(variable) %>%
            summarise(
                median = median(value),
                hdi_95_lower = hdi(value, .width = 0.95)[1],
                hdi_95_upper = hdi(value, .width = 0.95)[2]
                ) %>%
            mutate(trait = trait) %>%
            mutate(model_no = model_no) %>%
            pivot_longer(
                cols = c(median, hdi_95_lower, hdi_95_upper),
                names_to = "statistic",
                values_to = "value"
            )
        loo_estimates <- read_csv(str_c(models[model_no], '/loo_estimates.csv')) %>%
            filter(metric == 'looic') %>%
            mutate(trait = trait) %>%
            mutate(model_no = model_no) %>%
            mutate(variable = 'looic') %>%
            pivot_longer(
                cols = c('Estimate', 'SE'),
                names_to = "statistic",
                values_to = "value"
            ) %>%
            select(-metric)
        results <- bind_rows(results, local_results, loo_estimates)
    }
}

write_csv(results, 'sensitivity_analysis_results.csv')
