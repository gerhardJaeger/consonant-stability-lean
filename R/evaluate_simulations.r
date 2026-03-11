library(tidyverse)
library(stringr)

bbm_files <- list.files("data/simulations_fitted/", full.names = TRUE) %>%
  keep(~ str_detect(.x, "BBM"))

ou_files <- list.files("data/simulations_fitted/", full.names = TRUE) %>%
  keep(~ str_detect(.x, "OU"))

bbm_params <- lapply(
    bbm_files,
    function(fn) as.numeric(str_match(fn, ".*BBM_sig_(.+).summary.csv")[2])) %>%
    unlist() %>%
    tibble(sd = .) %>%
    mutate(idx = row_number(), .before = sd)

ou_params <- lapply(
    ou_files,
    function(fn) str_match(fn, ",*OU_sd_(.+)_thalf_(.+).summary.csv")[,-1] %>% as.numeric()
) %>%
    map_dfr(~ as_tibble_row(set_names(.x, c("sd", "t_half")))) %>%
    mutate(idx = row_number(), .before = sd)

results_bbm <- tibble(
  variable = character(),
  name = character(),
  value = numeric()
)

for (i in seq_along(bbm_files)) {
    read_csv(bbm_files[i], show_col_types = FALSE) %>%
        mutate(idx = i) %>%
        select(idx, variable, mean, median, sd, mad, q5, q95) %>%
        pivot_longer(
        cols = -c(idx, variable)
        ) %>%
        rbind(
        results_bbm,
        .
        ) -> results_bbm
}
bbm_params %>%
  inner_join(results_bbm, by = "idx") -> results_bbm

results_ou <- tibble(
  variable = character(),
  name = character(),
  value = numeric()
)

for (i in seq_along(ou_files)) {
    read_csv(ou_files[i], show_col_types = FALSE) %>%
        mutate(idx = i) %>%
        select(idx, variable, mean, median, sd, mad, q5, q95) %>%
        pivot_longer(
        cols = -c(idx, variable)
        ) %>%
        rbind(
        results_ou,
        .
        ) -> results_ou
}
ou_params %>%
  inner_join(results_ou, by = "idx") -> results_ou

results_ou %>%
    filter(variable == "t_half", name == "mean") %>%
    ggplot(aes(x = sd, y = value)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue")

p_recovery <- results_ou %>%
    filter(variable == "t_half", name == "mean") %>%
    ggplot(aes(x = t_half, y = value)) +
    geom_point() +
    geom_smooth(method = "gam", se = TRUE, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

p_recovery

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/simulation_recovery.svg", p_recovery, width = 10, height = 6, device = "svg")
