library(tidyverse)

sim_dir <- 'data/simulations_fitted/'

files <- list.files(sim_dir, pattern = '^OU_sd_.*\\.summary\\.csv$', full.names = TRUE)

sim_data <- map_dfr(files, function(f) {
  fname  <- basename(f)
  sd_val <- as.numeric(str_match(fname, 'sd_([0-9.]+)')[,2])
  th_val <- as.numeric(str_match(fname, 'thalf_([0-9.]+)')[,2])
  read_csv(f, show_col_types = FALSE) %>%
    filter(variable == 't_half') %>%
    transmute(true_thalf = th_val, sigma = sd_val,
              est_median = median, est_q5 = q5, est_q95 = q95)
})

# empirical range (from main paper)
emp_lo <- 13.3
emp_hi <- 26.1

# reference data for y=x line
ref <- tibble(x = c(0.4, 110))

# highlight sigma = 1 (used in supplementary Table 3); show others faintly
sim_data <- sim_data %>%
  mutate(highlight = sigma == 1,
         sigma_fac = factor(sigma))

p <- ggplot(sim_data, aes(x = true_thalf, y = est_median, group = sigma_fac)) +

  # empirical range shading
  annotate('rect', xmin = emp_lo, xmax = emp_hi,
           ymin = -Inf, ymax = Inf,
           fill = '#4dac26', alpha = 0.10) +
  annotate('segment', x = emp_lo, xend = emp_lo,
           y = 0.3, yend = 115,
           linetype = 'dashed', color = '#4dac26', linewidth = 0.4) +
  annotate('segment', x = emp_hi, xend = emp_hi,
           y = 0.3, yend = 115,
           linetype = 'dashed', color = '#4dac26', linewidth = 0.4) +
  annotate('text', x = sqrt(emp_lo * emp_hi), y = 110,
           label = 'empirical\nrange', color = '#4dac26',
           size = 2.8, vjust = 1, hjust = 0.5, lineheight = 0.9) +

  # y = x reference
  geom_line(data = ref, aes(x = x, y = x), inherit.aes = FALSE,
            color = 'grey40', linetype = 'solid', linewidth = 0.5) +

  # faint lines for all sigma
  geom_line(data = filter(sim_data, !highlight),
            color = 'grey70', linewidth = 0.35) +

  # error bars for sigma = 1
  geom_errorbar(data = filter(sim_data, highlight),
                aes(ymin = est_q5, ymax = est_q95),
                width = 0.06, color = '#2166ac', linewidth = 0.5, alpha = 0.7) +

  # line for sigma = 1
  geom_line(data = filter(sim_data, highlight),
            color = '#2166ac', linewidth = 0.8) +

  # points for sigma = 1
  geom_point(data = filter(sim_data, highlight),
             color = '#2166ac', size = 2) +

  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c('0.5', '1', '2', '5', '10', '20', '50', '100')
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c('0.5', '1', '2', '5', '10', '20', '50', '100')
  ) +
  coord_cartesian(xlim = c(0.4, 120), ylim = c(0.4, 120)) +
  labs(
    x = expression('True'~italic(t)[half]~'(kya)'),
    y = expression('Estimated'~italic(t)[half]~'(kya, posterior median)')
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = 'none'
  )

# annotation for reference line and sigma=1 line
p <- p +
  annotate('text', x = 0.42, y = 100, label = 'y = x',
           color = 'grey40', size = 3, hjust = 0, vjust = -0.3, fontface = 'italic') +
  annotate('text', x = 0.42, y = 65,
           label = expression(sigma == 1~'(90% CI)'),
           color = '#2166ac', size = 3, hjust = 0) +
  annotate('text', x = 0.42, y = 42,
           label = 'other \u03c3 values',
           color = 'grey60', size = 3, hjust = 0)

ggsave('figures/simulation_recovery.pdf', p,
       width = 5, height = 4.5, device = cairo_pdf)
ggsave('figures/simulation_recovery.svg', p,
       width = 5, height = 4.5)
