library(tidyverse)

results <- read_csv('sensitivity_analysis_results.csv')

best_transformation <- c(
  plosive   = 'log',
  fricative = 'identity',
  affricate = 'log',
  series    = 'log'
)

t_half_data <- results %>%
  filter(variable == 't_half') %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(
    transformation = if_else(model_no <= 6, 'log', 'identity'),
    is_best = transformation == best_transformation[trait],
    is_baseline = (model_no == 1 & transformation == 'log') |
                  (model_no == 7 & transformation == 'identity'),
    trait = factor(trait,
                   levels = c('plosive', 'fricative', 'affricate', 'series'),
                   labels = c('Plosive', 'Fricative', 'Affricate', 'Series')),
    model_label = factor(str_c('M', model_no), levels = str_c('M', 12:1))
  )

# baseline medians for reference lines (best-transformation model 1 or 7)
baselines <- t_half_data %>%
  filter(is_baseline, is_best) %>%
  select(trait, median)

cols <- c('log' = '#2166ac', 'identity' = '#d6604d')

p <- ggplot(t_half_data,
            aes(x = median, y = model_label,
                color = transformation, alpha = is_best)) +
  geom_vline(data = baselines, aes(xintercept = median),
             linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
  geom_errorbarh(aes(xmin = hdi_95_lower, xmax = hdi_95_upper),
                 height = 0.35, linewidth = 0.5) +
  geom_point(aes(shape = is_baseline, size = is_baseline)) +
  scale_color_manual(
    values = cols,
    name   = 'Transformation',
    labels = c('identity' = 'Identity (models 7\u201312)',
               'log'      = 'Log (models 1\u20136)')
  ) +
  scale_alpha_manual(values = c('TRUE' = 1, 'FALSE' = 0.45), guide = 'none') +
  scale_shape_manual(values = c('TRUE' = 18, 'FALSE' = 16), guide = 'none') +
  scale_size_manual(values  = c('TRUE' = 3.5, 'FALSE' = 2), guide = 'none') +
  facet_wrap(~trait, scales = 'free_x', ncol = 4) +
  labs(
    x = expression(italic(t)[half] ~ '(kya)'),
    y = NULL
  ) +
  geom_hline(yintercept = 6.5, linetype = 'dotted',
             color = 'grey70', linewidth = 0.4) +
  theme_bw(base_size = 10) +
  theme(
    legend.position   = 'bottom',
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    strip.background  = element_blank(),
    strip.text        = element_text(face = 'bold', size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 7.5)
  )

ggsave('figures/sensitivity_analysis.pdf', p,
       width = 8, height = 4.5, device = cairo_pdf)
ggsave('figures/sensitivity_analysis.svg', p,
       width = 8, height = 4.5)
