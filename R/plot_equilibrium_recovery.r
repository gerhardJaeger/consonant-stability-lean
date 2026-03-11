library(tidyverse)
library(patchwork)

sim_dir <- 'data/simulations_fitted/'

# ---- load all simulation files ------------------------------------------
read_sim <- function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  list(
    mu    = filter(d, variable == 'mu'),
    sigma = filter(d, variable == 'sigma'),
    thalf = filter(d, variable == 't_half')
  )
}

# OU files
ou_files <- list.files(sim_dir, pattern = '^OU_sd_.*\\.summary\\.csv$',
                       full.names = TRUE)

ou_data <- map_dfr(ou_files, function(f) {
  fname  <- basename(f)
  sd_val <- as.numeric(str_match(fname, 'sd_([0-9.]+)')[,2])
  th_val <- as.numeric(str_match(fname, 'thalf_([0-9.]+)')[,2])
  s <- read_sim(f)
  tibble(
    true_sd    = sd_val,
    true_V     = sd_val^2,
    true_thalf = th_val,
    mu_med     = s$mu$median,
    sigma_med  = s$sigma$median,
    thalf_med  = s$thalf$median,
    V_est      = s$sigma$median^2 * s$thalf$median / (2 * log(2))
  )
})

# BBM files
bbm_files <- list.files(sim_dir, pattern = '^BBM_sig_.*\\.summary\\.csv$',
                        full.names = TRUE)

bbm_data <- map_dfr(bbm_files, function(f) {
  s <- read_sim(f)
  tibble(
    mu_med    = s$mu$median,
    sigma_med = s$sigma$median,
    thalf_med = s$thalf$median,
    V_est     = s$sigma$median^2 * s$thalf$median / (2 * log(2))
  )
}) %>% mutate(id = row_number())

# true mu and V: use the most tightly-estimated case (sd=1, small t_half)
ref_row <- ou_data %>% filter(true_sd == 1, true_thalf == 2)
true_mu <- ref_row$mu_med
true_sd_val <- 1.0   # V_true = sd^2 = 1

# =========================================================================
# PANEL A: OU simulations — true vs estimated equilibrium density
# =========================================================================
x_vals <- seq(-0.5, 11, length.out = 600)

# True equilibrium: N(true_mu, true_V) for sd=1 → V=1
true_curve <- tibble(
  x = x_vals,
  y = dnorm(x_vals, mean = true_mu, sd = 1),
  type = 'True equilibrium'
)

# Estimated equilibria for sd=1, selected t_half values
selected_th <- c(2, 10, 20, 50, 100)
palette_ou  <- colorRampPalette(c('#9ecae1', '#08306b'))(length(selected_th))

ou_curves <- map2_dfr(selected_th, palette_ou, function(th, col) {
  row <- filter(ou_data, true_sd == 1, true_thalf == th)
  tibble(
    x     = x_vals,
    y     = dnorm(x_vals, mean = row$mu_med, sd = sqrt(row$V_est)),
    label = factor(str_c(th, ' kya'), levels = str_c(selected_th, ' kya')),
    col   = col
  )
})

p_ou <- ggplot() +
  geom_line(data = ou_curves,
            aes(x = x, y = y, color = label, group = label),
            linewidth = 0.8, alpha = 0.5) +
  geom_function(fun = dnorm, args = list(mean = true_mu, sd = true_sd_val),
                color = '#c0392b', linetype = 'dashed', linewidth = 1.3) +
  annotate('text', x = true_mu + 1.6,
           y = dnorm(true_mu + 1.6, true_mu, true_sd_val) + 0.03,
           label = 'true equilibrium', hjust = 0, size = 3,
           color = '#c0392b', fontface = 'italic') +
  scale_color_manual(
    values = setNames(palette_ou, str_c(selected_th, ' kya')),
    name   = expression(italic(t)[half])
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, NA)) +
  labs(
    title = 'OU simulations',
    x     = 'Obstruent series count',
    y     = 'Equilibrium density'
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title        = element_text(face = 'bold', size = 10),
    panel.grid.minor  = element_blank(),
    legend.position   = 'inside',
    legend.position.inside = c(0.98, 0.98),
    legend.justification   = c(1, 1),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7.5),
    legend.background = element_rect(fill = alpha('white', 0.8), color = NA)
  )

# =========================================================================
# PANEL B: BBM simulations — inferred vs theoretical uniform
# =========================================================================
B_max <- 9   # max obstruent series count

bbm_curves <- map_dfr(seq_len(nrow(bbm_data)), function(i) {
  row <- bbm_data[i, ]
  tibble(
    x   = x_vals,
    y   = dnorm(x_vals, mean = row$mu_med, sd = sqrt(row$V_est)),
    id  = i
  )
})

unif_curve <- tibble(
  x = c(0, 0, B_max, B_max),
  y = c(0, 1/B_max, 1/B_max, 0)
)

p_bbm <- ggplot() +
  geom_line(data = bbm_curves,
            aes(x = x, y = y, group = id),
            color = '#d6604d', linewidth = 0.6, alpha = 0.6) +
  geom_line(data = unif_curve,
            aes(x = x, y = y),
            color = 'black', linetype = 'dashed', linewidth = 1) +
  annotate('text', x = B_max / 2, y = 1/B_max + 0.01,
           label = 'theoretical uniform', hjust = 0.5, size = 3,
           color = 'black', fontface = 'italic') +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, NA)) +
  labs(
    title = 'BBM simulations',
    x     = 'Obstruent series count',
    y     = 'Equilibrium density'
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = 'bold', size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = 'none'
  )

# =========================================================================
# Combine and save
# =========================================================================
p_combined <- (p_ou | p_bbm) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 11))

ggsave('figures/equilibrium_recovery.pdf', p_combined,
       width = 8.5, height = 4, device = cairo_pdf)
ggsave('figures/equilibrium_recovery.svg', p_combined,
       width = 8.5, height = 4)
