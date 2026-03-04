## Generate figures

library(tidyverse)
library(sf)
options(repr.plot.width = 20, repr.plot.height = 12.4)
library(svglite)
library(patchwork)
library(rnaturalearth)

d <- read_csv("../data/data_pruned.csv")
# tree <- read.nexus("../data/tree_pruned.nex")

p1 <- d %>%
    ggplot(aes(x = series_markedness_fullness)) + 
    geom_histogram(bins = 50, aes(y = after_stat(density)), fill = "grey", color = "black") +
    geom_density(fill = "#CC79A7", linewidth=1, alpha = 0.4) +
    labs(x = "series markedness fullness") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +  # Ensure integer x-ticks
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(text = element_text(size = 24))

p2 <- d %>%
    ggplot(aes(x = plosive_markedness_fullness)) + 
    geom_histogram(bins = 50, aes(y = after_stat(density)), fill = "grey", color = "black") +
    geom_density(fill = "#4C72B0", linewidth=1, alpha = 0.4) +
    labs(x = "plosive markedness fullness") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(text = element_text(size = 24))

p3 <- d %>%
    ggplot(aes(x = fricative_markedness_fullness)) + 
    geom_histogram(bins = 50, aes(y = after_stat(density)), fill = "grey", color = "black") +
    geom_density(fill = "#55A868", linewidth=1, alpha = 0.4) + # Fixed color
    labs(x = "fricative markedness fullness") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(text = element_text(size = 24)) + 
    theme(axis.title.y = element_blank())

p4 <- d %>%
    ggplot(aes(x = affricate_markedness_fullness)) + 
    geom_histogram(bins = 50, aes(y = after_stat(density)), fill = "grey", color = "black") +
    geom_density(fill = "#E69F00", linewidth=1, alpha = 0.4) +
    labs(x = "affricate markedness fullness") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(text = element_text(size = 24)) + 
    theme(axis.title.y = element_blank())


combined_plot <- (p2 | p3 | p4) / p1
ggsave("../figures/histograms.svg", combined_plot, width = 20, height = 12.4)
combined_plot

create_violin_plot <- function(data, y_var, title, use_pseudo_log = FALSE) {
  plot <- data %>%
    filter(!is.na(.data[[y_var]])) %>% # Filter out missing values
    mutate(Macroarea = fct_reorder(Macroarea, !!sym(y_var), .fun = median)) %>%
    ggplot(aes(x = Macroarea, y = .data[[y_var]])) +
    geom_violin(aes(fill = Macroarea)) +
    scale_fill_viridis_d(option = "plasma") +
    geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
    stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "black") +
    # theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      text = element_text(size = 18),
      panel.grid.major.y = element_line(color = "grey80")
    ) +
    labs(x = NULL, y = title)
  
  # Apply pseudo-log transformation or keep default linear scale
  if (use_pseudo_log) {
    plot <- plot + scale_y_continuous(trans = "pseudo_log")
  }
  
  return(plot)
}

# Create individual plots without y-axis limits
plot_plosive <- create_violin_plot(d, "plosive_markedness_fullness", "Plosive Markedness Fullness", use_pseudo_log = TRUE) +
  theme(legend.position = "none")

plot_fricative <- create_violin_plot(d, "fricative_markedness_fullness", "Fricative Markedness Fullness", use_pseudo_log = TRUE) +
  theme(legend.position = "none")

plot_affricate <- create_violin_plot(d, "affricate_markedness_fullness", "Affricate Markedness Fullness", use_pseudo_log = TRUE) +
  theme(legend.position = "none")

# Larger plot for the bottom row without a legend
plot_series <- create_violin_plot(d, "series_markedness_fullness", "Series Markedness Fullness", use_pseudo_log = TRUE) +
  theme(legend.position = "none") +
  labs(x = "Macro Area")

# Combine plots with adjusted row heights
combined_plot <- (plot_plosive + plot_fricative + plot_affricate) / plot_series +
  plot_layout(heights = c(1, 1.2))

# Display the combined plot
combined_plot

ggsave("../figures/macro_areas.svg", plot = combined_plot, width = 18, height = 12, device = "svg")

d_geo <- d %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Change to Equal Earth projection centered at 160°E
target_crs <- st_crs("+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Define cutting polygon
offset <- 180 - 160
polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

# Cut and transform
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>% 
    st_make_valid() %>% 
    st_difference(polygon) %>% 
    st_transform(crs = target_crs)

ocean <- st_polygon(
    list(cbind(c(
        seq(-20, 339, len = 100),
        rep(339, 100),
        seq(338, -19, len = 100),
        rep(-20, 100)),
        c(rep(-90, 100),
        seq(-89, 89, len = 100),
        rep(90, 100),
        seq(89, -90, len = 100)
    )))
  ) |>
  st_sfc(crs = "WGS84") |>
  st_as_sf()

sf::st_transform(ocean, crs = target_crs)

# Convert d to an sf object
d_geo <- st_as_sf(d, coords = c("Longitude", "Latitude"), crs = 4326)

pmf_map <- ggplot(data = worldMap) +
  geom_sf(data = ocean, fill = "#d0e6f7") +
  geom_sf(color = NA, fill = "gray70") +  # Light gray background for the world map
  geom_sf(data = d_geo, aes(fill = plosive_markedness_fullness, size = plosive_markedness_fullness), 
          shape = 21, color = "black", stroke = 0.2, alpha = 0.8) +  # Add black borders to circles and use fill colors
  scale_fill_viridis_c(
    trans = "sqrt",
    option = "inferno",
    direction = -1,  # Reversed inferno for high values in dark colors
    name = NULL,
    #breaks = c(1, 3, 5, 7, 10, 12),  # Key intervals for better legend readability
    guide = guide_colorbar(barwidth = 1, barheight = 15)
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +  # Variable point size based on value
  theme_minimal() +  # Light background theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  ggtitle("Plosive Markedness Fullness") +  # Add the main title
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5) 
  )
pmf_map

ggsave("../figures/plosives_map.svg", plot = pmf_map, width = 18, height = 12, device = "svg")

fmf_map <- ggplot(data = worldMap) +
  geom_sf(data = ocean, fill = "#d0e6f7") +
  geom_sf(color = NA, fill = "gray70") +  # Light gray background for the world map
  geom_sf(data = d_geo, aes(fill = fricative_markedness_fullness, size = fricative_markedness_fullness), 
          shape = 21, color = "black", stroke = 0.2, alpha = 0.8) +  # Add black borders to circles and use fill colors
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,  # Reversed inferno for high values in dark colors
    name = NULL,
    #breaks = c(1, 3, 5, 7, 10, 12),  # Key intervals for better legend readability
    guide = guide_colorbar(barwidth = 1, barheight = 15)
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +  # Variable point size based on value
  theme_minimal() +  # Light background theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  ggtitle("Fricative Markedness Fullness") +  # Add the main title
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5) 
  )
fmf_map

ggsave("../figures/fricatives_map.svg", plot = fmf_map, width = 18, height = 12, device = "svg")

amf_map <- ggplot(data = worldMap) +
  geom_sf(data = ocean, fill = "#d0e6f7") +
  geom_sf(color = NA, fill = "gray70") +  # Light gray background for the world map
  geom_sf(data = d_geo, aes(fill = affricate_markedness_fullness, size = affricate_markedness_fullness), 
          shape = 21, color = "black", stroke = 0.2, alpha = 0.8) +  # Add black borders to circles and use fill colors
  scale_fill_viridis_c(
    trans = "sqrt",
    option = "inferno",
    direction = -1,  # Reversed inferno for high values in dark colors
    name = NULL,
    #breaks = c(1, 3, 5, 7, 10, 12),  # Key intervals for better legend readability
    guide = guide_colorbar(barwidth = 1, barheight = 15)
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +  # Variable point size based on value
  theme_minimal() +  # Light background theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  ggtitle("Affricate Markedness Fullness") +  # Add the main title
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5) 
  )
amf_map

ggsave("../figures/affricates_map.svg", plot = amf_map, width = 18, height = 12, device = "svg")

library(ggplot2)

smf_map <- ggplot(data = worldMap) +
  geom_sf(data = ocean, fill = "#d0e6f7") +
  geom_sf(color = NA, fill = "gray70") +  # Light gray background for the world map
  geom_sf(data = d_geo, aes(fill = series_markedness_fullness, size = series_markedness_fullness), 
          shape = 21, color = "black", stroke = 0.2, alpha = 0.8) +  # Data points
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,
    name = NULL
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +  # Adjust point sizes
  theme_minimal() +  # Light background theme
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  ggtitle("Series Markedness Fullness") +  # Add the main title
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5) 
  )

smf_map

ggsave("../figures/series_map.svg", plot = smf_map, width = 18, height = 12, device = "svg")

