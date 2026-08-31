############################################################
## Publication-quality two-panel figure for mFDR-MDR trade-off
##
## Two-group Gaussian mixture model:
##   theta ~ Bernoulli(pi)
##   X | theta = 0 ~ N(0, 1)
##   X | theta = 1 ~ N(mu, 1)
##
## Threshold rule:
##   delta_c(X) = 1{X >= c}
############################################################

rm(list = ls())

library(ggplot2)
library(dplyr)
library(patchwork)

############################################################
## Basic functions
############################################################

mFDR_normal <- function(c, pi, mu) {
  fp <- (1 - pi) * pnorm(c, lower.tail = FALSE)
  tp <- pi * pnorm(c - mu, lower.tail = FALSE)
  fp / (fp + tp)
}

mdr_normal <- function(c, mu) {
  pnorm(c - mu)
}

solve_c_for_mFDR <- function(q, pi, mu) {
  if (q >= 1 - pi) {
    return(-Inf)
  }
  
  target <- q * pi / ((1 - q) * (1 - pi))
  
  f <- function(c) {
    pnorm(c, lower.tail = FALSE, log.p = TRUE) -
      pnorm(c - mu, lower.tail = FALSE, log.p = TRUE) -
      log(target)
  }
  
  lower <- -30
  upper <- max(30, mu + 30)
  
  while (f(upper) > 0 && upper < 300) {
    upper <- upper + 30
  }
  
  while (f(lower) < 0 && lower > -300) {
    lower <- lower - 30
  }
  
  uniroot(f, lower = lower, upper = upper)$root
}

mdr_min_under_mFDR <- function(q, pi, mu) {
  c_q <- solve_c_for_mFDR(q = q, pi = pi, mu = mu)
  
  if (is.infinite(c_q) && c_q < 0) {
    return(0)
  }
  
  mdr_normal(c = c_q, mu = mu)
}

############################################################
## Aesthetic settings
############################################################

theme_panel <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      
      axis.title.x = element_text(
        color = "black",
        size = base_size + 1,
        margin = margin(t = 3)
      ),
      axis.title.y = element_text(
        color = "black",
        size = base_size + 1,
        margin = margin(r = 3)
      ),
      
      plot.title = element_text(
        face = "bold",
        color = "black",
        size = base_size + 2,
        hjust = 0.5,
        margin = margin(t = 0, b = -40)
      ),
      
      plot.subtitle = element_blank(),
      legend.title = element_text(color = "black", size = base_size),
      legend.text = element_text(color = "black", size = base_size - 1),
      plot.margin = margin(12, 6, 2, 6)
    )
}

############################################################
## Panel A: mFDR-MDR trade-off curve with fixed pi
############################################################

pi_tradeoff <- 0.05

mu_values <- c(1, 2, 3, 4, 5)

c_grid <- seq(-8, 10, length.out = 3000)

df_tradeoff <- expand.grid(
  c = c_grid,
  mu = mu_values
) %>%
  as_tibble() %>%
  mutate(
    pi = pi_tradeoff,
    mFDR = mFDR_normal(c = c, pi = pi, mu = mu),
    MDR = mdr_normal(c = c, mu = mu),
    mu_label = factor(
      paste0("mu = ", mu),
      levels = paste0("mu = ", mu_values)
    )
  ) %>%
  filter(
    is.finite(mFDR),
    is.finite(MDR)
  )

label_targets <- tibble(
  mu = c(1, 2, 3, 4, 5),
  target_mFDR = c(0.64, 0.38, 0.15, 0.060, 0.035),
  nudge_x = c(0.000, 0.000, 0.000, 0.035, 0.050),
  nudge_y = c(0.000, 0.000, 0.000, 0.015, 0.018)
)

df_curve_labels <- df_tradeoff %>%
  inner_join(label_targets, by = "mu") %>%
  group_by(mu) %>%
  slice_min(abs(mFDR - target_mFDR), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    ## Important: use plotmath-compatible label, not Unicode mu.
    label = paste0("mu == ", mu),
    label_x = pmin(mFDR + nudge_x, 0.97),
    label_y = pmin(MDR + nudge_y, 0.98),
    mu_label = factor(
      paste0("mu = ", mu),
      levels = paste0("mu = ", mu_values)
    )
  )

p_tradeoff <- ggplot(
  df_tradeoff,
  aes(x = mFDR, y = MDR, color = mu_label)
) +
  geom_line(linewidth = 1.18, lineend = "round") +
  geom_label(
    data = df_curve_labels,
    aes(x = label_x, y = label_y, label = label),
    parse = TRUE,
    size = 3.1,
    label.size = 0.15,
    label.padding = unit(0.12, "lines"),
    fill = "white",
    alpha = 0.93,
    show.legend = FALSE
  ) +
  scale_color_viridis_d(
    option = "viridis",
    direction = 1,
    begin = 0.05,
    end = 0.90
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.015))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.015))
  ) +
  labs(
    title = bquote("mFDR-MDR trade-off (" * pi * " = " * .(pi_tradeoff) * ")"),
    x = "mFDR",
    y = "MDR"
  ) +
  theme_panel(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5,
      margin = margin(t = 0, b = -40)
    )
  )

############################################################
## Panel B: pi-mu difficulty map under an mFDR constraint
############################################################

q_level <- 0.10

pi_min_map <- 0.001
pi_max_map <- 0.20

pi_grid <- seq(pi_min_map, pi_max_map, length.out = 180)

mu_grid <- seq(1.0, 6.0, length.out = 180)

df_map <- expand.grid(
  pi = pi_grid,
  mu = mu_grid
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    MDR_min = mdr_min_under_mFDR(q = q_level, pi = pi, mu = mu),
    kappa = mu^2 / 2 - log((1 - pi) / pi)
  ) %>%
  ungroup()

df_boundary <- tibble(
  pi = pi_grid
) %>%
  mutate(
    mu = sqrt(2 * log((1 - pi) / pi))
  ) %>%
  filter(mu >= min(mu_grid), mu <= max(mu_grid))

df_kappa_label <- df_boundary %>%
  slice_min(abs(pi - 0.060), n = 1, with_ties = FALSE) %>%
  mutate(
    ## Important: use plotmath-compatible label, not Unicode kappa.
    label = "kappa == 0",
    mu_label_pos = mu + 0.30,
    pi_label_pos = pi - 0.012
  )

p_map <- ggplot(df_map, aes(x = mu, y = pi)) +
  geom_tile(aes(fill = MDR_min)) +
  geom_contour(
    aes(z = MDR_min),
    breaks = c(0.1, 0.25, 0.5, 0.75, 0.9),
    color = "white",
    linewidth = 0.35,
    alpha = 0.95
  ) +
  geom_line(
    data = df_boundary,
    aes(x = mu, y = pi),
    inherit.aes = FALSE,
    linewidth = 1,
    linetype = "longdash",
    color = "red4"
  ) +
  geom_text(
    data = df_kappa_label,
    aes(x = mu_label_pos, y = pi_label_pos, label = label),
    parse = TRUE,
    inherit.aes = FALSE,
    size = 4.2,
    angle = -62,
    color = "red4",
    fontface = "bold"
  ) +
  scale_x_continuous(
    limits = c(1.0, 6),
    breaks = seq(1, 6, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(pi_min_map, pi_max_map),
    breaks = c(0.001, 0.05, 0.10, 0.15, 0.20),
    labels = c("0.001", "0.05", "0.10", "0.15", "0.20"),
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    limits = c(0, 1),
    name = expression(MDR[mFDR]^"*"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(36, "mm"),
      barwidth = grid::unit(4.0, "mm")
    )
  ) +
  labs(
    title = sprintf("Minimum MDR under mFDR <= %.2f", q_level),
    x = expression("Signal strength " * mu),
    y = expression("Signal proportion " * pi)
  ) +
  theme_panel(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    aspect.ratio = 1,
    plot.title = element_text(
      face = "plain",
      size = 13,
      hjust = 0.5,
      margin = margin(t = 0, b = -40)
    ),
    axis.text.y = element_text(size = 10)
  )

############################################################
## Combine two panels
############################################################

p_combined <- p_tradeoff + p_map +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(2, 2, 2, 2)
    )
  )

if (interactive()) print(p_combined)