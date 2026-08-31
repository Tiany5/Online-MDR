# ==============================================================================
# Combine Settings 1-3 into a single 3x2 figure (rows = settings, columns =
# MDR / FDR). This reads the cached numerical results (Setting{1,2,3}_results.rds)
# produced by the individual scripts, so NO experiment is rerun. Labels use the
# SCT (Static calibration threshold) and RTK (rolling top-k) abbreviations.
#
# After swapping: Setting 1 = Skewed (Exp/Gamma), Setting 2 = Gaussian shift,
# Setting 3 = High-dimensional.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

code_dir <- file.path(Sys.getenv("OMDRC_ROOT", unset = getwd()), "code-semi")
fig_dir <- file.path(code_dir, "setting1,2,3")

m <- seq(from = 100, to = 1000, by = 50)
alpha <- 0.1

out1 <- readRDS(file.path(fig_dir, "Setting1_results.rds"))
out2 <- readRDS(file.path(fig_dir, "Setting2_results.rds"))
out3 <- readRDS(file.path(fig_dir, "Setting3_results.rds"))

my_colors <- c(
  "OMDRC.OR" = "#F8766D",
  "OMDRC.DD" = "#7CAE00",
  "SCT" = "#619CFF",
  "RTK" = "#FF61C3",
  "Adj-SAFFRON" = "#C77CFF"
)
my_shapes <- c(
  "OMDRC.OR" = 16, "OMDRC.DD" = 17,
  "SCT" = 18, "RTK" = 8, "Adj-SAFFRON" = 3
)

custom_theme <- theme_bw() +
  theme(
    plot.subtitle = element_text(size = 22, hjust = 0.5, margin = margin(b = 5)),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.key.size = unit(1.4, "lines"),
    axis.text = element_text(size = 19, colour = "black"),
    axis.title = element_text(size = 22),
    panel.grid.major = element_line(colour = "grey92", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)
  )

method_suffixes <- c("or", "dd", "static", "topk", "lord")
method_labels <- c("OMDRC.OR", "OMDRC.DD", "SCT", "RTK", "Adj-SAFFRON")

prepare_tidy_data <- function(res_list, metric_prefix, time_points,
                              suffixes, labels) {
  pieces <- lapply(seq_along(suffixes), function(j) {
    mat <- res_list[[paste0(metric_prefix, ".", suffixes[j])]]
    estimate <- colMeans(mat, na.rm = TRUE)
    se <- apply(mat, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(mat))
    data.frame(
      value = estimate,
      lower = pmax(0, estimate - 1.96 * se),
      upper = pmin(1, estimate + 1.96 * se),
      type = labels[j],
      t = time_points
    )
  })
  df <- bind_rows(pieces)
  df$type <- factor(df$type, levels = labels)
  df
}

build_mdr <- function(out, sub) {
  df <- prepare_tidy_data(out, "mdr", m, method_suffixes, method_labels)
  ggplot(df, aes(x = t, y = value, color = type, shape = type, group = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(size = 1.0) +
    geom_point(size = 2.4) +
    geom_hline(yintercept = alpha, linetype = "dashed", size = 1.0,
               colour = "black") +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_shape_manual(values = my_shapes) +
    labs(subtitle = sub, x = "Time (t)", y = "MDR") +
    custom_theme +
    coord_cartesian(ylim = c(0, NA)) +
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
}

build_fdr <- function(out, sub) {
  df <- prepare_tidy_data(out, "fdr", m, method_suffixes, method_labels)
  ggplot(df, aes(x = t, y = value, color = type, shape = type, group = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(size = 1.0) +
    geom_point(size = 2.4) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_shape_manual(values = my_shapes) +
    labs(subtitle = sub, x = "Time (t)", y = "FDR") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    custom_theme
}

g_3x2 <- (
  build_mdr(out1, "(a.1)") + build_fdr(out1, "(a.2)")
) / (
  build_mdr(out2, "(b.1)") + build_fdr(out2, "(b.2)")
) / (
  build_mdr(out3, "(c.1)") + build_fdr(out3, "(c.2)")
) +
  plot_layout(guides = "collect") &
  # The five-entry legend must fit inside the 10in figure width, otherwise the
  # last label (Adj-SAFFRON) and the right column are clipped: keep the keys and
  # the inter-entry spacing compact. The extra right plot margin keeps the
  # overhanging last x tick label ("1000") from being cut as well.
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "lines"),
        legend.spacing.x = unit(0.2, "lines"),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
        plot.margin = margin(t = 5.5, r = 14, b = 5.5, l = 5.5))

if (interactive()) print(g_3x2)
ggsave(file.path(fig_dir, "Setting123_3x2.pdf"), plot = g_3x2,
       width = 10, height = 13)
ggsave(file.path(fig_dir, "Setting123_3x2.png"), plot = g_3x2,
       width = 10, height = 13, dpi = 300)
cat("Saved 3x2 figure to", file.path(fig_dir, "Setting123_3x2.{pdf,png}"), "\n")
