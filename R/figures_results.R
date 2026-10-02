## Benchmark results figures.
##
## Sample size is a discrete design factor (30, 50, 100, 200, 500), so no
## figure draws a line between sample sizes. Each figure is built around the
## comparison it is meant to communicate:
##   Fig 1  held out AUC: dot and 95% interval per model at each n (dodged)
##   Fig 2  paired contrasts as forest plots: rows are n, whiskers are 95% CI,
##          filled markers where the interval excludes zero, marker size is the
##          share of replicates with a positive difference
##   Fig 3  optimism as dumbbells: random and spatial block estimates for the
##          same cell joined by a segment, so the (small) effect of blocking is
##          the length of the segment
##   Fig 4  coverage effect as a forest plot (dispersed minus clustered), and
##          failure rates as a labelled tile matrix
##   Fig 5  calibration slope: median and interquartile range per model at
##          each n, against the reference slope of 1
## Inputs are the tables written by R/analysis.R (results/analysis_out/T_*.csv)
## and results/results_all.rds (for calibration quantiles). Fig 1 is Figure S5 and
## Fig 2b is Figure S6 of the Supporting Information.
## The five model colours are the same as in every other results figure.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

model_lv  <- c("rf", "rf_spatial", "gam", "gam_spatial", "gmrf")
model_cols <- c(rf = "#E69F00", rf_spatial = "#D55E00", gam = "#56B4E9", gam_spatial = "#0072B2", gmrf = "#009E73")
model_labs <- c(rf = "RF", rf_spatial = "RF + coordinates", gam = "GAM", gam_spatial = "GAM + s(X,Y)", gmrf = "GMRF (sdmTMB)")
n_lv <- c(30, 50, 100, 200, 500)
facet_lab <- labeller(coverage = c(dispersed = "Dispersed subsample", clustered = "Clustered subsample"),
                      target = c(conditional = "Conditional target", marginal = "Marginal target"))
fx <- function(d) {
  if ("model" %in% names(d))    d$model    <- factor(d$model, levels = model_lv)
  if ("n" %in% names(d))        d$n_f      <- factor(d$n, levels = n_lv)
  if ("coverage" %in% names(d)) d$coverage <- factor(d$coverage, levels = c("dispersed", "clustered"))
  if ("target" %in% names(d))   d$target   <- factor(d$target, levels = c("conditional", "marginal"))
  d
}
base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey70"), legend.position = "bottom",
        strip.text = element_text(size = 10))
pd <- position_dodge(width = 0.65)

## ---- Fig 1: held out AUC, dot and interval at each n --------------------------
cell <- read.csv("results/analysis_out/T_cell_summaries.csv") |> fx() |>
  mutate(se = auc_sd / sqrt(reps), lo = auc - qt(0.975, pmax(reps - 1, 1)) * se, hi = auc + qt(0.975, pmax(reps - 1, 1)) * se)
f1 <- ggplot(cell, aes(n_f, auc, colour = model)) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 2.2) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  labs(x = "Training records", y = "Held out AUC\n(mean and 95% interval over replicates)") +
  base_theme
ggsave("results/analysis_out/Fig1_v2_auc_dots.png", f1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig1_v2_auc_dots.pdf", f1, width = 8, height = 6, bg = "white")

## ---- Fig 2: paired contrasts as forest plots -----------------------------------
ctr <- read.csv("results/analysis_out/T_paired_contrasts_auc.csv") |>
  mutate(key = case_when(grepl("^rf_spatial - rf", contrast) ~ "rf_spatial",
                         grepl("^gam_spatial - gam", contrast) ~ "gam_spatial",
                         grepl("^gmrf - gam_spatial", contrast) ~ "gmrf",
                         grepl("^rf - gam", contrast) ~ "rf",
                         grepl("^rf_spatial - gam_spatial", contrast) ~ "rf_spatial_vs_gam_spatial"),
         group = ifelse(key %in% c("rf_spatial", "gam_spatial", "gmrf"), "Adding spatial structure", "Algorithm family"),
         excl0 = lo > 0 | hi < 0,
         key = factor(key, levels = c("rf_spatial", "gam_spatial", "gmrf", "rf", "rf_spatial_vs_gam_spatial")),
         n_f = factor(n, levels = rev(n_lv)),
         coverage = factor(coverage, levels = c("dispersed", "clustered")),
         target = factor(target, levels = c("conditional", "marginal")))
ctr_labs <- c(rf_spatial = "RF + coordinates minus RF", gam_spatial = "GAM + s(X,Y) minus GAM",
              gmrf = "GMRF minus GAM + s(X,Y)", rf = "RF minus GAM (non spatial)",
              rf_spatial_vs_gam_spatial = "RF + coordinates minus GAM + s(X,Y)")
ctr_cols <- c(rf_spatial = "#D55E00", gam_spatial = "#0072B2", gmrf = "#009E73", rf = "#E69F00",
              rf_spatial_vs_gam_spatial = "#CC79A7")
forest <- function(d, cols, labs, xlab) {
  ggplot(d, aes(mean, n_f, colour = key)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), linewidth = 0.6) +
    geom_point(aes(shape = excl0, size = sign_pos), position = position_dodge(width = 0.7), fill = "white", stroke = 0.9) +
    facet_grid(target ~ coverage, labeller = facet_lab) +
    scale_colour_manual(values = cols, labels = labs, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "excludes zero", `FALSE` = "includes zero"),
                       name = "95% interval") +
    scale_size_continuous(range = c(1.2, 3.6), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                          name = "Share of replicates with positive difference") +
    labs(x = xlab, y = "Training records") +
    base_theme + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey90"),
                       legend.box = "vertical", legend.spacing.y = unit(0, "pt")) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 3)),
           shape = guide_legend(order = 2, override.aes = list(size = 3)),
           size = guide_legend(order = 3))
}
f2a <- forest(ctr |> filter(group == "Adding spatial structure"), ctr_cols, ctr_labs,
              "Paired difference in held out AUC (adding spatial structure within a family)")
f2b <- forest(ctr |> filter(group == "Algorithm family"), ctr_cols, ctr_labs,
              "Paired difference in held out AUC (random forest minus GAM)")
ggsave("results/analysis_out/Fig2_v2_spatial_forest.png", f2a, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig2_v2_spatial_forest.pdf", f2a, width = 8, height = 7, bg = "white")
ggsave("results/analysis_out/Fig2b_v2_algorithm_forest.png", f2b, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig2b_v2_algorithm_forest.pdf", f2b, width = 8, height = 7, bg = "white")

## ---- Fig 3: optimism as dumbbells (random vs spatial block) --------------------
opt <- read.csv("results/analysis_out/T_optimism_headline.csv") |>
  mutate(lo = mean - 1.96 * se, hi = mean + 1.96 * se,
         cv = factor(cv, levels = c("random", "spatial"), labels = c("Random 5 fold", "Spatial block 5 fold")),
         n_f = factor(n, levels = rev(n_lv)),
         coverage = factor(coverage, levels = c("dispersed", "clustered")),
         target = factor(target, levels = c("conditional", "marginal")))
seg <- opt |> select(n_f, coverage, target, cv, mean) |> pivot_wider(names_from = cv, values_from = mean)
f3 <- ggplot(opt, aes(y = n_f)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(data = seg, aes(x = `Random 5 fold`, xend = `Spatial block 5 fold`, yend = n_f), colour = "grey45", linewidth = 0.8) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = cv), position = position_dodge(width = 0.5), linewidth = 0.45, alpha = 0.55) +
  geom_point(aes(x = mean, colour = cv, shape = cv), size = 2.8, fill = "white", stroke = 1) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_colour_manual(values = c("Random 5 fold" = "grey15", "Spatial block 5 fold" = "grey15"), name = NULL) +
  scale_shape_manual(values = c("Random 5 fold" = 16, "Spatial block 5 fold" = 24), name = NULL) +
  annotate("text", x = -0.13, y = 5.45, label = "pessimistic", size = 3, colour = "grey30", hjust = 0) +
  annotate("text", x = 0.01, y = 5.45, label = "optimistic", size = 3, colour = "grey30", hjust = 0) +
  labs(x = "Optimism: cross validation AUC minus held out AUC (models pooled, 95% interval)", y = "Training records") +
  base_theme + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey90"))
ggsave("results/analysis_out/Fig3_v2_optimism_dumbbell.png", f3, width = 8, height = 6.2, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig3_v2_optimism_dumbbell.pdf", f3, width = 8, height = 6.2, bg = "white")

## ---- Fig 4: coverage effect forest and failure tile matrix ---------------------
cov <- read.csv("results/analysis_out/T_coverage_effect_auc.csv") |> fx() |>
  mutate(n_f = factor(n, levels = rev(n_lv)), excl0 = lo > 0 | hi < 0)
f4a <- ggplot(cov, aes(mean, n_f, colour = model)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), linewidth = 0.6) +
  geom_point(aes(shape = excl0), position = position_dodge(width = 0.7), size = 2.4, fill = "white", stroke = 0.9) +
  facet_wrap(~ target, labeller = facet_lab) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "excludes zero", `FALSE` = "includes zero"),
                     name = "95% interval") +
  labs(x = "Dispersed minus clustered, paired difference in held out AUC", y = "Training records") +
  base_theme + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey90"),
                     legend.box = "vertical", legend.spacing.y = unit(0, "pt")) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)), shape = guide_legend(order = 2, override.aes = list(size = 3)))
fail <- read.csv("results/analysis_out/T_failure_rates.csv") |> fx()
f4b <- ggplot(fail, aes(n_f, model, fill = unusable_rate)) +
  geom_tile(colour = "grey80", linewidth = 0.5) +
  geom_text(aes(label = ifelse(unusable_rate == 0, "0", sprintf("%.0f%%", 100 * unusable_rate)),
                colour = unusable_rate > 0.25), size = 3.2) +
  facet_wrap(~ coverage, labeller = facet_lab) +
  scale_fill_gradient(low = "white", high = "#4a4a4a", limits = c(0, 0.45), labels = scales::percent, name = "Unusable fits") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  scale_y_discrete(limits = rev(model_lv), labels = model_labs) +
  labs(x = "Training records", y = NULL) +
  base_theme + theme(panel.grid = element_blank(), legend.position = "right")
f4 <- f4a / f4b + plot_layout(heights = c(1.4, 1))
ggsave("results/analysis_out/Fig4_v2_coverage_failures.png", f4, width = 8, height = 8.2, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig4_v2_coverage_failures.pdf", f4, width = 8, height = 8.2, bg = "white")

## ---- Fig 5: calibration slope, median and IQR ---------------------------------
res <- readRDS("results/results_all.rds") |> as_tibble() |>
  filter(is.na(error), !(degenerate %in% TRUE), is.finite(truth_auc), is.finite(truth_cal)) |> fx()
cal <- res |> group_by(model, n_f, coverage, target) |>
  summarise(med = median(truth_cal), q1 = quantile(truth_cal, 0.25), q3 = quantile(truth_cal, 0.75), .groups = "drop")
f5 <- ggplot(cal, aes(n_f, med, colour = model)) +
  geom_hline(yintercept = 1, colour = "grey50") +
  geom_linerange(aes(ymin = q1, ymax = q3), position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 2.2) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  coord_cartesian(ylim = c(0, 1.6)) +
  labs(x = "Training records", y = "Calibration slope on the held out target\n(median and interquartile range)") +
  base_theme
ggsave("results/analysis_out/Fig5_v2_calibration.png", f5, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/analysis_out/Fig5_v2_calibration.pdf", f5, width = 8, height = 6, bg = "white")
cat("Wrote Fig1_v2 to Fig5_v2 in results/analysis_out/\n")
