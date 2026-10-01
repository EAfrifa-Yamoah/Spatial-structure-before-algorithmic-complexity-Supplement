## Route 2 benchmark: analysis of results_all.rds
##
## Inputs : results_all.rds (one row per species x target x n x coverage x rep x model)
## Outputs: tables (CSV) and figures (PNG, PDF) in ./analysis_out/
##
## Analytical choices
##  * Primary estimand for each contrast is the paired difference within a
##    replicate (same species, rep, n, coverage, target, hence the same training
##    subsample and the same evaluation set). Pairing removes the between
##    subsample variance that dominates at small n. A cell mean of paired
##    differences is reported with a t based 95% CI over replicates (species
##    pooled) and with the proportion of replicates in which the difference is
##    positive ("sign consistency"), which is the robustness statistic the
##    Methods promised.
##  * A mixed model (lmer) with species/rep random intercepts is fitted as a
##    robustness check on the same conclusions. With three species the species
##    variance component is weakly identified; treat its output as secondary.
##  * Optimism = cross validation estimate minus the held out truth on the
##    matched target. Positive values overstate transferable skill.
##  * Runs with an error, non convergence (gmrf) or a degenerate fit
##    (gam_spatial) are excluded from performance summaries and reported
##    separately as failure rates. Excluding failures makes the performance
##    summaries conditional on a usable fit, which is what a practitioner sees.
##
## Packages: dplyr, tidyr, ggplot2, patchwork, lme4, emmeans

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
dir.create("results/analysis_out", showWarnings = FALSE, recursive = TRUE)
res_file <- if (file.exists("results/results_all.rds")) "results/results_all.rds" else stop("results/results_all.rds not found")
res <- readRDS(res_file) |> as_tibble()

## ---- 0. Bookkeeping ---------------------------------------------------------
res <- res |>
  mutate(model    = factor(model, levels = c("rf", "rf_spatial", "gam", "gam_spatial", "gmrf")),
         n        = as.integer(n),
         coverage = factor(coverage, levels = c("dispersed", "clustered")),
         target   = factor(target, levels = c("conditional", "marginal")),
         ## unusable: an error, no finite prediction, or a degenerate smooth surface.
         ## The sdmTMB sanity flag (converged) is kept as a separate diagnostic and
         ## examined in a sensitivity analysis rather than used for exclusion, because
         ## the scaffold stores only the combined hessian/nlminb/gradient flag and
         ## the gradient check alone fails many otherwise usable fits.
         failed   = !is.na(error) | (degenerate %in% TRUE) | !is.finite(truth_auc),
         opt_random  = cvr_auc - truth_auc,
         opt_spatial = cvs_auc - truth_auc)

n_runs <- res |> distinct(species, target, n, coverage, rep) |> nrow()
cat(sprintf("Replicate runs: %d (%d species, %d reps per cell); model rows: %d\n",
            n_runs, n_distinct(res$species), max(res$rep), nrow(res)))
cat("Species:", paste(unique(res$species), collapse = "; "), "\n")

## ---- 1. Failure rates -------------------------------------------------------
fail <- res |>
  group_by(model, n, coverage) |>
  summarise(runs = n(),
            error_rate      = mean(!is.na(error)),
            nonconv_rate    = mean(!converged, na.rm = TRUE),
            degenerate_rate = mean(degenerate %in% TRUE),
            unusable_rate   = mean(failed),
            .groups = "drop")
write.csv(fail, "results/analysis_out/T_failure_rates.csv", row.names = FALSE)

ok <- res |> filter(!failed)

## ---- 2. Cell summaries of truth metrics and optimism ------------------------
cell <- ok |>
  group_by(model, n, coverage, target) |>
  summarise(reps = n(),
            auc = mean(truth_auc), auc_sd = sd(truth_auc),
            brier = mean(truth_brier), cal = median(truth_cal, na.rm = TRUE),
            opt_r = mean(opt_random, na.rm = TRUE), opt_r_sd = sd(opt_random, na.rm = TRUE),
            opt_s = mean(opt_spatial, na.rm = TRUE), opt_s_sd = sd(opt_spatial, na.rm = TRUE),
            .groups = "drop")
write.csv(cell, "results/analysis_out/T_cell_summaries.csv", row.names = FALSE)

## Main text Table 4: mean held out AUC over all usable fits by model and n, with SD and count
t4 <- ok |> group_by(model, n) |>
  summarise(k = n(), mean = mean(truth_auc), sd = sd(truth_auc), se = sd / sqrt(k), .groups = "drop")
write.csv(t4, "results/analysis_out/T_table4_auc_by_model_n.csv", row.names = FALSE)
## Supplementary Table S12: AUC by species under dispersed conditional coverage
t_sp <- ok |> filter(coverage == "dispersed", target == "conditional") |>
  group_by(species, model, n) |> summarise(auc = mean(truth_auc), .groups = "drop") |>
  pivot_wider(names_from = n, values_from = auc)
write.csv(t_sp, "results/analysis_out/T_species_auc_dispersed_conditional.csv", row.names = FALSE)

## ---- 3. Paired contrasts ----------------------------------------------------
wide <- ok |>
  select(species, rep, n, coverage, target, model, truth_auc) |>
  pivot_wider(names_from = model, values_from = truth_auc)

paired_summary <- function(d, diff, label, by = c("n", "coverage", "target")) {
  d |>
    mutate(diff = {{ diff }}) |>
    filter(is.finite(diff)) |>
    group_by(across(all_of(by))) |>
    summarise(contrast = label, reps = n(),
              mean = mean(diff), se = sd(diff) / sqrt(n()),
              lo = if (n() > 1) mean - qt(0.975, n() - 1) * se else NA_real_,
              hi = if (n() > 1) mean + qt(0.975, n() - 1) * se else NA_real_,
              sign_pos = mean(diff > 0),
              .groups = "drop")
}
contrasts <- bind_rows(
  paired_summary(wide, rf_spatial - rf,          "rf_spatial - rf (spatial effect, RF)"),
  paired_summary(wide, gam_spatial - gam,        "gam_spatial - gam (spatial effect, GAM)"),
  paired_summary(wide, gmrf - gam_spatial,       "gmrf - gam_spatial (GMRF vs smooth surface)"),
  paired_summary(wide, rf - gam,                 "rf - gam (algorithm, non spatial)"),
  paired_summary(wide, rf_spatial - gam_spatial, "rf_spatial - gam_spatial (algorithm, spatial)")
)
write.csv(contrasts, "results/analysis_out/T_paired_contrasts_auc.csv", row.names = FALSE)

## Coverage effect: dispersed minus clustered, paired on species, rep, n, target, model
cov_wide <- ok |>
  select(species, rep, n, coverage, target, model, truth_auc) |>
  pivot_wider(names_from = coverage, values_from = truth_auc)
coverage_effect <- paired_summary(cov_wide, dispersed - clustered,
                                  "dispersed - clustered", by = c("n", "target", "model"))
write.csv(coverage_effect, "results/analysis_out/T_coverage_effect_auc.csv", row.names = FALSE)

## Optimism summaries with CI (species pooled)
opt <- ok |>
  pivot_longer(c(opt_random, opt_spatial), names_to = "cv", values_to = "optimism") |>
  mutate(cv = recode(cv, opt_random = "random 5 fold", opt_spatial = "spatial block 5 fold")) |>
  filter(is.finite(optimism)) |>
  group_by(model, n, coverage, target, cv) |>
  summarise(reps = n(), mean = mean(optimism), se = sd(optimism) / sqrt(n()),
            .groups = "drop") |>
  mutate(lo = ifelse(reps > 1, mean - qt(0.975, pmax(reps - 1, 1)) * se, NA_real_),
         hi = ifelse(reps > 1, mean + qt(0.975, pmax(reps - 1, 1)) * se, NA_real_))
write.csv(opt, "results/analysis_out/T_optimism.csv", row.names = FALSE)

## Optimism averaged over models, by n x coverage x target x cv (headline table)
opt_head <- ok |>
  pivot_longer(c(opt_random, opt_spatial), names_to = "cv", values_to = "optimism") |>
  filter(is.finite(optimism)) |>
  group_by(n, coverage, target, cv) |>
  summarise(reps = n(), mean = mean(optimism), se = sd(optimism) / sqrt(n()), .groups = "drop") |>
  mutate(cv = recode(cv, opt_random = "random", opt_spatial = "spatial"))
write.csv(opt_head, "results/analysis_out/T_optimism_headline.csv", row.names = FALSE)

## ---- 4. Mixed model robustness check ----------------------------------------
## Sensitivity: gmrf held out AUC with and without the strict sanity flag
gmrf_sens <- ok |> filter(model == "gmrf") |>
  group_by(n, coverage, converged) |>
  summarise(reps = n(), auc = mean(truth_auc), .groups = "drop") |>
  pivot_wider(names_from = converged, values_from = c(reps, auc), names_prefix = "sanity_")
write.csv(gmrf_sens, "results/analysis_out/T_gmrf_sanity_sensitivity.csv", row.names = FALSE)
strict <- res |> group_by(model, n, coverage) |> summarise(sanity_pass = mean(converged, na.rm = TRUE), .groups = "drop")
write.csv(strict, "results/analysis_out/T_sanity_pass_rates.csv", row.names = FALSE)

if (n_distinct(res$species) > 1 && requireNamespace("lme4", quietly = TRUE) && requireNamespace("emmeans", quietly = TRUE)) {
  library(lme4); library(emmeans)
  ok_mm <- ok |> mutate(n_f = factor(n), sp_rep = interaction(species, rep))
  m_auc <- lmer(truth_auc ~ model * n_f * coverage * target + (1 | species) + (1 | sp_rep),
                data = ok_mm, REML = TRUE)
  em <- emmeans(m_auc, ~ model | n_f + coverage + target)
  mm_spatial <- contrast(em, method = list("rf_spatial - rf" = c(-1, 1, 0, 0, 0),
                                           "gam_spatial - gam" = c(0, 0, -1, 1, 0),
                                           "gmrf - gam_spatial" = c(0, 0, 0, -1, 1),
                                           "rf - gam" = c(1, 0, -1, 0, 0)),
                         adjust = "none") |> as.data.frame()
  write.csv(mm_spatial, "results/analysis_out/T_mixed_model_contrasts.csv", row.names = FALSE)
  vc <- as.data.frame(VarCorr(m_auc)); write.csv(vc, "results/analysis_out/T_mixed_model_varcomp.csv", row.names = FALSE)
  m_opt <- lmer(opt_random ~ model * n_f * coverage * target + (1 | species) + (1 | sp_rep),
                data = ok_mm |> filter(is.finite(opt_random)))
  em_opt <- emmeans(m_opt, ~ coverage * target | n_f) |> as.data.frame()
  write.csv(em_opt, "results/analysis_out/T_mixed_model_optimism_by_cell.csv", row.names = FALSE)
}

## ---- 5. Figures -------------------------------------------------------------
## One categorical palette for the five models, used in every figure (Okabe Ito).
model_cols <- c(rf = "#E69F00", rf_spatial = "#D55E00", gam = "#56B4E9",
                gam_spatial = "#0072B2", gmrf = "#009E73")
model_labs <- c(rf = "RF", rf_spatial = "RF + coordinates", gam = "GAM",
                gam_spatial = "GAM + s(X,Y)", gmrf = "GMRF (sdmTMB)")
base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(size = 10), panel.border = element_rect(fill = NA, colour = "grey70"))
facet_lab <- labeller(coverage = c(dispersed = "Dispersed subsample", clustered = "Clustered subsample"),
                      target = c(conditional = "Conditional target", marginal = "Marginal target"))

## Figure A: held out AUC against n
fa <- ggplot(cell, aes(n, auc, colour = model, group = model)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  labs(x = "Training records (log scale)", y = "Held out AUC (mean over replicates)") +
  base_theme
ggsave("results/analysis_out/FigA_auc_vs_n.png", fa, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigA_auc_vs_n.pdf", fa, width = 8, height = 6, bg = "white")

## Figure B: spatial effect within family (paired), with CI and sign consistency
cb <- contrasts |>
  filter(contrast %in% c("rf_spatial - rf (spatial effect, RF)", "gam_spatial - gam (spatial effect, GAM)",
                         "gmrf - gam_spatial (GMRF vs smooth surface)")) |>
  mutate(family = recode(contrast,
                         "rf_spatial - rf (spatial effect, RF)" = "rf_spatial",
                         "gam_spatial - gam (spatial effect, GAM)" = "gam_spatial",
                         "gmrf - gam_spatial (GMRF vs smooth surface)" = "gmrf"),
         family = factor(family, levels = names(model_cols)))
fb <- ggplot(cb, aes(n, mean, colour = family, group = family)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = family), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(aes(size = sign_pos)) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_colour_manual(values = model_cols, name = NULL,
                      labels = c(rf_spatial = "RF + coordinates minus RF",
                                 gam_spatial = "GAM + s(X,Y) minus GAM",
                                 gmrf = "GMRF minus GAM + s(X,Y)")) +
  scale_fill_manual(values = model_cols, guide = "none") +
  scale_size_continuous(range = c(1, 4), limits = c(0, 1), name = "Share of replicates with positive difference") +
  labs(x = "Training records (log scale)", y = "Paired difference in held out AUC (95% CI)") +
  base_theme + theme(legend.box = "vertical")
ggsave("results/analysis_out/FigB_spatial_effect.png", fb, width = 8, height = 6.5, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigB_spatial_effect.pdf", fb, width = 8, height = 6.5, bg = "white")

## Figure C: optimism of random and spatial cross validation
fc <- ggplot(opt, aes(n, mean, colour = model, linetype = cv, shape = cv, group = interaction(model, cv))) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  scale_linetype_manual(values = c("random 5 fold" = "solid", "spatial block 5 fold" = "22"), name = "Cross validation") +
  scale_shape_manual(values = c("random 5 fold" = 16, "spatial block 5 fold" = 17), name = "Cross validation") +
  labs(x = "Training records (log scale)", y = "Optimism: cross validation AUC minus held out AUC") +
  base_theme + theme(legend.box = "vertical")
ggsave("results/analysis_out/FigC_optimism.png", fc, width = 8, height = 6.5, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigC_optimism.pdf", fc, width = 8, height = 6.5, bg = "white")


## Figure C2: optimism pooled over models, with 95% CI (headline view)
opt_pool <- ok |>
  pivot_longer(c(opt_random, opt_spatial), names_to = "cv", values_to = "optimism") |>
  mutate(cv = recode(cv, opt_random = "random 5 fold", opt_spatial = "spatial block 5 fold")) |>
  filter(is.finite(optimism)) |>
  group_by(n, coverage, target, cv) |>
  summarise(reps = n(), mean = mean(optimism), se = sd(optimism) / sqrt(n()),
            lo = mean - qt(0.975, n() - 1) * se, hi = mean + qt(0.975, n() - 1) * se, .groups = "drop")
fc2 <- ggplot(opt_pool, aes(n, mean, linetype = cv, shape = cv, group = cv)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = cv), fill = "grey70", alpha = 0.3) +
  geom_line(linewidth = 0.7, colour = "grey20") + geom_point(size = 2.2, colour = "grey20") +
  facet_grid(target ~ coverage, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_linetype_manual(values = c("random 5 fold" = "solid", "spatial block 5 fold" = "22"), name = "Cross validation") +
  scale_shape_manual(values = c("random 5 fold" = 16, "spatial block 5 fold" = 17), name = "Cross validation") +
  labs(x = "Training records (log scale)", y = "Optimism, models pooled: cross validation AUC minus held out AUC (95% CI)") +
  base_theme
ggsave("results/analysis_out/FigC2_optimism_pooled.png", fc2, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigC2_optimism_pooled.pdf", fc2, width = 8, height = 6, bg = "white")

## Figure D: coverage effect (dispersed minus clustered) and failure rates
fd1 <- ggplot(coverage_effect, aes(n, mean, colour = model, group = model)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~ target, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  scale_fill_manual(values = model_cols, guide = "none") +
  labs(x = NULL, y = "Dispersed minus clustered,\nheld out AUC (95% CI)") +
  base_theme
fd2 <- ggplot(fail, aes(n, unusable_rate, colour = model, group = model)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~ coverage, labeller = facet_lab) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 500)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = model_cols, labels = model_labs, name = NULL) +
  labs(x = "Training records (log scale)", y = "Unusable fits (error, no prediction\nor degenerate surface)") +
  base_theme
fd <- fd1 / fd2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", axis.title.y = element_text(size = 10, lineheight = 0.9, margin = margin(r = 6)))
ggsave("results/analysis_out/FigD_coverage_and_failures.png", fd, width = 8, height = 8, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigD_coverage_and_failures.pdf", fd, width = 8, height = 8, bg = "white")

## ---- 6. Console digest ------------------------------------------------------
options(width = 160)
cat("\n== Failure rates (unusable fits) by model x n, averaged over coverage ==\n")
print(fail |> group_by(model, n) |> summarise(unusable = mean(unusable_rate), .groups = "drop") |>
        pivot_wider(names_from = n, values_from = unusable) |> as.data.frame(), digits = 2)
cat("\n== gmrf sanity flag sensitivity (held out AUC, usable fits) ==\n")
print(as.data.frame(gmrf_sens), digits = 3)
cat("\n== Held out AUC by model x n, averaged over coverage and target ==\n")
print(cell |> group_by(model, n) |> summarise(auc = mean(auc), .groups = "drop") |>
        pivot_wider(names_from = n, values_from = auc) |> as.data.frame(), digits = 3)
cat("\n== Paired contrasts (AUC), by cell ==\n")
print(as.data.frame(contrasts |> mutate(across(c(mean, lo, hi), ~ round(.x, 3)), sign_pos = round(sign_pos, 2))),
      row.names = FALSE)
cat("\n== Coverage effect (dispersed minus clustered) ==\n")
print(as.data.frame(coverage_effect |> mutate(across(c(mean, lo, hi), ~ round(.x, 3)), sign_pos = round(sign_pos, 2))),
      row.names = FALSE)
cat("\n== Optimism headline (models pooled) ==\n")
print(as.data.frame(opt_head |> mutate(across(c(mean, se), ~ round(.x, 3)))), row.names = FALSE)
cat("\nWrote tables and figures to analysis_out/\n")
