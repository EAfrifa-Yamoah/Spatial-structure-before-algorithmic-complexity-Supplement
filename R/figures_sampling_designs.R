## Maps of the Route 2 sampling designs, drawn from the same functions the
## benchmark uses (route2_benchmark_scaffold.R, sections 1 to 8), so what is
## shown is exactly what the benchmark draws.
##
## Main text
##   FigS_design_1_sampling      targets (rows) x evaluation set and example
##                               subsamples at n = 30 and 200 (columns)
##   FigS_design_2_cv_folds      random and spatial block folds, n = 200
## Supplementary (all scenarios)
##   FigS_design_S1_conditional  conditional target: coverage (rows) x n (columns)
##   FigS_design_S2_marginal     marginal target: coverage (rows) x n (columns)
##   FigS_design_S3_folds_dispersed   folds for dispersed subsamples, all n
##   FigS_design_S4_folds_clustered   folds for clustered subsamples, all n
##
## Colour meanings are fixed across all of these figures and do not reuse the
## five model colours of the results figures:
##   grey        survey hauls in the training pool that were not drawn
##   black       training subsample
##   #CC79A7     evaluation set (held out target)
## Fold identity uses viridis with distinct shapes (folds only).
##
## Packages: those of the scaffold, plus ggplot2, patchwork, sf, rnaturalearth,
## rnaturalearthdata (1:50m land polygons).

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(sf); library(dplyr) })
dir.create("results/analysis_out", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load the scaffold up to the driver (no smoke test, no full run) ----
src <- readLines("R/benchmark_scaffold.R")
cut <- grep("^## ---- 9\\. Driver", src)
eval(parse(text = src[seq_len(cut - 1)]))

FIG_SEED <- 20260925          # figure specific seed; does not touch the benchmark
set.seed(FIG_SEED)
sp <- species[1]
d <- make_response(hauls, sp)
parts <- make_partitions(d)
N_LEVELS <- c(30, 50, 100, 200, 500)

## ---- 2. Basemap and shared styling -----------------------------------------
land <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
bb <- st_bbox(c(xmin = min(hauls$lon) - 1.6, xmax = max(hauls$lon) + 0.9,
                ymin = min(hauls$lat) - 0.6, ymax = max(hauls$lat) + 0.5), crs = 4326)
land <- suppressWarnings(st_crop(st_make_valid(land), bb))

role_cols <- c(pool = "grey72", training = "black", evaluation = "#CC79A7")
role_shp  <- c(pool = 16, training = 17, evaluation = 15)
role_labs <- c(pool = "Survey hauls in the training pool", training = "Training subsample",
               evaluation = "Evaluation set (held out target)")
## degree labels through plotmath so they render in any locale
lab_lon <- function(x) parse(text = paste0(abs(x), "*degree*W"))
lab_lat <- function(x) parse(text = paste0(x, "*degree*N"))

map_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.4),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(size = 9.5),
        axis.text = element_text(size = 7.5, colour = "grey30"),
        legend.position = "bottom", panel.spacing = unit(4, "pt"))

## scale bar (200 km), exact at its latitude; drawn in one named panel only
scale_bar_layers <- function(panel_df, lat0 = 33.0, lon0 = -125.6, km = 200) {
  dlon <- km / (111.32 * cos(lat0 * pi / 180))
  sb <- cbind(panel_df, x = lon0, xend = lon0 + dlon, y = lat0, ylab = lat0 + 0.35, lab = paste(km, "km"))
  list(geom_segment(data = sb, aes(x = x, xend = xend, y = y, yend = y), linewidth = 1.2, colour = "grey20"),
       geom_text(data = sb, aes(x = (x + xend) / 2, y = ylab, label = lab), size = 2.6, colour = "grey20"))
}

base_layers <- function() list(
  geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.25),
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE),
  scale_x_continuous(breaks = c(-124, -120), labels = lab_lon),
  scale_y_continuous(breaks = c(35, 40, 45), labels = lab_lat),
  labs(x = NULL, y = NULL))

role_scales <- function() list(
  scale_colour_manual(values = role_cols, breaks = names(role_cols), labels = role_labs, name = NULL),
  scale_shape_manual(values = role_shp, breaks = names(role_cols), labels = role_labs, name = NULL),
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1)),
         shape  = guide_legend(override.aes = list(size = 2.6, alpha = 1))))

## ---- 3. Draw one example subsample per scenario ----------------------------
draw_scenario <- function(tg, cv, nn) {
  pool <- parts[[tg]]$pool
  sub <- subsample(pool, nn, cv)
  list(sub = sub, extent = attr(sub, "extent_frac"),
       band = if (cv == "clustered") range(sub$lat) else c(NA, NA))
}
scen <- expand.grid(target = c("conditional", "marginal"), coverage = c("dispersed", "clustered"),
                    n = N_LEVELS, stringsAsFactors = FALSE)
draws <- lapply(seq_len(nrow(scen)), function(i) draw_scenario(scen$target[i], scen$coverage[i], scen$n[i]))

panel_data <- function(i, with_eval = FALSE) {
  s <- scen[i, ]; dr <- draws[[i]]
  pool <- parts[[s$target]]$pool; ev <- parts[[s$target]]$eval
  out <- bind_rows(pool |> mutate(role = "pool"),
                   if (with_eval) ev |> mutate(role = "evaluation"),
                   dr$sub |> mutate(role = "training")) |>
    mutate(target = s$target, coverage = s$coverage, n = s$n)
  attr(out, "band") <- dr$band; attr(out, "extent") <- dr$extent
  out
}
nice_target <- c(conditional = "Conditional target", marginal = "Marginal target")
nice_cov <- c(dispersed = "Dispersed", clustered = "Clustered")

## ---- 4. Main figure 1: targets x (evaluation set, n = 30, n = 200) ----------
main_cells <- list()
for (tg in c("conditional", "marginal")) {
  pool <- parts[[tg]]$pool; ev <- parts[[tg]]$eval
  main_cells[[length(main_cells) + 1]] <- bind_rows(pool |> mutate(role = "pool"), ev |> mutate(role = "evaluation")) |>
    mutate(target = tg, panel = "Evaluation set", extent = NA_real_, band_lo = NA_real_, band_hi = NA_real_)
  for (cv in c("dispersed", "clustered")) for (nn in c(30, 200)) {
    i <- which(scen$target == tg & scen$coverage == cv & scen$n == nn)
    pd <- panel_data(i)
    main_cells[[length(main_cells) + 1]] <- pd |>
      mutate(panel = sprintf("%s, n = %d", nice_cov[cv], nn), extent = attr(pd, "extent"),
             band_lo = attr(pd, "band")[1], band_hi = attr(pd, "band")[2])
  }
}
dd <- bind_rows(main_cells) |>
  mutate(panel = factor(panel, levels = c("Evaluation set", "Dispersed, n = 30", "Dispersed, n = 200",
                                          "Clustered, n = 30", "Clustered, n = 200")),
         target = factor(nice_target[target], levels = nice_target),
         role = factor(role, levels = names(role_cols)))
bands <- dd |> filter(!is.na(band_lo)) |> distinct(target, panel, band_lo, band_hi)
ext_lab <- dd |> filter(!is.na(extent)) |> distinct(target, panel, extent) |>
  mutate(lab = sprintf("extent %.0f%%", 100 * extent))
if (any(ext_lab$extent < 0.05)) ext_lab$lab[ext_lab$extent < 0.05] <- sprintf("extent %.1f%%", 100 * ext_lab$extent[ext_lab$extent < 0.05])

draw_map <- function(dd, bands, ext_lab, facets, eval_size = 0.45) {
  ggplot() +
    base_layers()[[1]] +
    geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = band_lo, ymax = band_hi), fill = "grey40", alpha = 0.12) +
    geom_hline(data = bands, aes(yintercept = band_lo), colour = "grey35", linetype = "22", linewidth = 0.3) +
    geom_hline(data = bands, aes(yintercept = band_hi), colour = "grey35", linetype = "22", linewidth = 0.3) +
    geom_point(data = dd |> filter(role == "pool"), aes(lon, lat, colour = role, shape = role), size = 0.25, alpha = 0.55) +
    geom_point(data = dd |> filter(role == "evaluation"), aes(lon, lat, colour = role, shape = role), size = eval_size, alpha = 0.8) +
    geom_point(data = dd |> filter(role == "training"), aes(lon, lat, colour = role, shape = role), size = 1.15) +
    geom_text(data = ext_lab, aes(x = bb["xmin"] + 0.3, y = bb["ymin"] + 0.5, label = lab),
              hjust = 0, vjust = 0, size = 2.6, colour = "grey20") +
    facets + base_layers()[-1] + role_scales() + map_theme
}
sb_panel <- data.frame(target = factor("Conditional target", levels = nice_target),
                       panel = factor("Evaluation set", levels = levels(dd$panel)))
f1 <- draw_map(dd, bands, ext_lab, facet_grid(target ~ panel)) + scale_bar_layers(sb_panel)
ggsave("results/analysis_out/FigS_design_1_sampling.png", f1, width = 9, height = 7.4, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigS_design_1_sampling.pdf", f1, width = 9, height = 7.4, bg = "white")

## ---- 5. Supplementary S1 and S2: every scenario, one figure per target ------
for (tg in c("conditional", "marginal")) {
  cells <- lapply(which(scen$target == tg), function(i) {
    ## the conditional evaluation set is a random 30 percent interleaved with the
    ## pool and would hide it; it is shown in the main figure instead
    pd <- panel_data(i, with_eval = (tg == "marginal"))
    pd |> mutate(extent = attr(pd, "extent"), band_lo = attr(pd, "band")[1], band_hi = attr(pd, "band")[2])
  })
  ds <- bind_rows(cells) |>
    mutate(coverage = factor(nice_cov[coverage], levels = nice_cov),
           n_lab = factor(sprintf("n = %d", n), levels = sprintf("n = %d", N_LEVELS)),
           role = factor(role, levels = names(role_cols)))
  bands_s <- ds |> filter(!is.na(band_lo)) |> distinct(coverage, n_lab, band_lo, band_hi)
  ext_s <- ds |> distinct(coverage, n_lab, extent) |>
    mutate(lab = ifelse(extent < 0.05, sprintf("extent %.1f%%", 100 * extent), sprintf("extent %.0f%%", 100 * extent)))
  sb_s <- data.frame(coverage = factor("Dispersed", levels = nice_cov), n_lab = factor("n = 30", levels = levels(ds$n_lab)))
  fs <- draw_map(ds, bands_s, ext_s, facet_grid(coverage ~ n_lab), eval_size = 0.3) + scale_bar_layers(sb_s)
  fn <- if (tg == "conditional") "FigS_design_S1_conditional" else "FigS_design_S2_marginal"
  ggsave(sprintf("results/analysis_out/%s.png", fn), fs, width = 9.5, height = 7.2, dpi = 300, bg = "white")
  ggsave(sprintf("results/analysis_out/%s.pdf", fn), fs, width = 9.5, height = 7.2, bg = "white")
}

## ---- 6. Cross validation folds ----------------------------------------------
fold_panel <- function(sub, folds, label, xl, yl, xb, yb, note = NULL) {
  d2 <- sub |> mutate(fold = factor(folds), panel = label)
  p <- ggplot() +
    geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.25) +
    geom_point(data = d2, aes(lon, lat, colour = fold, shape = fold), size = 1.6) +
    facet_wrap(~ panel) +
    coord_sf(xlim = xl, ylim = yl, expand = FALSE) +
    scale_x_continuous(breaks = xb, labels = lab_lon) +
    scale_y_continuous(breaks = yb, labels = lab_lat) +
    scale_colour_viridis_d(name = "Fold", end = 0.92, drop = FALSE) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8), name = "Fold", drop = FALSE) +
    labs(x = NULL, y = NULL) + map_theme + theme(legend.position = "right")
  if (!is.null(note)) p <- p + annotate("label", x = xl[1] + 0.05 * diff(xl), y = yl[2] - 0.03 * diff(yl),
                                        label = note, hjust = 0, vjust = 1, size = 2.5, label.size = 0, fill = "white", alpha = 0.85)
  p
}
fold_set <- function(sub, fs) {
  note <- sprintf("block %.0f km%s", fs$size / 1000, if (isTRUE(fs$capped)) " (capped)" else "")
  if (!identical(fs$method, "block")) note <- paste0(note, "; ", fs$method)
  note
}
lims_for <- function(sub, coverage) {
  if (coverage == "clustered") {
    yl <- range(sub$lat) + c(-0.35, 0.35)
    xl <- mean(range(sub$lon)) + c(-1, 1) * max(diff(range(sub$lon)) / 2 + 0.4, diff(yl) * 0.75)
    xb <- pretty(xl, n = 3); xb <- xb[xb > xl[1] + 0.25 & xb < xl[2] - 0.25]
    yb <- pretty(yl, n = 3); yb <- yb[yb > yl[1] & yb < yl[2]]
  } else {
    xl <- c(bb["xmin"], bb["xmax"]); yl <- c(bb["ymin"], bb["ymax"]); xb <- c(-124, -120); yb <- c(35, 40, 45)
  }
  list(xl = xl, yl = yl, xb = xb, yb = yb)
}

## Main figure 2: n = 200, both coverages, both schemes (from the same draws as above)
f2_rows <- lapply(c("dispersed", "clustered"), function(cv) {
  i <- which(scen$target == "conditional" & scen$coverage == cv & scen$n == 200)
  sub <- draws[[i]]$sub
  fr <- cv_folds_random(sub); fsp <- cv_folds_spatial(sub)
  L <- lims_for(sub, cv)
  pr <- fold_panel(sub, fr, sprintf("%s, n = 200\nrandom 5 fold", nice_cov[cv]), L$xl, L$yl, L$xb, L$yb)
  ps <- fold_panel(sub, fsp$folds, sprintf("%s, n = 200\nspatial block 5 fold", nice_cov[cv]), L$xl, L$yl, L$xb, L$yb,
                   note = fold_set(sub, fsp))
  pr | ps
})
f2 <- (f2_rows[[1]] / f2_rows[[2]]) + plot_layout(guides = "collect", heights = c(1.7, 1)) & theme(legend.position = "right")
ggsave("results/analysis_out/FigS_design_2_cv_folds.png", f2, width = 7.5, height = 9.5, dpi = 300, bg = "white")
ggsave("results/analysis_out/FigS_design_2_cv_folds.pdf", f2, width = 7.5, height = 9.5, bg = "white")

## Supplementary S3 and S4: folds at every n, one figure per coverage
for (cv in c("dispersed", "clustered")) {
  panels <- list()
  for (nn in N_LEVELS) {
    i <- which(scen$target == "conditional" & scen$coverage == cv & scen$n == nn)
    sub <- draws[[i]]$sub
    fr <- cv_folds_random(sub); fsp <- cv_folds_spatial(sub)
    L <- lims_for(sub, cv)
    panels[[length(panels) + 1]] <- fold_panel(sub, fr, sprintf("n = %d: random", nn), L$xl, L$yl, L$xb, L$yb)
    panels[[length(panels) + 1]] <- fold_panel(sub, fsp$folds, sprintf("n = %d: spatial block", nn), L$xl, L$yl, L$xb, L$yb,
                                               note = fold_set(sub, fsp))
  }
  ## arrange: columns = n, rows = scheme
  ord <- c(seq(1, 9, 2), seq(2, 10, 2))
  fsx <- wrap_plots(panels[ord], ncol = 5, byrow = TRUE) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  fn <- if (cv == "dispersed") "FigS_design_S3_folds_dispersed" else "FigS_design_S4_folds_clustered"
  h <- if (cv == "dispersed") 8.5 else 5.2
  ggsave(sprintf("results/analysis_out/%s.png", fn), fsx, width = 11, height = h, dpi = 300, bg = "white")
  ggsave(sprintf("results/analysis_out/%s.pdf", fn), fsx, width = 11, height = h, bg = "white")
}

## ---- 7. Record what was drawn ----------------------------------------------
drawn <- scen |> mutate(species = sp, seed = FIG_SEED,
                        extent_frac = sapply(draws, function(x) x$extent),
                        redraws = sapply(draws, function(x) attr(x$sub, "redraws")),
                        presences = sapply(draws, function(x) sum(x$sub$y)))
write.csv(drawn, "results/analysis_out/T_design_figure_draws.csv", row.names = FALSE)
cat("Species drawn:", sp, "; seed", FIG_SEED, "\n")
print(drawn)
