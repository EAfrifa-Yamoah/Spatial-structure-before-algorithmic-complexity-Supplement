## Figure 2. Evidence map as a matrix, drawn from Appendix S3 (coding sheet).
##
## Rows    : the contrasts of Table 1.
## Columns : sample size bands (below 100; 100 to 500; above 500; structural or NA).
## Cell    : number of marine sources (single hue sequential ramp, colour blind safe);
##           non marine, simulated and mixed system sources as grey outlined markers.
##
## Rules applied (state these in the caption or Appendix S3):
##  * A study contributes to every band its reported sample size range spans
##    (n_min to n_max). A study with n_min and n_max both NA, or a methods paper
##    coded "structural guidance", falls in the "structural or NA" band.
##  * A study can contribute to more than one row, as in Table 1.
##  * "Marine" means marine == "yes". Studies coded "partial" (mixed systems)
##    are counted with the non marine markers.
##  * The counts are generated from the coding sheet, not from Table 1. Any
##    difference between the two is a signal to reconcile the sheet and Table 1
##    before submission; the script prints both for comparison.
##
## Packages: readr, dplyr, tidyr, stringr, ggplot2, ragg (optional, for PNG)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)
})

sheet_path <- "data/Appendix_S3_evidence_map_coding_sheet.csv"
out_png    <- "results/analysis_out/Figure2_evidence_map.png"
out_pdf    <- "results/analysis_out/Figure2_evidence_map.pdf"

cs <- read_csv(sheet_path, show_col_types = FALSE) |>
  filter(!is.na(year), !str_detect(row_status, "add")) |>
  mutate(n_min = suppressWarnings(as.numeric(n_min)),
         n_max = suppressWarnings(as.numeric(n_max)),
         n_max = coalesce(n_max, n_min))

## ---- 1. Row membership (one row of the matrix per contrast) -----------------
crossed   <- cs$contrasts_crossed == "yes"
is_app    <- cs$contrast_spatial_vs_nonspatial == "no" &
             cs$contrast_algorithm_family == "no" &
             cs$contrast_coverage_or_n == "no" &
             cs$contrast_validation_design == "no" &
             cs$metric_family %in% c("application", "projection outcome")

rows <- list(
  "Spatial against non spatial\nwithin a model family" =
      cs$contrast_spatial_vs_nonspatial == "yes",
  "Algorithm family, spatial\ntreatment held constant" =
      cs$contrast_algorithm_family == "yes",
  "Both contrasts crossed\non the same data" = crossed,
  "Sampling coverage against\nrecord count" =
      cs$contrast_coverage_or_n == "yes",
  "Random against blocked\nvalidation" =
      cs$contrast_validation_design == "yes",
  "Minimum sample size\nexperiment or requirement" =
      str_detect(cs$finding_supported, "minimum"),
  "Applied demonstrations" = is_app
)

## ---- 2. Band membership (a study spans every band its n range covers) -------
bands <- c("Below 100", "100 to 500", "Above 500", "Structural or\nnot applicable")
band_member <- function(n_min, n_max, structural) {
  if (structural || is.na(n_min)) return(bands[4])
  b <- character(0)
  if (n_min < 100)               b <- c(b, bands[1])
  if (n_max >= 100 & n_min <= 500) b <- c(b, bands[2])
  if (n_max > 500)               b <- c(b, bands[3])
  b
}
structural <- cs$contrast_coverage_or_n == "structural guidance"

long <- bind_rows(lapply(names(rows), function(r) {
  idx <- which(rows[[r]] %in% TRUE)
  bind_rows(lapply(idx, function(i) {
    tibble(contrast = r,
           band     = band_member(cs$n_min[i], cs$n_max[i], structural[i]),
           study    = cs$short_ref[i],
           marine   = cs$marine[i] == "yes")
  }))
}))

cells <- long |>
  group_by(contrast, band) |>
  summarise(non_marine = sum(!marine),
            marine     = sum(marine),
            studies    = paste(study, collapse = "; "),
            .groups = "drop") |>
  complete(contrast = names(rows), band = bands,
           fill = list(marine = 0L, non_marine = 0L, studies = "")) |>
  mutate(contrast = factor(contrast, levels = rev(names(rows))),
         band     = factor(band, levels = bands))

## ---- 3. Comparison against Table 1 (printed, not plotted) -------------------
by_row <- long |>
  group_by(contrast) |>
  summarise(sources        = n_distinct(study),
            marine         = n_distinct(study[marine]),
            below_100      = n_distinct(study[band == bands[1]]),
            marine_below_100 = n_distinct(study[marine & band == bands[1]]),
            .groups = "drop") |>
  mutate(contrast = str_replace_all(contrast, "\n", " "))
cat("\nCounts derived from the coding sheet (compare with Table 1):\n")
print(as.data.frame(by_row), row.names = FALSE)
cat("\nCell contents (for checking):\n")
print(as.data.frame(cells |> filter(studies != "") |>
                      arrange(desc(contrast), band) |>
                      mutate(contrast = str_replace_all(contrast, "\n", " "),
                             band = str_replace_all(band, "\n", " "))),
      row.names = FALSE)

## ---- 4. Draw ----------------------------------------------------------------
## Single hue sequential ramp (blue), zero drawn as white so empty cells read as empty.
ramp <- c("#deebf7", "#9ecae1", "#4292c6", "#08519c")   # light to dark blue
max_m <- max(cells$marine)

p <- ggplot(cells, aes(band, contrast)) +
  geom_tile(aes(fill = marine), colour = "grey60", linewidth = 0.4, width = 0.96, height = 0.96) +
  ## marine count printed in the cell centre
  geom_text(data = filter(cells, marine > 0),
            aes(label = marine, colour = marine > max_m * 0.6),
            size = 5, fontface = "bold", show.legend = FALSE) +
  ## non marine, simulated and mixed sources: grey outlined marker, lower right
  geom_point(data = filter(cells, non_marine > 0),
             aes(x = as.numeric(band) + 0.30, y = as.numeric(contrast) - 0.28,
                 shape = "Non marine, simulated or\nmixed system sources (count)"),
             size = 7, stroke = 0.9, colour = "grey30", fill = "white") +
  scale_shape_manual(values = c(21), name = NULL) +
  geom_text(data = filter(cells, non_marine > 0),
            aes(x = as.numeric(band) + 0.30, y = as.numeric(contrast) - 0.28,
                label = non_marine),
            size = 3.2, colour = "grey30") +
  ## dashed frame around the below 100 column
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = length(rows) + 0.5,
           fill = NA, colour = "grey20", linetype = "22", linewidth = 0.6) +
  scale_fill_gradientn(colours = c("white", ramp), values = c(0, 0.001, 0.33, 0.66, 1),
                       limits = c(0, max_m), breaks = 0:max_m,
                       name = "Marine sources") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 0.7, clip = "off") +
  labs(x = "Sample size band", y = NULL) +
  guides(fill = guide_colourbar(barheight = unit(3.2, "cm"), barwidth = unit(0.4, "cm"),
                                frame.colour = "grey60", ticks.colour = "grey60",
                                order = 1),
         shape = guide_legend(order = 2, override.aes = list(size = 6))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 10, lineheight = 0.9),
        axis.text.y = element_text(size = 10, lineheight = 0.9, hjust = 1),
        axis.title.x.top = element_text(size = 10, margin = margin(b = 6)),
        legend.position = "right",
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 9, lineheight = 0.9),
        legend.spacing.y = unit(0.6, "cm"),
        plot.margin = margin(6, 6, 6, 6))

fig <- p

w <- 8.2; h <- 5.6
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(out_png, width = w, height = h, units = "in", res = 300)
  print(fig); dev.off()
} else {
  ggsave(out_png, fig, width = w, height = h, dpi = 300, bg = "white")
}
ggsave(out_pdf, fig, width = w, height = h, bg = "white")
cat("\nWrote", out_png, "and", out_pdf, "\n")
