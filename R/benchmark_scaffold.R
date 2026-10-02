## =====================================================================
## Benchmark scaffold: spatial structure or algorithm choice?
## Structured subsamples of a harmonised marine trawl survey.
##
## The loader is written against the documented surveyjoin interface
## (cache_data, load_sql_data, get_data, surv_db). Run from the repository
## root. Without BENCH_RUN_FULL=TRUE the script runs the smoke test in
## Section 11; with it, the full design in Section 9.
##
## Design (manuscript Section 2):
##   sample sizes   n  in {30, 50, 100, 200, 500}
##   coverage       in {"dispersed", "clustered"}
##   targets        in {"conditional", "marginal"}
##   models         RF, spatial RF, GAM, spatial GAM, sdmTMB (GMRF)
##   replicates     BENCH_REPS per cell (10 for the committed results)
##   metrics        AUC, Brier score, calibration slope
##   validation     truth on held out set; random 5 fold; spatial block 5 fold
## =====================================================================

## ---- 0. Packages and reproducibility ---------------------------------
pkgs <- c("dplyr", "tidyr", "sf", "ranger", "mgcv", "sdmTMB", "blockCV",
          "future.apply", "surveyjoin", "DBI", "rnaturalearth")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "),
  "\n  surveyjoin: remotes::install_github('DFO-NOAA-Pacific/surveyjoin')",
  "\n  rnaturalearth uses the 1:50m coastline from rnaturalearthdata")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(sf); library(ranger); library(mgcv)
  library(sdmTMB); library(blockCV); library(future.apply); library(surveyjoin)
})
## Package versions used for the committed results: results/sessioninfo_smoke_test.txt
set.seed(20260919)

## ---- 1. Data ----------------------------------------------------------
## Settings. One survey keeps sampling design and gear constant, which is
## what the benchmark needs. NWFSC.Combo (US West Coast groundfish bottom
## trawl, annual from 2003) is the default; it is the system of Thorson
## et al. (2015). Other options from get_survey_names(): "SYN QCS",
## "SYN WCVI", "Gulf of Alaska", "eastern Bering Sea".
SURVEY   <- "NWFSC.Combo"
YEARS    <- 2003:2023            # adjust to what the cached data holds
UTM_EPSG <- 32610                # zone 10N for the US West Coast and BC;
                                 # use 32603 to 32609 for Alaska surveys
## Candidate species by common name as used in get_species(). Prevalence
## is computed below and the final choice is made from that table, not
## from these guesses. Aim: one common (> 40%), one or two intermediate
## (10 to 40%), one rare (< 10%).
CANDIDATES <- c("dover sole", "sablefish", "north pacific hake", "lingcod",
                "petrale sole", "canary rockfish")

## surveyjoin workflow as documented: cache_data() downloads the files,
## load_sql_data() builds the local SQLite database, get_data() queries it.
## Run the first two once; they are slow and idempotent.
if (!dir.exists(get_cache_folder())) cache_data()
load_sql_data()

## get_data() returns one row per haul by species with columns including
## survey_name, event_id, date, year, lat_start, lon_start, depth_m,
## effort, effort_units, performance, bottom_temp_c, catch_numbers,
## catch_weight, common_name. Whether hauls with zero catch are included is
## not stated in the documentation, so we check, and if absences are
## missing we rebuild them from the haul table through surv_db().
load_hauls <- function(species_names, survey = SURVEY, years = YEARS) {
  d <- get_data(common = species_names, surveys = survey, years = years)
  d <- as.data.frame(d)
  stopifnot(all(c("event_id", "lat_start", "lon_start", "depth_m", "effort",
                  "year", "catch_weight", "common_name") %in% names(d)))
  d$common_name <- tolower(d$common_name)

  ## Keep usable hauls
  if ("performance" %in% names(d)) {
    d <- d[is.na(d$performance) | tolower(d$performance) == "satisfactory", ]
  }
  d <- d[!is.na(d$effort) & d$effort > 0 & !is.na(d$lat_start) &
         !is.na(d$lon_start) & !is.na(d$depth_m), ]
  d$bottom_temp <- if ("bottom_temp_c" %in% names(d)) d$bottom_temp_c else NA_real_

  ## Wide: one row per haul, one catch column per species
  hauls <- d %>%
    distinct(event_id, year, lat_start, lon_start, depth_m, effort, bottom_temp)
  catch <- d %>%
    group_by(event_id, common_name) %>%
    summarise(cw = sum(catch_weight, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = common_name, values_from = cw, values_fill = 0)
  hauls <- left_join(hauls, catch, by = "event_id")
  for (s in species_names) if (!s %in% names(hauls)) hauls[[s]] <- 0

  ## Absence check without database access. If get_data() omitted zero
  ## catch hauls, the number of hauls returned would depend on which species
  ## were requested: a rare species alone would return far fewer hauls than
  ## the candidate set. If the count is the same, absences are included.
  n_all <- nrow(hauls)
  n_rare <- tryCatch({
    rare <- names(which.min(sapply(species_names, function(s) mean(hauls[[s]] > 0))))
    dr <- as.data.frame(get_data(common = rare, surveys = survey, years = years))
    length(unique(dr$event_id))
  }, error = function(e) NA_integer_)
  message(n_all, " hauls for the candidate set; ", ifelse(is.na(n_rare), "unknown", n_rare),
          " hauls when requesting the rarest candidate alone")
  if (!is.na(n_rare) && n_rare < 0.95 * n_all) {
    stop("get_data() appears to omit zero catch hauls (", n_rare, " vs ", n_all,
         "). Absences must be rebuilt from the haul table; see fill_absences_from_db().")
  }
  if (!is.na(n_rare)) message("Haul count is invariant to the species requested: absences are included.")

  hauls <- hauls %>%
    rename(haul_id = event_id, lon = lon_start, lat = lat_start, depth = depth_m) %>%
    mutate(across(all_of(species_names), ~ replace_na(.x, 0)))
  hauls
}

## Database helpers. Table and column names inside the SQLite database are
## not documented on the package site; both helpers discover them by name
## pattern and report what they found, so the lookups can be fixed once.
find_haul_table <- function(con) {
  tabs <- DBI::dbListTables(con)
  ht <- tabs[grepl("haul", tabs, ignore.case = TRUE)]
  if (!length(ht)) {
    message("Tables in the surveyjoin database: ", paste(tabs, collapse = ", "))
    return(NA_character_)
  }
  ht[1]
}

get_surv_db <- function() {
  ## surv_db() is documented online but not exported in every installed
  ## version; fall back to the internal function if present.
  if (exists("surv_db", where = asNamespace("surveyjoin"), inherits = FALSE))
    return(get("surv_db", envir = asNamespace("surveyjoin"))())
  stop("surv_db() not available in this surveyjoin version")
}

count_db_hauls <- function(survey, years) {
  n <- tryCatch({
    con <- get_surv_db(); on.exit(DBI::dbDisconnect(con), add = TRUE)
    ht <- find_haul_table(con); if (is.na(ht)) return(NA_integer_)
    flds <- DBI::dbListFields(con, ht)
    survey_col <- flds[grepl("survey", flds, ignore.case = TRUE)][1]
    year_col   <- flds[grepl("^year$", flds, ignore.case = TRUE)][1]
    if (is.na(survey_col) || is.na(year_col)) {
      message("Columns in ", ht, ": ", paste(flds, collapse = ", "))
      return(NA_integer_)
    }
    q <- sprintf("SELECT COUNT(*) AS n FROM %s WHERE %s = '%s' AND %s BETWEEN %d AND %d",
                 ht, survey_col, survey, year_col, min(years), max(years))
    DBI::dbGetQuery(con, q)$n
  }, error = function(e) { message("count_db_hauls: ", conditionMessage(e)); NA_integer_ })
  as.integer(n)
}

fill_absences_from_db <- function(hauls, species_names, survey, years) {
  con <- get_surv_db(); on.exit(DBI::dbDisconnect(con), add = TRUE)
  ht <- find_haul_table(con)
  if (is.na(ht)) stop("No haul table found; inspect DBI::dbListTables(surv_db())")
  flds <- DBI::dbListFields(con, ht)
  message("Columns in ", ht, ": ", paste(flds, collapse = ", "))
  survey_col <- flds[grepl("survey", flds, ignore.case = TRUE)][1]
  all_hauls <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE %s = '%s'", ht, survey_col, survey))
  all_hauls <- as.data.frame(all_hauls)
  ## Map database column names onto the names used above. Edit here if the
  ## printed column list differs.
  need <- c(event_id = "event_id", year = "year", lat_start = "lat_start",
            lon_start = "lon_start", depth_m = "depth_m", effort = "effort")
  miss <- need[!need %in% names(all_hauls)]
  if (length(miss)) stop("Haul table lacks: ", paste(miss, collapse = ", "),
                         ". Edit the `need` mapping in fill_absences_from_db().")
  all_hauls$bottom_temp <- if ("bottom_temp_c" %in% names(all_hauls)) all_hauls$bottom_temp_c else NA_real_
  all_hauls <- all_hauls[all_hauls$year %in% years &
                         !is.na(all_hauls$effort) & all_hauls$effort > 0 &
                         !is.na(all_hauls$lat_start) & !is.na(all_hauls$lon_start) &
                         !is.na(all_hauls$depth_m), ]
  out <- left_join(select(all_hauls, event_id, year, lat_start, lon_start, depth_m, effort, bottom_temp),
                   select(hauls, event_id, all_of(species_names)), by = "event_id")
  out <- mutate(out, across(all_of(species_names), ~ replace_na(.x, 0)))
  message(nrow(out), " hauls after filling absences")
  out
}

hauls <- load_hauls(CANDIDATES)
n_before <- nrow(hauls)
hauls <- hauls[!is.na(hauls$bottom_temp), ]
message(nrow(hauls), " hauls loaded from ", SURVEY, " for ", length(YEARS), " years (",
        n_before - nrow(hauls), " dropped for missing bottom temperature)")

## Project to km for spatial models and distances
to_km <- function(df, epsg = UTM_EPSG) {
  pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  xy  <- st_coordinates(st_transform(pts, crs = epsg)) / 1000
  df$X <- xy[, 1]; df$Y <- xy[, 2]
  df
}
hauls <- to_km(hauls)

## Ecologically meaningful coordinate axes for tree based learners.
## Distance to shore from the Natural Earth coastline. Distance to the
## shelf break needs a 200 m isobath; if you have one as an sf LINESTRING
## pass it as `shelf`, otherwise the axis is dropped and depth carries the
## cross shelf information.
add_axes <- function(df, epsg = UTM_EPSG, shelf = NULL) {
  coast <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")  # 1:50m from rnaturalearthdata; scale = 10 needs rnaturalearthhires
  coast <- st_transform(coast, epsg)
  pts   <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  pts   <- st_transform(pts, epsg)
  df$dist_shore <- as.numeric(st_distance(pts, st_union(coast))) / 1000
  if (!is.null(shelf)) {
    df$dist_shelf <- as.numeric(st_distance(pts, st_transform(shelf, epsg))) / 1000
  }
  df
}
hauls <- add_axes(hauls)
axes <- intersect(c("X", "Y", "dist_shore", "dist_shelf"), names(hauls))

## ---- 2. Species -------------------------------------------------------
prevalence <- sort(sapply(CANDIDATES, function(s) mean(hauls[[s]] > 0)), decreasing = TRUE)
print(round(prevalence, 3))
## Choose from the printed table before continuing: one common, one or two
## intermediate, one rare. Record the choice and the prevalence in Table 2.
species <- c(
  names(prevalence)[which(prevalence > 0.40)][1],
  names(prevalence)[which(prevalence > 0.10 & prevalence <= 0.40)][1],
  names(prevalence)[which(prevalence <= 0.10 & prevalence > 0.02)][1]
)
species <- species[!is.na(species)]
if (length(species) < 3) warning("Fewer than three prevalence classes found among CANDIDATES; add candidates.")
message("Species selected: ", paste(species, collapse = "; "))

## Presence absence response
make_response <- function(df, sp) { df$y <- as.integer(df[[sp]] > 0); df }

## ---- 3. Evaluation partitions (once per species) ----------------------
## Conditional target: random 30% of hauls held out.
## Marginal target: contiguous region, the northern 30% of hauls by Y
##   (the survey's dominant gradient runs north to south).
make_partitions <- function(df) {
  n <- nrow(df)
  cond_idx <- sample(n, size = round(0.3 * n))
  ycut <- quantile(df$Y, 0.7)
  marg_idx <- which(df$Y > ycut)
  list(
    conditional = list(pool = df[-cond_idx, ], eval = df[cond_idx, ]),
    marginal    = list(pool = df[-marg_idx, ], eval = df[marg_idx, ])
  )
}

## ---- 4. Subsampling ---------------------------------------------------
## Dispersed: stratified random across a regular grid over the pool.
## Clustered: random seed haul and its n nearest neighbours.
CLUSTER_FRAC <- 0.10   # clustered coverage: band width as a fraction of the pool's north to south extent

subsample <- function(pool, n, coverage, grid_km = 25, min_pres = 5, cluster_frac = CLUSTER_FRAC) {
  for (attempt in 1:50) {
    if (coverage == "dispersed") {
      ## Stratified across a regular grid so the subsample spans the pool
      gx <- floor(pool$X / grid_km); gy <- floor(pool$Y / grid_km)
      cell <- paste(gx, gy)
      cells <- unique(cell)
      take_cells <- cells[sample.int(length(cells), min(n, length(cells)))]
      ## One haul per chosen cell. Index into the vector of candidates rather
      ## than calling sample() on it: sample(x, 1) with a single candidate x
      ## draws from 1:x instead of returning x, which would pick an arbitrary
      ## haul (and possibly a duplicate) whenever a grid cell holds one haul.
      pick <- vapply(take_cells, function(cc) {
        idx <- which(cell == cc)
        idx[sample.int(length(idx), 1)]
      }, integer(1))
      if (length(pick) < n) {
        rest <- setdiff(seq_len(nrow(pool)), pick)
        pick <- c(pick, rest[sample.int(length(rest), n - length(pick))])
      }
    } else {
      ## A band of fixed width (cluster_frac of the north to south extent),
      ## placed at random, with n hauls sampled at random inside it. Extent
      ## is therefore controlled independently of n, unlike nearest
      ## neighbour clustering, which shrinks to a single station at small n.
      yr <- range(pool$Y); width <- cluster_frac * diff(yr)
      y0 <- runif(1, yr[1], yr[2] - width)
      inband <- which(pool$Y >= y0 & pool$Y < y0 + width)
      if (length(inband) < n) next
      pick <- sample(inband, n)
    }
    sub <- pool[pick, ]
    if (sum(sub$y) >= min_pres && sum(sub$y == 0) >= min_pres) {
      attr(sub, "redraws") <- attempt - 1
      ## realised extent as fraction of pool extent (operationalises coverage)
      attr(sub, "extent_frac") <- (diff(range(sub$X)) * diff(range(sub$Y))) /
                                  (diff(range(pool$X)) * diff(range(pool$Y)))
      return(sub)
    }
  }
  NULL
}

## ---- 5. Models --------------------------------------------------------
## Covariates: kept small and fixed. Hauls lacking bottom_temp are dropped
## in Section 1 (258 of 12672 in the committed run; see results/pilot_run_log.txt).
covs <- c("depth", "bottom_temp", "year")   # year numeric; hauls lacking bottom_temp dropped below

k_for_n <- function(n) max(3, min(10, floor(n / 10)))   # basis dimension cap
mesh_cutoff_for <- function(sub) {
  ## cutoff so that nodes stay well below n; tune once on n = 500 and fix
  ext <- max(diff(range(sub$X)), diff(range(sub$Y)))
  max(ext / 8, 5)
}

fit_predict <- function(model, train, test) {
  n <- nrow(train); k <- k_for_n(n)
  out <- list(p = rep(NA_real_, nrow(test)), converged = TRUE, sanity_ok = NA,
              hessian_ok = NA, nlminb_ok = NA, gradients_ok = NA, fit_error = NA_character_)
  tryCatch({
    if (model == "rf") {
      f <- ranger(as.factor(y) ~ ., data = train[, c("y", covs)],
                  probability = TRUE, num.trees = 500, min.node.size = 5)
      out$p <- predict(f, data = test[, covs, drop = FALSE])$predictions[, "1"]
    } else if (model == "rf_spatial") {
      f <- ranger(as.factor(y) ~ ., data = train[, c("y", covs, axes)],
                  probability = TRUE, num.trees = 500, min.node.size = 5)
      out$p <- predict(f, data = test[, c(covs, axes), drop = FALSE])$predictions[, "1"]
    } else if (model == "gam") {
      form <- as.formula(paste("y ~", paste(sprintf("s(%s, k = %d)", covs, k), collapse = " + ")))
      f <- gam(form, data = train, family = binomial(), method = "REML")
      out$p <- as.numeric(predict(f, newdata = test, type = "response"))
    } else if (model == "gam_spatial") {
      form <- as.formula(paste("y ~", paste(sprintf("s(%s, k = %d)", covs, k), collapse = " + "),
                               sprintf("+ s(X, Y, k = %d)", max(6, min(30, floor(n / 5))))))
      f <- gam(form, data = train, family = binomial(), method = "REML")
      out$p <- as.numeric(predict(f, newdata = test, type = "response"))
    } else if (model == "gmrf") {
      mesh <- make_mesh(train, xy_cols = c("X", "Y"), cutoff = mesh_cutoff_for(train))
      form <- as.formula(paste("y ~", paste(sprintf("s(%s, k = %d)", covs, k), collapse = " + ")))
      f <- sdmTMB(form, data = train, mesh = mesh, family = binomial(), spatial = "on")
      s <- sanity(f, silent = TRUE)
      out$hessian_ok <- isTRUE(s$hessian_ok); out$nlminb_ok <- isTRUE(s$nlminb_ok)
      out$gradients_ok <- isTRUE(s$gradients_ok)
      out$converged <- out$hessian_ok && out$nlminb_ok && out$gradients_ok
      out$sanity_ok <- isTRUE(s$all_ok)      # also range, sigma and SE checks
      pr <- predict(f, newdata = test)
      out$p <- plogis(pr$est)
    }
  }, error = function(e) { out$converged <<- FALSE; out$fit_error <<- conditionMessage(e) })
  out
}

## ---- 6. Metrics -------------------------------------------------------
auc <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
brier <- function(y, p) mean((p - y)^2)
cal_slope <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  if (length(unique(y)) < 2) return(NA_real_)
  co <- tryCatch(coef(glm(y ~ qlogis(p), family = binomial()))[2], error = function(e) NA_real_)
  as.numeric(co)
}
metrics <- function(y, p) c(auc = auc(y, p), brier = brier(y, p), cal_slope = cal_slope(y, p))

## ---- 7. Cross validation within a subsample ---------------------------
cv_folds_random <- function(sub, k = 5) sample(rep(seq_len(k), length.out = nrow(sub)))

cv_folds_spatial <- function(sub, k = 5, epsg = UTM_EPSG) {
  ## Points in a projected CRS in metres, built from lon/lat so that block
  ## sizes and ranges are in consistent units.
  pts <- st_transform(st_as_sf(sub, coords = c("lon", "lat"), crs = 4326), epsg)
  bb  <- st_bbox(pts)
  ext_m <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])

  ## Residual autocorrelation range from a preliminary spatial GAM; fall
  ## back to the response if that fails. Units: metres.
  rng_m <- tryCatch({
    g <- gam(y ~ s(depth, k = k_for_n(nrow(sub))) + s(X, Y, k = min(10, max(4, floor(nrow(sub) / 6)))),
             data = sub, family = binomial(), method = "REML")
    pts$res <- residuals(g, type = "pearson")
    cv_spatial_autocor(x = pts, column = "res", plot = FALSE, progress = FALSE)$range
  }, error = function(e) NA_real_)
  if (is.na(rng_m) || !is.finite(rng_m) || rng_m <= 0) {
    rng_m <- tryCatch(cv_spatial_autocor(x = pts, column = "y", plot = FALSE, progress = FALSE)$range,
                      error = function(e) NA_real_)
  }
  if (is.na(rng_m) || !is.finite(rng_m) || rng_m <= 0) rng_m <- ext_m / 5

  ## Block size: the range, capped so that at least about k blocks fit the
  ## subsample extent. When the cap binds, blocks are smaller than the
  ## range and some leakage remains; that is recorded, because it is a real
  ## constraint on small clustered samples and belongs in the results.
  size_m <- min(rng_m, ext_m / ceiling(sqrt(k)))
  capped <- size_m < rng_m

  folds <- tryCatch(
    cv_spatial(x = pts, column = "y", k = k, size = size_m, selection = "random",
               iteration = 50, progress = FALSE, report = FALSE, plot = FALSE)$folds_ids,
    error = function(e) NULL)
  method <- "block"
  if (is.null(folds) || length(unique(folds)) < 2) {
    ## Spatial clustering always yields k folds; used only when blocking fails
    folds <- cv_cluster(x = pts, k = k)$folds_ids
    method <- "cluster"
  }
  list(folds = folds, range = rng_m, size = size_m, capped = capped, method = method)
}

cv_estimate <- function(model, sub, folds) {
  p <- rep(NA_real_, nrow(sub))
  for (f in sort(unique(folds))) {
    tr <- sub[folds != f, ]; te <- sub[folds == f, ]
    if (length(unique(tr$y)) < 2) next
    p[folds == f] <- fit_predict(model, tr, te)$p
  }
  ok <- !is.na(p)
  metrics(sub$y[ok], p[ok])
}

## ---- 8. One cell by replicate -----------------------------------------
run_one <- function(sp, target, n, coverage, rep, parts) {
  tryCatch(run_one_inner(sp, target, n, coverage, rep, parts),
           error = function(e) data.frame(species = sp, target = target, n = n,
             coverage = coverage, rep = rep, model = NA_character_, converged = FALSE,
             error = conditionMessage(e), stringsAsFactors = FALSE))
}

run_one_inner <- function(sp, target, n, coverage, rep, parts) {
  t0 <- proc.time()[["elapsed"]]
  pool <- parts[[target]]$pool; eval <- parts[[target]]$eval
  sub <- subsample(pool, n, coverage)
  if (is.null(sub)) return(data.frame(species = sp, target = target, n = n, coverage = coverage,
                                      rep = rep, model = NA_character_, converged = FALSE,
                                      error = "subsample: presence minimum not met in 50 attempts",
                                      stringsAsFactors = FALSE))
  fr <- cv_folds_random(sub)
  fs <- cv_folds_spatial(sub)
  rows <- list()
  for (m in c("rf", "rf_spatial", "gam", "gam_spatial", "gmrf")) {
    fp <- fit_predict(m, sub, eval)
    truth <- if (all(is.na(fp$p))) c(auc = NA, brier = NA, cal_slope = NA) else metrics(eval$y, fp$p)
    est_r <- tryCatch(cv_estimate(m, sub, fr), error = function(e) c(auc = NA, brier = NA, cal_slope = NA))
    est_s <- tryCatch(cv_estimate(m, sub, fs$folds), error = function(e) c(auc = NA, brier = NA, cal_slope = NA))
    rows[[m]] <- data.frame(
      species = sp, target = target, n = n, coverage = coverage, rep = rep, model = m,
      converged = fp$converged, sanity_ok = fp$sanity_ok,
      hessian_ok = fp$hessian_ok, nlminb_ok = fp$nlminb_ok, gradients_ok = fp$gradients_ok,
      fit_error = fp$fit_error, extent_frac = attr(sub, "extent_frac"),
      redraws = attr(sub, "redraws"), block_range_m = fs$range, block_size_m = fs$size,
      block_capped = fs$capped, block_method = fs$method, error = NA_character_,
      elapsed_s = NA_real_,
      truth_auc = truth["auc"], truth_brier = truth["brier"], truth_cal = truth["cal_slope"],
      cvr_auc = est_r["auc"], cvr_brier = est_r["brier"],
      cvs_auc = est_s["auc"], cvs_brier = est_s["brier"],
      degenerate = isTRUE(sd(fp$p, na.rm = TRUE) < 1e-3),
      row.names = NULL)
  }
  out <- do.call(rbind, rows)
  out$elapsed_s <- proc.time()[["elapsed"]] - t0   # whole replicate, all models
  out
}

## ---- 9. Driver --------------------------------------------------------
## RUN_FULL = FALSE runs only the smoke test in Section 11. Set the
## environment variable BENCH_RUN_FULL=TRUE (or RUN_FULL <- TRUE) for the
## full design on the HPC. BENCH_WORKERS sets the number of cores.
RUN_FULL <- toupper(Sys.getenv("BENCH_RUN_FULL", "FALSE")) == "TRUE"
N_WORKERS <- as.integer(Sys.getenv("BENCH_WORKERS", max(1, parallel::detectCores() - 1)))
N_REPS <- as.integer(Sys.getenv("BENCH_REPS", "10"))

design <- expand.grid(n = c(30, 50, 100, 200, 500),
                      coverage = c("dispersed", "clustered"),
                      target = c("conditional", "marginal"),
                      rep = seq_len(N_REPS), stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE)
if (RUN_FULL) {
  plan(multisession, workers = N_WORKERS)
  results <- list()
  for (sp in species) {
    d <- make_response(hauls, sp)
    parts <- make_partitions(d)
    res_sp <- future_lapply(seq_len(nrow(design)), function(i) {
      r <- design[i, ]
      run_one(sp, r$target, r$n, r$coverage, r$rep, parts)
    }, future.seed = TRUE)
    results[[sp]] <- dplyr::bind_rows(res_sp)
    saveRDS(results[[sp]], sprintf("results/results_%s.rds", gsub(" ", "_", sp)))
  }
  res <- dplyr::bind_rows(results)
  message(sum(!is.na(res$error)), " replicate level errors recorded; see res$error")
  res$opt_random  <- res$cvr_auc - res$truth_auc     # optimism: positive overstates
  res$opt_spatial <- res$cvs_auc - res$truth_auc
  saveRDS(res, "results/results_all.rds")
}

## ---- 10. Analysis (sketch; implemented in R/analysis.R) ----------------
## library(lme4); library(emmeans)
## m_auc <- lmer(truth_auc ~ model * factor(n) * coverage * target + (1 | species/rep), data = res)
## spatial effect within family: emmeans contrasts rf_spatial - rf, gam_spatial - gam, by n, coverage
## algorithm effect: rf - gam, rf_spatial - gam_spatial, by n
## coverage effect: dispersed - clustered by n, model, target
## optimism: lmer(opt_random ~ model * factor(n) * coverage * target + (1 | species/rep))
## robustness: proportion of replicates with the contrast sign, by cell
## failure modes: aggregate(!converged | degenerate ~ model + n + coverage, res, mean)

## ---- 11. Smoke test (runs when RUN_FULL is FALSE) ---------------------
if (!RUN_FULL) {
  d <- make_response(hauls, species[1]); parts <- make_partitions(d)
  message("Smoke test on ", species[1], ": n = 500 dispersed conditional, then n = 30 clustered marginal")
  one500 <- run_one(species[1], "conditional", 500, "dispersed", 1, parts)
  one30  <- run_one(species[1], "marginal",     30, "clustered", 1, parts)
  smoke <- dplyr::bind_rows(one500, one30)
  print(as.data.frame(smoke[, intersect(c("n","coverage","target","model","converged","sanity_ok","truth_auc","truth_brier",
        "truth_cal","cvr_auc","cvs_auc","block_range_m","block_size_m","block_capped","block_method",
        "extent_frac","degenerate","elapsed_s","error"), names(smoke))]), digits = 3)
  hrs <- mean(c(one500$elapsed_s[1], one30$elapsed_s[1])) * nrow(design) * length(species) / 3600
  message(sprintf("Projected single core run time for the full design: %.1f hours (%d species, %d replicates per cell)",
                  hrs, length(species), N_REPS))
  saveRDS(smoke, "results/smoke_test.rds")
}
## Check: all five models return finite truth_auc at n = 500; gmrf converged;
## block_range_m and block_size_m are in metres and plausible (tens of km);
## block_capped is TRUE for clustered n = 30 and FALSE for dispersed n = 500;
## block_method is "block" in most cells; extent_frac is smaller for
## clustered than dispersed; the error column is NA.
