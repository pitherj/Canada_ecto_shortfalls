# =============================================================================
# Linnean Shortfall — Per-sample Sequence-based Rarefaction/Extrapolation
# =============================================================================
# Sample-specific abundance-based rarefaction + extrapolation of EcM SH
# richness for every Canadian GlobalFungi sample, following the iNEXT
# (analytical) approach used by van Galen and colleagues at SPUN. Each
# GlobalFungi sample is treated as a separate community and its read counts
# (per UNITE SH code) are passed to iNEXT::iNEXT() with q = 0,
# datatype = "abundance", and endpoint = 2 × sample sequencing depth (the
# package default; matches the SPUN published methods). For each sample we
# extract: observed SH richness (S.obs), Chao1 asymptotic richness with 95 %
# confidence interval, sample completeness Ĉ at observed depth, and the
# extrapolated richness at 2 × depth with its 95 % CI.
#
# Caveats (relevant for interpretation, not for the computation):
#   - Reads are not biologically independent (PCR amplification bias and rRNA
#     copy-number variation), so per-sample Chao1 from reads tends to
#     over-estimate per-sample SH richness. The standardization is what makes
#     patterns comparable across samples.
#   - The "asymptote at 2 × depth" is iNEXT's default extrapolation endpoint,
#     not the curve's true plateau. Samples far from saturation will show
#     wide CIs and the endpoint estimate is closer to an extrapolated count
#     than an asymptotic richness.
#   - Per-sample estimates do not aggregate to a Canada-wide richness by
#     summing or averaging — that is what the site-based Chao2 in
#     09_linnean.R Step 5 addresses (in incidence space).
#
# Pipeline structure (each step independently sentinel-guarded):
#   Step 1.  Build per-sample EcM read-count vectors (cheap; always runs)
#   Step 2.  Run iNEXT per sample → cache to data_derived/checkpoints/
#            linnean_inext_per_sample.rds  (slow; sentinel-guarded by RDS)
#   Step 3.  Build per-sample CSV from the RDS  (cheap; always rebuilds)
#   Step 4.  Build distribution-summary CSV     (cheap; always rebuilds)
#   Step 5.  Build sample-location richness map (sentinel-guarded by PNG).
#            All samples with valid coordinates and a successful iNEXT
#            estimate are mapped (no Ĉ filter). The per-sample CSV retains
#            the same set.
#
# Column names in linnean_inext_per_sample.csv:
#   ecm_reads       — sum of reads in this sample mapped to EcM-genus SHs
#                     (NOT total sequencing depth; see comment in Step 1)
#   sh_obs          — observed unique SH count from the input vector
#   s_obs_inext     — observed richness as reported by iNEXT (≡ sh_obs)
#   s_est_chao1     — Chao1 asymptotic estimator
#   s_est_chao1_lci — lower 95% CI on Chao1
#   s_est_chao1_uci — upper 95% CI on Chao1
#   coverage_obs    — sample completeness Ĉ at observed depth (Chao & Jost 2012)
#   s_ext_2x        — richness at the 2× extrapolation endpoint
#   s_ext_2x_lci    — lower 95% CI at 2× endpoint
#   s_ext_2x_uci    — upper 95% CI at 2× endpoint
#   coverage_2x     — Ĉ at the 2× endpoint
#   s_est_over_obs  — ratio S.est / S.obs
#   inext_status    — "ok" or error message
#
# To force any step to re-run, delete its sentinel file and re-source this
# script. The map step in particular is independent of the iNEXT loop, so
# tweaks to the figure can be made without redoing the per-sample fits.
#
# Outputs:
#   paths$linnean_inext_rds         — data_derived/checkpoints/linnean_inext_per_sample.rds
#   paths$linnean_inext_per_sample  — data_derived/linnean/linnean_inext_per_sample.csv
#   paths$linnean_inext_summary     — data_derived/linnean/linnean_inext_summary.csv
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(iNEXT)
library(sf)
library(terra)
library(scales)

sf::sf_use_s2(FALSE)

MAP_DPI <- 300L

# Completeness reference level retained as a diagnostic metric in
# inext_summary (counts samples meeting Ĉ ≥ 0.5); no longer used to
# filter the map — all samples with valid coordinates and a successful
# iNEXT estimate are mapped.
MAP_MIN_COVERAGE <- 0.5

# ---- Step 1: Build per-sample EcM read-count vectors ------------------------

ts("Step 1: Building per-sample EcM read-count vectors (GlobalFungi)...")

gf <- emf |>
  dplyr::filter(source == "GlobalFungi",
                !is.na(sample_ID), !is.na(sh_code), !is.na(abundance),
                abundance > 0) |>
  dplyr::select(sample_ID, sh_code, abundance, lat, lon, coord_in_canada)

ts(sprintf("  GlobalFungi non-zero records: %d  |  unique samples: %d  |  unique SHs: %d",
           nrow(gf),
           dplyr::n_distinct(gf$sample_ID),
           dplyr::n_distinct(gf$sh_code)))

# Per-sample read-count vectors (named by SH for traceability)
sample_vecs <- split(gf, gf$sample_ID) |>
  lapply(function(d) setNames(as.integer(d$abundance), d$sh_code))

# Per-sample summary: post-EcM-filter read sum and observed SH count, plus
# the sample's coordinates pulled from the long-format records. NOTE:
# `ecm_reads` is the SUM of reads in this sample that mapped to UNITE SHs
# whose genus is flagged ectomycorrhizal in FungalTraits — it is NOT total
# sequencing depth. Variation across samples therefore reflects (a) the
# original sequencing-depth differences across GlobalFungi studies, (b)
# substrate-driven variation in the EcM proportion of the fungal community
# (root tip / mineral horizon EcM-forest samples >> grassland or peatland
# samples), and (c) primer / amplification effects. iNEXT is run on the
# per-sample EcM read vector, so this is the correct quantity to feed in,
# but downstream interpretation should treat the wide range as an artefact
# of pooled studies rather than a uniform sampling design. We pull lat/lon
# via a defensive as.numeric() because read_csv occasionally infers
# character types for all-NA leading columns.
sample_meta <- gf |>
  dplyr::group_by(sample_ID) |>
  dplyr::summarise(
    ecm_reads      = sum(abundance, na.rm = TRUE),
    sh_obs         = dplyr::n_distinct(sh_code),
    lat            = suppressWarnings(as.numeric(dplyr::first(lat))),
    lon            = suppressWarnings(as.numeric(dplyr::first(lon))),
    coord_in_canada = dplyr::first(coord_in_canada),
    .groups        = "drop"
  )

# Retain all samples with at least one EcM SH code detected; samples where
# iNEXT cannot compute an estimate (e.g. f2 = 0) are handled gracefully by
# inext_one() and flagged via inext_status.
sample_meta <- dplyr::filter(sample_meta, sh_obs >= 1L)
sample_vecs <- sample_vecs[as.character(sample_meta$sample_ID)]

ts(sprintf("  Samples retained for iNEXT: %d", nrow(sample_meta)))
ts(sprintf("  Of which have validated coordinates (coord_in_canada == TRUE): %d",
           sum(sample_meta$coord_in_canada == TRUE, na.rm = TRUE)))

# ---- Step 2: Run iNEXT per sample (sentinel-guarded by RDS checkpoint) ------
# Defensive AsyEst / size_based extraction: iNEXT column conventions have
# shifted across versions. We try multiple recognised patterns and fall back
# to "first row" semantics when q = 0 was specified (in which case AsyEst
# should have exactly one row per assemblage anyway).

inext_one <- function(x_vec, sample_id) {
  out <- tibble::tibble(
    sample_ID    = sample_id,
    s_obs_inext  = NA_real_,
    s_est_chao1  = NA_real_,
    s_est_chao1_lci = NA_real_,
    s_est_chao1_uci = NA_real_,
    coverage_obs = NA_real_,
    s_ext_2x     = NA_real_,
    s_ext_2x_lci = NA_real_,
    s_ext_2x_uci = NA_real_,
    coverage_2x  = NA_real_,
    inext_status = "ok"
  )

  res <- tryCatch(
    suppressWarnings(iNEXT::iNEXT(x_vec, q = 0, datatype = "abundance",
                                  nboot = 50, conf = 0.95)),
    error = function(e) {
      out$inext_status <<- paste0("error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res)) return(out)

  # ---- Asymptotic estimator (AsyEst) -----------------------------------------
  asy <- res$AsyEst
  if (!is.null(asy) && nrow(asy) > 0L) {
    # Filter to q = 0 across the various iNEXT versions
    pick <- rep(TRUE, nrow(asy))
    if ("Diversity" %in% names(asy)) {
      pick <- grepl("species\\s*richness", as.character(asy$Diversity),
                    ignore.case = TRUE) |
              as.character(asy$Diversity) %in% c("0", "q = 0", "q=0")
    } else if ("Order.q" %in% names(asy)) {
      pick <- asy$Order.q == 0
    }
    asy0 <- asy[pick, , drop = FALSE]
    # Defensive fallback: with q = 0 specified, AsyEst should have exactly
    # one row per assemblage; if the filter caught nothing, take row 1.
    if (nrow(asy0) == 0L) asy0 <- asy[1L, , drop = FALSE]

    # Column-name flexibility: Observed/Estimator have stable names across
    # versions; LCL/UCL are sometimes prefixed with "qD." in development
    # builds.
    pick_col <- function(df, candidates) {
      hit <- candidates[candidates %in% names(df)][1L]
      if (is.na(hit)) NA_real_ else df[[hit]][1L]
    }
    out$s_obs_inext     <- pick_col(asy0, c("Observed", "S.obs", "qD.obs"))
    out$s_est_chao1     <- pick_col(asy0, c("Estimator", "S.est", "qD.est"))
    out$s_est_chao1_lci <- pick_col(asy0, c("95% Lower", "LCL", "qD.LCL"))
    out$s_est_chao1_uci <- pick_col(asy0, c("95% Upper", "UCL", "qD.UCL"))
  }

  # ---- Size-based curve (observed-depth Ĉ + 2× endpoint) ---------------------
  sb <- res$iNextEst$size_based
  # Some iNEXT versions stash the size-based table directly under iNextEst
  if (is.null(sb) && is.data.frame(res$iNextEst))      sb <- res$iNextEst
  if (is.null(sb) && !is.null(res$iNextEst[[1L]]))     sb <- res$iNextEst[[1L]]

  if (!is.null(sb) && nrow(sb) > 0L) {
    if ("Order.q" %in% names(sb)) sb <- sb[sb$Order.q == 0, , drop = FALSE]
    if (nrow(sb) > 0L) {
      method_col <- if ("Method" %in% names(sb)) sb$Method
                    else if ("method" %in% names(sb)) sb$method
                    else NULL
      if (!is.null(method_col)) {
        ref <- sb[grepl("^observed",     method_col, ignore.case = TRUE), , drop = FALSE]
        ext <- sb[grepl("^extrapolation", method_col, ignore.case = TRUE), , drop = FALSE]
      } else {
        ref <- sb[1L, , drop = FALSE]
        ext <- sb[nrow(sb), , drop = FALSE]
      }

      sc_col <- intersect(c("SC", "Coverage", "qSC"), names(sb))[1L]
      qd_col <- intersect(c("qD", "qD.obs", "Estimator"), names(sb))[1L]

      if (nrow(ref) >= 1L && !is.na(sc_col)) out$coverage_obs <- ref[[sc_col]][1L]
      if (nrow(ext) >= 1L) {
        last <- ext[nrow(ext), , drop = FALSE]
        if (!is.na(qd_col)) out$s_ext_2x  <- last[[qd_col]][1L]
        if ("qD.LCL" %in% names(last)) out$s_ext_2x_lci <- last[["qD.LCL"]][1L]
        if ("qD.UCL" %in% names(last)) out$s_ext_2x_uci <- last[["qD.UCL"]][1L]
        if (!is.na(sc_col)) out$coverage_2x <- last[[sc_col]][1L]
      }
    }
  }
  out
}

if (file.exists(paths$linnean_inext_rds)) {

  ts("Step 2: Loading cached iNEXT results from RDS checkpoint...")
  inext_results <- readRDS(paths$linnean_inext_rds)
  ts(sprintf("  Loaded %d rows from %s",
             nrow(inext_results), basename(paths$linnean_inext_rds)))

} else {

  ts(sprintf("Step 2: Running iNEXT per sample (q = 0, datatype = 'abundance', endpoint = 2× depth, nboot = 50)..."))
  ts(sprintf("  Iterating over %d samples...", length(sample_vecs)))
  # Seed set once, immediately before the per-sample bootstrap loop, so the
  # Chao1 CIs (nboot = 50, via iNEXT::iNEXT() inside inext_one()) are
  # reproducible across runs of this checkpointed RDS.
  set.seed(3492)
  res_list <- vector("list", length(sample_vecs))
  ids <- names(sample_vecs)
  t0 <- Sys.time()
  for (i in seq_along(sample_vecs)) {
    res_list[[i]] <- inext_one(sample_vecs[[i]], ids[i])
    if (i %% 100L == 0L || i == length(sample_vecs)) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      ts(sprintf("    %d / %d samples complete  (%.1f s elapsed)",
                 i, length(sample_vecs), elapsed))
    }
  }
  inext_results <- dplyr::bind_rows(res_list)
  saveRDS(inext_results, paths$linnean_inext_rds)
  ts(sprintf("  Cached %d-row RDS checkpoint -> %s",
             nrow(inext_results), basename(paths$linnean_inext_rds)))

}

n_failed <- sum(inext_results$inext_status != "ok", na.rm = TRUE)
n_est_na <- sum(is.na(inext_results$s_est_chao1))
ts(sprintf("  iNEXT errors: %d  |  s_est_chao1 NA: %d  |  with estimate: %d",
           n_failed, n_est_na, sum(!is.na(inext_results$s_est_chao1))))

# ---- Step 3: Per-sample CSV (always rebuilds from RDS) ----------------------

ts("Step 3: Assembling per-sample CSV...")

per_sample <- dplyr::left_join(sample_meta, inext_results, by = "sample_ID") |>
  dplyr::mutate(
    s_est_over_obs = ifelse(is.na(s_est_chao1) | is.na(sh_obs) | sh_obs == 0,
                            NA_real_, s_est_chao1 / sh_obs)
  ) |>
  dplyr::arrange(dplyr::desc(s_est_chao1))

readr::write_csv(per_sample, paths$linnean_inext_per_sample)
ts(sprintf("  Saved %s (%d rows)",
           basename(paths$linnean_inext_per_sample), nrow(per_sample)))

# ---- Step 4: Distribution summary -------------------------------------------

ts("Step 4: Computing distribution summary...")

inext_summary <- tibble::tibble(
  metric = c(
    "Samples processed",
    "Samples with successful iNEXT estimation",
    "Median observed SH richness per sample (S.obs)",
    "Median Chao1 asymptotic SH richness per sample (S.est)",
    "Median S.est / S.obs ratio per sample",
    "Median sample completeness Ĉ at observed depth",
    "Samples below 90% completeness (Ĉ < 0.90)",
    "Samples below 90% completeness (% of estimated samples)",
    sprintf("Samples with Ĉ ≥ %.2f at observed depth (diagnostic)", MAP_MIN_COVERAGE),
    "Median EcM read sum per sample (post-EcM-filter; not total sequencing depth)"
  ),
  value = c(
    nrow(per_sample),
    sum(per_sample$inext_status == "ok", na.rm = TRUE),
    round(stats::median(per_sample$sh_obs,         na.rm = TRUE), 1),
    round(stats::median(per_sample$s_est_chao1,    na.rm = TRUE), 1),
    round(stats::median(per_sample$s_est_over_obs, na.rm = TRUE), 2),
    round(stats::median(per_sample$coverage_obs,   na.rm = TRUE), 3),
    sum(!is.na(per_sample$coverage_obs) & per_sample$coverage_obs < 0.9),
    round(100 * sum(!is.na(per_sample$coverage_obs) &
                      per_sample$coverage_obs < 0.9) /
            max(sum(!is.na(per_sample$coverage_obs)), 1L), 1),
    sum(!is.na(per_sample$coverage_obs) &
          per_sample$coverage_obs >= MAP_MIN_COVERAGE),
    round(stats::median(per_sample$ecm_reads,      na.rm = TRUE), 0)
  )
)

readr::write_csv(inext_summary, paths$linnean_inext_summary)
ts(sprintf("  Saved %s", basename(paths$linnean_inext_summary)))


ts("10_linnean_inext.R complete.")
