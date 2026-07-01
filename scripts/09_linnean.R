# =============================================================================
# Linnean Shortfall (Updated): Taxonomic Richness
# =============================================================================
# How many EcM fungal taxa are recorded in Canada, at each taxonomic level?
# How does this compare to the global inventory of known EcM species and genera?
#
# Updates from original 01_linnean.R:
#   - iNEXT accumulation curves removed
#   - GBIF physical specimen records added (rgbif, checkpointed)
#   - Species-level assignment rates added:
#       GlobalFungi: proportion of SH codes with a species-level UNITE name
#       GenBank:     proportion of sequences with BLAST identity >= 97%
#   - Pooled-abundance Chao1 (Canada-wide GF) removed: read counts are a poor
#     proxy for abundance in metabarcoding data, and the estimator is already
#     complemented by the site-based Chao2 (Step 5) and the per-sample Chao1
#     (10_linnean_inext.R).
#
# Workflow:
#   1.  Taxonomic richness counts (SH codes, genera, species, lineages)
#   2.  GlobalFungi singleton SH codes (total abundance = 1)
#   3.  Compare observed genera to global EcM genus inventory (FungalTraits)
#   4.  Species-level assignment rates by source
#   5.  GBIF physical specimen records for Canadian EcM fungi
#
# Raw data downloaded:
#   data_raw/gbif/<key>.zip           — GBIF occurrence download zip (preserved/living
#                                      specimens, Canada, Kingdom Fungi)
#
# Checkpoint files (data_derived/temp/):
#   gbif_ecm_canada_raw.csv          — parsed GBIF records (pre-EcM filter)
#
# Outputs:
#   data_derived/linnean_summary.csv
#   data_derived/linnean_genus_coverage.csv
#   data_derived/linnean_gbif_ecm_canada.csv  — GBIF EcM specimen records
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(rgbif)

gbif_ckpt <- file.path(paths$temp_dir, "gbif_ecm_canada_raw.csv")

# ---- Step 1: Taxonomic richness counts ---------------------------------------

ts("Step 1: Counting unique taxa...")

n_sh      <- dplyr::n_distinct(emf$sh_code)
n_genus   <- dplyr::n_distinct(emf$genus)

# Named species only (UNITE species field; entries ending _sp are unresolved)
n_species <- dplyr::n_distinct(
  emf$species[!is.na(emf$species) & !grepl("_sp$", emf$species)]
)
n_lineage <- dplyr::n_distinct(
  emf$ectomycorrhiza_lineage[!is.na(emf$ectomycorrhiza_lineage)]
)

ts(sprintf("  Unique SH codes: %d", n_sh))
ts(sprintf("  Unique genera:   %d", n_genus))
ts(sprintf("  Unique species:  %d", n_species))
ts(sprintf("  Unique EcM lineages: %d", n_lineage))

# Coordinate-filtered equivalents (records validated within the Canada boundary)
emf_coords <- dplyr::filter(emf, coord_in_canada == TRUE)

n_sh_coords          <- dplyr::n_distinct(emf_coords$sh_code)
n_named_sh_coords    <- emf_coords |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(sh_code) |>
  dplyr::n_distinct()
n_genus_coords       <- dplyr::n_distinct(emf_coords$genus)
n_species_coords     <- dplyr::n_distinct(
  emf_coords$species[!is.na(emf_coords$species) & !grepl("_sp$", emf_coords$species)]
)

ts(sprintf("  Unique SH codes (with coords):         %d", n_sh_coords))
ts(sprintf("  Named-species SH codes (with coords):  %d", n_named_sh_coords))
ts(sprintf("  Unique genera (with coords):           %d", n_genus_coords))
ts(sprintf("  Unique species (with coords):          %d", n_species_coords))

# ---- Step 2: GlobalFungi singleton SH codes ----------------------------------

ts("Step 2: Identifying GlobalFungi singleton SH codes (total abundance = 1)...")

gf_all <- dplyr::filter(emf, source == "GlobalFungi")

sh_abundance <- gf_all |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(total_abundance = sum(abundance, na.rm = TRUE),
                   .groups = "drop")

singleton_sh     <- dplyr::filter(sh_abundance, total_abundance == 1L)$sh_code
n_sh_singletons  <- length(singleton_sh)
n_sh_gf_total    <- dplyr::n_distinct(gf_all$sh_code)
n_sh_nonsing     <- n_sh_gf_total - n_sh_singletons
pct_singletons   <- round(100 * n_sh_singletons / n_sh_gf_total, 1)

ts(sprintf("  GlobalFungi SH codes: %d total | %d singletons (%.1f%%) | %d non-singleton",
           n_sh_gf_total, n_sh_singletons, pct_singletons, n_sh_nonsing))

# ---- Step 3: Compare observed genera to global EcM inventory ----------------

ts("Step 3: Loading FungalTraits for global EcM genus inventory...")

ft <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
  dplyr::rename_with(tolower)

global_ecm_genera <- ft |>
  dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
  dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
  dplyr::distinct(genus_lower, .keep_all = TRUE)

n_global_ecm_genera <- nrow(global_ecm_genera)
ts(sprintf("  Known EcM genera globally (FungalTraits): %d", n_global_ecm_genera))

our_genera_lower <- tolower(trimws(unique(emf$genus)))

genus_coverage <- global_ecm_genera |>
  dplyr::mutate(observed_in_canada = genus_lower %in% our_genera_lower) |>
  dplyr::select(genus, genus_lower, observed_in_canada,
                ectomycorrhiza_lineage_template,
                ectomycorrhiza_exploration_type_template)

n_observed_genera   <- sum(genus_coverage$observed_in_canada)
pct_observed_genera <- round(100 * n_observed_genera / n_global_ecm_genera, 1)

ts(sprintf("  EcM genera observed in Canada: %d (%.1f%% of global EcM genera)",
           n_observed_genera, pct_observed_genera))

our_genera_not_in_ft <- our_genera_lower[
  !our_genera_lower %in% global_ecm_genera$genus_lower
]
ts(sprintf("  Our genera absent from FungalTraits EcM list: %d",
           length(our_genera_not_in_ft)))
if (length(our_genera_not_in_ft) > 0) {
  ts(sprintf("    Examples: %s", paste(head(our_genera_not_in_ft, 5), collapse = ", ")))
}

readr::write_csv(genus_coverage,
                 file.path(paths$out_linnean, "linnean_genus_coverage.csv"))
ts("  Saved linnean_genus_coverage.csv")

# ---- Step 4: Species-level assignment rates by source -----------------------

ts("Step 4: Calculating species-level assignment rates by source...")

# --- 4a. GlobalFungi: proportion of SH codes with species-level UNITE name ---
# A SH code has a species-level name when the 'species' field does NOT end in
# '_sp' (which indicates UNITE could not assign a species epithet).

gf_species_rate <- gf_all |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(
    has_species = any(!is.na(species) & !grepl("_sp$", species)),
    .groups = "drop"
  )

n_gf_sh_with_sp  <- sum(gf_species_rate$has_species)
n_gf_sh_total    <- nrow(gf_species_rate)
pct_gf_sp        <- round(100 * n_gf_sh_with_sp / n_gf_sh_total, 1)

ts(sprintf("  GlobalFungi SH codes with species-level name: %d / %d (%.1f%%)",
           n_gf_sh_with_sp, n_gf_sh_total, pct_gf_sp))

# --- 4b. GenBank: proportion of sequences with BLAST identity >= 97% ----------
# The 'identity' column is the pairwise % identity from the BLAST search
# against UNITE. A threshold of 97% is a commonly used species-level cutoff.

IDENTITY_THRESHOLD <- 97

gb_all <- dplyr::filter(emf, source == "GenBank")
n_gb_total     <- nrow(gb_all)
n_gb_above_97  <- sum(!is.na(gb_all$identity) & gb_all$identity >= IDENTITY_THRESHOLD)
pct_gb_above_97 <- round(100 * n_gb_above_97 / n_gb_total, 1)

ts(sprintf("  GenBank sequences with identity >= %d%%: %d / %d (%.1f%%)",
           IDENTITY_THRESHOLD, n_gb_above_97, n_gb_total, pct_gb_above_97))

# Distribution of identity values (summary)
if (any(!is.na(gb_all$identity))) {
  id_summary <- summary(gb_all$identity[!is.na(gb_all$identity)])
  ts(sprintf("  GenBank identity range: %.1f – %.1f%% (median %.1f%%)",
             id_summary["Min."], id_summary["Max."], id_summary["Median"]))
}

# GenBank SH codes with species-level UNITE name (analogous to GlobalFungi metric)
gb_sh_total      <- dplyr::n_distinct(gb_all$sh_code)
gb_sh_with_sp    <- gb_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  nrow()
pct_gb_sh_sp <- round(100 * gb_sh_with_sp / gb_sh_total, 1)

ts(sprintf("  GenBank SH codes (all): %d", gb_sh_total))
ts(sprintf("  GenBank SH codes with species-level UNITE name: %d / %d (%.1f%%)",
           gb_sh_with_sp, gb_sh_total, pct_gb_sh_sp))

# --- 4b.5 Overlap of ALL SH codes (named + unnamed) between GF and GenBank ----
# Used to compute the dark fraction for source-exclusive and shared subsets
# (Table S4 in the SI). gf_all and gb_all are already defined above.
gf_all_sh <- unique(gf_all$sh_code)
gb_all_sh <- unique(gb_all$sh_code)
n_sh_all_shared  <- length(intersect(gf_all_sh, gb_all_sh))
n_sh_all_gf_only <- length(setdiff(gf_all_sh, gb_all_sh))
n_sh_all_gb_only <- length(setdiff(gb_all_sh, gf_all_sh))

ts(sprintf("  All SH codes: shared: %d | GF only: %d | GB only: %d",
           n_sh_all_shared, n_sh_all_gf_only, n_sh_all_gb_only))

# --- 4c. Overlap of named-species SH codes between GlobalFungi and GenBank ----
gf_sp_sh <- gf_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(sh_code) |>
  unique()

gb_sp_sh <- gb_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  dplyr::pull(sh_code) |>
  unique()

n_sp_sh_shared  <- length(intersect(gf_sp_sh, gb_sp_sh))
n_sp_sh_gf_only <- length(setdiff(gf_sp_sh, gb_sp_sh))
n_sp_sh_gb_only <- length(setdiff(gb_sp_sh, gf_sp_sh))
n_sp_sh_total   <- n_sp_sh_shared + n_sp_sh_gf_only + n_sp_sh_gb_only

ts(sprintf("  Named-species SH codes: %d total | shared: %d | GF only: %d | GB only: %d",
           n_sp_sh_total, n_sp_sh_shared, n_sp_sh_gf_only, n_sp_sh_gb_only))

# --- 4d. Overlap of unique named species (epithet-level) between GF and GenBank ----
# Unique named species are the distinct UNITE species epithets carried by the
# named-species SH codes above. Multiple SH codes can share an epithet (see
# Table S2 caption), so these counts are lower than and not simply derivable
# from the SH-level counts in 4c.
gf_sp_names <- gf_all |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(species) |>
  unique()

gb_sp_names <- gb_all |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  dplyr::pull(species) |>
  unique()

n_sp_shared   <- length(intersect(gf_sp_names, gb_sp_names))
n_sp_gf_only  <- length(setdiff(gf_sp_names, gb_sp_names))
n_sp_gb_only  <- length(setdiff(gb_sp_names, gf_sp_names))
n_sp_total    <- n_sp_shared + n_sp_gf_only + n_sp_gb_only

ts(sprintf("  Unique named species: %d total | shared: %d | GF only: %d | GB only: %d",
           n_sp_total, n_sp_shared, n_sp_gf_only, n_sp_gb_only))

# =============================================================================
# Step 5: Coverage diagnostics, asymptotic richness estimators, and triangulation
# =============================================================================
# This section addresses the bottom-up dimension of the Linnean shortfall:
# given the SH codes detected in the Canadian sequence dataset, how close is
# the sample to saturation, and what is a defensible lower bound on the
# unobserved richness?
#
# Design choices (see project notes):
#   1. SH codes are the primary unit (consistent with the rest of the pipeline)
#   2. GlobalFungi is the primary source for Chao-style estimators because GF
#      has a more uniform sampling design than GenBank; GenBank is opportunistic
#      and host-targeted, which violates the iid assumption Chao requires.
#      The combined GF + GenBank dataset is run as a sensitivity comparison.
#   3. Stratification is by Canadian ecozone where ≥30 GF samples are available
#      (the n threshold below which Chao becomes unreliable). Ecozones below
#      threshold are reported with Sobs and n only.
#   4. Coverage diagnostics (sample completeness Ĉ, Q1/Q2, accumulation curves)
#      are the main outcomes. Asymptotic point estimates are reported
#      alongside as supporting context, framed explicitly as lower bounds.
#   5. Singletons are retained (not excluded as a sensitivity check): the
#      dataset's singletons are very likely real detections rather than
#      sequencing artefacts.
#   6. Top-down anchors triangulate the bottom-up estimate against (a) a
#      genus-multiplier extrapolation from FungalTraits, and (b) the van Galen
#      et al. (2025) continental dark-taxa estimate for North America (80%).
#
# References:
#   Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation:
#     standardizing samples by completeness rather than size. Ecology 93:
#     2533–2547.
#   Colwell, R. K. et al. (2012) Models and estimators linking individual-
#     based and sample-based rarefaction, extrapolation and comparison of
#     assemblages. Journal of Plant Ecology 5: 3–21.
# =============================================================================

ts("Step 5: Coverage diagnostics, asymptotic estimators, and triangulation...")

library(sf)
library(terra)
library(iNEXT)
library(patchwork)

ECOZONE_N_THRESHOLD   <- 30L  # Minimum unique sites per ecozone for rarefaction/extrapolation curves
ECOZONE_COV_THRESHOLD <- 10L  # Minimum unique sites per ecozone for coverage (Ĉ) diagnostics only
sf::sf_use_s2(FALSE)

# ---- 5a. Build site × SH incidence matrices ---------------------------------
# Sampling units throughout this section are "sites" as defined by the shared
# add_site_id() helper in 00_setup.R: lat/lon rounded to 3 decimal places
# (~100 m at temperate latitudes) and concatenated into a stable composite
# key. This matches the project convention used by build_site_sh_matrix() in
# the Pinus banksiana accumulation curve (20_sampling_maps.R) and the distance-
# decay turnover analysis. Collapsing samples to one row per site removes
# within-site autocorrelation and GPS-precision artefacts, so coverage Ĉ and
# Chao2 reflect spatial sampling completeness rather than within-site read
# replication. Trade-off: discards within-site detection-frequency information
# and reduces n, so per-ecozone estimators are noisier — accepted as the more
# transparent framing.

ts("  5a. Building site × SH incidence matrices (3-decimal site binning)...")

# Filter to records with non-missing SH code and parsed coordinates
emf_sh  <- dplyr::filter(emf, !is.na(sh_code))
gf_for_sites <- emf_sh |>
  dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE)
combined_for_sites <- emf_sh |>
  dplyr::filter(coord_in_canada == TRUE)

# Build site × SH incidence matrices via the shared helper. min_records = 1L
# preserves single-record sites (singletons are needed for Chao2).
inc_gf       <- build_site_sh_matrix(gf_for_sites,       min_records = 1L)
inc_combined <- build_site_sh_matrix(combined_for_sites, min_records = 1L)
storage.mode(inc_gf)       <- "integer"
storage.mode(inc_combined) <- "integer"
ts(sprintf("    GF site matrix:       %d sites × %d SHs",
           nrow(inc_gf), ncol(inc_gf)))
ts(sprintf("    Combined site matrix: %d sites × %d SHs",
           nrow(inc_combined), ncol(inc_combined)))

# Site → coordinate lookup (used by the spatial join below). Built from the
# combined-source data so it covers every site that appears in either matrix.
site_coords <- combined_for_sites |>
  add_site_id() |>
  dplyr::distinct(site, site_lat, site_lon)

# ---- 5b. Spatial join: unique sites to ecozones -----------------------------
ts("  5b. Joining unique sites to Canadian ecozones...")

ecoregions_raw     <- sf::st_read(paths$ecoregions_processed, quiet = TRUE)
ecozone_names_tbl  <- readr::read_csv(paths$ecozone_names, show_col_types = FALSE)
ecoregions_named   <- dplyr::left_join(ecoregions_raw, ecozone_names_tbl, by = "ECOZONE") |>
  dplyr::mutate(NAME_EN = dplyr::if_else(is.na(NAME_EN),
                                          paste("Ecozone", ECOZONE),
                                          NAME_EN))

# One spatial join per unique site (rounded coordinate); downstream per-ecozone
# subsetting filters the GF site matrix by ecozone label.
site_pts <- sf::st_as_sf(site_coords, coords = c("site_lon", "site_lat"),
                         crs = 4326) |>
  sf::st_transform(sf::st_crs(ecoregions_named))

site_ecozone <- sf::st_join(site_pts, ecoregions_named, join = sf::st_within) |>
  sf::st_drop_geometry() |>
  dplyr::select(site, ecozone = NAME_EN) |>
  dplyr::distinct()

ts(sprintf("    Sites mapped to ecozone: %d / %d (unmapped: %d)",
           sum(!is.na(site_ecozone$ecozone)),
           nrow(site_ecozone),
           sum(is.na(site_ecozone$ecozone))))

# ---- 5c. Helper functions for richness estimation ----------------------------

# Convert a site × SH binary incidence matrix to an iNEXT incidence-frequency
# vector: first element = number of sites (T), remaining elements = per-SH
# detection frequency across sites (Q_i values).
# unname() is essential: c(nrow(M), colSums(M)) produces a named vector whose
# first element has a blank name; some iNEXT versions misparse this. as.numeric()
# ensures a plain double vector regardless of storage.mode of M.
to_inext_freq <- function(M) unname(as.numeric(c(nrow(M), colSums(M))))

# Run iNEXT on a single stratum and return a list with:
#   $summary  — one-row tibble of diagnostics + Chao2 asymptotic estimate
#   $curve    — tibble of the size-based rarefaction/extrapolation curve
#               columns: sites, richness, richness_lci, richness_uci, coverage
# When nrow(M) < n_threshold the estimator and curve fields are NULL/NA.
# nboot = 50 is sufficient for SE/CI estimation; raise for publication figures.
run_inext_stratum <- function(M, label, n_threshold = ECOZONE_N_THRESHOLD,
                              nboot = 50L) {
  n_sites  <- nrow(M)
  n_sh_obs <- ncol(M)
  Q1       <- sum(colSums(M) == 1L)
  Q2       <- sum(colSums(M) == 2L)

  # coverage_hat: iNEXT::DataInfo() returns SC at observed sample size via a
  # fast deterministic call (no bootstrapping). Called for all strata so that
  # below-threshold ecozones also get a coverage value.
  di      <- tryCatch(
    iNEXT::DataInfo(to_inext_freq(M), datatype = "incidence_freq"),
    error = function(e) NULL
  )
  cov_hat <- if (!is.null(di) && "SC" %in% names(di)) di$SC[1L] else NA_real_

  base <- tibble::tibble(
    stratum      = label,
    n_sites      = n_sites,
    n_sh_obs     = n_sh_obs,
    Q1           = Q1,
    Q2           = Q2,
    coverage_hat = cov_hat,
    chao2        = NA_real_,
    chao2_se     = NA_real_,
    chao2_lci    = NA_real_,
    chao2_uci    = NA_real_
  )

  if (n_sites < n_threshold || n_sh_obs == 0L) {
    return(list(summary = base, curve = NULL))
  }

  res <- tryCatch(
    suppressWarnings(
      iNEXT::iNEXT(to_inext_freq(M), q = 0,
                   datatype = "incidence_freq",
                   nboot    = nboot,
                   conf     = 0.95)
    ),
    error = function(e) {
      warning(sprintf("run_inext_stratum: iNEXT failed for '%s': %s",
                      label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) return(list(summary = base, curve = NULL))

  # ---- Asymptotic estimate (Chao2 at q = 0) ----------------------------------
  asy <- res$AsyEst
  if (!is.null(asy) && nrow(asy) > 0L) {
    pick <- rep(TRUE, nrow(asy))
    if ("Diversity" %in% names(asy)) {
      pick <- grepl("species\\s*richness", as.character(asy$Diversity),
                    ignore.case = TRUE)
    } else if ("Order.q" %in% names(asy)) {
      pick <- asy$Order.q == 0L
    }
    a0 <- if (any(pick)) asy[pick, , drop = FALSE] else asy[1L, , drop = FALSE]

    pick_col <- function(df, cands) {
      hit <- cands[cands %in% names(df)][1L]
      if (is.na(hit)) NA_real_ else df[[hit]][1L]
    }
    base$chao2     <- pick_col(a0, c("Estimator", "S.est"))
    base$chao2_se  <- pick_col(a0, c("Est_s.e.", "s.e.", "SE", "Std.Error"))
    base$chao2_lci <- pick_col(a0, c("95% Lower", "LCL", "qD.LCL"))
    base$chao2_uci <- pick_col(a0, c("95% Upper", "UCL", "qD.UCL"))
  }

  # ---- Size-based rarefaction / extrapolation curve --------------------------
  sb <- res$iNextEst$size_based
  if (is.null(sb) && is.data.frame(res$iNextEst)) sb <- res$iNextEst
  if (is.null(sb) && !is.null(res$iNextEst[[1L]]))  sb <- res$iNextEst[[1L]]

  curve <- NULL
  if (!is.null(sb) && nrow(sb) > 0L) {
    if ("Order.q" %in% names(sb)) sb <- sb[sb$Order.q == 0L, , drop = FALSE]
    t_col   <- intersect(c("t", "m", "x"),        names(sb))[1L]
    qd_col  <- intersect(c("qD", "Richness"),      names(sb))[1L]
    lci_col <- intersect(c("qD.LCL", "LCL"),       names(sb))[1L]
    uci_col <- intersect(c("qD.UCL", "UCL"),       names(sb))[1L]
    sc_col  <- intersect(c("SC", "Coverage"),      names(sb))[1L]
    met_col <- intersect(c("Method", "method"),    names(sb))[1L]
    curve <- tibble::tibble(
      sites        = if (!is.na(t_col))   sb[[t_col]]   else seq_len(nrow(sb)),
      richness     = if (!is.na(qd_col))  sb[[qd_col]]  else NA_real_,
      richness_lci = if (!is.na(lci_col)) sb[[lci_col]] else NA_real_,
      richness_uci = if (!is.na(uci_col)) sb[[uci_col]] else NA_real_,
      coverage     = if (!is.na(sc_col))  sb[[sc_col]]  else NA_real_,
      method       = if (!is.na(met_col)) sb[[met_col]] else NA_character_
    )
  }

  list(summary = base, curve = curve)
}

# ---- 5d. Per-ecozone and Canada-wide estimators -----------------------------

ts("  5c. Computing iNEXT incidence-based coverage diagnostics and asymptotic estimators...")
ts("       (datatype = 'incidence_freq', q = 0, nboot = 50)")

# Build per-ecozone GF site-incidence sub-matrices
ecozone_levels <- sort(unique(stats::na.omit(site_ecozone$ecozone)))
gf_site_ids    <- rownames(inc_gf)
ecozone_mats   <- lapply(ecozone_levels, function(ez) {
  ids <- site_ecozone$site[
    !is.na(site_ecozone$ecozone) & site_ecozone$ecozone == ez
  ]
  ids <- intersect(ids, gf_site_ids)
  if (length(ids) == 0L) return(NULL)
  m <- inc_gf[ids, , drop = FALSE]
  m[, colSums(m) > 0L, drop = FALSE]
})
names(ecozone_mats) <- ecozone_levels
ecozone_mats <- ecozone_mats[!vapply(ecozone_mats, is.null, logical(1L))]

# Compute iNEXT results for every stratum (singletons included).
# Accumulation curves are extracted here and stored for the RDS (5e below).
# Seed set once, immediately before the first bootstrap call, so the Chao2
# CIs (nboot = 50, via iNEXT::iNEXT()) are reproducible across runs. Without
# this, S_Chao2 itself is deterministic but Chao2_95_LCI/UCI are not.
set.seed(3492)
ts("    Canada (GF)...")
res_canada_gf       <- run_inext_stratum(inc_gf,       "Canada (GF)")
ts("    Canada (GF + GenBank)...")
res_canada_combined <- run_inext_stratum(inc_combined, "Canada (GF + GenBank)")
ts("    Per-ecozone strata...")
res_ecozones <- lapply(names(ecozone_mats), function(ez) {
  ts(sprintf("      %s...", ez))
  run_inext_stratum(ecozone_mats[[ez]], ez)
})
names(res_ecozones) <- names(ecozone_mats)

# Bind summary rows → estimators_full
# (singletons are retained in the data: the dataset's singletons are very
# likely real detections rather than artefacts, so no exclusion sensitivity
# run is performed. The `singletons` column is kept at a constant value of
# "included" for backward compatibility with existing downstream filters.)
estimators_full <- dplyr::bind_rows(
  res_canada_gf$summary,
  res_canada_combined$summary,
  dplyr::bind_rows(lapply(res_ecozones, `[[`, "summary"))
) |> dplyr::mutate(singletons = "included")

estimators <- estimators_full

readr::write_csv(estimators,
                 file.path(paths$out_linnean, "linnean_extrapolation_estimators.csv"))
ts(sprintf("    Saved linnean_extrapolation_estimators.csv (%d rows)", nrow(estimators)))

# Coverage-only diagnostic table (per stratum, GF only)
coverage_tbl <- estimators_full |>
  dplyr::select(stratum, n_sites, n_sh_obs, Q1, Q2, coverage_hat) |>
  dplyr::arrange(dplyr::desc(n_sites))
readr::write_csv(coverage_tbl,
                 file.path(paths$out_linnean, "linnean_extrapolation_coverage.csv"))
ts("    Saved linnean_extrapolation_coverage.csv")

# ---- 5e. Rarefaction / extrapolation curves (iNEXT size-based) ---------------
# Curves were already computed as part of run_inext_stratum() above; here we
# simply collect, save, and plot them. The RDS stores per-stratum tibbles with
# columns: sites, richness, richness_lci, richness_uci, coverage, method.
# 'method' is the iNEXT Method field: "Rarefaction", "Observed", "Extrapolation".

ts("  5d. Saving iNEXT rarefaction/extrapolation curves...")

# Collect viable ecozone curves (strata with ≥ ECOZONE_N_THRESHOLD sites)
acc_ecozone_curves <- lapply(names(res_ecozones), function(ez) {
  r <- res_ecozones[[ez]]
  if (is.null(r$curve)) return(NULL)
  if ((r$summary$n_sites) < ECOZONE_N_THRESHOLD) return(NULL)
  r$curve
})
names(acc_ecozone_curves) <- names(res_ecozones)
acc_ecozone_curves <- acc_ecozone_curves[!vapply(acc_ecozone_curves, is.null, logical(1L))]

saveRDS(
  list(canada_gf       = res_canada_gf$curve,
       canada_combined = res_canada_combined$curve,
       ecozone         = acc_ecozone_curves),
  file.path(paths$out_linnean, "linnean_accumulation.rds")
)
ts("    Saved linnean_accumulation.rds (iNEXT size-based curves; tibble format)")

# The accumulation curves are saved above as linnean_accumulation.rds. The
# corresponding diagnostic panel figure is not part of the manuscript, so it is
# not written to figures/.

# Compute the viable-ecozone count here (before freeing ecozone_mats below)
n_viable_ecozones <- sum(
  vapply(ecozone_mats, function(m) nrow(m) >= ECOZONE_N_THRESHOLD, logical(1L))
)

# Free large objects before continuing to GBIF step
rm(inc_gf, inc_combined, ecozone_mats, ecoregions_named,
   site_pts, site_ecozone, site_coords)

# ---- Step 6: GBIF physical specimen records ----------------------------------
# Query GBIF for preserved/living fungal specimens in Canada, then filter to
# EcM genera using our observed genus list. Checkpointed to avoid re-querying.
#
# Prerequisites: GBIF credentials in environment variables:
#   gbif_user, gbif_pwd, gbif_email
# Set these with usethis::edit_r_environ() or Sys.setenv().

if (file.exists(gbif_ckpt)) {
  ts("Step 5: Loading checkpointed GBIF specimen records...")
  gbif_raw <- readr::read_csv(
    gbif_ckpt,
    show_col_types = FALSE,
    col_types = readr::cols(
      decimalLatitude  = readr::col_double(),
      decimalLongitude = readr::col_double(),
      .default         = readr::col_character()
    )
  )
  ts(sprintf("  Records loaded from checkpoint: %d", nrow(gbif_raw)))
} else {

  gbif_user  <- Sys.getenv("GBIF_USER")
  gbif_pwd   <- Sys.getenv("GBIF_PWD")
  gbif_email <- Sys.getenv("GBIF_EMAIL")

  if (any(c(gbif_user, gbif_pwd, gbif_email) == "")) {
    warning(
      "GBIF credentials not found in environment variables ",
      "(GBIF_USER, GBIF_PWD, GBIF_EMAIL).\n",
      "Set these with usethis::edit_r_environ() and restart R.\n",
      "Skipping GBIF step."
    )
    gbif_raw <- NULL
  } else {

    ts(sprintf(
      "Step 5: Submitting GBIF download request (Canada, Fungi, physical specimens)..."
    ))

    gbif_dl <- tryCatch(
      rgbif::occ_download(
        rgbif::pred("country",               "CA"),
        rgbif::pred("HAS_GEOSPATIAL_ISSUE",  FALSE),
        rgbif::pred("taxonKey",              5),    # Kingdom Fungi
        rgbif::pred_in("BASIS_OF_RECORD",
                       c("PRESERVED_SPECIMEN", "LIVING_SPECIMEN")),
        format = "SIMPLE_CSV",
        user   = gbif_user,
        pwd    = gbif_pwd,
        email  = gbif_email
      ),
      error = function(e) {
        message("GBIF download submission failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gbif_dl)) {
      ts("  Waiting for GBIF download to complete (may take several minutes)...")
      rgbif::occ_download_wait(gbif_dl)

      ts("  Importing GBIF download...")
      gbif_zip_dir <- here::here("data_raw", "gbif")
      dir.create(gbif_zip_dir, showWarnings = FALSE, recursive = TRUE)
      gbif_zip <- rgbif::occ_download_get(gbif_dl,
                                           path      = gbif_zip_dir,
                                           overwrite = TRUE)
      gbif_raw <- rgbif::occ_download_import(gbif_zip) |>
        as.data.frame()

      ts(sprintf("  Raw GBIF records: %d", nrow(gbif_raw)))
      readr::write_csv(gbif_raw, gbif_ckpt)
      ts(sprintf("  Saved GBIF checkpoint -> %s", basename(gbif_ckpt)))

      ts("  GBIF citation:")
      print(rgbif::gbif_citation(gbif_zip))
    } else {
      gbif_raw <- NULL
    }
  }
}

# Filter GBIF records to EcM genera and species-level records
if (!is.null(gbif_raw) && nrow(gbif_raw) > 0) {

  ts("  Filtering GBIF records to species rank and EcM genera...")

  # Use our observed EcM genera as the filter list
  ecm_genera_lower <- tolower(trimws(unique(emf$genus)))

  # GBIF 'genus' column should be available in SIMPLE_CSV
  gbif_sp <- dplyr::filter(gbif_raw, taxonRank == "SPECIES")

  gbif_ecm <- gbif_sp |>
    dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
    dplyr::filter(genus_lower %in% ecm_genera_lower)

  n_gbif_ecm         <- nrow(gbif_ecm)
  n_gbif_ecm_species <- dplyr::n_distinct(gbif_ecm$species)
  n_gbif_ecm_genera  <- dplyr::n_distinct(gbif_ecm$genus)

  ts(sprintf("  GBIF EcM specimen records: %d", n_gbif_ecm))
  ts(sprintf("  GBIF EcM species: %d  |  genera: %d",
             n_gbif_ecm_species, n_gbif_ecm_genera))

  readr::write_csv(gbif_ecm, paths$gbif_ecm)
  ts(sprintf("  Saved -> %s", basename(paths$gbif_ecm)))

  # Filter to EcM genera with NO sequence data in Canada
  # Uses the FungalTraits EcM genus list from Step 3 (global_ecm_genera),
  # selecting genera flagged as not observed in our sequence dataset.
  ecm_genera_nosequence <- genus_coverage$genus_lower[
    !genus_coverage$observed_in_canada
  ]

  gbif_ecm_nosequence <- gbif_sp |>
    dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
    dplyr::filter(genus_lower %in% ecm_genera_nosequence)

  n_gbif_ecm_nosequence         <- nrow(gbif_ecm_nosequence)
  n_gbif_ecm_nosequence_species <- dplyr::n_distinct(gbif_ecm_nosequence$species)
  n_gbif_ecm_nosequence_genera  <- dplyr::n_distinct(gbif_ecm_nosequence$genus)

  ts(sprintf("  GBIF EcM (genera without sequence data) records: %d",
             n_gbif_ecm_nosequence))
  ts(sprintf("  GBIF EcM (no sequence) species: %d  |  genera: %d",
             n_gbif_ecm_nosequence_species, n_gbif_ecm_nosequence_genera))

  readr::write_csv(gbif_ecm_nosequence, paths$gbif_ecm_nosequence)
  ts(sprintf("  Saved -> %s", basename(paths$gbif_ecm_nosequence)))

} else {
  ts("  GBIF data not available — skipping EcM filter.")
  n_gbif_ecm                    <- NA_integer_
  n_gbif_ecm_species            <- NA_integer_
  n_gbif_ecm_genera             <- NA_integer_
  n_gbif_ecm_nosequence         <- NA_integer_
  n_gbif_ecm_nosequence_species <- NA_integer_
  n_gbif_ecm_nosequence_genera  <- NA_integer_
}

# ---- Step 7: Save summary table ----------------------------------------------

# Pull a small set of headline numbers from the new Step 5 outputs to include
# in the summary, so the master tibble alongside `linnean_summary.csv` reflects
# the bottom-up extrapolation results without duplicating the per-stratum
# detail (which lives in linnean_extrapolation_*.csv).
canada_full      <- estimators_full[estimators_full$stratum == "Canada (GF)", ]
chao_canada      <- canada_full$chao2
chao_canada_se   <- canada_full$chao2_se
chao_canada_lci  <- canada_full$chao2_lci
chao_canada_uci  <- canada_full$chao2_uci
cov_canada_gf    <- canada_full$coverage_hat
# n_viable_ecozones was computed in Step 5 before the rm() of ecozone_mats

linnean_summary <- tibble::tibble(
  metric = c(
    "Unique UNITE v10 SH codes (combined dataset, all records)",
    "Unique UNITE v10 SH codes (records with coordinates only)",
    "Named-species SH codes (records with coordinates only)",
    "Unique EcM genera (combined dataset, FungalTraits-filtered)",
    "Unique EcM genera (records with coordinates only)",
    "Unique named species (UNITE taxonomy, excl. _sp; combined dataset)",
    "Unique named species (records with coordinates only)",
    "All SH codes: shared between GlobalFungi and GenBank",
    "All SH codes: GlobalFungi only (not in GenBank)",
    "All SH codes: GenBank only (not in GlobalFungi)",
    "Named-species SH codes: total unique across GF + GenBank (regardless of coords)",
    "Named-species SH codes: shared between GlobalFungi and GenBank",
    "Named-species SH codes: GlobalFungi only",
    "Named-species SH codes: GenBank only",
    "Unique named species: shared between GlobalFungi and GenBank",
    "Unique named species: GlobalFungi only",
    "Unique named species: GenBank only",
    "Unique EcM lineages observed in Canada (lineage definitions from FungalTraits)",
    "Known EcM genera globally (FungalTraits)",
    "EcM genera observed in Canada",
    "% of global EcM genera observed in Canada",
    "Canadian EcM genera absent from FungalTraits EcM list",
    "GlobalFungi: unique SH codes (Canadian dataset)",
    "GlobalFungi: singleton SH codes (total abundance = 1 across Canadian samples)",
    "GlobalFungi: singleton SH codes (% of GF Canadian SHs)",
    "GlobalFungi: non-singleton SH codes",
    "GlobalFungi: SH codes with species-level UNITE name",
    "GlobalFungi: SH codes with species-level name (% of GF Canadian SHs)",
    sprintf("GenBank: sequences with BLAST identity >= %d%% (sequence count, not SH count)", IDENTITY_THRESHOLD),
    sprintf("GenBank: sequences with BLAST identity >= %d%% (%% of GenBank Canadian sequences)", IDENTITY_THRESHOLD),
    "GenBank: unique SH codes (Canadian dataset)",
    "GenBank: SH codes with species-level UNITE name",
    "GenBank: SH codes with species-level name (% of GenBank Canadian SHs)",
    "Site coverage Ĉ (GlobalFungi, Canada-wide; 3-decimal lat/lon binning)",
    "Chao2 SH richness, Canada-wide GF site × SH matrix (lower bound)",
    "Chao2 SH richness bootstrap SE, Canada-wide GF site × SH matrix",
    "Chao2 SH richness bootstrap 95% LCI, Canada-wide GF site × SH matrix",
    "Chao2 SH richness bootstrap 95% UCI, Canada-wide GF site × SH matrix",
    sprintf("Viable ecozones for coverage diagnostics (≥ %d unique GF sites)", ECOZONE_COV_THRESHOLD),
    "GBIF: physical EcM specimen records — genera WITH sequence data (Canada)",
    "GBIF: EcM species represented (genera with sequence data)",
    "GBIF: EcM genera represented (genera with sequence data)",
    "GBIF: physical EcM specimen records — genera WITHOUT sequence data (Canada)",
    "GBIF: EcM species represented (genera without sequence data)",
    "GBIF: EcM genera represented (genera without sequence data)"
  ),
  value = c(
    n_sh, n_sh_coords, n_named_sh_coords, n_genus, n_genus_coords, n_species, n_species_coords,
    n_sh_all_shared, n_sh_all_gf_only, n_sh_all_gb_only,
    n_sp_sh_total, n_sp_sh_shared, n_sp_sh_gf_only, n_sp_sh_gb_only,
    n_sp_shared, n_sp_gf_only, n_sp_gb_only,
    n_lineage,
    n_global_ecm_genera, n_observed_genera, pct_observed_genera,
    length(our_genera_not_in_ft),
    n_sh_gf_total, n_sh_singletons, pct_singletons, n_sh_nonsing,
    n_gf_sh_with_sp, pct_gf_sp,
    n_gb_above_97, pct_gb_above_97,
    gb_sh_total, gb_sh_with_sp, pct_gb_sh_sp,
    round(cov_canada_gf, 3),
    round(chao_canada, 0),
    round(chao_canada_se, 1),
    round(chao_canada_lci, 0),
    round(chao_canada_uci, 0),
    n_viable_ecozones,
    n_gbif_ecm, n_gbif_ecm_species, n_gbif_ecm_genera,
    n_gbif_ecm_nosequence, n_gbif_ecm_nosequence_species, n_gbif_ecm_nosequence_genera
  )
)

readr::write_csv(linnean_summary,
                 file.path(paths$out_linnean, "linnean_summary.csv"))
ts("Saved linnean_summary.csv")
ts("09_linnean.R complete.")
