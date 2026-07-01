# =============================================================================
# 00_setup.R  —  Shared setup for the ECM_manuscript pipeline
# =============================================================================
# WHAT THIS FILE IS
#   Every numbered script in scripts/ begins with:
#       source(here::here("scripts", "00_setup.R"))
#   Sourcing it (a) loads the handful of always-on packages, (b) defines the
#   coordinate reference systems, small helper functions, and the single
#   canonical `paths` list used everywhere, and (c) loads the primary EcM
#   dataset into the object `emf` if it has already been produced.
#
#   It writes nothing except creating empty output sub-directories under
#   data_derived/ the first time it runs.
#
# WHY A SINGLE SETUP FILE
#   Keeping paths, constants, and helpers in one place means an analysis script
#   never hard-codes a file path or re-implements a shared rule (e.g. how a
#   "site" is defined). To point the pipeline at data held elsewhere, edit the
#   `paths` list below and nothing else.
#
# CONVENTIONS USED THROUGHOUT THE PIPELINE
#   * File paths are always built with here::here() from the project root
#     (the folder containing ECM_manuscript.Rproj). No absolute paths.
#   * Package functions are called with explicit namespaces (dplyr::filter(),
#     sf::st_read(), ...). The five library() calls below are the only
#     exceptions, kept for readability of the very common verbs.
# =============================================================================

# ---- Always-on packages -----------------------------------------------------
library(here)     # here::here() builds paths from the project root
library(dplyr)    # data manipulation verbs
library(tidyr)    # pivoting / reshaping
library(readr)    # fast CSV read/write
library(ggplot2)  # figures

# ---- Timestamped console logging --------------------------------------------
# ts("message") prints "[HH:MM:SS] message" so long-running scripts leave a
# readable progress trail.
ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# ---- Coordinate reference systems -------------------------------------------
# crs_wgs84  : lat/lon, used for all coordinate-based operations.
# crs_albers : Canada Albers Equal Area Conic, used for mapping (equal-area so
#              areas and point densities are not distorted).
crs_wgs84  <- "EPSG:4326"
crs_albers <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# =============================================================================
# Helper functions
# =============================================================================

# ---- parse_gb_latlon(): parse a GenBank lat_lon string -----------------------
# GenBank stores coordinates as a single string like "45.1234 N 75.4321 W".
# Returns a one-row data.frame with numeric columns lat, lon (NA if unparseable).
parse_gb_latlon <- function(x) {
  result <- data.frame(lat = NA_real_, lon = NA_real_)
  if (is.na(x) || !nchar(trimws(x))) return(result)
  m <- regmatches(x, regexpr("([0-9.]+)\\s+([NS])\\s+([0-9.]+)\\s+([EW])", x))
  if (!length(m)) return(result)
  parts <- strsplit(trimws(m), "\\s+")[[1]]
  lat <- as.numeric(parts[1]) * ifelse(parts[2] == "S", -1, 1)
  lon <- as.numeric(parts[3]) * ifelse(parts[4] == "W", -1, 1)
  data.frame(lat = lat, lon = lon)
}

# ---- add_site_id(): define sampling "sites" by binning coordinates -----------
# A "site" is the bin produced by rounding lat/lon to 3 decimal places
# (~100 m at temperate latitudes). This collapses GPS-precision noise across
# repeat readings at one physical location without merging distinct sites.
# Adds columns site_lat, site_lon, and site (a "lat_lon" key).
# This is the canonical site rule; pass digits = 4 only if a finer bin is
# genuinely required.
add_site_id <- function(df, digits = 3L) {
  dplyr::mutate(df,
    site_lat = round(lat, digits),
    site_lon = round(lon, digits),
    site     = paste(site_lat, site_lon, sep = "_"))
}

# ---- snap_30arcsec(): snap a coordinate to the 30 arc-second (~1 km) grid -----
# 30 arc-seconds = 1/120 degree. Used by the Wallacean occupancy counts so that
# occurrences aggregate into ~1 km^2 grid cells (matching common SDM practice).
snap_30arcsec <- function(x) round(x * 120) / 120

# ---- build_site_sh_matrix(): long records -> binary site x SH matrix ----------
# Turns a long table (columns lat, lon, sh_code) into a presence/absence
# site x SH-code matrix for incidence-based richness estimation (iNEXT / Chao2).
# min_records drops sites with fewer than that many distinct SH detections;
# keep it at 1 for richness extrapolation, which needs singletons.
build_site_sh_matrix <- function(df, min_records = 1L) {
  df |>
    add_site_id() |>
    dplyr::group_by(site) |>
    dplyr::filter(dplyr::n() >= min_records) |>
    dplyr::ungroup() |>
    dplyr::distinct(site, sh_code) |>
    dplyr::mutate(present = 1L) |>
    tidyr::pivot_wider(names_from = sh_code, values_from = present,
                       values_fill = 0L) |>
    tibble::column_to_rownames("site") |>
    as.matrix()
}

# ---- canonicalize_host(): normalize a raw host-plant name --------------------
# Raw host strings come from the GenBank `host` qualifier and `isolation_source`
# field, and from GlobalFungi's `dominant_plant_species` / `other_plant_species`
# fields. This returns a clean "Genus species" (or "Genus") string, or NA when
# the input cannot be resolved to a Latin binomial. Steps, in order:
#   1. strip surrounding quotes/brackets and outer punctuation/whitespace
#   2. collapse internal whitespace
#   3. keep only the first two whitespace-separated tokens
#   4. strip residual per-token punctuation
#   5. apply the curated typo-correction lookup `typo_fixes`
#   6. case-normalize: Titlecase genus, lowercase species epithet (so that
#      "pinus banksiana" and "PINUS BANKSIANA" are both accepted)
#   7. reject anything that is not a valid Latin name pattern (common names such
#      as "pine", or entries with digits, become NA)
# Extend `typo_fixes` as new misspellings surface.
canonicalize_host <- function(x) {
  if (length(x) == 0L) return(character(0L))
  out <- trimws(as.character(x))
  out <- gsub('^["‘’“”\\(\\[]+', "", out, perl = TRUE)
  out <- gsub('["‘’“”\\)\\]]+$',  "", out, perl = TRUE)
  out <- gsub("^[[:punct:][:space:]]+|[[:punct:][:space:]]+$", "", out)
  out <- gsub("\\s+", " ", out)
  out <- trimws(out)
  empty_idx <- !is.na(out) & !nzchar(out)
  out[empty_idx] <- NA_character_

  toks <- strsplit(out, "\\s+")
  out <- vapply(seq_along(toks), function(i) {
    if (is.na(out[i])) return(NA_character_)
    p <- toks[[i]]
    if (length(p) == 0L) return(NA_character_)
    p <- gsub("^[[:punct:]]+|[[:punct:]]+$", "", p)
    p <- p[nzchar(p)]
    if (length(p) >= 2L) paste(p[1L], p[2L])
    else if (length(p) == 1L) p[1L]
    else NA_character_
  }, character(1L))

  typo_fixes <- c("Abies balamifera" = "Abies balsamea")
  ix <- match(out, names(typo_fixes))
  out[!is.na(ix)] <- unname(typo_fixes[ix[!is.na(ix)]])

  toks2 <- strsplit(out, " ", fixed = TRUE)
  out <- vapply(seq_along(toks2), function(i) {
    if (is.na(out[i])) return(NA_character_)
    p <- toks2[[i]]
    if (length(p) == 0L) return(NA_character_)
    g <- paste0(toupper(substring(p[1L], 1L, 1L)), tolower(substring(p[1L], 2L)))
    if (length(p) >= 2L) paste(g, tolower(p[2L])) else g
  }, character(1L))

  valid <- grepl("^[A-Z][a-z]+( [a-z][a-z-]*)?$", out)
  out[!valid] <- NA_character_
  out
}

# ---- load_bien2_range(): read one BIEN2 range polygon, clipped to Canada -----
# Reads the modelled range shapefile for `species` from data_raw/bien2_ranges/,
# repairs geometry, dissolves to a single polygon, reprojects to WGS84, and
# intersects with the Canada boundary. Returns an sf polygon, or NULL (with a
# warning) if no shapefile exists for that species.
load_bien2_range <- function(species, canada_wgs84) {
  sp_us  <- gsub(" ", "_", species)
  sp_dir <- file.path(paths$bien2_ranges_dir, sp_us)
  shp    <- list.files(sp_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp) == 0) {
    warning(sprintf("No BIEN2 shapefile found for %s", species))
    return(NULL)
  }
  sf::st_read(shp[[1]], quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_union() |>
    sf::st_transform("EPSG:4326") |>
    sf::st_intersection(suppressWarnings(sf::st_buffer(sf::st_make_valid(canada_wgs84), 0)))
}

# =============================================================================
# Canonical paths
# =============================================================================
# Two data roots:
#   data_raw/      read-only inputs (obtained from external sources; see
#                  data_raw/DATA-DICTIONARY.md). Never written to by the pipeline.
#   data_derived/  everything the pipeline produces. Starts empty; running the
#                  scripts in order recreates it (see data_derived/DATA-DICTIONARY.md).
# Figures go to figures/ with manuscript file names (Figure-xx_description.png).
paths <- list(

  # ---- directory roots -------------------------------------------------------
  data_raw     = here::here("data_raw"),
  data_derived = here::here("data_derived"),
  figures      = here::here("figures"),
  checkpoints  = here::here("data_derived", "checkpoints"),
  temp_dir     = here::here("data_derived", "checkpoints"),

  # ---- per-shortfall output directories --------------------------------------
  out_linnean       = here::here("data_derived", "linnean"),
  out_wallacean     = here::here("data_derived", "wallacean"),
  out_prestonian    = here::here("data_derived", "prestonian"),
  out_darwinian     = here::here("data_derived", "darwinian"),
  out_raunkiaeran   = here::here("data_derived", "raunkiaeran"),
  out_hutchinsonian = here::here("data_derived", "hutchinsonian"),
  out_eltonian      = here::here("data_derived", "eltonian"),

  # ---- raw: spatial ----------------------------------------------------------
  canada_bound_raw = here::here("data_raw", "admin_boundaries", "canada_bound_raw.gpkg"),
  usa_bound_raw    = here::here("data_raw", "admin_boundaries", "usa_bound_raw.gpkg"),
  mexico_bound_raw = here::here("data_raw", "admin_boundaries", "mexico_bound_raw.gpkg"),
  ecoregions_raw   = here::here("data_raw", "ecoregions", "Ecoregions", "ecoregions.shp"),
  ne_canada_raw    = here::here("data_raw", "natural_earth", "canada_ne.gpkg"),
  ne_lakes_raw     = here::here("data_raw", "natural_earth", "lakes_ne.gpkg"),

  # ---- raw: sequence sources -------------------------------------------------
  gf_metadata     = here::here("data_raw", "GlobalFungi", "GlobalFungi_5_sample_metadata.txt"),
  gf_sh_abundance = here::here("data_raw", "GlobalFungi", "GlobalFungi_5_SH_abundance_ITS1_ITS2.txt"),
  # UNITE reference FASTA — pinned to the 2024-04-04 build (see
  # data_raw/DATA-DICTIONARY.md for why this build, not the newest, is used).
  # GlobalFungi's pre-assigned SH codes must be decoded against this same build,
  # and GenBank SH codes are assigned in-house by vsearch against it, so both
  # roles deliberately point to one file.
  unite_fasta = here::here("data_raw", "UNITE", "sh_general_release_dynamic_04.04.2024_dev.fasta"),

  # ---- raw: trait / genome / climate reference -------------------------------
  fungaltraits   = here::here("data_raw", "fungaltraits", "FungalTraits_1-2.csv"),
  mycocosm_list  = here::here("data_raw", "mycocosm", "mycocosm_organism_list.csv"),
  climate_raster = here::here("data_raw", "climate", "wc2.1_country", "CAN_wc2.1_30s_bio.tif"),
  biotime_db     = here::here("data_raw", "biotime"),
  gsmc_metadata  = here::here("data_raw", "GSMc_02-09-2021",
                              "Tedersoo L, Mikryukov V, Anslan S et al. Fungi_GSMc_sample_metadata.txt"),
  # van Galen et al. (2025) per-sample proportion of unassigned EcM OTUs.
  # Only this small CSV is needed (for the Linnean per-sample dark-fraction
  # comparison); the large dark-taxa GeoTIFFs are NOT part of this project.
  van_galen_per_sample = here::here("data_raw", "van_Galen_per_sample",
                                    "GFv5_EcM_unassigned_per_sample.csv"),

  # ---- raw: BIEN2 modelled range shapefiles (from 06_bien2_ranges.R) ---------
  bien2_ranges_dir = here::here("data_raw", "bien2_ranges"),

  # ---- derived: processed spatial (from 01_spatial_data.R) -------------------
  canada_bound         = here::here("data_derived", "spatial", "canada_simple.gpkg"),
  canada_albers        = here::here("data_derived", "spatial", "canada_ne_albers.gpkg"),
  lakes_albers         = here::here("data_derived", "spatial", "lakes_canada_albers.gpkg"),
  ecoregions_processed = here::here("data_derived", "spatial", "ecoregions_processed.gpkg"),
  ecozone_names        = here::here("data_derived", "spatial", "ecozone_names.csv"),

  # ---- derived: sequence pipeline checkpoints & outputs ----------------------
  gf_ids_out       = here::here("data_derived", "checkpoints", "globalfungi_canada_ids.txt"),
  gf_sh_subset_out = here::here("data_derived", "checkpoints", "globalfungi_canada_SH_abundance.txt"),
  gf_meta_out      = here::here("data_derived", "checkpoints", "globalfungi_canada_metadata.csv"),
  unite_taxonomy   = here::here("data_derived", "checkpoints", "unite_sh_taxonomy.csv"),
  gf_sh_unmatched  = here::here("data_derived", "checkpoints", "globalfungi_sh_unmatched.csv"),
  gb_fasta_out     = here::here("data_derived", "checkpoints", "genbank_emf_canada.fasta"),
  gf_long_out      = here::here("data_derived", "globalfungi_canada_long.csv"),
  gb_long_out      = here::here("data_derived", "genbank_emf_canada_long.csv"),
  emf_combined     = here::here("data_derived", "emf_canada_combined.csv"),
  emf_data         = here::here("data_derived", "emf_canada_em_only.csv"),

  # ---- derived: host reference (from 05_fungalroot_hosts.R) ------------------
  fungalroot_sp = here::here("data_derived", "clean_fungalroot_species.csv"),
  host_species  = here::here("data_derived", "ecm_native_canada_host_species.csv"),

  # ---- derived: BIEN2 host rasters (from 07_host_rasters.R) -------------------
  bien_ranges        = here::here("data_derived", "checkpoints", "bien2_ecm_host_ranges.gpkg"),
  bien_richness      = here::here("data_derived", "spatial", "bien_host_richness_0.5deg.tif"),
  bien_species_stack = here::here("data_derived", "spatial", "bien_host_species_stack.tif"),
  bien_data_rich     = here::here("data_derived", "spatial", "bien_host_data_richness_0.5deg.tif"),
  bien_proportion    = here::here("data_derived", "spatial", "bien_host_data_proportion_0.5deg.tif"),
  bien_ecoregions    = here::here("data_derived", "spatial", "bien_ecoregions_with_host_habitat.gpkg"),

  # ---- derived: cross-referenced shortfall outputs ---------------------------
  gbif_ecm            = here::here("data_derived", "linnean", "linnean_gbif_ecm_canada.csv"),
  gbif_ecm_nosequence = here::here("data_derived", "linnean", "linnean_gbif_ecm_nosequence_canada.csv"),
  linnean_inext_rds        = here::here("data_derived", "checkpoints", "linnean_inext_per_sample.rds"),
  linnean_inext_per_sample = here::here("data_derived", "linnean", "linnean_inext_per_sample.csv"),
  linnean_inext_summary    = here::here("data_derived", "linnean", "linnean_inext_summary.csv"),
  eltonian_global        = here::here("data_derived", "eltonian", "eltonian_global_host_associations.csv"),
  eltonian_host_matching = here::here("data_derived", "eltonian", "eltonian_host_matching.csv"),
  prestonian_out      = here::here("data_derived", "prestonian", "prestonian_biotime_matches.csv"),
  darwinian_out       = here::here("data_derived", "darwinian", "darwinian_mycocosm_matches.csv"),

  # ---- figures (manuscript file names) ---------------------------------------
  # Main text
  fig_sampling_map     = here::here("figures", "Figure-01_sampling_map.png"),
  fig_wallacean_occ    = here::here("figures", "Figure-02_wallacean_occupancy.png"),
  fig_climate_gap      = here::here("figures", "Figure-03_climate_gap.png"),
  fig_host_bivariate   = here::here("figures", "Figure-04_host_bivariate_map.png"),
  # Supplemental
  fig_density_world    = here::here("figures", "Figure-S1_gf_sampling_density_world.png"),
  fig_ecozone_sampling = here::here("figures", "Figure-S2_ecozone_sampling_map.png"),
  fig_gbif_specimens   = here::here("figures", "Figure-S3_gbif_specimens.png")
)

# ---- Create empty data_derived/ output sub-directories the first time --------
for (d in c(paths$checkpoints, paths$figures, paths$biotime_db,
            paths$out_linnean, paths$out_wallacean, paths$out_prestonian,
            paths$out_darwinian, paths$out_raunkiaeran, paths$out_hutchinsonian,
            paths$out_eltonian, here::here("data_derived", "spatial"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Load the primary EcM dataset if it exists ------------------------------
# `emf` is the combined, EcM-filtered GlobalFungi + GenBank record table
# produced by 04_combine_ecm_dataset.R. Group 0 scripts that build it run before
# it exists, so loading is conditional. The na and col_types arguments below
# suppress spurious readr parsing warnings (literal "NULL" year strings from
# GenBank, and the rarely-populated secondary_lifestyle column).
if (file.exists(paths$emf_data)) {
  ts("Loading EcM dataset (GlobalFungi + GenBank)...")
  emf <- readr::read_csv(paths$emf_data, show_col_types = FALSE,
                         na = c("", "NA", "NULL"),
                         col_types = readr::cols(
                           secondary_lifestyle = readr::col_character(),
                           .default = readr::col_guess()
                         )) |>
    dplyr::filter(source %in% c("GenBank", "GlobalFungi"))
  ts(sprintf("  Records: %d  |  GlobalFungi: %d  |  GenBank: %d",
             nrow(emf), sum(emf$source == "GlobalFungi"), sum(emf$source == "GenBank")))
}
