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

# ---- Figure 5 grey-background colour -----------------------------------------
# Figure 5 (manuscript.qmd) is a hand-assembled schematic composite that reuses
# panels from Figures 1, 3, 4, S1, S4, and S5 on a light-grey backdrop. The
# scripts producing those figures (10, 17, 18, and 20) each save TWO versions:
# the normal white-background copy used in the manuscript/SI (paths$fig_*),
# and an additional copy on this grey background (paths$fig_*_grey) for use
# only when manually assembling Figure 5.
fig5_grey_bg <- "#F2F2F2"

# ---- Multi-format figure export (journal submission) -------------------------
# save_fig_formats() wraps ggplot2::ggsave() to write the SAME figure in several
# file formats in one call. The pipeline names every figure with a .png path
# (paths$fig_*). FACETS (Canadian Science Publishing) does not accept PNG and
# instead asks for .tif / .jpg (among others), so each manuscript figure is
# written as PNG (unchanged) plus JPG and TIFF alongside it.
#
#   path    : canonical .png target (e.g. paths$fig_sampling_map). Its extension
#             is swapped to build the other formats in the same directory.
#   plot    : the ggplot / patchwork object to save.
#   formats : which formats to write. Default = PNG (keeps existing pipeline
#             behaviour) plus JPG and TIFF for submission.
#   ...     : forwarded verbatim to ggplot2::ggsave() (width, height, dpi, bg,
#             units, ...), so each call keeps the sizing it already specified.
#
# Format-specific defaults (overridable through ...):
#   JPG  -> quality = 100        (FACETS: max-quality JPEG, >= 300 dpi)
#   TIFF -> compression = "none" (FACETS: uncompressed TIFF)
save_fig_formats <- function(path, plot,
                             formats = c("png", "jpg", "tif"), ...) {
  stem <- tools::file_path_sans_ext(path)
  dots <- list(...)
  for (fmt in formats) {
    ext <- switch(fmt, jpeg = "jpg", tiff = "tif", fmt)           # normalize
    dev <- switch(ext, png = "png", jpg = "jpeg", tif = "tiff",
                  stop("save_fig_formats(): unsupported format '", fmt, "'"))
    args <- c(list(filename = paste0(stem, ".", ext),
                   plot = plot, device = dev), dots)
    if (ext == "jpg" && is.null(args$quality))     args$quality     <- 100
    if (ext == "tif" && is.null(args$compression)) args$compression <- "none"
    do.call(ggplot2::ggsave, args)
  }
  invisible(path)
}

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

# ---- read_big_tsv_subset(): read selected columns from a very large TSV -------
# data.table::fread(..., select = ) still has to stream and tokenize the whole
# file, which overruns R's 2^31-1-byte single-string limit on multi-GB inputs
# such as the ~13 GB GlobalFungi SH abundance matrix. This helper instead streams
# the file once through `awk`, writing ONLY the requested columns to a temporary
# TSV, then fread()s that much smaller file.
#   path : path to the large tab-separated file (with a header row)
#   cols : column NAMES to keep (first is typically the id column, e.g. sample_ID);
#          names not present in the file's header are dropped
# Returns a data.table with the kept columns, in the order given.
read_big_tsv_subset <- function(path, cols) {
  header <- names(data.table::fread(path, sep = "\t", quote = "", nrows = 0L))
  keep   <- cols[cols %in% header]
  if (length(keep) == 0L)
    stop("None of the requested columns are present in ", basename(path))
  idx <- match(keep, header)                          # 1-based column positions
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  awk_prog <- paste0(
    'BEGIN{FS=OFS="\\t"; n=split(cols,a,",")} ',
    '{out=$(a[1]); for(i=2;i<=n;i++) out=out OFS $(a[i]); print out}'
  )
  cmd <- sprintf("awk -v cols=%s %s %s > %s",
                 shQuote(paste(idx, collapse = ",")),
                 shQuote(awk_prog), shQuote(path), shQuote(tmp))
  if (system(cmd) != 0L || !file.exists(tmp))
    stop("awk column subset failed for ", basename(path))
  data.table::fread(tmp, sep = "\t", quote = "")
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
  # FungalTraits v1.2 (Põlme et al. 2020): pinned reference, as with UNITE above.
  # primary_lifestyle == "ectomycorrhizal" defines EcM dataset membership, so the
  # reference is held fixed rather than tracked to a live source. To re-pin,
  # replace the file, update data_raw/fungaltraits/fungaltraits_version.txt, and
  # re-validate the EcM genus count.
  fungaltraits   = here::here("data_raw", "fungaltraits", "FungalTraits_1-2.csv"),
  mycocosm_list  = here::here("data_raw", "mycocosm", "mycocosm_organism_list.csv"),
  climate_raster = here::here("data_raw", "climate", "wc2.1_country", "CAN_wc2.1_30s_bio.tif"),
  biotime_db     = here::here("data_raw", "biotime"),
  # van Galen et al. (2025) per-sample proportion of unassigned EcM OTUs.
  # Only this small CSV is needed (for the Linnean per-sample dark-fraction
  # comparison); the large dark-taxa GeoTIFFs are NOT part of this project.
  van_galen_per_sample = here::here("data_raw", "van_Galen_per_sample",
                                    "GFv5_EcM_unassigned_per_sample.csv"),
  # van Galen et al. (2025) dark-taxa richness GeoTIFF (band percentage_dark_taxa
  # is used for the dark-diversity map, Figure S1).
  van_galen_tif = here::here("data_raw", "van_Galen_et_al_dark_taxa_code_and_data",
                             "4.Dark_EcM_taxa_richness_maps",
                             "Dark_taxa_geospatial_layers.tif"),

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
  fungalroot_sp     = here::here("data_derived", "clean_fungalroot_species.csv"),
  fungalroot_genera = here::here("data_derived", "clean_fungalroot_genera_table_s2.csv"),
  host_species      = here::here("data_derived", "ecm_native_canada_host_species.csv"),

  # ---- derived: BIEN2 host rasters (from 07_host_rasters.R) -------------------
  bien_ranges        = here::here("data_derived", "checkpoints", "bien2_ecm_host_ranges.gpkg"),
  bien_richness      = here::here("data_derived", "spatial", "bien_host_richness_0.5deg.tif"),
  bien_species_stack = here::here("data_derived", "spatial", "bien_host_species_stack.tif"),
  bien_data_rich     = here::here("data_derived", "spatial", "bien_host_data_richness_0.5deg.tif"),
  bien_proportion    = here::here("data_derived", "spatial", "bien_host_data_proportion_0.5deg.tif"),

  # ---- derived: cross-referenced shortfall outputs ---------------------------
  gbif_ecm            = here::here("data_derived", "linnean", "linnean_gbif_ecm_canada.csv"),
  gbif_ecm_nosequence = here::here("data_derived", "linnean", "linnean_gbif_ecm_nosequence_canada.csv"),
  gbif_plotted_counts = here::here("data_derived", "linnean", "linnean_gbif_plotted_counts.csv"),
  eltonian_global        = here::here("data_derived", "eltonian", "eltonian_global_host_associations.csv"),
  eltonian_host_matching = here::here("data_derived", "eltonian", "eltonian_host_matching.csv"),
  prestonian_out      = here::here("data_derived", "prestonian", "prestonian_biotime_matches.csv"),
  darwinian_out       = here::here("data_derived", "darwinian", "darwinian_mycocosm_matches.csv"),

  # ---- figures (manuscript file names) ---------------------------------------
  # Main text
  fig_sampling_map      = here::here("figures", "Figure-01_sampling_map.png"),
  fig_sampling_map_grey = here::here("figures", "Figure-01_sampling_map_grey.png"),
  fig_wallacean_occ    = here::here("figures", "Figure-02_wallacean_occupancy.png"),
  fig_climate_gap      = here::here("figures", "Figure-03_climate_gap.png"),
  fig_climate_gap_grey = here::here("figures", "Figure-03_climate_gap_grey.png"),
  fig_host_bivariate      = here::here("figures", "Figure-04_host_bivariate_map.png"),
  fig_host_bivariate_grey = here::here("figures", "Figure-04_host_bivariate_map_grey.png"),
  fig_shortfalls_summary  = here::here("figures", "Figure-05_shortfalls_summary.png"),
  # Supplemental
  fig_dark_diversity      = here::here("figures", "Figure-S1_dark_diversity.png"),
  fig_dark_diversity_grey = here::here("figures", "Figure-S1_dark_diversity_grey.png"),
  fig_density_world    = here::here("figures", "Figure-S2_gf_sampling_density_world.png"),
  fig_depth_discard    = here::here("figures", "Figure-S3_depth_discard.png"),
  fig_ecozone_sampling      = here::here("figures", "Figure-S4_ecozone_sampling_map.png"),
  fig_ecozone_sampling_grey = here::here("figures", "Figure-S4_ecozone_sampling_map_grey.png"),
  fig_gbif_specimens      = here::here("figures", "Figure-S5_gbif_specimens.png"),
  fig_gbif_specimens_grey = here::here("figures", "Figure-S5_gbif_specimens_grey.png")
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
