# =============================================================================
# Prestonian Shortfall: Temporal Data Coverage
# =============================================================================
# Do any of the EcM fungal taxa recorded in our Canadian dataset appear in
# BioTIME — a global database of time-series biodiversity data? Matching
# our taxa to BioTIME quantifies the Prestonian shortfall: the degree to
# which temporal change in EcM fungal communities is undocumented.
#
# BioTIME is downloaded directly from:
#   RDS:      https://biotime.st-andrews.ac.uk/dl_request.php?dl=raw_rds
#   Metadata: https://biotime.st-andrews.ac.uk/dl_request.php?dl=metadata_csv
#
# The RDS contains a data frame with columns including:
#   STUDY_ID, ID_ALL_RAW_DATA, YEAR, taxon, genus, species, valid_name, ...
# We filter to Fungi records and match valid_name against our EcM species list.
# At genus level we also identify any BioTIME species from our EcM genera.
#
# For matched taxa we use STUDY_ID and study metadata to characterise whether
# those records represent true time-series (duration > 1 year), and what
# realm/biome they come from.
#
# Workflow:
#   1.  Load or download BioTIME RDS and metadata CSV
#   2.  Filter metadata to fungi studies (taxa + organisms both contain "fungi")
#   3.  Filter main BioTIME data to those study IDs
#   4.  Match valid_name against our EcM species (exact) and genera (genus-level)
#   5.  Spatial-temporal assessment: group by (valid_name, latitude, longitude),
#       count distinct years per location — n_years > 1 = genuine time-series
#   6.  Save summary and matched records
#
# Raw data files (data_raw/biotime/):
#   *.rds          — BioTIME raw data (scanned by extension; downloaded as
#                    biotime_raw.rds if absent); metadata CSV is REQUIRED
#   *.csv          — BioTIME study metadata (scanned by extension; downloaded
#                    as biotime_metadata.csv if absent)
#
# Outputs:
#   data_derived/prestonian_biotime_matches.csv      — BioTIME records for matched EcM taxa
#   data_derived/prestonian_taxon_summary.csv        — per-taxon summary: match level,
#                                                 n_records, year range, n_locations,
#                                                 max_years_at_location, has_repeat_sampling
#   data_derived/prestonian_location_timeseries.csv  — per-taxon x location x year detail
#   data_derived/prestonian_study_metadata.csv       — metadata for studies with EcM matches
#   data_derived/prestonian_summary.csv              — aggregate statistics
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(httr)

biotime_dir      <- here::here("data_raw", "biotime")
biotime_rds_path <- file.path(biotime_dir, "biotime_raw.rds")
biotime_csv_path <- file.path(biotime_dir, "biotime_metadata.csv")

BIOTIME_RDS_URL  <- "https://biotime.st-andrews.ac.uk/dl_request.php?dl=raw_rds"
BIOTIME_META_URL <- "https://biotime.st-andrews.ac.uk/dl_request.php?dl=metadata_csv"

dir.create(biotime_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Step 1: Load or download BioTIME RDS ------------------------------------
# Scan data_raw/biotime/ for any .rds file (catches user-downloaded files with
# names differing from the default, e.g. biotime_v2_data_raw_2025.rds).

existing_rds <- list.files(biotime_dir, pattern = "\\.rds$",
                            full.names = TRUE, ignore.case = TRUE)

if (length(existing_rds) > 0L) {
  biotime_data <- readRDS(existing_rds[1L])
} else {
  dl <- tryCatch(
    httr::GET(
      BIOTIME_RDS_URL,
      httr::write_disk(biotime_rds_path, overwrite = TRUE),
      httr::progress(),
      httr::config(timeout = 900)
    ),
    error = function(e) {
      warning("  Download failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dl) || httr::http_error(dl)) {
    stop(
      "Failed to download BioTIME RDS.\n",
      "Download manually from https://biotime.st-andrews.ac.uk\n",
      "and save to: data_raw/biotime/"
    )
  }
  biotime_data <- readRDS(biotime_rds_path)
}

# ---- Load or download BioTIME metadata CSV -----------------------------------

existing_csv <- list.files(biotime_dir, pattern = "\\.csv$",
                            full.names = TRUE, ignore.case = TRUE)

if (length(existing_csv) > 0L) {
  biotime_meta <- readr::read_csv(existing_csv[1L], show_col_types = FALSE)
} else {
  dl_meta <- tryCatch(
    httr::GET(
      BIOTIME_META_URL,
      httr::write_disk(biotime_csv_path, overwrite = TRUE),
      httr::config(timeout = 300)
    ),
    error = function(e) {
      warning("  Metadata download failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(dl_meta) && !httr::http_error(dl_meta)) {
    biotime_meta <- readr::read_csv(biotime_csv_path, show_col_types = FALSE)
  } else {
    biotime_meta <- NULL
  }
}

# Normalise column names to lower case
names(biotime_data) <- tolower(names(biotime_data))
if (!is.null(biotime_meta)) names(biotime_meta) <- tolower(names(biotime_meta))

# ---- Step 2: Identify fungi studies via metadata ----------------------------
# Filter metadata to studies where both 'taxa' and 'organisms' contain "fungi".
# This is the primary filter — more reliable than filtering the main data by a
# taxon column, and avoids loading all non-fungal records into the analysis.

if (is.null(biotime_meta)) {
  stop(
    "BioTIME metadata is required for this analysis.\n",
    "Download from https://biotime.st-andrews.ac.uk and save to data_raw/biotime/"
  )
}

fungi_studies <- biotime_meta |>
  dplyr::filter(
    grepl("fungi", taxa,      ignore.case = TRUE),
    grepl("fungi", organisms, ignore.case = TRUE)
  )

fungi_study_ids <- unique(fungi_studies$study_id)

# Retain useful metadata columns for later output
meta_cols <- intersect(
  c("study_id", "taxa", "organisms", "realm", "biome_map",
    "cent_lat", "cent_long", "start_year", "end_year", "duration",
    "number_of_species", "contact_1"),
  names(fungi_studies)
)
fungi_study_meta <- dplyr::select(fungi_studies, dplyr::all_of(meta_cols))

# ---- Step 3: Filter main BioTIME data to fungi studies ----------------------

bt_fungi <- dplyr::filter(biotime_data, study_id %in% fungi_study_ids)

# Build valid_name from genus + species if column is absent
if (!"valid_name" %in% names(bt_fungi)) {
  bt_fungi <- bt_fungi |>
    dplyr::mutate(valid_name = dplyr::if_else(
      !is.na(genus) & !is.na(species),
      paste(trimws(genus), trimws(species)),
      NA_character_
    ))
}

bt_fungi_species <- bt_fungi |>
  dplyr::filter(!is.na(valid_name), nzchar(valid_name)) |>
  dplyr::distinct(valid_name) |>
  dplyr::mutate(bt_genus = sub("^(\\S+).*", "\\1", valid_name))

# ---- Step 4: Match against our Canadian EcM taxa ----------------------------

# UNITE species uses "Genus_epithet" format; convert to "Genus epithet"
our_species <- emf |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::mutate(species_clean = trimws(gsub("_", " ", species))) |>
  dplyr::distinct(species_clean) |>
  dplyr::pull(species_clean)

our_genera <- unique(trimws(emf$genus))

# Exact species-level match
bt_ecm_sp <- dplyr::filter(bt_fungi_species, valid_name %in% our_species)
n_species_matched <- nrow(bt_ecm_sp)

# Genus-level match (catches species not named to our exact UNITE epithets)
bt_ecm_gn <- dplyr::filter(bt_fungi_species, bt_genus %in% our_genera)
n_genus_matched <- nrow(bt_ecm_gn)

# All matched valid_names (union of both match levels)
matched_names <- unique(c(bt_ecm_sp$valid_name, bt_ecm_gn$valid_name))

bt_ecm_records <- dplyr::filter(bt_fungi, valid_name %in% matched_names)
n_matched_records <- nrow(bt_ecm_records)

# ---- Step 5: Spatial-temporal assessment ------------------------------------
# For each matched taxon, determine whether there are repeat observations at
# the same location across multiple years. Group by (valid_name, latitude,
# longitude) and count distinct years — locations with n_years > 1 represent
# genuine time-series sampling.

if (n_matched_records > 0) {

  bt_ecm_records <- bt_ecm_records |>
    dplyr::mutate(
      match_level = dplyr::case_when(
        valid_name %in% bt_ecm_sp$valid_name ~ "species",
        TRUE                                 ~ "genus"
      )
    )

  # Per-location, per-taxon time series assessment
  location_timeseries <- bt_ecm_records |>
    dplyr::filter(!is.na(latitude), !is.na(longitude), !is.na(year)) |>
    dplyr::group_by(valid_name, match_level, study_id, latitude, longitude) |>
    dplyr::summarise(
      n_years         = dplyr::n_distinct(year),
      years           = paste(sort(unique(year)), collapse = "; "),
      n_records       = dplyr::n(),
      total_abundance = sum(abundance, na.rm = TRUE),
      total_biomass   = sum(biomass,   na.rm = TRUE),
      .groups         = "drop"
    ) |>
    dplyr::arrange(valid_name, dplyr::desc(n_years))

  n_locations_multiyear <- sum(location_timeseries$n_years > 1L, na.rm = TRUE)

  # Per-taxon summary across all locations
  taxon_summary <- bt_ecm_records |>
    dplyr::group_by(valid_name, match_level) |>
    dplyr::summarise(
      study_ids         = paste(sort(unique(study_id)), collapse = "; "),
      n_studies         = dplyr::n_distinct(study_id),
      n_records         = dplyr::n(),
      year_min          = min(year, na.rm = TRUE),
      year_max          = max(year, na.rm = TRUE),
      n_locations       = dplyr::n_distinct(
                            paste(latitude, longitude), na.rm = TRUE),
      max_years_at_location = max(
        location_timeseries$n_years[
          location_timeseries$valid_name == dplyr::cur_group()$valid_name
        ], na.rm = TRUE),
      has_repeat_sampling = max_years_at_location > 1L,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_records))

  readr::write_csv(taxon_summary,
                   file.path(paths$out_prestonian, "prestonian_taxon_summary.csv"))

  # Attach study metadata to the EcM-matched studies
  study_ids_with_ecm <- unique(bt_ecm_records$study_id)
  study_meta_ecm <- dplyr::filter(fungi_study_meta,
                                   study_id %in% study_ids_with_ecm)

  n_studies_with_ecm  <- length(study_ids_with_ecm)
  n_repeat_taxa       <- sum(taxon_summary$has_repeat_sampling, na.rm = TRUE)

  readr::write_csv(bt_ecm_records,
                   file.path(paths$out_prestonian, "prestonian_biotime_matches.csv"))

  readr::write_csv(location_timeseries,
                   file.path(paths$out_prestonian, "prestonian_location_timeseries.csv"))

  readr::write_csv(study_meta_ecm,
                   file.path(paths$out_prestonian, "prestonian_study_metadata.csv"))

} else {
  study_ids_with_ecm  <- character(0L)
  n_studies_with_ecm  <- 0L
  n_repeat_taxa       <- 0L
  n_locations_multiyear <- 0L
  readr::write_csv(dplyr::slice(bt_fungi, 0L),
                   file.path(paths$out_prestonian, "prestonian_biotime_matches.csv"))
}

# ---- Step 6: Save summary table ----------------------------------------------

our_genera_in_bt    <- our_genera[our_genera %in% bt_fungi_species$bt_genus]
our_genera_not_in_bt <- setdiff(our_genera, our_genera_in_bt)

# Named-taxon (species)-level absence: species-level exact match only (genus-
# level matches do not establish that the *named species* itself is present).
our_species_in_bt    <- unique(bt_ecm_sp$valid_name)
our_species_not_in_bt <- setdiff(our_species, our_species_in_bt)

prestonian_summary <- tibble::tibble(
  metric = c(
    "EcM genera in our Canadian dataset",
    "EcM species (UNITE, named) in our Canadian dataset",
    "BioTIME fungi studies (taxa + organisms filter)",
    "BioTIME species names in fungi studies",
    "BioTIME species matching our EcM species (exact)",
    "BioTIME species within our EcM genera (genus-level)",
    "BioTIME records for matched EcM taxa",
    "BioTIME studies containing matched EcM taxa",
    "Matched EcM taxa with repeat sampling (same location, >1 year)",
    "Location x taxon combinations with >1 year of data",
    "EcM genera with NO BioTIME record",
    "% of our EcM genera absent from BioTIME",
    "EcM named taxa (species) with NO BioTIME record",
    "% of our EcM named taxa absent from BioTIME"
  ),
  value = c(
    dplyr::n_distinct(emf$genus),
    length(our_species),
    nrow(fungi_studies),
    nrow(bt_fungi_species),
    n_species_matched,
    n_genus_matched,
    n_matched_records,
    n_studies_with_ecm,
    n_repeat_taxa,
    n_locations_multiyear,
    length(our_genera_not_in_bt),
    round(100 * length(our_genera_not_in_bt) / dplyr::n_distinct(emf$genus), 1),
    length(our_species_not_in_bt),
    round(100 * length(our_species_not_in_bt) / length(our_species), 1)
  )
)

readr::write_csv(prestonian_summary,
                 file.path(paths$out_prestonian, "prestonian_summary.csv"))
