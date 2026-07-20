# =============================================================================
# Linnean Shortfall: Taxonomic Richness
# =============================================================================
# How many EcM fungal taxa are recorded in Canada, at each taxonomic level?
# How does this compare to the global inventory of known EcM species and genera?
#
# Species-level assignment rates are reported per source:
#   GlobalFungi: proportion of SH codes with a species-level UNITE name
#   GenBank:     proportion of SH codes with a species-level UNITE name
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

n_sh      <- dplyr::n_distinct(emf$sh_code, na.rm = TRUE)  # na.rm: genus-resolved GenBank rows carry sh_code = NA
n_genus   <- dplyr::n_distinct(emf$genus)

# Named species only (UNITE species field; entries ending _sp are unresolved)
n_species <- dplyr::n_distinct(
  emf$species[!is.na(emf$species) & !grepl("_sp$", emf$species)]
)
n_lineage <- dplyr::n_distinct(
  emf$ectomycorrhiza_lineage[!is.na(emf$ectomycorrhiza_lineage)]
)

# Coordinate-filtered equivalents (records validated within the Canada boundary)
emf_coords <- dplyr::filter(emf, coord_in_canada == TRUE)

n_sh_coords          <- dplyr::n_distinct(emf_coords$sh_code, na.rm = TRUE)
n_named_sh_coords    <- emf_coords |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(sh_code) |>
  dplyr::n_distinct()
n_genus_coords       <- dplyr::n_distinct(emf_coords$genus)
n_species_coords     <- dplyr::n_distinct(
  emf_coords$species[!is.na(emf_coords$species) & !grepl("_sp$", emf_coords$species)]
)

# ---- Step 2: GlobalFungi singleton SH codes ----------------------------------

gf_all <- dplyr::filter(emf, source == "GlobalFungi")

sh_abundance <- gf_all |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(total_abundance = sum(abundance, na.rm = TRUE),
                   .groups = "drop")

singleton_sh     <- dplyr::filter(sh_abundance, total_abundance == 1L)$sh_code
n_sh_singletons  <- length(singleton_sh)
n_sh_gf_total    <- dplyr::n_distinct(gf_all$sh_code, na.rm = TRUE)
n_sh_nonsing     <- n_sh_gf_total - n_sh_singletons
pct_singletons   <- round(100 * n_sh_singletons / n_sh_gf_total, 1)

# ---- Step 3: Compare observed genera to global EcM inventory ----------------

ft <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
  dplyr::rename_with(tolower)

global_ecm_genera <- ft |>
  dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
  dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
  dplyr::distinct(genus_lower, .keep_all = TRUE)

n_global_ecm_genera <- nrow(global_ecm_genera)

our_genera_lower <- tolower(trimws(unique(emf$genus)))

genus_coverage <- global_ecm_genera |>
  dplyr::mutate(observed_in_canada = genus_lower %in% our_genera_lower) |>
  dplyr::select(genus, genus_lower, observed_in_canada,
                ectomycorrhiza_lineage_template,
                ectomycorrhiza_exploration_type_template)

n_observed_genera   <- sum(genus_coverage$observed_in_canada)
pct_observed_genera <- round(100 * n_observed_genera / n_global_ecm_genera, 1)

our_genera_not_in_ft <- our_genera_lower[
  !our_genera_lower %in% global_ecm_genera$genus_lower
]

readr::write_csv(genus_coverage,
                 file.path(paths$out_linnean, "linnean_genus_coverage.csv"))

# ---- Step 4: Species-level assignment rates by source -----------------------

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

# --- 4b. GenBank: total EcM sequence records ---------------------------------
# Every retained GenBank record cleared the 98.5% vsearch identity threshold by
# construction (see 03_genbank.R), so no additional identity filter is applied
# here — we simply report the total record count. The `identity` column
# (vsearch alignment identity) is summarised below for reference only.

gb_all <- dplyr::filter(emf, source == "GenBank")
n_gb_total <- nrow(gb_all)

# Distribution of identity values (summary)

# GenBank SH codes with species-level UNITE name (analogous to GlobalFungi metric)
gb_sh_total      <- dplyr::n_distinct(gb_all$sh_code, na.rm = TRUE)
gb_sh_with_sp    <- gb_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  nrow()
pct_gb_sh_sp <- round(100 * gb_sh_with_sp / gb_sh_total, 1)

# --- 4b.5 Overlap of ALL SH codes (named + unnamed) between GF and GenBank ----
# Used to compute the dark fraction for source-exclusive and shared subsets
# (Table S4 in the SI). gf_all and gb_all are already defined above.
# Exclude NA sh_code (genus-resolved GenBank rows) so it is not counted as a
# spurious "GenBank-only" SH in the set operations below.
gf_all_sh <- unique(gf_all$sh_code[!is.na(gf_all$sh_code)])
gb_all_sh <- unique(gb_all$sh_code[!is.na(gb_all$sh_code)])
n_sh_all_shared  <- length(intersect(gf_all_sh, gb_all_sh))
n_sh_all_gf_only <- length(setdiff(gf_all_sh, gb_all_sh))
n_sh_all_gb_only <- length(setdiff(gb_all_sh, gf_all_sh))

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

# ---- Step 6: GBIF physical specimen records ----------------------------------
# Query GBIF for preserved/living fungal specimens in Canada, then filter to
# EcM genera using our observed genus list. Checkpointed to avoid re-querying.
#
# Prerequisites: GBIF credentials in environment variables:
#   gbif_user, gbif_pwd, gbif_email
# Set these with usethis::edit_r_environ() or Sys.setenv().

if (file.exists(gbif_ckpt)) {
  gbif_raw <- readr::read_csv(
    gbif_ckpt,
    show_col_types = FALSE,
    col_types = readr::cols(
      decimalLatitude  = readr::col_double(),
      decimalLongitude = readr::col_double(),
      .default         = readr::col_character()
    )
  )

} else if (length(list.files(here::here("data_raw", "gbif"), pattern = "\\.zip$")) > 0) {

  # A GBIF download ZIP is already provided in data_raw/gbif/. Import it directly
  # instead of submitting a fresh (slow, credential-gated) download.
  # occ_download_import(key, path) reads a local "<key>.zip" from `path`, where
  # the key is the file name without its extension.
  gbif_zip_dir  <- here::here("data_raw", "gbif")
  gbif_zip_file <- list.files(gbif_zip_dir, pattern = "\\.zip$", full.names = TRUE)[1]
  gbif_key      <- sub("\\.zip$", "", basename(gbif_zip_file))
  gbif_raw <- rgbif::occ_download_import(key = gbif_key, path = gbif_zip_dir) |>
    as.data.frame()
  readr::write_csv(gbif_raw, gbif_ckpt)

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
        warning("GBIF download submission failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gbif_dl)) {
      rgbif::occ_download_wait(gbif_dl)

      gbif_zip_dir <- here::here("data_raw", "gbif")
      dir.create(gbif_zip_dir, showWarnings = FALSE, recursive = TRUE)
      gbif_zip <- rgbif::occ_download_get(gbif_dl,
                                           path      = gbif_zip_dir,
                                           overwrite = TRUE)
      gbif_raw <- rgbif::occ_download_import(gbif_zip) |>
        as.data.frame()

      readr::write_csv(gbif_raw, gbif_ckpt)

    } else {
      gbif_raw <- NULL
    }
  }
}

# Filter GBIF records to EcM genera and species-level records
if (!is.null(gbif_raw) && nrow(gbif_raw) > 0) {

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

  readr::write_csv(gbif_ecm, paths$gbif_ecm)

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

  readr::write_csv(gbif_ecm_nosequence, paths$gbif_ecm_nosequence)

} else {
  n_gbif_ecm                    <- NA_integer_
  n_gbif_ecm_species            <- NA_integer_
  n_gbif_ecm_genera             <- NA_integer_
  n_gbif_ecm_nosequence         <- NA_integer_
  n_gbif_ecm_nosequence_species <- NA_integer_
  n_gbif_ecm_nosequence_genera  <- NA_integer_
}

# ---- Step 7: Save summary table ----------------------------------------------

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
    "GenBank: total EcM sequence records",
    "GenBank: unique SH codes (Canadian dataset)",
    "GenBank: SH codes with species-level UNITE name",
    "GenBank: SH codes with species-level name (% of GenBank Canadian SHs)",
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
    n_gb_total,
    gb_sh_total, gb_sh_with_sp, pct_gb_sh_sp,
    n_gbif_ecm, n_gbif_ecm_species, n_gbif_ecm_genera,
    n_gbif_ecm_nosequence, n_gbif_ecm_nosequence_species, n_gbif_ecm_nosequence_genera
  )
)

readr::write_csv(linnean_summary,
                 file.path(paths$out_linnean, "linnean_summary.csv"))
