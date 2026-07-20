# =============================================================================
# Prepare FungalRoot Data: Download, Parse, Classify, Standardize
# =============================================================================
# Produces a computationally reproducible table of EcM host plant species
# "demonstrated" by the FungalRoot database, used to determine which plant
# genera (and species) count as EcM hosts for the Eltonian and Hutchinsonian
# shortfall analyses.
#
# Sources (two independent sources):
#   1. GBIF DwC-A archive of the FungalRoot occurrence data (species rule)
#     GBIF dataset key: 744edc21-8dd2-474e-8a0b-b8c3d56a3c2d
#     Download URL pattern:
#       https://api.gbif.org/v1/occurrence/download/request/<KEY>.zip
#     OR download interactively from:
#       https://www.gbif.org/occurrence/download?dataset_key=744edc21-...
#   2. FungalRoot's own published genus-level recommendation, Supplementary
#      Table S2 (Soudzilovskaia et al. 2020 New Phytologist, supporting
#      information), sheet "Table S2" — genus rule. Scanned by file
#      extension from data_raw/fungalroot/ (one .xlsx expected).
#
# DwC-A structure:
#   occurrences.csv: core file; column 0 = ID (primary key), column 44 = scientificName
#   measurements.csv: extension; coreid links to occurrences.csv ID;
#                     columns: coreid, measurementType, measurementValue
#
# Table S2 structure:
#   Row 1: title; row 2: blank; row 3: header ("Genus", "Mycorrhizal type");
#   row 4 onward: one row per genus. Read with `skip = 2`.
#
# Host-determination rule (two independent routes, combined with OR):
#
#   1. SPECIES RULE. A species is
#      "ecm_demonstrated" if at least one raw FungalRoot occurrence record
#      assigns it one of the following unambiguous EcM-positive
#      MycorrhizalType labels:
#        "EcM", "EcM, no AM colonization", "ErM,EcM"
#      Two labels are explicitly EXCLUDED as ambiguous:
#        - "EcM, AM undetermined" — does not distinguish a genuine EcM record
#          from a non-EcM record where AM status simply wasn't assessed (e.g.
#          ferns such as Dryopteris filix-mas, whose only EcM-flavoured record
#          carries this label and is contradicted by separate AM /
#          non-mycorrhizal records for the same species).
#        - "EcM,AM" (dual-positive call) — asserts both types were found in
#          the same record. Admitting it as qualifying evidence is a route to
#          false positives: a substantial set of genera would qualify ONLY via
#          this label (including Acer, Allium, Fraxinus, Juglans, Rubus and
#          Sambucus), with no unambiguous EcM evidence anywhere in the genus.
#
#   2. GENUS RULE. A genus counts as an EcM host genus if FungalRoot's own
#      published Table S2 recommendation (Soudzilovskaia et al. 2020 supporting
#      information; genus-level calls made at >=67% consistency of species
#      diagnosis) assigns it "EcM" or "EcM-AM" (`GENUS_EM_QUALIFYING_TYPES`).
#      "EcM-AM" captures genuine dual-mycorrhizal genera (e.g. Salix, Populus,
#      Eucalyptus) that the binary species rule cannot express.
#
#   3. COMBINATION. The two routes are independent and combined with OR: a
#      species/genus qualifies as an EcM host if it satisfies the species
#      rule directly, OR its genus qualifies under the Table S2 genus rule
#      (even if that particular species has no occurrence-level record of its
#      own). Note that Bistorta/Persicaria are called NM-AM in Table S2, so
#      Bistorta vivipara qualifies only via its own occurrence record — hence
#      the routes are combined with OR rather than the genus rule replacing
#      the species rule. This script does NOT perform the OR merge itself: it
#      writes two separate, clearly labelled outputs (see "Output" below), each
#      tagged with its `evidence_source`. 06_host_species.R performs the merge.
#
# Table S2 data quality:
#   A small number of Table S2 rows carry messy free-text or apparent
#   data-entry-error MycorrhizalTypeGenus values (e.g. "species-specific: EcM
#   or NM-AM", a value of "Thysanothus" for genus Thysanotus, "nM" for genus
#   Winklera). None of these match `GENUS_EM_QUALIFYING_TYPES` and so are
#   excluded by default (see Step 3b); none affect the qualifying genus count
#   either way.
#
# Name standardization:
#   Uses rgbif::name_backbone_checklist() to resolve taxonomic synonyms (this
#   is what lets "Persicaria vivipara" resolve to the same accepted name as
#   "Bistorta vivipara"). Only species with at least one EcM-positive raw
#   record are sent to the backbone lookup (species with no positive record
#   can never affect the species rule, so resolving their synonymy is
#   unnecessary and would needlessly inflate the number of API calls).
#   Batched at 200 names per request. Table S2 genera are NOT sent through
#   this backbone lookup (see Step 3b) — Table S2 is genus-level only, and
#   genus-level synonymy is not resolved by this script.
#
# Raw data files (data_raw/fungalroot/):
#   744edc21-...zip — GBIF DwC-A archive (scanned by prefix; downloaded via
#                     rgbif::occ_download() if absent)
#   *.xlsx          — FungalRoot Supplementary Tables S1-S4 workbook (scanned
#                     by extension; must contain exactly one .xlsx file);
#                     "Table S2" sheet is used here.
#
# Checkpoints (data_derived/checkpoints/):
#   fungalroot_dwca_extracted/        — extracted DwC-A directory
#   fungalroot_species_raw.csv        — parsed (PlantBinomial, MycorrhizalType)
#                                        before classification, one row per
#                                        distinct raw occurrence/measurement pair
#   fungalroot_species_backbone.csv   — rgbif backbone results for the
#                                        EcM-demonstrated candidate species
#   fungalroot_table_s2.csv           — parsed Table S2 (Genus,
#                                        MycorrhizalTypeGenus), pre-filter
#
# Output:
#   data_derived/clean_fungalroot_species.csv  (species rule)
#     UpdatedPlantBinomial — backbone-resolved accepted name
#     UpdatedGenus         — genus parsed from UpdatedPlantBinomial; retained
#                            for transparency, not used for genus-level
#                            qualification (that comes from Table S2)
#     ecm_demonstrated     — always TRUE in this table (only demonstrated
#                            species are retained; kept as an explicit column
#                            for self-documentation and so downstream code
#                            does not have to infer it)
#     evidence             — the qualifying raw MycorrhizalType label(s)
#                            backing the call, semicolon-separated
#     PlantBinomial         — the raw FungalRoot name(s) (pre-backbone),
#                            semicolon-separated when more than one raw
#                            synonym resolves to the same accepted name
#     evidence_source      — always "occurrence" in this table; distinguishes
#                            this route from the Table S2 genus route below
#
#   data_derived/clean_fungalroot_genera_table_s2.csv  (genus rule)
#     Genus                 — genus name as given in Table S2
#     mycorrhizal_type      — Table S2's MycorrhizalTypeGenus call ("EcM" or
#                            "EcM-AM" only; this table is pre-filtered to
#                            qualifying genera)
#     evidence_source      — always "table_s2"
#
#   NOTE — this script writes TWO output tables. A Canadian native species can
#   qualify as an EcM host via either route (its own row in the species table,
#   OR its genus appearing in the genus table) without itself appearing in the
#   species table. 06_host_species.R reads both and performs the OR merge.
#
# Runtime notes:
#   - rgbif backbone queries: a few minutes, scales with the number of
#     EcM-demonstrated candidate species (a small fraction of all ~17,000
#     raw FungalRoot occurrence/measurement rows).
#   - Table S2 read (Step 3b): near-instant; no API calls involved.
#   - Delete checkpoint files to force re-run of individual steps.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(rgbif)
library(readxl)

# ---- Paths -------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
fr_dir   <- here("data_raw", "fungalroot")
out_dir  <- here("data_derived")
temp_dir <- paths$temp_dir  # alias for data_derived/checkpoints (see 00_setup.R)

DWCA_ZIP_PATTERN <- "744edc21"  # prefix to locate existing zip

OUT_SPECIES <- file.path(out_dir, "clean_fungalroot_species.csv")
OUT_GENERA  <- file.path(out_dir, "clean_fungalroot_genera_table_s2.csv")

ckpt_species_raw      <- file.path(temp_dir, "fungalroot_species_raw.csv")
ckpt_species_backbone <- file.path(temp_dir, "fungalroot_species_backbone.csv")
ckpt_table_s2         <- file.path(temp_dir, "fungalroot_table_s2.csv")

# The three unambiguous EcM-positive raw MycorrhizalType labels (see header
# comment above for the rationale, including why "EcM, AM undetermined" is
# deliberately excluded).
EM_POSITIVE_LABELS <- c("EcM", "EcM, no AM colonization", "ErM,EcM")

# Table S2 genus-level calls that qualify a genus as an EcM host genus (see
# header comment above for rationale; excludes "uncertain" and the handful of
# messy free-text values, which are logged but not treated as qualifying).
GENUS_EM_QUALIFYING_TYPES <- c("EcM", "EcM-AM")

# ---- Helper: rgbif backbone lookup, batched ----------------------------------

backbone_standardize <- function(names_vec, batch_size = 200L) {
  names_vec <- as.character(names_vec)
  n_batches <- ceiling(length(names_vec) / batch_size)
  results   <- vector("list", n_batches)
  for (i in seq_len(n_batches)) {
    idx <- seq((i - 1L) * batch_size + 1L,
               min(i * batch_size, length(names_vec)))
    results[[i]] <- tryCatch(
      rgbif::name_backbone_checklist(name_data = data.frame(scientificName = names_vec[idx])),
      error = function(e) {
        warning("  Batch ", i, " failed: ", conditionMessage(e))
        data.frame(verbatim_name = names_vec[idx],
                   canonicalName  = NA_character_,
                   matchType      = "NONE",
                   stringsAsFactors = FALSE)
      }
    )
    if (i < n_batches) Sys.sleep(2)
  }
  dplyr::bind_rows(results)
}

# Extract best updated name from backbone result
# Prefers canonicalName for ACCEPTED/SYNONYM; falls back to verbatim_name
extract_updated_name <- function(backbone_df) {
  dplyr::mutate(
    backbone_df,
    UpdatedName = dplyr::case_when(
      matchType != "NONE" & !is.na(canonicalName) & nzchar(canonicalName) ~ canonicalName,
      TRUE ~ verbatim_name
    )
  ) |>
    dplyr::select(verbatim_name, UpdatedName)
}

# =============================================================================
# Step 1: Locate or download DwC-A zip
# =============================================================================

existing_zips <- list.files(fr_dir, pattern = DWCA_ZIP_PATTERN, full.names = TRUE)
existing_zip  <- existing_zips[grepl("\\.zip$", existing_zips)]

if (length(existing_zip) == 0L) {

  gbif_user  <- Sys.getenv("GBIF_USER")
  gbif_pwd   <- Sys.getenv("GBIF_PWD")
  gbif_email <- Sys.getenv("GBIF_EMAIL")

  if (!nzchar(gbif_user) || !nzchar(gbif_pwd) || !nzchar(gbif_email)) {
    stop(
      "GBIF credentials not found in environment.\n",
      "Add to .Renviron:\n",
      "  GBIF_USER=your_username\n",
      "  GBIF_PWD=your_password\n",
      "  GBIF_EMAIL=your_email\n",
      "Then re-run, or download manually from:\n",
      "  https://www.gbif.org/occurrence/download?dataset_key=744edc21-8dd2-474e-8a0b-b8c3d56a3c2d\n",
      "and save the zip to: ", fr_dir
    )
  }

  dl_key <- rgbif::occ_download(
    rgbif::pred("datasetKey", "744edc21-8dd2-474e-8a0b-b8c3d56a3c2d"),
    format = "DWCA",
    user   = gbif_user,
    pwd    = gbif_pwd,
    email  = gbif_email
  )
  rgbif::occ_download_wait(dl_key, status_ping = 30L)

  zip_dest <- file.path(fr_dir, paste0(dl_key, ".zip"))
  rgbif::occ_download_get(dl_key, path = fr_dir, overwrite = TRUE)
  existing_zip <- zip_dest
} else {
  existing_zip <- existing_zip[1L]
}

# =============================================================================
# Step 2: Extract DwC-A zip
# =============================================================================

extract_dir <- file.path(temp_dir, "fungalroot_dwca_extracted")

# Identify which files are needed (avoid re-extracting if present)
need_extract <- !all(file.exists(
  file.path(extract_dir, c("occurrences.csv", "measurements.csv"))
))

if (need_extract) {
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  unzip(existing_zip, exdir = extract_dir, overwrite = TRUE)
}

occ_file  <- file.path(extract_dir, "occurrences.csv")
meas_file <- file.path(extract_dir, "measurements.csv")

if (!file.exists(occ_file) || !file.exists(meas_file)) {
  stop("Expected files not found after extraction.\n",
       "Check zip contents: ", existing_zip)
}

# =============================================================================
# Step 3: Parse raw species-level MycorrhizalType records from the DwC-A
# =============================================================================

if (file.exists(ckpt_species_raw)) {
  species_raw <- readr::read_csv(ckpt_species_raw, show_col_types = FALSE)
} else {

  # occurrences.csv: ID and scientificName (read only the two needed columns)
  occ <- readr::read_csv(
    occ_file,
    col_select  = c("ID", "scientificName"),
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  # measurements.csv: coreid, measurementType, measurementValue
  # Filter to "Mycorrhiza type" only — one row per coreid after this filter
  meas <- readr::read_csv(
    meas_file,
    col_types = readr::cols(.default = "c"),
    name_repair = "minimal"
  )
  colnames(meas) <- c("coreid", "measurementType", "measurementValue")

  meas_mtype <- dplyr::filter(meas, measurementType == "Mycorrhiza type") |>
    dplyr::select(coreid, MycorrhizalType = measurementValue)
  rm(meas)

  # Join: measurements.coreid = occurrences.ID → get scientificName
  species_raw <- dplyr::inner_join(
    meas_mtype,
    dplyr::rename(occ, coreid = ID),
    by = "coreid"
  ) |>
    dplyr::select(PlantBinomial = scientificName, MycorrhizalType) |>
    dplyr::filter(!is.na(PlantBinomial), !is.na(MycorrhizalType),
                  nzchar(PlantBinomial), nzchar(MycorrhizalType)) |>
    dplyr::distinct()

  readr::write_csv(species_raw, ckpt_species_raw)
}

# =============================================================================
# Step 3b: Read FungalRoot Table S2 (genus-level recommendation)
# =============================================================================
# Structurally distinct source from the DwC-A occurrence data above: a single
# hand-curated genus-level recommendation table, published as supplementary
# information (see header comment for what this does and does not replace).

if (file.exists(ckpt_table_s2)) {
  table_s2 <- readr::read_csv(ckpt_table_s2, show_col_types = FALSE)
} else {

  s2_file <- list.files(fr_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(s2_file) != 1L) {
    stop("Expected exactly one .xlsx file in ", fr_dir, "; found ",
         length(s2_file), ".")
  }

  table_s2_raw <- readxl::read_excel(
    s2_file, sheet = "Table S2", skip = 2, col_names = TRUE
  )
  colnames(table_s2_raw)[1:2] <- c("Genus", "MycorrhizalTypeGenus")

  table_s2 <- table_s2_raw |>
    dplyr::select(Genus, MycorrhizalTypeGenus) |>
    dplyr::filter(!is.na(Genus), nzchar(trimws(Genus)))

  readr::write_csv(table_s2, ckpt_table_s2)
}

# Inspect distinct MycorrhizalTypeGenus values before trusting the qualifying
# filter below -- the published table has a handful of messy free-text rows
# (e.g. "species-specific: EcM or NM-AM", "NM-AM, rarely EcM") and at least
# one apparent data-entry error (genus Thysanotus carries the value
# "Thysanothus", not a valid mycorrhizal-type code). These are logged and
# EXCLUDED by default because they do not match GENUS_EM_QUALIFYING_TYPES
# via exact-match filtering; this is the intended, conservative default
# (silently including a free-text/error row as a qualifying genus would be
# worse than silently excluding one that never had a clean EcM/EcM-AM call).

genus_qualifying_s2 <- table_s2 |>
  dplyr::filter(MycorrhizalTypeGenus %in% GENUS_EM_QUALIFYING_TYPES) |>
  dplyr::distinct(Genus) |>
  dplyr::pull(Genus)

# =============================================================================
# Step 4: Apply the species rule — identify EcM-demonstrated species
# =============================================================================
# A species is "ecm_demonstrated" if at least one raw occurrence record
# carries one of the three unambiguous EM_POSITIVE_LABELS defined above. No
# other recoding or reconciliation of MycorrhizalType labels is performed:
# every other label (ambiguous, AM-only, non-mycorrhizal, etc.) is simply
# irrelevant to this binary call.

# GBIF scientificName strings include author citations (e.g. "Abies alba
# Mill."); take the first two whitespace-delimited tokens as the working
# binomial and discard anything that doesn't reduce to a clean two-word name.
species_raw <- species_raw |>
  dplyr::mutate(
    PlantBinomial = trimws(PlantBinomial),
    BinomialClean = sub("^(\\S+\\s+\\S+).*$", "\\1", PlantBinomial)
  ) |>
  dplyr::filter(grepl("^\\S+\\s+\\S+$", BinomialClean))

positive_raw <- dplyr::filter(species_raw, MycorrhizalType %in% EM_POSITIVE_LABELS)

demonstrated_binomials <- sort(unique(positive_raw$BinomialClean))

# =============================================================================
# Step 5: Standardize the demonstrated candidate names via GBIF backbone
# =============================================================================
# Only species with >=1 positive record are sent to the backbone — species
# with none can never satisfy the species rule (under any synonym), so
# resolving their synonymy is unnecessary. (The Table S2 genus rule does not
# depend on this backbone step at all — see Step 3b.)

if (file.exists(ckpt_species_backbone)) {
  backbone_sp <- readr::read_csv(ckpt_species_backbone, show_col_types = FALSE)
} else {
  backbone_sp <- backbone_standardize(demonstrated_binomials)
  readr::write_csv(backbone_sp, ckpt_species_backbone)
}

name_map_sp <- extract_updated_name(backbone_sp)

# =============================================================================
# Step 6: Resolve synonyms and build the final species table
# =============================================================================
# Multiple raw binomials (e.g. a basionym and its current combination) can
# resolve to the same backbone-accepted name; this step collapses those onto
# one row per accepted name, retaining the raw name(s) and qualifying
# evidence label(s) for transparency.

species_resolved <- positive_raw |>
  dplyr::select(BinomialClean, MycorrhizalType) |>
  dplyr::left_join(
    dplyr::rename(name_map_sp, BinomialClean = verbatim_name),
    by = "BinomialClean"
  ) |>
  dplyr::mutate(
    UpdatedPlantBinomial = dplyr::if_else(
      is.na(UpdatedName) | !nzchar(UpdatedName),
      BinomialClean,
      UpdatedName
    )
  ) |>
  dplyr::select(PlantBinomial = BinomialClean, MycorrhizalType, UpdatedPlantBinomial)

species_clean_final <- species_resolved |>
  dplyr::group_by(UpdatedPlantBinomial) |>
  dplyr::summarise(
    PlantBinomial = paste(sort(unique(PlantBinomial)), collapse = "; "),
    evidence      = paste(sort(unique(MycorrhizalType)), collapse = "; "),
    .groups       = "drop"
  ) |>
  dplyr::mutate(
    UpdatedGenus     = sub("^(\\S+).*$", "\\1", UpdatedPlantBinomial),
    ecm_demonstrated = TRUE,
    evidence_source  = "occurrence"  # tags this route now that a second
                                      # (Table S2 genus) route exists
  ) |>
  dplyr::arrange(UpdatedPlantBinomial) |>
  dplyr::select(UpdatedPlantBinomial, UpdatedGenus, ecm_demonstrated, evidence,
                PlantBinomial, evidence_source)

# =============================================================================
# Step 7: Save species-rule output
# =============================================================================

readr::write_csv(species_clean_final, OUT_SPECIES)

# =============================================================================
# Step 7b: Save Table S2 genus-rule output
# =============================================================================
# Separate output table: qualifies GENERA, not species. A Canadian native
# species can qualify as an EcM host via this table even if it has no row of
# its own in clean_fungalroot_species.csv (see header "Output" section).

genus_table_out <- table_s2 |>
  dplyr::filter(Genus %in% genus_qualifying_s2) |>
  dplyr::rename(mycorrhizal_type = MycorrhizalTypeGenus) |>
  dplyr::mutate(evidence_source = "table_s2")

readr::write_csv(genus_table_out, OUT_GENERA)

# =============================================================================
# Summary
# =============================================================================

