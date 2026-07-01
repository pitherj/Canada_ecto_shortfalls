# =============================================================================
# Prepare FungalRoot Data: Download, Parse, Classify, Standardize
# =============================================================================
# Produces a computationally reproducible table of EcM host plant species
# "demonstrated" by the FungalRoot database, used to determine which plant
# genera (and species) count as EcM hosts for the Eltonian and Hutchinsonian
# shortfall analyses.
#
# Source (single source as of 2026-06-27 — see "Design history" below):
#   GBIF DwC-A archive of the FungalRoot dataset
#     GBIF dataset key: 744edc21-8dd2-474e-8a0b-b8c3d56a3c2d
#     Download URL pattern:
#       https://api.gbif.org/v1/occurrence/download/request/<KEY>.zip
#     OR download interactively from:
#       https://www.gbif.org/occurrence/download?dataset_key=744edc21-...
#
# DwC-A structure:
#   occurrences.csv: core file; column 0 = ID (primary key), column 44 = scientificName
#   measurements.csv: extension; coreid links to occurrences.csv ID;
#                     columns: coreid, measurementType, measurementValue
#
# Host-determination rule (replaces the Table S2 genus-spreadsheet approach
# used prior to 2026-06-27 — see "Design history" below):
#
#   1. SPECIES RULE. A species is "ecm_demonstrated" if at least one raw
#      FungalRoot occurrence record assigns it one of the following
#      unambiguous EcM-positive MycorrhizalType labels:
#        "EcM", "EcM, no AM colonization", "ErM,EcM"
#      Two labels are explicitly EXCLUDED as ambiguous:
#        - "EcM, AM undetermined" — does not distinguish a genuine EcM record
#          from a non-EcM record where AM status simply wasn't assessed (e.g.
#          ferns such as Dryopteris filix-mas, whose only EcM-flavoured record
#          carries this label and is contradicted by separate AM /
#          non-mycorrhizal records for the same species).
#        - "EcM,AM" (dual-positive call; excluded 2026-06-27, see "Design
#          history") — asserts both types were found in the same record, which
#          in practice behaved as a second route to false positives: 38% of
#          genera that qualified under the original rule (including the
#          Canadian list's Acer, Allium, Fraxinus, Juglans, Rubus, Sambucus)
#          did so ONLY via this label, with no unambiguous EcM evidence
#          anywhere in the genus.
#
#   2. GENUS RULE. A genus counts as an EcM host genus if AT LEAST ONE of its
#      species ANYWHERE in the global FungalRoot data (no geographic
#      restriction) satisfies the species rule above. There is no
#      proportion/majority threshold and no minimum sample size — this is a
#      deliberate choice (see "Design history") to avoid the asymmetric risk
#      of erroneously excluding a genuine host genus (e.g. Bistorta, via
#      B. vivipara) because most of its congeners are non-EcM and dilute a
#      proportion-based score below an arbitrary cutoff. The cost is the
#      converse error — accepting some genera (e.g. Dryopteris, a fern genus)
#      that may not be genuine EcM hosts. Per Jason's explicit guidance, this
#      trade-off is accepted: false negatives (excluding a real host) are
#      considered worse than false positives (including a doubtful one).
#
# Design history (for the audit trail):
#   Until 2026-06-27, genus-level host status was read directly from
#   Soudzilovskaia et al. (2020) Supplementary Table S2, a hand-curated
#   genus-level spreadsheet. That table (a) is no longer used at all, (b)
#   implicitly restricted hosts in places to trees/shrubs, and (c) could not
#   resolve species-level synonym cases (e.g. Persicaria vivipara / Bistorta
#   vivipara) correctly. The current approach instead derives genus status
#   directly and exclusively from the same occurrence-level DwC-A species
#   data used for the species table, via the two rules above. Several
#   alternative genus rules (simple majority vote; 2/3 global proportion;
#   2/3 Canada-restricted proportion; a minimum-n backoff to a global
#   proportion) were evaluated and rejected — each either wrongly excluded
#   Bistorta vivipara, miscalibrated against biogeographically diluted genera
#   (e.g. Salix, Betula), or reintroduced an arbitrary free parameter. The
#   "any demonstrated species anywhere" rule was chosen because it requires
#   no free parameters and passes both worked test cases (correctly excludes
#   Dryopteris's 33 species, none of which has unambiguous EcM evidence;
#   correctly includes Bistorta via B. vivipara).
#
#   Refinement (2026-06-27, same day): an initial version of the species rule
#   also accepted the dual-positive label "EcM,AM" as qualifying evidence. On
#   first full run this produced 698 native Canadian host species — a 4.7x
#   jump from the prior Table-S2-derived figure of 148 — and inspection showed
#   two distinct causes. (1) 38% of qualifying genera (118 of 698 Canadian
#   species, across 20 genera including Acer, Allium, Fraxinus, Juglans,
#   Rubus, Sambucus) qualified ONLY via "EcM,AM", with no unambiguous EcM
#   evidence anywhere in the genus — this is the cause addressed by dropping
#   the label, as done here. (2) Separately, a handful of mega-diverse genera
#   with no other EcM signal (e.g. Potentilla, Saxifraga, Polygonum,
#   Pedicularis) qualified via a single unambiguous "EcM, no AM colonization"
#   record in just 1-3 species globally, then propagated via the genus rule to
#   dozens of Canadian congeners (155 species from 7 raw records, combined).
#   Cause (2) is accepted as the deliberate, known cost of the no-threshold
#   genus rule (see above) and is NOT addressed by this label change — it
#   remains an open design trade-off, flagged to Jason 2026-06-27.
#
# Name standardization:
#   Uses rgbif::name_backbone_checklist() to resolve taxonomic synonyms (this
#   is what lets "Persicaria vivipara" resolve to the same accepted name as
#   "Bistorta vivipara"). Only species with at least one EcM-positive raw
#   record are sent to the backbone lookup (species with no positive record
#   can never affect either the species or genus rule, so resolving their
#   synonymy is unnecessary and would needlessly inflate the number of API
#   calls). Batched at 200 names per request.
#
# Raw data files (data_raw/fungalroot/):
#   744edc21-...zip — GBIF DwC-A archive (scanned by prefix; downloaded via
#                     rgbif::occ_download() if absent)
#
# Checkpoints (data_derived/checkpoints/):
#   fungalroot_dwca_extracted/        — extracted DwC-A directory
#   fungalroot_species_raw.csv        — parsed (PlantBinomial, MycorrhizalType)
#                                        before classification, one row per
#                                        distinct raw occurrence/measurement pair
#   fungalroot_species_backbone.csv   — rgbif backbone results for the
#                                        EcM-demonstrated candidate species
#
# Output:
#   data_derived/clean_fungalroot_species.csv
#     UpdatedPlantBinomial — backbone-resolved accepted name
#     UpdatedGenus         — genus parsed from UpdatedPlantBinomial
#     ecm_demonstrated     — always TRUE in this table (only demonstrated
#                            species are retained; kept as an explicit column
#                            for self-documentation and so downstream code
#                            does not have to infer it)
#     evidence             — the qualifying raw MycorrhizalType label(s)
#                            backing the call, semicolon-separated
#     PlantBinomial         — the raw FungalRoot name(s) (pre-backbone),
#                            semicolon-separated when more than one raw
#                            synonym resolves to the same accepted name
#
#   This table is the sole input to 06_host_species.R, which derives
#   the qualifying genus list (em_genera) from `UpdatedGenus`, and the
#   per-species `host_demonstrated` flag from `UpdatedPlantBinomial`.
#
# Runtime notes:
#   - rgbif backbone queries: a few minutes, scales with the number of
#     EcM-demonstrated candidate species (a small fraction of all ~17,000
#     raw FungalRoot occurrence/measurement rows).
#   - Delete checkpoint files to force re-run of individual steps.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(rgbif)

# ---- Paths -------------------------------------------------------------------

fr_dir   <- here("data_raw", "fungalroot")
out_dir  <- here("data_derived")
temp_dir <- paths$temp_dir  # alias for data_derived/checkpoints (see 00_setup.R)

DWCA_ZIP_PATTERN <- "744edc21"  # prefix to locate existing zip

OUT_SPECIES <- file.path(out_dir, "clean_fungalroot_species.csv")

ckpt_species_raw      <- file.path(temp_dir, "fungalroot_species_raw.csv")
ckpt_species_backbone <- file.path(temp_dir, "fungalroot_species_backbone.csv")

# The four unambiguous EcM-positive raw MycorrhizalType labels (see header
# comment above for the rationale, including why "EcM, AM undetermined" is
# deliberately excluded).
EM_POSITIVE_LABELS <- c("EcM", "EcM, no AM colonization", "ErM,EcM")

# Simple timestamped message helper
ts <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

# ---- Helper: rgbif backbone lookup, batched ----------------------------------

backbone_standardize <- function(names_vec, batch_size = 200L) {
  names_vec <- as.character(names_vec)
  n_batches <- ceiling(length(names_vec) / batch_size)
  results   <- vector("list", n_batches)
  for (i in seq_len(n_batches)) {
    idx <- seq((i - 1L) * batch_size + 1L,
               min(i * batch_size, length(names_vec)))
    ts(sprintf("  Backbone batch %d / %d (%d names)...", i, n_batches, length(idx)))
    results[[i]] <- tryCatch(
      rgbif::name_backbone_checklist(name_data = data.frame(scientificName = names_vec[idx])),
      error = function(e) {
        message("  Batch ", i, " failed: ", conditionMessage(e))
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

ts("Step 1: Locating FungalRoot DwC-A zip...")

existing_zips <- list.files(fr_dir, pattern = DWCA_ZIP_PATTERN, full.names = TRUE)
existing_zip  <- existing_zips[grepl("\\.zip$", existing_zips)]

if (length(existing_zip) == 0L) {
  ts("  No local zip found. Attempting GBIF API download...")
  ts("  Note: GBIF requires a registered download. Initiating via rgbif::occ_download()...")
  ts("  This requires GBIF credentials in .Renviron:")
  ts("    GBIF_USER=your_username")
  ts("    GBIF_PWD=your_password")
  ts("    GBIF_EMAIL=your_email")

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
  ts(sprintf("  Download key: %s  (GBIF will email when ready)", dl_key))
  ts("  Waiting for download to complete (this may take 5–20 min)...")
  rgbif::occ_download_wait(dl_key, status_ping = 30L)

  dl_info  <- rgbif::occ_download_meta(dl_key)
  zip_dest <- file.path(fr_dir, paste0(dl_key, ".zip"))
  rgbif::occ_download_get(dl_key, path = fr_dir, overwrite = TRUE)
  existing_zip <- zip_dest
  ts(sprintf("  Downloaded: %s", basename(existing_zip)))
} else {
  existing_zip <- existing_zip[1L]
  ts(sprintf("  Using existing zip: %s", basename(existing_zip)))
}

# =============================================================================
# Step 2: Extract DwC-A zip
# =============================================================================

ts("Step 2: Extracting DwC-A zip...")

extract_dir <- file.path(temp_dir, "fungalroot_dwca_extracted")

# Identify which files are needed (avoid re-extracting if present)
need_extract <- !all(file.exists(
  file.path(extract_dir, c("occurrences.csv", "measurements.csv"))
))

if (need_extract) {
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  unzip(existing_zip, exdir = extract_dir, overwrite = TRUE)
  ts(sprintf("  Extracted to: %s", extract_dir))
} else {
  ts("  Extracted files already present; skipping unzip.")
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
  ts("Step 3: Loading checkpointed raw species data...")
  species_raw <- readr::read_csv(ckpt_species_raw, show_col_types = FALSE)
  ts(sprintf("  Rows: %d", nrow(species_raw)))
} else {
  ts("Step 3: Parsing DwC-A species data...")

  # occurrences.csv: ID and scientificName (read only the two needed columns)
  occ <- readr::read_csv(
    occ_file,
    col_select  = c("ID", "scientificName"),
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  ts(sprintf("  Occurrence records: %d", nrow(occ)))

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
  ts(sprintf("  Mycorrhiza type measurement rows: %d", nrow(meas_mtype)))
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

  ts(sprintf("  Distinct (PlantBinomial, MycorrhizalType) rows: %d", nrow(species_raw)))
  readr::write_csv(species_raw, ckpt_species_raw)
  ts("  Saved raw species checkpoint.")
}

# =============================================================================
# Step 4: Apply the species rule — identify EcM-demonstrated species
# =============================================================================
# A species is "ecm_demonstrated" if at least one raw occurrence record
# carries one of the three unambiguous EM_POSITIVE_LABELS defined above. No
# other recoding or reconciliation of MycorrhizalType labels is performed:
# every other label (ambiguous, AM-only, non-mycorrhizal, etc.) is simply
# irrelevant to this binary call.

ts("Step 4: Applying species rule (>=1 unambiguous EcM-positive record)...")

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
ts(sprintf("  Distinct raw binomials with >=1 positive record: %d",
           length(demonstrated_binomials)))

# =============================================================================
# Step 5: Standardize the demonstrated candidate names via GBIF backbone
# =============================================================================
# Only species with >=1 positive record are sent to the backbone — species
# with none can never satisfy the species rule (under any synonym) nor
# contribute to the genus rule, so resolving their synonymy is unnecessary.

ts(sprintf("Step 5: Standardizing %d candidate species names via GBIF backbone...",
           length(demonstrated_binomials)))

if (file.exists(ckpt_species_backbone)) {
  ts("  Loading checkpointed backbone results...")
  backbone_sp <- readr::read_csv(ckpt_species_backbone, show_col_types = FALSE)
} else {
  backbone_sp <- backbone_standardize(demonstrated_binomials)
  readr::write_csv(backbone_sp, ckpt_species_backbone)
  ts("  Saved backbone checkpoint.")
}

name_map_sp <- extract_updated_name(backbone_sp)

match_summary <- table(backbone_sp$matchType, useNA = "ifany")
ts(sprintf("  Match types: %s",
           paste(names(match_summary), as.integer(match_summary),
                 sep = " = ", collapse = ", ")))

# =============================================================================
# Step 6: Resolve synonyms and build the final species table
# =============================================================================
# Multiple raw binomials (e.g. a basionym and its current combination) can
# resolve to the same backbone-accepted name; this step collapses those onto
# one row per accepted name, retaining the raw name(s) and qualifying
# evidence label(s) for transparency.

ts("Step 6: Building final species table...")

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
    ecm_demonstrated = TRUE
  ) |>
  dplyr::arrange(UpdatedPlantBinomial) |>
  dplyr::select(UpdatedPlantBinomial, UpdatedGenus, ecm_demonstrated, evidence, PlantBinomial)

ts(sprintf("  Final demonstrated species: %d  |  qualifying genera: %d",
           nrow(species_clean_final),
           dplyr::n_distinct(species_clean_final$UpdatedGenus)))

# =============================================================================
# Step 7: Save output
# =============================================================================

readr::write_csv(species_clean_final, OUT_SPECIES)
ts(sprintf("Saved -> %s", basename(OUT_SPECIES)))

# =============================================================================
# Summary
# =============================================================================

ts("=== 05_prepare_fungalroot.R complete ===")
ts(sprintf("Species output: %d EcM-demonstrated species (%d genera) -> %s",
           nrow(species_clean_final),
           dplyr::n_distinct(species_clean_final$UpdatedGenus),
           OUT_SPECIES))
ts("Downstream: 06_host_species.R selects Canadian host species")
ts("directly from UpdatedPlantBinomial (species-level only, as of 2026-06-27;")
ts("UpdatedGenus is retained here for transparency but no longer used to")
ts("select species downstream).")
