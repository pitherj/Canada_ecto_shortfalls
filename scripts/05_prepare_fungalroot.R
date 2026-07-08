# =============================================================================
# Prepare FungalRoot Data: Download, Parse, Classify, Standardize
# =============================================================================
# Produces a computationally reproducible table of EcM host plant species
# "demonstrated" by the FungalRoot database, used to determine which plant
# genera (and species) count as EcM hosts for the Eltonian and Hutchinsonian
# shortfall analyses.
#
# Sources (two independent sources as of 2026-07-05 — see "Design history"
# below):
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
# Host-determination rule (two independent routes, combined with OR; see
# "Design history" below for how this rule has evolved):
#
#   1. SPECIES RULE (unchanged since 2026-06-27). A species is
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
#        - "EcM,AM" (dual-positive call; excluded 2026-06-27, see "Design
#          history") — asserts both types were found in the same record, which
#          in practice behaved as a second route to false positives: 38% of
#          genera that qualified under the original rule (including the
#          Canadian list's Acer, Allium, Fraxinus, Juglans, Rubus, Sambucus)
#          did so ONLY via this label, with no unambiguous EcM evidence
#          anywhere in the genus.
#
#   2. GENUS RULE (replaced 2026-07-05; see "Design history" below). A genus
#      counts as an EcM host genus if FungalRoot's own published Table S2
#      recommendation (Soudzilovskaia et al. 2020 supporting information;
#      genus-level calls made at >=67% consistency of species diagnosis)
#      assigns it "EcM" or "EcM-AM" (`GENUS_EM_QUALIFYING_TYPES`). This
#      replaces the previous ad hoc rule ("any species anywhere in FungalRoot
#      satisfies the species rule"), which was a non-standard, home-grown
#      construction not used by other FungalRoot-based studies.
#
#   3. COMBINATION. The two routes are independent and combined with OR: a
#      species/genus qualifies as an EcM host if it satisfies the species
#      rule directly, OR its genus qualifies under the Table S2 genus rule
#      (even if that particular species has no occurrence-level record of its
#      own). This script does NOT perform that OR merge itself — it produces
#      two separate, clearly labelled outputs (see "Output" below), each
#      tagged with its `evidence_source`. The merge is downstream code's
#      responsibility; see the closing summary below for the required
#      06_host_species.R follow-up this creates.
#
# Design history (for the audit trail):
#   Until 2026-06-27, genus-level host status was read directly from
#   Soudzilovskaia et al. (2020) Supplementary Table S2, a hand-curated
#   genus-level spreadsheet. On 2026-06-27 that table was dropped entirely in
#   favour of deriving genus status directly from the occurrence-level DwC-A
#   species data (the "any demonstrated species anywhere" rule), because Table
#   S2 (a) implicitly restricted hosts in places to trees/shrubs, and (b)
#   could not resolve species-level synonym cases (e.g. Persicaria vivipara /
#   Bistorta vivipara) correctly. Several alternative genus rules (simple
#   majority vote; 2/3 global proportion; 2/3 Canada-restricted proportion; a
#   minimum-n backoff to a global proportion) were evaluated and rejected at
#   that time — each either wrongly excluded Bistorta vivipara, miscalibrated
#   against biogeographically diluted genera (e.g. Salix, Betula), or
#   reintroduced an arbitrary free parameter.
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
#   Cause (2) was accepted at the time as the known cost of the no-threshold
#   occurrence-based genus rule.
#
#   Reintroduction of Table S2, genus-level only (2026-07-05): the
#   occurrence-based genus rule above was never used downstream in practice —
#   06_host_species.R has selected Canadian host species from
#   `UpdatedPlantBinomial` (species-level only) since 2026-06-27, with
#   `UpdatedGenus` retained only for transparency. It also had no way to
#   express genuine dual-mycorrhizal-type genera, since it inherits a binary
#   EcM/not-EcM call from the occurrence species rule. FungalRoot's own
#   published Table S2 genus-level recommendation is reinstated here, this
#   time restricted to genus-level qualification only (the species-level rule
#   above, and the synonym problem that caused Table S2 to be dropped for
#   species-level use, are both unaffected). Verified in this project's copy
#   of Table S2 (`data_raw/fungalroot/nph16569-sup-0002-tabless1-s4.xlsx`):
#     - Acer, Fraxinus, and Juglans are all called "AM" in Table S2,
#       consistent with (not contradicting) their exclusion under the
#       occurrence-based species rule — Table S2 does not newly include them.
#     - Salix, Populus, and Eucalyptus are called "EcM-AM" (genuine dual-type
#       genera) — information the binary occurrence-based rule structurally
#       cannot express. These are the concrete new cases Table S2 adds.
#     - Bistorta and Persicaria (both spellings present as separate genus
#       rows) are called "NM-AM" (non-EcM) by Table S2. Bistorta vivipara
#       therefore remains EcM-demonstrated ONLY via its own direct
#       occurrence-level record (species rule, unchanged) — Table S2 gives no
#       genus-level support for Bistorta/Persicaria, which is exactly why the
#       two rules are combined with OR rather than the genus rule replacing
#       the species rule.
#     - Regarding Dryopteris specifically: the ad hoc occurrence-based genus
#       rule this replaces would have wrongly included Dryopteris under an
#       earlier version of the species rule (before "EcM, AM undetermined"
#       was excluded — see the species rule above). Under the CURRENT
#       (already-refined) species rule, Dryopteris has zero species with an
#       unambiguous EcM-positive record, so it was already excluded from the
#       implicit occurrence-based genus set before this change; Table S2's
#       independent "AM" call for Dryopteris corroborates that exclusion but
#       was not itself required to fix it. Flagged to Jason 2026-07-05 so the
#       rationale here is accurate rather than presuming a bug that the
#       species-rule refinement had already resolved.
#   Table S2 contains a small number of messy free-text / apparent
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
#   data_derived/clean_fungalroot_species.csv  (species rule; UNCHANGED
#     structure, plus one new column)
#     UpdatedPlantBinomial — backbone-resolved accepted name
#     UpdatedGenus         — genus parsed from UpdatedPlantBinomial (no longer
#                            the basis of any genus-level rule — see Design
#                            history above; retained for transparency only)
#     ecm_demonstrated     — always TRUE in this table (only demonstrated
#                            species are retained; kept as an explicit column
#                            for self-documentation and so downstream code
#                            does not have to infer it)
#     evidence             — the qualifying raw MycorrhizalType label(s)
#                            backing the call, semicolon-separated
#     PlantBinomial         — the raw FungalRoot name(s) (pre-backbone),
#                            semicolon-separated when more than one raw
#                            synonym resolves to the same accepted name
#     evidence_source      — always "occurrence" in this table (new column;
#                            distinguishes this route from the Table S2 route
#                            below now that there are two)
#
#   data_derived/clean_fungalroot_genera_table_s2.csv  (genus rule; NEW output)
#     Genus                 — genus name as given in Table S2
#     mycorrhizal_type      — Table S2's MycorrhizalTypeGenus call ("EcM" or
#                            "EcM-AM" only; this table is pre-filtered to
#                            qualifying genera)
#     evidence_source      — always "table_s2"
#
#   IMPORTANT — this changes the script's output contract: there are now TWO
#   output tables, not one, and a Canadian native species can qualify as an
#   EcM host via either (its own row in the species table, OR its genus
#   appearing in the genus table) without itself appearing in the species
#   table. 06_host_species.R currently reads only the species table and does
#   NOT yet consult the genus table — this is a required follow-up, flagged
#   explicitly rather than fixed silently here (see closing summary below).
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
# Step 3b: Read FungalRoot Table S2 (genus-level recommendation)
# =============================================================================
# Structurally distinct source from the DwC-A occurrence data above: a single
# hand-curated genus-level recommendation table, published as supplementary
# information (see header comment for what this does and does not replace).

if (file.exists(ckpt_table_s2)) {
  ts("Step 3b: Loading checkpointed Table S2...")
  table_s2 <- readr::read_csv(ckpt_table_s2, show_col_types = FALSE)
  ts(sprintf("  Rows: %d", nrow(table_s2)))
} else {
  ts("Step 3b: Reading FungalRoot Table S2 (genus-level recommendation)...")

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

  ts(sprintf("  Rows read: %d", nrow(table_s2)))
  readr::write_csv(table_s2, ckpt_table_s2)
  ts("  Saved Table S2 checkpoint.")
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
unusual_values <- table_s2 |>
  dplyr::filter(!MycorrhizalTypeGenus %in%
                  c("AM", "NM-AM", "NM", "OM", "uncertain", "EcM-AM", "EcM", "ErM")) |>
  dplyr::distinct(MycorrhizalTypeGenus)
if (nrow(unusual_values) > 0L) {
  ts(sprintf("  NOTE: %d unusual MycorrhizalTypeGenus value(s) found and excluded by default: %s",
             nrow(unusual_values), paste(unusual_values$MycorrhizalTypeGenus, collapse = " | ")))
}

genus_qualifying_s2 <- table_s2 |>
  dplyr::filter(MycorrhizalTypeGenus %in% GENUS_EM_QUALIFYING_TYPES) |>
  dplyr::distinct(Genus) |>
  dplyr::pull(Genus)

ts(sprintf("  Qualifying genera (EcM or EcM-AM in Table S2): %d",
           length(genus_qualifying_s2)))

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
# with none can never satisfy the species rule (under any synonym), so
# resolving their synonymy is unnecessary. (The Table S2 genus rule does not
# depend on this backbone step at all — see Step 3b.)

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
    ecm_demonstrated = TRUE,
    evidence_source  = "occurrence"  # tags this route now that a second
                                      # (Table S2 genus) route exists
  ) |>
  dplyr::arrange(UpdatedPlantBinomial) |>
  dplyr::select(UpdatedPlantBinomial, UpdatedGenus, ecm_demonstrated, evidence,
                PlantBinomial, evidence_source)

ts(sprintf("  Final demonstrated species: %d  |  qualifying genera: %d",
           nrow(species_clean_final),
           dplyr::n_distinct(species_clean_final$UpdatedGenus)))

# =============================================================================
# Step 7: Save species-rule output
# =============================================================================

readr::write_csv(species_clean_final, OUT_SPECIES)
ts(sprintf("Saved -> %s", basename(OUT_SPECIES)))

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
ts(sprintf("Saved -> %s (%d qualifying genera)", basename(OUT_GENERA), nrow(genus_table_out)))

# =============================================================================
# Summary
# =============================================================================

genera_species_route <- unique(species_clean_final$UpdatedGenus)
genera_table_s2_route <- unique(genus_table_out$Genus)
genera_union <- union(genera_species_route, genera_table_s2_route)

ts("=== 05_prepare_fungalroot.R complete ===")
ts(sprintf("Species output: %d EcM-demonstrated species (%d genera) -> %s",
           nrow(species_clean_final),
           length(genera_species_route),
           OUT_SPECIES))
ts(sprintf("Genus output:   %d Table S2 qualifying genera -> %s",
           length(genera_table_s2_route), OUT_GENERA))
ts(sprintf("Genera covered by EITHER route (union): %d  (overlap between routes: %d)",
           length(genera_union),
           length(intersect(genera_species_route, genera_table_s2_route))))
ts("")
ts("REQUIRED FOLLOW-UP (not made here — see header 'Output' section):")
ts("06_host_species.R currently selects Canadian host species directly from")
ts("UpdatedPlantBinomial (species-level only, as of 2026-06-27) and does NOT")
ts("consult clean_fungalroot_genera_table_s2.csv at all. It needs a companion")
ts("change to also flag any Canadian-native species whose genus appears in")
ts("clean_fungalroot_genera_table_s2.csv as host_demonstrated, even when that")
ts("species has no row of its own in clean_fungalroot_species.csv. This has")
ts("deliberately NOT been implemented as part of this script's changes.")
