# =============================================================================
# EcM Native Host Species for Canada (BIEN-based)
# =============================================================================
# Assemble a list of EcM host plant species that are (a) directly demonstrated
# as EcM hosts by the FungalRoot species-level table (see NOTE below),
# (b) native to Canada according to the Native Status Resolver (NSR), and
# (c) present in the BIEN range database. No growth-form restriction is
# applied at this stage — herbaceous species are included on the same
# footing as trees and shrubs.
#
# In addition, BIEN trait data are used to assign each species a growth form
# (e.g., tree, shrub, herb), purely as descriptive metadata for downstream
# figures/tables. Where multiple trait values exist for a species, the most
# frequently reported value is retained.
#
# Workflow:
#   1.  Load FungalRoot species-level table and derive the demonstrated-
#       species name set (see NOTE)
#   2.  Query BIEN for all plant species with distributions in Canada
#   3.  Filter to native Canadian species using NSR (batched, 500 spp/request)
#   4.  Intersect with the demonstrated-species set (species-level, not genus)
#   5.  Retrieve growth forms from BIEN trait database
#   6.  Join growth forms to host species table
#   7a. Impute missing growth_form from GIFT trait 1.2.1
#   7b. Congener-modal fallback for any growth_form still missing
#   7c. Curated manual overrides for any growth_form still missing after 7a/7b
#   8.  Save (with host_demonstrated flag — see NOTE)
#
# NOTE on EcM-host determination (redesigned 2026-06-27; refined again same
# day — see "Design history" below):
# Host status is derived exclusively from the species-level FungalRoot table
# (data_derived/clean_fungalroot_species.csv, produced by 05_prepare_fungalroot.R
# directly from GBIF DwC-A occurrence-level "Mycorrhiza type" records). A
# species is "ecm_demonstrated" if FungalRoot carries >=1 unambiguous
# EcM-positive occurrence record for it (see EM_POSITIVE_LABELS in
# 05_prepare_fungalroot.R). Selection into this table requires the SPECIES
# ITSELF to satisfy that rule — there is no genus-level fallback. Every row
# therefore has host_demonstrated == TRUE by construction; the column is kept
# (rather than omitted) purely for self-documentation, mirroring the
# `ecm_demonstrated` convention in clean_fungalroot_species.csv.
#
# Design history (for the audit trail):
#   05_prepare_fungalroot.R also defines a GENUS RULE (a genus qualifies as
#   an EcM host genus, "host_broad", if >=1 of its species ANYWHERE in
#   FungalRoot satisfies the species rule, with no proportion/majority
#   threshold). This script originally (2026-06-27, first pass) selected
#   Canadian host species by genus membership in that broader set, on the
#   reasoning that a species-only rule would wrongly exclude a genuine host
#   like Bistorta vivipara if its congeners diluted a proportion-based score.
#   That version returned 698 species (initially), or 580 after a same-day
#   refinement excluding the ambiguous "EcM,AM" label (see 1-0's header) —
#   versus 148 species under the prior Table-S2-derived approach, a ~4x jump.
#   Inspection showed the increase was structurally dominated by a small
#   number of mega-diverse genera (e.g. Potentilla, Saxifraga, Polygonum,
#   Pedicularis) qualifying via a single unambiguous record in 1-3 species
#   globally, then propagating via the genus rule to dozens of Canadian
#   congeners with no other EcM signal (155 species from 7 raw records,
#   combined). Jason judged the direct species-level count (145 species,
#   intersected against the full BIEN Canada-native flora) the more credible
#   figure and asked that the pipeline select on it directly. The genus rule
#   itself is unaffected and still computed in 05_prepare_fungalroot.R
#   (`UpdatedGenus` column of clean_fungalroot_species.csv) for transparency
#   and possible future use, but as of this revision it is NOT used to select
#   species into this script's output. Bistorta vivipara remains correctly
#   included because it has its own direct, unambiguous FungalRoot record —
#   the genus-level fallback was never actually load-bearing for that case.
#
# Checkpoint files (data_derived/checkpoints/):
#   bien_ecm_canada_species.csv     — species list after NSR native filter
#                                      and demonstrated-species intersection
#   bien_ecm_growthforms.csv        — growth form trait data from BIEN
#   gift_growthforms.csv            — growth form trait data from GIFT
#                                     (used only for species lacking BIEN data)
#
# Output:
#   data_derived/ecm_native_canada_host_species.csv
#     Columns: species, host_demonstrated (always TRUE), growth_form
#
# Runtime notes:
#   - Step 2 (BIEN_list_country) ~5 min
#   - Step 3 (NSR batched)       ~5–15 min (depends on server)
#   - Step 5 (BIEN traits)       ~10–30 min
#   All slow steps are checkpointed; delete checkpoint files to force re-run.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(BIEN)
library(NSR)
library(GIFT)

sf::sf_use_s2(FALSE)

# Note: an earlier copy of clean_host_name() lived here but was never called;
# host-name cleaning is now centralized in canonicalize_host() (00_setup.R).

species_ckpt    <- file.path(paths$temp_dir, "bien_ecm_canada_species.csv")
growthform_ckpt <- file.path(paths$temp_dir, "bien_ecm_growthforms.csv")
gift_gf_ckpt    <- file.path(paths$temp_dir, "gift_growthforms.csv")

# ---- Step 1: FungalRoot-derived demonstrated species ------------------------

ts("Step 1: Loading FungalRoot species table; deriving demonstrated species...")
ft_species <- readr::read_csv(paths$fungalroot_sp, show_col_types = FALSE)
demonstrated_species <- sort(unique(ft_species$UpdatedPlantBinomial))
ts(sprintf("  EcM-demonstrated species (host_demonstrated): %d",
           length(demonstrated_species)))

# ---- Steps 2–4: BIEN query + NSR native filter + FungalRoot intersection ----

if (file.exists(species_ckpt)) {
  ts("Steps 2-4: Loading checkpointed native EcM host species list...")
  em_canada_species <- readr::read_csv(species_ckpt, show_col_types = FALSE)$species
  ts(sprintf("  EcM host species native to Canada (from checkpoint): %d",
             length(em_canada_species)))

  # Validate against the current FungalRoot demonstrated-species list.
  # Stale checkpoints can contain species that were demonstrated in an older
  # FungalRoot version (or selected under the now-superseded genus rule) but
  # no longer qualify directly; remove them now.
  invalid_spp <- em_canada_species[!em_canada_species %in% demonstrated_species]
  if (length(invalid_spp) > 0L) {
    ts(sprintf("  Warning: %d checkpoint species no longer in current FungalRoot demonstrated-species list — removing:",
               length(invalid_spp)))
    for (sp in invalid_spp) ts(sprintf("    - %s", sp))
    em_canada_species <- em_canada_species[em_canada_species %in% demonstrated_species]
    ts(sprintf("  Validated species count: %d", length(em_canada_species)))
  }
} else {

  ts("Step 2: Querying BIEN for all plant species in Canada...")
  bien_canada <- BIEN::BIEN_list_country("Canada", new.world = TRUE,
                                          cultivated = FALSE)
  bien_canada_species <- unique(bien_canada$scrubbed_species_binomial)
  ts(sprintf("  BIEN species with distributions in Canada: %d",
             length(bien_canada_species)))

  ts("Step 3: Checking native status via NSR (batched, 500 species/request)...")
  NSR_BATCH <- 500L
  nsr_batches <- split(bien_canada_species,
                       ceiling(seq_along(bien_canada_species) / NSR_BATCH))
  nsr_list <- vector("list", length(nsr_batches))
  for (i in seq_along(nsr_batches)) {
    sp_batch <- nsr_batches[[i]]
    ts(sprintf("  NSR batch %d / %d (%d species)...",
               i, length(nsr_batches), length(sp_batch)))
    nsr_list[[i]] <- tryCatch(
      NSR::NSR_simple(species = sp_batch,
                      country = rep("Canada", length(sp_batch))),
      error = function(e) {
        message(sprintf("  NSR batch %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }
  nsr_status <- dplyr::bind_rows(Filter(Negate(is.null), nsr_list))

  status_tbl <- table(nsr_status$native_status, useNA = "ifany")
  ts(sprintf("  NSR status breakdown: %s",
             paste(names(status_tbl), as.integer(status_tbl),
                   sep = " = ", collapse = ", ")))

  native_species <- nsr_status$species[
    !is.na(nsr_status$native_status) & nsr_status$native_status == "N"
  ]
  ts(sprintf("  Species native to Canada (NSR status == 'N'): %d",
             length(native_species)))

  ts("Step 4: Intersecting with FungalRoot demonstrated species...")
  em_canada_species <- sort(intersect(native_species, demonstrated_species))
  ts(sprintf("  EcM host species native to Canada (species-level FungalRoot match): %d",
             length(em_canada_species)))

  readr::write_csv(data.frame(species = em_canada_species), species_ckpt)
  ts(sprintf("  Saved species checkpoint -> %s", basename(species_ckpt)))
}

# ---- Step 5: BIEN growth form trait data ------------------------------------
# BIEN_trait_traitbyspecies() returns one row per trait record; multiple
# values may exist per species. We retain the modal value per species.

if (file.exists(growthform_ckpt)) {
  ts("Step 5: Loading checkpointed growth form data...")
  growthforms_raw <- readr::read_csv(growthform_ckpt, show_col_types = FALSE)
} else {
  ts(sprintf("Step 5: Querying BIEN growth forms for %d species...",
             length(em_canada_species)))
  ts("  (This may take 10-30 min depending on BIEN server load)")

  # Batch to avoid timeouts; BIEN_trait_traitbyspecies handles a vector but
  # large vectors can time out.
  TRAIT_BATCH <- 200L
  trait_batches <- split(em_canada_species,
                         ceiling(seq_along(em_canada_species) / TRAIT_BATCH))
  trait_list <- vector("list", length(trait_batches))
  for (i in seq_along(trait_batches)) {
    ts(sprintf("  Trait batch %d / %d...", i, length(trait_batches)))
    trait_list[[i]] <- tryCatch(
      BIEN::BIEN_trait_traitbyspecies(
        species = trait_batches[[i]],
        trait   = "whole plant growth form"
      ) |>
        dplyr::select(scrubbed_species_binomial, trait_value) |>
        dplyr::mutate(trait_value = tolower(trait_value)),
      error = function(e) {
        message(sprintf("  Trait batch %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }
  growthforms_raw <- dplyr::bind_rows(Filter(Negate(is.null), trait_list))
  readr::write_csv(growthforms_raw, growthform_ckpt)
  ts(sprintf("  Saved growth form checkpoint -> %s", basename(growthform_ckpt)))
}

ts(sprintf("  Growth form records retrieved: %d", nrow(growthforms_raw)))

# ---- Step 6: Summarise to one growth form per species -----------------------
# Use the most commonly reported trait value; ties broken arbitrarily by
# slice_max (first alphabetically among ties).

ecm_species_with_growthform <- growthforms_raw |>
  dplyr::count(scrubbed_species_binomial, trait_value) |>
  dplyr::group_by(scrubbed_species_binomial) |>
  dplyr::slice_max(n, n = 1L, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(species = scrubbed_species_binomial, growth_form = trait_value)

ts(sprintf("  Species with growth form data: %d / %d",
           dplyr::n_distinct(ecm_species_with_growthform$species),
           length(em_canada_species)))

# ---- Step 7: Join to host species table -------------------------------------
# host_demonstrated is TRUE for every row by construction (selection in
# Step 4 already required direct membership in demonstrated_species); it is
# computed explicitly here rather than just hard-coded, as a self-documenting
# sanity check (mirrors the `ecm_demonstrated` convention in
# clean_fungalroot_species.csv).

host_species_table <- data.frame(species = em_canada_species) |>
  dplyr::mutate(host_demonstrated = species %in% demonstrated_species) |>
  dplyr::left_join(ecm_species_with_growthform, by = "species")

ts(sprintf("  host_demonstrated FALSE (should be 0): %d  |  growth_form NA: %d",
           sum(!host_species_table$host_demonstrated),
           sum(is.na(host_species_table$growth_form))))

# ---- Step 7a: Impute missing growth_form from GIFT trait 1.2.1 --------------
# BIEN traits did not cover all species. GIFT (trait_ID 1.2.1 = plant growth
# form) is used as a second source for remaining NAs. Checkpointed.

n_na_gf <- sum(is.na(host_species_table$growth_form))
if (n_na_gf > 0L) {
  ts(sprintf("Step 7a: Imputing growth_form for %d species from GIFT...", n_na_gf))

  if (file.exists(gift_gf_ckpt)) {
    ts("  Loading checkpointed GIFT growth form data...")
    gift_gf <- readr::read_csv(gift_gf_ckpt, show_col_types = FALSE)
  } else {
    ts("  Querying GIFT (trait 1.2.1 — plant growth form)...")
    gift_raw <- tryCatch(
      GIFT::GIFT_traits(trait_IDs  = "1.2.1",
                        agreement  = 0.66,
                        bias_ref   = FALSE,
                        bias_deriv = FALSE),
      error = function(e) {
        message("  GIFT query failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gift_raw)) {
      # Column containing the trait value is named after the trait ID
      gf_col <- grep("1\\.2\\.1", names(gift_raw), value = TRUE)[1L]
      gift_gf <- gift_raw |>
        dplyr::select(species = work_species,
                      growth_form_gift = dplyr::all_of(gf_col)) |>
        dplyr::filter(!is.na(growth_form_gift)) |>
        dplyr::mutate(growth_form_gift = tolower(trimws(growth_form_gift)))
      readr::write_csv(gift_gf, gift_gf_ckpt)
      ts(sprintf("  Saved GIFT checkpoint -> %s  (%d records)",
                 basename(gift_gf_ckpt), nrow(gift_gf)))
    } else {
      gift_gf <- NULL
    }
  }

  if (!is.null(gift_gf) && nrow(gift_gf) > 0L) {
    host_species_table <- host_species_table |>
      dplyr::left_join(gift_gf, by = "species") |>
      dplyr::mutate(
        growth_form = dplyr::if_else(is.na(growth_form), growth_form_gift, growth_form)
      ) |>
      dplyr::select(-growth_form_gift)

    n_imputed_gf <- n_na_gf - sum(is.na(host_species_table$growth_form))
    ts(sprintf("  Imputed: %d  |  Still NA: %d", n_imputed_gf,
               sum(is.na(host_species_table$growth_form))))
  }
}

# ---- Step 7b: Congener modal fallback for remaining NA growth_form ----------
# For any species still lacking a growth form, adopt the most common growth
# form among other species of the same genus already in the table.

n_na_gf_remaining <- sum(is.na(host_species_table$growth_form))
if (n_na_gf_remaining > 0L) {
  ts(sprintf("Step 7b: Congener-based growth form fallback for %d species...",
             n_na_gf_remaining))

  genus_modal_gf <- host_species_table |>
    dplyr::filter(!is.na(growth_form)) |>
    dplyr::mutate(genus = sub(" .*", "", species)) |>
    dplyr::count(genus, growth_form) |>
    dplyr::group_by(genus) |>
    dplyr::slice_max(n, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(genus, growth_form_genus = growth_form)

  host_species_table <- host_species_table |>
    dplyr::mutate(genus = sub(" .*", "", species)) |>
    dplyr::left_join(genus_modal_gf, by = "genus") |>
    dplyr::mutate(
      growth_form = dplyr::if_else(is.na(growth_form), growth_form_genus, growth_form)
    ) |>
    dplyr::select(-genus, -growth_form_genus)

  n_imputed_congener <- n_na_gf_remaining - sum(is.na(host_species_table$growth_form))
  ts(sprintf("  Imputed: %d  |  Still NA: %d", n_imputed_congener,
             sum(is.na(host_species_table$growth_form))))
}

# ---- Step 7c: Manual growth_form overrides for species unresolved by all ----
#               three automated sources (BIEN, GIFT, congener-modal fallback)
# Curated, species-level corrections for cases where growth_form remained NA
# because BIEN trait data, GIFT trait 1.2.1, and the genus-modal fallback all
# failed to resolve a value (e.g. no other Canadian congener in the table had
# a known growth form). Add new entries here as they surface; mirrors the
# `typo_fixes` lookup convention in canonicalize_host() (00_setup.R).
#   Saxifraga oppositifolia (purple saxifrage) — low cushion/mat-forming
#   arctic-alpine herb; unambiguous from species account, but no Saxifraga
#   congener in this host list carried a resolved growth form for the
#   genus-modal fallback to draw on (flagged 2026-06-28).
manual_growth_form <- c(
  "Saxifraga oppositifolia" = "herb"
)
n_na_gf_manual <- sum(is.na(host_species_table$growth_form))
if (n_na_gf_manual > 0L) {
  ix <- match(host_species_table$species, names(manual_growth_form))
  matched <- !is.na(ix) & is.na(host_species_table$growth_form)
  if (any(matched)) {
    host_species_table$growth_form[matched] <- unname(manual_growth_form[ix[matched]])
    ts(sprintf("Step 7c: Applied %d manual growth_form override(s). Still NA: %d",
               sum(matched), sum(is.na(host_species_table$growth_form))))
  }
}

# ---- Step 8: Save -----------------------------------------------------------

gf_summary <- sort(table(host_species_table$growth_form), decreasing = TRUE)
ts(sprintf("  Growth form breakdown: %s",
           paste(names(gf_summary), as.integer(gf_summary),
                 sep = " = ", collapse = ", ")))

readr::write_csv(host_species_table, paths$host_species)
ts(sprintf("Saved -> %s", basename(paths$host_species)))
ts(sprintf("Final host species count: %d", nrow(host_species_table)))
ts("06_host_species.R complete.")
