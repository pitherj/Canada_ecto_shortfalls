# =============================================================================
# EcM Native Host Species for Canada (BIEN-based)
# =============================================================================
# Assemble a list of EcM host plant species that are (a) demonstrated as EcM
# hosts by FungalRoot via EITHER of two independent routes (see NOTE below),
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
#   1.  Load the FungalRoot species-level table (occurrence route) AND the
#       Table S2 genus-level table (genus route) — see NOTE
#   2.  Query BIEN for all plant species with distributions in Canada
#   3.  Filter to native Canadian species using NSR (batched, 500 spp/request)
#   4.  Select hosts: native species satisfying the species rule directly,
#       OR whose genus qualifies under the Table S2 genus rule (always
#       recomputed fresh — see NOTE)
#   5.  Retrieve growth forms from BIEN trait database
#   6.  Join growth forms to host species table
#   7a. Impute missing growth_form from GIFT trait 1.2.1
#   7b. Congener-modal fallback for any growth_form still missing
#   7c. Curated manual overrides for any growth_form still missing after 7a/7b
#   8.  Save (with host_demonstrated flag and evidence_source — see NOTE)
#
# NOTE on EcM-host determination:
# Host status is now derived from TWO independent routes, combined with OR:
#   (a) SPECIES ROUTE: data_derived/clean_fungalroot_species.csv (produced by
#       05_prepare_fungalroot.R from GBIF DwC-A occurrence-level "Mycorrhiza
#       type" records). A species qualifies directly if FungalRoot carries
#       >=1 unambiguous EcM-positive occurrence record for it (see
#       EM_POSITIVE_LABELS in 05_prepare_fungalroot.R).
#   (b) GENUS ROUTE: data_derived/clean_fungalroot_genera_table_s2.csv
#       (produced by 05_prepare_fungalroot.R from FungalRoot's own published
#       Supplementary Table S2 genus-level recommendation). A species
#       qualifies if its genus is in this table (called "EcM" or "EcM-AM"),
#       even if that species has no occurrence-level record of its own.
# A native Canadian species is selected into this script's output if EITHER
# route applies. `evidence_source` records which route(s) actually applied
# ("occurrence", "table_s2", or "both"); `host_demonstrated` is TRUE for
# every row by construction (kept, rather than omitted, purely for
# self-documentation, mirroring the `ecm_demonstrated` convention in
# clean_fungalroot_species.csv).
#
# Checkpoint files (data_derived/checkpoints/):
#   bien_nsr_native_species.csv     — all BIEN Canada species with NSR
#                                      native_status == "N" (expensive;
#                                      independent of FungalRoot)
#   bien_ecm_growthforms.csv        — growth form trait data from BIEN
#   gift_growthforms.csv            — growth form trait data from GIFT
#                                     (used only for species lacking BIEN data)
#
# Output:
#   data_derived/ecm_native_canada_host_species.csv
#     Columns: species, host_demonstrated (always TRUE), evidence_source
#     ("occurrence" | "table_s2" | "both"), growth_form, growth_form_source
#     ("bien" | "gift" | "congener" | "manual" — which stage of the Step 5-7c
#     chain supplied the growth form)
#
# Runtime notes:
#   - Step 2 (BIEN_list_country) ~5 min
#   - Step 3 (NSR batched)       ~5–15 min (depends on server)
#   - Step 4 (FungalRoot select) seconds — always recomputed, not checkpointed
#   - Step 5 (BIEN traits)       ~10–30 min
#   Steps 2-3 and 5 are checkpointed; delete checkpoint files to force re-run.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(BIEN)
library(NSR)
library(GIFT)

sf::sf_use_s2(FALSE)

# Note: an earlier copy of clean_host_name() lived here but was never called;
# host-name cleaning is now centralized in canonicalize_host() (00_setup.R).

native_ckpt     <- file.path(paths$temp_dir, "bien_nsr_native_species.csv")
growthform_ckpt <- file.path(paths$temp_dir, "bien_ecm_growthforms.csv")
gift_gf_ckpt    <- file.path(paths$temp_dir, "gift_growthforms.csv")

# ---- Step 1: FungalRoot species table (occurrence route) + Table S2 genus
#              table (genus route) ---------------------------------------

ft_species <- readr::read_csv(paths$fungalroot_sp, show_col_types = FALSE)
demonstrated_species <- sort(unique(ft_species$UpdatedPlantBinomial))

ft_genera <- readr::read_csv(paths$fungalroot_genera, show_col_types = FALSE)
genus_qualifying_s2 <- sort(unique(ft_genera$Genus))

# ---- Steps 2-3: BIEN query + NSR native filter (expensive; checkpointed) ----
# Only the BIEN/NSR native-flora universe is checkpointed here. It is
# independent of FungalRoot, so it does not need to be recomputed just
# because 05_prepare_fungalroot.R's outputs change.

if (file.exists(native_ckpt)) {
  native_species <- readr::read_csv(native_ckpt, show_col_types = FALSE)$species
} else {

  bien_canada <- BIEN::BIEN_list_country("Canada", new.world = TRUE,
                                          cultivated = FALSE)
  bien_canada_species <- unique(bien_canada$scrubbed_species_binomial)

  NSR_BATCH <- 500L
  nsr_batches <- split(bien_canada_species,
                       ceiling(seq_along(bien_canada_species) / NSR_BATCH))
  nsr_list <- vector("list", length(nsr_batches))
  for (i in seq_along(nsr_batches)) {
    sp_batch <- nsr_batches[[i]]
    nsr_list[[i]] <- tryCatch(
      NSR::NSR_simple(species = sp_batch,
                      country = rep("Canada", length(sp_batch))),
      error = function(e) {
        warning(sprintf("  NSR batch %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }
  nsr_status <- dplyr::bind_rows(Filter(Negate(is.null), nsr_list))

  native_species <- nsr_status$species[
    !is.na(nsr_status$native_status) & nsr_status$native_status == "N"
  ]

  readr::write_csv(data.frame(species = native_species), native_ckpt)
}

# ---- Step 4: Select EcM hosts — species rule OR Table S2 genus rule --------
# Always recomputed fresh (cheap; no API calls) so the selection reflects
# whatever clean_fungalroot_species.csv and clean_fungalroot_genera_table_s2.csv
# currently contain, even when native_ckpt above is loaded from an old cache.

native_genus <- sub(" .*", "", native_species)

species_via_occurrence <- intersect(native_species, demonstrated_species)
species_via_genus      <- native_species[native_genus %in% genus_qualifying_s2]

em_canada_species <- sort(union(species_via_occurrence, species_via_genus))

evidence_lookup <- data.frame(species = em_canada_species) |>
  dplyr::mutate(
    via_occurrence = species %in% species_via_occurrence,
    via_genus      = species %in% species_via_genus,
    evidence_source = dplyr::case_when(
      via_occurrence & via_genus ~ "both",
      via_occurrence             ~ "occurrence",
      TRUE                       ~ "table_s2"
    )
  ) |>
  dplyr::select(species, evidence_source)

# ---- Step 5: BIEN growth form trait data ------------------------------------
# BIEN_trait_traitbyspecies() returns one row per trait record; multiple
# values may exist per species. We retain the modal value per species.

if (file.exists(growthform_ckpt)) {
  growthforms_raw <- readr::read_csv(growthform_ckpt, show_col_types = FALSE)
} else {

  # Batch to avoid timeouts; BIEN_trait_traitbyspecies handles a vector but
  # large vectors can time out.
  TRAIT_BATCH <- 200L
  trait_batches <- split(em_canada_species,
                         ceiling(seq_along(em_canada_species) / TRAIT_BATCH))
  trait_list <- vector("list", length(trait_batches))
  for (i in seq_along(trait_batches)) {
    trait_list[[i]] <- tryCatch(
      BIEN::BIEN_trait_traitbyspecies(
        species = trait_batches[[i]],
        trait   = "whole plant growth form"
      ) |>
        dplyr::select(scrubbed_species_binomial, trait_value) |>
        dplyr::mutate(trait_value = tolower(trait_value)),
      error = function(e) {
        warning(sprintf("  Trait batch %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }
  growthforms_raw <- dplyr::bind_rows(Filter(Negate(is.null), trait_list))
  readr::write_csv(growthforms_raw, growthform_ckpt)
}

# ---- Step 6: Summarise to one growth form per species -----------------------
# Use the most commonly reported trait value; ties broken arbitrarily by
# slice_max (first alphabetically among ties).

ecm_species_with_growthform <- growthforms_raw |>
  dplyr::count(scrubbed_species_binomial, trait_value) |>
  dplyr::group_by(scrubbed_species_binomial) |>
  dplyr::slice_max(n, n = 1L, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(species = scrubbed_species_binomial, growth_form = trait_value)

# ---- Step 7: Join to host species table -------------------------------------
# host_demonstrated is TRUE for every row by construction (selection in
# Step 4 already required the species rule OR the Table S2 genus rule to
# apply); it is computed explicitly here rather than just hard-coded, as a
# self-documenting sanity check (mirrors the `ecm_demonstrated` convention in
# clean_fungalroot_species.csv). evidence_source (from Step 4) records which
# route(s) actually applied for each species.

host_species_table <- data.frame(species = em_canada_species) |>
  dplyr::mutate(host_demonstrated = species %in% em_canada_species) |>
  dplyr::left_join(evidence_lookup, by = "species") |>
  dplyr::left_join(ecm_species_with_growthform, by = "species") |>
  # growth_form_source records WHICH of the four stages below supplied the
  # growth form, so the provenance of every value is auditable downstream.
  dplyr::mutate(growth_form_source = dplyr::if_else(is.na(growth_form),
                                                    NA_character_, "bien"))

# ---- Step 7a: Impute missing growth_form from GIFT trait 1.2.1 --------------
# BIEN traits did not cover all species. GIFT (trait_ID 1.2.1 = plant growth
# form) is used as a second source for remaining NAs. Checkpointed.

n_na_gf <- sum(is.na(host_species_table$growth_form))
if (n_na_gf > 0L) {

  if (file.exists(gift_gf_ckpt)) {
    gift_gf <- readr::read_csv(gift_gf_ckpt, show_col_types = FALSE)
  } else {
    gift_raw <- tryCatch(
      GIFT::GIFT_traits(trait_IDs  = "1.2.1",
                        agreement  = 0.66,
                        bias_ref   = FALSE,
                        bias_deriv = FALSE),
      error = function(e) {
        warning("  GIFT query failed: ", conditionMessage(e))
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
    } else {
      gift_gf <- NULL
    }
  }

  if (!is.null(gift_gf) && nrow(gift_gf) > 0L) {
    host_species_table <- host_species_table |>
      dplyr::left_join(gift_gf, by = "species") |>
      dplyr::mutate(
        growth_form_source = dplyr::if_else(
          is.na(growth_form) & !is.na(growth_form_gift), "gift", growth_form_source),
        growth_form = dplyr::if_else(is.na(growth_form), growth_form_gift, growth_form)
      ) |>
      dplyr::select(-growth_form_gift)

  }
}

# ---- Step 7b: Congener modal fallback for remaining NA growth_form ----------
# For any species still lacking a growth form, adopt the most common growth
# form among other species of the same genus already in the table.

n_na_gf_remaining <- sum(is.na(host_species_table$growth_form))
if (n_na_gf_remaining > 0L) {

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
      growth_form_source = dplyr::if_else(
        is.na(growth_form) & !is.na(growth_form_genus), "congener", growth_form_source),
      growth_form = dplyr::if_else(is.na(growth_form), growth_form_genus, growth_form)
    ) |>
    dplyr::select(-genus, -growth_form_genus)

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
#   genus-modal fallback to draw on.
manual_growth_form <- c(
  "Saxifraga oppositifolia" = "herb"
)
n_na_gf_manual <- sum(is.na(host_species_table$growth_form))
if (n_na_gf_manual > 0L) {
  ix <- match(host_species_table$species, names(manual_growth_form))
  matched <- !is.na(ix) & is.na(host_species_table$growth_form)
  if (any(matched)) {
    host_species_table$growth_form[matched] <- unname(manual_growth_form[ix[matched]])
    host_species_table$growth_form_source[matched] <- "manual"
  }
}

# ---- Step 8: Save -----------------------------------------------------------

readr::write_csv(host_species_table, paths$host_species)
