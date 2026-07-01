# =============================================================================
# Eltonian Shortfall (Updated): Host–Fungus Interaction Knowledge
# =============================================================================
#
# PART A — Canadian scope
#   What proportion of potential EcM host–fungus interactions in Canada have
#   been observed in our dataset?
#   Host species: BIEN-based native Canadian EcM hosts (06_host_species.R)
#   Fungal taxa:  EcM sequences from GenBank and GlobalFungi in our dataset
#
# PART B — Global scope (three new questions)
#   Q1: Among EcM host species native to Canada, how many have associated EcM
#       sequence records ANYWHERE in the world?
#   Q2: For each Canadian EcM host species, which and how many EcM fungal
#       species have been associated with it globally?
#   Q3: For each EcM fungal species found in Canada, which and how many EcM
#       host species has it been associated with globally?
#
#   Global data sources:
#     GlobalFungi: full sample metadata (78 MB) — host info from root samples
#                  full SH abundance matrix (13 GB) — fungal species per sample
#     GenBank:     global query for our EcM genera (no Canada filter)
#
# IMPORTANT NOTE on dominant_plant_species:
#   GlobalFungi samples have a 'dominant_plant_species' field recording which
#   plant species dominates the sample site. This can only be interpreted as
#   the likely mycorrhizal host when sample_type == "root" (i.e., the sample
#   is a direct root sample rather than a bulk soil or other substrate).
#   The same logic applies to 'other_plant_species'. This constraint is applied
#   throughout both Part A and Part B.
#
# Prerequisite files:
#   data_derived/ecm_native_canada_host_species.csv  (06_host_species.R)
#   data_raw/GlobalFungi/GlobalFungi_5_sample_metadata.txt    (full global)
#   data_raw/GlobalFungi/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt (full, 13 GB)
#   data_derived/unite_sh_taxonomy.csv               (01_globalfungi_prep.R)
#
# Outputs:
#   data_derived/eltonian_host_list.csv
#   data_derived/eltonian_host_matching.csv
#   data_derived/eltonian_matrix_genus.csv
#   data_derived/eltonian_matrix_sh.csv
#   data_derived/eltonian_matrix_species.csv         — host x named-EcM-species
#       binary matrix, trimmed to taxa with >= 1 observed pair (same
#       convention as the genus/SH matrices above)
#   data_derived/eltonian_species_occurrence_counts.csv — per host-species x
#       named-fungal-species pair, count of supporting samples/records
#       (occurrences); backs the singleton-association statistic in
#       eltonian_summary.csv
#   data_derived/eltonian_matrix_genus_genus.csv     — host-genus x fungal-genus
#       binary matrix (both axes collapsed to genus from the species-exact
#       `matched` interactions), trimmed to taxa with >= 1 observed pair,
#       same convention as the SH/species/genus matrices above
#   data_derived/eltonian_genus_occurrence_counts.csv — per host-genus x
#       fungal-genus pair, count of supporting samples/records (occurrences);
#       backs the genus-level singleton-association statistic in
#       eltonian_summary.csv
#   data_derived/eltonian_summary.csv
#   data_derived/eltonian_global_host_coverage.csv    — Q1: Canadian hosts with global data
#   data_derived/eltonian_global_host_to_fungi.csv    — Q2: per-host EcM fungal associations
#   data_derived/eltonian_global_fungi_to_hosts.csv   — Q3: per-fungal-species host associations
#   data_derived/eltonian/eltonian_sample_type_tally_canada.csv — diagnostic: sample_type
#       composition of EcM-positive Canadian GF samples with a
#       dominant_plant_species entry
#   data_derived/eltonian/eltonian_sample_type_tally_global.csv — diagnostic: sample_type
#       composition of EcM-positive global GF samples (>=1 of our Canadian
#       EcM SH codes detected) with a dominant_plant_species entry
#   data_derived/eltonian/eltonian_genbank_tissue_tally_canada.csv — diagnostic: tissue-type
#       composition of Canadian GenBank EcM records with host information
#       (host_taxon), keyword-binned from the free-text isolation_src field
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(rentrez)
library(data.table)

# Host-name cleaning is handled by the shared canonicalize_host() helper in
# 00_setup.R. GlobalFungi `dominant_plant_species` and `other_plant_species`
# values use underscores as token separators in some records, so we convert
# underscores to spaces before passing into canonicalize_host().
clean_host_name <- function(x) canonicalize_host(gsub("_", " ", x))

# =============================================================================
# PART A: Canadian scope
# =============================================================================

# ---- A1. Load Canadian EcM host species list ---------------------------------

ts("Part A — Step 1: Loading BIEN-based Canadian EcM host species list...")
if (!file.exists(paths$host_species)) {
  stop("Host species file not found: ", paths$host_species,
       "\nRun 06_host_species.R first.")
}
host_tbl <- readr::read_csv(paths$host_species, show_col_types = FALSE)
n_host_species <- nrow(host_tbl)
n_host_genera  <- dplyr::n_distinct(sub("^(\\S+).*", "\\1", host_tbl$species))
ts(sprintf("  Canadian EcM host species: %d  |  genera: %d",
           n_host_species, n_host_genera))
readr::write_csv(host_tbl, file.path(paths$out_eltonian, "eltonian_host_list.csv"))

# ---- A2. Extract and clean host strings from EMF dataset (Canadian records) --
# IMPORTANT: dominant_plant_species and other_plant_species are only reliable
# host indicators for root samples (sample_type == "root"). For soil and other
# sample types these fields record the dominant plant at the site but do NOT
# imply a direct mycorrhizal association. The host_taxon field (from GenBank
# host or isolation_source parsing) is used unconditionally.

ts("Part A — Step 2: Extracting host associations from Canadian EMF dataset...")

# Root samples only for GlobalFungi dominant/other plant species fields
emf_gf_root <- dplyr::filter(emf,
                               source == "GlobalFungi",
                               sample_type == "root")

ts(sprintf("  GlobalFungi root samples: %d of %d total GF records",
           nrow(emf_gf_root), sum(emf$source == "GlobalFungi")))

host_long <- dplyr::bind_rows(

  # dominant_plant_species: root samples only
  emf_gf_root |>
    dplyr::filter(!is.na(dominant_plant_species)) |>
    dplyr::select(sh_code, genus, species, dominant_plant_species) |>
    dplyr::rename(host_raw = dominant_plant_species) |>
    dplyr::mutate(host_field = "dominant_plant_species"),

  # other_plant_species: root samples only
  emf_gf_root |>
    dplyr::filter(!is.na(other_plant_species)) |>
    dplyr::select(sh_code, genus, species, other_plant_species) |>
    dplyr::rename(host_raw = other_plant_species) |>
    dplyr::mutate(host_field = "other_plant_species"),

  # host_taxon: GenBank field, usable unconditionally
  emf |>
    dplyr::filter(source == "GenBank", !is.na(host_taxon)) |>
    dplyr::select(sh_code, genus, species, host_taxon) |>
    dplyr::rename(host_raw = host_taxon) |>
    dplyr::mutate(host_field = "host_taxon")
)

ts(sprintf("  Total host strings extracted: %d", nrow(host_long)))

host_long <- host_long |>
  dplyr::mutate(
    host_clean  = clean_host_name(host_raw),
    host_genus  = sub("^(\\S+).*", "\\1", host_clean),
    matched     = host_clean %in% host_tbl$species,
    match_genus = host_genus %in% sub("^(\\S+).*", "\\1", host_tbl$species)
  )

n_matched       <- sum(host_long$matched,     na.rm = TRUE)
n_genus_matched <- sum(host_long$match_genus, na.rm = TRUE)
ts(sprintf("  Matched at species level: %d (%.1f%%)",
           n_matched, 100 * n_matched / nrow(host_long)))
ts(sprintf("  Matched at genus level:   %d (%.1f%%)",
           n_genus_matched, 100 * n_genus_matched / nrow(host_long)))

readr::write_csv(host_long,
                 file.path(paths$out_eltonian, "eltonian_host_matching.csv"))

# ---- A2b. Diagnostic: sample_type tally among Canadian GF samples with a -----
#           dominant_plant_species entry
# This does NOT change the root-only filtering used above (A2) or anywhere
# else in the script -- it is purely descriptive, to quantify how much
# dominant_plant_species data exist for sample types other than "root" (e.g.
# soil) that are excluded from host matching because, per the note at the top
# of this script, dominant_plant_species/other_plant_species are only
# interpretable as a likely host for direct root samples.
#
# Counted at the distinct-sample level (sample_ID), not row level, because
# `emf` carries one row per (sample x SH code): a sample with several
# co-occurring EcM detections would otherwise be tallied once per detection.

ts("Part A — Step 2b: Tallying sample_type among Canadian GF samples with a dominant_plant_species entry...")

sample_type_tally_canada <- emf |>
  dplyr::filter(source == "GlobalFungi", !is.na(dominant_plant_species)) |>
  dplyr::distinct(sample_ID, sample_type) |>
  dplyr::count(sample_type, name = "n_samples") |>
  dplyr::arrange(dplyr::desc(n_samples))

ts("  Canada sample_type tally (distinct samples with dominant_plant_species):")
print(as.data.frame(sample_type_tally_canada))

readr::write_csv(sample_type_tally_canada,
                 file.path(paths$out_eltonian, "eltonian_sample_type_tally_canada.csv"))
ts("  Saved eltonian_sample_type_tally_canada.csv")

# ---- A2c. Diagnostic: tissue-type tally among Canadian GenBank records with -
#           host information (host_taxon)
# Mirrors A2b above, but for GenBank. Unlike GlobalFungi, GenBank has no
# controlled-vocabulary sample_type field -- the closest analogue is the
# free-text 'isolation_src' field (NCBI /isolation_source qualifier), which is
# heterogeneous (56 distinct strings in this dataset, e.g. "root system",
# "ectomycorrhiza", "soil adjacent to Pinus albicaulis", "40 year-old
# Douglas-fir stand"). Per the note at the top of this script, host_taxon is
# used UNCONDITIONALLY for GenBank host matching (no tissue-type filter,
# unlike the root-only restriction applied to GlobalFungi). This diagnostic
# quantifies what tissue type actually underlies that unconditional host
# attribution, among records that have host information.
#
# Tissue categories are assigned by keyword search (case-insensitive) on
# isolation_src, in this priority order:
#   Root / ectomycorrhizal tissue  - matches root|mycorrhiz|ecm (whole word)
#                                     and does NOT also match a soil keyword
#   Soil / rhizosphere             - matches soil|rhizosphere, no root keyword
#   Mixed root + soil              - matches BOTH a root and a soil keyword
#                                     (e.g. "pine roots from urban soil" --
#                                     note this mechanical rule will also tag
#                                     some root-only collections that merely
#                                     mention the soil habitat; treat as an
#                                     approximate, documented limitation)
#   Other / non-tissue description - text present but no root/soil keyword
#                                     (e.g. stand age, forest type, locality)
#   Not recorded                   - isolation_src missing entirely
#
# Counted at the record level (one row = one GenBank accession; confirmed
# 1:1 in this dataset, unlike the long-format GF rows in A2b).

ts("Part A — Step 2c: Tallying tissue-type composition of Canadian GenBank records with host_taxon...")

genbank_canada <- dplyr::filter(emf, source == "GenBank")
n_gb_total <- nrow(genbank_canada)

genbank_with_host <- dplyr::filter(genbank_canada, !is.na(host_taxon))
n_gb_with_host <- nrow(genbank_with_host)

root_kw <- "root|mycorrhiz|\\becm\\b"
soil_kw <- "soil|rhizosphere"

genbank_tissue <- genbank_with_host |>
  dplyr::mutate(
    has_root = grepl(root_kw, isolation_src, ignore.case = TRUE),
    has_soil = grepl(soil_kw, isolation_src, ignore.case = TRUE),
    tissue_category = dplyr::case_when(
      is.na(isolation_src)        ~ "Not recorded",
      has_root & has_soil         ~ "Mixed root + soil",
      has_root                    ~ "Root / ectomycorrhizal tissue",
      has_soil                    ~ "Soil / rhizosphere",
      TRUE                        ~ "Other / non-tissue description"
    )
  )

genbank_tissue_tally_canada <- dplyr::bind_rows(
  tibble::tibble(category = "Total EcM fungal records",      n = n_gb_total),
  tibble::tibble(category = "Records with host information", n = n_gb_with_host),
  genbank_tissue |>
    dplyr::count(tissue_category, name = "n") |>
    dplyr::rename(category = tissue_category) |>
    dplyr::arrange(dplyr::desc(n))
)

ts("  GenBank tissue-type tally (Canada):")
print(as.data.frame(genbank_tissue_tally_canada))

readr::write_csv(genbank_tissue_tally_canada,
                 file.path(paths$out_eltonian, "eltonian_genbank_tissue_tally_canada.csv"))
ts("  Saved eltonian_genbank_tissue_tally_canada.csv")

# ---- A3. Build interaction matrices ------------------------------------------

ts("Part A — Step 3: Building interaction matrices...")

matched_interactions <- dplyr::filter(host_long, matched)

genus_pairs <- matched_interactions |>
  dplyr::distinct(host_clean, genus) |>
  dplyr::rename(host_species = host_clean, fungal_genus = genus)

genus_matrix <- genus_pairs |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_species, names_from = fungal_genus,
                     values_from = present, values_fill = 0L)
readr::write_csv(genus_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_genus.csv"))

sh_pairs <- matched_interactions |>
  dplyr::distinct(host_clean, sh_code) |>
  dplyr::rename(host_species = host_clean)

sh_matrix <- sh_pairs |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_species, names_from = sh_code,
                     values_from = present, values_fill = 0L)
readr::write_csv(sh_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_sh.csv"))
ts("  Saved eltonian_matrix_genus.csv and eltonian_matrix_sh.csv")

# ---- A3b. Species-level Canada pairs ----------------------------------------
# Join sh_pairs with UNITE taxonomy to resolve SH codes → named species
ts("Part A — Step 3b: Building species-level Canada pairs (SH → named species)...")
sh_lookup <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
species_pairs_canada <- sh_pairs |>
  dplyr::left_join(
    sh_lookup |> dplyr::select(sh_code, species),
    by = "sh_code"
  ) |>
  dplyr::mutate(
    fungal_species = dplyr::if_else(
      !is.na(species) & !grepl("_sp$", species),
      trimws(gsub("_", " ", species)),
      NA_character_
    )
  ) |>
  dplyr::filter(!is.na(fungal_species)) |>
  dplyr::distinct(host_species, fungal_species)
ts(sprintf("  Species-level Canada pairs: %d (host × named EcM species)",
           nrow(species_pairs_canada)))

# ---- A3c. Species-level matrix, fill statistics, and occurrence counts -------
# Trimmed host x named-EcM-species matrix, following the same convention as
# eltonian_matrix_genus.csv / eltonian_matrix_sh.csv above: rows and columns
# only for taxa with >= 1 observed pair, NOT the full host x species
# denominator. (Confirmed with Jason 2026-06-27: the full 147 x 1079 grid is
# reported as a numeric fill statistic below rather than materialized as a
# mostly-empty CSV.)

ts("Part A — Step 3c: Building species-level matrix and occurrence/fill statistics...")

species_matrix <- species_pairs_canada |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_species, names_from = fungal_species,
                     values_from = present, values_fill = 0L)
readr::write_csv(species_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_species.csv"))
ts("  Saved eltonian_matrix_species.csv")

# Full-matrix fill statistics. Unlike the matrix CSV above (trimmed to taxa
# with data), the fill rate is calculated against the FULL potential host x
# named-fungal-species grid: all BIEN-based Canadian EcM host species
# (n_host_species) x all named EcM fungal species detected anywhere in the
# Canadian dataset (species not ending "_sp" in `emf$species`), regardless of
# whether a given fungal species has any host information at all.
all_named_fungal_species <- emf |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::mutate(fungal_species = trimws(gsub("_", " ", species))) |>
  dplyr::distinct(fungal_species) |>
  dplyr::pull(fungal_species)
n_named_fungal_species <- length(all_named_fungal_species)

n_matrix_cells   <- n_host_species * n_named_fungal_species
n_cells_filled   <- nrow(species_pairs_canada)
n_cells_empty    <- n_matrix_cells - n_cells_filled
pct_cells_filled <- round(100 * n_cells_filled / n_matrix_cells, 3)
pct_cells_empty  <- round(100 * n_cells_empty  / n_matrix_cells, 3)

ts(sprintf("  Full host x named-species matrix: %d x %d = %s cells",
           n_host_species, n_named_fungal_species,
           format(n_matrix_cells, big.mark = ",")))
ts(sprintf("  Filled: %s (%.3f%%)  |  Empty: %s (%.3f%%)",
           format(n_cells_filled, big.mark = ","), pct_cells_filled,
           format(n_cells_empty,  big.mark = ","), pct_cells_empty))

# Per-pair occurrence counts: how many supporting samples/records (rows of
# matched_interactions, i.e. distinct GlobalFungi sample x SH-detection rows
# or GenBank accession rows) underlie each observed host-species x
# fungal-species pair. Built from matched_interactions directly -- NOT from
# the already-deduplicated sh_pairs / species_pairs_canada -- so occurrence
# multiplicity is preserved. Uses the same sh_code -> fungal_species
# resolution (via sh_lookup) as species_pairs_canada above, for consistency.
species_occurrence_counts <- matched_interactions |>
  dplyr::select(-species) |>
  dplyr::left_join(sh_lookup |> dplyr::select(sh_code, species), by = "sh_code") |>
  dplyr::mutate(
    fungal_species = dplyr::if_else(
      !is.na(species) & !grepl("_sp$", species),
      trimws(gsub("_", " ", species)),
      NA_character_
    )
  ) |>
  dplyr::filter(!is.na(fungal_species)) |>
  dplyr::rename(host_species = host_clean) |>
  dplyr::count(host_species, fungal_species, name = "n_occurrences") |>
  dplyr::arrange(dplyr::desc(n_occurrences))

# Sanity check: this should resolve to exactly the same set of distinct pairs
# as species_pairs_canada (same sh_code -> fungal_species join, same source
# rows collapsed to the same keys). A mismatch would indicate the two code
# paths have diverged.
if (nrow(species_occurrence_counts) != nrow(species_pairs_canada)) {
  warning(sprintf(
    paste("species_occurrence_counts (%d pairs) does not match",
          "species_pairs_canada (%d pairs) -- check join logic."),
    nrow(species_occurrence_counts), nrow(species_pairs_canada)))
}

readr::write_csv(species_occurrence_counts,
                 file.path(paths$out_eltonian, "eltonian_species_occurrence_counts.csv"))
ts("  Saved eltonian_species_occurrence_counts.csv")

n_pairs_total       <- nrow(species_occurrence_counts)
n_pairs_singleton   <- sum(species_occurrence_counts$n_occurrences == 1L)
pct_pairs_singleton <- round(100 * n_pairs_singleton / n_pairs_total, 1)

ts(sprintf("  Host x species pairs supported by a single occurrence: %d / %d (%.1f%%)",
           n_pairs_singleton, n_pairs_total, pct_pairs_singleton))

# ---- A3d. Genus-level matrix, fill statistics, and occurrence counts --------
# Genus x genus analogue of the host x named-species block above (A3c).
# Both axes are collapsed to genus: the host axis from host_clean (species)
# to host_genus, the fungal axis is already genus-resolution (emf$genus).
# Uses matched_interactions (i.e. the SAME species-exact host match, `matched`,
# used to build genus_pairs/sh_pairs/species_pairs_canada above) collapsed to
# host_genus -- NOT the looser `match_genus` flag in host_long -- so the
# matching criterion stays consistent across the SH-, species-, and
# genus-resolution blocks. (Confirmed with Jason 2026-06-28: the looser
# match_genus criterion gives slightly different counts -- 12 host genera x
# 69 fungal genera, 263 filled cells -- because it admits host records whose
# species name didn't resolve to the BIEN list but whose genus did; the
# stricter, already-established `matched` criterion is used here instead.)
#
# As with the species-level matrix, the saved CSV is trimmed to taxa with
# >= 1 observed pair; the full-grid fill rate is reported as a numeric
# statistic in eltonian_summary.csv rather than materialized as a mostly-
# empty CSV.

ts("Part A — Step 3d: Building genus-level matrix and occurrence/fill statistics...")

genus_genus_pairs_canada <- matched_interactions |>
  dplyr::distinct(host_genus, genus) |>
  dplyr::rename(fungal_genus = genus)
ts(sprintf("  Genus-level Canada pairs: %d (host genus x fungal genus)",
           nrow(genus_genus_pairs_canada)))

genus_genus_matrix <- genus_genus_pairs_canada |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_genus, names_from = fungal_genus,
                     values_from = present, values_fill = 0L)
readr::write_csv(genus_genus_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_genus_genus.csv"))
ts("  Saved eltonian_matrix_genus_genus.csv")

# Full-matrix fill statistics, against the FULL potential host-genus x
# fungal-genus grid: all BIEN-based Canadian EcM host genera (n_host_genera)
# x all EcM fungal genera detected anywhere in the Canadian dataset
# (n_fungal_genera_total), regardless of whether a given fungal genus has
# any host information at all.
n_fungal_genera_total <- dplyr::n_distinct(emf$genus)

n_matrix_cells_genus   <- n_host_genera * n_fungal_genera_total
n_cells_filled_genus   <- nrow(genus_genus_pairs_canada)
n_cells_empty_genus    <- n_matrix_cells_genus - n_cells_filled_genus
pct_cells_filled_genus <- round(100 * n_cells_filled_genus / n_matrix_cells_genus, 2)
pct_cells_empty_genus  <- round(100 * n_cells_empty_genus  / n_matrix_cells_genus, 2)

ts(sprintf("  Full host-genus x fungal-genus matrix: %d x %d = %s cells",
           n_host_genera, n_fungal_genera_total,
           format(n_matrix_cells_genus, big.mark = ",")))
ts(sprintf("  Filled: %s (%.2f%%)  |  Empty: %s (%.2f%%)",
           format(n_cells_filled_genus, big.mark = ","), pct_cells_filled_genus,
           format(n_cells_empty_genus,  big.mark = ","), pct_cells_empty_genus))

# Per-pair occurrence counts: how many supporting samples/records underlie
# each observed host-genus x fungal-genus pair. Built from matched_interactions
# directly (preserves occurrence multiplicity), mirroring
# eltonian_species_occurrence_counts.csv above.
genus_occurrence_counts <- matched_interactions |>
  dplyr::count(host_genus, genus, name = "n_occurrences") |>
  dplyr::rename(fungal_genus = genus) |>
  dplyr::arrange(dplyr::desc(n_occurrences))

readr::write_csv(genus_occurrence_counts,
                 file.path(paths$out_eltonian, "eltonian_genus_occurrence_counts.csv"))
ts("  Saved eltonian_genus_occurrence_counts.csv")

n_pairs_total_genus       <- nrow(genus_occurrence_counts)
n_pairs_singleton_genus   <- sum(genus_occurrence_counts$n_occurrences == 1L)
pct_pairs_singleton_genus <- round(100 * n_pairs_singleton_genus / n_pairs_total_genus, 1)

ts(sprintf("  Host-genus x fungal-genus pairs supported by a single occurrence: %d / %d (%.1f%%)",
           n_pairs_singleton_genus, n_pairs_total_genus, pct_pairs_singleton_genus))

# ---- A4. Canadian-scope coverage statistics ----------------------------------

ts("Part A — Step 4: Calculating Canadian-scope coverage statistics...")

n_obs_genus_pairs       <- nrow(genus_pairs)
n_obs_sh_pairs          <- nrow(sh_pairs)
n_hosts_with_genus_data <- dplyr::n_distinct(genus_pairs$host_species)
n_hosts_with_sh_data    <- dplyr::n_distinct(sh_pairs$host_species)
n_genera_with_host      <- dplyr::n_distinct(genus_pairs$fungal_genus)
n_sh_with_host          <- dplyr::n_distinct(sh_pairs$sh_code)
n_potential_genus       <- n_host_species * dplyr::n_distinct(emf$genus)
n_potential_sh          <- n_host_species * dplyr::n_distinct(emf$sh_code)
genus_per_host          <- dplyr::count(genus_pairs, host_species, name = "n_genera")
host_per_genus          <- dplyr::count(genus_pairs, fungal_genus, name = "n_hosts")

# =============================================================================
# PART B: Global scope
# =============================================================================
# Paths to full GlobalFungi data (not just Canadian subset)
gf_meta_path <- file.path(here::here("data_raw"), "GlobalFungi",
                            "GlobalFungi_5_sample_metadata.txt")
gf_sh_path   <- file.path(here::here("data_raw"), "GlobalFungi",
                            "GlobalFungi_5_SH_abundance_ITS1_ITS2.txt")
unite_lookup_path <- paths$unite_taxonomy

# Checkpoints for slow operations
global_meta_ckpt <- file.path(paths$temp_dir, "gf_global_root_metadata.csv")
global_sh_ckpt   <- file.path(paths$temp_dir, "gf_global_ecm_sh_subset.rds")
global_gb_ckpt   <- file.path(paths$temp_dir, "genbank_global_ecm_meta.csv")

# ---- B1. GlobalFungi: global root-sample metadata ---------------------------
# Q1 (metadata-only approach): which Canadian host species appear as
# dominant_plant_species in root samples anywhere in the full GlobalFungi DB?

ts("Part B — Step 1: Reading GlobalFungi global root sample metadata...")

if (!file.exists(gf_meta_path)) {
  ts(sprintf("  GlobalFungi metadata not found at: %s", gf_meta_path))
  ts("  Skipping GlobalFungi global analysis.")
  gf_global_root <- NULL
} else if (file.exists(global_meta_ckpt)) {
  ts("  Loading checkpointed GlobalFungi root sample metadata...")
  gf_global_root <- readr::read_csv(global_meta_ckpt, show_col_types = FALSE)
  ts(sprintf("  Root samples loaded: %d", nrow(gf_global_root)))
} else {
  ts("  Reading GlobalFungi global metadata (~78 MB)...")
  # Use data.table::fread for speed on large file
  gf_meta_all <- data.table::fread(
    gf_meta_path,
    sep            = "\t",
    quote          = "",
    select         = c("sample_ID", "country", "sample_type",
                       "dominant_plant_species", "other_plant_species"),
    showProgress   = TRUE
  )
  ts(sprintf("  Total GlobalFungi samples: %d", nrow(gf_meta_all)))

  # Filter to root samples only
  gf_global_root <- dplyr::filter(gf_meta_all, sample_type == "root")
  ts(sprintf("  Root samples globally: %d (%.1f%%)",
             nrow(gf_global_root),
             100 * nrow(gf_global_root) / nrow(gf_meta_all)))
  rm(gf_meta_all)

  readr::write_csv(gf_global_root, global_meta_ckpt)
  ts(sprintf("  Saved root metadata checkpoint -> %s",
             basename(global_meta_ckpt)))
}

# ---- B2. Q1: Which Canadian host species have global GF root records? --------

if (!is.null(gf_global_root) && nrow(gf_global_root) > 0) {
  ts("Part B — Step 2 (Q1): Matching Canadian hosts against global GF root samples...")

  # Extract all host strings from global root samples (dominant + other)
  global_hosts <- c(
    clean_host_name(gf_global_root$dominant_plant_species[
      !is.na(gf_global_root$dominant_plant_species)]),
    clean_host_name(gf_global_root$other_plant_species[
      !is.na(gf_global_root$other_plant_species)])
  )
  global_hosts <- unique(global_hosts[!is.na(global_hosts)])

  # Match against our Canadian host list
  canada_hosts_with_global_gf <- intersect(host_tbl$species, global_hosts)
  n_canada_hosts_with_gf <- length(canada_hosts_with_global_gf)
  pct_with_gf <- round(100 * n_canada_hosts_with_gf / n_host_species, 1)
  ts(sprintf("  Canadian hosts with GF root records globally: %d / %d (%.1f%%)",
             n_canada_hosts_with_gf, n_host_species, pct_with_gf))

  # Build sample-host lookup for the next steps
  # (sample_ID → cleaned dominant_plant_species) for root samples only
  gf_root_hosts <- dplyr::bind_rows(
    gf_global_root |>
      dplyr::filter(!is.na(dominant_plant_species)) |>
      dplyr::transmute(sample_ID,
                       host_clean = clean_host_name(dominant_plant_species)),
    gf_global_root |>
      dplyr::filter(!is.na(other_plant_species)) |>
      dplyr::transmute(sample_ID,
                       host_clean = clean_host_name(other_plant_species))
  ) |>
    dplyr::filter(!is.na(host_clean)) |>
    dplyr::distinct()

  # Filter to sample_IDs where host is a Canadian EcM host species
  gf_root_canadian_hosts <- dplyr::filter(gf_root_hosts,
                                           host_clean %in% host_tbl$species)
  ts(sprintf("  Global root sample-host pairs (Canadian hosts): %d",
             nrow(gf_root_canadian_hosts)))

  # Save Q1 result
  q1_host_coverage <- data.frame(
    species         = host_tbl$species,
    has_global_gf   = host_tbl$species %in% canada_hosts_with_global_gf
  )
  readr::write_csv(q1_host_coverage,
                   file.path(paths$out_eltonian, "eltonian_global_host_coverage.csv"))
  ts("  Saved eltonian_global_host_coverage.csv")

} else {
  ts("  GlobalFungi global metadata not available; skipping Q1.")
  gf_root_canadian_hosts <- NULL
  n_canada_hosts_with_gf <- NA_integer_
  pct_with_gf <- NA_real_
}

# ---- B3. Q2 and Q3: GlobalFungi SH abundance (global, 13 GB) ----------------
# This step requires reading the full 13 GB SH abundance matrix.
# Strategy: select only the SH code columns corresponding to our Canadian EcM
# taxa, filter to root samples, join with metadata for host info.
# This approach efficiently answers Q3 (Canadian EcM fungi → global hosts).
# Q2 (Canadian hosts → global EcM fungi) is derived from the same join.
#
# Memory note: loading just our EcM SH codes (a few hundred columns) from a
# very wide file is much more efficient than loading all 77k+ SH columns.

ts("Part B — Step 3 (Q2 + Q3): GlobalFungi global SH abundance analysis...")

if (!file.exists(gf_sh_path)) {
  ts(sprintf("  GlobalFungi SH abundance file not found: %s", gf_sh_path))
  ts("  Skipping Q2/Q3 GlobalFungi analysis.")
  gf_sh_ecm <- NULL

} else if (file.exists(global_sh_ckpt)) {
  ts("  Loading checkpointed GF global EcM SH subset...")
  gf_sh_ecm <- readRDS(global_sh_ckpt)
  ts(sprintf("  Loaded: %d rows", nrow(gf_sh_ecm)))

} else {
  # Identify our Canadian EcM SH codes (from emf dataset)
  our_sh_codes <- unique(emf$sh_code[!is.na(emf$sh_code)])
  ts(sprintf("  Our Canadian EcM SH codes: %d", length(our_sh_codes)))

  # Always check the file header first to find the intersection of our SH codes
  # with the file's columns. Our codes use a specific UNITE version suffix
  # (e.g., .10FU) which may differ from the version used in GlobalFungi.
  ts("  Reading GlobalFungi file header to check SH code versions...")
  header_cols <- names(data.table::fread(gf_sh_path, sep = "\t", quote = "",
                                          nrows = 0L))
  present_sh <- intersect(our_sh_codes, header_cols)
  ts(sprintf("  SH codes matched in GlobalFungi file: %d / %d",
             length(present_sh), length(our_sh_codes)))

  if (length(present_sh) == 0L) {
    # Version mismatch: strip version suffix and try prefix matching
    # e.g., "SH1052460.10FU" -> "SH1052460"
    our_sh_prefix   <- sub("\\.[0-9]+FU$", "", our_sh_codes)
    header_prefix   <- sub("\\.[0-9]+FU$", "", header_cols)
    matched_idx     <- match(our_sh_prefix, header_prefix)
    matched_header  <- header_cols[matched_idx[!is.na(matched_idx)]]
    present_sh      <- matched_header
    ts(sprintf("  Version-stripped SH matches: %d / %d",
               length(present_sh), length(our_sh_codes)))
  }

  if (length(present_sh) == 0L) {
    ts("  No SH codes matched in GlobalFungi file — skipping Q2/Q3 (GF).")
    gf_sh_ecm_wide <- NULL
  } else {
    ts("  Reading GlobalFungi SH abundance matrix (this may take 10-30 min)...")
    ts("  Selecting only matched EcM SH columns from the full 13 GB file...")
    gf_sh_ecm_wide <- data.table::fread(
      gf_sh_path,
      sep          = "\t",
      quote        = "",
      select       = c("sample_ID", present_sh),
      showProgress = TRUE
    )
  }

  ts(sprintf("  Read: %d samples x %d EcM SH columns",
             nrow(gf_sh_ecm_wide), ncol(gf_sh_ecm_wide) - 1))

  # Melt to long format, filter to non-zero abundances
  ts("  Melting to long format...")
  gf_sh_ecm <- data.table::melt(
    gf_sh_ecm_wide,
    id.vars       = "sample_ID",
    variable.name = "sh_code",
    value.name    = "abundance",
    variable.factor = FALSE
  ) |>
    dplyr::filter(abundance > 0L) |>
    as.data.frame()
  rm(gf_sh_ecm_wide)

  ts(sprintf("  Non-zero abundance records: %d", nrow(gf_sh_ecm)))
  saveRDS(gf_sh_ecm, global_sh_ckpt)
  ts(sprintf("  Saved SH subset checkpoint -> %s", basename(global_sh_ckpt)))
}

# ---- B3b. Diagnostic: sample_type tally among EcM-positive global GF -------
#           samples with a dominant_plant_species entry
# Mirrors A2b above (Canadian scope), but for the global scope. The Canadian
# tally is implicitly EcM-only because it is built from `emf`, which is
# already filtered to ectomycorrhizal fungal records (FungalTraits
# primary_lifestyle == "ectomycorrhizal"; see emf_canada_em_only.csv). For
# the two tallies to be directly comparable, this global tally is restricted
# the same way: to global samples carrying >=1 non-zero detection of one of
# our Canadian EcM SH codes. `gf_sh_ecm` (built just above in Part B Step 3)
# already provides exactly that sample set, so this step reuses it rather
# than triggering a second full 13 GB matrix scan. Sample-type and
# dominant_plant_species values still require a fresh (cheap, ~78 MB) read
# of the metadata file, joined here on sample_ID.

global_sampletype_tally_out <- file.path(paths$out_eltonian,
                                         "eltonian_sample_type_tally_global.csv")

ts("Part B — Step 3b: Tallying sample_type among EcM-positive global GF samples with a dominant_plant_species entry...")

if (file.exists(global_sampletype_tally_out)) {
  ts("  Loading existing global sample_type tally...")
  sample_type_tally_global <- readr::read_csv(global_sampletype_tally_out,
                                              show_col_types = FALSE)
} else if (is.null(gf_sh_ecm) || !file.exists(gf_meta_path)) {
  ts("  Global EcM SH subset or metadata not available; skipping global sample_type tally.")
  sample_type_tally_global <- NULL
} else {
  ecm_positive_samples <- unique(gf_sh_ecm$sample_ID)
  ts(sprintf("  Global samples with >=1 of our Canadian EcM SH codes detected: %d",
             length(ecm_positive_samples)))

  ts("  Reading sample_ID + sample_type + dominant_plant_species from full GlobalFungi metadata...")
  gf_meta_minimal <- data.table::fread(
    gf_meta_path,
    sep    = "\t",
    quote  = "",
    select = c("sample_ID", "sample_type", "dominant_plant_species")
  )

  sample_type_tally_global <- gf_meta_minimal |>
    dplyr::filter(sample_ID %in% ecm_positive_samples,
                  !is.na(dominant_plant_species)) |>
    dplyr::distinct(sample_ID, sample_type) |>
    dplyr::count(sample_type, name = "n_samples") |>
    dplyr::arrange(dplyr::desc(n_samples))
  rm(gf_meta_minimal)

  readr::write_csv(sample_type_tally_global, global_sampletype_tally_out)
  ts(sprintf("  Saved %s", basename(global_sampletype_tally_out)))
}

if (!is.null(sample_type_tally_global)) {
  ts("  Global sample_type tally (EcM-positive samples with dominant_plant_species):")
  print(as.data.frame(sample_type_tally_global))
}

# Q2 and Q3 from GlobalFungi
if (!is.null(gf_sh_ecm) && !is.null(gf_root_canadian_hosts)) {

  ts("  Computing Q2 and Q3 from GlobalFungi global data...")
  # (sh_lookup already loaded in Part A — Step 3b)

  # Attach UNITE species to each SH code observation
  gf_sh_ecm_tax <- dplyr::left_join(gf_sh_ecm, sh_lookup, by = "sh_code") |>
    dplyr::mutate(
      fungal_species = dplyr::if_else(!is.na(species) & !grepl("_sp$", species),
                                       trimws(gsub("_", " ", species)),
                                       NA_character_)
    )

  # Q3: For each Canadian EcM fungal species, which hosts globally (root samples)?
  # NOTE on relationship = "many-to-many": sample_ID is legitimately duplicated
  # on both sides of every join below against gf_root_canadian_hosts.
  # gf_sh_ecm_tax / gf_sh_ecm have one row per (sample_ID, sh_code) with
  # non-zero abundance, so a sample with several co-occurring EcM fungi
  # contributes several rows. gf_root_canadian_hosts can also have >1 row per
  # sample_ID when both dominant_plant_species and other_plant_species
  # resolve to (different) Canadian EcM host species. Each join below
  # therefore enumerates every fungus x host pair co-occurring at a sample —
  # the intended co-occurrence design for this shortfall — and is always
  # followed immediately by distinct(), so the many-to-many expansion never
  # inflates a downstream count.
  q3 <- gf_sh_ecm_tax |>
    dplyr::filter(!is.na(fungal_species)) |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::distinct(fungal_species, host_clean) |>
    dplyr::rename(host_species = host_clean) |>
    dplyr::group_by(fungal_species) |>
    dplyr::summarise(
      n_host_species = dplyr::n_distinct(host_species),
      host_species   = paste(sort(unique(host_species)), collapse = "; "),
      .groups        = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_host_species))

  readr::write_csv(q3,
                   file.path(paths$out_eltonian, "eltonian_global_fungi_to_hosts.csv"))
  ts(sprintf("  Q3: EcM species with global host associations: %d", nrow(q3)))
  ts("  Saved eltonian_global_fungi_to_hosts.csv")

  # Q2: For each Canadian host species, which EcM fungal species globally (root)?
  # many-to-many expected here too (see NOTE above Q3) — distinct() follows.
  q2 <- gf_sh_ecm_tax |>
    dplyr::filter(!is.na(fungal_species)) |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::distinct(host_clean, fungal_species) |>
    dplyr::rename(host_species = host_clean) |>
    dplyr::filter(host_species %in% host_tbl$species) |>
    dplyr::group_by(host_species) |>
    dplyr::summarise(
      n_ecm_species  = dplyr::n_distinct(fungal_species),
      ecm_species    = paste(sort(unique(fungal_species)), collapse = "; "),
      .groups        = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_ecm_species))

  readr::write_csv(q2,
                   file.path(paths$out_eltonian, "eltonian_global_host_to_fungi.csv"))
  ts(sprintf("  Q2: Canadian hosts with global EcM species data: %d", nrow(q2)))
  ts("  Saved eltonian_global_host_to_fungi.csv")

  # ---- Long-form pair data frames for the interactions Rds (Part C) ---------
  ts("  Building long-form global pair data frames (GF)...")

  # Host → SH code (global, GF root samples, Canadian hosts only)
  # many-to-many expected here too (see NOTE above Q3) — distinct() follows.
  host_sh_global_gf <- gf_sh_ecm |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::filter(host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, sh_code)
  ts(sprintf("  host_sh_global_gf:      %d host–SH pairs", nrow(host_sh_global_gf)))

  # Host → named EcM species (global, GF)
  # many-to-many expected here too (see NOTE above Q3) — distinct() follows.
  host_species_global_gf <- gf_sh_ecm_tax |>
    dplyr::filter(!is.na(fungal_species)) |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::filter(host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, fungal_species)
  ts(sprintf("  host_species_global_gf: %d host–species pairs",
             nrow(host_species_global_gf)))

  # Host → genus (global, GF, derived from species)
  host_genus_global_gf <- host_species_global_gf |>
    dplyr::mutate(fungal_genus = sub(" .*", "", fungal_species)) |>
    dplyr::distinct(host_species, fungal_genus)
  ts(sprintf("  host_genus_global_gf:   %d host–genus pairs",
             nrow(host_genus_global_gf)))

} else {
  ts("  GlobalFungi global SH data not available; skipping Q2/Q3 (GF).")
  q2 <- NULL
  q3 <- NULL
  host_sh_global_gf      <- NULL
  host_species_global_gf <- NULL
  host_genus_global_gf   <- NULL
}

# ---- B4. Global GenBank: EcM genera queried without Canada filter -----------
# For each of our EcM genera, query GenBank globally (no country constraint)
# and retrieve host_taxon metadata. This answers Q3 from the GenBank side
# (and partially Q1/Q2 by the host_taxon field).
#
# Each genus is searched separately and results are combined.
# The search mirrors 03_genbank_fetch.R but removes AND "Canada"[Country].

ts("Part B — Step 4 (GenBank global): Querying GenBank for our EcM genera globally...")

if (file.exists(global_gb_ckpt)) {
  ts("  Loading checkpointed GenBank global metadata...")
  gb_global_meta <- readr::read_csv(global_gb_ckpt, show_col_types = FALSE)
  ts(sprintf("  Records loaded: %d", nrow(gb_global_meta)))
} else {

  if (Sys.getenv("ENTREZ_KEY") == "") {
    ts("  ENTREZ_KEY not set in .Renviron — GenBank global step will be rate-limited.")
    ts("  Set ENTREZ_KEY for faster fetching.")
  }

  our_genera_list <- unique(trimws(emf$genus))
  ts(sprintf("  Querying GenBank globally for %d EcM genera...",
             length(our_genera_list)))

  # Helper: query one genus globally, return host_taxon metadata
  query_genus_globally <- function(genus_name) {
    search_term <- paste0(
      '"', genus_name, '"[Organism]',
      ' AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])',
      ' AND "Fungi"[Organism]'
    )
    tryCatch({
      search <- rentrez::entrez_search(
        db          = "nuccore",
        term        = search_term,
        retmax      = 0L,
        use_history = TRUE
      )
      if (search$count == 0L) return(NULL)

      # Cap at 2000 records per genus to keep runtime manageable
      n_fetch    <- min(search$count, 2000L)
      batch_size <- 200L
      n_batches  <- ceiling(n_fetch / batch_size)
      all_ids    <- character(0L)

      for (i in seq_len(n_batches)) {
        batch <- tryCatch(
          rentrez::entrez_search(
            db          = "nuccore",
            term        = search_term,
            retmax      = batch_size,
            retstart    = (i - 1L) * batch_size,
            use_history = FALSE
          ),
          error = function(e) NULL
        )
        if (!is.null(batch)) all_ids <- c(all_ids, batch$ids)
        Sys.sleep(0.15)
      }

      if (length(all_ids) == 0L) return(NULL)

      # Fetch metadata (esummary) for these IDs
      meta_batches <- split(all_ids,
                            ceiling(seq_along(all_ids) / 200L))
      meta_list    <- vector("list", length(meta_batches))

      for (j in seq_along(meta_batches)) {
        summ <- tryCatch(
          rentrez::entrez_summary(db = "nuccore", id = meta_batches[[j]]),
          error = function(e) NULL
        )
        if (!is.null(summ)) {
          if (inherits(summ, "esummary")) summ <- list(summ)
          meta_list[[j]] <- dplyr::bind_rows(lapply(summ, function(s) {
            subtype <- strsplit(if (!is.null(s$subtype)) s$subtype else "", "\\|")[[1L]]
            subname <- strsplit(if (!is.null(s$subname)) s$subname else "", "\\|")[[1L]]
            get_sub <- function(key) {
              idx <- which(subtype == key)
              if (length(idx) > 0L) subname[idx[1L]] else NA_character_
            }
            tibble::tibble(
              genus_queried   = genus_name,
              accession       = s$caption,
              organism        = s$organism,
              country_gb      = get_sub("country"),
              host_gb         = get_sub("host"),
              isolation_src   = get_sub("isolation_source")
            )
          }))
        }
        Sys.sleep(0.15)
      }
      dplyr::bind_rows(meta_list)
    }, error = function(e) {
      message(sprintf("  Genus '%s' query failed: %s", genus_name, conditionMessage(e)))
      NULL
    })
  }

  ts("  (This step may take 30-120 min depending on genus list size.)")
  ts("  Results are checkpointed; delete checkpoint to re-run.")

  genus_results <- vector("list", length(our_genera_list))
  for (i in seq_along(our_genera_list)) {
    ts(sprintf("  Genus %d/%d: %s",
               i, length(our_genera_list), our_genera_list[i]))
    genus_results[[i]] <- query_genus_globally(our_genera_list[i])
    Sys.sleep(0.5)
  }

  gb_global_meta <- dplyr::bind_rows(Filter(Negate(is.null), genus_results))
  ts(sprintf("  GenBank global records retrieved: %d", nrow(gb_global_meta)))

  readr::write_csv(gb_global_meta, global_gb_ckpt)
  ts(sprintf("  Saved GenBank global checkpoint -> %s", basename(global_gb_ckpt)))
}

# Extract host_taxon from GenBank global records
parse_host_from_text <- function(x) {
  pattern   <- "(?:of|on)\\s+([A-Z][a-z]+(?:\\s+[a-z]+)?)"
  has_match <- grepl(pattern, x, perl = TRUE)
  result    <- rep(NA_character_, length(x))
  result[has_match] <- sub(paste0(".*", pattern, ".*"), "\\1",
                            x[has_match], perl = TRUE)
  result
}

if (nrow(gb_global_meta) > 0) {
  gb_global_meta <- gb_global_meta |>
    dplyr::mutate(
      host_taxon_gb = dplyr::coalesce(
        dplyr::if_else(!is.na(host_gb) & nzchar(trimws(host_gb)),
                       trimws(host_gb), NA_character_),
        parse_host_from_text(isolation_src)
      ),
      host_clean = clean_host_name(host_taxon_gb)
    )

  n_gb_with_host <- sum(!is.na(gb_global_meta$host_clean))
  ts(sprintf("  GenBank global records with host info: %d / %d",
             n_gb_with_host, nrow(gb_global_meta)))

  # Q3 (GenBank): per fungal genus, which Canadian host species globally?
  gb_q3 <- gb_global_meta |>
    dplyr::filter(!is.na(host_clean),
                  host_clean %in% host_tbl$species) |>
    dplyr::distinct(genus_queried, host_clean) |>
    dplyr::rename(fungal_genus = genus_queried, host_species = host_clean) |>
    dplyr::group_by(fungal_genus) |>
    dplyr::summarise(
      n_host_species_gb = dplyr::n_distinct(host_species),
      host_species_gb   = paste(sort(unique(host_species)), collapse = "; "),
      .groups           = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_host_species_gb))

  readr::write_csv(gb_q3,
                   file.path(paths$out_eltonian, "eltonian_global_gb_fungi_to_hosts.csv"))
  ts(sprintf("  Q3 (GenBank): EcM genera with host associations: %d", nrow(gb_q3)))
  ts("  Saved eltonian_global_gb_fungi_to_hosts.csv")

  # Q1 (GenBank): which Canadian host species appear as hosts globally in GB?
  gb_q1_hosts <- unique(gb_global_meta$host_clean[
    !is.na(gb_global_meta$host_clean) &
      gb_global_meta$host_clean %in% host_tbl$species
  ])
  ts(sprintf("  Q1 (GenBank): Canadian hosts with global GenBank records: %d / %d",
             length(gb_q1_hosts), n_host_species))

  # Long-form GenBank global pairs for the interactions Rds (Part C)
  # Host → genus (global, GenBank, from genus_queried field)
  host_genus_global_gb <- gb_global_meta |>
    dplyr::filter(!is.na(host_clean), host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, fungal_genus = genus_queried)
  ts(sprintf("  host_genus_global_gb:   %d host–genus pairs (GenBank)",
             nrow(host_genus_global_gb)))

  # Host → named EcM species (global, GenBank, from organism field)
  host_species_global_gb <- gb_global_meta |>
    dplyr::filter(!is.na(host_clean), host_clean %in% host_tbl$species,
                  !is.na(organism), nzchar(trimws(organism))) |>
    dplyr::mutate(fungal_species = trimws(organism)) |>
    dplyr::filter(!grepl("^(Fungi|Basidiomycota|Ascomycota|uncultured|sp\\.|NA)$",
                         fungal_species, ignore.case = TRUE),
                  nchar(fungal_species) > 3L) |>
    dplyr::distinct(host_species = host_clean, fungal_species)
  ts(sprintf("  host_species_global_gb: %d host–species pairs (GenBank)",
             nrow(host_species_global_gb)))

} else {
  ts("  No GenBank global records retrieved.")
  gb_q1_hosts            <- character(0L)
  gb_q3                  <- NULL
  host_genus_global_gb   <- NULL
  host_species_global_gb <- NULL
}

# =============================================================================
# PART C: Assemble and save interaction list
# =============================================================================
# Stores 12 long-form pair data frames as a named list in one Rds file.
# Each element has two columns: focal taxon and partner taxon.
# "host_" prefix = host species as focal; reversed = fungal taxon as focal.
# Canada scope: interactions observed in our Canadian EMF dataset.
# Global scope: interactions from GF root samples + GenBank globally.
# =============================================================================

ts("Part C: Assembling eltonian_interactions.rds...")

# ---------- Canada scope, host as focal --------------------------------------
host_sh_canada      <- sh_pairs             # cols: host_species, sh_code
host_genus_canada   <- genus_pairs          # cols: host_species, fungal_genus
host_species_canada <- species_pairs_canada # cols: host_species, fungal_species

# ---------- Canada scope, fungal taxon as focal (same rows, roles reversed) --
sh_host_canada      <- dplyr::select(host_sh_canada,      sh_code,       host_species)
genus_host_canada   <- dplyr::select(host_genus_canada,   fungal_genus,  host_species)
species_host_canada <- dplyr::select(host_species_canada, fungal_species, host_species)

# ---------- Global scope, host as focal — combine GF + GenBank ---------------
host_sh_global <- if (!is.null(host_sh_global_gf) && nrow(host_sh_global_gf) > 0) {
  host_sh_global_gf
} else {
  data.frame(host_species = character(0L), sh_code = character(0L))
}

host_species_global <- dplyr::bind_rows(
  host_species_global_gf,
  host_species_global_gb
) |> dplyr::distinct()

host_genus_global <- dplyr::bind_rows(
  host_genus_global_gf,
  host_genus_global_gb
) |> dplyr::distinct()

# ---------- Global scope, fungal taxon as focal (same rows, roles reversed) --
sh_host_global <- if (nrow(host_sh_global) > 0) {
  dplyr::select(host_sh_global, sh_code, host_species)
} else {
  data.frame(sh_code = character(0L), host_species = character(0L))
}

species_host_global <- if (nrow(host_species_global) > 0) {
  dplyr::select(host_species_global, fungal_species, host_species)
} else {
  data.frame(fungal_species = character(0L), host_species = character(0L))
}

genus_host_global <- if (nrow(host_genus_global) > 0) {
  dplyr::select(host_genus_global, fungal_genus, host_species)
} else {
  data.frame(fungal_genus = character(0L), host_species = character(0L))
}

# ---------- Assemble and save ------------------------------------------------
eltonian_interactions <- list(
  # Canada scope — host species as focal
  host_sh_canada      = host_sh_canada,
  host_species_canada = host_species_canada,
  host_genus_canada   = host_genus_canada,
  # Canada scope — fungal taxon as focal
  sh_host_canada      = sh_host_canada,
  species_host_canada = species_host_canada,
  genus_host_canada   = genus_host_canada,
  # Global scope — host species as focal (GF root + GenBank)
  host_sh_global      = host_sh_global,
  host_species_global = host_species_global,
  host_genus_global   = host_genus_global,
  # Global scope — fungal taxon as focal
  sh_host_global      = sh_host_global,
  species_host_global = species_host_global,
  genus_host_global   = genus_host_global
)

saveRDS(eltonian_interactions,
        file.path(paths$out_eltonian, "eltonian_interactions.rds"))
ts("  Saved eltonian_interactions.rds")
for (nm in names(eltonian_interactions)) {
  ts(sprintf("    $%-25s  %d rows", nm, nrow(eltonian_interactions[[nm]])))
}

# =============================================================================
# Summary tables
# =============================================================================

ts("Saving summary tables...")

eltonian_summary <- tibble::tibble(
  metric = c(
    # Canadian-scope denominators and basic interaction counts
    "Canadian EcM host plant species (denominator: BIEN-based native list)",
    "Canadian EcM host genera (denominator: BIEN-based native list)",
    "Canadian host species with >= 1 observed EcM-fungus association in Canada (paired with fungal genus)",
    "Canadian host species with >= 1 observed EcM-fungus association in Canada (paired with fungal SH code)",
    "% of Canadian EcM host species with any documented fungal association in Canada (paired with fungal genus)",
    "% of Canadian EcM host species with any documented fungal association in Canada (paired with fungal SH code)",
    "Canadian EcM fungal genera with >= 1 documented host species in Canada",
    "Canadian EcM SH codes with >= 1 documented host species in Canada",
    "Potential host x EcM-genus interaction pairs (n_hosts × n_EcM_genera; full Canadian scope)",
    "Observed host x EcM-genus interaction pairs (Canadian scope)",
    "% of potential host x EcM-genus interaction pairs observed (Canadian scope)",
    "Mean EcM genera documented per Canadian host species (averaged over hosts with >= 1 association)",
    "Max EcM genera documented for any Canadian host species",
    "Mean Canadian host species per EcM genus (averaged over genera with >= 1 host)",
    "Max Canadian host species for any EcM genus",
    # Species-level interaction matrix: fill rate and occurrence support
    "Named EcM fungal species detected anywhere in the Canadian dataset (full column denominator)",
    "Full host x named-species matrix size (n_host_species x n_named_fungal_species)",
    "Filled cells in full host x named-species matrix (observed host-species x fungal-species pairs)",
    "Empty cells in full host x named-species matrix (no observed pair)",
    "% of full host x named-species matrix cells filled",
    "% of full host x named-species matrix cells empty",
    "Observed host x named-species pairs supported by exactly 1 occurrence (sample/record)",
    "% of observed host x named-species pairs supported by exactly 1 occurrence",
    # Genus-level interaction matrix: fill rate and occurrence support
    "EcM fungal genera detected anywhere in the Canadian dataset (full genus-matrix column denominator)",
    "Full host-genus x fungal-genus matrix size (n_host_genera x n_fungal_genera_total)",
    "Filled cells in full host-genus x fungal-genus matrix (observed host-genus x fungal-genus pairs)",
    "Empty cells in full host-genus x fungal-genus matrix (no observed pair)",
    "% of full host-genus x fungal-genus matrix cells filled",
    "% of full host-genus x fungal-genus matrix cells empty",
    "Observed host-genus x fungal-genus pairs supported by exactly 1 occurrence (sample/record)",
    "% of observed host-genus x fungal-genus pairs supported by exactly 1 occurrence",
    # Global-scope: how well are Canadian hosts and Canadian EcM fungi documented worldwide?
    "Canadian host species recorded as hosts in GlobalFungi root samples anywhere in the world",
    "% of Canadian host species with any global GlobalFungi root record",
    "Canadian host species recorded as hosts in GenBank EcM records anywhere in the world",
    "Canadian host species with documented global EcM-species associations (GlobalFungi root samples)",
    "Canadian EcM species with documented global host associations (GlobalFungi root samples)"
  ),
  value = c(
    n_host_species,
    n_host_genera,
    n_hosts_with_genus_data <- dplyr::n_distinct(genus_pairs$host_species),
    n_hosts_with_sh_data    <- dplyr::n_distinct(sh_pairs$host_species),
    round(100 * n_hosts_with_genus_data / n_host_species, 1),
    round(100 * n_hosts_with_sh_data    / n_host_species, 1),
    n_genera_with_host,
    n_sh_with_host,
    n_potential_genus,
    n_obs_genus_pairs,
    round(100 * n_obs_genus_pairs / n_potential_genus, 2),
    round(mean(genus_per_host$n_genera), 1),
    max(genus_per_host$n_genera),
    round(mean(host_per_genus$n_hosts), 1),
    max(host_per_genus$n_hosts),
    n_named_fungal_species,
    n_matrix_cells,
    n_cells_filled,
    n_cells_empty,
    pct_cells_filled,
    pct_cells_empty,
    n_pairs_singleton,
    pct_pairs_singleton,
    n_fungal_genera_total,
    n_matrix_cells_genus,
    n_cells_filled_genus,
    n_cells_empty_genus,
    pct_cells_filled_genus,
    pct_cells_empty_genus,
    n_pairs_singleton_genus,
    pct_pairs_singleton_genus,
    ifelse(!is.na(n_canada_hosts_with_gf), n_canada_hosts_with_gf, NA_real_),
    ifelse(!is.na(pct_with_gf), pct_with_gf, NA_real_),
    length(gb_q1_hosts),
    ifelse(!is.null(q2), nrow(q2), NA_real_),
    ifelse(!is.null(q3), nrow(q3), NA_real_)
  )
)

readr::write_csv(eltonian_summary,
                 file.path(paths$out_eltonian, "eltonian_summary.csv"))
ts("Saved eltonian_summary.csv")
print(as.data.frame(eltonian_summary))
ts("19_eltonian.R complete.")
