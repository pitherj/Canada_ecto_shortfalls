# =============================================================================
# Eltonian Shortfall: Host–Fungus Interaction Knowledge
# =============================================================================
#
# PART A — Canadian scope
#   What proportion of potential EcM host–fungus interactions in Canada have
#   been observed in our dataset?
#   Host species: BIEN-based native Canadian EcM hosts (06_host_species.R)
#   Fungal taxa:  EcM sequences from GenBank and GlobalFungi in our dataset
#
# PART B — Global scope
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
#   Part A also maps this knowledge onto geographic space (steps A5-A9): which
#   parts of Canada hold EcM host habitat, and what proportion of the host
#   species present in each 0.5-degree grid cell have EcM sequence data? The
#   host species counted as "having data" are exactly those flagged `matched`
#   in the Part A host-matching table, so the maps and the interaction
#   analyses share one definition of a documented host.
#
# Prerequisite files:
#   data_derived/ecm_native_canada_host_species.csv  (06_host_species.R)
#   data_derived/spatial/canada_simple.gpkg          (01_spatial_data.R)
#   data_derived/spatial/bien_host_richness_0.5deg.tif  (08_host_rasters.R)
#   data_derived/spatial/bien_host_species_stack.tif    (08_host_rasters.R)
#   data_raw/GlobalFungi/GlobalFungi_5_sample_metadata.txt    (full global)
#   data_raw/GlobalFungi/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt (full, 13 GB)
#   data_derived/unite_sh_taxonomy.csv               (02_globalfungi.R)
#
# Outputs:
#   data_derived/eltonian/eltonian_host_list.csv
#   data_derived/eltonian/eltonian_host_matching.csv
#   data_derived/eltonian/eltonian_matrix_genus.csv
#   data_derived/eltonian/eltonian_matrix_sh.csv
#   data_derived/eltonian/eltonian_matrix_species.csv — host x named-EcM-species
#       binary matrix, trimmed to taxa with >= 1 observed pair (same
#       convention as the genus/SH matrices above)
#   data_derived/eltonian/eltonian_species_occurrence_counts.csv — per host-species x
#       named-fungal-species pair, count of supporting samples/records
#       (occurrences); backs the singleton-association statistic in
#       eltonian_summary.csv
#   data_derived/eltonian/eltonian_matrix_genus_genus.csv — host-genus x fungal-genus
#       binary matrix (both axes collapsed to genus from the species-exact
#       `matched` interactions), trimmed to taxa with >= 1 observed pair,
#       same convention as the SH/species/genus matrices above
#   data_derived/eltonian/eltonian_genus_occurrence_counts.csv — per host-genus x
#       fungal-genus pair, count of supporting samples/records (occurrences);
#       backs the genus-level singleton-association statistic in
#       eltonian_summary.csv
#   data_derived/eltonian/eltonian_interactions.rds  — assembled interaction
#       objects retained for downstream reuse
#   data_derived/eltonian/eltonian_summary.csv
#   data_derived/eltonian/eltonian_global_host_coverage.csv    — Q1: Canadian hosts with global data
#   data_derived/eltonian/eltonian_global_host_to_fungi.csv    — Q2: per-host EcM fungal associations
#   data_derived/eltonian/eltonian_global_fungi_to_hosts.csv   — Q3: per-fungal-species host associations
#   data_derived/eltonian/eltonian_global_gb_fungi_to_hosts.csv — Q3, GenBank scope
#   data_derived/eltonian/eltonian_sample_type_tally_canada.csv — diagnostic: sample_type
#       composition of EcM-positive Canadian GF samples with a
#       dominant_plant_species entry
#   data_derived/eltonian/eltonian_sample_type_tally_global.csv — diagnostic: sample_type
#       composition of EcM-positive global GF samples (>=1 of our Canadian
#       EcM SH codes detected) with a dominant_plant_species entry
#   data_derived/eltonian/eltonian_genbank_tissue_tally_canada.csv — diagnostic: tissue-type
#       composition of Canadian GenBank EcM records with host information
#       (host_taxon), keyword-binned from the free-text isolation_src field
#   data_derived/eltonian/eltonian_host_raster_summary.csv — 0.5-degree grid-cell
#       summary of host-habitat coverage (Part A, step A8)
#   data_derived/spatial/bien_host_data_richness_0.5deg.tif   — host spp. with EcM data
#   data_derived/spatial/bien_host_data_proportion_0.5deg.tif — proportion with EcM data
#   figures/Figure-04_host_bivariate_map.png      (paths$fig_host_bivariate)      -- white bg, used in manuscript
#   figures/Figure-04_host_bivariate_map_grey.png (paths$fig_host_bivariate_grey) -- #F2F2F2 bg, Figure 5 panel source
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(rentrez)
library(data.table)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

# The GADM-derived Canada polygon has minor topology issues that the S2
# spherical geometry engine rejects; switch sf to planar geometry for the
# spatial steps below.
sf::sf_use_s2(FALSE)

# Host-name cleaning is handled by the shared canonicalize_host() helper in
# 00_setup.R. GlobalFungi `dominant_plant_species` and `other_plant_species`
# values use underscores as token separators in some records, so we convert
# underscores to spaces before passing into canonicalize_host().
clean_host_name <- function(x) canonicalize_host(gsub("_", " ", x))

# =============================================================================
# PART A: Canadian scope
# =============================================================================

# ---- A1. Load Canadian EcM host species list ---------------------------------

if (!file.exists(paths$host_species)) {
  stop("Host species file not found: ", paths$host_species,
       "\nRun 06_host_species.R first.")
}
host_tbl <- readr::read_csv(paths$host_species, show_col_types = FALSE)
n_host_species <- nrow(host_tbl)
n_host_genera  <- dplyr::n_distinct(sub("^(\\S+).*", "\\1", host_tbl$species))
readr::write_csv(host_tbl, file.path(paths$out_eltonian, "eltonian_host_list.csv"))

# ---- A2. Extract and clean host strings from EMF dataset (Canadian records) --
# IMPORTANT: dominant_plant_species and other_plant_species are only reliable
# host indicators for root samples (sample_type == "root"). For soil and other
# sample types these fields record the dominant plant at the site but do NOT
# imply a direct mycorrhizal association. The host_taxon field (from GenBank
# host or isolation_source parsing) is used unconditionally.

# Root samples only for GlobalFungi dominant/other plant species fields
emf_gf_root <- dplyr::filter(emf,
                               source == "GlobalFungi",
                               sample_type == "root")

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

host_long <- host_long |>
  dplyr::mutate(
    host_clean  = clean_host_name(host_raw),
    host_genus  = sub("^(\\S+).*", "\\1", host_clean),
    matched     = host_clean %in% host_tbl$species,
    match_genus = host_genus %in% sub("^(\\S+).*", "\\1", host_tbl$species)
  )

n_matched       <- sum(host_long$matched,     na.rm = TRUE)
n_genus_matched <- sum(host_long$match_genus, na.rm = TRUE)

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

sample_type_tally_canada <- emf |>
  dplyr::filter(source == "GlobalFungi", !is.na(dominant_plant_species)) |>
  dplyr::distinct(sample_ID, sample_type) |>
  dplyr::count(sample_type, name = "n_samples") |>
  dplyr::arrange(dplyr::desc(n_samples))

readr::write_csv(sample_type_tally_canada,
                 file.path(paths$out_eltonian, "eltonian_sample_type_tally_canada.csv"))

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

readr::write_csv(genbank_tissue_tally_canada,
                 file.path(paths$out_eltonian, "eltonian_genbank_tissue_tally_canada.csv"))

# ---- A3. Build interaction matrices ------------------------------------------

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
  dplyr::filter(!is.na(sh_code)) |>   # genus-resolved GenBank rows carry sh_code = NA
  dplyr::distinct(host_clean, sh_code) |>
  dplyr::rename(host_species = host_clean)

sh_matrix <- sh_pairs |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_species, names_from = sh_code,
                     values_from = present, values_fill = 0L)
readr::write_csv(sh_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_sh.csv"))

# ---- A3b. Species-level Canada pairs ----------------------------------------
# Join sh_pairs with UNITE taxonomy to resolve SH codes → named species
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

# ---- A3c. Species-level matrix, fill statistics, and occurrence counts -------
# Trimmed host x named-EcM-species matrix, following the same convention as
# eltonian_matrix_genus.csv / eltonian_matrix_sh.csv above: rows and columns
# only for taxa with >= 1 observed pair, NOT the full host x species
# denominator. We report the full 147 x 1079 grid as a numeric fill statistic
# below rather than materializing it as a mostly-empty CSV.

species_matrix <- species_pairs_canada |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_species, names_from = fungal_species,
                     values_from = present, values_fill = 0L)
readr::write_csv(species_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_species.csv"))

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

n_pairs_total       <- nrow(species_occurrence_counts)
n_pairs_singleton   <- sum(species_occurrence_counts$n_occurrences == 1L)
pct_pairs_singleton <- round(100 * n_pairs_singleton / n_pairs_total, 1)

# ---- A3d. Genus-level matrix, fill statistics, and occurrence counts --------
# Genus x genus analogue of the host x named-species block above (A3c).
# Both axes are collapsed to genus: the host axis from host_clean (species)
# to host_genus, the fungal axis is already genus-resolution (emf$genus).
# Uses matched_interactions (i.e. the SAME species-exact host match, `matched`,
# used to build genus_pairs/sh_pairs/species_pairs_canada above) collapsed to
# host_genus -- NOT the looser `match_genus` flag in host_long -- so the
# matching criterion stays consistent across the SH-, species-, and
# genus-resolution blocks. The looser match_genus criterion would give
# slightly different counts -- 12 host genera x 69 fungal genera, 263 filled
# cells -- because it admits host records whose species name did not resolve to
# the BIEN list but whose genus did; we use the stricter `matched` criterion
# here instead.
#
# As with the species-level matrix, the saved CSV is trimmed to taxa with
# >= 1 observed pair; the full-grid fill rate is reported as a numeric
# statistic in eltonian_summary.csv rather than materialized as a mostly-
# empty CSV.

genus_genus_pairs_canada <- matched_interactions |>
  dplyr::distinct(host_genus, genus) |>
  dplyr::rename(fungal_genus = genus)

genus_genus_matrix <- genus_genus_pairs_canada |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(id_cols = host_genus, names_from = fungal_genus,
                     values_from = present, values_fill = 0L)
readr::write_csv(genus_genus_matrix,
                 file.path(paths$out_eltonian, "eltonian_matrix_genus_genus.csv"))

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

n_pairs_total_genus       <- nrow(genus_occurrence_counts)
n_pairs_singleton_genus   <- sum(genus_occurrence_counts$n_occurrences == 1L)
pct_pairs_singleton_genus <- round(100 * n_pairs_singleton_genus / n_pairs_total_genus, 1)

# ---- A4. Canadian-scope coverage statistics ----------------------------------

n_obs_genus_pairs       <- nrow(genus_pairs)
n_hosts_with_genus_data <- dplyr::n_distinct(genus_pairs$host_species)
n_hosts_with_sh_data    <- dplyr::n_distinct(sh_pairs$host_species)
n_genera_with_host      <- dplyr::n_distinct(genus_pairs$fungal_genus)
n_sh_with_host          <- dplyr::n_distinct(sh_pairs$sh_code, na.rm = TRUE)
n_potential_genus       <- n_host_species * dplyr::n_distinct(emf$genus)
genus_per_host          <- dplyr::count(genus_pairs, host_species, name = "n_genera")
host_per_genus          <- dplyr::count(genus_pairs, fungal_genus, name = "n_hosts")

# ---- A5. Spatial prerequisites for the host-coverage rasters -----------------
# The 0.5-degree host-richness raster and the matching per-species binary stack
# are produced by 08_host_rasters.R. The Canada boundary is used both to
# convert host layers to a Canada-only extent and to draw the maps in Step A9.

canada_bound <- sf::st_read(paths$canada_bound, quiet = TRUE)

if (!file.exists(paths$bien_richness)) {
  stop("Host richness raster not found: ", paths$bien_richness,
       "\nRun 08_host_rasters.R first.")
}
richness_wgs84 <- terra::rast(paths$bien_richness)
names(richness_wgs84) <- "richness"

# ---- A6. Identify host species with EcM data in our dataset ------------------
# `host_long` (Step A2) is the host-matching table: one row per host string
# extracted from the EcM dataset, with `host_clean` canonicalized and `matched`
# flagging strings that resolve to a species on the BIEN-based Canadian EcM
# host list. The species that matched are exactly the host species for which we
# hold at least one EcM sequence record.

em_species_with_data <- unique(
  host_long$host_clean[host_long$matched & !is.na(host_long$host_clean)]
)

# Tree / non-tree host species subsets (used for bivariate maps)
em_canada_tree_species    <- host_tbl$species[host_tbl$growth_form %in% "tree"]
em_canada_nontree_species <- host_tbl$species[!host_tbl$growth_form %in% "tree"]
em_species_with_data_tree    <- em_species_with_data[em_species_with_data %in% em_canada_tree_species]
em_species_with_data_nontree <- em_species_with_data[em_species_with_data %in% em_canada_nontree_species]

# ---- A7. Load species stack and derive all richness rasters ------------------
# The per-species binary stack (one layer per host species, named by species)
# is produced by 08_host_rasters.R alongside bien_host_richness_0.5deg.tif.
# All subset richness rasters are derived by filtering layers + summing, with
# no range re-rasterization needed here.

if (!file.exists(paths$bien_species_stack)) {
  stop("Species stack not found: ", paths$bien_species_stack,
       "\nRun 08_host_rasters.R first.", call. = FALSE)
}

species_stack <- terra::rast(paths$bien_species_stack)
canada_vect   <- terra::vect(sf::st_transform(canada_bound, 4326))

# Normalise layer names: BIEN stores species with underscores ("Abies_amabilis")
# but all other species lists in this project use spaces ("Abies amabilis").
names(species_stack) <- gsub("_", " ", names(species_stack))

# Helper: sum selected layers; uses na.rm = FALSE so outside-Canada cells
# (all NA in the masked stack) remain NA rather than becoming 0.
# zero_to_na = TRUE matches the richness_wgs84 convention (no-habitat cells
# are NA); FALSE preserves 0 for data richness (0 = habitat exists, no data).
sum_species_layers <- function(stack, species_names, zero_to_na = FALSE) {
  idx <- which(names(stack) %in% species_names)
  if (length(idx) == 0L) {
    warning("No matching layers found in species stack for provided species names.")
    return(NULL)
  }
  r <- terra::app(stack[[idx]], fun = "sum", na.rm = FALSE)
  if (zero_to_na) r <- terra::ifel(r == 0, NA, r)
  r
}

# All-species data richness (hosts present in our EcM dataset)
data_richness_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data, zero_to_na = FALSE
)
names(data_richness_wgs84) <- "data_richness"
terra::writeRaster(data_richness_wgs84, paths$bien_data_rich, overwrite = TRUE)

# Tree-only richness (all tree host species)
richness_tree_wgs84 <- sum_species_layers(
  species_stack, em_canada_tree_species, zero_to_na = TRUE
)
names(richness_tree_wgs84) <- "richness"

# Tree-only data richness (tree host species in our EcM dataset)
data_richness_tree_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data_tree, zero_to_na = FALSE
)
names(data_richness_tree_wgs84) <- "data_richness"

# Non-tree richness and data richness
richness_nontree_wgs84 <- sum_species_layers(
  species_stack, em_canada_nontree_species, zero_to_na = TRUE
)
names(richness_nontree_wgs84) <- "richness"

data_richness_nontree_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data_nontree, zero_to_na = FALSE
)
names(data_richness_nontree_wgs84) <- "data_richness"

# ---- A8. Proportion raster and raster summary --------------------------------

proportion_wgs84 <- terra::clamp(
  data_richness_wgs84 / richness_wgs84,
  lower = 0, upper = 1
)
# Cells where richness is NA (outside Canada) remain NA
names(proportion_wgs84) <- "proportion"

terra::writeRaster(proportion_wgs84, paths$bien_proportion, overwrite = TRUE)

# Summary statistics
prop_vals <- terra::values(proportion_wgs84)
prop_vals <- prop_vals[!is.na(prop_vals)]

raster_summary <- tibble::tibble(
  metric = c(
    "Grid cells (0.5°) with EcM host habitat",
    "Cells with EcM host habitat but no sequence records",
    "Cells with EcM host habitat and >=1 sequence record",
    "Mean proportion of host spp. with sequence data (per cell)",
    "Median proportion of host spp. with sequence data (per cell)",
    "Max proportion of host spp. with sequence data (per cell)"
  ),
  value = c(
    length(prop_vals),
    sum(prop_vals == 0),
    sum(prop_vals > 0),
    round(mean(prop_vals), 3),
    round(median(prop_vals), 3),
    round(max(prop_vals), 3)
  )
)
readr::write_csv(raster_summary,
                 file.path(paths$out_eltonian, "eltonian_host_raster_summary.csv"))

# ---- A9. Bivariate map (host richness x proportion with data) ----------------
# Two panels: (1) all EcM host species, (2) tree host species only.
# Breaks are computed independently for each panel so each uses its own tertiles.

canada_albers <- sf::st_transform(canada_bound, crs_albers)

bivar_colors <- c(
  "1-1" = "#e8e8e8", "2-1" = "#ace4e4", "3-1" = "#5ac8c8",
  "1-2" = "#dfb0d6", "2-2" = "#a5add3", "3-2" = "#5698b9",
  "1-3" = "#d272aa", "2-3" = "#ad6aad", "3-3" = "#7759a1"
)

# Helper: build one bivariate map panel + inset legend
# richness_r, data_r are WGS84 SpatRasters; title_label is a string; bg sets
# the panel/plot background colour (white for the manuscript version, grey
# for the Figure 5 schematic source version -- see build_bivariate_figure()).
make_bivar_panel <- function(richness_r, data_r, title_label, bg = "white") {

  prop_r <- terra::clamp(data_r / richness_r, lower = 0, upper = 1)
  names(prop_r) <- "proportion"

  rich_alb <- terra::project(richness_r, crs_albers, method = "near")
  data_alb <- terra::project(data_r,     crs_albers, method = "near")
  prop_alb <- terra::project(prop_r,     crs_albers, method = "bilinear")
  names(rich_alb) <- "richness"
  names(data_alb) <- "n_with_data"
  names(prop_alb) <- "proportion"

  # Clip to the Canada boundary (study-area restriction, not a projection fix).
  # The 0.5-deg BIEN host-richness grid carries ~565 richness > 0 cells outside
  # Canada (host ranges crossing into the contiguous US, plus half-cell straddle
  # of the coarse cells along the border). This is a Canada-scoped figure, so
  # those out-of-country cells are masked out.
  canada_v_b <- terra::vect(canada_albers)
  rich_alb   <- terra::mask(rich_alb, canada_v_b)
  data_alb   <- terra::mask(data_alb, canada_v_b)
  prop_alb   <- terra::mask(prop_alb, canada_v_b)

  # Breaks computed from this panel's values only
  rv <- terra::values(rich_alb)
  rv <- rv[!is.na(rv) & rv > 0]
  rich_breaks <- c(-Inf, stats::quantile(rv, probs = c(1/3, 2/3), names = FALSE), Inf)
  prop_breaks <- c(-Inf, 1/3, 2/3, Inf)

  bdf <- as.data.frame(c(rich_alb, data_alb, prop_alb), xy = TRUE) |>
    stats::setNames(c("x", "y", "richness", "n_with_data", "proportion")) |>
    dplyr::filter(!is.na(richness), richness > 0) |>
    dplyr::mutate(
      proportion = dplyr::if_else(is.na(proportion), 0, proportion),
      rich_class = as.integer(cut(richness,   breaks = rich_breaks, labels = 1:3)),
      prop_class = as.integer(cut(proportion, breaks = prop_breaks, labels = 1:3)),
      bi_class   = paste0(rich_class, "-", prop_class)
    )

  rich_q      <- round(stats::quantile(rv, probs = c(0, 1/3, 2/3, 1), names = FALSE))
  rich_labels <- paste0(rich_q[1:3], "\u2013", rich_q[2:4])
  prop_labels <- c("0\u201333%", "33\u201367%", "67\u2013100%")

  legend_df <- expand.grid(x = factor(1:3), y = factor(1:3)) |>
    dplyr::mutate(bi_class = paste0(as.integer(x), "-", as.integer(y)))

  p_leg <- ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = y, fill = bi_class)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::scale_fill_manual(values = bivar_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = rich_labels) +
    ggplot2::scale_y_discrete(labels = prop_labels) +
    ggplot2::labs(x = "Host richness", y = "Prop. with data") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text  = ggplot2::element_text(colour = "black"),
                   axis.title = ggplot2::element_text(colour = "black", size = 7))

  p_map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = canada_albers, fill = "grey95", colour = NA) +
    ggplot2::geom_raster(data = bdf, ggplot2::aes(x = x, y = y, fill = bi_class)) +
    ggplot2::geom_sf(data = canada_albers, fill = NA, colour = "grey40",
                     linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = bivar_colors, guide = "none") +
    ggplot2::coord_sf(crs = crs_albers) +
    ggplot2::labs(title = title_label) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.background = ggplot2::element_rect(fill = bg, colour = NA)
    ) +
    patchwork::inset_element(p_leg,
                             left = 0.72, bottom = 0.62, right = 0.99, top = 0.87)
  p_map
}

# build_bivariate_figure(bg) assembles the two-panel (non-tree / tree)
# composite at the given background colour.
build_bivariate_figure <- function(bg) {
  p_all  <- make_bivar_panel(richness_nontree_wgs84, data_richness_nontree_wgs84,
                             "Non-tree EcM host species", bg = bg)
  p_tree <- make_bivar_panel(richness_tree_wgs84, data_richness_tree_wgs84,
                             "Tree EcM host species", bg = bg)

  (p_all / p_tree) +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.background = ggplot2::element_rect(fill = bg, colour = NA))
    )
}

save_fig_formats(paths$fig_host_bivariate, build_bivariate_figure("white"),
                 width = 10, height = 14, dpi = 300, bg = "white")

# Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
# schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
ggplot2::ggsave(paths$fig_host_bivariate_grey, build_bivariate_figure(fig5_grey_bg),
                width = 10, height = 14, dpi = 300, bg = fig5_grey_bg)

# =============================================================================
# PART B: Global scope
# =============================================================================
# Paths to full GlobalFungi data (not just Canadian subset)
gf_meta_path <- file.path(here::here("data_raw"), "GlobalFungi",
                            "GlobalFungi_5_sample_metadata.txt")
gf_sh_path   <- file.path(here::here("data_raw"), "GlobalFungi",
                            "GlobalFungi_5_SH_abundance_ITS1_ITS2.txt")

# Checkpoints for slow operations
global_meta_ckpt <- file.path(paths$temp_dir, "gf_global_root_metadata.csv")
global_sh_ckpt   <- file.path(paths$temp_dir, "gf_global_ecm_sh_subset.rds")
global_gb_ckpt   <- file.path(paths$temp_dir, "genbank_global_ecm_meta.csv")

# ---- B1. GlobalFungi: global root-sample metadata ---------------------------
# Q1 (metadata-only approach): which Canadian host species appear as
# dominant_plant_species in root samples anywhere in the full GlobalFungi DB?

if (!file.exists(gf_meta_path)) {
  gf_global_root <- NULL
} else if (file.exists(global_meta_ckpt)) {
  gf_global_root <- readr::read_csv(global_meta_ckpt, show_col_types = FALSE)
} else {
  # Use data.table::fread for speed on large file
  gf_meta_all <- data.table::fread(
    gf_meta_path,
    sep            = "\t",
    quote          = "",
    select         = c("sample_ID", "country", "sample_type",
                       "dominant_plant_species", "other_plant_species"),
    showProgress   = TRUE
  )

  # Filter to root samples only
  gf_global_root <- dplyr::filter(gf_meta_all, sample_type == "root")
  rm(gf_meta_all)

  readr::write_csv(gf_global_root, global_meta_ckpt)
}

# ---- B2. Q1: Which Canadian host species have global GF root records? --------

if (!is.null(gf_global_root) && nrow(gf_global_root) > 0) {

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

  # Save Q1 result
  q1_host_coverage <- data.frame(
    species         = host_tbl$species,
    has_global_gf   = host_tbl$species %in% canada_hosts_with_global_gf
  )
  readr::write_csv(q1_host_coverage,
                   file.path(paths$out_eltonian, "eltonian_global_host_coverage.csv"))

} else {
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

if (!file.exists(gf_sh_path)) {
  gf_sh_ecm <- NULL

} else if (file.exists(global_sh_ckpt)) {
  gf_sh_ecm <- readRDS(global_sh_ckpt)

} else {
  # Identify our Canadian EcM SH codes (from emf dataset)
  our_sh_codes <- unique(emf$sh_code[!is.na(emf$sh_code)])

  # Always check the file header first to find the intersection of our SH codes
  # with the file's columns. Our codes use a specific UNITE version suffix
  # (e.g., .10FU) which may differ from the version used in GlobalFungi.
  header_cols <- names(data.table::fread(gf_sh_path, sep = "\t", quote = "",
                                          nrows = 0L))
  present_sh <- intersect(our_sh_codes, header_cols)

  if (length(present_sh) == 0L) {
    # Version mismatch: strip version suffix and try prefix matching
    # e.g., "SH1052460.10FU" -> "SH1052460"
    our_sh_prefix   <- sub("\\.[0-9]+FU$", "", our_sh_codes)
    header_prefix   <- sub("\\.[0-9]+FU$", "", header_cols)
    matched_idx     <- match(our_sh_prefix, header_prefix)
    matched_header  <- header_cols[matched_idx[!is.na(matched_idx)]]
    present_sh      <- matched_header
  }

  if (length(present_sh) == 0L) {
    gf_sh_ecm_wide <- NULL
  } else {
    # awk-streaming subset helper (00_setup.R) avoids fread()'s 2^31-byte string
    # limit on the ~13 GB matrix.
    gf_sh_ecm_wide <- read_big_tsv_subset(gf_sh_path, c("sample_ID", present_sh))
  }

  # Melt to long format, filter to non-zero abundances
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

  saveRDS(gf_sh_ecm, global_sh_ckpt)
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

if (file.exists(global_sampletype_tally_out)) {
  sample_type_tally_global <- readr::read_csv(global_sampletype_tally_out,
                                              show_col_types = FALSE)
} else if (is.null(gf_sh_ecm) || !file.exists(gf_meta_path)) {
  sample_type_tally_global <- NULL
} else {
  ecm_positive_samples <- unique(gf_sh_ecm$sample_ID)

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
}

# Q2 and Q3 from GlobalFungi
if (!is.null(gf_sh_ecm) && !is.null(gf_root_canadian_hosts)) {

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

  # ---- Long-form pair data frames for the interactions Rds (Part C) ---------

  # Host → SH code (global, GF root samples, Canadian hosts only)
  # many-to-many expected here too (see NOTE above Q3) — distinct() follows.
  host_sh_global_gf <- gf_sh_ecm |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::filter(host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, sh_code)

  # Host → named EcM species (global, GF)
  # many-to-many expected here too (see NOTE above Q3) — distinct() follows.
  host_species_global_gf <- gf_sh_ecm_tax |>
    dplyr::filter(!is.na(fungal_species)) |>
    dplyr::inner_join(gf_root_canadian_hosts, by = "sample_ID",
                      relationship = "many-to-many") |>
    dplyr::filter(host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, fungal_species)

  # Host → genus (global, GF, derived from species)
  host_genus_global_gf <- host_species_global_gf |>
    dplyr::mutate(fungal_genus = sub(" .*", "", fungal_species)) |>
    dplyr::distinct(host_species, fungal_genus)

} else {
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
# The search mirrors 03_genbank.R but removes AND "Canada"[Country].

if (file.exists(global_gb_ckpt)) {
  gb_global_meta <- readr::read_csv(global_gb_ckpt, show_col_types = FALSE)
} else {

  if (Sys.getenv("ENTREZ_KEY") == "")
    warning("ENTREZ_KEY not set in .Renviron — the GenBank global step will be ",
            "rate-limited. Set ENTREZ_KEY for faster fetching.", call. = FALSE)

  our_genera_list <- unique(trimws(emf$genus))

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

      # Retrieve metadata for ALL matching records (no per-genus cap). Pagination
      # goes through the NCBI history server (web_history), which — unlike a
      # historyless retstart search — has no ~10,000-record ceiling, so even the
      # most heavily sequenced genera are retrieved in full.
      batch_size <- 200L
      starts     <- seq(0L, search$count - 1L, by = batch_size)
      meta_list  <- vector("list", length(starts))

      for (k in seq_along(starts)) {
        summ <- tryCatch(
          rentrez::entrez_summary(
            db          = "nuccore",
            web_history = search$web_history,
            retstart    = starts[k],
            retmax      = batch_size
          ),
          error = function(e) NULL
        )
        if (!is.null(summ)) {
          if (inherits(summ, "esummary")) summ <- list(summ)
          meta_list[[k]] <- dplyr::bind_rows(lapply(summ, function(s) {
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
      warning(sprintf("  Genus '%s' query failed: %s", genus_name, conditionMessage(e)))
      NULL
    })
  }

  genus_results <- vector("list", length(our_genera_list))
  for (i in seq_along(our_genera_list)) {
    genus_results[[i]] <- query_genus_globally(our_genera_list[i])
    Sys.sleep(0.5)
  }

  gb_global_meta <- dplyr::bind_rows(Filter(Negate(is.null), genus_results))

  readr::write_csv(gb_global_meta, global_gb_ckpt)
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

  # Q1 (GenBank): which Canadian host species appear as hosts globally in GB?
  gb_q1_hosts <- unique(gb_global_meta$host_clean[
    !is.na(gb_global_meta$host_clean) &
      gb_global_meta$host_clean %in% host_tbl$species
  ])

  # Long-form GenBank global pairs for the interactions Rds (Part C)
  # Host → genus (global, GenBank, from genus_queried field)
  host_genus_global_gb <- gb_global_meta |>
    dplyr::filter(!is.na(host_clean), host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, fungal_genus = genus_queried)

  # Host → named EcM species (global, GenBank, from organism field)
  host_species_global_gb <- gb_global_meta |>
    dplyr::filter(!is.na(host_clean), host_clean %in% host_tbl$species,
                  !is.na(organism), nzchar(trimws(organism))) |>
    dplyr::mutate(fungal_species = trimws(organism)) |>
    dplyr::filter(!grepl("^(Fungi|Basidiomycota|Ascomycota|uncultured|sp\\.|NA)$",
                         fungal_species, ignore.case = TRUE),
                  nchar(fungal_species) > 3L) |>
    dplyr::distinct(host_species = host_clean, fungal_species)

} else {
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

# =============================================================================
# Summary tables
# =============================================================================

# Mycobiont-focal association counts (named fungal species -> host species),
# for the Eltonian "mycobiont side" paragraph. Canada scope from the Canadian
# EMF dataset; global scope from GlobalFungi root samples + GenBank worldwide
# (species_host_canada / species_host_global are assembled in Part C).
myco_counts <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(c(n = 0L, mx = 0L))
  cts <- dplyr::count(df, fungal_species, name = "n")
  c(n = nrow(cts), mx = max(cts$n))
}
myco_can  <- myco_counts(species_host_canada)
myco_glob <- myco_counts(species_host_global)

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
    "Canadian EcM species with documented global host associations (GlobalFungi root samples)",
    # Mycobiont-focal: named fungal species with >= 1 documented host species, by scope
    "Named EcM fungal species with >= 1 documented host species (Canada scope)",
    "Max host species documented for any named EcM fungal species (Canada scope)",
    "Named EcM fungal species with >= 1 documented host species (global scope; GlobalFungi root + GenBank worldwide)",
    "Max host species documented for any named EcM fungal species (global scope)"
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
    ifelse(!is.null(q3), nrow(q3), NA_real_),
    unname(myco_can["n"]),
    unname(myco_can["mx"]),
    unname(myco_glob["n"]),
    unname(myco_glob["mx"])
  )
)

readr::write_csv(eltonian_summary,
                 file.path(paths$out_eltonian, "eltonian_summary.csv"))
