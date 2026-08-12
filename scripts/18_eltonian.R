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
# IMPORTANT NOTE on the root-evidence rule:
#   Every fungus-host association reported by this script rests on root-derived
#   evidence, by two different mechanisms:
#     GlobalFungi — samples have a 'dominant_plant_species' field recording
#       which plant species dominates the sample site. This can only be
#       interpreted as the likely mycorrhizal host when sample_type == "root"
#       (i.e. a direct root sample rather than bulk soil or another substrate).
#       The same logic applies to 'other_plant_species'.
#     GenBank — there is no controlled sample-type field, so we require the
#       free-text isolation_source to name root or ectomycorrhizal material
#       (keyword rule root_kw, Step A2-0). Records from soil, rhizosphere,
#       non-root tissue, a habitat description, or a blank field are excluded
#       from all association analyses; they remain in the taxonomic and
#       geographic assessments carried out by the other scripts.
#   This constraint is applied throughout both Part A and Part B.
#
#   REVISION HISTORY. Until 2026-08 the GenBank host field was used with no
#   tissue test at all, while GlobalFungi was root-only. The two databases were
#   therefore held to different evidential standards, and because both feed
#   `host_long`, the asymmetry propagated to every downstream statistic. The
#   same revision fixed a silent data loss in the host-name cleaner, which kept
#   only the first plant named in a multi-plant label (see Step A2b).
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
#   data_derived/eltonian/eltonian_host_matching.csv — the ANALYSIS table: one
#       row per (record x named plant), root-derived evidence only
#   data_derived/eltonian/eltonian_host_names_metadata.csv — the DESCRIPTIVE
#       inventory: the same rows without the root-evidence restriction, plus a
#       `root_evidence` flag. Backs the host-name quality tables (S2, S3)
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
#   data_derived/eltonian/eltonian_pair_confidence.csv — how many observed
#       host x fungus pairs rest on samples that named a single candidate plant
#       ("unambiguous") versus several ("ambiguous"), for both the SH-code and
#       the named-species pair sets
#   data_derived/eltonian/eltonian_sensitivity_confidence.csv — the headline
#       Canadian numbers recomputed from unambiguous (single-host-sample)
#       evidence only, as a sensitivity check on the pair set above
#   data_derived/eltonian/eltonian_sample_type_tally_canada.csv — diagnostic: sample_type
#       composition of EcM-positive Canadian GF samples with a
#       dominant_plant_species entry
#   data_derived/eltonian/eltonian_sample_type_tally_canada_all.csv — diagnostic:
#       sample_type composition of ALL EcM-positive Canadian GF samples (no
#       dominant_plant_species restriction), for comparison with the worldwide
#       GlobalFungi tally written by 02_globalfungi.R
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
# imply a direct mycorrhizal association.
#
# REVISION 2026-08 (root-evidence rule). The GenBank host_taxon field was
# previously used UNCONDITIONALLY, i.e. with no tissue test at all. That held
# the two databases to different evidential standards: GlobalFungi host links
# were root-samples-only while GenBank host links came from soil, rhizosphere,
# habitat descriptions or a blank tissue field on the same footing as a root
# tip. Both branches feed `host_long`, and `matched_interactions` (Step A3) is
# the sole ancestor of every downstream host-association statistic, so the
# asymmetry propagated everywhere. GenBank records are now restricted to those
# whose free-text isolation_source names root or ectomycorrhizal material,
# using the SAME keyword rule as the tissue tally in Step A2c below. A fungus
# sequenced from soil beside a tree is not evidence that the tree was its
# partner.

# ---- A2-0. Tissue keyword rules (shared by Step A2 and the tally in A2c) -----
# Defined here, ABOVE their first use, so the host-record restriction in A2 and
# the descriptive tally in A2c can never drift apart. Word boundaries (\\b) are
# used where a bare substring would produce false positives -- most importantly
# "stem", which is a substring of "root system" (184 records).
#
#   root_kw   root / ectomycorrhizal material
#   soil_kw   soil, rhizosphere, or forest-floor duff
#   tiss_kw   plainly non-root fungal or plant tissue (sporocarps, leaves,
#             needles, bark, stems, algal thallus)
root_kw <- "root|mycorrhiz|\\becm\\b"
soil_kw <- "soil|rhizosphere|\\bduff\\b"
tiss_kw <- paste0("basidiocarp|sporocarp|\\bleaf\\b|\\bleaves\\b|\\bneedle",
                  "s?\\b|\\bbark\\b|\\bstem\\b|thallus|macroalga")

# How many plant names does a free-text label contain? Labels separate plants
# with commas or semicolons. A missing label contributes no candidate plants.
n_listed <- function(x) {
  ifelse(is.na(x), 0L, lengths(strsplit(x, "[,;]")))
}

# Root samples only for GlobalFungi dominant/other plant species fields.
# n_hosts_in_record counts ALL candidate plants named for the sample, i.e. the
# dominant_plant_species and other_plant_species fields added together. It is
# the sample-level count, not the per-field count, and it is what the
# confidence flag below keys on -- see the note in Step A2b.
emf_gf_root <- dplyr::filter(emf,
                               source == "GlobalFungi",
                               sample_type == "root") |>
  dplyr::mutate(n_hosts_in_record = n_listed(dominant_plant_species) +
                                    n_listed(other_plant_species))

# host_long_all is the FULL inventory of host strings carried by the metadata:
# every GlobalFungi root-sample plant field, and every GenBank host field
# whatever the tissue. It is kept because the paper asks two different
# questions of these strings:
#
#   * The host-NAME inventory -- how many distinct plant names appear at all,
#     and how many are ornamentals, misspellings, or otherwise unusable -- is a
#     statement about metadata quality, and should see every recorded name.
#     Supplemental Tables S2 and S3 are built from this.
#   * The host-ASSOCIATION analysis is a statement about evidence, and takes
#     only root-derived records (the `root_evidence` flag below).
#
# Building one table and filtering it guarantees the two can never drift apart.
host_long_all <- dplyr::bind_rows(

  # dominant_plant_species: root samples only
  emf_gf_root |>
    dplyr::filter(!is.na(dominant_plant_species)) |>
    dplyr::select(sh_code, genus, species, n_hosts_in_record,
                  dominant_plant_species) |>
    dplyr::rename(host_raw = dominant_plant_species) |>
    dplyr::mutate(host_field = "dominant_plant_species",
                  root_evidence = TRUE),

  # other_plant_species: root samples only
  emf_gf_root |>
    dplyr::filter(!is.na(other_plant_species)) |>
    dplyr::select(sh_code, genus, species, n_hosts_in_record,
                  other_plant_species) |>
    dplyr::rename(host_raw = other_plant_species) |>
    dplyr::mutate(host_field = "other_plant_species",
                  root_evidence = TRUE),

  # host_taxon: the GenBank host field, for every record that has one. The
  # root-evidence test is RECORDED here rather than applied here, so that the
  # descriptive inventory keeps the full set (see the note above). Only records
  # whose isolation_source names root or ectomycorrhizal material carry
  # root_evidence = TRUE and so reach the association analyses; the rest (soil,
  # rhizosphere, habitat descriptions, non-root tissue, or a blank field) are
  # excluded from every host-association result, though they remain in the
  # taxonomic and geographic assessments elsewhere in the paper.
  # A GenBank record has only one host field, so its candidate-plant count is
  # simply the number of names in that field (almost always exactly one).
  emf |>
    dplyr::filter(source == "GenBank", !is.na(host_taxon)) |>
    dplyr::mutate(n_hosts_in_record = n_listed(host_taxon),
                  root_evidence = grepl(root_kw, isolation_src,
                                        ignore.case = TRUE)) |>
    dplyr::select(sh_code, genus, species, n_hosts_in_record, root_evidence,
                  host_taxon) |>
    dplyr::rename(host_raw = host_taxon) |>
    dplyr::mutate(host_field = "host_taxon")
)

# ---- A2b. Split multi-species host labels ------------------------------------
# Some metadata fields name several plants for a single sample, e.g.
# "Picea mariana, Picea glauca, Pinus banksiana". canonicalize_host() keeps only
# the first two words of a string, so before this revision every plant after the
# first was silently discarded -- Picea glauca, for instance, was named in 36 of
# the 44 Canadian root samples and captured from none of them.
#
# We now split each label on commas and semicolons and treat every named plant
# as its own candidate host. Because a fungus found in a sample that listed
# eight plants is much weaker evidence than one found where only a single plant
# was named, we record how many plants were named and derive a confidence
# label. Nothing is discarded; the weaker records are flagged so they can be
# reported separately or excluded in a sensitivity check.
#
#   host_items / host_item  the individual plant names cut out of host_raw
#   host_rank               position within this field's label (1 = first)
#   n_hosts_in_field        how many plants THIS field named
#   n_hosts_in_record       how many candidate plants the whole source record
#                           named (see below)
#   pair_confidence         "unambiguous" if the source record named exactly
#                           one candidate plant, "ambiguous" if several
#
# WHY THE CONFIDENCE FLAG USES n_hosts_in_record, NOT n_hosts_in_field.
# A GlobalFungi root sample carries two host fields. Reading them separately
# would call a sample "unambiguous" whenever dominant_plant_species named a
# single plant -- even when other_plant_species listed seven more at the same
# sample. That is exactly the ambiguity the flag exists to record, so the count
# is taken over the whole record. In this dataset 8 of the 44 Canadian root
# samples name a single plant; 31 name eight and 5 name eleven.
#
# RECOVERING THE PRE-REVISION BEHAVIOUR: filter(host_rank == 1 &
# n_hosts_in_field == 1) on host_long returns exactly the rows the script
# produced before the splitting step was added.

host_long_all <- host_long_all |>
  dplyr::mutate(
    host_items       = strsplit(host_raw, "[,;]"),
    n_hosts_in_field = lengths(host_items)
  ) |>
  tidyr::unnest_longer(host_items,
                       values_to  = "host_item",
                       indices_to = "host_rank") |>
  dplyr::mutate(host_item = trimws(host_item)) |>
  dplyr::filter(nzchar(host_item))

host_long_all <- host_long_all |>
  dplyr::mutate(
    host_clean      = clean_host_name(host_item),
    host_genus      = sub("^(\\S+).*", "\\1", host_clean),
    matched         = host_clean %in% host_tbl$species,
    match_genus     = host_genus %in% sub("^(\\S+).*", "\\1", host_tbl$species),
    pair_confidence = dplyr::if_else(n_hosts_in_record == 1L,
                                     "unambiguous", "ambiguous")
  )

# The descriptive inventory: every host name the metadata carries, with the
# root_evidence flag showing which ones reach the association analyses.
# Supplemental Tables S2 and S3 are built from this file.
readr::write_csv(host_long_all,
                 file.path(paths$out_eltonian,
                           "eltonian_host_names_metadata.csv"))

# The analysis table: root-derived evidence only. Everything downstream of this
# point uses host_long, so the root-evidence rule applies to every reported
# association statistic, map, and matrix.
host_long <- dplyr::filter(host_long_all, root_evidence)

n_matched       <- sum(host_long$matched,     na.rm = TRUE)
n_genus_matched <- sum(host_long$match_genus, na.rm = TRUE)

readr::write_csv(host_long,
                 file.path(paths$out_eltonian, "eltonian_host_matching.csv"))

# ---- A2c. Diagnostic: sample_type tally among Canadian GF samples with a -----
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

# Denominator for the diagnostic tally and for the host-information-coverage
# statistic reported in the manuscript (Eltonian intro paragraph): all
# distinct Canadian GlobalFungi samples with >= 1 EcM fungal SH code detected
# (i.e. every sample_type, not just those with dominant_plant_species filled
# in -- that restriction is applied below).
n_gf_samples_total <- emf |>
  dplyr::filter(source == "GlobalFungi") |>
  dplyr::distinct(sample_ID) |>
  nrow()

sample_type_tally_canada <- emf |>
  dplyr::filter(source == "GlobalFungi", !is.na(dominant_plant_species)) |>
  dplyr::distinct(sample_ID, sample_type) |>
  dplyr::count(sample_type, name = "n_samples") |>
  dplyr::arrange(dplyr::desc(n_samples))

readr::write_csv(sample_type_tally_canada,
                 file.path(paths$out_eltonian, "eltonian_sample_type_tally_canada.csv"))

# Companion tally over ALL EcM-positive Canadian GlobalFungi samples, i.e.
# without the "has a dominant_plant_species entry" restriction applied above.
# Added 2026-08 so the manuscript can state what proportion of the Canadian
# EcM-bearing samples are root samples (3%) against the worldwide GlobalFungi
# figure (19%, tallied in 02_globalfungi.R) -- the point being that thin root
# coverage is a Canadian sampling gap rather than a property of the database.
sample_type_tally_canada_all <- emf |>
  dplyr::filter(source == "GlobalFungi") |>
  dplyr::distinct(sample_ID, sample_type) |>
  dplyr::count(sample_type, name = "n_samples") |>
  dplyr::mutate(pct = round(100 * n_samples / sum(n_samples), 1)) |>
  dplyr::arrange(dplyr::desc(n_samples))

readr::write_csv(sample_type_tally_canada_all,
                 file.path(paths$out_eltonian,
                           "eltonian_sample_type_tally_canada_all.csv"))

# ---- A2d. Diagnostic: tissue-type tally among Canadian GenBank records with -
#           host information (host_taxon)
# Mirrors A2c above, but for GenBank. Unlike GlobalFungi, GenBank has no
# controlled-vocabulary sample_type field -- the closest analogue is the
# free-text 'isolation_src' field (NCBI /isolation_source qualifier), which is
# heterogeneous (65 distinct strings in this dataset, e.g. "root system",
# "ectomycorrhiza", "soil adjacent to Pinus albicaulis", "40 year-old
# Douglas-fir stand").
#
# Since the 2026-08 revision this tally is no longer purely descriptive: Step
# A2 applies exactly the same root_kw test to decide which GenBank records may
# contribute a host association. The categories below therefore also document
# what was kept and what was excluded. Records counted under "Root /
# ectomycorrhizal tissue" and "Root and soil both mentioned" are the ones
# retained (they match root_kw); the other three categories are excluded.
#
# Tissue categories are assigned by keyword search (case-insensitive) on
# isolation_src, in this priority order:
#   Root / ectomycorrhizal tissue  - matches root_kw and no soil keyword
#   Root and soil both mentioned   - matches BOTH root_kw and soil_kw (e.g.
#                                     "mycorrhizal root tips of Pinus
#                                     albicaulis ... and soil adjacent to
#                                     Pinus albicaulis"). Renamed in 2026-08
#                                     from "Mixed root + soil", which implied
#                                     a mixed sample; in this dataset all such
#                                     strings describe root material collected
#                                     from a site where soil was also sampled,
#                                     so they are treated as root evidence.
#   Soil / rhizosphere / duff      - matches soil_kw, no root keyword
#   Non-root tissue                - plainly non-root fungal or plant tissue
#                                     (sporocarp, leaf, needle, bark, stem,
#                                     algal thallus); matches tiss_kw only
#   Habitat description only       - text present but naming no material at
#                                     all (stand age, forest type, locality),
#                                     so the tissue sampled is unknown rather
#                                     than known to be non-root
#   Not recorded                   - isolation_src missing entirely
#
# Counted at the record level (one row = one GenBank accession; confirmed
# 1:1 in this dataset, unlike the long-format GF rows in A2c).

genbank_canada <- dplyr::filter(emf, source == "GenBank")
n_gb_total <- nrow(genbank_canada)

genbank_with_host <- dplyr::filter(genbank_canada, !is.na(host_taxon))
n_gb_with_host <- nrow(genbank_with_host)

genbank_tissue <- genbank_with_host |>
  dplyr::mutate(
    has_root   = grepl(root_kw, isolation_src, ignore.case = TRUE),
    has_soil   = grepl(soil_kw, isolation_src, ignore.case = TRUE),
    has_tissue = grepl(tiss_kw, isolation_src, ignore.case = TRUE),
    tissue_category = dplyr::case_when(
      is.na(isolation_src)        ~ "Not recorded",
      has_root & has_soil         ~ "Root and soil both mentioned",
      has_root                    ~ "Root / ectomycorrhizal tissue",
      has_soil                    ~ "Soil / rhizosphere / duff",
      has_tissue                  ~ "Non-root tissue",
      TRUE                        ~ "Habitat description only"
    )
  )

# Three roll-up rows are appended so the manuscript can quote the grouped
# figures without re-deriving them from the six categories:
#   "Root-derived (retained ...)"  = the records that pass the root_kw test in
#                                    Step A2 and so contribute host associations
#   "Tissue provenance unknown"    = blank field + habitat-only descriptions;
#                                    these say nothing about the material
#   "Demonstrably non-root"        = soil / rhizosphere / duff + non-root tissue
# The second and third are the correction to the submitted text, which lumped
# them together as "likely derived from samples other than ectomycorrhiza or
# host roots" -- a claim the habitat-only strings do not support.
.tcat <- function(...) sum(genbank_tissue$tissue_category %in% c(...))

genbank_tissue_tally_canada <- dplyr::bind_rows(
  tibble::tibble(category = "Total EcM fungal records",      n = n_gb_total),
  tibble::tibble(category = "Records with host information", n = n_gb_with_host),
  genbank_tissue |>
    dplyr::count(tissue_category, name = "n") |>
    dplyr::rename(category = tissue_category) |>
    dplyr::arrange(dplyr::desc(n)),
  tibble::tibble(
    category = c("Root-derived (retained for host associations)",
                 "Tissue provenance unknown (blank or habitat only)",
                 "Demonstrably non-root (soil, rhizosphere, duff, other tissue)"),
    n = c(.tcat("Root / ectomycorrhizal tissue", "Root and soil both mentioned"),
          .tcat("Not recorded", "Habitat description only"),
          .tcat("Soil / rhizosphere / duff", "Non-root tissue")))
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

# A given host x fungus pair may be supported by several records of differing
# quality (see the confidence flag added in Step A2b). We take the best
# available evidence: if ANY supporting record came from a sample naming a
# single plant, the pair is "unambiguous"; if the pair rests entirely on
# samples that named several candidate plants, it is "ambiguous".
sh_pairs <- matched_interactions |>
  dplyr::filter(!is.na(sh_code)) |>   # genus-resolved GenBank rows carry sh_code = NA
  dplyr::group_by(host_species = host_clean, sh_code) |>
  dplyr::summarise(
    pair_confidence = dplyr::if_else(any(pair_confidence == "unambiguous"),
                                     "unambiguous", "ambiguous"),
    .groups = "drop"
  )

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
  # Same best-available-evidence rule as sh_pairs above, applied one level up:
  # several SH codes can resolve to the same named species, so a host x species
  # pair is unambiguous if any of its supporting SH-level pairs was.
  dplyr::group_by(host_species, fungal_species) |>
  dplyr::summarise(
    pair_confidence = dplyr::if_else(any(pair_confidence == "unambiguous"),
                                     "unambiguous", "ambiguous"),
    .groups = "drop"
  )

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
  # Alongside the occurrence count we carry the confidence flag (added 2026-08),
  # applying the same best-available-evidence rule used for sh_pairs and
  # species_pairs_canada. n_records_unambiguous records how many of the
  # supporting records came from a single-host sample, so a reader can see how
  # thin the evidence is for a pair that is flagged unambiguous on one record.
  dplyr::group_by(host_species, fungal_species) |>
  dplyr::summarise(
    n_occurrences        = dplyr::n(),
    n_records_unambiguous = sum(pair_confidence == "unambiguous"),
    pair_confidence      = dplyr::if_else(n_records_unambiguous > 0L,
                                          "unambiguous", "ambiguous"),
    .groups = "drop"
  ) |>
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

# ---- A3c-bis. Confidence composition of the Canadian pair set ---------------
# Reported in the manuscript so readers can see how much of the fungus-host
# evidence rests on samples where several candidate plants were listed. Both
# pair sets are tallied: the SH-code-level set (the finest fungal resolution)
# and the named-species-level set (the one whose fill rate the manuscript
# reports as the headline interaction matrix).
pair_conf_tally <- dplyr::bind_rows(
  sh_pairs |>
    dplyr::count(pair_confidence, name = "n_pairs") |>
    dplyr::mutate(pair_set = "host x SH code"),
  species_occurrence_counts |>
    dplyr::count(pair_confidence, name = "n_pairs") |>
    dplyr::mutate(pair_set = "host x named fungal species")
) |>
  dplyr::group_by(pair_set) |>
  dplyr::mutate(pct = round(100 * n_pairs / sum(n_pairs), 1)) |>
  dplyr::ungroup() |>
  dplyr::select(pair_set, pair_confidence, n_pairs, pct)

readr::write_csv(pair_conf_tally,
                 file.path(paths$out_eltonian, "eltonian_pair_confidence.csv"))

# ---- A3c-ter. Sensitivity analysis: unambiguous evidence only ---------------
# Reviewers will reasonably ask what happens if the associations that rest
# entirely on multi-plant samples are simply dropped. This block recomputes the
# headline Canadian numbers from single-host records only. It writes its own
# small CSV and changes nothing above: the reported analysis remains the full
# set, with this as the stated sensitivity check.
#
# "Unambiguous" here means the source record named exactly one candidate plant:
# every retained GenBank root record, plus the 8 Canadian GlobalFungi root
# samples whose two plant fields between them name a single species.
sens_interactions <- matched_interactions |>
  dplyr::filter(n_hosts_in_record == 1L)

sens_species_pairs <- sens_interactions |>
  dplyr::filter(!is.na(sh_code)) |>
  dplyr::select(-species) |>
  dplyr::left_join(sh_lookup |> dplyr::select(sh_code, species), by = "sh_code") |>
  dplyr::mutate(
    fungal_species = dplyr::if_else(
      !is.na(species) & !grepl("_sp$", species),
      trimws(gsub("_", " ", species)), NA_character_)
  ) |>
  dplyr::filter(!is.na(fungal_species)) |>
  dplyr::distinct(host_species = host_clean, fungal_species)

n_hosts_sens   <- dplyr::n_distinct(sens_interactions$host_clean)
n_myco_sens    <- dplyr::n_distinct(sens_species_pairs$fungal_species)
n_pairs_sens   <- nrow(sens_species_pairs)

eltonian_sensitivity <- tibble::tibble(
  metric = c(
    "Canadian host species with >= 1 documented association (all evidence)",
    "Canadian host species with >= 1 documented association (unambiguous evidence only)",
    "Eltonian host shortfall, hosts lacking any association (all evidence)",
    "Eltonian host shortfall, hosts lacking any association (unambiguous evidence only)",
    "Named EcM fungal species with >= 1 documented host (all evidence)",
    "Named EcM fungal species with >= 1 documented host (unambiguous evidence only)",
    "Host x named-species pairs (all evidence)",
    "Host x named-species pairs (unambiguous evidence only)"
  ),
  value = c(
    dplyr::n_distinct(matched_interactions$host_clean),
    n_hosts_sens,
    n_host_species - dplyr::n_distinct(matched_interactions$host_clean),
    n_host_species - n_hosts_sens,
    dplyr::n_distinct(species_pairs_canada$fungal_species),
    n_myco_sens,
    n_pairs_total,
    n_pairs_sens
  )
)

readr::write_csv(eltonian_sensitivity,
                 file.path(paths$out_eltonian, "eltonian_sensitivity_confidence.csv"))

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

  # REVISION 2026-08 (root-evidence rule, global scope). The Canadian GenBank
  # branch in Step A2 keeps only records whose isolation_source names root or
  # ectomycorrhizal material. The same test is applied here, for the same
  # reason and with the same keyword pattern. Without it the two scopes would
  # be defined differently from one another -- the Canadian figures strict and
  # the global ones mixed -- which is worse than applying no filter at all.
  gb_global_root <- gb_global_meta |>
    dplyr::filter(grepl(root_kw, isolation_src, ignore.case = TRUE))

  n_gb_global_with_host <- sum(!is.na(gb_global_root$host_clean))

  # Q3 (GenBank): per fungal genus, which Canadian host species globally?
  gb_q3 <- gb_global_root |>
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
  gb_q1_hosts <- unique(gb_global_root$host_clean[
    !is.na(gb_global_root$host_clean) &
      gb_global_root$host_clean %in% host_tbl$species
  ])

  # Long-form GenBank global pairs for the interactions Rds (Part C)
  # Host → genus (global, GenBank, from genus_queried field)
  host_genus_global_gb <- gb_global_root |>
    dplyr::filter(!is.na(host_clean), host_clean %in% host_tbl$species) |>
    dplyr::distinct(host_species = host_clean, fungal_genus = genus_queried)

  # Host → named EcM species (global, GenBank, from organism field)
  host_species_global_gb <- gb_global_root |>
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
    # Sample-level host-information coverage (Eltonian intro paragraph)
    "GlobalFungi samples in Canada with >= 1 EcM fungal SH code detected",
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
    n_gf_samples_total,
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
