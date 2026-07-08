# =============================================================================
# Wallacean Shortfall: Geographic Sampling Coverage
# =============================================================================
# How well distributed across Canada are the EcM fungal sampling locations,
# and how many locations has each taxon been recorded from?
#
# Low mean locations per taxon and a high proportion of taxa known from a
# single location are the primary indicators of the Wallacean shortfall.
#
# The analysis is run independently for three datasets to allow direct
# comparison:
#   GlobalFungi — all GF records with decimal-degree coordinates
#   GenBank     — GenBank records with parsed coordinates (canada_basis %in%
#                 c("both", "coordinates_only"); ~40.7% of GenBank records)
#   Combined    — union of the above two, deduplicated on
#                 (sh_code, genus, species, lat, lon) so co-located records
#                 from both sources are counted once
#
# The geographic sampling map is produced by 17_hutchinsonian.R
# (hutchinsonian_sampling_map.png) and 2-7_sampling_maps.R
# (sampling_map_gbif.png); it is not duplicated here.
#
# Outputs (data_derived/wallacean/):
#   wallacean_location_summary.csv     — unique location counts (3 datasets)
#   wallacean_sampling_intensity.csv   — locations-per-taxon stats (3 datasets)
#   wallacean_locs_per_sh.csv          — per-SH location counts (dataset col)
#   wallacean_locs_per_genus.csv       — per-genus location counts (dataset col)
#   wallacean_locs_per_species.csv     — per-species location counts (dataset col)
#   wallacean_global_gf_locs_per_species.csv — global GF grid-cell occupancy per
#                                       named species (30 arc-second cells), for
#                                       the GlobalFungi-only occupancy analyses
#                                       in EMF_shortfalls_SI.Rmd §S3 (Figure M,
#                                       Table NN, Figure YY)
#   wallacean_global_gf_locs_per_sh.csv — global GF grid-cell occupancy per SH
#                                       code (named + unnamed), for the same
#                                       §S3 analyses.
#
#   Both global-occupancy files are produced by a single Step 5, which reads
#   the full ~13 GB global SH abundance matrix once (not the smaller species
#   matrix) and is long-running — run separately on demand (same pattern as
#   the Table S1 global comparator stages). Species-level occupancy is
#   derived from our own UNITE SH-code taxonomy (paths$unite_taxonomy), not
#   from GlobalFungi's bundled species-level abundance matrix — see the Step 5
#   header comment below for why.
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- 0. Build three analysis datasets ---------------------------------------

ts("Building analysis datasets...")

gf <- dplyr::filter(emf, source == "GlobalFungi", coord_in_canada == TRUE)

# GenBank: same filter — coord_in_canada == TRUE ensures coordinates were
# validated within the GADM Canada boundary at dataset assembly.
gb_geo <- dplyr::filter(emf, source == "GenBank", coord_in_canada == TRUE)

combined <- dplyr::bind_rows(
  dplyr::select(gf,     sh_code, genus, species, lat, lon),
  dplyr::select(gb_geo, sh_code, genus, species, lat, lon)
) |> dplyr::distinct()

n_gf_no_coords <- sum(emf$source == "GlobalFungi" & is.na(emf$lat))
n_gb_total     <- sum(emf$source == "GenBank")
n_gb_no_coords <- sum(emf$source == "GenBank" & (is.na(emf$lat) | is.na(emf$lon)))

ts(sprintf("  GlobalFungi records w/ coords:  %d", nrow(gf)))
ts(sprintf("  GlobalFungi records w/o coords: %d (excluded)", n_gf_no_coords))
ts(sprintf("  GenBank records w/ coords: %d / %d (%.1f%%)",
           nrow(gb_geo), n_gb_total, 100 * nrow(gb_geo) / n_gb_total))
ts(sprintf("  GenBank records w/o coords:     %d (excluded)", n_gb_no_coords))
ts(sprintf("  Combined (deduplicated rows):   %d", nrow(combined)))

datasets      <- list(GlobalFungi = gf, GenBank = gb_geo, Combined = combined)
dataset_names <- names(datasets)

# ---- 1. Unique sampling locations -------------------------------------------

ts("Counting unique sampling locations per dataset...")

location_summary <- dplyr::bind_rows(lapply(dataset_names, function(nm) {
  d <- datasets[[nm]]
  tibble::tibble(
    dataset  = nm,
    n_records_with_coords = nrow(d),
    n_unique_locations    = dplyr::n_distinct(d$lat, d$lon)
  )
}))

readr::write_csv(location_summary,
                 file.path(paths$out_wallacean, "wallacean_location_summary.csv"))
ts("  Saved wallacean_location_summary.csv")

# ---- 2. Locations per taxon -------------------------------------------------

ts("Calculating locations per taxon across three datasets...")

# Helper: distinct (taxon, lat, lon) → count of unique locations per taxon
count_locs_per_taxon <- function(df, taxon_col) {
  df |>
    dplyr::filter(!is.na(.data[[taxon_col]])) |>
    dplyr::distinct(.data[[taxon_col]], lat, lon) |>
    dplyr::count(.data[[taxon_col]], name = "n_locations") |>
    dplyr::arrange(dplyr::desc(n_locations))
}

summarise_locs <- function(df, level_name) {
  tibble::tibble(
    taxonomic_level     = level_name,
    n_taxa              = nrow(df),
    mean_locs           = round(mean(df$n_locations), 2),
    median_locs         = median(df$n_locations),
    max_locs            = max(df$n_locations),
    min_locs            = min(df$n_locations),
    n_single_location   = sum(df$n_locations == 1L),
    pct_single_location = round(100 * mean(df$n_locations == 1L), 1)
  )
}

all_sh      <- list()
all_genus   <- list()
all_species <- list()
intensity   <- list()

for (nm in dataset_names) {
  d  <- datasets[[nm]]
  sh <- count_locs_per_taxon(d, "sh_code")
  gn <- count_locs_per_taxon(d, "genus")
  sp <- count_locs_per_taxon(
    dplyr::filter(d, !is.na(species), !grepl("_sp$", species)), "species"
  )
  ts(sprintf("  %-12s — SH: %d  Genera: %d  Species: %d",
             nm, nrow(sh), nrow(gn), nrow(sp)))

  all_sh[[nm]]    <- dplyr::mutate(sh, dataset = nm)
  all_genus[[nm]] <- dplyr::mutate(gn, dataset = nm)
  all_species[[nm]] <- dplyr::mutate(sp, dataset = nm)

  intensity[[nm]] <- dplyr::bind_rows(
    summarise_locs(sh, "SH code"),
    summarise_locs(gn, "Genus"),
    summarise_locs(sp, "Species")
  ) |>
    dplyr::mutate(dataset = nm) |>
    dplyr::relocate(dataset)
}

locs_sh      <- dplyr::bind_rows(all_sh)
locs_genus   <- dplyr::bind_rows(all_genus)
locs_species <- dplyr::bind_rows(all_species)
sampling_intensity <- dplyr::bind_rows(intensity)

readr::write_csv(locs_sh,
                 file.path(paths$out_wallacean, "wallacean_locs_per_sh.csv"))
readr::write_csv(locs_genus,
                 file.path(paths$out_wallacean, "wallacean_locs_per_genus.csv"))
readr::write_csv(locs_species,
                 file.path(paths$out_wallacean, "wallacean_locs_per_species.csv"))
readr::write_csv(sampling_intensity,
                 file.path(paths$out_wallacean, "wallacean_sampling_intensity.csv"))
ts("  Saved wallacean_locs_per_*.csv and wallacean_sampling_intensity.csv")

# NOTE: Per-taxon locations-per-taxon histograms (formerly Sections 3 and 4,
# producing figures/wallacean_histograms.png and
# figures/wallacean_histograms_combined.png) have been removed. The streamlined
# Wallacean section of EMF_shortfalls_SI.Rmd (§S3) builds its own GlobalFungi-
# only, grid-cell-based occupancy figures (Figure M, Figure YY) directly from
# the CSVs above and from the global-occupancy files written in Step 5 below.

# ---- 5. Global GF locations per Canadian SH code and named species ----------
#
# Former design (pre-2026-06-28): two separate steps. Step 5 matched named
# species by STRING NAME against GlobalFungi's own bundled species-level
# abundance matrix (`GlobalFungi_5_species_abundance_ITS1_ITS2.txt`) — i.e.
# GlobalFungi's internal SH-to-taxonomy pipeline, not a join through our
# pinned UNITE build. Step 6 matched SH codes EXACTLY by code against the
# SH-level matrix, going through our own UNITE lookup as everywhere else in
# this pipeline does.
#
# That split was an inconsistency: every other taxonomy join in this pipeline
# (0-3, 0-4, and Step 5's own SH-code population below) decodes through our
# pinned UNITE build (`paths$unite_taxonomy`) specifically so that GlobalFungi
# and GenBank SH codes are interpreted against an identical reference — see
# the SH-build-mismatch discussion in 02_globalfungi.R. Trusting
# GlobalFungi's bundled species names instead breaks that invariant. An
# empirical check (2026-06-28) found 98.6% name agreement with GlobalFungi's
# species columns and no detectable build-mismatch signature in the residual
# — but there was no validation gate, so a future GlobalFungi release built
# against a different UNITE version could silently degrade the species-level
# output (more species falsely scored as zero global locations) without
# raising an error.
#
# Fix: derive the named-species -> SH-code crosswalk from our own UNITE
# lookup table (Step 5b) instead of from GlobalFungi's species columns, then
# query the SH-level matrix ONCE (Step 5c) for the union of:
#   (a) Canadian SH codes (Step 5a) — at least one georeferenced Canadian
#       record (combined GlobalFungi + GenBank); an SH code with none can
#       never contribute a non-zero count to a Canada-side location tally.
#   (b) every globally-known SH code carrying a named Canadian species'
#       UNITE epithet (Step 5b) — NOT restricted to SH codes actually
#       detected in Canada, because a species' global range can include
#       SH-code variants we never sampled here; restricting to Canada-
#       detected codes would understate global occupancy relative to the
#       former GlobalFungi-species-table approach, which pooled across all
#       of a species' SH-code variants worldwide.
# This also removes the need for a second, separate read of the ~3.8 GB
# species matrix — both outputs now come from one pass over the ~13 GB
# SH-level matrix.
#
# Outputs (names/columns unchanged from the former two-step design):
#   data_derived/wallacean/wallacean_global_gf_locs_per_species.csv
#     species       — UNITE species epithet (underscores, matching emf$species)
#     n_locs_global — distinct 30 arc-second grid cells (~1 km²) globally,
#                     unioned across all of that species' SH-code variants
#   data_derived/wallacean/wallacean_global_gf_locs_per_sh.csv
#     sh_code       — UNITE SH code (named + unnamed; e.g. "SH0969952.10FU")
#     n_locs_global — distinct 30 arc-second grid cells globally
#
# This step reads ~2,000-3,000 columns x ~hundreds of thousands of samples in
# the global SH matrix. fread must scan every byte of the ~13 GB file
# regardless of column selection, so this is long-running — guarded by the
# same sentinel-file pattern as every other expensive step in this pipeline,
# and intended to be run once, on demand, rather than as part of routine
# pipeline re-runs (consistent with the "run separately by the user" pattern
# already used for the Table S1 global comparator stages; see
# scripts/13_wallacean_global_comparator.R). Because both output files are
# produced from a single shared read, the sentinel guard regenerates both
# together if either is missing.

out_global_gf_sp <- file.path(paths$out_wallacean,
                               "wallacean_global_gf_locs_per_species.csv")
out_global_gf_sh <- file.path(paths$out_wallacean,
                               "wallacean_global_gf_locs_per_sh.csv")

if (!file.exists(out_global_gf_sp) || file.size(out_global_gf_sp) == 0L ||
    !file.exists(out_global_gf_sh) || file.size(out_global_gf_sh) == 0L) {

  ts("Step 5: Counting global GF locations per Canadian SH code and named species...")

  # -- 5a. Per-SH-code population: Canadian SH codes (named + unnamed) with
  # >= 1 georeferenced Canadian record. See the header comment above for why
  # this population is restricted to coord_in_canada == TRUE.
  canada_sh <- emf |>
    dplyr::filter(!is.na(sh_code), coord_in_canada == TRUE) |>
    dplyr::distinct(sh_code) |>
    dplyr::pull(sh_code)
  ts(sprintf("  %d Canadian SH codes (per-SH population)", length(canada_sh)))

  # -- 5b. Per-species crosswalk: every named Canadian species mapped to ALL
  # globally-known SH codes carrying that UNITE species epithet, read from
  # our own UNITE lookup (built in 0-3 Step 1) rather than from GlobalFungi's
  # bundled species table. sh_code is unique per row in this lookup (one
  # representative sequence per SH), so this is a clean many-to-one join key.
  named_species <- emf |>
    dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
    dplyr::distinct(species) |>
    dplyr::pull(species)
  ts(sprintf("  %d named Canadian species (per-species population)",
             length(named_species)))

  unite_lookup <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
  species_sh_xwalk <- unite_lookup |>
    dplyr::filter(species %in% named_species) |>
    dplyr::distinct(species, sh_code)
  ts(sprintf("  %d named species resolved to %d global SH codes via UNITE lookup",
             dplyr::n_distinct(species_sh_xwalk$species),
             dplyr::n_distinct(species_sh_xwalk$sh_code)))

  species_not_in_unite <- setdiff(named_species, species_sh_xwalk$species)
  if (length(species_not_in_unite) > 0L) {
    ts(sprintf("  WARNING: %d named species not found in UNITE lookup (will be zero-filled below):",
               length(species_not_in_unite)))
    for (sp in species_not_in_unite) ts(sprintf("    - %s", sp))
  }

  # -- 5c. Single column-selective read of the global SH abundance matrix ----
  target_sh <- union(canada_sh, species_sh_xwalk$sh_code)

  ts("  Reading GF SH abundance matrix header...")
  sh_header <- data.table::fread(
    paths$gf_sh_abundance, nrows = 0L, sep = "\t"
  )
  gf_sh_col_names <- names(sh_header)
  ts(sprintf("  GF SH matrix: %d SH-code columns", length(gf_sh_col_names) - 1L))

  matched_sh   <- intersect(target_sh, gf_sh_col_names)
  unmatched_sh <- setdiff(target_sh, gf_sh_col_names)
  ts(sprintf("  Matched: %d / %d target SH codes (unmatched, zero-filled below: %d)",
             length(matched_sh), length(target_sh), length(unmatched_sh)))

  if (length(matched_sh) == 0L) {
    stop("No SH codes matched GF SH matrix column names — check sh_code formatting.")
  }

  # Read only the relevant columns from the ~13 GB SH abundance matrix, via the
  # awk-streaming helper in 00_setup.R (avoids fread()'s 2^31-byte string limit).
  ts("  Extracting matched columns from GF SH abundance matrix (long-running)...")
  sh_mat <- read_big_tsv_subset(paths$gf_sh_abundance, c("sample_ID", matched_sh))
  ts(sprintf("  Read %d samples x %d SH codes", nrow(sh_mat), length(matched_sh)))

  # Pivot to long; keep only presence rows.
  ts("  Pivoting to long format and filtering to presence...")
  sh_long <- data.table::melt(
    sh_mat,
    id.vars         = "sample_ID",
    variable.name   = "sh_code",
    value.name      = "abundance",
    variable.factor = FALSE
  )
  sh_long <- sh_long[abundance > 0L, .(sample_ID, sh_code)]
  ts(sprintf("  %d sample x SH-code presence records", nrow(sh_long)))

  # Load GF metadata; apply the same quality filters used throughout.
  ts("  Loading GF metadata and applying quality filters...")
  gf_meta_all <- data.table::fread(
    paths$gf_metadata,
    sep    = "\t",
    select = c("sample_ID", "latitude", "longitude",
               "barcoding_region", "manipulated", "sample_type")
  )
  gf_meta_filt <- gf_meta_all[
    barcoding_region %in% c("ITS2", "ITSboth") &
    manipulated == "NO" &
    !sample_type %in% c("shoot", "air", "water", "sediment") &
    !is.na(latitude) & !is.na(longitude)
  ]
  ts(sprintf("  %d / %d global samples pass quality filters and have coordinates",
             nrow(gf_meta_filt), nrow(gf_meta_all)))

  # Join presence records to filtered metadata.
  sh_geo <- merge(
    sh_long, gf_meta_filt[, .(sample_ID, latitude, longitude)],
    by = "sample_ID", all = FALSE
  )
  ts(sprintf("  %d presence records after metadata join", nrow(sh_geo)))

  # 30 arc-second grid cell key (~1 km²), via the shared snap_30arcsec()
  # helper from 00_setup.R — consistent with the grid used in the SI Rmd's
  # SDM sufficiency tables.
  sh_geo[, grid_key := paste(snap_30arcsec(latitude), snap_30arcsec(longitude))]

  # -- 5d. Per-SH-code output (former Step 6) ---------------------------------
  global_locs_sh <- sh_geo[
    sh_code %in% canada_sh,
    .(n_locs_global = data.table::uniqueN(grid_key)),
    by = sh_code
  ]
  data.table::setorder(global_locs_sh, -n_locs_global)

  # Include Canadian SH codes with no match in the GF SH matrix, or matched
  # but with zero presence records after quality filtering (n_locs_global = 0).
  zero_sh <- tibble::tibble(
    sh_code       = setdiff(canada_sh, global_locs_sh$sh_code),
    n_locs_global = 0L
  )
  global_locs_sh_full <- dplyr::bind_rows(as.data.frame(global_locs_sh), zero_sh) |>
    dplyr::arrange(dplyr::desc(n_locs_global))

  ts(sprintf("  SH codes with >= 1 global location: %d",
             sum(global_locs_sh_full$n_locs_global >= 1L)))
  ts(sprintf("  SH codes with >= 30 global locations: %d",
             sum(global_locs_sh_full$n_locs_global >= 30L)))

  readr::write_csv(global_locs_sh_full, out_global_gf_sh)
  ts(sprintf("  Saved -> wallacean/wallacean_global_gf_locs_per_sh.csv (%d rows)",
             nrow(global_locs_sh_full)))

  # -- 5e. Per-species output (former Step 5) ---------------------------------
  # Attach species labels via the crosswalk BEFORE de-duplicating grid cells,
  # so that two different SH-code variants of the same species detected at
  # the same physical location count as one location, not two. sh_code is
  # unique per row in species_sh_xwalk, so this is a many-to-one join (one
  # species label per presence record); no cartesian expansion risk.
  sp_geo <- merge(sh_geo, species_sh_xwalk, by = "sh_code", all = FALSE)
  global_locs_sp <- sp_geo[
    ,
    .(n_locs_global = data.table::uniqueN(grid_key)),
    by = species
  ]
  data.table::setorder(global_locs_sp, -n_locs_global)

  # Include named species with no matched SH code present in the GF SH
  # matrix, or matched but with zero presence records (n_locs_global = 0).
  zero_sp <- tibble::tibble(
    species       = setdiff(named_species, global_locs_sp$species),
    n_locs_global = 0L
  )
  global_locs_sp_full <- dplyr::bind_rows(as.data.frame(global_locs_sp), zero_sp) |>
    dplyr::arrange(dplyr::desc(n_locs_global))

  ts(sprintf("  Species with >= 1 global location: %d",
             sum(global_locs_sp_full$n_locs_global >= 1L)))
  ts(sprintf("  Species with >= 30 global locations: %d",
             sum(global_locs_sp_full$n_locs_global >= 30L)))

  readr::write_csv(global_locs_sp_full, out_global_gf_sp)
  ts(sprintf("  Saved -> wallacean/wallacean_global_gf_locs_per_species.csv (%d rows)",
             nrow(global_locs_sp_full)))

} else {
  ts("Step 5: Global GF location files already exist — skipping.")
  n_ge30_sh <- sum(
    readr::read_csv(out_global_gf_sh, show_col_types = FALSE)$n_locs_global >= 30L
  )
  n_ge30_sp <- sum(
    readr::read_csv(out_global_gf_sp, show_col_types = FALSE)$n_locs_global >= 30L
  )
  ts(sprintf("  SH codes with >= 30 global GF locations: %d", n_ge30_sh))
  ts(sprintf("  Species with >= 30 global GF locations: %d", n_ge30_sp))
}

# =============================================================================
# Figure 2 — grid-cell occupancy histograms (GlobalFungi only)
# =============================================================================
# Four panels: named species and SH codes, each at Canada and global scope.
# Occupancy = number of distinct 30 arc-second grid cells (~1 km^2) a taxon
# occupies. Requires the global occupancy CSVs written in Step 5 above.
# =============================================================================
if (!file.exists(paths$fig_wallacean_occ) &&
    file.exists(out_global_gf_sp) && file.exists(out_global_gf_sh)) {

  ts("Figure 2: building occupancy histograms...")
  library(patchwork)   # provides the | and / plot-composition operators

  # Canada-scope occupancy: distinct 30 arc-second cells per taxon
  gf_can <- emf |>
    dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE,
                  !is.na(lat), !is.na(lon)) |>
    dplyr::mutate(grid_lat = snap_30arcsec(lat), grid_lon = snap_30arcsec(lon))

  occ_sh_canada <- gf_can |>
    dplyr::distinct(sh_code, grid_lat, grid_lon) |>
    dplyr::count(sh_code, name = "n_canada")
  occ_sp_canada <- gf_can |>
    dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
    dplyr::distinct(species, grid_lat, grid_lon) |>
    dplyr::count(species, name = "n_canada")

  # Global-scope occupancy: pre-computed per taxon in Step 5
  global_sp <- readr::read_csv(out_global_gf_sp, show_col_types = FALSE)
  global_sh <- readr::read_csv(out_global_gf_sh, show_col_types = FALSE)

  occ_sp <- occ_sp_canada |>
    dplyr::left_join(global_sp, by = "species") |>
    dplyr::mutate(n_global = tidyr::replace_na(n_locs_global, 0L))
  occ_sh <- occ_sh_canada |>
    dplyr::left_join(global_sh, by = "sh_code") |>
    dplyr::mutate(n_global = tidyr::replace_na(n_locs_global, 0L))

  make_occ_hist <- function(values, xlab, ylab, title) {
    ggplot2::ggplot(tibble::tibble(v = values), ggplot2::aes(x = v)) +
      ggplot2::geom_histogram(bins = 40, fill = "#4393C3", colour = "white",
                              linewidth = 0.2) +
      ggplot2::labs(x = xlab, y = ylab, title = title) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                     panel.grid.minor = ggplot2::element_blank())
  }
  xlab_can  <- "No. 30 arc-second grid cells occupied in Canada"
  xlab_glob <- "No. 30 arc-second grid cells occupied globally"

  fig2 <- (make_occ_hist(occ_sp$n_canada, xlab_can,  "Number of species",  "(a) Named species — Canada") |
           make_occ_hist(occ_sh$n_canada, xlab_can,  "Number of SH codes", "(b) SH codes — Canada")) /
          (make_occ_hist(occ_sp$n_global, xlab_glob, "Number of species",  "(c) Named species — global") |
           make_occ_hist(occ_sh$n_global, xlab_glob, "Number of SH codes", "(d) SH codes — global"))

  ggplot2::ggsave(paths$fig_wallacean_occ, fig2, width = 11, height = 8, dpi = 300)
  ts(sprintf("  Saved %s", basename(paths$fig_wallacean_occ)))
} else {
  ts("Figure 2: skipping (already exists, or global occupancy CSVs not yet built).")
}

ts("11_wallacean.R complete.")
