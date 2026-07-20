# =============================================================================
# Hutchinsonian Shortfall: Ecozone and Climate-Space Sampling Coverage
# =============================================================================
# PURPOSE
#   How completely has Canada's ecological space been sampled for
#   ectomycorrhizal (EcM) fungi? Two complementary partitions of that space are
#   used, and both ask the same question of it:
#     Steps 1-5  Ecozones  — Canada's 15 named terrestrial ecozones: which are
#                sampled at all, which meet the minimum-sample thresholds used
#                elsewhere in the SI (>=10, >=30), and at what density?
#     Steps 6-9  Climate   — Canada's two-dimensional temperature x
#                precipitation space: what proportion of the places sharing
#                each climate have been sampled, and where in geographic space
#                do the well- and poorly-covered climates fall?
#
#   Sampling points use both GenBank and GlobalFungi coordinates, restricted to
#   records whose coordinates were validated inside the GADM Canada boundary at
#   dataset assembly (coord_in_canada == TRUE, 04_combine_ecm_dataset.R).
#
# CLIMATE METHOD (Steps 6-9)
#   The representation is adapted from Figure 6 in the following preprint:
# Blitz the Gap: a nation-wide effort to guide citizen science toward the needs
# of biodiversity science. K. Hebert et al.
# https://doi.org/10.32942/X2T09G
#
#   1. Define Canada's climate space from two WorldClim v2.1 layers:
#        BIO1  = mean annual temperature (MAT, deg C)
#        BIO12 = mean annual precipitation (MAP, mm)
#      The native 30 arc-second raster is AGGREGATED to ~10 arc-minutes
#      (AGG_FACTOR = 20) so the spatial grain matches the colleagues' workflow.
#   2. Bin the two-dimensional climate space into an N_CLIMATE_BINS x
#      N_CLIMATE_BINS grid (equal-width bins per axis). Each occupied bin is a
#      "climate zone". count_canada = number of grid cells per zone.
#   3. A grid cell is "sampled" if >= 1 EcM sample falls in it.
#      coverage(zone) = (sampled cells in zone) / (Canada cells in zone),
#      i.e. the proportion of the places sharing a climate that have been
#      sampled at all. coverage lies in [0, 1]; 0 = climate present but never
#      sampled.
#   4. Produce a 6-panel figure laid out as 2 rows x 3 columns:
#        Columns: Available climate frequency | GlobalFungi | GF + GenBank
#        Row 1 (a, b, c): geographic maps
#        Row 2 (d, e, f): the same quantity in climate space
#      The "Available climate" column shows climate frequency (count_canada =
#      Canadian cells per zone, i.e. how common each climate is); this is the
#      denominator of coverage and is analogous to the colleagues' Fig S5 a/c.
#      The GlobalFungi and GF + GenBank columns show sampling coverage.
#      In panel d, zones never sampled by GlobalFungi + GenBank combined (the
#      grey, coverage == 0 zones in panel f) are additionally outlined in
#      white, so the unsampled climates are visible directly in panel d
#      without needing to cross-reference panel f.
#
# CLIMATE SCOPES (Steps 6-9)
#   - "GlobalFungi": GF samples with >= 1 EcM SH code (named or not). Because the
#     downstream dataset (emf) is already EcM-only, every GF sample qualifies.
#   - "GlobalFungi + GenBank": the GF samples above plus GenBank EcM records.
#   In both scopes the sampling unit for coverage is the climate grid cell
#   (a cell counts once regardless of how many samples fall in it).
#
# INPUTS
#   emf                                                    (auto-loaded by 00_setup.R)
#   data_derived/spatial/canada_simple.gpkg                (paths$canada_bound)
#   data_derived/spatial/ecoregions_processed.gpkg         (paths$ecoregions_processed)
#   data_derived/spatial/lakes_canada_albers.gpkg          (paths$lakes_albers)
#   data_derived/spatial/ecozone_names.csv                 (paths$ecozone_names)
#   data_raw/climate/wc2.1_country/CAN_wc2.1_30s_bio.tif   (paths$climate_raster)
#
# OUTPUTS
#   data_derived/hutchinsonian/hutchinsonian_ecozone_summary.csv       — ecozone
#       sampling-threshold counts (Table S16b)
#   data_derived/hutchinsonian/hutchinsonian_ecozone_sample_counts.csv — per-ecozone
#       counts, areas and sampling densities (Table S1)
#   figures/Figure-S4_ecozone_sampling_map.png      (paths$fig_ecozone_sampling)      -- white bg, used in SI
#   figures/Figure-S4_ecozone_sampling_map_grey.png (paths$fig_ecozone_sampling_grey) -- #F2F2F2 bg, Figure 5 panel source
#   figures/Figure-03_climate_gap.png               (paths$fig_climate_gap)           -- white bg, used in manuscript
#   figures/Figure-03_climate_gap_grey.png          (paths$fig_climate_gap_grey)      -- #F2F2F2 bg, Figure 5 panel source
#
# NOTES
#   - Deterministic: no random sampling, so no seed is required.
#   - Sentinel-file guard: the climate-gap figure (Steps 6-9) is rebuilt only if
#     BOTH versions do not yet exist. Delete one or both to force regeneration.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

sf::sf_use_s2(FALSE)

# ---- Step 1: Prepare sampling points ----------------------------------------
# coord_in_canada == TRUE guarantees the coordinates were validated against the
# GADM Canada boundary when the dataset was assembled (04_combine_ecm_dataset.R).

# GenBank records with validated Canadian coordinates, reduced to unique
# locations. The GlobalFungi locations are added alongside these in Step 2.
gb_parsed <- emf |>
  dplyr::filter(source == "GenBank", coord_in_canada == TRUE) |>
  dplyr::distinct(lat, lon)

# Canada boundary (used by the map steps below)
canada_bound <- sf::st_read(paths$canada_bound, quiet = TRUE)

# ---- Step 2: Ecozone setup and sampling-locations map (Figure S4) ------------

canada_albers <- sf::st_transform(canada_bound, crs_albers)

# -- Shared ecozone / ecoregion setup (used by sampling map and ecozone climate space plot) ----

ecoregions_raw <- sf::st_read(paths$ecoregions_processed, quiet = TRUE) |>
  sf::st_transform(crs_albers) |>
  sf::st_make_valid()

ecozone_names_tbl <- readr::read_csv(paths$ecozone_names, show_col_types = FALSE)

ecoregions_named <- dplyr::left_join(
  ecoregions_raw, ecozone_names_tbl, by = "ECOZONE"
) |>
  dplyr::mutate(
    NAME_EN = dplyr::if_else(is.na(NAME_EN),
                             paste("Ecozone", ECOZONE),
                             NAME_EN)
  )

# Clip ecoregions to Canada boundary
# st_buffer(0) forces topology rebuild and resolves side-location conflicts
# that st_make_valid() alone cannot fix
ecoregions_clipped <- sf::st_intersection(
  sf::st_buffer(sf::st_make_valid(ecoregions_named), 0),
  sf::st_buffer(sf::st_make_valid(canada_albers),    0)
) |> sf::st_make_valid()

# Official Environment Canada ecozone colours, extracted from the vector fills
# of the published NRCan/EC terrestrial ecozones map PDF.
# Names ordered alphabetically to match sort(unique(...)).
unique_ecozones <- sort(unique(ecoregions_clipped$NAME_EN))
n_ecozones      <- length(unique_ecozones)
ecozone_colors  <- c(
  "Arctic Cordillera"  = "#E9F6FD",
  "Atlantic Maritime"  = "#B2DFDB",
  "Boreal Cordillera"  = "#8BBBE4",
  "Boreal Plains"      = "#EDF2C3",
  "Boreal Shield"      = "#99CD75",
  "Hudson Plains"      = "#F1EE8B",
  "Mixedwood Plains"   = "#FFE8A2",
  "Montane Cordillera" = "#CBE19A",
  "Northern Arctic"    = "#C1E8FA",
  "Pacific Maritime"   = "#44B763",
  "Prairies"           = "#FFF799",
  "Southern Arctic"    = "#6DCFF6",
  "Taiga Cordillera"   = "#E7D7EA",
  "Taiga Plains"       = "#ECEDEE",
  "Taiga Shield"       = "#BDC0C2"
)

# Major lakes
lakes_sf <- tryCatch(
  sf::st_read(paths$lakes_albers, quiet = TRUE) |>
    sf::st_transform(crs_albers) |>
    sf::st_make_valid() |>
    sf::st_buffer(0) |>
    sf::st_intersection(sf::st_buffer(sf::st_make_valid(canada_albers), 0)),
  error = function(e) { warning("Lakes file not found — skipping.", call. = FALSE); NULL }
)

# Sampling locations with source label (used by both sampling map and ecozone climate space plot).
# GenBank is bound first so it is drawn underneath; GlobalFungi is bound second so it
# is drawn on top.
# Each source is deduplicated within itself (gb_parsed already is; GF deduped below);
# cross-source deduplication is intentionally omitted so that co-located points from
# both sources are retained, with the GF point rendered on top.
locs_src_sf <- dplyr::bind_rows(
  dplyr::select(gb_parsed, lat, lon) |>
    dplyr::mutate(source = "GenBank"),
  dplyr::filter(emf, source == "GlobalFungi", coord_in_canada == TRUE) |>
    dplyr::distinct(lat, lon) |>
    dplyr::mutate(source = "GlobalFungi")
) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  sf::st_transform(crs_albers)

# -- Sampling locations map (ecozone style) ------------------------------------

# Pad the plotted extent beyond the Canada boundary bbox so the outer map
# border doesn't sit flush against the figure edge. Buffer is a fraction of
# each dimension's range (5%), applied symmetrically on both sides.
sampling_map_bbox   <- sf::st_bbox(canada_albers)
sampling_map_buffer <- 0.05
sampling_map_xrange  <- sampling_map_bbox[3] - sampling_map_bbox[1]
sampling_map_yrange  <- sampling_map_bbox[4] - sampling_map_bbox[2]
sampling_map_xlim <- c(sampling_map_bbox[1] - sampling_map_buffer * sampling_map_xrange,
                       sampling_map_bbox[3] + sampling_map_buffer * sampling_map_xrange)
sampling_map_ylim <- c(sampling_map_bbox[2] - sampling_map_buffer * sampling_map_yrange,
                       sampling_map_bbox[4] + sampling_map_buffer * sampling_map_yrange)

# build_sampling_map(bg) returns the map with panel/plot background colour
# `bg`; called once for the white (SI) version and once for the grey
# (#F2F2F2, Figure 5 schematic panel) version. The sampling-point marker fill
# ("white", inside geom_sf/guides below) is a data-encoding colour, not the
# plot background, and is intentionally left unchanged in both versions.
build_sampling_map <- function(bg) {
  ggplot2::ggplot() +
    # Ecoregions: rainbow fill by ecozone, grey outlines show ecoregion boundaries
    ggplot2::geom_sf(data = ecoregions_clipped,
                     ggplot2::aes(fill = NAME_EN),
                    colour = "grey30", linewidth = 0, alpha = 0.8) +
    # Major lakes
    (if (!is.null(lakes_sf))
      ggplot2::geom_sf(data = lakes_sf,
                       fill = "dodgerblue", colour = "darkblue",
                       linewidth = 0.2, alpha = 0.9)
    else NULL) +
    # Canada outer boundary
    ggplot2::geom_sf(data = canada_albers,
                     fill = NA, colour = "black", linewidth = 0.2) +
    # Sampling locations: circle (GlobalFungi) or triangle (GenBank),
    # white fill, colour-blind-safe outline
    ggplot2::geom_sf(data = locs_src_sf,
                     ggplot2::aes(colour = source, shape = source),
                     fill = "white", size = 2.2, stroke = 1.1, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = ecozone_colors, name = "Ecozones") +
    ggplot2::scale_colour_manual(
      breaks = c("GlobalFungi", "GenBank"),
      values = c("GlobalFungi" = "#648FFF", "GenBank" = "#E69F00"),
      name   = "Sample source"
    ) +
    ggplot2::scale_shape_manual(
      breaks = c("GlobalFungi", "GenBank"),
      values = c("GlobalFungi" = 21, "GenBank" = 24),
      name   = "Sample source"
    ) +
    ggplot2::guides(
      fill   = ggplot2::guide_legend(order = 1,
                                     override.aes = list(alpha = 0.8, colour = NA)),
      colour = ggplot2::guide_legend(order = 2,
                                     override.aes = list(shape = c(21L, 24L),
                                                         fill  = "white",
                                                         size  = 4)),
      shape  = "none"
    ) +
    ggplot2::coord_sf(
      xlim   = sampling_map_xlim,
      ylim   = sampling_map_ylim,
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 7),
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),
      plot.background  = ggplot2::element_rect(fill = bg, colour = NA)
    )
}

save_fig_formats(paths$fig_ecozone_sampling, build_sampling_map("white"),
                 width = 12, height = 9, dpi = 300, bg = "white")

# Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
# schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
ggplot2::ggsave(paths$fig_ecozone_sampling_grey, build_sampling_map(fig5_grey_bg),
                width = 12, height = 9, dpi = 300, bg = fig5_grey_bg)

# ---- Step 3: Ecozone-level sampling summary (Table S16b) -------------------
# Companion to the ecoregion-level summary (Table S16): of Canada's
# n_ecozones named ecozones (from ecoregions_clipped, Step 2), how many have
# *any* EcM fungal sampling, and how many meet the minimum-sample thresholds
# already used elsewhere in the SI (>=10 sites: Figure S8 coverage plot;
# >=30 sites: Figure S9 accumulation curves)?
#
# Three sample definitions are reported side by side, because "sample" is
# ambiguous across the SI:
#   1. GlobalFungi + GenBank combined, raw unique sampling locations
#      (locs_src_sf, Step 2 -- same data underlying Table S16 and the
#      sampling map).
#   2. GlobalFungi only, project-standard 3-decimal-binned "sites"
#      (add_site_id(), 00_setup.R) -- same scope/unit as the >=10 / >=30
#      thresholds already used for Figures S8 and S9.
#   3. GlobalFungi + GenBank combined, 3-decimal-binned "sites" -- same
#      combined scope as (1), but using the project's canonical site unit
#      instead of raw coordinates.

# Helper: assign a set of points (sf, any CRS) to an ecozone and count points
# per ecozone.
#
# IMPORTANT: the join is against the detailed, UNCLIPPED ecoregion layer
# (`ecoregions_named`), NOT `ecoregions_clipped`. `ecoregions_clipped` is the
# ecoregions re-intersected with `canada_simple`, which is built with
# simplifyGeom(tolerance = 0.1 deg). That leaves the 49th-parallel border as a
# few sparse vertices with very long straight edges; a segment straight in
# lon/lat bows up to ~46 km northward when reprojected to Albers, so the clip
# trims a wide border strip and silently drops genuine near-border points (all
# GenBank) from the join. Joining against the unclipped ecoregions removes that
# discrepancy and matches how 09_linnean.R assigns ecozones; `ecoregions_clipped`
# is retained only for drawing the map fills.
#
# `snap_tol_m`: a point falling just offshore of the ecoregion coastline (inside
# Canada only because the coarse boundary bulges over water) matches no polygon
# under st_within. Such points are snapped to the nearest ecozone, but only
# within snap_tol_m, so this fallback can never silently absorb a genuinely
# mislocated coordinate. Distances are in metres because ecoregions_named is in
# the (metre-based) Albers projection.
count_points_per_ecozone <- function(pts_sf, snap_tol_m = 5000) {
  pts_sf <- sf::st_transform(pts_sf, sf::st_crs(ecoregions_named))
  joined <- sf::st_join(pts_sf, dplyr::select(ecoregions_named, NAME_EN),
                        join = sf::st_within)
  # Bounded nearest-ecozone fallback for points inside no polygon.
  miss <- which(is.na(joined$NAME_EN))
  if (length(miss) > 0L) {
    near_idx  <- sf::st_nearest_feature(pts_sf[miss, ], ecoregions_named)
    near_dist <- as.numeric(sf::st_distance(pts_sf[miss, ],
                                            ecoregions_named[near_idx, ],
                                            by_element = TRUE))
    ok <- near_dist <= snap_tol_m
    joined$NAME_EN[miss[ok]] <- ecoregions_named$NAME_EN[near_idx[ok]]
  }
  joined <- sf::st_drop_geometry(joined)
  dplyr::count(dplyr::filter(joined, !is.na(NAME_EN)), NAME_EN, name = "n")
}

ez_threshold_row <- function(counts, def_label) {
  tibble::tibble(
    definition    = def_label,
    n_ecozones    = n_ecozones,
    ecozones_ge1  = sum(counts$n >= 1L),
    ecozones_ge10 = sum(counts$n >= 10L),
    ecozones_ge30 = sum(counts$n >= 30L)
  )
}

# Definition 1: GF + GenBank combined, raw unique sampling locations
counts_combined_raw <- count_points_per_ecozone(locs_src_sf)

# Definition 2: GlobalFungi only, 3-decimal-binned sites
sites_gf_only <- emf |>
  dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE) |>
  add_site_id() |>
  dplyr::distinct(site, site_lat, site_lon)
pts_sites_gf_only <- sf::st_as_sf(sites_gf_only,
                                  coords = c("site_lon", "site_lat"), crs = 4326)
counts_gf_sites <- count_points_per_ecozone(pts_sites_gf_only)

# Definition 3: GF + GenBank combined, 3-decimal-binned sites
sites_combined <- emf |>
  dplyr::filter(coord_in_canada == TRUE) |>
  add_site_id() |>
  dplyr::distinct(site, site_lat, site_lon)
pts_sites_combined <- sf::st_as_sf(sites_combined,
                                   coords = c("site_lon", "site_lat"), crs = 4326)
counts_combined_sites <- count_points_per_ecozone(pts_sites_combined)

ecozone_threshold_summary <- dplyr::bind_rows(
  ez_threshold_row(counts_combined_raw,   "GF + GenBank combined, raw unique locations"),
  ez_threshold_row(counts_gf_sites,       "GlobalFungi only, 3-decimal-binned sites"),
  ez_threshold_row(counts_combined_sites, "GF + GenBank combined, 3-decimal-binned sites")
)

for (i in seq_len(nrow(ecozone_threshold_summary))) {
  r <- ecozone_threshold_summary[i, ]
}

readr::write_csv(ecozone_threshold_summary,
                 file.path(paths$out_hutchinsonian, "hutchinsonian_ecozone_summary.csv"))

# ---- Step 4: Per-ecozone sample counts by source (Table S1) -----------------
# One row per named ecozone with raw unique locations and 3-decimal-binned sites,
# split by source. "Total" columns sum the two source columns (a location/site
# shared by both sources counts once per source). Same within-polygon +
# nearest-snap join (count_points_per_ecozone()) as Step 3.

counts_gf_locs_ez <- count_points_per_ecozone(dplyr::filter(locs_src_sf, source == "GlobalFungi"))
counts_gb_locs_ez <- count_points_per_ecozone(dplyr::filter(locs_src_sf, source == "GenBank"))
counts_gf_sites_ez <- count_points_per_ecozone(pts_sites_gf_only)

sites_gb_only <- emf |>
  dplyr::filter(source == "GenBank", coord_in_canada == TRUE) |>
  add_site_id() |>
  dplyr::distinct(site, site_lat, site_lon)
pts_sites_gb_only <- sf::st_as_sf(sites_gb_only,
                                  coords = c("site_lon", "site_lat"), crs = 4326)
counts_gb_sites_ez <- count_points_per_ecozone(pts_sites_gb_only)

ecozone_sample_counts <- tibble::tibble(ecozone = unique_ecozones) |>
  dplyr::left_join(dplyr::rename(counts_gf_locs_ez,  gf_locations = n), by = c("ecozone" = "NAME_EN")) |>
  dplyr::left_join(dplyr::rename(counts_gb_locs_ez,  gb_locations = n), by = c("ecozone" = "NAME_EN")) |>
  dplyr::left_join(dplyr::rename(counts_gf_sites_ez, gf_sites = n),     by = c("ecozone" = "NAME_EN")) |>
  dplyr::left_join(dplyr::rename(counts_gb_sites_ez, gb_sites = n),     by = c("ecozone" = "NAME_EN")) |>
  dplyr::mutate(dplyr::across(c(gf_locations, gb_locations, gf_sites, gb_sites),
                              ~ tidyr::replace_na(as.integer(.), 0L))) |>
  dplyr::mutate(total_locations = gf_locations + gb_locations,
                total_sites     = gf_sites + gb_sites) |>
  dplyr::arrange(ecozone)

# ---- Step 5: Ecozone areas -> sampling density (Table S1) -------------------
# Area of each named ecozone within Canada, from the SAME unclipped, Albers
# ecoregion polygons used to assign the sample counts (so counts and areas share
# one polygon basis). Ecozone area = sum of its constituent ecoregion areas.

ecozone_areas <- ecoregions_named |>
  sf::st_drop_geometry() |>
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(ecoregions_named)) / 1e6) |>
  dplyr::group_by(ecozone = NAME_EN) |>
  dplyr::summarise(area_km2 = sum(area_km2), .groups = "drop")

ecozone_sample_counts <- ecozone_sample_counts |>
  dplyr::left_join(ecozone_areas, by = "ecozone") |>
  dplyr::mutate(
    density_total = total_locations / area_km2 * 10000,   # locations per 10,000 km^2
    density_gf    = gf_locations   / area_km2 * 10000
  )

readr::write_csv(ecozone_sample_counts,
                 file.path(paths$out_hutchinsonian, "hutchinsonian_ecozone_sample_counts.csv"))

# =============================================================================
# Climate-space sampling coverage (Steps 6-9)
# =============================================================================
# Steps 1-5 above measure sampling coverage across Canada's ecozones. Steps 6-9
# ask the same question of Canada's *climate* space, and show where in
# geographic space the well- and poorly-covered climates fall. See the METHOD
# and SCOPES notes in the header block at the top of this script.
# =============================================================================

# ---- Tunable parameters for the climate-space analysis ----------------------
AGG_FACTOR      <- 20L   # 30 arc-sec -> ~10 arc-min (20 * 30" = 600" = 10')
N_CLIMATE_BINS  <- 50L   # bins per climate axis (50 x 50 grid of climate zones)

# Colour palettes match the colleagues' figures (via the colorspace package):
#   - "Batlow" sequential  -> available climate frequency (their Figs S4/S5)
#   - "Temps"  divergingx   -> sampling coverage, reversed so poor = red and
#                              full = teal (their Fig 6)
FREQ_PALETTE         <- "Batlow"
COVERAGE_PALETTE     <- "Temps"
NEVER_SAMPLED_COLOUR <- "grey85"   # climate present in Canada but never sampled

fig_out      <- paths$fig_climate_gap
fig_out_grey <- paths$fig_climate_gap_grey

# Sentinel guard: skip Steps 6-9 entirely only if BOTH versions of the
# climate-gap figure already exist. Delete one or both to force regeneration.
if (!(file.exists(fig_out) && file.exists(fig_out_grey))) {

  # ---- Step 6: Build the aggregated, Canada-masked climate grid ------------

  if (!file.exists(paths$climate_raster))
    stop("Climate raster not found: ", paths$climate_raster, call. = FALSE)

  clim_full <- terra::rast(paths$climate_raster)

  # WorldClim layer naming varies by download method; match BIO1 / BIO12 by the
  # numeric suffix rather than assuming a fixed band order.
  layer_names <- names(clim_full)
  mat_idx <- grep("bio_?1$",  layer_names)[1]
  map_idx <- grep("bio_?12$", layer_names)[1]
  if (is.na(mat_idx) || is.na(map_idx))
    stop("Could not locate BIO1 and BIO12 layers in: ",
         paste(layer_names, collapse = ", "), call. = FALSE)

  clim <- clim_full[[c(mat_idx, map_idx)]]
  names(clim) <- c("MAT", "MAP")

  # Some WorldClim distributions store MAT x10 as integers; convert if detected.
  mat_rng <- terra::values(clim[["MAT"]])
  if (!all(is.na(mat_rng)) && max(abs(mat_rng), na.rm = TRUE) > 100) {
    clim[["MAT"]] <- clim[["MAT"]] / 10
  }

  # Aggregate 30" -> ~10' (mean of contributing fine cells), then mask to the
  # project's canonical Canada boundary so the climate cloud and the coverage
  # denominator use Canadian land cells only.
  clim <- terra::aggregate(clim, fact = AGG_FACTOR, fun = "mean", na.rm = TRUE)

  # This step rebuilds the Canada polygon rather than reusing `canada_bound`
  # from Step 1: it needs an explicitly WGS84, st_make_valid()-repaired
  # version for terra::mask() and for the neatline difference below, which is
  # not how `canada_bound` / `canada_albers` are constructed.
  canada_wgs84 <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
    sf::st_transform(4326) |>
    sf::st_make_valid()
  canada_vect <- terra::vect(canada_wgs84)
  clim <- terra::mask(terra::crop(clim, canada_vect), canada_vect)

  # Per-cell climate values, aligned to cell indices (NA outside Canada).
  cell_vals <- terra::values(clim)            # matrix: columns MAT, MAP
  cell_ok   <- !is.na(cell_vals[, "MAT"]) & !is.na(cell_vals[, "MAP"])

  # ---- Step 7: Bin climate space into an N x N grid of "climate zones" -----

  mat_ok <- cell_vals[cell_ok, "MAT"]
  map_ok <- cell_vals[cell_ok, "MAP"]

  # Equal-width breaks spanning the Canadian climate range on each axis.
  mat_breaks <- seq(min(mat_ok), max(mat_ok), length.out = N_CLIMATE_BINS + 1)
  map_breaks <- seq(min(map_ok), max(map_ok), length.out = N_CLIMATE_BINS + 1)
  mat_width  <- diff(mat_breaks)[1]
  map_width  <- diff(map_breaks)[1]

  # Assign each Canadian cell to a bin on each axis. all.inside = TRUE folds the
  # uppermost edge value into the last bin (1..N_CLIMATE_BINS).
  cell_xbin <- rep(NA_integer_, length(cell_ok))
  cell_ybin <- rep(NA_integer_, length(cell_ok))
  cell_xbin[cell_ok] <- findInterval(cell_vals[cell_ok, "MAT"], mat_breaks,
                                     rightmost.closed = TRUE, all.inside = TRUE)
  cell_ybin[cell_ok] <- findInterval(cell_vals[cell_ok, "MAP"], map_breaks,
                                     rightmost.closed = TRUE, all.inside = TRUE)
  cell_bin <- ifelse(cell_ok, paste(cell_xbin, cell_ybin, sep = "_"), NA)

  # count_canada: number of Canadian cells per climate zone.
  zone_canada <- tibble::as_tibble(table(bin = cell_bin[cell_ok]),
                                   .name_repair = "minimal")
  names(zone_canada) <- c("bin", "count_canada")
  zone_canada$count_canada <- as.integer(zone_canada$count_canada)

  # ---- Step 8: Coverage per scope ------------------------------------------
  # For a set of sample coordinates, find the unique Canadian climate cells they
  # fall in, tally sampled cells per zone, and divide by count_canada.
  # Returns a per-zone tibble joined to zone geometry (bin centres) plus a
  # coverage SpatRaster for the geographic map.
  compute_coverage <- function(pts_df, scope_label) {

    xy    <- as.matrix(pts_df[, c("lon", "lat")])
    cells <- terra::cellFromXY(clim, xy)
    cells <- unique(cells[!is.na(cells)])
    cells <- cells[cell_ok[cells]]          # keep only cells with climate data

    samp_bins <- cell_bin[cells]
    zone_samp <- tibble::as_tibble(table(bin = samp_bins), .name_repair = "minimal")
    names(zone_samp) <- c("bin", "count_sampled")
    zone_samp$count_sampled <- as.integer(zone_samp$count_sampled)

    # Per-zone coverage table (all occupied zones; unsampled -> 0).
    zones <- zone_canada |>
      dplyr::left_join(zone_samp, by = "bin") |>
      dplyr::mutate(
        count_sampled = dplyr::coalesce(count_sampled, 0L),
        coverage      = count_sampled / count_canada,
        xbin = as.integer(sub("_.*$", "", bin)),
        ybin = as.integer(sub("^.*_", "", bin)),
        # Bin centres in real climate units (for geom_tile placement).
        mat_centre = mat_breaks[xbin] + mat_width / 2,
        map_centre = map_breaks[ybin] + map_width / 2,
        scope      = scope_label
      )

    # Coverage raster for the map: each Canadian cell takes its zone's coverage.
    cov_lookup <- stats::setNames(zones$coverage, zones$bin)
    cov_vec    <- unname(cov_lookup[cell_bin])   # NA outside Canada
    cov_rast   <- clim[["MAT"]]
    terra::values(cov_rast) <- cov_vec
    names(cov_rast) <- "coverage"
    # Cells whose zone is never sampled (coverage 0) -> NA so the grey Canada
    # base layer shows through on the map.
    cov_rast <- terra::ifel(cov_rast == 0, NA, cov_rast)

    list(zones = zones, cov_rast = cov_rast)
  }

  # GlobalFungi: every EcM-only GF sample qualifies (>= 1 EcM SH by construction).
  gf_pts <- emf |>
    dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE) |>
    dplyr::distinct(lat, lon)

  # GenBank EcM records with validated Canadian coordinates. Constructed the
  # same way as `gb_parsed` in Step 1, but kept separate so this climate
  # analysis reads as a self-contained unit alongside gf_pts / comb_pts.
  gb_pts <- emf |>
    dplyr::filter(source == "GenBank", coord_in_canada == TRUE) |>
    dplyr::distinct(lat, lon)

  # Combined scope: union of distinct coordinates.
  comb_pts <- dplyr::bind_rows(gf_pts, gb_pts) |>
    dplyr::distinct(lat, lon)

  cov_gf   <- compute_coverage(gf_pts,   "GlobalFungi")
  cov_comb <- compute_coverage(comb_pts, "GlobalFungi + GenBank")

  # ---- Available climate frequency (the "denominator" of coverage) ----------
  # count_canada = number of Canadian grid cells per climate zone, i.e. how
  # COMMON each climate is across Canada. This is the "available climate space"
  # shown in the colleagues' Figures S4/S5 (a = climate space, c = geographic).
  zone_freq <- zone_canada |>
    dplyr::mutate(
      xbin = as.integer(sub("_.*$", "", bin)),
      ybin = as.integer(sub("^.*_", "", bin)),
      mat_centre = mat_breaks[xbin] + mat_width / 2,
      map_centre = map_breaks[ybin] + map_width / 2
    ) |>
    # Flag zones never sampled by GlobalFungi + GenBank combined (i.e. the
    # zones shown grey in panel f) so panel d can outline them in white.
    dplyr::left_join(dplyr::select(cov_comb$zones, bin, coverage), by = "bin") |>
    dplyr::mutate(never_sampled_comb = coverage == 0)

  # Climate-frequency raster: each Canadian cell takes its zone's count_canada.
  freq_lookup <- stats::setNames(zone_freq$count_canada, zone_freq$bin)
  freq_rast   <- clim[["MAT"]]
  terra::values(freq_rast) <- unname(freq_lookup[cell_bin])
  names(freq_rast) <- "frequency"

  # ---- Step 9: Build the six panels (2 rows x 3 columns) -------------------
  #   Columns: Available climate frequency | GlobalFungi | GF + GenBank
  #   Rows:    geographic maps (top)       | climate space (bottom)

  canada_albers_clim <- sf::st_transform(canada_wgs84, crs_albers)

  # Tight Canada bounding box, used by coord_sf so each map fills its panel.
  bbox_can <- sf::st_bbox(canada_albers_clim)
  xlim_can <- c(bbox_can[["xmin"]], bbox_can[["xmax"]])
  ylim_can <- c(bbox_can[["ymin"]], bbox_can[["ymax"]])

  # Neatline / "outside-Canada" mask: bounding box MINUS Canada. Drawn (white)
  # on top of the raster but under the boundary line, it hides the part of any
  # ~10' grid cell that straddles the coastline/border. This clips the *display*
  # to Canada WITHOUT discarding any data -- every Canadian cell is still shown
  # in full; only the over-the-border fringe of edge cells is covered.
  canada_exterior <- sf::st_difference(
    sf::st_as_sfc(bbox_can),
    # st_buffer(0) after st_make_valid rebuilds topology and avoids the
    # "side location conflict" GEOS throws on the GADM-derived polygon.
    sf::st_union(sf::st_buffer(sf::st_make_valid(canada_albers_clim), 0))
  )

  # -- Two fill scales (colorspace palettes, matching the colleagues) ---------
  # Compact horizontal colourbars: each is drawn ONCE, inside a climate-space
  # panel (frequency -> panel d, coverage -> panel e); the other panels that
  # share the scale suppress it with guides(fill = "none").
  cbar <- ggplot2::guide_colourbar(barwidth  = grid::unit(3.0, "cm"),
                                   barheight = grid::unit(0.3, "cm"),
                                   title.position = "top")

  # Coverage scale: "Temps" divergingx, reversed (poor = red, full = teal), as
  # in their Fig 6. Fixed 0-1 limits so the GF and GF + GenBank panels are
  # directly comparable; mid = 0.5 centres the light tone of the ramp.
  cov_fill_scale <- colorspace::scale_fill_continuous_divergingx(
    palette  = COVERAGE_PALETTE, rev = TRUE, mid = 0.5,
    limits   = c(0, 1),
    breaks   = c(0, 0.25, 0.5, 0.75, 1),
    labels   = scales::percent_format(accuracy = 1),
    na.value = "transparent",
    name     = "Climate coverage (% of zone's cells sampled)",
    guide    = cbar
  )

  # Frequency scale: "Batlow" sequential (rare = dark, common = light), as in
  # their Figs S4/S5. Square-root transformed because count_canada is strongly
  # right-skewed (median 10, max ~1500).
  freq_fill_scale <- colorspace::scale_fill_continuous_sequential(
    palette  = FREQ_PALETTE, rev = FALSE, trans = "sqrt",
    na.value = "transparent",
    name     = "Climate frequency (Canadian cells per zone)",
    breaks   = c(1, 100, 400, 900, 1400),
    guide    = cbar
  )

  # Theme placing a legend inside the upper-left (empty cold/wet corner) of a
  # climate-space panel, on a semi-transparent background so tiles show through.
  legend_inside <- ggplot2::theme(
    legend.position        = "inside",
    legend.position.inside = c(0.02, 0.90),
    legend.justification.inside = c(0, 1),
    legend.direction       = "horizontal",
    legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.65),
                                              colour = NA),
    legend.title = ggplot2::element_text(size = 8.5),
    legend.text  = ggplot2::element_text(size = 7.5)
  )

  # Shared theming helpers (keep the six panels visually consistent).
  # Both take a `bg` argument (panel/plot background colour) so the whole
  # composite can be rendered twice -- once white (manuscript), once
  # #F2F2F2 (Figure 5 schematic source panel) -- from the same panel code.
  clim_theme_fn <- function(bg) {
    ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = bg, colour = NA),
        plot.background  = ggplot2::element_rect(fill = bg, colour = NA)
      )
  }
  map_theme_fn <- function(bg) {
    ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.margin      = ggplot2::margin(2, 2, 2, 2),
        panel.background = ggplot2::element_rect(fill = bg, colour = NA),
        plot.background  = ggplot2::element_rect(fill = bg, colour = NA)
      )
  }
  clim_labs <- ggplot2::labs(x = "Mean annual temperature (\u00b0C)",
                             y = "Mean annual precipitation (mm)")

  # Helper: project a Canada-grid raster to Albers and return a tidy data frame
  # for geom_raster(). The input raster is already restricted to Canada (masked
  # to the Canada boundary in Step 6), so NO further masking is done here -- all
  # Canadian cells are shown. Drawing geom_raster() under coord_sf(crs =
  # crs_albers) (below) keeps the raster and the geom_sf boundary in the SAME
  # projection; without that explicit crs the two layers misregister and the
  # raster appears to spill past the boundary line.
  #
  # na.rm = FALSE keeps the FULL regular grid (NA cells included). Dropping NA
  # rows would leave gaps in the column positions, which makes geom_raster warn
  # that "pixels are placed at uneven horizontal intervals"; NA cells render
  # transparent (scale na.value) so keeping them is harmless.
  rast_to_df <- function(r, valname) {
    r <- terra::project(r, crs_albers, method = "near")
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    names(df) <- c("x", "y", valname)
    df
  }

  # -- Coverage climate-space panel (tiles in MAT x MAP space) ----------------
  build_cov_clim_panel <- function(zones, show_legend = FALSE, bg = "white") {
    p <- ggplot2::ggplot() +
      # Climate zones present in Canada but never sampled: solid grey.
      ggplot2::geom_tile(
        data = dplyr::filter(zones, coverage == 0),
        ggplot2::aes(x = mat_centre, y = map_centre),
        width = mat_width, height = map_width, fill = NEVER_SAMPLED_COLOUR
      ) +
      # Sampled zones: coloured by coverage.
      ggplot2::geom_tile(
        data = dplyr::filter(zones, coverage > 0),
        ggplot2::aes(x = mat_centre, y = map_centre, fill = coverage),
        width = mat_width, height = map_width
      ) +
      cov_fill_scale + clim_labs + clim_theme_fn(bg)
    if (show_legend) p + legend_inside else p + ggplot2::guides(fill = "none")
  }

  # -- Coverage geographic-map panel (no legend; no title) --------------------
  build_cov_map_panel <- function(cov_rast, bg = "white") {
    df <- rast_to_df(cov_rast, "coverage")
    ggplot2::ggplot() +
      # Grey Canada landmass shows through where coverage is 0 / NA.
      ggplot2::geom_sf(data = canada_albers_clim, fill = NEVER_SAMPLED_COLOUR, colour = NA) +
      ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = coverage)) +
      # Neatline mask: hide raster cell fringes that straddle the border.
      # Filled with `bg` (not a fixed "white") so it blends into the panel
      # background in both the white and grey versions of this figure.
      ggplot2::geom_sf(data = canada_exterior, fill = bg, colour = NA) +
      ggplot2::geom_sf(data = canada_albers_clim, fill = NA, colour = "grey40",
                       linewidth = 0.2) +
      cov_fill_scale +
      ggplot2::coord_sf(crs = crs_albers, xlim = xlim_can, ylim = ylim_can,
                        expand = FALSE) +
      map_theme_fn(bg) + ggplot2::guides(fill = "none")
  }

  # -- Available-climate climate-space panel (every zone coloured) ------------
  # Zones never sampled by GlobalFungi + GenBank combined (coverage == 0 in
  # panel f, shown there as solid grey) are additionally outlined in white
  # here, so the reader can see directly -- without cross-referencing panel f
  # -- which climates remain completely unsampled.
  build_freq_clim_panel <- function(show_legend = FALSE, bg = "white") {
    p <- ggplot2::ggplot(zone_freq) +
      ggplot2::geom_tile(
        ggplot2::aes(x = mat_centre, y = map_centre, fill = count_canada),
        width = mat_width, height = map_width
      ) +
      ggplot2::geom_tile(
        data = dplyr::filter(zone_freq, never_sampled_comb),
        ggplot2::aes(x = mat_centre, y = map_centre),
        width = mat_width, height = map_width,
        fill = NA, colour = "white", linewidth = 0.15
      ) +
      freq_fill_scale + clim_labs + clim_theme_fn(bg)
    if (show_legend) p + legend_inside else p + ggplot2::guides(fill = "none")
  }

  # -- Available-climate geographic-map panel (no legend; no title) -----------
  build_freq_map_panel <- function(bg = "white") {
    df <- rast_to_df(freq_rast, "frequency")
    ggplot2::ggplot() +
      # Canvas layer under the raster; filled with `bg` for the same reason
      # as the neatline mask below.
      ggplot2::geom_sf(data = canada_albers_clim, fill = bg, colour = NA) +
      ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = frequency)) +
      # Neatline mask: hide raster cell fringes that straddle the border.
      ggplot2::geom_sf(data = canada_exterior, fill = bg, colour = NA) +
      ggplot2::geom_sf(data = canada_albers_clim, fill = NA, colour = "grey40",
                       linewidth = 0.2) +
      freq_fill_scale +
      ggplot2::coord_sf(crs = crs_albers, xlim = xlim_can, ylim = ylim_can,
                        expand = FALSE) +
      map_theme_fn(bg) + ggplot2::guides(fill = "none")
  }

  # Build the six panels and assemble the composite for a given background
  # colour. Legends are drawn only inside panels d (frequency) and e
  # (coverage); titles are omitted (described in the figure caption). Row 1 =
  # maps, row 2 = climate space; top row given slightly less height so the
  # maps fill their panels; tags a-f read row-major.
  build_combined <- function(bg) {
    p_freq_map  <- build_freq_map_panel(bg)
    p_gf_map    <- build_cov_map_panel(cov_gf$cov_rast, bg)
    p_cmb_map   <- build_cov_map_panel(cov_comb$cov_rast, bg)
    p_freq_clim <- build_freq_clim_panel(show_legend = TRUE, bg = bg)            # panel d
    p_gf_clim   <- build_cov_clim_panel(cov_gf$zones, show_legend = TRUE, bg = bg)   # panel e
    p_cmb_clim  <- build_cov_clim_panel(cov_comb$zones, show_legend = FALSE, bg = bg) # panel f

    patchwork::wrap_plots(
      p_freq_map,  p_gf_map,  p_cmb_map,
      p_freq_clim, p_gf_clim, p_cmb_clim,
      ncol = 3, nrow = 2, heights = c(0.85, 1)
    ) +
      patchwork::plot_annotation(
        tag_levels = "a",
        theme = ggplot2::theme(plot.background = ggplot2::element_rect(fill = bg, colour = NA))
      )
  }

  save_fig_formats(fig_out, build_combined("white"),
                   width = 15, height = 10, dpi = 300, bg = "white")

  # Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
  # schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
  ggplot2::ggsave(fig_out_grey, build_combined(fig5_grey_bg),
                  width = 15, height = 10, dpi = 300, bg = fig5_grey_bg)

}
