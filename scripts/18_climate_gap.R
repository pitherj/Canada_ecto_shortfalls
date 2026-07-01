# =============================================================================
# 18_climate_gap.R — Hutchinsonian Shortfall: climate-space sampling coverage
# =============================================================================
# PURPOSE
#   Quantify and visualize how completely Canada's *climate space* has been
#   sampled for ectomycorrhizal (EcM) fungi, and where in geographic space the
#   well- vs. poorly-covered climates fall.
#
#   The representation is adapted from Figure 6 in the following preprint:
# Blitz the Gap: a nation-wide effort to guide citizen science toward the needs
# of biodiversity science. K. Hebert et al.
# https://doi.org/10.32942/X2T09G
#
# METHOD (matches the colleagues' coverage definition)
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
#
# SCOPES
#   - "GlobalFungi": GF samples with >= 1 EcM SH code (named or not). Because the
#     downstream dataset (emf) is already EcM-only, every GF sample qualifies.
#   - "GlobalFungi + GenBank": the GF samples above plus GenBank EcM records.
#   In both scopes the sampling unit for coverage is the climate grid cell
#   (a cell counts once regardless of how many samples fall in it).
#
# INPUTS
#   emf                                        (auto-loaded by 00_setup.R)
#   data_raw/climate/wc2.1_country/CAN_wc2.1_30s_bio.tif   (paths$climate_raster)
#   data_derived/spatial/canada_simple.gpkg                   (paths$canada_bound)
#
# OUTPUT
#   figures/hutchinsonian_climate_gap.png      (replaces Figure S11 in the SI)
#
# NOTES
#   - Deterministic: no random sampling, so no seed is required.
#   - Sentinel-file guard: the figure is rebuilt only if it does not yet exist.
#     Delete it to force regeneration.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

# ---- Tunable parameters (exposed at top, mirroring 10_linnean_inext.R) -----
AGG_FACTOR      <- 20L   # 30 arc-sec -> ~10 arc-min (20 * 30" = 600" = 10')
N_CLIMATE_BINS  <- 50L   # bins per climate axis (50 x 50 grid of climate zones)

# Colour palettes match the colleagues' figures (via the colorspace package):
#   - "Batlow" sequential  -> available climate frequency (their Figs S4/S5)
#   - "Temps"  divergingx   -> sampling coverage, reversed so poor = red and
#                              full = teal (their Fig 6)
FREQ_PALETTE         <- "Batlow"
COVERAGE_PALETTE     <- "Temps"
NEVER_SAMPLED_COLOUR <- "grey85"   # climate present in Canada but never sampled

fig_out <- paths$fig_climate_gap

# Sentinel guard: skip the whole script if the figure already exists.
if (file.exists(fig_out)) {
  ts(sprintf("Figure already exists, skipping: %s", basename(fig_out)))
} else {

  sf::sf_use_s2(FALSE)  # GADM-derived Canada polygon has minor topology issues

  # ===========================================================================
  # Step 1 — Build the aggregated, Canada-masked climate grid
  # ===========================================================================
  ts("Step 1: Loading and aggregating WorldClim climate layers...")

  if (!file.exists(paths$climate_raster))
    stop("Climate raster not found: ", paths$climate_raster, call. = FALSE)

  clim_full <- terra::rast(paths$climate_raster)

  # WorldClim layer naming varies by download method; match BIO1 / BIO12 by the
  # numeric suffix rather than assuming a fixed band order (same defensive trick
  # used in 17_hutchinsonian.R).
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
    ts("  MAT divided by 10 (detected x10 integer encoding)")
  }

  # Aggregate 30" -> ~10' (mean of contributing fine cells), then mask to the
  # project's canonical Canada boundary so the climate cloud and the coverage
  # denominator use Canadian land cells only.
  ts(sprintf("  Aggregating by factor %d (30\" -> ~10')...", AGG_FACTOR))
  clim <- terra::aggregate(clim, fact = AGG_FACTOR, fun = "mean", na.rm = TRUE)

  canada_wgs84 <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
    sf::st_transform(4326) |>
    sf::st_make_valid()
  canada_vect <- terra::vect(canada_wgs84)
  clim <- terra::mask(terra::crop(clim, canada_vect), canada_vect)

  # Per-cell climate values, aligned to cell indices (NA outside Canada).
  cell_vals <- terra::values(clim)            # matrix: columns MAT, MAP
  cell_ok   <- !is.na(cell_vals[, "MAT"]) & !is.na(cell_vals[, "MAP"])
  ts(sprintf("  Canadian climate cells (~10' grain): %d", sum(cell_ok)))

  # ===========================================================================
  # Step 2 — Bin climate space into an N x N grid of "climate zones"
  # ===========================================================================
  ts(sprintf("Step 2: Binning climate space into %d x %d zones...",
             N_CLIMATE_BINS, N_CLIMATE_BINS))

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
  ts(sprintf("  Occupied climate zones: %d", nrow(zone_canada)))

  # ===========================================================================
  # Step 3 — Coverage per scope
  # ===========================================================================
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

    ts(sprintf("  [%s] sampled cells: %d | zones with coverage > 0: %d / %d",
               scope_label, length(cells), sum(zones$coverage > 0),
               nrow(zones)))

    list(zones = zones, cov_rast = cov_rast)
  }

  ts("Step 3: Extracting sample coordinates and computing coverage...")

  # GlobalFungi: every EcM-only GF sample qualifies (>= 1 EcM SH by construction).
  gf_pts <- emf |>
    dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE) |>
    dplyr::distinct(lat, lon)

  # GenBank EcM records with validated Canadian coordinates.
  gb_pts <- emf |>
    dplyr::filter(source == "GenBank", coord_in_canada == TRUE) |>
    dplyr::distinct(lat, lon)

  # Combined scope: union of distinct coordinates.
  comb_pts <- dplyr::bind_rows(gf_pts, gb_pts) |>
    dplyr::distinct(lat, lon)

  ts(sprintf("  GF coords: %d | GB coords: %d | combined: %d",
             nrow(gf_pts), nrow(gb_pts), nrow(comb_pts)))

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
    )

  # Climate-frequency raster: each Canadian cell takes its zone's count_canada.
  freq_lookup <- stats::setNames(zone_freq$count_canada, zone_freq$bin)
  freq_rast   <- clim[["MAT"]]
  terra::values(freq_rast) <- unname(freq_lookup[cell_bin])
  names(freq_rast) <- "frequency"

  # ===========================================================================
  # Step 4 -- Build the six panels (2 rows x 3 columns)
  # ===========================================================================
  #   Columns: Available climate frequency | GlobalFungi | GF + GenBank
  #   Rows:    geographic maps (top)       | climate space (bottom)
  ts("Step 4: Building figure panels...")

  canada_albers <- sf::st_transform(canada_wgs84, crs_albers)

  # Tight Canada bounding box, used by coord_sf so each map fills its panel.
  bbox_can <- sf::st_bbox(canada_albers)
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
    sf::st_union(sf::st_buffer(sf::st_make_valid(canada_albers), 0))
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
  clim_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  map_theme <- ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.margin      = ggplot2::margin(2, 2, 2, 2),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
  clim_labs <- ggplot2::labs(x = "Mean annual temperature (\u00b0C)",
                             y = "Mean annual precipitation (mm)")

  # Helper: project a Canada-grid raster to Albers and return a tidy data frame
  # for geom_raster(). The input raster is already restricted to Canada (masked
  # to the Canada boundary in Step 1), so NO further masking is done here -- all
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
  build_cov_clim_panel <- function(zones, show_legend = FALSE) {
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
      cov_fill_scale + clim_labs + clim_theme
    if (show_legend) p + legend_inside else p + ggplot2::guides(fill = "none")
  }

  # -- Coverage geographic-map panel (no legend; no title) --------------------
  build_cov_map_panel <- function(cov_rast) {
    df <- rast_to_df(cov_rast, "coverage")
    ggplot2::ggplot() +
      # Grey Canada landmass shows through where coverage is 0 / NA.
      ggplot2::geom_sf(data = canada_albers, fill = NEVER_SAMPLED_COLOUR, colour = NA) +
      ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = coverage)) +
      # Neatline mask: hide raster cell fringes that straddle the border.
      ggplot2::geom_sf(data = canada_exterior, fill = "white", colour = NA) +
      ggplot2::geom_sf(data = canada_albers, fill = NA, colour = "grey40",
                       linewidth = 0.2) +
      cov_fill_scale +
      ggplot2::coord_sf(crs = crs_albers, xlim = xlim_can, ylim = ylim_can,
                        expand = FALSE) +
      map_theme + ggplot2::guides(fill = "none")
  }

  # -- Available-climate climate-space panel (every zone coloured) ------------
  build_freq_clim_panel <- function(show_legend = FALSE) {
    p <- ggplot2::ggplot(zone_freq) +
      ggplot2::geom_tile(
        ggplot2::aes(x = mat_centre, y = map_centre, fill = count_canada),
        width = mat_width, height = map_width
      ) +
      freq_fill_scale + clim_labs + clim_theme
    if (show_legend) p + legend_inside else p + ggplot2::guides(fill = "none")
  }

  # -- Available-climate geographic-map panel (no legend; no title) -----------
  build_freq_map_panel <- function() {
    df <- rast_to_df(freq_rast, "frequency")
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = canada_albers, fill = "white", colour = NA) +
      ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = frequency)) +
      # Neatline mask: hide raster cell fringes that straddle the border.
      ggplot2::geom_sf(data = canada_exterior, fill = "white", colour = NA) +
      ggplot2::geom_sf(data = canada_albers, fill = NA, colour = "grey40",
                       linewidth = 0.2) +
      freq_fill_scale +
      ggplot2::coord_sf(crs = crs_albers, xlim = xlim_can, ylim = ylim_can,
                        expand = FALSE) +
      map_theme + ggplot2::guides(fill = "none")
  }

  # Build the six panels. Legends are drawn only inside panels d (frequency)
  # and e (coverage); titles are omitted (described in the figure caption).
  p_freq_map  <- build_freq_map_panel()
  p_gf_map    <- build_cov_map_panel(cov_gf$cov_rast)
  p_cmb_map   <- build_cov_map_panel(cov_comb$cov_rast)
  p_freq_clim <- build_freq_clim_panel(show_legend = TRUE)            # panel d
  p_gf_clim   <- build_cov_clim_panel(cov_gf$zones, show_legend = TRUE)   # panel e
  p_cmb_clim  <- build_cov_clim_panel(cov_comb$zones, show_legend = FALSE) # panel f

  # Assemble 2 x 3 (row 1 = maps, row 2 = climate space). Top row given slightly
  # less height (heights) so the maps fill their panels; tags a-f read row-major.
  combined <- patchwork::wrap_plots(
    p_freq_map,  p_gf_map,  p_cmb_map,
    p_freq_clim, p_gf_clim, p_cmb_clim,
    ncol = 3, nrow = 2, heights = c(0.85, 1)
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  ggplot2::ggsave(fig_out, combined, width = 15, height = 10, dpi = 300)
  ts(sprintf("Saved -> %s", basename(fig_out)))

  ts("18_climate_gap.R complete.")
}
