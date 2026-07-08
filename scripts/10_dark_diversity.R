# =============================================================================
# 10_dark_diversity.R  —  Dark ('undescribed') EcM fungal taxa map (Figure S1)
# =============================================================================
# PURPOSE
#   Map the estimated percentage of undescribed ("dark") ectomycorrhizal fungal
#   taxa across Canada, from the van Galen et al. (2025) dark-taxa raster. This
#   is Figure S1 of the manuscript's Supplemental Materials.
#
# INPUTS
#   paths$van_galen_tif   van Galen et al. (2025) dark-taxa GeoTIFF (data_raw/;
#                         band `percentage_dark_taxa`, 0-100 %)
#   paths$canada_bound,
#   paths$lakes_albers    basemap layers (01_spatial_data.R)
#
# OUTPUT (figures/)
#   Figure-S1_dark_diversity.png        (paths$fig_dark_diversity)      -- white bg, used in SI
#   Figure-S1_dark_diversity_grey.png   (paths$fig_dark_diversity_grey) -- #F2F2F2 bg, Figure 5 panel source
#
# SOURCE
#   van Galen LG et al. (2025) The biogeography and conservation of Earth's
#   'dark' ectomycorrhizal fungi. Current Biology 35: R563-R574.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)

sf::sf_use_s2(FALSE)

DARK_AGGREG_FACTOR <- 100   # aggregate the raster to ~100 km before reprojection
MAP_DPI            <- 300

if (file.exists(paths$fig_dark_diversity) && file.exists(paths$fig_dark_diversity_grey)) {
  ts("Figure S1 (dark diversity, both versions) already exists — skipping.")
} else if (!file.exists(paths$van_galen_tif)) {
  stop("van Galen raster not found: ", paths$van_galen_tif)
} else {

  ts("Figure S1: building dark-diversity map...")

  # ---- Basemap layers --------------------------------------------------------
  canada_bound  <- sf::st_read(paths$canada_bound, quiet = TRUE) |> sf::st_make_valid()
  canada_albers <- sf::st_transform(canada_bound, crs_albers)
  canada_v      <- terra::vect(canada_albers)

  lakes_sf <- tryCatch({
    lakes_raw <- sf::st_read(paths$lakes_albers, quiet = TRUE) |>
      sf::st_transform(crs_albers) |> sf::st_make_valid() |> sf::st_buffer(0)
    sf::st_agr(lakes_raw) <- "constant"
    sf::st_intersection(lakes_raw, sf::st_buffer(sf::st_make_valid(canada_albers), 0))
  }, error = function(e) { ts("  Lakes not found — skipping."); NULL })

  bbox_can <- sf::st_bbox(canada_albers)
  xlim_can <- c(bbox_can[["xmin"]] - 150000, bbox_can[["xmax"]] + 150000)
  ylim_can <- c(bbox_can[["ymin"]] - 150000, bbox_can[["ymax"]] + 150000)

  # ---- Load and prepare the dark-taxa raster ---------------------------------
  # The GeoTIFF has three bands; we use `percentage_dark_taxa` (0-100 %).
  vg_raw    <- terra::rast(paths$van_galen_tif)[["percentage_dark_taxa"]]
  vg_agg    <- terra::aggregate(vg_raw, fact = DARK_AGGREG_FACTOR, fun = "mean", na.rm = TRUE)
  vg_albers <- terra::project(vg_agg, crs_albers, method = "bilinear")
  vg_canada <- terra::mask(vg_albers, canada_v)

  # Keep the full regular grid (na.rm = FALSE) so geom_raster() sees evenly
  # spaced pixels; NA cells render transparent via na.value below.
  vg_df <- terra::as.data.frame(vg_canada, xy = TRUE, na.rm = FALSE) |>
    dplyr::rename(dark_pct = percentage_dark_taxa)

  vg_min    <- floor(min(vg_df$dark_pct, na.rm = TRUE))
  vg_max    <- ceiling(max(vg_df$dark_pct, na.rm = TRUE))
  vg_breaks <- pretty(c(vg_min, vg_max), n = 5)

  # ---- Plot --------------------------------------------------------------
  # build_plot(bg) returns the map with panel/plot background colour `bg`;
  # called once for the white (manuscript/SI) version and once for the grey
  # (#F2F2F2, Figure 5 schematic panel) version.
  build_plot <- function(bg) {
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = canada_albers, fill = "grey95", colour = NA) +
      ggplot2::geom_raster(data = vg_df, ggplot2::aes(x = x, y = y, fill = dark_pct)) +
      (if (!is.null(lakes_sf))
        ggplot2::geom_sf(data = lakes_sf, fill = "dodgerblue", colour = "darkblue",
                         linewidth = 0.2, alpha = 0.9) else NULL) +
      ggplot2::geom_sf(data = canada_albers, fill = NA, colour = "black", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c(
        option = "F", direction = -1, na.value = "transparent",
        name = "Dark EcM fungal taxa (%)",
        limits = c(vg_min, vg_max), breaks = vg_breaks,
        guide = ggplot2::guide_colorbar(barwidth = 8, barheight = 0.5,
                                        title.position = "top", title.hjust = 0.5)) +
      ggplot2::coord_sf(crs = crs_albers, xlim = xlim_can, ylim = ylim_can, expand = FALSE) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        legend.position  = "bottom",
        panel.background = ggplot2::element_rect(fill = bg, colour = NA),
        plot.background  = ggplot2::element_rect(fill = bg, colour = NA))
  }

  ggplot2::ggsave(paths$fig_dark_diversity, build_plot("white"),
                  width = 12, height = 9, dpi = MAP_DPI, bg = "white")
  ts(sprintf("  Saved %s", basename(paths$fig_dark_diversity)))

  # Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
  # schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
  ggplot2::ggsave(paths$fig_dark_diversity_grey, build_plot(fig5_grey_bg),
                  width = 12, height = 9, dpi = MAP_DPI, bg = fig5_grey_bg)
  ts(sprintf("  Saved %s", basename(paths$fig_dark_diversity_grey)))
}

ts("10_dark_diversity.R complete.")
