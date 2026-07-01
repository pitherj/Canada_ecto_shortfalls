# =============================================================================
# 20_sampling_maps.R  —  Sampling-location overview maps
# =============================================================================
# PURPOSE
#   Produce the two whole-Canada sampling maps used in the manuscript:
#     * Figure 1  — EcM sequence sampling locations (GlobalFungi + GenBank) on a
#                   clean administrative-boundary basemap.
#     * Figure S3 — GBIF physical EcM specimen records in Canada, split by
#                   whether the genus also has sequence data in Canada.
#
# INPUTS
#   paths$emf_data              combined EcM dataset (04_combine_ecm_dataset.R)
#   paths$gbif_ecm              GBIF EcM specimen records         (09_linnean.R)
#   paths$gbif_ecm_nosequence   GBIF records for genera lacking sequence data
#   paths$canada_bound,
#   paths$canada_bound_raw,
#   paths$lakes_albers          spatial basemap layers            (01_spatial_data.R)
#
# OUTPUTS (figures/)
#   Figure-01_sampling_map.png     (paths$fig_sampling_map)
#   Figure-S3_gbif_specimens.png   (paths$fig_gbif_specimens)
#
# NOTE
#   Coordinates are validated against the Canada boundary once, at dataset
#   assembly (04_combine_ecm_dataset.R), so here we simply filter on
#   coord_in_canada == TRUE. Each map is skipped if its file already exists.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)

sf::sf_use_s2(FALSE)

MAP_DPI <- 300

# ---- Which maps need building? ----------------------------------------------
need_g0   <- !file.exists(paths$fig_sampling_map)
need_gbif <- !file.exists(paths$fig_gbif_specimens)

if (!(need_g0 || need_gbif)) {
  ts("Both sampling maps already exist — nothing to do.")
}

# =============================================================================
# Shared basemap layers
# =============================================================================
if (need_g0 || need_gbif) {
  ts("Loading base spatial layers...")

  canada_bound  <- sf::st_read(paths$canada_bound, quiet = TRUE) |> sf::st_make_valid()
  canada_albers <- sf::st_transform(canada_bound, crs_albers)

  # Provinces/territories (GADM level 1) for the white administrative basemap
  provinces_albers <- terra::vect(paths$canada_bound_raw) |>
    terra::project(crs_albers) |>
    sf::st_as_sf()

  # Major lakes, clipped to Canada (optional — skipped if the file is missing)
  # st_agr(x) <- "constant" declares that lake attributes belong to each polygon
  # and are not redistributed by the clip; it only silences an sf warning.
  lakes_sf <- tryCatch({
    lakes_raw <- sf::st_read(paths$lakes_albers, quiet = TRUE) |>
      sf::st_transform(crs_albers) |> sf::st_make_valid() |> sf::st_buffer(0)
    sf::st_agr(lakes_raw) <- "constant"
    sf::st_intersection(lakes_raw, sf::st_buffer(sf::st_make_valid(canada_albers), 0))
  }, error = function(e) { ts("  Lakes not found — skipping."); NULL })

  # Common plot extent (Canada bbox + a small margin)
  bbox_can <- sf::st_bbox(canada_albers)
  xlim_can <- c(bbox_can[["xmin"]] - 150000, bbox_can[["xmax"]] + 150000)
  ylim_can <- c(bbox_can[["ymin"]] - 150000, bbox_can[["ymax"]] + 150000)
}

# A reusable ggplot theme for both maps
map_theme <- ggplot2::theme_void() +
  ggplot2::theme(
    legend.justification = c("right", "top"),
    legend.background = ggplot2::element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
    legend.margin     = ggplot2::margin(5, 6, 5, 6),
    legend.text       = ggplot2::element_text(size = 10),
    legend.title      = ggplot2::element_text(size = 11, face = "bold"),
    panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background   = ggplot2::element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# Figure 1 — EcM sequence sampling locations (GlobalFungi + GenBank)
# =============================================================================
ts("Figure 1: sequence sampling-location map...")
if (!need_g0) {
  ts("  Already exists — skipping.")
} else {
  # One point per unique coordinate per source (coordinates pre-validated)
  locs_g0 <- dplyr::bind_rows(
    dplyr::filter(emf, source == "GenBank",     coord_in_canada == TRUE) |>
      dplyr::distinct(lat, lon) |> dplyr::mutate(source = "GenBank"),
    dplyr::filter(emf, source == "GlobalFungi", coord_in_canada == TRUE) |>
      dplyr::distinct(lat, lon) |> dplyr::mutate(source = "GlobalFungi")
  ) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(crs_albers)

  p_g0 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = provinces_albers, fill = "white", colour = "grey55", linewidth = 0.3) +
    (if (!is.null(lakes_sf))
      ggplot2::geom_sf(data = lakes_sf, fill = "lightsteelblue2", colour = "steelblue4",
                       linewidth = 0.15, alpha = 0.9) else NULL) +
    ggplot2::geom_sf(data = locs_g0, ggplot2::aes(colour = source, shape = source),
                     fill = "white", size = 2.5, stroke = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(breaks = c("GlobalFungi", "GenBank"),
      values = c("GlobalFungi" = "#0072B2", "GenBank" = "#E69F00"), name = "Sample source") +
    ggplot2::scale_shape_manual(breaks = c("GlobalFungi", "GenBank"),
      values = c("GlobalFungi" = 21L, "GenBank" = 24L), name = "Sample source") +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(shape = c(21L, 24L), fill = "white", size = 4)), shape = "none") +
    ggplot2::coord_sf(xlim = xlim_can, ylim = ylim_can, expand = FALSE) +
    map_theme + ggplot2::theme(legend.position = c(0.85, 0.75))

  ggplot2::ggsave(paths$fig_sampling_map, p_g0, width = 12, height = 9, dpi = MAP_DPI)
  ts(sprintf("  Saved %s", basename(paths$fig_sampling_map)))
}

# =============================================================================
# Figure S3 — GBIF physical EcM specimen records
# =============================================================================
# Blue circles: genera that also have sequence data in Canada.
# Orange triangles: genera without sequence data in Canada.
# =============================================================================
ts("Figure S3: GBIF physical specimen map...")
if (!need_gbif) {
  ts("  Already exists — skipping.")
} else if (!file.exists(paths$gbif_ecm) || !file.exists(paths$gbif_ecm_nosequence)) {
  ts("  GBIF specimen files not found — run 09_linnean.R first.")
} else {
  gbif_labels <- c("Genera with sequence data in Canada",
                   "Genera without sequence data in Canada")

  gbif_combined <- dplyr::bind_rows(
    readr::read_csv(paths$gbif_ecm,            show_col_types = FALSE) |>
      dplyr::mutate(category = gbif_labels[1]),
    readr::read_csv(paths$gbif_ecm_nosequence, show_col_types = FALSE) |>
      dplyr::mutate(category = gbif_labels[2])
  ) |>
    dplyr::filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
    dplyr::distinct(decimalLatitude, decimalLongitude, category, .keep_all = TRUE) |>
    dplyr::mutate(category = factor(category, levels = gbif_labels))

  ts(sprintf("  Georeferenced GBIF records — with data: %d  |  without: %d",
             sum(gbif_combined$category == gbif_labels[1]),
             sum(gbif_combined$category == gbif_labels[2])))

  gbif_sf <- sf::st_as_sf(gbif_combined, coords = c("decimalLongitude", "decimalLatitude"),
                          crs = 4326) |>
    sf::st_filter(canada_bound) |>
    sf::st_transform(crs_albers)

  p_gbif <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = provinces_albers, fill = "white", colour = "grey55", linewidth = 0.3) +
    (if (!is.null(lakes_sf))
      ggplot2::geom_sf(data = lakes_sf, fill = "lightsteelblue2", colour = "steelblue4",
                       linewidth = 0.15, alpha = 0.9) else NULL) +
    ggplot2::geom_sf(data = gbif_sf, ggplot2::aes(colour = category, shape = category),
                     fill = "white", size = 2.5, stroke = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(breaks = gbif_labels, values = c("#0072B2", "#E69F00"),
                                 name = "GBIF genus category") +
    ggplot2::scale_shape_manual(breaks = gbif_labels, values = c(21L, 24L),
                                name = "GBIF genus category") +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(shape = c(21L, 24L), fill = "white", size = 4)), shape = "none") +
    ggplot2::coord_sf(xlim = xlim_can, ylim = ylim_can, expand = FALSE) +
    map_theme + ggplot2::theme(legend.position = c(0.85, 0.85))

  ggplot2::ggsave(paths$fig_gbif_specimens, p_gbif, width = 12, height = 9, dpi = MAP_DPI)
  ts(sprintf("  Saved %s", basename(paths$fig_gbif_specimens)))
}

ts("20_sampling_maps.R complete.")
