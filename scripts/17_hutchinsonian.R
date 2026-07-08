# =============================================================================
# Hutchinsonian Shortfall (Updated): Environmental and Ecozone Coverage
# =============================================================================
# What proportion of Canadian ecological space with EcM host habitat has
# been sampled for EcM fungi?
#
# Updates from original 04_hutchinsonian.R:
#   - BIEN-based host raster (0.5°, from 08_host_rasters.R) replaces
#     USTreeAtlas 1° raster
#   - Both GenBank and GlobalFungi coordinates used for sampling points
#     (GenBank lat/lon parsed from lat_lon_gb where available)
#   - ecozone climate space analysis added using WorldClim data
#
# Design history: this script previously also computed an ecoregion-level
# habitat-coverage metric (hutchinsonian_ecoregion_summary.csv, using the
# ecoregion habitat layer from 08_host_rasters.R). It was never reported in
# the manuscript or supplemental materials -- all reported Hutchinsonian
# coverage is ecozone-level -- so it was removed 2026-07-05 (pruned by
# Jason's request) along with its upstream input in 08_host_rasters.R.
#
# Prerequisite files:
#   data_derived/bien_host_richness_0.5deg.tif       (08_host_rasters.R)
#   data_derived/bien_host_species_stack.tif         (08_host_rasters.R)
#   data_derived/ecm_native_canada_host_species.csv  (06_host_species.R)
#   data_derived/eltonian_host_matching.csv          (19_eltonian.R)
#   data_raw/climate/wc2.1_country/CAN_wc2.1_30s_bio.tif
#
# Outputs:
#   data_derived/bien_host_data_richness_0.5deg.tif    — host spp. with EcM data
#   data_derived/bien_host_data_proportion_0.5deg.tif  — proportion with EcM data
#   data_derived/hutchinsonian_raster_summary.csv
#   data_derived/hutchinsonian_ecozone_summary.csv     — ecozone sampling-threshold
#                                                    counts (Table S16b)
#   figures/Figure-04_host_bivariate_map.png        (paths$fig_host_bivariate)      -- white bg, used in manuscript
#   figures/Figure-04_host_bivariate_map_grey.png   (paths$fig_host_bivariate_grey) -- #F2F2F2 bg, Figure 5 panel source
#   figures/Figure-S4_ecozone_sampling_map.png      (paths$fig_ecozone_sampling)      -- white bg, used in SI
#   figures/Figure-S4_ecozone_sampling_map_grey.png (paths$fig_ecozone_sampling_grey) -- #F2F2F2 bg, Figure 5 panel source
#   figures/hutchinsonian_climate_space.png (currently commented out — see Step 8)
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

sf::sf_use_s2(FALSE)

# ---- Step 1: Prepare sampling points ----------------------------------------

ts("Step 1: Extracting sampling locations from EMF dataset...")

# GlobalFungi: use coord_in_canada == TRUE (coordinates validated within GADM Canada)
gf_locs <- emf |>
  dplyr::filter(source == "GlobalFungi", coord_in_canada == TRUE) |>
  dplyr::distinct(sample_ID, lat, lon)

# GenBank: same filter
gb_parsed <- emf |>
  dplyr::filter(source == "GenBank", coord_in_canada == TRUE) |>
  dplyr::distinct(lat, lon)

all_locs <- dplyr::bind_rows(gf_locs, gb_parsed) |>
  dplyr::distinct(lat, lon)

ts(sprintf("  GlobalFungi locations: %d", nrow(gf_locs)))
ts(sprintf("  GenBank locations (stored coords): %d", nrow(gb_parsed)))
ts(sprintf("  Unique locations combined: %d", nrow(all_locs)))

# Canada boundary (needed for map steps below)
canada_bound <- sf::st_read(paths$canada_bound, quiet = TRUE)

# coord_in_canada == TRUE guarantees coordinates were validated within the GADM
# Canada boundary at dataset assembly (04_combine_ecm_dataset.R).
locs_sf <- sf::st_as_sf(all_locs, coords = c("lon", "lat"), crs = 4326)
ts(sprintf("  Unique locations as sf: %d", nrow(locs_sf)))

# ---- Step 2: Load spatial prerequisites -------------------------------------

ts("Step 2: Loading spatial prerequisites...")

if (!file.exists(paths$bien_richness)) {
  stop("Host richness raster not found: ", paths$bien_richness,
       "\nRun 08_host_rasters.R first.")
}
richness_wgs84 <- terra::rast(paths$bien_richness)
names(richness_wgs84) <- "richness"
ts(sprintf("  Richness raster: %d x %d cells", nrow(richness_wgs84), ncol(richness_wgs84)))

# ---- Step 3: Identify host species with EcM data in our dataset -------------

ts("Step 3: Identifying host species with EcM data in dataset...")

host_tbl <- readr::read_csv(paths$host_species, show_col_types = FALSE)
em_canada_species <- host_tbl$species

# Use eltonian host matching output if available; else build from host fields
matching_file <- paths$eltonian_host_matching

# Host-name cleaning is handled by the shared canonicalize_host() helper in
# 00_setup.R. Underscores in GlobalFungi `dominant_plant_species` /
# `other_plant_species` values are converted to spaces before tokenization.
clean_host_name <- function(x) canonicalize_host(gsub("_", " ", x))

if (file.exists(matching_file)) {
  ts("  Using eltonian_host_matching.csv for host species with data...")
  host_matching   <- readr::read_csv(matching_file, show_col_types = FALSE)
  em_species_with_data <- unique(
    host_matching$host_clean[host_matching$matched & !is.na(host_matching$host_clean)]
  )
} else {
  ts("  eltonian_host_matching.csv not found; building from host fields...")
  # dominant_plant_species / other_plant_species: root samples only (GlobalFungi)
  emf_gf_root <- dplyr::filter(emf, source == "GlobalFungi", sample_type == "root")
  all_hosts <- unique(c(
    clean_host_name(emf_gf_root$dominant_plant_species[!is.na(emf_gf_root$dominant_plant_species)]),
    clean_host_name(emf_gf_root$other_plant_species[!is.na(emf_gf_root$other_plant_species)]),
    clean_host_name(emf$host_taxon[!is.na(emf$host_taxon)])
  ))
  em_species_with_data <- unique(all_hosts[!is.na(all_hosts) & all_hosts %in% em_canada_species])
}

ts(sprintf("  Host species represented in our dataset: %d / %d",
           length(em_species_with_data), length(em_canada_species)))

# Tree / non-tree host species subsets (used for bivariate maps)
em_canada_tree_species    <- host_tbl$species[host_tbl$growth_form %in% "tree"]
em_canada_nontree_species <- host_tbl$species[!host_tbl$growth_form %in% "tree"]
em_species_with_data_tree    <- em_species_with_data[em_species_with_data %in% em_canada_tree_species]
em_species_with_data_nontree <- em_species_with_data[em_species_with_data %in% em_canada_nontree_species]
ts(sprintf("  Tree host species:     %d total / %d with EcM data",
           length(em_canada_tree_species), length(em_species_with_data_tree)))
ts(sprintf("  Non-tree host species: %d total / %d with EcM data",
           length(em_canada_nontree_species), length(em_species_with_data_nontree)))

# ---- Step 4: Load species stack and derive all richness rasters --------------
# The per-species binary stack (one layer per host species, named by species)
# is produced by 08_host_rasters.R alongside bien_richness_0.5deg.tif.
# All subset richness rasters are derived by filtering layers + summing, with
# no range re-rasterization needed here.

ts("Step 4: Loading species stack...")

if (!file.exists(paths$bien_species_stack)) {
  stop("Species stack not found: ", paths$bien_species_stack,
       "\nRun 08_host_rasters.R first.", call. = FALSE)
}

species_stack <- terra::rast(paths$bien_species_stack)
canada_vect   <- terra::vect(sf::st_transform(canada_bound, 4326))

# Normalise layer names: BIEN stores species with underscores ("Abies_amabilis")
# but all other species lists in this project use spaces ("Abies amabilis").
names(species_stack) <- gsub("_", " ", names(species_stack))

ts(sprintf("  Stack loaded: %d layers (%d species)",
           terra::nlyr(species_stack), terra::nlyr(species_stack)))

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
ts(sprintf("  Deriving all-species data richness (%d spp. with data)...",
           length(em_species_with_data)))
data_richness_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data, zero_to_na = FALSE
)
names(data_richness_wgs84) <- "data_richness"
terra::writeRaster(data_richness_wgs84, paths$bien_data_rich, overwrite = TRUE)
ts(sprintf("  Saved -> %s", basename(paths$bien_data_rich)))

# Tree-only richness (all tree host species)
ts(sprintf("  Deriving tree richness (%d tree spp.)...",
           length(em_canada_tree_species)))
richness_tree_wgs84 <- sum_species_layers(
  species_stack, em_canada_tree_species, zero_to_na = TRUE
)
names(richness_tree_wgs84) <- "richness"

# Tree-only data richness (tree host species in our EcM dataset)
ts(sprintf("  Deriving tree data richness (%d tree spp. with data)...",
           length(em_species_with_data_tree)))
data_richness_tree_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data_tree, zero_to_na = FALSE
)
names(data_richness_tree_wgs84) <- "data_richness"

# Non-tree richness and data richness
ts(sprintf("  Deriving non-tree richness (%d non-tree spp.)...",
           length(em_canada_nontree_species)))
richness_nontree_wgs84 <- sum_species_layers(
  species_stack, em_canada_nontree_species, zero_to_na = TRUE
)
names(richness_nontree_wgs84) <- "richness"

ts(sprintf("  Deriving non-tree data richness (%d non-tree spp. with data)...",
           length(em_species_with_data_nontree)))
data_richness_nontree_wgs84 <- sum_species_layers(
  species_stack, em_species_with_data_nontree, zero_to_na = FALSE
)
names(data_richness_nontree_wgs84) <- "data_richness"

# ---- Step 5: Proportion raster ----------------------------------------------

ts("Step 5: Computing proportion raster...")

proportion_wgs84 <- terra::clamp(
  data_richness_wgs84 / richness_wgs84,
  lower = 0, upper = 1
)
# Cells where richness is NA (outside Canada) remain NA
names(proportion_wgs84) <- "proportion"

terra::writeRaster(proportion_wgs84, paths$bien_proportion, overwrite = TRUE)
ts(sprintf("  Saved proportion raster -> %s", basename(paths$bien_proportion)))

# Summary statistics
prop_vals <- terra::values(proportion_wgs84)
prop_vals <- prop_vals[!is.na(prop_vals)]
ts(sprintf("  Cells with host habitat: %d  |  cells with any data: %d (%.1f%%)",
           length(prop_vals),
           sum(prop_vals > 0),
           100 * sum(prop_vals > 0) / length(prop_vals)))

# ---- Step 6: Raster summary --------------------------------------------------

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
                 file.path(paths$out_hutchinsonian, "hutchinsonian_raster_summary.csv"))
ts("Saved hutchinsonian_raster_summary.csv")

# ---- Step 7: Maps ------------------------------------------------------------

ts("Step 7: Setting up spatial layers for map production...")

canada_albers <- sf::st_transform(canada_bound, crs_albers)

# -- Shared ecozone / ecoregion setup (used by sampling map and ecozone climate space plot) ----
ts("  Loading ecozone and ecoregion layers...")

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
  error = function(e) { ts("  Lakes file not found — skipping."); NULL }
)

# Sampling locations with source label (used by both sampling map and ecozone climate space plot).
# GenBank is bound first so it is drawn underneath; GlobalFungi is bound second so it
# is drawn on top — matching the layer order used in sampling_map_group0.png.
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
ts("  Creating sampling locations map...")

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
    # white fill, colour-blind-safe outline (matching sampling_map_group0 style)
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
ts(sprintf("  Saved %s", basename(paths$fig_ecozone_sampling)))

# Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
# schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
ggplot2::ggsave(paths$fig_ecozone_sampling_grey, build_sampling_map(fig5_grey_bg),
                width = 12, height = 9, dpi = 300, bg = fig5_grey_bg)
ts(sprintf("  Saved %s", basename(paths$fig_ecozone_sampling_grey)))

# -- Bivariate map (host richness x proportion with data) ---------------------
# Two panels: (1) all EcM host species, (2) tree host species only.
# Breaks are computed independently for each panel so each uses its own tertiles.
ts("  Creating bivariate maps...")

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
ts(sprintf("  Saved %s", basename(paths$fig_host_bivariate)))

# Grey (#F2F2F2) version: source panel for the hand-assembled Figure 5
# schematic (see fig5_grey_bg in 00_setup.R); not used elsewhere.
ggplot2::ggsave(paths$fig_host_bivariate_grey, build_bivariate_figure(fig5_grey_bg),
                width = 10, height = 14, dpi = 300, bg = fig5_grey_bg)
ts(sprintf("  Saved %s", basename(paths$fig_host_bivariate_grey)))

# ---- Step 8: ecozone climate space analysis -------------------------
# Superseded by the new MAT/MAP sampling maps. Retained for reference.
# if (FALSE) {  # ----- BEGIN COMMENTED-OUT BLOCK -----

# ts("Step 9: ecozone climate space analysis...")
# 
# if (!file.exists(paths$climate_raster)) {
#   ts(sprintf("  Climate raster not found: %s", paths$climate_raster))
#   ts("  Skipping climate space analysis.")
# } else {
#   ts("  Loading WorldClim raster...")
#   clim_full <- terra::rast(paths$climate_raster)
# 
#   # Layer names vary by download method; match bio_1 and bio_12 by pattern
#   layer_names <- names(clim_full)
#   mat_idx <- grep("bio_?1$",  layer_names, value = FALSE)[1]
#   map_idx <- grep("bio_?12$", layer_names, value = FALSE)[1]
# 
#   if (is.na(mat_idx) || is.na(map_idx)) {
#     ts(sprintf("  Could not find bio_1 and bio_12 in layers: %s",
#                paste(layer_names, collapse = ", ")))
#     ts("  Skipping climate space analysis.")
#   } else {
#     clim_mat <- clim_full[[mat_idx]]
#     clim_map <- clim_full[[map_idx]]
# 
#     # Convert MAT from ×10 if values suggest ×10 encoding
#     mat_vals_raw <- terra::values(clim_mat)
#     if (!all(is.na(mat_vals_raw)) && max(abs(mat_vals_raw), na.rm = TRUE) > 100) {
#       clim_mat <- clim_mat / 10
#       ts("  MAT values divided by 10 (detected ×10 encoding)")
#     }
# 
#     # ---- Random climate sample per ecozone (balanced representation) --------
#     # Dissolve ecoregions to ecozones in WGS84, then spatSample n points each.
#     # This ensures each ecozone is equally represented in climate space
#     # regardless of its geographic area.
#     ts("  Sampling climate per ecozone...")
# 
#     n_per_ecozone <- 300L
# 
#     ecozones_agg <- ecoregions_clipped |>
#       dplyr::group_by(NAME_EN) |>
#       dplyr::summarise(.groups = "drop") |>
#       sf::st_make_valid() |>
#       sf::st_transform(4326)
# 
#     clim_stack        <- c(clim_mat, clim_map)
#     names(clim_stack) <- c("MAT", "MAP")
# 
#     grid_clim_df <- dplyr::bind_rows(lapply(unique_ecozones, function(ez) {
#       ez_vect  <- terra::vect(ecozones_agg[ecozones_agg$NAME_EN == ez, ])
#       clim_ez  <- terra::mask(terra::crop(clim_stack, ez_vect), ez_vect)
#       n_cells  <- sum(!is.na(terra::values(clim_ez[[1]])))
#       if (n_cells == 0L) return(NULL)
#       pts <- terra::spatSample(clim_ez,
#                                size   = min(n_per_ecozone, n_cells),
#                                method = "random",
#                                na.rm  = TRUE)
#       if (nrow(pts) == 0L) return(NULL)
#       # Access by position: terra may not preserve "MAT"/"MAP" layer names
#       # through crop/mask/arithmetic, so column names are not reliable
#       tibble::tibble(NAME_EN = ez, MAT = pts[[1]], MAP = pts[[2]])
#     })) |>
#       dplyr::filter(!is.na(MAT), !is.na(MAP))
# 
#     # Remove per-ecozone MAT outliers using Tukey fences (Q1/Q3 ± 1.5×IQR).
#     # Boundary cells at ecozone edges can carry erroneous or atypical values
#     # (e.g. a single Arctic Cordillera cell near the coastline showing 0 °C).
#     n_before <- nrow(grid_clim_df)
#     grid_clim_df <- grid_clim_df |>
#       dplyr::group_by(NAME_EN) |>
#       dplyr::filter(
#         MAT >= stats::quantile(MAT, 0.25) - 1.5 * stats::IQR(MAT),
#         MAT <= stats::quantile(MAT, 0.75) + 1.5 * stats::IQR(MAT)
#       ) |>
#       dplyr::ungroup()
#     n_removed <- n_before - nrow(grid_clim_df)
#     if (n_removed > 0L)
#       ts(sprintf("  Removed %d MAT outlier(s) via per-ecozone Tukey filter", n_removed))
# 
#     ts(sprintf("  Background climate points: %d across %d ecozones",
#                nrow(grid_clim_df), dplyr::n_distinct(grid_clim_df$NAME_EN)))
# 
#     # ---- Climate at actual EcM sampling locations ---------------------------
#     ts("  Extracting climate at sampling locations...")
# 
#     locs_wgs84 <- sf::st_transform(locs_src_sf, 4326)
#     locs_vect  <- terra::vect(locs_wgs84)
#     mat_samp   <- terra::extract(clim_mat, locs_vect)[, 2]
#     map_samp   <- terra::extract(clim_map, locs_vect)[, 2]
# 
#     samp_clim_df <- tibble::tibble(
#       source = locs_src_sf$source,
#       MAT    = mat_samp,
#       MAP    = map_samp
#     ) |>
#       dplyr::filter(!is.na(MAT), !is.na(MAP))
# 
#     ts(sprintf("  Sampling points with climate data: %d", nrow(samp_clim_df)))
# 
#     # ---- ecozone climate space plot -----------------------------------------------------
#     # Background uses colour aesthetic (ecozone_colors, same as sampling map).
#     # alpha = 0.2 matches the visual weight of the alpha = 0.2 ecozone fills on
#     # the sampling map, so both plots read as the same palette.
#     # Sampling points use a mapped shape aesthetic (no mapped colour) so that
#     # a second legend ("Sample source") can appear without conflicting with the
#     # ecozone colour scale.  override.aes in the shape guide injects the correct
#     # fill/colour/stroke so the legend symbols are identical to the sampling map.
#     p_climate_space <- ggplot2::ggplot() +
#       # Background: solid circles coloured by ecozone, muted to match map fills
#       ggplot2::geom_point(data = grid_clim_df,
#                           ggplot2::aes(x = MAT, y = MAP, colour = NAME_EN),
#                           shape = 19, size = 1.5, alpha = 0.2) +
#       # Sampling — GlobalFungi: open circle, black outline, white fill
#       ggplot2::geom_point(
#         data = dplyr::filter(samp_clim_df, source == "GlobalFungi"),
#         ggplot2::aes(x = MAT, y = MAP, shape = source),
#         fill = "white", colour = "black",
#         size = 3, stroke = 0.8, alpha = 0.7
#       ) +
#       # Sampling — GenBank: open triangle, red outline, white fill
#       ggplot2::geom_point(
#         data = dplyr::filter(samp_clim_df, source == "GenBank"),
#         ggplot2::aes(x = MAT, y = MAP, shape = source),
#         fill = "white", colour = "red",
#         size = 3, stroke = 0.8, alpha = 0.7
#       ) +
#       ggplot2::scale_colour_manual(
#         values = ecozone_colors,
#         name   = "Ecozone",
#         guide  = ggplot2::guide_legend(order = 1)
#       ) +
#       ggplot2::scale_shape_manual(
#         values = c("GlobalFungi" = 21L, "GenBank" = 24L),
#         name   = "Sample source",
#         guide  = ggplot2::guide_legend(
#           order        = 2,
#           override.aes = list(
#             fill   = "white",
#             colour = c("black", "red"),
#             size   = 3,
#             stroke = 0.8
#           )
#         )
#       ) +
#       ggplot2::labs(
#         x = "Mean Annual Temperature (\u00b0C)",
#         y = "Mean Annual Precipitation (mm)"
#       ) +
#       ggplot2::theme_bw() +
#       ggplot2::theme(
#         legend.position = "right",
#         legend.text     = ggplot2::element_text(size = 7),
#         legend.title    = ggplot2::element_text(size = 9, face = "bold")
#       )
# 
#     ggplot2::ggsave(file.path(paths$figures, "hutchinsonian_climate_space.png"),
#                     p_climate_space, width = 10, height = 7, dpi = 300)
#     ts("  Saved hutchinsonian_climate_space.png")
#   }
# }  # ----- END COMMENTED-OUT BLOCK -----
# }  # end if (FALSE)

# ---- Step 9: Ecozone-level sampling summary (Table S16b) -------------------
# Companion to the ecoregion-level summary in Step 6 (Table S16): of Canada's
# n_ecozones named ecozones (from ecoregions_clipped, Step 8), how many have
# *any* EcM fungal sampling, and how many meet the minimum-sample thresholds
# already used elsewhere in the SI (>=10 sites: Figure S8 coverage plot;
# >=30 sites: Figure S9 accumulation curves)?
#
# Three sample definitions are reported side by side, because "sample" is
# ambiguous across the SI:
#   1. GlobalFungi + GenBank combined, raw unique sampling locations
#      (locs_src_sf, Step 8 -- same data underlying Table S16 and the
#      sampling map).
#   2. GlobalFungi only, project-standard 3-decimal-binned "sites"
#      (add_site_id(), 00_setup.R) -- same scope/unit as the >=10 / >=30
#      thresholds already used for Figures S8 and S9.
#   3. GlobalFungi + GenBank combined, 3-decimal-binned "sites" -- same
#      combined scope as (1), but using the project's canonical site unit
#      instead of raw coordinates.

ts("Step 9: Ecozone-level sampling summary...")

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
  ts(sprintf("  [%s] total=%d  >=1: %d  >=10: %d  >=30: %d",
             r$definition, r$n_ecozones, r$ecozones_ge1, r$ecozones_ge10, r$ecozones_ge30))
}

readr::write_csv(ecozone_threshold_summary,
                 file.path(paths$out_hutchinsonian, "hutchinsonian_ecozone_summary.csv"))
ts("  Saved hutchinsonian_ecozone_summary.csv")

# ---- Step 10: Per-ecozone sample counts by source (Table S1) ----------------
# One row per named ecozone with raw unique locations and 3-decimal-binned sites,
# split by source. "Total" columns sum the two source columns (a location/site
# shared by both sources counts once per source). Same within-polygon +
# nearest-snap join as Step 10.
ts("Step 10: Per-ecozone sample counts by source...")

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

# ---- Step 11: Ecozone areas -> sampling density (Table S1) -------------------
# Area of each named ecozone within Canada, from the SAME unclipped, Albers
# ecoregion polygons used to assign the sample counts (so counts and areas share
# one polygon basis). Ecozone area = sum of its constituent ecoregion areas.
ts("Step 11: Computing ecozone areas for sampling density...")

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
ts(sprintf("  Saved hutchinsonian_ecozone_sample_counts.csv (%d ecozones; total locations = %d)",
           nrow(ecozone_sample_counts), sum(ecozone_sample_counts$total_locations)))

ts("17_hutchinsonian.R complete.")
