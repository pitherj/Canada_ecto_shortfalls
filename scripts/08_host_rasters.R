# =============================================================================
# Build BIEN2 EcM Host Rasters and Ecoregion Habitat Layer
# =============================================================================
# Loads BIEN2 modelled range shapefiles for all EcM host plant species native
# to Canada (from 06_host_species.R), rasterizes them at 0.5°
# resolution, and identifies which ecoregions contain EcM host habitat.
#
# Range shapefiles must already be present in data_raw/bien2_ranges/ — run
# scripts/07_bien2_ranges.R first.
#
# Workflow:
#   1.  Load host species list (output of 06_host_species.R)
#   2.  Load BIEN2 range shapefiles from data_raw/bien2_ranges/
#   3.  Filter ranges to those overlapping Canada
#   4.  Rasterize: count host species richness per 0.5° cell
#   5.  Identify ecoregions with EcM host habitat (spatial join)
#
# Checkpoint files (data_derived/checkpoints/):
#   bien2_ecm_host_ranges.gpkg      — BIEN2 range polygons (filtered to Canada)
#
# Outputs:
#   data_derived/spatial/bien_host_richness_0.5deg.tif    — host species richness
#   data_derived/spatial/bien_host_species_stack.tif      — per-species binary raster
#                                                      stack (one layer per host
#                                                      species; named by species)
#   data_derived/spatial/bien_ecoregions_with_host_habitat.gpkg — ecoregion polygons
#                                                      with EcM host presence
#
# Runtime notes:
#   - Step 2 (shapefile loading) ~1–5 min depending on species count
#   - Step 4 (rasterization loop) ~5–20 min
#   Both steps are checkpointed; delete checkpoint files to force re-run.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)
library(terra)

sf::sf_use_s2(FALSE)  # Use GEOS rather than S2 for intersection operations

# ---- Step 1: Load host species list ------------------------------------------

ts("Step 1: Loading EcM host species list...")
if (!file.exists(paths$host_species)) {
  stop("Host species file not found: ", paths$host_species,
       "\nRun 06_host_species.R first.")
}
host_tbl        <- readr::read_csv(paths$host_species, show_col_types = FALSE)
em_canada_species <- host_tbl$species
ts(sprintf("  Host species: %d", length(em_canada_species)))

# ---- Step 2: Load BIEN2 range shapefiles (checkpointed) ---------------------

if (file.exists(paths$bien_ranges)) {
  ts("Step 2: Loading checkpointed BIEN2 range polygons...")
  ranges_canada <- sf::st_read(paths$bien_ranges, quiet = TRUE)
  ts(sprintf("  Range polygons loaded: %d features (%d unique species)",
             nrow(ranges_canada),
             dplyr::n_distinct(ranges_canada$species)))
} else {

  if (!dir.exists(paths$bien2_ranges_dir)) {
    stop("BIEN2 ranges directory not found: ", paths$bien2_ranges_dir,
         "\nRun scripts/07_bien2_ranges.R first.")
  }

  ts(sprintf("Step 2: Loading BIEN2 shapefiles for %d species from disk...",
             length(em_canada_species)))

  range_list <- lapply(em_canada_species, function(sp) {
    sp_us  <- gsub(" ", "_", sp)
    sp_dir <- file.path(paths$bien2_ranges_dir, sp_us)
    shp    <- list.files(sp_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp) == 0) return(NULL)
    tryCatch(
      sf::st_read(shp[[1]], quiet = TRUE) |>
        sf::st_make_valid() |>
        sf::st_transform("EPSG:4326") |>
        dplyr::transmute(species = sp),   # retain only species + geometry
      error = function(e) {
        message(sprintf("  Skipping %s: %s", sp, conditionMessage(e)))
        NULL
      }
    )
  })

  n_found <- sum(!vapply(range_list, is.null, logical(1)))
  ts(sprintf("  Shapefiles read: %d / %d species", n_found, length(em_canada_species)))
  if (n_found < length(em_canada_species))
    message(sprintf("  (%d species have no shapefile — not available from biendata.org)",
                    length(em_canada_species) - n_found))

  ranges_raw <- dplyr::bind_rows(Filter(Negate(is.null), range_list))

  # Filter ranges to those overlapping Canada
  ts("Step 3: Filtering ranges to Canada boundary...")
  canada_wgs84 <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
    sf::st_transform(4326)

  ranges_canada <- sf::st_filter(ranges_raw, canada_wgs84)
  ts(sprintf("  Range polygons overlapping Canada: %d features (%d unique species)",
             nrow(ranges_canada),
             dplyr::n_distinct(ranges_canada$species)))

  sf::st_write(ranges_canada, paths$bien_ranges,
               delete_dsn = TRUE, quiet = TRUE)
  ts(sprintf("  Saved range checkpoint -> %s", basename(paths$bien_ranges)))
}

# ---- Spatial setup (unconditional) ------------------------------------------

canada_wgs84 <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
  sf::st_transform(4326)

canada_bbox <- sf::st_bbox(canada_wgs84)

# 0.5° raster template covering Canada bounding box
rast_template <- terra::rast(
  xmin = floor(canada_bbox["xmin"]),
  xmax = ceiling(canada_bbox["xmax"]),
  ymin = floor(canada_bbox["ymin"]),
  ymax = ceiling(canada_bbox["ymax"]),
  resolution = 0.5,
  crs = "EPSG:4326"
)

# Canada mask raster (for masking output)
canada_vect <- terra::vect(canada_wgs84)

# Vectorise ranges once for rasterization
ranges_vect <- terra::vect(ranges_canada)

# ---- Step 4: Rasterize host species richness + per-species stack (checkpointed)

# Both outputs are always written together; require both to skip re-run.
if (file.exists(paths$bien_richness) && file.exists(paths$bien_species_stack)) {
  ts("Step 4: Loading checkpointed host richness raster and species stack...")
  richness_raster <- terra::rast(paths$bien_richness)
  species_stack   <- terra::rast(paths$bien_species_stack)
  ts(sprintf("  Richness raster: %d x %d cells  |  Stack layers: %d",
             nrow(richness_raster), ncol(richness_raster),
             terra::nlyr(species_stack)))
} else {

  sp_list <- unique(ranges_canada$species)
  n_sp    <- length(sp_list)
  ts(sprintf("Step 4: Rasterizing %d host species into per-species stack...", n_sp))

  # Collect individual binary rasters in a list, then assemble stack once.
  # This avoids a running-sum pattern that discards per-species information.
  sp_layers <- vector("list", n_sp)

  for (i in seq_along(sp_list)) {
    sp_poly      <- ranges_vect[ranges_vect$species == sp_list[i], ]
    sp_rast      <- terra::rasterize(sp_poly, rast_template, field = 1)
    sp_layers[[i]] <- terra::ifel(is.na(sp_rast), 0L, 1L)
    if (i %% 20 == 0 || i == n_sp)
      ts(sprintf("    Processed %d / %d species", i, n_sp))
  }

  # Assemble stack; layer names = species names for downstream filtering.
  # Normalise to spaces ("Abies amabilis") to match host_tbl and all other
  # species lists in this project (BIEN range$species uses underscores).
  species_stack <- terra::rast(sp_layers)
  names(species_stack) <- gsub("_", " ", sp_list)

  # Mask outside Canada (sets all layers to NA for out-of-Canada cells)
  species_stack <- terra::mask(species_stack, canada_vect)

  # Richness = sum across layers (na.rm = FALSE → outside Canada stays NA)
  richness_raster <- terra::app(species_stack, fun = "sum", na.rm = FALSE)
  richness_raster <- terra::ifel(richness_raster == 0, NA, richness_raster)
  names(richness_raster) <- "host_richness"

  terra::writeRaster(richness_raster, paths$bien_richness,      overwrite = TRUE)
  terra::writeRaster(species_stack,   paths$bien_species_stack, overwrite = TRUE)
  ts(sprintf("  Saved richness raster -> %s", basename(paths$bien_richness)))
  ts(sprintf("  Saved species stack   -> %s", basename(paths$bien_species_stack)))

  # Summary
  rich_vals <- terra::values(richness_raster)
  rich_vals <- rich_vals[!is.na(rich_vals)]
  ts(sprintf("  Cells with host habitat: %d  |  median richness: %.0f  |  max: %d",
             length(rich_vals), median(rich_vals), max(rich_vals)))
}

# ---- Step 5: Ecoregion habitat layer -----------------------------------------

ts("Step 5: Building ecoregion host habitat layer...")

if (!file.exists(paths$ecoregions_processed)) {
  warning("Ecoregion file not found: ", paths$ecoregions_processed,
          "\nSkipping ecoregion habitat layer.")
} else {
  ecoregions <- sf::st_read(paths$ecoregions_processed, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)

  ts(sprintf("  Ecoregions loaded: %d", nrow(ecoregions)))

  # Determine which ecoregions intersect any host range polygon
  # Use st_intersects for a sparse logical matrix
  intersects_mat <- sf::st_intersects(ecoregions, ranges_canada,
                                       sparse = TRUE)
  ecoregions$has_ecm_host_habitat <- lengths(intersects_mat) > 0

  n_habitat_ecor  <- sum(ecoregions$has_ecm_host_habitat)
  n_total_ecor    <- nrow(ecoregions)
  ts(sprintf("  Ecoregions with EcM host habitat: %d / %d (%.1f%%)",
             n_habitat_ecor, n_total_ecor,
             100 * n_habitat_ecor / n_total_ecor))

  sf::st_write(ecoregions, paths$bien_ecoregions,
               delete_dsn = TRUE, quiet = TRUE)
  ts(sprintf("  Saved -> %s", basename(paths$bien_ecoregions)))
}

ts("08_host_rasters.R complete.")
