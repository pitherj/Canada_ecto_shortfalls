# =============================================================================
# 01_spatial_data.R  —  Download and process the base spatial layers
# =============================================================================
# PURPOSE
#   Obtain and prepare every spatial layer the rest of the pipeline needs:
#   national/administrative boundaries, Canadian ecoregions (with ecozone
#   names), and a Natural Earth Canada boundary plus major lakes for mapping.
#
# INPUTS
#   Downloaded automatically from public servers on first run:
#     * GADM administrative boundaries (Canada, level 1) via {geodata}
#     * Canadian ecoregions + ecozone-name table (agr.gc.ca)
#     * Natural Earth Canada boundary + lakes via {rnaturalearth}
#   Downloaded copies are written to data_raw/ (admin_boundaries/, ecoregions/,
#   natural_earth/). If those files are already present, download is skipped.
#
# OUTPUTS (data_derived/spatial/)
#   canada_simple.gpkg          dissolved, simplified Canada boundary (WGS84)
#   ecoregions_processed.gpkg   Canadian ecoregions, reprojected to WGS84
#   ecozone_names.csv           ecozone ID -> English/French name lookup
#   lakes_canada_albers.gpkg    major lakes clipped to Canada (Albers)
#
# NOTE
#   Every step checks whether its output already exists and skips if so, so the
#   script is safe to re-run. Delete an output file to force it to regenerate.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(geodata)
library(terra)
library(sf)
library(rnaturalearth)
library(foreign)

options(timeout = max(6000, getOption("timeout")))
sf::sf_use_s2(FALSE)   # GADM polygons have minor topology GEOS tolerates but s2 rejects

spatial_dir <- here::here("data_derived", "spatial")
dir.create(spatial_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1 — Download raw spatial layers (skipped if already in data_raw/)
# =============================================================================

# ---- 1a. Administrative boundaries (GADM level 1) ---------------------------
admin_dir <- file.path(paths$data_raw, "admin_boundaries")
dir.create(admin_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(paths$canada_bound_raw)) {
  temp <- file.path(paths$temp_dir, "gadm_download")
  dir.create(temp, showWarnings = FALSE, recursive = TRUE)
  terra::writeVector(geodata::gadm("CAN", 1, path = temp), paths$canada_bound_raw, overwrite = TRUE)
}

# ---- 1b. Ecoregions + ecozone names (Canadian National Ecological Framework) -
ecor_dir <- file.path(paths$data_raw, "ecoregions")
dir.create(ecor_dir, showWarnings = FALSE, recursive = TRUE)

if (!(length(list.files(ecor_dir, pattern = "\\.shp$", recursive = TRUE)) > 0)) {
  zip_path <- file.path(paths$temp_dir, "ecoregions.zip")
  download.file("https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
                destfile = zip_path, mode = "wb")
  unzip(zip_path, exdir = ecor_dir, overwrite = TRUE)
  shp_files <- list.files(ecor_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  if (length(shp_files) == 0) stop("No shapefile found after extracting ecoregions archive.")
  # small ecozone-name lookup table (DBF), same server
  ecoz_dbf <- file.path(dirname(shp_files[1]), "ecozone_names.dbf")
  if (!file.exists(ecoz_dbf))
    download.file("https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/zn_names.dbf",
                  destfile = ecoz_dbf, mode = "wb")
}

# ---- 1c. Natural Earth Canada boundary + lakes ------------------------------
ne_dir <- file.path(paths$data_raw, "natural_earth")
dir.create(ne_dir, showWarnings = FALSE, recursive = TRUE)

if (!(all(file.exists(paths$ne_canada_raw, paths$ne_lakes_raw)))) {
  canada_sf <- rnaturalearth::ne_countries(scale = "medium", country = "Canada", returnclass = "sf")
  terra::writeVector(terra::vect(canada_sf), paths$ne_canada_raw, overwrite = TRUE)
  lakes_sf <- rnaturalearth::ne_download(scale = "medium", type = "lakes",
                                         category = "physical", returnclass = "sf")
  terra::writeVector(terra::vect(lakes_sf), paths$ne_lakes_raw, overwrite = TRUE)
}

# =============================================================================
# PART 2 — Process layers into the derived forms the pipeline consumes
# =============================================================================

# ---- 2a. Dissolve + simplify the Canada boundary ----------------------------
# The provincial/territorial polygons are dissolved into a single national
# outline, then simplified (0.1-degree tolerance) to keep file size and
# plotting cost down. This boundary sets `coord_in_canada` in
# 04_combine_ecm_dataset.R and is the basemap for every map in the pipeline.

if (!(file.exists(paths$canada_bound))) {
  canada <- terra::vect(paths$canada_bound_raw)

  canada_union  <- terra::aggregate(canada, dissolve = TRUE)
  canada_simple <- terra::simplifyGeom(canada_union, tolerance = 0.1)

  terra::writeVector(canada_simple, paths$canada_bound, overwrite = TRUE)
}

# ---- 2b. Ecoregions (reproject to WGS84) + ecozone name lookup --------------
if (!(file.exists(paths$ecoregions_processed))) {
  shp_files <- list.files(ecor_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  if (length(shp_files) == 0) stop("Ecoregions shapefile not found. Re-run Part 1.")
  ecor_shp <- shp_files[1]

  # Read via sf to preserve French accents in region names
  ecoregions <- sf::st_read(ecor_shp, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  sf::st_write(ecoregions, paths$ecoregions_processed, delete_dsn = TRUE, quiet = TRUE)

  ecoz_dbf <- file.path(dirname(ecor_shp), "ecozone_names.dbf")
  if (file.exists(ecoz_dbf) && !file.exists(paths$ecozone_names)) {
    ecozone_names <- foreign::read.dbf(ecoz_dbf, as.is = TRUE)
    write.csv(ecozone_names, paths$ecozone_names, row.names = FALSE)
  }
}

# ---- 2c. Natural Earth lakes -> Albers, clipped to Canada -------------------
# The Natural Earth country outline is used only as the clipping mask here, so
# it is projected in memory and not written out.

if (!file.exists(paths$lakes_albers)) {
  canada_ne <- terra::vect(paths$ne_canada_raw)
  lakes_ne  <- terra::vect(paths$ne_lakes_raw)

  canada_albers_v <- terra::project(canada_ne, crs_albers)
  lakes_albers_v  <- terra::project(lakes_ne,  crs_albers)
  lakes_canada    <- terra::intersect(lakes_albers_v, canada_albers_v)

  terra::writeVector(lakes_canada, paths$lakes_albers, overwrite = TRUE)
}

