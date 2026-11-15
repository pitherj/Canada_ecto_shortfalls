# =============================================================================
# 12_wallacean_density_map.R  —  Global GlobalFungi EcM sampling density
# =============================================================================
# PURPOSE
#   Map the number of GlobalFungi samples with at least one EcM SH code, per
#   100 x 100 km grid cell, worldwide (Mollweide equal-area projection).
#   EcM samples are identified by SH presence:
#     1. UNITE taxonomy for every SH code (paths$unite_taxonomy).
#     2. EcM SH codes flagged with the same definition as the Canadian
#        dataset (classify_ecm() in 00_setup.R: FungalTraits EcM genus, or
#        the two-thirds family rule).
#     3. Only those columns are read from the 13 GB abundance matrix.
#     4. A sample is kept if it has at least one EcM read.
#     5. The same quality filters as the Canadian dataset (02): ITS2 or
#        ITSboth; not manipulated; sample_type in GF_SAMPLE_TYPES_KEPT (soil,
#        root); valid coordinates.
#
# INPUTS
#   paths$unite_taxonomy (02); paths$fungaltraits; GlobalFungi v5 metadata and
#   SH abundance matrix (data_raw)
#
# OUTPUTS
#   figures/Figure-S2_gf_sampling_density_world.png
#   data_derived/checkpoints/gf_global_ecm_sample_ids.csv  sample IDs and
#       coordinates passing all filters (also read by 13)
#
# MANUSCRIPT ITEMS IT FEEDS
#   Figure S2; the global sampling-location count in Table 1 (via 13).
#
# RE-RUN WHEN
#   GlobalFungi, the UNITE build, FungalTraits, the EcM rules or the sample
#   filters change. Delete the checkpoint and the figure to rebuild. Step 3
#   reads the 13 GB matrix: 5-15 minutes.
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Packages ----------------------------------------------------------------
# data.table  : fast large-file reading (fread with column selection)
# sf          : vector geometry (points, grid, spatial joins)
# ggplot2     : plotting (loaded via 00_setup.R, re-stated for clarity)
# rnaturalearth / rnaturalearthdata : land polygons for map backgrounds
# scales      : log-scale colour formatting

library(data.table)

# ---- Output paths ------------------------------------------------------------

out_sample_ids <- here::here("data_derived", "checkpoints", "gf_global_ecm_sample_ids.csv")
out_fig_world  <- paths$fig_density_world

# ---- Projection strings ------------------------------------------------------
# World map: Mollweide equal-area (PROJ.4)
crs_moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

# Canada map: Albers Equal Area Conic (defined in 00_setup.R as crs_albers)

# ---- Grid cell size ----------------------------------------------------------
CELL_M <- 100000L   # 100 km in metres (both projections use metres)

# =============================================================================
# STEP 1: Identify EcM SH codes via UNITE taxonomy × FungalTraits
# =============================================================================
# Skip the expensive Step 2 (reading 13 GB file) if the checkpoint already
# exists. Both steps are guarded independently.

# 1a. Load UNITE SH taxonomy
# Columns: sh_code, kingdom, phylum, class, order, family, genus, species
unite_tax <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)

# 1b. Which UNITE SH codes are EcM?
#     The genus column in unite_sh_taxonomy.csv stores the genus as extracted
#     from the UNITE taxonomy string (e.g., "g__Russula" or just "Russula").
#     Strip the "g__" prefix if present, then lowercase for matching.
unite_tax <- dplyr::mutate(
  unite_tax,
  genus_clean = tolower(sub("^g__", "", genus))
)

# EcM SH codes: same definition as the Canadian dataset (classify_ecm() in
# 00_setup.R: FungalTraits EcM genus, or the two-thirds family rule).
ecm_ref <- ecm_reference()
ecm_sh_codes <- unite_tax |>
  dplyr::filter(classify_ecm(genus, family, ecm_ref) %in% c("genus", "family")) |>
  dplyr::pull(sh_code) |>
  unique()

if (length(ecm_sh_codes) == 0L) {
  stop(
    "No EcM SH codes identified. ",
    "Check that genus names in unite_sh_taxonomy.csv match FungalTraits GENUS column."
  )
}

# =============================================================================
# STEP 2: Read EcM SH columns from the 13 GB abundance matrix and identify
#         samples with at least one EcM SH present (row sum > 0)
# =============================================================================

if (file.exists(out_sample_ids)) {

  ecm_samples <- readr::read_csv(out_sample_ids, show_col_types = FALSE)

} else {

  # 2a. Read the file header to find which of our EcM SH codes are present
  #     as columns in the abundance matrix. This avoids requesting columns
  #     that don't exist, which would cause fread() to error.
  gf_header_cols <- names(
    data.table::fread(paths$gf_sh_abundance, sep = "\t", quote = "", nrows = 0L)
  )

  # The first column is sample_ID; the rest are SH codes.
  # Match our EcM SH codes against the matrix columns (exact match first).
  ecm_cols_present <- intersect(ecm_sh_codes, gf_header_cols)

  # If exact matches are few, try prefix matching (strips UNITE version suffix,
  # e.g., "SH1052460.10FU" → "SH1052460") — handles version mismatches.
  if (length(ecm_cols_present) < 0.1 * length(ecm_sh_codes)) {
    ecm_prefix   <- sub("\\.[0-9]+FU$", "", ecm_sh_codes)
    col_prefix   <- sub("\\.[0-9]+FU$", "", gf_header_cols)
    matched_idx  <- match(ecm_prefix, col_prefix)
    ecm_cols_present <- gf_header_cols[matched_idx[!is.na(matched_idx)]]
  }

  if (length(ecm_cols_present) == 0L) {
    stop(
      "No EcM SH code columns found in the abundance matrix. ",
      "Check that paths$gf_sh_abundance points to the correct file and that ",
      "the UNITE taxonomy version matches the SH codes in the matrix."
    )
  }

  # 2b/2c. Identify samples with >= 1 EcM SH detection WITHOUT loading the matrix
  #        into R. Reading ~13,000 columns as a wide matrix overruns R's 2^31-byte
  #        string limit (and would be a ~2 GB object). We only need a per-row
  #        presence test, so we stream the matrix once through awk: for each
  #        sample row, print sample_ID if ANY EcM column has a read count > 0.
  #        The result is just a short list of sample IDs.
  scol    <- match("sample_ID", gf_header_cols)          # sample_ID column position
  ecm_idx <- match(ecm_cols_present, gf_header_cols)      # EcM column positions
  ids_tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(ids_tmp), add = TRUE)
  awk_prog <- paste0(
    'BEGIN{FS="\\t"; n=split(cols,a,",")} ',
    'NR==1{next} ',                                       # skip header row
    '{for(i=1;i<=n;i++){ if($(a[i])+0>0){ print $scol; break } }}'
  )
  awk_cmd <- sprintf("awk -v scol=%d -v cols=%s %s %s > %s",
                     scol, shQuote(paste(ecm_idx, collapse = ",")),
                     shQuote(awk_prog), shQuote(paths$gf_sh_abundance),
                     shQuote(ids_tmp))
  if (system(awk_cmd) != 0L || !file.exists(ids_tmp))
    stop("awk EcM-detection scan of the GF SH abundance matrix failed.")
  ecm_sample_ids <- readLines(ids_tmp)
  unlink(ids_tmp)

  # 2d. Load sample metadata; apply quality filters; retain only EcM samples
  gf_meta <- data.table::fread(
    paths$gf_metadata,
    sep    = "\t",
    quote  = "",
    select = c("sample_ID", "latitude", "longitude",
               "barcoding_region", "manipulated", "sample_type")
  )

  # Apply the same quality filters as the Canadian pipeline (02_globalfungi.R)
  gf_meta_filt <- gf_meta[
    barcoding_region %in% c("ITS2", "ITSboth") &
    manipulated == "NO" &
    sample_type %in% GF_SAMPLE_TYPES_KEPT &
    !is.na(latitude)  & !is.na(longitude) &
    latitude  >= -90  & latitude  <= 90 &
    longitude >= -180 & longitude <= 180
  ]

  # Retain only EcM-confirmed samples
  ecm_samples <- gf_meta_filt[sample_ID %in% ecm_sample_ids,
                               .(sample_ID, latitude, longitude)]

  # Save checkpoint
  readr::write_csv(as.data.frame(ecm_samples), out_sample_ids)

}

# =============================================================================
# STEP 3: World sampling density map (Mollweide, 100 km grid)
# =============================================================================

if (!file.exists(out_fig_world)) {

  # 3a. Convert EcM sample points to sf, reproject to Mollweide
  pts_world <- sf::st_as_sf(
    ecm_samples,
    coords = c("longitude", "latitude"),
    crs    = crs_wgs84
  ) |>
    sf::st_transform(crs_moll)

  # 3b. Build a global 100 km × 100 km grid covering the Mollweide extent.
  #     The Mollweide world bounds are approximately ±18,000 km W–E and
  #     ±9,000 km N–S. We add a small buffer and let st_make_grid handle edges.
  world_bbox <- sf::st_bbox(
    c(xmin = -18040096, xmax = 18040096,
      ymin = -9020048,  ymax = 9020048),
    crs = sf::st_crs(crs_moll)
  )
  grid_world <- sf::st_make_grid(
    sf::st_as_sfc(world_bbox),
    cellsize = CELL_M,
    what     = "polygons"
  ) |> sf::st_sf()

  # 3c. Count samples per grid cell via spatial join.
  #     Assign a cell ID to the grid first, then join points → grid cells,
  #     then count by cell ID. This correctly handles multiple points per cell.
  grid_world$cell_id <- seq_len(nrow(grid_world))

  pts_with_cell_world <- sf::st_join(
    pts_world,
    grid_world[, "cell_id"],
    join = sf::st_within
  )

  cell_counts_world <- pts_with_cell_world |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(cell_id)) |>
    dplyr::count(cell_id, name = "n_samples")

  grid_world <- dplyr::left_join(grid_world, cell_counts_world, by = "cell_id") |>
    dplyr::mutate(n_samples = dplyr::coalesce(n_samples, 0L))

  # Keep only cells with at least 1 sample for plotting
  grid_filled_world <- dplyr::filter(grid_world, n_samples > 0L)

  # 3d. Land outline for background (low-resolution)
  land_world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf") |>
    sf::st_transform(crs_moll)

  # 3e. Plot
  p_world <- ggplot2::ggplot() +
    ggplot2::geom_sf(data  = land_world,
                     fill  = "grey85",
                     colour = "grey60",
                     linewidth = 0.15) +
    ggplot2::geom_sf(data = grid_filled_world,
                     ggplot2::aes(fill = n_samples),
                     colour = NA) +
    ggplot2::scale_fill_viridis_c(
      option    = "plasma",
      trans     = "log10",
      name      = "No. samples",
      labels    = scales::label_comma(),
      guide     = ggplot2::guide_colourbar(
        barwidth  = 0.8,
        barheight = 8,
        title.position = "top"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = c(0.95, 0.75)
    )

  save_fig_formats(
    out_fig_world,
    p_world,
    width    = 14, height = 7, dpi = 300, units = "in"
  )

}

# The Canada-only sampling density map is not used in the manuscript, so it is
# not produced here. This script's sole figure output is the global density map
# (Figure S2) built in Step 3 above.

