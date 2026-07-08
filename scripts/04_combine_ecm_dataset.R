# =============================================================================
# Combine Sources and Filter for EcM Fungi
# =============================================================================
# Step 1: Stack GlobalFungi and GenBank long-format tables, harmonise columns,
#         validate ALL coordinates against the GADM Canada boundary (nullifying
#         out-of-Canada lat/lon so downstream scripts need no spatial filtering).
#         Output: data_derived/emf_canada_combined.csv
#
# Step 2: Filter combined dataset to ectomycorrhizal taxa only using
#         FungalTraits v1.2 genus-level primary_lifestyle (pinned reference).
#         Output: data_derived/emf_canada_em_only.csv  (primary dataset for all
#                 downstream analyses)
#
# Prerequisites:
#   02_globalfungi.R → data_derived/globalfungi_canada_long.csv
#   03_genbank.R     → data_derived/genbank_emf_canada_long.csv
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(sf)

# ---- Step 1: Combine GlobalFungi and GenBank --------------------------------

if (file.exists(paths$emf_combined)) {
  ts("Step 1: emf_canada_combined.csv already exists — skipping.")
} else {

  ts("Step 1: Combining GlobalFungi and GenBank outputs...")

  gf <- readr::read_csv(paths$gf_long_out, show_col_types = FALSE)
  ts(sprintf("  GlobalFungi: %d rows, %d unique SHs, %d unique samples",
             nrow(gf), dplyr::n_distinct(gf$sh_code, na.rm = TRUE), dplyr::n_distinct(gf$sample_ID)))

  gb <- readr::read_csv(paths$gb_long_out, show_col_types = FALSE)
  ts(sprintf("  GenBank:     %d rows, %d unique SHs, %d unique accessions",
             nrow(gb), dplyr::n_distinct(gb$sh_code, na.rm = TRUE), dplyr::n_distinct(gb$accession)))

  # Harmonise lat/lon column names
  if (all(c("latitude", "longitude") %in% names(gf)))
    gf <- dplyr::rename(gf, lat = latitude, lon = longitude)

  # Parse GenBank lat_lon_gb string → numeric lat/lon
  # parse_gb_latlon() is defined in 00_setup.R and operates on a single value;
  # apply it row-wise here to build lat/lon columns.
  gb_ll <- dplyr::bind_rows(lapply(gb$lat_lon_gb, parse_gb_latlon))
  gb    <- dplyr::mutate(gb, lat = gb_ll$lat, lon = gb_ll$lon)

  gf <- dplyr::mutate(gf, source = "GlobalFungi")

  combined <- dplyr::bind_rows(gf, gb)
  ts(sprintf("  Combined (pre-spatial filter): %d rows total", nrow(combined)))

  # SH overlap summary
  ts(sprintf("  SH overlap — GF: %d, GB: %d, shared: %d, either: %d",
             dplyr::n_distinct(gf$sh_code, na.rm = TRUE),
             dplyr::n_distinct(gb$sh_code, na.rm = TRUE),
             length(intersect(unique(gf$sh_code), unique(gb$sh_code))),
             dplyr::n_distinct(combined$sh_code, na.rm = TRUE)))

  # Spatial containment check — applied to ALL records (GlobalFungi + GenBank)
  # with parsed coordinates. A logical column `coord_in_canada` is added:
  #   TRUE  — coordinates present and within the (buffered) GADM Canada boundary
  #   FALSE — coordinates present but outside the buffered boundary
  #   NA    — no coordinates in the source data
  #
  # Original lat/lon values are preserved in all cases. Downstream spatial
  # analyses filter on coord_in_canada == TRUE rather than !is.na(lat), so
  # records with out-of-boundary coordinates are excluded from spatial work
  # but retain their coordinates for enumeration and traceability.
  #
  # The boundary is buffered by 2 km (in the equal-area CRS `crs_albers`, so
  # the buffer distance is a true 2,000 m regardless of latitude) before the
  # containment test. This absorbs GADM coastline simplification / topology
  # artefacts: a handful of GlobalFungi records with country == "Canada"
  # metadata but coordinates on the coast or on small islands (e.g. Gulf of
  # St. Lawrence) fell just outside the unbuffered polygon, which silently
  # dropped them from every downstream analysis gated on
  # coord_in_canada == TRUE (this caused the Table 1 vs Table 3 discrepancy
  # in the manuscript: 2 named species / 4 SH codes; see Supp. Methods,
  # Decisions Log E2/E3). Any records still flagged coord_in_canada == FALSE
  # after the 2 km buffer are far enough outside the boundary that their
  # country == "Canada" metadata (if present) should be treated as
  # authoritative only for non-spatial enumeration, not spatial analyses.
  # For GenBank records, canada_basis is also updated to
  # "coordinates_outside_canada" for traceability.
  ts("  Checking coordinate containment within Canada boundary (GADM, +2 km buffer)...")
  canada_bound_sf <- sf::st_read(paths$canada_bound, quiet = TRUE)
  canada_bound_buffered <- canada_bound_sf |>
    sf::st_transform(crs_albers) |>
    sf::st_buffer(dist = 2000) |>   # 2 km, in metres (crs_albers is equal-area, units = m)
    sf::st_transform(4326)
  coord_idx <- which(!is.na(combined$lat) & !is.na(combined$lon))
  ts(sprintf("  Records with coordinates: %d", length(coord_idx)))

  # Initialise flag: NA for records with no coordinates
  combined$coord_in_canada <- NA

  if (length(coord_idx) > 0) {
    pts <- sf::st_as_sf(combined[coord_idx, ], coords = c("lon", "lat"), crs = 4326)
    sf::sf_use_s2(FALSE)
    in_can <- lengths(sf::st_intersects(pts, canada_bound_buffered)) > 0
    sf::sf_use_s2(TRUE)
    n_outside <- sum(!in_can)
    ts(sprintf("  In Canada: %d | outside boundary (coord_in_canada = FALSE): %d",
               sum(in_can), n_outside))
    combined$coord_in_canada[coord_idx] <- in_can
    # Update canada_basis for GenBank records whose coordinates fall outside Canada
    if (n_outside > 0) {
      outside_idx <- coord_idx[!in_can]
      gb_outside <- outside_idx[combined$source[outside_idx] == "GenBank"]
      if (length(gb_outside) > 0)
        combined$canada_basis[gb_outside] <- "coordinates_outside_canada"
    }
  }

  readr::write_csv(combined, paths$emf_combined)
  ts("  Saved emf_canada_combined.csv")
}

# ---- Step 2: Filter for EcM using FungalTraits ------------------------------

if (file.exists(paths$emf_data)) {
  ts("Step 2: emf_canada_em_only.csv already exists — skipping.")
} else {

  ts("Step 2: Filtering for ectomycorrhizal taxa using FungalTraits...")

  # FungalTraits v1.2 (Põlme et al. 2020): pinned reference file (see 00_setup.R
  # and data_raw/fungaltraits/fungaltraits_version.txt for provenance).
  ft <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
    dplyr::rename_with(tolower)

  em_genera <- ft |>
    dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
    dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
    dplyr::distinct(genus_lower)
  ts(sprintf("  EcM genera in FungalTraits: %d", nrow(em_genera)))

  ft_join <- ft |>
    dplyr::select(genus, primary_lifestyle, secondary_lifestyle,
                  dplyr::matches("lineage_template"),
                  dplyr::matches("exploration_type_template")) |>
    dplyr::rename_with(~ sub("_template$", "", .x)) |>
    dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
    dplyr::select(-genus) |>
    dplyr::distinct(genus_lower, .keep_all = TRUE)

  combined <- readr::read_csv(paths$emf_combined, show_col_types = FALSE) |>
    dplyr::mutate(genus_lower = tolower(trimws(genus)))
  combined  <- dplyr::left_join(combined, ft_join, by = "genus_lower")

  em_only <- combined |>
    dplyr::filter(genus_lower %in% em_genera$genus_lower) |>
    dplyr::select(-genus_lower)

  ts(sprintf("  EM records: %d | SHs: %d | genera: %d",
             nrow(em_only), dplyr::n_distinct(em_only$sh_code, na.rm = TRUE),
             dplyr::n_distinct(em_only$genus)))
  ts("  By source:"); print(dplyr::count(em_only, source))

  readr::write_csv(em_only, paths$emf_data)
  ts("  Saved emf_canada_em_only.csv")
}

ts("04_combine_ecm_dataset.R complete.")
