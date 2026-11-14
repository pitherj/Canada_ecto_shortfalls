# =============================================================================
# Combine Sources and Filter for EcM Fungi
# =============================================================================
# Step 1: Stack GlobalFungi and GenBank long-format tables, harmonise columns,
#         check ALL coordinates against the Canada boundary (coord_in_canada),
#         and keep only GenBank records whose own location evidence places
#         them in Canada (the GenBank location rule, 1b).
#         Outputs: data_derived/emf_canada_combined.csv
#                  data_derived/genbank_excluded_non_canadian.csv (dropped
#                  GenBank records, with the reason)
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

if (!(file.exists(paths$emf_combined))) {

  gf <- readr::read_csv(paths$gf_long_out, show_col_types = FALSE)

  gb <- readr::read_csv(paths$gb_long_out, show_col_types = FALSE)

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

  # ---- 1a. Coordinate check: is each point inside Canada? -------------------
  # Applied to ALL records (GlobalFungi + GenBank) that have coordinates. Adds
  # a logical column `coord_in_canada`:
  #   TRUE  - coordinates inside the Canada boundary (buffered by 2 km)
  #   FALSE - coordinates outside it
  #   NA    - no coordinates in the source data
  # Coordinates are always kept, so records stay traceable. Spatial analyses
  # downstream use only records with coord_in_canada == TRUE.
  #
  # The boundary (paths$canada_bound, built in 01_spatial_data.R) is simplified
  # at a 1 km tolerance. Buffering it by 2 km (in the equal-area projection
  # `crs_albers`, so 2 km everywhere) absorbs that simplification and small
  # coordinate imprecision at coastal and island sites.
  canada_bound_buffered <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
    sf::st_transform(crs_albers) |>
    sf::st_buffer(dist = 2000) |>   # 2 km, in metres
    sf::st_transform(4326)
  coord_idx <- which(!is.na(combined$lat) & !is.na(combined$lon))
  combined$coord_in_canada <- NA
  if (length(coord_idx) > 0) {
    pts <- sf::st_as_sf(combined[coord_idx, ], coords = c("lon", "lat"), crs = 4326)
    sf::sf_use_s2(FALSE)
    combined$coord_in_canada[coord_idx] <-
      lengths(sf::st_intersects(pts, canada_bound_buffered)) > 0
    sf::sf_use_s2(TRUE)
  }

  # ---- 1b. GenBank location rule: which GenBank records are Canadian? -------
  # The GenBank search matches "Canada" anywhere in a record (03_genbank.R), so
  # each record's own location evidence decides whether it is kept:
  #   - its location field (country_gb, GenBank's /geo_loc_name) begins with
  #     "Canada"                                           -> kept;
  #   - its location field names another country           -> dropped;
  #   - it has no location field, and its coordinates are
  #     inside Canada                                      -> kept;
  #   - it has no location field, and its coordinates are
  #     outside Canada                                     -> dropped;
  #   - it has neither a location field nor coordinates    -> dropped.
  # Where the field and the coordinates disagree, the field decides. Such a
  # record is kept, but its coordinates fail the check in 1a, so it is left out
  # of spatial analyses. GlobalFungi records were already restricted to
  # Canada by their own country field in 02_globalfungi.R.
  #
  # `location_evidence` records which case applies to each GenBank record.
  # Dropped records are written, with the reason, to paths$genbank_excluded,
  # so the exclusion can be audited.
  combined <- combined |>
    dplyr::select(-dplyr::any_of("canada_basis")) |>   # older label from 03
    dplyr::mutate(location_evidence = dplyr::case_when(
      source != "GenBank"                                 ~ NA_character_,
      grepl("^Canada", country_gb)                        ~ "location field: Canada",
      !is.na(country_gb)                                  ~ "location field: another country",
      coord_in_canada %in% TRUE                           ~ "no location field; coordinates inside Canada",
      coord_in_canada %in% FALSE                          ~ "no location field; coordinates outside Canada",
      TRUE                                                ~ "no location field or coordinates"))

  kept_evidence <- c("location field: Canada",
                     "no location field; coordinates inside Canada")
  is_excluded <- combined$source == "GenBank" &
                 !combined$location_evidence %in% kept_evidence

  combined |>
    dplyr::filter(is_excluded) |>
    dplyr::select(accession, organism, sh_code, genus, species, country_gb,
                  lat_lon_gb, coord_in_canada, location_evidence) |>
    readr::write_csv(paths$genbank_excluded)

  combined <- dplyr::filter(combined, !is_excluded)

  readr::write_csv(combined, paths$emf_combined)
}

# ---- Step 2: Filter for EcM using FungalTraits ------------------------------

if (!(file.exists(paths$emf_data))) {

  # FungalTraits v1.2 (Põlme et al. 2020): pinned reference file (see 00_setup.R
  # and data_raw/fungaltraits/fungaltraits_version.txt for provenance).
  ft <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
    dplyr::rename_with(tolower)

  em_genera <- ft |>
    dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
    dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
    dplyr::distinct(genus_lower)

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

  readr::write_csv(em_only, paths$emf_data)
}

