# =============================================================================
# 04_combine_ecm_dataset.R  —  Combine sources; keep Canadian EcM records
# =============================================================================
# PURPOSE
#   Step 1: stack the GlobalFungi and GenBank records, check every coordinate
#           against the Canada boundary (coord_in_canada), and keep only
#           GenBank records whose own location evidence places them in Canada
#           (the GenBank location rule, 1b).
#   Step 2: keep ectomycorrhizal (EcM) records only: genera listed as EcM in
#           FungalTraits v1.2, plus records FungalTraits cannot classify whose
#           family is at least two-thirds EcM in FungalTraits (the family rule).
#   The result, emf_canada_em_only.csv, is the dataset every later script uses
#   (loaded as `emf` by 00_setup.R).
#
# INPUTS
#   paths$gf_long_out    GlobalFungi records (02_globalfungi.R)
#   paths$gb_long_out    GenBank records (03_genbank.R)
#   paths$canada_bound   Canada boundary (01_spatial_data.R)
#   paths$fungaltraits   FungalTraits v1.2 (data_raw, pinned)
#
# OUTPUTS (data_derived/)
#   emf_canada_combined.csv            all guilds, after the location rule
#   genbank_excluded_non_canadian.csv  GenBank records dropped, with the reason
#   emf_canada_em_only.csv             the EcM dataset
#   ecm_family_rule_table.csv          families considered by the family rule
#
# MANUSCRIPT ITEMS IT FEEDS
#   Every result: all shortfall analyses read the EcM dataset.
#
# RE-RUN WHEN
#   02 or 03 outputs change, the boundary (01) changes, or FungalTraits is
#   repinned. Each step is skipped when its output exists: delete
#   emf_canada_combined.csv (Step 1) and/or emf_canada_em_only.csv (Step 2).
#   Then re-run 05 onward.
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
  # coords_in_canada() (00_setup.R) tests each point against the Canada
  # boundary (simplified at 1 km in 01_spatial_data.R) buffered by 2 km. The
  # buffer absorbs that simplification and small coordinate imprecision at
  # coastal and island sites.
  combined$coord_in_canada <- coords_in_canada(combined$lat, combined$lon)

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
# A record is treated as ectomycorrhizal (EcM) in one of two ways.
#
#   1. Genus rule. Its UNITE genus is listed as EcM (primary_lifestyle) in
#      FungalTraits v1.2 (Põlme et al. 2020). This covers almost all records.
#
#   2. Family rule. FungalTraits cannot classify two kinds of record:
#        - records whose UNITE genus is absent from FungalTraits altogether,
#          mostly genera named after FungalTraits was compiled (e.g. the
#          genera split from Cortinarius in 2022, such as Phlegmacium);
#        - records whose SH has no genus label in UNITE, only a family
#          (e.g. "Boletaceae_gen_Incertae_sedis").
#      Such a record is treated as EcM if at least FAMILY_ECM_MIN_FRAC (two
#      thirds) of the genera FungalTraits lists in the record's UNITE family
#      are EcM. Two thirds admits families that are EcM in practice, such as
#      Cortinariaceae (8 of its 12 FungalTraits genera are EcM; the other four
#      are saprotrophs), and refuses mixed families such as Thelephoraceae
#      (8 of 14) and Pezizaceae (25 of 48).
#
# A record whose genus is in FungalTraits but NOT listed as EcM is never
# admitted: FungalTraits has classified it. The rules are implemented in
# ecm_reference() and classify_ecm() in 00_setup.R (FAMILY_ECM_MIN_FRAC is set
# there), and shared with 09_linnean.R.
#
# Records admitted by the family rule have no genus-level traits. They get the
# family's EcM lineage if every EcM genus in that family shares one lineage,
# and no exploration type. For SHs with no genus label, `genus` is set to NA
# (so they are never counted as a genus) and UNITE's label is kept in
# `genus_label_unite`. Column `ecm_basis` records which rule admitted each
# record ("genus" or "family").
#
# Output: paths$emf_data (the EcM dataset used by every later script), and
# paths$ecm_family_rule_table: every family the family rule considered, with
# its FungalTraits counts, whether it passed, and how many SHs and records it
# admitted or refused.

if (!(file.exists(paths$emf_data))) {

  # FungalTraits v1.2: pinned reference file (see 00_setup.R and
  # data_raw/fungaltraits/fungaltraits_version.txt for provenance).
  # ecm_reference() and classify_ecm() are defined in 00_setup.R, so this
  # script and 09_linnean.R use one identical EcM definition.
  ref <- ecm_reference()

  # Genus-level traits joined to records admitted by the genus rule.
  ft_join <- ref$ft |>
    dplyr::select(genus_lower, primary_lifestyle, secondary_lifestyle,
                  dplyr::matches("lineage_template"),
                  dplyr::matches("exploration_type_template")) |>
    dplyr::rename_with(~ sub("_template$", "", .x)) |>
    dplyr::distinct(genus_lower, .keep_all = TRUE)

  combined <- readr::read_csv(paths$emf_combined, show_col_types = FALSE,
                              guess_max = 1e6) |>
    dplyr::mutate(genus_lower  = tolower(trimws(genus)),
                  family_clean = sub("_.*$", "", family),
                  no_genus_label = !has_genus_label(genus),
                  ecm_basis    = classify_ecm(genus, family, ref))

  # Record of the family rule, for the supplement and for checking: every
  # family the rule considered that has at least one EcM genus in FungalTraits.
  combined |>
    dplyr::filter(ecm_basis %in% c("family", "family rule: refused")) |>
    dplyr::mutate(reason = dplyr::if_else(no_genus_label, "no genus label",
                                          "genus not in FungalTraits")) |>
    dplyr::inner_join(ref$family, by = c("family_clean" = "family")) |>
    dplyr::filter(ft_ecm_genera > 0) |>
    dplyr::group_by(family = family_clean, reason, ft_genera, ft_ecm_genera,
                    frac_ecm, family_passes) |>
    dplyr::summarise(n_sh = dplyr::n_distinct(sh_code, na.rm = TRUE),
                     n_records = dplyr::n(),
                     genus_labels = paste(sort(unique(genus)), collapse = "; "),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(family_passes), dplyr::desc(frac_ecm), dplyr::desc(n_sh)) |>
    readr::write_csv(paths$ecm_family_rule_table)

  admitted_genus <- combined |>
    dplyr::filter(ecm_basis == "genus") |>
    dplyr::left_join(ft_join, by = "genus_lower") |>
    dplyr::mutate(genus_label_unite = genus)

  admitted_family <- combined |>
    dplyr::filter(ecm_basis == "family") |>
    dplyr::left_join(dplyr::select(ref$family, family, family_lineage),
                     by = c("family_clean" = "family")) |>
    dplyr::mutate(genus_label_unite = genus,
                  genus = dplyr::if_else(no_genus_label, NA_character_, genus),
                  primary_lifestyle = "ectomycorrhizal",
                  ectomycorrhiza_lineage = family_lineage) |>
    dplyr::select(-family_lineage)

  em_only <- dplyr::bind_rows(admitted_genus, admitted_family) |>
    dplyr::select(-genus_lower, -family_clean, -no_genus_label)

  readr::write_csv(em_only, paths$emf_data)
}
