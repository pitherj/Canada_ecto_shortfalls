# =============================================================================
# Diagnostic: what the GenBank location rule kept and dropped
# =============================================================================
# PURPOSE
#   04_combine_ecm_dataset.R keeps a GenBank record only if its own location
#   evidence places it in Canada (location field begins "Canada", or no field
#   and coordinates inside Canada). This check summarizes the result: records
#   kept and dropped by reason, the countries named by dropped records, the
#   kept records whose coordinates disagree with their "Canada" field, and the
#   EcM taxa that the dropped records alone would have contributed.
#
# INPUTS
#   paths$emf_combined      combined dataset after the rule (04)
#   paths$genbank_excluded  GenBank records dropped by the rule (04)
#   paths$fungaltraits      for the EcM definition (classify_ecm(), 00_setup.R)
#   emf                     the EcM dataset (04)
#
# OUTPUTS (paths$out_diagnostics/)
#   genbank_location_summary.csv          metric, value
#   genbank_excluded_countries.csv        dropped records by country named
#   genbank_kept_coordinate_conflicts.csv kept records ("Canada" field) whose
#                                         coordinates fall outside Canada
#
# MANUSCRIPT ITEMS IT FEEDS
#   SM1 text on the GenBank location rule; reviewer responses on the search.
#
# RE-RUN WHEN
#   04_combine_ecm_dataset.R is re-run.
#   Run from the project root:  Rscript scripts/diagnostics/genbank_location_summary.R
# =============================================================================

source(here::here("scripts", "00_setup.R"))

ct <- readr::cols(.default = "c")
kept <- readr::read_csv(paths$emf_combined, col_types = ct) |>
  dplyr::filter(source == "GenBank")
dropped <- readr::read_csv(paths$genbank_excluded, col_types = ct)

# EcM taxa that only the dropped records carry (same EcM definition as 04).
ref <- ecm_reference()
dropped_ecm <- dropped |>
  dplyr::mutate(ecm_basis = classify_ecm(genus, NA_character_, ref)) |>
  dplyr::filter(ecm_basis == "genus")
is_named <- function(sp) !is.na(sp) & !grepl("_sp$", sp)
lost_sh      <- setdiff(unique(stats::na.omit(dropped_ecm$sh_code)), unique(emf$sh_code))
lost_species <- setdiff(unique(dropped_ecm$species[is_named(dropped_ecm$species)]),
                        unique(emf$species))
lost_genera  <- setdiff(unique(stats::na.omit(dropped_ecm$genus)),
                        unique(stats::na.omit(emf$genus)))

conflicts <- kept |>
  dplyr::filter(location_evidence == "location field: Canada",
                coord_in_canada == "FALSE") |>
  dplyr::select(accession, organism, sh_code, genus, species, country_gb, lat_lon_gb)

by_reason <- dplyr::bind_rows(
  dplyr::count(kept, location_evidence) |> dplyr::mutate(status = "kept"),
  dplyr::count(dropped, location_evidence) |> dplyr::mutate(status = "dropped"))

summary_tbl <- tibble::tibble(
  metric = c(paste0("GenBank records ", by_reason$status, ": ", by_reason$location_evidence),
             "Kept records whose coordinates fall outside Canada (field says Canada)",
             "EcM SH codes carried only by dropped records",
             "EcM named species carried only by dropped records",
             "EcM genera carried only by dropped records"),
  value = c(by_reason$n, nrow(conflicts), length(lost_sh), length(lost_species),
            length(lost_genera)))

readr::write_csv(summary_tbl, file.path(paths$out_diagnostics, "genbank_location_summary.csv"))
dropped |>
  dplyr::filter(location_evidence == "location field: another country") |>
  dplyr::mutate(country = sub(":.*$", "", country_gb)) |>
  dplyr::count(country, sort = TRUE, name = "n_records") |>
  readr::write_csv(file.path(paths$out_diagnostics, "genbank_excluded_countries.csv"))
readr::write_csv(conflicts, file.path(paths$out_diagnostics, "genbank_kept_coordinate_conflicts.csv"))

print(summary_tbl, n = Inf)
