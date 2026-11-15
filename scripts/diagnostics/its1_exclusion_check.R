# =============================================================================
# Diagnostic: what would ITS1-only GlobalFungi samples add?
# =============================================================================
# PURPOSE
#   The pipeline keeps only GlobalFungi samples sequenced for ITS2 (alone or
#   with ITS1), because GlobalFungi assigns ITS1 and ITS2 sequences to SHs
#   independently, and because it lets GenBank records be processed the same
#   way (ITS2 only). This check measures what that choice leaves out: it takes
#   the Canadian samples sequenced for ITS1 only, applies every other filter
#   used in 02_globalfungi.R (unmanipulated; soil or root), keeps their EcM
#   records (same EcM definition as 04, classify_ecm() in 00_setup.R), and
#   counts what they would add to the EcM dataset.
#
# INPUTS
#   paths$gf_metadata      GlobalFungi v5 sample metadata (data_raw)
#   paths$gf_sh_abundance  GlobalFungi v5 SH abundance matrix (~13 GB, data_raw)
#   paths$unite_taxonomy   SH -> taxonomy lookup (02_globalfungi.R)
#   emf                    the EcM dataset (04_combine_ecm_dataset.R)
#
# OUTPUTS
#   paths$out_diagnostics/its1_exclusion_check.csv   metric, value
#   paths$out_diagnostics/its1_only_added_taxa.csv   the SH codes that ITS1-only
#                                                    samples would add, with
#                                                    their UNITE taxonomy
#
# MANUSCRIPT ITEMS IT FEEDS
#   The SM1 paragraph and reviewer response on excluding ITS1-only samples.
#
# RE-RUN WHEN
#   The EcM dataset changes (any re-run of 02-04), or GlobalFungi is updated.
#   Run from the project root:  Rscript scripts/diagnostics/its1_exclusion_check.R
#   Takes several minutes (one pass through the 13 GB matrix).
# =============================================================================

source(here::here("scripts", "00_setup.R"))

if (!exists("emf")) stop("The EcM dataset is missing: run 01-04 first.")

# ---- Step 1: Canadian ITS1-only samples passing the other filters ------------
# Same filters as Step 2 of 02_globalfungi.R, except the barcode region.
meta <- data.table::fread(paths$gf_metadata, sep = "\t", quote = "") |>
  as.data.frame()
its1_meta <- dplyr::filter(meta,
                           country          == "Canada",
                           barcoding_region == "ITS1",
                           manipulated      == "NO",
                           sample_type      %in% c("soil", "root"))
rm(meta)

# ---- Step 2: their rows of the SH abundance matrix ---------------------------
# Same awk approach as Step 3 of 02_globalfungi.R: keep the header line and
# every row whose first field (sample_ID) is in the ID list.
ids_file <- tempfile(fileext = ".txt")
sub_file <- tempfile(fileext = ".tsv")
writeLines(its1_meta$sample_ID, ids_file)
cmd <- paste("awk -F'\\t'",
             "'NR==FNR{ids[$1]=1; next} FNR==1{print; next} $1 in ids {print}'",
             shQuote(ids_file), shQuote(paths$gf_sh_abundance), ">", shQuote(sub_file))
if (system(cmd) != 0L) stop("awk extraction failed")

its1_long <- data.table::fread(sub_file, sep = "\t", quote = "") |>
  data.table::melt(id.vars = "sample_ID", variable.name = "sh_code",
                   value.name = "abundance", variable.factor = FALSE) |>
  dplyr::filter(abundance > 0) |>
  dplyr::left_join(readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE),
                   by = "sh_code") |>
  dplyr::left_join(dplyr::select(its1_meta, sample_ID, latitude, longitude),
                   by = "sample_ID")
unlink(c(ids_file, sub_file))

# ---- Step 3: keep EcM records (same definition as 04) -------------------------
ref <- ecm_reference()
its1_ecm <- its1_long |>
  dplyr::mutate(ecm_basis = classify_ecm(genus, family, ref)) |>
  dplyr::filter(ecm_basis %in% c("genus", "family"))

# ---- Step 4: what would they add? --------------------------------------------
is_named <- function(sp) !is.na(sp) & !grepl("_sp$", sp)
emf_sh      <- unique(emf$sh_code[!is.na(emf$sh_code)])
emf_species <- unique(emf$species[is_named(emf$species)])
emf_cells   <- emf |>
  dplyr::filter(coord_in_canada %in% TRUE) |>
  dplyr::distinct(cell = paste(snap_30arcsec(lat), snap_30arcsec(lon))) |>
  dplyr::pull(cell)

its1_ecm <- its1_ecm |>
  dplyr::mutate(in_canada = coords_in_canada(latitude, longitude) %in% TRUE,
                cell = paste(snap_30arcsec(latitude), snap_30arcsec(longitude)))

new_sh      <- setdiff(unique(its1_ecm$sh_code), emf_sh)
new_species <- setdiff(unique(its1_ecm$species[is_named(its1_ecm$species)]), emf_species)
new_cells   <- setdiff(unique(its1_ecm$cell[its1_ecm$in_canada]), emf_cells)

summary_tbl <- tibble::tibble(
  metric = c("Canadian ITS1-only samples (unmanipulated, soil or root)",
             "EcM records in those samples",
             "EcM SH codes in those samples",
             "EcM SH codes not already in the EcM dataset",
             "SH codes in the EcM dataset",
             "% increase in SH codes if ITS1-only samples were added",
             "Named species not already in the EcM dataset",
             "Named species in the EcM dataset",
             "% increase in named species if ITS1-only samples were added",
             "30 arc-second grid cells in Canada not already sampled"),
  value = c(nrow(its1_meta), nrow(its1_ecm), dplyr::n_distinct(its1_ecm$sh_code),
            length(new_sh), length(emf_sh),
            round(100 * length(new_sh) / length(emf_sh), 1),
            length(new_species), length(emf_species),
            round(100 * length(new_species) / length(emf_species), 1),
            length(new_cells)))

readr::write_csv(summary_tbl, file.path(paths$out_diagnostics, "its1_exclusion_check.csv"))
its1_ecm |>
  dplyr::filter(sh_code %in% new_sh) |>
  dplyr::group_by(sh_code, family, genus, species) |>
  dplyr::summarise(n_samples = dplyr::n_distinct(sample_ID), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_samples)) |>
  readr::write_csv(file.path(paths$out_diagnostics, "its1_only_added_taxa.csv"))

print(summary_tbl, n = Inf)
