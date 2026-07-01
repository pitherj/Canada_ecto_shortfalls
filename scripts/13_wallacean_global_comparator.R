# =============================================================================
# Wallacean Shortfall: Global GlobalFungi Comparator Metrics
# =============================================================================
# Computes summary metrics describing the FULL global GlobalFungi v5 database
# for use as a comparator against the Canada-specific dataset in Table S1 of
# EMF_shortfalls_SI.Rmd ("EcM Fungal Records Dataset").
#
# Two independently-guarded stages, following the per-step sentinel pattern
# used in 10_linnean_inext.R:
#
#   Stage 1 (cheap — header/taxonomy/metadata only; seconds):
#     - Total global GlobalFungi v5 samples (raw, all taxa)
#     - Global samples passing the same quality filters applied to the
#       Canadian dataset in 02_globalfungi.R
#     - Global EcM SH codes / genera / named species detected as columns in
#       the abundance matrix. Identification logic is identical to
#       12_wallacean_density_map.R Step 1 (UNITE x FungalTraits genus
#       cross-reference, with the same exact-match-then-prefix-fallback rule
#       for matching SH codes to matrix column names).
#     - Global unique EcM sampling locations, derived from the existing
#       checkpoint data_derived/checkpoints/gf_global_ecm_sample_ids.csv
#       (produced by 12_wallacean_density_map.R).
#     Output: data_derived/checkpoints/gf_global_comparator_cheap.csv
#
#   Stage 2 (expensive — full 13 GB matrix scan via data.table::fread(select=...);
#   5-15+ minutes on a modern laptop, per the runtime note in 2-6a):
#     - Total global EcM detection records (sample x SH-code pairs with
#       abundance > 0), restricted to quality-filtered samples
#     - Total global EcM raw read abundance (sum of abundance values)
#     - Total global quality-filtered samples with >=1 EcM detection
#       (the global analogue of the Canadian "GlobalFungi samples,
#       quality-filtered, with EcM fungal records" count in Table S1;
#       derived from the same in-memory EcM-column matrix subset used for
#       the two metrics above, so no additional matrix read is required)
#     Output: data_derived/checkpoints/gf_global_comparator_volume.csv
#
# Prerequisites:
#   data_derived/checkpoints/unite_sh_taxonomy.csv         (from 0-3 / 0-4)
#   data_derived/checkpoints/gf_global_ecm_sample_ids.csv  (from 2-6a — run that
#                                                       script first if absent)
#
# Note on Stage 2: this is a NEW, expensive pipeline stage. It is split out
# from Stage 1 so that Table S1 can render with the cheap comparator numbers
# even before Stage 2 has been run. Delete the Stage 2 output file to force
# re-computation.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(data.table)

out_cheap  <- here::here("data_derived", "checkpoints", "gf_global_comparator_cheap.csv")
out_volume <- here::here("data_derived", "checkpoints", "gf_global_comparator_volume.csv")

# A small helper used by both stages to identify which UNITE SH codes belong
# to EcM genera, and which of those are present as columns in the abundance
# matrix. Kept as a local function (not added to 00_setup.R) since it is
# only used by this script and 12_wallacean_density_map.R, which
# implements the same logic inline rather than via a shared helper.
identify_ecm_columns <- function() {
  unite_tax <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE) |>
    dplyr::mutate(genus_clean = tolower(sub("^g__", "", genus)))

  ft <- data.table::fread(paths$fungaltraits, select = c("GENUS", "primary_lifestyle"))
  ft[, genus_lower := tolower(GENUS)]
  ecm_genera_lower <- ft[primary_lifestyle == "ectomycorrhizal", tolower(GENUS)]

  ecm_sh_codes <- unite_tax |>
    dplyr::filter(genus_clean %in% ecm_genera_lower) |>
    dplyr::pull(sh_code) |>
    unique()

  gf_header_cols <- names(
    data.table::fread(paths$gf_sh_abundance, sep = "\t", quote = "", nrows = 0L)
  )
  gf_header_cols <- gf_header_cols[-1]  # drop sample_ID

  ecm_cols_present <- intersect(ecm_sh_codes, gf_header_cols)
  if (length(ecm_cols_present) < 0.1 * length(ecm_sh_codes)) {
    ecm_prefix   <- sub("\\.[0-9]+FU$", "", ecm_sh_codes)
    col_prefix   <- sub("\\.[0-9]+FU$", "", gf_header_cols)
    matched_idx  <- match(ecm_prefix, col_prefix)
    ecm_cols_present <- unique(gf_header_cols[matched_idx[!is.na(matched_idx)]])
  }

  list(unite_tax = unite_tax, ecm_genera_lower = ecm_genera_lower,
       ecm_cols_present = ecm_cols_present)
}

# A second helper: the quality-filtered global sample_ID set, using the same
# filters as 02_globalfungi.R / 2-6a.
quality_filtered_sample_ids <- function() {
  gf_meta <- data.table::fread(
    paths$gf_metadata, sep = "\t", quote = "",
    select = c("sample_ID", "latitude", "longitude",
               "barcoding_region", "manipulated", "sample_type")
  )
  list(
    meta = gf_meta,
    filtered_ids = gf_meta[
      barcoding_region %in% c("ITS2", "ITSboth") &
      manipulated == "NO" &
      !sample_type %in% c("shoot", "air", "water", "sediment") &
      !is.na(latitude) & !is.na(longitude) &
      latitude  >= -90  & latitude  <= 90 &
      longitude >= -180 & longitude <= 180,
      sample_ID
    ]
  )
}

# =============================================================================
# STAGE 1: cheap header / taxonomy / metadata-based metrics
# =============================================================================

if (file.exists(out_cheap)) {
  ts("Stage 1: gf_global_comparator_cheap.csv already exists — skipping.")
} else {

  ts("Stage 1: Computing cheap global GlobalFungi comparator metrics...")

  ecm <- identify_ecm_columns()
  ts(sprintf("  Global EcM SH codes detected in matrix: %d", length(ecm$ecm_cols_present)))

  ecm_tax <- dplyr::filter(ecm$unite_tax, sh_code %in% ecm$ecm_cols_present)
  n_ecm_sh_global    <- length(ecm$ecm_cols_present)
  n_ecm_genus_global <- dplyr::n_distinct(ecm_tax$genus_clean)
  n_ecm_sp_global    <- dplyr::n_distinct(
    ecm_tax$species[!is.na(ecm_tax$species) & !grepl("_sp$", ecm_tax$species)]
  )

  filt <- quality_filtered_sample_ids()
  n_samples_global_raw      <- nrow(filt$meta)
  n_samples_global_filtered <- length(filt$filtered_ids)

  ecm_sample_ids_path <- here::here("data_derived", "checkpoints", "gf_global_ecm_sample_ids.csv")
  if (!file.exists(ecm_sample_ids_path)) {
    stop(
      "Prerequisite checkpoint missing: data_derived/checkpoints/gf_global_ecm_sample_ids.csv. ",
      "Run scripts/12_wallacean_density_map.R first."
    )
  }
  n_ecm_locs_global <- readr::read_csv(ecm_sample_ids_path, show_col_types = FALSE) |>
    dplyr::distinct(latitude, longitude) |>
    nrow()

  out <- tibble::tibble(
    metric = c(
      "GlobalFungi v5 samples worldwide (raw, all taxa)",
      "GlobalFungi v5 samples worldwide (quality-filtered, all taxa)",
      "EcM SH codes detected (GlobalFungi-wide, matrix presence)",
      "EcM genera detected (GlobalFungi-wide, matrix presence)",
      "EcM named species detected (GlobalFungi-wide, matrix presence)",
      "EcM unique sampling locations (GlobalFungi-wide)"
    ),
    value = c(
      n_samples_global_raw,
      n_samples_global_filtered,
      n_ecm_sh_global,
      n_ecm_genus_global,
      n_ecm_sp_global,
      n_ecm_locs_global
    )
  )
  readr::write_csv(out, out_cheap)
  ts("  Saved gf_global_comparator_cheap.csv")
}

# =============================================================================
# STAGE 2: expensive full-matrix scan (global EcM detection volume)
# =============================================================================

if (file.exists(out_volume)) {
  ts("Stage 2: gf_global_comparator_volume.csv already exists — skipping.")
} else {

  ts("Stage 2: Scanning full abundance matrix for global EcM detection volume...")
  ts("  (Reads several hundred EcM SH columns from ~13 GB; expect 5-15+ minutes.)")

  ecm  <- identify_ecm_columns()
  filt <- quality_filtered_sample_ids()

  gf_ecm_mat <- data.table::fread(
    paths$gf_sh_abundance, sep = "\t", quote = "",
    select = c("sample_ID", ecm$ecm_cols_present)
  )
  gf_ecm_mat <- gf_ecm_mat[sample_ID %in% filt$filtered_ids]
  ts(sprintf("  Matrix read and filtered: %d quality-filtered samples x %d EcM SH columns",
             nrow(gf_ecm_mat), length(ecm$ecm_cols_present)))

  vals <- gf_ecm_mat[, -"sample_ID", with = FALSE]
  n_ecm_records_global <- sum(vals > 0L)
  n_ecm_reads_global   <- sum(as.numeric(unlist(vals, use.names = FALSE)))
  # Samples with at least one positive-abundance EcM detection (row-wise any()
  # over the EcM columns already loaded above for n_ecm_records_global — this
  # is the global analogue of the Canadian "quality-filtered, with EcM fungal
  # records" sample count, and requires no additional matrix read).
  n_samples_with_ecm_global <- sum(rowSums(vals > 0L) > 0L)

  out <- tibble::tibble(
    metric = c(
      "EcM detection records (GlobalFungi-wide, quality-filtered samples)",
      "EcM total read abundance, Sigma reads (GlobalFungi-wide, quality-filtered samples)",
      "GlobalFungi v5 samples worldwide (quality-filtered, with EcM detection)"
    ),
    value = c(n_ecm_records_global, n_ecm_reads_global, n_samples_with_ecm_global)
  )
  readr::write_csv(out, out_volume)
  ts("  Saved gf_global_comparator_volume.csv")
}

ts("13_wallacean_global_comparator.R complete.")
