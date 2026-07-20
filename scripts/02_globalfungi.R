# =============================================================================
# GlobalFungi Sequence Extraction
# =============================================================================
# Extracts Canadian EcM fungal sequence records from GlobalFungi v5.
# Each step is checkpointed — safe to re-run; skips completed steps.
#
# Manual prerequisites (no public API — download manually):
#   GlobalFungi v5:  https://globalfungi.com  → Downloads
#     data_raw/GlobalFungi/GlobalFungi_5_sample_metadata.txt       (~75 MB)
#     data_raw/GlobalFungi/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt (~13 GB)
#     data_raw/GlobalFungi/GlobalFungi_5_species_abundance_ITS1_ITS2.txt
#
#   UNITE general FASTA (_dev variant):  https://unite.ut.ee/repository.php
#     The reference build is PINNED, not auto-detected (paths$unite_fasta in
#     00_setup.R). UNITE SH codes are renumbered/merged across builds, so a
#     build other than the one GlobalFungi used to pre-assign its SH codes
#     silently breaks this join. Do not repin without re-validating coverage
#     (see the gate at Step 4 below) against the new build first.
#
#   To update either dataset: delete the corresponding data_derived/ checkpoint
#   files so the relevant step regenerates on the next run.
#
# Outputs (all in data_derived/):
#   temp/unite_sh_taxonomy.csv              — SH code → taxonomy lookup
#   temp/globalfungi_canada_metadata.csv    — Canada sample metadata
#   temp/globalfungi_canada_ids.txt         — Canada sample ID list (for awk)
#   temp/globalfungi_canada_SH_abundance.txt — Canada rows from 13 GB matrix
#   globalfungi_canada_long.csv             — final long-format table
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(data.table)

# Maximum tolerated fraction of GlobalFungi SH codes left unmatched after the
# UNITE taxonomy join (Step 4 below). The pinned 2024-04-04 build leaves an
# irreducible residual of 0.93% (723 / 77,793 codes). This threshold sits just
# above that residual: it tolerates the expected gap, but trips if a future
# UNITE build (or other upstream change) introduces a large-scale mismatch,
# which manifests as the great majority of codes failing to join.
SH_MAX_UNMATCHED_FRAC <- 0.02

# ---- Step 0: Locate UNITE FASTA (pinned reference, not auto-detected) -------
#
# The UNITE build is pinned explicitly via paths$unite_fasta in
# 00_setup.R (currently 2024-04-04, _dev variant). Pinning matters because a
# newer UNITE build renumbers SH codes, which breaks the Step 4 taxonomy join.
# Repinning to a different build requires re-validating coverage at the Step 4
# gate below first.

unite_fasta_path <- paths$unite_fasta

if (!file.exists(unite_fasta_path)) {
  stop(
    "Pinned UNITE FASTA not found: ", unite_fasta_path, "\n",
    "Either restore this file to data_raw/UNITE/, or — if intentionally\n",
    "switching builds — update paths$unite_fasta (and\n",
    "paths$unite_fasta, which must match it) in 00_setup.R, then\n",
    "re-validate match coverage against GlobalFungi's SH codes before\n",
    "trusting the new build."
  )
}

# ---- Step 1: Parse UNITE FASTA headers → SH taxonomy lookup -----------------

if (!file.exists(paths$unite_taxonomy)) {

  headers_raw <- system(paste("grep '^>'", shQuote(unite_fasta_path)), intern = TRUE)
  headers_raw <- sub("^>", "", headers_raw)

  # Format: name|accession|SH_code|refs_type|taxonomy_string
  parts     <- strsplit(headers_raw, "|", fixed = TRUE)
  sh_code   <- vapply(parts, `[[`, character(1L), 3L)
  taxonomy  <- vapply(parts, function(x) x[length(x)], character(1L))

  pull_level <- function(tax, prefix) {
    has_match         <- grepl(prefix, tax, fixed = TRUE)
    result            <- rep(NA_character_, length(tax))
    result[has_match] <- sub(paste0(".*", prefix, "([^;|]*).*"), "\\1", tax[has_match])
    result
  }

  dt <- tibble::tibble(
    sh_code = sh_code,
    kingdom = pull_level(taxonomy, "k__"),
    phylum  = pull_level(taxonomy, "p__"),
    class   = pull_level(taxonomy, "c__"),
    order   = pull_level(taxonomy, "o__"),
    family  = pull_level(taxonomy, "f__"),
    genus   = pull_level(taxonomy, "g__"),
    species = pull_level(taxonomy, "s__")
  )

  sh_tax <- dplyr::distinct(dt, sh_code, kingdom, phylum, class, order, family, genus)

  # Species consensus: most common non-_sp name if it accounts for >= 50% of named seqs
  sh_species <- dt |>
    dplyr::group_by(sh_code) |>
    dplyr::summarise(
      species = {
        all_sp <- species[!is.na(species)]
        named  <- all_sp[!grepl("_sp$", all_sp, ignore.case = TRUE)]
        if (length(named) == 0L) {
          all_sp[1L]
        } else {
          tbl      <- sort(table(named), decreasing = TRUE)
          top_frac <- tbl[[1L]] / length(named)
          if (top_frac >= 0.5) names(tbl)[[1L]]
          else paste0(sub("_.*", "", names(tbl)[[1L]]), "_sp")
        }
      },
      .groups = "drop"
    )

  sh_lookup <- dplyr::left_join(sh_tax, sh_species, by = "sh_code")
  readr::write_csv(sh_lookup, paths$unite_taxonomy)

}

# ---- Step 2: Filter GlobalFungi metadata for Canada samples -----------------

if (!file.exists(paths$gf_meta_out)) {

  meta <- data.table::fread(paths$gf_metadata, sep = "\t", quote = "")

  canada_meta <- dplyr::filter(
    as.data.frame(meta),
    country          == "Canada",
    barcoding_region %in% c("ITS2", "ITSboth"),
    manipulated      == "NO",
    !sample_type     %in% c("shoot", "air", "water", "sediment")
  )

  readr::write_csv(canada_meta, paths$gf_meta_out)
  writeLines(canada_meta$sample_ID, paths$gf_ids_out)

}

# ---- Step 3: Extract Canada rows from 13 GB SH abundance matrix (awk) -------

if (!file.exists(paths$gf_sh_subset_out)) {

  cmd <- paste(
    "awk -F'\\t'",
    "'NR==FNR{ids[$1]=1; next} FNR==1{print; next} $1 in ids {print}'",
    shQuote(paths$gf_ids_out),
    shQuote(paths$gf_sh_abundance),
    ">",
    shQuote(paths$gf_sh_subset_out)
  )

  ret <- system(cmd)
  if (ret != 0L) stop("awk extraction failed with exit code ", ret)
  sz <- round(file.info(paths$gf_sh_subset_out)$size / 1e6, 1)

} else {
  sz <- round(file.info(paths$gf_sh_subset_out)$size / 1e6, 1)
}

# ---- Step 4: Pivot to long format, join taxonomy + metadata -----------------

if (!file.exists(paths$gf_long_out)) {

  canada_sh <- data.table::fread(paths$gf_sh_subset_out, sep = "\t", quote = "")

  canada_long <- data.table::melt(
    canada_sh,
    id.vars         = "sample_ID",
    variable.name   = "sh_code",
    value.name      = "abundance",
    variable.factor = FALSE
  )
  canada_long <- dplyr::filter(canada_long, abundance > 0)

  sh_lookup   <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
  canada_long <- dplyr::left_join(canada_long, sh_lookup, by = "sh_code")

  # ---- Coverage-check gate -------------------------------------------------
  # Fail fast rather than silently propagating NA taxonomy downstream: a UNITE
  # build mismatch can leave the great majority of rows unmatched with no error
  # raised. Reported at both the
  # row level (consistent with other diagnostics in this script) and the
  # unique-SH-code level for diagnosability; unmatched codes are persisted
  # to a checkpoint rather than just logged.
  unmatched_rows <- is.na(canada_long$kingdom)
  pct_unmatched  <- mean(unmatched_rows)
  unmatched_sh   <- unique(canada_long$sh_code[unmatched_rows])

  if (length(unmatched_sh) > 0L) {
    readr::write_csv(tibble::tibble(sh_code = unmatched_sh), paths$gf_sh_unmatched)
  }

  if (pct_unmatched > SH_MAX_UNMATCHED_FRAC) {
    stop(sprintf(
      paste0(
        "UNITE taxonomy join left %.2f%% of GlobalFungi rows unmatched ",
        "(threshold: %.2f%%). This usually means the pinned UNITE build ",
        "(paths$unite_fasta in 00_setup.R) no longer matches ",
        "the build GlobalFungi used to assign its SH codes. Re-validate ",
        "coverage before changing the pin, and ",
        "inspect %s for the unmatched codes."
      ),
      100 * pct_unmatched, 100 * SH_MAX_UNMATCHED_FRAC, paths$gf_sh_unmatched
    ))
  }

  canada_meta <- readr::read_csv(paths$gf_meta_out, show_col_types = FALSE)
  meta_cols <- c(
    "sample_ID", "latitude", "longitude", "country",
    "barcoding_region", "sample_type", "environment_type",
    "ecosystem_classification", "year_of_sampling_from",
    "dominant_plant_species", "other_plant_species",
    "paper_ID", "paper_title", "paper_authors", "paper_doi"
  )
  canada_long <- dplyr::left_join(
    as.data.frame(canada_long),
    dplyr::select(canada_meta, dplyr::any_of(meta_cols)),
    by = "sample_ID"
  )

  readr::write_csv(canada_long, paths$gf_long_out)

}

