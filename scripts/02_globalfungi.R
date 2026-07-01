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
#     The reference build is PINNED, not auto-detected — see
#     paths$unite_fasta in 00_setup.R and
#     docs/unite_sh_code_mismatch_memo.md for why. UNITE SH codes are
#     renumbered/merged across builds, so picking a different build than the
#     one GlobalFungi used to pre-assign its SH codes silently breaks this
#     join. Do not repin without re-validating coverage (see the gate at
#     Step 4 below) against the new build first.
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
# UNITE taxonomy join (Step 4 below). The validated, pinned 2024-04-04 build
# leaves a documented, investigated residual of 0.93% (723 / 77,793 codes;
# see docs/unite_sh_code_mismatch_memo.md §3.8) that could not be resolved
# further. This threshold sits just above that known residual: it tolerates
# the documented gap but will trip if a future UNITE build (or other
# upstream change) reintroduces a large-scale mismatch like the one this
# pin was put in place to prevent (~97% unmatched, see memo §2–§3).
SH_MAX_UNMATCHED_FRAC <- 0.02

# ---- Step 0: Locate UNITE FASTA (pinned reference, not auto-detected) -------
#
# The UNITE build is pinned explicitly via paths$unite_fasta in
# 00_setup.R (currently 2024-04-04, _dev variant). This used to be
# auto-detected as "most recent by filename date," which silently broke the
# Step 4 taxonomy join when a newer UNITE build renumbered SH codes — see
# docs/unite_sh_code_mismatch_memo.md for the full investigation. Repinning
# to a different build requires re-validating coverage at the Step 4 gate
# below first.

unite_fasta_path <- paths$unite_fasta

if (!file.exists(unite_fasta_path)) {
  stop(
    "Pinned UNITE FASTA not found: ", unite_fasta_path, "\n",
    "Either restore this file to data_raw/UNITE/, or — if intentionally\n",
    "switching builds — update paths$unite_fasta (and\n",
    "paths$unite_fasta, which must match it) in 00_setup.R, then\n",
    "re-validate match coverage against GlobalFungi's SH codes before\n",
    "trusting the new build. See docs/unite_sh_code_mismatch_memo.md."
  )
}

cur_date_str <- regmatches(
  basename(unite_fasta_path),
  regexpr("\\d{2}\\.\\d{2}\\.\\d{4}", basename(unite_fasta_path))
)
cur_date_fmt <- if (length(cur_date_str) == 1L) {
  format(as.Date(cur_date_str, "%d.%m.%Y"), "%Y-%m-%d")
} else {
  "unknown date"
}
ts(sprintf("Step 0: UNITE release (pinned): %s (%s)", basename(unite_fasta_path), cur_date_fmt))

# ---- Step 1: Parse UNITE FASTA headers → SH taxonomy lookup -----------------

if (!file.exists(paths$unite_taxonomy)) {
  ts("Step 1: Parsing UNITE FASTA headers...")

  headers_raw <- system(paste("grep '^>'", shQuote(unite_fasta_path)), intern = TRUE)
  ts("  Read", length(headers_raw), "FASTA header lines")
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
  ts("  Unique SH codes:", nrow(sh_lookup))
  readr::write_csv(sh_lookup, paths$unite_taxonomy)
  ts("  Saved -> temp/unite_sh_taxonomy.csv")

} else {
  ts("Step 1: UNITE taxonomy lookup already exists — skipping.")
}

# ---- Step 2: Filter GlobalFungi metadata for Canada samples -----------------

if (!file.exists(paths$gf_meta_out)) {
  ts("Step 2: Filtering GlobalFungi metadata for Canada samples...")

  meta <- data.table::fread(paths$gf_metadata, sep = "\t", quote = "")
  ts("  Loaded", nrow(meta), "total samples")

  canada_meta <- dplyr::filter(
    as.data.frame(meta),
    country          == "Canada",
    barcoding_region %in% c("ITS2", "ITSboth"),
    manipulated      == "NO",
    !sample_type     %in% c("shoot", "air", "water", "sediment")
  )

  ts("  Retained", nrow(canada_meta), "Canada samples")
  ts("  By barcoding_region:"); print(dplyr::count(canada_meta, barcoding_region))
  ts("  By sample_type:");     print(dplyr::count(canada_meta, sample_type))

  readr::write_csv(canada_meta, paths$gf_meta_out)
  writeLines(canada_meta$sample_ID, paths$gf_ids_out)
  ts("  Saved -> temp/globalfungi_canada_metadata.csv + _ids.txt")

} else {
  ts("Step 2: Canada metadata already exists — skipping.")
  n_ids <- length(readLines(paths$gf_ids_out))
  ts(sprintf("  %d Canada sample IDs available.", n_ids))
}

# ---- Step 3: Extract Canada rows from 13 GB SH abundance matrix (awk) -------

if (!file.exists(paths$gf_sh_subset_out)) {
  ts("Step 3: Extracting Canada rows from 13 GB SH abundance matrix (awk)...")
  ts("  This may take a few minutes.")

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
  ts(sprintf("  Done. Output: %s MB -> temp/globalfungi_canada_SH_abundance.txt", sz))

} else {
  sz <- round(file.info(paths$gf_sh_subset_out)$size / 1e6, 1)
  ts(sprintf("Step 3: Canada SH abundance file already exists (%s MB) — skipping.", sz))
}

# ---- Step 4: Pivot to long format, join taxonomy + metadata -----------------

if (!file.exists(paths$gf_long_out)) {
  ts("Step 4: Pivoting SH abundance matrix to long format...")

  canada_sh <- data.table::fread(paths$gf_sh_subset_out, sep = "\t", quote = "")
  ts(sprintf("  Dimensions: %d samples x %d SH columns",
             nrow(canada_sh), ncol(canada_sh) - 1L))

  canada_long <- data.table::melt(
    canada_sh,
    id.vars         = "sample_ID",
    variable.name   = "sh_code",
    value.name      = "abundance",
    variable.factor = FALSE
  )
  canada_long <- dplyr::filter(canada_long, abundance > 0)
  ts("  Non-zero records:", nrow(canada_long))

  ts("  Joining UNITE taxonomy...")
  sh_lookup   <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
  canada_long <- dplyr::left_join(canada_long, sh_lookup, by = "sh_code")

  # ---- Coverage-check gate -------------------------------------------------
  # Fail fast rather than silently propagating NA taxonomy downstream — this
  # is the exact failure mode documented in
  # docs/unite_sh_code_mismatch_memo.md (a UNITE build mismatch previously
  # left ~97% of rows unmatched with no error raised). Reported at both the
  # row level (consistent with other diagnostics in this script) and the
  # unique-SH-code level for diagnosability; unmatched codes are persisted
  # to a checkpoint rather than just logged.
  unmatched_rows <- is.na(canada_long$kingdom)
  pct_unmatched  <- mean(unmatched_rows)
  unmatched_sh   <- unique(canada_long$sh_code[unmatched_rows])
  ts(sprintf("  UNITE taxonomy coverage: %.2f%% of rows unmatched (%d unique SH codes)",
             100 * pct_unmatched, length(unmatched_sh)))

  if (length(unmatched_sh) > 0L) {
    readr::write_csv(tibble::tibble(sh_code = unmatched_sh), paths$gf_sh_unmatched)
    ts("  Unmatched SH codes written -> ", paths$gf_sh_unmatched)
  }

  if (pct_unmatched > SH_MAX_UNMATCHED_FRAC) {
    stop(sprintf(
      paste0(
        "UNITE taxonomy join left %.2f%% of GlobalFungi rows unmatched ",
        "(threshold: %.2f%%). This usually means the pinned UNITE build ",
        "(paths$unite_fasta in 00_setup.R) no longer matches ",
        "the build GlobalFungi used to assign its SH codes. See ",
        "docs/unite_sh_code_mismatch_memo.md before changing the pin, and ",
        "inspect %s for the unmatched codes."
      ),
      100 * pct_unmatched, 100 * SH_MAX_UNMATCHED_FRAC, paths$gf_sh_unmatched
    ))
  }

  ts("  Joining sample metadata...")
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
  ts(sprintf("  Saved globalfungi_canada_long.csv  (%d rows, %d unique SHs, %d unique samples)",
             nrow(canada_long),
             dplyr::n_distinct(canada_long$sh_code),
             dplyr::n_distinct(canada_long$sample_ID)))

} else {
  ts("Step 4: globalfungi_canada_long.csv already exists — skipping.")
}

ts("02_globalfungi.R complete.  Run 03_genbank.R next.")
