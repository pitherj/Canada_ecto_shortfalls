# =============================================================================
# 00_setup.R  —  Shared setup for the ECM_manuscript pipeline
# =============================================================================
# WHAT THIS FILE IS
#   Every numbered script in scripts/ begins with:
#       source(here::here("scripts", "00_setup.R"))
#   Sourcing it (a) loads the handful of always-on packages, (b) defines the
#   coordinate reference systems, small helper functions, and the single
#   canonical `paths` list used everywhere, and (c) loads the primary EcM
#   dataset into the object `emf` if it has already been produced.
#
#   It writes nothing except creating empty output sub-directories under
#   data_derived/ the first time it runs.
#
# WHY A SINGLE SETUP FILE
#   Keeping paths, constants, and helpers in one place means an analysis script
#   never hard-codes a file path or re-implements a shared rule (e.g. how a
#   "site" is defined). To point the pipeline at data held elsewhere, edit the
#   `paths` list below and nothing else.
#
# CONVENTIONS USED THROUGHOUT THE PIPELINE
#   * File paths are always built with here::here() from the project root
#     (the folder containing ECM_manuscript.Rproj). No absolute paths.
#   * Package functions are called with explicit namespaces (dplyr::filter(),
#     sf::st_read(), ...). The five library() calls below are the only
#     exceptions, kept for readability of the very common verbs.
# =============================================================================

# ---- Always-on packages -----------------------------------------------------
library(here)     # here::here() builds paths from the project root
library(dplyr)    # data manipulation verbs
library(tidyr)    # pivoting / reshaping
library(readr)    # fast CSV read/write
library(ggplot2)  # figures

# ---- Coordinate reference systems -------------------------------------------
# crs_wgs84  : lat/lon, used for all coordinate-based operations.
# crs_albers : Canada Albers Equal Area Conic, used for mapping (equal-area so
#              areas and point densities are not distorted).
crs_wgs84  <- "EPSG:4326"
crs_albers <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# ---- Figure 5 grey-background colour -----------------------------------------
# Figure 5 (manuscript.qmd) is a hand-assembled schematic composite that reuses
# panels from Figures 1, 3, 4, S1, S4, and S5 on a light-grey backdrop. The
# scripts producing those figures (10, 17, 18, and 20) each save TWO versions:
# the normal white-background copy used in the manuscript/SI (paths$fig_*),
# and an additional copy on this grey background (paths$fig_*_grey) for use
# only when manually assembling Figure 5.
fig5_grey_bg <- "#F2F2F2"

# ---- Multi-format figure export (journal submission) -------------------------
# save_fig_formats() wraps ggplot2::ggsave() to write the SAME figure in several
# file formats in one call. The pipeline names every figure with a .png path
# (paths$fig_*). FACETS (Canadian Science Publishing) does not accept PNG and
# instead asks for .tif / .jpg (among others), so each manuscript figure is
# written as PNG (unchanged) plus JPG and TIFF alongside it.
#
#   path    : canonical .png target (e.g. paths$fig_sampling_map). Its extension
#             is swapped to build the other formats in the same directory.
#   plot    : the ggplot / patchwork object to save.
#   formats : which formats to write. Default = PNG (keeps existing pipeline
#             behaviour) plus JPG and TIFF for submission.
#   ...     : forwarded verbatim to ggplot2::ggsave() (width, height, dpi, bg,
#             units, ...), so each call keeps the sizing it already specified.
#
# Format-specific defaults (overridable through ...):
#   JPG  -> quality = 100        (FACETS: max-quality JPEG, >= 300 dpi)
#   TIFF -> compression = "none" (FACETS: uncompressed TIFF)
save_fig_formats <- function(path, plot,
                             formats = c("png", "jpg", "tif"), ...) {
  stem <- tools::file_path_sans_ext(path)
  dots <- list(...)
  for (fmt in formats) {
    ext <- switch(fmt, jpeg = "jpg", tiff = "tif", fmt)           # normalize
    dev <- switch(ext, png = "png", jpg = "jpeg", tif = "tiff",
                  stop("save_fig_formats(): unsupported format '", fmt, "'"))
    args <- c(list(filename = paste0(stem, ".", ext),
                   plot = plot, device = dev), dots)
    if (ext == "jpg" && is.null(args$quality))     args$quality     <- 100
    if (ext == "tif" && is.null(args$compression)) args$compression <- "none"
    do.call(ggplot2::ggsave, args)
  }
  invisible(path)
}

# =============================================================================
# Helper functions
# =============================================================================

# ---- parse_gb_latlon(): parse a GenBank lat_lon string -----------------------
# GenBank stores coordinates as a single string like "45.1234 N 75.4321 W".
# Returns a one-row data.frame with numeric columns lat, lon (NA if unparseable).
parse_gb_latlon <- function(x) {
  result <- data.frame(lat = NA_real_, lon = NA_real_)
  if (is.na(x) || !nchar(trimws(x))) return(result)
  m <- regmatches(x, regexpr("([0-9.]+)\\s+([NS])\\s+([0-9.]+)\\s+([EW])", x))
  if (!length(m)) return(result)
  parts <- strsplit(trimws(m), "\\s+")[[1]]
  lat <- as.numeric(parts[1]) * ifelse(parts[2] == "S", -1, 1)
  lon <- as.numeric(parts[3]) * ifelse(parts[4] == "W", -1, 1)
  data.frame(lat = lat, lon = lon)
}

# ---- add_site_id(): define sampling "sites" by binning coordinates -----------
# A "site" is the bin produced by rounding lat/lon to 3 decimal places
# (~100 m at temperate latitudes). This collapses GPS-precision noise across
# repeat readings at one physical location without merging distinct sites.
# Adds columns site_lat, site_lon, and site (a "lat_lon" key).
# This is the canonical site rule; pass digits = 4 only if a finer bin is
# genuinely required.
add_site_id <- function(df, digits = 3L) {
  dplyr::mutate(df,
    site_lat = round(lat, digits),
    site_lon = round(lon, digits),
    site     = paste(site_lat, site_lon, sep = "_"))
}

# ---- snap_30arcsec(): snap a coordinate to the 30 arc-second (~1 km) grid -----
# 30 arc-seconds = 1/120 degree. Used by the Wallacean occupancy counts so that
# occurrences aggregate into ~1 km^2 grid cells (matching common SDM practice).
snap_30arcsec <- function(x) round(x * 120) / 120

# ---- canonicalize_host(): normalize a raw host-plant name --------------------
# Raw host strings come from the GenBank `host` qualifier and `isolation_source`
# field, and from GlobalFungi's `dominant_plant_species` / `other_plant_species`
# fields. This returns a clean "Genus species" (or "Genus") string, or NA when
# the input cannot be resolved to a Latin binomial. Steps, in order:
#   1. strip surrounding quotes/brackets and outer punctuation/whitespace
#   2. collapse internal whitespace
#   3. keep only the first two whitespace-separated tokens
#   4. strip residual per-token punctuation
#   5. apply the curated typo-correction lookup `typo_fixes`
#   6. case-normalize: Titlecase genus, lowercase species epithet (so that
#      "pinus banksiana" and "PINUS BANKSIANA" are both accepted)
#   7. reject anything that is not a valid Latin name pattern (common names such
#      as "pine", or entries with digits, become NA)
# Extend `typo_fixes` as new misspellings surface.
canonicalize_host <- function(x) {
  if (length(x) == 0L) return(character(0L))
  out <- trimws(as.character(x))
  out <- gsub('^["‘’“”\\(\\[]+', "", out, perl = TRUE)
  out <- gsub('["‘’“”\\)\\]]+$',  "", out, perl = TRUE)
  out <- gsub("^[[:punct:][:space:]]+|[[:punct:][:space:]]+$", "", out)
  out <- gsub("\\s+", " ", out)
  out <- trimws(out)
  empty_idx <- !is.na(out) & !nzchar(out)
  out[empty_idx] <- NA_character_

  toks <- strsplit(out, "\\s+")
  out <- vapply(seq_along(toks), function(i) {
    if (is.na(out[i])) return(NA_character_)
    p <- toks[[i]]
    if (length(p) == 0L) return(NA_character_)
    p <- gsub("^[[:punct:]]+|[[:punct:]]+$", "", p)
    p <- p[nzchar(p)]
    if (length(p) >= 2L) paste(p[1L], p[2L])
    else if (length(p) == 1L) p[1L]
    else NA_character_
  }, character(1L))

  typo_fixes <- c("Abies balamifera" = "Abies balsamea")
  ix <- match(out, names(typo_fixes))
  out[!is.na(ix)] <- unname(typo_fixes[ix[!is.na(ix)]])

  toks2 <- strsplit(out, " ", fixed = TRUE)
  out <- vapply(seq_along(toks2), function(i) {
    if (is.na(out[i])) return(NA_character_)
    p <- toks2[[i]]
    if (length(p) == 0L) return(NA_character_)
    g <- paste0(toupper(substring(p[1L], 1L, 1L)), tolower(substring(p[1L], 2L)))
    if (length(p) >= 2L) paste(g, tolower(p[2L])) else g
  }, character(1L))

  valid <- grepl("^[A-Z][a-z]+( [a-z][a-z-]*)?$", out)
  out[!valid] <- NA_character_
  out
}

# =============================================================================
# Batch-download safety helpers
# =============================================================================
# Several scripts retrieve data from external services (NCBI, BIEN, NSR, GBIF,
# GIFT) in batches, because a single request for thousands of records times out.
# The original pattern was:
#
#     tryCatch(fetch(batch), error = function(e) { warning(...); NULL })
#
# which SKIPS a failed batch. Combined with a checkpoint written regardless and
# an if (!file.exists(...)) guard on the step, that turns a transient network
# error into permanent, invisible data loss: nothing compares what came back
# against what was asked for, the truncated file looks complete, and every later
# run reuses it. That is exactly how 3,001 of 60,911 GenBank records went missing
# (see repro/genbank_fetch_gap_prefix_state.md).
#
# The two helpers below are the shared remedy. Scripts use them together:
# fetch_with_retry() gives a flaky batch several chances, and
# assert_fetch_complete() refuses to let the step continue if records are still
# missing afterwards.

# ---- fetch_with_retry(): run a fetch, retrying transient failures ------------
# Calls fetch_fun() up to `max_attempts` times, pausing between tries (the pause
# is multiplied by the attempt number, so it backs off politely if the server is
# rate-limiting rather than merely glitching).
#
#   fetch_fun    : a zero-argument function that performs the request
#   what         : short label for messages, e.g. "NSR native status"
#   batch_i      : which batch this is (for messages); NULL if not batched
#   n_batches    : how many batches in total (for messages)
#
# Returns whatever fetch_fun() returned on the first success, or NULL if every
# attempt failed. NULL means "this batch is lost" — the CALLER must count those
# and pass the count to assert_fetch_complete(). Never treat NULL as "no data".
fetch_with_retry <- function(fetch_fun, what, batch_i = NULL, n_batches = NULL,
                             max_attempts = 5L, pause = 3) {
  label <- if (is.null(batch_i)) what else
    sprintf("%s batch %d/%d", what, batch_i, n_batches)
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(fetch_fun(), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (attempt < max_attempts) {
      message(sprintf("  %s failed (attempt %d of %d): %s\n    retrying in %g s",
                      label, attempt, max_attempts,
                      conditionMessage(result), pause * attempt))
      Sys.sleep(pause * attempt)
    } else {
      warning(sprintf("%s failed on all %d attempts: %s", label, max_attempts,
                      conditionMessage(result)), call. = FALSE)
    }
  }
  NULL
}

# ---- assert_fetch_complete(): stop unless the step got what it asked for -----
# Compares records actually returned against records requested and stops with an
# informative error if the step came up short. Two independent trip conditions:
#   (a) any batch exhausted its retries — a known, located loss;
#   (b) the shortfall exceeds max_missing_frac — catches partial returns that
#       raised no error at all.
#
# max_missing_frac is not zero because a few records can legitimately disappear
# between listing and retrieval. Set it far below the size of one batch, so that
# losing even a single whole batch trips the guard.
#
# Call this BEFORE writing the checkpoint, so a short result is never persisted.
assert_fetch_complete <- function(n_returned, n_requested, n_failed_batches,
                                  what, max_missing_frac = 0.001,
                                  partial_path = NULL) {
  shortfall <- n_requested - n_returned
  frac      <- if (n_requested == 0L) 0 else shortfall / n_requested
  if (n_failed_batches == 0L && frac <= max_missing_frac) {
    return(invisible(TRUE))
  }
  stop(sprintf(
    paste0(
      "%s is incomplete: requested %s, received %s (shortfall %s, %.3f%%; ",
      "tolerance %.3f%%). Batches that failed every retry: %d.\n",
      "  Refusing to write a truncated checkpoint — a short file here would be ",
      "reused by every later run and the loss would become permanent.\n",
      "  This is usually a transient service outage: wait, then re-run.%s"
    ),
    what, format(n_requested, big.mark = ","), format(n_returned, big.mark = ","),
    format(shortfall, big.mark = ","), 100 * frac, 100 * max_missing_frac,
    n_failed_batches,
    if (is.null(partial_path)) "" else
      paste0("\n  Partial output left at: ", partial_path)
  ), call. = FALSE)
}

# ---- write_checkpoint_atomically(): write, verify, then promote --------------
# Writes to a ".partial" sibling and renames it over the target only after
# `write_fun` has returned without error. A run that dies part-way therefore
# leaves the ".partial" file (obviously not a checkpoint) rather than a
# half-written file that the next run's file.exists() test accepts as complete.
#   path      : final checkpoint path
#   write_fun : function(tmp_path) that writes the content to tmp_path
write_checkpoint_atomically <- function(path, write_fun) {
  tmp <- paste0(path, ".partial")
  write_fun(tmp)
  if (!file.exists(tmp))
    stop("write_checkpoint_atomically(): nothing was written to ", tmp,
         call. = FALSE)
  if (!file.rename(tmp, path))
    stop("write_checkpoint_atomically(): could not move ", tmp, " to ", path,
         call. = FALSE)
  invisible(path)
}

# ---- EcM definition: classify_ecm() and its reference tables ---------------
# One definition of "ectomycorrhizal" (EcM), used by 04_combine_ecm_dataset.R
# for our sequence records and by 09_linnean.R for GBIF specimen genera.
#   1. Genus rule: the genus is listed as EcM (primary_lifestyle) in
#      FungalTraits v1.2.
#   2. Family rule: if the genus is absent from FungalTraits, or there is no
#      genus label, the record is EcM when at least FAMILY_ECM_MIN_FRAC (two
#      thirds) of the genera FungalTraits lists in its family are EcM.
# A genus that FungalTraits lists as NOT EcM is never admitted.
FAMILY_ECM_MIN_FRAC <- 2 / 3

# ecm_reference(): read FungalTraits once and return what classify_ecm() needs:
# the EcM genera, all FungalTraits genera, and a per-family summary (number of
# genera, number EcM, whether the family passes, and the shared EcM lineage
# if all EcM genera in the family have the same one).
ecm_reference <- function(ft_path = paths$fungaltraits) {
  ft <- readr::read_csv(ft_path, show_col_types = FALSE) |>
    dplyr::rename_with(tolower) |>
    dplyr::mutate(genus_lower = tolower(trimws(genus)),
                  is_ecm      = primary_lifestyle %in% "ectomycorrhizal")
  family <- ft |>
    dplyr::filter(!is.na(family)) |>
    dplyr::group_by(family) |>
    dplyr::summarise(
      ft_genera      = dplyr::n_distinct(genus_lower),
      ft_ecm_genera  = dplyr::n_distinct(genus_lower[is_ecm]),
      family_lineage = {
        l <- unique(stats::na.omit(ectomycorrhiza_lineage_template[is_ecm]))
        if (length(l) == 1L) l else NA_character_
      },
      .groups = "drop") |>
    dplyr::mutate(frac_ecm      = ft_ecm_genera / ft_genera,
                  family_passes = frac_ecm >= FAMILY_ECM_MIN_FRAC)
  list(ft = ft,
       em_genera = unique(ft$genus_lower[ft$is_ecm]),
       ft_genera = unique(ft$genus_lower),
       family    = family)
}

# has_genus_label(): FALSE for a missing genus or a UNITE placeholder such as
# "Boletaceae_gen_Incertae_sedis" (an SH placed only to family).
has_genus_label <- function(genus) {
  !is.na(genus) & nzchar(trimws(genus)) & !grepl("Incertae_sedis", genus)
}

# classify_ecm(): for vectors of genus and family names, return one of
#   "genus"                  EcM by the genus rule
#   "family"                 EcM by the family rule
#   "FungalTraits: not EcM"  genus listed in FungalTraits, not as EcM
#   "family rule: refused"   family rule applied, family below the threshold
#                            or not in FungalTraits
classify_ecm <- function(genus, family, ref) {
  g   <- tolower(trimws(genus))
  fam <- sub("_.*$", "", family)               # "Boletaceae_fam..." -> "Boletaceae"
  lab <- has_genus_label(genus)
  passes <- fam %in% ref$family$family[ref$family$family_passes]
  dplyr::case_when(
    lab & g %in% ref$em_genera ~ "genus",
    lab & g %in% ref$ft_genera ~ "FungalTraits: not EcM",
    passes                     ~ "family",
    TRUE                       ~ "family rule: refused")
}

# ---- coords_in_canada(): is each point inside Canada? ------------------------
# TRUE if the point lies inside the Canada boundary (paths$canada_bound, built
# in 01_spatial_data.R and simplified at 1 km) buffered by 2 km; FALSE if
# outside; NA if lat or lon is missing. The buffer is applied in the
# equal-area projection, so it is 2 km everywhere. Used by
# 04_combine_ecm_dataset.R (coord_in_canada) and 09_linnean.R (Table S4).
coords_in_canada <- function(lat, lon, buffer_m = 2000) {
  out <- rep(NA, length(lat))
  ok  <- which(!is.na(lat) & !is.na(lon))
  if (length(ok) == 0L) return(out)
  boundary <- sf::st_read(paths$canada_bound, quiet = TRUE) |>
    sf::st_transform(crs_albers) |>
    sf::st_buffer(dist = buffer_m) |>
    sf::st_transform(4326)
  pts <- sf::st_as_sf(data.frame(lon = lon[ok], lat = lat[ok]),
                      coords = c("lon", "lat"), crs = 4326)
  old <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old), add = TRUE)
  out[ok] <- lengths(sf::st_intersects(pts, boundary)) > 0
  out
}

# ---- read_fasta_tbl() / write_fasta_tbl(): FASTA <-> table --------------------
# read: one row per sequence, `header` (without ">") and `seq` (joined across
# wrapped lines); an empty table if the file is missing or has no sequences.
read_fasta_tbl <- function(path) {
  empty <- tibble::tibble(header = character(0), seq = character(0))
  if (!file.exists(path)) return(empty)
  lines <- readr::read_lines(path)
  idx <- which(startsWith(lines, ">"))
  if (length(idx) == 0L) return(empty)
  ends <- c(idx[-1] - 1L, length(lines))
  tibble::tibble(
    header = sub("^>", "", lines[idx]),
    seq    = vapply(seq_along(idx), function(i)
      if (ends[i] > idx[i]) paste(lines[(idx[i] + 1L):ends[i]], collapse = "") else "",
      character(1)))
}
write_fasta_tbl <- function(tbl, path) {
  readr::write_lines(paste0(">", tbl$header, "\n", tbl$seq), path)
}

# ---- assign_its2_to_sh(): ITSx + vsearch SH assignment -------------------------
# The rule used for our own GenBank data (03_genbank.R) and for the genome
# sequences (15b): extract ITS2 with ITSx, drop fragments shorter than
# `min_len`, match against the pinned UNITE FASTA with vsearch at `identity`,
# and keep every reference tied for the best match.
# Ties are resolved as in 03: each tied hit is reduced to a "tie identity" (the
# species name for a named species, the SH code for an un-named "_sp" SH). If
# all tied hits share one identity, the sequence is credited to every SH code
# among them; otherwise it is credited to none ("ambiguous_tie").
# `seqs` needs columns `key` (a short unique label) and `seq`. Returns one row
# per key: outcome ("sh_assigned", "ambiguous_tie", "no_unite_hit",
# "no_its2"), `sh_codes` ("; "-separated; only for "sh_assigned"),
# `unite_names` (UNITE names of the tied best hits) and best identity.
assign_its2_to_sh <- function(seqs, work_prefix, identity = 0.985,
                              min_len = 100L, threads = 4L) {
  itsx_bin    <- Sys.which("ITSx")
  vsearch_bin <- Sys.which("vsearch")
  if (!nzchar(itsx_bin))    stop("ITSx not found on PATH (see 03_genbank.R header).")
  if (!nzchar(vsearch_bin)) stop("vsearch not found on PATH (see 03_genbank.R header).")
  input  <- paste0(work_prefix, "_input.fasta")
  prefix <- paste0(work_prefix, "_itsx")
  query  <- paste0(work_prefix, "_vsearch_query.fasta")
  hitsf  <- paste0(work_prefix, "_vsearch_hits.txt")
  # remove outputs of any earlier run, so stale files are never read
  unlink(c(input, Sys.glob(paste0(prefix, "*")), query, hitsf))

  write_fasta_tbl(dplyr::transmute(seqs, header = key, seq), input)
  ret <- system(paste(shQuote(itsx_bin), "-i", shQuote(input), "-o", shQuote(prefix),
                      "-t F", "--save_regions ITS2", "--cpu", threads,
                      "--graphical F", "--silent T"))
  if (ret != 0L) stop("ITSx failed with exit code ", ret)
  its2 <- read_fasta_tbl(paste0(prefix, ".ITS2.fasta")) |>
    dplyr::transmute(key = sub("\\|.*$", "", header),
                     seq = gsub("[^ACGTUNRYSWKMBDHV]", "N", toupper(seq))) |>
    dplyr::filter(nchar(seq) >= min_len)

  hits <- tibble::tibble(key = character(0), target_id = character(0),
                         pident = numeric(0))
  if (nrow(its2) > 0L) {
    write_fasta_tbl(dplyr::rename(its2, header = key), query)
    ret <- system(paste(shQuote(vsearch_bin), "--usearch_global", shQuote(query),
                        "--db", shQuote(paths$unite_fasta), "--id", identity,
                        "--strand both", "--maxaccepts 0", "--maxrejects 0",
                        "--blast6out", shQuote(hitsf), "--top_hits_only",
                        "--threads", threads, "--quiet"))
    if (ret != 0L) stop("vsearch failed with exit code ", ret)
    if (file.exists(hitsf) && file.size(hitsf) > 0)
      hits <- readr::read_tsv(hitsf, col_names = FALSE, show_col_types = FALSE,
                              col_types = readr::cols(.default = "c")) |>
        dplyr::transmute(key = X1, target_id = X2, pident = as.numeric(X3))
  }
  per_key <- hits |>
    dplyr::mutate(sh_code  = sub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", target_id),
                  unite_sp = sub("\\|.*$", "", target_id),
                  tie_id   = dplyr::if_else(grepl("_sp$", unite_sp), sh_code, unite_sp)) |>
    dplyr::group_by(key) |>
    dplyr::summarise(best_identity = max(pident),
                     n_tie_ids     = dplyr::n_distinct(tie_id),
                     sh_codes      = paste(sort(unique(sh_code)), collapse = "; "),
                     unite_names   = paste(sort(unique(unite_sp)), collapse = "; "),
                     .groups = "drop")
  dplyr::select(seqs, key) |>
    dplyr::left_join(per_key, by = "key") |>
    dplyr::mutate(
      outcome = dplyr::case_when(
        !key %in% its2$key ~ "no_its2",
        is.na(n_tie_ids)   ~ "no_unite_hit",
        n_tie_ids > 1L     ~ "ambiguous_tie",
        TRUE               ~ "sh_assigned"),
      sh_codes = dplyr::if_else(outcome == "sh_assigned", sh_codes, NA_character_)) |>
    dplyr::select(key, outcome, sh_codes, unite_names, best_identity)
}

# ---- read_big_tsv_subset(): read selected columns from a very large TSV -------
# data.table::fread(..., select = ) still has to stream and tokenize the whole
# file, which overruns R's 2^31-1-byte single-string limit on multi-GB inputs
# such as the ~13 GB GlobalFungi SH abundance matrix. This helper instead streams
# the file once through `awk`, writing ONLY the requested columns to a temporary
# TSV, then fread()s that much smaller file.
#   path : path to the large tab-separated file (with a header row)
#   cols : column NAMES to keep (first is typically the id column, e.g. sample_ID);
#          names not present in the file's header are dropped
# Returns a data.table with the kept columns, in the order given.
read_big_tsv_subset <- function(path, cols) {
  header <- names(data.table::fread(path, sep = "\t", quote = "", nrows = 0L))
  keep   <- cols[cols %in% header]
  if (length(keep) == 0L)
    stop("None of the requested columns are present in ", basename(path))
  idx <- match(keep, header)                          # 1-based column positions
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  awk_prog <- paste0(
    'BEGIN{FS=OFS="\\t"; n=split(cols,a,",")} ',
    '{out=$(a[1]); for(i=2;i<=n;i++) out=out OFS $(a[i]); print out}'
  )
  cmd <- sprintf("awk -v cols=%s %s %s > %s",
                 shQuote(paste(idx, collapse = ",")),
                 shQuote(awk_prog), shQuote(path), shQuote(tmp))
  if (system(cmd) != 0L || !file.exists(tmp))
    stop("awk column subset failed for ", basename(path))
  data.table::fread(tmp, sep = "\t", quote = "")
}

# =============================================================================
# Canonical paths
# =============================================================================
# Two data roots:
#   data_raw/      read-only inputs (obtained from external sources; see
#                  data_raw/DATA-DICTIONARY.md). Never written to by the pipeline.
#   data_derived/  everything the pipeline produces. Starts empty; running the
#                  scripts in order recreates it (see data_derived/DATA-DICTIONARY.md).
# Figures go to figures/ with manuscript file names (Figure-xx_description.png).
paths <- list(

  # ---- directory roots -------------------------------------------------------
  data_raw     = here::here("data_raw"),
  data_derived = here::here("data_derived"),
  figures      = here::here("figures"),
  checkpoints  = here::here("data_derived", "checkpoints"),
  # ---- raw: live-query snapshots ---------------------------------------------
  # Results of queries to services with no dated releases (GenBank, GBIF name
  # matching, BIEN/NSR, GIFT). They cannot be downloaded again in identical
  # form, so they are raw data: kept in data_raw/snapshots/ and deposited.
  # Each is written by its script on a first run (when absent) and only read
  # afterwards. To refresh one, delete it and re-run the script that writes it.
  snapshots           = here::here("data_raw", "snapshots"),
  snap_gb_ids         = here::here("data_raw", "snapshots", "genbank_emf_canada_ids.txt"),           # 03
  snap_gb_fasta       = here::here("data_raw", "snapshots", "genbank_emf_canada.fasta"),             # 03
  snap_gb_meta        = here::here("data_raw", "snapshots", "genbank_emf_canada_metadata.csv"),      # 03
  snap_gb_nosequence  = here::here("data_raw", "snapshots", "genbank_no_sequence_available.csv"),    # 03
  snap_gb_fetch_log   = here::here("data_raw", "snapshots", "genbank_fetch_log.txt"),                # 03
  snap_fr_backbone    = here::here("data_raw", "snapshots", "fungalroot_species_backbone.csv"),      # 05
  snap_bien_native    = here::here("data_raw", "snapshots", "bien_nsr_native_species.csv"),          # 06
  snap_bien_growth    = here::here("data_raw", "snapshots", "bien_ecm_growthforms.csv"),             # 06
  snap_bien_growth_queried = here::here("data_raw", "snapshots", "bien_ecm_growthforms_queried.csv"),# 06
  snap_gift_growth    = here::here("data_raw", "snapshots", "gift_growthforms.csv"),                 # 06
  snap_gb_global_meta = here::here("data_raw", "snapshots", "genbank_global_ecm_meta.csv"),          # 18
  snap_gb_global_queried = here::here("data_raw", "snapshots", "genbank_global_ecm_genera_queried.csv"), # 18
  # MycoCosm genome -> same-strain GenBank ITS link (15b). Also live-query
  # results, but kept in their own writable folder because 15b adds rows for
  # newly listed genomes.
  mycocosm_its_dir    = here::here("data_raw", "mycocosm_its"),
  mycocosm_its_lookup = here::here("data_raw", "mycocosm_its", "mycocosm_its_lookup.csv"),
  mycocosm_its_seqs   = here::here("data_raw", "mycocosm_its", "mycocosm_its_sequences.fasta"),
  temp_dir     = here::here("data_derived", "checkpoints"),

  # ---- per-shortfall output directories --------------------------------------
  out_linnean       = here::here("data_derived", "linnean"),
  out_wallacean     = here::here("data_derived", "wallacean"),
  out_prestonian    = here::here("data_derived", "prestonian"),
  out_darwinian     = here::here("data_derived", "darwinian"),
  # SH-level genome coverage (15b_darwinian_sh_match.R)
  darwinian_genome_sh   = here::here("data_derived", "darwinian", "darwinian_genome_sh_assignments.csv"),
  darwinian_sh_coverage = here::here("data_derived", "darwinian", "darwinian_sh_coverage.csv"),
  darwinian_sh_summary  = here::here("data_derived", "darwinian", "darwinian_sh_summary.csv"),
  darwinian_not_assessed = here::here("data_derived", "darwinian", "darwinian_genomes_not_assessed.csv"),
  out_raunkiaeran   = here::here("data_derived", "raunkiaeran"),
  out_hutchinsonian = here::here("data_derived", "hutchinsonian"),
  out_eltonian      = here::here("data_derived", "eltonian"),
  # Geographic-possibility sensitivity analysis (21_geographic_possibility.R),
  # a supplementary check on the Eltonian interaction-matrix denominator.
  out_geo_possibility = here::here("data_derived", "geo_possibility"),

  # ---- raw: spatial ----------------------------------------------------------
  canada_bound_raw = here::here("data_raw", "admin_boundaries", "canada_bound_raw.gpkg"),
  ecoregions_raw   = here::here("data_raw", "ecoregions", "Ecoregions", "ecoregions.shp"),
  ne_canada_raw    = here::here("data_raw", "natural_earth", "canada_ne.gpkg"),
  ne_lakes_raw     = here::here("data_raw", "natural_earth", "lakes_ne.gpkg"),

  # ---- raw: sequence sources -------------------------------------------------
  gf_metadata     = here::here("data_raw", "GlobalFungi", "GlobalFungi_5_sample_metadata.txt"),
  gf_sh_abundance = here::here("data_raw", "GlobalFungi", "GlobalFungi_5_SH_abundance_ITS1_ITS2.txt"),
  # UNITE reference FASTA — pinned to the 2024-04-04 build (see
  # data_raw/DATA-DICTIONARY.md for why this build, not the newest, is used).
  # GlobalFungi's pre-assigned SH codes must be decoded against this same build,
  # and GenBank SH codes are assigned in-house by vsearch against it, so both
  # roles deliberately point to one file.
  unite_fasta = here::here("data_raw", "UNITE", "sh_general_release_dynamic_04.04.2024_dev.fasta"),

  # ---- raw: trait / genome / climate reference -------------------------------
  # FungalTraits v1.2 (Põlme et al. 2020): pinned reference, as with UNITE above.
  # primary_lifestyle == "ectomycorrhizal" defines EcM dataset membership, so the
  # reference is held fixed rather than tracked to a live source. To re-pin,
  # replace the file, update data_raw/fungaltraits/fungaltraits_version.txt, and
  # re-validate the EcM genus count.
  fungaltraits   = here::here("data_raw", "fungaltraits", "FungalTraits_1-2.csv"),
  mycocosm_list  = here::here("data_raw", "mycocosm", "mycocosm_organism_list.csv"),
  climate_raster = here::here("data_raw", "climate", "wc2.1_country", "CAN_wc2.1_30s_bio.tif"),
  biotime_db     = here::here("data_raw", "biotime"),
  # van Galen et al. (2025) per-sample proportion of unassigned EcM OTUs.
  # Only this small CSV is needed (for the Linnean per-sample dark-fraction
  # comparison); the large dark-taxa GeoTIFFs are NOT part of this project.
  van_galen_per_sample = here::here("data_raw", "van_Galen_per_sample",
                                    "GFv5_EcM_unassigned_per_sample.csv"),
  # van Galen et al. (2025) dark-taxa richness GeoTIFF (band percentage_dark_taxa
  # is used for the dark-diversity map, Figure S1).
  van_galen_tif = here::here("data_raw", "van_Galen_et_al_dark_taxa_code_and_data",
                             "4.Dark_EcM_taxa_richness_maps",
                             "Dark_taxa_geospatial_layers.tif"),

  # ---- raw: BIEN2 modelled range shapefiles (from 07_bien2_ranges.R) ---------
  bien2_ranges_dir = here::here("data_raw", "bien2_ranges"),

  # ---- derived: processed spatial (from 01_spatial_data.R) -------------------
  canada_bound         = here::here("data_derived", "spatial", "canada_simple.gpkg"),
  lakes_albers         = here::here("data_derived", "spatial", "lakes_canada_albers.gpkg"),
  ecoregions_processed = here::here("data_derived", "spatial", "ecoregions_processed.gpkg"),
  ecozone_names        = here::here("data_derived", "spatial", "ecozone_names.csv"),

  # ---- derived: sequence pipeline checkpoints & outputs ----------------------
  gf_ids_out       = here::here("data_derived", "checkpoints", "globalfungi_canada_ids.txt"),
  gf_sh_subset_out = here::here("data_derived", "checkpoints", "globalfungi_canada_SH_abundance.txt"),
  gf_meta_out      = here::here("data_derived", "checkpoints", "globalfungi_canada_metadata.csv"),
  unite_taxonomy   = here::here("data_derived", "checkpoints", "unite_sh_taxonomy.csv"),
  gf_sh_unmatched  = here::here("data_derived", "checkpoints", "globalfungi_sh_unmatched.csv"),
  # Sample-type composition of GlobalFungi v5, worldwide and Canada-wide
  # (02_globalfungi.R, Step 2b). Read by the manuscript Methods.
  gf_sample_type_tally = here::here("data_derived", "globalfungi_sample_type_tally.csv"),
  # Canadian samples by sample type, kept or dropped by the sample-type
  # inclusion list (02_globalfungi.R, Step 2).
  gf_sample_type_filter = here::here("data_derived", "globalfungi_canada_sample_type_filter.csv"),
  gb_fasta_out    = here::here("data_raw", "snapshots", "genbank_emf_canada.fasta"),
  gf_long_out      = here::here("data_derived", "globalfungi_canada_long.csv"),
  gb_long_out      = here::here("data_derived", "genbank_emf_canada_long.csv"),
  emf_combined     = here::here("data_derived", "emf_canada_combined.csv"),
  # GenBank records dropped by the location rule in 04_combine_ecm_dataset.R,
  # with the reason for each.
  genbank_excluded = here::here("data_derived", "genbank_excluded_non_canadian.csv"),
  # Families considered by the EcM family rule in 04_combine_ecm_dataset.R.
  ecm_family_rule_table = here::here("data_derived", "ecm_family_rule_table.csv"),
  emf_data         = here::here("data_derived", "emf_canada_em_only.csv"),

  # ---- derived: host reference (from 05_prepare_fungalroot.R) ----------------
  fungalroot_sp     = here::here("data_derived", "clean_fungalroot_species.csv"),
  fungalroot_genera = here::here("data_derived", "clean_fungalroot_genera_table_s2.csv"),
  host_species      = here::here("data_derived", "ecm_native_canada_host_species.csv"),

  # ---- derived: BIEN2 host rasters (from 08_host_rasters.R) -------------------
  bien_ranges        = here::here("data_derived", "checkpoints", "bien2_ecm_host_ranges.gpkg"),
  bien_richness      = here::here("data_derived", "spatial", "bien_host_richness_0.5deg.tif"),
  bien_species_stack = here::here("data_derived", "spatial", "bien_host_species_stack.tif"),
  bien_data_rich     = here::here("data_derived", "spatial", "bien_host_data_richness_0.5deg.tif"),
  bien_proportion    = here::here("data_derived", "spatial", "bien_host_data_proportion_0.5deg.tif"),

  # ---- derived: cross-referenced shortfall outputs ---------------------------
  gbif_ecm            = here::here("data_derived", "linnean", "linnean_gbif_ecm_canada.csv"),
  gbif_ecm_nosequence = here::here("data_derived", "linnean", "linnean_gbif_ecm_nosequence_canada.csv"),
  gbif_plotted_counts = here::here("data_derived", "linnean", "linnean_gbif_plotted_counts.csv"),
  # Table S4 by sequence matching (09_linnean.R, Step 6b). Raw GenBank
  # downloads for each candidate genus are cached in data_raw/ (they cannot be
  # reproduced exactly later); the per-genus result table feeds Table S4.
  table_s4_genbank_dir = here::here("data_raw", "genbank_genus_queries"),
  table_s4             = here::here("data_derived", "linnean", "linnean_table_s4.csv"),
  table_s4_sequences   = here::here("data_derived", "linnean", "linnean_table_s4_sequence_matches.csv"),
  # SUPERSEDED (09b_table_s4_check.R, retired 2026-11-14); kept until the SM1
  # Table S4 code is rewritten in Phase 5.
  # Table S4 check (09b_table_s4_check.R): hand-curated synonym table and the
  # cached GenBank query results (raw inputs), plus the per-genus result table
  # and the GBIF records of genera confirmed to lack Canadian ITS sequences.
  table_s4_synonyms       = here::here("data_raw", "table_s4_check", "table_s4_genus_synonyms.csv"),
  table_s4_ncbi_dir       = here::here("data_raw", "table_s4_check", "ncbi_genus_queries"),
  table_s4_reconciliation = here::here("data_derived", "linnean", "linnean_table_s4_reconciliation.csv"),
  gbif_ecm_nosequence_confirmed = here::here("data_derived", "linnean", "linnean_gbif_ecm_nosequence_confirmed_canada.csv"),
  eltonian_global        = here::here("data_derived", "eltonian", "eltonian_global_host_associations.csv"),
  eltonian_host_matching = here::here("data_derived", "eltonian", "eltonian_host_matching.csv"),
  prestonian_out      = here::here("data_derived", "prestonian", "prestonian_biotime_matches.csv"),
  darwinian_out       = here::here("data_derived", "darwinian", "darwinian_mycocosm_matches.csv"),

  # ---- figures (manuscript file names) ---------------------------------------
  # Main text
  fig_sampling_map      = here::here("figures", "Figure-01_sampling_map.png"),
  fig_sampling_map_grey = here::here("figures", "Figure-01_sampling_map_grey.png"),
  fig_wallacean_occ    = here::here("figures", "Figure-02_wallacean_occupancy.png"),
  fig_climate_gap      = here::here("figures", "Figure-03_climate_gap.png"),
  fig_climate_gap_grey = here::here("figures", "Figure-03_climate_gap_grey.png"),
  fig_host_bivariate      = here::here("figures", "Figure-04_host_bivariate_map.png"),
  fig_host_bivariate_grey = here::here("figures", "Figure-04_host_bivariate_map_grey.png"),
  fig_shortfalls_summary  = here::here("figures", "Figure-05_shortfalls_summary.png"),
  # Supplemental
  fig_dark_diversity      = here::here("figures", "Figure-S1_dark_diversity.png"),
  fig_dark_diversity_grey = here::here("figures", "Figure-S1_dark_diversity_grey.png"),
  fig_density_world    = here::here("figures", "Figure-S2_gf_sampling_density_world.png"),
  fig_depth_discard    = here::here("figures", "Figure-S3_depth_discard.png"),
  fig_ecozone_sampling      = here::here("figures", "Figure-S4_ecozone_sampling_map.png"),
  fig_ecozone_sampling_grey = here::here("figures", "Figure-S4_ecozone_sampling_map_grey.png"),
  fig_gbif_specimens      = here::here("figures", "Figure-S5_gbif_specimens.png"),
  fig_gbif_specimens_grey = here::here("figures", "Figure-S5_gbif_specimens_grey.png")
)

# ---- Create empty data_derived/ output sub-directories the first time --------
for (d in c(paths$checkpoints, paths$figures, paths$biotime_db,
            paths$out_linnean, paths$out_wallacean, paths$out_prestonian,
            paths$out_darwinian, paths$out_raunkiaeran, paths$out_hutchinsonian,
            paths$out_eltonian, paths$out_geo_possibility,
            here::here("data_derived", "spatial"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Load the primary EcM dataset if it exists ------------------------------
# `emf` is the combined, EcM-filtered GlobalFungi + GenBank record table
# produced by 04_combine_ecm_dataset.R. Group 0 scripts that build it run before
# it exists, so loading is conditional. The na and col_types arguments below
# suppress spurious readr parsing warnings (literal "NULL" year strings from
# GenBank, and the rarely-populated secondary_lifestyle column).
if (file.exists(paths$emf_data)) {
  emf <- readr::read_csv(paths$emf_data, show_col_types = FALSE,
                         na = c("", "NA", "NULL"),
                         col_types = readr::cols(
                           secondary_lifestyle = readr::col_character(),
                           .default = readr::col_guess()
                         )) |>
    dplyr::filter(source %in% c("GenBank", "GlobalFungi"))
}
