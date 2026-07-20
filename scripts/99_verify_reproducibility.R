# =============================================================================
# 99_verify_reproducibility.R  —  Verify the pipeline reproduces its own results
# =============================================================================
# WHAT THIS DOES
#   Captures a fingerprint of every pipeline output (file hashes + the headline
#   numeric metrics that the manuscript reports), then — after you re-run the
#   pipeline from scratch — compares the new outputs against that fingerprint
#   and reports what changed.
#
#   This tests COMPUTATIONAL reproducibility: that your code, run on your
#   archived raw data, regenerates your reported results. It does NOT test data
#   acquisition (the network-dependent steps still call GBIF/BIEN/NSR/Entrez,
#   which can legitimately return different data on a different day).
#
# HOW TO RUN
#   Step A — before deleting anything, capture the baseline:
#       Rscript scripts/99_verify_reproducibility.R baseline
#   Step B — re-run the pipeline (see WORKFLOW below), then:
#       Rscript scripts/99_verify_reproducibility.R verify
#
# WORKFLOW (the full check)
#   1. Rscript scripts/99_verify_reproducibility.R baseline
#   2. Back up data_derived/ and figures/ (belt and braces)
#   3. Delete data_derived/ and figures/  <- the cold-cache run is the real test
#   4. Rscript scripts/run_all.R           (redirect stderr: 2> repro/run_stderr.log)
#   5. Rscript scripts/99_verify_reproducibility.R verify
#   6. Re-render the manuscript + supplements and diff the text (see README note)
#
# WHERE OUTPUT GOES
#   Everything is written to repro/ at the PROJECT ROOT — deliberately NOT
#   inside data_derived/, because step 3 above deletes data_derived/ and would
#   otherwise destroy the baseline you are comparing against.
#
#     repro/baseline_files.csv    — path, bytes, md5 for every tracked output
#     repro/baseline_metrics.csv  — headline numeric metrics (see below)
#     repro/baseline_env.csv      — R version, platform, key package versions
#     repro/verify_files.csv      — same three, captured at verify time
#     repro/verify_metrics.csv
#     repro/verify_env.csv
#     repro/verify_report.txt     — the human-readable comparison
#
# WHAT COUNTS AS A METRIC
#   Every (metric, value) pair in the per-shortfall *_summary.csv files — these
#   are the numbers the manuscript quotes — plus a few structural counts (row
#   counts of the main tables, the host-list evidence_source breakdown, and the
#   per-source unique sampling locations behind Figure 1). Metrics are captured
#   automatically, so this file needs no editing when a summary gains a row.
#
# OPTIONS (environment variables)
#   ECM_REPRO_CHECKPOINTS=true   also fingerprint data_derived/checkpoints/
#                                (~920 MB of API/parse caches). Off by default:
#                                caches are regenerated from live services, so
#                                they differ legitimately and add noise.
# =============================================================================

# ---- Configuration -----------------------------------------------------------

REPRO_DIR <- here::here("repro")

# Directories whose contents are treated as pipeline OUTPUT.
TRACK_DIRS <- c(here::here("data_derived"), here::here("figures"))

# Caches, run records and this tool's own output are not results.
INCLUDE_CHECKPOINTS <- tolower(Sys.getenv("ECM_REPRO_CHECKPOINTS", "false")) %in%
  c("true", "1", "yes")

# Byte-identical comparison is not meaningful for these: image and container
# formats embed creation timestamps, so they differ on every write even when
# the content is identical. They are compared on size and flagged for eyeball.
BINARY_EXT <- c("png", "jpg", "jpeg", "tif", "tiff", "pdf", "gpkg", "rds")

# ---- Helpers -----------------------------------------------------------------

# Collect every tracked output file, with its size and MD5 hash.
fingerprint_files <- function() {
  files <- unlist(lapply(TRACK_DIRS, function(d)
    if (dir.exists(d)) list.files(d, recursive = TRUE, full.names = TRUE) else character(0)))

  drop <- grepl("(^|/)\\.DS_Store$", files) |
          grepl("(^|/)run_log\\.csv$", files) |
          grepl("/repro/", files)
  if (!INCLUDE_CHECKPOINTS) drop <- drop | grepl("/checkpoints/", files)
  files <- files[!drop]

  if (length(files) == 0L)
    return(data.frame(path = character(0), bytes = numeric(0), md5 = character(0)))

  rel <- sub(paste0("^", here::here(), "/"), "", files)
  data.frame(path  = rel,
             bytes = file.size(files),
             md5   = unname(tools::md5sum(files)),
             stringsAsFactors = FALSE)[order(rel), ]
}

# Pull the headline numbers: every (metric, value) row in any *_summary.csv,
# plus a few structural counts. Returns a two-column key/value table.
fingerprint_metrics <- function() {
  out <- list()

  # (a) all metric/value summary tables, namespaced by their file name
  dd <- here::here("data_derived")
  summaries <- if (dir.exists(dd))
    list.files(dd, pattern = "_summary\\.csv$", recursive = TRUE, full.names = TRUE) else character(0)
  for (f in summaries) {
    tab <- try(utils::read.csv(f, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tab, "try-error")) next
    if (!all(c("metric", "value") %in% names(tab))) next   # different shape; skip
    tag <- tools::file_path_sans_ext(basename(f))
    out[[length(out) + 1L]] <- data.frame(
      key   = paste0(tag, "::", tab$metric),
      value = as.character(tab$value),
      stringsAsFactors = FALSE)
  }

  # (b) structural counts: row counts of the main derived tables
  main_tables <- c("emf_canada_em_only.csv", "emf_canada_combined.csv",
                   "ecm_native_canada_host_species.csv",
                   "clean_fungalroot_species.csv",
                   "clean_fungalroot_genera_table_s2.csv",
                   "genbank_emf_canada_long.csv", "globalfungi_canada_long.csv")
  for (nm in main_tables) {
    p <- file.path(dd, nm)
    if (!file.exists(p)) next
    tab <- try(utils::read.csv(p, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tab, "try-error")) next
    out[[length(out) + 1L]] <- data.frame(
      key = paste0("nrow::", nm), value = as.character(nrow(tab)),
      stringsAsFactors = FALSE)
  }

  # (c) host-list evidence_source split — the two-route (species OR genus) rule
  hp <- file.path(dd, "ecm_native_canada_host_species.csv")
  if (file.exists(hp)) {
    h <- try(utils::read.csv(hp, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(h, "try-error") && "evidence_source" %in% names(h)) {
      tb <- table(h$evidence_source)
      out[[length(out) + 1L]] <- data.frame(
        key = paste0("host_evidence_source::", names(tb)),
        value = as.character(as.integer(tb)), stringsAsFactors = FALSE)
    }
    if (!inherits(h, "try-error") && "growth_form" %in% names(h)) {
      tb <- table(h$growth_form)
      out[[length(out) + 1L]] <- data.frame(
        key = paste0("host_growth_form::", names(tb)),
        value = as.character(as.integer(tb)), stringsAsFactors = FALSE)
    }
  }

  # (d) per-source unique sampling locations (the Figure 1 / Table 1 numbers)
  wp <- file.path(dd, "wallacean", "wallacean_location_summary.csv")
  if (file.exists(wp)) {
    w <- try(utils::read.csv(wp, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(w, "try-error") && all(c("dataset", "n_unique_locations") %in% names(w)))
      out[[length(out) + 1L]] <- data.frame(
        key = paste0("unique_locations::", w$dataset),
        value = as.character(w$n_unique_locations), stringsAsFactors = FALSE)
  }

  res <- if (length(out)) do.call(rbind, out) else
    data.frame(key = character(0), value = character(0), stringsAsFactors = FALSE)
  res[order(res$key), ]
}

# Record the software environment, so an unexplained diff can be attributed.
fingerprint_env <- function() {
  pkgs <- c("dplyr", "readr", "tidyr", "sf", "terra", "ggplot2",
            "data.table", "rgbif", "BIEN", "rentrez")
  vers <- vapply(pkgs, function(p) {
    v <- try(as.character(utils::packageVersion(p)), silent = TRUE)
    if (inherits(v, "try-error")) NA_character_ else v
  }, character(1))
  data.frame(key   = c("R.version", "platform", paste0("pkg::", pkgs)),
             value = c(R.version.string, R.version$platform, unname(vers)),
             stringsAsFactors = FALSE)
}

write_fingerprint <- function(prefix) {
  dir.create(REPRO_DIR, showWarnings = FALSE, recursive = TRUE)
  f <- fingerprint_files(); m <- fingerprint_metrics(); e <- fingerprint_env()
  utils::write.csv(f, file.path(REPRO_DIR, paste0(prefix, "_files.csv")),   row.names = FALSE)
  utils::write.csv(m, file.path(REPRO_DIR, paste0(prefix, "_metrics.csv")), row.names = FALSE)
  utils::write.csv(e, file.path(REPRO_DIR, paste0(prefix, "_env.csv")),     row.names = FALSE)
  list(files = f, metrics = m, env = e)
}

# ---- Mode --------------------------------------------------------------------

.args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(.args)) .args[1] else
  if (exists("VERIFY_MODE")) get("VERIFY_MODE") else "baseline"
if (!MODE %in% c("baseline", "verify"))
  stop("Usage: Rscript scripts/99_verify_reproducibility.R [baseline|verify]", call. = FALSE)

# ---- baseline: capture the reference state -----------------------------------

if (MODE == "baseline") {
  fp <- write_fingerprint("baseline")
  cat(sprintf("Baseline captured: %d files, %d metrics -> %s\n",
              nrow(fp$files), nrow(fp$metrics), REPRO_DIR))
  cat("Next: back up data_derived/ and figures/, delete them, then run scripts/run_all.R\n")
}

# ---- verify: recapture and compare -------------------------------------------

if (MODE == "verify") {
  base_f <- file.path(REPRO_DIR, "baseline_files.csv")
  if (!file.exists(base_f))
    stop("No baseline found in ", REPRO_DIR,
         ". Run with 'baseline' BEFORE re-running the pipeline.", call. = FALSE)

  b_files   <- utils::read.csv(base_f, stringsAsFactors = FALSE)
  b_metrics <- utils::read.csv(file.path(REPRO_DIR, "baseline_metrics.csv"), stringsAsFactors = FALSE)
  b_env     <- utils::read.csv(file.path(REPRO_DIR, "baseline_env.csv"),     stringsAsFactors = FALSE)
  now       <- write_fingerprint("verify")

  rpt <- c(sprintf("Reproducibility verification — %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "")

  # --- 1. Metrics: the numbers the manuscript reports. Any change here matters.
  mm <- merge(b_metrics, now$metrics, by = "key", all = TRUE,
              suffixes = c("_baseline", "_now"))
  mm$state <- ifelse(is.na(mm$value_baseline), "NEW",
              ifelse(is.na(mm$value_now),      "MISSING",
              ifelse(mm$value_baseline == mm$value_now, "same", "CHANGED")))
  metric_problems <- mm[mm$state != "same", ]

  rpt <- c(rpt, "== METRICS (headline numbers) ==",
           sprintf("  tracked: %d | unchanged: %d | changed: %d | new: %d | missing: %d",
                   nrow(mm), sum(mm$state == "same"), sum(mm$state == "CHANGED"),
                   sum(mm$state == "NEW"), sum(mm$state == "MISSING")))
  if (nrow(metric_problems)) {
    rpt <- c(rpt, "", "  --- differences ---")
    for (i in seq_len(nrow(metric_problems)))
      rpt <- c(rpt, sprintf("  [%s] %s: baseline=%s  now=%s",
                            metric_problems$state[i], metric_problems$key[i],
                            metric_problems$value_baseline[i], metric_problems$value_now[i]))
  }

  # --- 2. Files. Split text (byte-comparable) from binary (timestamped).
  ff <- merge(b_files, now$files, by = "path", all = TRUE, suffixes = c("_baseline", "_now"))
  ext <- tolower(tools::file_ext(ff$path))
  ff$kind <- ifelse(ext %in% BINARY_EXT, "binary", "text")
  ff$state <- ifelse(is.na(ff$md5_baseline), "NEW",
              ifelse(is.na(ff$md5_now),      "MISSING",
              ifelse(ff$md5_baseline == ff$md5_now, "identical", "differs")))

  text_diff   <- ff[ff$kind == "text"   & ff$state == "differs", ]
  bin_diff    <- ff[ff$kind == "binary" & ff$state == "differs", ]
  missing_new <- ff[ff$state %in% c("NEW", "MISSING"), ]

  rpt <- c(rpt, "", "== FILES ==",
           sprintf("  tracked: %d | identical: %d | text differs: %d | binary differs: %d | new/missing: %d",
                   nrow(ff), sum(ff$state == "identical"), nrow(text_diff),
                   nrow(bin_diff), nrow(missing_new)))

  if (nrow(text_diff)) {
    rpt <- c(rpt, "", "  --- text outputs whose CONTENT changed (investigate) ---")
    for (i in seq_len(nrow(text_diff)))
      rpt <- c(rpt, sprintf("  %s  (%s -> %s bytes)", text_diff$path[i],
                            text_diff$bytes_baseline[i], text_diff$bytes_now[i]))
  }
  if (nrow(bin_diff)) {
    rpt <- c(rpt, "", "  --- binary/image outputs: byte differences are EXPECTED",
             "      (embedded timestamps). Size change is the signal worth checking. ---")
    for (i in seq_len(nrow(bin_diff))) {
      db <- bin_diff$bytes_now[i] - bin_diff$bytes_baseline[i]
      rpt <- c(rpt, sprintf("  %s  (%+d bytes%s)", bin_diff$path[i], db,
                            if (db == 0) " — size identical, almost certainly fine" else ""))
    }
  }
  if (nrow(missing_new)) {
    rpt <- c(rpt, "", "  --- new / missing outputs ---")
    for (i in seq_len(nrow(missing_new)))
      rpt <- c(rpt, sprintf("  [%s] %s", missing_new$state[i], missing_new$path[i]))
  }

  # --- 3. Environment drift is informational, not a failure.
  ee <- merge(b_env, now$env, by = "key", all = TRUE, suffixes = c("_baseline", "_now"))
  env_diff <- ee[is.na(ee$value_baseline) | is.na(ee$value_now) |
                 ee$value_baseline != ee$value_now, ]
  rpt <- c(rpt, "", "== ENVIRONMENT ==")
  if (nrow(env_diff)) {
    rpt <- c(rpt, "  (differences below may explain numeric drift)")
    for (i in seq_len(nrow(env_diff)))
      rpt <- c(rpt, sprintf("  %s: baseline=%s  now=%s", env_diff$key[i],
                            env_diff$value_baseline[i], env_diff$value_now[i]))
  } else {
    rpt <- c(rpt, "  identical to baseline")
  }

  # --- 4. Verdict. Metrics and text outputs are what must hold.
  fail <- nrow(metric_problems) > 0 || nrow(text_diff) > 0 || nrow(missing_new) > 0
  rpt <- c(rpt, "", "== VERDICT ==",
           if (fail) "  FAIL — see differences above"
           else      "  PASS — all headline metrics and text outputs reproduce",
           "",
           "  Note: this does not check the manuscript itself. Re-render it and",
           "  diff the text against the submitted .docx to confirm the reported",
           "  numbers are unchanged (see the README note).")

  writeLines(rpt, file.path(REPRO_DIR, "verify_report.txt"))
  cat(paste(rpt, collapse = "\n"), "\n")

  if (fail && !interactive()) quit(status = 1L)
}
