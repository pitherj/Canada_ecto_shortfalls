# =============================================================================
# Download BIEN2 Range Shapefiles from biendata.org
# =============================================================================
# Downloads modelled range shapefiles for all native Canadian EcM host species
# from biendata.org (Moulatlet et al. 2025, PNAS), then unzips each archive
# into a per-species subdirectory under data_raw/bien2_ranges/.
#
# This script is part of Group 0 (data acquisition) and must be run before
# 08_host_rasters.R.  It depends on the host species list produced by
# 06_host_species.R.
#
# Citation:
#   Moulatlet GM et al. (2025) General laws of biodiversity: Climatic niches
#   predict plant range size and ecological dominance globally. PNAS.
#   Data portal: https://www.biendata.org
#
# Prerequisites:
#   data_derived/ecm_native_canada_host_species.csv  (06_host_species.R)
#
# Outputs (per species):
#   data_raw/bien2_ranges/<Genus_species>.zip          — downloaded archive
#   data_raw/bien2_ranges/<Genus_species>/<files>      — unzipped shapefile
# Summary log:
#   data_raw/bien2_ranges/download_log.csv             — per-species status log
#
# Re-run behaviour:
#   - Species with an existing unzipped .shp are pre-filtered out before the
#     loop (no wasted iterations).
#   - Species logged as not_available / unzip_error / network_error are skipped.
#   - Species logged as rate_limited ARE retried on the next run.
#   - Species with an existing .zip but no .shp are unzipped without
#     re-downloading.
#   - HTTP 429 triggers a wait of RATE_LIMIT_WAIT_SEC then retries (up to
#     MAX_RETRIES times); if all retries fail the species is logged as
#     rate_limited and the run continues.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(httr)

# ---- Paths -------------------------------------------------------------------

out_dir  <- paths$bien2_ranges_dir
log_path <- file.path(out_dir, "download_log.csv")

# ---- Parameters --------------------------------------------------------------

BASE_URL            <- "https://biendata.org/api/download/range"
SLEEP_SEC           <- 1.5       # polite pause between requests
TIMEOUT_SEC         <- 60        # per-request timeout
RATE_LIMIT_WAIT_SEC <- 16 * 60   # wait on HTTP 429 (site says 15 min)
MAX_RETRIES         <- 3         # retry attempts per species on 429

# ---- Check prerequisite ------------------------------------------------------

if (!file.exists(paths$host_species)) {
  stop("Host species list not found:\n  ", paths$host_species,
       "\nRun 06_host_species.R first.")
}

hosts        <- readr::read_csv(paths$host_species, show_col_types = FALSE)
species_list <- sort(hosts$species)
n_total      <- length(species_list)

ts(sprintf("[bien2] %d species in host list", n_total))

# ---- Create output directory -------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: check for .shp in a directory -----------------------------------

shp_exists <- function(sp_dir) {
  length(list.files(sp_dir, pattern = "\\.shp$", full.names = FALSE)) > 0
}

# ---- Pre-filter 1: species with shapefiles already on disk ------------------
# These are skipped before the loop; no iteration needed.

already_shp <- vapply(species_list, function(sp) {
  sp_dir <- file.path(out_dir, gsub(" ", "_", sp))
  dir.exists(sp_dir) && shp_exists(sp_dir)
}, logical(1))

skip_results <- lapply(species_list[already_shp], function(sp) {
  sp_dir <- file.path(out_dir, gsub(" ", "_", sp))
  shp    <- list.files(sp_dir, pattern = "\\.shp$", full.names = TRUE)[[1]]
  data.frame(species = sp, status = "skipped",
             shp = shp, note = NA_character_,
             stringsAsFactors = FALSE)
})

to_process <- species_list[!already_shp]

# ---- Pre-filter 2: species already settled in the log -----------------------
# Statuses considered settled (will NOT be retried): not_available, unzip_error,
# network_error.  rate_limited IS retried.

settled_results <- list()

if (file.exists(log_path)) {
  log_prev      <- readr::read_csv(log_path, show_col_types = FALSE)
  settled_stats <- c("not_available", "unzip_error", "network_error")
  settled_rows  <- log_prev[log_prev$species %in% to_process &
                              log_prev$status %in% settled_stats, ]
  if (nrow(settled_rows) > 0) {
    settled_results <- lapply(seq_len(nrow(settled_rows)), function(j) {
      data.frame(species = settled_rows$species[[j]],
                 status  = settled_rows$status[[j]],
                 shp     = NA_character_,
                 note    = settled_rows$note[[j]],
                 stringsAsFactors = FALSE)
    })
    to_process <- to_process[!to_process %in% settled_rows$species]
  }
}

n_skip_shp <- sum(already_shp)
n_skip_log <- length(settled_results)
n_todo     <- length(to_process)

ts(sprintf("[bien2] %d species skipped (shapefile on disk).", n_skip_shp))
if (n_skip_log > 0)
  ts(sprintf("[bien2] %d species skipped (settled in log).", n_skip_log))
ts(sprintf("[bien2] %d species to download.", n_todo))

# ---- Download + unzip loop ---------------------------------------------------

if (n_todo > 0) {

  dl_results <- lapply(seq_along(to_process), function(i) {

    sp     <- to_process[[i]]
    sp_us  <- gsub(" ", "_", sp)
    sp_zip <- file.path(out_dir, paste0(sp_us, ".zip"))
    sp_dir <- file.path(out_dir, sp_us)

    # ---- Download if zip not already on disk ---------------------------------
    if (!file.exists(sp_zip) || file.size(sp_zip) <= 500) {

      Sys.sleep(SLEEP_SEC)

      # Retry loop — handles HTTP 429
      dl_result <- NULL
      for (attempt in seq_len(MAX_RETRIES)) {

        resp <- tryCatch(
          httr::GET(paste0(BASE_URL, "?species=", sp_us),
                    httr::timeout(TIMEOUT_SEC)),
          error = function(e) {
            message(sprintf("  [%d/%d] network error: %s — %s",
                            i, n_todo, sp, conditionMessage(e)))
            NULL
          }
        )

        if (is.null(resp)) {
          dl_result <- data.frame(species = sp, status = "network_error",
                                  shp = NA_character_, note = "request failed",
                                  stringsAsFactors = FALSE)
          break
        }

        code <- httr::status_code(resp)
        ct   <- httr::headers(resp)[["content-type"]]
        ct   <- if (is.null(ct)) "" else ct

        # ---- Rate-limited: wait and retry ------------------------------------
        if (code == 429) {
          if (attempt < MAX_RETRIES) {
            ts(sprintf("  [%d/%d] HTTP 429 (attempt %d/%d) — waiting %.0f min: %s",
                       i, n_todo, attempt, MAX_RETRIES,
                       RATE_LIMIT_WAIT_SEC / 60, sp))
            Sys.sleep(RATE_LIMIT_WAIT_SEC)
            next
          } else {
            ts(sprintf("  [%d/%d] HTTP 429 — exhausted %d retries, logging as rate_limited: %s",
                       i, n_todo, MAX_RETRIES, sp))
            dl_result <- data.frame(species = sp, status = "rate_limited",
                                    shp = NA_character_,
                                    note = paste0("HTTP 429 after ", MAX_RETRIES, " attempts"),
                                    stringsAsFactors = FALSE)
            break
          }
        }

        # ---- Server returned JSON → no range available ----------------------
        if (!(code == 200 &&
              grepl("zip|octet-stream|shapefile|binary", ct, ignore.case = TRUE))) {
          server_msg <- tryCatch(
            httr::content(resp, "parsed", encoding = "UTF-8")[["message"]],
            error = function(e) NULL
          )
          note <- if (!is.null(server_msg)) server_msg else
            paste0("HTTP ", code, " ct=", ct)
          message(sprintf("  [%d/%d] not available: %s  (%s)", i, n_todo, sp, note))
          dl_result <- data.frame(species = sp, status = "not_available",
                                  shp = NA_character_, note = note,
                                  stringsAsFactors = FALSE)
          break
        }

        # ---- Success ---------------------------------------------------------
        writeBin(httr::content(resp, "raw"), sp_zip)
        ts(sprintf("  [%d/%d] downloaded: %s  (%.0f KB)",
                   i, n_todo, sp, file.size(sp_zip) / 1024))
        dl_result <- "downloaded"
        break

      }  # end retry loop

      if (!identical(dl_result, "downloaded")) return(dl_result)

    } else {
      message(sprintf("  [%d/%d] zip exists, unzipping: %s", i, n_todo, sp))
    }

    # ---- Unzip into per-species subdirectory ---------------------------------
    dir.create(sp_dir, showWarnings = FALSE)
    unzip_ok <- tryCatch({
      utils::unzip(sp_zip, exdir = sp_dir)
      TRUE
    }, error = function(e) {
      message(sprintf("  [%d/%d] unzip error: %s — %s",
                      i, n_todo, sp, conditionMessage(e)))
      FALSE
    })

    if (!unzip_ok || !shp_exists(sp_dir)) {
      return(data.frame(species = sp, status = "unzip_error",
                        shp = NA_character_, note = "unzip failed or no .shp found",
                        stringsAsFactors = FALSE))
    }

    shp <- list.files(sp_dir, pattern = "\\.shp$", full.names = TRUE)[[1]]
    return(data.frame(species = sp, status = "success",
                      shp = shp, note = NA_character_,
                      stringsAsFactors = FALSE))
  })

} else {
  dl_results <- list()
}

# ---- Assemble and write log --------------------------------------------------

download_log <- dplyr::bind_rows(skip_results, settled_results, dl_results)
readr::write_csv(download_log, log_path)

n_ok   <- sum(download_log$status %in% c("success", "skipped"))
n_new  <- sum(download_log$status == "success")
n_skip <- sum(download_log$status == "skipped")
n_miss <- sum(download_log$status == "not_available")
n_rl   <- sum(download_log$status == "rate_limited")
n_err  <- sum(download_log$status %in% c("network_error", "unzip_error"))
n_tot  <- nrow(download_log)

ts("")
ts(sprintf("[bien2] ── Complete ─────────────────────────────────────────"))
ts(sprintf("  Available    : %d / %d  (%.0f%%)",
           n_ok, n_tot, 100 * n_ok / n_tot))
ts(sprintf("  Newly fetched: %d  |  skipped: %d", n_new, n_skip))
ts(sprintf("  Not available: %d  |  errors: %d", n_miss, n_err))
if (n_rl > 0)
  ts(sprintf("  Rate-limited : %d  (re-run after cooldown to retry)", n_rl))
ts(sprintf("  Log -> %s", basename(log_path)))

if (n_miss > 0) {
  message("")
  message("[bien2] Species not found on biendata.org:")
  message(paste0("  ",
                 download_log$species[download_log$status == "not_available"],
                 collapse = "\n"))
}
if (n_rl > 0) {
  message("")
  message("[bien2] Species to retry (rate-limited — re-run this script):")
  message(paste0("  ",
                 download_log$species[download_log$status == "rate_limited"],
                 collapse = "\n"))
}

ts("07_bien2_ranges.R complete.")
