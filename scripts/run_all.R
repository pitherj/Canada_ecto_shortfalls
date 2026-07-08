# =============================================================================
# run_all.R  —  Run the whole ECM_manuscript pipeline in order
# =============================================================================
# WHAT THIS DOES
#   Sources scripts 01–21 in numeric order (which is the dependency order),
#   printing a banner and elapsed time for each. Every script is checkpoint-
#   guarded to varying degrees, so re-running skips work whose output already
#   exists; delete a specific output in data_derived/ or figures/ to force that
#   step to recompute.
#
# HOW TO RUN
#   From the project root (the folder containing ECM_manuscript.Rproj):
#       source(here::here("scripts", "run_all.R"))
#   or from a shell:
#       Rscript scripts/run_all.R
#
# OPTIONS (edit the two settings below, or set the environment variables)
#   SKIP_HEAVY  — if TRUE, skip the two on-demand steps that each scan the full
#                 ~13 GB GlobalFungi matrix (12_wallacean_density_map.R and
#                 13_wallacean_global_comparator.R). Handy for a routine rebuild
#                 that doesn't need the global comparators refreshed.
#                 Env override: ECM_SKIP_HEAVY=true / false
#   STOP_ON_ERROR — if TRUE (default), abort as soon as any script errors; if
#                 FALSE, log the failure and continue to the next script.
#                 Env override: ECM_STOP_ON_ERROR=true / false
# =============================================================================

# ---- Settings ---------------------------------------------------------------
SKIP_HEAVY    <- tolower(Sys.getenv("ECM_SKIP_HEAVY",    "false")) %in% c("true", "1", "yes")
STOP_ON_ERROR <- !(tolower(Sys.getenv("ECM_STOP_ON_ERROR", "true")) %in% c("false", "0", "no"))

# The two heavy, on-demand global-matrix steps.
HEAVY_SCRIPTS <- c("12_wallacean_density_map.R", "13_wallacean_global_comparator.R")

# ---- Discover the ordered script list ---------------------------------------
# here::here() anchors paths to the project root regardless of the working
# directory. We match only the NN_*.R analysis scripts (00_setup.R is sourced
# by each of them, and run_all.R excludes itself).
script_dir <- here::here("scripts")
scripts <- sort(list.files(script_dir, pattern = "^[0-9]{2}_.*\\.R$"))
# 00_setup.R is not an analysis step — it is sourced by each numbered script
# itself, so exclude it from the run list (it would otherwise appear as a
# no-op "[1/N]" step).
scripts <- setdiff(scripts, "00_setup.R")

if (SKIP_HEAVY) scripts <- setdiff(scripts, HEAVY_SCRIPTS)

banner <- function(txt) {
  cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
}

banner(sprintf("ECM_manuscript pipeline — %d scripts%s", length(scripts),
               if (SKIP_HEAVY) " (heavy global-matrix steps skipped)" else ""))

# ---- Run each script, timing and logging as we go ---------------------------
# Bookkeeping variables use a ".ra_" prefix so a sourced script can't clobber
# them, and each script is sourced into its OWN fresh environment (local = ...)
# so the objects a script creates (paths, emf, results, i, ...) stay isolated
# and never overwrite this runner's state. Scripts communicate through files in
# data_derived/, not through shared in-memory objects, so isolation is safe.
.ra_log <- data.frame(script = scripts, status = NA_character_,
                      seconds = NA_real_, stringsAsFactors = FALSE)

for (.ra_i in seq_along(scripts)) {
  .ra_s <- scripts[.ra_i]
  banner(sprintf("[%d/%d] %s", .ra_i, length(scripts), .ra_s))
  .ra_t0 <- Sys.time()
  .ra_ok <- tryCatch({
    source(file.path(script_dir, .ra_s), local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    message(sprintf("  !! ERROR in %s: %s", .ra_s, conditionMessage(e)))
    FALSE
  })
  .ra_log$seconds[.ra_i] <- round(as.numeric(difftime(Sys.time(), .ra_t0, units = "secs")), 1)
  .ra_log$status[.ra_i]  <- if (.ra_ok) "ok" else "ERROR"
  cat(sprintf("  -> %s in %.1f s\n", .ra_log$status[.ra_i], .ra_log$seconds[.ra_i]))
  if (!.ra_ok && STOP_ON_ERROR) {
    banner(sprintf("Stopped: %s failed. Fix it and re-run (completed steps are skipped).", .ra_s))
    stop(sprintf("Pipeline halted at %s", .ra_s), call. = FALSE)
  }
}

# ---- Summary ----------------------------------------------------------------
banner("Pipeline summary")
print(.ra_log, row.names = FALSE)
cat(sprintf("\nTotal: %.1f s across %d scripts (%d ok, %d error).\n",
            sum(.ra_log$seconds, na.rm = TRUE), nrow(.ra_log),
            sum(.ra_log$status == "ok"), sum(.ra_log$status == "ERROR")))
cat("\nAll data_derived/ outputs and figures/ regenerated.\n")
