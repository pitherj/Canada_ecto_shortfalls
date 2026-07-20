# =============================================================================
# run_all.R  —  Run the whole ECM_manuscript pipeline in order
# =============================================================================
# WHAT THIS DOES
#   Sources scripts 01–21 in numeric order (which is the dependency order),
#   timing each one and writing the run record to data_derived/run_log.csv
#   (script, status, seconds, started). Every script is checkpoint-
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
# no-op "[1/N]" step). 99_verify_reproducibility.R is a stand-alone verification
# tool, run manually before and after a pipeline run, not a pipeline step.
scripts <- setdiff(scripts, c("00_setup.R", "99_verify_reproducibility.R"))

if (SKIP_HEAVY) scripts <- setdiff(scripts, HEAVY_SCRIPTS)

# ---- Run each script, timing and logging as we go ---------------------------
# Bookkeeping variables use a ".ra_" prefix so a sourced script can't clobber
# them, and each script is sourced into its OWN fresh environment (local = ...)
# so the objects a script creates (paths, emf, results, i, ...) stay isolated
# and never overwrite this runner's state. Scripts communicate through files in
# data_derived/, not through shared in-memory objects, so isolation is safe.
#
# The log is rewritten after every script, so a halted run still leaves a
# complete record of everything that ran up to the point of failure. Errors are
# raised as warnings (and, under STOP_ON_ERROR, a stop()) so that an unattended
# run still surfaces them on stderr.
.ra_logfile <- here::here("data_derived", "run_log.csv")
dir.create(dirname(.ra_logfile), showWarnings = FALSE, recursive = TRUE)

.ra_log <- data.frame(script = scripts, status = NA_character_,
                      seconds = NA_real_, started = NA_character_,
                      stringsAsFactors = FALSE)

for (.ra_i in seq_along(scripts)) {
  .ra_s  <- scripts[.ra_i]
  .ra_t0 <- Sys.time()
  .ra_log$started[.ra_i] <- format(.ra_t0, "%Y-%m-%d %H:%M:%S")
  .ra_ok <- tryCatch({
    source(file.path(script_dir, .ra_s), local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    warning(sprintf("ERROR in %s: %s", .ra_s, conditionMessage(e)), call. = FALSE)
    FALSE
  })
  .ra_log$seconds[.ra_i] <- round(as.numeric(difftime(Sys.time(), .ra_t0, units = "secs")), 1)
  .ra_log$status[.ra_i]  <- if (.ra_ok) "ok" else "ERROR"
  utils::write.csv(.ra_log, .ra_logfile, row.names = FALSE)
  if (!.ra_ok && STOP_ON_ERROR)
    stop(sprintf("Pipeline halted at %s - see %s", .ra_s, .ra_logfile), call. = FALSE)
}

utils::write.csv(.ra_log, .ra_logfile, row.names = FALSE)
