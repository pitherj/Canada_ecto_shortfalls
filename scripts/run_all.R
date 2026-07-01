# =============================================================================
# run_all.R  —  Run the whole ECM_manuscript pipeline in order
# =============================================================================
# WHAT THIS DOES
#   Sources scripts 01–20 in numeric order (which is the dependency order),
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

if (SKIP_HEAVY) scripts <- setdiff(scripts, HEAVY_SCRIPTS)

banner <- function(txt) {
  cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
}

banner(sprintf("ECM_manuscript pipeline — %d scripts%s", length(scripts),
               if (SKIP_HEAVY) " (heavy global-matrix steps skipped)" else ""))

# ---- Run each script, timing and logging as we go ---------------------------
results <- data.frame(script = scripts, status = NA_character_,
                      seconds = NA_real_, stringsAsFactors = FALSE)

for (i in seq_along(scripts)) {
  s <- scripts[i]
  banner(sprintf("[%d/%d] %s", i, length(scripts), s))
  t0 <- Sys.time()
  ok <- tryCatch({
    source(file.path(script_dir, s))
    TRUE
  }, error = function(e) {
    message(sprintf("  !! ERROR in %s: %s", s, conditionMessage(e)))
    FALSE
  })
  results$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results$status[i]  <- if (ok) "ok" else "ERROR"
  cat(sprintf("  -> %s in %.1f s\n", results$status[i], results$seconds[i]))
  if (!ok && STOP_ON_ERROR) {
    banner(sprintf("Stopped: %s failed. Fix it and re-run (completed steps are skipped).", s))
    stop(sprintf("Pipeline halted at %s", s), call. = FALSE)
  }
}

# ---- Summary ----------------------------------------------------------------
banner("Pipeline summary")
print(results, row.names = FALSE)
cat(sprintf("\nTotal: %.1f s across %d scripts (%d ok, %d error).\n",
            sum(results$seconds, na.rm = TRUE), nrow(results),
            sum(results$status == "ok"), sum(results$status == "ERROR")))
cat("\nNext: render the Supplemental Materials with\n",
    '  quarto::quarto_render(here::here("supplemental_materials.qmd"))\n', sep = "")
