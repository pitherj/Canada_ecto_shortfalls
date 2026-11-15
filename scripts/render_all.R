# =============================================================================
# render_all.R  —  Render the manuscript and both Supplemental Materials
# =============================================================================
# PURPOSE
#   Render manuscript/manuscript.qmd (PDF and Word, as set in its YAML),
#   manuscript/supplemental_SM1.qmd and manuscript/supplemental_SM2.qmd (HTML)
#   in one step, then check every rendered document for a missing value.
#   The documents read their numbers from data_derived/ as they render, so a
#   number that cannot be found appears as "NA". Tables 1-4 already stop the
#   render if they contain NA (assert_no_na() in manuscript.qmd); this script
#   extends the check to the running text and to SM1 and SM2. It stops, and
#   lists each offending line, if any rendered document contains a standalone
#   "NA" outside its reference list (where author initials such as "Hynson NA"
#   are legitimate).
#
# INPUTS
#   manuscript/manuscript.qmd, supplemental_SM1.qmd, supplemental_SM2.qmd
#   data_derived/ and figures/ (the full pipeline must have been run)
#   Quarto on PATH (https://quarto.org); a LaTeX installation for the PDF
#   (quarto install tinytex)
#
# OUTPUTS
#   manuscript/_output/<document>/  rendered documents (not tracked in git)
#
# MANUSCRIPT ITEMS IT FEEDS
#   All of them: this produces the documents.
#
# RE-RUN WHEN
#   Any pipeline output or any of the three QMD files changes.
#   Run from the project root:  Rscript scripts/render_all.R
# =============================================================================

source(here::here("scripts", "00_setup.R"))

ms_dir  <- here::here("manuscript")
out_dir <- file.path(ms_dir, "_output")
dir.create(out_dir, showWarnings = FALSE)

docs <- c("manuscript.qmd", "supplemental_SM1.qmd", "supplemental_SM2.qmd")

if (!nzchar(Sys.which("quarto"))) stop("Quarto not found on PATH (https://quarto.org).")

# ---- 1. Render each document -------------------------------------------------
# Quarto is run from inside manuscript/, so each document finds its
# references.bib and facets.csl beside it. Each document gets its own folder in
# _output/, because Quarto empties --output-dir at every render. A failed
# render stops the script.
for (d in docs) {
  message("Rendering ", d, " ...")
  doc_dir <- file.path(out_dir, tools::file_path_sans_ext(d))
  status <- withr::with_dir(ms_dir,
    system2("quarto", c("render", shQuote(d), "--output-dir", shQuote(doc_dir))))
  if (status != 0L) stop("Quarto failed to render ", d, " (exit status ", status, ").")
}

# ---- 2. Check the rendered text for missing values ---------------------------
# Each Word and HTML output is converted to Markdown with Quarto's bundled
# pandoc, which keeps code formatting as `backticks`. Text in code formatting
# is removed before checking, because the documents deliberately write NA as
# a value there (e.g. SM2: "sh_code = NA"). Lines after the reference-list
# heading are not checked.
rendered <- list.files(out_dir, pattern = "\\.(docx|html)$", full.names = TRUE,
                       recursive = TRUE)
na_pattern <- "(^|[^[:alpha:]])NA([^[:alpha:]]|$)"

problems <- purrr::map_dfr(rendered, function(f) {
  txt <- system2("quarto", c("pandoc", shQuote(f), "-t", "gfm", "--wrap=none"),
                 stdout = TRUE)
  txt <- gsub("`[^`]*`", "", txt)                 # drop code-formatted text
  txt <- gsub("<[^>]+>", " ", txt)                # drop any leftover HTML tags
  # drop embedded image data (long runs of base64 text with no spaces)
  txt <- txt[!grepl("^\\S{100,}$", trimws(txt))]
  ref_start <- grep("^#*\\s*(References|REFERENCES|Literature cited)\\s*(\\{.*\\})?\\s*$", txt)
  if (length(ref_start)) txt <- txt[seq_len(ref_start[1] - 1L)]
  hits <- grep(na_pattern, txt)
  tibble::tibble(file = rep(basename(f), length(hits)),
                 line = hits,
                 text = substr(txt[hits], 1, 160))
})

if (nrow(problems) > 0L) {
  print(as.data.frame(problems), right = FALSE)
  stop(nrow(problems), " line(s) in the rendered documents contain \"NA\". ",
       "A value could not be found in data_derived/; see the lines above.",
       call. = FALSE)
}

message("Rendered ", length(docs), " documents to ", out_dir,
        "; no missing values found in ", length(rendered), " Word/HTML outputs.")
