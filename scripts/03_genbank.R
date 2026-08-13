# =============================================================================
# GenBank Sequence Retrieval and SH Assignment
# =============================================================================
# Retrieves fungal ITS sequences from GenBank for Canadian records, extracts
# the ITS2 sub-region with ITSx (the GenBank-side analogue of GlobalFungi's
# ITS2/ITSboth restriction), and assigns each ITS2 fragment to UNITE Species
# Hypotheses via vsearch at 98.5% identity (matches GlobalFungi's documented
# BLASTn SH-assignment criterion). Ties are resolved with dark-taxa awareness
# and a genus-level fallback (Step 6). Each step is checkpointed — safe to
# re-run; skips completed steps. See the GenBank methods in
# FACETS/supplemental_materials_SM1_FACETS.qmd for the full rationale and
# citations.
#
# Download completeness (Steps 1-3b)
#   The retrieval steps compare records returned against records requested and
#   stop if they come up short; failed batches are retried before being declared
#   lost; and a checkpoint is only put in place once its completeness check has
#   passed, so a failed run cannot leave a truncated file for a later run to
#   adopt. If a checkpoint exists but is INCOMPLETE, Steps 2/3/3b fetch only the
#   records it is missing and append them, rather than skipping the step (which
#   would make the gap permanent) or re-downloading everything (which would mix
#   retrieval dates). See the GB_FETCH_* / GB_MAX_MISSING_FRAC constants below.
#   This behaviour was added in 2026-08 after a silent loss of 3,001 of 60,911
#   records; repro/genbank_fetch_gap_prefix_state.md documents that incident.
#
# Prerequisites:
#   02_globalfungi.R must have been run (provides UNITE lookup)
#   ITSx (>= 1.1) + HMMER (3.x) installed on PATH:
#     HMMER:  brew install hmmer  (macOS)  /  sudo apt install hmmer  (Ubuntu)
#     ITSx:   https://microbiology.se/software/itsx/ — put the 'ITSx' script
#             (and its bundled HMM-profile directory) on PATH
#   vsearch installed on PATH:
#     macOS:  brew install vsearch
#     Ubuntu: sudo apt install vsearch
#     Other:  https://github.com/torognes/vsearch/releases
#   NCBI API key in .Renviron:  ENTREZ_KEY=<your key>
#     Free key: https://www.ncbi.nlm.nih.gov/account/
#
# Search strategy:
#   "Canada[Country] AND Fungi[Organism]
#    AND (internal transcribed spacer OR ITS)"
#
# Outputs (all in data_derived/ or data_derived/temp/):
#   temp/genbank_emf_canada_ids.txt          — GenBank UIDs
#   temp/genbank_emf_canada.fasta            — retrieved FASTA sequences
#   temp/genbank_emf_canada_metadata.csv     — esummary metadata
#   temp/genbank_no_sequence_available.csv   — UIDs individually confirmed to
#                                              have no retrievable sequence
#                                              (WGS/TLS "master" records), and
#                                              therefore excluded from the
#                                              Step 3b completeness requirement
#   temp/genbank_itsx.ITS2.fasta             — ITSx-extracted ITS2 regions
#   temp/genbank_itsx.ITS1.fasta             — ITSx-extracted ITS1 (its_region flag only)
#   temp/genbank_vsearch_query.fasta         — ITS2 fragments (>=100 bp, sanitized)
#   temp/genbank_vsearch_hits.txt            — vsearch blast6 output (tied hits)
#   temp/genbank_ambiguous_species_excluded.csv — records excluded (genus-ambiguous ties)
#   temp/genbank_its1_resolution_crosstab.csv    — ITS1-discard diagnostic (reporting only)
#   genbank_emf_canada_long.csv              — final annotated table; carries
#                                              its_region and taxonomic_resolution
#                                              columns, and may contain rows with
#                                              sh_code = NA (genus-resolved)
#
# Host name handling:
#   The final table carries two host columns:
#     host_taxon_raw  — the value extracted from the structured `host`
#                       qualifier (or, as fallback, regex-parsed from the
#                       `isolation_source` field). Preserved for traceability.
#     host_taxon      — canonicalized form used by all downstream analyses.
#                       See canonicalize_host() defined below for the rules.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(rentrez)

if (Sys.getenv("ENTREZ_KEY") == "") {
  stop("NCBI API key not found. Set ENTREZ_KEY in your .Renviron file.\n",
       "  Free key: https://www.ncbi.nlm.nih.gov/account/")
}

# UNITE FASTA: pinned reference, not auto-detected (see 00_setup.R). This MUST match
# paths$unite_fasta used by 02_globalfungi.R —
# if the two roles pointed to different UNITE builds, combined "GF + GB"
# SH counts would be inflated by spurious cross-build mismatches, since SH
# numbering is not stable across builds.
unite_fasta_path <- paths$unite_fasta
if (!file.exists(unite_fasta_path)) {
  stop(
    "Pinned UNITE FASTA not found: ", unite_fasta_path, "\n",
    "Either restore this file to data_raw/UNITE/, or — if intentionally\n",
    "switching builds — update paths$unite_fasta (and\n",
    "paths$unite_fasta, which must match it) in 00_setup.R."
  )
}

# Minimum length (bp) an ITSx-extracted ITS2 fragment must reach to be used
# for SH assignment. ITSx occasionally calls a spuriously short "ITS2 region"
# (as short as a few bp) on low-quality input; matching such a fragment
# against UNITE at 98.5% would trivially tie with thousands of unrelated SH
# codes, so fragments below this floor are excluded before assignment.
ITS2_MIN_LENGTH <- 100L

# ITSx's HMMER backend errors on any single input sequence longer than
# ~100,000 bp. A handful of whole-genome/scaffold GenBank entries swept in by
# the loose text search exceed this, and are excluded before ITSx.
ITSX_MAX_INPUT_LENGTH <- 50000L

# Maximum tolerated fraction of vsearch hits left unmatched after the UNITE
# taxonomy join (Step 5 below). Unlike the GlobalFungi join in 02_globalfungi.R
# (where SH
# codes are pre-assigned externally and a residual mismatch is expected and
# documented), GenBank's SH codes are assigned in-house via vsearch against
# this same file, so the join is self-consistent by construction and should
# be ~100% matched. Any non-trivial unmatched fraction here indicates a
# broken invariant (e.g. unite_fasta_genbank pointing to a different file
# than the one vsearch actually ran against) rather than an expected gap.
SH_MAX_UNMATCHED_FRAC <- 0.01

# ---- Batch-download completeness guards (Steps 1, 2, 3) ---------------------
# Steps 1-3 retrieve records from NCBI in batches. Until 2026-08 a batch that
# hit a transient network or server error was skipped with only a warning: the
# affected records were left out of the checkpoint, nothing compared records
# returned against records requested, and the checkpoint was written as though
# complete. Because those steps are guarded by if (!file.exists(...)), every
# later run then reused the truncated checkpoint and the loss became permanent.
# That is how 3,001 of 60,911 records went missing from the original pull. The
# three constants below, together with the guards and the repair path in Steps
# 2, 3 and 3b, close that hole.
#
# GB_FETCH_MAX_ATTEMPTS — how many times a single batch is tried before it is
#   declared failed. NCBI's eutils endpoints time out sporadically under load
#   and the same request usually succeeds on a repeat, so retrying (rather than
#   skipping) is the correct first response to a failed batch.
# GB_FETCH_RETRY_PAUSE — seconds to wait before retrying. The wait is
#   multiplied by the attempt number (3s, 6s, 9s, 12s), so the script backs off
#   politely if NCBI is rate-limiting rather than merely glitching.
# GB_MAX_MISSING_FRAC — the fraction of REQUESTED records a step may still be
#   short of, after retries, before it is treated as broken. The tolerance is
#   not zero because a handful of UIDs can legitimately vanish between the
#   search and the fetch (NCBI occasionally withdraws or suppresses a record).
#   It is deliberately set far below the size of one batch: a single lost
#   200-record batch is 0.33% of a 60,911-record pull, comfortably above 0.1%,
#   so losing even one whole batch trips the guard. A batch that exhausts its
#   retries is fatal regardless of this fraction (see gb_assert_fetch_complete).
GB_FETCH_MAX_ATTEMPTS <- 5L
GB_FETCH_RETRY_PAUSE  <- 3
GB_MAX_MISSING_FRAC   <- 0.001

# Additional path variables
gb_ids_path      <- file.path(paths$temp_dir, "genbank_emf_canada_ids.txt")
gb_fasta_path    <- file.path(paths$temp_dir, "genbank_emf_canada.fasta")
gb_meta_path     <- file.path(paths$temp_dir, "genbank_emf_canada_metadata.csv")
gb_itsx_prefix   <- file.path(paths$temp_dir, "genbank_itsx")
gb_itsx_its1     <- paste0(gb_itsx_prefix, ".ITS1.fasta")
gb_itsx_its2     <- paste0(gb_itsx_prefix, ".ITS2.fasta")
gb_vsearch_query <- file.path(paths$temp_dir, "genbank_vsearch_query.fasta")
gb_vsearch_path  <- file.path(paths$temp_dir, "genbank_vsearch_hits.txt")
gb_ambiguous_path <- file.path(paths$temp_dir, "genbank_ambiguous_species_excluded.csv")
gb_its1_diag_path <- file.path(paths$temp_dir, "genbank_its1_resolution_crosstab.csv")
gb_fetch_log     <- file.path(paths$temp_dir, "genbank_fetch_log.txt")
gb_nosequence_path <- file.path(paths$temp_dir, "genbank_no_sequence_available.csv")

# ---- Small local helper: read a (possibly multi-line) FASTA into a tibble --
# ITSx and vsearch both consume/produce plain FASTA files. Base R has no
# built-in FASTA reader, and this project does not otherwise depend on a
# sequence-I/O package, so a small dependency-free helper is defined here.
# Returns one row per sequence: `header` (without the leading ">") and `seq`
# (concatenated across any wrapped lines).
read_fasta_tbl <- function(path) {
  lines <- readr::read_lines(path)
  is_header <- startsWith(lines, ">")
  header_idx <- which(is_header)
  if (length(header_idx) == 0L) {
    return(tibble::tibble(header = character(0), seq = character(0)))
  }
  seq_end <- c(header_idx[-1] - 1L, length(lines))
  tibble::tibble(
    header = sub("^>", "", lines[header_idx]),
    seq = purrr::map2_chr(header_idx + 1L, seq_end, function(s, e) {
      if (s > e) return("")
      paste0(lines[s:e], collapse = "")
    })
  )
}

# ---- Small local helper: list the accessions already present in a FASTA -----
# Sequences are requested from NCBI by UID but arrive keyed by ACCESSION, so
# "which records do I already have a sequence for?" is answered in accession
# space. entrez_fetch writes headers of the form ">ACCESSION.VERSION descr...",
# so the accession is everything up to the first "." or space. Returns an empty
# vector if the file does not exist yet.
gb_fasta_accessions <- function(path) {
  if (!file.exists(path)) return(character(0L))
  headers <- grep("^>", readr::read_lines(path), value = TRUE)
  sub("^>([^ .]+).*$", "\\1", headers)
}

# ---- Small local helper: fetch one batch, retrying transient failures -------
# Runs fetch_fun(ids) up to GB_FETCH_MAX_ATTEMPTS times, pausing between tries.
# Returns whatever fetch_fun returned on the first success, or NULL if every
# attempt failed. NULL means "this batch is lost" — callers count those and hand
# the count to gb_assert_fetch_complete(), which stops the script. A batch is
# never silently dropped.
gb_fetch_batch <- function(fetch_fun, ids, what, batch_i, n_batches) {
  for (attempt in seq_len(GB_FETCH_MAX_ATTEMPTS)) {
    result <- tryCatch(fetch_fun(ids), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (attempt < GB_FETCH_MAX_ATTEMPTS) {
      message(sprintf(
        "  %s batch %d/%d failed (attempt %d of %d): %s\n    retrying in %g s",
        what, batch_i, n_batches, attempt, GB_FETCH_MAX_ATTEMPTS,
        conditionMessage(result), GB_FETCH_RETRY_PAUSE * attempt))
      Sys.sleep(GB_FETCH_RETRY_PAUSE * attempt)
    } else {
      warning(sprintf("%s batch %d/%d failed on all %d attempts: %s",
                      what, batch_i, n_batches, GB_FETCH_MAX_ATTEMPTS,
                      conditionMessage(result)), call. = FALSE)
    }
  }
  NULL
}

# ---- Small local helper: the completeness guard ------------------------------
# Compares records actually returned against records requested and stops with an
# informative error if the step came up short. Modelled on the
# SH_MAX_UNMATCHED_FRAC coverage gate in Step 5. Two independent trip
# conditions:
#   (a) any batch exhausted its retries — a known, located loss;
#   (b) the shortfall exceeds GB_MAX_MISSING_FRAC — catches partial returns that
#       raised no error at all.
# `partial_path`, when supplied, is named in the error message so the operator
# knows where the half-written file was left.
gb_assert_fetch_complete <- function(n_returned, n_requested, n_failed_batches,
                                     what, partial_path = NULL) {
  shortfall <- n_requested - n_returned
  frac      <- if (n_requested == 0L) 0 else shortfall / n_requested
  if (n_failed_batches == 0L && frac <= GB_MAX_MISSING_FRAC) {
    return(invisible(TRUE))
  }
  stop(sprintf(
    paste0(
      "GenBank %s fetch is incomplete: requested %s records, received %s ",
      "(shortfall %s, %.3f%%; tolerance %.3f%%). Batches that failed every ",
      "one of their %d attempts: %d.\n",
      "  This is the failure that must NOT be allowed to pass silently: a ",
      "truncated checkpoint written here would be reused by every later run.\n",
      "  Nothing has been written to the checkpoint. Re-run the script — the ",
      "repair path fetches only the records that are still missing.%s"
    ),
    what, format(n_requested, big.mark = ","), format(n_returned, big.mark = ","),
    format(shortfall, big.mark = ","), 100 * frac, 100 * GB_MAX_MISSING_FRAC,
    GB_FETCH_MAX_ATTEMPTS, n_failed_batches,
    if (is.null(partial_path)) "" else
      paste0("\n  Partial output left at: ", partial_path)
  ), call. = FALSE)
}

# ---- Small local helper: fetch FASTA sequences for a set of UIDs -------------
# Returns the raw FASTA text in per-batch chunks, the number of sequence records
# those chunks contain, and how many batches were lost outright. The caller
# writes the chunks and applies the completeness guard.
gb_fetch_sequences <- function(uids, batch_size = 200L) {
  if (length(uids) == 0L)
    return(list(chunks = character(0L), n_seq = 0L, n_failed_batches = 0L))
  batches <- split(uids, ceiling(seq_along(uids) / batch_size))
  chunks  <- character(0L)
  n_seq   <- 0L
  n_failed_batches <- 0L
  for (i in seq_along(batches)) {
    fasta <- gb_fetch_batch(
      function(ids) rentrez::entrez_fetch(db = "nuccore", id = ids,
                                          rettype = "fasta", retmode = "text"),
      batches[[i]], "sequence", i, length(batches)
    )
    if (is.null(fasta)) {
      n_failed_batches <- n_failed_batches + 1L
    } else if (nzchar(fasta)) {
      chunks <- c(chunks, fasta)
      n_seq  <- n_seq + sum(startsWith(strsplit(fasta, "\n")[[1L]], ">"))
    }
    Sys.sleep(0.4)
  }
  list(chunks = chunks, n_seq = n_seq, n_failed_batches = n_failed_batches)
}

# ---- Small local helper: confirm, one record at a time, that NCBI really has
# ---- no sequence for a UID ---------------------------------------------------
# A batch fetch that comes up short WITHOUT raising an error has two possible
# causes, and they must not be conflated:
#   (a) records were lost in transit — the defect this script exists to prevent;
#   (b) the record genuinely has no retrievable sequence. GenBank's WGS and TLS
#       "master" records are the usual case: an accession ending 00000000 is a
#       project-level container that indexes component sequences but holds none
#       itself, so a FASTA request for it correctly returns an empty body.
# This helper separates the two by re-requesting each UID ON ITS OWN. A record
# requested individually cannot be a casualty of a batch dropout, so an empty
# response to a single-record request is a property of the record. Anything that
# does come back is a genuine recovery and is returned for appending; anything
# that errors on every attempt is still fatal.
gb_confirm_no_sequence <- function(uids) {
  confirmed <- character(0L)
  chunks    <- character(0L)
  n_seq     <- 0L
  for (k in seq_along(uids)) {
    one <- gb_fetch_batch(
      function(ids) rentrez::entrez_fetch(db = "nuccore", id = ids,
                                          rettype = "fasta", retmode = "text"),
      uids[k], "single-record sequence", k, length(uids)
    )
    if (is.null(one))
      stop("UID ", uids[k], " errored on all ", GB_FETCH_MAX_ATTEMPTS,
           " attempts during single-record confirmation. This is a fetch ",
           "failure, not an absent sequence — re-run when NCBI is responsive.",
           call. = FALSE)
    if (!nzchar(trimws(one))) {
      confirmed <- c(confirmed, uids[k])
    } else {
      chunks <- c(chunks, one)
      n_seq  <- n_seq + sum(startsWith(strsplit(one, "\n")[[1L]], ">"))
    }
    Sys.sleep(0.4)
  }
  list(confirmed = confirmed, chunks = chunks, n_seq = n_seq)
}

# ---- Small local helper: turn one esummary record into a one-row tibble ------
# Split out of the Step 3 loop so that the first run and the repair path parse
# esummary output through exactly the same code.
gb_parse_esummary <- function(s) {
  subtype <- strsplit(if (!is.null(s$subtype)) s$subtype else "", "\\|")[[1L]]
  subname <- strsplit(if (!is.null(s$subname)) s$subname else "", "\\|")[[1L]]
  get_sub <- function(key) {
    idx <- which(subtype == key)
    if (length(idx)) subname[idx[1L]] else NA_character_
  }
  tibble::tibble(
    uid             = s$uid,          accession = s$caption,
    title           = s$title,        organism  = s$organism,
    taxid           = s$taxid,        seq_length = as.integer(s$slen),
    country_gb      = get_sub("country"),
    lat_lon_gb      = get_sub("lat_lon"),
    collection_date = get_sub("collection_date"),
    isolation_src   = get_sub("isolation_source"),
    host_gb         = get_sub("host")
  )
}

# ---- Small local helper: fetch esummary metadata for a set of UIDs -----------
# Same contract as gb_fetch_sequences(): returns the parsed rows plus the count
# of batches lost outright, and leaves the guard to the caller.
gb_fetch_metadata <- function(uids, batch_size = 200L) {
  if (length(uids) == 0L)
    return(list(meta = tibble::tibble(), n_failed_batches = 0L))
  batches   <- split(uids, ceiling(seq_along(uids) / batch_size))
  meta_list <- vector("list", length(batches))
  n_failed_batches <- 0L
  for (i in seq_along(batches)) {
    summ <- gb_fetch_batch(
      function(ids) rentrez::entrez_summary(db = "nuccore", id = ids),
      batches[[i]], "metadata", i, length(batches)
    )
    if (is.null(summ)) {
      n_failed_batches <- n_failed_batches + 1L
    } else {
      if (inherits(summ, "esummary")) summ <- list(summ)
      meta_list[[i]] <- dplyr::bind_rows(lapply(summ, gb_parse_esummary))
    }
    Sys.sleep(0.4)
  }
  list(meta = dplyr::bind_rows(meta_list), n_failed_batches = n_failed_batches)
}

# ---- Small local helper: the UID -> accession index of the metadata file -----
# Reads ONLY the two columns the completeness bookkeeping needs, and reads them
# as text so that UIDs compare reliably against gb_ids_path (which is read as
# lines of text; read as numbers they would be at the mercy of R's default
# 15-significant-digit formatting).
gb_meta_index <- function(path) {
  empty <- tibble::tibble(uid = character(0L), accession = character(0L))
  if (!file.exists(path)) return(empty)
  idx <- readr::read_csv(
    path, show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols_only(uid       = readr::col_character(),
                                 accession = readr::col_character())
  )
  if (nrow(idx) == 0L) empty else idx
}

# Note: canonicalize_host() is defined in 00_setup.R and applied at Step 6
# below to the GenBank `host_taxon` column. The same function is used at
# consumption time by 18_eltonian.R, 17_hutchinsonian.R, and
# 19_sampling_maps.R for the GlobalFungi host fields.

# ---- Report last fetch (if log exists) --------------------------------------

if (file.exists(gb_fetch_log)) {
  writeLines(paste(" ", readLines(gb_fetch_log)))
}

# ---- Step 1: Search GenBank and save UIDs -----------------------------------

search_query <- paste(
  #  'ectomyco*[All Fields]',
  #  'AND "Canada"[Country]',
  '"Canada"[Country]',
  'AND "Fungi"[Organism]',
  'AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])'
)

if (!file.exists(gb_ids_path)) {

  initial <- rentrez::entrez_search(
    db = "nuccore", term = search_query, retmax = 0L, use_history = TRUE
  )
  n_total <- initial$count

  all_ids   <- character(0L)
  batch_size <- 500L
  n_batches  <- ceiling(n_total / batch_size)

  n_failed_batches <- 0L

  for (i in seq_len(n_batches)) {
    batch <- gb_fetch_batch(
      function(ids) rentrez::entrez_search(
        db = "nuccore", term = search_query,
        retmax   = batch_size, retstart = (i - 1L) * batch_size
      ),
      NULL, "UID search", i, n_batches
    )
    if (is.null(batch)) n_failed_batches <- n_failed_batches + 1L
    else                all_ids <- c(all_ids, batch$ids)
    Sys.sleep(0.4)
  }

  all_ids <- unique(all_ids)

  # Completeness guard. The UID list is the reference every later completeness
  # check is measured against, so a short list here would quietly redefine what
  # "complete" means for the whole pipeline. Nothing is written until this
  # passes, so a failed search cannot leave a truncated gb_ids_path behind for a
  # later run to adopt as the frozen list.
  gb_assert_fetch_complete(length(all_ids), n_total, n_failed_batches,
                           "UID search")

  writeLines(all_ids, gb_ids_path)

  writeLines(c(
    paste("Fetch timestamp :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Search query    :", search_query),
    paste("Total hits      :", n_total),
    paste("Unique UIDs     :", length(all_ids))
  ), gb_fetch_log)

} else {
  all_ids <- readLines(gb_ids_path)
}

# Guard: a 0-hit search is almost always a malformed query (e.g. an
# unbalanced quote in search_query) rather than a genuine absence of
# records, and silently proceeding produces empty checkpoints at every
# downstream step until Step 5's coverage check fails on NaN. Fail loudly
# here instead. If GenBank really has 0 matching records, delete
# gb_ids_path and re-run to get past this guard.
if (length(all_ids) == 0L) {
  stop(
    "GenBank search returned 0 UIDs — likely a malformed search_query.\n",
    "  Query used: ", search_query, "\n",
    "  Check quoting/syntax above before re-running. If 0 hits is ",
    "genuinely correct, delete ", gb_ids_path, " and re-run to proceed."
  )
}

# UIDs confirmed by an EARLIER run to have no retrievable sequence (see
# gb_confirm_no_sequence() above). They are excluded from the "still missing"
# sets so that a permanently sequence-less record is not re-requested on every
# run. The ceiling check at the end of Step 3b still audits the whole list, so
# this cannot become a quiet dumping ground.
gb_known_no_seq <- if (file.exists(gb_nosequence_path))
  as.character(readr::read_csv(gb_nosequence_path, show_col_types = FALSE,
                               progress = FALSE)$uid) else character(0L)

# ---- Step 2: Batch-fetch FASTA sequences ------------------------------------
#
# Two paths, chosen automatically from what is already on disk:
#
#   FIRST RUN  (no FASTA checkpoint) — fetch a sequence for every UID.
#   REPAIR RUN (FASTA checkpoint present) — fetch ONLY the sequences missing
#     from it and APPEND them. Sequences already held are never re-requested.
#     That restraint is deliberate and matters for more than speed: re-fetching
#     a record risks silently picking up a sequence or annotation that its
#     submitter has revised since the original pull, which would turn the
#     dataset into a mixture of two retrieval dates and make any before/after
#     comparison uninterpretable.
#
# Sequences are requested by UID but arrive keyed by ACCESSION, so working out
# which UIDs still lack a sequence needs the UID -> accession map, which lives
# in the metadata checkpoint built by Step 3. On a first run that map does not
# exist yet; on a repair run it may itself be short. Step 2 therefore repairs
# every gap it can currently resolve and DEFERS any UID whose accession is not
# yet known. Step 3b below revisits the deferred set once metadata is complete.

if (!file.exists(gb_fasta_path)) {

  message("Step 2: no sequence checkpoint — fetching all ",
          format(length(all_ids), big.mark = ","), " records.")
  seq_res <- gb_fetch_sequences(all_ids)

  # Write to a ".partial" sibling and rename only after the guard has passed.
  # This is what stops a genuinely failed run from leaving a truncated file
  # that the if (!file.exists(...)) test above would accept as complete on the
  # next run — the exact mechanism by which the original loss became permanent.
  fasta_partial <- paste0(gb_fasta_path, ".partial")
  con <- file(fasta_partial, open = "w")
  for (chunk in seq_res$chunks) cat(chunk, file = con)
  close(con)

  # A batch lost outright is always fatal, and the .partial file is left behind
  # rather than promoted to a checkpoint. An error-free shortfall is NOT fatal
  # here: it is left to Step 3b, which re-requests each residual record on its
  # own to distinguish a record that has no sequence from one that was dropped,
  # and then applies the hard completeness gate.
  if (seq_res$n_failed_batches > 0L)
    gb_assert_fetch_complete(seq_res$n_seq, length(all_ids),
                             seq_res$n_failed_batches, "sequence", fasta_partial)
  file.rename(fasta_partial, gb_fasta_path)

} else {

  # ---- Repair path ---------------------------------------------------------
  have_acc      <- gb_fasta_accessions(gb_fasta_path)
  meta_idx      <- gb_meta_index(gb_meta_path)
  resolvable    <- meta_idx[meta_idx$uid %in% all_ids, ]
  need_seq_uids <- setdiff(resolvable$uid[!resolvable$accession %in% have_acc],
                           gb_known_no_seq)
  deferred_uids <- setdiff(all_ids, meta_idx$uid)

  if (length(need_seq_uids) > 0L || length(deferred_uids) > 0L) {
    message("Step 2: sequence checkpoint holds ",
            format(length(have_acc), big.mark = ","), " of ",
            format(length(all_ids), big.mark = ","), " records.")
    message("  missing sequences to fetch now : ",
            format(length(need_seq_uids), big.mark = ","))
    message("  deferred (accession not yet known, revisited at Step 3b): ",
            format(length(deferred_uids), big.mark = ","))
  }

  if (length(need_seq_uids) > 0L) {
    seq_res <- gb_fetch_sequences(need_seq_uids)
    # Guard BEFORE writing anything, and only on a batch lost outright (see the
    # note in the first-run branch above). An error-free shortfall falls through
    # to Step 3b for per-record confirmation.
    if (seq_res$n_failed_batches > 0L)
      gb_assert_fetch_complete(seq_res$n_seq, length(need_seq_uids),
                               seq_res$n_failed_batches, "sequence top-up")

    # Append to a copy, then swap it in. The existing sequences are copied
    # byte-for-byte and never rewritten, and a failure part-way through leaves
    # the original checkpoint untouched.
    fasta_partial <- paste0(gb_fasta_path, ".partial")
    if (!file.copy(gb_fasta_path, fasta_partial, overwrite = TRUE))
      stop("Could not stage ", fasta_partial, " for the sequence top-up.",
           call. = FALSE)
    con <- file(fasta_partial, open = "a")
    for (chunk in seq_res$chunks) cat(chunk, file = con)
    close(con)
    file.rename(fasta_partial, gb_fasta_path)

    if (seq_res$n_seq > 0L)
      message("  appended ", format(seq_res$n_seq, big.mark = ","), " sequences.")
  }
}

# ---- Step 3: Batch-fetch metadata via esummary ------------------------------
#
# Same two paths, and the same restraint, as Step 2: on a repair run only the
# UIDs with no row in the metadata checkpoint are fetched, and their rows are
# appended. Rows already present are left exactly as they were.

if (!file.exists(gb_meta_path)) {

  message("Step 3: no metadata checkpoint — fetching all ",
          format(length(all_ids), big.mark = ","), " records.")
  meta_res <- gb_fetch_metadata(all_ids)
  # Fatal on a batch lost outright; an error-free shortfall is judged by the
  # completeness gate at the end of Step 3b, which tolerates up to
  # GB_MAX_MISSING_FRAC for records withdrawn from GenBank since the search.
  if (meta_res$n_failed_batches > 0L)
    gb_assert_fetch_complete(nrow(meta_res$meta), length(all_ids),
                             meta_res$n_failed_batches, "metadata")

  meta_partial <- paste0(gb_meta_path, ".partial")
  readr::write_csv(meta_res$meta, meta_partial)
  file.rename(meta_partial, gb_meta_path)

} else {

  # ---- Repair path ---------------------------------------------------------
  have_uid       <- gb_meta_index(gb_meta_path)$uid
  need_meta_uids <- setdiff(all_ids, have_uid)

  if (length(need_meta_uids) > 0L) {
    message("Step 3: metadata checkpoint holds ",
            format(length(have_uid), big.mark = ","), " of ",
            format(length(all_ids), big.mark = ","), " records.")
    message("  missing metadata rows to fetch : ",
            format(length(need_meta_uids), big.mark = ","))

    meta_res <- gb_fetch_metadata(need_meta_uids)
    if (meta_res$n_failed_batches > 0L)
      gb_assert_fetch_complete(nrow(meta_res$meta), length(need_meta_uids),
                               meta_res$n_failed_batches, "metadata top-up")

    # gb_parse_esummary() emits the columns in the same order as the existing
    # header, so the new rows can be appended without a header and without
    # re-serializing (and thereby risking a change to) the rows already there.
    meta_partial <- paste0(gb_meta_path, ".partial")
    if (!file.copy(gb_meta_path, meta_partial, overwrite = TRUE))
      stop("Could not stage ", meta_partial, " for the metadata top-up.",
           call. = FALSE)
    readr::write_csv(meta_res$meta, meta_partial, append = TRUE)
    file.rename(meta_partial, gb_meta_path)

    message("  appended ", format(nrow(meta_res$meta), big.mark = ","), " rows.")
  }
}

# ---- Step 3b: Reconcile the two checkpoints and assert completeness ---------
#
# Runs on every invocation, and is the step that makes the pair of checkpoints
# self-healing. Two jobs:
#
#   1. Now that the metadata is complete, every UID's accession is known, so any
#      sequence gap Step 2 had to defer can finally be identified and filled.
#   2. Assert that both checkpoints now account for the full frozen UID list,
#      with no duplicates. These are the checks whose absence let 3,001 records
#      go missing unnoticed in the original pull: a checkpoint that does not
#      reach the UID list is not allowed to be used.
#
# Cost when everything is already complete is a couple of seconds (reading the
# two checkpoints), which is a fair price for the guarantee.

meta_idx      <- gb_meta_index(gb_meta_path)
have_acc      <- gb_fasta_accessions(gb_fasta_path)
still_missing <- setdiff(
  meta_idx$uid[meta_idx$uid %in% all_ids & !meta_idx$accession %in% have_acc],
  gb_known_no_seq)

gb_append_sequences <- function(chunks) {
  if (length(chunks) == 0L) return(invisible(NULL))
  fasta_partial <- paste0(gb_fasta_path, ".partial")
  if (!file.copy(gb_fasta_path, fasta_partial, overwrite = TRUE))
    stop("Could not stage ", fasta_partial, " for the top-up.", call. = FALSE)
  con <- file(fasta_partial, open = "a")
  for (chunk in chunks) cat(chunk, file = con)
  close(con)
  file.rename(fasta_partial, gb_fasta_path)
}

if (length(still_missing) > 0L) {
  message("Step 3b: ", format(length(still_missing), big.mark = ","),
          " sequences still missing once metadata was complete — fetching.")
  seq_res <- gb_fetch_sequences(still_missing)

  # A batch that errored on every attempt is fatal, full stop — there is no
  # legitimate reading of that, and tolerating it is what caused the original
  # loss. Only an error-free shortfall is allowed to proceed to confirmation.
  if (seq_res$n_failed_batches > 0L)
    gb_assert_fetch_complete(seq_res$n_seq, length(still_missing),
                             seq_res$n_failed_batches, "deferred sequence top-up")

  gb_append_sequences(seq_res$chunks)
  have_acc <- gb_fasta_accessions(gb_fasta_path)
  if (seq_res$n_seq > 0L)
    message("  appended ", format(seq_res$n_seq, big.mark = ","), " sequences.")

  # Whatever is STILL absent is now re-requested one record at a time, so that
  # "GenBank has no sequence for this record" can be told apart from "the fetch
  # dropped it". See gb_confirm_no_sequence() above.
  residual <- setdiff(
    meta_idx$uid[meta_idx$uid %in% all_ids & !meta_idx$accession %in% have_acc],
    gb_known_no_seq)
  if (length(residual) > 0L) {
    message("  ", format(length(residual), big.mark = ","),
            " record(s) returned no sequence and no error — confirming each ",
            "individually.")
    conf <- gb_confirm_no_sequence(residual)

    gb_append_sequences(conf$chunks)
    have_acc <- gb_fasta_accessions(gb_fasta_path)
    if (conf$n_seq > 0L)
      message("  recovered ", format(conf$n_seq, big.mark = ","),
              " on the individual retry.")

    if (length(conf$confirmed) > 0L) {
      # Persist rather than merely warn, so the exclusion is auditable — the
      # same treatment given to genus-ambiguous ties (gb_ambiguous_path) and to
      # unmatched SH codes in 02_globalfungi.R.
      no_seq_new <- gb_meta_index(gb_meta_path) |>
        dplyr::filter(uid %in% conf$confirmed) |>
        dplyr::left_join(
          readr::read_csv(gb_meta_path, show_col_types = FALSE,
                          progress = FALSE,
                          col_types = readr::cols_only(
                            uid      = readr::col_character(),
                            organism = readr::col_character(),
                            title    = readr::col_character())),
          by = "uid"
        ) |>
        dplyr::mutate(confirmed_on = format(Sys.Date()))

      # Accumulate: a record confirmed on an earlier run must stay in the file,
      # otherwise it would drop out of gb_known_no_seq and be re-requested (and
      # re-confirmed) forever.
      no_seq_tbl <- if (file.exists(gb_nosequence_path))
        dplyr::bind_rows(
          readr::read_csv(gb_nosequence_path, show_col_types = FALSE,
                          progress = FALSE,
                          col_types = readr::cols(.default = readr::col_character())),
          no_seq_new
        ) |> dplyr::distinct(uid, .keep_all = TRUE) else no_seq_new

      readr::write_csv(no_seq_tbl, gb_nosequence_path)
      message("  ", nrow(no_seq_new), " record(s) confirmed to have NO ",
              "retrievable sequence; logged to ", basename(gb_nosequence_path),
              " (", nrow(no_seq_tbl), " total)")
    }
  }
}

# Records confirmed to have no sequence are excluded from the completeness
# requirement below — but only after the individual confirmation above, and only
# up to a hard ceiling. If a large number of records suddenly claim to have no
# sequence, that is a systemic problem, not a handful of master records, and it
# must not be laundered through this exemption.
n_no_sequence <- if (file.exists(gb_nosequence_path))
  nrow(readr::read_csv(gb_nosequence_path, show_col_types = FALSE,
                       progress = FALSE)) else 0L

if (n_no_sequence > length(all_ids) * GB_MAX_MISSING_FRAC)
  stop(sprintf(
    paste0(
      "%s records claim to have no retrievable sequence (ceiling: %s, %.3f%% ",
      "of the frozen UID list). A handful of WGS/TLS master records is normal; ",
      "this many indicates a systemic retrieval problem. Inspect %s."
    ),
    format(n_no_sequence, big.mark = ","),
    format(floor(length(all_ids) * GB_MAX_MISSING_FRAC), big.mark = ","),
    100 * GB_MAX_MISSING_FRAC, gb_nosequence_path), call. = FALSE)

# --- Final completeness assertions on both checkpoints ----------------------
gb_report_gap <- function(n_have, what) sprintf(
  "  %-28s %s of %s (%s short)", what,
  format(n_have, big.mark = ","), format(length(all_ids), big.mark = ","),
  format(length(all_ids) - n_have, big.mark = ",")
)

seq_accounted <- length(have_acc) + n_no_sequence

if (seq_accounted < length(all_ids) * (1 - GB_MAX_MISSING_FRAC) ||
    nrow(meta_idx) < length(all_ids) * (1 - GB_MAX_MISSING_FRAC)) {
  stop(paste0(
    "GenBank checkpoints do not account for the frozen UID list in\n  ",
    gb_ids_path, "\n",
    gb_report_gap(seq_accounted,  "sequences accounted for"), "\n",
    gb_report_gap(nrow(meta_idx), "rows in metadata CSV"), "\n",
    "  Tolerance is ", sprintf("%.3f%%", 100 * GB_MAX_MISSING_FRAC),
    " of the UID list.\n",
    "  Re-run this script: Steps 2, 3 and 3b fetch only what is still missing.\n",
    "  If some UIDs have genuinely been withdrawn from GenBank, confirm that\n",
    "  individually before relaxing GB_MAX_MISSING_FRAC."
  ), call. = FALSE)
}

if (anyDuplicated(have_acc) > 0L)
  stop("Duplicate accessions in ", gb_fasta_path, " (",
       sum(duplicated(have_acc)), " duplicated). A top-up has been applied ",
       "twice, or a batch overlapped. Restore the checkpoint and re-run.",
       call. = FALSE)

if (anyDuplicated(meta_idx$uid) > 0L || anyDuplicated(meta_idx$accession) > 0L)
  stop("Duplicate records in ", gb_meta_path, " (",
       sum(duplicated(meta_idx$uid)), " duplicated uid, ",
       sum(duplicated(meta_idx$accession)), " duplicated accession). ",
       "Restore the checkpoint and re-run.", call. = FALSE)

message("Step 3b: checkpoints complete — ",
        format(length(have_acc), big.mark = ","), " sequences",
        if (n_no_sequence > 0L)
          paste0(" (+ ", n_no_sequence, " confirmed to have none)") else "",
        ", ", format(nrow(meta_idx), big.mark = ","), " metadata rows, ",
        "against ", format(length(all_ids), big.mark = ","), " frozen UIDs.")

# Loaded once, after reconciliation, so Step 6 always sees the complete table.
gb_meta <- readr::read_csv(gb_meta_path, show_col_types = FALSE)

# ---- Step 4: Extract the ITS2 region with ITSx ------------------------------
#
# GenBank sequences carry no marker-subregion metadata, so — unlike GlobalFungi,
# whose records are already restricted to the ITS2 (or ITSboth) barcoding
# region — the text search sweeps in ITS1-only, ITS2-only and full-ITS records
# together. A length-based proxy cannot reliably separate ITS1-only from ITS2
# sequences, because fungal ITS1 is itself commonly > 200 bp. ITSx
# (Bengtsson-Palme et al. 2013) is therefore used to detect and extract the ITS2
# sub-region directly (Fungi HMM profiles) — the GenBank-side analogue of
# GlobalFungi's ITS2/ITSboth restriction. ITSx is run once against the whole
# FASTA (faster than looping in R). See the GenBank methods in
# FACETS/supplemental_materials_SM1_FACETS.qmd ("GenBank records") for the full
# rationale and citations.

if (!file.exists(gb_itsx_its2) || !file.exists(gb_itsx_its1)) {

  itsx_bin <- Sys.which("ITSx")
  if (!nzchar(itsx_bin))
    stop("ITSx not found on PATH. ITSx requires HMMER 3.x.\n",
         "  Install HMMER:  brew install hmmer  (macOS)  or  sudo apt install hmmer  (Ubuntu)\n",
         "  Then download ITSx from https://microbiology.se/software/itsx/ and place the\n",
         "  'ITSx' script (and its bundled HMM-profile directory) on your PATH.")

  # Drop any input sequence longer than ITSX_MAX_INPUT_LENGTH before ITSx
  # (see the constant's comment above for why this is necessary and harmless).
  raw_seqs <- read_fasta_tbl(gb_fasta_path)
  itsx_input <- raw_seqs |> dplyr::filter(nchar(seq) <= ITSX_MAX_INPUT_LENGTH)
  itsx_input_path <- file.path(paths$temp_dir, "genbank_itsx_input.fasta")
  readr::write_lines(paste0(">", itsx_input$header, "\n", itsx_input$seq), itsx_input_path)

  # -t F           : restrict HMM profile search to the Fungi profile set.
  # --save_regions : write the extracted ITS1 and ITS2 fragments as separate
  #                  FASTA files. ITS1 is kept only for the `its_region`
  #                  provenance flag (Step 6) — it is NOT used for SH matching.
  # --cpu 4        : matches the --threads 4 used for vsearch below.
  cmd <- paste(
    shQuote(itsx_bin),
    "-i", shQuote(itsx_input_path),
    "-o", shQuote(gb_itsx_prefix),
    "-t F",
    "--save_regions ITS1,ITS2",
    "--cpu 4",
    "--graphical F",
    "--silent T"
  )
  ret <- system(cmd)
  if (ret != 0L) stop("ITSx failed with exit code ", ret)

}

# ---- Step 4b: Apply the ITS2 minimum-length floor and sanitize before vsearch
# Builds the vsearch query FASTA. Every eligible record is matched on its ITS2
# fragment ALONE, regardless of whether ITS1 was also detected. GlobalFungi
# classifies ITS1 and ITS2 as two separate, independently BLASTn-classified
# pools and never concatenates them (Vetrovsky et al. 2020), so ITS2 fragments
# are matched alone rather than concatenated with ITS1. Note that ITS1
# fragments are NOT independently classified here — a documented simplification
# relative to GlobalFungi (see the GenBank methods in
# FACETS/supplemental_materials_SM1_FACETS.qmd).

if (!file.exists(gb_vsearch_query)) {

  its2_all   <- read_fasta_tbl(gb_itsx_its2)
  its2_pass  <- its2_all |>
    dplyr::filter(nchar(seq) >= ITS2_MIN_LENGTH)

  # Sanitize before writing: vsearch's FASTA reader accepts only standard
  # nucleotide/IUPAC codes and aborts the ENTIRE run (not just the offending
  # record) on any other character. A few raw GenBank deposits carry gap
  # characters ("-", likely alignment-derived consensus sequences), which
  # ITSx passes through unchanged. Any non-IUPAC character is replaced with
  # "N" rather than dropped, so sequence positions/length are preserved.
  its2_pass <- dplyr::mutate(its2_pass,
    seq = gsub("[^ACGTUNRYSWKMBDHV]", "N", toupper(seq))
  )

  readr::write_lines(paste0(">", its2_pass$header, "\n", its2_pass$seq), gb_vsearch_query)

}

# ---- Step 5: UNITE SH assignment via vsearch (98.5% identity) -------------
#
# Threshold matches GlobalFungi's documented SH-assignment criterion (BLASTn,
# >= 98.5% identity; globalfungi.com methods).
# --strand both                 : also search the reverse-complement strand;
#                                 a plus-strand-only search silently drops
#                                 queries whose true best match is on the
#                                 minus strand.
# --maxaccepts 0 --maxrejects 0 : exhaustive search — vsearch's defaults stop
#                                 at the first hit clearing --id, which is not
#                                 guaranteed to be the single best match.
#                                 Slower, but worth the one-off cost.
# --top_hits_only               : report every reference tied for the best
#                                 identity (enables the Step 6 tie resolution).

if (!file.exists(gb_vsearch_path)) {

  vsearch_bin <- Sys.which("vsearch")
  if (!nzchar(vsearch_bin))
    stop("vsearch not found on PATH. Install: brew install vsearch (macOS) or",
         " sudo apt install vsearch (Ubuntu)")

  cmd <- paste(
    shQuote(vsearch_bin),
    "--usearch_global", shQuote(gb_vsearch_query),
    "--db",            shQuote(unite_fasta_path),
    "--id 0.985",
    "--strand both",
    "--maxaccepts 0",
    "--maxrejects 0",
    "--blast6out",     shQuote(gb_vsearch_path),
    "--top_hits_only",
    "--threads 4"
  )
  ret <- system(cmd)
  if (ret != 0L) stop("vsearch failed with exit code ", ret)

}

# ---- Step 6: Parse, resolve ties, filter, annotate and save -----------------

if (!file.exists(paths$gb_long_out)) {

  blast6_cols <- c("query_id", "target_id", "identity", "aln_length",
                   "mismatches", "gap_opens", "q_start", "q_end",
                   "s_start", "s_end", "evalue", "bitscore")
  hits <- readr::read_tsv(gb_vsearch_path, col_names = blast6_cols, show_col_types = FALSE)

  hits <- dplyr::mutate(hits,
    sh_code   = sub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*", "\\1", target_id),
    # query_id looks like "ACCESSION.VERSION|F|ITS2" (ITSx's header suffix);
    # strip both the ITSx suffix and the GenBank version number.
    accession = sub("\\.\\d+\\|.*$", "", query_id)
  )

  sh_lookup <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
  hits      <- dplyr::left_join(hits, sh_lookup, by = "sh_code")

  # ---- Coverage-check gate -------------------------------------------------
  # This join should be ~100% matched by construction: vsearch assigned
  # these SH codes against the same pinned UNITE file that
  # unite_sh_taxonomy.csv was built from. A non-trivial unmatched fraction
  # means that invariant is broken (e.g. the two path entries have drifted
  # apart, or unite_sh_taxonomy.csv is stale relative to the file vsearch
  # actually ran against) — see 00_setup.R.
  pct_unmatched_gb <- mean(is.na(hits$kingdom))
  if (pct_unmatched_gb > SH_MAX_UNMATCHED_FRAC) {
    stop(sprintf(
      paste0(
        "GenBank vsearch hits joined to UNITE taxonomy with %.2f%% ",
        "unmatched (threshold: %.2f%%). This join should be ~100%% ",
        "matched by construction — check that paths$unite_fasta ",
        "and paths$unite_taxonomy (00_setup.R) are both built from the ",
        "same UNITE file."
      ),
      100 * pct_unmatched_gb, 100 * SH_MAX_UNMATCHED_FRAC
    ))
  }

  # Override species/genus with target_id-specific annotation
  parse_tax_field <- function(target_id, prefix) {
    has_match <- grepl(prefix, target_id, fixed = TRUE)
    result    <- rep(NA_character_, length(target_id))
    result[has_match] <- sub(paste0(".*", prefix, "([^;|]*).*"), "\\1",
                             target_id[has_match])
    result
  }
  hits <- hits |>
    dplyr::mutate(
      .sp  = parse_tax_field(target_id, "s__"),
      .gn  = parse_tax_field(target_id, "g__"),
      species = dplyr::if_else(!is.na(.sp) & nzchar(.sp), .sp, species),
      genus   = dplyr::if_else(!is.na(.gn) & nzchar(.gn), .gn, genus)
    ) |>
    dplyr::select(-.sp, -.gn)

  # ---- Species-level tie resolution, with dark-taxa fix and genus fallback --
  # Exhaustive search + --top_hits_only reports every UNITE reference tied for
  # the best identity for a query (several rows per query is common). Each tied
  # hit is reduced to a "tie identity" unit:
  #   - NAMED species (not ending in "_sp"): the species string. Two tied hits
  #     with the same species but different SH codes are NOT ambiguous
  #     (intraspecific ITS variation can split a species across SH clusters).
  #   - DARK ("_sp") species: the SH code itself, per UNITE's convention
  #     (Nilsson et al. 2019; Ryberg & Nilsson 2018) that dark-taxon identity is
  #     defined by SH cluster membership, not by the uninformative "_sp"
  #     placeholder — two tied "_sp" hits from different SH codes ARE different
  #     dark taxa and must not be collapsed together.
  # A record is retained at SH-level resolution if all tied hits share one tie
  # identity. If they disagree there but every tied hit shares a genus, the
  # record is retained at GENUS-only resolution (sh_code and species set to NA)
  # rather than dropped — "climb to the lowest unambiguous rank" (Heeger et al.
  # 2019). Records ambiguous even at genus level are excluded and logged to
  # gb_ambiguous_path. taxonomic_resolution ("sh"/"genus") records which case
  # applies. See the GenBank methods in
  # FACETS/supplemental_materials_SM1_FACETS.qmd for the full rationale.
  tie_summary <- hits |>
    dplyr::mutate(
      is_dark      = grepl("_sp$", species) | is.na(species),
      tie_identity = dplyr::if_else(is_dark, sh_code, species)
    ) |>
    dplyr::group_by(query_id) |>
    dplyr::summarise(n_tied_hits          = dplyr::n(),
                      n_distinct_identity = dplyr::n_distinct(tie_identity),
                      n_distinct_genus    = dplyr::n_distinct(genus),
                      .groups = "drop")

  sh_resolved_ids    <- tie_summary$query_id[tie_summary$n_distinct_identity == 1L]
  genus_resolved_ids <- tie_summary$query_id[tie_summary$n_distinct_identity > 1L &
                                               tie_summary$n_distinct_genus == 1L]
  ambiguous_ids      <- tie_summary$query_id[tie_summary$n_distinct_genus > 1L]

  hits_ambiguous <- hits |> dplyr::filter(query_id %in% ambiguous_ids)
  readr::write_csv(hits_ambiguous, gb_ambiguous_path)

  # ---- Diagnostic: does ITS2-only matching (discarding ITS1) cost resolution?
  # Reporting-only cross-tabulation of taxonomic_resolution (sh/genus/excluded)
  # against whether the accession also had a genuine ITS1 detection. If
  # genus/excluded rates are similar regardless of ITS1 presence, the ITS2-only
  # simplification is not measurably costing resolution; no pipeline behaviour
  # depends on this file.
  its1_accessions <- unique(sub("\\.\\d+\\|.*$", "", read_fasta_tbl(gb_itsx_its1)$header))
  query_accession_lookup <- dplyr::distinct(hits, query_id, accession)
  resolution_diag <- query_accession_lookup |>
    dplyr::mutate(
      resolution = dplyr::case_when(
        query_id %in% sh_resolved_ids    ~ "sh",
        query_id %in% genus_resolved_ids ~ "genus",
        query_id %in% ambiguous_ids      ~ "excluded",
        TRUE ~ NA_character_
      ),
      had_its1 = accession %in% its1_accessions
    ) |>
    dplyr::count(had_its1, resolution, name = "n") |>
    dplyr::group_by(had_its1) |>
    dplyr::mutate(pct_within_its1_group = 100 * n / sum(n)) |>
    dplyr::ungroup()
  readr::write_csv(resolution_diag, gb_its1_diag_path)

  # SH-resolved: keep one representative row per query (lowest SH code, purely
  # for a deterministic pick; all rows in a group share identity and tie_identity).
  hits_sh <- hits |>
    dplyr::filter(query_id %in% sh_resolved_ids) |>
    dplyr::group_by(query_id) |>
    dplyr::slice_min(sh_code, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(taxonomic_resolution = "sh")

  # Genus-resolved: keep one representative row for traceability, but the SH-
  # and species-level calls are genuinely unresolved, so sh_code and species
  # are set to NA. genus is retained since every tied hit agreed on it.
  hits_genus <- hits |>
    dplyr::filter(query_id %in% genus_resolved_ids) |>
    dplyr::group_by(query_id) |>
    dplyr::slice_min(sh_code, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(sh_code = NA_character_, species = NA_character_,
                  taxonomic_resolution = "genus")

  hits <- dplyr::bind_rows(hits_sh, hits_genus)

  if (!exists("gb_meta")) gb_meta <- readr::read_csv(gb_meta_path, show_col_types = FALSE)
  hits <- dplyr::left_join(hits, gb_meta, by = "accession")

  # ITS-region provenance flag, analogous to GlobalFungi's barcoding_region
  # (ITS2 vs ITSboth). Purely informational: every record is matched against
  # UNITE on its ITS2 fragment alone regardless of this flag (see Step 4b). It
  # records whether ITSx also found a genuine standalone ITS1 region in the
  # same source sequence.
  hits <- dplyr::mutate(hits,
    its_region = dplyr::if_else(accession %in% its1_accessions, "ITSboth", "ITS2")
  )

  # Provenance column
  hits <- dplyr::mutate(hits,
    canada_basis = dplyr::case_when(
      !is.na(country_gb) & !is.na(lat_lon_gb) ~ "both",
      !is.na(country_gb) &  is.na(lat_lon_gb) ~ "country_only",
       is.na(country_gb) & !is.na(lat_lon_gb) ~ "coordinates_only",
      TRUE                                     ~ "search_only"
    )
  )

  # Host taxon (structured field first, then regex on isolation_src).
  # The raw extracted value is preserved as `host_taxon_raw` for traceability;
  # `host_taxon` carries the canonicalized form used by all downstream scripts
  # (18_eltonian.R, 17_hutchinsonian.R, 19_sampling_maps.R). See
  # canonicalize_host() in 00_setup.R for the cleaning rules.
  parse_host_from_text <- function(x) {
    pattern   <- "(?:of|on)\\s+([A-Z][a-z]+(?:\\s+[a-z]+)?)"
    has_match <- grepl(pattern, x, perl = TRUE)
    result    <- rep(NA_character_, length(x))
    result[has_match] <- sub(paste0(".*", pattern, ".*"), "\\1", x[has_match], perl = TRUE)
    result
  }
  hits <- dplyr::mutate(hits,
    host_taxon_raw = dplyr::coalesce(
      dplyr::if_else(!is.na(host_gb) & nzchar(trimws(host_gb)), trimws(host_gb), NA_character_),
      parse_host_from_text(isolation_src)
    ),
    host_taxon = canonicalize_host(host_taxon_raw),
    source = "GenBank"
  )

  # NOTE: genbank_emf_canada_long.csv can contain rows with
  # sh_code = NA (genus-resolved rows). Any downstream code that groups, joins
  # or summarises by sh_code must handle NA explicitly — n_distinct(sh_code)
  # here uses na.rm = TRUE so the NA is not counted as its own "distinct SH".
  readr::write_csv(hits, paths$gb_long_out)
}

