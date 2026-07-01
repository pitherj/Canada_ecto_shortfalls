# =============================================================================
# GenBank Sequence Retrieval and SH Assignment
# =============================================================================
# Retrieves EcM ITS sequences from GenBank for Canadian records, assigns them
# to UNITE Species Hypotheses via vsearch at 98.5% identity (matches
# GlobalFungi's documented BLASTn SH-assignment criterion; see Step 4 below).
# Each step is checkpointed — safe to re-run; skips completed steps.
#
# Prerequisites:
#   02_globalfungi.R must have been run (provides UNITE lookup)
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
#   temp/genbank_emf_canada_ids.txt       — GenBank UIDs
#   temp/genbank_emf_canada.fasta         — retrieved FASTA sequences
#   temp/genbank_emf_canada_metadata.csv  — esummary metadata
#   temp/genbank_vsearch_hits.txt         — vsearch blast6 output
#   genbank_emf_canada_long.csv           — final annotated table
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

# UNITE FASTA: pinned reference, not auto-detected (see 00_setup.R and
# docs/unite_sh_code_mismatch_memo.md). This MUST match
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
    "paths$unite_fasta, which must match it) in 00_setup.R.\n",
    "See docs/unite_sh_code_mismatch_memo.md."
  )
}
ts(sprintf("UNITE FASTA (pinned): %s", basename(unite_fasta_path)))

# Maximum tolerated fraction of vsearch hits left unmatched after the UNITE
# taxonomy join (Step 5 below). Unlike the GlobalFungi join in 0-3 (where SH
# codes are pre-assigned externally and a residual mismatch is expected and
# documented), GenBank's SH codes are assigned in-house via vsearch against
# this same file, so the join is self-consistent by construction and should
# be ~100% matched. Any non-trivial unmatched fraction here indicates a
# broken invariant (e.g. unite_fasta_genbank pointing to a different file
# than the one vsearch actually ran against) rather than an expected gap.
SH_MAX_UNMATCHED_FRAC <- 0.01

# Additional path variables
gb_ids_path     <- file.path(paths$temp_dir, "genbank_emf_canada_ids.txt")
gb_fasta_path   <- file.path(paths$temp_dir, "genbank_emf_canada.fasta")
gb_meta_path    <- file.path(paths$temp_dir, "genbank_emf_canada_metadata.csv")
gb_vsearch_path <- file.path(paths$temp_dir, "genbank_vsearch_hits.txt")
gb_fetch_log    <- file.path(paths$temp_dir, "genbank_fetch_log.txt")

# Note: canonicalize_host() is defined in 00_setup.R and applied at Step 5
# below to the GenBank `host_taxon` column. The same function is used at
# consumption time by 19_eltonian.R, 17_hutchinsonian.R, and
# 20_sampling_maps.R for the GlobalFungi host fields.

# ---- Report last fetch (if log exists) --------------------------------------

if (file.exists(gb_fetch_log)) {
  ts("Previous fetch log:")
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
  ts("Step 1: Searching GenBank (nuccore)...")
  ts("  Query:", search_query)

  initial <- rentrez::entrez_search(
    db = "nuccore", term = search_query, retmax = 0L, use_history = TRUE
  )
  n_total <- initial$count
  ts("  Total hits:", n_total)

  all_ids   <- character(0L)
  batch_size <- 500L
  n_batches  <- ceiling(n_total / batch_size)

  for (i in seq_len(n_batches)) {
    ts(sprintf("  Retrieving IDs: batch %d / %d...", i, n_batches))
    batch <- tryCatch(
      rentrez::entrez_search(
        db = "nuccore", term = search_query,
        retmax   = batch_size, retstart = (i - 1L) * batch_size
      ),
      error = function(e) { ts("  WARNING:", conditionMessage(e)); NULL }
    )
    if (!is.null(batch)) all_ids <- c(all_ids, batch$ids)
    Sys.sleep(0.4)
  }

  all_ids <- unique(all_ids)
  writeLines(all_ids, gb_ids_path)
  ts("  Unique UIDs:", length(all_ids), "-> saved to temp/")

  writeLines(c(
    paste("Fetch timestamp :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Search query    :", search_query),
    paste("Total hits      :", n_total),
    paste("Unique UIDs     :", length(all_ids))
  ), gb_fetch_log)

} else {
  all_ids <- readLines(gb_ids_path)
  ts("Step 1: IDs already exist (", length(all_ids), "UIDs) — skipping.")
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

# ---- Step 2: Batch-fetch FASTA sequences ------------------------------------

if (!file.exists(gb_fasta_path)) {
  ts("Step 2: Fetching FASTA sequences (", length(all_ids), "sequences)...")
  batch_size <- 200L
  batches    <- split(all_ids, ceiling(seq_along(all_ids) / batch_size))
  n_fetched  <- 0L
  con        <- file(gb_fasta_path, open = "w")

  for (i in seq_along(batches)) {
    ts(sprintf("  Batch %d / %d...", i, length(batches)))
    fasta <- tryCatch(
      rentrez::entrez_fetch(db = "nuccore", id = batches[[i]],
                            rettype = "fasta", retmode = "text"),
      error = function(e) { ts("  WARNING:", conditionMessage(e)); NULL }
    )
    if (!is.null(fasta) && nzchar(fasta)) {
      cat(fasta, file = con)
      n_fetched <- n_fetched + sum(startsWith(strsplit(fasta, "\n")[[1L]], ">"))
    }
    Sys.sleep(0.4)
  }
  close(con)
  ts("  Sequences written:", n_fetched)

} else {
  n_seqs <- as.integer(system(paste("grep -c '^>'", shQuote(gb_fasta_path)), intern = TRUE))
  ts("Step 2: FASTA already exists (", n_seqs, "sequences) — skipping.")
}

# ---- Step 3: Batch-fetch metadata via esummary ------------------------------

if (!file.exists(gb_meta_path)) {
  ts("Step 3: Fetching metadata via esummary...")
  batch_size <- 200L
  batches    <- split(all_ids, ceiling(seq_along(all_ids) / batch_size))
  meta_list  <- vector("list", length(batches))

  for (i in seq_along(batches)) {
    ts(sprintf("  Batch %d / %d...", i, length(batches)))
    summ <- tryCatch(
      rentrez::entrez_summary(db = "nuccore", id = batches[[i]]),
      error = function(e) { ts("  WARNING:", conditionMessage(e)); NULL }
    )
    if (!is.null(summ)) {
      if (inherits(summ, "esummary")) summ <- list(summ)
      meta_list[[i]] <- dplyr::bind_rows(lapply(summ, function(s) {
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
      }))
    }
    Sys.sleep(0.4)
  }

  gb_meta <- dplyr::bind_rows(meta_list)
  readr::write_csv(gb_meta, gb_meta_path)
  ts("  Saved metadata:", nrow(gb_meta), "records.")

} else {
  gb_meta <- readr::read_csv(gb_meta_path, show_col_types = FALSE)
  ts("Step 3: Metadata already exists (", nrow(gb_meta), "records) — skipping.")
}

# ---- Step 4: UNITE SH assignment via vsearch (98.5% identity) -------------
#
# Threshold matches GlobalFungi's documented SH-assignment criterion ("Each
# extracted ITS1 or ITS2 sequence was assigned to the closest UNITE 10.0
# Species Hypothesis using BLASTn ... if sequence similarity was >= 98.5%.",
# globalfungi.com methods). Previously this was 97% (a common general ITS
# OTU-clustering threshold), which is more permissive than GlobalFungi's own
# criterion and meant the two sources could apply different inclusion rules
# for SH membership. Raised to 98.5% for cross-source consistency.

if (!file.exists(gb_vsearch_path)) {
  ts("Step 4: Assigning sequences to UNITE SHs via vsearch...")

  vsearch_bin <- Sys.which("vsearch")
  if (!nzchar(vsearch_bin))
    stop("vsearch not found on PATH. Install: brew install vsearch (macOS) or",
         " sudo apt install vsearch (Ubuntu)")

  cmd <- paste(
    shQuote(vsearch_bin),
    "--usearch_global", shQuote(gb_fasta_path),
    "--db",            shQuote(unite_fasta_path),
    "--id 0.985",
    "--blast6out",     shQuote(gb_vsearch_path),
    "--top_hits_only",
    "--threads 4"
  )
  ret <- system(cmd)
  if (ret != 0L) stop("vsearch failed with exit code ", ret)
  ts("  vsearch complete:", length(readLines(gb_vsearch_path)), "hits.")

} else {
  ts("Step 4: vsearch output already exists (", length(readLines(gb_vsearch_path)), "hits) — skipping.")
}

# ---- Step 5: Parse, filter, annotate and save --------------------------------

if (!file.exists(paths$gb_long_out)) {
  ts("Step 5: Joining vsearch hits with taxonomy and metadata...")

  blast6_cols <- c("query_id", "target_id", "identity", "aln_length",
                   "mismatches", "gap_opens", "q_start", "q_end",
                   "s_start", "s_end", "evalue", "bitscore")
  hits <- readr::read_tsv(gb_vsearch_path, col_names = blast6_cols, show_col_types = FALSE)

  hits <- dplyr::mutate(hits,
    sh_code   = sub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*", "\\1", target_id),
    accession = sub("\\.\\d+$", "", query_id)
  )

  sh_lookup <- readr::read_csv(paths$unite_taxonomy, show_col_types = FALSE)
  hits      <- dplyr::left_join(hits, sh_lookup, by = "sh_code")

  # ---- Coverage-check gate -------------------------------------------------
  # This join should be ~100% matched by construction: vsearch assigned
  # these SH codes against the same pinned UNITE file (paths$unite_fasta_
  # genbank) that unite_sh_taxonomy.csv was built from. A non-trivial
  # unmatched fraction here means that invariant is broken (e.g. the two
  # path entries have drifted apart, or unite_sh_taxonomy.csv is stale
  # relative to the file vsearch actually ran against) — see 00_setup.R
  # and docs/unite_sh_code_mismatch_memo.md.
  pct_unmatched_gb <- mean(is.na(hits$kingdom))
  ts(sprintf("  UNITE taxonomy coverage (GenBank/vsearch join): %.2f%% unmatched",
             100 * pct_unmatched_gb))
  if (pct_unmatched_gb > SH_MAX_UNMATCHED_FRAC) {
    stop(sprintf(
      paste0(
        "GenBank vsearch hits joined to UNITE taxonomy with %.2f%% ",
        "unmatched (threshold: %.2f%%). This join should be ~100%% ",
        "matched by construction — check that paths$unite_fasta ",
        "and paths$unite_taxonomy (00_setup.R) are both built from the ",
        "same UNITE file. See docs/unite_sh_code_mismatch_memo.md."
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

  if (!exists("gb_meta")) gb_meta <- readr::read_csv(gb_meta_path, show_col_types = FALSE)
  hits <- dplyr::left_join(hits, gb_meta, by = "accession")

  # Quality filters
  hits <- dplyr::filter(hits, !is.na(seq_length), seq_length >= 200L, aln_length >= 200L)
  ts("  After quality filters (seq_length & aln_length >= 200 bp):", nrow(hits), "records")

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
  # (Eltonian, Hutchinsonian, 3-0 host-accumulation). See canonicalize_host()
  # at the top of this script for the cleaning rules.
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

  # Cleanup diagnostics — log totals and example transformations
  n_raw      <- sum(!is.na(hits$host_taxon_raw))
  n_clean    <- sum(!is.na(hits$host_taxon))
  n_dropped  <- sum(!is.na(hits$host_taxon_raw) & is.na(hits$host_taxon))
  n_modified <- sum(!is.na(hits$host_taxon_raw) & !is.na(hits$host_taxon) &
                    hits$host_taxon != hits$host_taxon_raw)
  ts(sprintf("  Host taxon canonicalization: %d raw -> %d valid (%d dropped, %d modified)",
             n_raw, n_clean, n_dropped, n_modified))
  if (n_modified > 0L) {
    examples <- hits |>
      dplyr::filter(!is.na(host_taxon_raw), !is.na(host_taxon),
                    host_taxon != host_taxon_raw) |>
      dplyr::distinct(host_taxon_raw, host_taxon) |>
      head(5L)
    ts("  Example modifications (raw -> canonical):")
    for (i in seq_len(nrow(examples)))
      ts(sprintf("    %s  ->  %s", examples$host_taxon_raw[i], examples$host_taxon[i]))
  }
  if (n_dropped > 0L) {
    dropped_examples <- hits |>
      dplyr::filter(!is.na(host_taxon_raw), is.na(host_taxon)) |>
      dplyr::count(host_taxon_raw, sort = TRUE) |>
      head(5L)
    ts("  Top dropped raw values (set to NA):")
    for (i in seq_len(nrow(dropped_examples)))
      ts(sprintf("    %-40s  (n = %d)",
                 dropped_examples$host_taxon_raw[i], dropped_examples$n[i]))
  }

  readr::write_csv(hits, paths$gb_long_out)
  ts(sprintf("  Saved genbank_emf_canada_long.csv (%d records, %d SHs, %d genera)",
             nrow(hits), dplyr::n_distinct(hits$sh_code), dplyr::n_distinct(hits$genus)))
} else {
  ts("Step 5: genbank_emf_canada_long.csv already exists — skipping.")
}

ts("03_genbank.R complete.  Run 04_combine_ecm_dataset.R next.")
