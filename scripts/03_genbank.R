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

  for (i in seq_len(n_batches)) {
    batch <- tryCatch(
      rentrez::entrez_search(
        db = "nuccore", term = search_query,
        retmax   = batch_size, retstart = (i - 1L) * batch_size
      ),
      error = function(e) { warning(conditionMessage(e), call. = FALSE); NULL }
    )
    if (!is.null(batch)) all_ids <- c(all_ids, batch$ids)
    Sys.sleep(0.4)
  }

  all_ids <- unique(all_ids)
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

# ---- Step 2: Batch-fetch FASTA sequences ------------------------------------

if (!file.exists(gb_fasta_path)) {
  batch_size <- 200L
  batches    <- split(all_ids, ceiling(seq_along(all_ids) / batch_size))
  n_fetched  <- 0L
  con        <- file(gb_fasta_path, open = "w")

  for (i in seq_along(batches)) {
    fasta <- tryCatch(
      rentrez::entrez_fetch(db = "nuccore", id = batches[[i]],
                            rettype = "fasta", retmode = "text"),
      error = function(e) { warning(conditionMessage(e), call. = FALSE); NULL }
    )
    if (!is.null(fasta) && nzchar(fasta)) {
      cat(fasta, file = con)
      n_fetched <- n_fetched + sum(startsWith(strsplit(fasta, "\n")[[1L]], ">"))
    }
    Sys.sleep(0.4)
  }
  close(con)

}

# ---- Step 3: Batch-fetch metadata via esummary ------------------------------

if (!file.exists(gb_meta_path)) {
  batch_size <- 200L
  batches    <- split(all_ids, ceiling(seq_along(all_ids) / batch_size))
  meta_list  <- vector("list", length(batches))

  for (i in seq_along(batches)) {
    summ <- tryCatch(
      rentrez::entrez_summary(db = "nuccore", id = batches[[i]]),
      error = function(e) { warning(conditionMessage(e), call. = FALSE); NULL }
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

} else {
  gb_meta <- readr::read_csv(gb_meta_path, show_col_types = FALSE)
}

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

