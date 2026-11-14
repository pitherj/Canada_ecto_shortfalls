# =============================================================================
# Table S4 check: do the "no Canadian sequence data" genera really lack it?
# =============================================================================
# WHY THIS SCRIPT EXISTS
#   09_linnean.R lists FungalTraits EcM genera that have GBIF physical specimens
#   in Canada but do not appear in our Canadian sequence dataset (`emf`). That
#   test compares genus NAMES across three independent naming systems: the GBIF
#   backbone (GBIF records), FungalTraits (the EcM genus list) and UNITE (the
#   genus label each of our sequences receives from its Species Hypothesis).
#   When the names disagree, a genus can look absent even though Canadian
#   sequences exist. A reviewer flagged this (e.g. Bankera = Phellodon;
#   Boletinus = Suillus). This script re-examines each of those genera
#   (32 in the 6 March 2026 download) against three independent lines of
#   evidence, and sorts each genus into one category.
#
# EVIDENCE USED, PER GENUS
#   1. Synonymy. A hand-curated table (data_raw/table_s4_check/
#      table_s4_genus_synonyms.csv) records the Index Fungorum current name of
#      each genus, with Index Fungorum record numbers. The raw Index Fungorum
#      responses it was built from (retrieved 2026-11-12) are archived beside it.
#   2. Our own GenBank download. The name the sequence's submitter gave
#      (GenBank `organism` field) is checked, rather than the UNITE label.
#   3. A targeted GenBank query for the genus, run twice: once restricted to
#      Canada, once worldwide. Search terms mirror 03_genbank.R. GenBank's
#      [Organism] field uses NCBI's taxonomy, which folds in many synonyms; the
#      query translation NCBI applied is saved so this can be inspected.
#
# CATEGORIES (assigned in this order; the first that applies wins)
#   a. "Synonym of a genus in our data"      — current name is a genus in emf
#   b. "Canadian ITS sequences exist"         — evidence 2 or 3 finds at least
#                                               one record whose GenBank country
#                                               field begins with "Canada"
#   c. "ITS sequences exist, none from Canada"— worldwide query > 0
#   d. "No ITS sequences found"               — worldwide query = 0
#   Genera in categories c and d are the genuine "no Canadian sequence" set.
#
# INPUTS
#   paths$gbif_ecm_nosequence  (from 09_linnean.R)
#   paths$gb_long_out          (from 03_genbank.R)
#   emf                        (loaded by 00_setup.R)
#   paths$table_s4_synonyms    (hand-curated, data_raw)
#
# RAW DATA DOWNLOADED (first run only; reused afterwards)
#   paths$table_s4_ncbi_dir/ncbi_query_log.csv  — one row per query: the exact
#       search string, NCBI's translation of it, the hit count, retrieval date
#   paths$table_s4_ncbi_dir/ncbi_query_uids.csv — every GenBank UID returned,
#       with the record's collection-country field for the Canada queries
#   Delete this folder to re-query GenBank (results will then reflect the new
#   retrieval date). Requires ENTREZ_KEY in .Renviron.
#
# OUTPUTS
#   paths$table_s4_reconciliation        — one row per genus, all evidence +
#                                          category (feeds Table S4)
#   paths$gbif_ecm_nosequence_confirmed  — GBIF records for genera in
#                                          categories c and d (feeds Figure S5)
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Step 1: Read inputs -----------------------------------------------------
# The GBIF records 09_linnean.R flagged as "no sequence data", and the list of
# genera present in our sequence dataset (UNITE labels).

gbif_noseq <- readr::read_csv(paths$gbif_ecm_nosequence, show_col_types = FALSE)
genera     <- sort(unique(gbif_noseq$genus))
emf_genera <- unique(emf$genus)

synonyms <- readr::read_csv(paths$table_s4_synonyms, show_col_types = FALSE)

# Every genus under review must have a row in the synonym table, so that no
# genus silently skips the synonymy check. Stop if any is missing.
missing_syn <- setdiff(genera, synonyms$genus_gbif)
if (length(missing_syn) > 0) {
  stop("These genera have no row in ", paths$table_s4_synonyms, ": ",
       paste(missing_syn, collapse = ", "),
       "\n  Look each up in Index Fungorum and add a row before re-running.")
}

# ---- Step 2: Evidence from our own GenBank download --------------------------
# Count Canadian GenBank records whose SUBMITTER named them in each genus (the
# first word of the `organism` field). These records are already in our
# download; they were filed under a different genus by UNITE, which is why the
# genus looked absent. The UNITE labels they received are kept for reporting.

# As in Step 3, a record counts as Canadian only if its own GenBank country
# field (country_gb) begins with "Canada".
gb <- readr::read_csv(paths$gb_long_out, show_col_types = FALSE,
                      col_select = c("accession", "organism", "genus", "country_gb")) |>
  dplyr::mutate(submitter_genus = sub(" .*$", "", organism))

own_gb <- gb |>
  dplyr::filter(submitter_genus %in% genera, grepl("^Canada", country_gb)) |>
  dplyr::group_by(genus_gbif = submitter_genus) |>
  dplyr::summarise(
    n_own_genbank_canada = dplyr::n_distinct(accession),
    unite_labels_given   = paste(sort(unique(genus)), collapse = "; "),
    .groups = "drop")

# ---- Step 3: Targeted GenBank queries (cached) -------------------------------
# For each genus, two queries: Canada only, and worldwide. Terms mirror the
# main search in 03_genbank.R, with the genus added as an [Organism] term.
# "Fungi"[Organism] is kept so that animal or plant genera sharing a name
# (e.g. the dragonfly genus Gomphus) are excluded.

ncbi_log_path  <- file.path(paths$table_s4_ncbi_dir, "ncbi_query_log.csv")
ncbi_uids_path <- file.path(paths$table_s4_ncbi_dir, "ncbi_query_uids.csv")

build_query <- function(genus, scope) {
  q <- paste0('"', genus, '"[Organism] AND "Fungi"[Organism] ',
              'AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])')
  if (scope == "canada") q <- paste(q, 'AND "Canada"[Country]')
  q
}

# run_query(): runs one search and pages through all UIDs (NCBI returns at
# most 10,000 per request). Stops if the number of UIDs collected does not
# equal the hit count NCBI reported, so a partial result is never saved.
run_query <- function(genus, scope) {
  q   <- build_query(genus, scope)
  res <- rentrez::entrez_search(db = "nuccore", term = q, retmax = 0)
  n   <- as.integer(res$count)
  uids <- character(0)
  start <- 0L
  while (start < n) {
    page <- rentrez::entrez_search(db = "nuccore", term = q,
                                   retstart = start, retmax = 10000L)
    uids  <- c(uids, page$ids)
    start <- start + 10000L
    Sys.sleep(0.15)  # stay under NCBI's request-rate limit
  }
  uids <- unique(uids)
  if (length(uids) != n) {
    stop("GenBank query for ", genus, " (", scope, ") returned ", length(uids),
         " UIDs but reported ", n, ". Not saving a partial result.")
  }
  # For the Canada query, also record each hit's collection-country field.
  # GenBank's [Country] search index is broader than the record's own country
  # field: it returns some records collected elsewhere (e.g. a USA holotype)
  # and some with no country recorded. Only hits whose country field begins
  # with "Canada" are later counted as Canadian.
  country <- rep(NA_character_, n)
  if (scope == "canada" && n > 0) {
    country <- get_country(uids)
  }
  list(
    log  = tibble::tibble(genus_gbif = genus, scope = scope, query = q,
                          query_translation = res$QueryTranslation,
                          n_hits = n, retrieved_on = as.character(Sys.Date())),
    uids = tibble::tibble(genus_gbif = rep(genus, n), scope = rep(scope, n),
                          uid = uids, country = country)
  )
}

# get_country(): fetches GenBank summary records (in batches of 200) and pulls
# out the `country` qualifier (collection locality, e.g. "Canada: Quebec").
# Returns NA where the record has no country qualifier. Same parsing approach
# as 03_genbank.R.
get_country <- function(uids) {
  out <- character(0)
  for (b in split(uids, ceiling(seq_along(uids) / 200))) {
    s <- rentrez::entrez_summary(db = "nuccore", id = b, always_return_list = TRUE)
    out <- c(out, vapply(b, function(id) {
      rec <- s[[id]]
      if (is.null(rec)) stop("No GenBank summary returned for UID ", id)
      st <- strsplit(if (!is.null(rec$subtype)) rec$subtype else "", "|", fixed = TRUE)[[1]]
      sn <- strsplit(if (!is.null(rec$subname)) rec$subname else "", "|", fixed = TRUE)[[1]]
      i  <- which(st == "country")
      if (length(i) == 0 || i[1] > length(sn)) NA_character_ else sn[i[1]]
    }, character(1)))
    Sys.sleep(0.15)
  }
  out
}

if (file.exists(ncbi_log_path) && file.exists(ncbi_uids_path)) {
  ncbi_log  <- readr::read_csv(ncbi_log_path,  show_col_types = FALSE)
  ncbi_uids <- readr::read_csv(ncbi_uids_path, show_col_types = FALSE,
                               col_types = readr::cols(.default = "c"))
  # The cache must cover every genus under review in both scopes.
  if (!all(genera %in% ncbi_log$genus_gbif)) {
    stop("Cached GenBank queries in ", paths$table_s4_ncbi_dir,
         " do not cover every genus. Delete the folder to re-query.")
  }
} else {
  if (Sys.getenv("ENTREZ_KEY") == "") {
    stop("NCBI API key not found. Set ENTREZ_KEY in your .Renviron file.")
  }
  rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
  results <- list()
  for (g in genera) for (s in c("canada", "world")) {
    results[[paste(g, s)]] <- run_query(g, s)
  }
  ncbi_log  <- dplyr::bind_rows(lapply(results, `[[`, "log"))
  ncbi_uids <- dplyr::bind_rows(lapply(results, `[[`, "uids"))
  dir.create(paths$table_s4_ncbi_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(ncbi_uids, ncbi_uids_path)
  readr::write_csv(ncbi_log,  ncbi_log_path)   # written last: marks a complete cache
}

# Three counts per genus:
#   n_genbank_its_world            — all fungal ITS records for the genus
#   n_genbank_its_canada_search    — records returned by the Canada search
#   n_genbank_its_canada_confirmed — of those, records whose own country
#                                    field begins with "Canada" (used below)
ncbi_counts <- ncbi_log |>
  dplyr::select(genus_gbif, scope, n_hits) |>
  tidyr::pivot_wider(names_from = scope, values_from = n_hits,
                     names_prefix = "n_genbank_its_") |>
  dplyr::rename(n_genbank_its_canada_search = n_genbank_its_canada) |>
  dplyr::left_join(
    ncbi_uids |>
      dplyr::filter(scope == "canada") |>
      dplyr::group_by(genus_gbif) |>
      dplyr::summarise(n_genbank_its_canada_confirmed =
                         sum(grepl("^Canada", country)), .groups = "drop"),
    by = "genus_gbif") |>
  dplyr::mutate(n_genbank_its_canada_confirmed =
                  dplyr::coalesce(n_genbank_its_canada_confirmed, 0L))

# ---- Step 4: Combine the evidence and assign categories ----------------------

gbif_counts <- gbif_noseq |>
  dplyr::group_by(genus_gbif = genus) |>
  dplyr::summarise(
    gbif_records = dplyr::n(),
    gbif_species = dplyr::n_distinct(species[!is.na(species) & species != ""]),
    .groups = "drop")

category_levels <- c("Synonym of a genus in our data",
                     "Canadian ITS sequences exist",
                     "ITS sequences exist, none from Canada",
                     "No ITS sequences found")

reconciliation <- gbif_counts |>
  dplyr::left_join(dplyr::select(synonyms, genus_gbif, current_genus,
                                 synonym_applied, if_genus_record,
                                 if_current_record, source),
                   by = "genus_gbif") |>
  dplyr::left_join(own_gb,      by = "genus_gbif") |>
  dplyr::left_join(ncbi_counts, by = "genus_gbif") |>
  dplyr::mutate(
    n_own_genbank_canada    = dplyr::coalesce(n_own_genbank_canada, 0L),
    current_genus_in_our_data = current_genus %in% emf_genera,
    category = dplyr::case_when(
      synonym_applied & current_genus_in_our_data           ~ category_levels[1],
      n_own_genbank_canada > 0 |
        n_genbank_its_canada_confirmed > 0                   ~ category_levels[2],
      n_genbank_its_world > 0                                ~ category_levels[3],
      TRUE                                                   ~ category_levels[4]),
    category = factor(category, levels = category_levels)
  ) |>
  dplyr::arrange(category, dplyr::desc(gbif_records))

readr::write_csv(reconciliation, paths$table_s4_reconciliation)

# ---- Step 5: GBIF records for the genuinely "no Canadian sequence" genera ----
# Only categories c and d. Figure S5 plots these as "without sequence data";
# records in categories a and b are plotted with the "with sequence data" set.

confirmed_genera <- reconciliation$genus_gbif[
  reconciliation$category %in% category_levels[3:4]]

gbif_noseq |>
  dplyr::filter(genus %in% confirmed_genera) |>
  readr::write_csv(paths$gbif_ecm_nosequence_confirmed)

# ---- Step 6: Console summary -------------------------------------------------
message("Table S4 check: ", length(genera), " genera reviewed")
print(dplyr::count(reconciliation, category, name = "n_genera"))
