# =============================================================================
# Darwinian Shortfall, part 2: genome coverage by Species Hypothesis (SH)
# =============================================================================
# PURPOSE
#   15_darwinian.R matches MycoCosm genomes to our EcM species by NAME. That
#   cannot say anything about the ~1,190 Canadian SHs that UNITE leaves
#   un-named (labels such as "Cortinarius_sp"), because an un-named SH has no
#   species name to match. A reviewer asked how many of the un-named units have
#   genome data. Name matching can also miss a genome whose species UNITE files
#   under an un-named SH (e.g. Cenococcum geophilum has a MycoCosm genome, but
#   our Canadian data contain only "Cenococcum_sp").
#
#   This script matches by SEQUENCE instead. For each MycoCosm genome in an EcM
#   genus it finds an ITS sequence deposited in GenBank for the SAME STRAIN,
#   extracts the ITS2 region, and assigns it to a UNITE SH with exactly the
#   rule used for our GenBank data in 03_genbank.R (ITSx, then vsearch at
#   98.5% identity against the same pinned UNITE file). A Canadian SH then
#   "has a genome" if at least one genome's ITS2 falls into it. Named and
#   un-named SHs are therefore judged by one identical rule.
#
# HOW THE GENOME-TO-ITS LINK IS FOUND (and kept between runs)
#   The link is stored in a lookup table, one row per genome:
#     data_raw/mycocosm_its/mycocosm_its_lookup.csv
#   On each run, only genomes NOT already in the table are searched for in
#   GenBank. Rows already present are reused unchanged, so results are stable
#   between runs and a twice-yearly update only queries newly added genomes.
#
#   Automatic search, per genome:
#     1. Split the MycoCosm name into genus, epithet and strain code
#        ("Laccaria bicolor S238N v1.0" -> Laccaria / bicolor / S238N).
#     2. Query GenBank: organism + strain code + ITS terms.
#     3. Accept a hit ONLY if the record's own strain / isolate / voucher /
#        culture-collection field matches the strain code (see
#        strain_matches() for the exact rule). The organism name alone is never
#        enough: that would just be name matching again.
#     4. If several records pass, prefer those whose title says they include
#        ITS2, then the longest.
#
# GENOMES NOT ASSESSED
#   Genomes the automatic search cannot link to an SH (no same-strain ITS
#   record, no ITS2 region, no UNITE match, or a tie between SHs) are written,
#   with the reason, to paths$darwinian_not_assessed. They are reported as not
#   assessed; nothing is resolved by hand, so every update runs unattended.
#   To force a genome to be searched again, delete its row from the lookup.
#
# INPUTS
#   paths$mycocosm_list   — MycoCosm genome list (copied by hand; see 15_darwinian.R)
#   paths$unite_fasta     — pinned UNITE reference FASTA (same file as 03_genbank.R)
#   emf                   — Canadian EcM dataset (loaded by 00_setup.R)
#   ENTREZ_KEY in .Renviron (only needed when new genomes must be searched)
#   ITSx (+ HMMER) and vsearch on PATH (see 03_genbank.R header for install)
#
# RAW DATA WRITTEN (cached; reused on later runs)
#   data_raw/mycocosm_its/mycocosm_its_lookup.csv      — genome -> ITS accession
#   data_raw/mycocosm_its/mycocosm_its_sequences.fasta — downloaded ITS sequences
#
# OUTPUTS (data_derived/darwinian/)
#   darwinian_genome_sh_assignments.csv — one row per genome: ITS accession,
#                                         outcome, and the SH code(s) credited
#   darwinian_sh_coverage.csv           — one row per Canadian SH: named or
#                                         un-named, and whether it has a genome
#   darwinian_sh_summary.csv            — summary counts (feeds the manuscript)
#   darwinian_genomes_not_assessed.csv  — genomes not linked to an SH, with
#                                         the reason
#
# MANUSCRIPT ITEMS IT FEEDS
#   The SH-level Darwinian result (main text, beside the named-species
#   result) and the reviewer response on un-named SHs.
#
# RE-RUN WHEN
#   The EcM dataset changes (04), the MycoCosm list is refreshed (new genomes
#   are searched automatically), or the UNITE build is repinned. Outputs are
#   rewritten on every run.
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Settings ---------------------------------------------------------------
# These match 03_genbank.R so that genome ITS sequences are treated exactly as
# our Canadian GenBank sequences were.
ITS2_MIN_LENGTH <- 100L     # ITS2 fragments shorter than this are not used
SH_IDENTITY     <- 0.985    # vsearch identity threshold for SH assignment

# Automatic GenBank search limits.
#   STRAIN_MIN_CHARS: strain codes shorter than this (after removing spaces and
#     punctuation, e.g. "11") match too many unrelated records to be trusted;
#     those genomes are reported as not assessed.
#   MAX_HITS_PER_QUERY: if GenBank cannot find the strain phrase it may quietly
#     drop that term and return every ITS record for the species. A hit count
#     above this limit is treated as that failure, and the genome is reported
#     as not assessed rather than downloading thousands of summaries.
#   MAX_ITS_LENGTH: records longer than this are whole-genome or scaffold
#     records, not ITS deposits, and are ignored.
STRAIN_MIN_CHARS   <- 3L
MAX_HITS_PER_QUERY <- 500L
MAX_ITS_LENGTH     <- 5000L

# File locations (paths list in 00_setup.R). The lookup and sequences are
# raw data (live-query results), so they live in data_raw/, in their own
# folder because data_raw/mycocosm/ is write-protected.
its_dir        <- paths$mycocosm_its_dir
lookup_path    <- paths$mycocosm_its_lookup
seq_cache_path <- paths$mycocosm_its_seqs
dir.create(its_dir, showWarnings = FALSE, recursive = TRUE)

out_assign  <- paths$darwinian_genome_sh
out_cover   <- paths$darwinian_sh_coverage
out_summary <- paths$darwinian_sh_summary
out_review  <- paths$darwinian_not_assessed

# read_fasta_tbl() and write_fasta_tbl() are defined in 00_setup.R.

# ---- Small helper: normalize a strain code for comparison --------------------
# Strain codes are written inconsistently ("CBS 559.96", "CBS559.96",
# "cbs 559-96"). Lower-casing and removing everything except letters and
# digits makes these compare equal ("cbs55996").
norm_strain <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# ---- Small helper: does a GenBank record's strain field match? ---------------
# `fields` holds the record's strain, isolate, specimen_voucher and
# culture_collection values. A match is EITHER
#   (a) a field equal to the strain code after normalization, OR
#   (b) for codes of 5+ characters, a field that contains the code (or is
#       contained in it), e.g. "S238N" inside "S238N-H82".
# Rule (b) is not applied to short codes, where "109" would wrongly match
# "1094". The record's title is deliberately not searched, for the same reason.
strain_matches <- function(fields, strain) {
  s <- norm_strain(strain)
  f <- norm_strain(fields)
  f <- f[nzchar(f)]
  if (length(f) == 0L) return(FALSE)
  any(f == s) ||
    (nchar(s) >= 5L && any(grepl(s, f, fixed = TRUE) |
                             vapply(f, function(v) nchar(v) >= 5L &&
                                      grepl(v, s, fixed = TRUE), logical(1))))
}

# =============================================================================
# Step 1: MycoCosm genomes in EcM genera, with name parts split out
# =============================================================================
# The genus universe is the set of genera in our Canadian dataset, as in
# 15_darwinian.R, so both scripts start from the same 186-or-so genomes.
#
# Genomes are identified by their full MycoCosm name (taxon_name), not by the
# `id` column: `id` is only the row number from the manual copy-and-paste, and
# changes whenever the list is refreshed.

mc <- readr::read_csv(paths$mycocosm_list, show_col_types = FALSE)
names(mc) <- tolower(gsub("\\s+", "_", names(mc)))

our_genera <- unique(emf$genus[!is.na(emf$genus)])

genomes <- mc |>
  dplyr::transmute(
    # Remove invisible "zero-width" characters copied from the web page, then
    # collapse repeated spaces.
    taxon_name = gsub("\\s+", " ", trimws(gsub("[​‌‍﻿]", "",
                                               taxon_name)))
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    # Drop the assembly version at the end ("v1.0", "v2.0", "v1.1").
    base  = sub("\\s+v[0-9]+(\\.[0-9]+)*$", "", taxon_name),
    genus = sub("^(\\S+).*", "\\1", base),
    rest  = trimws(sub("^\\S+", "", base))
  ) |>
  dplyr::filter(genus %in% our_genera)

# Split the remainder into "epithet" and "strain". The epithet is the first
# word, plus a second word when the first is a qualifier (aff., cf.) or when a
# rank marker (var., ssp., subsp.) follows. "sp." / "sp. nov." genomes have no
# epithet: everything after them is the strain code.
split_name <- function(rest) {
  w <- strsplit(rest, " ", fixed = TRUE)[[1]]
  if (length(w) == 0L) return(c(epithet = NA, strain = NA))
  if (w[1] == "sp.") {
    k <- if (length(w) >= 2 && w[2] == "nov.") 2L else 1L
    return(c(epithet = NA, strain = paste(w[-seq_len(k)], collapse = " ")))
  }
  k <- 1L
  if (w[1] %in% c("aff.", "cf.") && length(w) >= 2) k <- 2L
  if (length(w) >= k + 2 && w[k + 1] %in% c("var.", "ssp.", "subsp.", "subspA")) {
    k <- k + 2L
  }
  c(epithet = paste(w[seq_len(k)], collapse = " "),
    strain  = paste(w[-seq_len(k)], collapse = " "))
}

parts <- t(vapply(genomes$rest, split_name, character(2)))
genomes <- genomes |>
  dplyr::mutate(
    epithet = unname(parts[, "epithet"]),
    strain  = dplyr::na_if(trimws(unname(parts[, "strain"])), "")
  ) |>
  dplyr::select(taxon_name, genus, epithet, strain)

message("Step 1: ", nrow(genomes), " MycoCosm genomes in EcM genera present in Canada.")

# =============================================================================
# Step 2: Find a same-strain ITS accession for each genome not yet in the lookup
# =============================================================================

lookup_cols <- c("taxon_name", "genus", "epithet", "strain", "query",
                 "query_translation", "n_hits", "its_accession", "its_length",
                 "status", "source", "retrieved_on", "notes")

if (file.exists(lookup_path)) {
  lookup <- readr::read_csv(lookup_path, show_col_types = FALSE,
                            col_types = readr::cols(.default = "c"))
} else {
  lookup <- tibble::as_tibble(setNames(
    replicate(length(lookup_cols), character(0), simplify = FALSE), lookup_cols))
}

# search_one(): runs the automatic search for one genome and returns one
# lookup row. `status` records the outcome:
#   matched          — a GenBank record with a matching strain field was found
#   no_strain_code   — the MycoCosm name carries no usable strain code
#   too_many_hits    — GenBank probably ignored the strain term (see settings)
#   no_verified_hit  — hits were found, but none had a matching strain field
search_one <- function(g) {
  row <- tibble::tibble(taxon_name = g$taxon_name, genus = g$genus,
                        epithet = g$epithet, strain = g$strain,
                        query = NA_character_, query_translation = NA_character_,
                        n_hits = NA_character_, its_accession = NA_character_,
                        its_length = NA_character_, status = NA_character_,
                        source = "auto", retrieved_on = as.character(Sys.Date()),
                        notes = NA_character_)

  if (is.na(g$strain) || nchar(norm_strain(g$strain)) < STRAIN_MIN_CHARS) {
    row$status <- "no_strain_code"
    return(row)
  }

  # Organism term: the binomial when there is a plain epithet, otherwise the
  # genus alone (for "sp.", "aff." and "cf." genomes, whose GenBank organism
  # name is unpredictable). The strain check in step 3 guards against the
  # broader search.
  plain_epithet <- !is.na(g$epithet) && grepl("^[a-z-]+$", g$epithet)
  organism <- if (plain_epithet) paste(g$genus, g$epithet) else g$genus
  q <- paste0('"', organism, '"[Organism] AND "', g$strain, '"[All Fields] ',
              'AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])')
  row$query <- q

  res <- rentrez::entrez_search(db = "nuccore", term = q, retmax = MAX_HITS_PER_QUERY)
  Sys.sleep(0.15)  # stay under NCBI's request-rate limit
  row$query_translation <- res$QueryTranslation
  row$n_hits <- as.character(res$count)

  if (res$count == 0) {
    row$status <- "no_verified_hit"
    return(row)
  }
  if (res$count > MAX_HITS_PER_QUERY) {
    row$status <- "too_many_hits"
    return(row)
  }

  # Fetch summaries and read each record's length and strain-type fields.
  s <- rentrez::entrez_summary(db = "nuccore", id = res$ids, always_return_list = TRUE)
  Sys.sleep(0.15)
  cand <- purrr::map_dfr(res$ids, function(id) {
    rec <- s[[id]]
    if (is.null(rec)) stop("No GenBank summary returned for UID ", id)
    st <- strsplit(if (!is.null(rec$subtype)) rec$subtype else "", "|", fixed = TRUE)[[1]]
    sn <- strsplit(if (!is.null(rec$subname)) rec$subname else "", "|", fixed = TRUE)[[1]]
    keep <- st %in% c("strain", "isolate", "specimen_voucher", "culture_collection")
    keep <- keep & seq_along(st) <= length(sn)
    tibble::tibble(accession = rec$accessionversion,
                   length    = as.integer(rec$slen),
                   title     = rec$title,
                   ok        = strain_matches(sn[keep], g$strain))
  })

  # Rank the strain-verified records. Records whose title says they include
  # ITS2 come first, because only ITS2 is used for SH assignment and a longer
  # ITS1-only record would otherwise win (this happened for Cenococcum
  # geophilum 1.58). Within that, the longest record wins; the accession
  # breaks any remaining tie so the choice is the same on every run.
  cand <- cand |>
    dplyr::filter(ok, length <= MAX_ITS_LENGTH) |>
    dplyr::mutate(has_its2_title = grepl(
      "internal transcribed spacer 2|ITS2|internal transcribed spacer region",
      title, ignore.case = TRUE)) |>
    dplyr::arrange(dplyr::desc(has_its2_title), dplyr::desc(length), accession)

  if (nrow(cand) == 0L) {
    row$status <- "no_verified_hit"
    return(row)
  }
  row$its_accession <- cand$accession[1]
  row$its_length    <- as.character(cand$length[1])
  row$status        <- "matched"
  if (nrow(cand) > 1L) {
    row$notes <- paste0(nrow(cand), " strain-verified records; longest kept. Others: ",
                        paste(cand$accession[-1], collapse = "; "))
  }
  row
}

to_search <- genomes |> dplyr::filter(!taxon_name %in% lookup$taxon_name)

if (nrow(to_search) > 0L) {
  if (Sys.getenv("ENTREZ_KEY") == "") {
    stop(nrow(to_search), " genome(s) are not yet in ", lookup_path,
         " and need a GenBank search, but no NCBI API key was found.\n",
         "  Set ENTREZ_KEY in your .Renviron file.")
  }
  rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
  message("Step 2: searching GenBank for ", nrow(to_search), " new genome(s)...")
  new_rows <- purrr::map_dfr(seq_len(nrow(to_search)),
                             function(i) search_one(to_search[i, ]))
  lookup <- dplyr::bind_rows(lookup, new_rows) |>
    dplyr::arrange(taxon_name)
  readr::write_csv(lookup, lookup_path, na = "")
} else {
  message("Step 2: every genome already in the lookup table; no GenBank search needed.")
}

# Work only with the genomes in the current MycoCosm list (the lookup can hold
# rows for genomes that have since dropped off the list; they are kept but not
# counted).
lookup_now <- lookup |> dplyr::filter(taxon_name %in% genomes$taxon_name)

# =============================================================================
# Step 3: Download (or reuse) the ITS sequences
# =============================================================================
# Sequences already in the cache file are reused; only missing accessions are
# downloaded. The download is checked: every requested accession must arrive.

wanted <- unique(stats::na.omit(lookup_now$its_accession[lookup_now$status == "matched"]))
cache  <- read_fasta_tbl(seq_cache_path) |>
  dplyr::mutate(accession = sub("\\s.*$", "", header))
missing_acc <- setdiff(wanted, cache$accession)

if (length(missing_acc) > 0L) {
  if (Sys.getenv("ENTREZ_KEY") == "") stop("ENTREZ_KEY needed to download ITS sequences.")
  rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
  fetched <- character(0)
  for (b in split(missing_acc, ceiling(seq_along(missing_acc) / 100))) {
    fetched <- c(fetched, rentrez::entrez_fetch(db = "nuccore", id = b,
                                                rettype = "fasta", retmode = "text"))
    Sys.sleep(0.15)
  }
  tmp <- tempfile(fileext = ".fasta")
  readr::write_lines(fetched, tmp)
  new_seqs <- read_fasta_tbl(tmp) |>
    dplyr::mutate(accession = sub("\\s.*$", "", header))
  still_missing <- setdiff(missing_acc, new_seqs$accession)
  if (length(still_missing) > 0L) {
    stop("GenBank did not return these ITS accessions: ",
         paste(still_missing, collapse = ", "), ". Re-run; nothing was saved.")
  }
  cache <- dplyr::bind_rows(cache, new_seqs)
  write_fasta_tbl(cache, seq_cache_path)
  message("Step 3: downloaded ", length(missing_acc), " ITS sequence(s).")
}

# =============================================================================
# Steps 4-5: Extract ITS2 (ITSx), assign to UNITE SHs (vsearch, 98.5%), and
#            decide which SH(s) each genome is credited to
# =============================================================================
# assign_its2_to_sh() (00_setup.R) applies the rule used for our own GenBank
# data in 03_genbank.R: ITS2 only, fragments of at least ITS2_MIN_LENGTH bp,
# vsearch at SH_IDENTITY against the pinned UNITE FASTA, every reference tied
# for the best match kept. Ties: each tied hit is reduced to a "tie identity"
# (the species name for a named species, the SH code for an un-named "_sp"
# SH). If all tied hits share one identity, the genome is credited to every SH
# code among them (a named species can be split across several SHs);
# otherwise it is credited to none ("ambiguous_tie").
#
# These steps re-run every time (seconds for ~150 sequences), so a changed
# lookup table can never leave stale results behind. Each sequence is
# labelled with a short key ("g1", "g2", ...) so results join back to the
# genome without parsing GenBank headers.

matched <- lookup_now |>
  dplyr::filter(status == "matched", !is.na(its_accession)) |>
  dplyr::mutate(key = paste0("g", dplyr::row_number())) |>
  dplyr::left_join(dplyr::select(cache, accession, seq),
                   by = c("its_accession" = "accession"))

per_genome <- if (nrow(matched) > 0L)
  assign_its2_to_sh(dplyr::select(matched, key, seq),
                    file.path(paths$temp_dir, "darwinian"),
                    identity = SH_IDENTITY, min_len = ITS2_MIN_LENGTH) else
  tibble::tibble(key = character(0), outcome = character(0), sh_codes = character(0),
                 unite_names = character(0), best_identity = numeric(0))

assignments <- lookup_now |>
  dplyr::left_join(dplyr::select(matched, taxon_name, key), by = "taxon_name") |>
  dplyr::left_join(per_genome, by = "key") |>
  dplyr::mutate(
    outcome = dplyr::case_when(
      status != "matched"         ~ paste0("no_its_accession (", status, ")"),
      outcome == "no_its2"        ~ "no_its2_extracted",
      outcome == "no_unite_hit"   ~ "no_unite_hit_at_98.5",
      TRUE                        ~ outcome
    ),
    credited_sh = dplyr::if_else(outcome == "sh_assigned", sh_codes, NA_character_)
  ) |>
  dplyr::select(taxon_name, genus, epithet, strain, its_accession, source,
                outcome, credited_sh, unite_names, best_identity, notes)

readr::write_csv(assignments, out_assign, na = "")

# Genomes not assessed: no accession found, or an accession that did not lead
# to an SH. Written with the reason; reported as not assessed.
review <- assignments |>
  dplyr::filter(outcome != "sh_assigned") |>
  dplyr::left_join(dplyr::select(lookup_now, taxon_name, query, n_hits), by = "taxon_name") |>
  dplyr::arrange(outcome, taxon_name)
readr::write_csv(review, out_review, na = "")

# =============================================================================
# Step 6: Coverage of Canadian SHs, named and un-named
# =============================================================================
# The Canadian SH universe is every SH code in `emf`. An SH counts as "named"
# using the same rule as 15_darwinian.R: its UNITE species label is present and
# does not end in a placeholder (sp, spp, cf, aff).

credited <- assignments |>
  dplyr::filter(!is.na(credited_sh)) |>
  tidyr::separate_rows(credited_sh, sep = "; ") |>
  dplyr::distinct(sh_code = credited_sh, taxon_name)

canada_sh <- emf |>
  dplyr::filter(!is.na(sh_code)) |>
  dplyr::distinct(sh_code, genus, species) |>
  dplyr::mutate(
    sp_space = tolower(gsub("_", " ", species)),
    named = !is.na(species) &
      !grepl("\\b(sp\\.?|spp\\.?|cf\\.?|aff\\.?)$", sp_space)
  ) |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(genus = dplyr::first(genus), species = dplyr::first(species),
                   named = any(named), .groups = "drop")

genome_genera <- unique(genomes$genus)

coverage <- canada_sh |>
  dplyr::left_join(
    credited |>
      dplyr::group_by(sh_code) |>
      dplyr::summarise(n_genomes = dplyr::n(),
                       genomes = paste(sort(taxon_name), collapse = "; "),
                       .groups = "drop"),
    by = "sh_code"
  ) |>
  dplyr::mutate(
    n_genomes        = dplyr::coalesce(n_genomes, 0L),
    has_genome       = n_genomes > 0L,
    # NA when the SH has no genus label (genus = NA, see 04): its genus, and
    # so whether that genus has a genome, is unknown.
    genus_has_genome = dplyr::if_else(is.na(genus), NA, genus %in% genome_genera)
  ) |>
  dplyr::arrange(sh_code)

readr::write_csv(coverage, out_cover, na = "")

# ---- Summary table ------------------------------------------------------------
pct <- function(a, b) if (b == 0) NA_real_ else round(100 * a / b, 1)

n_sh        <- nrow(coverage)
n_named     <- sum(coverage$named)
n_unnamed   <- sum(!coverage$named)
g_sh        <- sum(coverage$has_genome)
g_named     <- sum(coverage$has_genome & coverage$named)
g_unnamed   <- sum(coverage$has_genome & !coverage$named)
unnamed_no_genus_genome <- sum(!coverage$named & coverage$genus_has_genome %in% FALSE)

# Named species (UNITE labels) with a genome via the SH route, for comparison
# with the name-based count in 15_darwinian.R.
named_sp_with_genome <- coverage |>
  dplyr::filter(named, has_genome) |>
  dplyr::distinct(species) |>
  nrow()

outcome_counts <- table(factor(assignments$outcome))

summary_tbl <- tibble::tibble(
  metric = c(
    "MycoCosm genomes in EcM genera present in Canada",
    "Genomes with a same-strain ITS accession",
    "Genomes assigned to a UNITE SH (98.5% ITS2)",
    "Genomes not assigned to an SH (see genomes_not_assessed file)",
    "Canadian SHs (all)",
    "Canadian SHs with >= 1 genome",
    "% of Canadian SHs with >= 1 genome",
    "Named Canadian SHs",
    "Named SHs with >= 1 genome",
    "% of named SHs with >= 1 genome",
    "Un-named Canadian SHs",
    "Un-named SHs with >= 1 genome",
    "% of un-named SHs with >= 1 genome",
    "% of un-named SHs WITHOUT a genome",
    "Un-named SHs in a genus with no MycoCosm genome (lower bound on absence)",
    "Named species (UNITE label) with >= 1 genome via SH match"
  ),
  value = c(
    nrow(genomes),
    sum(lookup_now$status == "matched"),
    sum(assignments$outcome == "sh_assigned"),
    sum(assignments$outcome != "sh_assigned"),
    n_sh, g_sh, pct(g_sh, n_sh),
    n_named, g_named, pct(g_named, n_named),
    n_unnamed, g_unnamed, pct(g_unnamed, n_unnamed),
    pct(n_unnamed - g_unnamed, n_unnamed),
    unnamed_no_genus_genome,
    named_sp_with_genome
  )
)

readr::write_csv(summary_tbl, out_summary)

message("Step 6: ", g_sh, " of ", n_sh, " Canadian SHs have >= 1 genome (",
        g_named, " named, ", g_unnamed, " un-named). ",
        nrow(review), " genome(s) not assessed.")
print(outcome_counts)
