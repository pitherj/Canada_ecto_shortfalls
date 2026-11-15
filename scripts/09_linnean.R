# =============================================================================
# 09_linnean.R  —  Linnean shortfall: taxonomic richness; Table S4
# =============================================================================
# PURPOSE
#   How many EcM fungal taxa are recorded in Canada, at each taxonomic level,
#   and how does this compare with the global inventory of known EcM genera?
#   Species-level assignment rates are reported per source (the proportion of
#   SH codes with a species-level UNITE name, for GlobalFungi and GenBank).
#   Also builds Table S4: EcM genera documented in Canada by GBIF specimens
#   but absent from our dataset, each checked by sequence matching (Step 6b).
#
#   Workflow: 1. richness counts (SH codes, genera, species, lineages);
#   2. GlobalFungi singleton SH codes; 3. observed genera vs the global EcM
#   genus inventory (FungalTraits); 4. species-level assignment rates;
#   5-6a. GBIF specimen records; 6b. Table S4 by sequence matching.
#
# INPUTS
#   emf (04); paths$fungaltraits; paths$van_galen_per_sample
#   data_raw/gbif/<key>.zip  GBIF occurrence download (preserved and living
#                            specimens, Canada, Kingdom Fungi)
#   NCBI GenBank via Entrez (only for Table S4 candidate genera not yet
#   downloaded; ENTREZ_KEY in .Renviron); ITSx and vsearch on PATH
#
# OUTPUTS
#   data_raw/genbank_genus_queries/  GenBank ITS records of the Table S4
#       candidate genera (query log, records with location fields,
#       sequences); extended automatically for new genera
#   data_derived/checkpoints/gbif_ecm_canada_raw.csv  parsed GBIF records
#   data_derived/linnean/
#     linnean_summary.csv               richness and related metrics
#     linnean_genus_coverage.csv        FungalTraits EcM genera, observed or not
#     linnean_gbif_ecm_canada.csv       GBIF records of genera in our data
#     linnean_gbif_ecm_nosequence_canada.csv  GBIF records of Table S4
#                                       candidate genera
#     linnean_gbif_ecm_nosequence_confirmed_canada.csv  those confirmed to
#                                       lack Canadian ITS sequences (Figure S5)
#     linnean_table_s4.csv              one row per candidate genus (Table S4)
#     linnean_table_s4_sequence_matches.csv  one row per sequence checked
#
# MANUSCRIPT ITEMS IT FEEDS
#   Table 1, Table 2, Table S4, Figure S5 (via 19), in-text Linnean results.
#
# RE-RUN WHEN
#   The EcM dataset (04) changes. Table S4 downloads only genera not already
#   in data_raw/genbank_genus_queries/. Then re-run 19.
# =============================================================================

source(here::here("scripts", "00_setup.R"))
library(rgbif)

gbif_ckpt <- file.path(paths$temp_dir, "gbif_ecm_canada_raw.csv")

# ---- Step 1: Taxonomic richness counts ---------------------------------------

n_sh      <- dplyr::n_distinct(emf$sh_code, na.rm = TRUE)  # na.rm: genus-resolved GenBank rows carry sh_code = NA
n_genus   <- dplyr::n_distinct(emf$genus, na.rm = TRUE)  # na.rm: records admitted by the family rule without a genus label have genus = NA (04)

# Named species only (UNITE species field; entries ending _sp are unresolved)
n_species <- dplyr::n_distinct(
  emf$species[!is.na(emf$species) & !grepl("_sp$", emf$species)]
)
# Genera that the named species belong to (fewer than n_genus, because some
# genera are represented only by SHs without a species name).
n_genus_named <- dplyr::n_distinct(
  emf$genus[!is.na(emf$species) & !grepl("_sp$", emf$species)], na.rm = TRUE
)
n_lineage <- dplyr::n_distinct(
  emf$ectomycorrhiza_lineage[!is.na(emf$ectomycorrhiza_lineage)]
)

# Coordinate-filtered equivalents (records validated within the Canada boundary)
emf_coords <- dplyr::filter(emf, coord_in_canada == TRUE)

n_sh_coords          <- dplyr::n_distinct(emf_coords$sh_code, na.rm = TRUE)
n_named_sh_coords    <- emf_coords |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(sh_code) |>
  dplyr::n_distinct()
n_genus_coords       <- dplyr::n_distinct(emf_coords$genus, na.rm = TRUE)
n_species_coords     <- dplyr::n_distinct(
  emf_coords$species[!is.na(emf_coords$species) & !grepl("_sp$", emf_coords$species)]
)

# ---- Step 2: GlobalFungi singleton SH codes ----------------------------------

gf_all <- dplyr::filter(emf, source == "GlobalFungi")

sh_abundance <- gf_all |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(total_abundance = sum(abundance, na.rm = TRUE),
                   .groups = "drop")

singleton_sh     <- dplyr::filter(sh_abundance, total_abundance == 1L)$sh_code
n_sh_singletons  <- length(singleton_sh)
n_sh_gf_total    <- dplyr::n_distinct(gf_all$sh_code, na.rm = TRUE)
n_sh_nonsing     <- n_sh_gf_total - n_sh_singletons
pct_singletons   <- round(100 * n_sh_singletons / n_sh_gf_total, 1)

# ---- Step 3: Compare observed genera to global EcM inventory ----------------

ft <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
  dplyr::rename_with(tolower)

global_ecm_genera <- ft |>
  dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
  dplyr::mutate(genus_lower = tolower(trimws(genus))) |>
  dplyr::distinct(genus_lower, .keep_all = TRUE)

n_global_ecm_genera <- nrow(global_ecm_genera)

our_genera_lower <- tolower(trimws(unique(emf$genus[!is.na(emf$genus)])))

genus_coverage <- global_ecm_genera |>
  dplyr::mutate(observed_in_canada = genus_lower %in% our_genera_lower) |>
  dplyr::select(genus, genus_lower, observed_in_canada,
                ectomycorrhiza_lineage_template,
                ectomycorrhiza_exploration_type_template)

n_observed_genera   <- sum(genus_coverage$observed_in_canada)
pct_observed_genera <- round(100 * n_observed_genera / n_global_ecm_genera, 1)

our_genera_not_in_ft <- our_genera_lower[
  !our_genera_lower %in% global_ecm_genera$genus_lower
]

readr::write_csv(genus_coverage,
                 file.path(paths$out_linnean, "linnean_genus_coverage.csv"))

# ---- Step 4: Species-level assignment rates by source -----------------------

# --- 4a. GlobalFungi: proportion of SH codes with species-level UNITE name ---
# A SH code has a species-level name when the 'species' field does NOT end in
# '_sp' (which indicates UNITE could not assign a species epithet).

gf_species_rate <- gf_all |>
  dplyr::group_by(sh_code) |>
  dplyr::summarise(
    has_species = any(!is.na(species) & !grepl("_sp$", species)),
    .groups = "drop"
  )

n_gf_sh_with_sp  <- sum(gf_species_rate$has_species)
n_gf_sh_total    <- nrow(gf_species_rate)
pct_gf_sp        <- round(100 * n_gf_sh_with_sp / n_gf_sh_total, 1)

# --- 4b. GenBank: total EcM sequence records ---------------------------------
# Every retained GenBank record cleared the 98.5% vsearch identity threshold by
# construction (see 03_genbank.R), so no additional identity filter is applied
# here — we simply report the total record count. The `identity` column
# (vsearch alignment identity) is summarised below for reference only.

gb_all <- dplyr::filter(emf, source == "GenBank")
n_gb_total <- nrow(gb_all)

# Distribution of identity values (summary)

# GenBank SH codes with species-level UNITE name (analogous to GlobalFungi metric)
gb_sh_total      <- dplyr::n_distinct(gb_all$sh_code, na.rm = TRUE)
gb_sh_with_sp    <- gb_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  nrow()
pct_gb_sh_sp <- round(100 * gb_sh_with_sp / gb_sh_total, 1)

# --- 4b.5 Overlap of ALL SH codes (named + unnamed) between GF and GenBank ----
# Used to compute the dark fraction for source-exclusive and shared subsets
# (Table S4 in the SI). gf_all and gb_all are already defined above.
# Exclude NA sh_code (genus-resolved GenBank rows) so it is not counted as a
# spurious "GenBank-only" SH in the set operations below.
gf_all_sh <- unique(gf_all$sh_code[!is.na(gf_all$sh_code)])
gb_all_sh <- unique(gb_all$sh_code[!is.na(gb_all$sh_code)])
n_sh_all_shared  <- length(intersect(gf_all_sh, gb_all_sh))
n_sh_all_gf_only <- length(setdiff(gf_all_sh, gb_all_sh))
n_sh_all_gb_only <- length(setdiff(gb_all_sh, gf_all_sh))

# --- 4c. Overlap of named-species SH codes between GlobalFungi and GenBank ----
gf_sp_sh <- gf_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(sh_code) |>
  unique()

gb_sp_sh <- gb_all |>
  dplyr::distinct(sh_code, species) |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  dplyr::pull(sh_code) |>
  unique()

n_sp_sh_shared  <- length(intersect(gf_sp_sh, gb_sp_sh))
n_sp_sh_gf_only <- length(setdiff(gf_sp_sh, gb_sp_sh))
n_sp_sh_gb_only <- length(setdiff(gb_sp_sh, gf_sp_sh))
n_sp_sh_total   <- n_sp_sh_shared + n_sp_sh_gf_only + n_sp_sh_gb_only

# --- 4d. Overlap of unique named species (epithet-level) between GF and GenBank ----
# Unique named species are the distinct UNITE species epithets carried by the
# named-species SH codes above. Multiple SH codes can share an epithet (see
# Table S2 caption), so these counts are lower than and not simply derivable
# from the SH-level counts in 4c.
gf_sp_names <- gf_all |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::pull(species) |>
  unique()

gb_sp_names <- gb_all |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species, ignore.case = TRUE)) |>
  dplyr::pull(species) |>
  unique()

n_sp_shared   <- length(intersect(gf_sp_names, gb_sp_names))
n_sp_gf_only  <- length(setdiff(gf_sp_names, gb_sp_names))
n_sp_gb_only  <- length(setdiff(gb_sp_names, gf_sp_names))

# ---- Step 6: GBIF physical specimen records ----------------------------------
# Query GBIF for preserved/living fungal specimens in Canada, then filter to
# EcM genera using our observed genus list. Checkpointed to avoid re-querying.
#
# Prerequisites: GBIF credentials in environment variables:
#   gbif_user, gbif_pwd, gbif_email
# Set these with usethis::edit_r_environ() or Sys.setenv().

if (file.exists(gbif_ckpt)) {
  # na = "NA": only the literal "NA" written by write_csv() is missing. An empty
  # field stays an empty string, exactly as in the table imported from the ZIP,
  # so a run that reads this checkpoint writes the same files as one that
  # imports the ZIP.
  gbif_raw <- readr::read_csv(
    gbif_ckpt,
    na = "NA",
    show_col_types = FALSE,
    col_types = readr::cols(
      decimalLatitude  = readr::col_double(),
      decimalLongitude = readr::col_double(),
      .default         = readr::col_character()
    )
  )

} else if (length(list.files(here::here("data_raw", "gbif"), pattern = "\\.zip$")) > 0) {

  # A GBIF download ZIP is already provided in data_raw/gbif/. Import it directly
  # instead of submitting a fresh (slow, credential-gated) download.
  # occ_download_import(key, path) reads a local "<key>.zip" from `path`, where
  # the key is the file name without its extension.
  gbif_zip_dir  <- here::here("data_raw", "gbif")
  gbif_zip_file <- list.files(gbif_zip_dir, pattern = "\\.zip$", full.names = TRUE)[1]
  gbif_key      <- sub("\\.zip$", "", basename(gbif_zip_file))
  gbif_raw <- rgbif::occ_download_import(key = gbif_key, path = gbif_zip_dir) |>
    as.data.frame()
  readr::write_csv(gbif_raw, gbif_ckpt)

} else {

  gbif_user  <- Sys.getenv("GBIF_USER")
  gbif_pwd   <- Sys.getenv("GBIF_PWD")
  gbif_email <- Sys.getenv("GBIF_EMAIL")

  if (any(c(gbif_user, gbif_pwd, gbif_email) == "")) {
    warning(
      "GBIF credentials not found in environment variables ",
      "(GBIF_USER, GBIF_PWD, GBIF_EMAIL).\n",
      "Set these with usethis::edit_r_environ() and restart R.\n",
      "Skipping GBIF step."
    )
    gbif_raw <- NULL
  } else {

    gbif_dl <- tryCatch(
      rgbif::occ_download(
        rgbif::pred("country",               "CA"),
        rgbif::pred("HAS_GEOSPATIAL_ISSUE",  FALSE),
        rgbif::pred("taxonKey",              5),    # Kingdom Fungi
        rgbif::pred_in("BASIS_OF_RECORD",
                       c("PRESERVED_SPECIMEN", "LIVING_SPECIMEN")),
        format = "SIMPLE_CSV",
        user   = gbif_user,
        pwd    = gbif_pwd,
        email  = gbif_email
      ),
      error = function(e) {
        warning("GBIF download submission failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gbif_dl)) {
      rgbif::occ_download_wait(gbif_dl)

      gbif_zip_dir <- here::here("data_raw", "gbif")
      dir.create(gbif_zip_dir, showWarnings = FALSE, recursive = TRUE)
      gbif_zip <- rgbif::occ_download_get(gbif_dl,
                                           path      = gbif_zip_dir,
                                           overwrite = TRUE)
      gbif_raw <- rgbif::occ_download_import(gbif_zip) |>
        as.data.frame()

      readr::write_csv(gbif_raw, gbif_ckpt)

    } else {
      gbif_raw <- NULL
    }
  }
}

# ---- Step 6a: GBIF specimen records of EcM genera ----------------------------
# Species-level GBIF specimen records are split into two sets:
#   - records of genera present in our Canadian sequence dataset ("with
#     sequence data"; paths$gbif_ecm);
#   - records of EcM genera whose name does NOT appear in our dataset
#     (paths$gbif_ecm_nosequence). These genera are the Table S4 candidates.
# A GBIF genus counts as EcM only if FungalTraits lists it as EcM (the genus
# rule of classify_ecm() in 00_setup.R). The two-thirds family rule used for
# our sequence records is NOT applied to GBIF genera: a named specimen genus
# that FungalTraits does not list is usually either a synonym of a listed
# genus (e.g. Midotis for Wynnella) or not EcM at all (e.g. Galera, a
# saprotroph in Hymenogastraceae, a family that passes the rule). For a
# sequence with no genus label, by contrast, the family is the only
# information available.

if (!is.null(gbif_raw) && nrow(gbif_raw) > 0) {

  ref <- ecm_reference()
  our_genera_lower <- tolower(trimws(unique(emf$genus[!is.na(emf$genus)])))

  gbif_sp <- gbif_raw |>
    dplyr::filter(taxonRank == "SPECIES") |>
    dplyr::mutate(genus_lower = tolower(trimws(genus)),
                  ecm_basis   = classify_ecm(genus, family, ref))

  gbif_ecm <- dplyr::filter(gbif_sp, genus_lower %in% our_genera_lower)
  n_gbif_ecm         <- nrow(gbif_ecm)
  n_gbif_ecm_species <- dplyr::n_distinct(gbif_ecm$species)
  n_gbif_ecm_genera  <- dplyr::n_distinct(gbif_ecm$genus)
  readr::write_csv(gbif_ecm, paths$gbif_ecm)

  gbif_ecm_nosequence <- gbif_sp |>
    dplyr::filter(ecm_basis == "genus",
                  !genus_lower %in% our_genera_lower)
  n_gbif_ecm_nosequence         <- nrow(gbif_ecm_nosequence)
  n_gbif_ecm_nosequence_species <- dplyr::n_distinct(gbif_ecm_nosequence$species)
  n_gbif_ecm_nosequence_genera  <- dplyr::n_distinct(gbif_ecm_nosequence$genus)
  readr::write_csv(gbif_ecm_nosequence, paths$gbif_ecm_nosequence)

  candidate_genera <- sort(unique(gbif_ecm_nosequence$genus))

  # ---- Step 6b: Table S4 -- check each candidate genus by sequence ------------
  # A genus name can be missing from our data even though the fungus is there
  # under another name (a synonym, or a genus UNITE files differently). Names
  # cannot settle this automatically, so sequences do. For each candidate genus:
  #   1. Download its fungal ITS records from GenBank, worldwide ("Genus"
  #      [Organism] plus the ITS terms used in 03_genbank.R), with each
  #      record's location field and coordinates. Cached in
  #      paths$table_s4_genbank_dir; a genus not yet in the cache is downloaded
  #      automatically.
  #   2. Assign each sequence's ITS2 to a UNITE SH with the rule used for our
  #      own data (assign_its2_to_sh() in 00_setup.R: ITSx, vsearch at 98.5%).
  #   3. Sort the genus into the first category that applies:
  #      a. "In our data under another name": at least MIN_SEQS_IN_OUR_SHS of
  #         its sequences fall in SHs present in our Canadian dataset. Two, not
  #         one, so that a single misidentified GenBank record cannot decide it.
  #      b. "Canadian ITS sequences in GenBank, not in our dataset": at least one
  #         record passes the GenBank location rule of 04_combine_ecm_dataset.R
  #         (location field begins "Canada", or no field and coordinates
  #         inside Canada).
  #      c. "ITS sequences exist, none from Canada".
  #      d. "No ITS sequences found".
  # Categories b-d are the gaps listed in Table S4. Figure S5 plots c and d as
  # "without sequence data" (paths$gbif_ecm_nosequence_confirmed).

  MIN_SEQS_IN_OUR_SHS <- 2L
  # GenBank records downloaded per genus, at most. Candidate genera are absent
  # from our data and usually have few records; the cap only prevents a very
  # large genus from making the download impractical. Records longer than
  # MAX_ITS_LENGTH are whole genomes or scaffolds, not ITS deposits, and are
  # not downloaded.
  MAX_SEQS_PER_GENUS <- 1000L
  MAX_ITS_LENGTH     <- 5000L

  gdir     <- paths$table_s4_genbank_dir
  log_path <- file.path(gdir, "genus_query_log.csv")
  rec_path <- file.path(gdir, "genus_records.csv")
  seq_path <- file.path(gdir, "genus_sequences.fasta")
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)

  read_or_empty <- function(path, template) {
    if (file.exists(path))
      readr::read_csv(path, show_col_types = FALSE,
                      col_types = readr::cols(.default = "c"))
    else template
  }
  qlog <- read_or_empty(log_path, tibble::tibble(
    genus = character(0), query = character(0), query_translation = character(0),
    n_hits = character(0), n_downloaded = character(0), retrieved_on = character(0)))
  recs <- read_or_empty(rec_path, tibble::tibble(
    genus = character(0), uid = character(0), accession = character(0),
    slen = character(0), country = character(0), lat_lon = character(0)))

  # gb_summaries(): accession, length, location field and coordinates for a
  # set of GenBank UIDs, in batches of 200.
  gb_summaries <- function(uids) {
    batches <- split(uids, ceiling(seq_along(uids) / 200))
    purrr::map_dfr(seq_along(batches), function(i) {
      sm <- fetch_with_retry(function() rentrez::entrez_summary(
              db = "nuccore", id = batches[[i]], always_return_list = TRUE),
            "GenBank summaries", i, length(batches))
      if (is.null(sm)) stop("GenBank summaries could not be downloaded; nothing saved.")
      Sys.sleep(0.15)
      purrr::map_dfr(batches[[i]], function(id) {
        r  <- sm[[id]]
        st <- strsplit(if (!is.null(r$subtype)) r$subtype else "", "|", fixed = TRUE)[[1]]
        sn <- strsplit(if (!is.null(r$subname)) r$subname else "", "|", fixed = TRUE)[[1]]
        field <- function(f) { j <- which(st == f); if (length(j) && j[1] <= length(sn)) sn[j[1]] else NA_character_ }
        tibble::tibble(uid = id, accession = r$accessionversion,
                       slen = as.character(r$slen),
                       country = field("country"), lat_lon = field("lat_lon"))
      })
    })
  }

  to_query <- setdiff(candidate_genera, qlog$genus)
  if (length(to_query) > 0L) {
    if (Sys.getenv("ENTREZ_KEY") == "")
      stop("NCBI API key not found. Set ENTREZ_KEY in your .Renviron file.")
    rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))
    message("Table S4: downloading GenBank ITS records for ", length(to_query), " genera")
    seq_cache <- read_fasta_tbl(seq_path)
    for (g in to_query) {
      q <- paste0('"', g, '"[Organism] AND "Fungi"[Organism] ',
                  'AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])')
      res <- fetch_with_retry(function() rentrez::entrez_search(
               db = "nuccore", term = q, retmax = MAX_SEQS_PER_GENUS),
             paste("GenBank search for", g))
      if (is.null(res)) stop("GenBank search for ", g, " failed; nothing saved for it.")
      meta <- if (length(res$ids)) gb_summaries(res$ids) else
        tibble::tibble(uid = character(0), accession = character(0), slen = character(0),
                       country = character(0), lat_lon = character(0))
      # Drop whole-genome records (longer than MAX_ITS_LENGTH) and whole-genome
      # "master" records (accessions such as KJNK00000000.1), which are index
      # entries that carry no sequence of their own.
      meta <- dplyr::filter(meta, as.integer(slen) <= MAX_ITS_LENGTH,
                            !grepl("^[A-Z]{4,6}0{6,}\\.[0-9]+$", accession))
      if (nrow(meta) > 0L) {
        batches <- split(meta$accession, ceiling(seq_len(nrow(meta)) / 200))
        fasta <- character(0)
        for (i in seq_along(batches)) {
          f <- fetch_with_retry(function() rentrez::entrez_fetch(
                 db = "nuccore", id = batches[[i]], rettype = "fasta", retmode = "text"),
               paste("GenBank sequences for", g), i, length(batches))
          if (is.null(f)) stop("GenBank sequences for ", g, " could not be downloaded; nothing saved for it.")
          fasta <- c(fasta, f)
          Sys.sleep(0.15)
        }
        tmp <- tempfile(fileext = ".fasta")
        readr::write_lines(fasta, tmp)
        new_seqs <- read_fasta_tbl(tmp)
        got <- sub("\\s.*$", "", new_seqs$header)
        if (!all(meta$accession %in% got))
          stop("GenBank did not return every sequence for ", g, "; nothing saved for it.")
        seq_cache <- dplyr::bind_rows(seq_cache,
                                      new_seqs[!got %in% sub("\\s.*$", "", seq_cache$header), ])
      }
      # Save sequences and records first, the log row last: a genus counts as
      # downloaded only once all three are written.
      write_fasta_tbl(seq_cache, seq_path)
      recs <- dplyr::bind_rows(recs, dplyr::mutate(meta, genus = g))
      readr::write_csv(recs, rec_path, na = "")
      qlog <- dplyr::bind_rows(qlog, tibble::tibble(
        genus = g, query = q, query_translation = res$QueryTranslation,
        n_hits = as.character(res$count), n_downloaded = as.character(nrow(meta)),
        retrieved_on = as.character(Sys.Date())))
      readr::write_csv(qlog, log_path, na = "")
      Sys.sleep(0.3)
    }
  }

  # Sequences of the candidate genera, labelled with short keys.
  seq_tbl <- read_fasta_tbl(seq_path) |>
    dplyr::transmute(accession = sub("\\s.*$", "", header), seq)
  cand_recs <- recs |>
    dplyr::filter(genus %in% candidate_genera) |>
    dplyr::inner_join(seq_tbl, by = "accession") |>
    dplyr::mutate(key = paste0("s", dplyr::row_number()))

  sh_result <- if (nrow(cand_recs) > 0L)
    assign_its2_to_sh(dplyr::select(cand_recs, key, seq),
                      file.path(paths$temp_dir, "table_s4")) else
    tibble::tibble(key = character(0), outcome = character(0),
                   sh_codes = character(0), best_identity = numeric(0))

  # SHs in our Canadian dataset, with the genus label each carries there.
  our_sh <- emf |>
    dplyr::filter(!is.na(sh_code)) |>
    dplyr::group_by(sh_code) |>
    dplyr::summarise(our_label = dplyr::first(dplyr::coalesce(genus, genus_label_unite)),
                     .groups = "drop")

  rec_ll <- dplyr::bind_rows(lapply(cand_recs$lat_lon, parse_gb_latlon))
  seq_matches <- cand_recs |>
    dplyr::select(-seq) |>
    dplyr::left_join(sh_result, by = "key") |>
    dplyr::mutate(
      lat = rec_ll$lat, lon = rec_ll$lon,
      canadian = grepl("^Canada", country) |
                 (is.na(country) & coords_in_canada(lat, lon) %in% TRUE),
      our_shs = vapply(strsplit(dplyr::coalesce(sh_codes, ""), "; "),
                       function(x) paste(intersect(x, our_sh$sh_code), collapse = "; "),
                       character(1)),
      in_our_shs = nzchar(our_shs),
      our_labels = vapply(strsplit(our_shs, "; "), function(x)
        paste(sort(unique(our_sh$our_label[our_sh$sh_code %in% x])), collapse = "; "),
        character(1))) |>
    dplyr::select(genus, accession, country, lat_lon, canadian, outcome,
                  sh_codes, best_identity, in_our_shs, our_shs, our_labels)
  readr::write_csv(seq_matches, paths$table_s4_sequences, na = "")

  category_levels <- c("In our data under another name",
                       "Canadian ITS sequences in GenBank, not in our dataset",
                       "ITS sequences exist, none from Canada",
                       "No ITS sequences found")

  per_genus <- seq_matches |>
    dplyr::group_by(genus) |>
    dplyr::summarise(
      n_sequences    = dplyr::n(),
      n_sh_assigned  = sum(outcome == "sh_assigned", na.rm = TRUE),
      n_in_our_shs   = sum(in_our_shs),
      n_canadian     = sum(canadian),
      our_labels     = paste(sort(unique(unlist(strsplit(our_labels[nzchar(our_labels)], "; ")))),
                             collapse = "; "),
      .groups = "drop")

  table_s4 <- gbif_ecm_nosequence |>
    dplyr::group_by(genus, family) |>
    dplyr::summarise(
      gbif_records = dplyr::n(),
      gbif_species = dplyr::n_distinct(species[!is.na(species) & species != ""]),
      ecm_basis    = dplyr::first(ecm_basis),
      .groups = "drop") |>
    dplyr::left_join(dplyr::select(qlog, genus, n_genbank_hits = n_hits,
                                   n_downloaded, retrieved_on), by = "genus") |>
    dplyr::left_join(per_genus, by = "genus") |>
    dplyr::mutate(
      dplyr::across(c(n_sequences, n_sh_assigned, n_in_our_shs, n_canadian),
                    ~ dplyr::coalesce(.x, 0L)),
      n_genbank_hits = as.integer(n_genbank_hits),
      category = dplyr::case_when(
        n_in_our_shs >= MIN_SEQS_IN_OUR_SHS ~ category_levels[1],
        n_canadian > 0                      ~ category_levels[2],
        n_sequences > 0                     ~ category_levels[3],
        TRUE                                ~ category_levels[4]),
      category = factor(category, levels = category_levels)) |>
    dplyr::arrange(category, dplyr::desc(gbif_records))
  readr::write_csv(table_s4, paths$table_s4, na = "")

  # GBIF records of genera confirmed to lack Canadian ITS sequences (c, d):
  # plotted as "without sequence data" in Figure S5 (19_sampling_maps.R).
  gbif_ecm_nosequence |>
    dplyr::filter(genus %in% table_s4$genus[table_s4$category %in% category_levels[3:4]]) |>
    readr::write_csv(paths$gbif_ecm_nosequence_confirmed)

  message("Table S4: ", length(candidate_genera), " candidate genera")
  print(dplyr::count(table_s4, category, name = "n_genera"))

} else {
  n_gbif_ecm                    <- NA_integer_
  n_gbif_ecm_species            <- NA_integer_
  n_gbif_ecm_genera             <- NA_integer_
  n_gbif_ecm_nosequence         <- NA_integer_
  n_gbif_ecm_nosequence_species <- NA_integer_
  n_gbif_ecm_nosequence_genera  <- NA_integer_
}

# ---- Step 7: Save summary table ----------------------------------------------

linnean_summary <- tibble::tibble(
  metric = c(
    "Unique UNITE v10 SH codes (combined dataset, all records)",
    "Unique UNITE v10 SH codes (records with coordinates only)",
    "Named-species SH codes (records with coordinates only)",
    "Unique EcM genera (combined dataset)",
    "Unique EcM genera (records with coordinates only)",
    "Unique named species (UNITE taxonomy, excl. _sp; combined dataset)",
    "Unique EcM genera of the named species (combined dataset)",
    "Unique named species (records with coordinates only)",
    "All SH codes: shared between GlobalFungi and GenBank",
    "All SH codes: GlobalFungi only (not in GenBank)",
    "All SH codes: GenBank only (not in GlobalFungi)",
    "Named-species SH codes: total unique across GF + GenBank (regardless of coords)",
    "Named-species SH codes: shared between GlobalFungi and GenBank",
    "Named-species SH codes: GlobalFungi only",
    "Named-species SH codes: GenBank only",
    "Unique named species: shared between GlobalFungi and GenBank",
    "Unique named species: GlobalFungi only",
    "Unique named species: GenBank only",
    "Unique EcM lineages observed in Canada (lineage definitions from FungalTraits)",
    "Known EcM genera globally (FungalTraits)",
    "FungalTraits EcM genera observed in Canada",
    "% of FungalTraits EcM genera observed in Canada",
    "Canadian EcM genera absent from FungalTraits EcM list",
    "GlobalFungi: unique SH codes (Canadian dataset)",
    "GlobalFungi: singleton SH codes (total abundance = 1 across Canadian samples)",
    "GlobalFungi: singleton SH codes (% of GF Canadian SHs)",
    "GlobalFungi: non-singleton SH codes",
    "GlobalFungi: SH codes with species-level UNITE name",
    "GlobalFungi: SH codes with species-level name (% of GF Canadian SHs)",
    "GenBank: total EcM sequence records",
    "GenBank: unique SH codes (Canadian dataset)",
    "GenBank: SH codes with species-level UNITE name",
    "GenBank: SH codes with species-level name (% of GenBank Canadian SHs)",
    "GBIF: physical EcM specimen records — genera WITH sequence data (Canada)",
    "GBIF: EcM species represented (genera with sequence data)",
    "GBIF: EcM genera represented (genera with sequence data)",
    "GBIF: physical EcM specimen records — genera WITHOUT sequence data (Canada)",
    "GBIF: EcM species represented (genera without sequence data)",
    "GBIF: EcM genera represented (genera without sequence data)"
  ),
  value = c(
    n_sh, n_sh_coords, n_named_sh_coords, n_genus, n_genus_coords, n_species, n_genus_named, n_species_coords,
    n_sh_all_shared, n_sh_all_gf_only, n_sh_all_gb_only,
    n_sp_sh_total, n_sp_sh_shared, n_sp_sh_gf_only, n_sp_sh_gb_only,
    n_sp_shared, n_sp_gf_only, n_sp_gb_only,
    n_lineage,
    n_global_ecm_genera, n_observed_genera, pct_observed_genera,
    length(our_genera_not_in_ft),
    n_sh_gf_total, n_sh_singletons, pct_singletons, n_sh_nonsing,
    n_gf_sh_with_sp, pct_gf_sp,
    n_gb_total,
    gb_sh_total, gb_sh_with_sp, pct_gb_sh_sp,
    n_gbif_ecm, n_gbif_ecm_species, n_gbif_ecm_genera,
    n_gbif_ecm_nosequence, n_gbif_ecm_nosequence_species, n_gbif_ecm_nosequence_genera
  )
)

readr::write_csv(linnean_summary,
                 file.path(paths$out_linnean, "linnean_summary.csv"))
