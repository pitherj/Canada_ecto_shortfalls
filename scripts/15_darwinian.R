# =============================================================================
# Darwinian Shortfall: Genomic Data Coverage
# =============================================================================
# How many EcM fungal taxa in our Canadian dataset have genome sequences
# available in MycoCosm (the JGI fungal genomics portal)?
#
# MycoCosm data (data_raw/mycocosm_organism_list.csv):
#   The automated export from the MycoCosm portal requires JGI login and
#   the bulk-export link was not functioning at time of data collection.
#   Data were obtained manually:
#     1. Navigate to https://mycocosm.jgi.doe.gov/fungi/fungi.info.html
#        (the JGI "All Fungi" portal page, no login required)
#     2. The page displays a table of all fungal genomes; select all rows
#        and copy the table content
#     3. Paste into a spreadsheet and save as:
#          data_raw/mycocosm/mycocosm_organism_list.csv
#     4. Manually rename the column headers to:
#          id, taxon_name, assembly_length, num_genes, publication
#   Because the source page lists only fungi, no kingdom/phylum filtering
#   is required.
#
# Workflow:
#   1.  Load MycoCosm organism list
#   2.  Note record count (already fungi-only)
#   3.  Extract genus from taxon_name; match against our EcM genus list
#   4.  Report coverage statistics
#
# Outputs:
#   data_derived/darwinian_mycocosm_matches.csv     — all MycoCosm records whose genus matches an EcM genus
#   data_derived/darwinian_genus_summary.csv        — per-genus genome count
#   data_derived/darwinian_species_matches.csv      — EcM species with an exact species-level genome match
#   data_derived/darwinian_summary.csv              — summary statistics (genus- and species-level)
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Step 1: Load MycoCosm organism list ------------------------------------

ts("Step 1: Loading MycoCosm organism list...")

if (!file.exists(paths$mycocosm_list)) {
  readr::write_csv(
    tibble::tibble(genus = character(), n_genomes = integer(),
                   record_ids = character()),
    paths$darwinian_out
  )
  readr::write_csv(
    tibble::tibble(metric = character(), value = numeric()),
    file.path(paths$out_darwinian, "darwinian_summary.csv")
  )
  stop(
    "MycoCosm organism list not found: ", paths$mycocosm_list, "\n",
    "See script header for manual download procedure.\n",
    "Source: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html",
    call. = FALSE
  )
}

# Read — try comma, then tab, then semicolon delimited
mc_raw <- tryCatch(
  readr::read_csv(paths$mycocosm_list, show_col_types = FALSE),
  error = function(e) {
    tryCatch(
      readr::read_tsv(paths$mycocosm_list, show_col_types = FALSE),
      error = function(e2) {
        readr::read_delim(paths$mycocosm_list, delim = ";",
                          show_col_types = FALSE)
      }
    )
  }
)

# Normalise column names to lower case
names(mc_raw) <- tolower(gsub("\\s+", "_", names(mc_raw)))

ts(sprintf("  MycoCosm records: %d", nrow(mc_raw)))
ts(sprintf("  Columns: %s", paste(names(mc_raw), collapse = ", ")))

# ---- Step 2: Record count (file is already fungi-only) ----------------------

ts("Step 2: Data sourced from fungi-only portal page — no filtering required.")
mc_fungi <- mc_raw
ts(sprintf("  MycoCosm fungal records: %d", nrow(mc_fungi)))

# ---- Step 3: Extract genus from organism name if no genus column -----------

# Expected column is taxon_name; fall back to other common name columns
name_col <- intersect(c("taxon_name", "name", "organism_name",
                         "scientific_name", "portal_id"), names(mc_fungi))[1]

if (!"genus" %in% names(mc_fungi)) {
  if (!is.na(name_col)) {
    ts(sprintf("  Parsing genus from column '%s'...", name_col))
    mc_fungi <- mc_fungi |>
      dplyr::mutate(genus = sub("^([A-Za-z]+).*", "\\1",
                                 trimws(.data[[name_col]])))
  } else {
    ts("  Could not extract genus — no suitable name column found.")
    mc_fungi$genus <- NA_character_
  }
}

mc_fungi <- mc_fungi |>
  dplyr::mutate(
    genus_lower = tolower(trimws(genus)),
    # Extract binomial (first two words) from taxon_name for species matching
    species_lower = if (!is.na(name_col)) {
      tolower(trimws(
        sub("^([A-Za-z]+\\s+[a-z]+).*", "\\1", trimws(.data[[name_col]]))
      ))
    } else {
      NA_character_
    }
  )

# ---- Step 4: Match against our EcM species and genus lists -------------------

ts("Step 4: Matching MycoCosm taxa against our EcM species and genus lists...")

# emf$species uses underscores; normalise to spaces for comparison.
# Exclude placeholder epithets (sp, sp., spp, cf, aff) — these are genus-level
# records only and should not be treated as species-level genome matches.
our_species_lower <- tolower(gsub("_", " ", unique(emf$species[!is.na(emf$species)])))
our_species_lower <- our_species_lower[
  !grepl("\\b(sp\\.?|spp\\.?|cf\\.?|aff\\.?)$", our_species_lower)
]
our_genera_lower  <- tolower(trimws(unique(emf$genus)))
# Denominator must match the universe actually searched for matches above:
# our_species_lower already excludes NA and genus-level placeholders
# (sp/spp/cf/aff). Using raw dplyr::n_distinct(emf$species) here previously
# inflated the denominator and made the with/absent percentages not sum to
# 100 (see Prestonian script's analogous, correctly-filtered species count).
n_our_species     <- length(our_species_lower)
n_our_genera      <- dplyr::n_distinct(emf$genus)

all_mc_genera  <- unique(mc_fungi$genus_lower[!is.na(mc_fungi$genus_lower)])
all_mc_species <- unique(mc_fungi$species_lower[!is.na(mc_fungi$species_lower)])
ts(sprintf("  Unique genera  in MycoCosm: %d", length(all_mc_genera)))
ts(sprintf("  Unique species in MycoCosm: %d", length(all_mc_species)))

mc_ecm <- mc_fungi |>
  dplyr::filter(genus_lower %in% our_genera_lower) |>
  dplyr::mutate(
    species_match = species_lower %in% our_species_lower
  )

n_mc_ecm_records  <- nrow(mc_ecm)
n_mc_ecm_genera   <- dplyr::n_distinct(mc_ecm$genus_lower)
n_mc_ecm_species  <- sum(mc_ecm$species_match, na.rm = TRUE)
n_ecm_sp_w_genome <- dplyr::n_distinct(
  mc_ecm$species_lower[mc_ecm$species_match]
)
ts(sprintf("  MycoCosm records matching our EcM genera:   %d", n_mc_ecm_records))
ts(sprintf("  EcM genera  with genome data: %d / %d (%.1f%%)",
           n_mc_ecm_genera, n_our_genera,
           100 * n_mc_ecm_genera / n_our_genera))
ts(sprintf("  EcM species with genome data: %d / %d (%.1f%%)",
           n_ecm_sp_w_genome, n_our_species,
           100 * n_ecm_sp_w_genome / n_our_species))

# ---- Step 5: Per-genus summary -----------------------------------------------

# Columns to retain in output (actual columns from manual export)
keep_cols <- intersect(
  c("id", "taxon_name", "genus", "genus_lower", "species_lower",
    "species_match", "assembly_length", "num_genes", "publication"),
  names(mc_ecm)
)

mc_ecm_out <- dplyr::select(mc_ecm, dplyr::all_of(keep_cols))

genus_genome_summary <- mc_ecm_out |>
  dplyr::group_by(genus_lower) |>
  dplyr::summarise(
    n_genomes = dplyr::n(),
    record_ids = if ("id" %in% names(mc_ecm_out)) {
      paste(unique(id), collapse = "; ")
    } else {
      NA_character_
    },
    .groups = "drop"
  ) |>
  dplyr::rename(genus = genus_lower) |>
  dplyr::arrange(dplyr::desc(n_genomes))

# Full MycoCosm records for EcM-genus matches
readr::write_csv(mc_ecm_out, paths$darwinian_out)
ts(sprintf("  Saved -> %s", basename(paths$darwinian_out)))

# Per-genus summary: how many genomes per EcM genus
readr::write_csv(genus_genome_summary,
                 file.path(paths$out_darwinian, "darwinian_genus_summary.csv"))
ts("  Saved -> darwinian_genus_summary.csv")

# Species-level match list: which of our EcM species have a genome in MycoCosm
species_match_list <- mc_ecm |>
  dplyr::filter(species_match) |>
  dplyr::select(dplyr::any_of(c("species_lower", "taxon_name", "id",
                                  "assembly_length", "num_genes", "publication"))) |>
  dplyr::rename(emf_species = species_lower) |>
  dplyr::arrange(emf_species)

readr::write_csv(species_match_list,
                 file.path(paths$out_darwinian, "darwinian_species_matches.csv"))
ts(sprintf("  Saved -> darwinian_species_matches.csv  (%d records)",
           nrow(species_match_list)))

# Genera absent from MycoCosm
our_genera_not_in_mc <- our_genera_lower[!our_genera_lower %in% all_mc_genera]
pct_absent <- round(100 * length(our_genera_not_in_mc) / n_our_genera, 1)
ts(sprintf("  EcM genera absent from MycoCosm: %d (%.1f%%)",
           length(our_genera_not_in_mc), pct_absent))

# ---- Step 6: Save summary table ----------------------------------------------

our_species_not_in_mc <- our_species_lower[!our_species_lower %in% all_mc_species]
pct_sp_absent <- round(100 * length(our_species_not_in_mc) / n_our_species, 1)

darwinian_summary <- tibble::tibble(
  metric = c(
    "EcM genera in our Canadian dataset",
    "EcM species (UNITE, named) in our Canadian dataset",
    "MycoCosm fungal records (total)",
    "MycoCosm records matching our EcM genera",
    "EcM genera with genome data in MycoCosm",
    "% of our EcM genera with genome data",
    "EcM genera absent from MycoCosm",
    "% of our EcM genera absent from MycoCosm",
    "EcM species with genome data in MycoCosm",
    "% of our EcM species with genome data",
    "EcM species absent from MycoCosm",
    "% of our EcM species absent from MycoCosm"
  ),
  value = c(
    n_our_genera,
    n_our_species,
    nrow(mc_fungi),
    n_mc_ecm_records,
    n_mc_ecm_genera,
    round(100 * n_mc_ecm_genera / n_our_genera, 1),
    length(our_genera_not_in_mc),
    pct_absent,
    n_ecm_sp_w_genome,
    round(100 * n_ecm_sp_w_genome / n_our_species, 1),
    length(our_species_not_in_mc),
    pct_sp_absent
  )
)

readr::write_csv(darwinian_summary,
                 file.path(paths$out_darwinian, "darwinian_summary.csv"))
ts("Saved darwinian_summary.csv")
print(as.data.frame(darwinian_summary))
ts("15_darwinian.R complete.")
