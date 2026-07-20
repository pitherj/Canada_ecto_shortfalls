# =============================================================================
# Raunkiæran Shortfall: Functional Trait Coverage
# =============================================================================
# How well are the functional traits of Canadian EcM fungal genera documented?
#
# FungalTraits v1.2 (Põlme et al. 2020) is the trait source. It is genus-level.
#
# Trait set
# ---------
# Only SIX FungalTraits columns carry trait information that is meaningful for
# EcM fungi. We assess coverage for exactly these six:
#
#   Morphological / structural traits (derivable from taxonomic description):
#     fruitbody_type_template
#     hymenium_type_template
#     ectomycorrhiza_exploration_type_template
#
#   Ecological / interaction traits:
#     secondary_lifestyle
#     endophytic_interaction_capability_template
#     specific_hosts
#
# Every other FungalTraits column is excluded, for one of two reasons:
#   (a) inapplicable to EcM fungi by lifestyle (decay/saprotroph, lichen
#       photobiont, plant-pathogen, animal-biotroph, aquatic-habitat traits);
#   (b) uninformative: growth_form_template is "filamentous_mycelium" for all
#       327 EcM genera, so it carries no information about the shortfall.
#
# Analyses:
#   1. Per-trait coverage: how many of our Canadian EcM genera have a
#      documented (non-missing) value for each of the six traits.
#   2. Per-genus count of documented traits (of six).
#   3. Value distributions for each trait.
#   4. Genera with a documented host-specificity entry.
#   5. Genera (if any) absent from FungalTraits entirely.
#
# Outputs (all in data_derived/raunkiaeran/):
#   raunkiaeran_trait_coverage.csv      - per-trait coverage counts and %
#   raunkiaeran_genus_summary.csv       - per-genus count of documented traits
#   raunkiaeran_trait_distributions.csv - value tallies for each trait
#   raunkiaeran_specific_hosts.csv      - genera with documented host specificity
#   raunkiaeran_absent_genera.csv       - genera absent from FungalTraits entirely
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Analysis options -------------------------------------------------------
# Should the literal string "unknown" be treated as a MISSING (undocumented)
# value? For a knowledge-coverage metric, an explicit "unknown" means the trait
# state is not known, so the default is TRUE. Set to FALSE to count "unknown"
# as a documented value.
TREAT_UNKNOWN_AS_MISSING <- TRUE

# ---- 1. Load FungalTraits (full CSV) ----------------------------------------

ft_raw <- readr::read_csv(paths$fungaltraits, show_col_types = FALSE) |>
  dplyr::rename_with(tolower) |>
  # Standardize the genus column name and case for joining
  dplyr::mutate(genus_lower = tolower(trimws(genus)))

# ---- 2. Define the six EcM-relevant trait columns ---------------------------
# Column names are lower-cased above, so we list them lower-case here. Each is
# tagged with a human-readable label and a class used to group the output
# table (morphological/structural vs ecological/interaction).

trait_meta <- tibble::tribble(
  ~trait,                                        ~trait_label,                 ~trait_class,
  "fruitbody_type_template",                     "Fruitbody type",             "Morphological/structural",
  "hymenium_type_template",                      "Hymenium type",              "Morphological/structural",
  "ectomycorrhiza_exploration_type_template",    "Mycorrhiza exploration type","Morphological/structural",
  "secondary_lifestyle",                         "Secondary lifestyle",        "Ecological/interaction",
  "endophytic_interaction_capability_template",  "Endophytic interaction capability", "Ecological/interaction",
  "specific_hosts",                              "Specific hosts",             "Ecological/interaction"
)

trait_cols <- trait_meta$trait

# Fail loudly if FungalTraits ever drops/renames one of these columns.
missing_cols <- setdiff(trait_cols, names(ft_raw))
if (length(missing_cols) > 0) {
  stop("FungalTraits is missing expected trait column(s): ",
       paste(missing_cols, collapse = ", "))
}

# ---- 3. Subset FungalTraits to EcM genera -----------------------------------
# FungalTraits has one row per genus; deduplicate on the join key to be safe.

ft_ecm <- ft_raw |>
  dplyr::filter(primary_lifestyle == "ectomycorrhizal") |>
  dplyr::distinct(genus_lower, .keep_all = TRUE)

# ---- 4. Our Canadian EcM genera ---------------------------------------------

our_genera <- emf |>
  dplyr::distinct(genus) |>
  dplyr::mutate(genus_lower = tolower(trimws(genus)))

n_our_genera <- nrow(our_genera)

# Left-join our genera to FungalTraits (keep all our genera). Drop the
# redundant 'genus' column from ft_ecm so our_genera$genus stays authoritative.
our_ft <- dplyr::left_join(our_genera,
                           dplyr::select(ft_ecm, -genus),
                           by = "genus_lower")

# ---- 5. Genera absent from FungalTraits entirely ----------------------------

absent_from_ft <- our_ft |>
  dplyr::filter(is.na(primary_lifestyle)) |>
  dplyr::pull(genus)

absent_df <- tibble::tibble(genus = absent_from_ft)
readr::write_csv(absent_df,
                 file.path(paths$out_raunkiaeran, "raunkiaeran_absent_genera.csv"))

# Restrict the coverage assessment to genera present in FungalTraits.
our_ft_matched <- dplyr::filter(our_ft, !is.na(primary_lifestyle))
n_matched <- nrow(our_ft_matched)

# ---- 6. "Documented" predicate ----------------------------------------------
# A value is documented if it is non-NA, non-empty, and (per the option above)
# not the literal string "unknown".
is_documented <- function(x) {
  s <- trimws(as.character(x))
  ok <- !is.na(x) & s != ""
  if (TREAT_UNKNOWN_AS_MISSING) ok <- ok & tolower(s) != "unknown"
  ok
}

# ---- 7. Per-trait coverage --------------------------------------------------

trait_coverage <- purrr::map_dfr(trait_cols, function(col) {
  vals  <- our_ft_matched[[col]]
  n_doc <- sum(is_documented(vals))
  tibble::tibble(
    trait          = col,
    n_documented   = n_doc,
    n_total        = n_matched,
    pct_documented = round(100 * n_doc / n_matched, 1)
  )
}) |>
  # Attach labels/classes, then order: structural first, by descending coverage.
  dplyr::left_join(trait_meta, by = "trait") |>
  dplyr::select(trait, trait_label, trait_class,
                n_documented, n_total, pct_documented) |>
  dplyr::arrange(trait_class, dplyr::desc(pct_documented))

readr::write_csv(trait_coverage,
                 file.path(paths$out_raunkiaeran, "raunkiaeran_trait_coverage.csv"))

# ---- 8. Per-genus trait-count summary ---------------------------------------
# For each matched genus, count how many of the six traits are documented.

genus_trait_counts <- our_ft_matched |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_traits_documented = sum(is_documented(dplyr::c_across(dplyr::all_of(trait_cols))))
  ) |>
  dplyr::ungroup() |>
  dplyr::select(genus, genus_lower, n_traits_documented) |>
  dplyr::arrange(dplyr::desc(n_traits_documented), genus)

readr::write_csv(genus_trait_counts,
                 file.path(paths$out_raunkiaeran, "raunkiaeran_genus_summary.csv"))

# ---- 9. Trait value distributions -------------------------------------------
# Tally genera per value for each trait (including a "(undocumented)" bucket so
# the gap is visible). specific_hosts is free text but tallied the same way.

trait_distributions <- purrr::map_dfr(trait_cols, function(col) {
  our_ft_matched |>
    dplyr::mutate(
      value = dplyr::if_else(is_documented(.data[[col]]),
                             trimws(as.character(.data[[col]])),
                             "(undocumented)")
    ) |>
    dplyr::count(value) |>
    dplyr::mutate(
      trait      = col,
      pct_genera = round(100 * n / n_matched, 1)
    ) |>
    dplyr::rename(n_genera = n) |>
    dplyr::select(trait, value, n_genera, pct_genera) |>
    dplyr::arrange(dplyr::desc(n_genera))
})

readr::write_csv(trait_distributions,
                 file.path(paths$out_raunkiaeran, "raunkiaeran_trait_distributions.csv"))

# ---- 10. Genera with documented host specificity ----------------------------

specific_hosts_df <- our_ft_matched |>
  dplyr::filter(is_documented(specific_hosts)) |>
  dplyr::select(genus, specific_hosts) |>
  dplyr::arrange(genus)

readr::write_csv(specific_hosts_df,
                 file.path(paths$out_raunkiaeran, "raunkiaeran_specific_hosts.csv"))

