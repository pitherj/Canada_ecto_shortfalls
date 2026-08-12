# =============================================================================
# 21_geographic_possibility.R
#   Sensitivity analysis: how much of the "full potential interaction matrix"
#   reported in the Eltonian section is geographically possible at all?
# =============================================================================
# STATUS OF THIS SCRIPT
#   Written 2026-08-05 as exploratory work in temp/. Promoted into the pipeline
#   on 2026-08-11, when its result was taken up in the revision: the manuscript
#   now states that even under the most permissive geographic filter we could
#   construct, half the matrix remains and the proportion of empty cells barely
#   moves. A number that appears in the manuscript must be reproducible by
#   scripts/run_all.R like every other, hence the promotion.
#
#   The analysis is a SENSITIVITY CHECK, not a replacement denominator. The
#   manuscript keeps the full matrix as its headline (framed as a metaweb --
#   the complete set of conceivable pairings) and cites this script only for
#   the reassurance that the conclusion does not depend on that choice. The
#   reasons not to adopt the filtered denominator are in the KNOWN LIMITATIONS
#   block below and are stated in the supplemental materials.
#
#   It READS pipeline inputs and WRITES ONLY to data_derived/geo_possibility/.
#
# THE QUESTION
#   scripts/18_eltonian.R section A3c reports the Canadian host x named-EcM-
#   fungal-species matrix as
#       n_host_species  x  n_named_fungal_species  =  the full grid of cells
#   and states what percentage of those cells is empty. That framing
#   treats every host x fungus pair as biologically meaningful, which overstates
#   the denominator: many pairs could never occur. Two independent reasons:
#       (1) geographic mismatch  - the host and the fungus never co-occur;
#       (2) incompatibility      - host breadth / phylogenetic constraint.
#   THIS SCRIPT ADDRESSES (1) ONLY. (2) is out of scope.
#
# THE METHOD (and why it is shaped this way)
#   The natural way to filter on geography would be to intersect host ranges
#   with fungal ranges. But fungal ranges in Canada are essentially unknown --
#   that is the Wallacean shortfall the same manuscript documents. So we invert
#   the problem and use only data we actually have:
#
#       For each named EcM fungal species, take its georeferenced Canadian
#       occurrence records. A host species is GEOGRAPHICALLY POSSIBLE for that
#       fungus if at least one of those records falls within the host's
#       modelled (BIEN) range.
#
#   No fungal range model is required. What comes out is, per fungus, the set of
#   host species whose ranges overlap the places the fungus has actually been
#   recorded.
#
# KNOWN LIMITATIONS (these are to be REPORTED, not fixed)
#   * This is a STRICT LOWER BOUND on the possible-host set. A fungus known from
#     one location gets a host set from one location, so the estimate is limited
#     by sampling effort rather than by the fungus's true range. Output 3 is the
#     diagnostic for exactly this.
#   * The bias direction is unhelpful. Shrinking the denominator while the
#     numerator (the observed pairs) stays fixed INCREASES the apparent fill
#     rate, making the shortfall look LESS severe. This is not a conservative
#     correction, and it pulls against the "our estimates are conservative"
#     position taken elsewhere in the Eltonian section.
#   * BIEN modelled ranges carry their own error. A pair that is observed in the
#     data but flagged geographically impossible indicates a range-map (or
#     coordinate) failure, not an error in our interaction records. Output 4
#     quantifies that rate and so bounds the reliability of the whole exercise.
#
# OUTPUTS (all to data_derived/geo_possibility/)
#   output1_possible_pairs_grain_a.csv    Output 1 - the possibility matrix
#   output2_grain_sensitivity.csv         Output 2 - possibility by spatial grain
#   output3_effort_by_species.csv         Output 3 - set size vs sampling effort
#   output3_rarefaction.csv               Output 3 - rarefaction curves
#   output3_effort_limitation.png         Output 3 - figure
#   output4_observed_but_impossible.csv   Output 4 - sanity check vs observed
#   output5_fill_statistics.csv           Output 5 - recomputed Eltonian numbers
#   output6_summary.csv                   Output 6 - headline numbers as a
#                                         metric/value table, read directly by
#                                         the manuscript and SM1
#   run_log.txt                           console summary of this run
#
# HOW TO RUN
#   From the project root:  source(here::here("scripts", "21_geographic_possibility.R"))
#   or as part of scripts/run_all.R. Runtime is a few minutes, dominated by the
#   100 km buffered-distance grain. Must be run AFTER 18_eltonian.R: it reads
#   that script's outputs and stops if its own denominator does not reproduce
#   them exactly.
# =============================================================================


# =============================================================================
# SECTION 0.  Setup, configuration, and the random seed
# =============================================================================
# Sourcing 00_setup.R gives us the canonical `paths` list, the CRS constants,
# and the loaded `emf` occurrence table. It is read-only for our purposes: the
# only thing it writes is empty data_derived/ sub-directories, all of which
# already exist, so sourcing it is a no-op on disk.

source(here::here("scripts", "00_setup.R"))

suppressPackageStartupMessages({
  library(sf)
  library(terra)
})

# ---- Where our outputs go ---------------------------------------------------
out_dir <- paths$out_geo_possibility
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- The random seed --------------------------------------------------------
# Output 3 rarefies occurrence records, which is a random subsampling step.
# The seed is fixed here, printed to the run log, and stored as a column in
# output3_rarefaction.csv so the exact draws can always be reproduced.
SEED <- 20260805L
set.seed(SEED)

# ---- Rarefaction settings ---------------------------------------------------
RAREFY_K    <- c(1L, 2L, 3L, 5L, 10L, 20L)  # records subsampled per species
RAREFY_REPS <- 20L                          # independent draws at each k

# ---- Buffer distances for the coarse-grain sensitivity (metres) -------------
BUFFERS_M <- c(25000, 100000)               # 25 km and 100 km

# ---- A tiny logging helper --------------------------------------------------
# Everything we print to the console is also appended to run_log.txt, so the
# numbers quoted in FINDINGS.md can be traced back to a specific run.
log_file <- file.path(out_dir, "run_log.txt")
if (file.exists(log_file)) file.remove(log_file)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_file, append = TRUE)
}

say("21_geographic_possibility.R  |  run ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("random seed = ", SEED)
say(strrep("=", 78))


# =============================================================================
# SECTION 1.  Load the inputs
# =============================================================================
# Five inputs, all produced by the existing pipeline and all read-only here:
#   1. the host-species list         (denominator rows)
#   2. the emf occurrence table      (denominator columns + the fungal points)
#   3. BIEN modelled host ranges     (polygons, clipped to Canada)
#   4. the 0.5-degree host stack     (raster version of the same polygons)
#   5. ecoregions + ecozone names    (the two coarsest grains)
# Plus the two Eltonian outputs we must reproduce and compare against.

say("\n[1] Loading inputs")

# ---- 1a. Host species: the ROW denominator ----------------------------------
# This is exactly what 18_eltonian.R section A1 does: n_host_species is simply
# the number of rows in this file. No filtering.
host_tbl       <- readr::read_csv(paths$host_species, show_col_types = FALSE)
host_all       <- host_tbl$species
n_host_species <- length(host_all)

# ---- 1b. Named fungal species: the COLUMN denominator -----------------------
# Exactly what 18_eltonian.R section A3c does: every distinct value of
# emf$species that is not missing and does not end in "_sp" (i.e. is resolved to
# a binomial), with underscores converted to spaces. Note this is taken over ALL
# emf records, INCLUDING those with no coordinates -- a fungal species with no
# georeferenced record still occupies a column of the matrix.
named_fungal_species <- emf |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::mutate(fungal_species = trimws(gsub("_", " ", species))) |>
  dplyr::distinct(fungal_species) |>
  dplyr::pull(fungal_species)
n_named_fungal_species <- length(named_fungal_species)

n_matrix_cells <- n_host_species * n_named_fungal_species

# ---- 1c. The observed pairs and their occurrence counts ---------------------
# eltonian_species_occurrence_counts.csv is the long form of the observed
# host x fungus matrix: one row per observed pair, with the number of supporting
# occurrences. Its row count is the number of FILLED cells in 18_eltonian.R.
occ_counts <- readr::read_csv(
  file.path(paths$out_eltonian, "eltonian_species_occurrence_counts.csv"),
  show_col_types = FALSE)

n_cells_filled      <- nrow(occ_counts)
n_pairs_singleton   <- sum(occ_counts$n_occurrences == 1L)

# ---- 1d. Verify the denominators against the submitted numbers --------------
# The whole point of this exercise is a number that is directly comparable to
# what is already in the manuscript. If our denominator does not reproduce the
# submitted one, everything downstream is incomparable, so we STOP rather than
# quietly proceed with a different denominator.
eltonian_summary <- readr::read_csv(
  file.path(paths$out_eltonian, "eltonian_summary.csv"), show_col_types = FALSE)

get_metric <- function(pattern) {
  hit <- eltonian_summary$value[grepl(pattern, eltonian_summary$metric, fixed = TRUE)]
  if (length(hit) != 1L)
    stop("Could not uniquely locate metric '", pattern, "' in eltonian_summary.csv")
  as.numeric(hit)
}

submitted <- list(
  n_hosts   = get_metric("Canadian EcM host plant species (denominator"),
  n_fungi   = get_metric("Named EcM fungal species detected anywhere in the Canadian dataset (full column denominator)"),
  n_cells   = get_metric("Full host x named-species matrix size"),
  n_filled  = get_metric("Filled cells in full host x named-species matrix"),
  pct_empty = get_metric("% of full host x named-species matrix cells empty"),
  n_single  = get_metric("Observed host x named-species pairs supported by exactly 1 occurrence"),
  pct_single = get_metric("% of observed host x named-species pairs supported by exactly 1 occurrence")
)

mismatches <- c(
  if (n_host_species         != submitted$n_hosts)  sprintf("hosts: recomputed %d vs submitted %d",  n_host_species, submitted$n_hosts),
  if (n_named_fungal_species != submitted$n_fungi)  sprintf("fungi: recomputed %d vs submitted %d",  n_named_fungal_species, submitted$n_fungi),
  if (n_matrix_cells         != submitted$n_cells)  sprintf("cells: recomputed %d vs submitted %d",  n_matrix_cells, submitted$n_cells),
  if (n_cells_filled         != submitted$n_filled) sprintf("filled: recomputed %d vs submitted %d", n_cells_filled, submitted$n_filled),
  if (n_pairs_singleton      != submitted$n_single) sprintf("singletons: recomputed %d vs submitted %d", n_pairs_singleton, submitted$n_single)
)
if (length(mismatches)) {
  stop("Denominator does NOT match scripts/18_eltonian.R section A3c:\n  ",
       paste(mismatches, collapse = "\n  "),
       "\nStopping: results would not be comparable to the submitted numbers.")
}

say("    host species (matrix rows)      : ", n_host_species)
say("    named fungal species (columns)  : ", n_named_fungal_species)
say("    full matrix cells               : ", format(n_matrix_cells, big.mark = ","))
say("    observed (filled) cells         : ", n_cells_filled)
say("    denominator matches 18_eltonian.R section A3c: TRUE")

# ---- 1e. BIEN modelled host ranges (polygons) -------------------------------
# Produced by 08_host_rasters.R: one multipolygon per host species, in WGS84,
# already filtered to species whose range overlaps Canada.
ranges <- sf::st_read(paths$bien_ranges, quiet = TRUE)
hosts_with_range <- ranges$species
say("    host species with a BIEN range  : ", length(hosts_with_range),
    " of ", n_host_species,
    "  (", n_host_species - length(hosts_with_range), " have no range and are")
say("                                      structurally impossible for every fungus)")

# ---- 1f. The 0.5-degree per-species binary stack ----------------------------
# One layer per host species, layer names = species names with spaces.
host_stack <- terra::rast(paths$bien_species_stack)

# ---- 1g. Ecoregions and ecozones --------------------------------------------
ecoregions <- sf::st_read(paths$ecoregions_processed, quiet = TRUE)


# =============================================================================
# SECTION 2.  The fungal occurrence points
# =============================================================================
# We keep the records that can contribute evidence:
#   * a named fungal species (same rule as the column denominator),
#   * a coordinate that falls in Canada (emf$coord_in_canada is TRUE exactly
#     when lat/lon are present and inside the Canada boundary).
#
# NOTE ON THE RECORD UNIT. A row of `emf` is one (sample x SH-detection) record.
# That is the same unit 18_eltonian.R counts as an "occurrence" when it builds
# eltonian_species_occurrence_counts.csv, so we use it here too. Many rows share
# a coordinate: there are tens of thousands of records but only ~1,500 distinct
# locations. We therefore do all spatial work ONCE per distinct location and
# then map records onto locations. This is both faster and exactly equivalent.

say("\n[2] Assembling fungal occurrence points")

fungal_records <- emf |>
  dplyr::filter(!is.na(species), !grepl("_sp$", species)) |>
  dplyr::mutate(fungal_species = trimws(gsub("_", " ", species))) |>
  dplyr::filter(coord_in_canada %in% TRUE, !is.na(lat), !is.na(lon)) |>
  dplyr::select(fungal_species, lat, lon)

# Distinct locations, and the index of each record's location
locs <- fungal_records |>
  dplyr::distinct(lat, lon) |>
  dplyr::mutate(point_id = dplyr::row_number())

fungal_records <- dplyr::left_join(fungal_records, locs, by = c("lat", "lon"))

n_records_total    <- nrow(fungal_records)
n_points           <- nrow(locs)
n_fungi_georef     <- dplyr::n_distinct(fungal_records$fungal_species)
n_fungi_no_georef  <- n_named_fungal_species - n_fungi_georef

say("    georeferenced Canadian records  : ", format(n_records_total, big.mark = ","))
say("    distinct locations              : ", n_points)
say("    named fungal species with >=1 georeferenced Canadian record: ",
    n_fungi_georef, " of ", n_named_fungal_species)
say("    named fungal species with NO georeferenced record          : ",
    n_fungi_no_georef,
    "  (", round(100 * n_fungi_no_georef / n_named_fungal_species, 1),
    "% of columns are structurally zero)")

# An sf point layer of the distinct locations, in WGS84 and in Albers.
pts_wgs   <- sf::st_as_sf(locs, coords = c("lon", "lat"), crs = crs_wgs84,
                          remove = FALSE)

# The equal-area CRS for buffering. The brief specifies the Albers CRS used by
# lakes_canada_albers.gpkg; we read it from that file and check it agrees with
# the `crs_albers` constant defined in 00_setup.R.
albers_from_file <- sf::st_layers(paths$lakes_albers)$crs[[1]]
if (!isTRUE(sf::st_crs(albers_from_file) == sf::st_crs(crs_albers)))
  warning("Albers CRS in lakes_canada_albers.gpkg differs from crs_albers in ",
          "00_setup.R; using the CRS from the file, as specified.")
crs_equal_area <- sf::st_crs(albers_from_file)

pts_aea <- sf::st_transform(pts_wgs, crs_equal_area)


# =============================================================================
# SECTION 3.  Point x host possibility, at five spatial grains
# =============================================================================
# The core data structure is, for each grain, a LIST of length n_points. Element
# i holds the integer indices (into `host_all`) of the host species whose range
# is judged to include location i under that grain. Everything downstream --
# the matrix, the grain comparison, the rarefaction -- is set arithmetic on
# these lists, which keeps the spatial work in one place.
#
# The five grains, from finest to coarsest:
#   (a) point_in_range : the point falls inside the BIEN polygon
#   (b) buffer_25km    : the point is within 25 km of the polygon
#   (b) buffer_100km   : the point is within 100 km of the polygon
#   (c) grid_0.5deg    : the point's 0.5-degree cell is occupied by the host
#   (d) ecoregion      : the point's ecoregion is intersected by the host's range
#   (e) ecozone        : same, at ecozone resolution
#
# They are NOT strictly nested -- (c) uses a rasterization of the same polygons,
# which both adds cells (any polygon touching a cell fills it) and drops cells
# (the stack is masked to Canada) -- but they are ordered by how permissive the
# co-occurrence rule is.

say("\n[3] Computing point x host possibility at each grain")

# A helper: turn an sf sparse index (list of range-row indices per point) into
# indices into `host_all`.
range_idx_to_host_idx <- function(sparse_list, range_species) {
  host_pos <- match(range_species, host_all)   # position of each range in host_all
  lapply(sparse_list, function(ix) sort(unique(host_pos[ix])))
}

grain_lists <- list()

# ---- (a) Point in polygon ---------------------------------------------------
say("    (a) point in BIEN range ...")
sp_a <- sf::st_intersects(pts_wgs, sf::st_make_valid(ranges))
grain_lists[["a_point_in_range"]] <- range_idx_to_host_idx(sp_a, ranges$species)

# ---- (b) Point within a buffered range --------------------------------------
# Rather than materialising buffered polygons (slow and memory-hungry for 227
# complex multipolygons), we ask the equivalent question directly: is the point
# within D metres of the polygon? st_is_within_distance uses a spatial index and
# gives exactly the same answer as buffering, in an equal-area CRS.
ranges_aea <- sf::st_transform(sf::st_make_valid(ranges), crs_equal_area)
for (d in BUFFERS_M) {
  say("    (b) point within ", d / 1000, " km of BIEN range ...")
  sp_b <- sf::st_is_within_distance(pts_aea, ranges_aea, dist = d)
  grain_lists[[paste0("b_buffer_", d / 1000, "km")]] <-
    range_idx_to_host_idx(sp_b, ranges$species)
}

# ---- (c) Shared 0.5-degree grid cell ----------------------------------------
# terra::extract returns one row per point and one column per host layer, with
# 1 = host present in that cell, 0 = absent, NA = cell outside the Canada mask.
say("    (c) shared 0.5-degree grid cell ...")
stack_vals <- terra::extract(host_stack, terra::vect(pts_wgs), ID = FALSE)
stack_host_pos <- match(names(host_stack), host_all)
grain_lists[["c_grid_0.5deg"]] <- lapply(seq_len(nrow(stack_vals)), function(i) {
  v <- as.numeric(stack_vals[i, ])
  sort(stack_host_pos[which(!is.na(v) & v == 1)])
})

# ---- (d) and (e) Shared ecoregion / ecozone ---------------------------------
# Two lookups are needed: which ecoregion polygon each point falls in, and which
# ecoregion polygons each host's range intersects. A point on a boundary can
# match more than one polygon; a point over water or outside the ecoregion
# coverage matches none and can therefore support no host at this grain.
say("    (d,e) shared ecoregion / ecozone ...")
pt_eco   <- sf::st_intersects(pts_wgs, sf::st_make_valid(ecoregions))
rng_eco  <- sf::st_intersects(sf::st_make_valid(ranges), sf::st_make_valid(ecoregions))

n_pts_no_ecoregion <- sum(lengths(pt_eco) == 0L)
if (n_pts_no_ecoregion > 0)
  say("        note: ", n_pts_no_ecoregion, " of ", n_points,
      " locations fall outside the ecoregion coverage (e.g. over water)")

# For each ecoregion polygon, the set of hosts whose range touches it. Inverting
# the host -> ecoregion list once is far cheaper than testing every point.
host_pos_rng   <- match(ranges$species, host_all)
eco_to_hosts   <- vector("list", nrow(ecoregions))
for (i in seq_along(rng_eco)) {
  for (e in rng_eco[[i]]) eco_to_hosts[[e]] <- c(eco_to_hosts[[e]], host_pos_rng[i])
}
eco_to_hosts <- lapply(eco_to_hosts, function(x) sort(unique(x)))

grain_lists[["d_ecoregion"]] <- lapply(pt_eco, function(ix) {
  if (!length(ix)) return(integer(0))
  sort(unique(unlist(eco_to_hosts[ix], use.names = FALSE)))
})

# Ecozone: collapse ecoregion polygons to their ECOZONE code, then repeat.
zone_of_eco  <- ecoregions$ECOZONE
zone_levels  <- sort(unique(zone_of_eco))
zone_to_hosts <- lapply(zone_levels, function(z) {
  sort(unique(unlist(eco_to_hosts[which(zone_of_eco == z)], use.names = FALSE)))
})
grain_lists[["e_ecozone"]] <- lapply(pt_eco, function(ix) {
  if (!length(ix)) return(integer(0))
  z <- match(unique(zone_of_eco[ix]), zone_levels)
  sort(unique(unlist(zone_to_hosts[z], use.names = FALSE)))
})


# =============================================================================
# SECTION 4.  Fungus x host possibility
# =============================================================================
# For each fungal species, take the locations of its records and union the host
# sets of those locations. A fungus with no georeferenced record gets the empty
# set -- correctly, since we have no geographic evidence for it at all.

say("\n[4] Rolling locations up to fungal species")

# Point indices per fungus, WITH multiplicity (needed for record counts) and
# de-duplicated (needed for the union).
pts_by_fungus_dup <- split(fungal_records$point_id, fungal_records$fungal_species)
pts_by_fungus     <- lapply(pts_by_fungus_dup, unique)

possible_hosts_at_grain <- function(grain) {
  gl <- grain_lists[[grain]]
  lapply(pts_by_fungus, function(p) sort(unique(unlist(gl[p], use.names = FALSE))))
}

poss_by_grain <- lapply(setNames(nm = names(grain_lists)), possible_hosts_at_grain)


# =============================================================================
# SECTION 5.  OUTPUT 1 - the geographic-possibility matrix (grain a)
# =============================================================================
# Written in long form. To keep the file compact we write ONLY the cells that
# are possible (possible == 1); ANY host x fungus pair that does not appear in
# this file has possible == 0. The `possible` column is retained so the schema
# is self-describing. The full denominator remains n_host_species x
# n_named_fungal_species = 382,160 cells.
#
# n_records_supporting counts the fungus's georeferenced Canadian records (rows
# of `emf`, the same unit as the Eltonian occurrence counts) that fall inside
# that host's range -- i.e. how much evidence sits behind the "possible" call.

say("\n[5] Output 1: geographic-possibility matrix at grain (a)")

gl_a <- grain_lists[["a_point_in_range"]]

output1 <- purrr::imap_dfr(pts_by_fungus_dup, function(pids, fsp) {
  # tabulate() counts, for each host index, how many of this fungus's records
  # fall in that host's range. pids carries multiplicity, so a host inside the
  # range of 12 records scores 12.
  counts <- tabulate(unlist(gl_a[pids], use.names = FALSE), nbins = n_host_species)
  hit <- which(counts > 0L)
  if (!length(hit)) return(NULL)
  tibble::tibble(host_species = host_all[hit],
                 fungal_species = fsp,
                 possible = 1L,
                 n_records_supporting = counts[hit])
})
output1 <- dplyr::arrange(output1, fungal_species, host_species)

readr::write_csv(output1, file.path(out_dir, "output1_possible_pairs_grain_a.csv"))
say("    possible cells at grain (a): ", format(nrow(output1), big.mark = ","),
    " of ", format(n_matrix_cells, big.mark = ","),
    "  (", round(100 * nrow(output1) / n_matrix_cells, 2), "%)")


# =============================================================================
# SECTION 6.  OUTPUT 2 - grain sensitivity
# =============================================================================
# How much of the full matrix is geographically possible, as the definition of
# "co-occurrence" is relaxed from a point-in-polygon test to a shared ecozone?
# If the answer is highly grain-dependent, the filter is not a stable quantity
# and should not be used as a headline denominator.

say("\n[6] Output 2: grain sensitivity")

grain_labels <- c(
  a_point_in_range = "(a) point inside BIEN range",
  b_buffer_25km    = "(b) point within 25 km of BIEN range",
  b_buffer_100km   = "(b) point within 100 km of BIEN range",
  c_grid_0.5deg    = "(c) shared 0.5-degree grid cell",
  d_ecoregion      = "(d) shared ecoregion",
  e_ecozone        = "(e) shared ecozone"
)

output2 <- purrr::imap_dfr(poss_by_grain, function(poss, grain) {
  n_possible <- sum(lengths(poss))
  tibble::tibble(
    grain                 = grain,
    grain_description     = unname(grain_labels[grain]),
    n_matrix_cells        = n_matrix_cells,
    n_cells_possible      = n_possible,
    prop_cells_possible   = n_possible / n_matrix_cells,
    pct_cells_possible    = round(100 * n_possible / n_matrix_cells, 2),
    mean_possible_hosts_per_fungus =
      round(n_possible / n_named_fungal_species, 1),
    median_possible_hosts_per_georef_fungus =
      stats::median(lengths(poss))
  )
})
output2 <- output2[match(names(grain_labels), output2$grain), ]

readr::write_csv(output2, file.path(out_dir, "output2_grain_sensitivity.csv"))
for (i in seq_len(nrow(output2)))
  say(sprintf("    %-40s %9s cells possible (%5.2f%% of matrix)",
              output2$grain_description[i],
              format(output2$n_cells_possible[i], big.mark = ","),
              output2$pct_cells_possible[i]))


# =============================================================================
# SECTION 7.  OUTPUT 3 - the effort-limitation diagnostic
# =============================================================================
# This is the decisive test. If the size of a fungus's possible-host set keeps
# growing with the number of records available for it, then the filter is
# measuring SAMPLING EFFORT, not geography, and the whole approach is
# uninformative. If it saturates, the filter is measuring something real.
#
# Two pieces of evidence:
#   7a. the observed relationship between set size and record count;
#   7b. rarefaction -- subsample k records per species and see whether the
#       curve flattens. Rarefaction is the cleaner test because it holds the
#       species constant and varies only the effort.

say("\n[7] Output 3: effort-limitation diagnostic")

# ---- 7a. Observed set size vs effort ----------------------------------------
# n_records is the record count (the effort axis the brief specifies).
# n_locations is included alongside it because records at the same coordinate
# add no geographic information at all, so it is the effort that actually
# matters for this estimator.
poss_a <- poss_by_grain[["a_point_in_range"]]

output3_effort <- tibble::tibble(
  fungal_species = names(pts_by_fungus_dup),
  n_records      = lengths(pts_by_fungus_dup),
  n_locations    = lengths(pts_by_fungus),
  n_possible_hosts = lengths(poss_a)
) |>
  dplyr::arrange(dplyr::desc(n_records))

readr::write_csv(output3_effort, file.path(out_dir, "output3_effort_by_species.csv"))

# Spearman correlation: a rank correlation, because the relationship is very
# unlikely to be linear and record counts are heavily right-skewed.
rho_records   <- suppressWarnings(stats::cor(output3_effort$n_records,
                                             output3_effort$n_possible_hosts,
                                             method = "spearman"))
rho_locations <- suppressWarnings(stats::cor(output3_effort$n_locations,
                                             output3_effort$n_possible_hosts,
                                             method = "spearman"))
say("    Spearman rho (possible hosts vs n records)   : ", round(rho_records, 3))
say("    Spearman rho (possible hosts vs n locations) : ", round(rho_locations, 3))

# ---- 7b. Rarefaction --------------------------------------------------------
# For each species and each k, draw k of its records at random (without
# replacement), union the host sets of the sampled locations, and record the
# size. Repeat RAREFY_REPS times and average.
#
# Two species sets are reported, because they answer different questions:
#   "all_eligible" - every species with at least k records. The species set
#                    changes with k, so the curve confounds effort with which
#                    species are included.
#   "n_ge_20"      - only species with >= 20 records, i.e. the SAME species at
#                    every k. This is the clean saturation test, and the one to
#                    read.

set.seed(SEED)  # re-seed here so the rarefaction is reproducible on its own

rarefy_species <- function(pids, k, reps) {
  n <- length(pids)
  if (n < k) return(NA_real_)
  # If a species has exactly k records there is only one possible draw, so one
  # replicate is enough and repeating it would just waste time.
  r <- if (n == k) 1L else reps
  sizes <- vapply(seq_len(r), function(j) {
    s <- if (n == k) pids else sample(pids, k)
    length(unique(unlist(gl_a[s], use.names = FALSE)))
  }, numeric(1))
  mean(sizes)
}

rare_raw <- purrr::map_dfr(RAREFY_K, function(k) {
  vals <- vapply(pts_by_fungus_dup, rarefy_species, numeric(1), k = k,
                 reps = RAREFY_REPS)
  tibble::tibble(k = k,
                 fungal_species = names(pts_by_fungus_dup),
                 n_records = lengths(pts_by_fungus_dup),
                 mean_possible_hosts = vals)
})

summarize_rare <- function(df, label) {
  df |>
    dplyr::filter(!is.na(mean_possible_hosts)) |>
    dplyr::group_by(k) |>
    dplyr::summarise(
      species_set   = label,
      n_species     = dplyr::n(),
      mean_hosts    = mean(mean_possible_hosts),
      median_hosts  = stats::median(mean_possible_hosts),
      q25_hosts     = stats::quantile(mean_possible_hosts, 0.25),
      q75_hosts     = stats::quantile(mean_possible_hosts, 0.75),
      .groups = "drop")
}

output3_rare <- dplyr::bind_rows(
  summarize_rare(rare_raw, "all_eligible"),
  summarize_rare(dplyr::filter(rare_raw, n_records >= 20L), "n_ge_20")
) |>
  dplyr::mutate(seed = SEED, n_reps = RAREFY_REPS) |>
  dplyr::relocate(species_set, .before = k)

readr::write_csv(output3_rare, file.path(out_dir, "output3_rarefaction.csv"))

# Saturation summary: the proportional gain in mean set size from k=10 to k=20,
# for the constant species set. A curve that has saturated shows a small gain.
sat <- dplyr::filter(output3_rare, species_set == "n_ge_20")
gain_10_20 <- with(sat, mean_hosts[k == 20] / mean_hosts[k == 10] - 1)
full_mean_ge20 <- mean(output3_effort$n_possible_hosts[output3_effort$n_records >= 20L])
say("    rarefaction (species with >= 20 records, n = ",
    sat$n_species[sat$k == 20], "):")
for (i in seq_len(nrow(sat)))
  say(sprintf("        k = %2d  mean possible hosts = %5.1f", sat$k[i], sat$mean_hosts[i]))
say("        gain from k=10 to k=20: ", round(100 * gain_10_20, 1), "%")
say("        mean at FULL sampling for the same species: ", round(full_mean_ge20, 1))

# ---- 7c. Figure -------------------------------------------------------------
# Two panels:
#   left  - the raw relationship, one point per fungal species
#   right - the rarefaction curves, with the full-sampling mean as a reference
p_left <- ggplot2::ggplot(output3_effort,
                          ggplot2::aes(x = n_records, y = n_possible_hosts)) +
  ggplot2::geom_point(alpha = 0.25, size = 1.1) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                       colour = "#B2182B", linewidth = 0.8) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Georeferenced Canadian records per fungal species (log scale)",
    y = "Geographically possible host species",
    title = "(A) Possible-host set size vs sampling effort",
    subtitle = sprintf("Grain (a), point in BIEN range; Spearman rho = %.2f",
                       rho_records)) +
  ggplot2::theme_minimal(base_size = 11)

p_right <- ggplot2::ggplot(output3_rare,
                           ggplot2::aes(x = k, y = mean_hosts,
                                        colour = species_set, fill = species_set)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_hosts, ymax = q75_hosts),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = full_mean_ge20, linetype = "dashed",
                      colour = "grey35") +
  ggplot2::annotate("text", x = min(RAREFY_K), y = full_mean_ge20,
                    label = "full sampling (species with >= 20 records)",
                    hjust = 0, vjust = -0.6, size = 3, colour = "grey35") +
  ggplot2::scale_colour_manual(values = c(all_eligible = "#2166AC",
                                          n_ge_20 = "#B2182B"),
                               name = "species set") +
  ggplot2::scale_fill_manual(values = c(all_eligible = "#2166AC",
                                        n_ge_20 = "#B2182B"),
                             guide = "none") +
  ggplot2::scale_x_continuous(breaks = RAREFY_K) +
  ggplot2::labs(
    x = "Records subsampled per fungal species (k)",
    y = "Mean geographically possible host species",
    title = "(B) Rarefaction: does the estimate saturate?",
    subtitle = sprintf("%d replicates per species per k; seed = %d",
                       RAREFY_REPS, SEED)) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

fig <- patchwork::wrap_plots(p_left, p_right, nrow = 1)
ggplot2::ggsave(file.path(out_dir, "output3_effort_limitation.png"), fig,
                width = 11, height = 5, dpi = 200, bg = "white")


# =============================================================================
# SECTION 8.  OUTPUT 4 - sanity check against the observed pairs
# =============================================================================
# Every host x fungus pair actually observed in the data ought to be
# geographically possible: the fungus was recorded somewhere, and the host was
# recorded with it. A pair flagged impossible therefore indicates a failure
# somewhere OTHER than the interaction record -- most likely a BIEN range that
# does not cover where the fungus was found, or a coordinate error.
#
# We separate three distinct failure modes, because they mean different things:
#   no_bien_range        - the host has no modelled range at all, so the test
#                          cannot be run (a coverage gap, not a range error)
#   fungus_no_georef     - the fungus has no georeferenced Canadian record, so
#                          it has no geography to test against
#   range_or_coord_error - both sides testable, and the test failed. THIS is the
#                          rate that bounds the reliability of the exercise.

say("\n[8] Output 4: observed pairs that are flagged impossible")

possible_lookup <- output1 |>
  dplyr::select(host_species, fungal_species) |>
  dplyr::mutate(geo_possible = TRUE)

output4 <- occ_counts |>
  dplyr::left_join(possible_lookup, by = c("host_species", "fungal_species")) |>
  dplyr::mutate(
    geo_possible     = !is.na(geo_possible),
    host_has_range   = host_species %in% hosts_with_range,
    fungus_has_georef = fungal_species %in% names(pts_by_fungus_dup),
    reason = dplyr::case_when(
      geo_possible        ~ "possible",
      !host_has_range     ~ "no_bien_range",
      !fungus_has_georef  ~ "fungus_no_georef",
      TRUE                ~ "range_or_coord_error"
    ))

readr::write_csv(dplyr::filter(output4, !geo_possible),
                 file.path(out_dir, "output4_observed_but_impossible.csv"))

n_obs             <- nrow(output4)
n_obs_impossible  <- sum(!output4$geo_possible)
reason_tally      <- dplyr::count(dplyr::filter(output4, !geo_possible), reason)

say("    observed pairs                              : ", n_obs)
say("    observed but flagged impossible at grain (a): ", n_obs_impossible,
    "  (", round(100 * n_obs_impossible / n_obs, 1), "%)")
for (i in seq_len(nrow(reason_tally)))
  say(sprintf("        %-22s %4d  (%.1f%% of all observed pairs)",
              reason_tally$reason[i], reason_tally$n[i],
              100 * reason_tally$n[i] / n_obs))

# The testable subset is the honest denominator for a "range-map error rate".
testable <- dplyr::filter(output4, host_has_range, fungus_has_georef)
n_testable      <- nrow(testable)
n_testable_fail <- sum(!testable$geo_possible)
say("    testable pairs (host has a range AND fungus has coordinates): ", n_testable)
say("    of those, geographically impossible: ", n_testable_fail,
    "  (", round(100 * n_testable_fail / n_testable, 1),
    "% -- the BIEN range-map / coordinate failure rate)")


# =============================================================================
# SECTION 9.  OUTPUT 5 - recomputed fill statistics
# =============================================================================
# The headline Eltonian numbers, recomputed against the geographically-possible
# denominator at each grain, alongside the submitted full-matrix values.
#
# TWO NUMERATORS ARE REPORTED, and the difference matters:
#   numerator_all      - all 739 observed pairs, matching the submitted
#                        analysis. Some of these are NOT in the possible set
#                        (Output 4), so the fill rate can in principle exceed
#                        what the possible set alone would allow. Directly
#                        comparable to the manuscript.
#   numerator_possible - only the observed pairs that are geographically
#                        possible at that grain. Internally consistent, but no
#                        longer directly comparable.
#
# The singleton statistic is reported two ways for the same reason: as a share
# of observed pairs (this does NOT change with the denominator -- it is a
# property of the observed data) and as a share of denominator cells (which
# does).

say("\n[9] Output 5: recomputed fill statistics")

fill_row <- function(label, description, denom, numer_all, numer_poss,
                     single_all, single_poss) {
  tibble::tibble(
    grain               = label,
    grain_description   = description,
    n_denominator_cells = denom,
    n_filled_all        = numer_all,
    pct_filled_all      = round(100 * numer_all / denom, 3),
    pct_empty_all       = round(100 * (denom - numer_all) / denom, 3),
    n_filled_possible_only   = numer_poss,
    pct_filled_possible_only = round(100 * numer_poss / denom, 3),
    pct_empty_possible_only  = round(100 * (denom - numer_poss) / denom, 3),
    n_singleton_all          = single_all,
    pct_singleton_of_observed_all = round(100 * single_all / numer_all, 1),
    pct_singleton_of_denominator  = round(100 * single_all / denom, 4),
    n_singleton_possible_only     = single_poss,
    pct_singleton_of_observed_possible_only =
      if (numer_poss > 0) round(100 * single_poss / numer_poss, 1) else NA_real_
  )
}

# The submitted, full-matrix baseline.
output5 <- fill_row("full_matrix", "Full matrix (as submitted)",
                    n_matrix_cells, n_cells_filled, n_cells_filled,
                    n_pairs_singleton, n_pairs_singleton)

# One row per grain. For each grain we need to know which observed pairs are
# possible at that grain, so we rebuild a membership test from poss_by_grain.
host_pos_obs   <- match(occ_counts$host_species, host_all)
fungus_of_obs  <- occ_counts$fungal_species

for (g in names(poss_by_grain)) {
  poss <- poss_by_grain[[g]]
  denom <- sum(lengths(poss))
  is_poss <- mapply(function(h, f) {
    ph <- poss[[f]]
    !is.na(h) && !is.null(ph) && h %in% ph
  }, host_pos_obs, fungus_of_obs)
  output5 <- dplyr::bind_rows(
    output5,
    fill_row(g, unname(grain_labels[g]), denom,
             n_cells_filled, sum(is_poss),
             n_pairs_singleton, sum(is_poss & occ_counts$n_occurrences == 1L)))
}

readr::write_csv(output5, file.path(out_dir, "output5_fill_statistics.csv"))
for (i in seq_len(nrow(output5)))
  say(sprintf("    %-40s denom %9s   %% empty (all obs) %7.3f   %% empty (possible obs only) %7.3f",
              output5$grain_description[i],
              format(output5$n_denominator_cells[i], big.mark = ","),
              output5$pct_empty_all[i], output5$pct_empty_possible_only[i]))

# =============================================================================
# SECTION 10.  OUTPUT 6 - headline numbers, in one tidy file
# =============================================================================
# ADDED 2026-08-11 on promoting this script into the pipeline. Everything below
# is already printed to run_log.txt, but the manuscript and the supplemental
# materials need to READ these values rather than have them typed in, so they
# are collected here as a metric/value table like every other summary file in
# the project.

say("\n[10] Output 6: headline numbers for the manuscript and SM1")

# Cells removed at grain (a) purely because a whole row or column is empty --
# a fungus with no georeferenced record, or a host with no modelled range.
# These are missing data rather than geography, and the share matters because
# it bounds how much of the filter's effect is genuinely geographic.
n_hosts_no_range  <- n_host_species - length(hosts_with_range)
n_cells_removed_a <- n_matrix_cells - sum(lengths(poss_by_grain[["a_point_in_range"]]))
n_cells_structural <- n_fungi_no_georef * n_host_species +
  n_hosts_no_range * (n_named_fungal_species - n_fungi_no_georef)

output6 <- tibble::tibble(
  metric = c(
    "Named fungal species with no georeferenced Canadian record",
    "% of named fungal species with no georeferenced Canadian record",
    "Host species with no BIEN modelled range",
    "Matrix cells removed by the finest-grain filter",
    "Cells removed solely because a row or column is structurally empty",
    "% of removed cells attributable to structurally empty rows/columns",
    "Spearman rho, possible hosts vs number of occurrence records",
    "Spearman rho, possible hosts vs number of distinct locations",
    "Rarefaction: % gain in possible hosts from k = 10 to k = 20 records",
    "Observed pairs flagged geographically impossible at the finest grain",
    "% of testable observed pairs flagged geographically impossible"
  ),
  value = c(
    n_fungi_no_georef,
    round(100 * n_fungi_no_georef / n_named_fungal_species, 1),
    n_hosts_no_range,
    n_cells_removed_a,
    n_cells_structural,
    round(100 * n_cells_structural / n_cells_removed_a, 1),
    round(rho_records,   3),
    round(rho_locations, 3),
    round(100 * gain_10_20, 1),
    n_obs_impossible,
    round(100 * n_testable_fail / n_testable, 1)
  )
)

readr::write_csv(output6, file.path(out_dir, "output6_summary.csv"))
for (i in seq_len(nrow(output6)))
  say(sprintf("    %-70s %10s", output6$metric[i], format(output6$value[i])))

say("\n", strrep("=", 78))
say("Done. Outputs written to ", out_dir)
