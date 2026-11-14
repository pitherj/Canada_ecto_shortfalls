# `scripts/` — Pipeline Overview

The analysis runs as an ordered sequence of numbered R scripts. Every script
begins with `source(here::here("scripts", "00_setup.R"))`, which loads packages,
the canonical `paths` list, shared helper functions, and (once it exists) the
primary dataset `emf`. Run the scripts in numerical order; each is
checkpoint-guarded, so re-running skips any step whose output already exists
(delete an output to force regeneration).

**Checkpoints that are checked rather than simply trusted.** For most steps
"the output exists" is taken to mean "the step is done". Four are stricter,
because a checkpoint built from a batched download or keyed to the EcM dataset
can exist and still be wrong:

| Step | What is checked | What happens if the check fails |
|---|---|---|
| `03_genbank.R` | Sequence and metadata checkpoints must account for every UID in the frozen search list | Only the missing records are fetched and appended (Steps 2, 3 and 3b) |
| `06_host_species.R` | Whether the NSR and BIEN checkpoints record which species were queried | Older files without that record are used, with a warning that their completeness cannot be verified |
| `18_eltonian.R` (B2) | Whether the cached global SH subset still matches the EcM dataset's SH codes | Rebuilt from the GlobalFungi matrix — local and deterministic, so rebuilding is safe |
| `18_eltonian.R` (B4) | Whether the cached global GenBank metadata covers every EcM genus | Only the uncovered genera are queried and appended; cached genera keep their original retrieval date |

So for these four, deleting an output is not the only way to force work: a
checkpoint that is incomplete or out of date is repaired in place on the next
run. A batch that fails every retry stops the script rather than being skipped,
and no checkpoint is written until its check has passed. This behaviour was
added in 2026-08 after a silent download fault cost 3,001 GenBank records — see
`repro/silent_batch_loss_audit.md`.

To run the whole pipeline, use the master runner, which sources `01`–`20` in
order and records per-script status and timing to `data_derived/run_log.csv`:

```r
source(here::here("scripts", "run_all.R"))
```

`run_all.R` has two settings (editable at the top of the file, or via
environment variables): `SKIP_HEAVY` (`ECM_SKIP_HEAVY=true`) skips the two
on-demand steps that scan the ~13 GB GlobalFungi matrix
(`12_wallacean_density_map.R`, `13_wallacean_global_comparator.R`); and
`STOP_ON_ERROR` (`ECM_STOP_ON_ERROR=false`) controls whether a failing script
aborts the run or is logged so the run continues.

You can also run scripts individually in numeric order — each is
checkpoint-guarded, so completed steps are skipped on re-run (with the four
exceptions noted above, which verify their checkpoints and repair them if
necessary).

The manuscript and Supplemental Materials sources in `FACETS/` read from
`data_derived/` to report in-text statistics and tables, so the pipeline must
run before they can be rendered (see `FACETS/README.md`).

After a run, `scripts/99_verify_reproducibility.R` fingerprints the outputs and
compares a later rerun against that baseline; it is excluded from `run_all.R`
and invoked directly.

## External requirements

| Requirement | Needed by | Notes |
|---|---|---|
| `ITSx` (≥ 1.1; developed under v1.1.3) | `03_genbank.R`, `15b_darwinian_sh_match.R` | Extracts the ITS1/ITS2 sub-regions; the script stops with an install hint if it is not on `PATH` |
| `HMMER` (3.x) | `03_genbank.R`, `15b_darwinian_sh_match.R` | ITSx's backend — ITSx will not run without it |
| `vsearch` (≥ 2.x) | `03_genbank.R`, `15b_darwinian_sh_match.R` | SH assignment of GenBank sequences and of MycoCosm genome ITS sequences |
| `awk` | `00_setup.R` helper, `02_globalfungi.R`, `12_wallacean_density_map.R`, `13_wallacean_global_comparator.R` | column/row subsetting of the multi-gigabyte GlobalFungi matrices |
| Internet + NCBI Entrez | `03_genbank.R`, `09b_table_s4_check.R` (first run only), `15b_darwinian_sh_match.R` (only when the MycoCosm list has genomes not yet in its lookup table), `18_eltonian.R` | GenBank fetches (an NCBI API key is recommended) |
| Internet + GBIF credentials | `09_linnean.R` | only if the archived GBIF ZIP in `data_raw/gbif/` is absent (see root `README.md`) |
| Internet (biendata.org) | `07_bien2_ranges.R` | BIEN2 range downloads |
| Internet (BIEN) | `06_host_species.R` | native-flora query |

## Script sequence

Scripts 01–08 assemble the dataset and reference layers; 09–20 run the seven
shortfall analyses (in manuscript order) and the sampling maps.

| # | Script | Purpose | Key inputs | Key outputs | Manuscript items |
|---|---|---|---|---|---|
| 01 | `01_spatial_data.R` | Download + process boundaries, ecoregions, lakes | GADM, Natural Earth, ecoregions (auto) | `data_derived/spatial/*` | basemaps for all maps |
| 02 | `02_globalfungi.R` | Assemble Canadian GlobalFungi records | GlobalFungi metadata + SH matrix, UNITE FASTA | `globalfungi_canada_long.csv` (+ checkpoints, UNITE taxonomy) | dataset |
| 03 | `03_genbank.R` | Fetch + SH-assign Canadian GenBank ITS records | NCBI GenBank, UNITE FASTA | `genbank_emf_canada_long.csv` (+ checkpoints) | dataset, Table 2 |
| 04 | `04_combine_ecm_dataset.R` | Combine sources, validate coords, EcM-filter via FungalTraits | 02 + 03 outputs, FungalTraits | `emf_canada_em_only.csv`, `emf_canada_combined.csv` | primary dataset |
| 05 | `05_prepare_fungalroot.R` | Classify EcM-demonstrated plant species | FungalRoot DwC-A | `clean_fungalroot_species.csv` | host list input |
| 06 | `06_host_species.R` | Native Canadian EcM host species list | 05 output, BIEN native flora | `ecm_native_canada_host_species.csv` | Eltonian denominator |
| 07 | `07_bien2_ranges.R` | Download BIEN2 modelled host ranges | 06 output, biendata.org | `data_raw/bien2_ranges/*` | host rasters input |
| 08 | `08_host_rasters.R` | Rasterize host richness + data coverage | 06 + 07 outputs, `emf` | `bien_host_*` rasters | Figure 4, Eltonian |
| 09 | `09_linnean.R` | Taxonomic accounting, dark fraction, GBIF specimens | `emf`, FungalTraits, UNITE, GBIF, van Galen CSV | `linnean_*` outputs (incl. `linnean_gbif_ecm_nosequence_canada.csv`) | Table 1, Table 2, in-text Linnean |
| 09b | `09b_table_s4_check.R` | Checks each genus that 09 flagged as lacking Canadian sequence data, against Index Fungorum synonymy and targeted GenBank queries | 09 output, `genbank_emf_canada_long.csv`, `emf`, `data_raw/table_s4_check/` | `linnean_table_s4_reconciliation.csv`, `linnean_gbif_ecm_nosequence_confirmed_canada.csv` | Table S4, Figure S5 |
| 10 | `10_dark_diversity.R` | Dark ("undescribed") EcM taxa map | van Galen raster, Canada boundary | `Figure-S1_dark_diversity.png` (+ `..._grey.png`) | Figure S1 |
| 11 | `11_wallacean.R` | Occurrence/occupancy + SDM sufficiency | `emf`, global GF SH matrix | `wallacean_*` | Table 3, Figure 2 |
| 12 | `12_wallacean_density_map.R` | Global GF EcM sampling-density map | global GF matrix + metadata | `Figure-S2_gf_sampling_density_world.png` | Figure S2 |
| 13 | `13_wallacean_global_comparator.R` | Global GF comparator metrics | global GF SH matrix | `gf_global_comparator_*.csv` | Table 1 (GlobalFungi-wide column) |
| 14 | `14_prestonian.R` | Temporal coverage vs. BioTIME | BioTIME, `emf` | `prestonian_*` | in-text Prestonian |
| 15 | `15_darwinian.R` | Genome availability vs. MycoCosm, matched by species name | MycoCosm list, `emf` | `darwinian_summary.csv`, `darwinian_mycocosm_matches.csv`, `darwinian_species_matches.csv`, `darwinian_genus_summary.csv` | in-text Darwinian |
| 15b | `15b_darwinian_sh_match.R` | Genome availability per UNITE SH, matched by sequence: links each MycoCosm genome to a same-strain GenBank ITS sequence, then assigns it to an SH (ITSx + vsearch at 98.5%, as in 03). Covers un-named SHs, which name matching cannot | MycoCosm list, NCBI GenBank, UNITE FASTA, `emf`, `data_raw/mycocosm_its/` | `data_raw/mycocosm_its/*` (lookup + sequences), `darwinian_sh_summary.csv`, `darwinian_sh_coverage.csv`, `darwinian_genome_sh_assignments.csv`, `darwinian_its_needs_review.csv` | reviewer response on un-named SHs (not yet cited in the manuscript) |
| 16 | `16_raunkiaeran.R` | Trait coverage from FungalTraits | FungalTraits, `emf` | `raunkiaeran_*` | Table 4 |
| 17 | `17_hutchinsonian.R` | Climate-space coverage, per-ecozone coverage, ecozone sampling map | `emf`, WorldClim, ecoregions (used to build ecozone polygons) | `hutchinsonian_*` (incl. `hutchinsonian_ecozone_sample_counts.csv`), `Figure-03_climate_gap.png` (+ `..._grey.png`), `Figure-S4_ecozone_sampling_map.png` (+ `..._grey.png`) | Figure 3, Figure S4, Table S1 |
| 18 | `18_eltonian.R` | Host–fungus interaction coverage (Canada + global; root-derived evidence only), host habitat coverage, host bivariate map | `emf`, host list, host rasters (08), global GF/GenBank | `eltonian_*` (incl. `eltonian_host_raster_summary.csv`, `eltonian_pair_confidence.csv`, `eltonian_sensitivity_confidence.csv`, `eltonian_sample_type_tally_canada_all.csv`), `bien_host_data_richness_0.5deg.tif`, `bien_host_data_proportion_0.5deg.tif`, `Figure-04_host_bivariate_map.png` (+ `..._grey.png`) | Figure 4, Tables S2 and S3, in-text Eltonian |
| 19 | `19_sampling_maps.R` | Whole-Canada sampling maps | `emf`, GBIF specimens (09, 09b) | `Figure-01_sampling_map.png` (+ `..._grey.png`), `Figure-S5_gbif_specimens.png` (+ `..._grey.png`), `linnean_gbif_plotted_counts.csv` | Figure 1, Figure S5 |
| 20 | `20_depth_discard.R` | % sequencing depth discarded by EcM-genus filtering | GF metadata (02), `emf` | `Figure-S3_depth_discard.png` | Figure S3 |
| 21 | `21_geographic_possibility.R` | Sensitivity check on the Eltonian matrix denominator: how much of the host × fungus grid is geographically possible, at six spatial grains | `emf`, host list, BIEN ranges (07/08), ecoregions (01), Eltonian outputs (18) | `geo_possibility/output1`–`output5*`, `run_log.txt` | in-text Eltonian sensitivity sentence, SM1 |

## Cross-step dependencies to note

- `07_bien2_ranges.R` needs the host list from `06_host_species.R`; the strict
  order across the block is `05 → 06 → 07 → 08`.
- `12`/`13` (global GlobalFungi comparators) require a full scan of the ~13 GB
  SH matrix and are the most compute-intensive steps; run them deliberately.
- `19_sampling_maps.R` reads the GBIF specimen tables written by `09_linnean.R`
  and `09b_table_s4_check.R`. Delete the Figure S5 files and
  `linnean_gbif_plotted_counts.csv` to force Figure S5 to rebuild after
  either script changes.
- `run_all.R` picks up `09b_table_s4_check.R` automatically and runs it
  straight after `09_linnean.R`; likewise `15b_darwinian_sh_match.R` runs
  straight after `15_darwinian.R`.
- `15b_darwinian_sh_match.R` does not read any output of `15_darwinian.R`,
  but both start from the same MycoCosm list. After refreshing
  `data_raw/mycocosm/mycocosm_organism_list.csv`, re-run both. `15b` then
  searches GenBank only for the newly added genomes, and lists any it cannot
  resolve in `data_derived/darwinian/darwinian_its_needs_review.csv`. Those
  can be resolved by hand in `data_raw/mycocosm_its/mycocosm_its_lookup.csv`
  (see `data_raw/DATA-DICTIONARY.md`). After any hand edit, re-run `15b`.
- `21_geographic_possibility.R` reads `eltonian_summary.csv` and
  `eltonian_species_occurrence_counts.csv` and **stops** unless its own
  recomputed denominator reproduces them exactly, so it must run after
  `18_eltonian.R` and be re-run whenever `18_eltonian.R` changes.

## Conventions

- **Paths**: always via `here::here()` / the `paths` list in `00_setup.R` — no
  absolute or relative paths in analysis code.
- **Figure 5 grey-background variants**: Figure 5 (`Figure-05_shortfalls_summary.png`)
  is a hand-assembled schematic composite manually built from panels taken
  from Figures 1, 3, 4, S1, and S4, placed on a light-grey (`fig5_grey_bg`,
  `#F2F2F2` in `00_setup.R`) backdrop. Scripts `10`, `17`, `18`, and `19`
  therefore each save two versions of their figure(s): the normal
  white-background copy (`paths$fig_*`, used in the manuscript/SI) and an
  additional grey copy (`paths$fig_*_grey`, filename suffixed `_grey`) used
  only when manually assembling Figure 5. Both are written to `figures/`.
- **Namespacing**: package functions are called explicitly (`dplyr::filter()`,
  `sf::st_read()`, …); only `here`, `dplyr`, `tidyr`, `readr`, `ggplot2` are
  attached in `00_setup.R`.
- **Sites**: a "site" is lat/lon rounded to 3 decimals (~100 m), via
  `add_site_id()`.
- **Coordinate systems**: WGS84 (`crs_wgs84`) for coordinate operations, Canada
  Albers (`crs_albers`) for mapping.
