# `data_derived/` — Derived Data: Provenance Guide and Data Dictionary

Everything the pipeline produces. **This directory starts empty**: the goal of
this repository is to regenerate all outputs *de novo* by running the scripts in
`scripts/` in order (see `scripts/README.md`), which should reproduce the values
reported in the manuscript and its Supplemental Materials.

Files are organized by subdirectory. Each entry names the script that produces
it. Many steps are checkpoint-guarded — a script skips work whose output already
exists — so to force regeneration, delete the relevant file.

---

## Root-level files

| File | Produced by | Description |
|---|---|---|
| `emf_canada_em_only.csv` | `04_combine_ecm_dataset.R` | **Primary analysis dataset.** EcM fungal records from Canada (GlobalFungi + GenBank), filtered to ectomycorrhizal genera via FungalTraits, with UNITE v10 taxonomy and trait columns. One row per (sample × SH code). Auto-loaded as `emf` by `00_setup.R`. Key column `coord_in_canada` (logical): TRUE = coordinates validated inside the GADM Canada boundary; FALSE = coordinates present but outside; NA = no coordinates. |
| `emf_canada_combined.csv` | `04_combine_ecm_dataset.R` | Pre-EcM-filter combined dataset (all fungal guilds, both sources); carries `coord_in_canada`. Large (>200 MB). |
| `globalfungi_canada_long.csv` | `02_globalfungi.R` | Long-format Canadian GlobalFungi records with UNITE taxonomy. One row per (sample × SH code). Large (>150 MB). |
| `genbank_emf_canada_long.csv` | `03_genbank.R` | Long-format Canadian GenBank EcM records with UNITE taxonomy, metadata, and provenance. Carries `host_taxon_raw` (raw extracted host string) and `host_taxon` (canonicalized via `canonicalize_host()`). |
| `clean_fungalroot_species.csv` | `05_prepare_fungalroot.R` | EcM-demonstrated plant species from FungalRoot (GBIF-backbone-matched). Columns: `UpdatedPlantBinomial`, `UpdatedGenus`, `ecm_demonstrated` (always TRUE), `evidence` (qualifying raw mycorrhizal-type label(s)), `PlantBinomial` (raw name(s)). |
| `ecm_native_canada_host_species.csv` | `06_host_species.R` | Native Canadian EcM host species (species-level FungalRoot evidence ∩ BIEN native flora). Columns: `species`, `host_demonstrated` (always TRUE), `growth_form`. This list is the denominator of potential host interactions in the Eltonian and Hutchinsonian analyses. |

---

## `spatial/` — Processed spatial layers

| File | Produced by | Description |
|---|---|---|
| `canada_simple.gpkg` | `01_spatial_data.R` | Dissolved, simplified Canada national boundary (WGS84). |
| `north_america_simple.gpkg` | `01_spatial_data.R` | Canada + contiguous USA + Mexico (WGS84), for basemaps. |
| `canada_ne_albers.gpkg` | `01_spatial_data.R` | Natural Earth Canada boundary, Canada Albers projection. |
| `lakes_canada_albers.gpkg` | `01_spatial_data.R` | Major Canadian lakes, Canada Albers projection. |
| `ecoregions_processed.gpkg` | `01_spatial_data.R` | Canadian ecoregion polygons (WGS84) with ecozone names joined. |
| `ecozone_names.csv` | `01_spatial_data.R` | Ecozone code → English/French name lookup. |
| `bien_host_richness_0.5deg.tif` | `08_host_rasters.R` | Raster: number of native EcM host species with a BIEN2 range per 0.5° cell. |
| `bien_host_species_stack.tif` | `08_host_rasters.R` | Multi-layer raster: one binary presence layer per host species. |
| `bien_host_data_richness_0.5deg.tif` | `08_host_rasters.R` | Raster: number of host species with ≥1 Canadian EcM sequence record per cell. |
| `bien_host_data_proportion_0.5deg.tif` | `08_host_rasters.R` | Raster: proportion of locally present host species with EcM data (data richness / host richness). |
| `bien_ecoregions_with_host_habitat.gpkg` | `08_host_rasters.R` | Ecoregion polygons flagged for predicted EcM host habitat. |

---

## `checkpoints/` — Intermediate pipeline files

Retained to make re-runs cheap; delete a file to force its step to re-execute.
Includes: `globalfungi_canada_ids.txt`, `globalfungi_canada_SH_abundance.txt`,
`globalfungi_canada_metadata.csv`, `unite_sh_taxonomy.csv` (SH → taxonomy lookup
parsed from the UNITE FASTA), `globalfungi_sh_unmatched.csv` (`02_globalfungi.R`);
`genbank_emf_canada.fasta`, `genbank_emf_canada_metadata.csv`,
`genbank_vsearch_hits.txt`, `genbank_fetch_log.txt` (`03_genbank.R`);
`bien2_ecm_host_ranges.gpkg`, `bien_ecm_canada_species.csv`,
`bien_ecm_growthforms.csv`, `gift_growthforms.csv` (`06_host_species.R`);
`gbif_ecm_canada_raw.csv` (`09_linnean.R`); `linnean_inext_per_sample.rds`
(`10_linnean_inext.R`, caches the slow iNEXT loop);
`gf_global_ecm_sample_ids.csv` (`12_wallacean_density_map.R`);
`gf_global_comparator_cheap.csv`, `gf_global_comparator_volume.csv`
(`13_wallacean_global_comparator.R`).

---

## `linnean/` — Linnean shortfall

| File | Description |
|---|---|
| `linnean_summary.csv` | Taxonomic richness at SH/genus/species level; singleton statistics; per-source species-assignment rates; site-based completeness Ĉ; Chao2 SH richness with bootstrap SE and 95% CI; GBIF specimen counts. Backs the in-text Linnean numbers and Table 1/Table 2 in the manuscript. |
| `linnean_genus_coverage.csv` | All FungalTraits EcM genera flagged for whether observed in Canada. |
| `linnean_extrapolation_coverage.csv` | Per-stratum (Canada-wide and per-ecozone) site-based completeness diagnostics: `stratum`, `n_sites`, `n_sh_obs`, `Q1`, `Q2`, `coverage_hat`. **Feeds Table S1.** |
| `linnean_extrapolation_estimators.csv` | Per-stratum asymptotic SH richness estimators (Chao2 etc.) with SE and CI; `singletons` column ("included"). **Feeds Table S1.** |
| `linnean_gbif_ecm_canada.csv` | GBIF physical specimen records for Canadian EcM genera **with** sequence data in Canada. Feeds Figure S3. |
| `linnean_gbif_ecm_nosequence_canada.csv` | GBIF specimen records for EcM genera **without** sequence data in Canada. Feeds Figure S3. |
| `linnean_inext_per_sample.csv` | Per-sample abundance-based rarefaction/extrapolation (iNEXT, q = 0). One row per Canadian GlobalFungi sample: `sample_ID`, `ecm_reads`, `sh_obs`, `lat`, `lon`, `s_est_chao1` (+ 95% CI), `coverage_obs`, `s_ext_2x` (+ CI), `coverage_2x`, `s_est_over_obs`, `inext_status`. |
| `linnean_inext_summary.csv` | Distribution-level summary of the per-sample iNEXT results (median observed and Chao1 richness, median Ĉ, etc.). Backs the in-text per-sample completeness numbers. |
| `linnean_accumulation.rds` | `vegan::specaccum` objects (diagnostic; used for the internal accumulation panel figure, not a manuscript figure). |

---

## `wallacean/` — Wallacean shortfall

| File | Description |
|---|---|
| `wallacean_location_summary.csv` | Counts of unique GlobalFungi sampling locations and records lacking coordinates. |
| `wallacean_sampling_intensity.csv` | Per-taxonomic-level location statistics (mean/median/max, % single-location taxa). Backs Table 3 and Figure 2. |
| `wallacean_locs_per_sh.csv` / `wallacean_locs_per_genus.csv` / `wallacean_locs_per_species.csv` | Per-SH / per-genus / per-species counts of distinct sampling locations (30 arc-second cells). |
| `wallacean_global_gf_locs_per_species.csv` | Per named species: unique 30 arc-second cells globally (mapped via our UNITE SH lookup against the global GlobalFungi SH matrix). |
| `wallacean_global_gf_locs_per_sh.csv` | Per SH code: unique 30 arc-second cells globally. |

---

## `raunkiaeran/` — Raunkiæran shortfall

| File | Description |
|---|---|
| `raunkiaeran_trait_coverage.csv` | Per-trait coverage across the six EcM-relevant FungalTraits columns: `trait`, `trait_label`, `trait_class`, `n_documented`, `n_total`, `pct_documented`. **Backs Table 4.** |
| `raunkiaeran_genus_summary.csv` | Per-genus count of documented traits (of six). |
| `raunkiaeran_trait_distributions.csv` | Value tallies per trait, including an explicit `(undocumented)` bucket. |
| `raunkiaeran_specific_hosts.csv` | Genera with a documented `specific_hosts` entry and the recorded host taxon. |

---

## `eltonian/` — Eltonian shortfall

| File | Description |
|---|---|
| `eltonian_summary.csv` | Interaction-coverage statistics: host × named-species matrix fill rate, % of observed pairs supported by a single record, plus genus-level analogues. Backs the in-text Eltonian numbers. |
| `eltonian_host_matching.csv` | Host-name matching results: raw input, cleaned name (`host_clean`), match status against the native host list, and host-field provenance. **Feeds Table S2.** |
| `eltonian_matrix_sh.csv` / `eltonian_matrix_species.csv` / `eltonian_matrix_genus.csv` / `eltonian_matrix_genus_genus.csv` | Binary host × fungus co-occurrence matrices at SH, named-species, host-species × fungal-genus, and host-genus × fungal-genus resolution (trimmed to observed rows/columns). |
| `eltonian_species_occurrence_counts.csv` / `eltonian_genus_occurrence_counts.csv` | Supporting-record counts per observed host–fungus pair (species and genus level); back the singleton-association statistics. |
| `eltonian_global_host_associations.csv` and `eltonian_global_*` | Global-scope host–fungus associations and per-host / per-fungus coverage used to establish which Canadian hosts/fungi have any documented partner anywhere. |
| `eltonian_sample_type_tally_*.csv`, `eltonian_genbank_tissue_tally_canada.csv` | Diagnostic tallies of the sample/tissue types underlying the host information. |

---

## `hutchinsonian/` — Hutchinsonian shortfall

| File | Description |
|---|---|
| `hutchinsonian_ecoregion_summary.csv` | EcM sampling coverage by Canadian ecozone (area, host-habitat proportion, sampling status). |
| `hutchinsonian_raster_summary.csv` | Summary statistics of 0.5° grid-cell coverage (proportion of host species with EcM data). |
| `hutchinsonian_ecozone_summary.csv` | Ecozone-level sample-threshold counts under three sample definitions. |

---

## `prestonian/` — Prestonian shortfall

| File | Description |
|---|---|
| `prestonian_summary.csv` | Temporal-coverage summary (backs the in-text 99.9% figure). |
| `prestonian_biotime_matches.csv` | BioTIME studies matched to Canadian EcM sampling. |
| `prestonian_study_metadata.csv` / `prestonian_taxon_summary.csv` / `prestonian_location_timeseries.csv` | Matched-study metadata, per-taxon coverage, and time-series records. |

---

## `darwinian/` — Darwinian shortfall

| File | Description |
|---|---|
| `darwinian_summary.csv` | Genomic-resource coverage summary (backs the in-text 95.4% figure). |
| `darwinian_mycocosm_matches.csv` / `darwinian_species_matches.csv` / `darwinian_genus_summary.csv` | Canadian EcM genera/species matched to MycoCosm genome records, and per-genus genome counts. |

---

## Directory structure summary

```
data_derived/                     # starts EMPTY; regenerated by scripts/
├── DATA-DICTIONARY.md            # this file
├── emf_canada_em_only.csv        # primary dataset (04)
├── emf_canada_combined.csv       # pre-filter combined (04)
├── globalfungi_canada_long.csv   # (02)
├── genbank_emf_canada_long.csv   # (03)
├── clean_fungalroot_species.csv  # (05)
├── ecm_native_canada_host_species.csv  # (06)
├── spatial/                      # processed layers + host rasters (01, 08)
├── checkpoints/                  # intermediate/cached files
├── linnean/  wallacean/  raunkiaeran/
├── eltonian/  hutchinsonian/  prestonian/  darwinian/
```
