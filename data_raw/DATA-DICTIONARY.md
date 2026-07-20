# `data_raw/` — Raw Input Data: Acquisition Guide and Data Dictionary

All external input data for the ECM_manuscript pipeline. **This directory is
read-only**: the pipeline never writes here, and nothing in it is produced by
the project's own scripts. Everything the pipeline generates goes to
`data_derived/` instead.

Most files were obtained from external sources and copied into this directory.
A few small base layers (administrative boundaries, ecoregions, Natural Earth,
FungalTraits) are re-downloaded automatically by `01_spatial_data.R` /
`04_combine_ecm_dataset.R` if absent, but copies are already present here.

> **Integrity note.** Two large files dominate this directory: the GlobalFungi
> SH abundance matrix (~13 GB) and the WorldClim climate raster (~1.3 GB). If a
> download is truncated the pipeline can fail silently, so expected sizes are
> given for each entry below.

## Directory structure summary

```
data_raw/                                   # READ-ONLY inputs
├── DATA-DICTIONARY.md                       # this file
├── GlobalFungi/                             # GlobalFungi v5 (~13 GB)
│   ├── GlobalFungi_5_sample_metadata.txt
│   ├── GlobalFungi_5_SH_abundance_ITS1_ITS2.txt
│   └── Metadata_table_description.pdf
├── UNITE/
│   └── sh_general_release_dynamic_04.04.2024_dev.fasta   # pinned reference
├── fungaltraits/
├── fungalroot/
├── biotime/
├── mycocosm/
│   └── mycocosm_organism_list.csv
├── climate/                                 # WorldClim 2.1 (~1.3 GB)
│   └── wc2.1_country/CAN_wc2.1_30s_bio.tif
├── bien2_ranges/                            # BIEN2 modelled ranges (07_bien2_ranges.R)
│   ├── download_log.csv
│   └── <Genus_species>/<model_name>.{dbf,prj,shp,shx}
├── van_Galen_per_sample/
│   └── GFv5_EcM_unassigned_per_sample.csv
├── van_Galen_et_al_dark_taxa_code_and_data/ # dark-taxa raster (10_dark_diversity.R)
│   └── 4.Dark_EcM_taxa_richness_maps/Dark_taxa_geospatial_layers.tif
├── admin_boundaries/                        # GADM (auto via 01_spatial_data.R)
├── ecoregions/                              # auto via 01_spatial_data.R
├── natural_earth/                           # auto via 01_spatial_data.R
└── gbif/                                     # GBIF specimen ZIP (09_linnean.R)
```

---

**On the field tables below.** Each tabular file gets a table of its fields
(name, type, description, units/allowed values), following the recommended
content in UBC's data dictionary guidance
(https://ubc-library-rc.github.io/rdm/content/07_data_dictionary.html#recommended-content).
For a handful of large external reference layers whose schema is mostly
inherited wholesale from the source database (Natural Earth boundaries,
GADM boundaries, GBIF occurrence download, BioTIME raw data), only the fields
this pipeline actually reads are tabled in detail; the remainder of that
source's native schema is noted with a pointer to the source's own
documentation rather than reproduced field-by-field here.

---

## GlobalFungi v5

**Location**: `data_raw/GlobalFungi/`
**Source**: https://globalfungi.com → Downloads
**Citation**: Větrovský T et al. (2020) GlobalFungi. *Scientific Data* 7: 228.

| File | Size | Description |
|---|---|---|
| `GlobalFungi_5_sample_metadata.txt` | ~78 MB | Tab-delimited sample metadata (one row per sample; country, coordinates, barcoding region, sample type, dominant/other plant species, etc.) |
| `GlobalFungi_5_SH_abundance_ITS1_ITS2.txt` | ~13 GB | Wide-format SH abundance matrix: rows = samples, columns = UNITE v10 SH codes, values = read counts |
| `Metadata_table_description.pdf` | ~80 KB | Field descriptions for the metadata table (from GlobalFungi); reproduced in the table below |

**Filters applied by pipeline** — `02_globalfungi.R` retains only:

- `country == "Canada"`
- `barcoding region %in% c("ITS2", "ITSboth")`
- `manipulated == "NO"`
- sample substrate not in `c("shoot", "air", "water", "sediment")`

then extracts the matching columns/rows from the abundance matrix using `awk`
and pivots to long format, dropping zero-abundance records.

> The GlobalFungi species-level abundance matrix
> (`GlobalFungi_5_species_abundance_ITS1_ITS2.txt`, ~4 GB) is **not used** by
> this pipeline and is not included here — the global Wallacean analysis works
> from the SH-level matrix decoded through our own pinned UNITE build.

#### `GlobalFungi_5_sample_metadata.txt` — fields

Source: `Metadata_table_description.pdf` (GlobalFungi). `"NA"` marks missing data throughout.

| Field | Type | Description | Allowed values / units |
|---|---|---|---|
| `sample_ID` | character | Unique identifier of the sample | — |
| `paper_ID` | character | Unique identifier of the source paper | — |
| `paper_title` | character | Title of paper | — |
| `paper_year` | integer | Year paper was published | — |
| `paper_authors` | character | List of paper authors | — |
| `paper_journal` | character | Journal where paper was published | — |
| `paper_doi` | character | DOI of the paper | — |
| `paper_sample_name` | character | Sample name used in the source paper | — |
| `latitude` | numeric | Sample latitude (WGS84) | decimal degrees |
| `longitude` | numeric | Sample longitude (WGS84) | decimal degrees |
| `elevation_study` | numeric | Elevation of sample location relative to sea level | m; negative = below sea level |
| `continent` | character | Continent or ocean | controlled: Africa; Antarctica; Asia; Australia; Europe; North America; South America; Atlantic/Arctic/Indian/Pacific/Southern Ocean |
| `country` | character | Country (marine samples: nearest country) | — |
| `location` | character | Geographic description of sample location | free text, e.g. "Deschutes National Forest, Oregon" |
| `year_of_sampling_from` | integer | Start of sampling year range | — |
| `year_of_sampling_to` | integer | End of sampling year range | — |
| `month_of_sampling` | integer | Month sample was collected | 1–12 |
| `day_of_sampling` | integer | Day sample was collected | 1–31 |
| `sample_type` | character | Defined sample type | controlled: air; algae; coral; deadwood; dust; fungal sporocarp; glacial ice debris; lichen; litter; moss; rhizosphere soil; root; root + rhizosphere soil; sediment; shoot; soil; topsoil; water |
| `sample_type_specification` | character | Free-text detail on `sample_type` | — |
| `environment_type` | character | Coarse environment classification | controlled: aquatic; anthropogenic; cropland; desert; grassland; forest; mangrove; shrubland; tundra; wetland; woodland |
| `ecosystem_classification` | character | Ecosystem classification per ENVO ontology hierarchy | ENVO terms (ontobee.org) |
| `dominant_plant_species` | character | Dominant plant species / sampled organism at the site | — |
| `other_plant_species` | character | Additional vegetation info | — |
| `manipulated` | character | Whether sample was experimentally manipulated | YES / NO |
| `experimental_manipulation_type` | character | Nature of manipulation | controlled: temperature; precipitation; nitrogen; and pairwise/triple combinations |
| `experimental_manipulation_direction` | character | Direction of manipulation (e.g. increase/decrease) | — |
| `experimental_manipulation_vegetation` | character | Vegetation of manipulated site | — |
| `experimental_manipulation_duration` | character | Duration of manipulation before sampling | — |
| `experimental_manipulation_frequency` | character | Frequency of manipulation applied | — |
| `experimental_manipulation_application_detail` | character | Detail on how manipulation was achieved | — |
| `experimental_manipulation_intensity` | character | Detail on manipulation magnitude | — |
| `pH` | numeric | Sample pH, as estimated by source-paper authors | — |
| `pH_method` | character | Method used to measure pH | controlled: H2O; KOH; KCl; CaCl2; other extract |
| `organic_matter_content` | numeric | Organic matter content | % |
| `organic_C_content` | numeric | Organic carbon content | % |
| `total_N_content` | numeric | Total N content | % |
| `total_Ca` | numeric | Total Ca content | ppm |
| `total_P` | numeric | Total P content | ppm |
| `total_K` | numeric | Total K content | ppm |
| `MAT_study` | numeric | Mean annual temperature of sample location | °C |
| `MAP_study` | numeric | Mean annual precipitation of sample location | mm |
| `area_GPS` | numeric | Precision/uncertainty of GPS localization | m² |
| `area_sampled` | numeric | Area covered by sampling | m² |
| `number_of_subsamples` | integer | Number of subsamples pooled into the sample | — |
| `sampling_info` | character | Additional sampling-design detail | — |
| `sample_depth` | character | Depth (or depth range) from sample surface | cm |
| `sample_info` | character | Additional sample metadata | — |
| `DNA_extraction_sample_mass` | numeric | Mass of sample used for DNA extraction | g |
| `DNA_extraction_size` | character | Extraction input size when mass unavailable | — |
| `DNA_extraction_method` | character | DNA extraction method + reference | — |
| `barcoding_region` | character | Barcoding region targeted by sequencing | controlled: ITS1; ITS2; ITSboth — **filtered by pipeline** |
| `PCR_primers` | character | PCR primer identity | — |
| `PCR_primers_sequence` | character | Forward/reverse primer sequences (no barcodes) | — |
| `sequencing_platform` | character | Sequencing platform | controlled: 454Roche; DNBSEQ-G400; Illumina; IonTorrent; Oxford Nanopore; PacBio; SOLiD |
| `ITS1_extracted` | integer | Number of ITS1 sequences extracted | — |
| `ITS2_extracted` | integer | Number of ITS2 sequences extracted | — |
| `ITS_total` | integer | Sum of ITS1 + ITS2 extracted sequences | — |
| `date_added` | date | Date added to GlobalFungi database | — |
| `submitted_by` | character | Submitting person | — |

#### `GlobalFungi_5_SH_abundance_ITS1_ITS2.txt` — structure

Not a conventional field table: this is a wide sample × SH-code matrix.

| Column | Type | Description |
|---|---|---|
| `sample_ID` | character | First column; joins to `GlobalFungi_5_sample_metadata.txt$sample_ID` |
| `SH<code>.<version>FU` | integer | One column per UNITE Species Hypothesis code (e.g. `SH0737578.10FU`); cell value = read count of that SH in that sample (0 = absent) |

---

## UNITE General FASTA (SH reference)

**Location**: `data_raw/UNITE/`
**Source**: https://unite.ut.ee/repository.php → "General FASTA releases"
**Citation**: Nilsson RH et al. (2019) *Nucleic Acids Research* 47: D259–D264.

| File | Size | Description |
|---|---|---|
| `sh_general_release_dynamic_04.04.2024_dev.fasta` | ~76 MB | **Pinned / required.** UNITE v10.0 general FASTA, build 2024-04-04, `_dev` (full-length) variant, DOI `10.15156/BIO/2959332`. One representative sequence per SH code (93,085 SH codes). |

**Format**: standard FASTA. Each header line encodes
`>accession|species|SH_code|taxonomy` (pipe-delimited); the sequence sourced
from the header's SH_code and taxonomy string is what the pipeline parses into
`unite_sh_taxonomy.csv` (a checkpoint).

**Why this build is pinned (not the newest available).** UNITE SH identifiers
are renumbered between dated builds of the same major version. GlobalFungi's
pre-assigned SH codes were generated against the 2024-04-04 build (which
resolves >99% of them, versus ~3% for the 2025-02-19 build), and GenBank SH
codes are assigned in-house by `vsearch` against this same file, so both
sources decode consistently. `paths$unite_fasta` in `00_setup.R` points every
role at this one file. Re-pin to a newer build only after confirming SH-code
match coverage is still high for the GlobalFungi component.

---

## FungalTraits (genus trait table)

**Location**: `data_raw/fungaltraits/`
**Primary source**: Põlme et al. (2020), Supplementary material 4
(`13225_2020_466_MOESM4_ESM.xlsx`, sheet `data` = Table S1, "Traits of genera").
No versioned data DOI is issued for this table; the source-XLSX md5
(`9168998921adc01b1cf63d9cdca8342e`) is recorded in `fungaltraits_version.txt`
as the verifiable fingerprint, and the authoritative frozen copy is the file
archived with this project's own deposit.
**Mirror / derivation route**: globalbioticinteractions/fungaltraits (GitHub),
commit `5deaa2f` — a documented XLSX→CSV conversion of the same MOESM4 sheet.
**Citation**: Põlme S et al. (2020) *Fungal Diversity* 105: 1–16
(DOI 10.1007/s13225-020-00466-2).

| File | Description |
|---|---|
| `FungalTraits_1-2.csv` | **Pinned / required.** Genus-level trait table (the `data` sheet of MOESM4); the single reference read by the pipeline. SHA256 and full provenance in `fungaltraits_version.txt` |
| `fungaltraits_version.txt` | Provenance record: primary-source DOI + supplementary-file name, source-XLSX md5, local SHA256, mirror commit, EcM genus count |
| `data-dictionary-FungalTraits_1-2.csv` | Column dictionary (from the authors); category/importance metadata for the trait fields, reproduced below |

Used by `04_combine_ecm_dataset.R` to assign EcM guild at the genus level
(`primary_lifestyle == "ectomycorrhizal"`) and by `09_linnean.R`,
`12_wallacean_density_map.R`, `13_wallacean_global_comparator.R`, and
`16_raunkiaeran.R`. The reference is **pinned** (like the UNITE FASTA):
`paths$fungaltraits` in `00_setup.R` points every role at this one file, so all
scripts share an identical EcM genus definition. To re-pin, replace the file,
update `fungaltraits_version.txt`, and re-validate the EcM genus count.

#### `FungalTraits_1-2.csv` — fields

One row per fungal genus. Column names below are as they appear in the file
(several carry a `_template` suffix from the authors' original spreadsheet).

| Field | Type | Description | Allowed values |
|---|---|---|---|
| `jrk_template` | character | Author working-column artifact (row/edit marker) — not used by the pipeline | — |
| `Phylum` | character | Fungal phylum | — |
| `Class` | character | Fungal class | — |
| `Order` | character | Fungal order | — |
| `Family` | character | Fungal family | — |
| `GENUS` | character | Fungal genus (join key used by the pipeline) | — |
| `COMMENT on genus` | character | Free-text author comment on the genus call | — |
| `primary_lifestyle` | character | Primary functional guild — **filtered to `"ectomycorrhizal"` by the pipeline to define the EcM genus list** | one of 30 controlled guild terms |
| `Secondary_lifestyle` | character | Secondary functional guild | selection of 30 terms, or free text |
| `Comment_on_lifestyle_template` | character | Free-text comment on guild assignment | — |
| `Endophytic_interaction_capability_template` | character | Endophytic interaction capacity | selection of 7 terms |
| `Plant_pathogenic_capacity_template` | character | Plant-pathogenic capacity | selection of 8 terms |
| `Decay_substrate_template` | character | Substrate(s) the genus can decay | selection of 16 terms |
| `Decay_type_template` | character | Type of decay capacity | selection of 8 terms |
| `Aquatic_habitat_template` | character | Aquatic habitat capacity | selection of 7 terms |
| `Animal_biotrophic_capacity_template` | character | Animal-biotrophic (parasitic) capacity | selection of 19 terms |
| `Specific_hosts` | character | Documented specific host taxa, where known — **read by `16_raunkiaeran.R`** | free text |
| `Growth_form_template` | character | Fungal growth form | selection of 15 terms |
| `Fruitbody_type_template` | character | Fruiting body type | selection of 23 terms |
| `Hymenium_type_template` | character | Fruiting body hymenium type | selection of 7 terms |
| `Ectomycorrhiza_exploration_type_template` | character | EcM hyphal exploration type strategy | selection of 7 terms |
| `Ectomycorrhiza_lineage_template` | character | EcM phylogenetic lineage | selection of 87 terms |
| `primary_photobiont` | character | Primary lichen photobiont (lichenized genera only) | free text |
| `secondary_photobiont` | character | Secondary lichen photobiont | free text |

---

## FungalRoot database (host mycorrhizal status)

**Location**: `data_raw/fungalroot/`
**Source**: GBIF Darwin Core Archive, DOI 10.15468/a7ujmj (species-level route),
and Soudzilovskaia et al. (2020) supplementary Table S2 (genus-level route).
**Citation**: Soudzilovskaia NA et al. (2020) *New Phytologist* 227: 955–966.

| File | Description |
|---|---|
| `744edc21-8dd2-474e-8a0b-b8c3d56a3c2d.232.zip` | GBIF Darwin Core Archive (occurrence-level FungalRoot data). `05_prepare_fungalroot.R` extracts this fresh into `data_derived/checkpoints/fungalroot_dwca_extracted/` on each run rather than keeping a loose extracted copy here. |
| `nph16569-sup-0002-tabless1-s4.xlsx` | Soudzilovskaia et al. (2020) supplementary tables; sheet `Table S2` is the genus-level mycorrhizal-type recommendation table read by the pipeline. |

The DwC-A carries the standard Darwin Core `occurrences.csv`/`measurements.csv`
pair plus `media.csv`, `eml.xml`, `meta.xml`, `dna_data.csv`; only the fields
below are actually read.

#### Fields used by the pipeline

| Source file | Field | Type | Description |
|---|---|---|---|
| `occurrences.csv` | `ID` | character | Darwin Core record ID; joined to `measurements.csv$coreid` |
| `occurrences.csv` | `scientificName` | character | Raw plant binomial (becomes `PlantBinomial`) |
| `measurements.csv` | `coreid` | character | Foreign key to `occurrences.csv$ID` |
| `measurements.csv` | `measurementType` | character | Measurement label; filtered to `"Mycorrhiza type"` |
| `measurements.csv` | `measurementValue` | character | Mycorrhizal-type code for the matched occurrence (becomes `MycorrhizalType`) |
| `nph16569-...xlsx`, sheet `Table S2` | column 1 (renamed `Genus`) | character | Genus name |
| `nph16569-...xlsx`, sheet `Table S2` | column 2 (renamed `MycorrhizalTypeGenus`) | character | Genus-level mycorrhizal-type recommendation (e.g. `EcM`, `EcM-AM`, `AM`, `NM`) |

`occurrences.csv` carries ~50 additional standard Darwin Core occurrence
fields (locality, coordinates, taxonomy ranks, collector, etc.) not read by
this pipeline. See the GBIF Darwin Core Archive documentation
(https://ipt.gbif.org/manual/en/ipt/latest/dwca-guide) for the full term list.

---

## BioTIME database (temporal series)

**Location**: `data_raw/biotime/`
**Source**: https://biotime.st-andrews.ac.uk
**Citation**: Dornelas M et al. (2018) *Global Ecology and Biogeography* 27: 760–786.

Read by `14_prestonian.R` for time-series records with Canadian locations.

| File | Description |
|---|---|
| `biotime_v2_rawdata_2025.rds` | BioTIME's full long-format observation table (one row per taxon × plot × sampling event). Only the fields below are read; the file also carries BioTIME's standard abundance/biomass and sample-description columns. See BioTIME's own data structure documentation (biotime.st-andrews.ac.uk) for the complete schema. |
| `biotime_metadata.csv` | One row per BioTIME study; used to identify fungi studies and their spatial/temporal scope. |

#### Fields used from `biotime_v2_rawdata_2025.rds`

| Field (lower-cased by pipeline) | Type | Description |
|---|---|---|
| `study_id` | character/integer | BioTIME study identifier; joins to `biotime_metadata.csv$STUDY_ID` |
| `genus` | character | Taxon genus |
| `species` | character | Taxon specific epithet |
| `latitude` | numeric | Observation latitude (decimal degrees) |
| `longitude` | numeric | Observation longitude (decimal degrees) |
| `year` | integer | Sampling year |

#### `biotime_metadata.csv` — fields

| Field | Description |
|---|---|
| `STUDY_ID` | Unique study identifier (join key) |
| `REALM` | Terrestrial / Freshwater / Marine |
| `CLIMATE` | Climate zone of the study |
| `GENERAL_TREAT`, `TREATMENT`, `TREAT_COMMENTS`, `TREAT_DATE` | Experimental treatment info, if any |
| `CEN_LATITUDE`, `CEN_LONGITUDE`, `CENT_LAT`, `CENT_LONG` | Study centroid coordinates (decimal degrees) |
| `HABITAT` | Habitat description |
| `PROTECTED_AREA` | Whether the study site is in a protected area |
| `AREA`, `AREA_SQ_KM`, `GRAIN_SIZE_TEXT` | Spatial extent / sampling grain |
| `BIOME_MAP` | Biome classification |
| `TAXA`, `ORGANISMS` | Taxonomic scope of the study — **filtered for `"fungi"` by the pipeline** |
| `TITLE` | Study title |
| `AB_BIO` | Whether abundance and/or biomass data are recorded |
| `DATA_POINTS` | Number of sampling events |
| `START_YEAR`, `END_YEAR` | Temporal coverage |
| `NUMBER_OF_SPECIES`, `NUMBER_OF_SAMPLES`, `NUMBER_LAT_LONG`, `TOTAL` | Study-level summary counts |
| `CONTACT_1`, `CONTACT_2`, `CONT_1_MAIL`, `CONT_2_MAIL` | Data contributor contact info |
| `PERMISSIONS` | Data use permissions |
| `WEB_LINK` | Source URL |
| `DATA_SOURCE` | Provenance of the digitized data |
| `DATE_TO_DB` | Date added to BioTIME |
| `METHODS`, `SUMMARY_METHODS` | Sampling methodology |
| `LINK_ID` | Cross-reference identifier |
| `COMMENTS` | Free-text notes |
| `DATES_CHANGED`, `CURATOR`, `LOC_ADDED`, `DATE_STUDY_ADDED` | Curation metadata |
| `ABUNDANCE_TYPE`, `BIOMASS_TYPE` | Units/type of abundance and biomass measures |
| `SAMPLE_DESC_NAME` | Name of the sample-description field used in the raw data |

---

## JGI MycoCosm organism list (genomes)

**Location**: `data_raw/mycocosm/mycocosm_organism_list.csv`
**Source**: https://mycocosm.jgi.doe.gov/fungi/fungi.info.html
**Purpose**: cross-referenced by `15_darwinian.R` to assess genome availability
for Canadian EcM genera.

| Field | Type | Description |
|---|---|---|
| `id` | character | MycoCosm portal ID for the organism/genome record |
| `taxon_name` | character | Organism scientific name (parsed to genus/species by the pipeline) |
| `assembly_length` | integer | Genome assembly length (bp) |
| `num_genes` | integer | Number of predicted genes |
| `publication` | character | Associated publication reference, if any |

---

## WorldClim 2.1 climate rasters

**Location**: `data_raw/climate/wc2.1_country/CAN_wc2.1_30s_bio.tif`
**Size**: ~1.3 GB
**Source**: https://worldclim.org/data/worldclim21.html → Country data → Canada
**Citation**: Fick SE & Hijmans RJ (2017) *International Journal of Climatology* 37: 4302–4315.

Multi-band bioclimatic raster (30 arc-second, standard WorldClim BIO1–BIO19
band order). `17_hutchinsonian.R` uses **BIO1** (mean annual temperature) and
**BIO12** (mean annual precipitation) to define Canada's climate space; the
other 17 bands are present but not read by the pipeline.

| Band | Variable | Units |
|---|---|---|
| BIO1 | Mean annual temperature | °C — **used** |
| BIO2 | Mean diurnal range | °C |
| BIO3 | Isothermality (BIO2/BIO7 × 100) | % |
| BIO4 | Temperature seasonality (SD × 100) | — |
| BIO5 | Max temperature of warmest month | °C |
| BIO6 | Min temperature of coldest month | °C |
| BIO7 | Temperature annual range (BIO5–BIO6) | °C |
| BIO8 | Mean temperature of wettest quarter | °C |
| BIO9 | Mean temperature of driest quarter | °C |
| BIO10 | Mean temperature of warmest quarter | °C |
| BIO11 | Mean temperature of coldest quarter | °C |
| BIO12 | Annual precipitation | mm — **used** |
| BIO13 | Precipitation of wettest month | mm |
| BIO14 | Precipitation of driest month | mm |
| BIO15 | Precipitation seasonality (CV) | % |
| BIO16 | Precipitation of wettest quarter | mm |
| BIO17 | Precipitation of driest quarter | mm |
| BIO18 | Precipitation of warmest quarter | mm |
| BIO19 | Precipitation of coldest quarter | mm |

---

## BIEN2 modelled host-range shapefiles

**Location**: `data_raw/bien2_ranges/`
**Source**: biendata.org REST API
**Citation**: Moulatlet GM et al. (2025) *PNAS* 122: e2517585122.

SDM-based range polygons for the native Canadian EcM host species, downloaded
programmatically by `07_bien2_ranges.R` (one subdirectory per species,
`<Genus_species>/<model_name>.{dbf,prj,shp,shx}`, plus `download_log.csv`). CRS
is embedded in each `.prj`; the `load_bien2_range()` helper in `00_setup.R`
reads, reprojects to WGS84, and clips to Canada. Consumed by `08_host_rasters.R`.
Only the range polygon geometry is used from each per-species `.shp`/`.dbf`;
the shapefile attribute table (model metadata from the BIEN2 API) is not read.

| `download_log.csv` column | Type | Description |
|---|---|---|
| `species` | character | Binomial species name |
| `status` | character | `success`, `skipped`, `not_available`, `rate_limited`, `network_error`, `unzip_error` |
| `shp` | character | Path to the downloaded `.shp` (if available) |
| `note` | character | Server message or error detail |

---

## van Galen et al. (2025) — per-sample dark fraction

**Location**: `data_raw/van_Galen_per_sample/GFv5_EcM_unassigned_per_sample.csv`
**Source**: Figshare deposit DOI 10.6084/m9.figshare.28830371, folder
`1.Proportion_of_EcM_fungal_OTUs_unassigned/`
**Citation**: van Galen LG et al. (2025) *Current Biology* 35: R563–R574.
**Purpose**: per-sample proportion of EcM SH codes unassigned to species level
for all GlobalFungi v5 samples. Filtered to `Country == "Canada"` in
`09_linnean.R` to provide the Canada-specific dark-fraction benchmark reported
in the manuscript's Linnean section.

| Column | Type | Description |
|---|---|---|
| `Sample_ID` | character | GlobalFungi v5 sample identifier |
| `Prop_EcM_OTUs_unassigned` | numeric | Proportion of EcM SH codes with no UNITE species-level name |
| `Prop_EcM_OTUs_assigned` | numeric | Complement of the above (sums to 1) |
| `Latitude` | numeric (°N) | Sample latitude |
| `Longitude` | numeric (°E) | Sample longitude |
| `Biome` | character | Biome classification |
| `Continent` | character | Continent |
| `Country` | character | Country |
| `Year_of_sampling` | integer | Calendar year |

## van Galen et al. (2025) — dark-taxa richness raster

**Location**: `data_raw/van_Galen_et_al_dark_taxa_code_and_data/4.Dark_EcM_taxa_richness_maps/Dark_taxa_geospatial_layers.tif`
**Size**: ~2.4 GB
**Source**: Figshare deposit DOI 10.6084/m9.figshare.28830371
**Citation**: van Galen LG et al. (2025) *Current Biology* 35: R563–R574.
**Purpose**: multi-band global GeoTIFF; `10_dark_diversity.R` uses the
`percentage_dark_taxa` band (0–100 %) to map undescribed EcM fungal diversity
across Canada (Figure S1).

| Band | Description |
|---|---|
| `dark_taxa_EcM_richness` | estimated absolute richness of dark EcM taxa |
| `percentage_dark_taxa` | percentage of EcM OTUs undescribed (0–100 %) — the band used |
| `dark_taxa_research_priority_metric` | composite research-priority score |

---

## Administrative boundaries (GADM)

**Location**: `data_raw/admin_boundaries/`
**Source**: GADM (www.gadm.org), via the `geodata` R package.

| File | Description |
|---|---|
| `canada_bound_raw.gpkg` | Canada national + provincial/territorial boundaries (level 1; 13 features) |

Auto-downloaded by `01_spatial_data.R` if absent. Used for coordinate
validation (`04_combine_ecm_dataset.R`) and mapping — only the geometry
column is used; the attribute fields below are GADM's native schema and are
not filtered on by the pipeline.

#### Attribute fields

| Field | Type | Description |
|---|---|---|
| `fid` | integer | Feature ID |
| `geom` | geometry | Polygon/multipolygon boundary |
| `GID_1` | character | GADM level-1 unique identifier |
| `GID_0` | character | GADM country identifier (`CAN`) |
| `COUNTRY` | character | Country name |
| `NAME_1` | character | Province/state/territory name |
| `VARNAME_1` | character | Alternate name spellings, pipe-separated |
| `NL_NAME_1` | character | Native-language name, if different |
| `TYPE_1` | character | Local administrative-division type (e.g. "Province") |
| `ENGTYPE_1` | character | English administrative-division type |
| `CC_1` | character | Local administrative-division code |
| `HASC_1` | character | Hierarchical administrative subdivision code |
| `ISO_1` | character | ISO 3166-2 subdivision code |

---

## Ecoregions (National Ecological Framework for Canada)

**Location**: `data_raw/ecoregions/`
**Source**: agr.gc.ca, downloaded by `01_spatial_data.R` (with the ecozone-name
lookup DBF). Used for ecozone-level coverage analyses.

| File | Description |
|---|---|
| `ecoregions.zip` | Original downloaded archive (kept for provenance; the pipeline reads the extracted shapefile below, not this zip) |
| `Ecoregions/ecoregions.shp` (+ `.shx`/`.prj`/`.dbf`/`.sbn`/`.sbx`) | Ecoregion polygons — **used by the pipeline** (`paths$ecoregions_raw`) |
| `Ecoregions/ecozone_names.dbf` | Ecozone code → English/French name lookup table |

#### `ecoregions.dbf` — fields

| Field | Type | Description |
|---|---|---|
| `AREA` | numeric | Polygon area (native units of the source shapefile) |
| `PERIMETER` | numeric | Polygon perimeter |
| `REGION_` | integer | Internal region sequence number |
| `REGION_ID` | integer | Ecoregion identifier |
| `ECOREGION` | integer | Ecoregion code |
| `REGION_NAM` | character | Ecoregion name (English) |
| `REGION_NOM` | character | Ecoregion name (French) |
| `ECOZONE` | integer | Ecozone code — joins to `ecozone_names.dbf$ECOZONE` |

#### `ecozone_names.dbf` — fields

| Field | Type | Description |
|---|---|---|
| `ECOZONE` | numeric | Ecozone code (join key) |
| `NAME_EN` | character | Ecozone name (English), e.g. "Arctic Cordillera" |
| `NAME_FR` | character | Ecozone name (French), e.g. "Cordillère arctique" |

---

## Natural Earth base layers

**Location**: `data_raw/natural_earth/`
**Source**: Natural Earth (www.naturalearthdata.com), via `rnaturalearth`.

| File | Description |
|---|---|
| `canada_ne.gpkg` | Canada boundary (medium scale, Admin-0 countries layer) |
| `lakes_ne.gpkg` | Major lakes (physical lakes layer) |

Only the geometry is used from these files (basemap outlines); no attribute
column is filtered or read by `01_spatial_data.R`. Each carries Natural
Earth's full native attribute schema (~140 fields for `canada_ne.gpkg`,
including multilingual name variants, ISO codes, and population/GDP
estimates; ~39 fields for `lakes_ne.gpkg`) — see Natural Earth's own field
documentation (naturalearthdata.com) for definitions of these unused fields.

---

## GBIF physical specimen occurrence download

**Location**: `data_raw/gbif/`
**Source**: GBIF occurrence download API, via `rgbif::occ_download()` in `09_linnean.R`
**Citation**: GBIF.org GBIF Occurrence Download (Kingdom Fungi, Canada,
PRESERVED/LIVING specimens).

A GBIF Simple-CSV occurrence ZIP (~90 MB) is provided here. `09_linnean.R` reads
this ZIP directly; it only calls the download API (which requires GBIF
credentials — see `README.md`) if the ZIP and its checkpoint CSV are both
absent.

#### Fields (GBIF Simple CSV / Darwin Core terms)

| Field | Description |
|---|---|
| `gbifID` | GBIF's unique record identifier |
| `datasetKey` | Source dataset UUID |
| `occurrenceID` | Source-system occurrence identifier |
| `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` | Taxonomic hierarchy |
| `infraspecificEpithet` | Infraspecific name, if any |
| `taxonRank` | Rank of the identification |
| `scientificName` | Full scientific name as recorded |
| `verbatimScientificName`, `verbatimScientificNameAuthorship` | As originally recorded, unstandardized |
| `countryCode` | ISO country code |
| `locality`, `stateProvince` | Location description |
| `occurrenceStatus` | Present/absent |
| `individualCount` | Number of individuals recorded |
| `publishingOrgKey` | Publishing organization UUID |
| `decimalLatitude`, `decimalLongitude` | Coordinates (WGS84) |
| `coordinateUncertaintyInMeters`, `coordinatePrecision` | Location precision |
| `elevation`, `elevationAccuracy`, `depth`, `depthAccuracy` | Site physical measurements |
| `eventDate`, `day`, `month`, `year` | Collection date |
| `taxonKey`, `speciesKey` | GBIF backbone taxonomy identifiers |
| `basisOfRecord` | e.g. `PRESERVED_SPECIMEN`, `LIVING_SPECIMEN` |
| `institutionCode`, `collectionCode`, `catalogNumber`, `recordNumber` | Specimen/collection identifiers |
| `identifiedBy`, `dateIdentified` | Identification provenance |
| `license` | Data use license |
| `rightsHolder` | Rights holder |
| `recordedBy` | Collector(s) |
| `typeStatus` | Nomenclatural type status, if any |
| `establishmentMeans` | Native/introduced status |
| `lastInterpreted` | Date GBIF last processed the record |
| `mediaType` | Associated media type, if any |
| `issue` | GBIF data-quality flags |

---

## Approximate sizes

| Item | Size |
|---|---|
| GlobalFungi metadata + SH matrix | ~13 GB |
| van Galen dark-taxa raster | ~2.4 GB |
| WorldClim climate raster | ~1.3 GB |
| BioTIME | ~170 MB |
| admin boundaries | ~107 MB |
| GBIF ZIP | ~90 MB |
| BIEN2 ranges | ~85 MB |
| UNITE `_dev` FASTA | ~76 MB |
| everything else | < 15 MB combined |
| **Total** | **~17 GB** |
