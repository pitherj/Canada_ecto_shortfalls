# Archiving plan — Borealis deposit

The **Borealis** deposit (UBC's Dataverse instance) is the citable archive of
record for this project, with a minted DOI. This GitHub repository is the
working home for code and manuscript sources; the two are complementary, not
duplicates.

Collection: **UBC BLERF** — <https://borealisdata.ca/dataverse/UBC_BLERF>

---

## What goes where

| Component | Size | GitHub | Borealis | Rationale |
|---|---|---|---|---|
| `scripts/` | < 1 MB | ✅ | ✅ | Snapshot at submission fixes the analysis to a citable version |
| `FACETS/` sources (`.qmd`, `.bib`, `.csl`) | < 1 MB | ✅ | ✅ | Same |
| `figures/` PNG + JPG | 32 MB | ✅ | ✅ | Needed to render the manuscript from a clean clone |
| `figures/` TIFF | 332 MB | ❌ | ✅ | FACETS submission deliverables; regenerable, and poorly suited to git |
| `data_derived/` | ~1.4 GB | ❌ | ✅ | The analysis outputs; too large for git, and the substance of the deposit |
| `data_raw/` | ~16 GB | ❌ | ⚠️ selective | See below |
| `repro/baseline_*` | < 1 MB | ✅ | ✅ | Evidence that the deposited code reproduces the deposited data |
| Both `DATA-DICTIONARY.md` | < 1 MB | ✅ | ✅ | Documentation travels with the data |

---

## Raw inputs: cite, don't copy

Most source datasets are published, versioned, and retrievable by DOI or a
stable versioned URL. Duplicating ~16 GB of third-party data would add storage
cost without adding retrievability, so the deposit **records how to obtain each
input** rather than copying it. Every entry in `data_raw/DATA-DICTIONARY.md`
carries source, version, access date and licence — enough to rebuild
`data_raw/` from scratch.

Three categories are exceptions.

### 1. Deposited, because citation does not make them retrievable

| Input | Size | Why |
|---|---|---|
| GBIF occurrence download (`data_raw/gbif/*.zip`) | 89 MB | GBIF deletes download files six months after creation — this one on **6 September 2026**. The DOI then resolves to a record of the query, not the data. Re-running the query returns a different result set, because the GBIF index changes continuously. |
| Live-query snapshots (`data_derived/checkpoints/`) | 921 MB | Results from BIEN, NSR, GIFT, MycoCosm and GenBank queries. These services have no versioned releases, so the exact result set cannot be re-obtained. Without them the pipeline is not reproducible. |

### 2. Excluded, because their licences prohibit redistribution

| Input | Licence | Handling |
|---|---|---|
| GADM administrative boundaries | Academic / non-commercial; no redistribution without permission | Not deposited. `data_derived/spatial/canada_simple.gpkg` is GADM-derived and is **also excluded**; `01_spatial_data.R` regenerates it. Note the Canadian ecozone layer used for the ecozone analyses comes from AAFC under the Open Government Licence, not from GADM, and is unaffected. |
| WorldClim 2.1 climate rasters | Academic / non-commercial; no redistribution without permission | Not deposited. Derived outputs are aggregate climate-space statistics, not redistributions of the raster data, so they are unaffected. |

### 3. Cited only — published and stably retrievable

UNITE · GlobalFungi v5 · FungalTraits · FungalRoot · BIEN · BIEN2 ·
BioTIME · MycoCosm · van Galen et al. (2025) · Natural Earth ·
National Ecological Framework for Canada (AAFC ecoregions and ecozone names)

---

## Licensing

Code and data are deposited under **separate licences**, because they are
different kinds of work and software licences map poorly onto datasets.

| Component | Licence |
|---|---|
| Code (`scripts/`, `FACETS/` sources) | MIT |
| Data (`data_derived/`, `figures/`) | CC BY-NC 4.0 |

**The NonCommercial term is inherited, not chosen.** FungalRoot and the GBIF
occurrence download are both CC BY-NC 4.0; products derived from them cannot
carry a more permissive licence.

> **Do not accept the Dataverse CC0 default.** CC0 would misrepresent the
> rights held in this dataset. Set the licence explicitly to CC BY-NC 4.0.

No input imposes a NoDerivatives or ShareAlike term, so nothing else propagates
onto the derived data. Per-source terms for all 20 inputs are tabulated in
`data_raw/DATA-DICTIONARY.md` > Licence summary.

**Four sources state no licence** — GlobalFungi, FungalTraits, GIFT and
MycoCosm. The dictionary records that fact rather than inferring terms. All
four are published, openly distributed research resources and all four are
cited; if any later states terms incompatible with redistribution, the affected
derived files can be withdrawn from the deposit without disturbing the rest.

---

## Deposit checklist

**Before depositing**

- [ ] Run the full pipeline cold and capture a fresh reproducibility baseline
      (`scripts/99_verify_reproducibility.R baseline`)
- [ ] Re-render the manuscript and supplements; confirm reported values are unchanged
- [x] Export Figure 5 in JPG (hand-assembled; no generating script). No TIFF:
      it would exceed 40 MB, and FACETS accepts JPG at >= 300 dpi.

**At deposit**

- [ ] Create the dataset in the UBC BLERF collection
- [ ] Set the licence explicitly to **CC BY-NC 4.0** (not CC0)
- [ ] Upload in the structure below
- [ ] Reserve the DOI

**After the DOI is minted**

- [ ] Replace the `XXX` (code) and `YYY` (data) placeholders in
      `FACETS/manuscript_FACETS_final.qmd`
- [ ] Add the DOI to `README.md` > How to cite
- [ ] Tag the matching GitHub release so the two archives correspond

---

## Deposit structure

```
<Borealis dataset>
├── README.md                      # repo root README, describing the whole project
├── ARCHIVING.md                   # this file
├── code/
│   ├── scripts/                   # the full pipeline, incl. 99_verify_reproducibility.R
│   └── FACETS/                    # manuscript + supplement sources, references.bib, facets.csl
├── data_derived/
│   ├── DATA-DICTIONARY.md
│   ├── <analysis outputs>
│   └── checkpoints/               # live-query snapshots (see above)
├── data_raw/
│   ├── DATA-DICTIONARY.md         # acquisition guide for all inputs
│   └── gbif/                      # the only raw input deposited in full
├── figures/                       # PNG, JPG and TIFF
└── repro/                         # reproducibility baseline
```

`data_raw/DATA-DICTIONARY.md` is the key document for anyone rebuilding the
analysis: it is the acquisition manifest for the ~16 GB of inputs that are
cited rather than copied.
