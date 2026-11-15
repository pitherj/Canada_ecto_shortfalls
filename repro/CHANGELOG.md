# Change log

Why the pipeline is built the way it is. The scripts describe only what the
current code does; the history behind it is recorded here, newest first. Git
history (`git log`) holds the full record of every change.

---

## 2026-11 — data fixes and clean-up (second revision)

Changes that alter results (Phase 2 of the clean-up plan):

- **Canada boundary (`01`).** The outline was simplified at a 0.1-degree
  tolerance, which moved it by up to tens of kilometres (Manning Park and
  Osoyoos, BC, tested as 27 and 42 km outside Canada). It is now simplified at
  1 km in the equal-area projection. The outline decides `coord_in_canada`
  (`04`), masks the host rasters (`08`, `18`), defines Canada's climate space
  (`17`) and filters the Figure S5 points (`19`).
- **GlobalFungi sample types (`02`).** An inclusion list (soil, root)
  replaced an exclusion list (shoot, air, water, sediment), dropping 67 lichen
  samples.
- **GenBank location rule (`03`, `04`).** NCBI has no search field for a
  record's collection country: the former `"Canada"[Country]` query ran as
  `"Canada"[All Fields]`, matching submitter addresses among other things. The
  query now says `[All Fields]`, and `04` keeps a record only if its own
  location field begins with "Canada", or it has no location field and its
  coordinates are inside Canada. The earlier `canada_basis` label recorded only
  whether location evidence existed, not what it said; it was replaced by
  `location_evidence`. Audit: `repro/2026-11-12_genbank_search_only_check.md`
  (kept out of git while the manuscript is under review).
- **EcM family rule (`04`).** FungalTraits v1.2 (2020) predates the 2022 split
  of *Cortinarius* and the 2019 split of *Inocybe*, and UNITE labels some SHs
  only to family. Records FungalTraits cannot classify are now EcM if at least
  two thirds of the FungalTraits genera in their family are EcM. A strict
  "every genus" rule was rejected because FungalTraits places four
  saprotrophic genera in Cortinariaceae.
- **Table S4 (`09`).** Built by sequence matching, replacing
  `09b_table_s4_check.R` (12 November 2026), which relied on a hand-curated
  Index Fungorum synonym table and stopped for a person whenever a new genus
  appeared. Index Fungorum's web service now blocks scripts, and no scriptable
  taxonomic backbone reproduced its synonymy, so names were dropped as the
  basis of the check. GBIF specimen genera use the FungalTraits genus rule
  only, because the family rule admitted a saprotroph (*Galera*).
- **Genome coverage by SH (`15b`).** Added 14 November 2026 to answer a
  reviewer's question about un-named SHs. Hand edits to its lookup table were
  removed so updates run unattended. Its "genus with no genome" count now
  excludes SHs without a genus label.
- **Global comparators (`11` Step 5, `12`, `13`).** They still used the old
  GlobalFungi sample-type exclusion list and the FungalTraits genus list
  alone, unlike the Canadian data. They now use the shared soil-and-root list
  and EcM rule (`GF_SAMPLE_TYPES_KEPT`, `classify_ecm()`).
- **Host-name typos (`00_setup.R`, `18`).** A hand-maintained list (one entry,
  `Abies balamifera` → `Abies balsamea`) was replaced by an automatic rule:
  a Canadian host name not in the host list is corrected when its epithet is
  at most 2 letters, and fewer than a third of its length, from exactly one
  listed congener. Every unmatched name goes to
  `eltonian_host_names_unmatched.csv`. On the current data this corrects
  *Betula payrifera* and *Abies balsames* (4 record-by-plant rows), which the
  old list missed. *Abies balamifera* (4 letters off) is now listed as
  unmatched rather than corrected. The cut-off of 2: every unmatched name
  3 letters from a host species is a distinct species (*Quercus robur* vs
  *Q. rubra*).

Changes that do not alter results (Phase 3):

- **Growth-form override removed (`06`).** The one hand-entered value
  (*Saxifraga oppositifolia* = herb) was no longer used: BIEN now supplies
  it. Species no source resolves would be recorded as `unknown`.

- Live-query results moved from `data_derived/checkpoints/` to
  `data_raw/snapshots/`.

Changes made with the manuscript update (Phase 5), which do not alter results:

- **Genus labels in `linnean_summary.csv` (`09`).** Since the family rule,
  the EcM genera in our data are no longer all FungalTraits EcM genera (124 in
  all, 116 of them on the FungalTraits list). Two labels were renamed so each
  says what it counts: "Unique EcM genera (combined dataset,
  FungalTraits-filtered)" became "Unique EcM genera (combined dataset)", and
  "EcM genera observed in Canada" (with its percentage) became "FungalTraits
  EcM genera observed in Canada". A new metric, "Unique EcM genera of the
  named species (combined dataset)", counts the genera the named species
  belong to (111), which the manuscript had described with the all-records
  figure. The dashboard reads the old labels and must be updated at its
  refresh.
- **GBIF checkpoint read (`09`).** When `09` read its saved GBIF table
  instead of importing the GBIF ZIP, empty text fields became missing values
  and were written as `NA`, so the three GBIF output files differed in bytes
  (not in records or counts) between a from-scratch run and a re-run. The
  checkpoint is now read with only the literal `NA` treated as missing; the
  re-run output is byte-identical to the from-scratch output.
- Shared helpers in `00_setup.R` (EcM definition, coordinate test, FASTA
  input/output, SH assignment), used by `04`, `09` and `15b`.
- Diagnostics in `scripts/diagnostics/`. The ITS1 exclusion check cited in the
  reviewer response had never been saved into the project; it was rewritten.
- History moved from script comments into this file; script headers
  standardized.

---

## 2026-08 — first revision

- **Silent batch loss in downloads (`03`, `05`, `06`, `18`).** Batches that hit
  a transient network or server error were skipped with only a warning (or, in
  `18`, with no warning at all; in `05`, replaced by "no match" rows). Nothing
  compared records returned with records requested, and because the steps are
  checkpointed, every later run reused the incomplete file. The GenBank pull in
  `03` lost 3,001 of 60,911 records this way; failed NSR and BIEN batches in
  `06` cost about 500 and 200 species each; a failed GIFT query pushed species
  to the coarser congener fallback. Every download now retries, checks
  completeness, and writes its checkpoint only once complete; an incomplete
  checkpoint is topped up. Details: `repro/genbank_fetch_gap_prefix_state.md`,
  `repro/genbank_fetch_gap_outcome.md`, `repro/silent_batch_loss_audit.md`.
  - `06`: the NSR checkpoint was natives-only, with no record of which species
    were queried; it now records every species and its status. Files in the
    older format are still read, with a warning.
- **Root-evidence rule (`18`).** GenBank host links were used with no tissue
  test, while GlobalFungi host links were root-samples-only. GenBank records
  are now restricted to root or ectomycorrhizal material, in both the Canadian
  and the global branches.
- **Multi-plant host labels (`18`).** The host-name cleaner kept only the
  first two words of a label, so every plant after the first in labels such as
  "Picea mariana, Picea glauca, Pinus banksiana" was lost (*Picea glauca* was
  named in 36 of 44 Canadian root samples and captured from none). Labels are
  now split first. The tissue category "Mixed root + soil" was renamed "Root
  and soil both mentioned".
- **Stale caches (`18`).** The global GlobalFungi SH cache and the global
  GenBank genus cache were checked only for existence. The SH cache in use
  during July–August 2026 dated from 1 July and carried 124 SH codes no longer
  in the dataset. Both are now checked against the current dataset; the
  GenBank cache records which genera were asked about in a companion file, and
  appends new records without rewriting existing ones (a read/write round trip
  had altered 5 `country_gb` values).
- **Sample-type tallies (`02`, `18`).** Percentages describing GlobalFungi's
  composition had been typed into the manuscript; they are now computed.
- **Geographic-possibility sensitivity check (`21`).** Written on 5 August 2026
  as exploratory work in `temp/`, and added to the pipeline on 11 August when
  its result was cited in the revision.
