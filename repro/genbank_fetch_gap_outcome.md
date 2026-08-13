# Outcome — repair of the GenBank silent fetch gap

**Repair run:** 2026-08-12 22:51 → 2026-08-13 01:51 PDT, branch `fix/genbank-fetch-gap`.
Companion to [genbank_fetch_gap_prefix_state.md](genbank_fetch_gap_prefix_state.md),
which records the damaged state this run started from.

## What was done

Only the 3,001 records identified as missing were fetched, and they were appended
to the existing checkpoints. No new Entrez search was issued, the frozen UID list
was not regenerated, and none of the 57,910 already-complete records was
re-requested. The whole pipeline was then re-run over the enlarged sequence set.

## Records recovered

| Step | Identified missing | Recovered | Not recoverable |
|---|---:|---:|---:|
| Sequences (appended to FASTA) | 2,001 | **2,000** | 1 |
| Metadata rows (appended to CSV) | 1,000 | **1,000** | 0 |
| **Total** | **3,001** | **3,000** | **1** |

### The one record that could not be recovered

`KHUX00000000` (uid 2496718099), *"TLS: soil metagenome ribosomal RNA internal
transcribed spacer region, targeted locus study"*, organism **soil metagenome**.

An accession ending in `00000000` is a WGS/TLS **master record**: a project-level
container that indexes component sequences (`KHUX01000001`, `KHUX01000002`, …)
but holds no sequence of its own. NCBI returns an empty body with no error, and
does so reproducibly. There is genuinely nothing to fetch, and the record is a
soil metagenome rather than a fungal specimen, so it would not have contributed
to the dataset in any case.

It is recorded in `data_derived/checkpoints/genbank_no_sequence_available.csv`
and excluded from the completeness requirement — but only after being
re-requested **on its own**, which is what distinguishes "this record has no
sequence" from "this record was dropped by a failed batch".

## Checkpoint state after repair

| Quantity | Before | After | Target |
|---|---:|---:|---:|
| Frozen UID list | 60,911 | 60,911 | 60,911 |
| Metadata rows | 59,911 | **60,911** | 60,911 |
| FASTA sequences | 58,910 | **60,910** | — |
| + confirmed to have no sequence | — | 1 | — |
| = sequences accounted for | 58,910 | **60,911** | 60,911 |

Duplicate accessions in the FASTA: 0. Duplicate uid or accession in the
metadata: 0. UIDs with no metadata row: 0.

### Proof the original records were not altered

The first 56,506,801 bytes of the repaired FASTA hash identically to the whole of
the original file, and the first 18,550,058 bytes of the repaired metadata CSV
hash identically to the whole of the original CSV (SHA-256 values in the
companion document). The repaired files are therefore the original files with new
records appended — no existing record was rewritten, re-serialized or
re-downloaded. `genbank_emf_canada_ids.txt` and `genbank_fetch_log.txt` are
byte-for-byte unchanged.

## Output row counts

| File | Before | After | Change |
|---|---:|---:|---:|
| `genbank_emf_canada_long.csv` | 31,517 | 32,612 | **+1,095** |
| `emf_canada_combined.csv` | 364,785 | 365,880 | **+1,095** |
| `emf_canada_em_only.csv` | 49,016 | 49,270 | **+254** |
| `globalfungi_canada_long.csv` | 333,268 | 333,268 | **0** |

The GlobalFungi table is unchanged and the combined table grew by exactly the
GenBank increase, which confirms the repair is confined to the GenBank branch.

## The two worked cases

Accessions reaching each stage (of the full accession range):

| Range | Stage | Before | After |
|---|---|---:|---:|
| HQ650744–HQ650767 (24) | sequence in FASTA | 0 | **24** |
| | in `genbank_emf_canada_long.csv` | 0 | **15** |
| | in `emf_canada_em_only.csv` | 0 | **15** |
| OQ410728–OQ410963 (236) | sequence in FASTA | 156 | **236** |
| | in `genbank_emf_canada_long.csv` | 95 | **152** |
| | in `emf_canada_em_only.csv` | 73 | **121** |

All 24 HQ sequences were recovered and 15 survive the documented ITS2
minimum-length and 98.5% Species Hypothesis criteria into the final dataset. For
the OQ range, the 80 sequences lost to the defect were recovered and 57 more
records reach the final GenBank table. The records that still drop out do so on
the stated criteria, which were not changed.

## Effect on reported results

Of 140 tracked headline metrics, **95 are unchanged and 45 changed**. Every
change is an increase in coverage except where noted. Selected values:

| Metric | Before | After |
|---|---:|---:|
| GenBank: total EcM sequence records | 9,149 | 9,403 |
| GenBank: unique SH codes (Canadian dataset) | 1,616 | 1,644 |
| Unique UNITE SH codes (combined dataset) | 2,805 | 2,822 |
| Unique named species (combined dataset) | 1,405 | 1,411 |
| Unique named species: GenBank only | 403 | 409 |
| Unique named species: shared GF + GenBank | 653 | 658 |
| Unique named species: GlobalFungi only | 349 | 344 |
| All SH codes: GenBank only | 743 | 760 |
| All SH codes: shared | 873 | 884 |
| All SH codes: GlobalFungi only | 1,189 | 1,178 |
| Unique sampling locations: GenBank | 1,137 | 1,180 |
| Unique sampling locations: combined | 1,623 | 1,666 |
| EcM species absent from MycoCosm | 1,338 | 1,344 |
| EcM named taxa with no BioTIME record | 1,403 | 1,409 |
| Named fungal species with no georeferenced Canadian record | 249 | 243 |
| % of named fungal species with no georeferenced Canadian record | 17.7 | 17.2 |

The decreases in the "GlobalFungi only" rows are **reclassifications, not
losses**: 11 SH codes and 5 named species moved from "GlobalFungi only" to
"shared", because GenBank now also holds them. The totals are exactly conserved
(1,189 − 11 = 1,178 and 873 + 11 = 884).

The SH-code set is a **strict superset** of the pre-fix set: 17 codes added, zero
removed.

### One group of changes is NOT attributable to this repair

Two Eltonian global metrics decreased:

| Metric | Before | After |
|---|---:|---:|
| Canadian EcM species with documented global host associations (GlobalFungi root samples) | 489 | 487 |
| Named EcM fungal species with ≥ 1 documented host species (global scope) | 646 | 644 |

These are **not** caused by the GenBank repair. They are the correction of a
separate, pre-existing staleness. `gf_global_ecm_sh_subset.rds` is a cache keyed
to the EcM dataset's SH codes, guarded only by `file.exists()` with no validity
check. The copy in place was built on **2026-07-01** and had never been rebuilt,
so it still carried **124 SH codes that were not in the EcM dataset even before
this repair** (they had been removed by the 2026-07-19 rebuild and the August
"root-evidence standard" revision). Deleting it as part of this run forced a
rebuild, and the new cache contains **0** codes absent from the current dataset.

The post-repair values are therefore the more correct ones; the pre-repair values
were computed partly from a six-week-stale SH code list.

## Incidental findings (not acted on)

1. **`gf_global_ecm_sh_subset.rds` has no staleness check** (18_eltonian.R, Part
   B2). Documented above. It is now consistent, but nothing prevents it going
   stale again the next time the EcM dataset changes.
2. **`genbank_global_ecm_meta.csv` has the same weakness** (18_eltonian.R, Part
   B4), keyed to the EcM *genus* list and built 2026-07-05. It is missing 2
   genera (*Fuscoboletinus*, *Scutiger*) — but it was missing exactly the same 2
   before the repair, and the repair added **no** new genera, so before and after
   stand on equal footing and no comparison is distorted. It was deliberately not
   rebuilt, because doing so would re-query GenBank globally at today's date and
   destroy the date isolation this repair depends on.
3. **`figures/Figure-05_shortfalls_summary.png/.jpg` is hand-assembled**, not
   produced by any script, so it could not regenerate. It was restored from git
   unchanged. **It is now stale**: it reuses panels from Figures 1, 3, 4, S1, S4
   and S5, and Figures 1, 3 and S4 changed in this run. It needs manual
   reassembly before submission.
4. **The high-resolution `.tif` figures were regenerated** by the pipeline as
   normal; they are gitignored by design and so do not appear in this branch.

## Reproducibility check

`data_derived/run_log.csv` records all 20 pipeline scripts with status `ok`.

`scripts/99_verify_reproducibility.R verify` reports `FAIL`, which is the
expected and correct result here: the tool compares against the pre-fix baseline
committed at the head of this branch, and the whole purpose of this work was to
change those numbers. The differences it lists are the before/after record
summarized above. Once the change is accepted, re-run

```
Rscript scripts/99_verify_reproducibility.R baseline
```

to re-point the baseline at the repaired state.
