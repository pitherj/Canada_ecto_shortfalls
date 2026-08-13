# Pre-fix state — GenBank silent fetch gap

**Captured:** 2026-08-12, on branch `fix/genbank-fetch-gap`, at commit `c450d00`
(the last commit before any fix was applied).

## Why this file exists

`scripts/03_genbank.R` lost whole 200-record download batches without saying so.
This file is the in-repo record of the damaged state, captured *before* the
repair, so that a reader can verify what was broken and by how much.

The damaged data files themselves cannot be committed: this repository tracks
code only (`data_derived/` is excluded by `.gitignore`, and two of the affected
outputs exceed GitHub's 100 MB per-file limit). The checksums and counts below
identify them exactly instead. A full byte-level copy of the pre-fix
`data_derived/`, `figures/` and `repro/` trees was taken outside the repository
at capture time.

## The defect

Steps 2 (sequence download, `entrez_fetch`) and 3 (metadata download,
`entrez_summary`) fetched records from NCBI in batches of 200. Each download was
wrapped in `tryCatch(..., error = function(e) { warning(...); NULL })`. When a
batch failed — a transient NCBI timeout — the loop emitted a warning, wrote
nothing for those 200 records, and carried on to the next batch.

Three things then combined to make the loss permanent and invisible:

1. Nothing after either loop compared the number of records returned against the
   number requested, so a short result was never detected.
2. The checkpoint file was written at the end of the loop regardless, so an
   incomplete file looked exactly like a complete one.
3. Both steps are wrapped in `if (!file.exists(...))`, so every later run saw the
   checkpoint already present and skipped the step entirely. The gap could never
   heal itself.

## The frozen search list

The original NCBI search was run once, on 2026-06-30, and its result was frozen:

```
Fetch timestamp : 2026-06-30 23:48:30 PDT
Search query    : "Canada"[Country] AND "Fungi"[Organism] AND ("internal transcribed spacer"[All Fields] OR "ITS"[All Fields])
Total hits      : 60911
Unique UIDs     : 60911
```

`data_derived/checkpoints/genbank_emf_canada_ids.txt` holds those 60,911 UIDs and
is the reference against which completeness is judged.

## Measured shortfall, against the frozen 60,911-UID list

| State | Records |
|---|---:|
| Complete on disk (sequence **and** metadata) | 57,910 |
| Metadata present, sequence **missing** | 2,001 |
| Sequence present, metadata **missing** | 1,000 |
| Overlap between the two gap sets | 0 |
| **Total needing repair** | **3,001** |

Both shortfalls sit within one record of an exact multiple of the 200-record
batch size (2,001 = 10 batches + 1; 1,000 = 5 batches), which is the signature of
whole batches failing rather than individual records being filtered out.

## File-level counts, pre-fix

| File | Records |
|---|---:|
| `checkpoints/genbank_emf_canada_ids.txt` | 60,911 UIDs (all unique) |
| `checkpoints/genbank_emf_canada_metadata.csv` | 59,911 rows (all unique uid and accession) |
| `checkpoints/genbank_emf_canada.fasta` | 58,910 sequences (all unique accession) |

## Downstream output row counts, pre-fix

| File | Rows |
|---|---:|
| `data_derived/genbank_emf_canada_long.csv` | 31,517 |
| `data_derived/globalfungi_canada_long.csv` | 333,268 |
| `data_derived/emf_canada_combined.csv` | 364,785 |
| `data_derived/emf_canada_em_only.csv` | 49,016 |

## Two worked cases

**HQ650744–HQ650767** — 24 sequences from an ectomycorrhizal study in coastal
British Columbia. All 24 have metadata rows carrying valid Canadian country
strings, coordinates, and ITS-appropriate lengths (428–680 bp). **Zero** have a
sequence in the FASTA, so none reaches the final dataset. They met every stated
inclusion criterion; their sequences were simply never downloaded.

**OQ410728–OQ410963** — 236 accessions. All 236 have metadata rows; only **156**
have sequences. The 80 missing sequences are this defect. A further 61 records
are excluded later on legitimate, documented criteria (the ITS2 minimum-length
floor and the 98.5% Species Hypothesis identity threshold) — that is correct
behaviour and is not affected by the repair.

## SHA-256 checksums of the pre-fix files

```
7103b63775e8f733af1d5b61961d1ee989047b083e5e294f76ac4b44db80d29f  data_derived/checkpoints/genbank_emf_canada_ids.txt
8f6a0f13de397e1f922ad54c1d354fff17fb712e3fc799eb7feb785a8badccc8  data_derived/checkpoints/genbank_emf_canada.fasta
1222d0c7d6886b42e988447437dc14cfccf799d90033ed323bbf6ab2afcf06e3  data_derived/checkpoints/genbank_emf_canada_metadata.csv
6b04a93e052f65332095b952d895a35982bb5f785ef59e0ae7e1c4b2fefacbf6  data_derived/checkpoints/genbank_fetch_log.txt
0c3094456e7b416719c1bb16a7a4f108c4371443d6bc1b64c3865f29ede03e6d  data_derived/genbank_emf_canada_long.csv
e1cd8750daf2e354190b7ca636d07eb290f1a454f892fc6d50b7bbb666624a40  data_derived/emf_canada_combined.csv
85300167a0703c2bf460aa9de3cc813f5df150b7f42a3c53731134e4caa92922  data_derived/emf_canada_em_only.csv
7ffdf8919e407b8f53e3986c4b3def8ad400597e4306ed6d4a8a69cfa5dfe13c  data_derived/globalfungi_canada_long.csv
```

## Companion fingerprint files

`repro/baseline_files.csv`, `repro/baseline_metrics.csv` and
`repro/baseline_env.csv` were re-captured at the same moment as this file, from
the same pre-fix outputs (111 files, 140 headline metrics). They replace an
earlier baseline taken on 2026-07-19, which predated the August "Revision 01"
and Eltonian multi-host changes and would therefore have conflated those changes
with this repair. Running

```
Rscript scripts/99_verify_reproducibility.R verify
```

after the repaired pipeline completes compares against this baseline, so the
reported differences are attributable to the GenBank repair alone.

## How the repair was constrained

The manuscript is under peer review, so the repair had to isolate the effect of
this defect from the effect of GenBank having grown or been revised since
2026-06-30. Accordingly the repair fetched **only** the 3,001 missing records and
appended them. The frozen UID list was not regenerated, no new Entrez search was
issued, and none of the 57,910 already-complete records was re-downloaded or
altered.
