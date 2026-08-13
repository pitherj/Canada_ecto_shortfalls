# Audit — silent batch loss across the pipeline

**Date:** 2026-08-13. Branch `fix/silent-batch-loss-remaining`.

Companion to [genbank_fetch_gap_prefix_state.md](genbank_fetch_gap_prefix_state.md).
That document describes a defect in `03_genbank.R` that silently lost 3,001 of
60,911 GenBank records. This document records the sweep of the rest of `scripts/`
for the same pattern, and what was done about each instance.

## The pattern being looked for

```r
result <- tryCatch(fetch(batch), error = function(e) { warning(...); NULL })
# ... failed batches dropped, checkpoint written regardless ...
```

Three ingredients have to be present for a transient network error to become
permanent, invisible data loss:

1. a failed batch is **skipped** rather than retried;
2. nothing compares **records returned against records requested**, so the
   shortfall is never noticed;
3. the step is guarded by `if (!file.exists(...))`, so the truncated checkpoint
   is **reused by every later run** and can never heal.

## Findings

| Location | Batch | Severity | Status |
|---|---|---|---|
| `06_host_species.R:116` NSR native status | 500 spp | **All three ingredients** | Fixed |
| `06_host_species.R:173` BIEN growth form | 200 spp | **All three ingredients** | Fixed |
| `06_host_species.R:228` GIFT traits | single call | Degrades silently | Fixed |
| `05_prepare_fungalroot.R:183` GBIF backbone | 200 names | Degrades silently | Fixed |
| `18_eltonian.R:1349` GenBank global pagination | 200 recs | **All three, and no warning at all** | Fixed |
| `18_eltonian.R` B2 `gf_global_ecm_sh_subset.rds` | — | Stale cache, no validity check | Fixed |
| `18_eltonian.R` B4 `genbank_global_ecm_meta.csv` | — | Stale cache, no validity check | Fixed |
| `02_globalfungi.R:197` awk extraction | — | Partial file promotable on crash | Fixed |
| `07_bien2_ranges.R:130` BIEN range download | per species | **None — already correct** | No change |
| `09_linnean.R:239` GBIF download | single call | Deliberate: credential-gated optional step | No change |
| `14_prestonian.R:69` BioTIME download | single call | Already `stop()`s | No change |
| `15_darwinian.R:58` MycoCosm read | single call | Delimiter fallback, not error suppression | No change |
| `08_host_rasters.R:62` shapefile read | per species | Per-species skip, handled downstream | No change |
| `10/17/19` lakes layer | single call | Cosmetic map layer fallback | No change |

`02_globalfungi.R` was the obvious candidate but proved largely clean: it has no
network batch loops at all, and its `awk` step already stopped on a non-zero exit
status. Its only weakness was that the shell redirection wrote directly to the
final checkpoint path, so a crash mid-write could leave a partial file that a
later run would accept.

`07_bien2_ranges.R` is the model the fixes follow: it retries, records a
per-species status, and re-attempts failures on the next run.

## What was changed

Three shared helpers were added to `00_setup.R` (purely additive — no existing
behaviour changes):

- **`fetch_with_retry()`** — runs a request up to five times with a backing-off
  pause. Returns `NULL` only when every attempt failed, which callers must count
  and treat as a loss, never as "no data".
- **`assert_fetch_complete()`** — compares returned against requested and stops
  unless both hold: no batch exhausted its retries, and the shortfall is within a
  tolerance set far below one batch. Called **before** the checkpoint is written.
- **`write_checkpoint_atomically()`** — writes to a `.partial` sibling and
  renames it into place only on success, so a run that dies leaves an obviously
  incomplete file rather than a plausible-looking checkpoint.

Two checkpoints were made self-describing, so that their completeness can be
judged later rather than being taken on trust:

- `bien_nsr_native_species.csv` now stores **every species queried** with the
  status NSR returned, instead of natives only. Natives are derived on read.
- `bien_ecm_growthforms_queried.csv` (new) records the species actually sent to
  BIEN, because BIEN returns nothing for a species it has no data on — so
  without it, "no trait data" and "batch lost" are indistinguishable.

Files in the older formats are still accepted, with a warning that states plainly
that their completeness cannot be verified and how to force a checked rebuild.
Nothing is silently regenerated, because rebuilding re-queries live services and
could change the host list while the manuscript is under review.

## Verification

- All three helpers carry unit tests covering success, recovery on a later
  attempt, exhaustion, the two guard trip conditions, tolerance behaviour, and
  that a died-mid-write leaves no checkpoint. All pass.
- `06_host_species.R` reproduces `ecm_native_canada_host_species.csv`
  **byte-for-byte** against the existing checkpoints, confirming the changes are
  behaviour-preserving on unchanged inputs.
- `02_globalfungi.R` and `05_prepare_fungalroot.R` re-run cleanly with outputs
  unchanged.
- The two `18_eltonian.R` staleness checks were dry-run against the real
  checkpoints: B2 reports the SH cache consistent (0 orphaned codes, so it is
  reused rather than needlessly rebuilt), and B4 correctly identifies exactly 2
  uncovered genera while preserving the 122 already cached.

## Outcome of the B4 genus top-up

The B4 staleness check found that `genbank_global_ecm_meta.csv` did not cover two
EcM genera, **_Fuscoboletinus_** and **_Scutiger_**, and had not since it was
built on 2026-07-05. `18_eltonian.R` was then run to top them up.

**Both genera returned zero GenBank records worldwide.** The queries succeeded;
there is simply nothing deposited under either name. Both are older generic names
that the FungalTraits reference still carries but under which GenBank holds no
ITS sequence data.

Consequently:

- `genbank_global_ecm_meta.csv` is **byte-for-byte unchanged** (154,336 rows, 122
  genera with records). Nothing was appended because there was nothing to append.
- **Every one of the 22 Eltonian output files is unchanged**, including
  `eltonian_summary.csv`. No reported number moved.
- The new sidecar `genbank_global_ecm_genera_queried.csv` now records all **124**
  genera as asked about, so the two zero-record genera will not be re-queried on
  every future run.

The gap was therefore real but empty: the cache had never asked about those two
genera, and asking produced nothing. This is the reassuring result — it confirms
the Eltonian global figures were not understated by the missing coverage.

## Two bugs found in this branch's own code, and fixed

Running the top-up exposed two defects in the first draft of these changes. Both
are recorded here because they are instructive, and because the second is a trap
that the rest of this pipeline could fall into.

1. **`genus_results[[i]] <- NULL` deletes the element rather than storing it.**
   In R, assigning `NULL` into a list slot removes that slot, shortening the list
   and throwing every later index out of step. The loop crashed with "subscript
   out of bounds". Fixed by using `genus_results[i] <- list(res)`, which stores a
   `NULL` instead of deleting. (The original pre-fix code had the same
   construction but never read the element back, so the fault was latent.)

2. **`NULL` was doing double duty** for "the query failed" and "this genus has no
   records", so a genus with nothing deposited would have tripped the
   completeness guard as though the fetch had broken. Fixed by returning an empty
   table for the legitimate zero-record case and reserving `NULL` for failure.

A third issue was caught by verification rather than by a crash: rewriting the
whole cache file after a top-up **silently altered 5 cached values**, because
`readr::read_csv()` trims leading whitespace from unquoted fields, so a
read/write round-trip edits the data. The affected cells were malformed
`country_gb` entries that reach no result, but the principle matters — topping up
is supposed to leave cached records exactly as retrieved, bytes included. The
top-up now **appends** to a copy of the cache instead of rewriting it, and skips
the write entirely when there is nothing to add.
