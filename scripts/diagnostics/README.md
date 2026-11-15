# Diagnostics

Stand-alone checks that Supplemental Materials SM1 or the reviewer responses
report, but that the main pipeline does not need. `run_all.R` does **not** run
them. Each one reads the pipeline's outputs, so run it after the main pipeline,
from the project root:

```bash
Rscript scripts/diagnostics/<script>.R
```

Outputs go to `data_derived/diagnostics/`.

| Script | What it checks | Run after |
|---|---|---|
| `its1_exclusion_check.R` | What Canada's ITS1-only GlobalFungi samples would add (EcM SH codes, named species, grid cells) if they were not excluded. Several minutes (reads the 13 GB matrix). | `02`–`04` |
| `genbank_location_summary.R` | What the GenBank location rule in `04` kept and dropped, by reason; countries named by dropped records; kept records whose coordinates disagree with their "Canada" field; EcM taxa carried only by dropped records. | `04` |

The GlobalFungi sample-type tallies need no separate check: `02_globalfungi.R`
writes `globalfungi_sample_type_tally.csv` and
`globalfungi_canada_sample_type_filter.csv`, and `18_eltonian.R` writes
`eltonian_sample_type_tally_canada_all.csv`.
