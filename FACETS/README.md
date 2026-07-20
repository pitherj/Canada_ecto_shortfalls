# FACETS manuscript sources

Quarto sources for the **FACETS** (Canadian Science Publishing) submission of
*"An assessment of biodiversity data shortfalls for ectomycorrhizal fungi in
Canada."*

These documents read from `data_derived/` and `figures/`, so every reported
statistic, table and figure traces back to the pipeline in `scripts/` rather
than being typed in by hand.

## Contents

| File | Purpose |
|---|---|
| `manuscript_FACETS_final.qmd` | The manuscript, in FACETS style: double-spaced, continuous line numbers, letter paper, 1-inch margins, 12 pt, table captions above. Renders to PDF and Word. |
| `supplemental_materials_SM1_FACETS.qmd` | SM1 — detailed methods. Self-contained HTML. |
| `supplemental_materials_SM2_FACETS.qmd` | SM2 — decisions log. Self-contained HTML. |
| `references.bib` | Shared bibliography for all three documents. |
| `facets.csl` | Official FACETS author–date citation style. |

Rendered outputs (`.docx`, `.pdf`, `.html`, `.tex`) are git-ignored; the
sources are tracked.

## Rendering

Run from the **project root** — paths resolve via `here::here()`:

```bash
quarto render FACETS/manuscript_FACETS_final.qmd            # -> PDF + Word
quarto render FACETS/supplemental_materials_SM1_FACETS.qmd  # -> HTML
quarto render FACETS/supplemental_materials_SM2_FACETS.qmd  # -> HTML
```

Requirements:

- **`data_derived/`** must be present — the documents read it for every inline
  statistic and table. Run `scripts/run_all.R`, or obtain it from the Borealis
  archive (see the root `README.md` > Archiving).
- **`figures/`** is tracked in the repo, so figures need no regeneration unless
  the analysis changes.
- The PDF target needs a working LaTeX (TinyTeX: `quarto install tinytex`).

## FACETS requirements

**Manuscript format.** Set in the YAML: `linestretch: 2`, continuous line
numbers via LaTeX `lineno`, letter paper, 1-inch margins, 12 pt, and
`tbl-cap-location: top`.

**Figure formats.** FACETS does not accept PNG; it requires TIFF or JPG at
≥ 300 dpi. The `save_fig_formats()` helper in `scripts/00_setup.R` emits PNG,
JPG and TIFF for every manuscript figure, from scripts `10`, `11`, `12`, `17`,
`18`, `19` and `20`. The grey-background panels used to assemble Figure 5 are
PNG only.

The `.tif` versions (~332 MB) are git-ignored; they are included in the Borealis deposit.

## Open items before submission

1. **Data and code DOIs** — the `XXX` (code) and `YYY` (data) placeholders at
   lines 288, 297, 299 and 615–616 of `manuscript_FACETS_final.qmd` need the
   real GitHub and Borealis DOIs.
2. **Word double-spacing and line numbers** — the `docx` output does not carry
   these automatically. Enable in Word (Layout ▸ Line Numbers ▸ Continuous;
   select all ▸ line spacing 2.0), or supply a `reference-doc` template. The
   PDF output already complies.
3. **Figure 5 in TIFF/JPG** — it is a hand-assembled composite with no
   generating script, so `save_fig_formats()` does not produce it. Export TIFF
   and JPG versions from the source design file.
4. **Cover letter** — FACETS asks for a significance statement.

## Sources

- FACETS Instructions to Authors: <https://www.facetsjournal.com/for-authors/instructions-to-authors>
- CSP figure preparation (formats/resolution): via the Instructions to Authors, "Preparation of graphic files".
