# =============================================================================
# 20_depth_discard.R  —  Sequencing depth discarded by EcM-genus filtering (Fig S3)
# =============================================================================
# PURPOSE
#   Quantify how much of each Canadian GlobalFungi soil sample's total sequencing
#   depth is discarded by the EcM-genus filtering step (FungalTraits). For every
#   soil sample with at least one EcM-classified read, we compute the percentage
#   of total ITS reads that are NOT retained after filtering to EcM genera, and
#   plot the distribution. This is Figure S3 of the manuscript.
#
# INPUTS
#   paths$gf_meta_out   Canadian GlobalFungi sample metadata checkpoint
#                       (data_derived/checkpoints/globalfungi_canada_metadata.csv;
#                       provides sample_type and ITS_total = total fungal ITS
#                       read count per sample, i.e. depth before EcM filtering).
#                       Written by 02_globalfungi.R.
#   paths$emf_data      combined EcM dataset (auto-loaded as `emf`); GlobalFungi
#                       rows give per-sample EcM read sums.
#
# OUTPUT (figures/)
#   Figure-S3_depth_discard.png   (paths$fig_depth_discard)
# =============================================================================

source(here::here("scripts", "00_setup.R"))

if (file.exists(paths$fig_depth_discard)) {
} else if (!(!file.exists(paths$gf_meta_out))) {

  # ---- 1. Canadian soil-sample total sequencing depth (ITS_total) ------------
  soil_meta <- readr::read_csv(paths$gf_meta_out, show_col_types = FALSE) |>
    dplyr::filter(sample_type == "soil") |>
    dplyr::distinct(sample_ID, ITS_total)

  # ---- 2. Per-sample EcM read sum (post-filter depth) ------------------------
  ecm_per_sample <- emf |>
    dplyr::filter(source == "GlobalFungi", sample_ID %in% soil_meta$sample_ID,
                  !is.na(sh_code), !is.na(abundance), abundance > 0) |>
    dplyr::group_by(sample_ID) |>
    dplyr::summarise(ecm_reads = sum(abundance), .groups = "drop")

  # ---- 3. Percentage of total depth discarded by EcM-genus filtering ---------
  # Every soil sample with >= 1 EcM-classified read (no minimum on observed
  # richness, since no richness estimator is being fed).
  dilution_df <- soil_meta |>
    dplyr::left_join(ecm_per_sample, by = "sample_ID") |>
    dplyr::mutate(ecm_reads     = dplyr::coalesce(ecm_reads, 0),
                  pct_discarded = 100 * (1 - ecm_reads / ITS_total)) |>
    dplyr::filter(ecm_reads >= 1)

  discard_q <- stats::quantile(dilution_df$pct_discarded, probs = c(0, .25, .5, .75, 1))

  # ---- 4. Histogram ----------------------------------------------------------
  p_hist <- ggplot2::ggplot(dilution_df, ggplot2::aes(x = pct_discarded)) +
    ggplot2::geom_histogram(binwidth = 2, boundary = 0,
                            fill = "grey60", colour = "white", linewidth = 0.1) +
    ggplot2::geom_vline(xintercept = discard_q[["50%"]], linetype = "dashed",
                        colour = "firebrick", linewidth = 0.6) +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0.01, 0)) +
    ggplot2::labs(x = "% of total sequencing depth discarded by EcM-genus filtering",
                  y = "Number of samples") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  save_fig_formats(paths$fig_depth_discard, p_hist, width = 10, height = 6, dpi = 300)
}

