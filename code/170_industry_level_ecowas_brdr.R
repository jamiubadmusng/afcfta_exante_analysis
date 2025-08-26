
# ====================== Industry-level analysis for ECOWAS =======================

rm(list=ls())

library(dplyr)
library(fixest)
library(purrr)
library(broom)
library(ggplot2)
library(readr)
library(tidyr)

# ----- Paths -----
base_dir  <- "C:/Users/muham/EGEI Dissertation/afcfta/input"
grav_file <- file.path(base_dir, "dgd_2_1.rds")
out_dir   <- "C:/Users/muham/EGEI Dissertation/afcfta/output/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----- Country sets -----
afc_countries <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM",
  "COG","CIV","COD","DJI","EGY","GNQ","SWZ","ETH","GAB","GMB","GHA",
  "GIN","GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS",
  "MAR","MOZ","NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","SOM",
  "ZAF","SSD","SDN","TZA","TGO","TUN","UGA","ZMB","ZWE"
)
eco_countries <- c("BEN","BFA","CPV","CIV","GMB","GHA","GIN","GNB","LBR","MLI","NER","NGA","SEN","SLE","TGO")

# ----- Data prep (merge + ECOWAS border variables; ifelse for dummies, time-varying via as.integer(factor())) -----
prepare_data <- function(itpde_file) {
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)

  grav <- readRDS(grav_file)

  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year")) %>%
    mutate(
      ecowas = ifelse(exporter %in% eco_countries & importer %in% eco_countries & exporter != importer, 1, 0),
      ecowas_to_af = ifelse(exporter %in% eco_countries & importer %in% afc_countries & !(importer %in% eco_countries) & exporter != importer, 1, 0),
      af_to_ecowas = ifelse(importer %in% eco_countries & exporter %in% afc_countries & !(exporter %in% eco_countries) & exporter != importer, 1, 0),
      af_n_ecowas  = ifelse(exporter %in% afc_countries & importer %in% afc_countries & !(exporter %in% eco_countries) & !(importer %in% eco_countries) & exporter != importer, 1, 0),
      ecowas_to_row = ifelse(exporter %in% eco_countries & !(importer %in% afc_countries) & exporter != importer, 1, 0),
      af_n_ecowas_to_row = ifelse(exporter %in% afc_countries & !(exporter %in% eco_countries) & !(importer %in% afc_countries) & exporter != importer, 1, 0),
      rr_n_ecowas = ifelse(!(exporter %in% afc_countries) & !(importer %in% afc_countries) & exporter != importer, 1, 0)
    ) %>%
    mutate(
      ecowas            = as.integer(factor(paste(ecowas, year, sep = "_"))),
      ecowas_to_af      = as.integer(factor(paste(ecowas_to_af, year, sep = "_"))),
      af_to_ecowas      = as.integer(factor(paste(af_to_ecowas, year, sep = "_"))),
      af_n_ecowas       = as.integer(factor(paste(af_n_ecowas, year, sep = "_"))),
      ecowas_to_row     = as.integer(factor(paste(ecowas_to_row, year, sep = "_"))),
      af_n_ecowas_to_row= as.integer(factor(paste(af_n_ecowas_to_row, year, sep = "_"))),
      rr_n_ecowas       = as.integer(factor(paste(rr_n_ecowas, year, sep = "_"))),

      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )

  df
}

# ----- Estimation (ECOWAS model) -----
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + rr_n_ecowas +
      ecowas + ecowas_to_af + af_to_ecowas + af_n_ecowas + af_n_ecowas_to_row + ecowas_to_row |
      exp_year + imp_year,
    data = df,
    vcov = ~ pair_id
  )

  broom::tidy(m, conf.int = TRUE) %>%
    filter(term %in% c("ecowas_to_af", "af_to_ecowas")) %>%
    mutate(industry_id = industry_id,
           industry_descr = industry_descr) %>%
    select(industry_id, industry_descr, term, estimate, std.error, conf.low, conf.high, p.value)
}

# ----- Run across 170 industries -----
itpde_files <- file.path(base_dir, sprintf("itpder2_%d.rds", 1:170))

results <- map_dfr(seq_along(itpde_files), function(i) {
  df <- prepare_data(itpde_files[i])
  # Get industry description from the first row (assuming it's the same for all rows in each file)
  industry_descr <- df$industry_descr[1]
  estimate_model(df, i, industry_descr)
})

# Save combined results
write_csv(results, file.path("C:/Users/muham/EGEI Dissertation/afcfta/output/csv/ecowas_border_results.csv"))

# ----- Robust rank-plot helper -----

make_rank_plot <- function(res, var_term, out_name, ylab,
                           # significance filter
                           sig_only = FALSE, alpha = 0.10,
                           # CI/SE filters (percentile-based)
                           se_pctl_max = 0.95,
                           ciwidth_pctl_max = 0.95,
                           # optional absolute caps (set to Inf to disable)
                           se_abs_max = Inf,
                           ci_abs_max = Inf,
                           # relative CI width cap: (high-low)/|estimate|
                           ci_ratio_max = Inf,
                           # final trimming by estimate value
                           est_pctl_range = c(0.05, 0.95)) {
  
  df0 <- res %>%
    dplyr::filter(term == var_term) %>%
    dplyr::mutate(
      estimate  = as.numeric(estimate),
      conf.low  = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      std.error = as.numeric(std.error),
      p.value   = as.numeric(p.value),
      ci_width  = conf.high - conf.low,
      ci_ratio  = (conf.high - conf.low) / pmax(abs(estimate), 1e-12)
    ) %>%
    dplyr::filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high),
                  is.finite(std.error), is.finite(ci_width), is.finite(ci_ratio))
  
  if (sig_only) {
    df0 <- dplyr::filter(df0, !is.na(p.value) & p.value <= alpha)
  }
  if (nrow(df0) == 0L) {
    warning(sprintf("No usable estimates for %s after initial filters.", var_term))
    return(invisible(NULL))
  }
  
  # percentile thresholds
  se_cut  = min(as.numeric(stats::quantile(df0$std.error, probs = se_pctl_max, na.rm = TRUE)), se_abs_max)
  ciw_cut = min(as.numeric(stats::quantile(df0$ci_width,  probs = ciwidth_pctl_max, na.rm = TRUE)), ci_abs_max)
  
  # apply CI/SE filters (drop very wide CIs / very large SEs / very imprecise ratios)
  plot_df <- df0 %>%
    dplyr::filter(std.error <= se_cut,
                  ci_width  <= ciw_cut,
                  ci_ratio  <= ci_ratio_max) %>%
    dplyr::arrange(estimate)
  
  if (nrow(plot_df) == 0L) {
    warning(sprintf("All %s estimates dropped by CI/SE filters; relax thresholds.", var_term))
    return(invisible(NULL))
  }
  
  # final trimming on estimate values
  est_bounds <- as.numeric(stats::quantile(plot_df$estimate, probs = est_pctl_range, na.rm = TRUE))
  plot_df <- plot_df %>%
    dplyr::filter(estimate > est_bounds[1], estimate < est_bounds[2]) %>%
    dplyr::arrange(estimate)
  
  if (nrow(plot_df) == 0L) {
    warning(sprintf("All %s estimates dropped by estimate trimming; relax est_pctl_range.", var_term))
    return(invisible(NULL))
  }
  
  # rank AFTER all filters
  plot_trim <- plot_df %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::select(rank, estimate, conf.low, conf.high)
  
  plot_long <- tidyr::pivot_longer(plot_trim,
                                   cols = c(estimate, conf.low, conf.high),
                                   names_to = "series", values_to = "value") %>%
    dplyr::mutate(series = factor(series,
                                  levels = c("estimate","conf.low","conf.high"),
                                  labels = c("Estimate","Lower 95 CL","Upper 95 CL")))
  
  p <- ggplot2::ggplot(plot_long,
                       ggplot2::aes(x = rank, y = value, color = series, linetype = series, group = series)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = c("Estimate" = "blue",
                                           "Lower 95 CL" = "red",
                                           "Upper 95 CL" = "green")) +
    ggplot2::scale_linetype_manual(values = c("Estimate" = "solid",
                                              "Lower 95 CL" = "dashed",
                                              "Upper 95 CL" = "dashed")) +
    ggplot2::labs(x = "Rank of Estimate in Terms of Size",
                  y = ylab,
                  color = NULL, linetype = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.box = "horizontal")
  
  ggplot2::ggsave(filename = file.path(out_dir, out_name), plot = p, width = 8, height = 5)
}

# ----- Make the two plots -----
make_rank_plot(results, "ecowas_to_af", "graph_ecowas_to_af.pdf",
  ylab = "ECOWAS_to_AfCFTA Border Effect",
  sig_only = FALSE,      
  alpha = 0.10,
  se_pctl_max = 0.95,    
  ciwidth_pctl_max = 0.95,  
  se_abs_max = Inf,         
  ci_abs_max = Inf,         
  ci_ratio_max = 10,        
  est_pctl_range = c(0.05, 0.95))

make_rank_plot(results, "af_to_ecowas", "graph_af_to_ecowas.pdf",
               ylab = "AfCFTA_to_ECOWAS Border Effect",
               sig_only = FALSE,
               alpha = 0.10,
               se_pctl_max = 0.95,    
               ciwidth_pctl_max = 0.95,  
               se_abs_max = Inf,         
               ci_abs_max = Inf,         
               ci_ratio_max = 10,        
               est_pctl_range = c(0.05, 0.95))

# ====================== End of ECOWAS analysis =======================