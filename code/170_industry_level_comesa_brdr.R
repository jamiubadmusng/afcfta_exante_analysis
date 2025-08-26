# ====================== Industry-level analysis for COMESA =======================

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
com_countries <- c("BDI","COM","COD","DJI","EGY","SWZ","ETH","KEN","LYB","MDG","MWI","MUS","RWA","SYC","SOM","SDN","TUN","UGA","ZMB","ZWE")

# ----- Data prep (merge + COMESA border variables; ifelse for dummies, time-varying via as.integer(factor())) -----
prepare_data <- function(itpde_file) {
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)

  grav <- readRDS(grav_file)

  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year")) %>%
    mutate(
      comesa = ifelse(exporter %in% com_countries & importer %in% com_countries & exporter != importer, 1, 0),
      comesa_to_af = ifelse(exporter %in% com_countries & importer %in% afc_countries & !(importer %in% com_countries) & exporter != importer, 1, 0),
      af_to_comesa = ifelse(importer %in% com_countries & exporter %in% afc_countries & !(exporter %in% com_countries) & exporter != importer, 1, 0),
      af_n_comesa  = ifelse(exporter %in% afc_countries & importer %in% afc_countries & !(exporter %in% com_countries) & !(importer %in% com_countries) & exporter != importer, 1, 0),
      comesa_to_row = ifelse(exporter %in% com_countries & !(importer %in% afc_countries) & exporter != importer, 1, 0),
      af_n_comesa_to_row = ifelse(exporter %in% afc_countries & !(exporter %in% com_countries) & !(importer %in% afc_countries) & exporter != importer, 1, 0),
      rr_n_comesa = ifelse(!(exporter %in% afc_countries) & !(importer %in% afc_countries) & exporter != importer, 1, 0)
    ) %>%
    mutate(
      comesa            = as.integer(factor(paste(comesa, year, sep = "_"))),
      comesa_to_af      = as.integer(factor(paste(comesa_to_af, year, sep = "_"))),
      af_to_comesa      = as.integer(factor(paste(af_to_comesa, year, sep = "_"))),
      af_n_comesa       = as.integer(factor(paste(af_n_comesa, year, sep = "_"))),
      comesa_to_row     = as.integer(factor(paste(comesa_to_row, year, sep = "_"))),
      af_n_comesa_to_row= as.integer(factor(paste(af_n_comesa_to_row, year, sep = "_"))),
      rr_n_comesa       = as.integer(factor(paste(rr_n_comesa, year, sep = "_"))),

      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )

  df
}

# ----- Estimation (COMESA model) -----
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + rr_n_comesa +
      comesa + comesa_to_af + af_to_comesa + af_n_comesa + af_n_comesa_to_row + comesa_to_row |
      exp_year + imp_year,
    data = df,
    vcov = ~ pair_id
  )

  broom::tidy(m, conf.int = TRUE) %>%
    filter(term %in% c("comesa_to_af", "af_to_comesa")) %>%
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
write_csv(results, file.path("C:/Users/muham/EGEI Dissertation/afcfta/output/csv/comesa_border_results.csv"))

# ----- Robust rank-plot helper -----
make_rank_plot <- function(res, var_term, out_name, ylab,
                           sig_only = FALSE,
                           se_pctl_max = 0.95,
                           ciwidth_pctl_max = 0.95,
                           est_pctl_range = c(0.05, 0.95)) {

  df0 <- res %>%
    filter(term == var_term) %>%
    mutate(
      estimate  = as.numeric(estimate),
      conf.low  = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      std.error = as.numeric(std.error),
      ci_width  = conf.high - conf.low
    ) %>%
    filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high),
           is.finite(std.error), is.finite(ci_width))

  if (sig_only) df0 <- filter(df0, !is.na(p.value) & p.value <= 0.10)
  if (nrow(df0) == 0L) {
    warning(sprintf("No usable estimates for %s after filters.", var_term))
    return(invisible(NULL))
  }

  se_cut     <- as.numeric(quantile(df0$std.error, probs = se_pctl_max, na.rm = TRUE))
  ciw_cut    <- as.numeric(quantile(df0$ci_width,  probs = ciwidth_pctl_max, na.rm = TRUE))
  est_bounds <- as.numeric(quantile(df0$estimate,  probs = est_pctl_range,  na.rm = TRUE))

  plot_df <- df0 %>%
    filter(std.error <= se_cut, ci_width <= ciw_cut) %>%
    filter(estimate > est_bounds[1], estimate < est_bounds[2]) %>%
    arrange(estimate)

  if (nrow(plot_df) == 0L) {
    warning(sprintf("All %s estimates dropped by thresholds; relax filters.", var_term))
    return(invisible(NULL))
  }

  plot_trim <- plot_df %>%
    mutate(rank = row_number()) %>%
    select(rank, estimate, conf.low, conf.high)

  plot_long <- pivot_longer(plot_trim,
                            cols = c(estimate, conf.low, conf.high),
                            names_to = "series", values_to = "value") %>%
    mutate(series = factor(series,
                           levels = c("estimate","conf.low","conf.high"),
                           labels = c("Estimate","Lower 95 CL","Upper 95 CL")))

  p <- ggplot(plot_long, aes(x = rank, y = value, color = series, linetype = series, group = series)) +
    geom_line() +
    scale_color_manual(values = c("Estimate" = "blue",
                                  "Lower 95 CL" = "red",
                                  "Upper 95 CL" = "green")) +
    scale_linetype_manual(values = c("Estimate" = "solid",
                                     "Lower 95 CL" = "dashed",
                                     "Upper 95 CL" = "dashed")) +
    labs(x = "Rank of Estimate in Terms of Size",
         y = ylab,
         color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal")

  ggsave(filename = file.path(out_dir, out_name), plot = p, width = 8, height = 5)
}

# ----- Make the two plots -----
make_rank_plot(results, "comesa_to_af", "graph_comesa_to_af.pdf",
               ylab = "COMESA_to_AfCFTA Border Effect",
               sig_only = FALSE,
               se_pctl_max = 0.95,
               ciwidth_pctl_max = 0.95,
               est_pctl_range = c(0.05, 0.95))

make_rank_plot(results, "af_to_comesa", "graph_af_to_comesa.pdf",
               ylab = "AfCFTA_to_COMESA Border Effect",
               sig_only = FALSE,
               se_pctl_max = 0.95,
               ciwidth_pctl_max = 0.95,
               est_pctl_range = c(0.05, 0.95))


# ====================== End of COMESA analysis =======================