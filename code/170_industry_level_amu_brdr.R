# ====================== Industry-level analysis for AMU =======================

# ---- clear environment -------
rm(list=ls())

# --- Packages ---
library(dplyr)     # data wrangling
library(fixest)    # fepois (PPML with FEs)
library(purrr)     # map / loops
library(broom)     # tidy(), conf.int
library(ggplot2)   # plotting
library(readr)     # write_csv
library(tidyr)     # pivot_longer


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
amu_countries <- c("DZA", "LBY", "MAR", "MRT", "TUN")

# ----- Data prep (merge + AMU border variables; preserves your time-varying encoding) -----
prepare_data <- function(itpde_file) {
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)
  
  grav <- readRDS(grav_file)
  
  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year")) %>%
    mutate(
      amu = ifelse(exporter %in% amu_countries &
                     importer %in% amu_countries &
                     exporter != importer, 1, 0),
      
      amu_to_af = ifelse(exporter %in% amu_countries &
                           importer %in% afc_countries &
                           !(importer %in% amu_countries) &
                           exporter != importer, 1, 0),
      
      af_to_amu = ifelse(importer %in% amu_countries &
                           exporter %in% afc_countries &
                           !(exporter %in% amu_countries) &
                           exporter != importer, 1, 0),
      
      af_n_amu = ifelse(exporter %in% afc_countries &
                          importer %in% afc_countries &
                          !(exporter %in% amu_countries) &
                          !(importer %in% amu_countries) &
                          exporter != importer, 1, 0),
      
      amu_to_row = ifelse(exporter %in% amu_countries &
                            !(importer %in% afc_countries) &
                            exporter != importer, 1, 0),
      
      af_n_amu_to_row = ifelse(exporter %in% afc_countries &
                                 !(exporter %in% amu_countries) &
                                 !(importer %in% afc_countries) &
                                 exporter != importer, 1, 0),
      
      rr_n_amu = ifelse(!(exporter %in% afc_countries) &
                          !(importer %in% afc_countries) &
                          exporter != importer, 1, 0)
    ) %>%
    mutate(
      amu            = as.integer(factor(paste(amu, year, sep = "_"))),
      amu_to_af      = as.integer(factor(paste(amu_to_af, year, sep = "_"))),
      af_to_amu      = as.integer(factor(paste(af_to_amu, year, sep = "_"))),
      af_n_amu       = as.integer(factor(paste(af_n_amu, year, sep = "_"))),
      amu_to_row     = as.integer(factor(paste(amu_to_row, year, sep = "_"))),
      af_n_amu_to_row= as.integer(factor(paste(af_n_amu_to_row, year, sep = "_"))),
      rr_n_amu       = as.integer(factor(paste(rr_n_amu, year, sep = "_"))),
      
      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )
  
  df
}

# ----- Estimation (AMU model) -----
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + rr_n_amu +
      amu + amu_to_af + af_to_amu + af_n_amu + af_n_amu_to_row + amu_to_row |
      exp_year + imp_year,
    data = df,
    vcov = ~ pair_id
  )
  
  broom::tidy(m, conf.int = TRUE) %>%
    filter(term %in% c("amu_to_af", "af_to_amu")) %>%
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
write_csv(results, file.path("C:/Users/muham/EGEI Dissertation/afcfta/output/csv/amu_border_results.csv"))

# ----- Helper to make rank-plot for one variable (robust & consistent) -----
make_rank_plot <- function(res, var_term, out_name, ylab,
                           sig_only = FALSE,
                           se_pctl_max = 0.95,         # keep SE <= 95th pct
                           ciwidth_pctl_max = 0.95,    # keep CI width <= 95th pct
                           est_pctl_range = c(0.05, 0.95)) {  # keep middle 90% by estimate
  
  df0 <- res %>%
    filter(term == var_term) %>%
    mutate(
      estimate = as.numeric(estimate),
      conf.low = as.numeric(conf.low),
      conf.high = as.numeric(conf.high),
      std.error = as.numeric(std.error),
      ci_width = conf.high - conf.low
    ) %>%
    filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high),
                  is.finite(std.error), is.finite(ci_width))
  
  if (sig_only) df0 <- filter(df0, !is.na(p.value) & p.value <= 0.10)
  if (nrow(df0) == 0L) {
    warning(sprintf("No usable estimates for %s after filters.", var_term))
    return(invisible(NULL))
  }
  
  # thresholds
  se_cut     <- as.numeric(quantile(df0$std.error, probs = se_pctl_max, na.rm = TRUE))
  ciw_cut    <- as.numeric(quantile(df0$ci_width,  probs = ciwidth_pctl_max, na.rm = TRUE))
  est_bounds <- as.numeric(quantile(df0$estimate, probs = est_pctl_range, na.rm = TRUE))
  
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
  
  plot_long <- pivot_longer(
    plot_trim, cols = c(estimate, conf.low, conf.high),
    names_to = "series", values_to = "value"
  ) %>%
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

# ----- Produce the two figures with your requested y-axis labels -----
make_rank_plot(results, "amu_to_af", "graph_amu_to_af.pdf",
               ylab = "AMU_to_AfCFTA Border Effect",
               sig_only = FALSE,          # set TRUE to keep only significant ones
               se_pctl_max = 0.95,
               ciwidth_pctl_max = 0.95,
               est_pctl_range = c(0.05, 0.95))

make_rank_plot(results, "af_to_amu", "graph_af_to_amu.pdf",
               ylab = "AfCFTA_to_AMU Border Effect",
               sig_only = FALSE,
               se_pctl_max = 0.95,
               ciwidth_pctl_max = 0.95,
               est_pctl_range = c(0.05, 0.95))


# ====================== End of AMU analysis =======================