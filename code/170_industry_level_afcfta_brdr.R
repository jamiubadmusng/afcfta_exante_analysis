# Industry-level analysis for Aggregate AfCFTA border --------------------------

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

# --- Paths ---
base_dir  <- "C:/Users/muham/EGEI Dissertation/afcfta/input"
grav_file <- file.path(base_dir, "dgd_2_1.rds")   # <-- .rds gravity file

# --- AfCFTA members (excl. Eritrea) ---
afc_countries <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM",
  "COG","CIV","COD","DJI","EGY","GNQ","SWZ","ETH","GAB","GMB","GHA",
  "GIN","GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS",
  "MAR","MOZ","NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","SOM",
  "ZAF","SSD","SDN","TZA","TGO","TUN","UGA","ZMB","ZWE"
)

# --- Data prep function (merge + time-varying border vars) ---
prepare_data <- function(itpde_file) {
  
  # read ITPDE industry data (.rds), then rename to match gravity file
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)
  
  # read gravity covariates (.rds) with exporter/importer columns
  grav <- readRDS(grav_file)
  
  # merge
  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year"))
  
  # construct time-varying border variables (as provided)
  df <- df %>%
    mutate(
      # indicator dummies
      afcfta_brdr = ifelse(exporter %in% afc_countries &
                             importer %in% afc_countries &
                             exporter != importer, 1, 0),
      
      afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                                 !(importer %in% afc_countries) &
                                 exporter != importer, 1, 0),
      
      row_to_row = ifelse(!(exporter %in% afc_countries) &
                            !(importer %in% afc_countries) &
                            exporter != importer, 1, 0)
    ) %>%
    mutate(
      # make them time-varying via year-indexed factor codes (your original approach)
      afcfta_brdr     = as.integer(factor(paste(afcfta_brdr, year, sep = "_"))),
      afcfta_row_brdr = as.integer(factor(paste(afcfta_row_brdr, year, sep = "_"))),
      row_to_row      = as.integer(factor(paste(row_to_row, year, sep = "_"))),
      
      # ids for FE and clustering
      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )
  
  return(df)
}

# --- Estimation function (your exact model) ---
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + row_to_row +
      afcfta_brdr + afcfta_row_brdr | exp_year + imp_year,
    data = df,
    vcov = ~ pair_id
  )
  
  # grab only afcfta_brdr
  tidy(m, conf.int = TRUE) %>%
    filter(term == "afcfta_brdr") %>%
    mutate(industry_id = industry_id,
           industry_descr = industry_descr) %>%
    select(industry_id, industry_descr, term, estimate, std.error, conf.low, conf.high, p.value)
}

# --- Loop over 170 industries ---
itpde_files <- file.path(base_dir, sprintf("itpder2_%d.rds", 1:170))

results <- map_dfr(seq_along(itpde_files), function(i) {
  df <- prepare_data(itpde_files[i])
  # Get industry description from the first row (assuming it's the same for all rows in each file)
  industry_descr <- df$industry_descr[1]
  estimate_model(df, i, industry_descr)
})

# --- Save results (afcfta_brdr only) ---
out_csv <- file.path("C:/Users/muham/EGEI Dissertation/afcfta/output/csv/afcfta_brdr_industry_results.csv")
write_csv(results, out_csv)

# --- Rank plot for afcfta_brdr (trim 5–95th percentile like Stata) ---
# results must contain: estimate, conf.low, conf.high
# 1) clean + sort
plot_df <- results %>%
  mutate(across(c(estimate, conf.low, conf.high), as.numeric)) %>%
  filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
  arrange(estimate)

# 2) trim by 5th–95th pct and THEN compute contiguous ranks 1..N
qs <- quantile(plot_df$estimate, probs = c(0.05, 0.95), na.rm = TRUE)

plot_trim <- plot_df %>%
  filter(estimate > qs[1], estimate < qs[2]) %>%
  arrange(estimate) %>%
  mutate(rank = row_number()) %>%
  select(rank, estimate, conf.low, conf.high)

# 3) long form for legend; set factor order and labels
plot_long <- plot_trim %>%
  pivot_longer(cols = c(estimate, conf.low, conf.high),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("estimate", "conf.low", "conf.high"),
                         labels = c("Estimate", "Lower 95 CL", "Upper 95 CL")))

# 4) plot (estimate=blue solid; lower=red dashed; upper=green dashed)
p <- ggplot(plot_long, aes(x = rank, y = value, color = series, linetype = series)) +
  geom_line() +
  scale_color_manual(values = c("Estimate" = "blue",
                                "Lower 95 CL" = "red",
                                "Upper 95 CL" = "green")) +
  scale_linetype_manual(values = c("Estimate" = "solid",
                                   "Lower 95 CL" = "dashed",
                                   "Upper 95 CL" = "dashed")) +
  labs(x = "Rank of Estimate in Terms of Size",
       y = "AfCFTA Border Effect (PPML Estimate)",
       color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal")

# save
out_pdf <- file.path("C:/Users/muham/EGEI Dissertation/afcfta/output/figures/graph_afcfta_brdr.pdf")
ggsave(out_pdf, plot = p, width = 8, height = 5)

# ======== End of Script ================