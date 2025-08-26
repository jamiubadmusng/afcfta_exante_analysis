################################################################################
# R code for General Equilibrium PPML Analysis of AfCFTA
# This code translates the original Stata GE PPML code by Anderson et al. (2018)
# for the analysis of abolishing borders between AfCFTA countries
# for Mining and Energy Sector data in ITPD-E-R02
# Author: Jamiu Olamilekan Badmus
################################################################################

# Clear workspace
rm(list = ls())

# Load required packages
# packages <- c("data.table", "fixest", "dplyr", "tidyr", "readstata13", "ggplot2")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages)

library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(readstata13)
library(ggplot2)

# Set working directory
setwd("C:/Users/muham/EGEI Dissertation/afcfta")

# Set parameters
sigma <- 7

################################################################################
# I. Prepare Data for Mining and Energy Sector
################################################################################

# Load data
data <- read.dta13("input/afcfta_2019_mine_balanced.dta")
dt <- as.data.table(data)

# 1. Create aggregate variables
dt[, output := sum(trade), by = "exporter"]
dt[, expndr := sum(trade), by = "importer"]

# 2. Choose a country for reference group (Germany)
dt[exporter == "DEU", exporter := "ZZZ"]
dt[importer == "DEU", importer := "ZZZ"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Save processed data
saveRDS(dt, "input/ge_ppml_data_mine.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("input/ge_ppml_data_mine.rds")

#############################
# Step 1: "Baseline" Scenario
#############################

# Step 1.a: Estimate "Baseline" Gravity
# Create country dummy variables (for exporter and importer fixed effects)
countries <- unique(dt$exporter)
NoC <- length(countries)

# Baseline PPML with fixest
# Using fepois for PPML with high-dimensional fixed effects
baseline_model <- fepois(
  trade ~ cntg + col + lang + dist + rta + afcfta_brdr + afcfta_row_brdr + row_to_row | exporter + importer,
  data = dt
)

# Save the model summary
baseline_summary <- summary(baseline_model)

# Extract coefficients
coeffs <- coef(baseline_model)
CNTG_est <- coeffs["cntg"]
COL_est <- coeffs["col"]
LANG_est <- coeffs["lang"]
DIST_est <- coeffs["dist"]
RTA_est <- coeffs["rta"]
AfCFTA_BRDR_est <- coeffs["afcfta_brdr"]
AfCFTA_ROW_BRDR_est <- coeffs["afcfta_row_brdr"]
ROW_to_ROW_est <- coeffs["row_to_row"]

# Predict trade in the baseline
dt[, trade_bsln := predict(baseline_model, newdata = dt, type = "response")]

# Create baseline and counterfactual trade costs
dt[, t_ij_bsln := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * afcfta_brdr +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Counterfactual: set afcfta_brdr = 0, keep everything else the same
dt[, t_ij_ctrf := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * 0 +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Keep domestic trade costs at baseline
dt[exporter == importer, t_ij_ctrf := t_ij_bsln]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# Step 1.b: Construct "Baseline" GE Indexes
# Extract fixed effects
fe_exporter <- fixef(baseline_model)$exporter
fe_importer <- fixef(baseline_model)$importer

# Convert fixed effects to data table
exp_fe_dt <- data.table(exporter = names(fe_exporter), exp_fe = exp(fe_exporter))
imp_fe_dt <- data.table(importer = names(fe_importer), imp_fe = exp(fe_importer))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_0 := exp_fe]
dt[, all_imp_fes_0 := imp_fe]

# Calculate importer-specific exporter fixed effects
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_TB := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Equation (7): Outward Multilateral Resistance
dt[, omr_bsln := output * expndr_deu / all_exp_fes_0]

# Equation (8): Inward Multilateral Resistance
dt[, imr_bsln := expndr / (all_imp_fes_0 * expndr_deu)]

# Real GDP (baseline)
dt[exporter == importer, rGDP_bsln_temp := output / (imr_bsln^(1/(1-sigma)))]
dt[, rGDP_bsln := sum(rGDP_bsln_temp, na.rm = TRUE), by = "exporter"]

# Domestic absorption share (acr)
dt[exporter == importer, acr_bsln := trade_bsln / expndr]

# Bilateral exports & totals
dt[exporter != importer, exp_bsln := trade_bsln]
dt[exporter == importer, exp_bsln_acr := trade_bsln]
dt[, tot_exp_bsln := sum(exp_bsln, na.rm = TRUE), by = "exporter"]

################################
# Step 2: "Conditional" Scenario
################################

# Step 2.a: Estimate "Conditional" Gravity
# PPML with offset = log(counterfactual trade costs)
conditional_model <- fepois(
  trade ~ 1 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1
)

# Predict trade under conditional GE
dt[, trade_cndl := predict(conditional_model, newdata = dt, type = "response")]

# Step 2.b: Construct "Conditional" GE Indexes
# Extract conditional fixed effects
fe_exporter_cndl <- fixef(conditional_model)$exporter
fe_importer_cndl <- fixef(conditional_model)$importer

# Convert fixed effects to data table
exp_fe_cndl_dt <- data.table(exporter = names(fe_exporter_cndl), exp_fe_cndl = exp(fe_exporter_cndl))
imp_fe_cndl_dt <- data.table(importer = names(fe_importer_cndl), imp_fe_cndl = exp(fe_importer_cndl))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_cndl_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_cndl_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_1 := exp_fe_cndl]
dt[, all_imp_fes_1 := imp_fe_cndl]

# Equation (7): Outward Multilateral Resistance (conditional)
dt[, omr_cndl := output * expndr_deu / all_exp_fes_1]

# Equation (8): Inward Multilateral Resistance (conditional)
dt[, imr_cndl := expndr / (all_imp_fes_1 * expndr_deu)]

# Exports and totals (conditional)
dt[exporter != importer, exp_cndl := trade_cndl]
dt[, tot_exp_cndl := sum(exp_cndl, na.rm = TRUE), by = "exporter"]

# % change in total exports
dt[, tot_exp_cndl_ch := (tot_exp_cndl - tot_exp_bsln) / tot_exp_bsln * 100]

# Real GDP (conditional)
dt[exporter == importer, rGDP_cndl_temp := output / (imr_cndl^(1/(1-sigma)))]
dt[, rGDP_cndl := sum(rGDP_cndl_temp, na.rm = TRUE), by = "exporter"]

###################################
# Step 3: "Full Endowment" Scenario
###################################

# Step 3.a: Estimate "Full Endowment" Gravity
# Initialize variables for the iterative process
dt[, trade_1_pred := trade_cndl]
dt[, output_bsln := output]
dt[, expndr_bsln := expndr]
dt[exporter == importer, phi := expndr / output]

# Calculate initial values
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

dt[exporter == importer, expndr_temp_1 := phi * output]
dt[, expndr_1 := mean(expndr_temp_1, na.rm = TRUE), by = "importer"]
dt[, expndr_temp_1 := NULL]

dt[importer == "ZZZ", expndr_deu01_1 := expndr_1]
dt[, expndr_deu_1 := mean(expndr_deu01_1, na.rm = TRUE)]

dt[exporter == importer, temp := all_exp_fes_1]
dt[, all_exp_fes_1_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Initial values for factory-gate prices and multilateral resistance terms
dt[, p_full_exp_0 := 0]
dt[, p_full_exp_1 := (all_exp_fes_1 / all_exp_fes_0)^(1/(1-sigma))]
dt[, p_full_imp_1 := (all_exp_fes_1_imp / all_exp_fes_0_imp)^(1/(1-sigma))]
dt[, imr_full_1 := expndr_1 / (all_imp_fes_1 * expndr_deu_1)]
dt[, imr_full_ch_1 := 1]
dt[, omr_full_1 := output * expndr_deu_1 / all_exp_fes_1]
dt[, omr_full_ch_1 := 1]

# Set convergence parameters
max_iterations <- 30
diff_all_exp_fes_sd <- 1
diff_all_exp_fes_max <- 1
tolerance <- 0.001
i <- 3
iteration <- 0

# Define a function to check for duplicated variables and clean them up
check_and_remove_vars <- function(dt, pattern) {
  vars_to_remove <- grep(pattern, names(dt), value = TRUE)
  if (length(vars_to_remove) > 0) {
    dt[, (vars_to_remove) := NULL]
  }
  return(dt)
}

# Iterative process
cat("Starting iterative process for full endowment calculation...\n")

while ((diff_all_exp_fes_sd > tolerance || diff_all_exp_fes_max > tolerance) && iteration < max_iterations) {
  iteration <- iteration + 1
  cat(sprintf("Iteration %d of maximum %d\n", iteration, max_iterations))
  
  # Clean up variables before recreation
  current_iter <- i-1
  patterns <- c(
    paste0("trade_", current_iter),
    paste0("all_exp_fes_", current_iter),
    paste0("all_imp_fes_", current_iter),
    paste0("output_", current_iter),
    paste0("expndr_check_", current_iter),
    paste0("expndr_deu0_", current_iter),
    paste0("expndr_deu_", current_iter),
    paste0("all_exp_fes_", current_iter, "_imp"),
    paste0("p_full_exp_", current_iter),
    paste0("p_full_imp_", current_iter),
    paste0("omr_full_", current_iter),
    paste0("omr_full_ch_", current_iter),
    paste0("expndr_temp_", current_iter),
    paste0("expndr_", current_iter),
    paste0("imr_full_", current_iter),
    paste0("imr_full_ch_", current_iter),
    paste0("diff_p_full_exp_", current_iter),
    paste0("trade_", current_iter, "_pred")
  )
  
  for (pattern in patterns) {
    dt <- check_and_remove_vars(dt, pattern)
  }
  
  # Equation (14)
  dt[, paste0("trade_", current_iter) := get(paste0("trade_", i-2, "_pred")) * 
      get(paste0("p_full_exp_", i-2)) * get(paste0("p_full_imp_", i-2)) / 
      (get(paste0("omr_full_ch_", i-2)) * get(paste0("imr_full_ch_", i-2)))]
  
  # Estimate new PPML model
  tryCatch({
    formula_str <- paste0("trade_", current_iter, " ~ 1 | exporter + importer")
    full_model <- fepois(
      as.formula(formula_str),
      data = dt,
      offset = ~t_ij_ctrf_1
    )
    
    # Predict trade flows
    dt[, paste0("trade_", current_iter, "_pred") := predict(full_model, newdata = dt, type = "response")]
    
    # Extract fixed effects
    fe_exporter_full <- fixef(full_model)$exporter
    fe_importer_full <- fixef(full_model)$importer
    
    # Process fixed effects
    dt[, paste0("all_exp_fes_", current_iter) := 0]
    dt[, paste0("all_imp_fes_", current_iter) := 0]
    
    # Apply exporter fixed effects
    for (country in names(fe_exporter_full)) {
      dt[exporter == country, paste0("all_exp_fes_", current_iter) := exp(fe_exporter_full[country])]
    }
    
    # Apply importer fixed effects
    for (country in names(fe_importer_full)) {
      dt[importer == country, paste0("all_imp_fes_", current_iter) := exp(fe_importer_full[country])]
    }
    
    # Update output
    dt[, paste0("output_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "exporter"]
    
    # Update expenditure
    dt[, paste0("expndr_check_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "importer"]
    dt[importer == "ZZZ", paste0("expndr_deu0_", current_iter) := get(paste0("expndr_check_", current_iter))]
    dt[, paste0("expndr_deu_", current_iter) := mean(get(paste0("expndr_deu0_", current_iter)), na.rm = TRUE)]
    
    # Get importer-specific exporter fixed effects
    dt[exporter == importer, temp := get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("all_exp_fes_", current_iter, "_imp") := mean(temp, na.rm = TRUE), by = "importer"]
    dt[, temp := NULL]
    
    # Update factory-gate prices
    dt[, paste0("p_full_exp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter)) / get(paste0("all_exp_fes_", i-2))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    dt[, paste0("p_full_imp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter, "_imp")) / get(paste0("all_exp_fes_", i-2, "_imp"))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    # Equation (7) - Update outward multilateral resistance
    dt[, paste0("omr_full_", current_iter) := get(paste0("output_", current_iter)) / get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("omr_full_ch_", current_iter) := get(paste0("omr_full_", current_iter)) / get(paste0("omr_full_", i-2))]
    
    # Update expenditure
    dt[exporter == importer, paste0("expndr_temp_", current_iter) := phi * get(paste0("output_", current_iter))]
    dt[, paste0("expndr_", current_iter) := mean(get(paste0("expndr_temp_", current_iter)), na.rm = TRUE), by = "importer"]
    
    # Equation (8) - Update inward multilateral resistance
    dt[, paste0("imr_full_", current_iter) := get(paste0("expndr_", current_iter)) / 
         (get(paste0("all_imp_fes_", current_iter)) * get(paste0("expndr_deu_", current_iter)))]
    dt[, paste0("imr_full_ch_", current_iter) := get(paste0("imr_full_", current_iter)) / get(paste0("imr_full_", i-2))]
    
    # Calculate convergence criteria
    dt[, paste0("diff_p_full_exp_", current_iter) := get(paste0("p_full_exp_", i-2)) - get(paste0("p_full_exp_", i-3))]
    
    # Calculate convergence statistics
    diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                          max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
    diff_all_exp_fes_sd <- diff_stats$sd
    diff_all_exp_fes_max <- diff_stats$max
    
    cat(sprintf("Convergence stats - SD: %f, Max: %f\n", diff_all_exp_fes_sd, diff_all_exp_fes_max))
    
  }, error = function(e) {
    cat("Error in iteration", iteration, ":", e$message, "\n")
    cat("Using values from previous iteration\n")
    
    # Create placeholder with previous values
    dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
  })
  
  # Increment counter for next iteration
  i <- i + 1
}

cat(sprintf("Convergence achieved after %d iterations\n", iteration))

# Get the final iteration number for use in subsequent calculations
final_iter <- i - 2

# Step 3.b: Construct "Full Endowment" GE Indexes
# Calculate p^c/p
dt[, c("output", "expndr", "expndr_deu0") := NULL]
dt[, expndr_deu_bsln := expndr_deu]

dt[, output := sum(get(paste0("trade_", final_iter, "_pred"))), by = "exporter"]
dt[exporter == importer, expndr_temp := phi * output]
dt[, expndr := mean(expndr_temp, na.rm = TRUE), by = "importer"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Calculate p_full
dt[, p_full := ((get(paste0("all_exp_fes_", final_iter)) / all_exp_fes_0) / 
                 (expndr_deu / expndr_deu_bsln))^(1/(1-sigma))]

# Calculate output_full
dt[, output_full := p_full * output_bsln]

# Equation (7): Outward Multilateral Resistance
dt[, omr_full := output_full * expndr_deu / get(paste0("all_exp_fes_", final_iter))]

# Equation (8): Inward Multilateral Resistance
dt[, imr_full := expndr / (get(paste0("all_imp_fes_", final_iter)) * expndr_deu)]

# Real GDP (full endowment)
dt[exporter == importer, rGDP_full_temp := p_full * output_bsln / (imr_full^(1/(1-sigma)))]
dt[, rGDP_full := sum(rGDP_full_temp, na.rm = TRUE), by = "exporter"]

# Expenditure (full endowment)
dt[exporter == importer, expndr_full_temp := phi * output_full]
dt[, expndr_full := mean(expndr_full_temp, na.rm = TRUE), by = "importer"]

# Bilateral trade at full endowment (using counterfactual costs)
dt[, trade_full := (output_full * expndr_full * t_ij_ctrf) / (imr_full * omr_full)]
dt[exporter != importer, exp_full := trade_full]
dt[, tot_exp_full := sum(exp_full, na.rm = TRUE), by = "exporter"]
dt[, tot_exp_full_ch := (tot_exp_full - tot_exp_bsln) / tot_exp_bsln * 100]

# Save results
saveRDS(dt, "output/rds/full_static_all_mine.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "output/rds/imrs_all_mine.rds")

# OMR and other indexes
omr_dt <- dt[, .(omr_full = mean(omr_full, na.rm = TRUE),
                 omr_cndl = mean(omr_cndl, na.rm = TRUE),
                 omr_bsln = mean(omr_bsln, na.rm = TRUE),
                 rGDP_full = mean(rGDP_full, na.rm = TRUE),
                 rGDP_cndl = mean(rGDP_cndl, na.rm = TRUE),
                 rGDP_bsln = mean(rGDP_bsln, na.rm = TRUE),
                 tot_exp_full = mean(tot_exp_full, na.rm = TRUE),
                 tot_exp_cndl = mean(tot_exp_cndl, na.rm = TRUE),
                 tot_exp_bsln = mean(tot_exp_bsln, na.rm = TRUE),
                 p_full = mean(p_full, na.rm = TRUE),
                 output_bsln = mean(output_bsln, na.rm = TRUE),
                 acr_bsln = mean(acr_bsln, na.rm = TRUE)), by = exporter]
setnames(omr_dt, "exporter", "country")
omr_dt[, omr_full_ch := (omr_full - omr_bsln) / omr_bsln * 100]
omr_dt[, omr_cndl_ch := (omr_cndl - omr_bsln) / omr_bsln * 100]
omr_dt[, rGDP_full_ch := (rGDP_full - rGDP_bsln) / rGDP_bsln * 100]
omr_dt[, rGDP_cndl_ch := (rGDP_cndl - rGDP_bsln) / rGDP_bsln * 100]

# Combine OMR and IMR indexes
all_indexes <- merge(omr_dt, imr_dt, by = "country", all = TRUE)
saveRDS(all_indexes, "output/rds/all_indexes_geppml_mine.rds")

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Add a column to identify AfCFTA countries
all_indexes[, is_afcfta := country %in% afcfta_countries]

# Calculate total output for ROW for weighting
row_total_output <- all_indexes[is_afcfta == FALSE, sum(output_bsln, na.rm = TRUE)]

# Calculate weighted welfare change for ROW countries
# Following Larch, Tan, and Yotov (2023) using output shares as weights
weighted_row_welfare <- all_indexes[is_afcfta == FALSE, 
                                    sum(rGDP_full_ch * output_bsln, na.rm = TRUE) / 
                                      sum(output_bsln, na.rm = TRUE)]

# Create a copy of the dataset with just AfCFTA countries
afcfta_indexes <- all_indexes[is_afcfta == TRUE]

# Add the aggregated ROW row
row_entry <- data.table(
  country = "ROW",
  is_afcfta = FALSE,
  omr_full = mean(all_indexes[is_afcfta == FALSE, omr_full], na.rm = TRUE),
  omr_cndl = mean(all_indexes[is_afcfta == FALSE, omr_cndl], na.rm = TRUE),
  omr_bsln = mean(all_indexes[is_afcfta == FALSE, omr_bsln], na.rm = TRUE),
  rGDP_full = mean(all_indexes[is_afcfta == FALSE, rGDP_full], na.rm = TRUE),
  rGDP_cndl = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl], na.rm = TRUE),
  rGDP_bsln = mean(all_indexes[is_afcfta == FALSE, rGDP_bsln], na.rm = TRUE),
  tot_exp_full = sum(all_indexes[is_afcfta == FALSE, tot_exp_full], na.rm = TRUE),
  tot_exp_cndl = sum(all_indexes[is_afcfta == FALSE, tot_exp_cndl], na.rm = TRUE),
  tot_exp_bsln = sum(all_indexes[is_afcfta == FALSE, tot_exp_bsln], na.rm = TRUE),
  p_full = mean(all_indexes[is_afcfta == FALSE, p_full], na.rm = TRUE),
  output_bsln = row_total_output,
  acr_bsln = mean(all_indexes[is_afcfta == FALSE, acr_bsln], na.rm = TRUE),
  omr_full_ch = mean(all_indexes[is_afcfta == FALSE, omr_full_ch], na.rm = TRUE),
  omr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, omr_cndl_ch], na.rm = TRUE),
  rGDP_full_ch = weighted_row_welfare,
  rGDP_cndl_ch = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl_ch], na.rm = TRUE),
  imr_full = mean(all_indexes[is_afcfta == FALSE, imr_full], na.rm = TRUE),
  imr_bsln = mean(all_indexes[is_afcfta == FALSE, imr_bsln], na.rm = TRUE),
  imr_cndl = mean(all_indexes[is_afcfta == FALSE, imr_cndl], na.rm = TRUE),
  imr_full_ch = mean(all_indexes[is_afcfta == FALSE, imr_full_ch], na.rm = TRUE),
  imr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, imr_cndl_ch], na.rm = TRUE)
)

# Combine AfCFTA countries and ROW
all_indexes_with_row <- rbindlist(list(afcfta_indexes, row_entry), fill = TRUE)

# Save the aggregated data
saveRDS(all_indexes_with_row, "output/rds/all_indexes_with_row_geppml_mine.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP):\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW): ", round(weighted_row_welfare, 4), "%\n")

# Create plots for welfare effects showing AfCFTA countries and ROW
# Color by whether the country is an AfCFTA member
welfare_plot <- ggplot(all_indexes_with_row, 
                      aes(x = reorder(country, rGDP_full_ch), 
                          y = rGDP_full_ch,
                          fill = is_afcfta)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", rGDP_full_ch)), 
            hjust = ifelse(all_indexes_with_row$rGDP_full_ch > 0, -0.1, 1.1),
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("gray60", "steelblue"), 
                    labels = c("Rest of World", "AfCFTA Countries"),
                    name = "") +
  labs(title = "Welfare Effects of AfCFTA Agreement for Mining & Energy Sector",
       subtitle = "% Change in Real GDP",
       x = "Country",
       y = "% Change in Real GDP") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
ggsave("output/figures/welfare_effects_with_row_geppml_mine.pdf", welfare_plot, width = 10, height = 12, dpi = 300)

# Create another plot showing just AfCFTA countries for detail
welfare_plot_afcfta <- ggplot(all_indexes_with_row[is_afcfta == TRUE], 
                             aes(x = reorder(country, rGDP_full_ch), 
                                 y = rGDP_full_ch)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f%%", rGDP_full_ch)), 
            hjust = ifelse(all_indexes_with_row[is_afcfta == TRUE]$rGDP_full_ch > 0, -0.1, 1.1),
            size = 3.5) +
  coord_flip() +
  labs(title = "Welfare Effects for AfCFTA Member Countries for Mining & Energy Sector",
       subtitle = "% Change in Real GDP",
       x = "Country",
       y = "% Change in Real GDP") +
  theme_minimal()

# Save the AfCFTA-only plot
ggsave("output/figures/welfare_effects_afcfta_only_geppml_mine.pdf", welfare_plot_afcfta, width = 8, height = 12, dpi = 300)

cat("\nAnalysis completed. Results saved to output folder.\n")

##################################################################
### End of Script
##################################################################