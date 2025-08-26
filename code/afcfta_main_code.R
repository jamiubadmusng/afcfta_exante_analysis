################################################################################
####### R Code for the main analysis of the ex-ante analysis of AfCFTA #########
###################### by Jamiu Olamilekan Badmus ##############################
################################################################################

# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Set the working directory ----------------------------------------------------
setwd("C:/Users/muham/EGEI Dissertation/afcfta")


# Load the necessary libraries -------------------------------------------------
library(tidyverse)
library(here)
library(fixest)
library(flextable)
library(huxtable)
library(readr)
library(msm)
library(car)


# Load the datasets ------------------------------------------------------------
itpder2 <- readRDS("input/itpder2.rds") # Trade Dataset
dgd_2_1 <- readRDS("input/dgd_2_1.rds") # DGD Dataset

# Merge datasets
afcfta_data <- left_join(
  itpder2,
  dgd_2_1 %>%
    select(year, exporter, importer, cntg, col, lang, dist, rta),
  by = c("year", "exporter", "importer")
)


# construct exporter-time, importer-time, and asymmetric pair fixed effects
afcfta_data <- afcfta_data %>%
  mutate(
    exp_year = as.integer(factor(paste(exporter, year, sep = "_"))),
    imp_year = as.integer(factor(paste(importer, year, sep = "_"))),
    pair_id = as.integer(factor(paste(exporter, importer, sep = "_")))
  )


################################################################################
### Research Objective I for AfCFTA --------------------------------------------
################################################################################

# AfCFTA members excluding Eritrea
afc_countries <- c("DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
                   "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
                   "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
                   "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
                   "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE")


# create time-varying border variables between AfCFTA countries and ROW
afcfta_data <- afcfta_data %>%
  mutate(
    # indicator for intra-AfCFTA international trade
    afcfta_brdr = ifelse(exporter %in% afc_countries &
                           importer %in% afc_countries &
                           exporter != importer, 1, 0),
    
    # indicator for AfCFTA exports to the ROW
    afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    # indicator for ROW exports to ROW
    row_to_row = ifelse(!(exporter %in% afc_countries) &
                          !(importer %in% afc_countries) &
                          exporter != importer, 1, 0)
  ) %>%
  mutate(
    afcfta_brdr = as.integer(factor(paste(afcfta_brdr, year, sep = "_"))),
    afcfta_row_brdr = as.integer(factor(paste(afcfta_row_brdr, year, sep = "_"))),
    row_to_row = as.integer(factor(paste(row_to_row, year, sep = "_")))
  )


## Agriculture trade
agri_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_obj)

# Delta method for Agriculture
expr_afcfta <- paste0("(exp(afcfta_brdr) - 1) * 100")
delta_afcfta_agri <- deltaMethod(object = coef(agri_obj), vcov. = vcov(agri_obj), g = expr_afcfta, parameterNames = names(coef(agri_obj)))
print(delta_afcfta_agri)


## Manufacturing trade
manu_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_obj)

# Delta method for Manufacturing
delta_afcfta_manu <- deltaMethod(object = coef(manu_obj), vcov. = vcov(manu_obj), g = expr_afcfta, parameterNames = names(coef(manu_obj)))
print(delta_afcfta_manu)

## Mining and Energy Trade
mine_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_obj)

# Delta method for Mining and Energy
delta_afcfta_mine <- deltaMethod(object = coef(mine_obj), vcov. = vcov(mine_obj), g = expr_afcfta, parameterNames = names(coef(mine_obj)))
print(delta_afcfta_mine)

# Service trade
serv_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_obj)

# Delta method for Services
delta_afcfta_serv <- deltaMethod(object = coef(serv_obj), vcov. = vcov(serv_obj), g = expr_afcfta, parameterNames = names(coef(serv_obj)))
print(delta_afcfta_serv)

# AfCFTA Latex results ---------------------------------------------------------
tab_afcfta <- huxreg("Agri" = agri_obj,
                     "Manu" = manu_obj,
                     "Mine" = mine_obj,
                     "Serv" = serv_obj,
                     coefs = c("CNTG" = "cntg",
                               "COLN" = "col",
                               "DIST" = "dist",
                               "LANG" = "lang",
                               "RTA" = "rta",
                               "INTL_BRDR" = "row_to_row",
                               "BRDR_AFC" = "afcfta_brdr",
                               "BRDR_AFC_ROW" = "afcfta_row_brdr"),
                     stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                     note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on intra-African Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on Sectoral Intra-African Trade") %>%
  set_label("tab_afcfta")
width(tab_afcfta) <- 1.1
cat(to_latex(tab_afcfta), file = "output/tables/tab_afcfta.tex")



################################################################################
### Research Objective II ------------------------------------------------------
################################################################################

# Arab Maghreb Union (AMU) -----------------------------------------------------
amu_countries <- c("DZA", "LBY", "MAR", "MRT", "TUN")

# constructing time-varying border variables
afcfta_data <- afcfta_data %>%
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
    amu = as.integer(factor(paste(amu, year, sep = "_"))),
    amu_to_af = as.integer(factor(paste(amu_to_af, year, sep = "_"))),
    af_to_amu = as.integer(factor(paste(af_to_amu, year, sep = "_"))),
    af_n_amu = as.integer(factor(paste(af_n_amu, year, sep = "_"))),
    amu_to_row = as.integer(factor(paste(amu_to_row, year, sep = "_"))),
    af_n_amu_to_row = as.integer(factor(paste(af_n_amu_to_row, year, sep = "_"))),
    rr_n_amu = as.integer(factor(paste(rr_n_amu, year, sep = "_")))
  )

# Estimation for AMU
# Agriculture trade
agri_amu = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_amu
                  + amu + amu_to_af + af_to_amu 
                  + af_n_amu + af_n_amu_to_row + amu_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_amu)

# Delta method for Agriculture
expr_amu <- paste0("(exp(amu) - 1) * 100")
expr_amu_to_af <- paste0("(exp(amu_to_af) - 1) * 100")
expr_af_to_amu <- paste0("(exp(af_to_amu) - 1) * 100")
delta_amu_agri <- deltaMethod(object = coef(agri_amu), vcov. = vcov(agri_amu), g = expr_amu, parameterNames = names(coef(agri_amu)))
delta_amu_to_af_agri <- deltaMethod(object = coef(agri_amu), vcov. = vcov(agri_amu), g = expr_amu_to_af, parameterNames = names(coef(agri_amu)))
delta_af_to_amu_agri <- deltaMethod(object = coef(agri_amu), vcov. = vcov(agri_amu), g = expr_af_to_amu, parameterNames = names(coef(agri_amu)))
print(delta_amu_agri)
print(delta_amu_to_af_agri)
print(delta_af_to_amu_agri)

# Manufacturing trade
manu_amu = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_amu
                  + amu + amu_to_af + af_to_amu 
                  + af_n_amu + af_n_amu_to_row + amu_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_amu)

# Delta method for Manufacturing
delta_amu_manu <- deltaMethod(object = coef(manu_amu), vcov. = vcov(manu_amu), g = expr_amu, parameterNames = names(coef(manu_amu)))
delta_amu_to_af_manu <- deltaMethod(object = coef(manu_amu), vcov. = vcov(manu_amu), g = expr_amu_to_af, parameterNames = names(coef(manu_amu)))
delta_af_to_amu_manu <- deltaMethod(object = coef(manu_amu), vcov. = vcov(manu_amu), g = expr_af_to_amu, parameterNames = names(coef(manu_amu)))
print(delta_amu_manu)
print(delta_amu_to_af_manu)
print(delta_af_to_amu_manu)

# Mining and Energy Trade
mine_amu = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_amu
                  + amu + amu_to_af + af_to_amu 
                  + af_n_amu + af_n_amu_to_row + amu_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_amu)

# Delta method for Mining and Energy
delta_amu_mine <- deltaMethod(object = coef(mine_amu), vcov. = vcov(mine_amu), g = expr_amu, parameterNames = names(coef(mine_amu)))
delta_amu_to_af_mine <- deltaMethod(object = coef(mine_amu), vcov. = vcov(mine_amu), g = expr_amu_to_af, parameterNames = names(coef(mine_amu)))
delta_af_to_amu_mine <- deltaMethod(object = coef(mine_amu), vcov. = vcov(mine_amu), g = expr_af_to_amu, parameterNames = names(coef(mine_amu)))
print(delta_amu_mine)
print(delta_amu_to_af_mine)
print(delta_af_to_amu_mine)

# Service trade
serv_amu = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_amu
                  + amu + amu_to_af + af_to_amu 
                  + af_n_amu + af_n_amu_to_row + amu_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_amu)

# Delta method for Services
delta_amu_serv <- deltaMethod(object = coef(serv_amu), vcov. = vcov(serv_amu), g = expr_amu, parameterNames = names(coef(serv_amu)))
#delta_amu_to_af_serv <- deltaMethod(object = coef(serv_amu), vcov. = vcov(serv_amu), g = expr_amu_to_af, parameterNames = names(coef(serv_amu)))
#delta_af_to_amu_serv <- deltaMethod(object = coef(serv_amu), vcov. = vcov(serv_amu), g = expr_af_to_amu, parameterNames = names(coef(serv_amu)))
print(delta_amu_serv)
#print(delta_amu_to_af_serv)
#print(delta_af_to_amu_serv)

# AMU Latex results
tab_amu <- huxreg("Agri" = agri_amu,
                  "Manu" = manu_amu,
                  "Mine" = mine_amu,
                  "Serv" = serv_amu,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_amu",
                            "BRDR_AMU" = "amu",
                            "BRDR_AMU_AFC" = "amu_to_af",
                            "BRDR_AFC_AMU" = "af_to_amu",
                            "BRDR_AFC" = "af_n_amu",
                            "BRDR_AFC_ROW" = "af_n_amu_to_row",
                            "BRDR_AMU_ROW" = "amu_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Arab Maghreb Union (AMU) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on AMU Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_amu")
width(tab_amu) <- 1.1
cat(to_latex(tab_amu), file = "output/tables/tab_amu.tex")




# Community of Sahel–Saharan States (CEN-SAD) — excluding Eritrea --------------
cen_countries <- c("BEN", "BFA", "CAF", "TCD", "COM", "CIV", "DJI", "EGY", "GMB", "GHA", "GIN",
                   "GNB", "KEN", "LBY", "MLI", "MRT", "MAR", "NER", "NGA", "SEN", "SLE", "SDN",
                   "SOM", "TGO", "TUN")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    cen = ifelse(exporter %in% cen_countries &
                   importer %in% cen_countries &
                   exporter != importer, 1, 0),
    
    cen_to_af = ifelse(exporter %in% cen_countries &
                         importer %in% afc_countries &
                         !(importer %in% cen_countries) &
                         exporter != importer, 1, 0),
    
    af_to_cen = ifelse(importer %in% cen_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% cen_countries) &
                         exporter != importer, 1, 0),
    
    af_n_cen = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% cen_countries) &
                        !(importer %in% cen_countries) &
                        exporter != importer, 1, 0),
    
    cen_to_row = ifelse(exporter %in% cen_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_cen_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% cen_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_cen = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    cen = as.integer(factor(paste(cen, year, sep = "_"))),
    cen_to_af = as.integer(factor(paste(cen_to_af, year, sep = "_"))),
    af_to_cen = as.integer(factor(paste(af_to_cen, year, sep = "_"))),
    af_n_cen = as.integer(factor(paste(af_n_cen, year, sep = "_"))),
    cen_to_row = as.integer(factor(paste(cen_to_row, year, sep = "_"))),
    af_n_cen_to_row = as.integer(factor(paste(af_n_cen_to_row, year, sep = "_"))),
    rr_n_cen = as.integer(factor(paste(rr_n_cen, year, sep = "_")))
  )

# Estimation for CEN-SAD
# Agriculture trade
agri_cen = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_cen
                  + cen + cen_to_af + af_to_cen 
                  + af_n_cen + af_n_cen_to_row + cen_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_cen)

# Delta method for Agriculture
expr_cen <- paste0("(exp(cen) - 1) * 100")
expr_cen_to_af <- paste0("(exp(cen_to_af) - 1) * 100")
expr_af_to_cen <- paste0("(exp(af_to_cen) - 1) * 100")
delta_cen_agri <- deltaMethod(object = coef(agri_cen), vcov. = vcov(agri_cen), g = expr_cen, parameterNames = names(coef(agri_cen)))
delta_cen_to_af_agri <- deltaMethod(object = coef(agri_cen), vcov. = vcov(agri_cen), g = expr_cen_to_af, parameterNames = names(coef(agri_cen)))
delta_af_to_cen_agri <- deltaMethod(object = coef(agri_cen), vcov. = vcov(agri_cen), g = expr_af_to_cen, parameterNames = names(coef(agri_cen)))
print(delta_cen_agri)
print(delta_cen_to_af_agri)
print(delta_af_to_cen_agri)

# Manufacturing trade
manu_cen = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_cen
                  + cen + cen_to_af + af_to_cen 
                  + af_n_cen + af_n_cen_to_row + cen_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_cen)

# Delta method for Manufacturing
delta_cen_manu <- deltaMethod(object = coef(manu_cen), vcov. = vcov(manu_cen), g = expr_cen, parameterNames = names(coef(manu_cen)))
delta_cen_to_af_manu <- deltaMethod(object = coef(manu_cen), vcov. = vcov(manu_cen), g = expr_cen_to_af, parameterNames = names(coef(manu_cen)))
delta_af_to_cen_manu <- deltaMethod(object = coef(manu_cen), vcov. = vcov(manu_cen), g = expr_af_to_cen, parameterNames = names(coef(manu_cen)))
print(delta_cen_manu)
print(delta_cen_to_af_manu)
print(delta_af_to_cen_manu)

# Mining and Energy Trade
mine_cen = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_cen
                  + cen + cen_to_af + af_to_cen 
                  + af_n_cen + af_n_cen_to_row + cen_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_cen)

# Delta method for Mining and Energy
delta_cen_mine <- deltaMethod(object = coef(mine_cen), vcov. = vcov(mine_cen), g = expr_cen, parameterNames = names(coef(mine_cen)))
delta_cen_to_af_mine <- deltaMethod(object = coef(mine_cen), vcov. = vcov(mine_cen), g = expr_cen_to_af, parameterNames = names(coef(mine_cen)))
delta_af_to_cen_mine <- deltaMethod(object = coef(mine_cen), vcov. = vcov(mine_cen), g = expr_af_to_cen, parameterNames = names(coef(mine_cen)))
print(delta_cen_mine)
print(delta_cen_to_af_mine)
print(delta_af_to_cen_mine)

# Service trade
serv_cen = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_cen
                  + cen + cen_to_af + af_to_cen 
                  + af_n_cen + af_n_cen_to_row + cen_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_cen)

# Delta method for Services
delta_cen_serv <- deltaMethod(object = coef(serv_cen), vcov. = vcov(serv_cen), g = expr_cen, parameterNames = names(coef(serv_cen)))
delta_cen_to_af_serv <- deltaMethod(object = coef(serv_cen), vcov. = vcov(serv_cen), g = expr_cen_to_af, parameterNames = names(coef(serv_cen)))
delta_af_to_cen_serv <- deltaMethod(object = coef(serv_cen), vcov. = vcov(serv_cen), g = expr_af_to_cen, parameterNames = names(coef(serv_cen)))
print(delta_cen_serv)
print(delta_cen_to_af_serv)
print(delta_af_to_cen_serv)

# CEN-SAD Latex results
tab_cen <- huxreg("Agri" = agri_cen,
                  "Manu" = manu_cen,
                  "Mine" = mine_cen,
                  "Serv" = serv_cen,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_cen",
                            "BRDR_CEN" = "cen",
                            "BRDR_CEN_AFC" = "cen_to_af",
                            "BRDR_AFC_CEN" = "af_to_cen",
                            "BRDR_AFC" = "af_n_cen",
                            "BRDR_AFC_ROW" = "af_n_cen_to_row",
                            "BRDR_CEN_ROW" = "cen_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Community of Sahel-Saharan States (CEN-SAD) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on CEN-SAD Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_cen")
width(tab_cen) <- 1.1
cat(to_latex(tab_cen), file = "output/tables/tab_cen.tex")



# Common Market for Eastern and Southern Africa (COMESA) — excluding Eritrea ---
com_countries <- c("BDI", "COM", "COD", "DJI", "EGY", "SWZ", "ETH", "KEN", "LYB", "MDG", "MWI",
                   "MUS", "RWA", "SYC", "SOM", "SDN", "TUN", "UGA", "ZMB", "ZWE")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    com = ifelse(exporter %in% com_countries &
                   importer %in% com_countries &
                   exporter != importer, 1, 0),
    
    com_to_af = ifelse(exporter %in% com_countries &
                         importer %in% afc_countries &
                         !(importer %in% com_countries) &
                         exporter != importer, 1, 0),
    
    af_to_com = ifelse(importer %in% com_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% com_countries) &
                         exporter != importer, 1, 0),
    
    af_n_com = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% com_countries) &
                        !(importer %in% com_countries) &
                        exporter != importer, 1, 0),
    
    com_to_row = ifelse(exporter %in% com_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_com_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% com_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_com = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    com = as.integer(factor(paste(com, year, sep = "_"))),
    com_to_af = as.integer(factor(paste(com_to_af, year, sep = "_"))),
    af_to_com = as.integer(factor(paste(af_to_com, year, sep = "_"))),
    af_n_com = as.integer(factor(paste(af_n_com, year, sep = "_"))),
    com_to_row = as.integer(factor(paste(com_to_row, year, sep = "_"))),
    af_n_com_to_row = as.integer(factor(paste(af_n_com_to_row, year, sep = "_"))),
    rr_n_com = as.integer(factor(paste(rr_n_com, year, sep = "_")))
  )

# Estimation for COMESA
# Agriculture trade
agri_com = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_com
                  + com + com_to_af + af_to_com 
                  + af_n_com + af_n_com_to_row + com_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_com)

# Delta method for Agriculture
expr_com <- paste0("(exp(com) - 1) * 100")
expr_com_to_af <- paste0("(exp(com_to_af) - 1) * 100")
expr_af_to_com <- paste0("(exp(af_to_com) - 1) * 100")
delta_com_agri <- deltaMethod(object = coef(agri_com), vcov. = vcov(agri_com), g = expr_com, parameterNames = names(coef(agri_com)))
delta_com_to_af_agri <- deltaMethod(object = coef(agri_com), vcov. = vcov(agri_com), g = expr_com_to_af, parameterNames = names(coef(agri_com)))
delta_af_to_com_agri <- deltaMethod(object = coef(agri_com), vcov. = vcov(agri_com), g = expr_af_to_com, parameterNames = names(coef(agri_com)))
print(delta_com_agri)
print(delta_com_to_af_agri)
print(delta_af_to_com_agri)

# Manufacturing trade
manu_com = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_com
                  + com + com_to_af + af_to_com 
                  + af_n_com + af_n_com_to_row + com_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_com)

# Delta method for Manufacturing
delta_com_manu <- deltaMethod(object = coef(manu_com), vcov. = vcov(manu_com), g = expr_com, parameterNames = names(coef(manu_com)))
delta_com_to_af_manu <- deltaMethod(object = coef(manu_com), vcov. = vcov(manu_com), g = expr_com_to_af, parameterNames = names(coef(manu_com)))
delta_af_to_com_manu <- deltaMethod(object = coef(manu_com), vcov. = vcov(manu_com), g = expr_af_to_com, parameterNames = names(coef(manu_com)))
print(delta_com_manu)
print(delta_com_to_af_manu)
print(delta_af_to_com_manu)

# Mining and Energy Trade
mine_com = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_com
                  + com + com_to_af + af_to_com 
                  + af_n_com + af_n_com_to_row + com_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_com)

# Delta method for Mining and Energy
delta_com_mine <- deltaMethod(object = coef(mine_com), vcov. = vcov(mine_com), g = expr_com, parameterNames = names(coef(mine_com)))
delta_com_to_af_mine <- deltaMethod(object = coef(mine_com), vcov. = vcov(mine_com), g = expr_com_to_af, parameterNames = names(coef(mine_com)))
delta_af_to_com_mine <- deltaMethod(object = coef(mine_com), vcov. = vcov(mine_com), g = expr_af_to_com, parameterNames = names(coef(mine_com)))
print(delta_com_mine)
print(delta_com_to_af_mine)
print(delta_af_to_com_mine)

# Service trade
serv_com = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_com
                  + com + com_to_af + af_to_com 
                  + af_n_com + af_n_com_to_row + com_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_com)

# Delta method for Services
#delta_com_serv <- deltaMethod(object = coef(serv_com), vcov. = vcov(serv_com), g = expr_com, parameterNames = names(coef(serv_com)))
delta_com_to_af_serv <- deltaMethod(object = coef(serv_com), vcov. = vcov(serv_com), g = expr_com_to_af, parameterNames = names(coef(serv_com)))
delta_af_to_com_serv <- deltaMethod(object = coef(serv_com), vcov. = vcov(serv_com), g = expr_af_to_com, parameterNames = names(coef(serv_com)))
#print(delta_com_serv)
print(delta_com_to_af_serv)
print(delta_af_to_com_serv)

# COMESA Latex results
tab_com <- huxreg("Agri" = agri_com,
                  "Manu" = manu_com,
                  "Mine" = mine_com,
                  "Serv" = serv_com,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_com",
                            "BRDR_COM" = "com",
                            "BRDR_COM_AFC" = "com_to_af",
                            "BRDR_AFC_COM" = "af_to_com",
                            "BRDR_AFC" = "af_n_com",
                            "BRDR_AFC_ROW" = "af_n_com_to_row",
                            "BRDR_COM_ROW" = "com_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Common Market for Eastern and Southern Africa (COMESA) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on COMESA Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_com")
width(tab_com) <- 1.1
cat(to_latex(tab_com), file = "output/tables/tab_com.tex")



# East African Community (EAC) -------------------------------------------------
eac_countries <- c("BDI", "COD", "KEN", "RWA", "SSD", "TZA", "UGA", "SOM")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    eac = ifelse(exporter %in% eac_countries &
                   importer %in% eac_countries &
                   exporter != importer, 1, 0),
    
    eac_to_af = ifelse(exporter %in% eac_countries &
                         importer %in% afc_countries &
                         !(importer %in% eac_countries) &
                         exporter != importer, 1, 0),
    
    af_to_eac = ifelse(importer %in% eac_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% eac_countries) &
                         exporter != importer, 1, 0),
    
    af_n_eac = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% eac_countries) &
                        !(importer %in% eac_countries) &
                        exporter != importer, 1, 0),
    
    eac_to_row = ifelse(exporter %in% eac_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_eac_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% eac_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_eac = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    eac = as.integer(factor(paste(eac, year, sep = "_"))),
    eac_to_af = as.integer(factor(paste(eac_to_af, year, sep = "_"))),
    af_to_eac = as.integer(factor(paste(af_to_eac, year, sep = "_"))),
    af_n_eac = as.integer(factor(paste(af_n_eac, year, sep = "_"))),
    eac_to_row = as.integer(factor(paste(eac_to_row, year, sep = "_"))),
    af_n_eac_to_row = as.integer(factor(paste(af_n_eac_to_row, year, sep = "_"))),
    rr_n_eac = as.integer(factor(paste(rr_n_eac, year, sep = "_")))
  )


# Estimation for EAC
# Agriculture trade
agri_eac = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eac
                  + eac + eac_to_af + af_to_eac 
                  + af_n_eac + af_n_eac_to_row + eac_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_eac)

# Delta method for Agriculture
expr_eac <- paste0("(exp(eac) - 1) * 100")
expr_eac_to_af <- paste0("(exp(eac_to_af) - 1) * 100")
expr_af_to_eac <- paste0("(exp(af_to_eac) - 1) * 100")
delta_eac_agri <- deltaMethod(object = coef(agri_eac), vcov. = vcov(agri_eac), g = expr_eac, parameterNames = names(coef(agri_eac)))
delta_eac_to_af_agri <- deltaMethod(object = coef(agri_eac), vcov. = vcov(agri_eac), g = expr_eac_to_af, parameterNames = names(coef(agri_eac)))
delta_af_to_eac_agri <- deltaMethod(object = coef(agri_eac), vcov. = vcov(agri_eac), g = expr_af_to_eac, parameterNames = names(coef(agri_eac)))
print(delta_eac_agri)
print(delta_eac_to_af_agri)
print(delta_af_to_eac_agri)

# Manufacturing trade
manu_eac = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eac
                  + eac + eac_to_af + af_to_eac 
                  + af_n_eac + af_n_eac_to_row + eac_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_eac)

# Delta method for Manufacturing
delta_eac_manu <- deltaMethod(object = coef(manu_eac), vcov. = vcov(manu_eac), g = expr_eac, parameterNames = names(coef(manu_eac)))
delta_eac_to_af_manu <- deltaMethod(object = coef(manu_eac), vcov. = vcov(manu_eac), g = expr_eac_to_af, parameterNames = names(coef(manu_eac)))
delta_af_to_eac_manu <- deltaMethod(object = coef(manu_eac), vcov. = vcov(manu_eac), g = expr_af_to_eac, parameterNames = names(coef(manu_eac)))
print(delta_eac_manu)
print(delta_eac_to_af_manu)
print(delta_af_to_eac_manu)

# Mining and Energy Trade
mine_eac = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eac
                  + eac + eac_to_af + af_to_eac 
                  + af_n_eac + af_n_eac_to_row + eac_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_eac)

# Delta method for Mining and Energy
delta_eac_mine <- deltaMethod(object = coef(mine_eac), vcov. = vcov(mine_eac), g = expr_eac, parameterNames = names(coef(mine_eac)))
delta_eac_to_af_mine <- deltaMethod(object = coef(mine_eac), vcov. = vcov(mine_eac), g = expr_eac_to_af, parameterNames = names(coef(mine_eac)))
delta_af_to_eac_mine <- deltaMethod(object = coef(mine_eac), vcov. = vcov(mine_eac), g = expr_af_to_eac, parameterNames = names(coef(mine_eac)))
print(delta_eac_mine)
print(delta_eac_to_af_mine)
print(delta_af_to_eac_mine)

# Service trade
serv_eac = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eac
                  + eac + eac_to_af + af_to_eac 
                  + af_n_eac + af_n_eac_to_row + eac_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_eac)

# Delta method for Services
#delta_eac_serv <- deltaMethod(object = coef(serv_eac), vcov. = vcov(serv_eac), g = expr_eac, parameterNames = names(coef(serv_eac)))
#delta_eac_to_af_serv <- deltaMethod(object = coef(serv_eac), vcov. = vcov(serv_eac), g = expr_eac_to_af, parameterNames = names(coef(serv_eac)))
#delta_af_to_eac_serv <- deltaMethod(object = coef(serv_eac), vcov. = vcov(serv_eac), g = expr_af_to_eac, parameterNames = names(coef(serv_eac)))
#print(delta_eac_serv)
#print(delta_eac_to_af_serv)
#print(delta_af_to_eac_serv)

# EAC Latex results
tab_eac <- huxreg("Agri" = agri_eac,
                  "Manu" = manu_eac,
                  "Mine" = mine_eac,
                  "Serv" = serv_eac,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_eac",
                            "BRDR_EAC" = "eac",
                            "BRDR_EAC_AFC" = "eac_to_af",
                            "BRDR_AFC_EAC" = "af_to_eac",
                            "BRDR_AFC" = "af_n_eac",
                            "BRDR_AFC_ROW" = "af_n_eac_to_row",
                            "BRDR_EAC_ROW" = "eac_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the East African Community (EAC) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on EAC Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_eac")
width(tab_eac) <- 1.1
cat(to_latex(tab_eac), file = "output/tables/tab_eac.tex")



# Economic Community of Central African States (ECCAS) -------------------------
ecc_countries <- c("AGO", "BDI", "CMR", "CAF", "TCD", "COG", "COD", "GNQ", "GAB", "RWA", "STP")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    ecc = ifelse(exporter %in% ecc_countries &
                   importer %in% ecc_countries &
                   exporter != importer, 1, 0),
    
    ecc_to_af = ifelse(exporter %in% ecc_countries &
                         importer %in% afc_countries &
                         !(importer %in% ecc_countries) &
                         exporter != importer, 1, 0),
    
    af_to_ecc = ifelse(importer %in% ecc_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% ecc_countries) &
                         exporter != importer, 1, 0),
    
    af_n_ecc = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% ecc_countries) &
                        !(importer %in% ecc_countries) &
                        exporter != importer, 1, 0),
    
    ecc_to_row = ifelse(exporter %in% ecc_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_ecc_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% ecc_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_ecc = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    ecc = as.integer(factor(paste(ecc, year, sep = "_"))),
    ecc_to_af = as.integer(factor(paste(ecc_to_af, year, sep = "_"))),
    af_to_ecc = as.integer(factor(paste(af_to_ecc, year, sep = "_"))),
    af_n_ecc = as.integer(factor(paste(af_n_ecc, year, sep = "_"))),
    ecc_to_row = as.integer(factor(paste(ecc_to_row, year, sep = "_"))),
    af_n_ecc_to_row = as.integer(factor(paste(af_n_ecc_to_row, year, sep = "_"))),
    rr_n_ecc = as.integer(factor(paste(rr_n_ecc, year, sep = "_")))
  )


# Estimation for ECCAS
# Agriculture trade
agri_ecc = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_ecc
                  + ecc + ecc_to_af + af_to_ecc 
                  + af_n_ecc + af_n_ecc_to_row + ecc_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_ecc)

# Delta method for Agriculture
expr_ecc <- paste0("(exp(ecc) - 1) * 100")
expr_ecc_to_af <- paste0("(exp(ecc_to_af) - 1) * 100")
expr_af_to_ecc <- paste0("(exp(af_to_ecc) - 1) * 100")
delta_ecc_agri <- deltaMethod(object = coef(agri_ecc), vcov. = vcov(agri_ecc), g = expr_ecc, parameterNames = names(coef(agri_ecc)))
delta_ecc_to_af_agri <- deltaMethod(object = coef(agri_ecc), vcov. = vcov(agri_ecc), g = expr_ecc_to_af, parameterNames = names(coef(agri_ecc)))
delta_af_to_ecc_agri <- deltaMethod(object = coef(agri_ecc), vcov. = vcov(agri_ecc), g = expr_af_to_ecc, parameterNames = names(coef(agri_ecc)))
print(delta_ecc_agri)
print(delta_ecc_to_af_agri)
print(delta_af_to_ecc_agri)

# Manufacturing trade
manu_ecc = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_ecc
                  + ecc + ecc_to_af + af_to_ecc 
                  + af_n_ecc + af_n_ecc_to_row + ecc_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_ecc)

# Delta method for Manufacturing
delta_ecc_manu <- deltaMethod(object = coef(manu_ecc), vcov. = vcov(manu_ecc), g = expr_ecc, parameterNames = names(coef(manu_ecc)))
delta_ecc_to_af_manu <- deltaMethod(object = coef(manu_ecc), vcov. = vcov(manu_ecc), g = expr_ecc_to_af, parameterNames = names(coef(manu_ecc)))
delta_af_to_ecc_manu <- deltaMethod(object = coef(manu_ecc), vcov. = vcov(manu_ecc), g = expr_af_to_ecc, parameterNames = names(coef(manu_ecc)))
print(delta_ecc_manu)
print(delta_ecc_to_af_manu)
print(delta_af_to_ecc_manu)

# Mining and Energy Trade
mine_ecc = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_ecc
                  + ecc + ecc_to_af + af_to_ecc 
                  + af_n_ecc + af_n_ecc_to_row + ecc_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_ecc)

# Delta method for Mining and Energy
delta_ecc_mine <- deltaMethod(object = coef(mine_ecc), vcov. = vcov(mine_ecc), g = expr_ecc, parameterNames = names(coef(mine_ecc)))
delta_ecc_to_af_mine <- deltaMethod(object = coef(mine_ecc), vcov. = vcov(mine_ecc), g = expr_ecc_to_af, parameterNames = names(coef(mine_ecc)))
delta_af_to_ecc_mine <- deltaMethod(object = coef(mine_ecc), vcov. = vcov(mine_ecc), g = expr_af_to_ecc, parameterNames = names(coef(mine_ecc)))
print(delta_ecc_mine)
print(delta_ecc_to_af_mine)
print(delta_af_to_ecc_mine)

# Service trade
serv_ecc = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_ecc
                  + ecc + ecc_to_af + af_to_ecc 
                  + af_n_ecc + af_n_ecc_to_row + ecc_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_ecc)

# Delta method for Services
#delta_ecc_serv <- deltaMethod(object = coef(serv_ecc), vcov. = vcov(serv_ecc), g = expr_ecc, parameterNames = names(coef(serv_ecc)))
#delta_ecc_to_af_serv <- deltaMethod(object = coef(serv_ecc), vcov. = vcov(serv_ecc), g = expr_ecc_to_af, parameterNames = names(coef(serv_ecc)))
#delta_af_to_ecc_serv <- deltaMethod(object = coef(serv_ecc), vcov. = vcov(serv_ecc), g = expr_af_to_ecc, parameterNames = names(coef(serv_ecc)))
#print(delta_ecc_serv)
#print(delta_ecc_to_af_serv)
#print(delta_af_to_ecc_serv)

# ECC Latex results
tab_ecc <- huxreg("Agri" = agri_ecc,
                  "Manu" = manu_ecc,
                  "Mine" = mine_ecc,
                  "Serv" = serv_ecc,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_ecc",
                            "BRDR_ECC" = "ecc",
                            "BRDR_ECC_AFC" = "ecc_to_af",
                            "BRDR_AFC_ecc" = "af_to_ecc",
                            "BRDR_AFC" = "af_n_ecc",
                            "BRDR_AFC_ROW" = "af_n_ecc_to_row",
                            "BRDR_ECC_ROW" = "ecc_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Economic Community of Central African States (ECCAS) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on ECCAS Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_ecc")
width(tab_ecc) <- 1.1
cat(to_latex(tab_ecc), file = "output/tables/tab_ecc.tex")




# Economic Community of West African States (ECOWAS) ---------------------------
eco_countries <- c("BEN", "BFA", "CPV", "CIV", "GMB", "GHA", "GIN", "GNB", "LBR", "MLI", 
                   "NER", "NGA", "SEN", "SLE", "TGO")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    eco = ifelse(exporter %in% eco_countries &
                   importer %in% eco_countries &
                   exporter != importer, 1, 0),
    
    eco_to_af = ifelse(exporter %in% eco_countries &
                         importer %in% afc_countries &
                         !(importer %in% eco_countries) &
                         exporter != importer, 1, 0),
    
    af_to_eco = ifelse(importer %in% eco_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% eco_countries) &
                         exporter != importer, 1, 0),
    
    af_n_eco = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% eco_countries) &
                        !(importer %in% eco_countries) &
                        exporter != importer, 1, 0),
    
    eco_to_row = ifelse(exporter %in% eco_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_eco_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% eco_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_eco = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    eco = as.integer(factor(paste(eco, year, sep = "_"))),
    eco_to_af = as.integer(factor(paste(eco_to_af, year, sep = "_"))),
    af_to_eco = as.integer(factor(paste(af_to_eco, year, sep = "_"))),
    af_n_eco = as.integer(factor(paste(af_n_eco, year, sep = "_"))),
    eco_to_row = as.integer(factor(paste(eco_to_row, year, sep = "_"))),
    af_n_eco_to_row = as.integer(factor(paste(af_n_eco_to_row, year, sep = "_"))),
    rr_n_eco = as.integer(factor(paste(rr_n_eco, year, sep = "_")))
  )


# Estimation for ECOWAS
# Agriculture trade
agri_eco = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eco
                  + eco + eco_to_af + af_to_eco 
                  + af_n_eco + af_n_eco_to_row + eco_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_eco)

# Delta method for Agriculture
expr_eco <- paste0("(exp(eco) - 1) * 100")
expr_eco_to_af <- paste0("(exp(eco_to_af) - 1) * 100")
expr_af_to_eco <- paste0("(exp(af_to_eco) - 1) * 100")
delta_eco_agri <- deltaMethod(object = coef(agri_eco), vcov. = vcov(agri_eco), g = expr_eco, parameterNames = names(coef(agri_eco)))
delta_eco_to_af_agri <- deltaMethod(object = coef(agri_eco), vcov. = vcov(agri_eco), g = expr_eco_to_af, parameterNames = names(coef(agri_eco)))
delta_af_to_eco_agri <- deltaMethod(object = coef(agri_eco), vcov. = vcov(agri_eco), g = expr_af_to_eco, parameterNames = names(coef(agri_eco)))
print(delta_eco_agri)
print(delta_eco_to_af_agri)
print(delta_af_to_eco_agri)

# Manufacturing trade
manu_eco = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eco
                  + eco + eco_to_af + af_to_eco 
                  + af_n_eco + af_n_eco_to_row + eco_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_eco)

# Delta method for Manufacturing
delta_eco_manu <- deltaMethod(object = coef(manu_eco), vcov. = vcov(manu_eco), g = expr_eco, parameterNames = names(coef(manu_eco)))
delta_eco_to_af_manu <- deltaMethod(object = coef(manu_eco), vcov. = vcov(manu_eco), g = expr_eco_to_af, parameterNames = names(coef(manu_eco)))
delta_af_to_eco_manu <- deltaMethod(object = coef(manu_eco), vcov. = vcov(manu_eco), g = expr_af_to_eco, parameterNames = names(coef(manu_eco)))
print(delta_eco_manu)
print(delta_eco_to_af_manu)
print(delta_af_to_eco_manu)

# Mining and Energy Trade
mine_eco = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eco
                  + eco + eco_to_af + af_to_eco 
                  + af_n_eco + af_n_eco_to_row + eco_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_eco)

# Delta method for Mining and Energy
delta_eco_mine <- deltaMethod(object = coef(mine_eco), vcov. = vcov(mine_eco), g = expr_eco, parameterNames = names(coef(mine_eco)))
delta_eco_to_af_mine <- deltaMethod(object = coef(mine_eco), vcov. = vcov(mine_eco), g = expr_eco_to_af, parameterNames = names(coef(mine_eco)))
delta_af_to_eco_mine <- deltaMethod(object = coef(mine_eco), vcov. = vcov(mine_eco), g = expr_af_to_eco, parameterNames = names(coef(mine_eco)))
print(delta_eco_mine)
print(delta_eco_to_af_mine)
print(delta_af_to_eco_mine)

# Service trade
serv_eco = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_eco
                  + eco + eco_to_af + af_to_eco 
                  + af_n_eco + af_n_eco_to_row + eco_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_eco)

# Delta method for Services
#delta_eco_serv <- deltaMethod(object = coef(serv_eco), vcov. = vcov(serv_eco), g = expr_eco, parameterNames = names(coef(serv_eco)))
#delta_eco_to_af_serv <- deltaMethod(object = coef(serv_eco), vcov. = vcov(serv_eco), g = expr_eco_to_af, parameterNames = names(coef(serv_eco)))
#delta_af_to_eco_serv <- deltaMethod(object = coef(serv_eco), vcov. = vcov(serv_eco), g = expr_af_to_eco, parameterNames = names(coef(serv_eco)))
#print(delta_eco_serv)
#print(delta_eco_to_af_serv)
#print(delta_af_to_eco_serv)

# ECO Latex results
tab_eco <- huxreg("Agri" = agri_eco,
                  "Manu" = manu_eco,
                  "Mine" = mine_eco,
                  "Serv" = serv_eco,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_eco",
                            "BRDR_ECO" = "eco",
                            "BRDR_ECO_AFC" = "eco_to_af",
                            "BRDR_AFC_ECO" = "af_to_eco",
                            "BRDR_AFC" = "af_n_eco",
                            "BRDR_AFC_ROW" = "af_n_eco_to_row",
                            "BRDR_ECO_ROW" = "eco_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Economic Community of West African States (ECOWAS) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on ECOWAS Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_eco")
width(tab_eco) <- 1.1
cat(to_latex(tab_eco), file = "output/tables/tab_eco.tex")



# Intergovernmental Authority on Development (IGAD) — excluding Eritrea
iga_countries <- c("DJI", "ETH", "KEN", "SOM", "SSD", "SDN", "UGA")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    iga = ifelse(exporter %in% iga_countries &
                   importer %in% iga_countries &
                   exporter != importer, 1, 0),
    
    iga_to_af = ifelse(exporter %in% iga_countries &
                         importer %in% afc_countries &
                         !(importer %in% iga_countries) &
                         exporter != importer, 1, 0),
    
    af_to_iga = ifelse(importer %in% iga_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% iga_countries) &
                         exporter != importer, 1, 0),
    
    af_n_iga = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% iga_countries) &
                        !(importer %in% iga_countries) &
                        exporter != importer, 1, 0),
    
    iga_to_row = ifelse(exporter %in% iga_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_iga_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% iga_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_iga = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    iga = as.integer(factor(paste(iga, year, sep = "_"))),
    iga_to_af = as.integer(factor(paste(iga_to_af, year, sep = "_"))),
    af_to_iga = as.integer(factor(paste(af_to_iga, year, sep = "_"))),
    af_n_iga = as.integer(factor(paste(af_n_iga, year, sep = "_"))),
    iga_to_row = as.integer(factor(paste(iga_to_row, year, sep = "_"))),
    af_n_iga_to_row = as.integer(factor(paste(af_n_iga_to_row, year, sep = "_"))),
    rr_n_iga = as.integer(factor(paste(rr_n_iga, year, sep = "_")))
  )


# Estimation for IGAD
# Agriculture trade
agri_iga = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_iga
                  + iga + iga_to_af + af_to_iga 
                  + af_n_iga + af_n_iga_to_row + iga_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_iga)

# Delta method for Agriculture
expr_iga <- paste0("(exp(iga) - 1) * 100")
expr_iga_to_af <- paste0("(exp(iga_to_af) - 1) * 100")
expr_af_to_iga <- paste0("(exp(af_to_iga) - 1) * 100")
delta_iga_agri <- deltaMethod(object = coef(agri_iga), vcov. = vcov(agri_iga), g = expr_iga, parameterNames = names(coef(agri_iga)))
delta_iga_to_af_agri <- deltaMethod(object = coef(agri_iga), vcov. = vcov(agri_iga), g = expr_iga_to_af, parameterNames = names(coef(agri_iga)))
delta_af_to_iga_agri <- deltaMethod(object = coef(agri_iga), vcov. = vcov(agri_iga), g = expr_af_to_iga, parameterNames = names(coef(agri_iga)))
print(delta_iga_agri)
print(delta_iga_to_af_agri)
print(delta_af_to_iga_agri)

# Manufacturing trade
manu_iga = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_iga
                  + iga + iga_to_af + af_to_iga 
                  + af_n_iga + af_n_iga_to_row + iga_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_iga)

# Delta method for Manufacturing
delta_iga_manu <- deltaMethod(object = coef(manu_iga), vcov. = vcov(manu_iga), g = expr_iga, parameterNames = names(coef(manu_iga)))
delta_iga_to_af_manu <- deltaMethod(object = coef(manu_iga), vcov. = vcov(manu_iga), g = expr_iga_to_af, parameterNames = names(coef(manu_iga)))
delta_af_to_iga_manu <- deltaMethod(object = coef(manu_iga), vcov. = vcov(manu_iga), g = expr_af_to_iga, parameterNames = names(coef(manu_iga)))
print(delta_iga_manu)
print(delta_iga_to_af_manu)
print(delta_af_to_iga_manu)

# Mining and Energy Trade
mine_iga = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_iga
                  + iga + iga_to_af + af_to_iga 
                  + af_n_iga + af_n_iga_to_row + iga_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_iga)

# Delta method for Mining and Energy
delta_iga_mine <- deltaMethod(object = coef(mine_iga), vcov. = vcov(mine_iga), g = expr_iga, parameterNames = names(coef(mine_iga)))
delta_iga_to_af_mine <- deltaMethod(object = coef(mine_iga), vcov. = vcov(mine_iga), g = expr_iga_to_af, parameterNames = names(coef(mine_iga)))
delta_af_to_iga_mine <- deltaMethod(object = coef(mine_iga), vcov. = vcov(mine_iga), g = expr_af_to_iga, parameterNames = names(coef(mine_iga)))
print(delta_iga_mine)
print(delta_iga_to_af_mine)
print(delta_af_to_iga_mine)

# Service trade
serv_iga = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_iga
                  + iga + iga_to_af + af_to_iga 
                  + af_n_iga + af_n_iga_to_row + iga_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_iga)

# Delta method for Services
#delta_iga_serv <- deltaMethod(object = coef(serv_iga), vcov. = vcov(serv_iga), g = expr_iga, parameterNames = names(coef(serv_iga)))
#delta_iga_to_af_serv <- deltaMethod(object = coef(serv_iga), vcov. = vcov(serv_iga), g = expr_iga_to_af, parameterNames = names(coef(serv_iga)))
#delta_af_to_iga_serv <- deltaMethod(object = coef(serv_iga), vcov. = vcov(serv_iga), g = expr_af_to_iga, parameterNames = names(coef(serv_iga)))
#print(delta_iga_serv)
#print(delta_iga_to_af_serv)
#print(delta_af_to_iga_serv)

# IGAD Latex results
tab_iga <- huxreg("Agri" = agri_iga,
                  "Manu" = manu_iga,
                  "Mine" = mine_iga,
                  "Serv" = serv_iga,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_iga",
                            "BRDR_IGA" = "iga",
                            "BRDR_IGA_AFC" = "iga_to_af",
                            "BRDR_AFC_IGA" = "af_to_iga",
                            "BRDR_AFC" = "af_n_iga",
                            "BRDR_AFC_ROW" = "af_n_iga_to_row",
                            "BRDR_IGA_ROW" = "iga_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the # Intergovernmental Authority on Development (IGAD) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on IGAD Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_iga")
width(tab_iga) <- 1.1
cat(to_latex(tab_iga), file = "output/tables/tab_iga.tex")



# Southern African Development Community (SADC) --------------------------------
sad_countries <- c("AGO", "BWA", "COM", "COD", "SWZ", "LSO", "MDG", "MWI", "MUS", "MOZ", "NAM", 
                   "ZAF", "SYC", "TZA", "ZMB", "ZWE")

# construct time-varying border variables
afcfta_data <- afcfta_data %>%
  mutate(
    sad = ifelse(exporter %in% sad_countries &
                   importer %in% sad_countries &
                   exporter != importer, 1, 0),
    
    sad_to_af = ifelse(exporter %in% sad_countries &
                         importer %in% afc_countries &
                         !(importer %in% sad_countries) &
                         exporter != importer, 1, 0),
    
    af_to_sad = ifelse(importer %in% sad_countries &
                         exporter %in% afc_countries &
                         !(exporter %in% sad_countries) &
                         exporter != importer, 1, 0),
    
    af_n_sad = ifelse(exporter %in% afc_countries &
                        importer %in% afc_countries & 
                        !(exporter %in% sad_countries) &
                        !(importer %in% sad_countries) &
                        exporter != importer, 1, 0),
    
    sad_to_row = ifelse(exporter %in% sad_countries &
                          !(importer %in% afc_countries) & 
                          exporter != importer, 1, 0),
    
    af_n_sad_to_row = ifelse(exporter %in% afc_countries &
                               !(exporter %in% sad_countries) &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    rr_n_sad = ifelse(!(exporter %in% afc_countries) &
                        !(importer %in% afc_countries) &
                        exporter != importer, 1, 0)
  ) %>%
  mutate(
    sad = as.integer(factor(paste(sad, year, sep = "_"))),
    sad_to_af = as.integer(factor(paste(sad_to_af, year, sep = "_"))),
    af_to_sad = as.integer(factor(paste(af_to_sad, year, sep = "_"))),
    af_n_sad = as.integer(factor(paste(af_n_sad, year, sep = "_"))),
    sad_to_row = as.integer(factor(paste(sad_to_row, year, sep = "_"))),
    af_n_sad_to_row = as.integer(factor(paste(af_n_sad_to_row, year, sep = "_"))),
    rr_n_sad = as.integer(factor(paste(rr_n_sad, year, sep = "_")))
  )


# Estimation for SADC
# Agriculture trade
agri_sad = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_sad
                  + sad + sad_to_af + af_to_sad 
                  + af_n_sad + af_n_sad_to_row + sad_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_sad)

# Delta method for Agriculture
expr_sad <- paste0("(exp(sad) - 1) * 100")
expr_sad_to_af <- paste0("(exp(sad_to_af) - 1) * 100")
expr_af_to_sad <- paste0("(exp(af_to_sad) - 1) * 100")
delta_sad_agri <- deltaMethod(object = coef(agri_sad), vcov. = vcov(agri_sad), g = expr_sad, parameterNames = names(coef(agri_sad)))
delta_sad_to_af_agri <- deltaMethod(object = coef(agri_sad), vcov. = vcov(agri_sad), g = expr_sad_to_af, parameterNames = names(coef(agri_sad)))
delta_af_to_sad_agri <- deltaMethod(object = coef(agri_sad), vcov. = vcov(agri_sad), g = expr_af_to_sad, parameterNames = names(coef(agri_sad)))
print(delta_sad_agri)
print(delta_sad_to_af_agri)
print(delta_af_to_sad_agri)

# Manufacturing trade
manu_sad = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_sad
                  + sad + sad_to_af + af_to_sad 
                  + af_n_sad + af_n_sad_to_row + sad_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_sad)

# Delta method for Manufacturing
delta_sad_manu <- deltaMethod(object = coef(manu_sad), vcov. = vcov(manu_sad), g = expr_sad, parameterNames = names(coef(manu_sad)))
delta_sad_to_af_manu <- deltaMethod(object = coef(manu_sad), vcov. = vcov(manu_sad), g = expr_sad_to_af, parameterNames = names(coef(manu_sad)))
delta_af_to_sad_manu <- deltaMethod(object = coef(manu_sad), vcov. = vcov(manu_sad), g = expr_af_to_sad, parameterNames = names(coef(manu_sad)))
print(delta_sad_manu)
print(delta_sad_to_af_manu)
print(delta_af_to_sad_manu)

# Mining and Energy Trade
mine_sad = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_sad
                  + sad + sad_to_af + af_to_sad 
                  + af_n_sad + af_n_sad_to_row + sad_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_sad)

# Delta method for Mining and Energy
delta_sad_mine <- deltaMethod(object = coef(mine_sad), vcov. = vcov(mine_sad), g = expr_sad, parameterNames = names(coef(mine_sad)))
delta_sad_to_af_mine <- deltaMethod(object = coef(mine_sad), vcov. = vcov(mine_sad), g = expr_sad_to_af, parameterNames = names(coef(mine_sad)))
delta_af_to_sad_mine <- deltaMethod(object = coef(mine_sad), vcov. = vcov(mine_sad), g = expr_af_to_sad, parameterNames = names(coef(mine_sad)))
print(delta_sad_mine)
print(delta_sad_to_af_mine)
print(delta_af_to_sad_mine)

# Service trade
serv_sad = fepois(trade ~ cntg + col + dist + lang + rta + rr_n_sad
                  + sad + sad_to_af + af_to_sad 
                  + af_n_sad + af_n_sad_to_row + sad_to_row
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_sad)

# Delta method for Services
delta_sad_serv <- deltaMethod(object = coef(serv_sad), vcov. = vcov(serv_sad), g = expr_sad, parameterNames = names(coef(serv_sad)))
#delta_sad_to_af_serv <- deltaMethod(object = coef(serv_sad), vcov. = vcov(serv_sad), g = expr_sad_to_af, parameterNames = names(coef(serv_sad)))
#delta_af_to_sad_serv <- deltaMethod(object = coef(serv_sad), vcov. = vcov(serv_sad), g = expr_af_to_sad, parameterNames = names(coef(serv_sad)))
print(delta_sad_serv)
#print(delta_sad_to_af_serv)
#print(delta_af_to_sad_serv)

# SADC Latex results
tab_sad <- huxreg("Agri" = agri_sad,
                  "Manu" = manu_sad,
                  "Mine" = mine_sad,
                  "Serv" = serv_sad,
                  coefs = c("CNTG" = "cntg",
                            "COLN" = "col",
                            "DIST" = "dist",
                            "LANG" = "lang",
                            "RTA" = "rta",
                            "INTL_BRDR" = "rr_n_sad",
                            "BRDR_SAD" = "sad",
                            "BRDR_SAD_AFC" = "sad_to_af",
                            "BRDR_AFC_SAD" = "af_to_sad",
                            "BRDR_AFC" = "af_n_sad",
                            "BRDR_AFC_ROW" = "af_n_sad_to_row",
                            "BRDR_SAD_ROW" = "sad_to_row"),
                  stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                  note = "Notes: Statistics based on author's calculations. The table presents the gravity estimates of ex-ante impacts of the AfCFTA agreement on the Southern African Development Community (SADC) members' trade with the rest of the AfCFTA members for Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Ex-ante Impacts of the AfCFTA Agreement on SADC Sectoral Trade with the Rest of AfCFTA members") %>%
  set_label("tab_sad")
width(tab_sad) <- 1.1
cat(to_latex(tab_sad), file = "output/tables/tab_sad.tex")


################################################################################
### End of script --------------------------------------------------------------
################################################################################