################################################################################
### R Code for the descriptive analysis the ex-ante quantification of AfCFTA ###
###################### by Jamiu Olamilekan Badmus ##############################
################################################################################

# clear the workspace ----------------------------------------------------------
rm(list = ls())

# Load the necessary libraries -------------------------------------------------
library(tidyverse)
library(readr)
library(dplyr)
library(knitr)
library(kableExtra)
library(scales)

# Set the working directory ----------------------------------------------------
setwd("C:/Users/muham/EGEI Dissertation/afcfta")

# Load the cleaned ITPDER2 dataset ---------------------------------------------
itpde <- readRDS("input/itpder2.rds") # Trade Dataset


################################################################################
## Summary Statistics ----------------------------------------------------------
################################################################################

# List of African countries
afc_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Arab Maghreb Union (AMU)
amu_countries <- c("DZA", "LBY", "MAR", "MRT", "TUN")

# Community of Sahel–Saharan States (CEN-SAD) — excluding Eritrea
cen_countries <- c("BEN", "BFA", "CAF", "TCD", "COM", "CIV", "DJI", "EGY", "GMB", "GHA", "GIN",
                   "GNB", "KEN", "LBY", "MLI", "MRT", "MAR", "NER", "NGA", "SEN", "SLE", "SDN",
                   "SOM", "TGO", "TUN")

# Common Market for Eastern and Southern Africa (COMESA) — excluding Eritrea
com_countries <- c("BDI", "COM", "COD", "DJI", "EGY", "SWZ", "ETH", "KEN", "LYB", "MDG", "MWI",
                   "MUS", "RWA", "SYC", "SOM", "SDN", "TUN", "UGA", "ZMB", "ZWE")

# East African Community (EAC)
eac_countries <- c("BDI", "COD", "KEN", "RWA", "SSD", "TZA", "UGA", "SOM")

# Economic Community of Central African States (ECCAS)
ecc_countries <- c("AGO", "BDI", "CMR", "CAF", "TCD", "COG", "COD", "GNQ", "GAB", "RWA", "STP")

# Economic Community of West African States (ECOWAS)
eco_countries <- c("BEN", "BFA", "CPV", "CIV", "GMB", "GHA", "GIN", "GNB", "LBR", "MLI", 
                   "NER", "NGA", "SEN", "SLE", "TGO")

# Intergovernmental Authority on Development (IGAD) — excluding Eritrea
iga_countries <- c("DJI", "ETH", "KEN", "SOM", "SSD", "SDN", "UGA")

# Southern African Development Community (SADC)
sad_countries <- c("AGO", "BWA", "COM", "COD", "SWZ", "LSO", "MDG", "MWI", "MUS", "MOZ", "NAM", 
                   "ZAF", "SYC", "TZA", "ZMB", "ZWE")


# Create each subset for summary statistics (include both intra-and inter-national flows)
intra_africa <- itpde %>%
  filter(exporter %in% afc_countries, importer %in% afc_countries, exporter != importer)

amu_to_afc <- itpde %>%
  filter(exporter %in% amu_countries,
         importer %in% afc_countries,
         !(importer %in% amu_countries), exporter != importer)

afc_to_amu <- itpde %>%
  filter(importer %in% amu_countries,
         exporter %in% afc_countries,
         !(exporter %in% amu_countries), exporter != importer)

cen_to_afc <- itpde %>%
  filter(exporter %in% cen_countries,
         importer %in% afc_countries,
         !(importer %in% cen_countries), exporter != importer)

afc_to_cen <- itpde %>%
  filter(importer %in% cen_countries,
         exporter %in% afc_countries,
         !(exporter %in% cen_countries), exporter != importer)

com_to_afc <- itpde %>%
  filter(exporter %in% com_countries,
         importer %in% afc_countries,
         !(importer %in% com_countries), exporter != importer)

afc_to_com <- itpde %>%
  filter(importer %in% com_countries,
         exporter %in% afc_countries,
         !(exporter %in% com_countries), exporter != importer)

eac_to_afc <- itpde %>%
  filter(exporter %in% eac_countries,
         importer %in% afc_countries,
         !(importer %in% eac_countries), exporter != importer)

afc_to_eac <- itpde %>%
  filter(importer %in% eac_countries,
         exporter %in% afc_countries,
         !(exporter %in% eac_countries), exporter != importer)

ecc_to_afc <- itpde %>%
  filter(exporter %in% ecc_countries,
         importer %in% afc_countries,
         !(importer %in% ecc_countries), exporter != importer)

afc_to_ecc <- itpde %>%
  filter(importer %in% ecc_countries,
         exporter %in% afc_countries,
         !(exporter %in% ecc_countries), exporter != importer)

eco_to_afc <- itpde %>%
  filter(exporter %in% eco_countries,
         importer %in% afc_countries,
         !(importer %in% eco_countries), exporter != importer)

afc_to_eco <- itpde %>%
  filter(importer %in% eco_countries,
         exporter %in% afc_countries,
         !(exporter %in% eco_countries), exporter != importer)

iga_to_afc <- itpde %>%
  filter(exporter %in% iga_countries,
         importer %in% afc_countries,
         !(importer %in% iga_countries), exporter != importer)

afc_to_iga <- itpde %>%
  filter(importer %in% iga_countries,
         exporter %in% afc_countries,
         !(exporter %in% iga_countries), exporter != importer)

sad_to_afc <- itpde %>%
  filter(exporter %in% sad_countries,
         importer %in% afc_countries,
         !(importer %in% sad_countries), exporter != importer)

afc_to_sad <- itpde %>%
  filter(importer %in% sad_countries,
         exporter %in% afc_countries,
         !(exporter %in% sad_countries), exporter != importer)


# Function to generate summary
sector_summary <- function(data, label) {
  data %>%
    group_by(broad_sector) %>%
    summarise(
      count = n(),
      pairs = n_distinct(paste(exporter, importer, sep = "-")),
      mean = mean(trade, na.rm = TRUE),
      median = median(trade, na.rm = TRUE),
      sd = sd(trade, na.rm = TRUE),
      min = min(trade, na.rm = TRUE),
      max = max(trade, na.rm = TRUE)
    ) %>%
    mutate(dataset = label)
}

# Apply function to generate all summaries
intra_africa_sum <- sector_summary(intra_africa, "Intra_Africa")
amu_to_afc_sum <- sector_summary(amu_to_afc, "AMU_to_AFC")
afc_to_amu_sum <- sector_summary(afc_to_amu, "AFC_to_AMU")
cen_to_afc_sum <- sector_summary(cen_to_afc, "CEN_to_AFC")
afc_to_cen_sum <- sector_summary(afc_to_cen, "AFC_to_CEN")
com_to_afc_sum <- sector_summary(com_to_afc, "COM_to_AFC")
afc_to_com_sum <- sector_summary(afc_to_com, "AFC_to_COM")
eac_to_afc_sum <- sector_summary(eac_to_afc, "EAC_to_AFC")
afc_to_eac_sum <- sector_summary(afc_to_eac, "AFC_to_EAC")
ecc_to_afc_sum <- sector_summary(ecc_to_afc, "ECC_to_AFC")
afc_to_ecc_sum <- sector_summary(afc_to_ecc, "AFC_to_ECC")
eco_to_afc_sum <- sector_summary(eco_to_afc, "ECO_to_AFC")
afc_to_eco_sum <- sector_summary(afc_to_eco, "AFC_to_ECO")
iga_to_afc_sum <- sector_summary(iga_to_afc, "IGA_to_AFC")
afc_to_iga_sum <- sector_summary(afc_to_iga, "AFC_to_IGA")
sad_to_afc_sum <- sector_summary(sad_to_afc, "SAD_to_AFC")
afc_to_sad_sum <- sector_summary(afc_to_sad, "AFC_to_SAD")

# Combine and format
combined_summary <- bind_rows(
  intra_africa_sum,
  amu_to_afc_sum, afc_to_amu_sum,
  cen_to_afc_sum, afc_to_cen_sum,
  com_to_afc_sum, afc_to_com_sum,
  eac_to_afc_sum, afc_to_eac_sum,
  ecc_to_afc_sum, afc_to_ecc_sum,
  eco_to_afc_sum, afc_to_eco_sum,
  iga_to_afc_sum, afc_to_iga_sum,
  sad_to_afc_sum, afc_to_sad_sum
) %>%
  arrange(dataset, broad_sector) %>%
  select(dataset, broad_sector, count, pairs, mean, median, sd, min, max) %>%
  mutate(across(where(is.numeric), round, 2))

# LaTeX table
# First check the structure of combined_summary to get correct row counts
dataset_counts <- table(combined_summary$dataset)
print(dataset_counts)

# Generate the kable with shorter lines
summary_stat <- kable(combined_summary[, -1], 
                      format = "latex", 
                      booktabs = TRUE,
                      col.names = c("Broad Sector", "Observations", "Pairs", "Mean", 
                                    "Median", "Std. Dev.", "Min", "Max"),
                      caption = paste("Summary Statistics of Trade by Broad Sector",
                                      "and Country Pairs")) %>%
  kable_styling(latex_options = c("hold_position", "striped", "scale_down"))

# Apply grouping dynamically based on the actual data
start_row <- 1
for (ds in unique(combined_summary$dataset)) {
  # Get count of rows for this dataset
  ds_count <- sum(combined_summary$dataset == ds)
  if (ds_count > 0) {
    end_row <- start_row + ds_count - 1
    summary_stat <- summary_stat %>%
      group_rows(ds, start_row, end_row)
    start_row <- end_row + 1
  }
}

# Save to file
save_kable(summary_stat, file = "output/tables/summary_statistics.tex")



################################################################################
## Trade flows and tariff graphs: Intra-African Trade --------------------------
################################################################################

# Plot intra-African trade flows and tariff over time by broad sector ----------
# Summarize trade by year and broad sector
intra_africa_year <- intra_africa %>%
  group_by(year, broad_sector) %>%
  summarise(total_trade = sum(trade, na.rm = TRUE), .groups = "drop")

# Load intra-Africa_tariff data
intra_africa_tariff <- read_csv("input/intra_africa_tariff.csv")

# Filter AHS tariff data
intra_africa_tariff <- intra_africa_tariff %>%
  filter(DutyType == "AHS") %>%
  select('Reporter Name', Product, 'Product Name', 'Partner Name', 'Tariff Year', DutyType, 'Weighted Average')

# Subset the year to 2000 to 2019
intra_africa_tariff <- intra_africa_tariff %>%
  filter(`Tariff Year` >= 2000 & `Tariff Year` <= 2019)

# Prepare the ahs_data for sectoral mapping
intra_africa_tariff <- intra_africa_tariff %>%
  mutate(hs6 = substr(Product, 1, 6),
         hs2 = substr(hs6, 1, 2),
         sector = case_when(
           hs2 %in% sprintf("%02d", 1:24) ~ "Agriculture",
           hs2 %in% c("25", "26", "27") ~ "Mining and Energy",
           hs2 %in% sprintf("%02d", 28:99) ~ "Manufacturing",
           TRUE ~ NA_character_
         ))

# Compute weighted sectoral averages
intra_africa_tariff <- intra_africa_tariff %>%
  filter(sector %in% c("Agriculture", "Manufacturing", "Mining and Energy")) %>%
  group_by(`Tariff Year`, sector) %>%
  summarise(
    avg_ahs_tariff = mean(`Weighted Average`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(year = `Tariff Year`, broad_sector = sector)

# Merge trade and tariff data for sectors with tariff data
intra_africa_trade_tariff <- intra_africa_year %>%
  filter(broad_sector %in% c("Agriculture", "Manufacturing", "Mining and Energy")) %>%
  left_join(intra_africa_tariff, by = c("year", "broad_sector"))

# Function to create individual sector plots with dual y-axes
create_sector_plot <- function(data, sector_name, subplot_label) {
  # Filter data for specific sector
  sector_data <- data %>% filter(broad_sector == sector_name)
  
  # Calculate scaling factor for secondary axis
  trade_max <- max(sector_data$total_trade, na.rm = TRUE)
  tariff_max <- max(sector_data$avg_ahs_tariff, na.rm = TRUE)
  scale_factor <- trade_max / tariff_max
  
  # Create the plot
  ggplot(sector_data, aes(x = year)) +
    # Trade line - primary y-axis (left)
    geom_line(aes(y = total_trade, color = "Total Trade"), 
              size = 1.2, alpha = 0.8) +
    geom_point(aes(y = total_trade, color = "Total Trade"), 
               size = 2, alpha = 0.8) +
    
    # Tariff line - secondary y-axis (right)
    geom_line(aes(y = avg_ahs_tariff * scale_factor, color = "Average Tariff"), 
              size = 1.2, linetype = "dashed", alpha = 0.8) +
    geom_point(aes(y = avg_ahs_tariff * scale_factor, color = "Average Tariff"), 
               size = 2, alpha = 0.8) +
    
    # Primary y-axis for trade
    scale_y_continuous(
      name = "Total Trade (Million USD)",
      labels = comma_format(),
      # Secondary y-axis for tariffs
      sec.axis = sec_axis(~ . / scale_factor, 
                          name = "AHS Weighted Average Tariff (%)",
                          labels = function(x) sprintf("%.1f", x))
    ) +
    
    # Customize colors
    scale_color_manual(values = c("Total Trade" = "steelblue", 
                                  "Average Tariff" = "darkred")) +
    
    # Labels
    labs(
      title = paste0("(", subplot_label, ") ", sector_name),
      x = "Year",
      color = NULL
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      axis.title.y.left = element_text(color = "steelblue", size = 10),
      axis.title.y.right = element_text(color = "darkred", size = 10),
      axis.text.y.left = element_text(color = "steelblue"),
      axis.text.y.right = element_text(color = "darkred")
    ) +
    
    scale_x_continuous(breaks = seq(2000, 2019, 4))
}

# Function to create Services plot (trade only)
create_services_plot <- function(data, sector_name, subplot_label) {
  # Filter data for Services sector
  sector_data <- data %>% filter(broad_sector == sector_name)
  
  # Create the plot
  ggplot(sector_data, aes(x = year)) +
    # Trade line only
    geom_line(aes(y = total_trade, color = "Total Trade"), 
              size = 1.2, alpha = 0.8) +
    geom_point(aes(y = total_trade, color = "Total Trade"), 
               size = 2, alpha = 0.8) +
    
    # Y-axis for trade only
    scale_y_continuous(
      name = "Total Trade (Million USD)",
      labels = comma_format()
    ) +
    
    # Customize colors
    scale_color_manual(values = c("Total Trade" = "steelblue")) +
    
    # Labels
    labs(
      title = paste0("(", subplot_label, ") ", sector_name),
      x = "Year",
      color = NULL
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(color = "steelblue", size = 10),
      axis.text.y = element_text(color = "steelblue")
    ) +
    
    scale_x_continuous(breaks = seq(2000, 2019, 4))
}

# Create individual plots
library(gridExtra)  # Load gridExtra package for grid.arrange

# Plot (a): Agriculture
plot_agri <- create_sector_plot(intra_africa_trade_tariff, "Agriculture", "a")

# Plot (b): Manufacturing  
plot_manu <- create_sector_plot(intra_africa_trade_tariff, "Manufacturing", "b")

# Plot (c): Mining and Energy
plot_mine <- create_sector_plot(intra_africa_trade_tariff, "Mining and Energy", "c")

# Plot (d): Services (trade only)
services_data <- intra_africa_year %>% filter(broad_sector == "Services")
plot_serv <- create_services_plot(intra_africa_year, "Services", "d")

# Save individual plots as PDF
ggsave("output/figures/intra_africa_agriculture_trade_tariff.pdf", 
       plot = plot_agri, 
       width = 8, height = 6, dpi = 300)

ggsave("output/figures/intra_africa_manufacturing_trade_tariff.pdf", 
       plot = plot_manu, 
       width = 8, height = 6, dpi = 300)

ggsave("output/figures/intra_africa_mining_trade_tariff.pdf", 
       plot = plot_mine, 
       width = 8, height = 6, dpi = 300)

ggsave("output/figures/intra_africa_services_trade.pdf", 
       plot = plot_serv, 
       width = 8, height = 6, dpi = 300)

# Create combined plot with all four subplots
combined_plot <- grid.arrange(plot_agri, plot_manu, plot_mine, plot_serv, 
                              nrow = 2, ncol = 2)

# Save combined plot
ggsave("output/figures/intra_africa_trade_and_tariffs_combined.pdf", 
       plot = combined_plot, 
       width = 16, height = 12, dpi = 300)


################################################################################
## Trade flows: RECs Trade -----------------------------------------------------
################################################################################

# Function to create REC trade flow plots
create_rec_plots <- function(rec_to_afc_data, afc_to_rec_data, rec_name, rec_short) {
  
  # Summarize trade by year and broad sector for REC to AfCFTA
  rec_to_afc_year <- rec_to_afc_data %>%
    group_by(year, broad_sector) %>%
    summarise(total_trade = sum(trade, na.rm = TRUE), .groups = "drop")
  
  # Summarize trade by year and broad sector for AfCFTA to REC
  afc_to_rec_year <- afc_to_rec_data %>%
    group_by(year, broad_sector) %>%
    summarise(total_trade = sum(trade, na.rm = TRUE), .groups = "drop")
  
  # Plot (a): REC to AfCFTA\REC
  plot_a <- ggplot(rec_to_afc_year, aes(x = year, y = total_trade, color = broad_sector)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    scale_y_continuous(
      name = "Total Trade (Million USD)",
      labels = comma_format()
    ) +
    scale_color_manual(values = c("Agriculture" = "#2E8B57", 
                                  "Manufacturing" = "#4682B4", 
                                  "Mining and Energy" = "#FF6347",
                                  "Services" = "#9932CC")) +
    labs(
      title = paste0("(a) ", rec_name, " to AfCFTA\\", rec_name),
      x = "Year",
      color = "Broad Sector"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(breaks = seq(2000, 2019, 4))
  
  # Plot (b): AfCFTA\REC to REC
  plot_b <- ggplot(afc_to_rec_year, aes(x = year, y = total_trade, color = broad_sector)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    scale_y_continuous(
      name = "Total Trade (Million USD)",
      labels = comma_format()
    ) +
    scale_color_manual(values = c("Agriculture" = "#2E8B57", 
                                  "Manufacturing" = "#4682B4", 
                                  "Mining and Energy" = "#FF6347",
                                  "Services" = "#9932CC")) +
    labs(
      title = paste0("(b) AfCFTA\\", rec_name, " to ", rec_name),
      x = "Year",
      color = "Broad Sector"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(breaks = seq(2000, 2019, 4))
  
  # Combine plots side by side
  combined_rec_plot <- grid.arrange(plot_a, plot_b, ncol = 2)
  
  # Save combined plot
  ggsave(paste0("output/figures/rec_trade_flows_", tolower(rec_short), ".pdf"), 
         plot = combined_rec_plot, 
         width = 16, height = 8, dpi = 300)
  
  return(combined_rec_plot)
}


# AMU (Arab Maghreb Union)
amu_plot <- create_rec_plots(amu_to_afc, afc_to_amu, "AMU", "AMU")

# CEN-SAD (Community of Sahel–Saharan States)
cen_plot <- create_rec_plots(cen_to_afc, afc_to_cen, "CEN-SAD", "CENSAD")

# COMESA (Common Market for Eastern and Southern Africa)
com_plot <- create_rec_plots(com_to_afc, afc_to_com, "COMESA", "COMESA")

# EAC (East African Community)
eac_plot <- create_rec_plots(eac_to_afc, afc_to_eac, "EAC", "EAC")

# ECCAS (Economic Community of Central African States)
ecc_plot <- create_rec_plots(ecc_to_afc, afc_to_ecc, "ECCAS", "ECCAS")

# ECOWAS (Economic Community of West African States)
eco_plot <- create_rec_plots(eco_to_afc, afc_to_eco, "ECOWAS", "ECOWAS")

# IGAD (Intergovernmental Authority on Development)
iga_plot <- create_rec_plots(iga_to_afc, afc_to_iga, "IGAD", "IGAD")

# SADC (Southern African Development Community)
sad_plot <- create_rec_plots(sad_to_afc, afc_to_sad, "SADC", "SADC")

###################################################################################
# End of the script ---------------------------------------------------------------
###################################################################################