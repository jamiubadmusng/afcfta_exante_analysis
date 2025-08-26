################################################################################
# Combine Welfare Effects from Different Sectors
# This script combines welfare effects (rGDP_full_ch) from multiple sectors:
# Agriculture, Manufacturing, Mining and Energy, Services, and Structural Gravity
# Author: Jamiu Olamilekan Badmus
################################################################################

# clear workspace
rm(list = ls())

# Check and install required packages if they're not already installed
#required_packages <- c("data.table", "dplyr", "xtable", "countrycode")
#new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
#if(length(new_packages) > 0) {
#  cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
#  install.packages(new_packages, repos="https://cran.rstudio.com/")
#}

# Load required packages
library(data.table)
library(dplyr)
library(xtable)
library(countrycode)

# Set working directory to load RDS files
setwd("C:/Users/muham/EGEI Dissertation/afcfta/output/rds")

# Step 1: Load all datasets
agri_data <- readRDS("all_indexes_with_row_geppml_agri.rds")
manu_data <- readRDS("all_indexes_with_row_geppml_manu.rds")
mine_data <- readRDS("all_indexes_with_row_geppml_mine.rds")
serv_data <- readRDS("all_indexes_with_row_geppml_serv.rds")
struc_data <- readRDS("all_indexes_with_row_geppml_struc.rds")

# Step 2: Extract relevant columns from each dataset
agri_welfare <- agri_data[, .(country_iso_3 = country, rGDP_full_ch_agri = rGDP_full_ch)]
manu_welfare <- manu_data[, .(country_iso_3 = country, rGDP_full_ch_manu = rGDP_full_ch)]
mine_welfare <- mine_data[, .(country_iso_3 = country, rGDP_full_ch_mine = rGDP_full_ch)]
serv_welfare <- serv_data[, .(country_iso_3 = country, rGDP_full_ch_serv = rGDP_full_ch)]
struc_welfare <- struc_data[, .(country_iso_3 = country, rGDP_full_ch_struc = rGDP_full_ch)]

# Step 3: Merge all datasets
combined_welfare <- agri_welfare
combined_welfare <- merge(combined_welfare, manu_welfare, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, mine_welfare, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, serv_welfare, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, struc_welfare, by = "country_iso_3", all = TRUE)


#Step 4: Add country names
# Create a function to handle special cases (ROW and countries not in countrycode package)
# Also providing a fallback if countrycode package isn't available
get_country_name <- function(iso3) {
  if (iso3 == "ROW") {
    return("Rest of World")
  } else if (iso3 == "SSD") {
    return("South Sudan")
  } else if (iso3 == "STP") {
    return("São Tomé and Príncipe")
  } else if (iso3 == "COD") {
    return("Democratic Republic of the Congo")
  } else if (iso3 == "COG") {
    return("Republic of the Congo")
  } else {
    # Try using countrycode package, but provide fallback
    tryCatch({
      name <- countrycode(iso3, "iso3c", "country.name")
      if (is.na(name)) {
        return(iso3) # Return the ISO code if no name is found
      } else {
        return(name)
      }
    }, error = function(e) {
      # If countrycode package isn't working, use a manual lookup table for common countries
      country_lookup <- list(
        "DZA" = "Algeria", "AGO" = "Angola", "BEN" = "Benin", "BWA" = "Botswana",
        "BFA" = "Burkina Faso", "BDI" = "Burundi", "CPV" = "Cape Verde", "CMR" = "Cameroon",
        "CAF" = "Central African Republic", "TCD" = "Chad", "COM" = "Comoros", "COG" = "Congo",
        "CIV" = "Côte d'Ivoire", "COD" = "Democratic Republic of the Congo", "DJI" = "Djibouti",
        "EGY" = "Egypt", "GNQ" = "Equatorial Guinea", "SWZ" = "Eswatini", "ETH" = "Ethiopia",
        "GAB" = "Gabon", "GMB" = "Gambia", "GHA" = "Ghana", "GIN" = "Guinea", "GNB" = "Guinea-Bissau",
        "KEN" = "Kenya", "LSO" = "Lesotho", "LBR" = "Liberia", "LBY" = "Libya", "MDG" = "Madagascar",
        "MWI" = "Malawi", "MLI" = "Mali", "MRT" = "Mauritania", "MUS" = "Mauritius", "MAR" = "Morocco",
        "MOZ" = "Mozambique", "NAM" = "Namibia", "NER" = "Niger", "NGA" = "Nigeria", "RWA" = "Rwanda",
        "SEN" = "Senegal", "SYC" = "Seychelles", "SLE" = "Sierra Leone", "SOM" = "Somalia",
        "ZAF" = "South Africa", "SDN" = "Sudan", "TZA" = "Tanzania", "TGO" = "Togo", "TUN" = "Tunisia",
        "UGA" = "Uganda", "ZMB" = "Zambia", "ZWE" = "Zimbabwe", "DEU" = "Germany"
      )
      
      if (iso3 %in% names(country_lookup)) {
        return(country_lookup[[iso3]])
      } else {
        return(iso3) # Return the ISO code if not in our lookup
      }
    })
  }
}

# Apply the function to get country names
combined_welfare[, country_name := sapply(country_iso_3, get_country_name)]

# Step 5: Reorder columns as specified
combined_welfare <- combined_welfare[, .(
  country_iso_3,
  country_name,
  rGDP_full_ch_agri,
  rGDP_full_ch_manu,
  rGDP_full_ch_mine,
  rGDP_full_ch_serv,
  rGDP_full_ch_struc
)]

# Step 6: Sort by alphabetical order of country names, with ROW at the end
# First, create a temporary sorting variable
combined_welfare[, sort_order := ifelse(country_iso_3 == "ROW", 2, 1)]
# Sort by sort_order and then by country_name
setorder(combined_welfare, sort_order, country_name)
# Remove the temporary sorting variable
combined_welfare[, sort_order := NULL]

# Step 7: Format values for better display (rounded to 2 decimal places)
format_columns <- c("rGDP_full_ch_agri", "rGDP_full_ch_manu", "rGDP_full_ch_mine", 
                    "rGDP_full_ch_serv", "rGDP_full_ch_struc")

for (col in format_columns) {
  combined_welfare[, (col) := round(get(col), 2)]
}

# Save the combined dataset
saveRDS(combined_welfare, "combined_welfare_effects.rds")
write.csv(combined_welfare, "C:/Users/muham/EGEI Dissertation/afcfta/output/csv/combined_welfare_effects.csv", row.names = FALSE)

# Step 8: Generate LaTeX table
# Make a copy for formatting
latex_table <- copy(combined_welfare)

# Format numbers to include % sign and ensure proper decimal places
for (col in format_columns) {
  latex_table[, (col) := sprintf("%.2f\\%%", get(col))]
}

# Create xtable object with appropriate formatting
x_table <- xtable(latex_table, 
                  caption = "Welfare Effects of AfCFTA Agreement Across Different Sectors (\\% Change in Real GDP)",
                  label = "tab:welfare")

# Add a caption note
caption_note <- "Notes: This table shows the welfare effects (percentage change in real GDP) from the AfCFTA agreement across different sectors. AGRI = Agriculture sector, MANU = Manufacturing sector, MINE = Mining and Energy sector, SERV = Services sector, STRUC = Manufacturing from Structural Gravity Database. All estimates are based on ITPD-E-R02 data except for STRUC column. ROW = Rest of World (weighted average for non-AfCFTA countries)."

# Generate LaTeX code
latex_code <- print(x_table, 
                    include.rownames = FALSE,
                    sanitize.text.function = function(x) x,
                    caption.placement = "top",
                    hline.after = c(-1, 0, nrow(latex_table)),
                    add.to.row = list(pos = list(nrow(latex_table)),
                                     command = "\\hline "),
                    floating = TRUE,
                    table.placement = "htbp",
                    tabular.environment = "tabular",
                    size = "\\small")

# Add the caption note to the LaTeX code
latex_code_with_note <- gsub("\\\\end\\{table\\}", 
                            paste0("\\\\caption*{", caption_note, "}\n\\\\end{table}"), 
                            latex_code)

# Rename columns for better LaTeX presentation
latex_code_with_note <- gsub("country\\_iso\\_3", "ISO-3", latex_code_with_note)
latex_code_with_note <- gsub("country\\_name", "Country", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_agri", "AGRI", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_manu", "MANU", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_mine", "MINE", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_serv", "SERV", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_struc", "STRUC", latex_code_with_note)

# Write the LaTeX code to a file
writeLines(latex_code_with_note, "C:/Users/muham/EGEI Dissertation/afcfta/output/table/welfare_effects_table.tex")

# Print completion message
cat("\nAnalysis completed. Results saved to:\n")
cat("1. combined_welfare_effects.rds (R data format)\n")
cat("2. combined_welfare_effects.csv (CSV format)\n")
cat("3. welfare_effects_table.tex (LaTeX table)\n")

# Print a preview of the combined data
cat("\nPreview of combined welfare effects:\n")
print(head(combined_welfare, 10))


################################################################################
# End of Script ----------------------------------------------------------------
################################################################################