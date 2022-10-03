## TCRG_CDR3_Shared_Clonotypes.R
# This script identifies clonoltypes  shared between animals and timepoints for TCRG.
# The output is a table which gives the CDR3 clonoltype and the count for each animal_timepoint sample
# For example:
# CDR3             | NM11_Dn14 | NM11_D4 | NM89_Dn14 | ... | shared
# CALWEVQQFGRKVKLF | 11        | 10      | 4         | ... | 6

# Removes ALL objects from the global environment
rm(list = ls())

# Installs and loads Packages
install.packages("tidyverse")
install.packages("data.table")
library(tidyverse)
library(data.table)

# Imports data formatted by TCR_Data_Preparation_MiXCR.R Script
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()

# Generates a list of unique clonotypes for the entire dataset for TCRG CDR3 cells
unique_tcrg_cdr3 <- data.frame(unique(macTCRgd$TCRG.CDR3.aa))
colnames(unique_tcrg_cdr3) <- c("CDR3")

# Setup columns for each animal and timepoint to keep counts for each clonotype
unique_tcrg_cdr3["NM11_Dn14"] <- 0
unique_tcrg_cdr3["NM11_D4"] <- 0

unique_tcrg_cdr3["NM89_Dn14"] <- 0
unique_tcrg_cdr3["NM89_D4"] <- 0

unique_tcrg_cdr3["NM251_Dn14"] <- 0
unique_tcrg_cdr3["NM251_D4"] <- 0

unique_tcrg_cdr3["NM295_Dn14"] <- 0
unique_tcrg_cdr3["NM295_D4"] <- 0

unique_tcrg_cdr3["shared"] <- 0


# Identify shared values via brute force.
# For each row (identified clonotype) in macTCRgd this will:
#    - find the appropriate row on the unique clonotype dataframe
#    - tally the correct column, determined by animal/day
#    - determine how many animal/day samples share that clonotype

for (i in 1:nrow(macTCRgd)) {
  for (j in 1:nrow(unique_tcrg_cdr3)) {
    if (macTCRgd$TCRG.CDR3.aa[i] == unique_tcrg_cdr3$CDR3[j]) {

      # NM11 Count
      if (macTCRgd$Animal.ID[i] == "NM11") {
        if (macTCRgd$Day[i] == -14) {
          unique_tcrg_cdr3$NM11_Dn14[j] <- unique_tcrg_cdr3$NM11_Dn14[j] + 1
          if (unique_tcrg_cdr3$NM11_Dn14[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
  
        if (macTCRgd$Day[i] == 4) {
          unique_tcrg_cdr3$NM11_D4[j] <- unique_tcrg_cdr3$NM11_D4[j] + 1
          if (unique_tcrg_cdr3$NM11_D4[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
      }

      # NM89 Count
      if (macTCRgd$Animal.ID[i] == "NM89") {
        if (macTCRgd$Day[i] == -14) {
          unique_tcrg_cdr3$NM89_Dn14[j] <- unique_tcrg_cdr3$NM89_Dn14[j] + 1
          if (unique_tcrg_cdr3$NM89_Dn14[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
        if (macTCRgd$Day[i] == 4) {
          unique_tcrg_cdr3$NM89_D4[j] <- unique_tcrg_cdr3$NM89_D4[j] + 1
          if (unique_tcrg_cdr3$NM89_D4[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
      }

      # NM251 Count
      if (macTCRgd$Animal.ID[i] == "NM251") {
        if (macTCRgd$Day[i] == -14) {
          unique_tcrg_cdr3$NM251_Dn14[j] <- unique_tcrg_cdr3$NM251_Dn14[j] + 1
          if (unique_tcrg_cdr3$NM251_Dn14[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }

        if (macTCRgd$Day[i] == 4) {
          unique_tcrg_cdr3$NM251_D4[j] <- unique_tcrg_cdr3$NM251_D4[j] + 1
          if (unique_tcrg_cdr3$NM251_D4[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
      }

      # NM295 Count
      if (macTCRgd$Animal.ID[i] == "NM295") {
        if (macTCRgd$Day[i] == -14) {
          unique_tcrg_cdr3$NM295_Dn14[j] <- unique_tcrg_cdr3$NM295_Dn14[j] + 1
          if (unique_tcrg_cdr3$NM295_Dn14[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }

        if (macTCRgd$Day[i] == 4) {
          unique_tcrg_cdr3$NM295_D4[j] <- unique_tcrg_cdr3$NM295_D4[j] + 1
          if (unique_tcrg_cdr3$NM295_D4[j] == 1) {
            unique_tcrg_cdr3$shared[j] <- unique_tcrg_cdr3$shared[j] + 1
          }
        }
      }
    }
  }
}

# Sort the results by how many animal/day timepoints share that clonotype
unique_tcrg_cdr3 <- unique_tcrg_cdr3[order(-unique_tcrg_cdr3$shared, unique_tcrg_cdr3$CDR3), ]

# Set the working directory and save outputs
setwd("Outputs")
write.csv(unique_tcrg_cdr3, "TCRG_CDR3_Shared_Clonotypes.csv")
setwd("..")

