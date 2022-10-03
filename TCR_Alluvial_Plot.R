## TCR_Alluvial_Plot.R
# This creates Alluvial Plots for all animals for 2 of 3 timepoints, normalized to 100%.
# Final analysis did not include day 15 data. 
# Plots are created for TCRD, TRCG, and paired TCRD-TCRG datasets.
# Additionally it exports a CSV containing counts of each sequence that exists on all days.

# Removes ALL objects from the global environment
rm(list = ls())

# Loads and Installs Packages
install.packages("tidyverse")
install.packages("ggalluvial")
library(tidyverse)
library(ggalluvial)

# Imports data formatted by TCR_Data_Preparation_MiXCR.R Script
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()

setwd("Outputs")

# Animal Selection for Data Analysis
AnimalID <- c("NM251", "NM295", "NM89", "NM11")

# Repeat all steps for every animal
for (n in 1:length(AnimalID)) {

  macTCRgd_AnimalID <- filter(macTCRgd, Animal.ID == AnimalID[n])

  # Generates a data.frame of unique TCRD Junctions from the original macTCRgd data.frame
  CDR3_4 <- bind_rows(
                        data.frame(
                          c(
                            unique(macTCRgd_AnimalID$TCRD.CDR3.aa),
                            unique(macTCRgd_AnimalID$TCRG.CDR3.aa),
                            unique(macTCRgd_AnimalID$Combined.CDR3))))

  colnames(CDR3_4) <- c("Unique_CDR3")

  # Add columns to the TCRD_Junction data frame for TCRD junction count and day
  CDR3_4["delta_count"] <- 0
  CDR3_4["gamma_count"] <- 0
  CDR3_4["combined_count"] <- 0
  CDR3_n14 <- CDR3_4
  #CDR3_15 <- CDR3_4

  CDR3_4["day"] <- "4"
  CDR3_n14["day"] <- "-14"
  #CDR3_15['day'] <- "15"


  # Iterating to count all junctions
  for (i in 1:nrow(macTCRgd_AnimalID)) {

    # Day 4 and -14 are the same length, so this iteration will work for both.
    for (j in 1:nrow(CDR3_4)) {
      if (CDR3_4$`Unique_CDR3`[j] == macTCRgd_AnimalID$TCRD.CDR3.aa[i]) {
      
        # Handling Delta Junctions
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          CDR3_n14$`delta_count`[j] <- CDR3_n14$`delta_count`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4") {
          CDR3_4$`delta_count`[j] <- CDR3_4$`delta_count`[j] + 1
        }
        # if(macTCRgd_AnimalID$Day[i] == "15") {
        #   CDR3_15$`delta_count`[j] <- CDR3_15$`delta_count`[j] + 1
        # }
      }

      # Handling Gamma Junctions
      if (CDR3_4$`Unique_CDR3`[j] == macTCRgd_AnimalID$TCRG.CDR3.aa[i]) {
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          CDR3_n14$`gamma_count`[j] <- CDR3_n14$`gamma_count`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4")         {
          CDR3_4$`gamma_count`[j] <- CDR3_4$`gamma_count`[j] + 1
        }
        # if(macTCRgd_AnimalID$Day[i] == "15") {
        #   CDR3_15$`gamma_count`[j] <- CDR3_15$`gamma_count`[j] + 1
        # }
      }

      # Handling Delta-Gamma Junction
      if (CDR3_4$`Unique_CDR3`[j] == macTCRgd_AnimalID$Combined.CDR3[i]) {
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          CDR3_n14$`combined_count`[j] <- CDR3_n14$`combined_count`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4") {
          CDR3_4$`combined_count`[j] <- CDR3_4$`combined_count`[j] + 1
        }
        # if(macTCRgd_AnimalID$Day[i] == "15") {
        #   CDR3_15$`combined_count`[j] <- CDR3_15$`combined_count`[j] + 1
        # }
      }
    }
  }


  # Sum chain recoveries from TCRD
  # These will be the same for TCRG and TCRD/TCRG Combination
  sum_AnimalID <- c(sum(CDR3_n14$`delta_count`), sum(CDR3_4$`delta_count`)
                #, sum(CDR3_15$`delta_count`)
                )

  # Ensure no division by zero, if there is no data for a specific day.
  for (i in 1:length(sum_AnimalID)) {
    if (sum_AnimalID[i] == 0) {
      sum_AnimalID[i] <- 1
    }
  }

  # Scales the count to format for plots. 
  # Note that CDR3_n14 <- CDR3_4, so they are the same length, so there is no out of bounds access
  for (i in 1:nrow(CDR3_n14)) {
    CDR3_n14$delta_count[i] <- (CDR3_n14$delta_count[i] / sum_AnimalID[1])
    CDR3_4$delta_count[i] <- (CDR3_4$delta_count[i] / sum_AnimalID[2])
    #CDR3_15$delta_count[i] <- (CDR3_15$delta_count[i] / sum_AnimalID[3])

    CDR3_n14$gamma_count[i] <- (CDR3_n14$gamma_count[i] / sum_AnimalID[1])
    CDR3_4$gamma_count[i] <- (CDR3_4$gamma_count[i] / sum_AnimalID[2])
    #CDR3_15$gamma_count[i] <- (CDR3_15$gamma_count[i] / sum_AnimalID[3])

    CDR3_n14$combined_count[i] <- (CDR3_n14$combined_count[i] / sum_AnimalID[1])
    CDR3_4$combined_count[i] <- (CDR3_4$combined_count[i] / sum_AnimalID[2])
    #CDR3_15$combined_count[i] <- (CDR3_15$combined_count[i] / sum_AnimalID[3])
  }

  # Select if you want to bind day 15 or not:
  Unique_CDR3_Data_Single_Animal <- bind_rows(CDR3_n14, CDR3_4)
  #Unique_CDR3_Data_Single_Animal <- bind_rows(CDR3_n14, CDR3_4, CDR3_15)

  # Delta Plot
  ggplot(Unique_CDR3_Data_Single_Animal,
        aes(x = factor(day, level = c("-14","4"
                                    #, "15"
                                    )),
              y = delta_count,
              fill = Unique_CDR3,
              stratum = Unique_CDR3,
              alluvium = Unique_CDR3,
              label = Unique_CDR3)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_flow(color = "grey50") +
    geom_stratum(color = "black") +
    theme(legend.position = "none") +
    labs(x = "Day",
        y = "Percent",
        title = AnimalID[n],
        subtitle = "TCRD CDR3 Usage")

  ggsave(paste(AnimalID[n], "TCRD", "CDR3.tiff", sep = "_"),
        dpi = 1200,
        height = 18,
        width = 10,
        units = "cm",
        limitsize = TRUE)

  # Gamma Plot
  ggplot(Unique_CDR3_Data_Single_Animal,
        aes(x = factor(day, level = c("-14", "4"
                                  #, "15"
                        )),
              y = gamma_count,
              fill = Unique_CDR3,
              stratum = Unique_CDR3,
              alluvium = Unique_CDR3,
              label = Unique_CDR3)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_flow(color = "grey50") +
    geom_stratum(color = "black") +
    theme(legend.position = "none") +
    labs(x = "Day",
        y = "Percent",
        title = AnimalID[n],
        subtitle = "TCRG CDR3 Usage")

  ggsave(paste(AnimalID[n], "TCRG", "CDR3.tiff", sep = "_"),
        dpi = 1200,
        height = 18,
        width = 10,
        units = "cm",
        limitsize = TRUE)

  # Combined Plot
  ggplot(Unique_CDR3_Data_Single_Animal,
        aes(x = factor(day, level = c("-14", "4"
                                  #, "15"
                        )),
              y = combined_count,
              fill = Unique_CDR3,
              stratum = Unique_CDR3,
              alluvium = Unique_CDR3,
              label = Unique_CDR3)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_flow(color = "grey50") +
    geom_stratum(color = "black") +
    theme(legend.position = "none") +
    labs(x = "Day",
        y = "Percent",
        title = AnimalID[n],
        subtitle = "TCRD/TCRG CDR3 Usage")

  ggsave(paste(AnimalID[n], "TCRD-TCRG", "CDR3.tiff", sep = "_"),
        dpi = 1200,
        height = 18,
        width = 10,
        units = "cm",
        limitsize = TRUE)


  # This next section is used to find all chains that occur in both -14 and 4 timepoints
  # This allows manual annotation of the alluvial plots, only paired chains are annotated
  Paired_CDR3_Table <- data.frame(unique(CDR3_4$Unique_CDR3))
  colnames(Paired_CDR3_Table) <- c("Unique_CDR3")

  # Add columns to the TCRD_Junction data frame for TCRD junction count and day
  Paired_CDR3_Table["day_n14"] <- 0
  Paired_CDR3_Table["day_4"] <- 0
  Paired_CDR3_Table["day_15"] <- 0

  data.frame(Paired_CDR3_Table)

  # Iterating to tally counts for all junctions
  for (i in 1:nrow(macTCRgd_AnimalID)) {
    for (j in 1:nrow(Paired_CDR3_Table)) {

      # Handling Delta Junctions
      if (Paired_CDR3_Table$`Unique_CDR3`[j] == macTCRgd_AnimalID$TCRD.CDR3.aa[i]) {
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          Paired_CDR3_Table$`day_n14`[j] <- Paired_CDR3_Table$`day_n14`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4") {
          Paired_CDR3_Table$`day_4`[j] <- Paired_CDR3_Table$`day_4`[j] + 1
        }
        #if (macTCRgd_AnimalID$Day[i] == "15") {
        # Paired_CDR3_Table$`day_15`[j] <- Paired_CDR3_Table$`day_15`[j] + 1
        #}
      }

      # Handling Gamma Junctions
      if (Paired_CDR3_Table$`Unique_CDR3`[j] == macTCRgd_AnimalID$TCRG.CDR3.aa[i]) {
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          Paired_CDR3_Table$`day_n14`[j] <- Paired_CDR3_Table$`day_n14`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4") {
          Paired_CDR3_Table$`day_4`[j] <- Paired_CDR3_Table$`day_4`[j] + 1
        }
        #if (macTCRgd_AnimalID$Day[i] == "15") {
        # Paired_CDR3_Table$`day_15`[j] <- Paired_CDR3_Table$`day_15`[j] + 1
        #}
      }

      # Handling Delta-Gamma Junction
      if (Paired_CDR3_Table$`Unique_CDR3`[j] == macTCRgd_AnimalID$Combined.CDR3[i]) {
        if (macTCRgd_AnimalID$Day[i] == "-14") {
          Paired_CDR3_Table$`day_n14`[j] <- Paired_CDR3_Table$`day_n14`[j] + 1
        }
        if (macTCRgd_AnimalID$Day[i] == "4") {
          Paired_CDR3_Table$`day_4`[j] <- Paired_CDR3_Table$`day_4`[j] + 1
        }
        #if (macTCRgd_AnimalID$Day[i] == "15") {
        # Paired_CDR3_Table$`day_15`[j] <- Paired_CDR3_Table$`day_15`[j] + 1
        #}
      }
    }
  }

  # Removes chains that are not temporally paired
  Paired_CDR3_Table <- filter(Paired_CDR3_Table,
                day_n14 != 0,
                day_4 != 0)

  write.csv(Paired_CDR3_Table, file = paste(AnimalID[n],
                                            "Paired",
                                            "CDR3",
                                            "Table.csv",
                                            sep = "-"))
}

setwd("..")