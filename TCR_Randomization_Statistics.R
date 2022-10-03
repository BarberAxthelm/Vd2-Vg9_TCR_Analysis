## TCR_Randomization_Statistics.R
# This script needs manual edits by the user to run different animals and TCR chains
# This script inputs the data formatted by "TCR_Data_Preperation.R" and creates a 
# randomized list of data to run statistics 
# These statistics are based recommendations from:
#    https://www.sciencedirect.com/science/article/abs/pii/S0022175907000397?via%3Dihub
# Note that all stats are run from a function, which is run in parallel to reduce the time
# and how many instances is function of the cores on the machine it is run from.
# This file contains the function, and code used to run the function (below it).


# Removes ALL objects from the global environment
rm(list = ls())

# Installs and Loads Packages
install.packages("tidyverse")
install.packages("data.table")
install.packages("vegan")
install.packages("DescTools")
install.packages("dplyr")
install.packages("parallel")
library(tidyverse)
library(data.table)
library(vegan)
library(DescTools)
library(dplyr)
library(parallel)

##################################################################################

# Function Calculations
RandomizationTestStatistics <- function(k, df) {
  n <- 1

  # Setup data.frame of output statistics for the function
  TEST_RandStats_TCR <- as.data.frame(matrix(0, ncol = 21, nrow = n))
  colnames(TEST_RandStats_TCR) <- c("Random.Simpson.Dn14", #1
                                     "Random.Simpson.D4", #2
                                     "Random.Simpson.Difference", #3
                                     "Random.Richness.Dn14", #4
                                     "Random.Richness.D4", #5
                                     "Random.Richness.Difference", #6
                                     "Random.Shannon.Dn14", #7
                                     "Random.Shannon.D4", #8
                                     "Random.Shannon.Difference", #9
                                     "Random.invSimpson.Dn14", #10
                                     "Random.invSimpson.D4", #11
                                     "Random.invSimpson.Difference", #12
                                     "Random.unbiasSimpson.Dn14", #13
                                     "Random.unbiasSimpson.D4", #14
                                     "Random.unbiasSimpson.Difference", #15
                                     "Random.FisherAlpha.Dn14", #16
                                     "Random.FisherAlpha.D4", #17
                                     "Random.FisherAlpha.Difference", #18
                                     "Random.Pielou.Dn14", #19
                                     "Random.Pielou.D4", #20
                                     "Random.Pielou.Difference" #21
                                     )


  TEST_unique_TCR <- as.data.frame(matrix(0, ncol = length(c(unique(df[,2]))), nrow = 2))
  colnames(TEST_unique_TCR) <- c(unique(df[, 2]))
  rownames(TEST_unique_TCR) <- c("Day -14",
                                  "Day 4")

  # Tallying unique TCRs for each day
  for (i in 1:nrow(df)) {
    for (j in colnames(TEST_unique_TCR)) {
      if (df[i, k + 1] == j) {
        if (df[i, 1] == -14) {
          TEST_unique_TCR[[j]][1] <- TEST_unique_TCR[[j]][1] + 1
        }
        if (df[i, 1] == 4) {
          TEST_unique_TCR[[j]][2] <- TEST_unique_TCR[[j]][2] + 1
        }
      }
    }
  }

  k <- n

  TEST_RandStats_TCR$Random.Simpson.Dn14[k] <- diversity(TEST_unique_TCR[1,], "simpson")
  TEST_RandStats_TCR$Random.Simpson.D4[k] <- diversity(TEST_unique_TCR[2,], "simpson")
  TEST_RandStats_TCR$Random.Simpson.Difference[k] <- abs((TEST_RandStats_TCR$Random.Simpson.Dn14[k]) - (TEST_RandStats_TCR$Random.Simpson.D4[k]))

  TEST_RandStats_TCR$Random.Richness.Dn14[k] <- specnumber(TEST_unique_TCR[1,])
  TEST_RandStats_TCR$Random.Richness.D4[k] <- specnumber(TEST_unique_TCR[2,])
  TEST_RandStats_TCR$Random.Richness.Difference[k] <- abs((TEST_RandStats_TCR$Random.Richness.Dn14[k]) - (TEST_RandStats_TCR$Random.Richness.D4[k]))

  TEST_RandStats_TCR$Random.Shannon.Dn14[k] = diversity(TEST_unique_TCR[1,], "shannon")
  TEST_RandStats_TCR$Random.Shannon.D4[k] = diversity(TEST_unique_TCR[2,], "shannon")
  TEST_RandStats_TCR$Random.Shannon.Difference[k] = abs(TEST_RandStats_TCR$Random.Shannon.Dn14[k] - TEST_RandStats_TCR$Random.Shannon.D4[k] )

  TEST_RandStats_TCR$Random.invSimpson.Dn14[k] = diversity(TEST_unique_TCR[1,], "invsimpson")
  TEST_RandStats_TCR$Random.invSimpson.D4[k] = diversity(TEST_unique_TCR[2,], "invsimpson")
  TEST_RandStats_TCR$Random.invSimpson.Difference[k] = abs(TEST_RandStats_TCR$Random.invSimpson.Dn14[k] - TEST_RandStats_TCR$Random.invSimpson.D4[k])

  TEST_RandStats_TCR$Random.unbiasSimpson.Dn14[k] = simpson.unb(TEST_unique_TCR[1,])
  TEST_RandStats_TCR$Random.unbiasSimpson.D4[k] = simpson.unb(TEST_unique_TCR[2,])
  TEST_RandStats_TCR$Random.unbiasSimpson.Difference[k] = abs(TEST_RandStats_TCR$Random.unbiasSimpson.Dn14[k] - TEST_RandStats_TCR$Random.unbiasSimpson.D4[k])

  TEST_RandStats_TCR$Random.FisherAlpha.Dn14[k] = fisher.alpha(TEST_unique_TCR[1,])
  TEST_RandStats_TCR$Random.FisherAlpha.D4[k] = fisher.alpha(TEST_unique_TCR[2,])
  TEST_RandStats_TCR$Random.FisherAlpha.Difference[k] = abs(TEST_RandStats_TCR$Random.FisherAlpha.Dn14[k] - TEST_RandStats_TCR$Random.FisherAlpha.D4[k])

  TEST_RandStats_TCR$Random.Pielou.Dn14[k] = (diversity(TEST_unique_TCR[1,], "shannon")) / (log(specnumber(TEST_unique_TCR[1,])))
  TEST_RandStats_TCR$Random.Pielou.D4[k] = (diversity(TEST_unique_TCR[2,], "shannon")) / (log(specnumber(TEST_unique_TCR[2,])))
  TEST_RandStats_TCR$Random.Pielou.Difference[k] = abs(TEST_RandStats_TCR$Random.Pielou.Dn14[k] - TEST_RandStats_TCR$Random.Pielou.D4[k])

  return(TEST_RandStats_TCR)

}


##################################################################################
# Select which animal to do the analysis on here, by uncommenting appropriating animal
AnimalID <- "NM251"
#AnimalID <- "NM295"
#AnimalID <- "NM89"
#AnimalID <- "NM11"

# Look below the function, for "# !Change" Comments for manual changes
# Option to change between TCRD/TCRG
# Must change the filename!
##################################################################################

# Imports data formatted by TCR_Data_Preparation_MiXCR.R Script
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()



## Creates Randomization values between day -14 and day 4 ##

# Filter data to animal and remove day 15 samples, if present
NM_Data <- filter(macTCRgd, Animal.ID == AnimalID)
NM_Data <- filter(mNM_Data, Day != 15)

# Sets up dataframe for randomized
TEST_NM_TCR <- as.data.frame(matrix(0, ncol = 10001, nrow = nrow(NM_Data)))
colnames(TEST_NM_TCR) <- c(0:10000)
colnames(TEST_NM_TCR)[1] <- "Day"

TEST_NM_TCR$Day <- NM_Data$Day

# Sample function randomizes
for (i in 2:10001) {
  # !Change TCRG/TCRD here!
  TEST_NM_TCR[i] <- c(sample(NM_Data$TCRD.CDR3.aa))
}

## Running the functions in parallel to speed up analysis process. ##
start_time <- Sys.time()
test_results <- mclapply(1:10000, RandomizationTestStatistics, df = TEST_NM_TCR, mc.cores = (detectCores()-2), mc.preschedule = TRUE)
end_time <- Sys.time()

# Keeping track of elapsed time
elapsed_time <- end_time - start_time
print(elapsed_time)

bound_outputs <- do.call("rbind", test_results)


## Save Output to File ##

setwd("Outputs")
# !Change Filename here
write.csv(bound_outputs, "TEST_NM251_TCRD_Randomization_Analysis.csv")
setwd("..")
