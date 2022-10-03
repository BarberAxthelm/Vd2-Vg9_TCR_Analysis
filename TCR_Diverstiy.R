## TCR_Diversity.R
# This script inputs the data formatted by "TCR_Data_Preperation.R"
# It runs the diversity analyses on the formatted data.
# These are run on TCRD, TCRG, then combined chains.
# It exports all data into CSVs.
# It then counts J-chain usages for individual animals and timepoints and exports to CSV.

# Removes ALL objects from the global environment
rm(list = ls())

# Install and Load Packages
install.packages("tidyverse")
install.packages("data.table")
install.packages("vegan")
install.packages("DescTools")
library(tidyverse)
library(data.table)
library(vegan)
library(DescTools)


# Import Data
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()


## Data Manipulation for Vegan Diversity Analysis ##

# Generate a list of unique CDR3s across all timepoints

# TCRD 
unique_TCRD_names <- macTCRgd %>%
  .$TCRD.CDR3.aa %>%
  unique()

# TCRG
unique_TCRG_names <- macTCRgd %>%
  .$TCRG.CDR3.aa %>%
  unique()

# Combined TCRD/TCRG
unique_Combined_names <- macTCRgd %>%
  .$Combined.CDR3 %>%
  unique()

# Create a list of row names for each of the data.frames
Row_Names <- c("NM11 Day -14", 
               "NM11 Day 4", 
               "NM89 Day -14", 
               "NM89 Day 4",
               "NM251 Day -14", 
               "NM251 Day 4", 
               "NM295 Day -14", 
               "NM295 Day 4"
               )

# Create a data.frame of CDR3 sequences to index:

# TCRD:
unique_TCRD <- as.data.frame(matrix(0, ncol = length(unique_TCRD_names), nrow = length(Row_Names)))
colnames(unique_TCRD) <- unique_TCRD_names
rownames(unique_TCRD) <- Row_Names

# TCRG:
unique_TCRG <- as.data.frame(matrix(0, ncol = length(unique_TCRG_names), nrow = length(Row_Names)))
colnames(unique_TCRG) <- unique_TCRG_names
rownames(unique_TCRG) <- Row_Names

# Combined TCRD/TCRG:
unique_Combined <- as.data.frame(matrix(0, ncol = length(unique_Combined_names), nrow = length(Row_Names)))
colnames(unique_Combined) <- unique_Combined_names
rownames(unique_Combined) <- Row_Names

# Creating a unique list of animal IDs and days to index for rows:

# NM11:
NM11 <- filter(macTCRgd, Animal.ID == "NM11")
NM11_Dn14 <- filter(NM11, Day == -14)
NM11_D4 <- filter(NM11, Day == 4)

# NM89:
NM89 <- filter(macTCRgd, Animal.ID == "NM89")
NM89_Dn14 <- filter(NM89, Day == -14)
NM89_D4 <- filter(NM89, Day == 4)
NM89_D15 <- filter(NM89, Day == 15)

# NM251:
NM251 <- filter(macTCRgd, Animal.ID == "NM251")
NM251_Dn14 <- filter(NM251, Day == -14)
NM251_D4 <- filter(NM251, Day == 4)

# NM295:
NM295 <- filter(macTCRgd, Animal.ID == "NM295")
NM295_Dn14 <- filter(NM295, Day == -14)
NM295_D4 <- filter(NM295, Day == 4)

# Indexing each animalID/time for the unique TCRD

# NM11 Day -14
for (i in 1:nrow(NM11_Dn14)){
  for (j in colnames(unique_TCRD)){
    if (NM11_Dn14$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][1] = unique_TCRD[[j]][1] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][1] references column with name "J" and row 1, which is nm11_dn14
    } 
  }
}

# NM11 Day 4
for (i in 1:nrow(NM11_D4)){
  for (j in colnames(unique_TCRD)){
    if (NM11_D4$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][2] = unique_TCRD[[j]][2] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][2] references column with name "J" and row 2, which is nm11_d4
    } 
  }
}

# NM89 Day -14
for (i in 1:nrow(NM89_Dn14)){
  for (j in colnames(unique_TCRD)){
    if (NM89_Dn14$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][3] = unique_TCRD[[j]][3] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][3] references column with name "J" and row 3, which is nm89_dn14
    } 
  }
}

# NM89 Day 4
for (i in 1:nrow(NM89_D4)){
  for (j in colnames(unique_TCRD)){
    if (NM89_D4$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][4] = unique_TCRD[[j]][4] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][4] references column with name "J" and row 4, which is nm89_d4
    } 
  }
}

# NM251 Day -14
for (i in 1:nrow(NM251_Dn14)){
  for (j in colnames(unique_TCRD)){
    if (NM251_Dn14$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][6] = unique_TCRD[[j]][6] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][5] references column with name "J" and row 5, which is nm251_dn14
    } 
  }
}

# NM251 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(unique_TCRD)){
    if (NM251_D4$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][7] = unique_TCRD[[j]][7] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][6] references column with name "J" and row 6, which is nm251_d4
    } 
  }
}

# NM295 Day -14
for (i in 1:nrow(NM295_Dn14)){
  for (j in colnames(unique_TCRD)){
    if (NM295_Dn14$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][8] = unique_TCRD[[j]][8] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][7] references column with name "J" and row 7, which is nm295_dn14
    } 
  }
}

# NM295 Day 4
for (i in 1:nrow(NM295_D4)){
  for (j in colnames(unique_TCRD)){
    if (NM295_D4$TCRD.CDR3.aa[i] == j){
      unique_TCRD[[j]][9] = unique_TCRD[[j]][9] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][8] references column with name "J" and row 8, which is nm295_d4
    } 
  }
}

# Indexing each animalID/time for the unique TCRG

# NM11 Day -14
for (i in 1:nrow(NM11_Dn14)){
  for (j in colnames(unique_TCRG)){
    if (NM11_Dn14$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][1] = unique_TCRG[[j]][1] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][1] references column with name "J" and row 1, which is nm11_dn14
    } 
  }
}

# NM11 Day 4
for (i in 1:nrow(NM11_D4)){
  for (j in colnames(unique_TCRG)){
    if (NM11_D4$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][2] = unique_TCRG[[j]][2] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][2] references column with name "J" and row 2, which is nm11_d4
    } 
  }
}

# NM89 Day -14
for (i in 1:nrow(NM89_Dn14)){
  for (j in colnames(unique_TCRG)){
    if (NM89_Dn14$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][3] = unique_TCRG[[j]][3] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][3] references column with name "J" and row 3, which is nm89_dn14
    } 
  }
}

# NM89 Day 4
for (i in 1:nrow(NM89_D4)){
  for (j in colnames(unique_TCRG)){
    if (NM89_D4$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][4] = unique_TCRG[[j]][4] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][4] references column with name "J" and row 4, which is nm89_d4
    } 
  }
}

# NM251 Day -14
for (i in 1:nrow(NM251_Dn14)){
  for (j in colnames(unique_TCRG)){
    if (NM251_Dn14$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][6] = unique_TCRG[[j]][6] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][5] references column with name "J" and row 5, which is nm251_dn14
    } 
  }
}

# NM251 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(unique_TCRG)){
    if (NM251_D4$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][7] = unique_TCRG[[j]][7] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][6] references column with name "J" and row 6, which is nm251_d4
    } 
  }
}

# NM295 Day -14
for (i in 1:nrow(NM295_Dn14)){
  for (j in colnames(unique_TCRG)){
    if (NM295_Dn14$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][8] = unique_TCRG[[j]][8] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][7] references column with name "J" and row 7, which is nm295_dn14
    } 
  }
}

# NM295 Day 4
for (i in 1:nrow(NM295_D4)){
  for (j in colnames(unique_TCRG)){
    if (NM295_D4$TCRG.CDR3.aa[i] == j){
      unique_TCRG[[j]][9] = unique_TCRG[[j]][9] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][8] references column with name "J" and row 8, which is nm295_d4
    } 
  }
}

#Indexing each animalID/time for the unique Combined TCRD/TCRG

# NM11 Day -14
for (i in 1:nrow(NM11_Dn14)){
  for (j in colnames(unique_Combined)){
    if (NM11_Dn14$Combined.CDR3[i] == j){
      unique_Combined[[j]][1] = unique_Combined[[j]][1] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][1] references column with name "J" and row 1, which is nm11_dn14
    } 
  }
}

# NM11 Day 4
for (i in 1:nrow(NM11_D4)){
  for (j in colnames(unique_Combined)){
    if (NM11_D4$Combined.CDR3[i] == j){
      unique_Combined[[j]][2] = unique_Combined[[j]][2] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][2] references column with name "J" and row 2, which is nm11_d4
    } 
  }
}

# NM89 Day -14
for (i in 1:nrow(NM89_Dn14)){
  for (j in colnames(unique_Combined)){
    if (NM89_Dn14$Combined.CDR3[i] == j){
      unique_Combined[[j]][3] = unique_Combined[[j]][3] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][3] references column with name "J" and row 3, which is nm89_dn14
    } 
  }
}

# NM89 Day 4
for (i in 1:nrow(NM89_D4)){
  for (j in colnames(unique_Combined)){
    if (NM89_D4$Combined.CDR3[i] == j){
      unique_Combined[[j]][4] = unique_Combined[[j]][4] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][4] references column with name "J" and row 4, which is nm89_d4
    } 
  }
}

# NM251 Day -14
for (i in 1:nrow(NM251_Dn14)){
  for (j in colnames(unique_Combined)){
    if (NM251_Dn14$Combined.CDR3[i] == j){
      unique_Combined[[j]][6] = unique_Combined[[j]][6] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][5] references column with name "J" and row 5, which is nm251_dn14
    } 
  }
}

# NM251 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(unique_Combined)){
    if (NM251_D4$Combined.CDR3[i] == j){
      unique_Combined[[j]][7] = unique_Combined[[j]][7] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][6] references column with name "J" and row 6, which is nm251_d4
    } 
  }
}

# NM295 Day -14
for (i in 1:nrow(NM295_Dn14)){
  for (j in colnames(unique_Combined)){
    if (NM295_Dn14$Combined.CDR3[i] == j){
      unique_Combined[[j]][8] = unique_Combined[[j]][8] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][7] references column with name "J" and row 7, which is nm295_dn14
    } 
  }
}

# NM295 Day 4
for (i in 1:nrow(NM295_D4)){
  for (j in colnames(unique_Combined)){
    if (NM295_D4$Combined.CDR3[i] == j){
      unique_Combined[[j]][9] = unique_Combined[[j]][9] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][8] references column with name "J" and row 8, which is nm295_d4
    } 
  }
}

# Use to check each row sums, which should equal the #of observations for the corresponding animal ID/timepoint data.frame
unique_TCRD_RowSums <- rowSums(unique_TCRD)
unique_TCRG_RowSums <- rowSums(unique_TCRG)
unique_Combined_RowSums <- rowSums(unique_Combined)

# Simpson Diversity Index:
Simpson_TCRD <- as.data.frame(diversity(unique_TCRD, "simpson"))
Simpson_TCRG <- as.data.frame(diversity(unique_TCRG, "simpson"))
Simpson_Combined <- as.data.frame(diversity(unique_Combined, "simpson"))

# Shannon-Weaver Diversity Index:
Shannon_TCRD <- as.data.frame(diversity(unique_TCRD, "shannon"))  ##as.data.frame(diversity(unique_TCRD)) generates the same result
Shannon_TCRG <- as.data.frame(diversity(unique_TCRG, "shannon"))  ##as.data.frame(diversity(unique_TCRG)) generates the same result
Shannon_Combined <- as.data.frame(diversity(unique_Combined, "shannon"))  ##as.data.frame(diversity(unique_Combined)) generates the same result

# Inverse Simpson Diversity Index:
invsimpson_TCRD <- as.data.frame(diversity(unique_TCRD, "invsimpson"))
invsimpson_TCRG <- as.data.frame(diversity(unique_TCRG, "invsimpson"))
invsimpson_Combined <- as.data.frame(diversity(unique_Combined, "invsimpson"))

# Unbiased Simpson
unbias_simpson_TCRD <- as.data.frame(simpson.unb(unique_TCRD))
unbias_simpson_TCRG <- as.data.frame(simpson.unb(unique_TCRG))
unbias_simpson_Combined <- as.data.frame(simpson.unb(unique_Combined))

# Fisher alpha
alpha_TCRD <- as.data.frame(fisher.alpha(unique_TCRD))
alpha_TCRG <- as.data.frame(fisher.alpha(unique_TCRG))
alpha_Combined <- as.data.frame(fisher.alpha(unique_Combined))

# Species richness:
Richness_TCRD <- as.data.frame(specnumber(unique_TCRD)) ## rowSums(unique_TCRD > 0) does the same (per vegan)...
Richness_TCRG <- as.data.frame(specnumber(unique_TCRG)) ## rowSums(unique_TCRG > 0) does the same (per vegan)...
Richness_Combined <- as.data.frame(specnumber(unique_Combined)) ## rowSums(unique_TCRG > 0) does the same (per vegan)...

# Pielou's evenness:
Pielou_TCRD <- as.data.frame(Shannon_TCRD/log(Richness_TCRD))
Pielou_TCRG <- as.data.frame(Shannon_TCRG/log(Richness_TCRG))
Pielou_Combined <- as.data.frame(Shannon_Combined/log(Richness_Combined))

## DescTools ##

# Setup table to add values to for Gini coefficeint and Gini-Simpson coefficent, from Desctools
DescTools_TCRD <- as.data.frame(matrix(0, ncol = 4, nrow = length(Row_Names)))
colnames(DescTools_TCRD) <- c("Gini coefficient", "Gini coefficient lwr.ci", "Gini coefficient upr.ci", "Gini-Simpson coefficient")
rownames(DescTools_TCRD) <- Row_Names

DescTools_TCRG <- as.data.frame(matrix(0, ncol = 4, nrow = length(Row_Names)))
colnames(DescTools_TCRG) <- c("Gini coefficient", "Gini coefficient lwr.ci", "Gini coefficient upr.ci", "Gini-Simpson coefficient")
rownames(DescTools_TCRG) <- Row_Names

DescTools_Combined <- as.data.frame(matrix(0, ncol = 4, nrow = length(Row_Names)))
colnames(DescTools_Combined) <- c("Gini coefficient", "Gini coefficient lwr.ci", "Gini coefficient upr.ci", "Gini-Simpson coefficient")
rownames(DescTools_Combined) <- Row_Names

# Gini coefficient
# Useful resource for the Gini coefficient, particularly with biases: https://www.statsdirect.com/help/nonparametric_methods/gini_coefficient.htm

for (i in 1:nrow(DescTools_TCRD)){
  DescTools_TCRD[i, 1:3] = Gini(as.numeric(unique_TCRD[i,]), conf.level=0.95)
}

for (i in 1:nrow(DescTools_TCRG)){
  DescTools_TCRG[i, 1:3] = Gini(as.numeric(unique_TCRG[i,]), conf.level=0.95)
}

for (i in 1:nrow(DescTools_Combined)){
  DescTools_Combined[i, 1:3] = Gini(as.numeric(unique_Combined[i,]), conf.level=0.95)
}

# Gini-Simpson coefficient:
for (i in 1:nrow(DescTools_TCRD)){
  DescTools_TCRD[i, 4] = GiniSimpson(factor(as.numeric(unique_TCRD[i,])))
}

for (i in 1:nrow(DescTools_TCRG)){
  DescTools_TCRG[i, 4] = GiniSimpson(factor(as.numeric(unique_TCRG[i,])))
}

for (i in 1:nrow(DescTools_Combined)){
  DescTools_Combined[i, 4] = GiniSimpson(factor(as.numeric(unique_Combined[i,])))
}

##################################################################################

Sequence_Counts <- as.data.frame(matrix(0, ncol = 1, nrow = length(Row_Names)))
colnames(Sequence_Counts) <- "Sequence.counts"
rownames(Sequence_Counts) <- Row_Names

Sequence_Counts[1, 1] <- nrow(NM11_Dn14)
Sequence_Counts[2, 1] <- nrow(NM11_D4)
Sequence_Counts[3, 1] <- nrow(NM89_Dn14)
Sequence_Counts[4, 1] <- nrow(NM89_D4)
Sequence_Counts[5, 1] <- nrow(NM89_D15)
Sequence_Counts[6, 1] <- nrow(NM251_Dn14)
Sequence_Counts[7, 1] <- nrow(NM251_D4)
Sequence_Counts[8, 1] <- nrow(NM295_Dn14)
Sequence_Counts[9, 1] <- nrow(NM295_D4)


##################################################################################

### Combine Statistical results into one table to export ###

Stats_table_names <- c("Gini coefficient", 
                       "Gini coefficient- 95% lower CI", 
                       "Gini coefficient- 95% upper CI", 
                       "Gini-Simpson coefficient",  
                       "Simpson Diversity Index", 
                       "Shannon Diversity Index", 
                       "Inverse Simpson Diversity Index", 
                       "Unbiased Simpson Diversity Index", 
                       "Fishers alpha Diversity Index", 
                       "Richness", 
                       "Pielou's Evenness Index",
                       "Sequence count")

# TCRD:
Stats_TCRD <- DescTools_TCRD %>%
  as.data.frame() %>%
  add_column("Simpson Diversity Index" = Simpson_TCRD[, 1]) %>%
  add_column("Shannon Diversity Index" = Shannon_TCRD[, 1]) %>%
  add_column("Inverse Simpson Diversity Index" = invsimpson_TCRD[, 1]) %>%
  add_column("Unbiased Simpson Diversity Index" = unbias_simpson_TCRD[, 1]) %>%
  add_column("Fishers alpha Diversity Index" = alpha_TCRD[, 1]) %>%
  add_column("Richness" = Richness_TCRD[, 1]) %>%
  add_column("Pielou's Evenness Index" = Pielou_TCRD[, 1]) %>%
  add_column("Sequence count" = Sequence_Counts[, 1])
  

colnames(Stats_TCRD) <- Stats_table_names

# TCRG:
Stats_TCRG <- DescTools_TCRG %>%
  as.data.frame() %>%
  add_column("Simpson Diversity Index" = Simpson_TCRG[, 1]) %>%
  add_column("Shannon Diversity Index" = Shannon_TCRG[, 1]) %>%
  add_column("Inverse Simpson Diversity Index" = invsimpson_TCRG[, 1]) %>%
  add_column("Unbiased Simpson Diversity Index" = unbias_simpson_TCRG[, 1]) %>%
  add_column("Fishers alpha Diversity Index" = alpha_TCRG[, 1]) %>%
  add_column("Richness" = Richness_TCRG[, 1]) %>%
  add_column("Pielou's Evenness Index" = Pielou_TCRG[, 1]) %>%
  add_column("Sequence count" = Sequence_Counts[, 1])

colnames(Stats_TCRG) <- Stats_table_names

# Combined TCRD/TCRG:
Stats_Combined <- DescTools_Combined %>%
  as.data.frame() %>%
  add_column("Simpson Diversity Index" = Simpson_Combined[, 1]) %>%
  add_column("Shannon Diversity Index" = Shannon_Combined[, 1]) %>%
  add_column("Inverse Simpson Diversity Index" = invsimpson_Combined[, 1]) %>%
  add_column("Unbiased Simpson Diversity Index" = unbias_simpson_Combined[, 1]) %>%
  add_column("Fishers alpha Diversity Index" = alpha_Combined[, 1]) %>%
  add_column("Richness" = Richness_Combined[, 1]) %>%
  add_column("Pielou's Evenness Index" = Pielou_Combined[, 1]) %>%
  add_column("Sequence count" = Sequence_Counts[, 1])

colnames(Stats_Combined) <- Stats_table_names

#Export stats tables as CSVs to plot data:
setwd("Outputs")
write.csv(Stats_TCRD, "Diversity_Statistics_TCRD.csv")
write.csv(Stats_TCRG, "Diversity_Statistics_TCRG.csv")
write.csv(Stats_Combined, "Diversity_Statistics_Combined.csv")
setwd("..")

#########################################

### J chain usage within Vd2/Vg9 TCRs ###

#TRDJ Chain usage

unique_TRDJ <- c(unique(macTCRgd$TCRD.Best.J.Allele))

TRDJ_usage <- as.data.frame(matrix(0, ncol = length(unique_TRDJ), nrow = length(Row_Names)))
colnames(TRDJ_usage) <- unique_TRDJ
rownames(TRDJ_usage) <- Row_Names

# NM11 Day -14
for (i in 1:nrow(NM11_Dn14)){
  for (j in colnames(TRDJ_usage)){
    if (NM11_Dn14$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][1] = TRDJ_usage[[j]][1] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][1] references column with name "J" and row 1, which is nm11_dn14
    } 
  }
}

# NM11 Day 4
for (i in 1:nrow(NM11_D4)){
  for (j in colnames(TRDJ_usage)){
    if (NM11_D4$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][2] = TRDJ_usage[[j]][2] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][2] references column with name "J" and row 1, which is nm11_d4
    } 
  }
}

# NM89 Day -14
for (i in 1:nrow(NM89_Dn14)){
  for (j in colnames(TRDJ_usage)){
    if (NM89_Dn14$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][3] = TRDJ_usage[[j]][3] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][3] references column with name "J" and row 3, which is nm89_dn14
    } 
  }
}

# NM89 Day 4
for (i in 1:nrow(NM89_D4)){
  for (j in colnames(TRDJ_usage)){
    if (NM89_D4$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][4] = TRDJ_usage[[j]][4] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][4] references column with name "J" and row 4, which is nm89_d4
    } 
  }
}

# NM251 Day -14
for (i in 1:nrow(NM251_Dn14)){
  for (j in colnames(TRDJ_usage)){
    if (NM251_Dn14$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][6] = TRDJ_usage[[j]][6] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][6] references column with name "J" and row 6, which is nm251_dn14
    } 
  }
}

# NM251 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(TRDJ_usage)){
    if (NM251_D4$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][7] = TRDJ_usage[[j]][7] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][7] references column with name "J" and row 7, which is nm251_d4
    } 
  }
}

# NM295 Day -14
for (i in 1:nrow(NM295_Dn14)){
  for (j in colnames(TRDJ_usage)){
    if (NM295_Dn14$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][8] = TRDJ_usage[[j]][8] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][8] references column with name "J" and row 8, which is nm295_dn14
    } 
  }
}

# NM295 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(TRDJ_usage)){
    if (NM295_D4$TCRD.Best.J.Allele[i] == j){
      TRDJ_usage[[j]][9] = TRDJ_usage[[j]][9] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][9] references column with name "J" and row 9, which is nm295_d4
    } 
  }
}

# TRGJ Chain usage

unique_TRGJ <- c(unique(macTCRgd$TCRG.Best.J.Allele))

TRGJ_usage <- as.data.frame(matrix(0, ncol = length(unique_TRGJ), nrow = length(Row_Names)))
colnames(TRGJ_usage) <- unique_TRGJ
rownames(TRGJ_usage) <- Row_Names

# NM11 Day -14
for (i in 1:nrow(NM11_Dn14)){
  for (j in colnames(TRGJ_usage)){
    if (NM11_Dn14$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][1] = TRGJ_usage[[j]][1] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][1] references column with name "J" and row 1, which is nm11_dn14
    } 
  }
}

# NM11 Day 4
for (i in 1:nrow(NM11_D4)){
  for (j in colnames(TRGJ_usage)){
    if (NM11_D4$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][2] = TRGJ_usage[[j]][2] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][2] references column with name "J" and row 1, which is nm11_d4
    } 
  }
}

# NM89 Day -14
for (i in 1:nrow(NM89_Dn14)){
  for (j in colnames(TRGJ_usage)){
    if (NM89_Dn14$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][3] = TRGJ_usage[[j]][3] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][3] references column with name "J" and row 3, which is nm89_dn14
    } 
  }
}

# NM89 Day 4
for (i in 1:nrow(NM89_D4)){
  for (j in colnames(TRGJ_usage)){
    if (NM89_D4$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][4] = TRGJ_usage[[j]][4] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][4] references column with name "J" and row 4, which is nm89_d4
    } 
  }
}

# NM251 Day -14
for (i in 1:nrow(NM251_Dn14)){
  for (j in colnames(TRGJ_usage)){
    if (NM251_Dn14$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][6] = TRGJ_usage[[j]][6] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][6] references column with name "J" and row 6, which is nm251_dn14
    } 
  }
}

# NM251 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(TRGJ_usage)){
    if (NM251_D4$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][7] = TRGJ_usage[[j]][7] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][7] references column with name "J" and row 7, which is nm251_d4
    } 
  }
}

# NM295 Day -14
for (i in 1:nrow(NM295_Dn14)){
  for (j in colnames(TRGJ_usage)){
    if (NM295_Dn14$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][8] = TRGJ_usage[[j]][8] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][8] references column with name "J" and row 8, which is nm295_dn14
    } 
  }
}

# NM295 Day 4
for (i in 1:nrow(NM251_D4)){
  for (j in colnames(TRGJ_usage)){
    if (NM295_D4$TCRG.Best.J.Allele[i] == j){
      TRGJ_usage[[j]][9] = TRGJ_usage[[j]][9] + 1
      # dataframe[[J]] in this loop will give you the entire "J" column
      # Adding dataframe[[J]][9] references column with name "J" and row 9, which is nm295_d4
    } 
  }
}

### Combine TRDJ and TRGJ Chain usages for export

TRJ_usage <- TRDJ_usage %>%
  as.data.frame() %>%
  add_column("BLANK" = 0) 

TRJ_usage <- cbind(TRJ_usage, TRGJ_usage)

# Export stats tables as CSVs to plot data:
# TRJ usage:
setwd("Outputs")
write.csv(TRJ_usage, "TRJ_usage.csv")
setwd("..")
