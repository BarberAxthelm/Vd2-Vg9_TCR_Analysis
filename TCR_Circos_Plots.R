## TCR_Circos_Plots.R
# This script inputs the data formatted by "TCR_Data_Preperation.R"
# It formats the data to create circos plots and normalized circos plots.
# Plots are created for TCRD, TCRG, and TCRD-TCRG. 
# It is run for a single animal. To run for the next animal:
# update animal ID and re-run.

# Removes ALL objects from the global environment
rm(list = ls())

# Install and Loads Packages
install.packages("tidyverse")
install.packages("ggalluvial")
install.packages("circlize")
install.packages("data.table")
library(tidyverse)
library(ggalluvial)
library(circlize)
library(data.table)

# Import data
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()

### Start of Analysis ###
setwd("Outputs")

####################################################
# Animal Selection for Data Analysis
# Only uncomment 1 animal at a time, select animal here before re-running script

AnimalID <- "NM251"
#AnimalID <- "NM295"
#AnimalID <- "NM89"
#AnimalID <- "NM11"

macTCRgd_AnimalID <- filter(macTCRgd, Animal.ID == AnimalID)
####################################################

# Circos plots of CDR3 usage #

# Find all unique CDR3 for individual days
macTCRgd_AnimalID_4 <- filter(macTCRgd_AnimalID, Day == 4)
macTCRgd_AnimalID_n14 <- filter(macTCRgd_AnimalID, Day == -14)
#macTCRgd_AnimalID_15 <- filter(macTCRgd_AnimalID, Day == 15)

# TCRD:
TCRD_CDR3_4 <- data.frame(unique(macTCRgd_AnimalID_4$TCRD.CDR3.aa))
TCRD_CDR3_n14 <- data.frame(unique(macTCRgd_AnimalID_n14$TCRD.CDR3.aa))
#TCRD_CDR3_n15 <- data.frame(unique(macTCRgd_AnimalID_n15$TCRD.CDR3.aa))

# TCRG:
TCRG_CDR3_4 <- data.frame(unique(macTCRgd_AnimalID_4$TCRG.CDR3.aa))
TCRG_CDR3_n14 <- data.frame(unique(macTCRgd_AnimalID_n14$TCRG.CDR3.aa))
#TCRG_CDR3_15 <- data.frame(unique(macTCRgd_AnimalID_15$TCRG.CDR3.aa))

# Combined CDR3s:
Combined_CDR3_4 <- data.frame(unique(macTCRgd_AnimalID_4$Combined.CDR3))
Combined_CDR3_n14 <- data.frame(unique(macTCRgd_AnimalID_n14$Combined.CDR3))
#Combined_CDR3_15 <- data.frame(unique(macTCRgd_AnimalID_15$Combined.CDR3))

# TCRD:
colnames(TCRD_CDR3_4) <- c("CDR3")
colnames(TCRD_CDR3_n14) <- c("CDR3")
#colnames(TCRD_CDR3_n15) <- c("CDR3")

# TCRG:
colnames(TCRG_CDR3_4) <- c("CDR3")
colnames(TCRG_CDR3_n14) <- c("CDR3")
#colnames(TCRG_CDR3_15) <- c("CDR3")

# Combined CDR3s:
colnames(Combined_CDR3_4) <- c("CDR3")
colnames(Combined_CDR3_n14) <- c("CDR3")
#colnames(Combined_CDR3_15) <- c("CDR3")

# Separating shared CDR3 across days

# TCRD:
shared_TCRD_CDR3 <- data.frame(intersect(TCRD_CDR3_4$CDR3, TCRD_CDR3_n14$CDR3))
unique_TCRD_CDR3_day4 <- data.frame(setdiff(TCRD_CDR3_4$CDR3, TCRD_CDR3_n14$CDR3))
unique_TCRD_CDR3_dayn14 <- data.frame(setdiff(TCRD_CDR3_n14$CDR3, TCRD_CDR3_4$CDR3))
colnames(shared_TCRD_CDR3) <- c("CDR3")
colnames(unique_TCRD_CDR3_day4) <- c("CDR3")
colnames(unique_TCRD_CDR3_dayn14) <- c("CDR3")

# TCRG:
shared_TCRG_CDR3 <- data.frame(intersect(TCRG_CDR3_4$CDR3, TCRG_CDR3_n14$CDR3))
unique_TCRG_CDR3_day4 <- data.frame(setdiff(TCRG_CDR3_4$CDR3, TCRG_CDR3_n14$CDR3))
unique_TCRG_CDR3_dayn14 <- data.frame(setdiff(TCRG_CDR3_n14$CDR3, TCRG_CDR3_4$CDR3))
colnames(shared_TCRG_CDR3) <- c("CDR3")
colnames(unique_TCRG_CDR3_day4) <- c("CDR3")
colnames(unique_TCRG_CDR3_dayn14) <- c("CDR3")

# Combined CDR3s:
shared_Combined_CDR3 <- data.frame(intersect(Combined_CDR3_4$CDR3, Combined_CDR3_n14$CDR3))
unique_Combined_CDR3_day4 <- data.frame(setdiff(Combined_CDR3_4$CDR3, Combined_CDR3_n14$CDR3))
unique_Combined_CDR3_dayn14 <- data.frame(setdiff(Combined_CDR3_n14$CDR3, Combined_CDR3_4$CDR3))
colnames(shared_Combined_CDR3) <- c("CDR3")
colnames(unique_Combined_CDR3_day4) <- c("CDR3")
colnames(unique_Combined_CDR3_dayn14) <- c("CDR3")

# Adding to/from for the circos plot
# From Day -14, to Day 4

# TCRD
shared_TCRD_CDR3["group_from"] <- "Day_n14"
shared_TCRD_CDR3["group_to"] <- "Day_4"

unique_TCRD_CDR3_day4["group_from"] <- "Day_4"
unique_TCRD_CDR3_day4["group_to"] <- "Day_4"

unique_TCRD_CDR3_dayn14["group_from"] <- "Day_n14"
unique_TCRD_CDR3_dayn14["group_to"] <- "Day_n14"

# TCRG
shared_TCRG_CDR3["group_from"] <- "Day_n14"
shared_TCRG_CDR3["group_to"] <- "Day_4"

unique_TCRG_CDR3_day4["group_from"] <- "Day_4"
unique_TCRG_CDR3_day4["group_to"] <- "Day_4"

unique_TCRG_CDR3_dayn14["group_from"] <- "Day_n14"
unique_TCRG_CDR3_dayn14["group_to"] <- "Day_n14"

# Combined CDR3s :
shared_Combined_CDR3["group_from"] <- "Day_n14"
shared_Combined_CDR3["group_to"] <- "Day_4"

unique_Combined_CDR3_day4["group_from"] <- "Day_4"
unique_Combined_CDR3_day4["group_to"] <- "Day_4"

unique_Combined_CDR3_dayn14["group_from"] <- "Day_n14"
unique_Combined_CDR3_dayn14["group_to"] <- "Day_n14"

# Setup Counts to weight the plot

# TCRD
shared_TCRD_CDR3["count_from"] <- 0
shared_TCRD_CDR3["count_to"] <- 0
unique_TCRD_CDR3_day4["count_from"] <- 0
unique_TCRD_CDR3_day4["count_to"] <- 0
unique_TCRD_CDR3_dayn14["count_from"] <- 0
unique_TCRD_CDR3_dayn14["count_to"] <- 0

# TCRG
shared_TCRG_CDR3["count_from"] <- 0
shared_TCRG_CDR3["count_to"] <- 0
unique_TCRG_CDR3_day4["count_from"] <- 0
unique_TCRG_CDR3_day4["count_to"] <- 0
unique_TCRG_CDR3_dayn14["count_from"] <- 0
unique_TCRG_CDR3_dayn14["count_to"] <- 0

# Combined CDR3s :
shared_Combined_CDR3["count_from"] <- 0
shared_Combined_CDR3["count_to"] <- 0
unique_Combined_CDR3_day4["count_from"] <- 0
unique_Combined_CDR3_day4["count_to"] <- 0
unique_Combined_CDR3_dayn14["count_from"] <- 0
unique_Combined_CDR3_dayn14["count_to"] <- 0

# Rearranging CDR3 sequences into alphabetical order

# TCRD:
shared_TCRD_CDR3 <- shared_TCRD_CDR3[order(shared_TCRD_CDR3$CDR3),]
unique_TCRD_CDR3_day4 <- unique_TCRD_CDR3_day4[order(unique_TCRD_CDR3_day4$CDR3),]
unique_TCRD_CDR3_dayn14 <- unique_TCRD_CDR3_dayn14[order(unique_TCRD_CDR3_dayn14$CDR3),]

# TCRG:
shared_TCRG_CDR3 <- shared_TCRG_CDR3[order(shared_TCRG_CDR3$CDR3),]
unique_TCRG_CDR3_day4 <- unique_TCRG_CDR3_day4[order(unique_TCRG_CDR3_day4$CDR3),]
unique_TCRG_CDR3_dayn14 <- unique_TCRG_CDR3_dayn14[order(unique_TCRG_CDR3_dayn14$CDR3),]

# Combined CDR3s:
shared_Combined_CDR3 <- shared_Combined_CDR3[order(shared_Combined_CDR3$CDR3),]
unique_Combined_CDR3_day4 <- unique_Combined_CDR3_day4[order(unique_Combined_CDR3_day4$CDR3),]
unique_Combined_CDR3_dayn14 <- unique_Combined_CDR3_dayn14[order(unique_Combined_CDR3_dayn14$CDR3),]

# Using Sum to calculate occurrences of CDR3

# TCRD:
for (i in 1:nrow(shared_TCRD_CDR3)) {
  shared_TCRD_CDR3$count_from[i] = sum(macTCRgd_AnimalID_n14$TCRD.CDR3.aa == shared_TCRD_CDR3$CDR3[i])
  shared_TCRD_CDR3$count_to[i] = sum(macTCRgd_AnimalID_4$TCRD.CDR3.aa == shared_TCRD_CDR3$CDR3[i])
}
for (i in 1:nrow(unique_TCRD_CDR3_day4)) {
  unique_TCRD_CDR3_day4$count_from[i] = sum(macTCRgd_AnimalID_4$TCRD.CDR3.aa == unique_TCRD_CDR3_day4$CDR3[i])
  unique_TCRD_CDR3_day4$count_to[i] = unique_TCRD_CDR3_day4$count_from[i]
}
for (i in 1:nrow(unique_TCRD_CDR3_dayn14)) {
  unique_TCRD_CDR3_dayn14$count_from[i] = sum(macTCRgd_AnimalID_n14$TCRD.CDR3.aa == unique_TCRD_CDR3_dayn14$CDR3[i])
  unique_TCRD_CDR3_dayn14$count_to[i] = unique_TCRD_CDR3_dayn14$count_from[i]
}

# TCRG:
for (i in 1:nrow(shared_TCRG_CDR3)) {
  shared_TCRG_CDR3$count_from[i] = sum(macTCRgd_AnimalID_n14$TCRG.CDR3.aa == shared_TCRG_CDR3$CDR3[i])
  shared_TCRG_CDR3$count_to[i] = sum(macTCRgd_AnimalID_4$TCRG.CDR3.aa == shared_TCRG_CDR3$CDR3[i])
}
for (i in 1:nrow(unique_TCRG_CDR3_day4)) {
  unique_TCRG_CDR3_day4$count_from[i] = sum(macTCRgd_AnimalID_4$TCRG.CDR3.aa == unique_TCRG_CDR3_day4$CDR3[i])
  unique_TCRG_CDR3_day4$count_to[i] = unique_TCRG_CDR3_day4$count_from[i]
}
for (i in 1:nrow(unique_TCRG_CDR3_dayn14)) {
  unique_TCRG_CDR3_dayn14$count_from[i] = sum(macTCRgd_AnimalID_n14$TCRG.CDR3.aa == unique_TCRG_CDR3_dayn14$CDR3[i])
  unique_TCRG_CDR3_dayn14$count_to[i] = unique_TCRG_CDR3_dayn14$count_from[i]
}

# Combined CDR3s:
for (i in 1:nrow(shared_Combined_CDR3)) {
  # If you add day 15 these should change, this assumes count_from in shared is always day -14
  shared_Combined_CDR3$count_from[i] = sum(macTCRgd_AnimalID_n14$Combined.CDR3 == shared_Combined_CDR3$CDR3[i])
  shared_Combined_CDR3$count_to[i] = sum(macTCRgd_AnimalID_4$Combined.CDR3 == shared_Combined_CDR3$CDR3[i])
}
for (i in 1:nrow(unique_Combined_CDR3_day4)) {
  unique_Combined_CDR3_day4$count_from[i] = sum(macTCRgd_AnimalID_4$Combined.CDR3 == unique_Combined_CDR3_day4$CDR3[i])
  unique_Combined_CDR3_day4$count_to[i] = unique_Combined_CDR3_day4$count_from[i]
}
for (i in 1:nrow(unique_Combined_CDR3_dayn14)) {
  unique_Combined_CDR3_dayn14$count_from[i] = sum(macTCRgd_AnimalID_n14$Combined.CDR3 == unique_Combined_CDR3_dayn14$CDR3[i])
  unique_Combined_CDR3_dayn14$count_to[i] = unique_Combined_CDR3_dayn14$count_from[i]
}

# Bind the shared and unique CDR3 sequences together and set the order

# TCRD
TCRD_CDR3_bind <- rbind(shared_TCRD_CDR3, unique_TCRD_CDR3_day4, unique_TCRD_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to,count_from, count_to)

# TCRG
TCRG_CDR3_bind <- rbind(shared_TCRG_CDR3, unique_TCRG_CDR3_day4, unique_TCRG_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to,count_from, count_to)

# Combined CDR3s :
Combined_CDR3_bind <- rbind(shared_Combined_CDR3, unique_Combined_CDR3_day4, unique_Combined_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to,count_from, count_to)


# Coerces bound data to be a data.table 
# If this doesn't work, run:  run: library(data.table) 
setDT(TCRD_CDR3_bind) 
setDT(TCRG_CDR3_bind) 
setDT(Combined_CDR3_bind) 

# Binds the CDR3 sequences, to the associated group_from_f and group_to_f data
# Then removed the CDR3 column

# TCRD
TCRD_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
TCRD_CDR3_bind[, CDR3 := NULL]   

# TCRG
TCRG_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
TCRG_CDR3_bind[, CDR3 := NULL]   

# Combined CDR3s :
Combined_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
Combined_CDR3_bind[, CDR3 := NULL]   

## Colour ##
# Adds color column to the bound CDR3 data.table

# TCRD
TCRD_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs

# TCRG
TCRG_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs

# Combined CDR3s :
Combined_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs


## Preparing Sectors around Circle ##

# TCRD
TCRD_CDR3_sectors_f <- union(TCRD_CDR3_bind[, group_from_f], TCRD_CDR3_bind[, group_to_f]) %>%
  sort()

# TCRG
TCRG_CDR3_sectors_f <- union(TCRG_CDR3_bind[, group_from_f], TCRG_CDR3_bind[, group_to_f]) %>%
  sort()

# Combined CDR3s :
Combined_CDR3_sectors_f <- union(Combined_CDR3_bind[, group_from_f], Combined_CDR3_bind[, group_to_f]) %>%
  sort()


# TCRD
TCRD_CDR3_colour_lookup <-
  union(TCRD_CDR3_bind[, group_from], TCRD_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)
TCRD_CDR3_sector_colours <- str_replace_all(TCRD_CDR3_sectors_f, "\\..+", "") %>%
  TCRD_CDR3_colour_lookup[.]
print(TCRD_CDR3_sector_colours)

# TCRG
TCRG_CDR3_colour_lookup <-
  union(TCRG_CDR3_bind[, group_from], TCRG_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)

TCRG_CDR3_sector_colours <- str_replace_all(TCRG_CDR3_sectors_f, "\\..+", "") %>%
  TCRG_CDR3_colour_lookup[.]

# Combined CDR3s 
Combined_CDR3_colour_lookup <-
  union(Combined_CDR3_bind[, group_from], Combined_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)
Combined_CDR3_sector_colours <- str_replace_all(Combined_CDR3_sectors_f, "\\..+", "") %>%
  Combined_CDR3_colour_lookup[.]

## Set Gap Info ## 
# CDR3 gap size = Universal accross CDR3s
# Increase the first number for small gaps between cells (within each day)
# Increase the second number for big gaps between sectors (aka between days)
CDR3_gap_sizes <- c(0.2, 1.0)

# TCRD
TCRD_CDR3_gap_degree <-
  sapply(table(names(TCRD_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

# TCRG
TCRG_CDR3_gap_degree <-
  sapply(table(names(TCRG_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

# Combined CDR3s :
Combined_CDR3_gap_degree <-
  sapply(table(names(Combined_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

## Setting the CDR3 xlims ## 
# Xlims are width of each cell

# TCRD
TCRD_CDR3_xlims <- matrix(0:0, nrow=length(TCRD_CDR3_sectors_f), ncol = 2)
for (i in 1:length((TCRD_CDR3_sectors_f))) {
  for (j in 1:nrow(TCRD_CDR3_bind)) {
    if (TCRD_CDR3_bind$group_from_f[j] == TCRD_CDR3_sectors_f[i]) {
      TCRD_CDR3_xlims[i,2] = TCRD_CDR3_bind$count_from[j]
    }
    if (TCRD_CDR3_bind$group_to_f[j] == TCRD_CDR3_sectors_f[i]) {
      TCRD_CDR3_xlims[i,2] = TCRD_CDR3_bind$count_to[j]
    }
  }
}

# TCRG
TCRG_CDR3_xlims <- matrix(0:0, nrow=length(TCRG_CDR3_sectors_f), ncol = 2)
for (i in 1:length((TCRG_CDR3_sectors_f))) {
  for (j in 1:nrow(TCRG_CDR3_bind)) {
    if (TCRG_CDR3_bind$group_from_f[j] == TCRG_CDR3_sectors_f[i]) {
      TCRG_CDR3_xlims[i,2] = TCRG_CDR3_bind$count_from[j]
    }
    if (TCRG_CDR3_bind$group_to_f[j] == TCRG_CDR3_sectors_f[i]) {
      TCRG_CDR3_xlims[i,2] = TCRG_CDR3_bind$count_to[j]
    }
  }
}

# Combined CDR3s 
Combined_CDR3_xlims <- matrix(0:0, nrow=length(Combined_CDR3_sectors_f), ncol = 2)
for (i in 1:length((Combined_CDR3_sectors_f))) {
  for (j in 1:nrow(Combined_CDR3_bind)) {
    if (Combined_CDR3_bind$group_from_f[j] == Combined_CDR3_sectors_f[i]) {
      Combined_CDR3_xlims[i,2] = Combined_CDR3_bind$count_from[j]
    }
    if (Combined_CDR3_bind$group_to_f[j] == Combined_CDR3_sectors_f[i]) {
      Combined_CDR3_xlims[i,2] = Combined_CDR3_bind$count_to[j]
    }
  }
}


## Plot! ##

## TCRD ##

tiff(filename = paste(AnimalID,"TCRD", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = TCRD_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(TCRD_CDR3_sectors_f, xlim = TCRD_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = TCRD_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(TCRD_CDR3_bind)) {
  row_i <- TCRD_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

##TCRG##

tiff(filename = paste(AnimalID,"TCRG", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = TCRG_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(TCRG_CDR3_sectors_f, xlim = TCRG_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = TCRG_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(TCRG_CDR3_bind)) {
  row_i <- TCRG_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

## Combined CDR3s ##

tiff(filename = paste(AnimalID,"Combined", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = Combined_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(Combined_CDR3_sectors_f, xlim = Combined_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = Combined_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(Combined_CDR3_bind)) {
  row_i <- Combined_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

####################################################
# Normalized
####################################################

#Normalized circos plots of CDR3 usage (all the corresponding data will have "norm_" in front of the name)

# Get the number of chains recovered at each time point
sum_macTCRgd_AnimalID <- c(nrow(macTCRgd_AnimalID_n14),
                           nrow(macTCRgd_AnimalID_4))

# Bind the shared and unique CDR3 sequences together and set the order

# TCRD
norm_TCRD_CDR3_bind <- rbind(shared_TCRD_CDR3, unique_TCRD_CDR3_day4, unique_TCRD_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to, count_from, count_to)

# TCRG
norm_TCRG_CDR3_bind <- rbind(shared_TCRG_CDR3, unique_TCRG_CDR3_day4, unique_TCRG_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to, count_from, count_to)

# Combined CDR3s 
norm_Combined_CDR3_bind <- rbind(shared_Combined_CDR3, unique_Combined_CDR3_day4, unique_Combined_CDR3_dayn14) %>%
  setorder(group_from, CDR3, group_to, count_from, count_to)

# Normalize the chain counts to the sums of the chain recoveries at each day:

# TCRD
for (i in 1:nrow(norm_TCRD_CDR3_bind)) {
  if (norm_TCRD_CDR3_bind$group_from[i] == "Day_n14") {
    norm_TCRD_CDR3_bind$count_from[i] <- (norm_TCRD_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_TCRD_CDR3_bind$group_from[i] == "Day_4") {
    norm_TCRD_CDR3_bind$count_from[i] <- (norm_TCRD_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[2]
  }
  if (norm_TCRD_CDR3_bind$group_to[i] == "Day_n14") {
    norm_TCRD_CDR3_bind$count_to[i] <- (norm_TCRD_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_TCRD_CDR3_bind$group_to[i] == "Day_4") {
    norm_TCRD_CDR3_bind$count_to[i] <- (norm_TCRD_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[2]
  }
}

# TCRG Normalize Count:
for (i in 1:nrow(norm_TCRG_CDR3_bind)) {
  if (norm_TCRG_CDR3_bind$group_from[i] == "Day_n14") {
    norm_TCRG_CDR3_bind$count_from[i] <- (norm_TCRG_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_TCRG_CDR3_bind$group_from[i] == "Day_4") {
    norm_TCRG_CDR3_bind$count_from[i] <- (norm_TCRG_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[2]
  }
  if (norm_TCRG_CDR3_bind$group_to[i] == "Day_n14") {
    norm_TCRG_CDR3_bind$count_to[i] <- (norm_TCRG_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_TCRG_CDR3_bind$group_to[i] == "Day_4") {
    norm_TCRG_CDR3_bind$count_to[i] <- (norm_TCRG_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[2]
  }
}


# TCR Combined Normalize Count:
for (i in 1:nrow(norm_Combined_CDR3_bind)) {
  if (norm_Combined_CDR3_bind$group_from[i] == "Day_n14") {
    norm_Combined_CDR3_bind$count_from[i] <- (norm_Combined_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_Combined_CDR3_bind$group_from[i] == "Day_4") {
    norm_Combined_CDR3_bind$count_from[i] <- (norm_Combined_CDR3_bind$count_from[i]) / sum_macTCRgd_AnimalID[2]
  }
  if (norm_Combined_CDR3_bind$group_to[i] == "Day_n14") {
    norm_Combined_CDR3_bind$count_to[i] <- (norm_Combined_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[1]
  }
  if (norm_Combined_CDR3_bind$group_to[i] == "Day_4") {
    norm_Combined_CDR3_bind$count_to[i] <- (norm_Combined_CDR3_bind$count_to[i]) / sum_macTCRgd_AnimalID[2]
  }
}


# Coerces bound data to be a data.table 
# If this doesn't work, run:  run: library(data.table) 
setDT(norm_TCRD_CDR3_bind) #TCRD
setDT(norm_TCRG_CDR3_bind) #TCRG
setDT(norm_Combined_CDR3_bind) #Combined

# Binds the CDR3 sequences, to the associated group_from_f and group_to_f data
# Then removes the CDR3 column

# TCRD
norm_TCRD_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
norm_TCRD_CDR3_bind[, CDR3 := NULL]

# TCRG
norm_TCRG_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
norm_TCRG_CDR3_bind[, CDR3 := NULL]   

# Combined CDR3s 
norm_Combined_CDR3_bind[, ":="(
  group_from_f = paste(group_from, CDR3, sep = "."),
  group_to_f = paste(group_to, CDR3, sep = "."))]
norm_Combined_CDR3_bind[, CDR3 := NULL]   

## Colour ##
# Adds color column to the bound CDR3 data.table

# TCRD
norm_TCRD_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs

# TCRG
norm_TCRG_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs

# Combined CDR3s 
norm_Combined_CDR3_bind[, colour := ifelse(group_from_f == group_to_f, "#FFFFFF00", "#00000050")]  # Change first to #FF0000FF to show red blobs


# Preparing Sectors around Circle

# TCRD
norm_TCRD_CDR3_sectors_f <- union(norm_TCRD_CDR3_bind[, group_from_f], norm_TCRD_CDR3_bind[, group_to_f]) %>%
  sort()

# TCRG
norm_TCRG_CDR3_sectors_f <- union(norm_TCRG_CDR3_bind[, group_from_f], norm_TCRG_CDR3_bind[, group_to_f]) %>%
  sort()

# Combined CDR3s 
norm_Combined_CDR3_sectors_f <- union(norm_Combined_CDR3_bind[, group_from_f], norm_Combined_CDR3_bind[, group_to_f]) %>%
  sort()

## Color_lookup ## 

# TCRD
norm_TCRD_CDR3_colour_lookup <-
  union(norm_TCRD_CDR3_bind[, group_from], norm_TCRD_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)
norm_TCRD_CDR3_sector_colours <- str_replace_all(norm_TCRD_CDR3_sectors_f, "\\..+", "") %>%
  norm_TCRD_CDR3_colour_lookup[.]

# TCRG
norm_TCRG_CDR3_colour_lookup <-
  union(norm_TCRG_CDR3_bind[, group_from], norm_TCRG_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)
norm_TCRG_CDR3_sector_colours <- str_replace_all(norm_TCRG_CDR3_sectors_f, "\\..+", "") %>%
  norm_TCRG_CDR3_colour_lookup[.]

# Combined CDR3s 
norm_Combined_CDR3_colour_lookup <-
  union(norm_Combined_CDR3_bind[, group_from], norm_Combined_CDR3_bind[, group_to]) %>% sort() %>%
  structure(seq_along(.) + 1, names = .)
norm_Combined_CDR3_sector_colours <- str_replace_all(norm_Combined_CDR3_sectors_f, "\\..+", "") %>%
  norm_Combined_CDR3_colour_lookup[.]

## CDR3 Gap Info ##

# CDR3 gap size = Universal accross CDR3s
# Increase the first number for small gaps between cells (within each day)
# Increase the second number for big gaps between sectors (aka between days)
CDR3_gap_sizes <- c(0.2, 1.0)

# TCRD
norm_TCRD_CDR3_gap_degree <-
  sapply(table(names(norm_TCRD_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

# TCRG
norm_TCRG_CDR3_gap_degree <-
  sapply(table(names(norm_TCRG_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

# Combined CDR3s 
norm_Combined_CDR3_gap_degree <-
  sapply(table(names(norm_Combined_CDR3_sector_colours)), function(i) c(rep(CDR3_gap_sizes[1], i-1), CDR3_gap_sizes[2])) %>%
  unlist() %>% unname() %>% as.vector()

## Setting the CDR3 xlims ## 
# Xlims are width of each cell

# TCRD
norm_TCRD_CDR3_xlims <- matrix(0:0, nrow=length(norm_TCRD_CDR3_sectors_f), ncol = 2)

for (i in 1:length((norm_TCRD_CDR3_sectors_f))) {
  for (j in 1:nrow(norm_TCRD_CDR3_bind)) {
    if (norm_TCRD_CDR3_bind$group_from_f[j] == norm_TCRD_CDR3_sectors_f[i]) {
      norm_TCRD_CDR3_xlims[i,2] = norm_TCRD_CDR3_bind$count_from[j]
    }
    if (norm_TCRD_CDR3_bind$group_to_f[j] == norm_TCRD_CDR3_sectors_f[i]) {
      norm_TCRD_CDR3_xlims[i,2] = norm_TCRD_CDR3_bind$count_to[j]
    }
  }
}

# TCRG
norm_TCRG_CDR3_xlims <- matrix(0:0, nrow=length(norm_TCRG_CDR3_sectors_f), ncol = 2)
for (i in 1:length((norm_TCRG_CDR3_sectors_f))) {
  for (j in 1:nrow(norm_TCRG_CDR3_bind)) {
    if (norm_TCRG_CDR3_bind$group_from_f[j] == norm_TCRG_CDR3_sectors_f[i]) {
      norm_TCRG_CDR3_xlims[i,2] = norm_TCRG_CDR3_bind$count_from[j]
    }
    if (norm_TCRG_CDR3_bind$group_to_f[j] == norm_TCRG_CDR3_sectors_f[i]) {
      norm_TCRG_CDR3_xlims[i,2] = norm_TCRG_CDR3_bind$count_to[j]
    }
  }
}

# Combined CDR3s 
norm_Combined_CDR3_xlims <- matrix(0:0, nrow=length(norm_Combined_CDR3_sectors_f), ncol = 2)
for (i in 1:length((norm_Combined_CDR3_sectors_f))) {
  for (j in 1:nrow(norm_Combined_CDR3_bind)) {
    if (norm_Combined_CDR3_bind$group_from_f[j] == norm_Combined_CDR3_sectors_f[i]) {
      norm_Combined_CDR3_xlims[i,2] = norm_Combined_CDR3_bind$count_from[j]
    }
    if (norm_Combined_CDR3_bind$group_to_f[j] == norm_Combined_CDR3_sectors_f[i]) {
      norm_Combined_CDR3_xlims[i,2] = norm_Combined_CDR3_bind$count_to[j]
    }
  }
}


## Plot! ##

## TCRD ##

tiff(filename = paste(AnimalID,"TCRD_Normalized", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = norm_TCRD_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(norm_TCRD_CDR3_sectors_f, xlim = norm_TCRD_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = norm_TCRD_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(norm_TCRD_CDR3_bind)) {
  row_i <- norm_TCRD_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

## TCRG ##

tiff(filename = paste(AnimalID,"TCRG_Norm", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = norm_TCRG_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(norm_TCRG_CDR3_sectors_f, xlim = norm_TCRG_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = norm_TCRG_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(norm_TCRG_CDR3_bind)) {
  row_i <- norm_TCRG_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

## Combined CDR3s ##

tiff(filename = paste(AnimalID,"Combined_Norm", "Circos", "CDR3.tiff", sep = "_"),
     height = 18,
     width = 18,
     units = "cm",
     res = 1200)

# Each "sector" is a separate patient/cell/feature combination
circos.clear()
circos.par(gap.degree = norm_Combined_CDR3_gap_degree,
           cell.padding = c(0.02, 0.1, 0.02, 0.1))
circos.initialize(norm_Combined_CDR3_sectors_f, xlim = norm_Combined_CDR3_xlims)
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = norm_Combined_CDR3_sector_colours, bg.border = NA)

for (i in 1:nrow(norm_Combined_CDR3_bind)) {
  row_i <- norm_Combined_CDR3_bind[i, ]
  circos.link(
    row_i[["group_from_f"]], c(0, row_i[["count_from"]]),
    row_i[["group_to_f"]], c(0, row_i[["count_to"]]),
    border = NA, col = row_i[["colour"]]
  )
}

dev.off()

setwd("..")