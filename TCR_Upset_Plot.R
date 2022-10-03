## TCR_Upset_Plot.R
# This inputs the data formatted by "TCR_Data_Preperation.R"
# and formats it specifically to create 3 upset plots: TCRD, TCRG, TCRD-TCRG.
# Good video resource for UpSet plots: https://www.youtube.com/watch?v=n9MRCZxJOfk&ab_channel=mathetal

# Removes ALL objects from the global environment
rm(list = ls())

# Installs and Loads Packages
install.packages("tidyverse")
install.packages("data.table")
install.packages("UpSetR")
library(tidyverse)
library(data.table)
library(UpSetR)
library(turner)


# Imports data formatted by TCR_Data_Preparation_MiXCR.R Script
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()


# Create Lists of CDR3 sequence (TCRD, TCRG, TCRD-TCRG) for each animal at each timepoint

# Create filtered data.frames by Animal.ID and Day:
 NM11_Dn14 <- macTCRgd %>%
   filter(Animal.ID == "NM11") %>%
   filter(Day == -14)

 NM11_D4 <- macTCRgd %>%
   filter(Animal.ID == "NM11") %>%
   filter(Day == 4)

NM89_Dn14 <- macTCRgd %>%
  filter(Animal.ID == "NM89") %>%
  filter(Day == -14)

NM89_D4 <- macTCRgd %>%
  filter(Animal.ID == "NM89") %>%
  filter(Day == 4)

NM89_D15 <- macTCRgd %>%
  filter(Animal.ID == "NM89") %>%
  filter(Day == 15)

NM251_Dn14 <- macTCRgd %>%
  filter(Animal.ID == "NM251") %>%
  filter(Day == -14)

NM251_D4 <- macTCRgd %>%
  filter(Animal.ID == "NM251") %>%
  filter(Day == 4)

NM295_Dn14 <- macTCRgd %>%
  filter(Animal.ID == "NM295") %>%
  filter(Day == -14)

NM295_D4 <- macTCRgd %>%
  filter(Animal.ID == "NM295") %>%
  filter(Day == 4)


# Create vectors (values) for TCRD, TCRG, and TCRD-TCRG from the filtered data.frames, and concatentate for each day

# Day -14; TCRD:
NM11_Dn14_TCRD <- as.vector(NM11_Dn14$TCRD.CDR3.aa)
NM89_Dn14_TCRD <- as.vector(NM89_Dn14$TCRD.CDR3.aa)
NM251_Dn14_TCRD <- as.vector(NM251_Dn14$TCRD.CDR3.aa)
NM295_Dn14_TCRD <- as.vector(NM295_Dn14$TCRD.CDR3.aa)


# Day -14; TCRG:
NM11_Dn14_TCRG <- as.vector(NM11_Dn14$TCRG.CDR3.aa)
NM89_Dn14_TCRG <- as.vector(NM89_Dn14$TCRG.CDR3.aa)
NM251_Dn14_TCRG <- as.vector(NM251_Dn14$TCRG.CDR3.aa)
NM295_Dn14_TCRG <- as.vector(NM295_Dn14$TCRG.CDR3.aa)

# Day -14; Combined TCRD-TCRG:
NM11_Dn14_Combined <- as.vector(NM11_Dn14$Combined.CDR3)
NM89_Dn14_Combined <- as.vector(NM89_Dn14$Combined.CDR3)
NM251_Dn14_Combined <- as.vector(NM251_Dn14$Combined.CDR3)
NM295_Dn14_Combined <- as.vector(NM295_Dn14$Combined.CDR3)

# Day 4; TCRD:
NM11_D4_TCRD <- as.vector(NM11_D4$TCRD.CDR3.aa)
NM89_D4_TCRD <- as.vector(NM89_D4$TCRD.CDR3.aa)
NM251_D4_TCRD <- as.vector(NM251_D4$TCRD.CDR3.aa)
NM295_D4_TCRD <- as.vector(NM295_D4$TCRD.CDR3.aa)

# Day 4; TCRG:
NM11_D4_TCRG <- as.vector(NM11_D4$TCRG.CDR3.aa)
NM89_D4_TCRG <- as.vector(NM89_D4$TCRG.CDR3.aa)
NM251_D4_TCRG <- as.vector(NM251_D4$TCRG.CDR3.aa)
NM295_D4_TCRG <- as.vector(NM295_D4$TCRG.CDR3.aa)

# Day 4; Combined TCRD-TCRG:
NM11_D4_Combined <- as.vector(NM11_D4$Combined.CDR3)
NM89_D4_Combined <- as.vector(NM89_D4$Combined.CDR3)
NM251_D4_Combined <- as.vector(NM251_D4$Combined.CDR3)
NM295_D4_Combined <- as.vector(NM295_D4$Combined.CDR3)


# Concatenating chains across multiple timepoints, to look for common CDR3 usage

# TCRD, Day -14 and day 4:

All_TCRD <- list(
  NM11_Dn14 = NM11_Dn14_TCRD,
  NM11_D4 = NM11_D4_TCRD,
  NM89_Dn14 = NM89_Dn14_TCRD,
  NM89_D4 = NM89_D4_TCRD,
  NM251_Dn14 = NM251_Dn14_TCRD,
  NM251_D4 = NM251_D4_TCRD,
  NM295_Dn14 = NM295_Dn14_TCRD,
  NM295_D4 = NM295_D4_TCRD
  )

# Rename list for UpSet plot labelling
names(All_TCRD) <- c("NM11 Day -14",
                    "NM11 Day 4",
                    "NM89 Day -14",
                    "NM89 Day 4",
                    "NM251 Day -14",
                    "NM251 Day 4",
                    "NM295 Day -14",
                    "NM295 Day 4"
                    )

# TCRG, Day -14 and day 4:
All_TCRG <- list(
  NM11_Dn14 = NM11_Dn14_TCRG,
  NM11_D4 = NM11_D4_TCRG,
  NM89_Dn14 = NM89_Dn14_TCRG,
  NM89_D4 = NM89_D4_TCRG,
  NM251_Dn14 = NM251_Dn14_TCRG,
  NM251_D4 = NM251_D4_TCRG,
  NM295_Dn14 = NM295_Dn14_TCRG,
  NM295_D4 = NM295_D4_TCRG
)

# Rename list for UpSet plot labelling
names(All_TCRG) <- c("NM11 Day -14",
                    "NM11 Day 4",
                    "NM89 Day -14",
                    "NM89 Day 4",
                    "NM251 Day -14",
                    "NM251 Day 4",
                    "NM295 Day -14",
                    "NM295 Day 4"
)

# Combined, day -14 and day 4:
All_Combined <- list(
  NM11_Dn14 = NM11_Dn14_Combined,
  NM11_D4 = NM11_D4_Combined,
  NM89_Dn14 = NM89_Dn14_Combined,
  NM89_D4 = NM89_D4_Combined,
  NM251_Dn14 = NM251_Dn14_Combined,
  NM251_D4 = NM251_D4_Combined,
  NM295_Dn14 = NM295_Dn14_Combined,
  NM295_D4 = NM295_D4_Combined
)

# Rename list for UpSet plot labelling
names(All_Combined) <- c("NM11 Day -14",
                    "NM11 Day 4",
                    "NM89 Day -14",
                    "NM89 Day 4",
                    "NM251 Day -14",
                    "NM251 Day 4",
                    "NM295 Day -14",
                    "NM295 Day 4"
)


## Upset Plots ##

setwd("Outputs")

# TCRD:

tiff(filename = paste(paste("TCRD", "UpSet.tiff", sep = "_")),
     height = 15,
     width = 50,
     units = "cm",
     res = 1200
)

upset(fromList(All_TCRD),
      sets = c("NM295 Day 4",
               "NM295 Day -14",
               "NM251 Day 4",
               "NM251 Day -14",
               "NM89 Day 4",
               "NM89 Day -14",
               "NM11 Day 4",
               "NM11 Day -14"
               ),
      nintersects = 100,
      mainbar.y.label = "TCRD Intersections",
      sets.x.label = "Unique TCRD sequences",
      number.angles = 0,
      order.by = "freq",
      #empty.intersections = "on",
      set_size.show = TRUE,
      text.scale = c(1.5, 1.2, 1.5, 1, 1.2, 1.2),
      keep.order = TRUE
      )

dev.off()

# TCRG:

tiff(filename = paste(paste("TCRG", "UpSet.tiff", sep = "_")),
     height = 15,
     width = 50,
     units = "cm",
     res = 1200
)

upset(fromList(All_TCRG),
      sets = c("NM295 Day 4",
               "NM295 Day -14",
               "NM251 Day 4",
               "NM251 Day -14",
               "NM89 Day 4",
               "NM89 Day -14",
               "NM11 Day 4",
               "NM11 Day -14"
      ),
      nintersects = 100,
      mainbar.y.label = "TCRG CDR3 Intersections",
      sets.x.label = "Unique TCRG CDR3 sequences",
      number.angles = 0,
      order.by = "freq",
      #empty.intersections = "on",
      set_size.show = TRUE,
      text.scale = c(1.5, 1.2, 1.5, 1, 1.2, 1.2),
      keep.order = TRUE
)

dev.off()

# Combined TCRD/TCRG:

tiff(filename = paste(paste("Combined", "UpSet.tiff", sep = "_")),
     height = 15,
     width = 50,
     units = "cm",
     res = 1200
     )

upset(fromList(All_Combined),
      sets = c("NM295 Day 4",
               "NM295 Day -14",
               "NM251 Day 4",
               "NM251 Day -14",
               "NM89 Day 4",
               "NM89 Day -14",
               "NM11 Day 4",
               "NM11 Day -14"
      ),
      nintersects = 100,
      mainbar.y.label = "TCRD/TCRG CDR3 Intersections",
      sets.x.label = "Unique TCRD/TCRG CDR3 sequences",
      number.angles = 0,
      order.by = "freq",
      #empty.intersections = "on",
      set_size.show = TRUE,
      text.scale = c(1.5, 1.2, 1.5, 1, 1.2, 1.2),
      keep.order = TRUE
)

dev.off()

setwd("..")