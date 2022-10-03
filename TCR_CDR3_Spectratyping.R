## TCR_CDR3_Spectratyping.R
# This script outputs a spectratyping bar plot showing the length distribution 
# of the TCRD/TCRG CDR3s (aa) and saves statistics for CDR3 Lengths for 2 of 3 
# timepoints in a CSV.

# Removes ALL objects from the global environment
rm(list = ls())

# Installs and Loads Packages
install.packages("tidyverse")
# DescTools package
install.packages("e1071")
library(tidyverse)
library(e1071)

## Loading MiXCR Data ##
macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()

# Color blind safe palette with black- for reference w/ ggplot2: https://personal.sron.nl/~pault/
cbbPalette <- c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499", "#000000")

setwd("Outputs")

## Start of Analysis ##

# Filter by Animal.ID
macTCRgd_NM11 <- filter(macTCRgd, Animal.ID == "NM11")
macTCRgd_NM89 <- filter(macTCRgd, Animal.ID == "NM89")
macTCRgd_NM251 <- filter(macTCRgd, Animal.ID == "NM251")
macTCRgd_NM295 <- filter(macTCRgd, Animal.ID == "NM295")

# Filter the Animal ID by Day
macTCRgd_NM11_D4 <- filter(macTCRgd_NM11, Day == 4)
macTCRgd_NM11_Dn14 <- filter(macTCRgd_NM11, Day == -14)
#macTCRgd_NM11_D15 <- filter(macTCRgd_NM11, Day == 15)

macTCRgd_NM89_D4 <- filter(macTCRgd_NM89, Day == 4)
macTCRgd_NM89_Dn14 <- filter(macTCRgd_NM89, Day == -14)
#macTCRgd_NM89_D15 <- filter(macTCRgd_NM89, Day == 15)

macTCRgd_NM251_D4 <- filter(macTCRgd_NM251, Day == 4)
macTCRgd_NM251_Dn14 <- filter(macTCRgd_NM251, Day == -14)
#macTCRgd_NM251_D15 <- filter(macTCRgd_NM251, Day == 15)

macTCRgd_NM295_D4 <- filter(macTCRgd_NM295, Day == 4)
macTCRgd_NM295_Dn14 <- filter(macTCRgd_NM295, Day == -14)
#macTCRgd_NM295_D15 <- filter(macTCRgd_NM295, Day == 15)

## Analysis Overview:
# For each animal, timepoint, and TCRD/TCRG:
# Create a dataframe of unique CDR3 lengths (aa)
# For each animal, TCRD and TCRG data is plotted seperately.
# At the end, all data for all animals/timepoints/TCR chains is combined for statistics.

## NM11 ##

# NM11, Day -14, TCRD:
unique_CDR3length_NM11_Dn14_TCRD <- as.numeric(unique(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM11_Dn14_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM11_Dn14_TCRD)))
colnames(CDR3length_NM11_Dn14_TCRD) <- c("count")
rownames(CDR3length_NM11_Dn14_TCRD) <- unique_CDR3length_NM11_Dn14_TCRD

for (i in 1:nrow(macTCRgd_NM11_Dn14)) {
  for (j in rownames(CDR3length_NM11_Dn14_TCRD)) {
    if (macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length[i] == j) {
      CDR3length_NM11_Dn14_TCRD[j, "count"] = CDR3length_NM11_Dn14_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM11_Dn14_TCRD <- CDR3length_NM11_Dn14_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM11_Dn14_TCRD) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM11")

# NM11, Day 4, TCRD:
unique_CDR3length_NM11_D4_TCRD <- as.numeric(unique(macTCRgd_NM11_D4$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM11_D4_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM11_D4_TCRD)))
colnames(CDR3length_NM11_D4_TCRD) <- c("count")
rownames(CDR3length_NM11_D4_TCRD) <- unique_CDR3length_NM11_D4_TCRD

for (i in 1:nrow(macTCRgd_NM11_D4)){
  for (j in rownames(CDR3length_NM11_D4_TCRD)){
    if (macTCRgd_NM11_D4$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM11_D4_TCRD[j, "count"] = CDR3length_NM11_D4_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM11_D4_TCRD <- CDR3length_NM11_D4_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM11_D4_TCRD) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM11")

# NM11 Combine Day -14 and Day 4 data.frames
CDR3length_NM11_TCRD <- CDR3length_NM11_Dn14_TCRD %>%
  bind_rows(CDR3length_NM11_D4_TCRD)

# Plot Data and Save Plot

plot_CDR3length_NM11_TCRD <- ggplot(CDR3length_NM11_TCRD, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM11_TCRD

ggsave(paste("NM11", "TCRD", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# NM11, Day -14, TCRG:
unique_CDR3length_NM11_Dn14_TCRG <- as.numeric(unique(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM11_Dn14_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM11_Dn14_TCRG)))
colnames(CDR3length_NM11_Dn14_TCRG) <- c("count")
rownames(CDR3length_NM11_Dn14_TCRG) <- unique_CDR3length_NM11_Dn14_TCRG

for (i in 1:nrow(macTCRgd_NM11_Dn14)){
  for (j in rownames(CDR3length_NM11_Dn14_TCRG)){
    if (macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM11_Dn14_TCRG[j, "count"] = CDR3length_NM11_Dn14_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM11_Dn14_TCRG <- CDR3length_NM11_Dn14_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM11_Dn14_TCRG) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM11")

# NM11, Day 4, TCRG:
unique_CDR3length_NM11_D4_TCRG <- as.numeric(unique(macTCRgd_NM11_D4$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM11_D4_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM11_D4_TCRG)))
colnames(CDR3length_NM11_D4_TCRG) <- c("count")
rownames(CDR3length_NM11_D4_TCRG) <- unique_CDR3length_NM11_D4_TCRG

for (i in 1:nrow(macTCRgd_NM11_D4)){
  for (j in rownames(CDR3length_NM11_D4_TCRG)){
    if (macTCRgd_NM11_D4$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM11_D4_TCRG[j, "count"] = CDR3length_NM11_D4_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM11_D4_TCRG <- CDR3length_NM11_D4_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM11_D4_TCRG) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM11")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM11_TCRG <- CDR3length_NM11_Dn14_TCRG %>%
  bind_rows(CDR3length_NM11_D4_TCRG)

# Plot and Save Data
plot_CDR3length_NM11_TCRG <- ggplot(CDR3length_NM11_TCRG, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM11_TCRG

ggsave(paste("NM11", "TCRG", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

## NM89 ##

# NM89, Day -14, TCRD:
unique_CDR3length_NM89_Dn14_TCRD <- as.numeric(unique(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_Dn14_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_Dn14_TCRD)))
colnames(CDR3length_NM89_Dn14_TCRD) <- c("count")
rownames(CDR3length_NM89_Dn14_TCRD) <- unique_CDR3length_NM89_Dn14_TCRD

for (i in 1:nrow(macTCRgd_NM89_Dn14)){
  for (j in rownames(CDR3length_NM89_Dn14_TCRD)){
    if (macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM89_Dn14_TCRD[j, "count"] = CDR3length_NM89_Dn14_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_Dn14_TCRD <- CDR3length_NM89_Dn14_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_Dn14_TCRD) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM89")

# NM89, Day 4, TCRD:
unique_CDR3length_NM89_D4_TCRD <- as.numeric(unique(macTCRgd_NM89_D4$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_D4_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_D4_TCRD)))
colnames(CDR3length_NM89_D4_TCRD) <- c("count")
rownames(CDR3length_NM89_D4_TCRD) <- unique_CDR3length_NM89_D4_TCRD

for (i in 1:nrow(macTCRgd_NM89_D4)){
  for (j in rownames(CDR3length_NM89_D4_TCRD)){
    if (macTCRgd_NM89_D4$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM89_D4_TCRD[j, "count"] = CDR3length_NM89_D4_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_D4_TCRD <- CDR3length_NM89_D4_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_D4_TCRD) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM89")

# NM89, Day 15, TCRD:
unique_CDR3length_NM89_D15_TCRD <- as.numeric(unique(macTCRgd_NM89_D15$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_D15_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_D15_TCRD)))
colnames(CDR3length_NM89_D15_TCRD) <- c("count")
rownames(CDR3length_NM89_D15_TCRD) <- unique_CDR3length_NM89_D15_TCRD

for (i in 1:nrow(macTCRgd_NM89_D15)){
  for (j in rownames(CDR3length_NM89_D15_TCRD)){
    if (macTCRgd_NM89_D15$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM89_D15_TCRD[j, "count"] = CDR3length_NM89_D15_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_D15_TCRD <- CDR3length_NM89_D15_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_D15_TCRD) %>%
  add_column("Day" = "Day 15") %>%
  add_column("Animal.ID" = "NM89")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM89_TCRD <- CDR3length_NM89_Dn14_TCRD %>%
  bind_rows(CDR3length_NM89_D4_TCRD)

# Plot Data and Save
plot_CDR3length_NM89_TCRD <- ggplot(CDR3length_NM89_TCRD, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM89_TCRD

ggsave(paste("NM89", "TCRD", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# Combine Day -14 and Day 4 and day 15 data.frame (NM89 only) ###
CDR3length_NM89_TCRD_D15 <- CDR3length_NM89_Dn14_TCRD %>%
  bind_rows(CDR3length_NM89_D4_TCRD) %>% 
  bind_rows(CDR3length_NM89_D15_TCRD)

# Plot Data and Save
plot_CDR3length_NM89_TCRD_D15 <- ggplot(CDR3length_NM89_TCRD_D15, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM89_TCRD_D15

ggsave(paste("NM89", "TCRD", "CDR3_Spectratyping_D15.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# NM89, Day -14, TCRG:
unique_CDR3length_NM89_Dn14_TCRG <- as.numeric(unique(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_Dn14_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_Dn14_TCRG)))
colnames(CDR3length_NM89_Dn14_TCRG) <- c("count")
rownames(CDR3length_NM89_Dn14_TCRG) <- unique_CDR3length_NM89_Dn14_TCRG

for (i in 1:nrow(macTCRgd_NM89_Dn14)){
  for (j in rownames(CDR3length_NM89_Dn14_TCRG)){
    if (macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM89_Dn14_TCRG[j, "count"] = CDR3length_NM89_Dn14_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_Dn14_TCRG <- CDR3length_NM89_Dn14_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_Dn14_TCRG) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM89")

# NM89, Day 4, TCRG:
unique_CDR3length_NM89_D4_TCRG <- as.numeric(unique(macTCRgd_NM89_D4$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_D4_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_D4_TCRG)))
colnames(CDR3length_NM89_D4_TCRG) <- c("count")
rownames(CDR3length_NM89_D4_TCRG) <- unique_CDR3length_NM89_D4_TCRG

for (i in 1:nrow(macTCRgd_NM89_D4)){
  for (j in rownames(CDR3length_NM89_D4_TCRG)){
    if (macTCRgd_NM89_D4$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM89_D4_TCRG[j, "count"] = CDR3length_NM89_D4_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_D4_TCRG <- CDR3length_NM89_D4_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_D4_TCRG) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM89")

# NM89, Day 15, TCRG:
unique_CDR3length_NM89_D15_TCRG <- as.numeric(unique(macTCRgd_NM89_D15$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM89_D15_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM89_D15_TCRG)))
colnames(CDR3length_NM89_D15_TCRG) <- c("count")
rownames(CDR3length_NM89_D15_TCRG) <- unique_CDR3length_NM89_D15_TCRG

for (i in 1:nrow(macTCRgd_NM89_D15)){
  for (j in rownames(CDR3length_NM89_D15_TCRG)){
    if (macTCRgd_NM89_D15$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM89_D15_TCRG[j, "count"] = CDR3length_NM89_D15_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM89_D15_TCRG <- CDR3length_NM89_D15_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM89_D15_TCRG) %>%
  add_column("Day" = "Day 15") %>%
  add_column("Animal.ID" = "NM89")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM89_TCRG <- CDR3length_NM89_Dn14_TCRG %>%
  bind_rows(CDR3length_NM89_D4_TCRG)

# Plot Data and Save
plot_CDR3length_NM89_TCRG <- ggplot(CDR3length_NM89_TCRG, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM89_TCRG

ggsave(paste("NM89", "TCRG", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# Combine Day -14 and Day 4 (option for Day 15)
CDR3length_NM89_TCRG_D15 <- CDR3length_NM89_Dn14_TCRG %>%
  bind_rows(CDR3length_NM89_D4_TCRG) #%>%
  #bind_rows(CDR3length_NM89_D15_TCRG)

# Plot Data and Save
plot_CDR3length_NM89_TCRG_D15 <- ggplot(CDR3length_NM89_TCRG_D15, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM89_TCRG_D15

ggsave(paste("NM89", "TCRG", "CDR3_Spectratyping_D15.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

## NM251 ##

# NM251, Day -14, TCRD:
unique_CDR3length_NM251_Dn14_TCRD <- as.numeric(unique(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM251_Dn14_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM251_Dn14_TCRD)))
colnames(CDR3length_NM251_Dn14_TCRD) <- c("count")
rownames(CDR3length_NM251_Dn14_TCRD) <- unique_CDR3length_NM251_Dn14_TCRD

for (i in 1:nrow(macTCRgd_NM251_Dn14)){
  for (j in rownames(CDR3length_NM251_Dn14_TCRD)){
    if (macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM251_Dn14_TCRD[j, "count"] = CDR3length_NM251_Dn14_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM251_Dn14_TCRD <- CDR3length_NM251_Dn14_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM251_Dn14_TCRD) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM251")

# NM251, Day 4, TCRD:
unique_CDR3length_NM251_D4_TCRD <- as.numeric(unique(macTCRgd_NM251_D4$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM251_D4_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM251_D4_TCRD)))
colnames(CDR3length_NM251_D4_TCRD) <- c("count")
rownames(CDR3length_NM251_D4_TCRD) <- unique_CDR3length_NM251_D4_TCRD

for (i in 1:nrow(macTCRgd_NM251_D4)){
  for (j in rownames(CDR3length_NM251_D4_TCRD)){
    if (macTCRgd_NM251_D4$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM251_D4_TCRD[j, "count"] = CDR3length_NM251_D4_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM251_D4_TCRD <- CDR3length_NM251_D4_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM251_D4_TCRD) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM251")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM251_TCRD <- CDR3length_NM251_Dn14_TCRD %>%
  bind_rows(CDR3length_NM251_D4_TCRD)

# Plot Data and Save
plot_CDR3length_NM251_TCRD <- ggplot(CDR3length_NM251_TCRD, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM251_TCRD

ggsave(paste("NM251", "TCRD", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# NM251, Day -14, TCRG:
unique_CDR3length_NM251_Dn14_TCRG <- as.numeric(unique(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM251_Dn14_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM251_Dn14_TCRG)))
colnames(CDR3length_NM251_Dn14_TCRG) <- c("count")
rownames(CDR3length_NM251_Dn14_TCRG) <- unique_CDR3length_NM251_Dn14_TCRG

for (i in 1:nrow(macTCRgd_NM251_Dn14)){
  for (j in rownames(CDR3length_NM251_Dn14_TCRG)){
    if (macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM251_Dn14_TCRG[j, "count"] = CDR3length_NM251_Dn14_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM251_Dn14_TCRG <- CDR3length_NM251_Dn14_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM251_Dn14_TCRG) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM251")

# NM251, Day 4, TCRG:
unique_CDR3length_NM251_D4_TCRG <- as.numeric(unique(macTCRgd_NM251_D4$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM251_D4_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM251_D4_TCRG)))
colnames(CDR3length_NM251_D4_TCRG) <- c("count")
rownames(CDR3length_NM251_D4_TCRG) <- unique_CDR3length_NM251_D4_TCRG

for (i in 1:nrow(macTCRgd_NM251_D4)){
  for (j in rownames(CDR3length_NM251_D4_TCRG)){
    if (macTCRgd_NM251_D4$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM251_D4_TCRG[j, "count"] = CDR3length_NM251_D4_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM251_D4_TCRG <- CDR3length_NM251_D4_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM251_D4_TCRG) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM251")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM251_TCRG <- CDR3length_NM251_Dn14_TCRG %>%
  bind_rows(CDR3length_NM251_D4_TCRG)

# Plot Data and Save
plot_CDR3length_NM251_TCRG <- ggplot(CDR3length_NM251_TCRG, 
                                    aes(x = `CDR3 length aa`, 
                                        y = `count`, 
                                        fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM251_TCRG

ggsave(paste("NM251", "TCRG", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

## NM295 ##

# NM295, Day -14, TCRD:
unique_CDR3length_NM295_Dn14_TCRD <- as.numeric(unique(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM295_Dn14_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM295_Dn14_TCRD)))
colnames(CDR3length_NM295_Dn14_TCRD) <- c("count")
rownames(CDR3length_NM295_Dn14_TCRD) <- unique_CDR3length_NM295_Dn14_TCRD

for (i in 1:nrow(macTCRgd_NM295_Dn14)){
  for (j in rownames(CDR3length_NM295_Dn14_TCRD)){
    if (macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM295_Dn14_TCRD[j, "count"] = CDR3length_NM295_Dn14_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM295_Dn14_TCRD <- CDR3length_NM295_Dn14_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM295_Dn14_TCRD) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM295")

# NM295, Day 4, TCRD:
unique_CDR3length_NM295_D4_TCRD <- as.numeric(unique(macTCRgd_NM295_D4$TCRD.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM295_D4_TCRD <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM295_D4_TCRD)))
colnames(CDR3length_NM295_D4_TCRD) <- c("count")
rownames(CDR3length_NM295_D4_TCRD) <- unique_CDR3length_NM295_D4_TCRD

for (i in 1:nrow(macTCRgd_NM295_D4)){
  for (j in rownames(CDR3length_NM295_D4_TCRD)){
    if (macTCRgd_NM295_D4$TCRD.CDR3.aa.length[i] == j){
      CDR3length_NM295_D4_TCRD[j, "count"] = CDR3length_NM295_D4_TCRD[j, "count"] + 1
    } 
  }
}

CDR3length_NM295_D4_TCRD <- CDR3length_NM295_D4_TCRD %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM295_D4_TCRD) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM295")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM295_TCRD <- CDR3length_NM295_Dn14_TCRD %>%
  bind_rows(CDR3length_NM295_D4_TCRD)

# Plot Data and Save
plot_CDR3length_NM295_TCRD <- ggplot(CDR3length_NM295_TCRD, 
                                     aes(x = `CDR3 length aa`, 
                                         y = `count`, 
                                         fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM295_TCRD

ggsave(paste("NM295", "TCRD", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# NM295, Day -14, TCRG:
unique_CDR3length_NM295_Dn14_TCRG <- as.numeric(unique(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM295_Dn14_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM295_Dn14_TCRG)))
colnames(CDR3length_NM295_Dn14_TCRG) <- c("count")
rownames(CDR3length_NM295_Dn14_TCRG) <- unique_CDR3length_NM295_Dn14_TCRG

for (i in 1:nrow(macTCRgd_NM295_Dn14)){
  for (j in rownames(CDR3length_NM295_Dn14_TCRG)){
    if (macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM295_Dn14_TCRG[j, "count"] = CDR3length_NM295_Dn14_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM295_Dn14_TCRG <- CDR3length_NM295_Dn14_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM295_Dn14_TCRG) %>%
  add_column("Day" = "Day -14") %>%
  add_column("Animal.ID" = "NM295")

# NM295, Day 4, TCRG:
unique_CDR3length_NM295_D4_TCRG <- as.numeric(unique(macTCRgd_NM295_D4$TCRG.CDR3.aa.length)) %>%
  sort(decreasing = FALSE)

CDR3length_NM295_D4_TCRG <- as.data.frame(matrix(0, ncol = 1, nrow = length(unique_CDR3length_NM295_D4_TCRG)))
colnames(CDR3length_NM295_D4_TCRG) <- c("count")
rownames(CDR3length_NM295_D4_TCRG) <- unique_CDR3length_NM295_D4_TCRG

for (i in 1:nrow(macTCRgd_NM295_D4)){
  for (j in rownames(CDR3length_NM295_D4_TCRG)){
    if (macTCRgd_NM295_D4$TCRG.CDR3.aa.length[i] == j){
      CDR3length_NM295_D4_TCRG[j, "count"] = CDR3length_NM295_D4_TCRG[j, "count"] + 1
    } 
  }
}

CDR3length_NM295_D4_TCRG <- CDR3length_NM295_D4_TCRG %>%
  add_column("CDR3 length aa" = unique_CDR3length_NM295_D4_TCRG) %>%
  add_column("Day" = "Day 4") %>%
  add_column("Animal.ID" = "NM295")

# Combine Day -14 and Day 4 data.frames
CDR3length_NM295_TCRG <- CDR3length_NM295_Dn14_TCRG %>%
  bind_rows(CDR3length_NM295_D4_TCRG)

# Plot Data and Save
plot_CDR3length_NM295_TCRG <- ggplot(CDR3length_NM295_TCRG, 
                                     aes(x = `CDR3 length aa`, 
                                         y = `count`, 
                                         fill = `Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_NM295_TCRG

ggsave(paste("NM295", "TCRG", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)


## All Animals and Timepoints ##

# All, TCRD:
# Combine all TCRD data.frames together, and merge Animal.ID and Day in new column
CDR3length_All_TCRD <- CDR3length_NM11_Dn14_TCRD %>%
  bind_rows(CDR3length_NM11_D4_TCRD) %>%
  bind_rows(CDR3length_NM89_Dn14_TCRD) %>% 
  bind_rows(CDR3length_NM89_D4_TCRD) %>%
  bind_rows(CDR3length_NM251_Dn14_TCRD) %>%
  bind_rows(CDR3length_NM251_D4_TCRD) %>%
  bind_rows(CDR3length_NM295_Dn14_TCRD) %>%
  bind_rows(CDR3length_NM295_D4_TCRD)

# Separate out the mutate function
# It doesn't work in the pipe if the data.frame has not been created first  
CDR3length_All_TCRD <- CDR3length_All_TCRD %>%
  mutate("Animal.ID.Day" = paste(CDR3length_All_TCRD$Animal.ID, CDR3length_All_TCRD$Day, sep = " "))

# Plot Data and Save
plot_CDR3length_All_TCRD <- ggplot(CDR3length_All_TCRD, 
                                     aes(x = `CDR3 length aa`, 
                                         y = `count`, 
                                         fill = `Animal.ID.Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_All_TCRD

ggsave(paste("All", "TCRD", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# All, TCRG:
# Combine all TCRG data.frames together, and merge Animal.ID and Day in new column
CDR3length_All_TCRG <- CDR3length_NM11_Dn14_TCRG %>%
  bind_rows(CDR3length_NM11_D4_TCRG) %>%
  bind_rows(CDR3length_NM89_Dn14_TCRG) %>% 
  bind_rows(CDR3length_NM89_D4_TCRG) %>%
  bind_rows(CDR3length_NM251_Dn14_TCRG) %>%
  bind_rows(CDR3length_NM251_D4_TCRG) %>%
  bind_rows(CDR3length_NM295_Dn14_TCRG) %>%
  bind_rows(CDR3length_NM295_D4_TCRG)

# Separate out the mutate function
# It doesn't work in the pipe if the data.frame has not been created first  
CDR3length_All_TCRG <- CDR3length_All_TCRG %>%
  mutate("Animal.ID.Day" = paste(CDR3length_All_TCRG$Animal.ID, CDR3length_All_TCRG$Day, sep = " "))

# Plot Data and Save
plot_CDR3length_All_TCRG <- ggplot(CDR3length_All_TCRG, 
                                   aes(x = `CDR3 length aa`, 
                                       y = `count`, 
                                       fill = `Animal.ID.Day`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbbPalette)

plot_CDR3length_All_TCRG

ggsave(paste("All", "TCRG", "CDR3_Spectratyping.tiff", sep = "_"),
       dpi = 1200,
       height = 18,
       width = 10,
       units = "cm",
       limitsize = TRUE)

# Calculating statistics (mean, median, skewness, kurtosis, etc). 
# Skewness and kurtosis to evaluate normal distribution:

# Mean:
mean(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length))

# Median:
median(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length))

# Range:
range(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length))

# IQR:
IQR(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length))

# Kurtosis:
kurtosis(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), type = 2)

# Skewness:
skewness(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), type = 2)

# Number of sequences
length(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length))


# Create data.frames to export Stats:

# TCRD:
Stats_CDR3length_TCRD <- as.data.frame(matrix(0, ncol = 9, nrow = length(unique(CDR3length_All_TCRD$Animal.ID.Day))))
colnames(Stats_CDR3length_TCRD) <- c("Mean CDR3 length (aa)", 
                                     "Median CDR3 length (aa)", 
                                     "Range (smallest) CDR3 length (aa)", 
                                     "Range (largest) CDR3 length (aa)", 
                                     "1st quartile CDR3 length (aa)", 
                                     "3rd quartile CDR3 length (aa)", 
                                     "Kurtosis CDR3 length (aa)", 
                                     "Skewness CDR3 length (aa)",
                                     "Number of sequences"
                                     )
rownames(Stats_CDR3length_TCRD) <- c(unique(CDR3length_All_TCRD$Animal.ID.Day))

# Mean:
Stats_CDR3length_TCRD[, 1] = c(mean(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length)) 
                               )

# Median:
Stats_CDR3length_TCRD[, 2] = c(median(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length)), 
                               median(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length))
                               )

# Range:
Stats_CDR3length_TCRD[1, 3:4] = c(range(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[2, 3:4] = c(range(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[3, 3:4] = c(range(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[4, 3:4] = c(range(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[5, 3:4] = c(range(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[6, 3:4] = c(range(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[7, 3:4] =c(range(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length)))
Stats_CDR3length_TCRD[8, 3:4] =c(range(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length)))

# First quartile
Stats_CDR3length_TCRD[, 5] = c(quantile(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length), probs = 0.25)
                               )

# Third quartile
Stats_CDR3length_TCRD[, 6] = c(quantile(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length), probs = 0.75)
                               )

# Kurtosis:
Stats_CDR3length_TCRD[, 7] = c(kurtosis(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length), type = 2)
                               )

# Skewness:
Stats_CDR3length_TCRD[, 8] = c(skewness(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length), type = 2)
                               )

# Number of sequences:
Stats_CDR3length_TCRD[, 9] = c(length(c(macTCRgd_NM11_Dn14$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM11_D4$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM89_Dn14$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM89_D4$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM251_Dn14$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM251_D4$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM295_Dn14$TCRD.CDR3.aa.length)), 
                               length(c(macTCRgd_NM295_D4$TCRD.CDR3.aa.length))
                               )

# TCRG:
Stats_CDR3length_TCRG <- as.data.frame(matrix(0, ncol = 9, nrow = length(unique(CDR3length_All_TCRG$Animal.ID.Day))))
colnames(Stats_CDR3length_TCRG) <- c("Mean CDR3 length (aa)", 
                                     "Median CDR3 length (aa)", 
                                     "Range (smallest) CDR3 length (aa)", 
                                     "Range (largest) CDR3 length (aa)", 
                                     "1st quartile CDR3 length (aa)", 
                                     "3rd quartile CDR3 length (aa)", 
                                     "Kurtosis CDR3 length (aa)", 
                                     "Skewness CDR3 length (aa)",
                                     "Number of sequences"
                                     )

rownames(Stats_CDR3length_TCRG) <- c(unique(CDR3length_All_TCRG$Animal.ID.Day))

# Mean:
Stats_CDR3length_TCRG[, 1] = c(mean(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length)), 
                               mean(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length)) 
                               )

# Median:
Stats_CDR3length_TCRG[, 2] = c(median(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length)), 
                               median(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length))
                               )

# Range:
Stats_CDR3length_TCRG[1, 3:4] = c(range(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[2, 3:4] = c(range(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[3, 3:4] = c(range(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[4, 3:4] = c(range(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[5, 3:4] = c(range(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[6, 3:4] = c(range(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[7, 3:4] =c(range(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length)))
Stats_CDR3length_TCRG[8, 3:4] =c(range(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length)))

# First quartile
Stats_CDR3length_TCRG[, 5] = c(quantile(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length), probs = 0.25), 
                               quantile(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length), probs = 0.25)
                               )

# Third quartile
Stats_CDR3length_TCRG[, 6] = c(quantile(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length), probs = 0.75), 
                               quantile(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length), probs = 0.75)
                               )

# Kurtosis:
Stats_CDR3length_TCRG[, 7] = c(kurtosis(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               kurtosis(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length), type = 2)
                               )

# Skewness:
Stats_CDR3length_TCRG[, 8] = c(skewness(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length), type = 2), 
                               skewness(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length), type = 2)
                               )

# Number of sequences:
Stats_CDR3length_TCRG[, 9] = c(length(c(macTCRgd_NM11_Dn14$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM11_D4$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM89_Dn14$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM89_D4$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM251_Dn14$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM251_D4$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM295_Dn14$TCRG.CDR3.aa.length)), 
                               length(c(macTCRgd_NM295_D4$TCRG.CDR3.aa.length)))


# Stats_CDR3length_TCRD CSV export:
write.csv(Stats_CDR3length_TCRD, "Spectratypinng_Statistics_TCRD.csv")

# Stats_CDR3length_TCRG CSV export:
write.csv(Stats_CDR3length_TCRG, "Spectratypinng_Statistics_TCRG.csv")