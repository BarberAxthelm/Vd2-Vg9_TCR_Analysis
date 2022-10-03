## TCR_Data_Preparation_MiXCR.R
# This script is used to prepare the data for use by the other analysis scripts. 
# It imports a CSV file containing data formatted appropriately post-MiXCR analysis,
# as discussed in the readme, and a CSV containing the MiXCR list of functional clones.
# It then filters the data to only MiXCR assembled functional clones to ensure the 
# analysis only occurs with productive TCRDV2/TCRGV9 paired CDR3 sequences, removes 
# all unused columns, and exports the data.
# It exports two tables: 
#    (1) MiXCR Data, filtered to include only paired functional sequences and 
#                    related sample infomation and flow cytometery data.
#    (2) Counts of functional paired sequences per sample.

# Removes ALL objects from the global environment
rm(list = ls())

# Intalling and loading all packages
install.packages("tidyverse")
install.packages("data.table")
library(tidyverse)
library(data.table)

## Import CSV ##

# Import combined AA sequence information with FACS data (note data order must match column names)
# Also renames column titles approprietly and converts the table to a data.frame
macTCRgd <- read_csv(
  file.choose(new = FALSE),
  col_names = c("TCRD Sequence number",
                "TCRD Best V Gene",
                "TCRD Best D Gene",
                "TCRD Best J Gene",
                "TCRD Best C Gene",
                "TCRD All V Genes",
                "TCRD All D Genes",
                "TCRD All J Genes",
                "TCRD All C Gene",
                "TCRD Best V Allele",
                "TCRD Best D Allele",
                "TCRD Best J Allele",
                "TCRD Best C Allele",
                "TCRD Best V Allele score",
                "TCRD Best D Allele score",
                "TCRD Best J Allele score",
                "TCRD Best C Allele score",
                "TCRD All V Allele score",
                "TCRD All D Allele score",
                "TCRD All J Allele score",
                "TCRD All C Allele score",
                "TCRD FR1 nt",
                "TCRD CDR1 nt",
                "TCRD FR2 nt",
                "TCRD CDR2 nt",
                "TCRD FR3 nt",
                "TCRD CDR3 nt",
                "TCRD FR4 nt",
                "TCRD FR1 nt quality",
                "TCRD CDR1 nt quality",
                "TCRD FR2 nt quality",
                "TCRD CDR2 nt quality",
                "TCRD FR3 nt quality",
                "TCRD CDR3 nt quality",
                "TCRD FR4 nt quality",
                "TCRD FR1 aa",
                "TCRD CDR1 aa",
                "TCRD FR2 aa",
                "TCRD CDR2 aa",
                "TCRD FR3 aa",
                "TCRD CDR3 aa",
                "TCRD FR4 aa",
                "TCRD Reference Points",
                "TCRD Sequence nt",
                "TCRD Sequence nt quality",
                "TCRD V Identity Percentage",
                "TCRD D Identity Percentage",
                "TCRD J Identity Percentage",
                "TCRD C Identity Percentage",
                "TCRD V Best Identity Percentage",
                "TCRD D Best Identity Percentage",
                "TCRD J Best Identity Percentage",
                "TCRD C Best Identity Percentage",
                "TCRD Sequence ID",
                "TCRD Sequence Length nt",
                "TCRG Sequence number",
                "TCRG Best V Gene",
                "TCRG Best D Gene",
                "TCRG Best J Gene",
                "TCRG Best C Gene",
                "TCRG All V Genes",
                "TCRG All D Genes",
                "TCRG All J Genes",
                "TCRG All C Gene",
                "TCRG Best V Allele",
                "TCRG Best D Allele",
                "TCRG Best J Allele",
                "TCRG Best C Allele",
                "TCRG Best V Allele score",
                "TCRG Best D Allele score",
                "TCRG Best J Allele score",
                "TCRG Best C Allele score",
                "TCRG All V Allele score",
                "TCRG All D Allele score",
                "TCRG All J Allele score",
                "TCRG All C Allele score",
                "TCRG FR1 nt",
                "TCRG CDR1 nt",
                "TCRG FR2 nt",
                "TCRG CDR2 nt",
                "TCRG FR3 nt",
                "TCRG CDR3 nt",
                "TCRG FR4 nt",
                "TCRG FR1 nt quality",
                "TCRG CDR1 nt quality",
                "TCRG FR2 nt quality",
                "TCRG CDR2 nt quality",
                "TCRG FR3 nt quality",
                "TCRG CDR3 nt quality",
                "TCRG FR4 nt quality",
                "TCRG FR1 aa",
                "TCRG CDR1 aa",
                "TCRG FR2 aa",
                "TCRG CDR2 aa",
                "TCRG FR3 aa",
                "TCRG CDR3 aa",
                "TCRG FR4 aa",
                "TCRG Reference Points",
                "TCRG Sequence nt",
                "TCRG Sequence nt quality",
                "TCRG V Identity Percentage",
                "TCRG D Identity Percentage",
                "TCRG J Identity Percentage",
                "TCRG C Identity Percentage",
                "TCRG V Best Identity Percentage",
                "TCRG D Best Identity Percentage",
                "TCRG J Best Identity Percentage",
                "TCRG C Best Identity Percentage",
                "TCRG Sequence ID",
                "TCRG Sequence Length nt",
                "Animal ID",
                "Day",
                "Plate number",
                "Well position",
                "Vg9 FITC MFI",
                "Compensated Vg9 FITC MFI",
                "Compensated NKG2A APC MFI",
                "Compensated CD3 Alexa700 MFI",
                "Compensated CD69 APCFire750 MFI",
                "Compensated Dump MFI",
                "Compensated DNAM BV605 MFI",
                "Compensated CD28 BV711 MFI",
                "Compensated CCR6 BV786 MFI",
                "Compensated Vd2 PE MFI",
                "Compensated CD26 PETR MFI",
                "Compensated CD95 PECy7 MFI",
                "FSC-A",
                "FSC-H",
                "FSC-W",
                "IdxCol",
                "IdxRow",
                "NKG2A APC MFI",
                "CD3 Alexa700 MFI",
                "CD69 APCFire750 MFI",
                "SSC-A",
                "SSC-H",
                "SSC-W",
                "FACS Time",
                "Dump MFI",
                "DNAM BV605 MFI",
                "CD28 BV711 MFI",
                "CCR6 BV786 MFI",
                "Vd2 PE MFI",
                "CD26 PETR MFI",
                "CD95 PECy7 MFI",
                "CCR6pos",
                "CD69pos",
                "DNAMpos",
                "Q1_NKG2AnegCD26pos",
                "Q1_NKG2AnegCD26pos CCR6pos",
                "Q1_NKG2AnegCD26pos CD69pos",
                "Q1_NKG2AnegCD26pos DNAMpos",
                "Q1_NKG2AnegCD26pos CD28posCD95neg",
                "Q1_NKG2AnegCD26pos CD28posCD95pos",
                "Q1_NKG2AnegCD26pos CD28negCD95pos",
                "Q1_NKG2AnegCD26pos CD28negCD95neg",
                "Q2_NKG2AposCD26pos",
                "Q2_NKG2AposCD26pos CCR6pos",
                "Q2_NKG2AposCD26pos CD69pos",
                "Q2_NKG2AposCD26pos DNAMpos",
                "Q2_NKG2AposCD26pos CD28posCD95neg",
                "Q2_NKG2AposCD26pos CD28posCD95pos",
                "Q2_NKG2AposCD26pos CD28negCD95pos",
                "Q2_NKG2AposCD26pos CD28negCD95neg",
                "Q3_NKG2AposCD26neg",
                "Q3_NKG2AposCD26neg CCR6pos",
                "Q3_NKG2AposCD26neg CD69pos",
                "Q3_NKG2AposCD26neg DNAMpos",
                "Q3_NKG2AposCD26neg CD28posCD95neg",
                "Q3_NKG2AposCD26neg CD28posCD95pos",
                "Q3_NKG2AposCD26neg CD28negCD95pos",
                "Q3_NKG2AposCD26neg CD28negCD95neg",
                "Q4_NKG2AnegCD26neg",
                "Q4_NKG2AnegCD26neg CCR6pos",
                "Q4_NKG2AnegCD26neg CD69pos",
                "Q4_NKG2AnegCD26neg DNAMpos",
                "Q4_NKG2AnegCD26neg CD28posCD95neg",
                "Q4_NKG2AnegCD26neg CD28posCD95pos",
                "Q4_NKG2AnegCD26neg CD28negCD95pos",
                "Q4_NKG2AnegCD26neg CD28negCD95neg",
                "CD28posCD95neg",
                "CD28posCD95pos",
                "CD28negCD95pos",
                "CD28negCD95neg"
  )) %>%
  data.frame()


## Filter and Format Data ##

# Removes the first row from the data.frame (which contains the headers from the csv due to the import method)
macTCRgd <- macTCRgd[-1, ]

# Converts characters to numeric values (note: Plate number won't convert to numeric)
macTCRgd <- macTCRgd %>%
  type_convert()

# Removes the columns that aren't needed (ex. no data from the MiXCR analysis, non-usable FACS data)
macTCRgd <- subset(macTCRgd, select = -c(TCRD.Best.C.Gene,
                                         TCRD.All.C.Gene,
                                         TCRD.Best.C.Allele,
                                         TCRD.Best.C.Allele.score,
                                         TCRD.All.C.Allele.score,
                                         TCRD.FR1.nt,
                                         TCRD.CDR1.nt,
                                         TCRD.FR2.nt,
                                         TCRD.CDR2.nt,
                                         TCRD.FR1.nt.quality,
                                         TCRD.CDR1.nt.quality,
                                         TCRD.FR2.nt.quality,
                                         TCRD.CDR2.nt.quality,
                                         TCRD.FR1.aa,
                                         TCRD.CDR1.aa,
                                         TCRD.FR2.aa,
                                         TCRD.CDR2.aa,
                                         TCRD.C.Identity.Percentage,
                                         TCRD.C.Best.Identity.Percentage,
                                         TCRG.Best.D.Gene,
                                         TCRG.Best.C.Gene,
                                         TCRG.All.D.Genes,
                                         TCRG.All.C.Gene,
                                         TCRG.Best.D.Allele,
                                         TCRG.Best.C.Allele,
                                         TCRG.Best.D.Allele.score,
                                         TCRG.Best.C.Allele.score,
                                         TCRG.All.D.Allele.score,
                                         TCRG.All.C.Allele.score,
                                         TCRG.FR1.nt,
                                         TCRG.CDR1.nt,
                                         TCRG.FR2.nt,
                                         TCRG.FR1.nt.quality,
                                         TCRG.CDR1.nt.quality,
                                         TCRG.FR2.nt.quality,
                                         TCRG.FR1.aa,
                                         TCRG.CDR1.aa,
                                         TCRG.FR2.aa,
                                         TCRG.D.Identity.Percentage,
                                         TCRG.C.Identity.Percentage,
                                         TCRG.D.Best.Identity.Percentage,
                                         TCRG.C.Best.Identity.Percentage,
                                         IdxCol,
                                         IdxRow))

# Filters sequences which are not TRDV2 and TRGV9
macTCRgd <- filter(macTCRgd, TCRD.Best.V.Gene == "TRDV2",
                   TCRG.Best.V.Gene == "TRGV9")

# Filters unpaired TRDV2/TRGV9 sequences
macTCRgd <- macTCRgd %>%
  filter(!is.na(TCRD.Best.V.Gene)) %>%
  filter(!is.na(TCRG.Best.V.Gene))

# Filters blank TRDV2/TRGV9 CDR3 sequences
macTCRgd <- macTCRgd %>%
  filter(!is.na(TCRD.CDR3.aa)) %>%
  filter(!is.na(TCRG.CDR3.aa))


## Process to Find Functional Clones ##

# Columns to indicate which CD3R3s are functional
macTCRgd["TCRD.functional"] <- 0
macTCRgd["TCRG.functional"] <- 0

# Imports MiXCR assembled functional clones, a list of functional TCRD CDR3 clones. 
# (TCRG clones imported below)
# We will use this to filter data to only use productive TCRD CDR3 Sequences with downstream analysis
TCRD_functional_list <- read_csv(
  file.choose(new = FALSE),
  col_names = c("Clone ID",
                "Clone count",
                "Clone fraction",
                "Target sequences",
                "Target quality",
                "All V hits with score",
                "All D hits with score",
                "All J hits with score",
                "All C hits with score",
                "All V alignments",
                "All D alignments",
                "All J alignments",
                "All C alignments",
                "nSeqFR1",
                "minQualFR1",
                "nSeqCDR1",
                "minQualCDR1",
                "nSeqFR2",
                "minQualFR2",
                "nSeqCDR2",
                "minQualCDR2",
                "nSeqFR3",
                "minQualFR3",
                "nSeqCDR3",
                "minQualCDR3",
                "nSeqFR4",
                "minQualFR4",
                "aaSeqFR1",
                "aaSeqCDR1",
                "aaSeqFR2",
                "aaSeqCDR2",
                "aaSeqFR3",
                "aaSeqCDR3",
                "aaSeqFR4",
                "Ref points"
  )) %>%
  data.frame()

# Removes the first row from the data.frame (which is the headers from the csv due to the import method)
TCRD_functional_list <- TCRD_functional_list[-1, ]

# Converts characters to numeric values
TCRD_functional_list <- TCRD_functional_list %>%
  type_convert()

# Removes unused columns (ex. no data from the MiXCR analysis, non-usable FACS data)
TCRD_functional_list <- subset(TCRD_functional_list, select = -c(All.C.hits.with.score,
                                                                 All.C.alignments,
                                                                 nSeqFR1,
                                                                 minQualFR1,
                                                                 nSeqCDR1,
                                                                 minQualCDR1,
                                                                 nSeqFR2,
                                                                 minQualFR2,
                                                                 nSeqCDR2,
                                                                 minQualCDR2,
                                                                 nSeqFR3,
                                                                 minQualFR3,
                                                                 nSeqFR4,
                                                                 minQualFR4,
                                                                 aaSeqFR1,
                                                                 aaSeqCDR1,
                                                                 aaSeqFR2,
                                                                 aaSeqCDR2,
                                                                 aaSeqFR3,
                                                                 aaSeqFR4
))

# Steps through macTCRgd and TCRD_functional_list and flags which clones are functional
for (i in 1:nrow(macTCRgd)) {
  for (j in 1:nrow(TCRD_functional_list)) {
    if (macTCRgd$TCRD.CDR3.aa[i] == TCRD_functional_list$aaSeqCDR3[j]) {
      macTCRgd$TCRD.functional[i] <- 1
    }
  }
}

# Imports MiXCR assembled functional clones, a list of functional TCRG CDR3 clones.
# We will use this to filter data to only use productive TCRG CDR3 Sequences with downstream analysis
TCRG_functional_list <- read_csv(
  file.choose(new = FALSE),
  col_names = c("Clone ID",
                "Clone count",
                "Clone fraction",
                "Target sequences",
                "Target quality",
                "All V hits with score",
                "All D hits with score",
                "All J hits with score",
                "All C hits with score",
                "All V alignments",
                "All D alignments",
                "All J alignments",
                "All C alignments",
                "nSeqFR1",
                "minQualFR1",
                "nSeqCDR1",
                "minQualCDR1",
                "nSeqFR2",
                "minQualFR2",
                "nSeqCDR2",
                "minQualCDR2",
                "nSeqFR3",
                "minQualFR3",
                "nSeqCDR3",
                "minQualCDR3",
                "nSeqFR4",
                "minQualFR4",
                "aaSeqFR1",
                "aaSeqCDR1",
                "aaSeqFR2",
                "aaSeqCDR2",
                "aaSeqFR3",
                "aaSeqCDR3",
                "aaSeqFR4",
                "Ref points"
  )) %>%
  data.frame()

# Removes the first row from the data.frame (which is the headers from the csv due to the importation method)
TCRG_functional_list <- TCRG_functional_list[-1, ]

# Converts characters to numeric values
TCRG_functional_list <- TCRG_functional_list %>%
  type_convert()

# Removes unused columns (ex. no data from the MiXCR analysis, non-usable FACS data)
TCRG_functional_list <- subset(TCRG_functional_list, select = -c(All.D.hits.with.score,
                                                                 All.C.hits.with.score,
                                                                 All.D.alignments,
                                                                 All.C.alignments,
                                                                 nSeqFR1,
                                                                 minQualFR1,
                                                                 nSeqCDR1,
                                                                 minQualCDR1,
                                                                 nSeqFR2,
                                                                 minQualFR2,
                                                                 nSeqCDR2,
                                                                 minQualCDR2,
                                                                 nSeqFR3,
                                                                 minQualFR3,
                                                                 nSeqFR4,
                                                                 minQualFR4,
                                                                 aaSeqFR1,
                                                                 aaSeqCDR1,
                                                                 aaSeqFR2,
                                                                 aaSeqCDR2,
                                                                 aaSeqFR3,
                                                                 aaSeqFR4
))

# Steps through macTCRgd and TCRG_functional_list and flags which clones are functional
for (i in 1:nrow(macTCRgd)) {
  for (j in 1:nrow(TCRG_functional_list)) {
    if (macTCRgd$TCRG.CDR3.aa[i] == TCRG_functional_list$aaSeqCDR3[j]) {
      macTCRgd$TCRG.functional[i] <- 1
    }
  }
}


# Create lists of excluded TCRD and TCRG cells, and filter the main list to include only functional cells

macTCRgd_excludedTCRD <- filter(macTCRgd, TCRD.functional == 0)
macTCRgd_excludedTCRG <- filter(macTCRgd, TCRG.functional == 0)

macTCRgd <- filter(macTCRgd, TCRD.functional == 1,
                   TCRG.functional == 1)

# Calculate CDR3 nucleotide (nt) and amino acid (aa) lengths
macTCRgd <- macTCRgd %>%
  add_column("TCRD.CDR3.aa.length" = 0) %>%
  add_column("TCRG.CDR3.aa.length" = 0) %>%
  add_column("TCRD.CDR3.nt.length" = 0) %>%
  add_column("TCRG.CDR3.nt.length" = 0)

for (i in 1:nrow(macTCRgd)){
  macTCRgd$TCRD.CDR3.aa.length[i] <- nchar(macTCRgd$TCRD.CDR3.aa[i])
  macTCRgd$TCRG.CDR3.aa.length[i] <- nchar(macTCRgd$TCRG.CDR3.aa[i])
  macTCRgd$TCRD.CDR3.nt.length[i] <- nchar(macTCRgd$TCRD.CDR3.nt[i])
  macTCRgd$TCRG.CDR3.nt.length[i] <- nchar(macTCRgd$TCRG.CDR3.nt[i])
}

# Confirm and ensures all additions are numeric.
macTCRgd <- macTCRgd %>%
  type_convert()

# Concatentates the TCRD and TCRG Junctions for amino acids
macTCRgd <- macTCRgd %>%
  mutate(
    Combined.CDR3 = paste(TCRD.CDR3.aa, TCRG.CDR3.aa, sep = "-")
  )

# Creating counts for TCRG/TCRD chain recovery for reference:
# Note day 15 is unused.

# NM11:
NM11 <- filter(macTCRgd, Animal.ID == "NM11")
NM11_Dn14 <- filter(NM11, Day == -14)
NM11_D4 <- filter(NM11, Day == 4)
#NM11_D15 <- filter(NM11, Day == 15)

# NM89:
NM89 <- filter(macTCRgd, Animal.ID == "NM89")
NM89_Dn14 <- filter(NM89, Day == -14)
NM89_D4 <- filter(NM89, Day == 4)
#NM89_D15 <- filter(NM89, Day == 15)

# NM251:
NM251 <- filter(macTCRgd, Animal.ID == "NM251")
NM251_Dn14 <- filter(NM251, Day == -14)
NM251_D4 <- filter(NM251, Day == 4)
#NM251_D15 <- filter(NM251, Day == 15)

# NM295:
NM295 <- filter(macTCRgd, Animal.ID == "NM295")
NM295_Dn14 <- filter(NM295, Day == -14)
NM295_D4 <- filter(NM295, Day == 4)
#NM295_D15 <- filter(NM295, Day == 15)

# Create a table with chain counts for each timepoint to reference:
Chain_Recovery_Counts <- as.data.frame(matrix(0, ncol = 1, nrow = 12))
colnames(Chain_Recovery_Counts) <- c("counts")
rownames(Chain_Recovery_Counts) <- c("NM11 Day -14",
                                     "NM11 Day 4",
                                     "NM11 Day 15",
                                     "NM89 Day -14",
                                     "NM89 Day 4",
                                     "NM89 Day 15",
                                     "NM251 Day -14",
                                     "NM251 Day 4",
                                     "NM251 Day 15",
                                     "NM295 Day -14",
                                     "NM295 Day 4",
                                     "NM295 Day 15")

# Counts of functional paired sequences per channel
# This means, both the TCRDV2 and TCRGV9 are functional and paired
Chain_Recovery_Counts[1, 1] <- nrow(NM11_Dn14)
Chain_Recovery_Counts[2, 1] <- nrow(NM11_D4)
#Chain_Recovery_Counts[3, 1] <- nrow(NM11_D15)
Chain_Recovery_Counts[4, 1] <- nrow(NM89_Dn14)
Chain_Recovery_Counts[5, 1] <- nrow(NM89_D4)
#Chain_Recovery_Counts[6, 1] <- nrow(NM89_D15)
Chain_Recovery_Counts[7, 1] <- nrow(NM251_Dn14)
Chain_Recovery_Counts[8, 1] <- nrow(NM251_D4)
#Chain_Recovery_Counts[9, 1] <- nrow(NM251_D15)
Chain_Recovery_Counts[10, 1] <- nrow(NM295_Dn14)
Chain_Recovery_Counts[11, 1] <- nrow(NM295_D4)
#Chain_Recovery_Counts[12, 1] <- nrow(NM295_D15)

# Delete data.frames for individual animals/days.
# Ensures they don't interfere with downstream analysis for the individual pipeline
# This will export the data to be used elsewhere.
rm(NM11,
   NM11_Dn14,
   NM11_D4,
   NM11_D15,
   NM89,
   NM89_Dn14,
   NM89_D4,
   NM89_D15,
   NM251,
   NM251_Dn14,
   NM251_D4,
   NM251_D15,
   NM295,
   NM295_Dn14,
   NM295_D4, 
   NM295_D15
)


## Export Formatted Data ## 

# Set the working directory
setwd("Outputs")

# macTCRgd: filtered to include only paired functional sequences and their sample info, and flow cytometery data.
write.csv(macTCRgd, "macTCRgd_Filtered.csv")

# Chain_Recovery_Counts
write.csv(Chain_Recovery_Counts, "Chain_recovery_macTCRgd_Filtered.csv")

# Reset Working Directory
setwd("..")
