## TCR_MSA_Sequence_Logos.R
# This script imports the formatted data (output of TCR_Data_Preparation_MiXCR.R)
# and filters by animal/timepoint to perform a multiple sequence alignment (msa)
# and a sequence logo generation (aa and nt) for both the gamma and delta chains.  
# The output is in both fasta and pdf file formats and labeled by animal, timepoint, and chain.

# Removes ALL objects from the global environment
rm(list = ls())

# Set the working directory
setwd("Outputs")

# Installs the msa package via Bioconductor, as per demo in the msa user manual
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")


# Installs and Loads Packages.
install.packages("tinytex")
tinytex::install_tinytex()
# To uninstall TinyTeX, run:
# tinytex::uninstall_tinytex(force = TRUE)
tinytex::parse_install("texinfo")
install.packages("tidyverse")
# Shows where the texshade file is located on the computer, for potential troubleshooting
system.file("tex", "texshade.sty", package = "msa")
library(msa)
library("tinytex")
library(tidyverse)


## Loading MiXCR Data ##

macTCRgd <- read_csv(
  file.choose(new = FALSE)) %>%
  data.frame()

## Start of Analysis ##

# Filter by Animal.ID
macTCRgd_NM11 <- filter(macTCRgd, Animal.ID == "NM11")
macTCRgd_NM89 <- filter(macTCRgd, Animal.ID == "NM89")
macTCRgd_NM251 <- filter(macTCRgd, Animal.ID == "NM251")
macTCRgd_NM295 <- filter(macTCRgd, Animal.ID == "NM295")

# Filter the Animal ID's by Day
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


## Alignments and Sequence logo generation ##

## NM11 ##

# NM11 D4 TCRD Amino Acids
macTCRgd_NM11_D4_TCRD_aa <- msa(as.vector(macTCRgd_NM11_D4$TCRD.CDR3.aa),
                       type = "Protein",
                       method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_D4_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top", # Specifies if you want to show a sequence logo or not
               logoColors = "rasmol", # Sets the colors for the sequence logo
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_D4_TCRD_aa.pdf",
               alFile = "macTCRgd_NM11_D4_TCRD_aa.fasta",
               verbose = FALSE)

# NM11 D4 TCRD Nucleotides
macTCRgd_NM11_D4_TCRD_nt <- msa(as.vector(macTCRgd_NM11_D4$TCRD.CDR3.nt),
                           type = "dna",
                           method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_D4_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top", # Specifies if you want to show a sequence logo or not
               logoColors = "rasmol", # Sets the colors for the sequence logo
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_D4_TCRD_nt.pdf",
               alFile = "macTCRgd_NM11_D4_TCRD_nt.fasta",
               verbose = FALSE)

# NM11 D -14 TCRD Amino Acids
macTCRgd_NM11_Dn14_TCRD_aa <- msa(as.vector(macTCRgd_NM11_Dn14$TCRD.CDR3.aa),
                           type = "Protein",
                           method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_Dn14_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_TCRD_Dn14_aa.pdf",
               alFile = " macTCRgd_NM11_TCRD_Dn14_aa.fasta",
               verbose = FALSE)

# NM11 D -14 TCRD Nucleotides
macTCRgd_NM11_Dn14_TCRD_nt <- msa(as.vector(macTCRgd_NM11_Dn14$TCRD.CDR3.nt),
                           type = "dna",
                           method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_Dn14_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_TCRD_Dn14_nt.pdf",
               alFile = "macTCRgd_NM11_TCRD_Dn14_nt.fasta",
               verbose = FALSE)

# NM11 D11 TCRG Amino Acids
macTCRgd_NM11_D4_TCRG_aa <- msa(as.vector(macTCRgd_NM11_D4$TCRG.CDR3.aa),
                           type = "Protein",
                           method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_D4_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_D4_TCRG_aa.pdf",
               alFile = "macTCRgd_NM11_D4_TCRG_aa.fasta",
               verbose = FALSE)

# NM11 D4 TCRG Nucleotides
macTCRgd_NM11_D4_TCRG_nt <- msa(as.vector(macTCRgd_NM11_D4$TCRG.CDR3.nt),
                           type = "dna",
                           method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_D4_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_D4_TCRG_nt.pdf",
               alFile = "macTCRgd_NM11_D4_TCRG_nt.fasta",
               verbose = FALSE)

# NM11 D-14 TCRG Amino Acids
macTCRgd_NM11_Dn14_TCRG_aa <- msa(as.vector(macTCRgd_NM11_Dn14$TCRG.CDR3.aa),
                             type = "Protein",
                             method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_Dn14_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_Dn14_TCRG_aa.pdf",
               alFile = "macTCRgd_NM11_Dn14_TCRG_aa.fasta",
               verbose = FALSE)

# NM11 D-14 TCRG Nucleotides
macTCRgd_NM11_Dn14_TCRG_nt <- msa(as.vector(macTCRgd_NM11_Dn14$TCRD.CDR3.nt),
                             type = "dna",
                             method = "ClustalW")
msaPrettyPrint(macTCRgd_NM11_Dn14_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM11_Dn14_TCRG_nt.pdf",
               alFile = "macTCRgd_NM11_Dn14_TCRG_nt.fasta",
               verbose = FALSE)


## NM89 ##

# NM89 D4 TCRD Amino Acids
macTCRgd_NM89_D4_TCRD_aa <- msa(as.vector(macTCRgd_NM89_D4$TCRD.CDR3.aa),
                                type = "Protein",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_D4_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_D4_TCRD_aa.pdf",
               alFile = "macTCRgd_NM89_D4_TCRD_aa.fasta",
               verbose = FALSE)

# NM89 D4 TCRD Nucleotides
macTCRgd_NM89_D4_TCRD_nt <- msa(as.vector(macTCRgd_NM89_D4$TCRD.CDR3.nt),
                                type = "dna",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_D4_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_D4_TCRD_nt.pdf",
               alFile = "macTCRgd_NM89_D4_TCRD_nt.fasta",
               verbose = FALSE)

# NM89 D-14 TCRD Amino Acids
macTCRgd_NM89_Dn14_TCRD_aa <- msa(as.vector(macTCRgd_NM89_Dn14$TCRD.CDR3.aa),
                                  type = "Protein",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_Dn14_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_TCRD_Dn14_aa.pdf",
               alFile = "macTCRgd_NM89_TCRD_Dn14_aa.fasta",
               verbose = FALSE)

# NM89 D-14 TCRD Nucleotides
macTCRgd_NM89_Dn14_TCRD_nt <- msa(as.vector(macTCRgd_NM89_Dn14$TCRD.CDR3.nt),
                                  type = "dna",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_Dn14_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_TCRD_Dn14_nt.pdf",
               alFile = "macTCRgd_NM89_TCRD_Dn14_nt.fasta",
               verbose = FALSE)

# NM89 D4 TCRG Amino Acids
macTCRgd_NM89_D4_TCRG_aa <- msa(as.vector(macTCRgd_NM89_D4$TCRG.CDR3.aa),
                                type = "Protein",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_D4_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_D4_TCRG_aa.pdf",
               alFile = "macTCRgd_NM89_D4_TCRG_aa.fasta",
               verbose = FALSE)

# NM89 D4 TCRG Nucleotides
macTCRgd_NM89_D4_TCRG_nt <- msa(as.vector(macTCRgd_NM89_D4$TCRG.CDR3.nt),
                                type = "dna",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_D4_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_D4_TCRG_nt.pdf",
               alFile = "macTCRgd_NM89_D4_TCRG_nt.fasta",
               verbose = FALSE)

# NM89 D-14 TCRG Amino Acids
macTCRgd_NM89_Dn14_TCRG_aa <- msa(as.vector(macTCRgd_NM89_Dn14$TCRG.CDR3.aa),
                                  type = "Protein",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_Dn14_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_Dn14_TCRG_aa.pdf",
               alFile = "macTCRgd_NM89_Dn14_TCRG_aa.fasta",
               verbose = FALSE)

# NM89 D-14 TCRG Nucleotides
macTCRgd_NM89_Dn14_TCRG_nt <- msa(as.vector(macTCRgd_NM89_Dn14$TCRD.CDR3.nt),
                                  type = "dna",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM89_Dn14_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM89_Dn14_TCRG_nt.pdf",
               alFile = "macTCRgd_NM89_Dn14_TCRG_nt.fasta",
               verbose = FALSE)

## NM251 ##

# NM251 D4 TCRD Amino Acids
macTCRgd_NM251_D4_TCRD_aa <- msa(as.vector(macTCRgd_NM251_D4$TCRD.CDR3.aa),
                                type = "Protein",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_D4_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_D4_TCRD_aa.pdf",
               alFile = "macTCRgd_NM251_D4_TCRD_aa.fasta",
               verbose = FALSE)

# NM251 D4 TCRD Nucleotides
macTCRgd_NM251_D4_TCRD_nt <- msa(as.vector(macTCRgd_NM251_D4$TCRD.CDR3.nt),
                                type = "dna",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_D4_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_D4_TCRD_nt.pdf",
               alFile = "macTCRgd_NM251_D4_TCRD_nt.fasta",
               verbose = FALSE)

# NM251 D-14 TCRD Amino Acids
macTCRgd_NM251_Dn14_TCRD_aa <- msa(as.vector(macTCRgd_NM251_Dn14$TCRD.CDR3.aa),
                                  type = "Protein",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_Dn14_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_TCRD_Dn14_aa.pdf",
               alFile = "macTCRgd_NM251_TCRD_Dn14_aa.fasta",
               verbose = FALSE)

# NM251 D-14 TCRD Nucleotides
macTCRgd_NM251_Dn14_TCRD_nt <- msa(as.vector(macTCRgd_NM251_Dn14$TCRD.CDR3.nt),
                                  type = "dna",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_Dn14_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_TCRD_Dn14_nt.pdf",
               alFile = "macTCRgd_NM251_TCRD_Dn14_nt.fasta",
               verbose = FALSE)

# NM251 D4 TCRG Amino Acids
macTCRgd_NM251_D4_TCRG_aa <- msa(as.vector(macTCRgd_NM251_D4$TCRG.CDR3.aa),
                                type = "Protein",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_D4_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_D4_TCRG_aa.pdf",
               alFile = "macTCRgd_NM251_D4_TCRG_aa.fasta",
               verbose = FALSE)

# NM251 D4 TCRG Nucleotides
macTCRgd_NM251_D4_TCRG_nt <- msa(as.vector(macTCRgd_NM251_D4$TCRG.CDR3.nt),
                                type = "dna",
                                method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_D4_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_D4_TCRG_nt.pdf",
               alFile = "macTCRgd_NM251_D4_TCRG_nt.fasta",
               verbose = FALSE)

# NM251 D-14 TCRG Amino Acids
macTCRgd_NM251_Dn14_TCRG_aa <- msa(as.vector(macTCRgd_NM251_Dn14$TCRG.CDR3.aa),
                                  type = "Protein",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_Dn14_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_Dn14_TCRG_aa.pdf",
               alFile = "macTCRgd_NM251_Dn14_TCRG_aa.fasta",
               verbose = FALSE)

# NM251 D-14 TCRG Nucleotides
macTCRgd_NM251_Dn14_TCRG_nt <- msa(as.vector(macTCRgd_NM251_Dn14$TCRD.CDR3.nt),
                                  type = "dna",
                                  method = "ClustalW")
msaPrettyPrint(macTCRgd_NM251_Dn14_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM251_Dn14_TCRG_nt.pdf",
               alFile = "macTCRgd_NM251_Dn14_TCRG_nt.fasta",
               verbose = FALSE)


## NM295 ##

# NM295 D4 TCRD Amino Acids
macTCRgd_NM295_D4_TCRD_aa <- msa(as.vector(macTCRgd_NM295_D4$TCRD.CDR3.aa),
                                 type = "Protein",
                                 method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_D4_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_D4_TCRD_aa.pdf",
               alFile = "macTCRgd_NM295_D4_TCRD_aa.fasta",
               verbose = FALSE)

# NM295 D4 TCRD Nucleotides
macTCRgd_NM295_D4_TCRD_nt <- msa(as.vector(macTCRgd_NM295_D4$TCRD.CDR3.nt),
                                 type = "dna",
                                 method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_D4_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_D4_TCRD_nt.pdf",
               alFile = "macTCRgd_NM295_D4_TCRD_nt.fasta",
               verbose = FALSE)

# NM295 D-14 TCRD Amino Acids
macTCRgd_NM295_Dn14_TCRD_aa <- msa(as.vector(macTCRgd_NM295_Dn14$TCRD.CDR3.aa),
                                   type = "Protein",
                                   method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_Dn14_TCRD_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_TCRD_Dn14_aa.pdf",
               alFile = "macTCRgd_NM295_TCRD_Dn14_aa.fasta",
               verbose = FALSE)

# NM295 D-14 TCRD Nucleotides
macTCRgd_NM295_Dn14_TCRD_nt <- msa(as.vector(macTCRgd_NM295_Dn14$TCRD.CDR3.nt),
                                   type = "dna",
                                   method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_Dn14_TCRD_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_TCRD_Dn14_nt.pdf",
               alFile = "macTCRgd_NM295_TCRD_Dn14_nt.fasta",
               verbose = FALSE)

# NM295 D4 TCRG Amino Acids
macTCRgd_NM295_D4_TCRG_aa <- msa(as.vector(macTCRgd_NM295_D4$TCRG.CDR3.aa),
                                 type = "Protein",
                                 method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_D4_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_D4_TCRG_aa.pdf",
               alFile = "macTCRgd_NM295_D4_TCRG_aa.fasta",
               verbose = FALSE)

# NM295 D4 TCRG Nucleotides
macTCRgd_NM295_D4_TCRG_nt <- msa(as.vector(macTCRgd_NM295_D4$TCRG.CDR3.nt),
                                 type = "dna",
                                 method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_D4_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_D4_TCRG_nt.pdf",
               alFile = "macTCRgd_NM295_D4_TCRG_nt.fasta",
               verbose = FALSE)

# NM295 D-14 TCRG Amino Acids
macTCRgd_NM295_Dn14_TCRG_aa <- msa(as.vector(macTCRgd_NM295_Dn14$TCRG.CDR3.aa),
                                   type = "Protein",
                                   method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_Dn14_TCRG_aa, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_Dn14_TCRG_aa.pdf",
               alFile = "macTCRgd_NM295_Dn14_TCRG_aa.fasta",
               verbose = FALSE)

# NM295 D-14 TCRG Nucleotides
macTCRgd_NM295_Dn14_TCRG_nt <- msa(as.vector(macTCRgd_NM295_Dn14$TCRD.CDR3.nt),
                                   type = "dna",
                                   method = "ClustalW")
msaPrettyPrint(macTCRgd_NM295_Dn14_TCRG_nt, 
               output = "pdf",
               showNames = "none", 
               showLogo = "top",
               logoColors = "rasmol",
               shadingMode = "similar",
               showLegend = FALSE,
               askForOverwrite = FALSE,
               file = "macTCRgd_NM295_Dn14_TCRG_nt.pdf",
               alFile = "macTCRgd_NM295_Dn14_TCRG_nt.fasta",
               verbose = FALSE)

# Reset Working Directory
setwd("..")
