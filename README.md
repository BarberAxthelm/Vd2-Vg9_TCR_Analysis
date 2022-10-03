# Vd2-Vg9_TCR_Analysis

## Notes

- MSA script cannot have spaces in the working directory filepath. If you have spaces in your filepath, please update the working directory before running.

## Data Prep


## Script Descriptions

The brief script descriptions are below. TCR_Data_Preparation_MiXCR.R script must be run before any other script, as it formats the table to be used with referenced column names.

### TCR_Alluvial_Plot.R

This script is used to create normalized Alluvial Plots for all animals. 

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - 12 Alluvial plots, combo of:
        - All animals ("NM251", "NM295", "NM89", "NM11")
        - TCRG, TCRD, Combined TCRD/TCRG
    - 4 CSV files, each containing counts of each sequence for a specific animal

### TCR_CDR3_Spectratyping.R

This script outputs a spectratyping bar plot showing the length distribution  of the TCRD/TCRG CDR3s (aa) and saves statistics for CDR3 lengths in a CSV.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - 12 Spectratyping bar plots, combo of:
        - All animals ("NM251", "NM295", "NM89", "NM11")
        - TCRG, TCRD, Combined TCRD/TCRG
    - 2 CSV files, showing spectratyping stats for TCRD and TCRG chais

### TCR_Circos_Plots.R

This script creates circos plots for a single animal, both normalized and non-normalized for TCRD, TCRG, and TCRD-TCRG chains. It must be manually adjusted to run for the other animals.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - 3 circos plots (RCRD, TCRG, TCRD-TCRG)
    - 3 normalized circos plots (RCRD, TCRG, TCRD-TCRG)

### TCR_Data_Preparation.R

This script is used to prepare the data for use by the other analysis scripts. 

It imports:

    - a CSV file containing data formatted appropriately post-MiXCR analysis
    - a CSV containing the MiXCR list of functional clones

It exports two tables as CSV files: 

    - MiXCR Data, filtered to include only paired functional sequences and 
                    related sample infomation and flow cytometery data.
    - Counts of functional paired sequences per sample.

### TCR_Diversity.R

This runs diversity analytics on TCRD, TCRG, and combined chains. 

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - a CSV containing diversity statistics.
    - a CSV containing J-chain usages for individual animals and timepoints.

### TCR_MSA_Sequence_Logos.R

This script is used to generate the multiple sequence alignment (msa) and sequence logo generations for each animal, timepoint, and chain (TCRD/TCRG), both as a nucleotide sequence and amino acid sequence.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - a fasta and pdf file for every combination of animal, timepoint, and chain.

### TCR_Randomization_Statistics.R

This script creates randomized lists based on the data, to run statistics on. It can be used to compare to the output of TCR_Diversity.R. Manual edits need to be made to run for each animal. It is setup to run for one animal at a time.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - a CSV containing diversity statistics for a particular animal.

### TCR_Upset_Plot.R

This scripts creates 3 upset plots (for all animals): TCRD, TCRG, TCRD-TCRG.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - 3 Upset Plots: TCRD, TCRG, TCRD-TCRG.

### TCRG_CDR3_Shared_Clonotypes.R

This script identifies clonoltypes  shared between animals and timepoints for TCRG.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - CSV with CDR3 clonoltype and the count for each animal_timepoint sample and a column indicating how many animals/samples share that clonoltype.

For example:

CDR3             | NM11_Dn14 | NM11_D4 | NM89_Dn14 | ... | shared
CALWEVQQFGRKVKLF | 11        | 10      | 4         | ... | 6