# Vd2-Vg9_TCR_Analysis

## Notes

- MSA script cannot have spaces in the working directory filepath. If you have spaces in your filepath, please update the working directory before running.

## Data Prep


## Script Descriptions

The brief script descriptions are below. TCR_Data_Preparation_MiXCR.R script must be run before any other script, as it formats the table to be used with referenced column names.

### TCR_Data_Preparation_MiXCR.R

This script is used to prepare the data for use by the other analysis scripts. 

It imports:

    - a CSV file containing data formatted appropriately post-MiXCR analysis
    = a CSV containing the MiXCR list of functional clones

It exports two tables as CSV files: 

    - MiXCR Data, filtered to include only paired functional sequences and 
                    related sample infomation and flow cytometery data.
    - Counts of functional paired sequences per sample.


### TCR_MSA_Sequence_Logos_Pipeline.R

This script is used to generate the multiple sequence alignment (msa) and sequence logo generations for each animal, timepoint, and chain (TCRD/TCRG), both as a nucleotide sequence and amino acid sequence.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - a fasta and pdf file for every combination of animal, timepoint, and chain.


### TCRG_CDR3_Shared_Clonotypes.R

This script identifies clonoltypes  shared between animals and timepoints for TCRG.

It imports:

    - CSV MiXCR Data file formated by TCR_Data_Preparation_MiXCR.R

It exports: 

    - CSV with CDR3 clonoltype and the count for each animal_timepoint sample and a column indicating how many animals/samples share that clonoltype.

For example:

CDR3             | NM11_Dn14 | NM11_D4 | NM89_Dn14 | ... | shared
CALWEVQQFGRKVKLF | 11        | 10      | 4         | ... | 6