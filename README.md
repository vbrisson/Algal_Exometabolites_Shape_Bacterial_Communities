# Algal_Exometabolites_Shape_Bacterial_Communities
Analysis code for "Dynamic Phaeodactylum tricornutum Exometabolites Shape Surrounding Bacterial Communities"

## raw_sequence_analysis
  This folder contains R code for analyzing raw sequencing data with DADA2
  The sequencing data are available from the Sequence Read Archive under BioProject number PRJNA803592:
     https://www.ncbi.nlm.nih.gov/bioproject/PRJNA803592/
  The taxonomy reference databases used in this analysis were from:
     https://benjjneb.github.io/dada2/training.html

## ASV_analyses
  This folder contains R code for analyzing ASV data using the phyloseq and ANCOM-BC
  The analyses include data filtering, adjustment for data compositionality, ANOVA analysis of alpha (Shannon) 
  diversity, Principal Coordinate Analysis (PCoA) of ASV community composition, PERMANOVA of community composition,
  and differential abundance analyses.  

## plotting_and_further_analyses
  This folder contains Jupyter notebooks for generating plots, and for additional analysis, of the results from the
  ASV analyses above, and of the metabolomics results from the Metabolite Atlas (https://github.com/biorack/metatlas)
  of LC-MS/MS.
