# Middleton_miRNA

R code used for data manipulation, statistical analyses and figure generation for the manuscript "Rhizospheric miRNAs affect the plant microbiota"

If you intend to re-run the analyses, please download the entire folder structure (click on "Code", then "Download ZIP" )

The raw data needed to run the scripts can be found on Zenodo: <https://zenodo.org/doi/10.5281/zenodo.11105307>. It should be copied to the /data/raw folder.

The scripts should be run in order (01 to 06).

**01-LoadPackages.R**: Download, install and load necessary packages

**02-LoadRawDataNormalise.R**: Load raw data and arrange.

**03-miRNAAbundance.R**: Analysis and figure generation for plant, rhizospheric and bacterial miRNAs.

**04-Transcriptomics.R**: RNA-seq analyses for isolated bacteria (Variovorax and Bacillus). qPCR analyses of selected genes of Variovorax in the rhizosphere of Arabidopsis treated with a miPEP. Figure generation.

**05-Community.R**: Look at the effect of miRNAs on rhizospheric microbial communities: 1) mutant Arabidopsis, 2) Arabidopsis treated with miPEP, 3) simplified soil community confronted with synthetic miRNAs. Figure generation.

**06-FlowCytometry.R**: Analyse data and create graphs from flow cytometry.

**07-Figures.R**: Create multi-panel figures for publication.
