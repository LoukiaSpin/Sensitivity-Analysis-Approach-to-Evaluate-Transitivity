# Quantifying study dissimilarities and exploring study contributions aid in streamlining the transitivity evaluation in network meta-analysis 

## Description of the repository

The repository offers the typical structure of separate folders for R (scripts to replicate the main and supplementary Tables and Figures) and output (file _Figures_). 
* Each R script indicates which Table(s) and Figure(s) it produces.
* The dataset for the star network can be found in the R folder under 19821440_Singh 2009_Outcome data.RData. This folder contains extra .RData files that refer to 1) the network meta-analysis results under the standard NMA and the enriching-through-weighting approach using the three proposed weighting approaches, and 2) the results of the multiple network meta-regression analyses to produce the Supplementary Table S1.

After downloading/cloning the repo, the user can use the .Rproj file to source all code.

## Output 

Prerequisite R packages: 
[rnmamod](https://CRAN.R-project.org/package=rnmamod),
[tracenma](https://CRAN.R-project.org/package=tracenma),
[gtsummary](https://CRAN.R-project.org/package=gtsummary),
[plyr](https://CRAN.R-project.org/package=plyr),
[netmeta](https://CRAN.R-project.org/package=netmeta)