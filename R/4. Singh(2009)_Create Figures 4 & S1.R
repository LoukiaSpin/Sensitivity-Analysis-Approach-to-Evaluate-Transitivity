#*******************************************************************************
#*
#*
#*                 Creating Figure 4 & Supplementary Figure S1                                                                                                                                                                                                                           
#*         <Study percentage contributions against study similarities>                                                                                                                                                                                                  
#* 
#* Author: Loukia M. Spineli
#* Date: September 2024
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("rnmamod", "tracenma", "netmeta")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions
source("./R/convert_long_to_wide_function.R")


## Load data 
# From 'tracenma'
data_set <- 
  get.dataset(pmid = 19821440, show.index = FALSE, show.type = FALSE)$Dataset

# The NMA outcome data 
load("./data/19821440_Singh 2009_Outcome data.RData")


## Prepare outcome dataset
# Convert long format (nmadb) into wide format (rnmamod)
data_nma <- convert_long_to_wide(data_nma0)

# Remove rows with zero events and sample size
data_nma_fin <- subset(data_nma, n.1 > 0)
rownames(data_nma_fin) <- data_set$trial

# Treatment names
treat_names <- c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT")


## Prepare dataset with characteristics
# STEP 1: Turn into data.frame
dataset_new0 <- as.data.frame(data_set)

# STEP 2: Remove the columns with treatment names
dataset_new <- dataset_new0[, -c(4,5)]

# STEP 3: Turn the treatment ID columns from 'double' to 'character'
dataset_new[, 2:3] <- lapply(dataset_new[, 2:3], as.character)

# STEP 4: Turn the ‘character’ characteristics into ‘integer’
dataset_new[, -c(1:3)] <- 
  lapply(dataset_new[, -c(1:3)], 
         function(x) if (typeof(x) == "character") as.factor(x))


## Gower's dissimilarity for all study pairs
# ?comp_clustering
study_diss <- 
  comp_clustering(input = dataset_new, 
                  threshold = 0.13,
                  get_plots = TRUE)


## Get the weights under approaches (a) and (b) (?weight_defined)
# Between-comparison similarity
rms_wgt <- weight_defined(diss_res = list(study_diss, "fixed"))

# Relative similarity index
relative_wgt <- weight_defined(diss_res = list(study_diss, "index"))


## Get contrast-based results for each study
contrast_res <- pairwise(list(t.1, t.2),
                         list(r.1, r.2), 
                         list(n.1, n.2),
                         data = data_nma_fin,
                         sm = "OR")


## Get study contributions using the 'Between-comparison similarities'
# ?study_perc_contrib
contrib_rms <- 
  study_perc_contrib(study_name = contrast_res$studlab,
                     base_t = contrast_res$treat1, 
                     exp_t = contrast_res$treat2, 
                     ref_t = 1,
                     obs_se = contrast_res$seTE,
                     covar = rms_wgt$weights,
                     covar_assum = "no",
                     model = "RE",
                     tau = 0.58) 


## Covariate-contribution plot with 'Between-comparison similarities' (?covar_contribution_plot)
# Basic parameters
tiff("./Figures/Figure 4.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_rms, 
                        comparisons = "basic",
                        drug_names = treat_names,
                        name_x_axis = "Between-comparison similarity",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        seq_by = 0.05)
dev.off()

# Functional parameters
tiff("./Figures/Figure S1.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_rms, 
                        comparisons = "functional",
                        drug_names = treat_names,
                        upper_limit = 60,
                        name_x_axis = "Between-comparison similarity",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        seq_by = 0.10)
dev.off()


## Summary of contributions - Basic parameters
# Restrict to studies with similarities above 60%
contrib_rms$perc_contribute[c(1, 2, 6), 4] <- 0

# Summary of non-zero contributions
summary(contrib_rms$perc_contribute[, 4:9][contrib_rms$perc_contribute[, 4:9] > 0])


## Summary of contributions - Functional parameters
summary(contrib_rms$perc_contribute[, -c(1:9, 25)][contrib_rms$perc_contribute[, -c(1:9, 25)] > 0])


## Get study contributions using the 'Relative similarity index'
# ?study_perc_contrib
contrib_relative <- 
  study_perc_contrib(study_name = contrast_res$studlab,
                     base_t = contrast_res$treat1, 
                     exp_t = contrast_res$treat2, 
                     ref_t = 1,
                     obs_se = contrast_res$seTE,
                     covar = relative_wgt$weights,
                     covar_assum = "no",
                     model = "RE",
                     tau = 0.58) 


## Covariate-contribution plot with 'Relative similarity index' (?covar_contribution_plot)
# Basic parameters
tiff("./Figures/Figure S2.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_relative, 
                        comparisons = "basic",
                        drug_names = treat_names,
                        name_x_axis = "Relative similarity index",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        seq_by = 0.035)
dev.off()

# Functional parameters
tiff("./Figures/Figure S3.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_relative, 
                        comparisons = "functional",
                        drug_names = treat_names,
                        upper_limit = 60,
                        name_x_axis = "Relative similarity index",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        seq_by = 0.035)
dev.off()
