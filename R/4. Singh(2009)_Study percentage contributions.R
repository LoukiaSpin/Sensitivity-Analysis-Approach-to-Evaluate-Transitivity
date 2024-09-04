#*******************************************************************************
#*
#*
#*           Motivating example: Singh et al., 2009 (PMID: 19821440)                                                                
#*        <Transitivity evaluation using Study percentage contributions>                                                                                                                                      
#*
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("rnmamod", "tracenma",  "reshape2", "ggplot2", "ggrepel")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions
source("./30_Analysis/rnmamod functions/convert_long_to_wide_function.R")
source("./30_Analysis/contribution functions/design.matrices_function.R")
source("./30_Analysis/contribution functions/study.perc.contrib_function.R")
source("./30_Analysis/contribution functions/covar.contribution.plot_function.R")


## Load data 
# From 'tracenma'
data_set <- 
  get.dataset(pmid = 19821440, show.index = FALSE, show.type = FALSE)$Dataset

# The NMA outcome data 
load("./31_Datasets/2009_19821440_Singh.RData")

# Convert long format (nmadb) into wide format (rnmamod)
data_nma <- convert_long_to_wide(data_nma0)

# Remove rows with zero events and sample size
data_nma_fin <- subset(data_nma, n.1 > 0)
rownames(data_nma_fin) <- data_set$trial

# Treatment names
treat_names <- c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT")


## Gower dissimilarity
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


## Dissimilarities
study_diss <- 
  comp_clustering(input = dataset_new, 
                  threshold = 0.13,
                  get_plots = TRUE)


## Define the similarity weights
# Root mean square
rms_wgt <- weight_defined(diss_res = list(study_diss, "rms"))


## Get contrast-based results for each study
contrast_res <- 
  netmeta::pairwise(list(t.1, t.2),
                    list(r.1, r.2), 
                    list(n.1, n.2),
                    data = data_nma_fin,
                    sm = "OR")


## Get study contributions to network meta-regression results
contrib_res <- 
  study_perc_contrib(study.name = contrast_res$studlab,
                     base.t = contrast_res$treat1, 
                     exp.t = contrast_res$treat2, 
                     ref.t = 1,
                     obs.se = contrast_res$seTE,
                     covar = rms_wgt$weights,
                     covar.assum = "no",
                     model = "RE",
                     tau = 0.58) 


## Covariate-contribution plot
# Basic parameters
tiff("./30_Analysis/Figure 4_Study contributions_basic pars.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
covar_contribution_plot(contr_res = contrib_res, 
                        comparisons = "basic",
                        drug_names = treat_names,
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5)
dev.off()

# Functional parameters
tiff("./30_Analysis/Figure S1_Study contributions_functional pars.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
covar_contribution_plot(contr_res = contrib_res, 
                        comparisons = "functional",
                        drug_names = treat_names,
                        upper.limit = 60,
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        seq_by = 0.11)
dev.off()


## Summary of contributions - Basic parameters
# Restrict to studies with similarities above 60%
contrib_res$perc_contribute[c(1, 2, 6), 4] <- 0

# Summary of non-zero contributions
summary(contrib_res$perc_contribute[, 4:9][contrib_res$perc_contribute[, 4:9] > 0])


## Summary of contributions - Functional parameters
summary(contrib_res$perc_contribute[, -c(1:9, 25)][contrib_res$perc_contribute[, -c(1:9, 25)] > 0])
