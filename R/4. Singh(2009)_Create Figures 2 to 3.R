#*******************************************************************************
#*
#*
#*                           Creating Figures 2 and 3                                                                                                                                                                                                                                                                                                        
#*         <Study percentage contributions against study similarities>                                                                                                                                                                                                  
#* 
#* Author: Loukia M. Spineli
#* Date: December 2024
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

# Make the study names 'prettier'
data_set$trial <- sub("\\s+[^ ]+$", "", data_set$trial)

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
dataset_new <- dataset_new0[, -c(4, 5)]

# STEP 3: Turn the treatment ID columns from 'double' to 'character'
dataset_new[, 2:3] <- lapply(dataset_new[, 2:3], as.character)

# STEP 4: Turn the ‘character’ characteristics into ‘integer’
dataset_new[, -c(1:3)] <- 
  lapply(dataset_new[, -c(1:3)], 
         function(x) if (typeof(x) == "character") as.factor(x) else x)


## Gower's dissimilarity for all study pairs
# ?comp_clustering
study_diss <- 
  comp_clustering(input = dataset_new, 
                  drug_names = treat_names,
                  threshold = 0.13, 
                  get_plots = TRUE)


## Get the between-comparisons dissimilarities
# ?plot_study_dissimilarities
between_diss <- plot_study_dissimilarities(results = study_diss)$diss_values[, "between_multiarm"] 


## Get contrast-based results for each study
contrast_res <- pairwise(studlab = id,
                         treat = list(t.1, t.2),
                         event = list(r.1, r.2), 
                         n = list(n.1, n.2),
                         data = data_nma_fin,
                         sm = "OR")


## Run RE-NMA with consistency
# ?run_model
primary <- 
  run_model(data = data_nma_fin,
            measure = "OR",
            model = "RE",
            heter_prior = list("halfnormal", 0, 1),
            D = 1, # positive outcome
            ref = 1,
            n_chains = 3,
            n_iter = 300000,
            n_burnin = 200000,
            n_thin = 5)


## Get study contributions using the 'Between-comparison similarities'
# ?study_perc_contrib
contrib_rms <- 
  study_perc_contrib(study_name = contrast_res$studlab,
                     base_t = contrast_res$treat1, 
                     exp_t = contrast_res$treat2, 
                     ref_t = 1,
                     obs_se = contrast_res$seTE,
                     covar = between_diss,
                     covar_assum = "no",
                     model = "RE",
                     tau = primary$tau[5])


## Covariate-contribution plot with 'Between-comparison similarities' (?covar_contribution_plot)
# Basic parameters
tiff("./Figures/Figure 2.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_rms, 
                        comparisons = "basic",
                        drug_names = treat_names,
                        name_x_axis = "Between-comparisons dissimilarity",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        percentage = TRUE)
dev.off()

# Functional parameters
tiff("./Figures/Figure 3.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
covar_contribution_plot(contr_res = contrib_rms, 
                        comparisons = "functional",
                        drug_names = treat_names,
                        upper_limit = 60,
                        name_x_axis = "Between-comparisons dissimilarity",
                        axis_title_size = 16,
                        axis_text_size = 16,
                        strip_text_size = 16,
                        subtitle_size = 16,
                        label_size = 5,
                        percentage = TRUE)
dev.off()
