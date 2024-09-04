#*******************************************************************************
#*
#*
#*           Motivating example: Singh et al., 2009 (PMID: 19821440)                                                                
#*              <Transitivity evaluation using Weighting>                                                                                                                              
#*
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("rnmamod", "tracenma")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions
source("./30_Analysis/rnmamod functions/convert_long_to_wide_function.R")


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

# Uniform distribution
uniform_wgt <- weight_defined(diss_res = list(study_diss, "uniform"))


## Illustrate study similarities
tiff("./30_Analysis/Figure 1_Rainbow plot.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
rainbow_similarities(results = study_diss,
                     axis_title_size = 14,
                     axis_text_size = 14,
                     label_size = 4.3) 
dev.off()


## Run the primary model
primary <- 
  run_model(data = data_nma_fin,
            measure = "OR",
            model = "RE",
            heter_prior = list("halfnormal", 0, 1),
            D = 1,
            ref = 1,
            n_chains = 3,
            n_iter = 300000,
            n_burnin = 200000,
            n_thin = 5,
            adjust_wgt = NULL)
#save(primary, file = "./30_Analysis/Singh 2009_R scripts/19821440 Singh_primary.RData")
load("./30_Analysis/Singh 2009_Primary NMA results.RData")
mcmc_diagnostics(net = primary,
                 par = c("EM", "tau"))


## Run the adjusted model with *fixed weights* (rms, root mean square)
rms_adjust <- 
  run_model(data = data_nma_fin,
            measure = "OR",
            model = "RE",
            heter_prior = list("halfnormal", 0, 1),
            D = 1,
            ref = 1,
            n_chains = 3,
            n_iter = 300000,
            n_burnin = 200000,
            n_thin = 5,
            adjust_wgt = rms_wgt$weights) 
#save(rms_adjust, file = "./30_Analysis/Singh 2009_Fixed weights results.RData")
load("./30_Analysis/Singh 2009_Fixed weights results.RData")
mcmc_diagnostics(net = rms_adjust,
                 par = c("EM", "tau"))


## Run the adjusted model with 'uniform prior distribution'
unif_adjust <- 
  run_model(data = data_nma_fin,
            measure = "OR",
            model = "RE",
            heter_prior = list("halfnormal", 0, 1),
            D = 1,
            ref = 1,
            n_chains = 3,
            n_iter = 300000,
            n_burnin = 200000,
            n_thin = 5,
            adjust_wgt = uniform_wgt$weights)
#save(unif_adjust, file = "./30_Analysis/Singh 2009_Uniform prior results.RData")
load("./30_Analysis/Singh 2009_Uniform prior results.RData")
mcmc_diagnostics(net = unif_adjust,
                 par = c("EM", "tau"))


## Forest plot with all model results
# Treatment effects of basic parameters
tiff("./30_Analysis/Figure 2_Forest plot Effects.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
forestplot_juxtapose(results = list(primary, rms_adjust, unif_adjust), 
                     compar = "PBO", 
                     name = c("Standard NMA", "Fixed weights", "Uniform prior distribution"), 
                     drug_names = treat_names,
                     axis_title_size = 16,
                     axis_text_size = 16,
                     caption_text_size = 12,
                     label_size = 5.5)$forest_plots
dev.off()

# Between-study standard deviation
tiff("./30_Analysis/Figure 3_Forest plot tau.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
forestplot_juxtapose(results = list(primary, rms_adjust, unif_adjust), 
                     compar = "PBO", 
                     name = c("Standard NMA", "Fixed weights", "Uniform prior distribution"), 
                     drug_names = treat_names,
                     axis_title_size = 16,
                     axis_text_size = 16,
                     label_size = 5.5)$tau_plot
dev.off()
