#*******************************************************************************
#*
#*
#*                           Creating Figures 1 to 3                                                                                                                                                                                                          
#*            <Rainbow plot and Enriching-through-weighting approach>                                                                                                                                                                                         
#* 
#* Author: Loukia M. Spineli
#* Date: September 2024
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("rnmamod", "tracenma")
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


## Define the similarity weights (?weight_defined)
# Fixed weights
fixed_wgt <- weight_defined(diss_res = list(study_diss, "fixed"))

# Uniform distribution
uniform_wgt <- weight_defined(diss_res = list(study_diss, "uniform"))

# Relative index
index_wgt <- weight_defined(diss_res = list(study_diss, "index"))


## Create the rainbow plot 
# ?rainbow_similarities
tiff("./Figures/Figure 1.tiff", 
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


## Run the 'standard' NMA (random-effects model)
# ?run_model
set.seed(05092024)
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
            n_thin = 5)
#save(primary, file = "./data/Results standard NMA.RData")
#load("./data/Results standard NMA.RData")

# MCMC diagnostics (?mcmc_diagnostics)
mcmc_diagnostics(net = primary,
                 par = c("EM", "tau"))


## Run the enriching-through-weighting model with *between-comparison similarities* 
set.seed(05092024)
fixed_adjust <- 
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
            adjust_wgt = fixed_wgt$weights) 
#save(fixed_adjust, file = "./data/Results fixed weights.RData")
#load("./data/Results fixed weights.RData")

# MCMC diagnostics 
mcmc_diagnostics(net = fixed_adjust,
                 par = c("EM", "tau"))


## Run the enriching-through-weighting model with 'uniform prior distribution'
set.seed(05092024)
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
#save(unif_adjust, file = "./data/Results uniform prior.RData")
#load("./data/Results uniform prior.RData")

# MCMC diagnostics 
mcmc_diagnostics(net = unif_adjust,
                 par = c("EM", "tau"))


## Run the enriching-through-weighting model with 'relative similarity index'
set.seed(05092024)
rel_index <- 
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
            adjust_wgt = index_wgt$weights)
#save(rel_index, file = "./data/Results relative index.RData")
#load("./data/Results relative index.RData")

# MCMC diagnostics 
mcmc_diagnostics(net = rel_index,
                 par = c("EM", "tau"))


## Forest plot with all model results
# Treatment effects of basic parameters (?forestplot_juxtapose)
tiff("./Figures/Figure 2.tiff", 
     height = 30, 
     width = 55, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
forestplot_juxtapose(results = list(primary, fixed_adjust, unif_adjust, rel_index), 
                     compar = "PBO", 
                     name = c("Standard NMA", "Between-comparison similarity", "Uniform prior distribution", "Relative similarity index"), 
                     drug_names = treat_names,
                     axis_title_size = 16,
                     axis_text_size = 16,
                     caption_text_size = 12,
                     label_size = 5.5)$forest_plots
dev.off()

# Between-study standard deviation
tiff("./Figures/Figure 3.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 600)
forestplot_juxtapose(results = list(primary, fixed_adjust, unif_adjust, rel_index), 
                     compar = "PBO", 
                     name = c("Standard NMA", "Between-comparison similarity", "Uniform prior distribution", "Relative similarity index"), 
                     drug_names = treat_names,
                     axis_title_size = 16,
                     axis_text_size = 16,
                     label_size = 5.5)$tau_plot
dev.off()
