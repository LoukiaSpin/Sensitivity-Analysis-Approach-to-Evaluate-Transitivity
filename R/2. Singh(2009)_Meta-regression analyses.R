#*******************************************************************************
#*
#*
#*           Motivating example: Singh et al., 2009 (PMID: 19821440)                                                                
#*              <Transitivity evaluation using Meta-regression>                                                                                                                              
#*
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("tracenma", "rnmamod", "plyr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load functions
source("./30_Analysis/rnmamod functions/convert_long_to_wide_function.R")


## Load data 
# From 'tracenma'
data_set <- 
  as.data.frame(get.dataset(pmid = 19821440, show.index = FALSE, show.type = FALSE)$Dataset)

# Copy to another name to 'revalue' categorical variables
data_set_new <- data_set[, c("sample.size", "duration", "disease.duration", "MTX.use", "anti-TNF", 
                             "prior.failure.TNF", "comb.biologic", "biologic.naive")]

# Transform categorical into numerical
data_set_new$MTX.use <- as.numeric(revalue(data_set_new$MTX.use, c("yes" = "1", "no" = "0")))
data_set_new$`anti-TNF` <- as.numeric(revalue(data_set_new$`anti-TNF`, c("yes" = "1", "no" = "0")))
data_set_new$prior.failure.TNF <- as.numeric(revalue(data_set_new$prior.failure.TNF, c("yes" = "1", "no" = "0")))
data_set_new$comb.biologic <- as.numeric(revalue(data_set_new$comb.biologic, c("yes" = "1", "no" = "0")))
data_set_new$biologic.naive <- as.numeric(revalue(data_set_new$biologic.naive, c("yes" = "1", "no" = "0")))


## Load outcome data for NMA
# NMA primary outcome
load("./31_Datasets/2009_19821440_Singh.RData")

# Convert long format (nmadb) into wide format (rnmamod)
data_nma <- convert_long_to_wide(data_nma0)

# Remove rows with zero events and sample size
data_nma_fin <- subset(data_nma, n.1 > 0)
rownames(data_nma_fin) <- data_set$trial

# Treatment names
treat_names <- c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT")


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
            n_thin = 5)
#save(primary, file = "./30_Analysis/Singh 2009_Primary NMA results.RData")
load("./30_Analysis/Singh 2009_Primary NMA results.RData")

# MCMC diagnostics
mcmc_diagnostics(net = primary,
                 par = c("EM", "tau"))


## Run network meta-regression with common regression coefficients
# Set the combinations
scenarios <- expand.grid(x = 1:dim(data_set_new)[2], y = c("independent", "exchangeable", "common"))

# List of results for all combinations
list_res <- 
  lapply(1:dim(scenarios)[1], 
         function(x) run_metareg(full = primary,
                                 # If the covariate is metric, transform into the logarithmic scale
                                 covariate = if (length(unique(data_set_new[, scenarios[x, 1]])) == 2) data_set_new[, scenarios[x, 1]] else log(data_set_new[, scenarios[x, 1]]),
                                 covar_assumption = scenarios[x, 2],
                                 # If the covariate is metric, consider the mean of the logarithm; otherwise, the presence of the characteristic (yes == 1)
                                 cov_value = if (length(unique(data_set_new[, scenarios[x, 1]])) == 2) 1 else mean(log(data_set_new[, scenarios[x, 1]])),
                                 n_chains = 3,
                                 n_iter = 300000,
                                 n_burnin = 200000,
                                 n_thin = 5)$beta_all)

# Name the elements of the list by the combinations
names(list_res) <- apply(expand.grid(x = colnames(data_set_new), 
                                     y = c("independent", "exchangeable", "common")), 
                         1, paste, collapse = " ")
#save(list_res, file = "./30_Analysis/Singh 2009_NMR results.RData")
load("./30_Analysis/Singh 2009_NMR results.RData")


## Keep median and bounds after exponentiation
#lapply(list_res, function(x) format(exp(x[1:6, c(5, 3, 7)]), scientific = FALSE, digits = 2))
lapply(list_res, function(x) ifelse(round(exp(x[1:6, c(5, 3, 7)]), 3) > 5, 
                                    format(exp(x[1:6, c(5, 3, 7)]), scientific = TRUE, digits = 3), 
                                    round(exp(x[1:6, c(5, 3, 7)]), 3)))
