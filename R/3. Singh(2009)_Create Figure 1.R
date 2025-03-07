#*******************************************************************************
#*
#*
#*                               Creating Figure 1                                                                                                                                                                                                                                       
#*       <Between- and within-comparison study-specific dissimilarities>                                                                                                                                                                                                             
#* 
#* Author: Loukia M. Spineli
#* Date: December 2024
#*       
#*******************************************************************************


## Load the latest development version
#devtools::install_github("LoukiaSpin/rnmamod", force = TRUE)


## Load libraries
list.of.packages <- c("rnmamod", "tracenma")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load data 
# From 'tracenma'
data_set <- 
  get.dataset(pmid = 19821440, show.index = FALSE, show.type = FALSE)$Dataset

# Make the study names 'prettier'
data_set$trial <- sub("\\s+[^ ]+$", "", data_set$trial)


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
                  drug_names = c("PBO", "ABA", "ADA", "ANA", "ETA", "INF", "RIT"),
                  threshold = 0.13,
                  get_plots = TRUE)


## Create the rainbow plot 
# ?rainbow_similarities
tiff("./Figures/Figure 1.tiff", 
     height = 30, 
     width = 50, 
     units = "cm", 
     compression = "lzw", 
     res = 300)
plot_study_dissimilarities(results = study_diss,
                           axis_title_size = 14,
                           axis_text_size = 14,
                           label_size = 4.3) 
dev.off()
