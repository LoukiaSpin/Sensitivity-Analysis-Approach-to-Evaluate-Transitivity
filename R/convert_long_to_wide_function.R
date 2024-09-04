#*******************************************************************************
#*
#*                Function to turn a gemtc dataset (long format)                                
#*                     into an rnmamod dataset (wide format)                                                                                                                                                                    
#*
#* Author: Loukia M. Spineli
#* Date: July 2024
#* 
#*******************************************************************************

convert_long_to_wide <- function (dataset) {

 
  # Calculate the number of arms per study
  na <- table(dataset[, 1])
  
  # Add an indicator variable to facilitate reshaping
  dataset$numbers <- unlist(lapply(na, function(x) x:1))
  
  # Get the desire data format for 'rnmamod'
  dataset_fin0 <- reshape(dataset, 
                          idvar = "id", 
                          timevar = "numbers", 
                          direction = "wide")
  
  # Sort columns by column name
  dataset_fin <- dataset_fin0[ , order(names(dataset_fin0))]
  
  
  # Rename the columns as required by rnmamod
  #colnames(dataset_fin) <-
  #  if (any(!startsWith(colnames(dataset_fin), "r.")) == TRUE) {  # Continuous outcome
  #
  #    gsub("treatment.", "t",
  #         gsub("sampleSize.", "n",
  #              gsub("mean.", "y",
  #                   gsub("std.dev.", "sd", colnames(dataset_fin)))))
  #} else {
  #  
  #  gsub("treatment.", "t",
  #       gsub("sampleSize.", "n",
  #            gsub("responders.", "r", colnames(dataset_fin))))
  #}
  
  # Turn treatment names into numbers
  #dataset_fin[, startsWith(colnames(dataset_fin), "t")] <- 
  #  matrix(unlist(lapply(dataset_fin[, startsWith(colnames(dataset_fin), "t")], 
  #                       function(x) as.numeric(x))), 
  #         nrow = dim(dataset_fin)[1], ncol = max(na), byrow = FALSE)
  
  return(dataset_fin)
}
