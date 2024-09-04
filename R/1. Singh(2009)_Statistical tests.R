#*******************************************************************************
#*
#*
#*                               Creating Table 1                                                                                                                                                    
#*               <Characteristics summary and Statistical tests>                                                                                                                                                            
#* 
#*       
#*******************************************************************************


## Load libraries
list.of.packages <- c("tracenma", "gtsummary")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)


## Load data 
# From 'tracenma'
data_set <- 
  as.data.frame(get.dataset(pmid = 19821440, show.index = FALSE, show.type = FALSE)$Dataset)

# Add comparisons
data_set$comparison <- paste(data_set$treat1, data_set$treat2)

# Re-arrange the levels
data_set$MTX.use <- factor(data_set$MTX.use, levels = c("yes", "no"))
data_set$`anti-TNF` <- factor(data_set$`anti-TNF`, levels = c("yes", "no"))
data_set$prior.drugs.failed <- factor(data_set$prior.drugs.failed, 
                                      levels = c("biologic", "DMARDs", "both", "none"))
data_set$prior.failure.TNF <- factor(data_set$prior.failure.TNF, levels = c("yes", "no"))
data_set$comb.biologic <- factor(data_set$comb.biologic, levels = c("yes", "no"))
data_set$biologic.naive <- factor(data_set$biologic.naive, levels = c("yes", "no"))


## Descriptive statistics
subset(data_set, select = -c(trial, treat1, treat2, arm1, arm2)) %>% 
  tbl_summary(
    by = comparison,
    statistic = list(all_continuous() ~ "{median} ({min}, {max})",
                     all_categorical() ~ "{n} ({p}%)"),
    type = list(sample.size ~ "continuous",
                duration ~ "continuous",
                disease.duration ~ "continuous",
                MTX.use ~ "categorical",
                RA.duration ~ "categorical",
                `anti-TNF` ~ "categorical",
                prior.drugs.failed ~ "categorical",
                prior.failure.TNF ~ "categorical",
                comb.biologic ~ "categorical",
                biologic.naive ~ "categorical"),
    digits = list(all_categorical() ~ c(0, 1))) %>%
  italicize_levels() 


## Perform multiple statistical tests
# Restrict to quantitative characteristics
quant <- data_set[, c("comparison", "sample.size", "duration", "disease.duration")]

# F-test of one-way ANOVA
res_quant <- 
  lapply(2:dim(quant)[2], function(x) oneway.test(quant[, x] ~ quant[, 1], var.equal = FALSE)$p.value)
names(res_quant) <- colnames(quant)[-1]

# Restrict to quantitative characteristics
qual <- subset(data_set, select = -c(trial, treat1, treat2, arm1, arm2, sample.size, duration, disease.duration))

# Chi-squared test
res_qual <- 
  lapply(1:(dim(qual)[2] - 1), function(x) chisq.test(table(qual[, x], qual[, 8]))$p.value)
names(res_qual) <- colnames(qual)[-8]
