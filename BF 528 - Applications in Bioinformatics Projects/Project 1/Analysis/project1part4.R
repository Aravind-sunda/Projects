# Caroline Sullivan
# BF528
# Project 1 Part 4
# Role: Analyst 

library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ----------------------- read in the data -------------------------------------
# reading in the normalized microarray from part 3, saved in a .csv file
# data is a sample x gene tibble

# column name for probeset is empty so we need to create new set of column names
columns <- colnames(read.csv("corrected_for_batch_effects.csv")) 
columns <- c('probesets', columns[-1])

# I want to remove everything but the initial piece with the patient ID
columns <- gsub('_.*','',columns)

# now we can read data, providing our own header and skipping the original
data <- readr::read_csv("corrected_for_batch_effects.csv",
                          col_names = columns, skip=1,
                          show_col_types = FALSE)

# ----------------------- 4.1: Filter 1 ----------------------------------------
# Filter for genes expressed in at least 20% of samples (i.e. for each gene, 
# at least 20% of the gene-expression values must be > log2(15)).
filter_by_expression <- function(data){
data_filtered_by_expression <- data %>%
  dplyr::filter(rowSums(across(where(is.numeric)) > log2(15)) >= (ncol(data)-1)*0.2)
  
return(data_filtered_by_expression)
  
}

filter_expression <- filter_by_expression(data)

# ----------------------- 4.2: Filter 2 ----------------------------------------
# Filter for genes that have a variance significantly different from the median 
# variance of all probe sets using a threshold of p < 0.01
# could be different in either direction, so we want a two-tailed test
filter_by_variance <- function(data){
  # before starting, get num sampls and deg of freedom
  n_samples <- ncol(data) - 1
  deg_freedom <- n_samples - 1
  
  # get the variance for each probeset
  vars <- data %>%
    dplyr::select(-probesets) %>%
    apply(1, var) %>%
    as_tibble(rownames='probesets') %>%
    dplyr::rename('variance'=value) %>%
    dplyr::mutate(probesets=data$probesets)
  
  # get the median variance
  median_var <- median(vars$variance)
  
  # add variances to data
  data_filtered_by_variance <- dplyr::full_join(data, vars, by='probesets')
  
  # calculate the test statistic
  # formula: ((n-1)×Var(P) / Varmed) where n in number of samples
  data_filtered_by_variance <- dplyr::mutate(data_filtered_by_variance,
    test_statistic=(deg_freedom)*variance/median_var)
  
  # set alpha (threshold for significance)
  alpha_val <- 0.01
  
  # get critical values with chi-squared test
  upper_quantile <- qchisq(1-(alpha_val/2), deg_freedom)
  lower_quantile <- qchisq(alpha_val/2, deg_freedom)
  
  # now filter for those lower than lower quantile, 
  # higher than upper quantile
  data_filtered_by_variance <- dplyr::filter(data_filtered_by_variance, 
                                              test_statistic > upper_quantile | 
                                               test_statistic < lower_quantile
                                              )
  
  # now remove the variance, test_statistic columns that were added
  data_filtered_by_variance <- dplyr::select(data_filtered_by_variance,
                                             -c('variance', 'test_statistic'))
  
  return(data_filtered_by_variance)
  
}

filter_variance <- filter_by_variance(data)

# ----------------------- 4.3: Filter 3 ----------------------------------------
# Filter for genes that have a coefficient of variation > 0.186.
# the coefficient of variation is the ratio of stdev to mean
filter_by_cv <- function(data){
  # get the stdev for each probeset
  st_devs <- data %>%
    dplyr::select(-probesets) %>%
    apply(1, sd) %>%
    as_tibble(rownames='probesets') %>%
    dplyr::rename('st_dev'=value) %>%
    dplyr::mutate(probesets=data$probesets)
  
  # get the mean for each probeset
  means <- data %>%
    dplyr::select(-probesets) %>%
    apply(1, mean) %>%
    as_tibble(rownames='probesets') %>%
    dplyr::rename('mean'=value) %>%
    dplyr::mutate(probesets=data$probesets)
  
  # add stdevs and means to data
  data_filtered_by_cv <- dplyr::full_join(data, st_devs, by='probesets') %>%
                         dplyr::full_join(means, by='probesets')
  
  # now add a column for coefficient of variation
  data_filtered_by_cv <- dplyr::mutate(data_filtered_by_cv,
                                       cv = st_dev/mean)
  
  # last step: filter
  data_filtered_by_cv <- dplyr::filter(data_filtered_by_cv,
                                       cv > 0.186)
  
  # remove columns that we created
  data_filtered_by_cv <- dplyr::select(data_filtered_by_cv,
                                       -c('mean', 'st_dev', 'cv'))
  
  return(data_filtered_by_cv)
  
}

filter_cv <- filter_by_cv(data)

# ----------------------- 4: combine filters -----------------------------------
# This will return a tibble with probes that pass all 3 filters
filter_by_all <- function(data){
  # first, apply all 3 filters to the data individually
  expression_filter <- filter_by_expression(data)
  variance_filter <- filter_by_variance(data)
  cv_filter <- filter_by_cv(data)
  
  # now, we want to join the tibbles so we end up with only 
  # genes contained in ALL 3 filtered tibbles
  all_filters <- dplyr::inner_join(expression_filter,
                                   variance_filter,
                                   by = colnames(data)) %>%
                 dplyr::inner_join(cv_filter,
                                   by = colnames(data))
  
  return(all_filters)
  
}

# 1,531 genes meet all 3 criteria
combined_filters <- filter_by_all(data)

# writing a file with the combined filter result
readr::write_csv(combined_filters, 'fully_filtered_data.csv')

# writing a file with just the second filter results
readr::write_csv(filter_variance, '4_2_filter_results.csv')

