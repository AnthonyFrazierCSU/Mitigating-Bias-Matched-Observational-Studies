################################################################################
####           Mitigating Bias in Matched Observational Studies             ####
####                    with Continuous Treatments                          ####
####                        Simulation Code                                 ####
####                Siyu Heng, Anthony Frazier, Wen Zhou                    ####
################################################################################

#### Libraries -----------------------------------------------------------------
# This section lists libraries you will need to load / install in order to run
# the simulation and obtain all output shown in the manuscript.

library(tidyverse)                                      # dataframe manipulation
library(tidyr)
library(dplyr)
library(nbpMatching)                          # non-bipartite matching algorithm
library(distr)
library(RFCDE)                                  # Conditional Density Estimation
library(ggplot2)                                              # creating figures
library(kableExtra)                                            # creating tables
library(foreach)                                               # parallelization
library(doParallel)                                            
library(devtools)                                        # for installing linCDE
library(randomForest)                    # for centering method used with linCDE
devtools::install_github("ZijunGao/LinCDE", build_vignettes = TRUE)      #linCDE
#### Utility Functions ---------------------------------------------------------
# This section lists useful functions that aren't necessary to perform the 
# simulation, but improve quality of life
#### round_df: rounds all numeric values in a dataframe to a specified digit.

#### CI_maker: creates wald-type CI given a point estimate and variance.

#### coverage_check: checks if paramter of interest lies within CI. 

round_df <- function(x, digits) {
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
CI_maker <- function(estimate, alpha, variance){
  return(c(estimate - (qnorm(1 - (alpha/2)) * sqrt(variance)), estimate + (qnorm(1 - (alpha/2)) * sqrt(variance))))
}
coverage_check <- function(CI, parameter){
  return(ifelse(CI[1] < parameter & CI[2] > parameter, 1, 0))
}

#### Data Generation -----------------------------------------------------------
# This section lists functions that create inputs for our matching and 
# estimation functions. 

#### generate_data_dose: generates treatment (Z), covariates (X1,...,X5), and 
# outcome (Y) variables based on specifications described in the manuscript.

#### propensity_matrix_oracle: calculates the true (oracle) probability of
# observed treatment assignment for all pairs of observations in a dataset.

#### propensity_matrix_linCDE: estimates the probability of observed
# treatment assignment for all pairs of observations. Uses linCDE package
# to estimate the conditional density of treatment given covariates, and employs
# centering via conditional mean estimation using the randomForest package. 

#### propensity_matrix_RFCDE: estimates the probability of observed
# treatment assignment for all pairs of observations. Uses RFCDE package.

generate_data_dose <- function(N){
  # integer N: determines number of observations generated 
  
  # Generate covariates
  x1 = rnorm(N, 0, 1); x2 = rnorm(N, 0, 1); x3 = rnorm(N, 0, 1); x4 = rnorm(N, 0, 1); x5 = rnorm(N, 0, 1);
  
  # Generate treatment variable
  error1 = rnorm(N, 0, 1)
  z = x1 + x2^2 + abs(x3 * x4) + ifelse(x4 > 0, 1, 0) + log(1 + abs(x5)) + error1
  
  # Generate outcome variable
  error2 = rnorm(N, 0, 3)
  y = z + 0.3*x1*z + 0.2*x3^3*z + exp(abs(x2 - x4)) - sin(x5) + error2
  
  # Place variables in a dataframe
  data = data.frame("Z" = z, "Y" = y, "X1" = x1, "X2" = x2, 
                    "X3" = x3, "X4" = x4, "X5" = x5)
  
  return(data)
}
propensity_matrix_oracle <- function(data){
  # dataframe data: data generated from generate_data_dose function
  
  # vector holding probability of treatment given one's own covariates
  same_propensity_vector <- 1:nrow(data)
  
  # vector holding center of conditional density of treatment given covariates
  location <- 1:nrow(data)
  
  for(i in 1:nrow(data)){
    
    # calculate the center observation i's conditional density of treatment
    location[i] <- data[i, 3] + data[i, 4]^2 + abs(data[i, 5] * data[i, 6]) 
                   + ifelse(data[i, 6] > 0, 1, 0) + log(1 + abs(data[i, 7]))
    
    # calculate true (oracle) conditional probability of treatment
    same_propensity_vector[i] <- dnorm(data[i, 1], location[i], 1)
    
  }
  
  # holds joint density values related to observed treatment assignment
  observed_treatment_matrix <- same_propensity_vector %*% base::t(same_propensity_vector)
  
  # will hold output (matrix of probabilities of observed treatment assignment)
  p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # calculate probability of opposing treatment assignment for each pair
      opposing_treatment <- dnorm(data[i, 1], location[j], 1) * dnorm(data[j, 1], location[i], 1)
      
      # calculates probability of observed treatment assignment
      p_matrix[i, j] <- observed_treatment_matrix[i, j] / (observed_treatment_matrix[i, j] + opposing_treatment)
      
      # ensures p_matrix is symmetric
      p_matrix[j, i] <- p_matrix[i, j]
      
    }
  }
  
  return(p_matrix)
}
propensity_matrix_linCDE <- function(data){
  # dataframe data: data generated from generate_data_dose function  
  
  # estimating conditional mean of treatment given covariates, fitting linCDE on residuals
  conditional_mean <- randomForest(x = data[, 3:7], y = data$Z)
  residuals <- data$Z - conditional_mean$predicted
  model <- LinCDE::LinCDE.boost(X = data[, 3:7], y = residuals, terminalSize = 20, verbose = FALSE)
  
  # will hold conditional density estimates for every combination of treatment
  # and covariates
  propensity_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:nrow(propensity_matrix)){
    propensity_matrix[i, ] <- predict(model, X = data[i, 3:7], y = residuals)
  }
  
  # will hold output (matrix of probabilities of observed treatment assignment)
  p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # calculates observed treatment assignment density for each pair
      observed_treatment <- propensity_matrix[i, i] * propensity_matrix[j, j]
      
      # calculate opposing treatment assignment density for each pair
      opposing_treatment <- propensity_matrix[j, i] * propensity_matrix[i, j]
      
      # calculates probability of observed treatment assignment
      p_matrix[i, j] <- observed_treatment / (observed_treatment + opposing_treatment)
      
      # ensures p_matrix is symmetric
      p_matrix[j, i] <- p_matrix[i, j]
      
    }
  }  
  
  return(p_matrix)
}
propensity_matrix_RFCDE <- function(data){
  # dataframe data: data generated from generate_data_dose function  
  
  model <- RFCDE(data[, 3:7], data$Z)
  predicted <- predict(model, as.matrix(data[, 3:7]), "CDE", data$Z)
  same_propensity_vector <- 1:nrow(data)
  
  for(i in 1:nrow(data)){
    
    same_propensity_vector[i] <- predicted[i, i] # "given" is first arg, predicting the second arg
    
  }
  
  observed_treatment_matrix <- same_propensity_vector %*% base::t(same_propensity_vector)
  p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      opposing_treatment <- predicted[j, i] * predicted[i, j]
      p_matrix[i, j] <- observed_treatment_matrix[i, j] / (observed_treatment_matrix[i, j] + opposing_treatment)
      p_matrix[j, i] <- p_matrix[i, j]
      
    }
  }
  
  return(p_matrix)
}

#### Matching Methods ----------------------------------------------------------
# this section lists the different matching algorithsm we compare in our 
# simulation. 

#### nbp_balanced: creates a set of matched pairs without a caliper

#### nbp_caliper: creates a set of matched pairs using a caliper with a 
# cutoff of delta.

nbp_balanced <- function(data, cutoff){
  # dataframe data: data generated from generate_data_dose function
  # double cutoff: function returns NULL if the standardized difference in 
  # covariates is larger than the cutoff for at least one covariate
  
  # the next set of lines uses the nbpMatch library to create matched pairs
  # for data based on mahalanobis distance of covariates.
  test.dist <- gendistance(data[, c(1, 3:7)], idcol = 1)
  test.mdm <- distancematrix(test.dist)
  test.match <- nonbimatch(test.mdm)
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  # the next section of code checks the standardized covariate difference between
  # the low-dose and high-dose groups. If there is at least one coviarate that has
  # a standardized difference larger that cutoff, we will not return the matched
  # pairs generated.
  standard_differences <- 1:5
  sigmas <- 1:5
  
  for(i in 3:7){
    xbar_before <- sum(data[, i])/nrow(data)
    sigma_before <- sqrt(sum((data[, i] - xbar_before)^2)/(nrow(data) - 1))
    xbar_high <- sum(data[which(row.names(data) %in% high_groups), i]) / nrow(data[which(row.names(data) %in% high_groups), ])
    xbar_low <- sum(data[which(row.names(data) %in% low_groups), i]) / nrow(data[which(row.names(data) %in% low_groups), ])
    
    standard_differences[i-2] <- abs((xbar_high - xbar_low) / sigma_before)
    sigmas[i-2] <- sigma_before
  }
  
  standardization_fail <- ifelse(max(standard_differences) > cutoff, TRUE, FALSE)
  balanced_pairs <- as.data.frame(matrix(nrow = 0, ncol = 2))
  
  # if our matched pairs are sufficiently well balanced, we will return the set
  # of matched pairs. If the set of matched pairs are not balanced on covariates,
  # an empty dataframe will be returned
  if(standardization_fail == FALSE){
    
    high_data <- matches[which(matches$group == "H"), ]
    balanced_pairs <- high_data[, c(2, 4)]
    
  }
  
  # the matched pairs we return will always have the higher-dose observation
  # on the first column and the low-dose observation on the second column.
  names(balanced_pairs) <- c("High", "Low")
  output <- list("pairs" = balanced_pairs, "standardized_diffs" = standard_differences)
  
  return(output)
  
}
nbp_caliper <- function(data, p_matrix, delta){
  # dataframe data: The data we want to make matches for
  # matrix p_matrix: matrix containing probability of observed treatment
  # assignment for each pair of observations
  # double delta: used for the propensity caliper
  
  # matrix encoding our caliper penalty
  caliper_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # if probability of observed treatment assignment is extreme, enforce
      # infinite penalty on distance matrix
      caliper_matrix[i, j] <- ifelse(p_matrix[i, j] > (1-delta) | p_matrix[i, j] < delta, Inf, 0)
      caliper_matrix[j, i] <- caliper_matrix[i, j]
      
    }
  }
  
  # making mahalanobis distance matrix
  
  test.dist <- gendistance(data[, c(1, 3:7)], idcol = 1)
  test.mdm <- distancematrix(test.dist)
  
  # applying caliper to mahalnobis distance matrix
  
  caliper_dist_matrix <- distancematrix(test.mdm + caliper_matrix)
  test.match <- nonbimatch(caliper_dist_matrix)
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  # calculates standardized covariate differences, if we want to compare the 
  # standardized covariate difference between mahalanobis-only matching and
  # calipered-matching
  standard_differences <- 1:5
  sigmas <- 1:5
  
  for(i in 3:7){
    xbar_before <- sum(data[, i])/nrow(data)
    sigma_before <- sqrt(sum((data[, i] - xbar_before)^2)/(nrow(data) - 1))
    xbar_high <- sum(data[which(row.names(data) %in% high_groups), i]) / nrow(data[which(row.names(data) %in% high_groups), ])
    xbar_low <- sum(data[which(row.names(data) %in% low_groups), i]) / nrow(data[which(row.names(data) %in% low_groups), ])
    
    standard_differences[i-2] <- abs((xbar_high - xbar_low) / sigma_before)
    sigmas[i-2] <- sigma_before
  }
  
  balanced_pairs <- as.data.frame(matrix(nrow = 0, ncol = 2))
  high_data <- matches[which(matches$group == "H"), ]
  balanced_pairs <- high_data[, c(2, 4)]
  names(balanced_pairs) <- c("High", "Low")
  
  output <- list("pairs" = balanced_pairs, "standardized_diffs" = standard_differences)
  
  return(output)  
}

#### Estimation Methods --------------------------------------------------------
# This section lists the different difference-in-means estimators we compare
# in our simulation

#### true_lambda: calculates what the estimand of interest is. Used so we
# can calculate bias, variance, and MSE of our estimators. 

#### naive_estimator: the difference-in-means estimator that is not bias-
# corrected.

#### bias_corrected_estimator: the difference-in-means estimator that includes
# the bias-correction terms. 

true_lambda <- function(data, pairs){
  # dataframe data: the data we want to perform a difference-in-means estimator on
  # dataframe pairs: the set of matched pairs we're using to calculate our
  # difference-in-means.
  
  numerator <- 1:nrow(pairs)
  denominator <- 1:nrow(pairs)
  
  for(i in 1:nrow(pairs)){
    
    # within each matched pair, we calculate the potential outcome for each 
    # observation under each observed treatment dose (within that matched pair).
    
    high_cf_high <- data[pairs[i, 1], 1] + 0.3*data[pairs[i, 1], 3]*data[pairs[i, 1], 1] 
    + 0.2*data[pairs[i, 1], 5]^3*data[pairs[i, 1], 1] + exp(abs(data[pairs[i, 1], 4] 
                                                                - data[pairs[i, 2], 6])) - sin(data[pairs[i, 1], 7])
    
    high_cf_low <- data[pairs[i, 2], 1] + 0.3*data[pairs[i, 1], 3]*data[pairs[i, 2], 1] 
    + 0.2*data[pairs[i, 1], 5]^3*data[pairs[i, 2], 1] + exp(abs(data[pairs[i, 1], 4] 
                                                                - data[pairs[i, 1], 6])) - sin(data[pairs[i, 1], 7])
    
    low_cf_low <- data[pairs[i, 2], 1] + 0.3*data[pairs[i, 2], 3]*data[pairs[i, 2], 1] 
    + 0.2*data[pairs[i, 2], 5]^3*data[pairs[i, 2], 1] + exp(abs(data[pairs[i, 2], 4] 
                                                                - data[pairs[i, 2], 6])) - sin(data[pairs[i, 2], 7])
    
    low_cf_high <- data[pairs[i, 1], 1] + 0.3*data[pairs[i, 2], 3]*data[pairs[i, 1], 1] 
    + 0.2*data[pairs[i, 2], 5]^3*data[pairs[i, 1], 1] + exp(abs(data[pairs[i, 2], 4] 
                                                                - data[pairs[i, 2], 6])) - sin(data[pairs[i, 2], 7])
    
    # calculates the numerator and denominator portions of lambda
    numerator[i] <- high_cf_high - high_cf_low + low_cf_high - low_cf_low
    denominator[i] <- 2*(data[pairs[i, 1], 1] - data[pairs[i, 2], 1])
    
  }
  
  # returns lambda
  return(sum(numerator) / sum(denominator))  
  
}
naive_estimator <- function(data, pairs){
  # dataframe data: the data we want to perform a difference-in-means estimator on
  # dataframe pairs: the set of matched pairs we're using to calculate our
  # difference-in-means.
    
  numerator <- 1:nrow(pairs)
  denominator <- 1:nrow(pairs)
  
  for(i in 1:nrow(pairs)){
    
    # because the set of matched pairs are ordered such that observations in the
    # first column are always the higher dose observation, we calculate the
    # numerator of our estimator just by calculating the first term, associated
    # with the 1st index observation being the higher dose observation. 
    numerator[i] <- data[pairs[i, 1], 2] - data[pairs[i, 2], 2]
    
    # difference in dose
    denominator[i] <- data[pairs[i, 1], 1] - data[pairs[i, 2], 1]
  }
  
  return(sum(numerator) / sum(denominator))
  
}
bias_corrected_estimator <- function(data, pairs, p_matrix, cutoff){
  # dataframe data: the data we want to perform a difference-in-means estimator on
  # dataframe pairs: the set of matched pairs we're using to calculate our
  # difference-in-means.
  # matrix p_matrix: matrix containing probability of observed treatment
  # assignment for each pair of observations
  # cutoff: allows user to regularize p_ij
  
  numerator <- 1:nrow(pairs)
  denominator <- 1:nrow(pairs)
  
  for(i in 1:nrow(pairs)){
    
    p = p_matrix[pairs[i, 1], pairs[i, 2]]
    
    if(p < cutoff){
      p <- cutoff
    }
    
    if(p > (1-cutoff)){
      p <- (1-cutoff)
    }
    
    # because the set of matched pairs are ordered such that observations in the
    # first column are always the higher dose observation, we calculate the
    # numerator of our estimator just by calculating the first term, associated
    # with the 1st index observation being the higher dose observation. 
    numerator[i] <- p^(-1)*(data[pairs[i, 1], 2] - data[pairs[i, 2], 2])
    
    denominator[i] <- 2*(data[pairs[i, 1], 1] - data[pairs[i, 2], 1])
  }
  
  return(sum(numerator) / sum(denominator))
  
}

#### Variance Estimators -------------------------------------------------------
# This section lists the different conservative variance estimators we use.

#### variance_lambda: calculates the conservative variance estimate for the naive
# estimator. excluding use of probability of observed treatment assignment. 

#### variance_bias-corrected: calculates the conservative variance estimate for the
# bias-corrected estimator. 

variance_naive <- function(data, pairs, estimate){
  numerator <- sum(((data[pairs[, 1], 2] - data[pairs[, 2], 2]) - estimate)^2)
  denominator <- 2*(sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])^2)
  return(numerator/denominator)
}
variance_bias_corrected <- function(data, pairs, p_matrix, estimate, cutoff){
  
  numerator <- 1:nrow(pairs)
  
  for(i in 1:length(numerator)){
    
    p = p_matrix[pairs[i, 1], pairs[i, 2]]
    delta <- p
    
    if(p < cutoff){
      delta <- cutoff
    }
    
    if(p > (1-cutoff)){
      delta <- (1-cutoff)
    }
    
    numerator[i] <- max(p, (1-p))*(delta^(-1)*((data[pairs[i, 1], 2] - data[pairs[i, 2], 2]) - estimate))^2
    
  }
  
  denominator <- 2*sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])^2
  return(sum(numerator)/denominator)
  
}

#### Simulation (MAE, RMSE, MCV, CR) -------------------------------
set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0

while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_dose(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

# Parallelization setup
numcores <- detectCores()
registerDoParallel(numcores - 1)

# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

output <- foreach(i = 1:iter, .combine = rbind,
                  .packages = c("nbpMatching",
                                "RFCDE",
                                "randomForest")) %dopar%
  {
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    p_mat_linCDE <- propensity_matrix_linCDE(test)
    p_mat_RFCDE <- propensity_matrix_RFCDE(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    RFCDE_pairs <- nbp_caliper(test, p_mat_RFCDE, 0.1)$pairs
    linCDE_pairs <- nbp_caliper(test, p_mat_linCDE, 0.1)$pairs
    
    ######################################################################################################################
    
    # estimators with mahalanobis only pairs, naive estimator
    lambda <- true_lambda(test, mahalanobis_pairs)
    e_naive <- naive_estimator(test, mahalanobis_pairs)
    var_naive <- variance_naive(test, mahalanobis_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_naive, var_naive, CI_naive, coverage_naive, "Naive", "Mahalanobis"))
    
    # estimators with mahalanobis only pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, var_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "Mahalanobis"))
    
    # estimators with mahalanobis only pairs, RFCDE estimator
    e_linCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, var_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "Mahalanobis"))
    
    ######################################################################################################################
    
    # estimators with RFCDE pairs, naive estimator
    lambda <- true_lambda(test, RFCDE_pairs)
    e_naive <- naive_estimator(test, RFCDE_pairs)
    var_naive <- variance_naive(test, RFCDE_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_naive, var_naive, CI_naive, coverage_naive, "Naive", "RFCDE_Caliper"))
    
    # estimators with RFCDE pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, var_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "RFCDE_Caliper"))
    
    # estimators with RFCDE pairs, RFCDE estimator
    e_linCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, var_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "RFCDE_Caliper"))
    
    ######################################################################################################################
    
    # estimators with linCDE pairs, naive estimator
    lambda <- true_lambda(test, linCDE_pairs)
    e_naive <- naive_estimator(test, linCDE_pairs)
    var_naive <- variance_naive(test, linCDE_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_naive, var_naive, CI_naive, coverage_naive, "Naive", "linCDE_Caliper"))
    
    # estimators with linCDE pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, var_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "linCDE_Caliper"))
    
    # estimators with linCDE pairs, RFCDE estimator
    e_linCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, var_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "linCDE_Caliper"))
    
    names(temp_data) <- c("lambda", "estimate", "variance", "lower_CI", "upper_CI", "coverage", "estimation_method", "matching_method")
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

output[, -c(7, 8)] <- lapply(output[, -c(7, 8)], as.numeric)

#write.csv(output, "sim_output.csv")

matching_method_vector <- unique(output$matching_method)
estimation_method_vector <- unique(output$estimation_method)
summary_data <- data.frame(matrix(nrow = 0, ncol = 4))
for(i in 1:length(estimation_method_vector)){
  
  method_data <- output[which(output$matching_method == matching_method_vector[i]), ]
  
  for(m in 1:length(matching_method_vector)){
    
    subdata <- method_data[which(method_data$estimation_method == estimation_method_vector[m]), ]
    mean_abs_error <- mean(abs(subdata[, 1] - subdata[, 2]))
    root_MSE <- sqrt(mean((subdata[, 1] - subdata[, 2])^2))
    variance_average <- mean(subdata[, 3])
    coverage_rate <- mean(subdata[, 6])
    
    summary_data <- rbind(summary_data, c(matching_method_vector[i], estimation_method_vector[m], 
                                          mean_abs_error, root_MSE, variance_average, coverage_rate))
    
  } 
}

summary_data[, -c(1, 2)] <- lapply(summary_data[, -c(1, 2)], as.numeric)
names(summary_data) <- c("Matching_Method", "Estimation_Method","Mean_Abs_Error", "RMSE", "Average_Variance", "Coverage_Rate")

#write.csv(summary_data, "sim_summary.csv")

#### Recreating Table 1 (MAE, RMSE table) ----

summary_data_rounded <- round_df(summary_data, 3)
table_data <- summary_data_rounded[, 1:4]
table_data <- do.call("rbind", list(table_data[1:3, ], table_data[7:9, ], table_data[4:6, ]))
kbl(table_data, align = 'c', booktabs = T, format = "latex",
    caption = "Simulation Results (Bias, MSE)", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  collapse_rows(columns = 1, latex_hline = "custom",  valign = "middle")

#### Recreating Table 2 (MCV, CR table) ----

summary_data_rounded <- round_df(summary_data, 3)
table_data <- summary_data_rounded[, c(1, 2, 5, 6)]
table_data <- do.call("rbind", list(table_data[1:3, ], table_data[7:9, ], table_data[4:6, ]))
kbl(table_data, align = 'c', booktabs = T, format = "latex",
    caption = "Simulation Results (Bias, MSE)", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  collapse_rows(columns = 1, latex_hline = "custom",  valign = "middle")






