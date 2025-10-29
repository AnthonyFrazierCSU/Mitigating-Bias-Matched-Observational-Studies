################################################################################
####           Bias Mitigation in Matched Observational Studies             ####
####                    with Continuous Treatments                          ####
####                        Simulation Code                                 ####
####                Anthony Frazier, Siyu Heng, Wen Zhou                    ####
################################################################################

################################################################################
####                        Major Functions                                 ####
################################################################################

#### Libraries -----------------------------------------------------------------
# This section lists libraries you will need to load / install in order to run
# the simulation and obtain all output shown in the manuscript.

#library(optmatch)                                                   # for bipartite matching
library(tidyverse)                                                  # for dataframe manipulation
library(randomForest)
library(nbpMatching)                                                # non-bipartite matching algorithm
library(kableExtra)
library(tidyr)
library(dplyr)
library(earth)
library("foreach")
library("doParallel")
library(pak)
pak::pkg_install("ZijunGao/LinCDE")
remotes::install_github("lee-group-CMU/RFCDE/r")
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

# for main manuscript
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
    location[i] <- data[i, 3] + data[i, 4]^2 + abs(data[i, 5] * data[i, 6]) + ifelse(data[i, 6] > 0, 1, 0) + log(1 + abs(data[i, 7]))
    
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
propensity_matrix_oracle2 <- function(data){
  # dataframe data: data generated from generate_data_dose function
  
  # vector holding probability of treatment given one's own covariates
  same_propensity_vector <- 1:nrow(data)
  
  # vector holding center of conditional density of treatment given covariates
  location <- 1:nrow(data)
  
  for(i in 1:nrow(data)){
    
    # calculate the center observation i's conditional density of treatment
    location[i] <- ifelse(data[i, 3] == 1 & data[i, 4] > 5, 4, 0) - data[i, 3] + data[i, 5] + (sin(data[i, 7]))^2 + 2*log(1 +data[i, 6]) + 2*exp(-abs(data[i, 8]))
    
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
  model <- LinCDE::LinCDE.boost(y = data[, 1], X = data[, 3:ncol(data)], terminalSize = 20, verbose = FALSE, centering = TRUE, centeringMethod = "linearRegression")
  
  # will hold conditional density estimates for every combination of treatment
  # and covariates
  propensity_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:nrow(propensity_matrix)){
    propensity_matrix[i, ] <- predict(model, X = data[i, 3:ncol(data)], y = data[, 1])
  }
  
  # propensity_matrix_2 <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  # propensity_matrix_2 <- sapply(propensity_matrix, predict(model, X = data[, 3:ncol(data)], y = data[, 1]))
  
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
  
  # Fit RFCDE model
  model <- RFCDE(data[, 3:ncol(data)], data[, 1])
  
  # creates matrix of predicted CDE values of Z given covariates (X)
  predicted <- predict(model, as.matrix(data[, 3:ncol(data)]), "CDE", data[, 1])
  
  # holds f(Z_n|X_n) for all observations. 
  same_propensity_vector <- 1:nrow(data)
  
  for(i in 1:nrow(data)){
    
    same_propensity_vector[i] <- predicted[i, i] # "given" is first arg, predicting the second arg
    
  }
  
  # holds joint density values related to observed treatment assignment
  observed_treatment_matrix <- same_propensity_vector %*% base::t(same_propensity_vector)
  
  # will hold output (matrix of probabilities of observed treatment assignment)
  p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # calculate probability of opposing treatment assignment for each pair
      opposing_treatment <- predicted[j, i] * predicted[i, j]
      
      # calculates probability of observed treatment assignment
      p_matrix[i, j] <- observed_treatment_matrix[i, j] / (observed_treatment_matrix[i, j] + opposing_treatment)
      
      # ensures p_matrix is symmetric
      p_matrix[j, i] <- p_matrix[i, j]
      
    }
  }
  
  return(p_matrix)
}
propensity_matrix_model <- function(data){
  # dataframe data: data generated from generate_data_dose function
  
  # create estimates for E[Z|x] for all observations
  model <- earth(data[, 1] ~., data = data[, 3:ncol(data)], degree = 2)
  fitted <- model$fitted.values
  
  # create estimate for variance
  sd_est <- sd(fitted - data[, 1])
  
  # vector holding probability of treatment given one's own covariates
  same_propensity_vector <- 1:nrow(data)
  
  # vector holding center of conditional density of treatment given covariates
  location <- 1:nrow(data)
  
  for(i in 1:nrow(data)){

    # calculate true (oracle) conditional probability of treatment
    same_propensity_vector[i] <- dnorm(data[i, 1], fitted[i], sd_est)
    
  }
  
  # holds joint density values related to observed treatment assignment
  observed_treatment_matrix <- same_propensity_vector %*% base::t(same_propensity_vector)
  
  # will hold output (matrix of probabilities of observed treatment assignment)
  p_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # calculate probability of opposing treatment assignment for each pair
      opposing_treatment <- dnorm(data[i, 1], fitted[j], sd_est) * dnorm(data[j, 1], fitted[i], sd_est)
      
      # calculates probability of observed treatment assignment
      p_matrix[i, j] <- observed_treatment_matrix[i, j] / (observed_treatment_matrix[i, j] + opposing_treatment)
      
      # ensures p_matrix is symmetric
      p_matrix[j, i] <- p_matrix[i, j]
      
    }
  }
  
  return(p_matrix)
}

# for supplementary materials
generate_data_uniform <- function(N){
  # integer N: determines number of observations generated 
  
  # Generate covariates
  x1 = runif(N, -1, 1); x2 = runif(N, -1, 1); x3 = runif(N, -1, 1); x4 = runif(N, -1, 1); x5 = runif(N, -1, 1);
  
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
generate_data_logistic <- function(N){
  # integer N: determines number of observations generated 
  
  # Generate covariates
  x1 = rlogis(N); x2 = rlogis(N); x3 = rlogis(N); x4 = rlogis(N); x5 = rlogis(N);
  
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
generate_data_dose_2 <- function(N){
  # integer N: determines number of observations generated 
  
  # Generate covariates
  x1 = rbinom(N, 1, 0.5); x2 = rbinom(N, 10, 0.75); x3 = rpois(N, 1.5); x4 = runif(N, 0, 3); x5 = runif(N, -1, 1); x6 = runif(N, -5, 5);
  
  # Generate treatment variable
  error1 = rnorm(N)
  z = ifelse(x1 == 1 & x2 > 5, 4, 0) - x1 + x3 + (sin(x5))^2 + 2*log(1 + x4) + 2*exp(-abs(x6)) + error1
  
  # Generate outcome variable
  error2 = rnorm(N, 0, 3)
  
  y = z + 0.7*x6*z + 2*x5^3*z + ifelse(x1 == 1, x3*(2*x4 + 1), x3) + 0.5*x4*abs(x5) + x2 + error2
  
  # Place variables in a dataframe
  data = data.frame("Z" = z, "Y" = y, "X1" = x1, "X2" = x2, 
                    "X3" = x3, "X4" = x4, "X5" = x5, "X6" = x6)
  
  return(data)
}

#### Matching Methods ----------------------------------------------------------
# this section lists the different matching algorithms we compare in our 
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
  test.dist <- gendistance(data[, c(1, 3:ncol(data))], idcol = 1)
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
  standard_differences <- 1:(ncol(data) - 2)
  sigmas <- 1:(ncol(data) - 2)
  
  for(i in 3:ncol(data)){
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
  
  test.dist <- gendistance(data[, c(1, 3:ncol(data))], idcol = 1)
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
  standard_differences <- 1:(ncol(data) - 3)
  sigmas <- 1:(ncol(data) - 3)
  
  for(i in 3:ncol(data)){
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
    
    high_cf_high <- data[pairs[i, 1], 1] + 0.3*data[pairs[i, 1], 3]*data[pairs[i, 1], 1] + 0.2*data[pairs[i, 1], 5]^3*data[pairs[i, 1], 1] + exp(abs(data[pairs[i, 1], 4] - data[pairs[i, 2], 6])) - sin(data[pairs[i, 1], 7])
    
    high_cf_low <- data[pairs[i, 2], 1] + 0.3*data[pairs[i, 1], 3]*data[pairs[i, 2], 1] + 0.2*data[pairs[i, 1], 5]^3*data[pairs[i, 2], 1] + exp(abs(data[pairs[i, 1], 4] - data[pairs[i, 1], 6])) - sin(data[pairs[i, 1], 7])
    
    low_cf_low <- data[pairs[i, 2], 1] + 0.3*data[pairs[i, 2], 3]*data[pairs[i, 2], 1] + 0.2*data[pairs[i, 2], 5]^3*data[pairs[i, 2], 1] + exp(abs(data[pairs[i, 2], 4] - data[pairs[i, 2], 6])) - sin(data[pairs[i, 2], 7])
    
    low_cf_high <- data[pairs[i, 1], 1] + 0.3*data[pairs[i, 2], 3]*data[pairs[i, 1], 1] + 0.2*data[pairs[i, 2], 5]^3*data[pairs[i, 1], 1] + exp(abs(data[pairs[i, 2], 4] - data[pairs[i, 2], 6])) - sin(data[pairs[i, 2], 7])
    
    # calculates the numerator and denominator portions of lambda
    numerator[i] <- high_cf_high - high_cf_low + low_cf_high - low_cf_low
    denominator[i] <- 2*(data[pairs[i, 1], 1] - data[pairs[i, 2], 1])
    
  }
  
  # returns lambda
  return(sum(numerator) / sum(denominator))  
  
}
true_lambda_2 <- function(data, pairs){
  # dataframe data: the data we want to perform a difference-in-means estimator on
  # dataframe pairs: the set of matched pairs we're using to calculate our
  # difference-in-means.
  
  numerator <- 1:nrow(pairs)
  denominator <- 1:nrow(pairs)
  
  for(i in 1:nrow(pairs)){
    
    # within each matched pair, we calculate the potential outcome for each 
    # observation under each observed treatment dose (within that matched pair).
    
    high_cf_high <- data[pairs[i, 1], 1] + 0.7*data[pairs[i, 1], 8]*data[pairs[i, 1], 1] + 2*(data[pairs[i, 1], 7]^3)*data[pairs[i, 1], 1] + ifelse(data[pairs[i, 1], 3] == 1, data[pairs[i, 1], 5]*(2*data[pairs[i, 1], 6] + 1), data[pairs[i, 1], 5]) + 0.5*data[pairs[i, 1], 6]*abs(data[pairs[i, 1], 7]) + data[pairs[i, 1], 4]
    
    high_cf_low <- data[pairs[i, 2], 1] + 0.7*data[pairs[i, 1], 8]*data[pairs[i, 2], 1] + 2*(data[pairs[i, 1], 7]^3)*data[pairs[i, 2], 1] + ifelse(data[pairs[i, 1], 3] == 1, data[pairs[i, 1], 5]*(2*data[pairs[i, 1], 6] + 1), data[pairs[i, 1], 5]) + 0.5*data[pairs[i, 1], 6]*abs(data[pairs[i, 1], 7]) + data[pairs[i, 1], 4]
    
    low_cf_low <- data[pairs[i, 2], 1] + 0.7*data[pairs[i, 2], 8]*data[pairs[i, 2], 1] + 2*(data[pairs[i, 2], 7])^3*data[pairs[i, 2], 1] + ifelse(data[pairs[i, 2], 3] == 1, data[pairs[i, 2], 5]*(2*data[pairs[i, 2], 6] + 1), data[pairs[i, 2], 5]) + 0.5*data[pairs[i, 2], 6]*abs(data[pairs[i, 2], 7]) + data[pairs[i, 2], 4]
    
    low_cf_high <- data[pairs[i, 1], 1] + 0.7*data[pairs[i, 2], 8]*data[pairs[i, 1], 1] + 2*(data[pairs[i, 2], 7])^3*data[pairs[i, 1], 1] + ifelse(data[pairs[i, 2], 3] == 1, data[pairs[i, 2], 5]*(2*data[pairs[i, 2], 6] + 1), data[pairs[i, 2], 5]) + 0.5*data[pairs[i, 2], 6]*abs(data[pairs[i, 2], 7]) + data[pairs[i, 2], 4]
    
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

#### variance_naive: calculates the conservative variance estimate for the naive
# estimator. excluding use of probability of observed treatment assignment. 

#### variance_bias-corrected: calculates the conservative variance estimate for the
# bias-corrected estimator. 

variance_naive <- function(data, pairs, estimate){
  
  # calculates a, centering constant that makes variance estimate least conservative
  # (in expectation)
  a <- estimate * (sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])) / (nrow(pairs))
  
  # numerator term of variance estimate
  numerator <- sum(((data[pairs[, 1], 2] - data[pairs[, 2], 2]) - a)^2)
  
  # denominator term of variance estimate
  denominator <- sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])^2
  estimate <- (numerator/denominator)*(nrow(pairs)/(nrow(pairs) - 1))
  return(estimate)
}
variance_bias_corrected <- function(data, pairs, p_matrix, estimate, cutoff){
  
  # calculates centering terms a1, a2, a3 used in variance estimation of bias
  # corrected method. 
  
  a <- estimate * (sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])) / (nrow(pairs))
  
  # Calculates each of the three sums involved in the numerator term of the variance
  # estimate for the bias-corrected method
  low_num <- rep(0, nrow(pairs))
  mid_num <- rep(0, nrow(pairs))
  high_num <- rep(0, nrow(pairs))
  denominator <- 1:nrow(pairs)
  
  for(i in 1:nrow(pairs)){
    
    p = p_matrix[pairs[i, 1], pairs[i, 2]]
    
    if(p < cutoff){
      low_num[i] <- ((cutoff)^(-1)*(data[pairs[i, 1], 2] - data[pairs[i, 2], 2]) - a)^2
    }
    
    if(p <= (1-cutoff) & p >= cutoff){
      mid_num[i] <- ((p)^(-1)*(data[pairs[i, 1], 2] - data[pairs[i, 2], 2]) - a)^2
    }
    
    if(p > (1-cutoff)){
      high_num[i] <- ((1 - cutoff)^(-1)*(data[pairs[i, 1], 2] - data[pairs[i, 2], 2]) - a)^2
    }
    
  }
  
  # denominator term of variance estimate
  denominator <- 4*sum(data[pairs[, 1], 1] - data[pairs[, 2], 1])^2
  estimate <- ((sum(low_num) + sum(mid_num) + sum(high_num)) / denominator)*(nrow(pairs)/(nrow(pairs) - 1))
  return(estimate)
}


################################################################################
####                        Main Simulation                                 ####
################################################################################

# Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR) -------------------------------
set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance <- data.frame(matrix(0, nrow = iter, ncol = 5))

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
    covariate_balance[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance, "covariate_balance_sim1.csv")

# Parallelization setup
numcores <- detectCores()
output <- data.frame(matrix(nrow = 0, ncol = 8))
registerDoParallel(numcores - 6)
  
  # now, we create matched pairs and estimate treatment effects with each matching
  # and estimation method we're interested in.
  
  small_output <- foreach(i = 1:iter, .combine = rbind,
                    .packages = c("nbpMatching",
                                  "RFCDE",
                                  "randomForest",
                                  "earth")) %dopar%
    {
      # holds true lambda and estimates in dataframe
      temp_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
      
      # extracts ith dataframe from the data list
      test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
      
      # creates p matrices, used for matching and estimation methods
      p_mat_linCDE <- propensity_matrix_linCDE(test)
      p_mat_RFCDE <- propensity_matrix_RFCDE(test)
      p_mat_model <- propensity_matrix_model(test)
      
      # creates matched pairs 
      mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
      RFCDE_pairs <- nbp_caliper(test, p_mat_RFCDE, 0.1)$pairs
      linCDE_pairs <- nbp_caliper(test, p_mat_linCDE, 0.1)$pairs
      model_pairs <- nbp_caliper(test, p_mat_model, 0.1)$pairs
      
      ######################################################################################################################
      
      # estimators with mahalanobis only pairs, naive estimator
      lambda <- true_lambda(test, mahalanobis_pairs)
      e_naive <- naive_estimator(test, mahalanobis_pairs)
      var_naive <- variance_naive(test, mahalanobis_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, model estimator
      e_model <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "Mahalanobis"))    
      
      ######################################################################################################################
      
      # estimators with RFCDE pairs, naive estimator
      lambda <- true_lambda(test, RFCDE_pairs)
      e_naive <- naive_estimator(test, RFCDE_pairs)
      var_naive <- variance_naive(test, RFCDE_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "RFCDE_Caliper"))

      # estimators with RFCDE pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "RFCDE_Caliper"))
      
      # estimators with RFCDE pairs, RFCDE estimator
      e_linCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "RFCDE_Caliper"))
      
      # estimators with RFCDE pairs, linCDE estimator
      e_model <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, RFCDE_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "RFCDE_Caliper"))    
      
      ######################################################################################################################
      
      # estimators with linCDE pairs, naive estimator
      lambda <- true_lambda(test, linCDE_pairs)
      e_naive <- naive_estimator(test, linCDE_pairs)
      var_naive <- variance_naive(test, linCDE_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "linCDE_Caliper"))

      # estimators with linCDE pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "linCDE_Caliper"))
      
      # estimators with linCDE pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "linCDE_Caliper"))
      
      # estimators with linCDE pairs, model estimator
      e_model <- bias_corrected_estimator(test, linCDE_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, linCDE_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "linCDE_Caliper"))  
      
      ######################################################################################################################
      
      # estimators with model pairs, naive estimator
      lambda <- true_lambda(test, model_pairs)
      e_naive <- naive_estimator(test, model_pairs)
      var_naive <- variance_naive(test, model_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "model_Caliper"))

      # estimators with model pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, model_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, model_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "model_Caliper"))
      
      # estimators with model pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, model_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, model_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "model_Caliper"))
      
      # estimators with model pairs, model estimator
      e_model <- bias_corrected_estimator(test, model_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, model_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "model_Caliper"))     
      names(temp_data) <- c("lambda", "estimate", "MOE", "lower_CI", "upper_CI", "coverage", "estimation_method", "matching_method")
      
      return(temp_data)
      
    }
  
  stopImplicitCluster()
  small_output[, -c(7, 8)] <- lapply(small_output[, -c(7, 8)], as.numeric)
  output <- rbind(output, small_output)

# Now, we'll compile the results of our simulation

output[, -c(7, 8)] <- lapply(output[, -c(7, 8)], as.numeric)
write.csv(output, "sim_output.csv")

matching_method_vector <- unique(output$matching_method)
estimation_method_vector <- unique(output$estimation_method)
summary_data <- data.frame(matrix(nrow = 0, ncol = 4))
for(i in 1:length(matching_method_vector)){
  
  method_data <- output[which(output$matching_method == matching_method_vector[i]), ]
  
  for(m in 1:length(estimation_method_vector)){
    
    subdata <- method_data[which(method_data$estimation_method == estimation_method_vector[m]), ]
    mean_abs_error <- mean(subdata[, 2] - subdata[, 1])
    root_MSE <- sqrt(mean((subdata[, 1] - subdata[, 2])^2))
    MOE_average <- mean(subdata[, 3])
    coverage_rate <- mean(subdata[, 6])
    summary_data <- rbind(summary_data, c(matching_method_vector[i], estimation_method_vector[m], 
                                          mean_abs_error, root_MSE, MOE_average, coverage_rate))
  } 
}

summary_data[, -c(1, 2)] <- lapply(summary_data[, -c(1, 2)], as.numeric)
summary_data[, 5] <- 2*summary_data[, 5]
names(summary_data) <- c("Matching_Method", "Estimation_Method","Mean_Error", "RMSE", "Average_Length", "Coverage_Rate")

write.csv(summary_data, "sim_summary.csv")

# Simulation 1 P Distribution (supplementary materials) -----------------------------------------------------------

# Parallelization setup
numcores <- detectCores()
output <- data.frame(matrix(nrow = 0, ncol = 500))

for(k in 1:100){
  
  registerDoParallel(numcores - 6)
  
  # now, we create matched pairs and estimate treatment effects with each matching
  # and estimation method we're interested in.
  
  small_output <- foreach(i = ((10*(k-1) + 1):(10*k)), .combine = rbind,
                          .packages = c("nbpMatching",
                                        "RFCDE",
                                        "randomForest",
                                        "earth")) %dopar%
    {
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    p_mat_linCDE <- propensity_matrix_linCDE(test)
    p_mat_RFCDE <- propensity_matrix_RFCDE(test)
    p_mat_model <- propensity_matrix_model(test)
    
    # creates matched pairs 
    pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)

    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_linCDE[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)

    for(j in 1:nrow(pairs)){
    p_pairs[j] <- p_mat_RFCDE[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)

    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_model[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
    }
 
  stopImplicitCluster()
  output <- rbind(output, small_output)
  print(k)
   
}

average_histograms <- data.frame(matrix(nrow = 4, ncol = 500))

for(i in 1:4){
  subindex <- seq(i, 996 + i, by = 4)
  subdata <- output[subindex, ]
  average_histograms[i, ] <- colMeans(subdata)
}

average_hist <- colMeans(output)
write.csv(average_hist, "average_hist_oracle.csv")

write.csv(average_histograms[1, ], "average_hist_oracle.csv")
write.csv(average_histograms[2, ], "average_hist_linCDE.csv")
write.csv(average_histograms[3, ], "average_hist_RFCDE.csv")
write.csv(average_histograms[4, ], "average_hist_model.csv")

# Simulation 2 Estimation and Inference (supplementary materials) -------------------------------
set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance <- data.frame(matrix(0, nrow = iter, ncol = 6))

while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_dose_2(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

#write.csv(covariate_balance, "covariate_balance_sim2.csv")

# Parallelization setup
numcores <- detectCores()

# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

output <- data.frame(matrix(nrow = 0, ncol = 8))

for(k in 1:(iter/10)){
  
  registerDoParallel(numcores - 6)
  
  # now, we create matched pairs and estimate treatment effects with each matching
  # and estimation method we're interested in.
  
  small_output <- foreach(i = ((10*(k-1) + 1):(10*k)), .combine = rbind,
                          .packages = c("nbpMatching",
                                        "RFCDE",
                                        "randomForest",
                                        "earth")) %dopar%
    {
      # holds true lambda and estimates in dataframe
      temp_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
      
      # extracts ith dataframe from the data list
      test <- as.data.frame(data_list[(8*(i-1) + 1):(8*i)])
      
      # creates p matrices, used for matching and estimation methods
      p_mat_linCDE <- propensity_matrix_linCDE(test)
      p_mat_RFCDE <- propensity_matrix_RFCDE(test)
      p_mat_model <- propensity_matrix_model(test)
      
      # creates matched pairs 
      mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
      RFCDE_pairs <- nbp_caliper(test, p_mat_RFCDE, 0.1)$pairs
      linCDE_pairs <- nbp_caliper(test, p_mat_linCDE, 0.1)$pairs
      model_pairs <- nbp_caliper(test, p_mat_model, 0.1)$pairs
      
      ######################################################################################################################
      
      # estimators with mahalanobis only pairs, naive estimator
      lambda <- true_lambda_2(test, mahalanobis_pairs)
      e_naive <- naive_estimator(test, mahalanobis_pairs)
      var_naive <- variance_naive(test, mahalanobis_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "Mahalanobis"))
      
      # estimators with mahalanobis only pairs, linCDE estimator
      e_model <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "Mahalanobis"))    
      
      ######################################################################################################################
      
      # estimators with RFCDE pairs, naive estimator
      lambda <- true_lambda_2(test, RFCDE_pairs)
      e_naive <- naive_estimator(test, RFCDE_pairs)
      var_naive <- variance_naive(test, RFCDE_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "RFCDE_Caliper"))
      
      # estimators with RFCDE pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "RFCDE_Caliper"))
      
      # estimators with RFCDE pairs, RFCDE estimator
      e_linCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "RFCDE_Caliper"))
      
      # estimators with RFCDE pairs, linCDE estimator
      e_model <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, RFCDE_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "RFCDE_Caliper"))    
      
      ######################################################################################################################
      
      # estimators with linCDE pairs, naive estimator
      lambda <- true_lambda_2(test, linCDE_pairs)
      e_naive <- naive_estimator(test, linCDE_pairs)
      var_naive <- variance_naive(test, linCDE_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "linCDE_Caliper"))
      
      # estimators with linCDE pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "linCDE_Caliper"))
      
      # estimators with linCDE pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "linCDE_Caliper"))
      
      # estimators with linCDE pairs, model estimator
      e_model <- bias_corrected_estimator(test, linCDE_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, linCDE_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "linCDE_Caliper"))  
      
      ######################################################################################################################
      
      # estimators with model pairs, naive estimator
      lambda <- true_lambda_2(test, model_pairs)
      e_naive <- naive_estimator(test, model_pairs)
      var_naive <- variance_naive(test, model_pairs, e_naive)
      CI_naive <- CI_maker(e_naive, 0.05, var_naive)
      coverage_naive <- coverage_check(CI_naive, lambda)
      MOE_naive <- qnorm(0.975) * sqrt(var_naive)
      temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "model_Caliper"))
      
      # estimators with model pairs, RFCDE estimator
      e_RFCDE <- bias_corrected_estimator(test, model_pairs, p_mat_RFCDE, 0.1)
      var_RFCDE <- variance_bias_corrected(test, model_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
      CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
      coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
      MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
      temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "model_Caliper"))
      
      # estimators with model pairs, linCDE estimator
      e_linCDE <- bias_corrected_estimator(test, model_pairs, p_mat_linCDE, 0.1)
      var_linCDE <- variance_bias_corrected(test, model_pairs, p_mat_linCDE, e_linCDE, 0.1)
      CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
      coverage_linCDE <- coverage_check(CI_linCDE, lambda)
      MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
      temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "model_Caliper"))
      
      # estimators with model pairs, model estimator
      e_model <- bias_corrected_estimator(test, model_pairs, p_mat_model, 0.1)
      var_model <- variance_bias_corrected(test, model_pairs, p_mat_model, e_model, 0.1)
      CI_model <- CI_maker(e_model, 0.05, var_model)
      coverage_model <- coverage_check(CI_model, lambda)
      MOE_model <- qnorm(0.975) * sqrt(var_model)
      temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "model_Caliper"))     
      names(temp_data) <- c("lambda", "estimate", "MOE", "lower_CI", "upper_CI", "coverage", "estimation_method", "matching_method")
      
      return(temp_data)
      
    }
  
  stopImplicitCluster()
  small_output[, -c(7, 8)] <- lapply(small_output[, -c(7, 8)], as.numeric)
  output <- rbind(output, small_output)
  print(k)
}

output[, -c(7, 8)] <- lapply(output[, -c(7, 8)], as.numeric)

for(i in 1:(nrow(output)/5)){
  agg_CI <- data.frame("lower_CI" = output[(5*i - 2):(5*i), 2] - output[(5*i - 4), 3], "upper_CI" = output[(5*i - 2):(5*i), 2] + output[(5*i - 4), 3])
  agg_CI <- cbind(agg_CI, "coverage" = ifelse(agg_CI[, 1] < output[(5*i - 4), 1] & agg_CI[, 2] > output[(5*i - 4), 1], 1, 0))
  add_CI <- do.call("cbind", list(output[(5*i - 2):(5*i), 1:2], "MOE" = rep(output[(5*i - 4), 3], 3), agg_CI, output[(5*i - 2):(5*i), 7:8]))
  add_CI$estimation_method <- c("RFCDE_agg", "linCDE_agg", "model_agg")
  output <- rbind(output, add_CI)
}

write.csv(output, "sim_output2.csv")

matching_method_vector <- unique(output$matching_method)
estimation_method_vector <- unique(output$estimation_method)
summary_data <- data.frame(matrix(nrow = 0, ncol = 4))
for(i in 1:length(matching_method_vector)){
  
  method_data <- output[which(output$matching_method == matching_method_vector[i]), ]
  
  for(m in 1:length(estimation_method_vector)){
    
    subdata <- method_data[which(method_data$estimation_method == estimation_method_vector[m]), ]
    mean_abs_error <- mean(subdata[, 2] - subdata[, 1])
    root_MSE <- sqrt(mean((subdata[, 1] - subdata[, 2])^2))
    MOE_average <- mean(subdata[, 3])
    coverage_rate <- mean(subdata[, 6])
    
    summary_data <- rbind(summary_data, c(matching_method_vector[i], estimation_method_vector[m], 
                                          mean_abs_error, root_MSE, MOE_average, coverage_rate))
    
  } 
}

summary_data[, -c(1, 2)] <- lapply(summary_data[, -c(1, 2)], as.numeric)
summary_data[, 5] <- 2*summary_data[, 5]
names(summary_data) <- c("Matching_Method", "Estimation_Method","Mean_Error", "RMSE", "Average_Length", "Coverage_Rate")
summary_data

write.csv(summary_data, "sim_summary2.csv")

# Simulation 2 P Distribution (supplementary materials) -----------------------------------------------------------

# Parallelization setup
numcores <- detectCores()
output <- data.frame(matrix(nrow = 0, ncol = 500))
# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

for(k in 1:100){
  
  registerDoParallel(numcores - 6)
  
  # now, we create matched pairs and estimate treatment effects with each matching
  # and estimation method we're interested in.
  
  small_output <- foreach(i = ((10*(k-1) + 1):(10*k)), .combine = rbind,
                          .packages = c("nbpMatching",
                                        "RFCDE",
                                        "randomForest",
                                        "earth")) %dopar%
    {
      # extracts ith dataframe from the data list
      test <- as.data.frame(data_list[(8*(i-1) + 1):(8*i)])
      
      # creates p matrices, used for matching and estimation methods
      p_mat_oracle <- propensity_matrix_oracle2(test)
      p_mat_linCDE <- propensity_matrix_linCDE(test)
      p_mat_RFCDE <- propensity_matrix_RFCDE(test)
      p_mat_model <- propensity_matrix_model(test)
      
      # creates matched pairs 
      pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
      
      # holds true lambda and estimates in dataframe
      temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
      
      p_pairs <- 1:nrow(pairs)
      for(j in 1:nrow(pairs)){
        p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
      }
      temp_data <- rbind(temp_data, sort(p_pairs))
      names(temp_data) <- 1:nrow(pairs)
      
      for(j in 1:nrow(pairs)){
        p_pairs[j] <- p_mat_linCDE[pairs[j, 1], pairs[j, 2]]
      }
      temp_data <- rbind(temp_data, sort(p_pairs))
      names(temp_data) <- 1:nrow(pairs)
      
      for(j in 1:nrow(pairs)){
        p_pairs[j] <- p_mat_RFCDE[pairs[j, 1], pairs[j, 2]]
      }
      temp_data <- rbind(temp_data, sort(p_pairs))
      names(temp_data) <- 1:nrow(pairs)
      
      for(j in 1:nrow(pairs)){
        p_pairs[j] <- p_mat_model[pairs[j, 1], pairs[j, 2]]
      }
      temp_data <- rbind(temp_data, sort(p_pairs))
      names(temp_data) <- 1:nrow(pairs)
      
      return(temp_data)
      
    }
  
  stopImplicitCluster()
  output <- rbind(output, small_output)
  print(k)
  
}

average_histograms <- data.frame(matrix(nrow = 4, ncol = 500))

for(i in 1:4){
  subindex <- seq(i, 996 + i, by = 4)
  subdata <- output[subindex, ]
  average_histograms[i, ] <- colMeans(subdata)
}

write.csv(average_histograms[1, ], "average_hist_oracle_sim2.csv")
write.csv(average_histograms[2, ], "average_hist_linCDE_sim2.csv")
write.csv(average_histograms[3, ], "average_hist_RFCDE_sim2.csv")
write.csv(average_histograms[4, ], "average_hist_model_sim2.csv")


#### Simulation Tables and Plots:

#### Results
# Simulation Result Tables ----
# creates table showing MAE, RMSE, Mean length and coverage results

# summary_data <- read.csv("sim_summary.csv")[, -1]
summary_data <- read.csv("sim_summary2.csv")[, -1]
rownames(summary_data) <- NULL
table_data <- round_df(summary_data, 3)
table_data <- do.call("rbind", list(table_data[1:4, ], table_data[9:12, ], table_data[5:8, ], table_data[13:16, ]))
kbl(table_data, align = 'c', booktabs = T, format = "latex", row.names = FALSE,
    caption = "Simulation Results (Bias, MSE, ML, CR)", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  collapse_rows(columns = c(1), latex_hline = "custom",  valign = "middle")

# Covariate Balance Tables ----
# average covariate balance tables for both simulation settings. 

balance_tab <- read.csv("covariate_balance_sim1.csv")[, -1]
mean_balance <- data.frame(as.list(round_df(colMeans(balance_tab), 3)))
kbl(mean_balance, align = 'c', booktabs = T, format = "latex", row.names = FALSE,
    caption = "Balance Table, First Simulation", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria")

balance_tab <- read.csv("covariate_balance_sim2.csv")[, -1]
mean_balance <- data.frame(as.list(round_df(colMeans(balance_tab), 3)))
kbl(mean_balance, align = 'c', booktabs = T, format = "latex", row.names = FALSE,
    caption = "Balance Table, Second Simulation", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria")

# P Distribution Plots ----
# average distribution of treatment assignment probabilities

#plot <- as.numeric(read.csv("average_hist_oracle_sim2.csv")[, -1])
#plot <- as.numeric(read.csv("average_hist_linCDE_sim2.csv")[, -1])
#plot <- as.numeric(read.csv("average_hist_RFCDE_sim2.csv")[, -1])
plot <- as.numeric(read.csv("average_hist_model_sim2.csv")[, -1])

ggplot()+
  geom_histogram(aes(x = plot), binwidth = 0.05, color = "black", fill = "white") + 
  theme_bw() +
  #labs(x = "Probability of Observed Treatment Assignment", y = "Frequency") +
  labs(x = "", y = "") +
  theme(axis.text = element_text(size = 12)) + 
  scale_x_continuous(limits = c(0, 1))
ggsave("P_distribution_model_sim2.png", width = 4.5, height = 3)


################################################################################
####                     Supplemental Simulations                           ####
################################################################################

#### Covariate concentration simulation -------------------------------

# uniform ----------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_uniform(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_uniform.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_uniform <- colMeans(output)
hist(average_histogram_uniform)
write.csv(average_histogram_uniform, "average_hist_uniform.csv")

# logistic --------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_logistic(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_logistic.csv")

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
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    
    # creates matched pairs 
    pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_laplace <- colMeans(output)
hist(average_histogram_laplace)
write.csv(average_histogram_laplace, "average_hist_logistic.csv")

# Gaussian ---------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
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
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_gaussian.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_gaussian <- colMeans(output)
hist(average_histogram_gaussian)
write.csv(average_histogram_gaussian, "average_hist_gaussian.csv")

# Cauchy --------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_cauchy(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_cauchy.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_cauchy <- colMeans(output)
hist(average_histogram_cauchy)
write.csv(average_histogram_cauchy, "average_hist_cauchy.csv")


#### Error concentration simulation -------------------------------

# uniform ----------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_uniform_errors(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_uniform_error.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_uniform(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_uniform <- colMeans(output)
hist(average_histogram_uniform)
write.csv(average_histogram_uniform, "average_hist_uniform_error.csv")

# Gaussian ---------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
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
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_gaussian_error.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_oracle(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_gaussian <- colMeans(output)
hist(average_histogram_gaussian)
write.csv(average_histogram_gaussian, "average_hist_gaussian_error.csv")

# laplacian --------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_laplace_errors(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_laplace_error.csv")

# Parallelization setup
numcores <- detectCores()
registerDoParallel(numcores - 1)

# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

output <- foreach(i = 1:iter, .combine = rbind,
                  .packages = c("nbpMatching",
                                "RFCDE",
                                "randomForest",
                                "VGAM")) %dopar%
  {
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_laplace(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_laplace <- colMeans(output)
hist(average_histogram_laplace)
write.csv(average_histogram_laplace, "average_hist_laplace_error.csv")

# Cauchy --------------------------------------------------------------------

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance_df <- data.frame(matrix(0, nrow = iter, ncol = 5))
while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_cauchy_errors(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_balanced(test, 0.2)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance_df[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
}

write.csv(covariate_balance_df, "covariate_balance_cauchy_error.csv")

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
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = nrow(pairs)))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_oracle <- propensity_matrix_cauchy(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    
    p_pairs <- 1:nrow(pairs)
    for(j in 1:nrow(pairs)){
      p_pairs[j] <- p_mat_oracle[pairs[j, 1], pairs[j, 2]]
    }
    temp_data <- rbind(temp_data, sort(p_pairs))
    names(temp_data) <- 1:nrow(pairs)
    
    return(temp_data)
    
  }

stopImplicitCluster()

# Now, we'll compile the results of our simulation

average_histogram_cauchy <- colMeans(output)
hist(average_histogram_cauchy)
write.csv(average_histogram_cauchy, "average_hist_cauchy_error.csv")

#### Appendix E ----------------------------------------------------------------

# Simulation E.1 Estimation and Inference (MAE, RMSE, ML, CR) --------------------

#### nbp_weighted: weights the (i, j)th entry in the distance matrix by 
# (z_i - z_j)^2. Used for Appendix D.
#### nbp_weighted_caliper: also implements our caliper with the nbp_weighted
# matching scheme. 

nbp_weighted <- function(data, cutoff){
  # dataframe data: data generated from generate_data_dose function
  # double cutoff: function returns NULL if the standardized difference in 
  # covariates is larger than the cutoff for at least one covariate
  
  # the next set of lines uses the nbpMatch library to create matched pairs
  # for data based on mahalanobis distance of covariates.
  test.dist <- gendistance(data[, c(1, 3:ncol(data))], idcol = 1)
  test.mdm <- distancematrix(test.dist)
  
  z_weight <- matrix(1, ncol = nrow(data), nrow = nrow(data))
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      z_weight[i, j] <- (data[i, 1] - data[j, 1])^2
      z_weight[j, i] <- (data[i, 1] - data[j, 1])^2
    }
  }
  z_weight <- (0.01 + z_weight) / (0.01 + min(z_weight))
  
  test.match <- nonbimatch(distancematrix(test.mdm / z_weight))
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  # the next section of code checks the standardized covariate difference between
  # the low-dose and high-dose groups. If there is at least one coviarate that has
  # a standardized difference larger that cutoff, we will not return the matched
  # pairs generated.
  standard_differences <- 1:(ncol(data) - 2)
  sigmas <- 1:(ncol(data) - 2)
  
  for(i in 3:ncol(data)){
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
nbp_weighted_caliper <- function(data, p_matrix, delta){
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
  
  test.dist <- gendistance(data[, c(1, 3:ncol(data))], idcol = 1)
  test.mdm <- distancematrix(test.dist)
  
  # applying caliper to mahalnobis distance matrix
  
  z_weight <- matrix(1, ncol = nrow(data), nrow = nrow(data))
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      z_weight[i, j] <- (data[i, 1] - data[j, 1])^2
      z_weight[j, i] <- (data[i, 1] - data[j, 1])^2
    }
  }
  
  z_weight <- (0.01 + z_weight) / (0.01 + min(z_weight))
  
  caliper_dist_matrix <- distancematrix((test.mdm/z_weight) + caliper_matrix)
  test.match <- nonbimatch(caliper_dist_matrix)
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  # calculates standardized covariate differences, if we want to compare the 
  # standardized covariate difference between mahalanobis-only matching and
  # calipered-matching
  standard_differences <- 1:(ncol(data) - 3)
  sigmas <- 1:(ncol(data) - 3)
  
  for(i in 3:ncol(data)){
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


set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance <- data.frame(matrix(0, nrow = iter, ncol = 5))

while(successful_matching < iter){
  
  # generate data with specified sample size
  test <- generate_data_dose(N = sample_size)
  
  # create matched pairs based on mahalanobis-only matching
  matched_pairs <- nbp_weighted(test, 0.6)
  
  # nrow(matched_pairs$pairs) == 0 only if the set of matched pairs violate the
  # balancing criteria specified (here, we use a cutoff of 0.2)
  if(nrow(matched_pairs$pairs) > 0){
    
    # If our matched pairs are well balanced, we'll save this dataset for future
    # use, as well as the mahalanobis-only matched pairs we generated.
    data_list <- append(data_list, test)
    mahalanobis_pairs_list <- append(mahalanobis_pairs_list, matched_pairs$pairs)
    covariate_balance[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  #if(successful_matching %% 100 == 0){
  #  print(successful_matching)
  #}
  
  print(successful_matching)
  
}

write.csv(covariate_balance, "covariate_balance_sim1_weighted.csv")

# Parallelization setup
numcores <- detectCores()
output <- data.frame(matrix(nrow = 0, ncol = 8))
registerDoParallel(numcores - 6)

# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

output <- foreach(i = 1:iter, .combine = rbind,
                        .packages = c("nbpMatching",
                                      "RFCDE",
                                      "randomForest",
                                      "earth")) %dopar%
  {
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_linCDE <- propensity_matrix_linCDE(test)
    p_mat_RFCDE <- propensity_matrix_RFCDE(test)
    p_mat_model <- propensity_matrix_model(test)
    
    # creates matched pairs 
    mahalanobis_pairs <- as.data.frame(mahalanobis_pairs_list[((2*i) - 1):(2*i)])
    RFCDE_pairs <- nbp_weighted_caliper(test, p_mat_RFCDE, 0.1)$pairs
    linCDE_pairs <- nbp_weighted_caliper(test, p_mat_linCDE, 0.1)$pairs
    model_pairs <- nbp_weighted_caliper(test, p_mat_model, 0.1)$pairs
    
    ######################################################################################################################
    
    # estimators with mahalanobis only pairs, naive estimator
    lambda <- true_lambda(test, mahalanobis_pairs)
    e_naive <- naive_estimator(test, mahalanobis_pairs)
    var_naive <- variance_naive(test, mahalanobis_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "Mahalanobis"))
    
    
    # estimators with mahalanobis only pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "Mahalanobis"))
    
    # estimators with mahalanobis only pairs, linCDE estimator
    e_linCDE <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "Mahalanobis"))
    
    # estimators with mahalanobis only pairs, model estimator
    e_model <- bias_corrected_estimator(test, mahalanobis_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, mahalanobis_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "Mahalanobis"))    
    
    ######################################################################################################################
    
    # estimators with RFCDE pairs, naive estimator
    lambda <- true_lambda(test, RFCDE_pairs)
    e_naive <- naive_estimator(test, RFCDE_pairs)
    var_naive <- variance_naive(test, RFCDE_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "RFCDE_Caliper"))
    
    # estimators with RFCDE pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "RFCDE_Caliper"))
    
    # estimators with RFCDE pairs, RFCDE estimator
    e_linCDE <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, RFCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "RFCDE_Caliper"))
    
    # estimators with RFCDE pairs, linCDE estimator
    e_model <- bias_corrected_estimator(test, RFCDE_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, RFCDE_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "RFCDE_Caliper"))    
    
    ######################################################################################################################
    
    # estimators with linCDE pairs, naive estimator
    lambda <- true_lambda(test, linCDE_pairs)
    e_naive <- naive_estimator(test, linCDE_pairs)
    var_naive <- variance_naive(test, linCDE_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "linCDE_Caliper"))
    
    # estimators with linCDE pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "linCDE_Caliper"))
    
    # estimators with linCDE pairs, linCDE estimator
    e_linCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "linCDE_Caliper"))
    
    # estimators with linCDE pairs, model estimator
    e_model <- bias_corrected_estimator(test, linCDE_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, linCDE_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "linCDE_Caliper"))  
    
    ######################################################################################################################
    
    # estimators with model pairs, naive estimator
    lambda <- true_lambda(test, model_pairs)
    e_naive <- naive_estimator(test, model_pairs)
    var_naive <- variance_naive(test, model_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "model_Caliper"))
    
    # estimators with model pairs, RFCDE estimator
    e_RFCDE <- bias_corrected_estimator(test, model_pairs, p_mat_RFCDE, 0.1)
    var_RFCDE <- variance_bias_corrected(test, model_pairs, p_mat_RFCDE, e_RFCDE, 0.1)
    CI_RFCDE <- CI_maker(e_RFCDE, 0.05, var_RFCDE)
    coverage_RFCDE <- coverage_check(CI_RFCDE, lambda)
    MOE_RFCDE <- qnorm(0.975) * sqrt(var_RFCDE)
    temp_data <- rbind(temp_data, c(lambda, e_RFCDE, MOE_RFCDE, CI_RFCDE, coverage_RFCDE, "RFCDE", "model_Caliper"))
    
    # estimators with model pairs, linCDE estimator
    e_linCDE <- bias_corrected_estimator(test, model_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, model_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "model_Caliper"))
    
    # estimators with model pairs, model estimator
    e_model <- bias_corrected_estimator(test, model_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, model_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "model_Caliper"))     
    names(temp_data) <- c("lambda", "estimate", "MOE", "lower_CI", "upper_CI", "coverage", "estimation_method", "matching_method")
    
    return(temp_data)
    
  }

stopImplicitCluster()
output[, -c(7, 8)] <- lapply(output[, -c(7, 8)], as.numeric)

write.csv(output, "sim_weighted_output.csv")

matching_method_vector <- unique(output$matching_method)
estimation_method_vector <- unique(output$estimation_method)
summary_data <- data.frame(matrix(nrow = 0, ncol = 4))
for(i in 1:length(matching_method_vector)){
  
  method_data <- output[which(output$matching_method == matching_method_vector[i]), ]
  
  for(m in 1:length(estimation_method_vector)){
    
    subdata <- method_data[which(method_data$estimation_method == estimation_method_vector[m]), ]
    mean_abs_error <- mean(subdata[, 2] - subdata[, 1])
    root_MSE <- sqrt(mean((subdata[, 1] - subdata[, 2])^2))
    MOE_average <- mean(subdata[, 3])
    coverage_rate <- mean(subdata[, 6])
    
    summary_data <- rbind(summary_data, c(matching_method_vector[i], estimation_method_vector[m], 
                                          mean_abs_error, root_MSE, MOE_average, coverage_rate))
    
  } 
}

summary_data[, -c(1, 2)] <- lapply(summary_data[, -c(1, 2)], as.numeric)
summary_data[, 5] <- 2*summary_data[, 5]
names(summary_data) <- c("Matching_Method", "Estimation_Method","Mean_Error", "RMSE", "Average_Length", "Coverage_Rate")

write.csv(summary_data, "sim_weighted_summary.csv")


balance_tab <- read.csv("covariate_balance_sim1_weighted.csv")[, -1]
mean_balance <- data.frame(as.list(round_df(colMeans(balance_tab), 3)))
kbl(mean_balance, align = 'c', booktabs = T, format = "latex", row.names = FALSE,
    caption = "Balance Table, First Simulation", position = "H") %>%
  kable_classic(full_width = F, html_font = "Cambria")

#### Appendix (blank)
# Simulation E.2 Estimation and Inference (MAE, RMSE, ML, CR) --------------------

nbp_pn_weight <- function(data, p_matrix, M, cutoff){
  # dataframe data: data generated from generate_data_dose function
  # double cutoff: function returns NULL if the standardized difference in 
  # covariates is larger than the cutoff for at least one covariate
  
  # the next set of lines uses the nbpMatch library to create matched pairs
  # for data based on mahalanobis distance of covariates.
  
  pn_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      
      # if probability of observed treatment assignment is extreme, enforce
      # infinite penalty on distance matrix
      pn_matrix[i, j] <- M*(p_matrix[i, j] - 0.5)^2
      pn_matrix[j, i] <- pn_matrix[i, j]
      
    }
  }
  
  # making mahalanobis distance matrix
  
  test.dist <- gendistance(data[, c(1, 3:ncol(data))], idcol = 1)
  test.mdm <- distancematrix(test.dist)
  
  # applying caliper to mahalnobis distance matrix
  
  pn_dist_matrix <- distancematrix(test.mdm + pn_matrix)
  test.match <- nonbimatch(pn_dist_matrix)
  matches <- test.match$matches
  matches$group <- ifelse(as.numeric(matches$Group1.ID) < as.numeric(matches$Group2.ID), "L", "H")
  high_groups <- matches[which(matches$group == "H"), 2]
  low_groups <- matches[which(matches$group == "L"), 2]
  
  # the next section of code checks the standardized covariate difference between
  # the low-dose and high-dose groups. If there is at least one coviarate that has
  # a standardized difference larger that cutoff, we will not return the matched
  # pairs generated.
  standard_differences <- 1:(ncol(data) - 2)
  sigmas <- 1:(ncol(data) - 2)
  
  for(i in 3:ncol(data)){
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

set.seed(12345)
iter = 1000
sample_size <- 1000
data_list <- list()
mahalanobis_pairs_list <- list()
successful_matching <- 0
covariate_balance <- data.frame(matrix(0, nrow = iter, ncol = 5))

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
    #covariate_balance[successful_matching + 1, ] <- matched_pairs$standardized_diffs
    successful_matching <- successful_matching + 1
    
  }
  
  # Progress check - we'll be notified after every 100 datasets / matched pairs 
  # are successfully made (pass the balance criteria)
  if(successful_matching %% 100 == 0){
    print(successful_matching)
  }
  
  #print(successful_matching)
  
}

#write.csv(covariate_balance, "covariate_balance_sim1_weighted.csv")

# Parallelization setup
numcores <- detectCores()
output <- data.frame(matrix(nrow = 0, ncol = 8))
registerDoParallel(numcores - 6)

# now, we create matched pairs and estimate treatment effects with each matching
# and estimation method we're interested in.

output <- foreach(i = 1:iter, .combine = rbind,
                        .packages = c("nbpMatching",
                                      "randomForest",
                                      "earth")) %dopar%
  {
    # holds true lambda and estimates in dataframe
    temp_data <- as.data.frame(matrix(nrow = 0, ncol = 8))
    
    # extracts ith dataframe from the data list
    test <- as.data.frame(data_list[(7*(i-1) + 1):(7*i)])
    
    # creates p matrices, used for matching and estimation methods
    p_mat_linCDE <- propensity_matrix_linCDE(test)
    p_mat_model <- propensity_matrix_model(test)
    
    # creates matched pairs 
    linCDE_pairs <- nbp_pn_weight(test, p_mat_linCDE, 5, 5)$pairs
    model_pairs <- nbp_pn_weight(test,  p_mat_model, 5, 5)$pairs
    
    ######################################################################################################################
    
    # estimators with linCDE pairs, naive estimator
    lambda <- true_lambda(test, linCDE_pairs)
    e_naive <- naive_estimator(test, linCDE_pairs)
    var_naive <- variance_naive(test, linCDE_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "linCDE"))
    
    # estimators with linCDE pairs, linCDE estimator
    e_linCDE <- bias_corrected_estimator(test, linCDE_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, linCDE_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "linCDE"))
    
    # estimators with linCDE pairs, model estimator
    e_model <- bias_corrected_estimator(test, linCDE_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, linCDE_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "linCDE"))    
    
    ######################################################################################################################
    
    # estimators with model pairs, naive estimator
    lambda <- true_lambda(test, model_pairs)
    e_naive <- naive_estimator(test, model_pairs)
    var_naive <- variance_naive(test, model_pairs, e_naive)
    CI_naive <- CI_maker(e_naive, 0.05, var_naive)
    coverage_naive <- coverage_check(CI_naive, lambda)
    MOE_naive <- qnorm(0.975) * sqrt(var_naive)
    temp_data <- rbind(temp_data, c(lambda, e_naive, MOE_naive, CI_naive, coverage_naive, "Naive", "model"))
    
    # estimators with modelpairs, linCDE estimator
    e_linCDE <- bias_corrected_estimator(test, model_pairs, p_mat_linCDE, 0.1)
    var_linCDE <- variance_bias_corrected(test, model_pairs, p_mat_linCDE, e_linCDE, 0.1)
    CI_linCDE <- CI_maker(e_linCDE, 0.05, var_linCDE)
    coverage_linCDE <- coverage_check(CI_linCDE, lambda)
    MOE_linCDE <- qnorm(0.975) * sqrt(var_linCDE)
    temp_data <- rbind(temp_data, c(lambda, e_linCDE, MOE_linCDE, CI_linCDE, coverage_linCDE, "linCDE", "model"))
    
    # estimators with model pairs, model estimator
    e_model <- bias_corrected_estimator(test, model_pairs, p_mat_model, 0.1)
    var_model <- variance_bias_corrected(test, model_pairs, p_mat_model, e_model, 0.1)
    CI_model <- CI_maker(e_model, 0.05, var_model)
    coverage_model <- coverage_check(CI_model, lambda)
    MOE_model <- qnorm(0.975) * sqrt(var_model)
    temp_data <- rbind(temp_data, c(lambda, e_model, MOE_model, CI_model, coverage_model, "model", "model"))    
 
    names(temp_data) <- c("lambda", "estimate", "MOE", "lower_CI", "upper_CI", "coverage", "estimation_method", "matching_method")
    
    return(temp_data)
    
  }

stopImplicitCluster()
output[, -c(7, 8)] <- lapply(output[, -c(7, 8)], as.numeric)

write.csv(output, "sim_pn_weight_output.csv")

matching_method_vector <- unique(output$matching_method)
estimation_method_vector <- unique(output$estimation_method)
summary_data <- data.frame(matrix(nrow = 0, ncol = 4))
for(i in 1:length(matching_method_vector)){
  
  method_data <- output[which(output$matching_method == matching_method_vector[i]), ]
  
  for(m in 1:length(estimation_method_vector)){
    
    subdata <- method_data[which(method_data$estimation_method == estimation_method_vector[m]), ]
    mean_abs_error <- mean(subdata[, 2] - subdata[, 1])
    root_MSE <- sqrt(mean((subdata[, 1] - subdata[, 2])^2))
    MOE_average <- mean(subdata[, 3])
    coverage_rate <- mean(subdata[, 6])
    
    summary_data <- rbind(summary_data, c(matching_method_vector[i], estimation_method_vector[m], 
                                          mean_abs_error, root_MSE, MOE_average, coverage_rate))
    
  } 
}

summary_data[, -c(1, 2)] <- lapply(summary_data[, -c(1, 2)], as.numeric)
summary_data[, 5] <- 2*summary_data[, 5]
names(summary_data) <- c("Matching_Method", "Estimation_Method","Mean_Error", "RMSE", "Average_Length", "Coverage_Rate")

write.csv(summary_data, "sim_pn_weight_summary.csv")

