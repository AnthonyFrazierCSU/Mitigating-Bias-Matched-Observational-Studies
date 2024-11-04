################################################################################
####           Bias Mitigation in Matched Observational Studies             ####
####                    with Continuous Treatments                          ####
####                        Application Code                                ####
####                Anthony Frazier, Siyu Heng, Wen Zhou                    ####
################################################################################

#### libraries, data files, functions --------------------------------------------------
library(ggplot2)                                    # plot creation
library(tidyverse)                                  # dataframe manipulation
library(randomForest)                               # used in linCDE, if using RF for centering
library(nbpMatching)                                # non-bipartite matching algorithm
library(kableExtra)                                 # table creation
library(tidyr)                                      # dataframe manipulation
library(dplyr)                                      # dataframe manipulation
devtools::install_github("ZijunGao/LinCDE", build_vignettes = TRUE) # linCDE
set.seed(12345)                                     # linCDE reproducibility

# data files
lowdose_data <- read.csv("matched_data_lowdose.csv")
highdose_data <- read.csv("matched_data_highdose.csv")
lowdose_0.1caliper_data <- read.csv("matched_data_0.1caliper_lowdose.csv")
highdose_0.1caliper_data <- read.csv("matched_data_0.1caliper_highdose.csv")
caliper_p_dist <- read.csv("full_p_distribution_0.1caliper.csv")[, -1]

# reformatting / cleaning data
names(lowdose_data)[1] <- "ID"
names(highdose_data)[1] <- "ID"
all_postmatch_data <- rbind(lowdose_data, highdose_data)
names(lowdose_0.1caliper_data)[1] <- "ID"
names(highdose_0.1caliper_data)[1] <- "ID"
all_postmatch_0.1caliper_data <- rbind(lowdose_0.1caliper_data, highdose_0.1caliper_data)

# functions
round_df <- function(x, digits) {
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
quantile_analysis <- function(high_df, low_df, quants, pvec, xi){
  
  dose_diffs <- high_df$cumulative_dose - low_df$cumulative_dose
  output <- data.frame(matrix(nrow = 0, ncol = 5))
  for(i in 1:length(quants)){
    subindex <- which(dose_diffs > quants[i])
    high_df_sub <- high_df[subindex, ]
    low_df_sub <- low_df[subindex, ]
    p_sub <- pvec[subindex]
    quant <- as.numeric(quants[i])
    
    numerator <- 1:nrow(high_df_sub)
    denominator <- 1:nrow(low_df_sub)
    for(i in 1:nrow(high_df_sub)){
      numerator[i] <- (high_df_sub[i, "cases"] - low_df_sub[i, "cases"])
      denominator[i] <- (high_df_sub[i, "cumulative_dose"] - low_df_sub[i, "cumulative_dose"])
    }
    ATE_classic_cases <- sum(numerator)/sum(denominator)
    
    # conservative confidence interval for classic estimator
    
    a <- ATE_classic_cases * (sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])) / nrow(high_df_sub)
    numerator <- sum(((high_df_sub[, "cases"] - low_df_sub[, "cases"]) - a)^2)
    denominator <- sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2
    variance_ATE_classic_cases <- (numerator/denominator) * (nrow(high_df_sub)/(nrow(high_df_sub) - 1))
    CI_ATE_classic_cases <- c(ATE_classic_cases - (qnorm(0.975) * sqrt(variance_ATE_classic_cases)), ATE_classic_cases + (qnorm(0.975) * sqrt(variance_ATE_classic_cases)))
    
    output <- rbind(output, c(quant, "Classic", ATE_classic_cases, CI_ATE_classic_cases))
    
    # classic Neyman estimator with covariate-adjusted variance estimator
    
    covariates <- c(6:8, 10:15, 20:46, 54:55)
    Q <- as.matrix((high_df_sub[, covariates] - low_df_sub[, covariates])/2)
    Q <- Q[, -ncol(Q)]
    hatQ <- Q %*% solve(base::t(Q)%*%Q) %*% base::t(Q)
    y <- 2*(high_df_sub[, "cases"] - low_df_sub[, "cases"])/sqrt(1-diag(hatQ))
    cov_adj_var <- (y %*% (diag(1, nrow(hatQ)) - hatQ) %*% y)/(4 * (sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2))
    CI_classic_cov_adj <- c(ATE_classic_cases - (qnorm(0.975) * sqrt(cov_adj_var)), ATE_classic_cases + (qnorm(0.975) * sqrt(cov_adj_var)))
    
    output <- rbind(output, c(quant, "Classic_covariate_adjusted", ATE_classic_cases, CI_classic_cov_adj))
    
    # Bias-corrected estimation, not regulated
    numerator <- 1:nrow(high_df_sub)
    denominator <- 1:nrow(high_df_sub)
    for(i in 1:nrow(high_df_sub)){
      numerator[i] <- p_sub[i]^(-1) * (high_df_sub[i, "cases"] - low_df_sub[i, "cases"])
      denominator[i] <- 2*(high_df_sub[i, "cumulative_dose"] - low_df_sub[i, "cumulative_dose"])
    }
    ATE_bias_corrected_cases <- sum(numerator)/sum(denominator)
    
    # Conservative confidence interval for bias-corrected estimator, not regulated
    
    a <- ATE_bias_corrected_cases * ((2*(high_df_sub[i, "cumulative_dose"] - low_df_sub[i, "cumulative_dose"]))/(nrow(high_df_sub)))
    
    numerator <- 1:nrow(high_df_sub)
    for(i in 1:nrow(high_df_sub)){
      numerator[i] <- ((p_sub[i])^(-1)*(high_df_sub[i, "cases"] - low_df_sub[i, "cases"]) - a)^2
    }
    denominator <- 4 * (sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2)
    variance_ATE_bias_corrected_cases <- (sum(numerator)/sum(denominator)) * (nrow(high_df_sub)/(nrow(high_df_sub) - 1))
    CI_ATE_bias_corrected_cases <- c(ATE_bias_corrected_cases - (qnorm(0.975) * sqrt(variance_ATE_bias_corrected_cases)), ATE_bias_corrected_cases + (qnorm(0.975) * sqrt(variance_ATE_bias_corrected_cases)))
    
    output <- rbind(output, c(quant, "bias_Corrected", ATE_bias_corrected_cases, CI_ATE_bias_corrected_cases))
    
    # bias-corrected estimator (non-regulated) with covariate-adjusted variance estimator
    covariates <- c(6:8, 10:15, 20:46, 54:55)
    Q <- as.matrix((high_df_sub[, covariates] - low_df_sub[, covariates])/2)
    Q <- Q[, -ncol(Q)]
    hatQ <- Q %*% solve(base::t(Q)%*%Q) %*% base::t(Q)
    bc_y <- (p_sub^(-1))*(high_df_sub[, "cases"] - low_df_sub[, "cases"])/sqrt(1-diag(hatQ))
    cov_adj_var <- (bc_y %*% (diag(1, nrow(hatQ)) - hatQ) %*% bc_y)/(4 * (sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2))
    CI_bc_cov_adj <- c(ATE_bias_corrected_cases - (qnorm(0.975) * sqrt(cov_adj_var)), ATE_bias_corrected_cases + (qnorm(0.975) * sqrt(cov_adj_var)))
    
    output <- rbind(output, c(quant, "bc_covariate_adjusted", ATE_bias_corrected_cases, CI_bc_cov_adj))
    
    # Bias-corrected estimation, regulated
    
    numerator <- 1:nrow(high_df_sub)
    denominator <- 1:nrow(high_df_sub)
    for(i in 1:nrow(high_df_sub)){
      
      p <- p_sub[i]
      if(p < xi){
        p <- xi
      }
      if(p > 1-xi){
        p <- 1-xi
      }
      
      numerator[i] <- p^(-1) * (high_df_sub[i, "cases"] - low_df_sub[i, "cases"])
      denominator[i] <- 2*(high_df_sub[i, "cumulative_dose"] - low_df_sub[i, "cumulative_dose"])
    }
    ATE_bcReg_cases <- sum(numerator)/sum(denominator)
    
    # Conservative confidence interval for bias-corrected estimator, regulated (xi = 0.1)
    
    a <- ATE_bcReg_cases * ((2*(high_df_sub[i, "cumulative_dose"] - low_df_sub[i, "cumulative_dose"]))/(nrow(high_df_sub)))
    low_num <- rep(0, nrow(high_df_sub)); mid_num <- rep(0, nrow(high_df_sub)); high_num <- rep(0, nrow(high_df_sub))
    
    for(i in 1:nrow(high_df_sub)){
      
      p = p_sub[i]
      
      if(p < xi){
        low_num[i] <- ((xi)^(-1)*(high_df_sub[i, "cases"] - low_df_sub[i, "cases"]) - a)^2
      }
      
      if(p <= (1-xi) & p >= xi){
        mid_num[i] <- ((p)^(-1)*(high_df_sub[i, "cases"] - low_df_sub[i, "cases"]) - a)^2
      }
      
      if(p > (1-xi)){
        high_num[i] <-((1 - xi)^(-1)*(high_df_sub[i, "cases"] - low_df_sub[i, "cases"]) - a)^2
      }
      
    }
    
    denominator <- 4*sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2
    variance_ATE_bcReg_cases <- ((sum(low_num) + sum(mid_num) + sum(high_num)) / denominator) * (nrow(high_df_sub)/(nrow(high_df_sub) - 1))
    CI_ATE_bcReg_cases <- c(ATE_bcReg_cases - (qnorm(0.975) * sqrt(variance_ATE_bcReg_cases)), ATE_bcReg_cases + (qnorm(0.975) * sqrt(variance_ATE_bcReg_cases)))
    
    output <- rbind(output, c(quant, "bcReg", ATE_bcReg_cases, CI_ATE_bcReg_cases))
    
    # bias-corrected estimator (regulated) with covariate-adjusted variance estimator
    psub_reg <- p_sub
    psub_reg[p_sub < xi] <- xi
    psub_reg[psub_reg > (1-xi)] <- 1-xi
    bcReg_y <- (psub_reg^(-1))*(high_df_sub[, "cases"] - low_df_sub[, "cases"])/sqrt(1-diag(hatQ))
    cov_adj_var <- (bcReg_y %*% (diag(1, nrow(hatQ)) - hatQ) %*% bcReg_y)/(4 * (sum(high_df_sub[, "cumulative_dose"] - low_df_sub[, "cumulative_dose"])^2))
    CI_bcReg_cov_adj <- c(ATE_bcReg_cases - (qnorm(0.975) * sqrt(cov_adj_var)), ATE_bcReg_cases + (qnorm(0.975) * sqrt(cov_adj_var)))
    
    output <- rbind(output, c(quant, "bcReg_covariate_adjusted", ATE_bcReg_cases, CI_bcReg_cov_adj))
    
  }
  
  names(output) <- c("Quantile", "Method", "Estimate", "Lower_CI", "Upper_CI")
  
  output[, c(1, 3:5)] <- lapply(output[, c(1, 3:5)], as.numeric)
  return(output)
  
}

#### creating p-matrix (linCDE) ------------------------------------------------

covariates <- c(6:8, 10:15, 20:46, 54:55)

# estimating conditional mean of treatment given covariates, fitting linCDE on residuals
model <- LinCDE::LinCDE.boost(y = all_postmatch_data[, 50], X = all_postmatch_data[, covariates], terminalSize = 20, verbose = FALSE, centering = TRUE, centeringMethod = "linearRegression")

# will hold conditional density estimates for every combination of treatment
# and covariates
propensity_matrix <- matrix(0, nrow = nrow(all_postmatch_data), ncol = nrow(all_postmatch_data))

for(i in 1:nrow(propensity_matrix)){
  propensity_matrix[i, ] <- LinCDE::predict.LinCDE(model, X = all_postmatch_data[i, covariates], y = all_postmatch_data[, 50])
}

# will hold output (matrix of probabilities of observed treatment assignment)
p_matrix <- matrix(0, nrow = nrow(all_postmatch_data), ncol = nrow(all_postmatch_data))

for(i in 1:(nrow(all_postmatch_data)-1)){
  for(j in (i+1):nrow(all_postmatch_data)){
    observed_treatment <- propensity_matrix[i, i] * propensity_matrix[j, j]
    opposing_treatment <- propensity_matrix[j, i] * propensity_matrix[i, j]
    p_matrix[i, j] <- observed_treatment / (observed_treatment + opposing_treatment)
    p_matrix[j, i] <- p_matrix[i, j]
  }
}

p_matched <- 1:nrow(highdose_data)
for(i in 1:nrow(highdose_data)){
  p_matched[i] <- p_matrix[i, i+1211]
}

#### Plots and summary statistics of probability of observed treatment assignment ----

# histogram showing distribution of probability of observed treatment assignment for all possible pairs
ggplot() + 
  geom_histogram(aes(x = p_matrix[upper.tri(p_matrix)]), binwidth = 0.05, fill = "white", color = "black") + 
  theme_bw() + 
  labs(x = "Probability of Observed Treatment Assignment", y = "Frequency") + 
  theme(axis.text = element_text(size = 12))
ggsave("P_distribution_all.png", width = 4.5, height = 3)

# histogram showing distribution of probability of observed treatment assignment of matched pairs
ggplot() + 
  geom_histogram(aes(x = p_matched), binwidth = 0.05, fill = "white", color = "black") + 
  theme_bw() + 
  labs(x = "Probability of Observed Treatment Assignment", y = "Frequency") + 
  theme(axis.text = element_text(size = 12))
ggsave("P_distribution_matched.png", width = 4.5, height = 3)

# 5-number summary of probability of observed treatment assignment for matched pairs
quantile(p_matched, c(0, 0.25, 0.5, 0.75, 1.0))

#### Average treatment effect estimation for original set of matched pairs, Cases -----

main_analysis_data <- quantile_analysis(highdose_data, lowdose_data, c(0), p_matched, 0.1)
main_analysis_data <- main_analysis_data[, -1]
write.csv(main_analysis_data, "main_analysis_summary.csv")

#### 0.1 caliper analysis ----

# histogram showing distribution of probability of observed treatment assignment of matched pairs
ggplot() + 
  geom_histogram(aes(x = caliper_p_dist), binwidth = 0.05, fill = "white", color = "black") + 
  theme_bw() + 
  labs(x = "Probability of Observed Treatment Assignment", y = "Frequency") + 
  theme(axis.text = element_text(size = 12))
ggsave("P_distribution_0.1caliper.png", width = 4.5, height = 3)

# 5-number summary of probability of observed treatment assignment for matched pairs
quantile(caliper_p_dist, c(0, 0.25, 0.5, 0.75, 1.0))

caliper_analysis <- quantile_analysis(highdose_0.1caliper_data, lowdose_0.1caliper_data, c(0), caliper_p_dist, 0.1)
caliper_analysis <- caliper_analysis[, -1]
#caliper_analysis <- caliper_analysis[-3, ]
write.csv(caliper_analysis, "main_caliper_analysis.csv")


#### Balance table for original set of 1,211 matched pairs ---------------------

bt <- read.csv("balance_table.csv")
bt <- bt[-((nrow(bt)-1):(nrow(bt))), c(1, 4)]
names(bt) <- c("Variable", "standardized_difference")
bt$Variable <- c("Female", "Below 18", "Above 65", "Black", "Hispanic", "Black/Hispanic",
                 "Driving Alone to Work", "Smoking", "Flu Vaccine", "Some College",
                 "Social Association", "Traffic", "Non Metro", "Percent Poverty", 
                 "Population Density", "Population", "Cases 1 Week Before", 
                 "Cases 2 Weeks Before", "Cases 3 Weeks Before",
                 "Cases 4 Weeks Before", "Cases 5 Weeks Before",
                 "Cases 6 Weeks Before", "Cases 7 Weeks Before",
                 "Cases 8 Weeks Before", "Cases 9 Weeks Before",
                 "Cases 10 Weeks Before", "Cases 11 Weeks Before",
                 "Deaths 1 Week Before", "Deaths 2 Weeks Before", "Deaths 3 Weeks Before",
                 "Deaths 4 Week Before", "Deaths 5 Weeks Before", "Deaths 6 Weeks Before",
                 "Deaths 7 Week Before", "Deaths 8 Weeks Before", "Deaths 9 Weeks Before",
                 "Deaths 10 Week Before", "Deaths 11 Weeks Before")
balance_df <- cbind(bt$Variable, abs(bt$standardized_difference))
balance_df <- data.frame(balance_df)
names(balance_df) <- c("Variable", "Standardized_difference")
balance_df[, 2] <- as.numeric(balance_df[, 2])
balance_df$Variable <- factor(balance_df$Variable, labels = rev(bt$Variable))
#top10_balance <- balance_df[1:10, ]

ggplot(data = balance_df) + 
  geom_point(aes(x = Variable, y = Standardized_difference)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  geom_hline(yintercept = 0.2, linetype = "dashed") + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Covariates", 
       y = "Standardized Difference-in-means") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))
ggsave("Covariate Balance for Original Matched Pairs.png", height = 4.5, width = 5.5)


