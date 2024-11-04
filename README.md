# Bias Mitigation in Matched Observational Studies with Continuous Treatments

## Introduction

The code is for reproducing the main and supplementary simulations for the paper:

* Bias Mitigation in Matched Observational Studies with Continuous Treatments: Calipered Non-Bipartite Matching and Bias-Corrected Estimation and Inference
<br /><i>Anthony Frazier, Siyu Heng, and Wen Zhou</i><br>

## Guidelines for Replication

### Results in Section 1 (Introduction: Background, Motivating Example, and Main Contributions)

Download the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the `Libraries, Data Files, Functions` code chunk to load the necessary libraries and functions, then run the `Balance Table for Original Set of 1,211 Matched Pairs` code chunk to produce figure 1. 

### Results in Section 5 (Simulation Studies)
  
Download the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, load the necessary libraries and functions, then run the simulations in the Main Simulation section. Specifically, run the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` and the `Simulation 2 Estimation and Inference (MAE, RMSE, ML, CR)`chunks of code to replicate results shown in tables 1 and 2 of the main manuscript. You may run the `Simulation Result Tables` chunk of code to produce a LaTeX table showing the MAE, RMSE, ML and CR results combined. 

### Results in Section 6 (Data Application: Re-Analyzing the Effects of Social Distancing Practice On COVID-19 Case Counts)

1) **Figure 2:** In the `Social Distancing Data Analysis Code.R` file, run the `Creating p-matrix (linCDE)` and the `Plots and Summary Statistics of Probability of Observed Treatment Assignment` code chunks to produce the two plots shown in figure 2.

2) **Table 3:** In the `Social Distancing Data Analysis Code.R` file, run the `Creating p-matrix (linCDE)` and the `Average Treatment Effect Estimation for Original Set of Matched Pairs` code chunks to replicate the ATE results pertaining to the non-calipered set of matched pairs. Run the `Average Treatment Effect Estimation for Calipered (xi=0.1) Set of Matched Pairs` code chunk to replicate the ATE results pertaining to the calipered set of matched pairs. These code chunks also provide the results shown in Appendix C, table S3. 

### Results in Web Appendix B (Extra Simulations)

1) **Section B.1 and B.2:** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, after running the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` code chunk, run the `Simulation 1 P Distribution (supplementary materials)` code chunk to produce the .csv files containing the distribution of treatment assignment probabilities for simulation 1. To replicate the figures in section B.1, run the `P Distribution Plots` code chunk, which will create a histogram for whichever .csv file you chose to read into the plot object. Repeat these steps for the `Simulation 2 Estimation and Inference (MAE, RMSE, ML, CR)` and `Simulation 2 P Distribution (supplementary materials)` code chunks to replicate the results shown in figure S2

2) **Section B.3:** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, after running the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` and the `Simulation 2 Estimation and Inference (MAE, RMSE, ML, CR)` code chunks, run the `Covariate Balance Tables` code chunk to create the covariate balance tables for simulations 1 and 2 separately.

### Results in Web Appendix C (Real Data Application)

1) **Section C.2:** In the `Social Distancing Data Analysis Code.R` file, run the ggplot code in the `0.1 caliper analysis` code chunk to replicate figure 3. 

2) **Section C.3:** The `balance_table.csv` file contains the balance table information shown in section C.3.

3) **Section C.4:** In the `Social Distancing Data Analysis Code.R` file, run the `creating p-matrix (linCDE)` and the `Urban vs Rural Subgroup Analysis` code chunks to replicate table 1 and both plots shown in figure 4.

4) **Section C.5:** In the `Social Distancing Data Analysis Code.R` file, run the `creating p-matrix (linCDE)` and the `Subgroup analysis (by quantile of dose differences)` code chunks to replicate table 3

5) **Section C.6:** All sensitivity analysis results can be replicated in the `Sensitivity Analysis` section of the `Social Distancing Data Analysis Code.R` file.

## License
Copyright 2024 Anthony Frazier, Siyu Heng, and Wen Zhou

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
