# Bias Mitigation in Matched Observational Studies with Continuous Treatments

## Introduction

The code is for reproducing the main and supplementary simulations for the paper:

* Bias Mitigation in Matched Observational Studies with Continuous Treatments: Calipered Non-Bipartite Matching and Bias-Corrected Estimation and Inference
<br /><i>Anthony Frazier, Siyu Heng, and Wen Zhou</i><br>

## Guidelines for Replication

### Results in Section 1 (Introduction)

Download the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the `Libraries, Data Files, Functions` code chunk to load the necessary libraries and functions, then run the `Balance Table for Original Set of 1,211 Matched Pairs` code chunk to produce figure 1. 

### Results in Section 4 (Simulation Studies)
  
Download the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, load the necessary libraries and functions, then run the simulations in the Main Simulation section. Specifically, run the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` chunk of code to replicate results shown in table 2 of the main manuscript. You may run the `Simulation Result Tables` chunk of code to produce a LaTeX table showing the Bias, MCIL and CR results combined. 

### Results in Section 5 (Data Application: Re-Analyzing the Effects of Social Distancing Practice On COVID-19 Case Counts)

1) **Figure 2:** In the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the `Creating p-matrix (linCDE)` and the `Plots and Summary Statistics of Probability of Observed Treatment Assignment` code chunks to produce the two plots shown in figure 2.

2) **Table 3:** In the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the `Creating p-matrix (linCDE)` and the `Average Treatment Effect Estimation for Original Set of Matched Pairs` code chunks to replicate the ATE results pertaining to the non-calipered set of matched pairs. Run the `Average Treatment Effect Estimation for Calipered (xi=0.075) Set of Matched Pairs` code chunk to replicate the ATE results pertaining to the calipered set of matched pairs. 

### Results in Web Appendix C: More Details on the Simulations

**Table S1** Download the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, load the necessary libraries and functions, then run the `Appendix B` code chunk.

### Results in Web Appendix C: More Details on the Simulations

1) **Table S2** The results shown in Table S2 are attained by following the instructions for replication the results in Section 4 of the main manuscript.

2) **Table S3** Download the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, load the necessary libraries and functions, then run the simulations in the Main Simulation section. Specifically, run the `Simulation 2 Estimation and Inference (Supplementary Materials)` chunk of code to replicate results shown in table S3 of the supplementary materials. You may run the `Simulation Result Tables` chunk of code to produce a LaTeX table showing the Bias, MCIL and CR results combined.

3) **Table S4** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, after running the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` and the `Simulation 2 Estimation and Inference (MAE, RMSE, ML, CR)` code chunks, run the `Covariate Balance Tables` code chunk to create the covariate balance tables for simulations 1 and 2 separately.

1) **Figures S1 and S2** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R` file, after running the `Simulation 1 Estimation and Inference (MAE, RMSE, ML, CR)` code chunk, run the `Simulation 1 P Distribution (supplementary materials)` code chunk to produce the .csv files containing the distribution of treatment assignment probabilities for simulation 1. To replicate figure S1, run the `P Distribution Plots` code chunk, which will create a histogram for whichever .csv file you chose to read into the plot object. Repeat these steps for the `Simulation 2 Estimation and Inference (Supplementary Materials)` and `Simulation 2 P Distribution (supplementary materials)` code chunks to replicate the results shown in figure S2.

### Results in Web Appendix D: More Details on the Data Application

1) **Tables S5, S6 and S7** In the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the `Balance Table for Original Set of 1,211 Matched Pairs` code chunk to replicates the covariate balance table before non-bipartite matching (Table S5) and the covariate balance table with conventional non-bipartite matching (Table S6). Run the ``Balance Table, Calipered (xi = 0.075) Versus Original Set`` to replicate the covariate balance table for the matched pairs formed by dose assignment discrepancy caliper (Table S7).

2) **Figure S3** In the `Bias Mitigation in Matched Observational Studies Application Code.R` file, run the ggplot code in the `Average Treatment Effect Estimation for Calipered (xi=0.075) Set of Matched Pairs` code chunk to replicate figure S3.

3) **Figure S4** In the `Bias Mitigation in Matched Observational Studies Application Code.R` file, running the code chunk ``Balance Table, Calipered (xi = 0.075) Versus Original Set`` also allows you to replicate figure S4.

### Results in Web Appendix E: Effectiveness of the Bias Mitigation Framework Under Alternative Matching Schemes

1) **Table S8** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R`, run the `Simulation E.1 Estimation and Inference (MAE, RMSE, ML, CR)` code chunk to replicate the results shown in Appendix E.1 of the supplementary materials.

2) **Table S9** In the `Bias Reduction for Matching with Continuous Treatments Simulation code.R`, run the `Simulation E.2 Estimation and Inference (MAE, RMSE, ML, CR)` code chunk to replicate the results shown in Appendix E.2 of the supplementary materials.

### Results in Web Appendix G: More Details Regarding the Dose Assignment Discrepancy Caliper

**Table S10** the covariate balance tables for all 8 matching schemes considered in Appendix G.3 are contained in the `Balance Tables for Appendix G.3` folder. Table S10 simply compiles the standardized difference-in-means from each of these 8 covariate balance tables.

## License
Copyright 2025 Anthony Frazier, Siyu Heng, and Wen Zhou

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
