## Mediation_Surrogacy

Here I share the code supporting our work on methods for counterfactual mediation analysis with multistate model for surrogate and clinical time-to-event outcomes.

#### Abstract
> We introduce a counterfactual-based mediation analysis for surrogacy evaluation with time-to-event surrogate and clinical outcomes. Our approach accommodates censoring and competing risks. We use a multistate model for risk prediction to account for both transitions towards the clinical outcome and transitions through the surrogate outcome. We use the counterfactual framework to define the natural direct and indirect effects with a causal interpretation. Based on these measures, we define the proportion of the treatment effect on the clinical outcome mediated by the surrogate outcome. We estimate the proportion for both the cumulative risk and restricted mean time lost. We illustrate our approach in a simulation study and using 18-year follow-up data from the SPCG-4 randomized trial of radical prostatectomy for prostate cancer. We assess time to metastasis as a surrogate outcome for prostate cancer-specific mortality.

#### Code and data sharing
For reproducbility and implementation, I share the code for all analyses presented in our working paper. The data from the SPCG-4 randomized trial are not publically available. The simulated dataset and demo may help you to understand the methods. 

#### 1. R functions to compute mediation measures according to our methods
* [FUN_getMedMeas_risk.R](FUN_getMedMeas_risk.R)
* [FUN_getMedMeas_rmtl.R](FUN_getMedMeas_rmtl.R)
* [FUN_AUCiw.R](FUN_AUCiw.R)
* [FUN_getMedMeas_Cox.R](FUN_getMedMeas_Cox.R) - get mediation measures accoridng to binary covariate for transition to state 1 (surrogate)

#### 2. R code used for application of method in SPCG-4 Randomized Controlled Trial 
* [SPCGanalysis_illnessdeath.R](SPCGanalysis_illnessdeath.R)
* [SPCGanalysis_msm.R](SPCGanalysis_msm.R)

#### 3. Simulated data and R code implementing method
* [demodata.txt](demodata.txt)
* [analysis_demo.R](analysis_demo.R)

#### 4. Code to replicate simulation study 
* [GenerateData.do](GenerateData.do) - a stata program to generate individual level time-to-event data through illness death model
* [scenario_table.txt](scenarioTable.txt) - table of scenario parameters
* [simulation_analysis.R](simulation_analysis.R) - 

