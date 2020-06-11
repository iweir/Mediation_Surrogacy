## Mediation_Surrogacy

Here I share the code supporting our work on methods for counterfactual mediation analysis with multistate model for surrogate and clinical time-to-event outcomes.

#### Abstract
> We introduce a counterfactual-based mediation analysis for surrogacy evaluation with time-to-event surrogate and clinical outcomes. Our approach accommodates censoring and competing risks. We use a multistate model for risk prediction to account for both transitions towards the clinical outcome and transitions through the surrogate outcome. We use the counterfactual framework to define the natural direct and indirect effects with a causal interpretation. Based on these measures, we define the proportion of the treatment effect on the clinical outcome mediated by the surrogate outcome. We estimate the proportion for both the cumulative risk and restricted mean time lost. We illustrate our approach in a simulation study and using 18-year follow-up data from the SPCG-4 randomized trial of radical prostatectomy for prostate cancer. We assess time to metastasis as a surrogate outcome for prostate cancer-specific mortality.

#### Code and data sharing
For reproducbility and transparency, I share the code for all analyses presented in my working paper. The data from the SPCG-4 randomized trial are not publically available however data generated from the simulation study may help you to understand the methods. 

#### 1. R function to compute mediation measures according to our methods
* XX FUN_getMedMeas_risk.R - with option to recieve CI
* XX FUN_getMedMeas_rmtl.R
* XX FUN_auc.R

#### 2. R code used for application of method in SPCG-4 Randomized Controlled Trial 
* SPCGanalysis_illnessdeath.R
* SPCGanalysis_msm.R

#### 3. Code to replicate simulation study 
* XX GenerateData.do - a stata program to generate individual level time-to-event data through illness death model
* XX scenario_table.txt 
* XX getMedMeas_Cox.R - get mediation measures accoridng to binary covariate for transition to state 1 (surrogate)
* XX simulation_analysis.R

