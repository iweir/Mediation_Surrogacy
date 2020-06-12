## Mediation_Surrogacy

Here I share the code supporting our work on methods for counterfactual mediation analysis with multistate model for surrogate and clinical time-to-event outcomes.

#### Abstract
> We introduce a counterfactual-based mediation analysis for surrogacy evaluation with time-to-event surrogate and clinical outcomes. Our approach accommodates censoring and competing risks. We use a multistate model for risk prediction to account for both transitions towards the clinical outcome and transitions through the surrogate outcome. We use the counterfactual framework to define the natural direct and indirect effects with a causal interpretation. Based on these measures, we define the proportion of the treatment effect on the clinical outcome mediated by the surrogate outcome. We estimate the proportion for both the cumulative risk and restricted mean time lost. We illustrate our approach in a simulation study and using 18-year follow-up data from the SPCG-4 randomized trial of radical prostatectomy for prostate cancer. We assess time to metastasis as a surrogate outcome for prostate cancer-specific mortality.

#### Code and data sharing
For reproducbility and implementation, I share the code for all analyses presented in our working paper. The data from the SPCG-4 randomized trial are not publically available. The simulated dataset and demo may help you to understand the methods. 

#### 1. R functions to compute mediation measures according to our methods
* [FUN_getMedMeas_risk.R](FUN_getMedMeas_risk.R) - to get mediation measures for multistate models according to difference in cumulative risks
* [FUN_getMedMeas_rmtl.R](FUN_getMedMeas_rmtl.R) - to get mediation measures for multistate models according to difference in rmtl
* [FUN_AUCiw.R](FUN_AUCiw.R) - to compute the area under the curve (rmtl)
* [FUN_getMedMeas_cox.R](FUN_getMedMeas_cox.R) - to get mediation measures from Cox model (used in simulation study)

#### 2. Simulated data and R code implementing method
* [demodata.txt](demodata.txt) - a simulated dataset for use demonstrating methods 
* [analysis_demo.R](analysis_demo.R) - demonstration of methods using demodata.txt

#### 3. Code to replicate simulation study 
* [GenerateData.do](GenerateData.do) - a stata program to generate individual level time-to-event data through illness death model
* [scenarioTable.txt](scenarioTable.txt) - table of scenario parameters
* [simulation_analysis.R](simulation_analysis.R) - code to reproduce simulation study

