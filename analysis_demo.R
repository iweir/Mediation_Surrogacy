# step 1 - prepare work space
  path <- "C:\\Users\\iweir\\Dropbox (Harvard University)\\SPCG-4\\"
  
  library(mstate)
  library(survival)
  
  source(paste0(path, "4. Github\\FUN_AUCiw.R"))
  source(paste0(path, "4. Github\\FUN_getMedMeas_Cox.R"))
  source(paste0(path, "4. Github\\FUN_getMedMeas_risk.R"))
  source(paste0(path, "4. Github\\FUN_getMedMeas_rmtl.R"))
  
  demodata <- read.table(paste0(path, "4. Github\\simulated_data.txt"), header=TRUE, sep="\t")

# step 2 - prepare data with mstate structure
  trans	<- trans.illdeath()
  
  tt		<- matrix(c(rep(NA,nrow(demodata)), 
                  demodata$timeM,demodata$timeD), 		
                nrow(demodata), 3)
  
  st		<- matrix(c(rep(NA,nrow(demodata)), 
                  demodata$MetaInd,demodata$DeathInd), 	
                nrow(demodata), 3)
  
  demodata_long	<- msprep(	time   = tt, 
                            status = st, 
                            data   = demodata, 
                            trans  = trans, 
                            id     = demodata$ID, 
                            keep   = c("Arm", "z1", "z2", "z3", "z4", 
                                       "arm0", "arm1", "timeM")
                          )
  
  # expand covariates (into ex. Arm.1, Arm.2, Arm.3 corresponding to the 3 transitions)
  demodata_long_cov	<- expand.covs(demodata_long, c("Arm", "z1", "z2", "z3", "z4", 
                                                "arm0", "arm1", "timeM"))      

# step 3 - fit models: ----

  ## simple cox model (model i in simulation study)
  fit.cox		<- coxph(Surv(timeD, DeathInd)~ Arm + MetaInd + z1 + z2 + z3 + z4, data=demodata)
  
  ## illness death model (model ii in simulation study)
  fit.msm	<- coxph(Surv(Tstart, Tstop, status) ~
                     Arm.1	 +	Arm.2	+	Arm.3	+
                     z1.1	   +	z1.2	+	z1.3	+
                     z2.1	   +	z2.2	+	z2.3	+
                     z3.1	   +	z3.2	+	z3.3	+
                     z4.1	   +	z4.2	+	z4.3	+ 
                     strata(trans), 
                     demodata_long_cov, method="breslow", timefix = FALSE)
  
  ## illness death model with time of transition to surrogate (model iii in simulation study)
  fit.msm.timeM	<- coxph(Surv(Tstart, Tstop, status) ~
                           Arm.1	 +	Arm.2	+	Arm.3	+
                           z1.1	   +	z1.2	+	z1.3	+
                           z2.1	   +	z2.2	+	z2.3	+
                           z3.1	   +	z3.2	+	z3.3	+
                           z4.1	   +	z4.2	+	z4.3	+ 
                           timeM.3 +
                           strata(trans), 
                           demodata_long_cov, method="breslow", timefix = FALSE)

# step 4 - Finally compute mediation measures according to the 3 models ----
  ## model (i)
  CoxMedRes 	           <- getMedMeas_cox(widedat  = demodata, 
                                          tau       = tau,
                                          cox.model = fit.cox)
  
  ## model (ii)
  mat 		               <- trans.illdeath()
  RiskMedRes 	           <- getMedMeas_risk(longdat        = demodata_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm, 
                                            transition.mat = mat, 
                                            CI             = FALSE)
  
  RmtlMedRes             <- getMedMeas_rmtl(longdat        = demodata_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm, 
                                            transition.mat = mat, 
                                            CI             = FALSE)
  
  ## model (iii)
  mat 		               <- trans.illdeath()
  RiskMedRes_timeM 	     <- getMedMeas_risk(longdat        = demodata_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm.timeM,
                                            transition.mat = mat,
                                            CI             = FALSE)
  
  RmtlMedRes_timeM 	     <- getMedMeas_rmtl(longdat        = demodata_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm.timeM,
                                            transition.mat = mat,
                                            CI             = FALSE)
  

  
  
