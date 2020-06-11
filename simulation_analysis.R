# step 1 - prepare work space - specify path, load libraries: ---------------------------------------

path    <- "/your/path/here"

#  I run each row of the scenarioTable.txt in parallel so 'j' refers to the 
#  particular scenario (row) you want to run 
j <- 1

# load libraries and functions
library(peperr)
library(survival)
library(mstate)

source(paste0(path, "/FUN_getMedMeas_risk.R"))
source(paste0(path, "/FUN_getMedMeas_Cox.R"))

# read in scenarios 
scen <- read.table(paste0(path, "/scenarioTable.txt"), header=TRUE, sep="\t") 

start   <- scen$start[j]
stop    <- scen$stop[j]
uni_low <- scen$uni_low[j]
uni_up  <- scen$uni_up[j]
scenNum <- scen$scenNum[j]

# prepare results storage
Niter		          <- stop-start+1
WEIR_result       <- matrix(NA, Niter, 5)
WEIR_result_timeM <- matrix(NA, Niter, 5)
Cox_result	      <- matrix(NA, Niter, 5)
data_descr        <- matrix(NA, Niter, 45)

# step 2 - Begin analysis on each replicate:
i <- 0
for(iteration in start:stop){
  
  i <- i +1
  print(paste0("scenario = ", scenNum, " and current iteration = ", iteration))
  
  puredat <- read.table(paste0(path, "/datafiles/stata_data_", iteration, ".txt"), header=TRUE, sep="\t")
 
  tau              <- 8
  prepdat          <- puredat[,c("a", "z1", "z2", "z3", "z4", "event", "stime", "stime3")]
  prepdat$ID       <- 1:nrow(prepdat)
  prepdat$Arm      <- prepdat$a
  prepdat$arm0     <- 0
  prepdat$arm1     <- 1
  
  # step 2a - introduce censoring and determine timeM and timeD
  
  prepdat$ctime   <- runif(nrow(prepdat), uni_low, uni_up) 
  
  prepdat$mintime <- pmin(tau, prepdat$stime, prepdat$ctime)
  
  prepdat$timeM   <- prepdat$mintime
  prepdat$MetaInd <- ifelse(prepdat$event%in% c(0,2) | prepdat$timeM==prepdat$ctime | prepdat$timeM==tau, 0, 1)
  
  prepdat$timeD    <- ifelse(prepdat$event==1, pmin(tau, prepdat$stime+prepdat$stime3, prepdat$ctime), prepdat$timeM)
  prepdat$DeathInd <- ifelse(   prepdat$event==0 | 
                                  (prepdat$event %in% c(1,2) & prepdat$timeD==tau) | 
                                  (prepdat$event %in% c(1,2) & prepdat$timeD==prepdat$ctime),
                                0, 1)
  
  # step 2b - calculate proportion censored
  data_descr[i,1] <- iteration
  data_descr[i,2] <- nrow(prepdat[which(prepdat$timeM==prepdat$ctime),])/nrow(prepdat)
  data_descr[i,3] <- nrow(prepdat[which(prepdat$timeD==prepdat$ctime),])/nrow(prepdat)
  
  # proportion truncation censored at Tau 
  data_descr[i,4] <- nrow(prepdat[which(prepdat$timeM==tau),])/nrow(prepdat)
  data_descr[i,5] <- nrow(prepdat[which(prepdat$timeD==tau),])/nrow(prepdat)
  
  keepdat          <- prepdat[,c("ID", "Arm", 
                                 "z1", "z2", "z3", "z4", 
                                 "timeM", "timeD", "MetaInd", "DeathInd", 
                                 "ctime", "stime", "stime3", "event",
                                 "arm0", "arm1")]
  
  # step 2c - Restructure data with the mstate functions ----
  trans	<- trans.illdeath()
  
  tt		<- matrix(c(rep(NA,nrow(keepdat)), 
                  keepdat$timeM,keepdat$timeD), 		
                nrow(keepdat), 3)
  
  st		<- matrix(c(rep(NA,nrow(keepdat)), 
                  keepdat$MetaInd,keepdat$DeathInd), 	
                nrow(keepdat), 3)
  
  dat_long	<- msprep(	time   = tt, 
                      status = st, 
                      data   = keepdat, 
                      trans  = trans, 
                      id     = keepdat$ID, 
                      keep   = c("Arm", "z1", "z2", "z3", "z4", 
                                 "arm0", "arm1", "timeM")
  )
  
  # expand covariates (into ex. Arm.1, Arm.2, Arm.3 corresponding to the 3 transitions)
  dat_long_cov	<- expand.covs(dat_long, c("Arm", "z1", "z2", "z3", "z4", 
                                          "arm0", "arm1", "timeM"))      
  
  data_descr[i,6] <- events(dat_long_cov)$Proportions[1,2]
  data_descr[i,7] <- events(dat_long_cov)$Proportions[1,3]
  data_descr[i,8] <- events(dat_long_cov)$Proportions[2,3]
  

  # step 2d - fit models: ----
  ## model (i)
  fit.cox		<- coxph(Surv(timeD, DeathInd)~ Arm + MetaInd + z1 + z2 + z3 + z4, data=keepdat)
  
  ## model (ii)
  fit.msm	<- coxph(Surv(Tstart, Tstop, status) ~
                     Arm.1	 +	Arm.2	+	Arm.3	+
                     z1.1	   +	z1.2	+	z1.3	+
                     z2.1	   +	z2.2	+	z2.3	+
                     z3.1	   +	z3.2	+	z3.3	+
                     z4.1	   +	z4.2	+	z4.3	+ 
                     strata(trans), 
                   dat_long_cov, method="breslow", timefix = FALSE)
  
  ## model (iii)
  fit.msm.timeM	<- coxph(Surv(Tstart, Tstop, status) ~
                           Arm.1	 +	Arm.2	+	Arm.3	+
                           z1.1	   +	z1.2	+	z1.3	+
                           z2.1	   +	z2.2	+	z2.3	+
                           z3.1	   +	z3.2	+	z3.3	+
                           z4.1	   +	z4.2	+	z4.3	+ 
                           timeM.3 +
                           strata(trans), 
                           dat_long_cov, method="breslow", timefix = FALSE)

  # save model coefficients
  data_descr[i,9:23]  <- fit.msm$coef
  data_descr[i,24:39] <- fit.msm.timeM$coef
  data_descr[i,40:45] <- fit.cox$coef
  
  # step 2e - Finally compute mediation measures according to the 3 models ----
  ## model (i)
  CoxMedRes 	           <- getMedMeas_cox(widedat   = keepdat, 
                                           tau       = tau,
                                           cox.model = fit.cox)
  
  Cox_result[i,1]        <- iteration
  Cox_result[i,2]        <- CoxMedRes$direct
  Cox_result[i,3]        <- CoxMedRes$indirect
  Cox_result[i,4]        <- CoxMedRes$total
  Cox_result[i,5]        <- CoxMedRes$PM
  
  ## model (ii)
  mat 		               <- trans.illdeath()
  medRes 	               <- getMedMeas_risk(longdat        = dat_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm, 
                                            transition.mat = mat, 
                                            CI             = FALSE)
  
  WEIR_result[i,1]       <- iteration
  WEIR_result[i,2]       <- medRes$direct[1] 
  WEIR_result[i,3]       <- medRes$indirect[1] 
  WEIR_result[i,4]       <- medRes$total[1] 
  WEIR_result[i,5]       <- medRes$PM[1] 
  
  ## model (iii)
  mat 		               <- trans.illdeath()
  medRes_timeM 	         <- getMedMeas_risk(longdat        = dat_long_cov, 
                                            tau            = tau, 
                                            model          = fit.msm.timeM,
                                            transition.mat = mat,
                                            CI             = FALSE)
  
  WEIR_result_timeM[i,1] <- iteration
  WEIR_result_timeM[i,2] <- medRes_timeM$direct[1] 
  WEIR_result_timeM[i,3] <- medRes_timeM$indirect[1] 
  WEIR_result_timeM[i,4] <- medRes_timeM$total[1] 
  WEIR_result_timeM[i,5] <- medRes_timeM$PM[1] 
  
}

# step 3 - add helpful column names and save results -------------------------------------------------------------
colnames(WEIR_result)     	<- c(	"iteration",
                                 "direct",
                                 "indirect",
                                 "total",
                                 "PM")

colnames(WEIR_result_timeM) <- c(	"iteration",
                                  "direct",
                                  "indirect",
                                  "total",
                                  "PM")

colnames(Cox_result)		    <- c(	"iteration",
                                 "direct",
                                 "indirect",
                                 "total",
                                 "PM")

colnames(data_descr)      	<- c(	"iteration",
                                 "propCensM",
                                 "propCensD",
                                 "propCensMtau",
                                 "propCensDtau",
                                 "prop12",
                                 "prop13",
                                 "prop23",
                                 "W1_Arm.1", "W1_Arm.2", "W1_Arm.3",
                                 "W1_x1.1" , "W1_x1.2" , "W1_x1.3",
                                 "W1_x2.1" , "W1_x2.2" , "W1_x2.3",
                                 "W1_x3.1" , "W1_x3.2" , "W1_x3.3",
                                 "W1_x4.1" , "W1_x4.2" , "W1_x4.3",
                                 "W2_Arm.1", "W2_Arm.2", "W2_Arm.3",
                                 "W2_x1.1" , "W2_x1.2" , "W2_x1.3",
                                 "W2_x2.1" , "W2_x2.2" , "W2_x2.3",
                                 "W2_x3.1" , "W2_x3.2" , "W2_x3.3",
                                 "W2_x4.1" , "W2_x4.2" , "W2_x4.3",
                                 "W2_timeM.3",
                                 "V_Arm", "V_MetaInd", 
                                 "V_x1", "V_x2", "V_x3", "V_x4")

write.table(WEIR_result, file=paste0(path, "/res_WEIR_scenario", scenNum, "_", start, "_to_", stop, ".txt"), row.names=FALSE, sep="\t")
write.table(WEIR_result_timeM, file=paste0(path, "/res_WEIR_timeM_scenario", scenNum, "_", start, "_to_", stop, ".txt"), row.names=FALSE, sep="\t")
write.table(Cox_result, file=paste0(path, "/res_cox_scenario", scenNum, "_", start, "_to_", stop, ".txt"), row.names=FALSE, sep="\t")
write.table(data_descr, file=paste0(path, "/res_datadescr_scenario", scenNum, "_", start, "_to_", stop, ".txt"), row.names=FALSE, sep="\t")
