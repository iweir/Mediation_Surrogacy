getMedMeas_rmtl <- function(longdata, tau, model, transition.mat, CI=FALSE){

  rmtl_byID_0 <- matrix(NA, max(longdata$id), 2+(length(tau)) )
  rmtl_byID_1 <- matrix(NA, max(longdata$id), 2+(length(tau)) )
  
    i         <- 0
  
  for(id in unique(longdata$id)){ 
    
    i               <- i+1
    tempdat         <- longdata[which(longdata$id==id),]
    tempdat$strata  <- tempdat$trans
    
    #----- ARM = 0 -----#
    
    # set arm variables for arm=0
    tempdat$Arm.1   <- tempdat$arm0.1
    tempdat$Arm.2   <- tempdat$arm0.2
    tempdat$Arm.3   <- tempdat$arm0.3
    if(any(names(tempdat) =="Arm.4")){ 
      tempdat$Arm.4 <- tempdat$arm0.4
    }
    
    # Predicted cumulative hazards using msfit
    msfit.temp0     <- msfit(model, tempdat, trans = transition.mat)
    
    # Predicted transition probabilities for counterfactual arm=0:
    pt.temp         <- probtrans(msfit.temp0, 0, variance=FALSE)[[1]][,c("time", "pstate3")]
    
    # RMTL at time points for counterfactual arm = 0
    rmtl_arm0       <- lapply(tau,function(x){AUCiw(pt.temp$time, pt.temp$pstate3, to = max(x))})
    temp_rmtl_arm0  <- do.call(cbind,rmtl_arm0)
    
    #----- ARM = 1 -----#
    
    # set arm variables for arm=1
    tempdat$Arm.1   <- tempdat$arm1.1
    tempdat$Arm.2   <- tempdat$arm1.2
    tempdat$Arm.3   <- tempdat$arm1.3
    if(any(names(tempdat) =="Arm.4")){ 
      tempdat$Arm.4 <- tempdat$arm1.4
    }
    
    # Predicted cumulative hazards using msfit
    msfit.temp1     <- msfit(model, tempdat, trans = transition.mat)
    
    # Predicted transition probabilities for counterfactual arm=1:
    pt.temp         <- probtrans(msfit.temp1, 0, variance=FALSE)[[1]][,c("time", "pstate3")]
    
    
    # RMTL at time points for counterfactual arm = 1
    rmtl_arm1       <- lapply(tau,function(x){AUCiw(pt.temp$time, pt.temp$pstate3, to = max(x))})
    temp_rmtl_arm1  <- do.call(cbind,rmtl_arm1)
    
    rmtl_byID_0[i,1]                   <- tempdat$id[1]
    rmtl_byID_0[i,2]                   <- tempdat$Arm[1]
    rmtl_byID_0[i,3:ncol(rmtl_byID_0)] <- temp_rmtl_arm0
    
    rmtl_byID_1[i,1] 	                 <- tempdat$id[1]
    rmtl_byID_1[i,2]                   <- tempdat$Arm[1]
    rmtl_byID_1[i,3:ncol(rmtl_byID_1)] <- temp_rmtl_arm1
    
  }
  
  # obtain aggregate measures:
  res0 <- aggregate(rmtl_byID_0[,3:ncol(rmtl_byID_0)], list(trueArm=rmtl_byID_0[,2]), mean)
  res1 <- aggregate(rmtl_byID_1[,3:ncol(rmtl_byID_1)], list(trueArm=rmtl_byID_1[,2]), mean)
  
  # compute mediation measures according to RMTL difference for each time point
  
  all_results <- data.frame(time     = numeric(length(tau)), 
                            indirect = numeric(length(tau)),
                            direct   = numeric(length(tau)),
                            total    = numeric(length(tau)),
                            PM       = numeric(length(tau)), stringsAsFactors=FALSE)
  
  index <- 0
  for(t in tau){
    index                   <- index+1
    all_results$time[index] <- t
    
    ST1M1 	<- res1[2,index+1]
    ST0M0 	<- res0[1,index+1]
    ST1M0 	<- res1[1,index+1]
    
    all_results$indirect[index] <- ST1M1-ST1M0
    all_results$direct[index]   <- ST1M0-ST0M0
    all_results$total[index]    <- ST1M1-ST0M0
    all_results$PM[index]       <- (ST1M1-ST1M0)/(ST1M1-ST0M0)
  }
  
  # OPTION HERE TO COMPUTE CONFIDENCE INTERVALS WITH PERTURBATION RESAMPLING  
  
  if(CI==TRUE){  
    
    CIresults <- data.frame(time         = numeric(length(tau)), 
                            indirect_low = numeric(length(tau)),
                            indirect_up  = numeric(length(tau)),
                            direct_low   = numeric(length(tau)),
                            direct_up    = numeric(length(tau)),
                            total_low    = numeric(length(tau)),
                            total_up     = numeric(length(tau)),
                            pm_low       = numeric(length(tau)),
                            pm_up        = numeric(length(tau)), stringsAsFactors=FALSE)    
    index <- 0
    
    for(t in tau){
      index       <- index+1	  
      B 		      <- 1000
      b_direct 	  <- as.numeric()
      b_indirect 	<- as.numeric()
      b_total 	  <- as.numeric()
      b_pm 	    	<- as.numeric()
      
      for(b in 1:B){
        
        n_C		<- nrow(rmtl_byID_0[which(rmtl_byID_0[,2]==0),])
        n_E		<- nrow(rmtl_byID_0[which(rmtl_byID_0[,2]==1),])
        
        Uc 		<- runif(n_C,0,1)
        Ue 		<- runif(n_E,0,1)
        Vc 		<- -log(Uc)
        Ve 		<- -log(Ue)
        
        ViPi_S1M1 	<- Ve*rmtl_byID_1[which(rmtl_byID_1[,2]==1), index+2]
        ViPi_S0M0 	<- Vc*rmtl_byID_0[which(rmtl_byID_0[,2]==0), index+2]
        ViPi_S1M0 	<- Vc*rmtl_byID_1[which(rmtl_byID_1[,2]==0), index+2]
        
        S1M1 		<- (1/n_E)*sum(ViPi_S1M1)
        S0M0 		<- (1/n_C)*sum(ViPi_S0M0)
        S1M0 		<- (1/n_C)*sum(ViPi_S1M0)
        
        # calculate mediation measures
        d		<- S1M0-S0M0
        i		<- S1M1-S1M0
        tot	<- d+i
        p		<- i/tot
        
        b_direct	  <- c(b_direct, d)
        b_indirect	<- c(b_indirect, i)
        b_total   	<- c(b_total, tot)
        b_pm		    <- c(b_pm, p)
        
      }
      
      direct_low		<- quantile(b_direct, c(0.025, 0.975))[[1]]
      direct_up		  <- quantile(b_direct, c(0.025, 0.975))[[2]]
      
      indirect_low	<- quantile(b_indirect, c(0.025, 0.975))[[1]]
      indirect_up		<- quantile(b_indirect, c(0.025, 0.975))[[2]]
      
      total_low	  	<- quantile(b_total, c(0.025, 0.975))[[1]]
      total_up		  <- quantile(b_total, c(0.025, 0.975))[[2]]
      
      pm_low	    	<- quantile(b_pm, c(0.025, 0.975))[[1]]
      pm_up			    <- quantile(b_pm, c(0.025, 0.975))[[2]]
      
      CIresults$time[index]         <- t
      CIresults$direct_low[index]   <- direct_low
      CIresults$direct_up[index]    <- direct_up
      CIresults$indirect_low[index] <- indirect_low
      CIresults$indirect_up[index]  <- indirect_up
      CIresults$total_low[index]    <- total_low
      CIresults$total_up[index]     <- total_up
      CIresults$pm_low[index]       <- pm_low
      CIresults$pm_up[index]        <- pm_up
      
    }
   
  all_results <- merge(all_results, CIresults)   
  }
  
return(all_results)

}
