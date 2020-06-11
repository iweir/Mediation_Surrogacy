getMedMeas_cox <- function(widedat, tau, cox.model){

	new_0		<-  data.frame(	Arm      = 0, 
			                timeD 	 = widedat$timeD, 
					DeathInd = widedat$DeathInd, 
					MetaInd	 = widedat$MetaInd,
					z1	 = widedat$z1,
					z2       = widedat$z2,
					z3	 = widedat$z3,
					z4	 = widedat$z4)

	new_1		<-  data.frame(	Arm       = 1, 
					timeD 	  = widedat$timeD, 
					DeathInd  = widedat$DeathInd, 
					MetaInd   = widedat$MetaInd,
					z1	  = widedat$z1,
					z2	  = widedat$z2,
					z3	  = widedat$z3,
					z4	  = widedat$z4)

	#----- 1-P(T(1,0) > t) -----#
	pp_10       <- 1-predictProb(cox.model, Surv(widedat$timeD, widedat$DeathInd), x=new_1, times=tau)
	RMLE_10     <- mean(cbind(pp_10, widedat$Arm)[widedat$Arm==0,1])

	#----- 1-P(T(0,0) > t) -----#
	pp_00       <- 1-predictProb(cox.model, Surv(widedat$timeD, widedat$DeathInd), x=new_0, times=tau)
	RMLE_00     <- mean(cbind(pp_00, widedat$Arm)[widedat$Arm==0,1])

	#----- 1-P(T(1,1) > t) -----# 
	pp_11       <- 1-predictProb(cox.model, Surv(widedat$timeD, widedat$DeathInd), x=new_1, times=tau)
	RMLE_11     <- mean(cbind(pp_11, widedat$Arm)[widedat$Arm==1,1])

	# direct effect:
	direct      <- RMLE_10-RMLE_00

	# indirect effect:
	indirect    <- RMLE_11-RMLE_10

	# total effect:
	total       <- RMLE_11-RMLE_00

	# mediation proportion:
	PM     	    <- indirect/total

	all_results <-	list(direct   = direct,
                             indirect = indirect,
                             total    = total,
			     PM	      = PM,
		             RMLE_11  = RMLE_11,
		             RMLE_00  = RMLE_00,
                             RMLE_10  = RMLE_10)

return(all_results)

}
