* note that you will need to install the package st0075
set more off

forvalues i=1/500 {
  	clear
	set obs 500
	generate a = rbinomial(1,0.5)
	generate z1 = rnormal(0,1)
	generate z2 = rnormal(0,1)
	generate z3 = rnormal(0,1)
	generate z4 = rnormal(0,1)

	* first generate the event times/indicators for competing risk of 0->1, 0->2 transition:
	* note that here the two numbers in all parentheses are the log(HR) for transition 1 (0->1) and transition 2 (0->2) 
	survsim stime event, distribution(weibull) cr ncr(2) lambdas(0.22 0.13) gammas (1.16 0.92) covariates(a -1.61 -1.24 z1 -0.015 0.09 z2 0.02 0.01 z3 0.37 0.09 z4 -0.83 -0.06)
	
	replace event = 0 if stime>8
	stset stime, failure(event==1)
	stcompet ci1=ci, compet1(2) by(a)
	
	stset stime, failure(event==2)
	stcompet ci2=ci, compet1(1) by(a)

	* second generate event times/indicators for 1->2 transition (we will do this for all participants and LATER determine path through model by
	* seeing if participant actually went down this transition path 

	survsim stime3, distribution(weibull) lambdas(0.42) gammas (1.94) covariates(a -0.78 z1 0.03 z2 -0.01 z3 0.03 z4 0.83)
	generate died3 = stime3 <= 8
	stset stime3, failure(died3=1)
  
  export delimited using "<insert_file_path>\stata_data_`i'.txt", delimiter(tab)
  
}
