set more off
cd "...\SimulvsRE\"
capture log close
log using Si1000_TrueRE_N1000T4_beta_1.log, replace
clear ado
clear all
set matsize 10000
do "...\allToest_bCMLE.ado"
do "...\Toest_AMEfeT4.ado"
do "...\Toest_REMix.ado"


/* Simulate from a RE (with two unobserved types) dynamic model and compare the following estimators of the AME:
 - estimador of the AME based on the identification results in Aguirregabiria and Carro (RESTAT, 2024). In this filed it is denoted as CMLE because it uses as estimator of the beta paremeter, the Conditional Maximum Likelihood Estimator from Chamberlain.
 - RE with two unobserved types, which is the simulated distribution in the DGP here.
 - The MLE that assumes no unobserved heterogeneity.
*/


local nsimul=1000
global T=4 //This cannot be changed as the program called to obtain the CMLE estimate of the AME is written only for the T=4 case.
display "T=", $T
global N=1000
display "N=", $N
global nobs=$N * $T
set obs 1
gen kk=.
save Si1000_TrueRE_N1000T4_beta_1.dta, replace
clear

set seed 1810

forvalues i=1/`nsimul' {
display "simul", `i'

  set obs $nobs
  scalar beta=1 //Change this to -1 to obtain the other simulated case reported in the paper.
    //---- individual and time index
    gen  id=int( (_n-1)/$T ) + 1
    qbys id:gen t=_n
	tsset id t
	//---- time varying error 
	gen eps =rlogistic()
	//---- Permanent unobserved heterogeneity
	display "It is a mixture of two mass points"
	gen eta1=-1
	gen eta2=0.5
	gen ty1=(runiform()<0.3)
	gen etai=ty1*eta1+(1-ty1)*eta2
	
	qbys id: replace etai=etai[1]

	qbys id: gen y=(etai+eps>=0)   if _n==1
	qbys id: replace y=(beta*y[_n-1]+etai+eps>=0) if _n>1
	qbys id: gen ylag=y[_n-1]
	
	tab ty1 if t==1 & y==0
	tab ty1 if t==1 & y==1

	/*Calculate the true AME in this simulated data*/
	gen MEtrue=logistic(beta+etai)-logistic(etai)
	egen trueAME=mean(MEtrue)
	
	
/*FE CML Estimator of beta*/
	sort id t
	allToest_bCMLE, nI($N) nT($T)
	gen estb_CMLE=r(all_estb_CMLE)

/*FE Estimator of the AME (Aguirregabiria and Carro, RESTAT 2024)*/
	sort id t
	sum estb_CMLE
	local betae=r(mean)
    Toest_AMEfeT4, betae(`betae') nI($N) nT($T)
	gen estAME_cmle_bcmle=r(estAME_fe)
	
/*RE Estimator (with two discrete types)*/
	sum y if t==1
	global pry0h=r(mean)
	Toest_REMix, ntipos(2) nI($N) nT($T) pry0($pry0h)
	gen estb_REMix=r(estbeta_REMix)
	gen estAME_REMix=r(estAME_REMix)

/*Estimator without Unobserved Heterogeneity*/
	logit y i.ylag
	mat B=e(b)
	scalar estb_noUH=B[1,2]
	gen estb_noUH=estb_noUH
	margins, dydx(ylag)
	mat Bame=r(b)
	scalar estAME_noUH=Bame[1,2]
	gen estAME_noUH=estAME_noUH
	
	keep in 1
	gen beta=beta
    keep trueAME estAME_cmle_bcmle estAME_REMix estAME_noUH beta estb_CMLE estb_REMix estb_noUH
    append using Si1000_TrueRE_N1000T4_beta_1.dta
    capture drop kk
    save Si1000_TrueRE_N1000T4_beta_1.dta, replace
    clear 
}


use Si1000_TrueRE_N1000T4_beta_1.dta, clear
gen SE_cmle= (estAME_cmle_bcmle- trueAME)^2
egen RMSE_cmle=mean( SE_cmle)
replace RMSE_cmle=sqrt(RMSE_cmle)
gen SE_REMix=( estAME_REMix- trueAME)^2
egen RMSE_REMix=mean( SE_REMix)
replace RMSE_REMix=sqrt(RMSE_REMix)

save Si1000_TrueRE_N1000T4_beta_1.dta, replace

sum trueAME estAME_cmle_bcmle estAME_REMix estAME_noUH beta estb_CMLE estb_REMix estb_noUH
sum RMSE_cmle RMSE_REMix 
  
log close

	
