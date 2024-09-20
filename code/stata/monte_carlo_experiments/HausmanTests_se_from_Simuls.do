cd "D:\FEident_structural\SimulvsRE\"
capture log close
log using HausmanTests_se_from_Simuls.log, replace

set more off

/* 
This program performs the Hausman Tests based on beta and on AME.
I use as variance of the estimators the empirical variance from the 1000 simulations.
*/

clear all

*******True Model: RE Finite Mixture, beta=-1********
use Si1000_TrueRE_N1000T4_beta_m1.dta, clear

sum estb_CMLE
gen var_estb_CMLE=r(sd)^2
sum estAME_cmle_bcmle
gen var_estAME_cmle=r(sd)^2
sum estb_REMix
gen var_estb_REMix=r(sd)^2
sum estAME_REMix
gen var_estAME_REMix=r(sd)^2
sum estb_noUH
gen var_estb_noUH=r(sd)^2
sum estAME_noUH
gen var_estAME_noUH=r(sd)^2

	*To Figure 1: validity of NoUH model in DGP FinMix(-1)
	gen HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE-var_estb_noUH)
	*The following order should not cause any change:
		replace HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE) if var_estb_CMLE<var_estb_noUH
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle-var_estAME_noUH)
	*The following order should not cause any change:
		replace HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_noUH
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: NoUH vs. FE in DGP FinMix(-1)) legend(on)
	graph save Graph Figure1.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2


*******True Model: RE Finite Mixture, beta=1********
use Si1000_TrueRE_N1000T4_beta_1.dta, clear

sum estb_CMLE
gen var_estb_CMLE=r(sd)^2
sum estAME_cmle_bcmle
gen var_estAME_cmle=r(sd)^2
sum estb_REMix
gen var_estb_REMix=r(sd)^2
sum estAME_REMix
gen var_estAME_REMix=r(sd)^2
sum estb_noUH
gen var_estb_noUH=r(sd)^2
sum estAME_noUH
gen var_estAME_noUH=r(sd)^2

	*To Figure 2: validity of NoUH model in DGP FinMix(+1)
	gen HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE-var_estb_noUH)
	*The following order should not cause any change:
		replace HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE) if var_estb_CMLE<var_estb_noUH
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle-var_estAME_noUH)
	*The following order should not cause any change:
		replace HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_noUH
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: NoUH vs. FE in DGP FinMix(+1)) legend(on)
	graph save Graph Figure2.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2

*******True Model: FE Mixture of Normals, beta=-1********
use Si1000_TrueFE_N1000T4_beta_m1.dta, clear

sum estb_CMLE
gen var_estb_CMLE=r(sd)^2
sum estAME_cmle_bcmle
gen var_estAME_cmle=r(sd)^2
sum estb_REMix
gen var_estb_REMix=r(sd)^2
sum estAME_REMix
gen var_estAME_REMix=r(sd)^2
sum estb_noUH
gen var_estb_noUH=r(sd)^2
sum estAME_noUH
gen var_estAME_noUH=r(sd)^2

	*To Figure 3: validity of NoUH model in DGP MixNor(-1)
    gen HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE-var_estb_noUH)
	*The following order should not cause any change:
		replace HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE) if var_estb_CMLE<var_estb_noUH
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle-var_estAME_noUH)
	*The following order should not cause any change:
		replace HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_noUH
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: NoUH vs. FE in DGP MixNor(-1)) legend(on)
	graph save Graph Figure3.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2

	*To Figure 5: validity of RE finite mixture model in DGP MixNor(-1)
	gen HS1= (estb_CMLE-estb_REMix)^2/(var_estb_CMLE-var_estb_REMix)
	*The following order is just in case estimated var of the FE is smaller than est Var of the RE:
*		replace HS1= (estb_CMLE-estb_REMix)^2/(var_estb_CMLE) if (var_estb_CMLE<var_estb_REMix) 
		replace HS1= (estb_CMLE-estb_REMix)^2/(var_estb_CMLE) if (var_estb_CMLE<var_estb_REMix) | (var_estAME_cmle<var_estAME_REMix)
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_REMix)^2/(var_estAME_cmle-var_estAME_REMix)
	*The following order is just in case estimated var of the FE is smaller than est Var of the RE:
*		replace HS2= (estAME_cmle_bcmle-estAME_REMix)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_REMix
		replace HS2= (estAME_cmle_bcmle-estAME_REMix)^2/(var_estAME_cmle) if (var_estb_CMLE<var_estb_REMix) | (var_estAME_cmle<var_estAME_REMix)
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: RE Finite Mix vs. FE) subtitle(DGP MixNor(-1)) legend(on)
	graph save Graph Figure5.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2

*******True Model: FE Mixture of Normals, beta=1********
use Si1000_TrueFE_N1000T4_beta_1.dta, clear

sum estb_CMLE
gen var_estb_CMLE=r(sd)^2
sum estAME_cmle_bcmle
gen var_estAME_cmle=r(sd)^2
sum estb_REMix
gen var_estb_REMix=r(sd)^2
sum estAME_REMix
gen var_estAME_REMix=r(sd)^2
sum estb_noUH
gen var_estb_noUH=r(sd)^2
sum estAME_noUH
gen var_estAME_noUH=r(sd)^2

	*To Figure 4: validity of NoUH model in DGP MixNor(+1)
	gen HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE-var_estb_noUH)
	*The following order should not cause any change:
		replace HS1= (estb_CMLE-estb_noUH)^2/(var_estb_CMLE) if var_estb_CMLE<var_estb_noUH
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle-var_estAME_noUH)
	*The following order should not cause any change:
		replace HS2= (estAME_cmle_bcmle-estAME_noUH)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_noUH
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: NoUH vs. FE in DGP MixNor(+1)) legend(on)
	graph save Graph Figure4.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2

	*To Figure 6: validity of RE finite mixture model in DGP MixNor(+1)
	gen HS1= (estb_CMLE-estb_REMix)^2/(var_estb_CMLE-var_estb_REMix)
	*The following order should not cause any change:
		replace HS1= (estb_CMLE-estb_REMix)^2/(var_estb_CMLE) if (var_estb_CMLE<var_estb_REMix)
	***
	gen pvalorHS1=chi2tail(1,HS1)
	cumul pvalorHS1, gen(cumulHS1)
	label variable cumulHS1 "HS1: test on beta"
	gen HS2= (estAME_cmle_bcmle-estAME_REMix)^2/(var_estAME_cmle-var_estAME_REMix)
	*The following order should not cause any change:
		replace HS2= (estAME_cmle_bcmle-estAME_REMix)^2/(var_estAME_cmle) if var_estAME_cmle<var_estAME_REMix
	***
	gen pvalorHS2=chi2tail(1,HS2)
	cumul pvalorHS2, gen(cumulHS2)
	label variable cumulHS2 "HS2: test on AME"
	twoway (line cumulHS1 pvalorHS1, sort lpattern(dash)) (line cumulHS2 pvalorHS2, sort lpattern(solid)), ytitle(Cumulative distribution) xtitle(p-value) title(p-values for Hausman Tests: RE Finite Mix vs. FE in DGP MixNor(+1)) legend(on)
	graph save Graph Figure6.gph, replace
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2
	sum estb_CMLE estAME_cmle_bcmle estb_REMix estAME_REMix HS1 pvalorHS1 HS2 pvalorHS2 if pvalorHS2>0.1
	drop HS1 pvalorHS1 cumulHS1 HS2 pvalorHS2 cumulHS2

log close

	
