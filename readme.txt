README File for the
Replication of E mpirical Results in:
"Identification of Average Marginal Effects
in Fixed Effects Dynamic Discrete Choice models"
Victor Aguirregabiria, University of Toronto, CEPR
Jesús M. Carro, Universidad Carlos III de Madrid

September 17, 2024

1 Introduction

This README file outlines the steps to replicate the results presented in: Aguirregabiria,
Victor, and Jesús Carro (Forthcoming). "Identification of Average Marginal Effects in Fixed
Effects Dynamic Discrete Choice Models,"
Review of Economics and Statistics

. All the code and data, including this README file, are in the folder:

/replication_package

The paper provides three sets of empirical results:

1. Results from the Monte Carlo exp eriments in Section 4, as shown in Tab le 1 and Figure 1 in the paper.
2. Descriptive statistics from the empirical app lication in Section 5.1, rep orted in Table 2 in the paper.
3. Model estimation results from the empirical application in Section 5.2, detailed in Tables 3 and 4 in the paper.

Below, we describe the code and data required for each of these three sets of results.

2 Monte Carlo Exp eriments in Section 4

The .do and .ado files for replicating th e Monte Carlo exp eriments an d generating the
results for Table 1 and Figure 1 in the paper are lo cated in the following folder:

/replication_package/monte_carlo_experiments

It is written in Stata and was run using Stata version 14MP. It consists of the following
three .ado files and four .do files:

allToest_bCMLE.ado	-- This Stata ado file contains a program for the Conditional Maximum
Likelihood estimation of slope parameters in Fixed Effects model. This ado file is called by
the do files.

Toest_REMix.ado 	-- This Stata ado file contai ns a program for the Maximum Likelihood
estimation of slope parameters in Random Effects model. This ado file is called by the do
files.

Toest_AMEfeT4.ado	-- This Stata ado file contains a program for the estimation of Average
Marginal Effects. This ado file is called by the do files.

Simul_TrueNoUH_vsRE.do	-- This Stata do file:
	1. Simulates data from the DGP without p ermanent unobserved heterogeneity: DGPs NoUH(-1)
	and NoUH(+1)
	2. Estimates slope parameter of the Fixed Effect model by CMLE calling program allToest_bCMLE.ado
	3. Estimates slope parameter of the model without unobserved heterogeneity.
	4. Estimates slope parameter of the Random Effect model by MLE calling program Toest_REMix.ado
	5. Estimates AME for the three di˙erent models calling program Toest_AMEfeT4.ado
	6. Using the estimates from all simulated samples, it calculates the means and standard deviations 
	of the estimators, wh ich are presented in rows NoUH(-1) and NoUH(+1) in Table 1 of the paper.

Simul_TrueRE_vsRE.do	-- This Stata do file:
	1. Simulates data from the D GP with Finite Mixture unobserved heterogeneity: DGPs FinMIx(-1) and FinMIx(+1)
	2. Estimates slope parameter of the Fixed Effect model by CMLE calling program allToest_bCMLE.ado
	3. Estimates slope parameter of the model without unobserved heterogeneity.
	4. Estimates slope parameter of the Random Effect model by MLE calling program Toest_REMix.ado
	5. Estimates AME for the three different models calling program Toest_AMEfeT4.ado
	6. Using the estimates from all simulated samples, it calculates the means and standard deviations of the 
	estimators, which are presented in rows FinMIx(-1) and FinMIx(+1) in Table 1 of the paper.

Simul_TrueFE_vsRE.do	-- This Stata do file:
	1. Simulates data from the DGP with Mixture of N ormals unobserved heterogeneity: DGPs MixNor(-1) and MixNor(+1)
	2. Estimates slope parameter of the Fixed Effect model by CMLE calling program allToest_bCMLE.ado
	3. Estimates slope parameter of the model without unobserved heterogeneity.
	4. Estimates slope parameter of the Random Effect model by MLE calling program Toest_REMix.ado
	5. Estimates AME for the three different models calling program Toest_AMEfeT4.ado
	6. Using the estimates from all simulated samples, it calculates the means and standard deviations of the estimators, 
	which are presented in rows MixNor(-1) and MixNor(+1) in Table 1 of the paper.

HausmanTests_se_from_Simuls.do	-- This Stata do file reads th e vectors of Monte Carlo estimates and standard errors generated 
by the do files Simul_TrueNoUH_vsRE.do, Simul_TrueRE_vsRE.do, and Simul_TrueFE_vsRE.do. It then calculates the Hausman statistics 
along with the corresponding p-values and generates the graphs shown in Figure 1 of the paper. Naturally, this do file should be 
executed after running the aforementioned do files. 

Instructions for replicating results:
	- The provided versions of the do files implement the Monte Carlo exp eriments for a value of the slope parameter beta = 1. 
	To implement the corresp onding Monte Carlo exp eriment for a value of the slop e parameter beta = -1, the user only needs 
	to change code line 35 in these do files. Specifically, replace co de line:
		scalar beta=1
	with
		scalar beta=-1

	- Do file HausmanTests_se_from_Simuls.do generating Figure 1 should be executed after running the do files
	Simul_TrueNoUH_vsRE.do, Simul_TrueRE_vsRE.do, and Simul_TrueFE_vsRE.do

3 Descriptive Statistics from Empirical Application: Section 5.1, Table 2

The Stata datafile and the .do file for generating the transition probability matrix in Table 2 in the paper are located 
in the following folder:

/replication_package/descriptive_statistics

The folder contains two files:
	consumer_withprices.dta -- Stata datafile with Nielsen consumer scanner panel data provided by Susumu Imai.

	application_trans_matrix.do -- Stata .do file that reads the datafile consumer_withprices.dta and generates 
	descriptive statistics including the transition probability matrix in Table 2 in the paper.

Instructions for replicating results:
	- Line 25 in the code loads the data file using the following command:
		use consumer_withprices.dta
	Make sure to specify the full folder path where this data file is stored on your computer.

4 Estimation of Parameters and AMEs in Empirical Application: Section 5.2, Tables 3 and 4

The datafiles and code for generating the statistics in Tables 3 and 4 in the paper are lo cated in the following folder:

	/replication_package/application_estimates_beta_ame

The folder contains three files:

	generate_xls_datafile.do	-- Stata .do file that reads that Stata datafile consumer_withprices.dta
	and generates the .xls file consumer_feestimation.xls. This Excel file is the input of the GAUSS
	program for the Fixed Effects estimation of slope parameters and AMEs.

	consumer_feestimation.xls	-- Datafile Excel format created by generate_xls_datafile.do

	fe_dyn_logit_EIKdata.gss	-- GAUSS program for the Fixed Effects estimation of slope parameters and AMEs 
	in Tables 3 and 4 in the paper. This code is self-contained as it includes all the procedures and functions 
	called by the main program.

Instructions for replicating results:
	- The program fe_dyn_logit_EIKdata.gss has been executed using GAUSS version 23.0. However, it should also work well 
	in earlier versions at least up to GAUSS version 12.0

