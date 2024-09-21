clear

// ---------------------------------------------- 
// application_trans_matrix.do
//
//
// 	It constructs dataset with information on
// 	consumer choice histories calculates
// 	matrix of transition probabilities 
//	in Table 2 of the paper:
//		Aguirregabiria, Victor, and Jesus Carro 
//		"Identification of Average Marginal Effects
//		in Fixed Effects Dynamic Discrete Choice Models"
//		Review of Economics and Statistics
//
// ---------------------------------------------- 

capture log close
log using application_trans_matrix.log, replace

// --------------------------------
//	1.  Reading datafile 
// -------------------------------- 

use consumer_withprices.dta

// -----------------------------------------------------------------
//	2.  Variable selection and construction of new variables 
// -----------------------------------------------------------------

// 2.1. Selecting variables that will be used
keep idconsumer week store brand size price_ind price_store

// 2.2. Removing consumer-week observations with no purchases of ketchup
drop if brand==0

// 2.3. Relabelling values of the variable brand such that:
//		- Following the paper notation, brand values are 0,1,2,3 instead of 1,2,3,4.
//		- Brand values are sorted according to aggregate market shares
//			with 0 for the largest (Heinz) and 3 for the lowest (Store brand)
//
//		Brand name	Old-code	New-code
//		Heinz		3			0
//		Hunt's		4			1
//		Del Monte	2			2
//		Store brand	1			3

replace brand = 0*(brand==3) + 1*(brand==4) + 2*(brand==2) + 3*(brand==1)

// Define the labels
label define brand_labels 0 "Heinz" 1 "Hunt's" 2 "Del Monte" 3 "Store brand"

// Apply the labels to the variable 'brand'
label values brand brand_labels

// 2.4. Creating variable T_i with number of purchasing events per consumer
egen Ti = sum(1), by(idconsumer)
tab Ti

// 2.5. Creating variable conta with values 1, 2, ... Ti 
//		counting the order of the purchasing events of eahc consumer
sort idconsumer week
gen conta = .
replace conta = 1 if idconsumer~=idconsumer[_n-1]
replace conta = conta[_n-1] + 1 if idconsumer==idconsumer[_n-1]

// 2.6. Constructing variable with consumer brand-choice history 
sort idconsumer week
gen history = .
replace history = brand * (10^(Ti-1)) if idconsumer~=idconsumer[_n-1]
replace history = history[_n-1] + brand * (10^(Ti-conta)) if idconsumer==idconsumer[_n-1]

// 2.7. Constructing variable with consumer brand-choice history: only first 4 purchases
sort idconsumer week
gen buff = .
replace buff = brand * 1000 if idconsumer~=idconsumer[_n-1]
replace buff = buff[_n-1] + brand *(10^(4-conta)) if conta>=2 & conta<=4
egen hist4 = sum(buff * (conta==4)), by(idconsumer)
drop buff 

// ----------------------------------------------------- 
// 	3.  Generating some Descriptive Statistics
// ----------------------------------------------------- 

// 3.1. Descriptive statistics of Ti
sum Ti, detail

// 3.2. Tabulate of hist4
tab hist4

// 3.3. Transition probabilities and market shares of purchasing events

xtset idconsumer conta
xttrans brand

log close 



