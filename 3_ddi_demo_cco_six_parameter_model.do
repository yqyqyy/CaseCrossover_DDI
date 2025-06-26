/*=========================================================================
DO FILE NAME:		3_ddi_demo_cco_six_parameter_model
AUTHOR:			Qiuyan Yu
EMAIL:			yuqiuyan@connect.hku.hk
DATE:			19 June 2025

Aim: 			To create datasets for running case-crossover studies
Steps:
1. Set up hazard and referent periods for each individual
2. Identify exposure in hazard and referent windows
3. Generate variables for regression
4. Run Conditional logistic regression (6 parameter model)
*=============================================================================*/

capture log close

/*******************************************************************************
Identify file locations
*******************************************************************************/
* open log file - no need as fast tool will create log files
log using "$logname", text replace

/**********************************************************************
* Open the dataset
***********************************************************************/
use "$inputdataset", clear

/*1. create hazard window and control windows*/
gen period_end = eventdate
gen period_start = eventdate - $riskperiod + 1

*create 2 obs per patient and create a period variable
expand 2, generate(period)

replace period_end = eventdate - $riskperiod - $washoutperiod + 2 if period == 0
replace period_start = period_end - $riskperiod + 1  if period == 0

format period_end period_start %td

label drop _merge
label def _merge 0 "Control window" 1 "Hazard window"

sort patid period_start

*drop patients if their risk periods were not within the valid study period
gen exclude = 1 if period_start < cohortentry 
bysort patid: egen max_exclude = max(exclude)
drop if max_exclude == 1
drop exclude max_exclude


/*2. identify object and precipitant drug exposure in risk periods*/

* Object drug
joinby patid using "$object_dataset", unmatched(master)
gen object = 1 if period_start<=rxst & rxst<=period_end 
replace object = 1 if period_start<=rxen & rxen<=period_end 
replace object = 1 if rxst <=period_start & rxen>=period_end 
replace object = 0 if object == .

gsort patid period -object
bysort patid period: keep if _n == 1

drop rxst rxen _merge

* Precipitant drug
joinby patid using "$precipitant_dataset", unmatched(master)
gen precipitant = 1 if period_start<=rxst & rxst<=period_end 
replace precipitant = 1 if period_start<=rxen & rxen<=period_end 
replace precipitant = 1 if rxst <=period_start & rxen>=period_end 
replace precipitant = 0 if precipitant == .

gsort patid period -precipitant
bysort patid period: keep if _n == 1

drop rxst rxen _merge

* Generate variables in 6-parameter model
gen min_row_drug = min(object, precipitant)
gen max_row_drug = max(object, precipitant)

bysort patid: egen min_col_object = min(object)
bysort patid: egen max_col_object = max(object)

bysort patid: egen min_col_precipitant = min(precipitant)
bysort patid: egen max_col_precipitant = max(precipitant)

* Generate object drug only variable
gen objectonly = 1 if min_col_object == 0 & max_col_object == 1 & max_col_precipitant == 0
replace objectonly = 0 if objectonly == 1 & object == 0
replace objectonly = 0 if objectonly == .

* Generate precipitant drug only variable
gen preciponly = 1 if min_col_precipitant == 0 & max_col_precipitant == 1 & max_col_object == 0
replace preciponly = 0 if preciponly == 1 & precipitant == 0
replace preciponly = 0 if preciponly == .

* Generate Joint exposure variable
gen joint = 1 if min_row_drug == 1 & max_row_drug == 1 & min_col_precipitant == 0 & min_col_object == 0
replace joint = 0 if joint == .

* Generate object drug in precipitant drug exposed time variable
gen object_onprecip = 1 if min_col_precipitant == 1 & min_col_object == 0 & max_col_object == 1
replace object_onprecip = 0 if object_onprecip == 1 & object == 0
replace object_onprecip = 0 if object_onprecip == .

* Generate precipitant drug in object drug exposed time variable
gen precip_onobject = 1 if min_col_object == 1 & min_col_precipitant == 0 & max_col_precipitant == 1
replace precip_onobject = 0 if precip_onobject == 1 & precipitant == 0
replace precip_onobject = 0 if precip_onobject == .

* Generate switch variable
gen switch = 1 if min_col_object == 0 & max_col_object == 1 & ///
					min_col_precipitant == 0 & max_col_precipitant == 1 & ///
					min_row_drug == 0 & max_row_drug == 1
					
replace switch = 0 if switch == 1 & period == 0 & precipitant == 1
replace switch = 0 if switch == 1 & period == 1 & precipitant == 1
replace switch = 0 if switch == .

drop min_row_drug max_row_drug min_col_object max_col_object min_col_precipitant max_col_precipitant 

* Remove people who do not have discordant pairs of exposure between periods
bysort patid: egen max_object= max(object)
bysort patid: egen max_prep= max(precipitant)

bysort patid: egen min_object= min(object)
bysort patid: egen min_prep= min(precipitant)

drop if max_object== 0 & max_prep == 0
drop if min_object== 1 & min_prep == 1
drop if min_prep == 1 & max_object == 0
drop if min_object == 1 & max_prep == 0

drop max_object max_prep min_object min_prep


* Number of event in each exposure patterns
preserve
duplicates drop patid, force
count
local event = r(N)
restore

preserve
keep if objectonly == 1
by patid: keep if _n==1
count if period == 0
local objectonly_control = r(N)
count if period == 1
local objectonly_hazard = r(N)
restore

preserve
keep if preciponly == 1
by patid: keep if _n==1
count if period == 0
local preciponly_control = r(N)
count if period == 1
local preciponly_hazard = r(N)
restore

preserve
keep if joint == 1
by patid: keep if _n==1
count if period == 0
local joint_control = r(N)
count if period == 1
local joint_hazard = r(N)
restore

preserve
keep if object_onprecip == 1
by patid: keep if _n==1
count if period == 0
local object_onprecip_control = r(N)
count if period == 1
local object_onprecip_hazard = r(N)
restore

preserve
keep if precip_onobject == 1
by patid: keep if _n==1
count if period == 0
local precip_onobject_control = r(N)
count if period == 1
local precip_onobject_hazard = r(N)
restore

preserve
keep if switch == 1
by patid: keep if _n==1
count if period == 0
local switch_control = r(N)
count if period == 1
local switch_hazard = r(N)
restore


* Run 6-parameter model
clogit period i.objectonly i.preciponly i.joint i.object_onprecip i.precip_onobject i.switch, ///
 group(patid) or level(99)

* Output to text file
cap file close tablecontent
file open tablecontent using "$outputtable", write text replace

file write tablecontent _tab ("Object only") _tab _tab ///
							("Precipitant only") _tab _tab ///
							("Joint") _tab _tab ///
							("Object while on precipitant") _tab _tab ///
							("Precipitant while on object") _tab _tab ///
							("Switch") _n
						
file write tablecontent ("$exp_interest") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab /// 
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _n
cap {						
lincom 1.objectonly, eform level(99)				
file write tablecontent _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _tab 
}

cap {	
lincom 1.preciponly, eform level(99)			
file write tablecontent %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _tab 
}

cap {	
lincom 1.joint, eform level(99)
file write tablecontent %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _tab 
}

cap {	
lincom 1.object_onprecip, eform level(99)
file write tablecontent %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _tab
}

cap {	
lincom 1.precip_onobject, eform level(99)
file write tablecontent %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _tab
}

cap {	
lincom 1.switch, eform level(99)
file write tablecontent %4.2f (r(estimate)) _tab %4.2f (r(lb)) (" - ") %4.2f (r(ub)) _n
}

file write tablecontent _tab ("objectonly_hazard") _tab ("objectonly_control") _tab ("preciponly_hazard") ///
						_tab ("preciponly_control") _tab ("joint_hazard") _tab ("joint_control") ///
						_tab ("object_onprecip_hazard") _tab ("object_onprecip_control") _tab ("precip_onobject_hazard") ///
						_tab ("precip_onobject_control") _tab ("switch_hazard") _tab ("switch_control") _n
						
file write tablecontent _tab (`objectonly_hazard') _tab (`objectonly_control') _tab (`preciponly_hazard') ///
						_tab (`preciponly_control') _tab (`joint_hazard') _tab (`joint_control') ///
						_tab (`object_onprecip_hazard') _tab (`object_onprecip_control') _tab (`precip_onobject_hazard') ///
						_tab (`precip_onobject_control') _tab (`switch_hazard') _tab (`switch_control') _n


file write tablecontent _tab ("Total N") _n
file write tablecontent _tab (`event')

file write tablecontent _n
file close tablecontent


* Output for excel file for running meta-analysis
cap file close tablecontent
file open tablecontent using "$outputexcel", write text replace
						
file write tablecontent ("$exp_interest") _tab ///
						 ("OR") _tab ("99% CI") _n
cap {						
lincom 1.objectonly, eform level(99)				
file write tablecontent ("Object only") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n 
scalar objectonly_$object = r(estimate)
}

cap {	
lincom 1.preciponly, eform level(99)			
file write tablecontent ("Precipitant only") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n
scalar preciponly_$object = r(estimate)
}

cap {	
lincom 1.joint, eform level(99)
file write tablecontent ("Joint") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n
scalar joint_$object = r(estimate)
}

cap {	
lincom 1.object_onprecip, eform level(99)
file write tablecontent ("Object while on precipitant") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n
scalar onprecip_$object = r(estimate)
}

cap {	
lincom 1.precip_onobject, eform level(99)
file write tablecontent ("Precipitant while on object") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n
scalar onobject_$object = r(estimate)
}

cap {	
lincom 1.switch, eform level(99)
file write tablecontent ("Switch") _tab %4.2f (r(estimate)) _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _n
scalar switch_$object = r(estimate)
}

file write tablecontent _n
file close tablecontent


*Estimating SE for active comparator analysis

clogit period i.objectonly i.preciponly i.joint i.object_onprecip i.precip_onobject i.switch, group(patid) level(99)
cap {						
lincom 1.objectonly, level(99)
scalar se_objectonly_$object = r(se)
}

cap {	
lincom 1.preciponly, level(99)		
scalar se_preciponly_$object = r(se)
}

cap {	
lincom 1.joint, level(99)
scalar se_joint_$object = r(se)
}

cap {	
lincom 1.object_onprecip, level(99)
scalar se_onprecip_$object = r(se)
}

cap {	
lincom 1.precip_onobject, level(99)
scalar se_onobject_$object = r(se) 
}

cap {	
lincom 1.switch, level(99)
scalar se_switch_$object = r(se)
}


log close
