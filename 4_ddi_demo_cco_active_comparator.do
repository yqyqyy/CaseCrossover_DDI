/*=========================================================================
DO FILE NAME:		4_ddi_demo_cco_active_comparator
AUTHOR:			Qiuyan Yu
EMAIL: 			yuqiuyan@connect.hku.hk
DATE:			19 June 2025

Aim: 			Compute adjusted estimate for active comparator 
			design of 6 parameters model

Note:			Confidence interval (CI) is 99% with z=2.575, 
			adjust the z value for corresponding CI as needed.
*=============================================================================*/

capture log close

/*******************************************************************************
Identify file locations
*******************************************************************************/
* open log file - no need as fast tool will create log files
log using "$logname", text replace

* Output to text file
cap file close tablecontent
file open tablecontent using "$ac_outputtable", write text replace

*Simple ratio of point estimate
local objectonly = objectonly_exp1 / objectonly_exp2
local preciponly = preciponly_exp1 / preciponly_exp2
local joint = joint_exp1 / joint_exp2
local object_onprecip = onprecip_exp1 / onprecip_exp2
local precip_onobject = onobject_exp1 / onobject_exp2
local switch = switch_exp1 / switch_exp2

*Simple ratio: SE
scalar se_objectonly = sqrt(se_objectonly_exp1 ^2 + se_objectonly_exp2 ^2)
local objectonly_upper =  exp((log(`objectonly') + 2.575 * se_objectonly))
local objectonly_lower =  exp((log(`objectonly') - 2.575 * se_objectonly))

scalar se_preciponly = sqrt(se_preciponly_exp1 ^2 + se_preciponly_exp2 ^2)
local preciponly_upper = exp((log(`preciponly') + 2.575 * se_preciponly))
local preciponly_lower = exp((log(`preciponly') - 2.575 * se_preciponly))

scalar se_joint = sqrt(se_joint_exp1 ^2 + se_joint_exp2 ^2)
local joint_upper = exp((log(`joint') + 2.575 * se_joint))
local joint_lower = exp((log(`joint') - 2.575 * se_joint))

scalar se_onprecip = sqrt(se_onprecip_exp1 ^2 + se_onprecip_exp2 ^2)
local object_onprecip_upper = exp((log(`object_onprecip') + 2.575 * se_onprecip))
local object_onprecip_lower = exp((log(`object_onprecip') - 2.575 * se_onprecip))

scalar se_onobject = sqrt(se_onobject_exp1 ^2 + se_onobject_exp2 ^2)
local precip_onobject_upper = exp((log(`precip_onobject') + 2.575 * se_onobject))
local precip_onobject_lower = exp((log(`precip_onobject') - 2.575 * se_onobject))

scalar se_switch = sqrt(se_switch_exp1 ^2 + se_switch_exp2 ^2)
local switch_upper = exp((log(`switch') + 2.575 * se_switch))
local switch_lower = exp((log(`switch') - 2.575 * se_switch))

file write tablecontent _tab ("Object only") _tab _tab ///
							("Precipitant only") _tab _tab ///
							("Joint") _tab _tab ///
							("Object while on precipitant") _tab _tab ///
							("Precipitant while on object") _tab _tab ///
							("Switch") _n
						
file write tablecontent ("Active comparator") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab /// 
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _tab ///
						 ("OR") _tab ("99% CI") _n


file write tablecontent %4.2f _tab (`objectonly') _tab %4.2f (`objectonly_lower') (" - ") %4.2f (`objectonly_upper') _tab

file write tablecontent %4.2f (`preciponly') _tab %4.2f (`preciponly_lower') (" - ") %4.2f (`preciponly_upper') _tab				
file write tablecontent %4.2f (`joint') _tab %4.2f (`joint_lower') (" - ") %4.2f (`joint_upper') _tab

file write tablecontent %4.2f (`object_onprecip') _tab %4.2f (`object_onprecip_lower') (" - ") %4.2f (`object_onprecip_upper') _tab

file write tablecontent %4.2f (`precip_onobject') _tab %4.2f (`precip_onobject_lower') (" - ") %4.2f (`precip_onobject_upper') _tab

file write tablecontent %4.2f (`switch') _tab %4.2f (`switch_lower') (" - ") %4.2f (`switch_upper') _n

file write tablecontent _n
file close tablecontent


log close
