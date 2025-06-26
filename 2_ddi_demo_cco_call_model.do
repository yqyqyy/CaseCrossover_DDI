/*=========================================================================
DO FILE NAME:		2_ddi_demo_cco_master
AUTHOR:			Qiuyan Yu
EMAIL:			yuqiuyan@connect.hku.hk
DATE:			19 June 2025
					
Aim: 			Call the model to run case-crossover studies
*=============================================================================*/

* open log file - no need as fast tool will create log files
log using "${pathLogs}/demo_cco_master", text replace


/*********************************************************************
* Run the Case-crossover study do-file for 6 parameter-model*************
*********************************************************************/
* Assuming the outcome of interest is major bleed, which is identified in the "cco_major_bleed" dataset
* Object drug information is extracted as the "object_data"
* Precipitant drug information is extracted as the "precipitant_data"

foreach outcome in major_bleed {
	foreach riskperiod in 30 {
	
*object + precipitant
global washoutperiod 30 
global riskperiod `riskperiod'
global inputdataset "$pathIn/cco_`outcome'"
global object exp1
global object_dataset "$pathIn/object_data"
global precipitant_dataset "$pathIn/precipitant_data"
global logname $pathLogs/`riskperiod'd/cco_exp1_`outcome'
global outputtable $pathResults/`riskperiod'd/table_cco_exp1_`outcome'.txt
global outputexcel $pathResults/`riskperiod'd/cco_exp1_`outcome'_excel.txt

cap noi do "$pathDofiles/3_ddi_demo_cco_six_parameter_model.do"

}
}



/****************************************************************************
* Run the case-crossover study for 6 parameter-model with active comparator
****************************************************************************/
foreach outcome in major_bleed {
	foreach riskperiod in 60 {
	
*object + precipitant
global washoutperiod 60 
global riskperiod `riskperiod'
global inputdataset "$pathIn/cco_`outcome'"
global object exp1
global object_dataset "$pathIn/object_data"
global precipitant_dataset "$pathIn/precipitant_data"
global logname $pathLogs/`riskperiod'd/cco_exp1_`outcome'
global outputtable $pathResults/`riskperiod'd/table_cco_exp1_`outcome'.txt
global outputexcel $pathResults/`riskperiod'd/cco_exp1_`outcome'_excel.txt

cap noi do "$pathDofiles/3_ddi_demo_cco_six_parameter_model.do"


*object + negative precipitant
global washoutperiod 60
global riskperiod `riskperiod'
global inputdataset "$pathIn/cco_`outcome'"
global object exp2
global object_dataset "$pathIn/object_data"
global precipitant_dataset "$pathIn/negative_precipitant_data"
global logname $pathLogs/`riskperiod'd/cco_exp2_`outcome'
global outputtable $pathResults/`riskperiod'd/table_cco_exp2_`outcome'.txt
global outputexcel $pathResults/`riskperiod'd/cco_exp2_`outcome'_excel.txt

cap noi do "$pathDofiles/3_ddi_demo_cco_six_parameter_model.do"


*Active comparator
global logname $pathLogs/`riskperiod'd/cco_ac_`outcome'
global ac_outputtable $pathResults/`riskperiod'd/AC_`outcome'.txt

cap noi do "$pathDofiles/4_ddi_demo_cco_active_comparator.do" 


}
}



