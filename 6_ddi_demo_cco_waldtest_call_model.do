/*=========================================================================
DO FILE NAME:		6_ddi_demo_cco_waldtest_call_model
AUTHOR:			Qiuyan Yu
EMAIL:			yuqiuyan@connect.hku.hk
DATE:			19 June 2025

Aim: 			Call model to run the wald test 
*=============================================================================*/


/*********************************************************************
* Run the Case-crossover study Wald test do-file for 6 parameter-model*************
*********************************************************************/

foreach outcome in major_bleed {
	foreach riskperiod in 30 {
	
*object + precipitant
global washoutperiod 30
global riskperiod `riskperiod'
global inputdataset "$pathIn/cco_`outcome'"
global object exp1
global object_dataset "$pathIn/object_data"
global precipitant_dataset "$pathIn/precipitant_data"
global outputtable $pathResults/cco_test/`riskperiod'd/table_cco_exp1_`outcome'.txt
global outputexcel $pathResults/cco_test/`riskperiod'd/cco_exp1_`outcome'_excel.txt

cap noi do "$pathDofiles/5_ddi_demo_cco_six_parameter_model_test.do"


*object + negative precipitant
global washoutperiod 30
global riskperiod `riskperiod'
global inputdataset "$pathIn/cco_`outcome'"
global object exp2
global object_dataset "$pathIn/object_data"
global precipitant_dataset "$pathIn/negative_precipitant_data"
global outputtable $pathResults/cco_test/`riskperiod'd/table_cco_exp2_`outcome'.txt
global outputexcel $pathResults/cco_test/`riskperiod'd/cco_exp2_`outcome'_excel.txt

cap noi do "$pathDofiles/5_ddi_demo_cco_six_parameter_model_test.do"

}
}




