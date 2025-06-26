/*=========================================================================
DO FILE NAME:			1_ddi_demo_paths.do
AUTHOR:				Qiuyan Yu
EMAIL: 				yuqiuyan@connect.hku.hk
DATE: 				19 June 2025

DESCRIPTION OF FILE:		Global macros for file paths
*=========================================================================*/

/*******************************************************************************
# DO FILES
*******************************************************************************/
clear 
set max_memory .

cd
cd "Stata Project/Demo_DDI" 


* Posted dofiles
global pathDofiles	 		"dofiles"


/*******************************************************************************
# DATASETS
*******************************************************************************/
global pathOut			"output_datafiles"	
global pathIn 			"input_datafiles"


/******************************************************************************* 
# OUTPUT
*******************************************************************************/
* log files
global pathLogs "logfiles"

* results text files
global pathResults "writeup"


