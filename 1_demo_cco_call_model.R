##########################################################################
# Author:   Qiuyan Yu
# Email:    yuqiuyan@connect.hku.hk
# Date:     19 June 2025
# Purpose:  To create datasets for running 6-parameter model case-crossover 
#           study using conditional logistic regression 
##########################################################################


# 0. Packages loading ----
library(data.table)
library(dplyr)
library(survival)



# 1. Load the cases and drugs datasets ----
# Assuming you have extracted all the outcome events as the case data, with "eventdate" indicating the outcome occurrence date, "cohortentry" and "cohortend" indicating the study start and end date for each patient. 
case_data <- readRDS("case.rds") 

# Assuming you have extracted all the object drug and precipitant drug as the separate data, with "rxst", "rxen", "drug" indicating the prescription start date, prescription end date, and drugname. 
object_data <- readRDS("object.rds")
precipitant_data <- readRDS("precipitant.rds")
negative_precipitant_data <- readRDS("precipitant.rds")

# Convert to data.table for efficient processing
setDT(case_data)
setDT(object_data)
setDT(precipitant_data)
setDT(negative_precipitant_data)

# Make sure all the date variables are in date format
case_data[, `:=`(eventdate = as.Date(eventdate, format = "%d%b%Y"), cohortentry = as.Date(cohortentry, format = "%d%b%Y"), cohortend = as.Date(cohortend, format = "%d%b%Y"))]
object_data[, `:=`(rxst = as.Date(rxst, format = "%d%b%Y"), rxen = as.Date(rxen, format = "%d%b%Y"))]
precipitant_data[, `:=`(rxst = as.Date(rxst, format = "%d%b%Y"), rxen = as.Date(rxen, format = "%d%b%Y"))]
negative_precipitant_data[, `:=`(rxst = as.Date(rxst, format = "%d%b%Y"), rxen = as.Date(rxen, format = "%d%b%Y"))]



# 2. Run the function to generate the dataset for 6-parameter model ----
cco_6paramdl_fun <- function(case_dt=case_data, object_dt=object_data, precipitant_dt=precipitant_data, riskperiod=30, washoutperiod=30) {
  
  # 1. Create hazard and referent windows (2 periods per patient)
  case_dt[, period_end := eventdate]
  case_dt[, period_start := eventdate - riskperiod + 1]
  
  # Expand dataset to 2 rows per patient for hazard (period=1) and referent (period=0) window
  case_expanded <- case_dt[rep(1:.N, each=2)]
  case_expanded[, period := rep(c(0,1), times=.N/2)]
  
  # Modify referent window dates
  case_expanded[period == 0, period_end := eventdate - riskperiod - washoutperiod + 2]
  case_expanded[period == 0, period_start := period_end - riskperiod + 1]
  
  # Filter patients whose risk periods are within valid study period
  case_expanded <- case_expanded[period_start >= cohortentry]
  
  
  
  # 2. Identify object and precipitant drug exposure in risk periods
  # Function to flag exposure during period
  flag_exposure <- function(dt, exposure_data, exposure_name) {
    # Merge exposure data by patient id
    dt <- merge(dt, exposure_data, by="patid", all.x=TRUE, suffixes=c("", "_exp"))
    
    # Flag exposure if any overlap between rxst-rxen and period_start-period_end
    dt[, (exposure_name) := 0]
    dt[!is.na(rxst) & (
      (rxst >= period_start & rxst <= period_end) |
        (rxen >= period_start & rxen <= period_end) |
        (rxst <= period_start & rxen >= period_end)
    ), (exposure_name) := 1]
    
    # Keep only one row per "patid" and "period" (first exposure record)
    dt <- dt[order(patid, period, -get(exposure_name))]
    dt <- dt[, .SD[1], by=.(patid, period)]
    
    # Drop exposure prescription start and prescription end date columns to avoid confusion
    dt[, c("rxst", "rxen") := NULL]
    
    return(dt)
   }
  
  # Apply exposure flagging functions
  case_expanded <- flag_exposure(case_expanded, object_dt, "object")
  case_expanded <- flag_exposure(case_expanded, precipitant_dt, "precipitant")
  
  
  
  # 3. Generate variables for 6-parameter model
  case_expanded[, min_row_drug := pmin(object, precipitant)]
  case_expanded[, max_row_drug := pmax(object, precipitant)]
  
  case_expanded[, min_col_object := min(object), by=patid]
  case_expanded[, max_col_object := max(object), by=patid]
  
  case_expanded[, min_col_precipitant := min(precipitant), by=patid]
  case_expanded[, max_col_precipitant := max(precipitant), by=patid]
  
  
  # Define 6 exposure initiation pattern variables
  case_expanded[, objectonly := fifelse(min_col_object == 0 & max_col_object == 1 & max_col_precipitant == 0 & object == 1, 1, 0)]
  case_expanded[, preciponly := fifelse(min_col_precipitant == 0 & max_col_precipitant == 1 & max_col_object == 0 & precipitant == 1, 1, 0)]
  case_expanded[, joint := fifelse(min_row_drug == 1 & max_row_drug == 1 & min_col_precipitant == 0 & min_col_object == 0, 1, 0)]
  case_expanded[, object_onprecip := fifelse(min_col_precipitant == 1 & min_col_object == 0 & max_col_object == 1 & object == 1, 1, 0)]
  case_expanded[, precip_onobject := fifelse(min_col_object == 1 & min_col_precipitant == 0 & max_col_precipitant == 1 & precipitant == 1, 1, 0)]
  
  case_expanded[, switch := fifelse(
    min_col_object == 0 & max_col_object == 1 &
      min_col_precipitant == 0 & max_col_precipitant == 1 &
      min_row_drug == 0 & max_row_drug == 1, 1, 0)]
  case_expanded[period == 0 & precipitant == 1 & switch == 1, switch := 0]
  case_expanded[period == 1 & precipitant == 1 & switch == 1, switch := 0]
  
  
  # Drop patients without discordant exposure pairs between periods
  case_expanded[, max_object := max(object), by=patid]
  case_expanded[, max_prep := max(precipitant), by=patid]
  case_expanded[, min_object := min(object), by=patid]
  case_expanded[, min_prep := min(precipitant), by=patid]
  
  case_expanded <- case_expanded[!(max_object == 0 & max_prep == 0)]
  case_expanded <- case_expanded[!(min_object == 1 & min_prep == 1)]
  case_expanded <- case_expanded[!(min_prep == 1 & max_object == 0)]
  case_expanded <- case_expanded[!(min_object == 1 & max_prep == 0)]
  
  return(case_expanded)
  
}

# Adjust the riskperiod, washoutperiod as needed specific to the study design
case_strata_exp1 <- cco_6paramdl_fun(case_dt=case_data, object_dt=object_data, precipitant_dt=precipitant_data, riskperiod=30, washoutperiod=30)

# If an active comparator is applied, run the function for negative precipitant drug using same conditions as precipitant drug
case_strata_exp2 <- cco_6paramdl_fun(case_dt=case_data, object_dt=object_data, precipitant_dt=negative_precipitant_data, riskperiod=30, washoutperiod=30)



# 3. Number of patients in each stratum per patient per period ----
count_cases_fun <- function(data=case_strata_exp1) {
  
  counts <- data[, .(
    objectonly_hazard = sum(objectonly == 1 & period == 1),
    objectonly_referent = sum(objectonly == 1 & period == 0),
    
    preciponly_hazard = sum(preciponly == 1 & period == 1),
    preciponly_referent = sum(preciponly == 1 & period == 0),
    
    joint_hazard = sum(joint == 1 & period == 1),
    joint_referent = sum(joint == 1 & period == 0),
    
    object_onprecip_hazard = sum(object_onprecip == 1 & period == 1),
    object_onprecip_referent = sum(object_onprecip == 1 & period == 0),
    
    precip_onobject_hazard = sum(precip_onobject == 1 & period == 1),
    precip_onobject_referent = sum(precip_onobject == 1 & period == 0),
    
    switch_hazard = sum(switch == 1 & period == 1),
    switch_referent = sum(switch == 1 & period == 0),
    
    total_cases = uniqueN(patid[period == 1])
  )] 
  
  # Sum across patients
  total_counts <- counts[, lapply(.SD, sum)]
  return(total_counts)
  
}

case_count_exp1 <- count_cases_fun(case_strata_exp1)
case_count_exp2 <- count_cases_fun(case_strata_exp2) # active comparator



# 4. Run the function to export results for the case-crossover analysis ----
cco_clogit_fun <- function(data=case_strata_exp1, conflevel=0.99, output_path, exp="exp1") {
  
  # 1. Run conditional logistic regression for the 6-parameter model using "clogit" function in "survival" package
  # Convert period to binary outcome: 1 = hazard, 0 = referent
  data[, period := as.integer(period)]
  
  # Prepare formula for 6-parameter model
  formula <- as.formula("period ~ objectonly + preciponly + joint + object_onprecip + precip_onobject + switch + strata(patid)")
  
  # Fit model
  model <- clogit(formula, data=data)
  
  # Summary with exponentiated coefficients (odds ratios)
  summary(model)
  exp(coef(model))
  
  
  
  # 2. Extract estimates and confidence intervals 
  confint_int <- confint(model, level=conflevel) 
  or <- exp(coef(model)) # odds ratios 
  ci_lower <- exp(confint_int[,1]) # lower CI
  ci_upper <- exp(confint_int[,2]) # upper CI
  se <- summary(model)$coef[,3]
  
  results <- data.table(
    Exposure = c("objectonly", "preciponly", "joint", "object_onprecip", "precip_onobject", "switch"),
    OR = or,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    SE = se
  )

  results_fmt <- results[, .(Exposure, OR = round(OR, 2), CI_lower = round(CI_lower, 2), CI_upper = round(CI_upper, 2), SE = round(SE, 2))]

  
  
  # 3. Save results to file (adjust path as needed)
  write.table(results_fmt, file=file.path(output_path, paste0("/table_cco_", exp, ".txt")), 
              sep="\t", row.names=FALSE, quote=FALSE)
  
  
  return(results)
  
}

# Adjust the riskperiod, washoutperiod, conflevel, output_path as needed specific to the study design
final_rs_exp1 <- cco_clogit_fun(data=case_strata_exp1, conflevel=0.99, output_path="C:/CCO", exp="exp1")
final_rs_exp2 <- cco_clogit_fun(data=case_strata_exp2, conflevel=0.99, output_path="C:/CCO", exp="exp2") # active comparator



# 5. Wald test to test difference between strata ----
wald_test_fun <- function(data=case_strata_exp1, output_path) {
  # Fit model
  model <- clogit(period ~ objectonly + preciponly + joint + object_onprecip + precip_onobject + switch + strata(patid), data=data)
  
  # Using Wald test with linearHypothesis from car package
  library(car)
  
  test_precip <- linearHypothesis(model, "preciponly = precip_onobject")
  test_object <- linearHypothesis(model, "objectonly = object_onprecip")
  test_precijoint <- linearHypothesis(model, "preciponly = joint")
  test_objectjoint <- linearHypothesis(model, "objectonly = joint")
  
  # Extract p-values
  preci_test_p <- test_precip$`Pr(>Chisq)`[2]
  object_test_p <- test_object$`Pr(>Chisq)`[2]
  precijoint_test_p <- test_precijoint$`Pr(>Chisq)`[2]
  objectjoint_test_p <- test_objectjoint$`Pr(>Chisq)`[2]
  
  # Output to text file 
  fileConn <- file(output_path, "w")
  writeLines(c(
    paste("Test object only and object while on precipitant parameters:\t", format(object_test_p, digits=4)),
    paste("Test precipitant only and precipitant while on object parameters:\t", format(preci_test_p, digits=4)),
    paste("Test object only and joint parameters:\t", format(objectjoint_test_p, digits=4)),
    paste("Test precipitant only and joint parameters:\t", format(precijoint_test_p, digits=4)),
    ""
  ), fileConn)
  close(fileConn)
  
}

wald_test_fun(data=case_strata_exp1, output_path="C:/CCO/test_exp1.txt")
wald_test_fun(data=case_strata_exp2, output_path="C:/CCO/test_exp2.txt")



# -------------------------------------------------------------------------


