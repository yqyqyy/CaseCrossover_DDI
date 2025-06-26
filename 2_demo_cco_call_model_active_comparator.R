##########################################################################
# Author:  Qiuyan Yu
# Emial:   yuqiuyan@connect.hku.hk
# Date:    19 June 2025
# Purpose: To run the active comparator analysis for case-crossover study
##########################################################################


# 0. Packages loading ----
library(data.table)
library(dplyr)
library(survival)



# 1. Simple ratio of point estimate ----
simple_ratio_fun <- function(z=2.575, output_path) {
  
  # Now you will have estimate of OR and SE for each stratum of each exposure groups (object+precipitant and object+active comparator)
  # Assign estimate
  objectonly_exp1 <- final_rs_exp1[Exposure=="objectonly", OR]
  preciponly_exp1 <- final_rs_exp1[Exposure=="preciponly", OR]
  joint_exp1 <- final_rs_exp1[Exposure=="joint", OR]
  object_onprecip_exp1 <- final_rs_exp1[Exposure=="object_onprecip", OR]
  precip_onobject_exp1 <- final_rs_exp1[Exposure=="precip_onobject", OR]
  switch_exp1 <- final_rs_exp1[Exposure=="switch", OR]
  
  objectonly_exp2 <- final_rs_exp2[Exposure=="objectonly", OR]
  preciponly_exp2 <- final_rs_exp2[Exposure=="preciponly", OR]
  joint_exp2 <- final_rs_exp2[Exposure=="joint", OR]
  object_onprecip_exp2 <- final_rs_exp2[Exposure=="object_onprecip", OR]
  precip_onobject_exp2 <- final_rs_exp2[Exposure=="precip_onobject", OR]
  switch_exp2 <- final_rs_exp2[Exposure=="switch", OR]
  
  
  # Assign SE
  se_objectonly_exp1 <- final_rs_exp1[Exposure=="objectonly", SE]
  se_preciponly_exp1 <- final_rs_exp1[Exposure=="preciponly", SE]
  se_joint_exp1 <- final_rs_exp1[Exposure=="joint", SE]
  se_object_onprecip_exp1 <- final_rs_exp1[Exposure=="object_onprecip", SE]
  se_precip_onobject_exp1 <- final_rs_exp1[Exposure=="precip_onobject", SE]
  se_switch_exp1 <- final_rs_exp1[Exposure=="switch", SE]
  
  se_objectonly_exp2 <- final_rs_exp2[Exposure=="objectonly", SE]
  se_preciponly_exp2 <- final_rs_exp2[Exposure=="preciponly", SE]
  se_joint_exp2 <- final_rs_exp2[Exposure=="joint", SE]
  se_object_onprecip_exp2 <- final_rs_exp2[Exposure=="object_onprecip", SE]
  se_precip_onobject_exp2 <- final_rs_exp2[Exposure=="precip_onobject", SE]
  se_switch_exp2 <- final_rs_exp2[Exposure=="switch", SE]
  
  
  # Compute estimates for each stratum (ratios)
  objectonly <- objectonly_exp1 / objectonly_exp2
  preciponly <- preciponly_exp1 / preciponly_exp2
  joint <- joint_exp1 / joint_exp2
  object_onprecip <- object_onprecip_exp1 / object_onprecip_exp2
  precip_onobject <- precip_onobject_exp1 / precip_onobject_exp2
  switch <- switch_exp1 / switch_exp2
  
  
  # Calculate combined standard errors
  se_objectonly <- sqrt(se_objectonly_exp1^2 + se_objectonly_exp2^2)
  se_preciponly <- sqrt(se_preciponly_exp1^2 + se_preciponly_exp2^2)
  se_joint <- sqrt(se_joint_exp1^2 + se_joint_exp2^2)
  se_onprecip <- sqrt(se_object_onprecip_exp1^2 + se_object_onprecip_exp2^2)
  se_onobject <- sqrt(se_precip_onobject_exp1^2 + se_precip_onobject_exp2^2)
  se_switch <- sqrt(se_switch_exp1^2 + se_switch_exp2^2)
  
  
  # z: z-score for confidence interval (default 2.575 for 99%)
  
  # Calculate confidence intervals by exponentiating log scale +/- z*SE
  calc_ci <- function(est, se) {
    lower <- exp(log(est) - z * se)
    upper <- exp(log(est) + z * se)
    c(lower = lower, upper = upper)
  }
  
  ci_objectonly <- calc_ci(objectonly, se_objectonly)
  ci_preciponly <- calc_ci(preciponly, se_preciponly)
  ci_joint <- calc_ci(joint, se_joint)
  ci_object_onprecip <- calc_ci(object_onprecip, se_onprecip)
  ci_precip_onobject <- calc_ci(precip_onobject, se_onobject)
  ci_switch <- calc_ci(switch, se_switch)
  
  
  # Create a data table to store results
  results <- data.table(
    Exposure = c("Object only", "Precipitant only", "Joint", 
                 "Object while on precipitant", "Precipitant while on object", "Switch"),
    OR = c(objectonly, preciponly, joint, object_onprecip, precip_onobject, switch),
    CI_lower = c(ci_objectonly["lower"], ci_preciponly["lower"], ci_joint["lower"],
                 ci_object_onprecip["lower"], ci_precip_onobject["lower"], ci_switch["lower"]),
    CI_upper = c(ci_objectonly["upper"], ci_preciponly["upper"], ci_joint["upper"],
                 ci_object_onprecip["upper"], ci_precip_onobject["upper"], ci_switch["upper"])
  )

  results_fmt <- results[, .(Exposure, OR = round(OR, 2), CI_lower = round(CI_lower, 2), CI_upper = round(CI_upper, 2))]

  
  write.table(results_fmt, file=file.path(output_path, paste0("/table_cco_ac.txt")), sep="\t", row.names=FALSE, quote=FALSE)

  return(results)  
}

ac_rs <- simple_ratio_fun(z = 2.575, output_path = "C:/CCO")



# -------------------------------------------------------------------------


