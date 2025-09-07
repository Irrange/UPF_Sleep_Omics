# Outcome-Wide Analysis

library(survival)
library(dplyr)
library(multcomp)

perform_cox_regression <- function(data, exposure, time_var, status_var) {

  formula_str <- paste0("Surv(", time_var, ", ", status_var, ") ~ ", exposure, " + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                        edu + tdi + diet_score + phq2 + energy")
  
  cox_model <- coxph(as.formula(formula_str), data = data)
  
  results <- summary(cox_model)
  
  return(results)
}

perform_multiple_testing_correction <- function(results, method = "fdr") {

  p_values <- sapply(results, function(x) x$coefficients[1, 5])
  

  adjusted_p_values <- p.adjust(p_values, method = method)
  

  corrected_results <- mapply(function(result, adj_p) {
    result$coefficients <- cbind(result$coefficients, adj_p_value = adj_p)
    return(result)
  }, results, adjusted_p_values, SIMPLIFY = FALSE)
  
  return(corrected_results)
}


perform_upf_outcome_analysis <- function(data, outcomes) {

  results <- lapply(outcomes, function(outcome) {

    time_var <- paste0(outcome, "_time")
    status_var <- paste0(outcome, "_status")
    
    if (!(time_var %in% names(data)) || !(status_var %in% names(data))) {
      return(NULL)
    }
    
    cox_model <- coxph(as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ upf_intake + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                          edu + tdi + diet_score + phq2 + energy")), 
                       data = data)
    
    result <- summary(cox_model)
    
    return(list(outcome = outcome, result = result))
  })
  
  results <- results[!sapply(results, is.null)]
  
  corrected_results <- perform_multiple_testing_correction(results)
  
  return(corrected_results)
}

perform_omics_outcome_analysis <- function(data, proteomic_score, metabolomic_score, outcomes) {

  proteomic_results <- lapply(outcomes, function(outcome) {

    time_var <- paste0(outcome, "_time")
    status_var <- paste0(outcome, "_status")
    
    if (!(time_var %in% names(data)) || !(status_var %in% names(data))) {
      return(NULL)
    }
    
    cox_model <- coxph(as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ proteomic_score + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                          edu + tdi + diet_score + phq2 + energy")), 
                       data = data)
    
    result <- summary(cox_model)
    
    return(list(outcome = outcome, result = result))
  })
  
  proteomic_results <- proteomic_results[!sapply(proteomic_results, is.null)]
  
  corrected_proteomic_results <- perform_multiple_testing_correction(proteomic_results)
  
  metabolomic_results <- lapply(outcomes, function(outcome) {
    time_var <- paste0(outcome, "_time")
    status_var <- paste0(outcome, "_status")
    
    cox_model <- coxph(as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ metabolomic_score + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                          edu + tdi + diet_score + phq2 + energy")), 
                       data = data)
    
    result <- summary(cox_model)
    
    return(list(outcome = outcome, result = result))
  })
  
  metabolomic_results <- metabolomic_results[!sapply(metabolomic_results, is.null)]
  
  corrected_metabolomic_results <- perform_multiple_testing_correction(metabolomic_results)
  
  results <- list(
    proteomic_results = corrected_proteomic_results,
    metabolomic_results = corrected_metabolomic_results
  )
  
  return(results)
}

perform_mortality_analysis <- function(data) {
  all_cause_mortality <- perform_cox_regression(data, "upf_intake", "death_time", "death_status")
  
  cancer_mortality <- perform_cox_regression(data, "upf_intake", "cancer_death_time", "cancer_death_status")
  
  nervous_system_mortality <- perform_cox_regression(data, "upf_intake", "nervous_death_time", "nervous_death_status")
  
  respiratory_mortality <- perform_cox_regression(data, "upf_intake", "respiratory_death_time", "respiratory_death_status")
  
  circulatory_mortality <- perform_cox_regression(data, "upf_intake", "circulatory_death_time", "circulatory_death_status")
  
  proteomic_all_cause <- perform_cox_regression(data, "proteomic_score", "death_time", "death_status")
  metabolomic_all_cause <- perform_cox_regression(data, "metabolomic_score", "death_time", "death_status")
  
  results <- list(
    all_cause_mortality = all_cause_mortality,
    cancer_mortality = cancer_mortality,
    nervous_system_mortality = nervous_system_mortality,
    respiratory_mortality = respiratory_mortality,
    circulatory_mortality = circulatory_mortality,
    proteomic_all_cause = proteomic_all_cause,
    metabolomic_all_cause = metabolomic_all_cause
  )
  
  return(results)
}


main <- function() {
  processed_data <- read.csv("path/to/processed_data.csv")
  
  outcomes <- c("cardiovascular_disease", "diabetes", "cancer", "dementia", "depression", 
                "chronic_respiratory_disease", "chronic_kidney_disease", "liver_disease",
                "parkinson_disease", "multiple_sclerosis")
  
  upf_outcome_results <- perform_upf_outcome_analysis(processed_data, outcomes)
  
  load("path/to/omics_signature_analysis_results.RData")
  
  omics_outcome_results <- perform_omics_outcome_analysis(processed_data, 
                                                          proteomic_signature$proteomic_score, 
                                                          metabolomic_signature$metabolomic_score, 
                                                          outcomes)
  
  mortality_results <- perform_mortality_analysis(processed_data)
  
  save(upf_outcome_results, omics_outcome_results, mortality_results,
       file = "path/to/outcome_wide_analysis_results.RData")
}

main()