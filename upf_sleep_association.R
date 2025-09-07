# UPF and Sleep Quality Association Analysis
library(dplyr)
library(ggplot2)
library(splines)


perform_upf_sleep_association <- function(data) {
  model1 <- lm(psqi_total ~ upf_intake, data = data)

  model2 <- lm(psqi_total ~ upf_intake + age + sex, data = data)
  
  model3 <- lm(psqi_total ~ upf_intake + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                 edu + tdi + diet_score + phq2 + energy, data = data)
  
  data_quartiles <- data %>%
    mutate(upf_quartile = cut(upf_intake, breaks = quantile(upf_intake, probs = seq(0, 1, 0.25)), 
                              include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")))
  
  quartile_model <- lm(psqi_total ~ upf_quartile + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                         edu + tdi + diet_score + phq2 + energy, data = data_quartiles)

  results <- list(
    model1 = summary(model1),
    model2 = summary(model2),
    model3 = summary(model3),
    quartile_model = summary(quartile_model)
  )
  
  return(results)
}


perform_psqi_components_analysis <- function(data) {
 
  components <- c("psqi_quality", "psqi_latency", "psqi_duration", "psqi_efficiency", 
                  "psqi_disturbance", "psqi_meds", "psqi_dysfunction")
  

  component_results <- lapply(components, function(component) {
    model <- lm(as.formula(paste0(component, " ~ upf_intake + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                  edu + tdi + diet_score + phq2 + energy")), data = data)
    return(summary(model))
  })
  
  names(component_results) <- components
  
  return(component_results)
}


perform_hrs_validation <- function(hrs_data) {

  model <- lm(jss_total ~ upf_intake + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                edu + diet_score + phq2 + energy, data = hrs_data)
  

  hrs_data_quartiles <- hrs_data %>%
    mutate(upf_quartile = cut(upf_intake, breaks = quantile(upf_intake, probs = seq(0, 1, 0.25)), 
                              include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")))
  
  quartile_model <- lm(jss_total ~ upf_quartile + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                         edu + diet_score + phq2 + energy, data = hrs_data_quartiles)
  
  results <- list(
    main_model = summary(model),
    quartile_model = summary(quartile_model)
  )
  
  return(results)
}


main <- function() {

  ukb_data <- read.csv("path/to/ukb_processed_data.csv")
  hrs_data <- read.csv("path/to/hrs_processed_data.csv")
  

  ukb_association_results <- perform_upf_sleep_association(ukb_data)
  

  psqi_components_results <- perform_psqi_components_analysis(ukb_data)
  
  hrs_validation_results <- perform_hrs_validation(hrs_data)
  

  save(ukb_association_results, psqi_components_results, 
       hrs_validation_results, meta_analysis_results, 
       file = "path/to/upf_sleep_association_results.RData")
}


main()