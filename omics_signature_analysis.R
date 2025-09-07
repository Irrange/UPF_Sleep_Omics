# Multi-omics Signature Construction and Analysis

library(glmnet)
library(caret)
library(doParallel)
library(foreach)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)


perform_elastic_net <- function(x, y, nfolds = 10) {

  ctrl <- trainControl(method = "cv", number = nfolds)
  
  tune_grid <- expand.grid(
    alpha = seq(0, 1, by = 0.1),
    lambda = seq(0.001, 0.1, length.out = 10)
  )
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  elastic_net_model <- train(
    x = x,
    y = y,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = tune_grid,
    standardize = TRUE
  )
  
  stopCluster(cl)
  
  return(elastic_net_model)
}


perform_univariate_analysis <- function(omics_data, upf_intake, p_threshold = 0.05) {

  n_tests <- ncol(omics_data)
  bonferroni_threshold <- p_threshold / n_tests
  
  univariate_results <- lapply(names(omics_data), function(molecule) {
    model <- lm(as.formula(paste0(molecule, " ~ upf_intake")), data = data.frame(upf_intake, omics_data[molecule]))
    p_value <- summary(model)$coefficients[2, 4]
    beta <- summary(model)$coefficients[2, 1]
    return(list(molecule = molecule, p_value = p_value, beta = beta))
  })
  
  significant_results <- do.call(rbind, lapply(univariate_results, function(x) {
    if (x$p_value < bonferroni_threshold) {
      return(data.frame(molecule = x$molecule, p_value = x$p_value, beta = x$beta))
    }
  }))
  

  significant_results <- significant_results[!sapply(significant_results, is.null), ]
  
  return(significant_results)
}


construct_proteomic_signature <- function(protein_data, upf_intake) {

  univariate_results <- perform_univariate_analysis(protein_data, upf_intake, p_threshold = 0.05)
  
  significant_proteins <- protein_data[, univariate_results$molecule]
  

  x <- as.matrix(significant_proteins)
  y <- scale(upf_intake)
  
  elastic_net_model <- perform_elastic_net(x, y)
  

  coef <- coef(elastic_net_model$finalModel, s = elastic_net_model$bestTune$lambda)
  significant_vars <- names(coef@x)[coef@i + 1]
  
  correlation <- cor(proteomic_score, upf_intake, use = "complete.obs")
  
  results <- list(
    univariate_results = univariate_results,
    model = elastic_net_model,
    significant_vars = significant_vars,
    proteomic_score = as.numeric(proteomic_score),
    correlation = correlation
  )
  
  return(results)
}


construct_metabolomic_signature <- function(metabolite_data, upf_intake) {

  univariate_results <- perform_univariate_analysis(metabolite_data, upf_intake, p_threshold = 0.05)
  

  significant_metabolites <- metabolite_data[, univariate_results$molecule]
  

  x <- as.matrix(significant_metabolites)
  y <- scale(upf_intake)
  
  elastic_net_model <- perform_elastic_net(x, y)
  

  coef <- coef(elastic_net_model$finalModel, s = elastic_net_model$bestTune$lambda)
  significant_vars <- names(coef@x)[coef@i + 1]

  correlation <- cor(metabolomic_score, upf_intake, use = "complete.obs")
  
  results <- list(
    univariate_results = univariate_results,
    model = elastic_net_model,
    significant_vars = significant_vars,
    metabolomic_score = as.numeric(metabolomic_score),
    correlation = correlation
  )
  
  return(results)
}


perform_functional_enrichment <- function(significant_proteins) {

  gene_symbols <- significant_proteins
  
  kegg_enrichment <- enrichKEGG(gene = gene_symbols, organism = "hsa", pvalueCutoff = 0.05)
  
  go_enrichment <- enrichGO(gene = gene_symbols, 
                            OrgDb = org.Hs.eg.db, 
                            keyType = "SYMBOL", 
                            ont = "BP", 
                            pvalueCutoff = 0.05)
  
  results <- list(
    kegg_enrichment = kegg_enrichment,
    go_enrichment = go_enrichment
  )
  
  return(results)
}

perform_omics_sleep_association <- function(data, proteomic_score, metabolomic_score) {

  proteomic_model <- lm(psqi_total ~ proteomic_score + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                          edu + tdi + diet_score + phq2 + energy, data = data)
  

  data_quartiles <- data %>%
    mutate(proteomic_quartile = cut(proteomic_score, breaks = quantile(proteomic_score, probs = seq(0, 1, 0.25)), 
                                    include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")))
  
  proteomic_quartile_model <- lm(psqi_total ~ proteomic_quartile + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                   edu + tdi + diet_score + phq2 + energy, data = data_quartiles)
  

  proteomic_spline_model <- lm(psqi_total ~ ns(proteomic_score, df = 3) + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                 edu + tdi + diet_score + phq2 + energy, data = data)
  

  metabolomic_model <- lm(psqi_total ~ metabolomic_score + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                            edu + tdi + diet_score + phq2 + energy, data = data)
  

  data_quartiles$metabolomic_quartile <- cut(metabolomic_score, breaks = quantile(metabolomic_score, probs = seq(0, 1, 0.25)), 
                                             include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
  
  metabolomic_quartile_model <- lm(psqi_total ~ metabolomic_quartile + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                     edu + tdi + diet_score + phq2 + energy, data = data_quartiles)
  

  metabolomic_spline_model <- lm(psqi_total ~ ns(metabolomic_score, df = 3) + age + sex + race1 + drink1 + smoke1 + bmi1 + met_cat +
                                   edu + tdi + diet_score + phq2 + energy, data = data)
  
  results <- list(
    proteomic_model = summary(proteomic_model),
    proteomic_quartile_model = summary(proteomic_quartile_model),
    proteomic_spline_model = summary(proteomic_spline_model),
    metabolomic_model = summary(metabolomic_model),
    metabolomic_quartile_model = summary(metabolomic_quartile_model),
    metabolomic_spline_model = summary(metabolomic_spline_model)
  )
  
  return(results)
}


main <- function() {

  processed_data <- read.csv("path/to/processed_data.csv")
  protein_data <- read.csv("path/to/protein_data.csv")
  metabolite_data <- read.csv("path/to/metabolite_data.csv")
  

  proteomic_signature <- construct_proteomic_signature(protein_data, processed_data$upf_intake)
  

  metabolomic_signature <- construct_metabolomic_signature(metabolite_data, processed_data$upf_intake)
  

  functional_enrichment <- perform_functional_enrichment(proteomic_signature$significant_vars)
  

  omics_sleep_association <- perform_omics_sleep_association(processed_data, 
                                                             proteomic_signature$proteomic_score, 
                                                             metabolomic_signature$metabolomic_score)
  

  save(proteomic_signature, metabolomic_signature, functional_enrichment, omics_sleep_association,
       file = "path/to/omics_signature_analysis_results.RData")
}


main()