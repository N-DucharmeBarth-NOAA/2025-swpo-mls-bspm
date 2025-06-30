# HMC Samples Condensation Function
# Extracts only the parameters block variables for different Stan models

library(data.table)

#' Get parameters block variables based on Stan executable
get_stan_parameters <- function(exec_name) {
  if (grepl("bspm_estq_softdep_mvprior", exec_name)) {
    return(c("raw_mv_params", "raw_logqeff", "raw_qdev_params", 
             "raw_logsigmap", "raw_sigmao_add", "raw_epsp", 
             "raw_qdev_period", "raw_edev", "lp__"))
  } else if (grepl("bspm_estq_flex", exec_name)) {
    return(c("raw_mv_params", "raw_logqeff", "raw_qdev_params", 
             "raw_logsigmap", "raw_sigmao_add", "raw_epsp", 
             "raw_qdev_period", "raw_edev", "lp__"))
  } else if (grepl("bspm_estF", exec_name)) {
    return(c("raw_logK", "raw_logr", "raw_logshape", "raw_logsigmap", 
             "raw_sigmao_add", "raw_epsp", "raw_sigmaf", "raw_F", "lp__"))
  } else {
    return(c("raw_mv_params", "raw_logqeff", "raw_qdev_params", 
             "raw_logsigmap", "raw_sigmao_add", "raw_epsp", 
             "raw_qdev_period", "raw_edev", "raw_logK", "raw_logr", 
             "raw_logshape", "raw_sigmaf", "raw_F", "lp__"))
  }
}

#' Condense HMC samples by keeping only Stan parameters block variables
#' @param model_dir Directory containing hmc_samples.csv, fit_summary.csv, stan_data.csv
#' @return Condensed samples data.table
condense_hmc_samples <- function(model_dir) {
  
  # Check required files exist
  hmc_path <- file.path(model_dir, "hmc_samples.csv")
  summary_path <- file.path(model_dir, "fit_summary.csv")
  
  if (!all(file.exists(hmc_path, summary_path))) {
    stop("Missing required files in ", model_dir)
  }
  
  # Load data
  hmc_samples <- fread(hmc_path)
  fit_summary <- fread(summary_path)
  
  # Get parameters to keep
  exec_name <- fit_summary$exec[1]
  stan_parameters <- get_stan_parameters(exec_name)
  
  # Extract and save
  condensed_samples <- hmc_samples[name %in% stan_parameters]
  output_path <- file.path(model_dir, "hmc_samples_compact.rda")
  save(condensed_samples, file = output_path)
  
  return(condensed_samples)
}
