# Catch PPC Plot Functions for SSP Model Analysis
# Adapted from CPUE PPC functions for catch/removals data

# =============================================================================
# CATCH PPC FUNCTIONS
# =============================================================================

#' Generate Catch PPC Density Overlay Plot
#' @param model_dir Directory containing model output files
#' @param params List of parameters (uses defaults if NULL)
generate_catch_ppc_dens <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data using ssp function
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC (catch has no index, only row)
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create density overlay plot (no grouping for catch - single time series)
  p <- bayesplot::ppc_dens_overlay(y = y_vec, yrep = yrep_mat) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC ECDF Overlay Plot
generate_catch_ppc_ecdf <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create ECDF overlay plot
  p <- bayesplot::ppc_ecdf_overlay(y = y_vec, yrep = yrep_mat) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC PIT ECDF Plot
generate_catch_ppc_pit_ecdf <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create PIT ECDF plot
  p <- bayesplot::ppc_pit_ecdf(y = y_vec, yrep = yrep_mat) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC Test Statistics Plot
generate_catch_ppc_stat <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Validate stat parameter
  if (!(length(params$stat) %in% c(1, 2))) {
    stop("Must select at least 1 or at most 2 PPC statistics.")
  }
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create test statistics plot
  if (length(params$stat) == 1) {
    p <- bayesplot::ppc_stat(y = y_vec, yrep = yrep_mat, stat = params$stat)
  } else {
    p <- bayesplot::ppc_stat_2d(y = y_vec, yrep = yrep_mat, stat = params$stat)
  }
  
  p <- p + get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC LOO-PIT Plot
generate_catch_ppc_loo_pit <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load required libraries for LOO
  if (!requireNamespace("rstantools", quietly = TRUE)) {
    stop("rstantools package required for LOO-PIT plots")
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("loo package required for LOO-PIT plots")
  }
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  # Calculate catch likelihood - only for rows with predictions
  lik1_dt <- ssp_calc_catch_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "value"), c("row", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_catch"]$iter)] %>%
    .[row %in% unique(fit_dt[metric == "ppd_catch"]$row)] %>%  # Only rows with predictions
    .[, .(iter, row, ll)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_catch_dt[, .(row, y)], lik1_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Calculate LOO
  r_eff <- loo::relative_eff(exp(log_lik_mat), unique(data$samples[, .(iter, chain)])$chain, cores = 1)
  loo_result <- loo::loo(log_lik_mat, r_eff = r_eff, cores = 1, save_psis = TRUE)
  psis <- loo_result$psis_object
  lw <- weights(psis)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create LOO-PIT plot
  p <- bayesplot::ppc_loo_pit_overlay(y = y_vec, yrep = yrep_mat, lw = lw) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC LOO-PIT QQ Plot
generate_catch_ppc_loo_qq <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load required libraries for LOO
  if (!requireNamespace("rstantools", quietly = TRUE)) {
    stop("rstantools package required for LOO-PIT plots")
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("loo package required for LOO-PIT plots")
  }
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  # Calculate catch likelihood - only for rows with predictions
  lik1_dt <- ssp_calc_catch_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "value"), c("row", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_catch"]$iter)] %>%
    .[row %in% unique(fit_dt[metric == "ppd_catch"]$row)] %>%  # Only rows with predictions
    .[, .(iter, row, ll)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_catch_dt[, .(row, y)], lik1_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Calculate LOO
  r_eff <- loo::relative_eff(exp(log_lik_mat), unique(data$samples[, .(iter, chain)])$chain, cores = 1)
  loo_result <- loo::loo(log_lik_mat, r_eff = r_eff, cores = 1, save_psis = TRUE)
  psis <- loo_result$psis_object
  lw <- weights(psis)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create LOO-PIT QQ plot
  p <- bayesplot::ppc_loo_pit_qq(y = y_vec, yrep = yrep_mat, lw = lw, compare = params$qqdist) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch PPC LOO Interval Plot
generate_catch_ppc_loo_interval <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load required libraries for LOO
  if (!requireNamespace("rstantools", quietly = TRUE)) {
    stop("rstantools package required for LOO interval plots")
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("loo package required for LOO interval plots")
  }
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract catch fit data
  fit_dt <- ssp_extract_catch_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  # Account for terminal year - predicted catch has one fewer observation
  obs_catch_dt <- fit_dt[metric == "obs_catch", .(row, value)] %>%
    setnames(., "value", "y") %>%
    .[row <= max(fit_dt[metric == "ppd_catch"]$row)]  # Match ppd_catch row range
  
  ppd_dt <- fit_dt[metric == "ppd_catch", .(iter, row, value)] %>%
    dcast(., row ~ iter)
  
  # Calculate catch likelihood - only for rows with predictions
  lik1_dt <- ssp_calc_catch_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "value"), c("row", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_catch"]$iter)] %>%
    .[row %in% unique(fit_dt[metric == "ppd_catch"]$row)] %>%  # Only rows with predictions
    .[, .(iter, row, ll)] %>%
    dcast(., row ~ iter)
  
  yrep_dt <- merge(obs_catch_dt, ppd_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_catch_dt[, .(row, y)], lik1_dt, by = "row", all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Calculate LOO
  r_eff <- loo::relative_eff(exp(log_lik_mat), unique(data$samples[, .(iter, chain)])$chain, cores = 1)
  loo_result <- loo::loo(log_lik_mat, r_eff = r_eff, cores = 1, save_psis = TRUE)
  psis <- loo_result$psis_object
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create LOO interval plot
  p <- bayesplot::ppc_loo_intervals(y = y_vec, yrep = yrep_mat, psis_object = psis) +
    get_ssp_theme()
  
  return(p)
}

# =============================================================================
# HELPER FUNCTION FOR CATCH LIKELIHOOD (if not already implemented)
# =============================================================================

#' Calculate catch likelihood for LOO-based PPC plots
#' @param samples_dt Samples data table
#' @param stan_data Stan data table
ssp_calc_catch_likelihood <- function(samples_dt, stan_data) {
  # Extract relevant data
  obs_removals <- stan_data[name == "obs_removals"]$value
  sigmac <- stan_data[name == "sigmac"]$value
  
  # Get predicted removals
  removals_dt <- samples_dt[name == "removals", .(iter, row, value)]
  
  # Calculate log-likelihood for each observation
  # Note: removals has one fewer observation than obs_removals (no terminal year prediction)
  ll_dt <- removals_dt[, {
    if (row <= length(obs_removals)) {  # Ensure we don't exceed obs_removals length
      mu_catch <- log(value) - 0.5 * sigmac^2
      ll <- dlnorm(obs_removals[row], meanlog = mu_catch, sdlog = sigmac, log = TRUE)
      .(T = row, ll = ll)
    } else {
      .(T = integer(0), ll = numeric(0))  # Return empty if row exceeds observations
    }
  }, by = .(iter, row)]
  
  return(ll_dt[, .(iter, T, value = ll)])
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Generate catch PPC plots
# model_dir <- "./path/to/model1/"
# 
# # Use defaults
# p1 <- generate_catch_ppc_dens(model_dir)
# p2 <- generate_catch_ppc_ecdf(model_dir)
# p3 <- generate_catch_ppc_pit_ecdf(model_dir)
# 
# # Custom parameters
# custom_ppc_params <- list(
#   scheme = "viridis",
#   prop = 0.5,
#   stat = c("mean", "sd"),
#   qqdist = "normal"
# )
# p4 <- generate_catch_ppc_stat(model_dir, custom_ppc_params)
# p5 <- generate_catch_ppc_loo_pit(model_dir, custom_ppc_params)
