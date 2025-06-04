# Standalone Plot Functions for SSP Model Analysis
# Part 2: Posterior Predictive Check (PPC) Functions

# Note: This requires Part 1 to be loaded first for utility functions

# =============================================================================
# PPC FUNCTIONS
# =============================================================================

#' Generate PPC Density Overlay Plot
#' @param model_dir Directory containing model output files
#' @param params List of parameters (uses defaults if NULL)
generate_ppc_dens <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract CPUE fit data using ssp function
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    active = ifelse(params$active, "TRUE", "FALSE"),
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  if (!params$group) {
    ygroup_vec <- yrep_dt$index
  }
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create density overlay plot
  if (!params$group) {
    p <- bayesplot::ppc_dens_overlay_grouped(y = y_vec, yrep = yrep_mat, group = index_names[ygroup_vec])
  } else {
    p <- bayesplot::ppc_dens_overlay(y = y_vec, yrep = yrep_mat)
  }
  
  p <- p + get_ssp_theme()
  
  return(p)
}

#' Generate PPC ECDF Overlay Plot
generate_ppc_ecdf <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    active = ifelse(params$active, "TRUE", "FALSE"),
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  if (!params$group) {
    ygroup_vec <- yrep_dt$index
  }
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create ECDF overlay plot
  if (!params$group) {
    p <- bayesplot::ppc_ecdf_overlay_grouped(y = y_vec, yrep = yrep_mat, group = index_names[ygroup_vec])
  } else {
    p <- bayesplot::ppc_ecdf_overlay(y = y_vec, yrep = yrep_mat)
  }
  
  p <- p + get_ssp_theme()
  
  return(p)
}

#' Generate PPC PIT ECDF Plot
generate_ppc_pit_ecdf <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    active = ifelse(params$active, "TRUE", "FALSE"),
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  if (!params$group) {
    ygroup_vec <- yrep_dt$index
  }
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create PIT ECDF plot
  if (!params$group) {
    p <- bayesplot::ppc_pit_ecdf_grouped(y = y_vec, yrep = yrep_mat, group = index_names[ygroup_vec])
  } else {
    p <- bayesplot::ppc_pit_ecdf(y = y_vec, yrep = yrep_mat)
  }
  
  p <- p + get_ssp_theme()
  
  return(p)
}

#' Generate PPC Test Statistics Plot
generate_ppc_stat <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppc
  
  # Validate stat parameter
  if (!(length(params$stat) %in% c(1, 2))) {
    stop("Must select at least 1 or at most 2 PPC statistics.")
  }
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = params$prop,
    active = ifelse(params$active, "TRUE", "FALSE"),
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  if (!params$group) {
    ygroup_vec <- yrep_dt$index
  }
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create test statistics plot
  if (length(params$stat) == 1) {
    if (!params$group) {
      p <- bayesplot::ppc_stat_grouped(y = y_vec, yrep = yrep_mat, group = index_names[ygroup_vec], stat = params$stat)
    } else {
      p <- bayesplot::ppc_stat(y = y_vec, yrep = yrep_mat, stat = params$stat)
    }
  } else {
    p <- bayesplot::ppc_stat_2d(y = y_vec, yrep = yrep_mat, stat = params$stat)
  }
  
  p <- p + get_ssp_theme()
  
  return(p)
}

#' Generate PPC LOO-PIT Plot
generate_ppc_loo_pit <- function(model_dir, params = NULL) {
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
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    active = "TRUE",
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  # Calculate likelihood
  lik1_dt <- ssp_calc_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "I", "value"), c("row", "index", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_cpue"]$iter)] %>%
    .[, .(iter, row, index, ll)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_cpue_dt, lik1_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
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

#' Generate PPC LOO-PIT QQ Plot
generate_ppc_loo_qq <- function(model_dir, params = NULL) {
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
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    active = "TRUE",
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  # Calculate likelihood
  lik1_dt <- ssp_calc_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "I", "value"), c("row", "index", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_cpue"]$iter)] %>%
    .[, .(iter, row, index, ll)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_cpue_dt, lik1_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
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

#' Generate PPC LOO Interval Plot
generate_ppc_loo_interval <- function(model_dir, params = NULL) {
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
  
  # Extract CPUE fit data
  fit_dt <- ssp_extract_cpue_fit(
    ssp_summary = data$summary,
    samples_dt = data$samples,
    stan_data = data$stan_data,
    settings = data$settings,
    sub_sample_prop = 1,  # Use full sample for LOO
    active = "TRUE",
    calc_std = "FALSE"
  )
  
  # Prepare data for PPC
  obs_cpue_dt <- fit_dt[metric == "obs_cpue", .(row, index, value)] %>%
    setnames(., "value", "y")
  
  ppd_dt <- fit_dt[metric == "ppd_cpue", .(iter, row, index, value)] %>%
    dcast(., row + index ~ iter)
  
  # Calculate likelihood
  lik1_dt <- ssp_calc_likelihood(data$samples, data$stan_data) %>%
    setnames(., c("T", "I", "value"), c("row", "index", "ll")) %>%
    .[iter %in% unique(fit_dt[metric == "ppd_cpue"]$iter)] %>%
    .[, .(iter, row, index, ll)] %>%
    dcast(., row + index ~ iter)
  
  yrep_dt <- merge(obs_cpue_dt, ppd_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  lik_dt <- merge(obs_cpue_dt, lik1_dt, by = c("row", "index"), all.x = TRUE) %>%
    na.omit(.)
  
  y_vec <- yrep_dt$y
  
  yrep_mat <- yrep_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
    .[, y := NULL] %>%
    as.matrix(.) %>%
    t(.)
  
  log_lik_mat <- lik_dt %>%
    .[, row := NULL] %>%
    .[, index := NULL] %>%
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
# EXAMPLE USAGE
# =============================================================================

# # Generate PPC plots
# model_dir <- "./path/to/model1/"
# 
# # Use defaults
# p1 <- generate_ppc_dens(model_dir)
# p2 <- generate_ppc_ecdf(model_dir)
# p3 <- generate_ppc_pit_ecdf(model_dir)
# 
# # Custom parameters
# custom_ppc_params <- list(
#   scheme = "viridis",
#   prop = 0.5,
#   active = FALSE,
#   group = FALSE,
#   stat = c("mean", "sd"),
#   qqdist = "normal"
# )
# p4 <- generate_ppc_stat(model_dir, custom_ppc_params)
# p5 <- generate_ppc_loo_pit(model_dir, custom_ppc_params)
