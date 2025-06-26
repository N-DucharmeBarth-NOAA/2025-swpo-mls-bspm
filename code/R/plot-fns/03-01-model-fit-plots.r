# Standalone Plot Functions for SSP Model Analysis
# Part 3: Model Fit Functions

# Note: This requires Parts 1-2 to be loaded first

# =============================================================================
# MODEL FIT FUNCTIONS
# =============================================================================

#' Generate Index Fit Plot (can handle multiple models)
#' @param model_dirs Vector of model directories
#' @param params List of parameters (uses defaults if NULL)
generate_index_fit <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract CPUE fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    fit_dt_list[[i]] <- ssp_extract_cpue_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      active = ifelse(params$active, "TRUE", "FALSE"),
      calc_std = "FALSE"
    )
  }
  
  # Combine and process data
  plot_dt_a <- rbindlist(fit_dt_list) %>%
    merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
    .[, group_id := paste0(run_id, "-", index)] %>%
    .[metric %in% c("obs_cpue", "pred_cpue")] %>%
    na.omit(.) %>%
    .[, index := factor(index, levels = sort(unique(index)), labels = index_names[as.numeric(as.character(sort(unique(index))))])] %>%
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)]
  
  plot_dt_b <- rbindlist(fit_dt_list) %>%
    merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
    .[, group_id := paste0(run_id, "-", index)] %>%
    .[metric %in% c("sigmao")] %>%
    na.omit(.) %>%
    .[, index := factor(index, levels = sort(unique(index)), labels = index_names[as.numeric(as.character(sort(unique(index))))])] %>%
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)] %>%
    .[, .(value = median(value)), by = .(run_id, metric, row, index, run_label, group_id, year)] %>%
    .[, iter := 0] %>%
    .[, .(run_id, metric, iter, row, index, value, run_label, group_id, year)]
  
  plot_dt <- rbind(plot_dt_a, plot_dt_b)
  
  # Process observation error data
  obs_se_dt <- plot_dt[metric == "sigmao"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, value)] %>%
    .[, se := round(value, digits = 3)]
  
  obs_cpue_dt <- plot_dt[metric == "obs_cpue"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, value)] %>%
    merge(., obs_se_dt[, .(run_id, run_label, group_id, index, row, year, se)], by = c("run_id", "run_label", "group_id", "index", "row", "year")) %>%
    .[, obs := round(value, digits = 3)] %>%
    .[, id2 := as.numeric(as.factor(paste0(index, "_", obs)))] %>%
    .[, id3 := as.numeric(as.factor(paste0(index, "_", obs, "_", se)))]
  
  pred_cpue_dt <- plot_dt[metric == "pred_cpue"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, iter, value)]
  
  # Process prediction data based on type
  if (params$type == "Median") {
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, .(median = median(value)), by = .(run_label, index, row, year)]
  } else if (params$type == "Spaghetti") {
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, group_id := paste0(run_label, "-", index, "-", iter)]
  } else {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, .(median = median(value), upper = quantile(value, probs = 1 - obs_quant), lower = quantile(value, probs = obs_quant)), by = .(run_label, index, row, year)]
  }
  
  # Check for consistent observations across models (only for shared index-year combinations)
  if (length(model_dirs) > 1 && !params$obs) {
    shared_combinations <- obs_cpue_dt[, .N, by = .(index, year)][N > 1]
    if (nrow(shared_combinations) > 0) {
      inconsistent <- obs_cpue_dt[shared_combinations, on = .(index, year)][
        , .(unique_obs = uniqueN(obs)), by = .(index, year)
      ][unique_obs > 1]
      if (nrow(inconsistent) > 0) {
        stop("Observed CPUE not identical between models for shared index-year combinations.")
      }
    }
  }
  
  # Process observation data
  if (params$obs) {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    obs_cpue_dt <- obs_cpue_dt[, .(index, row, year, obs, se)] %>%
      .[, upper := qlnorm(1 - obs_quant, meanlog = log(obs), sdlog = se)] %>%
      .[, lower := qlnorm(obs_quant, meanlog = log(obs), sdlog = se)]
  } else {
    obs_cpue_dt <- obs_cpue_dt[, .(index, row, year, obs)]
  }
  
  if (nrow(plot_dt) == 0 || nrow(obs_cpue_dt) == 0 || nrow(pred_cpue_dt) == 0) {
    stop("No data available for plotting")
  }

  # Define plot dims
  facet_variable = uniqueN(pred_cpue_dt$index)
  plot_ncol = if(is.null(params$ncol)) min(c(3, facet_variable)) else params$ncol
  plot_nrow = ceiling(facet_variable/plot_ncol)
  
  # Create plot
  p <- pred_cpue_dt %>%
    ggplot() +
    ylab("Index") +
    xlab("Year") +
    geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~index, 
           ncol = plot_ncol,
           nrow = plot_nrow)
  
  # Add observation error bars if requested
  if (params$obs) {
    p <- p + geom_segment(data = obs_cpue_dt, aes(x = year, xend = year, y = lower, yend = upper), linewidth = 1.05)
  }
  
  # Add prediction lines/ribbons based on type
  if (params$type == "Median") {
    p <- p + geom_line(aes(x = year, y = median, color = run_label), linewidth = 1.1)
  } else if (params$type == "Spaghetti") {
    p <- p + geom_line(aes(x = year, y = value, color = run_label, group = group_id), alpha = 0.1)
  } else {
    p <- p + 
      geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = run_label), alpha = 0.4) +
      geom_line(aes(x = year, y = median, color = run_label), linewidth = 1.1)
  }
  
  # Add observed data points
  p <- p + geom_point(data = obs_cpue_dt, aes(x = year, y = obs), color = "white", fill = "black", shape = 21, size = 3)
  
  # Add styling
  p <- p + 
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Index Fit Posterior Predictive Plot
generate_index_fit_ppd <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract CPUE fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    fit_dt_list[[i]] <- ssp_extract_cpue_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      active = ifelse(params$active, "TRUE", "FALSE"),
      calc_std = "FALSE"
    )
  }
  
  # Process data (similar to index_fit but using ppd_cpue instead of pred_cpue)
  plot_dt <- rbindlist(fit_dt_list) %>%
    merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
    .[, group_id := paste0(run_id, "-", index)] %>%
    .[metric %in% c("obs_cpue", "sigmao", "ppd_cpue")] %>%
    na.omit(.) %>%
    .[, index := factor(index, levels = sort(unique(index)), labels = index_names[as.numeric(as.character(sort(unique(index))))])] %>%
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)]
  
  # Process observation error and observed data (same as index_fit)
  obs_se_dt <- plot_dt[metric == "sigmao"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, value)] %>%
    .[, se := round(value, digits = 3)]
  
  obs_cpue_dt <- plot_dt[metric == "obs_cpue"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, value)] %>%
    merge(., obs_se_dt[, .(run_id, run_label, group_id, index, row, year, se)], by = c("run_id", "run_label", "group_id", "index", "row", "year")) %>%
    .[, obs := round(value, digits = 3)] %>%
    .[, id2 := as.numeric(as.factor(paste0(index, "_", obs)))] %>%
    .[, id3 := as.numeric(as.factor(paste0(index, "_", obs, "_", se)))]
  
  # Use posterior predictive data instead of predicted data
  pred_cpue_dt <- plot_dt[metric == "ppd_cpue"] %>%
    .[, .(run_id, run_label, group_id, index, row, year, iter, value)]
  
  # Process prediction data based on type (same logic as index_fit)
  if (params$type == "Median") {
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, .(median = median(value)), by = .(run_label, index, row, year)]
  } else if (params$type == "Spaghetti") {
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, group_id := paste0(run_label, "-", index, "-", iter)]
  } else {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    pred_cpue_dt <- pred_cpue_dt[, .(run_label, index, row, year, iter, value)] %>%
      .[, .(median = median(value), upper = quantile(value, probs = 1 - obs_quant), lower = quantile(value, probs = obs_quant)), by = .(run_label, index, row, year)]
  }
  
  # Check for consistent observations across models (only for shared index-year combinations)
  if (length(model_dirs) > 1 && !params$obs) {
    shared_combinations <- obs_cpue_dt[, .N, by = .(index, year)][N > 1]
    if (nrow(shared_combinations) > 0) {
      inconsistent <- obs_cpue_dt[shared_combinations, on = .(index, year)][
        , .(unique_obs = uniqueN(obs)), by = .(index, year)
      ][unique_obs > 1]
      if (nrow(inconsistent) > 0) {
        stop("Observed CPUE not identical between models for shared index-year combinations.")
      }
    }
  }
  
  if (params$obs) {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    obs_cpue_dt <- obs_cpue_dt[, .(index, row, year, obs, se)] %>%
      .[, upper := qlnorm(1 - obs_quant, meanlog = log(obs), sdlog = se)] %>%
      .[, lower := qlnorm(obs_quant, meanlog = log(obs), sdlog = se)]
  } else {
    obs_cpue_dt <- obs_cpue_dt[, .(index, row, year, obs)]
  }
  
  if (nrow(plot_dt) == 0 || nrow(obs_cpue_dt) == 0 || nrow(pred_cpue_dt) == 0) {
    stop("No data available for plotting")
  }

  # Define plot dims
  facet_variable = uniqueN(pred_cpue_dt$index)
  plot_ncol = if(is.null(params$ncol)) min(c(3, facet_variable)) else params$ncol
  plot_nrow = ceiling(facet_variable/plot_ncol)
  
  # Create plot (same structure as index_fit)
  p <- pred_cpue_dt %>%
    ggplot() +
    ylab("Index") +
    xlab("Year") +
    geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~index, 
           ncol = plot_ncol,
           nrow = plot_nrow)
  
  if (params$obs) {
    p <- p + geom_segment(data = obs_cpue_dt, aes(x = year, xend = year, y = lower, yend = upper), linewidth = 1.05)
  }
  
  if (params$type == "Median") {
    p <- p + geom_line(aes(x = year, y = median, color = run_label), linewidth = 1.1)
  } else if (params$type == "Spaghetti") {
    p <- p + geom_line(aes(x = year, y = value, color = run_label, group = group_id), alpha = 0.1)
  } else {
    p <- p + 
      geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = run_label), alpha = 0.4) +
      geom_line(aes(x = year, y = median, color = run_label), linewidth = 1.1)
  }
  
  p <- p + geom_point(data = obs_cpue_dt, aes(x = year, y = obs), color = "white", fill = "black", shape = 21, size = 3)
  
  p <- p + 
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Index Fit Residuals Plot
generate_index_fit_residuals <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract CPUE fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    calc_std_flag <- ifelse(params$resid == "Standardized", "TRUE", "FALSE")
    fit_dt_list[[i]] <- ssp_extract_cpue_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      active = ifelse(params$active, "TRUE", "FALSE"),
      calc_std = calc_std_flag
    )
  }
  
  # Process residuals based on type
  if (params$resid == "Ordinary") {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-", index)] %>%
      .[metric %in% c("residual")] %>%
      .[, .(value = median(value)), by = .(run_id, run_label, group_id, metric, row, index)]
    ylab_txt <- "Ordinary residual"
  } else if (params$resid == "Standardized") {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-", index)] %>%
      .[metric %in% c("std_residual")] %>%
      .[, .(value = median(value)), by = .(run_id, run_label, group_id, metric, row, index)]
    ylab_txt <- "Standardized residual"
  } else {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-", index)] %>%
      .[metric %in% c("pit_residual")] %>%
      .[, .(run_id, run_label, group_id, metric, row, index, value)]
    ylab_txt <- "PIT residual"
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Add jitter for multiple models and prepare data
  step <- 1 / (length(model_dirs) + 1)
  jitter_seq <- seq(from = -0.5, to = 0.5, by = step)[-1]
  
  plot_dt <- plot_dt %>%
    na.omit(.) %>%
    .[, index := factor(index, levels = sort(unique(index)), labels = index_names[as.numeric(as.character(sort(unique(index))))])] %>%
    .[, run_label := factor(run_label, levels = sort(unique(run_label)))] %>%
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)] %>%
    .[, year := year + jitter_seq[as.numeric(run_label)]]
  
  # Calculate runs test for single model
  if (length(model_dirs) == 1) {
    if (!requireNamespace("randtests", quietly = TRUE)) {
      warning("randtests package required for runs test")
      runs_dt <- NULL
    } else {
      u_index <- unique(plot_dt$index)
      runs_vec <- rep(NA, length(u_index))
      
      if (params$resid == "PIT") {
        for (i in 1:length(u_index)) {
          runs_vec[i] <- randtests::runs.test(plot_dt[index == u_index[i]]$value - 0.5, threshold = 0, alternative = "left.sided")[["p.value"]]
        }
        runs_vec <- ifelse(runs_vec < 0.05, "Runs test: Fail", "Runs test: Pass")
        runs_dt <- data.table(index = u_index, runs_test = runs_vec) %>%
          .[, x := min(plot_dt$year) + 0.05 * diff(range(plot_dt$year))] %>%
          .[, y := 0.9]
      } else {
        for (i in 1:length(u_index)) {
          runs_vec[i] <- randtests::runs.test(plot_dt[index == u_index[i]]$value, threshold = 0, alternative = "left.sided")[["p.value"]]
        }
        runs_vec <- ifelse(runs_vec > 0.05, "Runs test: Fail", "Runs test: Pass")
        runs_dt <- data.table(index = u_index, runs_test = runs_vec) %>%
          .[, x := min(plot_dt$year) + 0.05 * diff(range(plot_dt$year))] %>%
          .[, y := 0.9 * max(plot_dt$value)]
      }
    }
  } else {
    runs_dt <- NULL
  }

  # Define plot dims
  facet_variable = uniqueN(plot_dt$index)
  plot_ncol = if(is.null(params$resid_ncol)) min(c(3, facet_variable)) else params$resid_ncol
  plot_nrow = ceiling(facet_variable/plot_ncol)
  
  # Create residual plot
  if (params$resid != "PIT") {
    p <- plot_dt %>%
      ggplot() +
      ylab(ylab_txt) +
      xlab("Year") +
      geom_hline(yintercept = 0) +
      facet_wrap(~index, 
           ncol = plot_ncol,
           nrow = plot_nrow) +
      geom_segment(aes(x = year, xend = year, y = 0, yend = value, color = run_label, group = group_id), linewidth = 1.05) +
      geom_point(aes(x = year, y = value, fill = run_label, group = group_id), color = "black", shape = 21, size = 3, alpha = 0.5) +
      viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
      viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
      get_ssp_theme()
  } else {
    p <- plot_dt %>%
      ggplot() +
      ylab(ylab_txt) +
      xlab("Year") +
      ylim(0, 1) +
      geom_hline(yintercept = 0.5) +
      facet_wrap(~index, 
           ncol = plot_ncol,
           nrow = plot_nrow) +
      geom_segment(aes(x = year, xend = year, y = 0.5, yend = value, color = run_label, group = group_id), linewidth = 1.05) +
      geom_point(aes(x = year, y = value, fill = run_label, group = group_id), color = "black", shape = 21, size = 3, alpha = 0.5) +
      viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
      viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
      get_ssp_theme()
  }
  
  # Add runs test results for single model
  if (length(model_dirs) == 1 && !is.null(runs_dt)) {
    p <- p + geom_label(data = runs_dt, aes(x = x, y = y, label = runs_test), hjust = "inward")
  }
  
  return(p)
}

#' Generate Prior-Posterior Parameter Comparison Plot
generate_ppp <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$ppp
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  parameter_map <- get_parameter_map()
  target_par <- params$leading_params
  
  # Get posterior data
  posterior_dt <- rbindlist(lapply(all_data, function(x) x$samples)) %>%
    setnames(., "value", "Posterior")
  
  # Get prior data
  prior_dt_list <- list()
  for (i in seq_along(all_data)) {
    prior_dt_list[[i]] <- ssp_prior_pushforward(
      ssp_summary = all_data[[i]]$summary,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings
    )
  }
  prior_dt <- rbindlist(prior_dt_list) %>%
    setnames(., "value", "Prior")
  
  # Merge and process data
  plot_dt <- merge(prior_dt, posterior_dt[, .(run_id, iter, chain, chain_iter, variable, name, row, col, Posterior)], 
                   by = c("run_id", "iter", "chain", "chain_iter", "variable", "name", "row", "col"))
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[, "input"]), "transformed"]
    plot_dt <- plot_dt %>%
      .[, .(run_id, iter, name, Prior, Posterior)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[, "transformed"])]
  } else {
    grab_parameters <- parameter_map[match(target_par, parameter_map[, "input"]), "raw"]
    plot_dt <- plot_dt %>%
      .[, .(run_id, iter, name, Prior, Posterior)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[, "raw"])]
  }
  
  # Reshape and process based on combine setting
  if (params$combine) {
    plot_dt <- merge(plot_dt, tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, .(run_label, iter, name, Prior, Posterior)] %>%
      melt(., id.vars = c("run_label", "iter", "name")) %>%
      .[, variable := factor(variable, levels = c("Prior", "Posterior"))] %>%
      .[variable == "Posterior", run_label := "Combined posterior"] %>%
      .[variable == "Prior", run_label := "Combined prior"] %>%
      .[, run_label := factor(run_label, levels = c("Combined prior", "Combined posterior"))]
  } else {
    plot_dt <- merge(plot_dt, tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, .(run_label, iter, name, Prior, Posterior)] %>%
      melt(., id.vars = c("run_label", "iter", "name")) %>%
      .[, variable := factor(variable, levels = c("Prior", "Posterior"))]
  }
  
  # Filter based on show parameter
  if (params$show == "Prior") {
    plot_dt <- plot_dt[variable == "Prior"]
  } else if (params$show == "Posterior") {
    plot_dt <- plot_dt[variable == "Posterior"]
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }

  # Define plot dims
  facet_variable = uniqueN(plot_dt$name)
  plot_ncol = if(is.null(params$ncol)) min(c(3, facet_variable)) else params$ncol
  plot_nrow = ceiling(facet_variable/plot_ncol)
  
  # Create plot
  p <- plot_dt %>%
    ggplot() +
    ylab("Density") +
    xlab("Parameter") +
    facet_wrap(~name, scales = "free_x", 
           ncol = plot_ncol,
           nrow = plot_nrow)
  
  if (params$show == "Prior" || params$show == "Both") {
    p <- p + geom_density(data = plot_dt[variable == "Prior"], 
                         aes(x = value, y = after_stat(scaled), color = run_label, linetype = variable), 
                         linewidth = 1.15)
  }
  if (params$show == "Posterior" || params$show == "Both") {
    p <- p + geom_density(data = plot_dt[variable == "Posterior"], 
                         aes(x = value, y = after_stat(scaled), color = run_label, fill = run_label, linetype = variable), 
                         alpha = 0.4, linewidth = 1.15)
  }
  
  p <- p + 
    scale_linetype_manual("Distribution type", values = c("dotted", "solid"), drop = FALSE) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Generate model fit plots
# model_dirs <- c("./path/to/model1/", "./path/to/model2/")
# 
# # Single model plots
# p1 <- generate_index_fit(model_dirs[1])
# p2 <- generate_index_fit_ppd(model_dirs[1])
# p3 <- generate_index_fit_residuals(model_dirs[1])
# 
# # Multi-model comparison
# p4 <- generate_index_fit(model_dirs)
# p5 <- generate_ppp(model_dirs)
# 
# # Custom parameters
# custom_fit_params <- list(
#   prop = 0.5,
#   active = FALSE,
#   obs = FALSE,
#   type = "Quantile",
#   quants = 90,
#   resid = "Standardized"
# )
# p6 <- generate_index_fit_residuals(model_dirs, custom_fit_params)
