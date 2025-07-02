# Catch Fit Plot Functions for SSP Model Analysis
# Adapted from index fit functions for catch/removals data

# =============================================================================
# CATCH FIT FUNCTIONS
# =============================================================================

#' Generate Catch Fit Plot (can handle multiple models)
#' @param model_dirs Vector of model directories
#' @param params List of parameters (uses defaults if NULL)
generate_catch_fit <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract catch fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    fit_dt_list[[i]] <- ssp_extract_catch_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      calc_std = "FALSE"
    )
  }

  yo_dt = tmp_summary[,.(run_label)] %>%
          unique(.) %>%
          .[,year_one:=extract_model_start_year(run_label)]
  
  # Combine and process data
  plot_dt_a <- rbindlist(fit_dt_list) %>%
    merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
    .[, group_id := paste0(run_id, "-catch")] %>%
    .[metric %in% c("obs_catch", "pred_catch","sigmac")] %>%
    na.omit(.) %>%
    merge(.,yo_dt,by="run_label") %>% 
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)]
  
  
  plot_dt <- rbind(plot_dt_a)
  
  # Process observation error data
  obs_se_dt <- plot_dt[metric == "sigmac"] %>%
    .[, .(run_id, run_label, group_id, row, year, value)] %>%
    .[, se := round(value, digits = 3)]
  
  obs_catch_dt <- plot_dt[metric == "obs_catch"] %>%
    .[, .(run_id, run_label, group_id, row, year, value)] %>%
    merge(., obs_se_dt[, .(run_id, run_label, group_id, row, year, se)], by = c("run_id", "run_label", "group_id", "row", "year")) %>%
    .[, obs := round(value, digits = 3)] %>%
    .[, id2 := as.numeric(as.factor(obs))] %>%
    .[, id3 := as.numeric(as.factor(paste0(obs, "_", se)))]
  
  pred_catch_dt <- plot_dt[metric == "pred_catch"] %>%
    .[, .(run_id, run_label, group_id, row, year, iter, value)]
  
  # Process prediction data based on type
  if (params$type == "Median") {
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, .(median = median(value)), by = .(run_label, row, year)]
  } else if (params$type == "Spaghetti") {
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, group_id := paste0(run_label, "-catch-", iter)]
  } else {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, .(median = median(value), upper = quantile(value, probs = 1 - obs_quant), lower = quantile(value, probs = obs_quant)), by = .(run_label, row, year)]
  }
  
  # Check for consistent observations across models
  if (length(model_dirs) > 1 && mean(table(obs_catch_dt$id2) %% uniqueN(obs_catch_dt$run_id) == 0) != 1 && !params$obs) {
    stop("Observed catch not identical between models. Choose models with the same catch or only a single model.")
  } else if (length(model_dirs) > 1 && mean(table(obs_catch_dt$id2) %% uniqueN(obs_catch_dt$run_id) == 0) != 1 && params$obs) {
    if (mean(table(obs_catch_dt$id3) %% uniqueN(obs_catch_dt$run_id) == 0) != 1) {
      stop("Observation error is selected to be plotted however it is not identical. Turn off observation error plotting, choose models with the same observation error or choose a single model.")
    }
  }
  
  # Process observation data
  if (params$obs) {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    obs_catch_dt <- obs_catch_dt[, .(row, year, obs, se)] %>%
      unique(.) %>%
      .[, upper := qlnorm(1 - obs_quant, meanlog = log(obs), sdlog = se)] %>%
      .[, lower := qlnorm(obs_quant, meanlog = log(obs), sdlog = se)]
  } else {
    obs_catch_dt <- obs_catch_dt[, .(row, year, obs)] %>%
      unique(.)
  }
  
  if (nrow(plot_dt) == 0 || nrow(obs_catch_dt) == 0 || nrow(pred_catch_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Create plot
  p <- pred_catch_dt %>%
    ggplot() +
    ylab("Catch") +
    xlab("Year")
  
  # Add observation error bars if requested
  if (params$obs) {
    p <- p + geom_segment(data = obs_catch_dt, aes(x = year, xend = year, y = lower, yend = upper), linewidth = 1.05)
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
  p <- p + geom_point(data = obs_catch_dt, aes(x = year, y = obs), color = "white", fill = "black", shape = 21, size = 3)
  
  # Add styling
  p <- p + 
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch Fit Posterior Predictive Plot
generate_catch_fit_ppd <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract catch fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    fit_dt_list[[i]] <- ssp_extract_catch_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      calc_std = "FALSE"
    )
  }

  yo_dt = tmp_summary[,.(run_label)] %>%
          unique(.) %>%
          .[,year_one:=extract_model_start_year(run_label)]

  # Process data (similar to catch_fit but using ppd_catch instead of pred_catch)
  plot_dt <- rbindlist(fit_dt_list) %>%
    merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
    .[, group_id := paste0(run_id, "-catch")] %>%
    .[metric %in% c("obs_catch", "sigmac", "ppd_catch")] %>%
    na.omit(.) %>%
    merge(.,yo_dt,by="run_label") %>%
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)]
  
  # Process observation error and observed data (same as catch_fit)
  obs_se_dt <- plot_dt[metric == "sigmac"] %>%
    .[, .(run_id, run_label, group_id, row, year, value)] %>%
    .[, se := round(value, digits = 3)]
  
  obs_catch_dt <- plot_dt[metric == "obs_catch"] %>%
    .[, .(run_id, run_label, group_id, row, year, value)] %>%
    merge(., obs_se_dt[, .(run_id, run_label, group_id, row, year, se)], by = c("run_id", "run_label", "group_id", "row", "year")) %>%
    .[, obs := round(value, digits = 3)] %>%
    .[, id2 := as.numeric(as.factor(obs))] %>%
    .[, id3 := as.numeric(as.factor(paste0(obs, "_", se)))]
  
  # Use posterior predictive data instead of predicted data
  pred_catch_dt <- plot_dt[metric == "ppd_catch"] %>%
    .[, .(run_id, run_label, group_id, row, year, iter, value)]
  
  # Process prediction data based on type (same logic as catch_fit)
  if (params$type == "Median") {
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, .(median = median(value)), by = .(run_label, row, year)]
  } else if (params$type == "Spaghetti") {
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, group_id := paste0(run_label, "-catch-", iter)]
  } else {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    pred_catch_dt <- pred_catch_dt[, .(run_label, row, year, iter, value)] %>%
      .[, .(median = median(value), upper = quantile(value, probs = 1 - obs_quant), lower = quantile(value, probs = obs_quant)), by = .(run_label, row, year)]
  }
  
  # Check consistency and process observations (same as catch_fit)
  if (length(model_dirs) > 1 && mean(table(obs_catch_dt$id2) %% uniqueN(obs_catch_dt$run_id) == 0) != 1 && !params$obs) {
    stop("Observed catch not identical between models. Choose models with the same catch or only a single model.")
  } else if (length(model_dirs) > 1 && mean(table(obs_catch_dt$id2) %% uniqueN(obs_catch_dt$run_id) == 0) != 1 && params$obs) {
    if (mean(table(obs_catch_dt$id3) %% uniqueN(obs_catch_dt$run_id) == 0) != 1) {
      stop("Observation error is selected to be plotted however it is not identical. Turn off observation error plotting, choose models with the same observation error or choose a single model.")
    }
  }
  
  if (params$obs) {
    obs_quant <- 0.5 * (1 - as.numeric(params$quants) / 100)
    obs_catch_dt <- obs_catch_dt[, .(row, year, obs, se)] %>%
      unique(.) %>%
      .[, upper := qlnorm(1 - obs_quant, meanlog = log(obs), sdlog = se)] %>%
      .[, lower := qlnorm(obs_quant, meanlog = log(obs), sdlog = se)]
  } else {
    obs_catch_dt <- obs_catch_dt[, .(row, year, obs)] %>%
      unique(.)
  }
  
  if (nrow(plot_dt) == 0 || nrow(obs_catch_dt) == 0 || nrow(pred_catch_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Create plot (same structure as catch_fit)
  p <- pred_catch_dt %>%
    ggplot() +
    ylab("Catch") +
    xlab("Year")
  
  if (params$obs) {
    p <- p + geom_segment(data = obs_catch_dt, aes(x = year, xend = year, y = lower, yend = upper), linewidth = 1.05)
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
  
  p <- p + geom_point(data = obs_catch_dt, aes(x = year, y = obs), color = "white", fill = "black", shape = 21, size = 3)
  
  p <- p + 
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate Catch Fit Residuals Plot
generate_catch_fit_residuals <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$fits
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Extract catch fit data for all models
  fit_dt_list <- list()
  for (i in seq_along(all_data)) {
    calc_std_flag <- ifelse(params$resid == "Standardized", "TRUE", "FALSE")
    fit_dt_list[[i]] <- ssp_extract_catch_fit(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      calc_std = calc_std_flag
    )
  }
  
  # Process residuals based on type
  if (params$resid == "Ordinary") {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-catch")] %>%
      .[metric %in% c("residual")] %>%
      .[, .(value = median(value)), by = .(run_id, run_label, group_id, metric, row)]
    ylab_txt <- "Ordinary residual"
  } else if (params$resid == "Standardized") {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-catch")] %>%
      .[metric %in% c("std_residual")] %>%
      .[, .(value = median(value)), by = .(run_id, run_label, group_id, metric, row)]
    ylab_txt <- "Standardized residual"
  } else {
    plot_dt <- rbindlist(fit_dt_list) %>%
      merge(., tmp_summary[, .(run_id, run_label)], by = "run_id") %>%
      .[, group_id := paste0(run_label, "-catch")] %>%
      .[metric %in% c("pit_residual")] %>%
      .[, .(run_id, run_label, group_id, metric, row, value)]
    ylab_txt <- "PIT residual"
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Add jitter for multiple models and prepare data
  step <- 1 / (length(model_dirs) + 1)
  jitter_seq <- seq(from = -0.5, to = 0.5, by = step)[-1]

  yo_dt = tmp_summary[,.(run_label)] %>%
          unique(.) %>%
          .[,year_one:=extract_model_start_year(run_label)]
  
  plot_dt <- plot_dt %>%
    na.omit(.) %>%
    merge(.,yo_dt,by="run_label") %>%     
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)] %>%
    .[, run_label := factor(run_label, levels = sort(unique(run_label)))] %>%
    .[, year := year + jitter_seq[as.numeric(run_label)]]
  
  # Calculate runs test for single model
  if (length(model_dirs) == 1) {
    if (!requireNamespace("randtests", quietly = TRUE)) {
      warning("randtests package required for runs test")
      runs_dt <- NULL
    } else {
      if (params$resid == "PIT") {
        runs_pval <- randtests::runs.test(plot_dt$value - 0.5, threshold = 0, alternative = "left.sided")[["p.value"]]
        runs_result <- ifelse(runs_pval < 0.05, "Runs test: Fail", "Runs test: Pass")
        runs_dt <- data.table(runs_test = runs_result) %>%
          .[, x := min(plot_dt$year) + 0.05 * diff(range(plot_dt$year))] %>%
          .[, y := 0.9]
      } else {
        runs_pval <- randtests::runs.test(plot_dt$value, threshold = 0, alternative = "left.sided")[["p.value"]]
        runs_result <- ifelse(runs_pval > 0.05, "Runs test: Fail", "Runs test: Pass")
        runs_dt <- data.table(runs_test = runs_result) %>%
          .[, x := min(plot_dt$year) + 0.05 * diff(range(plot_dt$year))] %>%
          .[, y := 0.9 * max(plot_dt$value)]
      }
    }
  } else {
    runs_dt <- NULL
  }
  
  # Create residual plot
  if (params$resid != "PIT") {
    p <- plot_dt %>%
      ggplot() +
      ylab(ylab_txt) +
      xlab("Year") +
      geom_hline(yintercept = 0) +
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

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Generate catch fit plots
# model_dirs <- c("./path/to/model1/", "./path/to/model2/")
# 
# # Single model plots
# p1 <- generate_catch_fit(model_dirs[1])
# p2 <- generate_catch_fit_ppd(model_dirs[1])
# p3 <- generate_catch_fit_residuals(model_dirs[1])
# 
# # Multi-model comparison
# p4 <- generate_catch_fit(model_dirs)
# 
# # Custom parameters
# custom_fit_params <- list(
#   prop = 0.5,
#   obs = FALSE,
#   type = "Quantile",
#   quants = 90,
#   resid = "Standardized"
# )
# p5 <- generate_catch_fit_residuals(model_dirs, custom_fit_params)
