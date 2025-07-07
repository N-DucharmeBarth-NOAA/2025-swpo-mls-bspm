# Standalone Plot Functions for SSP Model Analysis
# Part 4D: Forecasting Function

# Note: This requires Parts 1-3, 4A, 4B, and 4C to be loaded first

# =============================================================================
# FORECASTING FUNCTION
# =============================================================================

#' Generate Forecast Plot
#' @param model_dirs Vector of model directories
#' @param params List of parameters (uses defaults if NULL)
generate_fcast <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$forecasts
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)

  if(is.null(params$model_names)){
    short_plot_names = tmp_summary$run_label
    short_plot_names = sapply(short_plot_names,function(x)strsplit(x,"-")[[1]][1])  
  } else {
    if(length(params$model_names)!=nrow(tmp_summary)){
      stop("`model_names` does not have the correct length")
    } else if(uniqueN(params$model_names)!=length(params$model_names)){
      stop("`model_names` can not have duplicate names")
    } else {
      short_plot_names = params$model_names
    }
  }

  # Parameter mapping
  parameter_map <- cbind(
    c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Surplus production"),
    c("D", "P", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "removals", "dev", "raw_epsp", "surplus_production")
  )
  colnames(parameter_map) <- c("input", "grab")
  target_par <- params$var
  
  # Generate forecasts for each model
  posterior_dt_list <- list()
  for (i in seq_along(all_data)) {
    # Run forecasts
    tmp_forecast <- ssp_forecast(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = all_data[[i]]$samples,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = params$prop,
      forecast_years = params$nyears,
      forecast_type = params$type,
      avg_years = params$avg_year,
      scalar = params$scalar,
      resample_raw_epsp = ifelse(params$resample_epsp, "TRUE", "FALSE")
    )
    
    posterior_dt_list[[i]] <- ssp_derived_quants_ts(
      ssp_summary = all_data[[i]]$summary,
      samples_dt = tmp_forecast,
      stan_data = all_data[[i]]$stan_data,
      settings = all_data[[i]]$settings,
      sub_sample_prop = 1
    )
  }
  
  # Process data
  grab_parameters <- parameter_map[match(target_par, parameter_map[, "input"]), "grab"]
  posterior_dt <- rbindlist(posterior_dt_list) %>%
    .[name %in% grab_parameters] %>%
    .[, name := factor(name, levels = parameter_map[match(target_par, parameter_map[, "input"]), "grab"], 
                      labels = parameter_map[match(target_par, parameter_map[, "input"]), "input"])] %>%
    .[, type := "Posterior"] %>%
    .[, type := factor(type, levels = c("Prior", "Posterior"))]
  
  plot_dt <- posterior_dt %>%
    merge(., tmp_summary[, .(run_id, run_label)])
  
  # Handle combine option
  if (params$combine) {
    plot_dt <- plot_dt %>%
      .[type == "Posterior", run_label := "Combined posterior"]
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }

  yo_dt = tmp_summary[,.(run_label)] %>%
          unique(.) %>%
          .[,year_one:=extract_model_start_year(run_label)]
  
  # Summarize data
  obs_quant <- 0.5 * (1 - (as.numeric(params$quants) - 1e-1) / 100)
  plot_dt <- plot_dt %>%
    .[!is.na(value)] %>%
    .[, .(med = median(value), avg = mean(value), 
          lp = quantile(value, probs = obs_quant), 
          up = quantile(value, probs = 1 - obs_quant)), 
      by = .(run_label, type, name, row)] %>%
     merge(.,yo_dt,by="run_label") %>%   
    .[row >= 1, year := year_one + (row - 1)] %>%
    .[row < 1, year := year_one + (row - 1)] %>%
    .[name %in% c("Process error", "Process error (raw)") & row > 0, year := year + 1]
  
  # Define plot dims
  facet_variable = uniqueN(plot_dt$name)
  plot_ncol = if(is.null(params$ncol)) min(c(3, facet_variable)) else params$ncol
  plot_nrow = ceiling(facet_variable/plot_ncol)

  # Create forecast plot
  plot_dt = plot_dt %>% .[,run_label:=factor(run_label,levels=tmp_summary$run_label,labels=short_plot_names)]
  p <- plot_dt %>%
    ggplot() +
    ylab("Metric") +
    xlab("Year") +
    facet_wrap(~name, scales = "free_y", 
           ncol = plot_ncol,
           nrow = plot_nrow) +
    geom_ribbon(data = plot_dt[type == "Posterior"], 
               aes(x = year, ymin = lp, ymax = up, fill = run_label), 
               alpha = 0.4, linewidth = 1.15) +
    geom_path(data = plot_dt[type == "Posterior"], 
             aes(x = year, y = med, color = run_label), linewidth = 1.15) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = max(plot_dt$year) - as.numeric(params$nyears) + 0.5) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Generate forecast plot
# model_dirs <- c("./path/to/model1/", "./path/to/model2/")
# 
# # Use defaults
# p1 <- generate_fcast(model_dirs)
# 
# # Custom parameters
# custom_forecast_params <- list(
#   var = c("Depletion (D)", "F_Fmsy"),
#   combine = FALSE,
#   prop = 0.25,
#   quants = 95,
#   nyears = 10,
#   resample_epsp = FALSE,
#   type = "MSY",
#   avg_year = 5,
#   scalar = 0.8
# )
# p2 <- generate_fcast(model_dirs, custom_forecast_params)
