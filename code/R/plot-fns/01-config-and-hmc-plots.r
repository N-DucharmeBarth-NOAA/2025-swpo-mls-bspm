# Standalone Plot Functions for SSP Model Analysis
# Part 1: Configuration and HMC Diagnostic Functions

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(viridis)
  library(bayesplot)
  library(GGally)
  library(MASS)
  library(randtests)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default parameters extracted from UI
get_default_params <- function() {
  list(
    hmc = list(
      leading_params = c("logK","x0","r"),
      raw = TRUE,
      diag = "None",
      eps = TRUE,
      lags = 30,
      scheme = "brightblue"
    ),
    ppc = list(
      scheme = "brightblue",
      prop = 0.25,
      active = TRUE,
      group = TRUE,
      stat = "median",
      qqdist = "uniform"
    ),
    fits = list(
      prop = 0.25,
      active = TRUE,
      obs = TRUE,
      type = "Median",
      quants = 95,
      resid = "PIT",
      ncol = NULL,  
      resid_ncol = NULL,
      model_names = NULL
    ),
    ppp = list(
      leading_params = c("logK","x0","r"),
      raw = TRUE,
      show = "Both",
      combine = FALSE,
      ncol = NULL,
      model_names = NULL
    ),
    ppts = list(
      var = c("Depletion (D)","F_Fmsy","Removals","Process error (mult.)"),
      show = "Both",
      combine = FALSE,
      prop = 0.25,
      quants = 95,
      ncol = NULL,
      model_names = NULL
    ),
    kbmj = list(
      show = "Both",
      combine = FALSE,
      prop = 0.25,
      uncertainty = TRUE,
      quants = 95,
      resolution = 100,
      model_names = NULL
    ),
    forecasts = list(
      var = c("Depletion (D)","F_Fmsy","Removals","Process error"),
      combine = FALSE,
      prop = 0.25,
      quants = 95,
      nyears = 5,
      resample_epsp = TRUE,
      type = "Catch",
      avg_year = 3,
      scalar = 1,
      ncol = NULL,
      model_names = NULL
    )
  )
}

# Global configuration (to be set by user)
set_global_config <- function( 
                             index_names = c("Index 1", "Index 2"), 
                             model_stem = "./output/",
                             height_per_panel = 350) {
  assign("index_names", index_names, envir = .GlobalEnv)
  assign("model_stem", model_stem, envir = .GlobalEnv)
  assign("height_per_panel", height_per_panel, envir = .GlobalEnv)
}

extract_model_start_year <- function(model_string) {
  has_sy <- grepl("-sy[0-9]{4}", model_string)
  result <- rep(1952L, length(model_string))
  result[has_sy] <- as.integer(sub(".*-sy([0-9]{4}).*", "\\1", model_string[has_sy]))
  result
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Load model data
load_model_data <- function(model_dir) {
  # Ensure model_dir ends with "/"
  if (!endsWith(model_dir, "/")) {
    model_dir <- paste0(model_dir, "/")
  }
  
  # Core required files
  data_list <- list(
    summary = fread(paste0(model_dir, "fit_summary.csv")),
    samples = fread(paste0(model_dir, "hmc_samples.csv")),
    stan_data = fread(paste0(model_dir, "stan_data.csv")),
    settings = fread(paste0(model_dir, "settings.csv"))
  )
  
  # Optional diagnostic files (may not exist for all models)
  rhat_file <- paste0(model_dir, "rhat.csv")
  if (file.exists(rhat_file)) {
    data_list$rhat <- fread(rhat_file)
  } else {
    warning(paste("rhat.csv not found in", model_dir, "- R-hat plots will not be available"))
    data_list$rhat <- NULL
  }
  
  neff_file <- paste0(model_dir, "neff.csv")
  if (file.exists(neff_file)) {
    data_list$neff <- fread(neff_file)
  } else {
    warning(paste("neff.csv not found in", model_dir, "- N_eff plots will not be available"))
    data_list$neff <- NULL
  }
  
  return(data_list)
}

# Parameter mapping
get_parameter_map <- function() {
  map <- cbind(
    c("logK","x0","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","sigmaf","lp__"),
    c("logK","x0","r","sigmao_add","sigmap","n","qeff","rho","sigma_qdev","sigmaf","lp__"),
    c("raw_logK","raw_logx0","raw_logr","raw_sigmao_add","raw_logsigmap","raw_logshape","raw_logqeff","raw_rho","raw_sigma_qdev","raw_sigmaf","lp__")
  )
  colnames(map) <- c("input","transformed","raw")
  return(map)
}

get_parameter_map_extended <- function() {
  map <- cbind(
    c("logK","x0","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","sigmaf","epsp","qdev_period","edev","lp__"),
    c("logK","x0","r","sigmao_add","sigmap","n","qeff","rho","sigma_qdev","sigmaf","epsp","qdev_period","edev","lp__"),
    c("raw_logK","raw_logx0","raw_logr","raw_sigmao_add","raw_logsigmap","raw_logshape","raw_logqeff","raw_rho","raw_sigma_qdev","raw_sigmaf","raw_epsp","raw_qdev_period","raw_edev","lp__")
  )
  colnames(map) <- c("input","transformed","raw")
  return(map)
}

# Standard theme
get_ssp_theme <- function() {
  theme(
    text = element_text(size = 20),
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
    strip.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
}

# =============================================================================
# HMC DIAGNOSTIC FUNCTIONS
# =============================================================================

#' Generate HMC Parallel Coordinates Plot
#' @param model_dir Directory containing model output files
#' @param params List of parameters (uses defaults if NULL)
generate_hmc_parcoord <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  tmp_summary <- data$summary
  parameter_map <- get_parameter_map()
  
  # Check diagnostics
  if (tmp_summary$treedepth == 0 && tmp_summary$divergent == 0 && params$diag != "None") {
    warning("No diagnostic issues found, but diagnostic coloring requested")
  }
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "transformed"]
    plot_title <- "Parcoord (transformed)"
    plot_dt <- data$samples %>%
      .[, .(iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"transformed"])]
  } else {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "raw"]
    plot_title <- "Parcoord (original)"
    plot_dt <- data$samples %>%
      .[, .(iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"raw"])]
  }
  
  # Handle diagnostics
  if (params$diag == "None") {
    plot_dt <- plot_dt[, .(iter, name, value)] %>%
      .[, type := "None"]
  } else if (params$diag == "Max. treedepth") {
    if (tmp_summary$treedepth > 0) {
      max_tree <- max(plot_dt$treedepth)
      plot_dt <- plot_dt[, .(iter, name, value, treedepth)] %>%
        .[, type := "None"] %>%
        .[treedepth == max_tree, type := "Max treedepth"] %>%
        .[, type := factor(type, levels = c("None", "Max treedepth"))]
    } else {
      plot_dt <- plot_dt[, .(iter, name, value, treedepth)] %>%
        .[, type := "None"]
    }
  } else {
    if (tmp_summary$divergent > 0) {
      plot_dt <- plot_dt[, .(iter, name, value, divergent)] %>%
        .[, type := "None"] %>%
        .[divergent == 1, type := "Divergent transition"] %>%
        .[, type := factor(type, levels = c("None", "Divergent transition"))]
    } else {
      plot_dt <- plot_dt[, .(iter, name, value, divergent)] %>%
        .[, type := "None"]
    }
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Create plot
  p <- plot_dt %>%
    ggplot() +
    ggtitle(plot_title) +
    ylab("Value") +
    xlab("Variable") +
    geom_line(aes(x = name, y = value, group = iter, color = type), alpha = 0.25, linewidth = 1.15) + 
    viridis::scale_color_viridis("Error type", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Error type", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate HMC Pairs Plot
generate_hmc_pairs <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  tmp_summary <- data$summary
  parameter_map <- get_parameter_map()
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "transformed"]
    plot_dt <- data$samples %>%
      .[, .(iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"transformed"])]
  } else {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "raw"]
    plot_dt <- data$samples %>%
      .[, .(iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"raw"])]
  }
  
  # Handle diagnostics
  if (params$diag == "None") {
    plot_dt <- plot_dt[, .(iter, name, value)] %>%
      .[, type := "None"]
  } else if (params$diag == "Max. treedepth") {
    if (tmp_summary$treedepth > 0) {
      max_tree <- max(plot_dt$treedepth)
      plot_dt <- plot_dt[, .(iter, name, value, treedepth)] %>%
        .[, type := "None"] %>%
        .[treedepth == max_tree, type := "Max treedepth"] %>%
        .[, type := factor(type, levels = c("None", "Max treedepth"))]
    } else {
      plot_dt <- plot_dt[, .(iter, name, value, treedepth)] %>%
        .[, type := "None"]
    }
  } else {
    if (tmp_summary$divergent > 0) {
      plot_dt <- plot_dt[, .(iter, name, value, divergent)] %>%
        .[, type := "None"] %>%
        .[divergent == 1, type := "Divergent transition"] %>%
        .[, type := factor(type, levels = c("None", "Divergent transition"))]
    } else {
      plot_dt <- plot_dt[, .(iter, name, value, divergent)] %>%
        .[, type := "None"]
    }
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Reshape for pairs plot
  plot_dt <- plot_dt %>%
    .[, .(type, iter, name, value)] %>%
    dcast(., type + iter ~ name)
  
  # Create pairs plot
  p <- plot_dt %>%
    ggpairs(., columns = 3:ncol(plot_dt), aes(color = type, alpha = 0.4)) + 
    viridis::scale_color_viridis("Error type", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Error type", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate HMC Trace Plot
generate_hmc_trace <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  tmp_summary <- data$summary
  parameter_map <- get_parameter_map()
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "transformed"]
    plot_title <- "Trace (transformed)"
    plot_dt <- data$samples %>%
      .[, .(chain, chain_iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"transformed"])]
  } else {
    grab_parameters <- parameter_map[match(params$leading_params, parameter_map[,"input"]), "raw"]
    plot_title <- "Trace (original)"
    plot_dt <- data$samples %>%
      .[, .(chain, chain_iter, name, value, treedepth, divergent)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"raw"])]
  }
  
  # Handle diagnostics
  if (params$diag == "None") {
    plot_dt <- plot_dt[, .(chain, chain_iter, name, value)] %>%
      .[, type := "None"]
  } else if (params$diag == "Max. treedepth") {
    if (tmp_summary$treedepth > 0) {
      max_tree <- max(plot_dt$treedepth)
      plot_dt <- plot_dt[, .(chain, chain_iter, name, value, treedepth)] %>%
        .[, type := "None"] %>%
        .[treedepth == max_tree, type := "Max treedepth"] %>%
        .[, type := factor(type, levels = c("None", "Max treedepth"))]
    } else {
      plot_dt <- plot_dt[, .(chain, chain_iter, name, value)] %>%
        .[, type := "None"]
    }
  } else {
    if (tmp_summary$divergent > 0) {
      plot_dt <- plot_dt[, .(chain, chain_iter, name, value, divergent)] %>%
        .[, type := "None"] %>%
        .[divergent == 1, type := "Divergent transition"] %>%
        .[, type := factor(type, levels = c("None", "Divergent transition"))]
    } else {
      plot_dt <- plot_dt[, .(chain, chain_iter, name, value)] %>%
        .[, type := "None"]
    }
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Create trace plot
  p <- plot_dt %>%
    ggplot() +
    ggtitle(plot_title) +
    ylab("Value") +
    xlab("Iteration") +
    facet_wrap(~name, scales = "free_y", ncol = min(c(3, uniqueN(plot_dt$name)))) +
    geom_path(aes(x = chain_iter, y = value, color = as.character(chain), group = paste0(chain, ".", name)), linewidth = 1.15)
  
  if (params$diag != "None" & uniqueN(plot_dt$type) > 1) {
    p <- p + geom_point(data = plot_dt[type != "None"], aes(x = chain_iter, y = value, fill = type), shape = 21, size = 2)
  }
  
  p <- p + 
    scale_color_brewer("Chain", palette = "Blues") +
    viridis::scale_fill_viridis("Error type", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE) +
    get_ssp_theme()
  
  return(p)
}

#' Generate HMC R-hat Plot
generate_hmc_rhat <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Check if rhat data is available
  if (is.null(data$rhat)) {
    stop("R-hat data not available. Ensure rhat.csv exists in the model directory.")
  }
  
  # Handle eps parameter
  if (params$eps) {
    parameter_map <- get_parameter_map_extended()
    target_par <- c(params$leading_params, "epsp")
  } else {
    parameter_map <- get_parameter_map()
    target_par <- params$leading_params
  }
  
  # Use loaded rhat data
  plot_dt <- data$rhat
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "transformed"]
    plot_dt <- plot_dt %>%
      .[, .(variable, name, row, rhat)] %>%
      .[name %in% grab_parameters] %>%
      .[!is.na(row), name := paste0(name, "[", row, "]")]
  } else {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "raw"]
    plot_dt <- plot_dt %>%
      .[, .(variable, name, row, rhat)] %>%
      .[name %in% grab_parameters] %>%
      .[!is.na(row), name := paste0(name, "[", row, "]")]
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  rhats <- plot_dt$rhat
  names(rhats) <- plot_dt$name
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create R-hat plot
  p <- bayesplot::mcmc_rhat(rhats, size = 2) +
    get_ssp_theme() +
    theme(
      panel.grid.major.x = element_line(color = 'gray70', linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Generate HMC Effective Sample Size Plot
generate_hmc_neff <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Check if neff data is available
  if (is.null(data$neff)) {
    stop("N_eff data not available. Ensure neff.csv exists in the model directory.")
  }
  
  # Handle eps parameter
  if (params$eps) {
    parameter_map <- get_parameter_map_extended()
    target_par <- c(params$leading_params, "epsp")
  } else {
    parameter_map <- get_parameter_map()
    target_par <- params$leading_params
  }
  
  # Use loaded neff data
  plot_dt <- data$neff
  
  # Parameter mapping
  if (params$raw) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "transformed"]
    plot_dt <- plot_dt %>%
      .[, .(variable, name, row, neff)] %>%
      .[name %in% grab_parameters] %>%
      .[!is.na(row), name := paste0(name, "[", row, "]")]
  } else {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "raw"]
    plot_dt <- plot_dt %>%
      .[, .(variable, name, row, neff)] %>%
      .[name %in% grab_parameters] %>%
      .[!is.na(row), name := paste0(name, "[", row, "]")]
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  neff <- plot_dt$neff
  names(neff) <- plot_dt$name
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create effective sample size plot
  p <- bayesplot::mcmc_neff(neff, size = 2) +
    get_ssp_theme() +
    theme(
      panel.grid.major.x = element_line(color = 'gray70', linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

#' Generate HMC Autocorrelation Plot
generate_hmc_acf <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  parameter_map <- get_parameter_map()
  target_par <- params$leading_params
  
  # Parameter mapping and data preparation
  if (params$raw) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "transformed"]
    plot_df <- data$samples %>%
      .[, .(chain, chain_iter, name, value)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"transformed"])] %>%
      dcast(., chain + chain_iter ~ name) %>%
      .[, chain_iter := NULL] %>%
      setnames(., "chain", "Chain") %>%
      as.data.frame(.)
  } else {
    grab_parameters <- parameter_map[match(target_par, parameter_map[,"input"]), "raw"]
    plot_df <- data$samples %>%
      .[, .(chain, chain_iter, name, value)] %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[,"raw"])] %>%
      dcast(., chain + chain_iter ~ name) %>%
      .[, chain_iter := NULL] %>%
      setnames(., "chain", "Chain") %>%
      as.data.frame(.)
  }
  
  if (nrow(plot_df) == 0) {
    stop("No data available for plotting")
  }
  
  # Set bayesplot color scheme
  bayesplot::color_scheme_set(params$scheme)
  
  # Create autocorrelation plot
  p <- bayesplot::mcmc_acf_bar(plot_df, lags = params$lags) +
    get_ssp_theme() +
    theme(
      panel.grid.minor.y = element_blank()
    ) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "gray70")
  
  return(p)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Set up environment
# set_global_config(year_one = 1994, index_names = c("Index 1", "Index 2"))
# 
# # Generate HMC diagnostic plots
# model_dir <- "./path/to/model1/"
# 
# # Use defaults
# p1 <- generate_hmc_parcoord(model_dir)
# p2 <- generate_hmc_trace(model_dir)
# 
# # Custom parameters
# custom_params <- list(
#   leading_params = c("logK", "r", "sigmap"),
#   raw = FALSE,
#   diag = "Divergences"
# )
# p3 <- generate_hmc_parcoord(model_dir, custom_params)
