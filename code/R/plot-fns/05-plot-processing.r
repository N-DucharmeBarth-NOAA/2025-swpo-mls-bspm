# Standalone Plot Functions for SSP Model Analysis
# Part 5: Batch Processing and Orchestration
# Note: This requires Parts 1-4 to be loaded first
# Generates comprehensive diagnostic plots for CPUE and catch data

# =============================================================================
# BATCH PROCESSING FUNCTIONS
# =============================================================================

#' Generate All Single Model Plots
#' @param model_dir Directory containing model output files
#' @param output_dir Directory to save plots
#' @param params Custom parameters (optional)
#' @param save_plots Whether to save plots to files
#' @param plot_format Format for saved plots (png, pdf, svg)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for raster formats
generate_all_single_model_plots <- function(model_dir, 
                                           output_dir = "./plots/", 
                                           params = NULL,
                                           save_plots = TRUE,
                                           plot_format = "png",
                                           width = 12,
                                           height = 8,
                                           dpi = 300) {
  
  if (is.null(params)) params <- get_default_params()
  if (save_plots && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  model_name <- basename(gsub("/$", "", model_dir))
  cat("Generating plots for model:", model_name, "\n")
  
  # Helper function to generate and save plots
  make_plot <- function(func, name, params_arg, prefix = "") {
    tryCatch({
      plot <- func(model_dir, params_arg)
      if (save_plots) {
        filename <- paste0(model_name, "_", prefix, name, ".", plot_format)
        ggsave(file.path(output_dir, filename), plot, width = width, height = height, dpi = dpi)
      }
      return(plot)
    }, error = function(e) {
      cat("    Error generating", name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  plots <- list()
  
  # HMC plots
  cat("  Generating HMC plots...\n")
  hmc_plots <- list(
    parcoord = generate_hmc_parcoord, pairs = generate_hmc_pairs, trace = generate_hmc_trace,
    rhat = generate_hmc_rhat, neff = generate_hmc_neff, acf = generate_hmc_acf
  )
  for (name in names(hmc_plots)) {
    plots[[paste0("hmc_", name)]] <- make_plot(hmc_plots[[name]], paste0("hmc_", name), params$hmc)
  }
  
  # CPUE PPC plots
  cat("  Generating CPUE PPC plots...\n")
  cpue_ppc_plots <- list(
    dens = generate_ppc_dens, ecdf = generate_ppc_ecdf, pit_ecdf = generate_ppc_pit_ecdf,
    stat = generate_ppc_stat, loo_pit = generate_ppc_loo_pit, 
    loo_qq = generate_ppc_loo_qq, loo_interval = generate_ppc_loo_interval
  )
  for (name in names(cpue_ppc_plots)) {
    plots[[paste0("cpue_ppc_", name)]] <- make_plot(cpue_ppc_plots[[name]], paste0("ppc_", name), params$ppc, "cpue_")
  }
  
  # Catch PPC plots
  cat("  Generating Catch PPC plots...\n")
  catch_ppc_plots <- list(
    dens = generate_catch_ppc_dens, ecdf = generate_catch_ppc_ecdf, pit_ecdf = generate_catch_ppc_pit_ecdf,
    stat = generate_catch_ppc_stat, loo_pit = generate_catch_ppc_loo_pit,
    loo_qq = generate_catch_ppc_loo_qq, loo_interval = generate_catch_ppc_loo_interval
  )
  for (name in names(catch_ppc_plots)) {
    plots[[paste0("catch_ppc_", name)]] <- make_plot(catch_ppc_plots[[name]], paste0("ppc_", name), params$ppc, "catch_")
  }
  
  cat("  Completed plots for model:", model_name, "\n")
  return(plots)
}

#' Generate All Multi-Model Comparison Plots
#' @param model_dirs Vector of model directories
#' @param output_dir Directory to save plots
#' @param params Custom parameters (optional)
#' @param save_plots Whether to save plots to files
#' @param plot_format Format for saved plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for raster formats
generate_all_multi_model_plots <- function(model_dirs, 
                                          output_dir = "./plots/", 
                                          params = NULL,
                                          save_plots = TRUE,
                                          plot_format = "png",
                                          width = 12,
                                          height = 8,
                                          dpi = 300) {
  
  if (is.null(params)) params <- get_default_params()
  if (save_plots && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("Generating multi-model comparison plots...\n")
  
  # Helper function for multi-model plots
  make_comparison_plot <- function(func, name, params_arg) {
    tryCatch({
      plot <- func(model_dirs, params_arg)
      if (save_plots) {
        ggsave(file.path(output_dir, paste0("comparison_", name, ".", plot_format)), 
               plot, width = width, height = height, dpi = dpi)
      }
      return(plot)
    }, error = function(e) {
      cat("    Error generating", name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  plots <- list()
  
  # Model fit plots
  cat("  Generating model fit plots...\n")
  fit_plots <- list(
    index_fit = generate_index_fit, index_fit_ppd = generate_index_fit_ppd, 
    index_fit_residuals = generate_index_fit_residuals,
    catch_fit = generate_catch_fit, catch_fit_ppd = generate_catch_fit_ppd,
    catch_fit_residuals = generate_catch_fit_residuals
  )
  for (name in names(fit_plots)) {
    plots[[name]] <- make_comparison_plot(fit_plots[[name]], name, params$fits)
  }
  
  # Prior-posterior and management plots
  cat("  Generating other comparison plots...\n")
  other_plots <- list(
    prior_posterior_params = list(func = generate_ppp, params = params$ppp),
    prior_posterior_timeseries = list(func = generate_ppts, params = params$ppts),
    kobe = list(func = generate_kb, params = params$kbmj),
    majuro = list(func = generate_mj, params = params$kbmj),
    forecasts = list(func = generate_fcast, params = params$forecasts)
  )
  for (name in names(other_plots)) {
    plots[[name]] <- make_comparison_plot(other_plots[[name]]$func, name, other_plots[[name]]$params)
  }
  
  cat("  Completed multi-model comparison plots\n")
  return(plots)
}

#' Generate All Plots for Multiple Models
#' @param model_dirs Vector of model directories
#' @param output_dir Base output directory
#' @param params Custom parameters (optional)
#' @param save_plots Whether to save plots to files
#' @param plot_format Format for saved plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for raster formats
#' @param parallel Whether to use parallel processing for single model plots
#' @param n_cores Number of cores for parallel processing
#' @param comparison_only Whether to generate only multi-model comparison plots (skip single model plots)
generate_all_plots <- function(model_dirs, 
                              output_dir = "./plots/", 
                              params = NULL,
                              save_plots = TRUE,
                              plot_format = "png",
                              width = 12,
                              height = 8,
                              dpi = 300,
                              parallel = FALSE,
                              n_cores = 2,
                              comparison_only = FALSE) {
  
  if (is.null(params)) params <- get_default_params()
  
  # Create output directories
  if (save_plots) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create subdirectories for each model (only if not comparison_only)
    if (!comparison_only) {
      for (model_dir in model_dirs) {
        model_name <- basename(gsub("/$", "", model_dir))
        model_output_dir <- file.path(output_dir, model_name)
        if (!dir.exists(model_output_dir)) {
          dir.create(model_output_dir, recursive = TRUE)
        }
      }
    }
    
    # Create comparison directory
    comparison_dir <- file.path(output_dir, "comparisons")
    if (!dir.exists(comparison_dir)) {
      dir.create(comparison_dir, recursive = TRUE)
    }
  }
  
  cat("=== SSP Model Analysis Plot Generation ===\n")
  cat("Models to process:", length(model_dirs), "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Plot format:", plot_format, "\n")
  cat("Mode:", ifelse(comparison_only, "Comparison plots only", "All plots"), "\n")
  if (!comparison_only) {
    cat("Parallel processing:", ifelse(parallel, paste("Yes (", n_cores, "cores)"), "No"), "\n")
  }
  cat("\n")
  
  all_plots <- list()
  
  # Generate single model plots (only if not comparison_only)
  if (!comparison_only) {
    cat("PHASE 1: Single Model Diagnostics\n")
    cat("===================================\n")
    
    if (parallel && requireNamespace("parallel", quietly = TRUE)) {
      cat("Using parallel processing...\n")
      cl <- parallel::makeCluster(n_cores)
      
      # Export necessary objects to cluster
      parallel::clusterEvalQ(cl, {
        library(data.table)
        library(magrittr)
        library(ggplot2)
        library(viridis)
        library(bayesplot)
        library(GGally)
        library(MASS)
        library(randtests)
      })
      
      # Export all custom functions to cluster
      parallel::clusterExport(cl, c(
        # Utility functions
        "get_default_params", "load_model_data", "get_parameter_map", "get_ssp_theme",
        # HMC functions
        "generate_hmc_parcoord", "generate_hmc_pairs", "generate_hmc_trace", 
        "generate_hmc_rhat", "generate_hmc_neff", "generate_hmc_acf",
        # CPUE PPC functions
        "generate_ppc_dens", "generate_ppc_ecdf", "generate_ppc_pit_ecdf", 
        "generate_ppc_stat", "generate_ppc_loo_pit", "generate_ppc_loo_qq", "generate_ppc_loo_interval",
        # Catch PPC functions
        "generate_catch_ppc_dens", "generate_catch_ppc_ecdf", "generate_catch_ppc_pit_ecdf",
        "generate_catch_ppc_stat", "generate_catch_ppc_loo_pit", "generate_catch_ppc_loo_qq", 
        "generate_catch_ppc_loo_interval",
        # Single model batch function
        "generate_all_single_model_plots",
        # Helper functions (need to export these too)
        "ssp_extract_cpue_fit", "ssp_extract_catch_fit", "ssp_calc_likelihood", "ssp_calc_catch_likelihood",
        "ssp_calc_rmse", "ssp_derived_quants", "ssp_derived_quants_ts", "ssp_forecast", "ssp_prior_pushforward",
        # Global variables
        "year_one", "index_names", "model_stem", "height_per_panel"
      ), envir = environment())
      
      single_model_plots <- parallel::parLapply(cl, model_dirs, function(model_dir) {
        model_name <- basename(gsub("/$", "", model_dir))
        model_output_dir <- ifelse(save_plots, file.path(output_dir, model_name), NULL)
        
        generate_all_single_model_plots(
          model_dir = model_dir,
          output_dir = model_output_dir,
          params = params,
          save_plots = save_plots,
          plot_format = plot_format,
          width = width,
          height = height,
          dpi = dpi
        )
      })
      
      parallel::stopCluster(cl)
      names(single_model_plots) <- sapply(model_dirs, function(x) basename(gsub("/$", "", x)))
      
    } else {
      single_model_plots <- list()
      for (model_dir in model_dirs) {
        model_name <- basename(gsub("/$", "", model_dir))
        model_output_dir <- ifelse(save_plots, file.path(output_dir, model_name), NULL)
        
        single_model_plots[[model_name]] <- generate_all_single_model_plots(
          model_dir = model_dir,
          output_dir = model_output_dir,
          params = params,
          save_plots = save_plots,
          plot_format = plot_format,
          width = width,
          height = height,
          dpi = dpi
        )
      }
    }
    
    all_plots$single_model <- single_model_plots
  } else {
    cat("SKIPPING: Single Model Diagnostics (comparison_only = TRUE)\n")
    cat("========================================================\n")
  }
  
  # Generate multi-model comparison plots
  phase_num <- ifelse(comparison_only, "PHASE 1", "PHASE 2")
  cat("\n", phase_num, ": Multi-Model Comparisons\n")
  cat("=================================\n")
  
  comparison_output_dir <- ifelse(save_plots, file.path(output_dir, "comparisons"), NULL)
  
  multi_model_plots <- generate_all_multi_model_plots(
    model_dirs = model_dirs,
    output_dir = comparison_output_dir,
    params = params,
    save_plots = save_plots,
    plot_format = plot_format,
    width = width,
    height = height,
    dpi = dpi
  )
  
  all_plots$multi_model <- multi_model_plots
  
  # Summary
  cat("\n=== SUMMARY ===\n")
  if (!comparison_only) {
    cat("Single model plots generated for", length(model_dirs), "models\n")
  }
  cat("Multi-model comparison plots generated\n")
  if (save_plots) {
    cat("All plots saved to:", output_dir, "\n")
  }
  cat("Plot generation completed successfully!\n")
  
  return(all_plots)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Complete workflow example
# 
# # Set up environment
# set_global_config(
#   year_one = 1994, 
#   index_names = c("Longline CPUE", "Drift Net CPUE"),
#   model_stem = "./model_outputs/"
# )
# 
# # Define model directories
# model_dirs <- c(
#   "./model_outputs/model1_baseline/",
#   "./model_outputs/model2_alternative/"
# )
# 
# # Generate all plots with defaults (now includes catch plots)
# all_plots <- generate_all_plots(
#   model_dirs = model_dirs,
#   output_dir = "./analysis_plots/",
#   save_plots = TRUE,
#   plot_format = "png",
#   parallel = TRUE,
#   n_cores = 2
# )
# 
# # Custom parameters example
# custom_params <- get_default_params()
# custom_params$hmc$leading_params <- c("logK", "r", "sigmap")
# custom_params$hmc$diag <- "Divergences"
# custom_params$ppc$prop <- 0.5
# custom_params$ppc$stat <- c("mean", "sd")  # For both CPUE and catch PPC
# custom_params$fits$type <- "Quantile"
# custom_params$fits$quants <- 90
# custom_params$fits$resid <- "Standardized"  # Applies to both CPUE and catch
# custom_params$kbmj$uncertainty <- TRUE
# custom_params$forecasts$nyears <- 10
# custom_params$forecasts$type <- "MSY"
# 
# # Generate plots with custom parameters
# custom_plots <- generate_all_plots(
#   model_dirs = model_dirs,
#   output_dir = "./custom_analysis_plots/",
#   params = custom_params,
#   save_plots = TRUE,
#   plot_format = "pdf",
#   width = 14,
#   height = 10
# )
# 
# # Generate only single model plots for specific model
# single_plots <- generate_all_single_model_plots(
#   model_dir = "./model_outputs/model1_baseline/",
#   output_dir = "./single_model_diagnostics/",
#   save_plots = TRUE
# )
# 
# # Generate only multi-model comparisons
# comparison_plots <- generate_all_multi_model_plots(
#   model_dirs = model_dirs,
#   output_dir = "./comparisons_only/",
#   save_plots = TRUE
# )
