# Standalone Plot Functions for SSP Model Analysis
# Part 5: Batch Processing and Orchestration

# Note: This requires Parts 1-4 to be loaded first

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
  
  # Create output directory if it doesn't exist
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  model_name <- basename(gsub("/$", "", model_dir))
  cat("Generating plots for model:", model_name, "\n")
  
  plots <- list()
  
  # HMC Diagnostic plots
  cat("  Generating HMC diagnostic plots...\n")
  tryCatch({
    plots$hmc_parcoord <- generate_hmc_parcoord(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_parcoord.", plot_format)), 
             plots$hmc_parcoord, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC parcoord:", e$message, "\n"))
  
  tryCatch({
    plots$hmc_pairs <- generate_hmc_pairs(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_pairs.", plot_format)), 
             plots$hmc_pairs, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC pairs:", e$message, "\n"))
  
  tryCatch({
    plots$hmc_trace <- generate_hmc_trace(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_trace.", plot_format)), 
             plots$hmc_trace, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC trace:", e$message, "\n"))
  
  tryCatch({
    plots$hmc_rhat <- generate_hmc_rhat(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_rhat.", plot_format)), 
             plots$hmc_rhat, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC rhat:", e$message, "\n"))
  
  tryCatch({
    plots$hmc_neff <- generate_hmc_neff(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_neff.", plot_format)), 
             plots$hmc_neff, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC neff:", e$message, "\n"))
  
  tryCatch({
    plots$hmc_acf <- generate_hmc_acf(model_dir, params$hmc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_hmc_acf.", plot_format)), 
             plots$hmc_acf, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating HMC acf:", e$message, "\n"))
  
  # PPC plots
  cat("  Generating PPC plots...\n")
  tryCatch({
    plots$ppc_dens <- generate_ppc_dens(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_dens.", plot_format)), 
             plots$ppc_dens, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC dens:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_ecdf <- generate_ppc_ecdf(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_ecdf.", plot_format)), 
             plots$ppc_ecdf, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC ecdf:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_pit_ecdf <- generate_ppc_pit_ecdf(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_pit_ecdf.", plot_format)), 
             plots$ppc_pit_ecdf, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC PIT ECDF:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_stat <- generate_ppc_stat(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_stat.", plot_format)), 
             plots$ppc_stat, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC stat:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_loo_pit <- generate_ppc_loo_pit(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_loo_pit.", plot_format)), 
             plots$ppc_loo_pit, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC LOO PIT:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_loo_qq <- generate_ppc_loo_qq(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_loo_qq.", plot_format)), 
             plots$ppc_loo_qq, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC LOO QQ:", e$message, "\n"))
  
  tryCatch({
    plots$ppc_loo_interval <- generate_ppc_loo_interval(model_dir, params$ppc)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(model_name, "_ppc_loo_interval.", plot_format)), 
             plots$ppc_loo_interval, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating PPC LOO interval:", e$message, "\n"))
  
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
  
  # Create output directory if it doesn't exist
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Generating multi-model comparison plots...\n")
  
  plots <- list()
  
  # Model fit plots
  cat("  Generating model fit plots...\n")
  tryCatch({
    plots$index_fit <- generate_index_fit(model_dirs, params$fits)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_index_fit.", plot_format)), 
             plots$index_fit, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating index fit:", e$message, "\n"))
  
  tryCatch({
    plots$index_fit_ppd <- generate_index_fit_ppd(model_dirs, params$fits)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_index_fit_ppd.", plot_format)), 
             plots$index_fit_ppd, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating index fit PPD:", e$message, "\n"))
  
  tryCatch({
    plots$index_fit_residuals <- generate_index_fit_residuals(model_dirs, params$fits)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_index_fit_residuals.", plot_format)), 
             plots$index_fit_residuals, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating index fit residuals:", e$message, "\n"))
  
  # Prior-posterior comparisons
  cat("  Generating prior-posterior comparison plots...\n")
  tryCatch({
    plots$ppp <- generate_ppp(model_dirs, params$ppp)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_prior_posterior_params.", plot_format)), 
             plots$ppp, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating prior-posterior parameters:", e$message, "\n"))
  
  tryCatch({
    plots$ppts <- generate_ppts(model_dirs, params$ppts)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_prior_posterior_timeseries.", plot_format)), 
             plots$ppts, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating prior-posterior time series:", e$message, "\n"))
  
  # Management plots
  cat("  Generating management plots...\n")
  tryCatch({
    plots$kobe <- generate_kb(model_dirs, params$kbmj)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_kobe.", plot_format)), 
             plots$kobe, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating Kobe plot:", e$message, "\n"))
  
  tryCatch({
    plots$majuro <- generate_mj(model_dirs, params$kbmj)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_majuro.", plot_format)), 
             plots$majuro, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating Majuro plot:", e$message, "\n"))
  
  # Forecast plots
  cat("  Generating forecast plots...\n")
  tryCatch({
    plots$forecast <- generate_fcast(model_dirs, params$forecasts)
    if (save_plots) {
      ggsave(file.path(output_dir, paste0("comparison_forecasts.", plot_format)), 
             plots$forecast, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) cat("    Error generating forecast plot:", e$message, "\n"))
  
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
generate_all_plots <- function(model_dirs, 
                              output_dir = "./plots/", 
                              params = NULL,
                              save_plots = TRUE,
                              plot_format = "png",
                              width = 12,
                              height = 8,
                              dpi = 300,
                              parallel = FALSE,
                              n_cores = 2) {
  
  if (is.null(params)) params <- get_default_params()
  
  # Create output directories
  if (save_plots) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create subdirectories for each model
    for (model_dir in model_dirs) {
      model_name <- basename(gsub("/$", "", model_dir))
      model_output_dir <- file.path(output_dir, model_name)
      if (!dir.exists(model_output_dir)) {
        dir.create(model_output_dir, recursive = TRUE)
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
  cat("Parallel processing:", ifelse(parallel, paste("Yes (", n_cores, "cores)"), "No"), "\n")
  cat("\n")
  
  all_plots <- list()
  
  # Generate single model plots
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
      # PPC functions
      "generate_ppc_dens", "generate_ppc_ecdf", "generate_ppc_pit_ecdf", 
      "generate_ppc_stat", "generate_ppc_loo_pit", "generate_ppc_loo_qq", "generate_ppc_loo_interval",
      # Single model batch function
      "generate_all_single_model_plots",
      # Helper functions (need to export these too)
      "ssp_extract_cpue_fit", "ssp_calc_likelihood", "ssp_calc_rmse", "ssp_derived_quants", 
      "ssp_derived_quants_ts", "ssp_forecast", "ssp_prior_pushforward",
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
  
  # Generate multi-model comparison plots
  cat("\nPHASE 2: Multi-Model Comparisons\n")
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
  cat("Single model plots generated for", length(model_dirs), "models\n")
  cat("Multi-model comparison plots generated\n")
  if (save_plots) {
    cat("All plots saved to:", output_dir, "\n")
  }
  cat("Plot generation completed successfully!\n")
  
  return(all_plots)
}

#' Create Plot Summary Report
#' @param model_dirs Vector of model directories
#' @param output_dir Directory containing generated plots
#' @param report_file Path for the summary report
create_plot_report <- function(model_dirs, output_dir = "./plots/", report_file = NULL) {
  
  if (is.null(report_file)) {
    report_file <- file.path(output_dir, "plot_summary_report.md")
  }
  
  model_names <- sapply(model_dirs, function(x) basename(gsub("/$", "", x)))
  
  # Create markdown report
  report_content <- c(
    "# SSP Model Analysis Plot Summary Report",
    "",
    paste("Generated on:", Sys.time()),
    paste("Number of models analyzed:", length(model_dirs)),
    "",
    "## Models Analyzed",
    "",
    paste("-", model_names),
    "",
    "## Generated Plots",
    "",
    "### Single Model Diagnostics",
    "",
    "For each model, the following diagnostic plots were generated:",
    "",
    "#### HMC Diagnostics",
    "- **Parallel Coordinates**: Parameter relationships across iterations",
    "- **Pairs Plot**: Pairwise parameter correlations",
    "- **Trace Plot**: Parameter evolution across chains",
    "- **R-hat Plot**: Convergence diagnostic",
    "- **Effective Sample Size**: Sampling efficiency",
    "- **Autocorrelation**: Parameter autocorrelation across lags",
    "",
    "#### Posterior Predictive Checks",
    "- **Density Overlay**: Observed vs predicted distributions",
    "- **ECDF Overlay**: Empirical cumulative distribution comparisons",
    "- **PIT ECDF**: Probability integral transform diagnostics",
    "- **Test Statistics**: Summary statistic comparisons",
    "- **LOO-PIT**: Leave-one-out cross-validation diagnostics",
    "- **LOO QQ**: Leave-one-out quantile-quantile plots",
    "- **LOO Intervals**: Leave-one-out posterior intervals",
    "",
    "### Multi-Model Comparisons",
    "",
    "#### Model Fits",
    "- **Index Fit**: Model predictions vs observed indices",
    "- **Index Fit PPD**: Posterior predictive distributions",
    "- **Index Fit Residuals**: Model residual analysis",
    "",
    "#### Prior-Posterior Comparisons",
    "- **Parameter Distributions**: Prior vs posterior parameter densities",
    "- **Time Series**: Prior vs posterior derived quantities over time",
    "",
    "#### Management Plots",
    "- **Kobe Plot**: Stock status (P/P_MSY vs F/F_MSY)",
    "- **Majuro Plot**: Alternative stock status (Depletion vs F/F_MSY)",
    "",
    "#### Projections",
    "- **Forecasts**: Future projections under different scenarios",
    "",
    "## File Organization",
    "",
    "```",
    paste0(output_dir, "/"),
    paste0("├── ", model_names[1], "/          # Individual model plots"),
    paste0("│   ├── ", model_names[1], "_hmc_*.png"),
    paste0("│   └── ", model_names[1], "_ppc_*.png"),
    if (length(model_names) > 1) paste0("├── ", model_names[2], "/          # Individual model plots") else "",
    if (length(model_names) > 1) paste0("│   ├── ", model_names[2], "_hmc_*.png") else "",
    if (length(model_names) > 1) paste0("│   └── ", model_names[2], "_ppc_*.png") else "",
    "└── comparisons/       # Multi-model comparisons",
    "    ├── comparison_index_fit.png",
    "    ├── comparison_kobe.png",
    "    ├── comparison_majuro.png",
    "    └── comparison_forecasts.png",
    "```",
    "",
    "## Usage Notes",
    "",
    "- All plots use consistent color schemes and themes for easy comparison",
    "- Single model plots focus on model diagnostics and validation",
    "- Multi-model plots enable direct comparison between different model runs",
    "- Plots are saved in high resolution suitable for publication",
    "",
    "---",
    "*Report generated by SSP Model Analysis Plot Functions*"
  )
  
  # Remove empty strings for models that don't exist
  report_content <- report_content[report_content != ""]
  
  writeLines(report_content, report_file)
  cat("Plot summary report saved to:", report_file, "\n")
  
  return(report_file)
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
# # Generate all plots with defaults
# all_plots <- generate_all_plots(
#   model_dirs = model_dirs,
#   output_dir = "./analysis_plots/",
#   save_plots = TRUE,
#   plot_format = "png",
#   parallel = TRUE,
#   n_cores = 2
# )
# 
# # Create summary report
# create_plot_report(
#   model_dirs = model_dirs,
#   output_dir = "./analysis_plots/",
#   report_file = "./analysis_plots/model_comparison_report.md"
# )
# 
# # Custom parameters example
# custom_params <- get_default_params()
# custom_params$hmc$leading_params <- c("logK", "r", "sigmap")
# custom_params$hmc$diag <- "Divergences"
# custom_params$ppc$prop <- 0.5
# custom_params$fits$type <- "Quantile"
# custom_params$fits$quants <- 90
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
