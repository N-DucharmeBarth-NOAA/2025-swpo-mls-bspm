# global.R - Load packages and data that are used across UI and server

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(markdown)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggthemes)
library(fresh)
library(bayesplot)
library(GGally)
library(MASS)
library(randtests)
library(loo)
library(rstantools)
# library(DT)

# Source helper functions
helper_files <- list.files("R/helper-fns", full.names = TRUE, pattern = "\\.r$|\\.R$")
if(length(helper_files) > 0) {
  sapply(helper_files, source)
}

# Source plot functions
plot_files <- list.files("R/plot-fns", full.names = TRUE, pattern = "\\.r$|\\.R$")
if(length(plot_files) > 0) {
  sapply(plot_files, source)
}

# Load and format summary data
summary_dt <- fread("data/summary_dt.csv")

# Format numeric columns for display
if("max_rhat" %in% colnames(summary_dt)) {
  summary_dt$max_rhat <- round(summary_dt$max_rhat, digits = 4)
}
if("min_neff" %in% colnames(summary_dt)) {
  summary_dt$min_neff <- round(summary_dt$min_neff)
}

# Discover available model directories
model_stem <- "data/model_runs"
all_dirs <- list.files(model_stem, recursive = FALSE)

# Filter to directories that contain required files
valid_dirs <- all_dirs[sapply(all_dirs, function(x) {
  dir_path <- file.path(model_stem, x)
  file.exists(file.path(dir_path, "hmc_samples.csv")) && 
  file.exists(file.path(dir_path, "fit_summary.csv"))
})]

# Configure global settings for plot functions
if(exists("set_global_config")) {
  set_global_config(
    year_one = 1952,  # Adjust based on your data
    index_names = c("DWFN CPUE", "DWFN CPUE"),  # Adjust based on your indices
    model_stem = model_stem,
    height_per_panel = 350
  )
}

# Helper function to get model directory path from model ID
get_model_path <- function(model_id) {
  file.path(model_stem, model_id)
}
