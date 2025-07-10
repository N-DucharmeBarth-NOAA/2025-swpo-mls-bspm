# build-shiny.R - Build script for local Shiny deployment

#' Create Symlink with Cross-Platform Support for Files and Directories
#' 
#' Creates symbolic links with robust error handling and appropriate fallbacks
#' 
#' @param target Path to target file or directory (source)
#' @param link Path where symlink should be created (destination)
create_symlink <- function(target, link) {
  # Detect if target is file or directory
  is_file <- file.exists(target) && !dir.exists(target)
  is_dir <- dir.exists(target)
  
  if (!is_file && !is_dir) {
    warning("Target does not exist: ", target)
    return(FALSE)
  }
  
  # Remove existing link if it exists
  if (file.exists(link) || dir.exists(link)) {
    unlink(link, recursive = TRUE, force = TRUE)
  }
  
  # Ensure parent directory exists
  dir.create(dirname(link), recursive = TRUE, showWarnings = FALSE)
  
  if (.Platform$OS.type == "windows") {
    success <- FALSE
    
    # Choose appropriate mklink flag
    mklink_flag <- if (is_file) "" else "/J"  # No flag for files, /J for directories
    
    # Method 1: Try mklink through cmd
    tryCatch({
      args <- if (is_file) {
        c("/c", "mklink", shQuote(link), shQuote(target))
      } else {
        c("/c", "mklink", "/J", shQuote(link), shQuote(target))
      }
      
      result <- system2("cmd", args = args, stdout = TRUE, stderr = TRUE)
      if (attr(result, "status") == 0 || is.null(attr(result, "status"))) {
        success <- TRUE
        type_name <- if (is_file) "file link" else "junction"
        cat("  ✓ Windows", type_name, "created:", basename(link), "\n")
      }
    }, error = function(e) {
      # Continue to next method
    })
    
    # Method 2: Try with Sys.which
    if (!success) {
      tryCatch({
        cmd_path <- Sys.which("cmd")
        if (cmd_path != "") {
          args <- if (is_file) {
            c("/c", "mklink", shQuote(link), shQuote(target))
          } else {
            c("/c", "mklink", "/J", shQuote(link), shQuote(target))
          }
          
          result <- system2(cmd_path, args = args, stdout = TRUE, stderr = TRUE)
          if (attr(result, "status") == 0 || is.null(attr(result, "status"))) {
            success <- TRUE
            type_name <- if (is_file) "file link" else "junction"
            cat("  ✓ Windows", type_name, "created:", basename(link), "\n")
          }
        }
      }, error = function(e) {
        # Continue to fallback
      })
    }
    
    # Fallback: Copy file or directory
    if (!success) {
      type_name <- if (is_file) "file" else "directory"
      cat("  ! Symlink creation failed, copying", type_name, "instead:", basename(link), "\n")
      
      tryCatch({
        if (is_file) {
          file.copy(target, link)
        } else {
          file.copy(target, dirname(link), recursive = TRUE)
          if (basename(target) != basename(link)) {
            file.rename(file.path(dirname(link), basename(target)), link)
          }
        }
        success <- TRUE
        cat("  ✓", stringr::str_to_title(type_name), "copied:", basename(link), "\n")
      }, error = function(e) {
        warning("Failed to create symlink or copy ", type_name, ": ", e$message)
        return(FALSE)
      })
    }
    
    return(success)
    
  } else {
    # Unix/Linux/Mac: Use standard symlink (works for both files and directories)
    tryCatch({
      file.symlink(target, link)
      type_name <- if (is_file) "file" else "directory"
      cat("  ✓ Symlink created for", type_name, ":", basename(link), "\n")
      return(TRUE)
    }, error = function(e) {
      warning("Failed to create symlink: ", e$message)
      
      # Fallback: Copy
      type_name <- if (is_file) "file" else "directory"
      cat("  ! Symlink creation failed, copying", type_name, "instead:", basename(link), "\n")
      
      tryCatch({
        if (is_file) {
          file.copy(target, link)
        } else {
          file.copy(target, dirname(link), recursive = TRUE)
          if (basename(target) != basename(link)) {
            file.rename(file.path(dirname(link), basename(target)), link)
          }
        }
        cat("  ✓", stringr::str_to_title(type_name), "copied:", basename(link), "\n")
        return(TRUE)
      }, error = function(e2) {
        warning("Failed to copy ", type_name, ": ", e2$message)
        return(FALSE)
      })
    })
  }
}

#' Build Local Shiny Application
#' 
#' Creates a self-contained Shiny app directory with all necessary files
#' for local deployment using your existing plot functions and data
#' 
#' @param proj_dir Project root directory (default: current project)
#' @param shiny_dir Target shiny directory (default: "shiny")
#' @param force Whether to overwrite existing shiny directory
build_local_shiny <- function(proj_dir = this.path::this.proj(), 
                             shiny_dir = "shiny", 
                             force = FALSE) {
  
  library(data.table)
  library(magrittr)
  
  # Setup paths
  shiny_path <- file.path(proj_dir, shiny_dir)
  code_dir <- file.path(proj_dir, "code", "R")
  data_dir <- file.path(proj_dir, "data", "output", "model_runs")
  
  # Check if shiny directory exists
  if (dir.exists(shiny_path) && !force) {
    stop("Shiny directory already exists. Use force = TRUE to overwrite.")
  }
  
  # Create shiny directory structure
  if (dir.exists(shiny_path) && force) {
    unlink(shiny_path, recursive = TRUE)
  }
  
  dir.create(shiny_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(shiny_path, "R"), showWarnings = FALSE)
  dir.create(file.path(shiny_path, "data","output"), showWarnings = FALSE)
  
  cat("Creating Shiny application structure...\n")
  
  # Create symlinks to plot functions
  cat("Creating symlinks to plot functions...\n")
  if (dir.exists(file.path(code_dir, "plot-fns"))) {
    create_symlink(normalizePath(file.path(code_dir, "plot-fns")), 
                   file.path(shiny_path, "R", "plot-fns"))
  } else {
    warning("Plot functions directory not found at: ", file.path(code_dir, "plot-fns"))
  }
  
  # Create symlinks to helper functions
  cat("Creating symlinks to helper functions...\n")
  if (dir.exists(file.path(code_dir, "helper-fns"))) {
    create_symlink(normalizePath(file.path(code_dir, "helper-fns")), 
                   file.path(shiny_path, "R", "helper-fns"))
  } else {
    warning("Helper functions directory not found at: ", file.path(code_dir, "helper-fns"))
  }

  # Create symlinks to shiny functions
  cat("Creating symlinks to Shiny app code...\n")
  shiny_app_dir <- file.path(code_dir, "shiny-app")
  if (dir.exists(shiny_app_dir)) {
    # Get all files in the shiny-app directory
    shiny_files <- list.files(shiny_app_dir, full.names = TRUE, recursive = FALSE)
    # Filter to only files (not subdirectories)
    shiny_files <- shiny_files[file.info(shiny_files)$isdir == FALSE]
    
    if (length(shiny_files) > 0) {
      cat("  Found", length(shiny_files), "files to link\n")
      
      # Create symlinks for each file
      success_count <- 0
      for (file_path in shiny_files) {
        file_name <- basename(file_path)
        target_path <- normalizePath(file_path)
        link_path <- file.path(shiny_path, file_name)
        
        if (create_symlink(target_path, link_path)) {
          success_count <- success_count + 1
        }
      }
      
      cat("  Successfully linked", success_count, "of", length(shiny_files), "files\n")
    } else {
      warning("No files found in shiny-app directory")
    }
  } else {
    warning("Shiny app code directory not found at: ", shiny_app_dir)
  }
  
  # Create symlink to model outputs (local deployment)
  cat("Creating data symlinks...\n")
  if (dir.exists(data_dir)) {
    create_symlink(normalizePath(data_dir), 
                   file.path(shiny_path, "data","output", "model_runs"))
  } else {
    warning("Model runs directory not found at: ", data_dir)
  }
  
  # Generate summary data for the app
  cat("Generating summary data...\n")
  create_summary_data(proj_dir, shiny_path)
  
  cat("Shiny application built successfully at:", shiny_path, "\n")
  cat("To run the app, navigate to the shiny directory and run app.R\n")
  
  return(invisible(shiny_path))
}

#' Parse Run Number Components with Robust Pattern Matching
#' 
#' Extracts model configuration components from run_number using multiple
#' parsing strategies to handle different naming conventions
#' 
#' @param run_number Character string containing the run identifier
#' @return List with extracted components and parsing method used
parse_run_components <- function(run_number) {
  # Initialize result with defaults - EXPANDED for new fields
  result <- list(
    run_num = NA_character_,
    cpue_index = NA_character_,
    sigma_catch = NA_real_,
    sigma_edev = NA_real_,
    n_step = NA_integer_,
    exec = NA_character_,           # NEW: executable type (B, BF, BX, BFX)
    start_year = NA_integer_,       # NEW: start year
    catch_scenario = NA_character_, # NEW: catch scenario (e.g., "0.2flat")
    step_scenario = NA_character_,  # NEW: step scenario (e.g., "5reg")
    parsing_method = "unknown"
  )
  
  # Extract run number (always first numeric part)
  run_match <- regexpr("^([0-9]+)", run_number)
  if (run_match > 0) {
    result$run_num <- regmatches(run_number, run_match)
  }
  
  # Extract everything after first dash as the stem
  stem_match <- regexpr("^[0-9]+-(.*)$", run_number)
  if (stem_match < 0) return(result)
  
  stem <- sub("^[0-9]+-", "", run_number)
  
  # Strategy 1: NEW Complex scenario pattern (cpue-exe{X}-sy{YEAR}-cs{SCENARIO}-e{VAL}-ss{SCENARIO})
  complex_pattern <- "^([^-]+)-exe([A-Z]+)-sy([0-9]+)-cs([^-]+)-e([0-9.]+)-ss([^_]+)"
  if (grepl(complex_pattern, stem)) {
    matches <- regmatches(stem, regexec(complex_pattern, stem))[[1]]
    if (length(matches) == 7) {
      result$cpue_index <- matches[2]
      result$exec <- matches[3]
      result$start_year <- as.integer(matches[4])
      result$catch_scenario <- matches[5]
      result$sigma_edev <- as.numeric(matches[6])
      result$step_scenario <- matches[7]
      
      # Extract numeric part from catch scenario if present
      catch_num_match <- regexpr("^([0-9.]+)", result$catch_scenario)
      if (catch_num_match > 0) {
        result$sigma_catch <- as.numeric(regmatches(result$catch_scenario, catch_num_match))
      }
      
      # Extract numeric part from step scenario if present
      step_num_match <- regexpr("^([0-9]+)", result$step_scenario)
      if (step_num_match > 0) {
        result$n_step <- as.integer(regmatches(result$step_scenario, step_num_match))
      }
      
      result$parsing_method <- "complex_scenario"
      return(result)
    }
  }
  
  # Strategy 2: Modern systematic pattern (cpue-c{val}-e{val}-s{val})
  # Translate to new format using provided guidance
  systematic_pattern <- "^([^-]+)-c([0-9.]+)-e([0-9.]+)-s([0-9]+)"
  if (grepl(systematic_pattern, stem)) {
    matches <- regmatches(stem, regexec(systematic_pattern, stem))[[1]]
    if (length(matches) == 5) {
      result$cpue_index <- matches[2]
      result$sigma_catch <- as.numeric(matches[3])
      result$sigma_edev <- as.numeric(matches[4])
      result$n_step <- as.integer(matches[5])
      
      # Translate to new format structure
      result$exec <- "B"                                    # Default for old systematic
      result$start_year <- 1952L                           # Default for old systematic  
      result$catch_scenario <- paste0(result$sigma_catch, "flat")
      result$step_scenario <- paste0(result$n_step, "reg")
      
      result$parsing_method <- "systematic"
      return(result)
    }
  }
  
  # Strategy 3: Legacy prior pattern (e.g., "2024cpueFPrior_0")
  prior_pattern <- "^([0-9]+)cpue([A-Za-z]+)Prior"
  if (grepl(prior_pattern, stem)) {
    matches <- regmatches(stem, regexec(prior_pattern, stem))[[1]]
    if (length(matches) == 3) {
      result$cpue_index <- tolower(matches[3])  # Convert to lowercase for consistency
      result$parsing_method <- "legacy_prior"
      return(result)
    }
  }
  
  # Strategy 4: Effort pattern (e.g., "2024cpueEffortQeff_0")
  effort_pattern <- "^([0-9]+)cpue([A-Za-z]+)([A-Za-z]+)"
  if (grepl(effort_pattern, stem)) {
    matches <- regmatches(stem, regexec(effort_pattern, stem))[[1]]
    if (length(matches) == 4) {
      result$cpue_index <- tolower(matches[3])
      result$parsing_method <- "legacy_effort"
      return(result)
    }
  }
  
  # Strategy 5: Fallback - extract any recognizable components
  # Look for individual patterns anywhere in the string
  
  # Extract cpue index if it's a known value
  known_indices <- c("au", "nz", "obs", "dwfn", "effort", "obsNoPF", "obsPFonly")
  for (idx in known_indices) {
    if (grepl(paste0("\\b", idx, "\\b"), stem, ignore.case = TRUE)) {
      result$cpue_index <- idx
      break
    }
  }
  
  # Extract executive type if present
  exec_match <- regexpr("exe([A-Z]+)", stem)
  if (exec_match > 0) {
    result$exec <- sub("exe([A-Z]+).*", "\\1", regmatches(stem, exec_match))
  }
  
  # Extract start year if present
  sy_match <- regexpr("sy([0-9]+)", stem)
  if (sy_match > 0) {
    result$start_year <- as.integer(sub("sy([0-9]+).*", "\\1", regmatches(stem, sy_match)))
  }
  
  # Extract catch scenario if present
  cs_match <- regexpr("cs([^-]+)", stem)
  if (cs_match > 0) {
    result$catch_scenario <- sub("cs([^-]+).*", "\\1", regmatches(stem, cs_match))
    # Try to extract numeric part
    catch_num_match <- regexpr("^([0-9.]+)", result$catch_scenario)
    if (catch_num_match > 0) {
      result$sigma_catch <- as.numeric(regmatches(result$catch_scenario, catch_num_match))
    }
  }
  
  # Extract step scenario if present
  ss_match <- regexpr("ss([^_]+)", stem)
  if (ss_match > 0) {
    result$step_scenario <- sub("ss([^_]+).*", "\\1", regmatches(stem, ss_match))
    # Try to extract numeric part
    step_num_match <- regexpr("^([0-9]+)", result$step_scenario)
    if (step_num_match > 0) {
      result$n_step <- as.integer(regmatches(result$step_scenario, step_num_match))
    }
  }
  
  # Extract catch error if present (standalone)
  if (is.na(result$sigma_catch)) {
    catch_match <- regexpr("c([0-9.]+)", stem)
    if (catch_match > 0) {
      result$sigma_catch <- as.numeric(sub("c([0-9.]+).*", "\\1", 
                                         regmatches(stem, catch_match)))
    }
  }
  
  # Extract effort dev if present (standalone)
  if (is.na(result$sigma_edev)) {
    edev_match <- regexpr("e([0-9.]+)", stem)
    if (edev_match > 0) {
      result$sigma_edev <- as.numeric(sub("e([0-9.]+).*", "\\1", 
                                        regmatches(stem, edev_match)))
    }
  }
  
  # Extract step size if present (standalone)
  if (is.na(result$n_step)) {
    step_match <- regexpr("s([0-9]+)", stem)
    if (step_match > 0) {
      result$n_step <- as.integer(sub("s([0-9]+).*", "\\1", 
                                    regmatches(stem, step_match)))
    }
  }
  
  result$parsing_method <- "fallback"
  return(result)
}

#' Create Summary Data for Shiny App with Enhanced Parsing
#' 
#' Generates a summary CSV file from all model runs for the model selection interface
#' with robust parsing of run components and comprehensive error handling
#' 
#' @param proj_dir Project root directory
#' @param shiny_path Target shiny application directory
create_summary_data <- function(proj_dir, shiny_path) {
  
  # Path to model runs
  model_runs_path <- file.path(proj_dir, "data", "output", "model_runs")
  
  if (!dir.exists(model_runs_path)) {
    warning("Model runs directory not found. Creating empty summary.")
    summary_dt <- data.table(model_id = character(0))
    fwrite(summary_dt, file.path(shiny_path, "data", "summary_dt.csv"))
    return(invisible())
  }
  
  # Use proven directory discovery logic from 06-script
  all_dirs <- list.files(model_runs_path, recursive = TRUE)
  all_dirs <- all_dirs[grep("fit_summary.csv", all_dirs, fixed = TRUE)]
  all_dirs <- gsub("fit_summary.csv", "", all_dirs, fixed = TRUE)
  if (length(grep("-ppc", all_dirs, fixed = TRUE)) > 0) {
    all_dirs <- all_dirs[-grep("-ppc", all_dirs, fixed = TRUE)]
  }
  
  if (length(all_dirs) == 0) {
    warning("No valid model directories found. Creating empty summary.")
    summary_dt <- data.table(model_id = character(0))
    fwrite(summary_dt, file.path(shiny_path, "data", "summary_dt.csv"))
    return(invisible())
  }
  
  # Read and combine summaries with robust parsing
  summary_df.list <- lapply(all_dirs, function(x) {
    dt <- fread(file.path(model_runs_path, x, "fit_summary.csv"))
    dt[, model_id := run_label]
    
    # Parse run components for each row
    parsed_components <- lapply(dt$run_number, parse_run_components)
    
    # Convert to data.table and bind
    components_dt <- rbindlist(parsed_components, fill = TRUE)
    dt <- cbind(dt, components_dt)
    
    return(dt)
  })
  
  # Combine all summaries and filter
  summary_dt <- rbindlist(summary_df.list, fill = TRUE) %>%
    .[run_retro == 0]
  
  # Log parsing success/failure rates for monitoring
  if ("parsing_method" %in% colnames(summary_dt)) {
    parsing_summary <- summary_dt[, .N, by = parsing_method]
    cat("Parsing methods used:\n")
    print(parsing_summary)
    
    failed_parsing <- summary_dt[parsing_method == "unknown", run_number]
    if (length(failed_parsing) > 0) {
      cat("Failed to parse:", paste(failed_parsing, collapse = ", "), "\n")
    }
  }
  
  # Select columns with parsed components first
  select_cols <- c("run_label","run_id", 
                   "n_par", "low_bfmi", "divergent", "treedepth", 
                   "max_rhat", "min_neff", "median_catch_rmse", 
                   "index_rmse_1", "index_rmse_2", "index_rmse_3", "index_rmse_4", "index_rmse_5", "index_rmse_6")
  
  # Only select columns that actually exist (handles missing parsed columns gracefully)
  available_cols <- intersect(select_cols, colnames(summary_dt))
  summary_dt <- summary_dt[, .SD, .SDcols = available_cols]
  
  # Write summary to shiny data directory
  fwrite(summary_dt, file.path(shiny_path, "data", "summary_dt.csv"))
  
  cat("Summary data created with", nrow(summary_dt), "rows from", length(all_dirs), "models\n")
  cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
}

#' Validate Shiny Directory Structure
#' 
#' Checks that all required files are present for the Shiny app to run
validate_shiny_structure <- function(shiny_path) {
  
  required_files <- c(
    "app.R",
    "ui.R", 
    "server.R",
    "global.R",
    "introduction.md",
    "data/summary_dt.csv"
  )
  
  required_dirs <- c(
    "R",
    "data",
    "data/output/model_runs"
  )
  
  missing_files <- character(0)
  missing_dirs <- character(0)
  
  # Check files
  for (file in required_files) {
    if (!file.exists(file.path(shiny_path, file))) {
      missing_files <- c(missing_files, file)
    }
  }
  
  # Check directories
  for (dir in required_dirs) {
    if (!dir.exists(file.path(shiny_path, dir))) {
      missing_dirs <- c(missing_dirs, dir)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("Missing files:\n")
    cat(paste(" -", missing_files, collapse = "\n"), "\n")
  }
  
  if (length(missing_dirs) > 0) {
    cat("Missing directories:\n")
    cat(paste(" -", missing_dirs, collapse = "\n"), "\n")
  }
  
  if (length(missing_files) == 0 && length(missing_dirs) == 0) {
    cat("✓ All required files and directories present\n")
    return(TRUE)
  } else {
    cat("✗ Missing required files or directories\n")
    return(FALSE)
  }
}

#' Run Build Process
#' 
#' Convenience function to build and validate the Shiny app
build_and_validate <- function(proj_dir = this.path::this.proj(), force = FALSE) {
  
  # Build the app
  shiny_path <- build_local_shiny(proj_dir, force = force)
  
  # Validate structure
  cat("\nValidating Shiny app structure...\n")
  is_valid <- validate_shiny_structure(shiny_path)
  
  if (is_valid) {
    cat("\n✓ Shiny app ready to run!\n")
    cat("To start the app:\n")
    cat("1. Navigate to:", shiny_path, "\n")
    cat("2. Run: source('app.R') or click 'Run App' in RStudio\n")
  } else {
    cat("\n✗ Shiny app setup incomplete. Please check missing files.\n")
  }
  
  return(invisible(shiny_path))
}

# Example usage:
# build_and_validate()
# 
# Or step by step:
# build_local_shiny()
# validate_shiny_structure("shiny")
