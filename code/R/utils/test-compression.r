# Test script for HMC samples compression and expansion
# Tests the condense -> expand workflow to ensure fidelity

library(data.table)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    


#' Test HMC compression and expansion workflow
#' @param model_dir Directory containing hmc_samples.csv, fit_summary.csv, stan_data.csv
#' @param tolerance Numerical tolerance for comparing values (default 1e-10)
#' @return List with test results
test_hmc_workflow <- function(model_dir, tolerance = 1e-10) {
  
  cat("=== Testing HMC Compression/Expansion Workflow ===\n")
  cat("Model directory:", model_dir, "\n\n")
  
  # Check if required files exist
  required_files <- c("hmc_samples.csv", "fit_summary.csv", "stan_data.csv")
  file_paths <- file.path(model_dir, required_files)
  missing_files <- required_files[!file.exists(file_paths)]
  
  if (length(missing_files) > 0) {
    stop("Missing files: ", paste(missing_files, collapse = ", "))
  }
  
  # Load original data
  cat("1. Loading original hmc_samples.csv...\n")
  original_samples <- fread(file.path(model_dir, "hmc_samples.csv"))
  cat("   Original rows:", format(nrow(original_samples), big.mark = ","), "\n")
  cat("   Original variables:", length(unique(original_samples$name)), "\n")
  
  # Test compression
  cat("\n2. Testing compression...\n")
  start_time <- Sys.time()
  condensed_samples <- condense_hmc_samples(model_dir)
  compress_time <- as.numeric(Sys.time() - start_time)
  
  cat("   Condensed rows:", format(nrow(condensed_samples), big.mark = ","), "\n")
  cat("   Condensed variables:", length(unique(condensed_samples$name)), "\n")
  cat("   Compression time:", round(compress_time, 2), "seconds\n")
  
  # Check compression file was created
  compact_file <- file.path(model_dir, "hmc_samples_compact.rda")
  if (!file.exists(compact_file)) {
    stop("hmc_samples_compact.rda was not created")
  }
  
  # Test expansion
  cat("\n3. Testing expansion...\n")
  start_time <- Sys.time()
  expanded_samples <- expand_hmc_samples(model_dir)
  expand_time <- as.numeric(Sys.time() - start_time)
  
  cat("   Expanded rows:", format(nrow(expanded_samples), big.mark = ","), "\n")
  cat("   Expanded variables:", length(unique(expanded_samples$name)), "\n")
  cat("   Expansion time:", round(expand_time, 2), "seconds\n")
  
  # Compare dimensions
  cat("\n4. Comparing dimensions...\n")
  orig_vars <- sort(unique(original_samples$name))
  exp_vars <- sort(unique(expanded_samples$name))
  
  missing_vars <- setdiff(orig_vars, exp_vars)
  extra_vars <- setdiff(exp_vars, orig_vars)
  
  cat("   Original variables:", length(orig_vars), "\n")
  cat("   Expanded variables:", length(exp_vars), "\n")
  cat("   Missing variables:", length(missing_vars), "\n")
  cat("   Extra variables:", length(extra_vars), "\n")
  
  if (length(missing_vars) > 0) {
    cat("   Missing:", paste(missing_vars, collapse = ", "), "\n")
  }
  if (length(extra_vars) > 0) {
    cat("   Extra:", paste(extra_vars, collapse = ", "), "\n")
  }
  
  # Compare values for overlapping variables
  cat("\n5. Comparing values...\n")
  common_vars <- intersect(orig_vars, exp_vars)
  
  if (length(common_vars) == 0) {
    cat("   ERROR: No common variables to compare!\n")
    return(list(success = FALSE, error = "No common variables"))
  }
  
  # Test a few iterations for detailed comparison
  test_iters <- head(unique(original_samples$iter), 3)
  value_mismatches <- 0
  total_comparisons <- 0
  
  for (var in head(common_vars, 10)) { # Test first 10 variables
    for (iter_val in test_iters) {
      
      orig_subset <- original_samples[name == var & iter == iter_val]
      exp_subset <- expanded_samples[name == var & iter == iter_val]
      
      if (nrow(orig_subset) != nrow(exp_subset)) {
        cat("   Row count mismatch for", var, "iter", iter_val, 
            ":", nrow(orig_subset), "vs", nrow(exp_subset), "\n")
        value_mismatches <- value_mismatches + 1
        next
      }
      
      if (nrow(orig_subset) > 0) {
        # Sort both by row index for comparison
        orig_sorted <- orig_subset[order(row)]
        exp_sorted <- exp_subset[order(row)]
        
        # Compare values
        value_diffs <- abs(orig_sorted$value - exp_sorted$value)
        max_diff <- max(value_diffs, na.rm = TRUE)
        total_comparisons <- total_comparisons + length(value_diffs)
        
        if (max_diff > tolerance) {
          cat("   Value mismatch for", var, "iter", iter_val, 
              ": max diff =", scientific(max_diff), "\n")
          value_mismatches <- value_mismatches + 1
        }
      }
    }
  }
  
  cat("   Total value comparisons:", format(total_comparisons, big.mark = ","), "\n")
  cat("   Value mismatches:", value_mismatches, "\n")
  
  # Summary statistics
  cat("\n6. Summary statistics...\n")
  compression_ratio <- nrow(condensed_samples) / nrow(original_samples)
  space_savings <- 1 - compression_ratio
  
  cat("   Compression ratio:", round(compression_ratio * 100, 1), "%\n")
  cat("   Space savings:", round(space_savings * 100, 1), "%\n")
  cat("   Total processing time:", round(compress_time + expand_time, 2), "seconds\n")
  
  # Determine success
  success <- (
    length(missing_vars) == 0 &&
    length(extra_vars) == 0 &&
    value_mismatches == 0 &&
    nrow(expanded_samples) == nrow(original_samples)
  )
  
  cat("\n=== TEST RESULT:", ifelse(success, "PASSED", "FAILED"), "===\n")
  
  # Clean up test file
  if (file.exists(compact_file)) {
    file.remove(compact_file)
    cat("Cleaned up hmc_samples_compact.rda\n")
  }
  
  return(list(
    success = success,
    original_rows = nrow(original_samples),
    condensed_rows = nrow(condensed_samples),
    expanded_rows = nrow(expanded_samples),
    compression_ratio = compression_ratio,
    space_savings = space_savings,
    compress_time = compress_time,
    expand_time = expand_time,
    missing_vars = missing_vars,
    extra_vars = extra_vars,
    value_mismatches = value_mismatches,
    total_comparisons = total_comparisons
  ))
}

#' Quick test on a single iteration
#' @param model_dir Directory containing model files
#' @param test_iter Which iteration to test (default 1)
quick_test <- function(model_dir, test_iter = 1) {
  
  cat("Quick test on iteration", test_iter, "\n")
  
  # Compress and expand
  condensed <- condense_hmc_samples(model_dir)
  expanded <- expand_hmc_samples(model_dir)
  
  # Load original for comparison
  original <- fread(file.path(model_dir, "hmc_samples.csv"))
  
  # Compare one iteration
  orig_iter <- original[iter == test_iter]
  exp_iter <- expanded[iter == test_iter]
  
  cat("Original rows for iter", test_iter, ":", nrow(orig_iter), "\n")
  cat("Expanded rows for iter", test_iter, ":", nrow(exp_iter), "\n")
  
  # Check a few specific variables
  test_vars <- c("logK", "r", "x", "removals")
  for (var in test_vars) {
    orig_var <- orig_iter[name == var]
    exp_var <- exp_iter[name == var]
    
    if (nrow(orig_var) > 0 && nrow(exp_var) > 0) {
      orig_val <- orig_var$value[1]
      exp_val <- exp_var$value[1]
      diff <- abs(orig_val - exp_val)
      cat(var, ": orig =", round(orig_val, 6), ", exp =", round(exp_val, 6), 
          ", diff =", scientific(diff), "\n")
    } else {
      cat(var, ": missing in original or expanded\n")
    }
  }
  
  # Clean up
  file.remove(file.path(model_dir, "hmc_samples_compact.rda"))
}

# Example usage:
test_results <- test_hmc_workflow("data/output/model_runs/0007-au-c0.2-e0.3-s5_0")
quick_test("data/output/model_runs/0007-au-c0.2-e0.3-s5_0")
