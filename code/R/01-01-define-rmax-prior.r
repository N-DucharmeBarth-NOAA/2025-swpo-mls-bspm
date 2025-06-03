# Nicholas Ducharme-Barth
# 2025/05/29
# Define priors for BSPM - Parallel processing version
# rmax, K, dep_1986

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(parallel)

#_____________________________________________________________________________________________________________________________
# define paths
proj_dir = this.path::this.proj()
dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# Setup parallel processing
# Detect number of cores and use all but one (or specify manually)
n_cores = parallel::detectCores() - 1  # Leave one core free
cat("Setting up parallel processing with", n_cores, "cores...\n")

# Create cluster
cl = makeCluster(n_cores)

# Export necessary objects to workers
clusterEvalQ(cl, {
    library(data.table)
    library(magrittr)
})

# Export all functions and objects to workers
clusterExport(cl, ls(envir = globalenv()))

# Source helper functions on each worker
clusterEvalQ(cl, {
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)
})

#_____________________________________________________________________________________________________________________________
# Generate parameter sets
n_samples = 500000
cat("Generating", n_samples, "parameter sets...\n")
param_dt = generate_random_params(n_samples = n_samples, seed = 123)

#_____________________________________________________________________________________________________________________________
# Define processing function for a single parameter set
process_parameter_set = function(i, param_dt) {
    
    # Extract parameters for current iteration
    params = param_dt[i]
    
    # Initialize result list
    result = list(
        sample_id = NA_integer_,
        rmax = NA_real_,
        epr_unfished_rmax = NA_real_,
        alpha = NA_real_,
        generation_time = NA_real_,
        inflection_point = NA_real_,
        amat50 = NA_integer_
    )
    
    tryCatch({
        # 1. Calculate Rmax
        rmax_result = sim_rmax(
            id = params$id,
            max_age = params$max_age,
            M_ref = params$M_ref,
            L1 = params$L1,
            L2 = params$L2,
            vbk = params$vbk,
            age1 = params$age1,
            age2 = params$age2,
            cv_len = params$cv_len,
            maturity_a = params$maturity_a,
            l50 = params$l50,
            weight_a = params$weight_a,
            weight_b = params$weight_b,
            sex_ratio = params$sex_ratio,
            reproductive_cycle = params$reproductive_cycle,
            h = params$h
        )
        
        # Store Rmax results
        result$sample_id = rmax_result$id
        result$rmax = rmax_result$rmax
        result$epr_unfished_rmax = rmax_result$epr_unfished
        result$alpha = rmax_result$alpha
        result$generation_time = rmax_result$generation_time
        result$inflection_point = rmax_result$inflection_point
        result$amat50 = rmax_result$amat50
        
        
    }, error = function(e) {
        result$error_message <<- as.character(e$message)
    })
    
    return(result)
}

#_____________________________________________________________________________________________________________________________
# Process parameter sets in parallel
cat("Processing", nrow(param_dt), "parameter sets in parallel...\n")
start_time = Sys.time()

results_list = parLapply(cl, 1:nrow(param_dt), process_parameter_set, param_dt = param_dt)

end_time = Sys.time()
processing_time = end_time - start_time
cat("Parallel processing completed in:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")

# Stop the cluster
stopCluster(cl)

#_____________________________________________________________________________________________________________________________
# Combine results back with original parameters
cat("Combining results...\n")

# Convert results list to data.table
results_dt = rbindlist(results_list, fill = TRUE)

# Merge with original parameter data
param_dt[, sample_id := 1:nrow(param_dt)]  # Ensure sample_id matches
final_dt = merge(param_dt, results_dt, by = "sample_id", all.x = TRUE)

#_____________________________________________________________________________________________________________________________
# Summary and diagnostics
cat("\nProcessing complete!\n")
cat("Summary of results:\n")

# Summary statistics for key outputs
summary_cols = c("rmax", "h", "generation_time")
for(col in summary_cols) {
    if(col %in% names(final_dt)) {
        cat("\n", col, ":\n")
        print(summary(final_dt[[col]], na.rm = TRUE))
    }
}

# Check for failed runs
failed_runs = final_dt[is.na(rmax), .N]
cat("\nFailed runs:", failed_runs, "out of", nrow(final_dt), "\n")

#_____________________________________________________________________________________________________________________________
# Save results
output_file = file.path(proj_dir,"data", "output", "bspm_parameter_priors_parallel.csv")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
fwrite(final_dt, output_file)
cat("Results saved to:", output_file, "\n")

# Also save as RDS for faster loading in R
rds_file = file.path(proj_dir,"data", "output", "bspm_parameter_priors_parallel.rds")
saveRDS(final_dt, rds_file)
cat("Results also saved as RDS to:", rds_file, "\n")

cat("Total processing time:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")

#_____________________________________________________________________________________________________________________________
# Filter data for successful runs (rmax > 0) and calculate shape parameter
cat("\nFiltering data for successful runs (rmax > 0)...\n")

# Read in the generated data (in case it needs to be reloaded)
# Uncomment the next two lines if you need to reload the data
# final_dt = fread(output_file)
# # or: final_dt = readRDS(rds_file)

# Filter for successful runs where rmax > 0 & rmax < 1.5
# These are plausible values based on the upper 95% of teleost Rmax values from
# Gravel et al. 2014 & Hutchings et al. 2012
rmax_filter = final_dt[!is.na(rmax) & rmax > 0 & rmax<1.5]

cat("Original dataset:", nrow(final_dt), "rows\n")
cat("Filtered dataset (rmax > 0):", nrow(rmax_filter), "rows\n")
cat("Proportion of successful runs:", round(nrow(rmax_filter)/nrow(final_dt), 3), "\n")

#_____________________________________________________________________________________________________________________________
# Calculate shape parameter based on inflection point

# Define helper functions for shape parameter calculation
tmp_fn = function(par, target){
    n = exp(par)
    phi = (1/n)^(1/(n-1))
    return((target - phi)^2)
}

wrap_fn = function(input_target){
    return(exp(optim(log(2), tmp_fn, target = input_target, method = "Brent", 
                     lower = -10, upper = 10)$par))
}

cat("Calculating shape parameters...\n")
start_time_shape = Sys.time()

# Calculate shape parameter for each inflection point
rmax_filter[, shape := sapply(inflection_point, wrap_fn)]

end_time_shape = Sys.time()
shape_time = end_time_shape - start_time_shape
cat("Shape parameter calculation completed in:", round(as.numeric(shape_time, units = "secs"), 2), "seconds\n")

#_____________________________________________________________________________________________________________________________
# Summary statistics for filtered dataset
cat("\nSummary statistics for filtered dataset:\n")

# Summary statistics for key outputs including shape
summary_cols_filtered = c("rmax", "h", "generation_time", "inflection_point", "shape")

for(col in summary_cols_filtered) {
    if(col %in% names(rmax_filter)) {
        cat("\n", col, ":\n")
        print(summary(rmax_filter[[col]], na.rm = TRUE))
    }
}

#_____________________________________________________________________________________________________________________________
# Save filtered results
output_file_filtered = file.path(proj_dir, "data", "output", "bspm_parameter_priors_filtered.csv")
fwrite(rmax_filter, output_file_filtered)
cat("\nFiltered results saved to:", output_file_filtered, "\n")

# Also save as RDS for faster loading in R
rds_file_filtered = file.path(proj_dir, "data", "output", "bspm_parameter_priors_filtered.rds")
saveRDS(rmax_filter, rds_file_filtered)
cat("Filtered results also saved as RDS to:", rds_file_filtered, "\n")

#_____________________________________________________________________________________________________________________________
# Final summary
cat("\n" , rep("=", 60), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
cat("Total original samples:", nrow(final_dt), "\n")
cat("Successful runs (rmax > 0):", nrow(rmax_filter), "\n")
cat("Success rate:", round(nrow(rmax_filter)/nrow(final_dt) * 100, 1), "%\n")
cat("Total processing time:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")
cat("Shape calculation time:", round(as.numeric(shape_time, units = "secs"), 2), "seconds\n")

# Check for any issues with shape calculation
shape_issues = rmax_filter[is.na(shape) | !is.finite(shape), .N]
if(shape_issues > 0) {
    cat("WARNING: ", shape_issues, " rows had issues with shape parameter calculation\n")
}

cat("Processing complete!\n")
