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
n_samples = 400000
cat("Generating", n_samples, "parameter sets...\n")
param_dt = generate_random_params(n_samples = n_samples, seed = 123)

# Add target average weight (read this in)
nz_wt = fread(file.path(proj_dir,"data","input","nz-all-sport-club-weights.csv"))
        colnames(nz_wt) = c("date","club","species","tagged","weight","boat","locality","year")
        nz_wt = nz_wt %>%
                .[!is.na(weight)&!is.na(year)&tagged=="",.(date,year,weight)] %>%
                .[year<1988] %>%
                .[,dd:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][2]))] %>%
                .[,mm:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][1]))] %>%
                .[,yy:=sapply(date,function(x)as.numeric(strsplit(x,"/")[[1]][3]))] %>%
                .[yy>51&yy<100,yy:=yy+1900] %>%
                .[yy<25,yy:=yy+2000] %>%
                .[,year:=yy] %>%
                .[,month:=c(2,2,2,5,5,5,8,8,8,11,11,11)[mm]] %>%
                .[,.(year,month,weight)] %>%
                .[,ts:=paste0(year,"-",month)]

param_dt[, target_average_weight := mean(nz_wt[year %in% 1985:1987]$weight)]

#_____________________________________________________________________________________________________________________________
# Define processing function for a single parameter set
process_parameter_set = function(i, param_dt) {
    
    # Extract parameters for current iteration
    params = param_dt[i]
    
    # Initialize result list
    result = list(
        sample_id = params$sample_id,
        rmax = NA_real_,
        epr_unfished_rmax = NA_real_,
        alpha = NA_real_,
        h = NA_real_,
        generation_time = NA_real_,
        inflection_point = NA_real_,
        amat50 = NA_integer_,
        avg_weight_unfished = NA_real_,
        F_est = NA_real_,
        spr = NA_real_,
        dep = NA_real_,
        dep_sb = NA_real_,
        epr_fished = NA_real_,
        epr_unfished_spr = NA_real_,
        convergence = NA_integer_,
        error_message = NA_character_
    )
    
    tryCatch({
        # 1. Calculate Rmax
        rmax_result = sim_rmax(
            id = params$sample_id,
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
            reproductive_cycle = params$reproductive_cycle
        )
        
        # Store Rmax results
        result$rmax = rmax_result$rmax
        result$epr_unfished_rmax = rmax_result$epr_unfished
        result$alpha = rmax_result$alpha
        result$h = rmax_result$h
        result$generation_time = rmax_result$generation_time
        result$inflection_point = rmax_result$inflection_point
        result$amat50 = rmax_result$amat50
        
        # 2. Calculate average weight in unfished conditions
        result$avg_weight_unfished = calc_avg_capture_weight_at_age_unfished(
            max_age = params$max_age,
            L1 = params$L1,
            L2 = params$L2,
            vbk = params$vbk,
            age1 = params$age1,
            age2 = params$age2,
            cv_len = params$cv_len,
            M_ref = params$M_ref,
            selex_slope = params$selex_slope,
            selex_l50 = params$selex_l50,
            maturity_a = params$maturity_a,
            l50 = params$l50,
            weight_a = params$weight_a,
            weight_b = params$weight_b,
            selexNZ_l50 = params$selexNZ_l50,
            selexNZ_slope = params$selexNZ_slope
        )
        
        # 3. Optimize F to achieve target average weight
        optim_result = optim(
            par = params$M_ref,
            fn = objective_function,
            target_avg_weight = params$target_average_weight,
            max_age = params$max_age,
            L1 = params$L1,
            L2 = params$L2,
            vbk = params$vbk,
            age1 = params$age1,
            age2 = params$age2,
            cv_len = params$cv_len,
            M_ref = params$M_ref,
            selex_slope = params$selex_slope,
            selex_l50 = params$selex_l50,
            maturity_a = params$maturity_a,
            l50 = params$l50,
            weight_a = params$weight_a,
            weight_b = params$weight_b,
            selexNZ_l50 = params$selexNZ_l50,
            selexNZ_slope = params$selexNZ_slope,
            method = "Brent",
            lower = 0,
            upper = 5
        )
        
        result$F_est = optim_result$par
        result$convergence = optim_result$convergence
        
        # 4. Calculate SPR with estimated F
        spr_result = calc_spr(
            id = params$sample_id,
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
            F = result$F_est,
            selex_l50 = params$selex_l50,
            selex_slope = params$selex_slope
        )
        
        # Store SPR results
        result$spr = spr_result$spr
        result$dep = spr_result$dep
        result$dep_sb = spr_result$dep_sb
        result$epr_fished = spr_result$epr_fished
        result$epr_unfished_spr = spr_result$epr_unfished
        
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

# Check convergence
convergence_summary = final_dt[!is.na(convergence), .N, by = convergence]
cat("Convergence summary:\n")
print(convergence_summary)

# Check for errors
error_count = final_dt[!is.na(error_message), .N]
cat("Runs with errors:", error_count, "out of", nrow(final_dt), "\n")

if(error_count > 0) {
    cat("Sample error messages:\n")
    print(final_dt[!is.na(error_message), .(sample_id, error_message)][1:min(5, .N)])
}

# Summary statistics for key outputs
summary_cols = c("rmax", "h", "generation_time", "avg_weight_unfished", "F_est", "spr")
for(col in summary_cols) {
    if(col %in% names(final_dt)) {
        cat("\n", col, ":\n")
        print(summary(final_dt[[col]], na.rm = TRUE))
    }
}

# Check for failed runs
failed_runs = final_dt[is.na(rmax) | is.na(F_est) | is.na(spr), .N]
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

# Filter for successful runs where rmax > 0
rmax_filter = final_dt[!is.na(rmax) & rmax > 0]

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
summary_cols_filtered = c("rmax", "h", "generation_time", "avg_weight_unfished", 
                         "F_est", "spr", "dep", "inflection_point", "shape")

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
