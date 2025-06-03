# Nicholas Ducharme-Barth
# 2025/06/03
# Calculate fishing mortality rates and relative depletion for multiple biological parameter sets
# Based on Gedamke and Hoenig 2006 - Parallel processing version
# Uses change in mean size to estimate Z1 and Z2, then calculates SPR-based depletion

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(parallel)

#_____________________________________________________________________________________________________________________________
# define helper functions
transitional_mean_weight = function(W_inf, Wc, vbk, Z1, Z2, d) {
    # Equation 6 - modified to be in terms of weight
    # d time in years
    
    # Calculate each component separately for clarity
    term1 = W_inf
    
    numerator_fraction = Z1 * Z2 * (W_inf - Wc) * 
                        (Z1 + vbk + (Z2 - Z1) * exp(-(Z2 + vbk) * d))
    
    denominator_fraction = (Z1 + vbk) * (Z2 + vbk) * 
                            (Z1 + (Z2 - Z1) * exp(-Z2 * d))
    
    mean_weight = term1 - (numerator_fraction / denominator_fraction)
    
    return(mean_weight)
}

negative_log_likelihood = function(par, data, bio_params) {
    # Equation 8 implementation
    # define leading parameters
    Z1 = par[1]
    Z2 = par[2]

    # extract pars from bio_params
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    Lc = bio_params$Lc
    
    # define W_inf and Wc
    age_vector = seq(from=0, to=bio_params$max_age*2, length.out=50)
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    W_inf = weight_a*max(length_at_age)^weight_b
    Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

    # process data
    proc_dat = data[weight>Wc&year>=1952] %>%
               .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
               .[,d:=year-1952]

    # calculate d_vec
    d_vec = proc_dat$d

    # Calculate predicted lengths (using your transitional equation)
    predicted = transitional_mean_weight(W_inf, Wc, vbk, Z1, Z2, d_vec)
    
    log_likelihood = sum(dnorm(proc_dat$mean_wt, 
                         mean = predicted, 
                         sd = proc_dat$sd_wt/sqrt(proc_dat$N), 
                         log = TRUE))
    
    return(-log_likelihood)  # Return negative for minimization
}

rel_depletion = function(result, bio_params) {
    # Calculate relative depletion from Z1 to Z2
    # define parameters
    Z1 = result$par[1]
    Z2 = result$par[2]
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    cv_len = bio_params$cv_len
    maturity_a = bio_params$maturity_a
    l50 = bio_params$l50
    sex_ratio = bio_params$sex_ratio
    reproductive_cycle = bio_params$reproductive_cycle
    
    # age 
    age_vector = 1:max_age
    
    # length
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    len_lower = seq(from=0, to=ceiling(max(length_at_age*(1+cv_len))), by=1)
    len_upper = len_lower + 1
    length_vec = len_lower + 0.5

    # survival
    survival_at_age_Z1 = rep(NA, max_age)
    survival_at_age_Z1[1] = 1
    survival_at_age_Z2 = rep(NA, max_age)
    survival_at_age_Z2[1] = 1
    for(i in 2:length(age_vector)) {
        survival_at_age_Z1[i] = survival_at_age_Z1[i-1] * exp(-Z1)
        survival_at_age_Z2[i] = survival_at_age_Z2[i-1] * exp(-Z2)
    }
    survival_at_age_Z1[max_age] = survival_at_age_Z1[max_age-1] / (1 - exp(-Z1))
    survival_at_age_Z2[max_age] = survival_at_age_Z2[max_age-1] / (1 - exp(-Z2))

    # maturity
    maturity_b = -maturity_a/l50
    maturity_at_length = (exp(maturity_a+maturity_b*length_vec)) / (1+exp(maturity_a+maturity_b*length_vec))

    # calc maturity at age using PLA
    pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
    maturity_at_age = as.vector(matrix(maturity_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
    maturity_at_age = maturity_at_age/max(maturity_at_age)
    
    # weight
    weight_at_length = weight_a * length_vec ^ weight_b
    weight_at_age = as.vector(matrix(weight_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)

    # reproductive output per year (sex-ratio & reproductive cycle)
    reproduction_at_age = (maturity_at_age * weight_at_age * sex_ratio)/reproductive_cycle

    # spr calc
    # lifetime average eggs per recruit in fished and unfished conditions
    epr_Z1 = sum(reproduction_at_age*survival_at_age_Z1)
    epr_Z2 = sum(reproduction_at_age*survival_at_age_Z2)
    rel_dep_n = sum(survival_at_age_Z2)/sum(survival_at_age_Z1)
    rel_dep_ssb = epr_Z2/epr_Z1
    
    return(data.table(sample_id=bio_params$sample_id, rel_dep_n = rel_dep_n, rel_dep_ssb = rel_dep_ssb))
}

mean_wt_resid = function(result, bio_params, data) {
    # Calculate mean weight residuals
    # define parameters
    Z1 = result$par[1]
    Z2 = result$par[2]
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    Lc = bio_params$Lc

    # define W_inf and Wc
    age_vector = seq(from=0, to=bio_params$max_age*2, length.out=50)
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    W_inf = weight_a*max(length_at_age)^weight_b
    Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

    # process data
    proc_dat = data[weight>Wc&year>=1952] %>%
            .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
            .[,d:=year-1952]

    # calculate d_vec
    d_vec = proc_dat$d

    # Calculate predicted lengths (using your transitional equation)
    predicted = transitional_mean_weight(W_inf, Wc, vbk, Z1, Z2, d_vec)
    
    proc_dat$pred_wt = predicted
    proc_dat$resid = proc_dat$mean_wt - proc_dat$pred_wt
    proc_dat$sample_id = bio_params$sample_id

    return(proc_dat[,.(sample_id, year, d, N, mean_wt, sd_wt, pred_wt, resid)])
}

#_____________________________________________________________________________________________________________________________
# define paths
proj_dir = this.path::this.proj()
dir_helper_fns = file.path(proj_dir, "code", "R", "helper-fns")
model_run_dir = file.path(proj_dir, "data", "output")
dir.create(model_run_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
sapply(file.path(dir_helper_fns, (list.files(dir_helper_fns))), source)

#_____________________________________________________________________________________________________________________________
# Setup parallel processing
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
    sapply(file.path(dir_helper_fns, (list.files(dir_helper_fns))), source)
})

#_____________________________________________________________________________________________________________________________
# Load data
cat("Loading biological parameters and NZ weight data...\n")

# Read in biological parameter results
bio_params_dt = fread(file.path(proj_dir, "data", "output", "bspm_parameter_priors_filtered.csv"))

#_____________________________________________________________________________________________________________________________
# Subsampling option
# Set to NULL to process all data, or specify number of samples
n_subsample = 10000  # Change this to desired number (e.g., 1000) or leave as NULL for all data
subsample_seed = 456  # Set seed for reproducible subsampling

if(!is.null(n_subsample) && n_subsample < nrow(bio_params_dt)) {
    cat("Subsampling", n_subsample, "parameter sets from", nrow(bio_params_dt), "total sets...\n")
    set.seed(subsample_seed)
    subsample_indices = sample(1:nrow(bio_params_dt), size = n_subsample, replace = FALSE)
    bio_params_dt = bio_params_dt[subsample_indices]
    cat("Subsampling complete. Processing", nrow(bio_params_dt), "parameter sets.\n")
} else if(!is.null(n_subsample) && n_subsample >= nrow(bio_params_dt)) {
    cat("Requested subsample size (", n_subsample, ") >= total available (", nrow(bio_params_dt), "). Processing all data.\n")
} else {
    cat("No subsampling requested. Processing all", nrow(bio_params_dt), "parameter sets.\n")
}

# Load NZ recreational data
nz_wt = fread(file.path(proj_dir, "data", "input", "nz-all-sport-club-weights.csv"))
colnames(nz_wt) = c("date", "club", "species", "tagged", "weight", "boat", "locality", "year")
nz_wt = nz_wt %>%
            .[!is.na(weight)&!is.na(year)&tagged=="", .(date, year, weight)] %>%
            .[year<1988] %>%
            .[, dd:=sapply(date, function(x)as.numeric(strsplit(x, "/")[[1]][2]))] %>%
            .[, mm:=sapply(date, function(x)as.numeric(strsplit(x, "/")[[1]][1]))] %>%
            .[, yy:=sapply(date, function(x)as.numeric(strsplit(x, "/")[[1]][3]))] %>%
            .[yy>51&yy<100, yy:=yy+1900] %>%
            .[yy<25, yy:=yy+2000] %>%
            .[, year:=yy] %>%
            .[, month:=c(2,2,2,5,5,5,8,8,8,11,11,11)[mm]] %>%
            .[, .(year, month, weight)] %>%
            .[, ts:=paste0(year, "-", month)] %>%
            na.omit(.)

cat("Loaded", nrow(bio_params_dt), "biological parameter sets (after subsampling if applied)\n")
cat("Loaded", nrow(nz_wt), "weight observations\n")

# Export data to workers
clusterExport(cl, c("nz_wt", "bio_params_dt"))

#_____________________________________________________________________________________________________________________________
# Define processing function for a single parameter set
process_depletion_set = function(i, bio_params_dt, nz_wt) {
    
    # Extract parameters for current iteration
    params = bio_params_dt[i]
    
    # Initialize result list
    result = list(
        sample_id = params$sample_id,
        Z1 = NA_real_,
        Z2 = NA_real_,
        convergence = NA_integer_,
        log_likelihood = NA_real_,
        rel_dep_n = NA_real_,
        rel_dep_ssb = NA_real_,
        error_message = NA_character_
    )
    
    tryCatch({
        # Fit Z1 and Z2 using optimization
        opt_result = nlminb(start = c(params$M_ref + 0.1, params$M_ref + 0.2),
                           objective = negative_log_likelihood,
                           data = nz_wt,
                           bio_params = params,
                           lower = rep(params$M_ref, 2),
                           upper = rep(3, 2))
        
        # Store optimization results
        result$Z1 = opt_result$par[1]
        result$Z2 = opt_result$par[2]
        result$convergence = opt_result$convergence
        result$log_likelihood = opt_result$objective
        
        # Calculate relative depletion if optimization converged
        if(opt_result$convergence == 0) {
            depletion_result = rel_depletion(opt_result, params)
            result$rel_dep_n = depletion_result$rel_dep_n
            result$rel_dep_ssb = depletion_result$rel_dep_ssb
        }
        
    }, error = function(e) {
        result$error_message <<- as.character(e$message)
    })
    
    return(result)
}

# Define separate function to calculate residuals for successful runs
calculate_residuals_parallel = function(successful_dt, nz_wt) {
    
    cat("Calculating residuals for", nrow(successful_dt), "successful parameter sets...\n")
    
    # Setup parallel processing for residuals
    n_cores_resid = min(parallel::detectCores() - 1, nrow(successful_dt))
    cl_resid = makeCluster(n_cores_resid)
    
    # Export necessary objects to workers
    clusterEvalQ(cl_resid, {
        library(data.table)
        library(magrittr)
    })
    
    # Export functions and data to workers
    clusterExport(cl_resid, c("mean_wt_resid", "transitional_mean_weight", "nz_wt"))
    
    # Function to calculate residuals for one parameter set
    calc_residuals_single = function(i, successful_dt, nz_wt) {
        params = successful_dt[i]
        
        # Recreate opt_result structure from stored values
        opt_result = list(par = c(params$Z1, params$Z2))
        
        tryCatch({
            residuals = mean_wt_resid(opt_result, params, nz_wt)
            return(residuals)
        }, error = function(e) {
            return(data.table(sample_id = params$sample_id, error = as.character(e$message)))
        })
    }
    
    # Calculate residuals in parallel
    residuals_list = parLapply(cl_resid, 1:nrow(successful_dt), calc_residuals_single, 
                              successful_dt = successful_dt, nz_wt = nz_wt)
    
    # Stop the cluster
    stopCluster(cl_resid)
    
    # Combine residuals
    residuals_dt = rbindlist(residuals_list, fill = TRUE)
    
    return(residuals_dt)
}

#_____________________________________________________________________________________________________________________________
# Process parameter sets in parallel
cat("Processing", nrow(bio_params_dt), "parameter sets in parallel...\n")
start_time = Sys.time()

results_list = parLapply(cl, 1:nrow(bio_params_dt), process_depletion_set, 
                        bio_params_dt = bio_params_dt, nz_wt = nz_wt)

end_time = Sys.time()
processing_time = end_time - start_time
cat("Parallel processing completed in:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")

# Stop the cluster
stopCluster(cl)

#_____________________________________________________________________________________________________________________________
# Combine results
cat("Combining results...\n")

# Convert results list to data.table
results_dt = rbindlist(results_list, fill = TRUE)

# Merge with original parameter data
final_dt = merge(bio_params_dt, results_dt, by = "sample_id", all.x = TRUE)

#_____________________________________________________________________________________________________________________________
# Summary and diagnostics
cat("\nProcessing complete!\n")
cat("Summary of results:\n")

# Check convergence
converged_runs = final_dt[convergence == 0 & !is.na(Z1), .N]
failed_runs = final_dt[is.na(Z1) | convergence != 0, .N]

cat("Converged runs:", converged_runs, "out of", nrow(final_dt), "\n")
cat("Failed/non-converged runs:", failed_runs, "\n")
cat("Success rate:", round(converged_runs/nrow(final_dt) * 100, 1), "%\n")

# Summary statistics for key outputs
if(converged_runs > 0) {
    converged_dt = final_dt[convergence == 0 & !is.na(Z1)]
    
    summary_cols = c("Z1", "Z2", "rel_dep_n", "rel_dep_ssb")
    for(col in summary_cols) {
        if(col %in% names(converged_dt)) {
            cat("\n", col, ":\n")
            print(summary(converged_dt[[col]], na.rm = TRUE))
        }
    }
}

#_____________________________________________________________________________________________________________________________
# Save results
output_file = file.path(proj_dir, "data", "output", 
                       if(!is.null(n_subsample)) paste0("depletion_results_parallel_n", n_subsample, ".csv") 
                       else "depletion_results_parallel.csv")
fwrite(final_dt, output_file)
cat("\nResults saved to:", output_file, "\n")

# Also save as RDS for faster loading in R
rds_file = file.path(proj_dir, "data", "output", 
                    if(!is.null(n_subsample)) paste0("depletion_results_parallel_n", n_subsample, ".rds")
                    else "depletion_results_parallel.rds")
saveRDS(final_dt, rds_file)
cat("Results also saved as RDS to:", rds_file, "\n")

#_____________________________________________________________________________________________________________________________
# Create filtered dataset with only successful runs
if(converged_runs > 0) {
    successful_dt = final_dt[convergence == 0 & !is.na(Z1) & !is.na(rel_dep_ssb)]
    
    output_file_success = file.path(proj_dir, "data", "output", 
                                   if(!is.null(n_subsample)) paste0("depletion_results_successful_n", n_subsample, ".csv")
                                   else "depletion_results_successful.csv")
    fwrite(successful_dt, output_file_success)
    cat("Successful results saved to:", output_file_success, "\n")
    
    rds_file_success = file.path(proj_dir, "data", "output", 
                                if(!is.null(n_subsample)) paste0("depletion_results_successful_n", n_subsample, ".rds")
                                else "depletion_results_successful.rds")
    saveRDS(successful_dt, rds_file_success)
    cat("Successful results also saved as RDS to:", rds_file_success, "\n")
    
    # Calculate residuals for successful runs using separate parallel processing
    residuals_dt = calculate_residuals_parallel(successful_dt, nz_wt)
    
    # Save residuals data.table
    if(nrow(residuals_dt) > 0) {
        residuals_output_file = file.path(proj_dir, "data", "output", 
                                         if(!is.null(n_subsample)) paste0("mean_wt_residuals_successful_n", n_subsample, ".csv")
                                         else "mean_wt_residuals_successful.csv")
        fwrite(residuals_dt, residuals_output_file)
        cat("Residuals saved to:", residuals_output_file, "\n")
        
        residuals_rds_file = file.path(proj_dir, "data", "output", 
                                      if(!is.null(n_subsample)) paste0("mean_wt_residuals_successful_n", n_subsample, ".rds")
                                      else "mean_wt_residuals_successful.rds")
        saveRDS(residuals_dt, residuals_rds_file)
        cat("Residuals also saved as RDS to:", residuals_rds_file, "\n")
    }
}

#_____________________________________________________________________________________________________________________________
# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
if(!is.null(n_subsample)) {
    cat("Subsampled to:", n_subsample, "parameter sets\n")
}
cat("Total parameter sets processed:", nrow(final_dt), "\n")
cat("Successful optimizations:", converged_runs, "\n")
cat("Success rate:", round(converged_runs/nrow(final_dt) * 100, 1), "%\n")
cat("Total processing time:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")

if(converged_runs > 0 & exists("residuals_dt") && nrow(residuals_dt) > 0) {
    cat("Residuals data points:", nrow(residuals_dt), "\n")
    cat("Unique years in residuals:", length(unique(residuals_dt$year)), "\n")
}

if(failed_runs > 0) {
    error_summary = final_dt[!is.na(error_message), .N, by = error_message]
    if(nrow(error_summary) > 0) {
        cat("\nError summary:\n")
        print(error_summary)
    }
}

cat("Processing complete!\n")
