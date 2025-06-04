# Nicholas Ducharme-Barth
# 2025/06/03
# Calculate fishing mortality rates and relative depletion for multiple biological parameter sets
# Based on Gedamke and Hoenig 2006 - Parallel processing version
# Uses change in mean size to estimate Z1 and Z2, then calculates SPR-based depletion
# Enhanced with additional filtering and visualization

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(parallel)
library(ggplot2)

#_____________________________________________________________________________________________________________________________
# define paths
proj_dir = this.path::this.proj()
dir_helper_fns = file.path(proj_dir, "code", "R", "helper-fns")
model_run_dir = file.path(proj_dir, "data", "output")
plot_dir = file.path(proj_dir, "plots", "depletion")
dir.create(model_run_dir, recursive = TRUE)
dir.create(plot_dir, recursive = TRUE)

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
n_subsample = NULL  # Change this to desired number (e.g., 1000) or leave as NULL for all data
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
# Apply additional filter: rel_dep_ssb < 1
cat("\nApplying additional filter: rel_dep_ssb < 1...\n")

if(converged_runs > 0) {
    # Initial successful runs
    successful_dt = final_dt[convergence == 0 & !is.na(Z1) & !is.na(rel_dep_ssb)]
    
    # Apply depletion filter
    filtered_dt = successful_dt[rel_dep_ssb < 1]
    
    cat("Runs before depletion filter:", nrow(successful_dt), "\n")
    cat("Runs after depletion filter (rel_dep_ssb < 1):", nrow(filtered_dt), "\n")
    cat("Filter retention rate:", round(nrow(filtered_dt)/nrow(successful_dt) * 100, 1), "%\n")
    
    if(nrow(filtered_dt) > 0) {
        cat("\nSummary of filtered rel_dep_ssb:\n")
        print(summary(filtered_dt$rel_dep_ssb))
    }
} else {
    filtered_dt = data.table()
    cat("No converged runs to filter.\n")
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
# Save filtered dataset and calculate residuals
if(nrow(filtered_dt) > 0) {
    output_file_filtered = file.path(proj_dir, "data", "output", 
                                   if(!is.null(n_subsample)) paste0("depletion_results_filtered_n", n_subsample, ".csv")
                                   else "depletion_results_filtered.csv")
    fwrite(filtered_dt, output_file_filtered)
    cat("Filtered results saved to:", output_file_filtered, "\n")
    
    rds_file_filtered = file.path(proj_dir, "data", "output", 
                                if(!is.null(n_subsample)) paste0("depletion_results_filtered_n", n_subsample, ".rds")
                                else "depletion_results_filtered.rds")
    saveRDS(filtered_dt, rds_file_filtered)
    cat("Filtered results also saved as RDS to:", rds_file_filtered, "\n")
    
    # Calculate residuals for filtered runs using separate parallel processing
    residuals_dt = calculate_residuals_parallel(filtered_dt, nz_wt)
    
    # Save residuals data.table
    if(nrow(residuals_dt) > 0) {
        residuals_output_file = file.path(proj_dir, "data", "output", 
                                         if(!is.null(n_subsample)) paste0("mean_wt_residuals_filtered_n", n_subsample, ".csv")
                                         else "mean_wt_residuals_filtered.csv")
        fwrite(residuals_dt, residuals_output_file)
        cat("Residuals saved to:", residuals_output_file, "\n")
        
        residuals_rds_file = file.path(proj_dir, "data", "output", 
                                      if(!is.null(n_subsample)) paste0("mean_wt_residuals_filtered_n", n_subsample, ".rds")
                                      else "mean_wt_residuals_filtered.rds")
        saveRDS(residuals_dt, residuals_rds_file)
        cat("Residuals also saved as RDS to:", residuals_rds_file, "\n")
    }
}

#_____________________________________________________________________________________________________________________________
# Calculate observed mean weights by year
if(nrow(residuals_dt) > 0) {
    cat("\nCalculating observed mean weights and creating visualization...\n")
    
    # Calculate observed mean weight by year from the NZ data
    # Use the same processing as in mean_wt_resid function
    obs_mean_wt = nz_wt[year>=1952] %>%
                  .[,.(mean_wt_obs=mean(weight), sd_wt_obs=sd(weight), .N),by=year]
    
    # Merge residuals with filtered results to get rel_dep_ssb values
    residuals_with_depletion = merge(residuals_dt, filtered_dt[,.(sample_id, rel_dep_ssb)], by = "sample_id") %>%
                               .[sample_id %in% sample(unique(filtered_dt$sample_id),300)]
    
    # Create the mean weight plot
    p_meanwt = ggplot() +
        # Spaghetti lines for predicted weights with color based on rel_dep_ssb
        geom_line(data = residuals_with_depletion, 
                 aes(x = year, y = pred_wt, group = sample_id, color = rel_dep_ssb),
                 alpha = 0.4, size = 0.3,linewidth=0.8) +
        # Error bars for observed weights
        geom_errorbar(data = obs_mean_wt,
                     aes(x = year, ymin = mean_wt_obs - sd_wt_obs/sqrt(N), 
                         ymax = mean_wt_obs + sd_wt_obs/sqrt(N)),
                     color = "black", alpha = 0.7, width = 0.5) +
        # Observed mean weights as points
        geom_point(data = obs_mean_wt, 
                  aes(x = year, y = mean_wt_obs),
                  shape = 21, color = "black", fill = "white",size = 2, alpha = 0.9) +
        viridis::scale_color_viridis("Relative Depletion\n(SSB)", 
                                    begin = 0.1, end = 0.9, 
                                    direction = 1, option = "turbo", 
                                    discrete = FALSE) +
        labs(x = "Year", 
             y = "Mean Weight (kg)",
             title = "Observed vs. Predicted Mean Weights",
             subtitle = "Red points: observed data, Colored lines: model predictions") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
              legend.position = "right",
              panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
              strip.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white"))
    
    # Save the plot
    plot_file = file.path(plot_dir, 
                         if(!is.null(n_subsample)) paste0("mean_weight_comparison_n", n_subsample, ".png")
                         else "mean_weight_comparison.png")
    
    ggsave(filename = plot_file, plot = p_meanwt, device = "png",
           width = 12, height = 8, units = "in", dpi = 300)
    cat("Mean weight plot saved to:", plot_file, "\n")
}

#_____________________________________________________________________________________________________________________________
# Calculate lognormal prior for rel_dep_ssb and create distribution plot
if(nrow(filtered_dt) > 0) {
    cat("\nCalculating lognormal prior for rel_dep_ssb...\n")
    
    # Fit lognormal distribution to rel_dep_ssb
    rel_dep_values = filtered_dt$rel_dep_ssb[filtered_dt$rel_dep_ssb > 0]  # Ensure positive values
    
    if(length(rel_dep_values) > 0) {
        # Fit lognormal parameters using maximum likelihood
        rel_dep_fn = function(par){-sum(dlnorm(rel_dep_values, meanlog = par[1], sdlog = par[2], log = TRUE))}
        rel_dep_pars = nlminb(c(log(mean(rel_dep_values)), sd(log(rel_dep_values))), rel_dep_fn)$par
        
        # Save parameters
        rel_dep_output_file = file.path(model_run_dir, 
                                       if(!is.null(n_subsample)) paste0("rel_dep_ssb_pars_n", n_subsample, ".csv")
                                       else "rel_dep_ssb_pars.csv")
        write.csv(rel_dep_pars, file = rel_dep_output_file, row.names = FALSE)
        cat("Rel_dep_ssb parameters saved to:", rel_dep_output_file, "\n")
        
        # Create distribution plot
        plot_file_dist = file.path(plot_dir, 
                                  if(!is.null(n_subsample)) paste0("prior_rel_dep_ssb_n", n_subsample, ".png")
                                  else "prior_rel_dep_ssb.png")
        
        png(filename = plot_file_dist, width = 8, height = 6, units = "in", bg = "white", res = 300)
        
        # Create histogram
        hist(rel_dep_values, freq = FALSE, breaks = 50, 
             xlab = "Relative Depletion (SSB)", 
             main = "Prior Distribution: Relative Depletion (SSB)",
             col = "lightblue", border = "white")
        
        # Add fitted lognormal curve
        plot_x = seq(from = 0.001, to = max(rel_dep_values), length.out = 1000)
        plot_y = dlnorm(plot_x, rel_dep_pars[1], rel_dep_pars[2])
        lines(plot_x, plot_y, col = "red", lwd = 3)
        
        # Add density curve for comparison
        lines(density(rel_dep_values, adjust = 1.2), col = "blue", lwd = 2, lty = 2)
        
        # Add legend
        legend("topright", 
               c("Fitted Lognormal", "Kernel Density", "Histogram"), 
               col = c("red", "blue", "lightblue"), 
               lty = c(1, 2, 1), lwd = c(3, 2, 1), 
               fill = c(NA, NA, "lightblue"),
               border = c(NA, NA, "white"),
               bty = "n")
        
        # Add parameter text
        text(x = max(rel_dep_values) * 0.7, y = max(plot_y) * 0.8,
             labels = paste0("Log-normal parameters:\nMeanlog = ", round(rel_dep_pars[1], 3),
                           "\nSDlog = ", round(rel_dep_pars[2], 3),
                           "\nN = ", length(rel_dep_values)),
             adj = 0, cex = 0.9, 
             bg = "white", box.col = "gray")
        
        dev.off()
        cat("Rel_dep_ssb distribution plot saved to:", plot_file_dist, "\n")
        
        # Print summary
        cat("\nRel_dep_ssb lognormal parameters:\n")
        cat("Meanlog:", round(rel_dep_pars[1], 4), "\n")
        cat("SDlog:", round(rel_dep_pars[2], 4), "\n")
        cat("Mean (arithmetic):", round(exp(rel_dep_pars[1] + rel_dep_pars[2]^2/2), 4), "\n")
        cat("Median:", round(exp(rel_dep_pars[1]), 4), "\n")
    } else {
        cat("No positive rel_dep_ssb values found for lognormal fitting.\n")
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
if(converged_runs > 0) {
    cat("Runs passing depletion filter (rel_dep_ssb < 1):", nrow(filtered_dt), "\n")
    cat("Final success rate:", round(nrow(filtered_dt)/nrow(final_dt) * 100, 1), "%\n")
}
cat("Total processing time:", round(as.numeric(processing_time, units = "mins"), 2), "minutes\n")

if(nrow(filtered_dt) > 0 & exists("residuals_dt") && nrow(residuals_dt) > 0) {
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
