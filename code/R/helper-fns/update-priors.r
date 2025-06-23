# Update priors from prior pushforward check output
# Based on the filtering approach in 03-03-prior-pushforward-q.r
# Creates new priors for logK, log_r, log_shape, atanh_rho, log_sigma_qdev, log_qeff, logsigmap, depletion

# ppc_output_path = "data/output/model_runs/0006-2024cpueEffortQeff_0-ppc"
# original_stan_data = stan.data
# catch_effort_dt = catch_effort_annual
# verbose = TRUE

update_priors_from_ppc = function(ppc_output_path, original_stan_data = NULL, 
                                  catch_effort_dt = NULL, verbose = TRUE) {
    
    # Load required libraries
    library(data.table)
    library(MASS)
    library(mvtnorm)
    
    # Read HMC samples from prior pushforward check
    hmc_samples_path = file.path(ppc_output_path, "hmc_samples.csv")
    if (!file.exists(hmc_samples_path)) {
        stop("hmc_samples.csv not found in: ", ppc_output_path)
    }
    
    hmc_dt = fread(hmc_samples_path) %>% .[,value:=as.numeric(value)]
    
    # Exclude divergent transitions if they exist
    n_divergent = uniqueN(hmc_dt[divergent == 1]$iter)
    if (n_divergent > 0) {
        if (verbose) cat("Excluding", n_divergent, "samples with divergent transitions\n")
        hmc_dt = hmc_dt[divergent == 0]
    }
    
    # Create unique simulation identifier 
    hmc_dt[, seed := paste(run_id, iter, chain, sep = "_")]
    
    # Check what variables are available in the samples
    available_vars = unique(hmc_dt$variable)
    
    # Extract model predictions for filtering using name and row columns
    # Get x (relative population) - should always be available
    if (!"x" %in% unique(hmc_dt$name)) {
        stop("No x variables found in HMC samples")
    }
    x_samples = hmc_dt[name == "x", .(seed, time_idx = row, dep = value)]
    
    # Get F (fishing mortality) - check if available  
    if ("F" %in% unique(hmc_dt$name)) {
        f_samples = hmc_dt[name == "F", .(seed, time_idx = row, F = value)]
        if (verbose) cat("Found F variables for filtering\n")
    } else {
        if (verbose) cat("F not available - using constant F = 0.5 for filtering\n")
        f_samples = x_samples[, .(seed, time_idx, F = 0.5)]  # Conservative placeholder
    }
    
    # Get removals (predicted catch) - should be available for effort models
    if (!"removals" %in% unique(hmc_dt$name)) {
        stop("No removals variables found in HMC samples")
    }
    removals_samples = hmc_dt[name == "removals", .(seed, time_idx = row, predicted_catch = value)]
    
    # Get transformed parameters for filtering and prior updating
    logK_samples = hmc_dt[variable == "logK", .(seed, logK = value)]
    r_samples = hmc_dt[variable == "r", .(seed, r = value)]
    shape_samples = hmc_dt[variable == "shape", .(seed, shape = value)]
    
    # Check that key parameters exist
    if (!"logK" %in% available_vars) stop("logK not found in HMC samples")
    if (!"r" %in% available_vars) stop("r not found in HMC samples")  
    if (!"shape" %in% available_vars) stop("shape not found in HMC samples")
    
    # Merge all simulation data including biological parameters
    sim_dt = merge(x_samples, f_samples, by = c("seed", "time_idx"), all = TRUE) %>%
              merge(., removals_samples, by = c("seed", "time_idx"), all = TRUE) %>%
              merge(., logK_samples, by = "seed", all = TRUE) %>%
              merge(., shape_samples, by = "seed", all = TRUE) 
    
    # Calculate population numbers: N = x * exp(logK)
    sim_dt[, n := dep * exp(logK)]
    
    # Add time variable (assuming model starts at 1952 based on your code)
    sim_dt[, time := 1951 + time_idx]  # Adjust base year as needed
    
    # Calculate summary metrics per simulation (seed) - matching your filtering approach
    sim_dt[, maxF := max(F, na.rm = TRUE), by = seed]
    sim_dt[, avgF := mean(F, na.rm = TRUE), by = seed]
    sim_dt[, mindep := min(dep, na.rm = TRUE), by = seed]
    sim_dt[, minN := min(n, na.rm = TRUE), by = seed]
    sim_dt[, maxCatch := max(predicted_catch, na.rm = TRUE), by = seed]
    sim_dt[, minCatch := min(predicted_catch, na.rm = TRUE), by = seed]
    sim_dt[, avgCatch := mean(predicted_catch, na.rm = TRUE), by = seed]
    
    # Period-specific metrics
    sim_dt[, avg_catch_early := mean(predicted_catch[time < 1962], na.rm = TRUE), by = seed]
    
    # Define survival filters (basic biological constraints)
    seed_surv = unique(sim_dt[mindep > 0.01 & minN > 15000 & shape < 20 & !is.infinite(maxF)]$seed)
    
    # F-based filters 
    seed_f = unique(sim_dt[seed %in% seed_surv & maxF < 2.5 & avgF < 0.8]$seed)
    
    # Catch-based filters (requires observed data for comparison)
    if (!is.null(catch_effort_dt)) {
        # Calculate observed catch statistics
        obs_catch_early_avg = mean(catch_effort_dt[year < 1962]$total_catch, na.rm = TRUE)
        obs_catch_overall_avg = mean(catch_effort_dt$total_catch, na.rm = TRUE)
        obs_max_catch = max(catch_effort_dt$total_catch, na.rm = TRUE)
        obs_min_catch = min(catch_effort_dt$total_catch, na.rm = TRUE)
        
        # Apply catch filters
        seed_catch = unique(sim_dt[
            seed %in% seed_f & 
            
            # Overall bounds
            maxCatch > 0.5 * obs_max_catch &
            maxCatch < 3.0 * obs_max_catch &
            minCatch > 0.5 * obs_min_catch &
            avgCatch > 0.5 * obs_catch_overall_avg &
            avgCatch < 2.0 * obs_catch_overall_avg &
            
            # Early period pattern  
            avg_catch_early > 0.5 * obs_catch_early_avg &
            avg_catch_early < 3 * obs_catch_early_avg
        ]$seed)
        
        final_seeds = seed_catch
        
        if (verbose) {
            cat("Filtering summary:\n")
            cat("  Total simulations:", length(unique(sim_dt$seed)), "\n")
            cat("  Passed survival filters:", length(seed_surv), "\n")
            cat("  Passed F filters:", length(seed_f), "\n") 
            cat("  Passed catch filters:", length(seed_catch), "\n")
        }
        
    } else {
        warning("No catch_effort_dt provided. Applying only survival and F filters.")
        final_seeds = seed_f
        
        if (verbose) {
            cat("Filtering summary:\n")
            cat("  Total simulations:", length(unique(sim_dt$seed)), "\n")
            cat("  Passed survival filters:", length(seed_surv), "\n")
            cat("  Passed F filters:", length(seed_f), "\n")
        }
    }
    
    if (length(final_seeds) < 50) {
        warning("Very few simulations passed filters (", length(final_seeds), 
                "). Consider relaxing constraints.")
    }
    
    # Extract transformed parameters for filtered simulations
    filtered_hmc = hmc_dt[seed %in% final_seeds]
    
    # Get parameter samples for prior fitting (on natural scales)
    param_vars = c("logK", "r", "shape", "qeff", "rho", "sigma_qdev", "sigmap")
    
    # Check what's available and extract
    param_samples = filtered_hmc[variable %in% param_vars, .(seed, variable, value)] %>%
        dcast(., seed ~ variable, value.var = "value")
    
    # Check for required parameters
    missing_params = setdiff(param_vars, names(param_samples))
    if (length(missing_params) > 0) {
        if (verbose) cat("Missing parameters (will skip):", paste(missing_params, collapse = ", "), "\n")
    }
    
    if (verbose) cat("Retained", nrow(param_samples), "samples for prior updating\n")
    
    # MLE fitting functions
    # Multivariate normal MLE with Cholesky parameterization
    mv_mle_fit = function(par, data) {
        mu = par[1:3]
        # Use Cholesky decomposition parameterization for positive definiteness
        L_vec = par[4:9]  # 6 parameters for lower triangular Cholesky factor
        L = matrix(0, 3, 3)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)  # This ensures positive definiteness
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(data, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }
    
    # Univariate normal MLE
    univariate_mle_fit = function(par, data) {
        mu = par[1]
        sigma = par[2]
        if(sigma <= 0) return(1e10)
        tryCatch({
            -sum(dnorm(data, mu, sigma, log=TRUE))
        }, error = function(e) 1e10)
    }
    
    # Bivariate normal MLE with Cholesky parameterization
    bivariate_mle_fit = function(par, data) {
        mu = par[1:2]
        # Use Cholesky decomposition parameterization for positive definiteness
        L_vec = par[3:5]  # 3 parameters for lower triangular Cholesky factor
        L = matrix(0, 2, 2)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)  # This ensures positive definiteness
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(data, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }
    
    # Fit priors using MLE
    priors_list = list()
    
    # Create multivariate parameter matrix for logK, log_r, log_shape
    if (all(c("logK", "r", "shape") %in% names(param_samples))) {
        mv_params_filtered = cbind(
            logK = param_samples$logK,
            log_r = log(param_samples$r),
            log_shape = log(param_samples$shape)
        )
        
        # Initial values using sample moments
        init_cov_filtered = cov(mv_params_filtered)
        init_chol_filtered = chol(init_cov_filtered)  # Upper triangular
        init_chol_lower_filtered = t(init_chol_filtered)       # Convert to lower triangular
        
        init_filtered = c(colMeans(mv_params_filtered), 
                        init_chol_lower_filtered[lower.tri(init_chol_lower_filtered, diag=TRUE)])
        
        # Fit with MLE
        mv_fit_filtered = nlminb(init_filtered, mv_mle_fit, data = mv_params_filtered,
                                control = list(eval.max = 1000, iter.max = 500))
        
        # Check convergence and extract parameters
        if(mv_fit_filtered$convergence != 0 && mv_fit_filtered$convergence != 4) {
            warning("Multivariate MLE fitting may not have converged properly. Using sample moments instead.")
            mv_mean_new = colMeans(mv_params_filtered)
            mv_cov_new = cov(mv_params_filtered)
        } else {
            # Extract fitted parameters from Cholesky parameterization
            mv_mean_new = mv_fit_filtered$par[1:3]
            L_fitted_filtered = matrix(0, 3, 3)
            L_fitted_filtered[lower.tri(L_fitted_filtered, diag=TRUE)] = mv_fit_filtered$par[4:9]
            mv_cov_new = L_fitted_filtered %*% t(L_fitted_filtered)
        }
        
        mv_cor_new = cov2cor(mv_cov_new)
        
        # Add names
        names(mv_mean_new) = dimnames(mv_cov_new)[[1]] = dimnames(mv_cov_new)[[2]] = c("logK", "log_r", "log_shape")
        
    } else {
        if (verbose) cat("Cannot create multivariate prior - missing logK, r, or shape\n")
        mv_mean_new = mv_cov_new = mv_cor_new = NULL
    }
    
    # qeff prior (log-normal) - MLE on log scale
    if ("qeff" %in% names(param_samples)) {
        qeff_log_samples = log(param_samples$qeff)
        init_qeff = c(mean(qeff_log_samples), sd(qeff_log_samples))
        
        qeff_fit = nlminb(init_qeff, univariate_mle_fit, data = qeff_log_samples,
                         control = list(eval.max = 1000, iter.max = 500))
        
        if(qeff_fit$convergence != 0 && qeff_fit$convergence != 4) {
            warning("qeff MLE fitting may not have converged properly. Using sample moments instead.")
            priors_list$qeff_pars = c(mean(qeff_log_samples), sd(qeff_log_samples))
        } else {
            priors_list$qeff_pars = qeff_fit$par
        }
        names(priors_list$qeff_pars) = c("meanlog", "sdlog")
    }
    
    # rho and sigma_qdev bivariate prior - MLE
    if (all(c("rho", "sigma_qdev") %in% names(param_samples))) {
        qdev_params_filtered = cbind(
            atanh_rho = atanh(param_samples$rho),
            log_sigma_qdev = log(param_samples$sigma_qdev)
        )
        
        # Initial values using sample moments
        init_cov_qdev = cov(qdev_params_filtered)
        init_chol_qdev = chol(init_cov_qdev)  # Upper triangular
        init_chol_lower_qdev = t(init_chol_qdev)       # Convert to lower triangular
        
        init_qdev = c(colMeans(qdev_params_filtered),
                     init_chol_lower_qdev[lower.tri(init_chol_lower_qdev, diag=TRUE)])
        
        # Fit with MLE
        qdev_fit = nlminb(init_qdev, bivariate_mle_fit, data = qdev_params_filtered,
                         control = list(eval.max = 1000, iter.max = 500))
        
        # Check convergence and extract parameters
        if(qdev_fit$convergence != 0 && qdev_fit$convergence != 4) {
            warning("qdev MLE fitting may not have converged properly. Using sample moments instead.")
            qdev_mean_new = colMeans(qdev_params_filtered)
            qdev_cov_new = cov(qdev_params_filtered)
        } else {
            # Extract fitted parameters from Cholesky parameterization
            qdev_mean_new = qdev_fit$par[1:2]
            L_fitted_qdev = matrix(0, 2, 2)
            L_fitted_qdev[lower.tri(L_fitted_qdev, diag=TRUE)] = qdev_fit$par[3:5]
            qdev_cov_new = L_fitted_qdev %*% t(L_fitted_qdev)
        }
        
        qdev_cor_new = cov2cor(qdev_cov_new)
        
        priors_list$mv_qdev_mean = qdev_mean_new
        priors_list$mv_qdev_cov = qdev_cov_new
        priors_list$mv_qdev_cor = qdev_cor_new
    }
    
    # sigmap prior (log-normal) - MLE on log scale
    if ("sigmap" %in% names(param_samples)) {
        sigmap_log_samples = log(param_samples$sigmap)
        init_sigmap = c(mean(sigmap_log_samples), sd(sigmap_log_samples))
        
        sigmap_fit = nlminb(init_sigmap, univariate_mle_fit, data = sigmap_log_samples,
                           control = list(eval.max = 1000, iter.max = 500))
        
        if(sigmap_fit$convergence != 0 && sigmap_fit$convergence != 4) {
            warning("sigmap MLE fitting may not have converged properly. Using sample moments instead.")
            priors_list$sigmap_pars = c(mean(sigmap_log_samples), sd(sigmap_log_samples))
        } else {
            priors_list$sigmap_pars = sigmap_fit$par
        }
        names(priors_list$sigmap_pars) = c("meanlog", "sdlog")
    }
    
    # depletion prior (from x at specific time point, e.g., t_dep) - MLE on log scale
    if (length(grep("x[",available_vars,fixed=TRUE))>0) {
        # Find the depletion time point (usually t_dep = 37 for 1988)
        t_dep = 37  # Adjust as needed
        depl_var = paste0("x[", t_dep, "]")
        if (depl_var %in% available_vars) {
            depl_samples = filtered_hmc[variable == depl_var, value]
            if (length(depl_samples) > 0) {
                depl_log_samples = log(depl_samples)
                init_depl = c(mean(depl_log_samples), sd(depl_log_samples))
                
                depl_fit = nlminb(init_depl, univariate_mle_fit, data = depl_log_samples,
                                 control = list(eval.max = 1000, iter.max = 500))
                
                if(depl_fit$convergence != 0 && depl_fit$convergence != 4) {
                    warning("depletion MLE fitting may not have converged properly. Using sample moments instead.")
                    priors_list$depletion_pars = c(mean(depl_log_samples), sd(depl_log_samples))
                } else {
                    priors_list$depletion_pars = depl_fit$par
                }
                names(priors_list$depletion_pars) = c("meanlog", "sdlog")
            }
        }
    }
    
    # Create updated stan.data structure
    if (is.null(original_stan_data)) {
        stan_data_updated = list()
        
        # Add multivariate prior if available
        if (!is.null(mv_mean_new)) {
            stan_data_updated$mv_prior_mean = mv_mean_new
            stan_data_updated$mv_prior_sd = sqrt(diag(mv_cov_new))
            stan_data_updated$mv_prior_corr = mv_cor_new
        }
        
        # Add individual priors
        if ("qeff_pars" %in% names(priors_list)) {
            stan_data_updated$prior_qeff_meanlog = priors_list$qeff_pars[["meanlog"]]
            stan_data_updated$prior_qeff_sdlog = priors_list$qeff_pars[["sdlog"]]
        }
        
        if ("mv_qdev_mean" %in% names(priors_list)) {
            stan_data_updated$mv_qdev_prior_mean = priors_list$mv_qdev_mean
            stan_data_updated$mv_qdev_prior_sd = sqrt(diag(priors_list$mv_qdev_cov))
            stan_data_updated$mv_qdev_prior_corr = priors_list$mv_qdev_cor
        }
        
        if ("sigmap_pars" %in% names(priors_list)) {
            stan_data_updated$PriorMean_logsigmap = priors_list$sigmap_pars[["meanlog"]]
            stan_data_updated$PriorSD_logsigmap = priors_list$sigmap_pars[["sdlog"]]
        }
        
        if ("depletion_pars" %in% names(priors_list)) {
            stan_data_updated$prior_depletion_meanlog = priors_list$depletion_pars[["meanlog"]]
            stan_data_updated$prior_depletion_sdlog = priors_list$depletion_pars[["sdlog"]]
        }
        
        if (verbose) cat("\nUpdated prior parameters calculated. Add your observation data to complete stan.data list.\n")
        
    } else {
        # Update existing stan.data
        stan_data_updated = original_stan_data
        
        # Update multivariate prior
        if (!is.null(mv_mean_new)) {
            stan_data_updated$mv_prior_mean = mv_mean_new
            stan_data_updated$mv_prior_sd = sqrt(diag(mv_cov_new))
            stan_data_updated$mv_prior_corr = mv_cor_new
        }
        
        # Update individual priors
        if ("qeff_pars" %in% names(priors_list)) {
            stan_data_updated$prior_qeff_meanlog = priors_list$qeff_pars[["meanlog"]]
            stan_data_updated$prior_qeff_sdlog = priors_list$qeff_pars[["sdlog"]]
        }
        
        if ("mv_qdev_mean" %in% names(priors_list)) {
            stan_data_updated$mv_qdev_prior_mean = priors_list$mv_qdev_mean
            stan_data_updated$mv_qdev_prior_sd = sqrt(diag(priors_list$mv_qdev_cov))
            stan_data_updated$mv_qdev_prior_corr = priors_list$mv_qdev_cor
        }
        
        if ("sigmap_pars" %in% names(priors_list)) {
            stan_data_updated$PriorMean_logsigmap = priors_list$sigmap_pars[["meanlog"]]
            stan_data_updated$PriorSD_logsigmap = priors_list$sigmap_pars[["sdlog"]]
        }
        
        if ("depletion_pars" %in% names(priors_list)) {
            stan_data_updated$prior_depletion_meanlog = priors_list$depletion_pars[["meanlog"]]
            stan_data_updated$prior_depletion_sdlog = priors_list$depletion_pars[["sdlog"]]
        }
        
        if (verbose) cat("\nExisting stan.data list updated with new priors.\n")
    }
    
    # Print summary of changes with comparisons
    if (verbose) {
        cat("\n=== PRIOR UPDATE SUMMARY (MLE Fitted) ===\n")
        
        if (!is.null(mv_mean_new) && !is.null(original_stan_data)) {
            cat("Multivariate prior for (logK, log_r, log_shape):\n")
            cat("  ORIGINAL Mean:", round(original_stan_data$mv_prior_mean, 3), "\n")
            cat("  UPDATED  Mean:", round(mv_mean_new, 3), "\n")
            cat("  ORIGINAL SD:  ", round(original_stan_data$mv_prior_sd, 3), "\n")
            cat("  UPDATED  SD:  ", round(sqrt(diag(mv_cov_new)), 3), "\n")
            cat("  ORIGINAL Correlations:\n")
            print(round(original_stan_data$mv_prior_corr, 3))
            cat("  UPDATED Correlations:\n")
            print(round(mv_cor_new, 3))
        } else if (!is.null(mv_mean_new)) {
            cat("Multivariate prior for (logK, log_r, log_shape):\n")
            cat("  Mean:", round(mv_mean_new, 3), "\n")
            cat("  SD:", round(sqrt(diag(mv_cov_new)), 3), "\n")
            cat("  Correlations:\n")
            print(round(mv_cor_new, 3))
        }
        
        if ("qeff_pars" %in% names(priors_list)) {
            cat("\nqeff prior (lognormal, MLE fitted):\n")
            if (!is.null(original_stan_data) && !is.null(original_stan_data$prior_qeff_meanlog)) {
                cat("  ORIGINAL Mean log:", round(original_stan_data$prior_qeff_meanlog, 3), "\n")
                cat("  UPDATED  Mean log:", round(priors_list$qeff_pars[["meanlog"]], 3), "\n")
                cat("  ORIGINAL SD log:  ", round(original_stan_data$prior_qeff_sdlog, 3), "\n")
                cat("  UPDATED  SD log:  ", round(priors_list$qeff_pars[["sdlog"]], 3), "\n")
            } else {
                cat("  Mean log:", round(priors_list$qeff_pars[["meanlog"]], 3), "\n")
                cat("  SD log:", round(priors_list$qeff_pars[["sdlog"]], 3), "\n")
            }
        }
        
        if ("mv_qdev_mean" %in% names(priors_list)) {
            cat("\nBivariate qdev prior for (atanh_rho, log_sigma_qdev, MLE fitted):\n")
            if (!is.null(original_stan_data) && !is.null(original_stan_data$mv_qdev_prior_mean)) {
                cat("  ORIGINAL Mean:", round(original_stan_data$mv_qdev_prior_mean, 3), "\n")
                cat("  UPDATED  Mean:", round(priors_list$mv_qdev_mean, 3), "\n")
                cat("  ORIGINAL SD:  ", round(original_stan_data$mv_qdev_prior_sd, 3), "\n")
                cat("  UPDATED  SD:  ", round(sqrt(diag(priors_list$mv_qdev_cov)), 3), "\n")
                cat("  ORIGINAL Correlation:\n")
                print(round(original_stan_data$mv_qdev_prior_corr, 3))
                cat("  UPDATED Correlation:\n")
                print(round(priors_list$mv_qdev_cor, 3))
            } else {
                cat("  Mean:", round(priors_list$mv_qdev_mean, 3), "\n")
                cat("  SD:", round(sqrt(diag(priors_list$mv_qdev_cov)), 3), "\n")
                cat("  Correlation:\n")
                print(round(priors_list$mv_qdev_cor, 3))
            }
        }
        
        if ("sigmap_pars" %in% names(priors_list)) {
            cat("\nsigmap prior (lognormal, MLE fitted):\n")
            if (!is.null(original_stan_data) && !is.null(original_stan_data$PriorMean_logsigmap)) {
                cat("  ORIGINAL Mean log:", round(original_stan_data$PriorMean_logsigmap, 3), "\n")
                cat("  UPDATED  Mean log:", round(priors_list$sigmap_pars[["meanlog"]], 3), "\n")
                cat("  ORIGINAL SD log:  ", round(original_stan_data$PriorSD_logsigmap, 3), "\n")
                cat("  UPDATED  SD log:  ", round(priors_list$sigmap_pars[["sdlog"]], 3), "\n")
            } else {
                cat("  Mean log:", round(priors_list$sigmap_pars[["meanlog"]], 3), "\n")
                cat("  SD log:", round(priors_list$sigmap_pars[["sdlog"]], 3), "\n")
            }
        }
        
        if ("depletion_pars" %in% names(priors_list)) {
            cat("\nDepletion prior (lognormal, MLE fitted):\n")
            if (!is.null(original_stan_data) && !is.null(original_stan_data$prior_depletion_meanlog)) {
                cat("  ORIGINAL Mean log:", round(original_stan_data$prior_depletion_meanlog, 3), "\n")
                cat("  UPDATED  Mean log:", round(priors_list$depletion_pars[["meanlog"]], 3), "\n")
                cat("  ORIGINAL SD log:  ", round(original_stan_data$prior_depletion_sdlog, 3), "\n")
                cat("  UPDATED  SD log:  ", round(priors_list$depletion_pars[["sdlog"]], 3), "\n")
            } else {
                cat("  Mean log:", round(priors_list$depletion_pars[["meanlog"]], 3), "\n")
                cat("  SD log:", round(priors_list$depletion_pars[["sdlog"]], 3), "\n")
            }
        }


        # Create comparison plots
        cat("\n=== GENERATING COMPARISON PLOTS ===\n")
        library(ggplot2)
        library(GGally)
        
        # Determine which parameters to include in plots
        plot_vars = c()
        if (!is.null(mv_mean_new)) plot_vars = c(plot_vars, "logK", "log_r", "log_shape")
        if ("qeff_pars" %in% names(priors_list)) plot_vars = c(plot_vars, "log_qeff")
        if ("mv_qdev_mean" %in% names(priors_list)) plot_vars = c(plot_vars, "atanh_rho", "log_sigma_qdev")
        if ("sigmap_pars" %in% names(priors_list)) plot_vars = c(plot_vars, "log_sigmap")
        if ("depletion_pars" %in% names(priors_list)) plot_vars = c(plot_vars, "log_depletion")
        
        if (length(plot_vars) > 0) {
            set.seed(123)
            n_samples = 1000
            
            # 1. Original prior samples (if available)
            if (!is.null(original_stan_data)) {
                original_samples_list = list()
                
                # Multivariate samples
                if (!is.null(original_stan_data$mv_prior_mean)) {
                    original_mv = mvrnorm(n_samples, 
                                        original_stan_data$mv_prior_mean, 
                                        diag(original_stan_data$mv_prior_sd) %*% 
                                        original_stan_data$mv_prior_corr %*% 
                                        diag(original_stan_data$mv_prior_sd))
                    original_samples_list$logK = original_mv[,1]
                    original_samples_list$log_r = original_mv[,2] 
                    original_samples_list$log_shape = original_mv[,3]
                }
                
                # Individual parameter samples
                if (!is.null(original_stan_data$prior_qeff_meanlog)) {
                    original_samples_list$log_qeff = rnorm(n_samples, 
                                                          original_stan_data$prior_qeff_meanlog,
                                                          original_stan_data$prior_qeff_sdlog)
                }
                
                if (!is.null(original_stan_data$mv_qdev_prior_mean)) {
                    original_qdev = mvrnorm(n_samples,
                                          original_stan_data$mv_qdev_prior_mean,
                                          diag(original_stan_data$mv_qdev_prior_sd) %*% 
                                          original_stan_data$mv_qdev_prior_corr %*% 
                                          diag(original_stan_data$mv_qdev_prior_sd))
                    original_samples_list$atanh_rho = original_qdev[,1]
                    original_samples_list$log_sigma_qdev = original_qdev[,2]
                }
                
                if (!is.null(original_stan_data$PriorMean_logsigmap)) {
                    original_samples_list$log_sigmap = rnorm(n_samples,
                                                           original_stan_data$PriorMean_logsigmap,
                                                           original_stan_data$PriorSD_logsigmap)
                }
                
                if (!is.null(original_stan_data$prior_depletion_meanlog)) {
                    original_samples_list$log_depletion = rnorm(n_samples,
                                                              original_stan_data$prior_depletion_meanlog,
                                                              original_stan_data$prior_depletion_sdlog)
                }
                
                # Create data frame
                original_df = data.frame(original_samples_list[plot_vars])
                original_df$type = "Original"
            }
            
            # 2. Actual filtered samples
            filtered_indices = sample(1:nrow(param_samples), min(n_samples, nrow(param_samples)))
            filtered_samples_list = list()
            
            if ("logK" %in% plot_vars) filtered_samples_list$logK = param_samples$logK[filtered_indices]
            if ("log_r" %in% plot_vars) filtered_samples_list$log_r = log(param_samples$r[filtered_indices])
            if ("log_shape" %in% plot_vars) filtered_samples_list$log_shape = log(param_samples$shape[filtered_indices])
            if ("log_qeff" %in% plot_vars) filtered_samples_list$log_qeff = log(param_samples$qeff[filtered_indices])
            if ("atanh_rho" %in% plot_vars) filtered_samples_list$atanh_rho = atanh(param_samples$rho[filtered_indices])
            if ("log_sigma_qdev" %in% plot_vars) filtered_samples_list$log_sigma_qdev = log(param_samples$sigma_qdev[filtered_indices])
            if ("log_sigmap" %in% plot_vars) filtered_samples_list$log_sigmap = log(param_samples$sigmap[filtered_indices])
            
            # For depletion, extract from filtered_hmc
            if ("log_depletion" %in% plot_vars) {
                t_dep = 37
                depl_var = paste0("x[", t_dep, "]")
                if (depl_var %in% available_vars) {
                    depl_filtered_samples = filtered_hmc[variable == depl_var, value]
                    if (length(depl_filtered_samples) >= length(filtered_indices)) {
                        filtered_samples_list$log_depletion = log(depl_filtered_samples[filtered_indices])
                    }
                }
            }
            
            filtered_df = data.frame(filtered_samples_list[plot_vars])
            filtered_df$type = "Actual Filtered"
            
            # 3. New fitted prior samples
            fitted_samples_list = list()
            
            # Multivariate samples
            if (!is.null(mv_mean_new)) {
                fitted_mv = mvrnorm(n_samples, mv_mean_new, mv_cov_new)
                fitted_samples_list$logK = fitted_mv[,1]
                fitted_samples_list$log_r = fitted_mv[,2]
                fitted_samples_list$log_shape = fitted_mv[,3]
            }
            
            # Individual parameter samples
            if ("qeff_pars" %in% names(priors_list)) {
                fitted_samples_list$log_qeff = rnorm(n_samples,
                                                    priors_list$qeff_pars[["meanlog"]],
                                                    priors_list$qeff_pars[["sdlog"]])
            }
            
            if ("mv_qdev_mean" %in% names(priors_list)) {
                fitted_qdev = mvrnorm(n_samples, priors_list$mv_qdev_mean, priors_list$mv_qdev_cov)
                fitted_samples_list$atanh_rho = fitted_qdev[,1]
                fitted_samples_list$log_sigma_qdev = fitted_qdev[,2]
            }
            
            if ("sigmap_pars" %in% names(priors_list)) {
                fitted_samples_list$log_sigmap = rnorm(n_samples,
                                                      priors_list$sigmap_pars[["meanlog"]],
                                                      priors_list$sigmap_pars[["sdlog"]])
            }
            
            if ("depletion_pars" %in% names(priors_list)) {
                fitted_samples_list$log_depletion = rnorm(n_samples,
                                                         priors_list$depletion_pars[["meanlog"]],
                                                         priors_list$depletion_pars[["sdlog"]])
            }
            
            fitted_df = data.frame(fitted_samples_list[plot_vars])
            fitted_df$type = "Fitted Updated"
            
            # Combine data frames
            if (!is.null(original_stan_data) && exists("original_df")) {
                combined_df = rbind(original_df, filtered_df, fitted_df)
                combined_df$type = factor(combined_df$type, levels = c("Original", "Actual Filtered", "Fitted Updated"))
                colors = c("Original" = "gray", "Actual Filtered" = "blue", "Fitted Updated" = "darkred")
            } else {
                combined_df = rbind(filtered_df, fitted_df)  
                combined_df$type = factor(combined_df$type, levels = c("Actual Filtered", "Fitted Updated"))
                colors = c("Actual Filtered" = "blue", "Fitted Updated" = "darkred")
            }
            
            # Create ggpairs plot
            p_pairs = ggpairs(combined_df, columns=1:length(plot_vars), aes(color=type, alpha=0.7),
                            upper=list(continuous="cor"),
                            lower=list(continuous="points"), 
                            diag=list(continuous="densityDiag")) +
                    scale_color_manual(values = colors) +
                    scale_fill_manual(values = colors) +
                    theme(
                            text = element_text(size = 12),
                            panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                            panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                            panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                            strip.background = element_rect(fill = "white"),
                            legend.key = element_rect(fill = "white")
                        )
            
            print(p_pairs)
            cat("Plot generated showing comparison of original, filtered, and updated priors\n")
            
        } else {
            cat("No parameters available for plotting\n")
            p_pairs = "No parameters available for plotting"
        }

    }
    
    # Return both the updated data structure and filtering summary
    result = list(
        stan_data = stan_data_updated,
        filtering_summary = list(
            total_simulations = length(unique(sim_dt$seed)),
            survival_passed = length(seed_surv),
            f_passed = length(seed_f),
            final_passed = length(final_seeds),
            final_seeds = final_seeds
        ),
        priors_fitted = names(priors_list),
        plot = p_pairs
    )
    
    return(result)
}

# Example usage:
# result = update_priors_from_ppc("data/output/model_runs/0006-2024cpueEffortQeff_0-ppc", 
#                                 original_stan_data = stan.data,
#                                 catch_effort_dt = catch_effort_annual)
# updated_stan_data = result$stan_data
# filtering_info = result$filtering_summary
