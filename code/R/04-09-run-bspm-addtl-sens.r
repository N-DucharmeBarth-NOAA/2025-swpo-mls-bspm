# Nicholas Ducharme-Barth
# 2025/07/03
# Comprehensive sensitivity analysis for BSPM models

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(rstan)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for model outputs
    dir.create(file.path(proj_dir,"data","output","model_runs"), showWarnings = FALSE, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load inputs
    catch_dt = fread(file.path(proj_dir,"data","input","catch.csv"))
    
    # Load effort data 
    effort_dt = fread(file.path(proj_dir, "data", "input", "WCPFC_L_PUBLIC_BY_FLAG_YR.CSV")) %>%
                .[, lat_short_d := ifelse(grepl("S$", lat_short), 
                                        -as.numeric(gsub("[NS]", "", lat_short)), 
                                        as.numeric(gsub("[NS]", "", lat_short))) + 2.5] %>%
                .[, lon_short_d := ifelse(grepl("W$", lon_short), 
                                        360 - as.numeric(gsub("[EW]", "", lon_short)), 
                                            as.numeric(gsub("[EW]", "", lon_short))) + 2.5] %>%
                .[lat_short_d > -60 & lat_short_d < 0 & lon_short_d < 230] %>%
                .[!(lat_short_d > -5 & lon_short_d > 210)] %>%
                .[,.(effort_scaled = sum(hhooks)/1e6),by=yy] %>%
                .[yy %in% 1952:2022] %>%
                setnames(., "yy", "time")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load reference model data and extract prior parameters
    
    # Function to extract prior parameters from stan_data.csv files
    extract_prior_params = function(model_id) {
        file_path = file.path(proj_dir, "data", "output", "model_runs", model_id, "stan_data.csv")
        if(!file.exists(file_path)) {
            warning(paste("File not found:", file_path))
            return(NULL)
        }
        dt = fread(file_path)
        return(dt)
    }
    
    # Function to calculate mean F distribution from sigmaF samples (for 0005 model)
    calc_F = function(model_id) {
        file_path = file.path(proj_dir, "data", "output", "model_runs", model_id, "hmc_samples.csv")
        if(!file.exists(file_path)) {
            warning(paste("File not found:", file_path))
            return(list(mean_F_mean = NA, mean_F_sd = NA))
        }
        dt = fread(file_path)
        samples = dt[name == "F", .(iter,row,value)] %>%
                  .[,.(mean(value),sd(value)),by=iter]
        if(length(samples) == 0) {
            warning(paste("F parameter not found in", model_id))
            return(list(mean_F_mean = NA, mean_F_sd = NA))
        }
        # Each sigmaF sample represents the mean F for that iteration
        return(list(mean_F_mean = mean(samples$V1, na.rm = TRUE), 
                   mean_F_sd = mean(samples$V2, na.rm = TRUE)))
    }

    # Load reference model data
    ref_0003_data = extract_prior_params("0003-2024cpueFPrior_0")
    ref_0005_data = extract_prior_params("0005-2024cpueMVPrior_0") 
    ref_0025_data = extract_prior_params("0025-dwfn-c0.2-e0.3-s5_0")
    
    # Extract F-based statistics from reference models
    ref_0005_meanF_stats = calc_F("0005-2024cpueMVPrior_0")
    ref_0025_F_variability = calc_F("0025-dwfn-c0.2-e0.3-s5_0")
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare common data inputs (following 04-08 pattern)
    
    # Load updated prior parameters from pushforward analysis
    mv_mean_catch = read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_mean_catch.csv"))$x
    mv_cov_catch = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_cov_catch.csv"))[,-1])
    mv_cor_catch = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_cor_catch.csv"))[,-1])

    # add x0 to the mix for 4D priors
    mv_cor_catch_4d = matrix(c(
    1.0,   mv_cor_catch[1,2], mv_cor_catch[1,3], 0.0,    # logK row
    mv_cor_catch[1,2], 1.0,   mv_cor_catch[2,3], 0.0,    # log_r row  
    mv_cor_catch[1,3], mv_cor_catch[2,3], 1.0,   0.0,    # log_shape row
        0.0,     0.0,    0.0,   1.0      # log_x0 row (uncorrelated)
    ), nrow=4, ncol=4)

    mv_mean_catch_4d = c(mv_mean_catch, log(1))  # tight around 0.95 initial depletion
    mv_prior_sd_4d = c(sqrt(diag(mv_cov_catch)), 0.025)           # small SD for x0
    
    # Independent qeff prior
    qeff_pars_catch = read.csv(file.path(proj_dir,"data","output","pushforward","q","qeff_pars_catch.csv"))$x
    
    # Bivariate rho/sigma_qdev prior
    mv_qdev_mean = read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_mean.csv"))$x
    mv_qdev_cov = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_cov.csv"))[,-1])
    mv_qdev_cor = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_cor.csv"))[,-1])
    
    # Other priors
    dep_pars = read.csv(file.path(proj_dir,"data","output","rel_dep_ssb_pars.csv"))
    sigmaF_par = read.csv(file.path(proj_dir,"data","output","sigmaF_pars_catch.csv"))

    # Extract standard deviations from covariance matrix (3D now)
    mv_prior_sd = sqrt(diag(mv_cov_catch))
    
    # Extract standard deviations for qdev parameters (2D)
    mv_qdev_prior_sd = sqrt(diag(mv_qdev_cov))

    # Validate correlation matrices are positive definite
    if(min(eigen(mv_cor_catch)$values) < 1e-10) {
        warning("3D correlation matrix has very small eigenvalues, adding regularization")
        mv_cor_catch = mv_cor_catch + diag(1e-8, 3)
    }

    if(min(eigen(mv_cor_catch_4d)$values) < 1e-10) {
        warning("4D correlation matrix has very small eigenvalues, adding regularization")
        mv_cor_catch_4d = mv_cor_catch_4d + diag(1e-8, 4)
    }
    
    if(min(eigen(mv_qdev_cor)$values) < 1e-10) {
        warning("2D correlation matrix has very small eigenvalues, adding regularization")  
        mv_qdev_cor = mv_qdev_cor + diag(1e-8, 2)
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load additional cpue indices
    au_cpue_dt = fread(file.path(proj_dir,"data","input","AU_cpue.csv"))
    nz_cpue_dt = fread(file.path(proj_dir,"data","input","NZ_cpue.csv"))
    obs_cpue_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP.csv"))
    obs_cpue_no_PF_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_no_PF.csv"))
    obs_cpue_PF_only_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_PF_only.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare data matrices (following existing pattern)
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    # Merge with effort data
    catch_effort_annual = merge(catch_annual, effort_dt, by.x = "year", by.y = "time", all.x = TRUE)
    
    time_years = catch_effort_annual$year
    n_years = length(time_years)
    
    # Initialize matrices for multiple indices
    n_indices = 6
    index_mat = matrix(-999, nrow = n_years, ncol = n_indices)
    se_mat = matrix(-999, nrow = n_years, ncol = n_indices)
    mean_se = rep(NA, n_indices)
    
    # Process DWFN CPUE (Index 1)
    cpue_dt = fread(file.path(proj_dir,"data","input","cpue.csv")) %>% .[,time:=floor(Time)]
    mean_se[1] = mean(cpue_dt$SE,na.rm=TRUE)
    
    for(i in 1:nrow(cpue_dt)) {
        year_idx = which(time_years == cpue_dt$time[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 1] = cpue_dt$Obs[i]/mean(cpue_dt$Obs, na.rm = TRUE)
            se_mat[year_idx, 1] = cpue_dt$SE[i]/mean_se[1]
        }
    }
    
    # australia cpue
    mean_se[2] = mean(au_cpue_dt$se_log,na.rm=TRUE)
    for(i in 1:nrow(au_cpue_dt)) {
        year_idx = which(time_years == au_cpue_dt$year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 2] = au_cpue_dt$obs[i]/mean(au_cpue_dt$obs)
            se_mat[year_idx, 2] = au_cpue_dt$se_log[i]/mean_se[2]
        }
    }

    # new zealand cpue
    mean_se[3] = mean(nz_cpue_dt$se_log,na.rm=TRUE)
    for(i in 1:nrow(nz_cpue_dt)) {
        year_idx = which(time_years == nz_cpue_dt$year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 3] = nz_cpue_dt$ob[i]/mean(nz_cpue_dt$ob)
            se_mat[year_idx, 3] = nz_cpue_dt$se_log[i]/mean_se[3]
        }
    }

    # Observer CPUE
    # calc SE from quantiles
    obs_cpue_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[4] = mean(obs_cpue_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_dt)) {
        year_idx = which(time_years == obs_cpue_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 4] = obs_cpue_dt$Estimate[i]/mean(obs_cpue_dt$Estimate)
            se_mat[year_idx, 4] = obs_cpue_dt$se_from_quantiles[i]/mean_se[4]
        }
    }

    # calc SE from quantiles
    obs_cpue_no_PF_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[5] = mean(obs_cpue_no_PF_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_no_PF_dt)) {
        year_idx = which(time_years == obs_cpue_no_PF_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 5] = obs_cpue_no_PF_dt$Estimate[i]/mean(obs_cpue_no_PF_dt$Estimate)
            se_mat[year_idx, 5] = obs_cpue_no_PF_dt$se_from_quantiles[i]/mean_se[5]
        }
    }

    # calc SE from quantiles
    obs_cpue_PF_only_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[6] = mean(obs_cpue_PF_only_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_PF_only_dt)) {
        year_idx = which(time_years == obs_cpue_PF_only_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 6] = obs_cpue_PF_only_dt$Estimate[i]/mean(obs_cpue_PF_only_dt$Estimate)
            se_mat[year_idx, 6] = obs_cpue_PF_only_dt$se_from_quantiles[i]/mean_se[6]
        }
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name_vec = c("bspm_estq_softdep_mvprior","bspm_estq_flex","bspm_estq_softdep_mvprior_x0","bspm_estq_flex_x0","bspm_estF_mvprior")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop sensitivity model grid using structured approach

    # Define sensitivity configuration structure
    create_sensitivity_config = function(model_type, exec_code, start_year = 1952, cpue_index = "dwfn",
                                       catch_scenario = "0.2flat", sigma_edev = 0.3, step_scenario = "5reg",
                                       use_alt_qeff_prior = FALSE, use_alt_shape_prior = FALSE, 
                                       x0_prior_value = NULL, use_mv_baseline = FALSE, 
                                       use_alt_sigmaF_prior = FALSE, use_alt_qeff_from_meanF = FALSE,
                                       use_flex_variant = FALSE) {
        
        # Generate descriptive name
        name_parts = c()
        if(use_alt_qeff_prior) name_parts = c(name_parts, "qalt")
        if(use_alt_qeff_from_meanF) name_parts = c(name_parts, "qmf")
        if(use_alt_shape_prior) name_parts = c(name_parts, "salt")
        if(!is.null(x0_prior_value)) name_parts = c(name_parts, paste0("x", x0_prior_value))
        if(use_mv_baseline) name_parts = c(name_parts, "mvbase")
        if(use_alt_sigmaF_prior) name_parts = c(name_parts, "sfalt")
        if(use_flex_variant) name_parts = c(name_parts, "flex")
        
        sensitivity_name = if(length(name_parts) > 0) paste(name_parts, collapse = "-") else "baseline"
        
        return(list(
            model_type = model_type,
            exec_code = exec_code,
            start_year = start_year,
            cpue_index = cpue_index,
            catch_scenario = catch_scenario,
            sigma_edev = sigma_edev,
            step_scenario = step_scenario,
            use_alt_qeff_prior = use_alt_qeff_prior,
            use_alt_shape_prior = use_alt_shape_prior,
            x0_prior_value = x0_prior_value,
            use_mv_baseline = use_mv_baseline,
            use_alt_sigmaF_prior = use_alt_sigmaF_prior,
            use_alt_qeff_from_meanF = use_alt_qeff_from_meanF,
            use_flex_variant = use_flex_variant,
            sensitivity_name = sensitivity_name
        ))
    }
    
    # Define core sensitivity scenarios
    sensitivity_configs = list(
        # bspm_estF_mvprior sensitivities
        create_sensitivity_config("bspm_estF_mvprior", "F", use_mv_baseline = TRUE),
        create_sensitivity_config("bspm_estF_mvprior", "F", use_alt_sigmaF_prior = TRUE),
        create_sensitivity_config("bspm_estF_mvprior", "F", use_alt_shape_prior = TRUE),
        create_sensitivity_config("bspm_estF_mvprior", "F", use_alt_shape_prior = TRUE,use_alt_sigmaF_prior = TRUE),
        
        # bspm_estq_softdep_mvprior sensitivities
        create_sensitivity_config("bspm_estq_softdep_mvprior", "B", use_alt_qeff_from_meanF = TRUE),
        create_sensitivity_config("bspm_estq_softdep_mvprior", "B", use_alt_shape_prior = TRUE),
        create_sensitivity_config("bspm_estq_softdep_mvprior", "B", use_alt_qeff_from_meanF = TRUE,use_alt_shape_prior = TRUE),

        # bspm_estq_softdep_mvprior_x0 sensitivities (start 1988)
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.7),
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.55),
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.85),
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.7, use_alt_shape_prior = TRUE),
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.7, use_alt_qeff_from_meanF = TRUE),
        create_sensitivity_config("bspm_estq_softdep_mvprior_x0", "BX", start_year = 1988, x0_prior_value = 0.7, use_alt_qeff_from_meanF = TRUE,use_alt_shape_prior = TRUE)
    )

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs and run sensitivity analysis

for(i in 1:length(sensitivity_configs)){
        
        config = sensitivity_configs[[i]]
        
        # Generate run label and number
        run_label_stem = paste0(config$cpue_index, "-exe", config$exec_code, 
                               "-sy", config$start_year, "-cs", config$catch_scenario, 
                               "-e", config$sigma_edev, "-ss", config$step_scenario, 
                               "-", config$sensitivity_name,"_0")
        run_number = 47 + i  # Starting from 0048
        run_number = sprintf("%04d", run_number)
        
        cat("Running model:", run_number, "-", run_label_stem, "\n")
        
        # Determine which CPUE index to fit and set sigmao_input accordingly
        cpue_idx = switch(config$cpue_index,
                         "dwfn" = 1,
                         "au" = 2, 
                         "nz" = 3,
                         "obs" = 4,
                         "obs_no_pf" = 5,
                         "obs_pf_only" = 6,
                         1)  # default to dwfn
        
        sigmao_input_val = mean_se[cpue_idx]
        
        # Set lambdas vector based on which index is being fit
        lambdas_vec = rep(0, n_indices)
        lambdas_vec[cpue_idx] = 1
        
        # Prepare base stan data
        tmp_data = list(
            T = n_years,
            I = n_indices,
            index = index_mat,
            sigmao_mat = se_mat,
            lambdas = lambdas_vec,
            t_dep = which(time_years == 1987),
            use_depletion_prior = 0,
            fit_to_data = 1,
            epsilon = 1e-10,
            effort = catch_effort_annual$effort_scaled,
            n_step = 5,
            n_periods = ceiling((n_years-1) / 5),
            sigma_edev = config$sigma_edev,
            sigmao_input = sigmao_input_val,
            obs_removals = catch_effort_annual$total_catch
        )
        
        # Base prior setup
        tmp.data.priors = list(
            mv_prior_mean = mv_mean_catch,
            mv_prior_sd = mv_prior_sd,
            mv_prior_corr = mv_cor_catch,
            prior_qeff_meanlog = qeff_pars_catch[1],
            prior_qeff_sdlog = qeff_pars_catch[2],
            mv_qdev_prior_mean = mv_qdev_mean,
            mv_qdev_prior_sd = mv_qdev_prior_sd,
            mv_qdev_prior_corr = mv_qdev_cor,
            PriorMean_logsigmap = -2.9311037,
            PriorSD_logsigmap = 0.2661089,
            PriorSD_sigmao_add = 0.2,
            PriorSD_sigmaf = sigmaF_par$x[1],
            prior_depletion_meanlog = dep_pars$x[1],
            prior_depletion_sdlog = dep_pars$x[1]
        )
        
        # Apply sensitivity-specific modifications
        if(config$use_alt_qeff_from_meanF && !is.na(ref_0005_meanF_stats$mean_F_mean)) {
            # Calculate qeff prior from mean F: qeff = mean_F / mean_effort
            mean_effort = mean(catch_effort_annual$effort_scaled, na.rm = TRUE)
            target_qeff = ref_0005_meanF_stats$mean_F_mean / (12*mean_effort) # account for q_dev (e.g., ~exp(2.5))
            tmp.data.priors$prior_qeff_meanlog = log(target_qeff)
            # Scale SD appropriately
            tmp.data.priors$prior_qeff_sdlog = ref_0005_meanF_stats$mean_F_sd / ref_0005_meanF_stats$mean_F_mean
            cat("  Using qeff prior from mean F: target qeff =", target_qeff, "\n")
        }
        
        if(config$use_alt_shape_prior && !is.null(ref_0003_data)) {
            shape_mean = ref_0003_data[name == "logshape" & type == "PriorMean", value]
            shape_sd = ref_0003_data[name == "logshape" & type == "PriorSD", value]
            if(length(shape_mean) > 0) tmp.data.priors$mv_prior_mean[3] = shape_mean
            if(length(shape_sd) > 0) tmp.data.priors$mv_prior_sd[3] = shape_sd
            cat("  Using alternative shape prior: mean =", shape_mean, ", sd =", shape_sd, "\n")
        }
        
        if(!is.null(config$x0_prior_value)) {
            # Switch to 4D prior structure
            tmp.data.priors$mv_prior_mean = mv_mean_catch_4d
            tmp.data.priors$mv_prior_sd = mv_prior_sd_4d
            tmp.data.priors$mv_prior_corr = mv_cor_catch_4d
            tmp.data.priors$mv_prior_mean[4] = log(config$x0_prior_value)
            cat("  Using x0 prior value:", config$x0_prior_value, "\n")
        }
        
        if(config$use_mv_baseline && !is.null(ref_0005_data)) {
            # Extract multivariate prior from 0005 model
            mv_mean_ref = ref_0005_data[name == "mv_prior_mean", value]
            mv_sd_ref = ref_0005_data[name == "mv_prior_sd", value]
            if(length(mv_mean_ref) == 3) tmp.data.priors$mv_prior_mean = mv_mean_ref
            if(length(mv_sd_ref) == 3) tmp.data.priors$mv_prior_sd = mv_sd_ref
            cat("  Using MV baseline priors from 0005 model\n")
        }
        
        if(config$use_alt_sigmaF_prior && !is.na(ref_0025_F_variability$mean_F_sd)) {
            # Use half-normal prior with mean equal to mean SD of Fs from 0025 model
            tmp.data.priors$PriorSD_sigmaf = ref_0025_F_variability$mean_F_sd*sqrt(pi/2)
            cat("  Using sigmaF prior from 0025 F variability: sd =", ref_0025_F_variability$mean_F_sd*sqrt(pi/2), "\n")
        }
        
        # Combine data and priors
        stan.data = c(tmp_data, tmp.data.priors)

        # Apply scenario-specific modifications
        
        # Start year adjustments
        if(config$start_year != 1952){
            year_diff = abs(diff(c(config$start_year, 1952)))
            stan.data$T = stan.data$T - year_diff
            stan.data$index = stan.data$index[-seq(1, year_diff), , drop=FALSE]
            stan.data$sigmao_mat = stan.data$sigmao_mat[-seq(1, year_diff), , drop=FALSE]
            stan.data$effort = stan.data$effort[-seq(1, year_diff)]
            stan.data$obs_removals = stan.data$obs_removals[-seq(1, year_diff)]
            stan.data$t_dep = stan.data$t_dep - year_diff
        }

        # Step scenario adjustments
        if(config$step_scenario == "5reg" & config$exec_code %in% c("B","BX")){
            stan.data$n_step = 5
            stan.data$n_periods = ceiling((stan.data$T-1) / stan.data$n_step)
        } else if(config$step_scenario == "5reg" & config$exec_code %in% c("BF","BFX")){
            stan.data$n_periods = ceiling((stan.data$T-1) / 5)
            stan.data$q_period_index = pmin(ceiling(1:stan.data$T / 5), stan.data$n_periods)
        }
        
        # Catch scenario
        if(config$catch_scenario == "0.2flat" & config$exec_code %in% c("B","BX","F")){
            stan.data$sigmac = 0.2
        } else if(config$catch_scenario == "0.2flat" & config$exec_code %in% c("BF","BFX")){
            stan.data$sigmac = rep(0.2, stan.data$T)
        }

        # Select appropriate stan model
        stan_c = switch(config$exec_code,
                       "B" = stan_c.list[[1]],    # bspm_estq_softdep_mvprior
                       "BF" = stan_c.list[[2]],   # bspm_estq_flex
                       "BX" = stan_c.list[[3]],   # bspm_estq_softdep_mvprior_x0
                       "BFX" = stan_c.list[[4]],  # bspm_estq_flex_x0
                       "F" = stan_c.list[[5]])    # bspm_estF_mvprior

        # Run model
        fit = fit_rstan(stan.data,
                        stan_c,
                        run_label = paste0(run_number,"-",run_label_stem),
                        exec_name = config$model_type,
                        seed = 321,
                        chains = 5,
                        n_thin = 10,
                        iter_keep = 200,
                        burnin.prop = 0.5,
                        adapt_delta = 0.99,
                        max_treedepth = 12,
                        silent = FALSE,
                        stan_save_dir = file.path(proj_dir,"data","output","model_runs"),
                        n_cores = 5)

        # Print key results
        key_pars = c("logK", "r", "shape", "sigmap")
        
        # Add x0 parameter for models that estimate it
        if(!is.null(config$x0_prior_value)) {
            key_pars = c(key_pars, "x0")
        }
        
        # Add sigmaF for F-estimation models
        if(config$model_type == "bspm_estF_mvprior") {
            key_pars = c(key_pars, "sigmaf")
        }
        
        print(fit, pars = key_pars)
        
        quick_diagnostics(fit)
        
        # Save sensitivity-specific diagnostics
        t = as.data.table(summary(fit)$summary)
        t$name = rownames(summary(fit)$summary)
        convergence_issues = t[n_eff < 500 | Rhat > 1.01]
        if(nrow(convergence_issues) > 0) {
            cat("Convergence issues detected for run", run_number, ":\n")
            print(convergence_issues[, .(name, mean, n_eff, Rhat)])
        }
        
        # Print sensitivity-specific summaries
        if(config$use_alt_qeff_from_meanF && "qeff" %in% t$name) {
            qeff_summary = t[name == "qeff", .(mean, sd)]
            cat("qeff posterior: mean =", qeff_summary$mean, ", sd =", qeff_summary$sd, "\n")
        }
        
        if(!is.null(config$x0_prior_value) && "x0" %in% t$name) {
            x0_summary = t[name == "x0", .(mean, sd)]
            cat("x0 posterior: mean =", x0_summary$mean, ", sd =", x0_summary$sd, "\n")
            cat("x0 depletion: mean =", exp(x0_summary$mean), "\n")
        }
        
        cat("Completed run", run_number, "-", run_label_stem, "\n\n")
}

