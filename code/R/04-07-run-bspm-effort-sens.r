# Nicholas Ducharme-Barth
# 2025/06/19
# Run BSPM with effort-based fishing mortality and updated prior structure

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
    cpue_dt = fread(file.path(proj_dir,"data","input","cpue.csv"))
    
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
    
    # Load updated prior parameters from pushforward analysis
    mv_mean_catch = read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_mean_catch.csv"))$x
    mv_cov_catch = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_cov_catch.csv"))[,-1])
    mv_cor_catch = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_cor_catch.csv"))[,-1])
    
    # Independent qeff prior
    qeff_pars_catch = read.csv(file.path(proj_dir,"data","output","pushforward","q","qeff_pars_catch.csv"))$x
    
    # Bivariate rho/sigma_qdev prior
    mv_qdev_mean = read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_mean.csv"))$x
    mv_qdev_cov = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_cov.csv"))[,-1])
    mv_qdev_cor = as.matrix(read.csv(file.path(proj_dir,"data","output","pushforward","q","mv_qdev_cor.csv"))[,-1])
    
    # Other priors
    dep_pars = read.csv(file.path(proj_dir,"data","output","rel_dep_ssb_pars.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load additional cpue indices
    au_cpue_dt = fread(file.path(proj_dir,"data","input","AU_cpue.csv"))
    nz_cpue_dt = fread(file.path(proj_dir,"data","input","NZ_cpue.csv"))
    obs_cpue_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP.csv"))
    obs_cpue_no_PF_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_no_PF.csv"))
    obs_cpue_PF_only_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_PF_only.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare data matrices
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    # Merge with effort data
    catch_effort_annual = merge(catch_annual, effort_dt, by.x = "year", by.y = "time", all.x = TRUE)
    
    # Fill missing effort with mean (if any)
    if(any(is.na(catch_effort_annual$effort_scaled))) {
        catch_effort_annual[is.na(effort_scaled), effort_scaled := mean(effort_scaled, na.rm = TRUE)]
    }
    
    time_years = catch_effort_annual$year
    n_years = length(time_years)
    
    index_mat = matrix(-999, nrow = n_years, ncol = 6)
    se_mat = matrix(-999, nrow = n_years, ncol = 6)
    mean_se = rep(NA,6)

    # dwfn cpue
    cpue_years = floor(cpue_dt$Time)
    mean_se[1] = mean(cpue_dt$SE,na.rm=TRUE)
    for(i in 1:nrow(cpue_dt)) {
        year_idx = which(time_years == cpue_years[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 1] = cpue_dt$Obs[i]
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
# prepare multivariate prior parameters
    # Extract standard deviations from covariance matrix (3D now)
    mv_prior_sd = sqrt(diag(mv_cov_catch))
    
    # Extract standard deviations for qdev parameters (2D)
    mv_qdev_prior_sd = sqrt(diag(mv_qdev_cov))

    # Validate correlation matrices are positive definite
    if(min(eigen(mv_cor_catch)$values) < 1e-10) {
        warning("3D correlation matrix has very small eigenvalues, adding regularization")
        mv_cor_catch = mv_cor_catch + diag(1e-8, 3)
    }
    
    if(min(eigen(mv_qdev_cor)$values) < 1e-10) {
        warning("2D correlation matrix has very small eigenvalues, adding regularization")  
        mv_qdev_cor = mv_qdev_cor + diag(1e-8, 2)
    }
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name = "bspm_estq_softdep_mvprior" # bspm_estq_softdep_mvprior # bspm_estq_optimized
    stan_c = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name,".stan")), model_name = exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
    model_config_df = rbind( expand.grid(cpue=c("au","nz","obs","obsNoPF","obsPFonly"),
                                  sigma_catch = c(0.2),
                                  sigma_edev = c(0.3),
                                  n_step=5),
                            expand.grid(cpue=c("dwfn"),
                                  sigma_catch = c(0.1,0.2,0.4),
                                  sigma_edev = c(0.3,0.5,1,2),
                                  n_step=c(2,5))
    )

    model_config_df = unique(model_config_df)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

for(i in 1:nrow(model_config_df)){
        run_label_stem = paste0(model_config_df$cpue[i],"-c",model_config_df$sigma_catch[i],"-e",model_config_df$sigma_edev[i],"-s",model_config_df$n_step[i],"_0")
        run_number = 6 + i
        if(run_number<10){
            run_number = paste0("000",run_number)
        } else {
            run_number = paste0("00",run_number)
        }

        if(model_config_df$cpue[i]=="dwfn"){
            lambda_vec = c(1,0,0,0,0,0)
            sigmao_input = mean_se[1]
        } else if(model_config_df$cpue[i]=="au"){
            lambda_vec = c(0,1,0,0,0,0)
            sigmao_input = mean_se[2]
        } else if(model_config_df$cpue[i]=="nz"){
            lambda_vec = c(0,0,1,0,0,0)
            sigmao_input = mean_se[3]
        } else if(model_config_df$cpue[i]=="obs"){
            lambda_vec = c(0,0,0,1,0,0)
            sigmao_input = mean_se[4]
        } else if(model_config_df$cpue[i]=="obsNoPF"){
            lambda_vec = c(0,0,0,0,1,0)
            sigmao_input = mean_se[5]
        } else {
            lambda_vec = c(0,0,0,0,0,1)
            sigmao_input = mean_se[6]
        }

 # Define effort parameters
    n_step = model_config_df$n_step[i]  # years per period for catchability
    n_periods = ceiling((n_years-1) / n_step)
            
    # specify data and priors
    tmp.data.priors = list(
        # 3D multivariate prior (logK, log_r, log_shape)
        mv_prior_mean = mv_mean_catch,
        mv_prior_sd = mv_prior_sd,
        mv_prior_corr = mv_cor_catch,
        
        # Independent qeff prior
        prior_qeff_meanlog = qeff_pars_catch[1],
        prior_qeff_sdlog = qeff_pars_catch[2],
        
        # Bivariate rho/sigma_qdev prior  
        mv_qdev_prior_mean = mv_qdev_mean,
        mv_qdev_prior_sd = mv_qdev_prior_sd,
        mv_qdev_prior_corr = mv_qdev_cor,
        
        # Other priors
        PriorMean_logsigmap = -2.9311037,
        PriorSD_logsigmap = 0.2661089,
        PriorSD_sigmao_add = 0.2,
        prior_depletion_meanlog = dep_pars$x[1],
        prior_depletion_sdlog = dep_pars$x[2]
    )

    tmp_data = list(
        T = as.integer(n_years),
        I = as.integer(ncol(index_mat)),
        index = index_mat,
        sigmao_mat = se_mat,
        lambdas = lambda_vec,
        t_dep = 37, # model starts t=1 in 1952, 37 is 1988
        use_depletion_prior = 0L, # Set to 1L to use prior, 0L to skip it
        fit_to_data = 1L, # Set to 1L to fit to data, 0L to skip likelihoods
        epsilon = 1e-10,
        
        # Effort-based parameters
        effort = catch_effort_annual$effort_scaled,
        n_step = as.integer(n_step),
        n_periods = as.integer(n_periods),
        sigma_edev = model_config_df$sigma_edev[i],  # Prior SD for effort deviations
        
        # Observation error parameters
        sigmao_input = sigmao_input,
        obs_removals = catch_effort_annual$total_catch,
        sigmac = model_config_df$sigma_catch[i]
    )
                     
    stan.data = c(tmp_data, tmp.data.priors)
    # stan.data.ppc = stan.data
    # stan.data.ppc$fit_to_data = 0L

    # ppc = fit_rstan(stan.data.ppc,
    #                 stan_c,
    #                 run_label = paste0(run_number,"-",run_label_stem,"-ppc"),
    #                 exec_name = exec_name,
    #                 seed = 321,
    #                 chains = 5,
    #                 n_thin = 10,
    #                 iter_keep = 1000,
    #                 burnin.prop = 0.5,
    #                 adapt_delta = 0.99,
    #                 max_treedepth = 12,
    #                 silent = FALSE,
    #                 stan_save_dir = file.path(proj_dir,"data","output","model_runs"),
    #                 n_cores = 5)

    fit = fit_rstan(stan.data,
                    stan_c,
                    run_label = paste0(run_number,"-",run_label_stem),
                    exec_name = exec_name,
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

    print(fit, pars = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev", 
                       "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))
    print(stan.data$obs_removals[c(3,70)])
    
    quick_diagnostics(fit)
    # compare_marginals(ppc,fit,c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev", 
    #                    "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))

    t=as.data.table(summary(fit)$summary)
    t$name = rownames(summary(fit)$summary)
    t[n_eff<500|Rhat>1.01]
}
    
    
   