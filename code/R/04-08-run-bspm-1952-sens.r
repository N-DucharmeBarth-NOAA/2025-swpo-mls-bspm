# Nicholas Ducharme-Barth
# 2025/06/30
# Run alternative models to test sensitivity to 1952 catch
# Given the observed 1952 effort either:
# 1) Observed catch is wrong, should be ~10x higher
# OR
# 2) Catchability in 1952 is uniquely different (e.g., first year of fishery so fished differently from technical or spatial standpoint)

# This appears to be having substantial impact on model estimated initial conditions in terms of population scale, initial depletion, and process error

# How to address this?
# 1) Down-weight catch (super-high 1952 obs error); let the model estimate what catch should be given effort
# 2) Externally correct catch based on observed effort (test sensitivity to this correction)
# 3) Begin the model in 1953 and estimate model initial conditions (x0)

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

    # add x0 to the mix
    mv_cor_catch_4d = matrix(c(
    1.0,   mv_cor_catch[1,2], mv_cor_catch[1,3], 0.0,    # logK row
    mv_cor_catch[1,2], 1.0,   mv_cor_catch[2,3], 0.0,    # log_r row  
    mv_cor_catch[1,3], mv_cor_catch[2,3], 1.0,   0.0,    # log_shape row
        0.0,     0.0,    0.0,   1.0      # log_x0 row (uncorrelated)
    ), nrow=4, ncol=4)

    # Where corr_12, corr_13, corr_23 are your existing 3D correlations
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

    if(min(eigen(mv_cor_catch_4d)$values) < 1e-10) {
        warning("4D correlation matrix has very small eigenvalues, adding regularization")
        mv_cor_catch_4d = mv_cor_catch_4d + diag(1e-8, 4)
    }
    
    if(min(eigen(mv_qdev_cor)$values) < 1e-10) {
        warning("2D correlation matrix has very small eigenvalues, adding regularization")  
        mv_qdev_cor = mv_qdev_cor + diag(1e-8, 2)
    }
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name_vec = c("bspm_estq_softdep_mvprior","bspm_estq_flex","bspm_estq_softdep_mvprior_x0","bspm_estq_flex_x0")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
# baseline model: 1952, dwfn, catch 0.2, effort 0.3, n_step 5


    model_config_df_base = data.frame(exec = c("B","BF","BX","BFX"),
                                 start_year = c(1952,1952,1952,1952),
                                 cpue = c("dwfn","dwfn","dwfn","dwfn"),
                                 catch_scenario = rep("0.2flat",4),
                                 sigma_edev = rep(0.3,4),
                                 step_scenario = rep("5reg",4))

    model_config_df_sens = data.frame(exec = c("BF","BF","B","BX","BX","BX","BX","BX"),
                                 start_year = c(1952,1952,1952,1953,1952,1952,1952,1988),
                                 cpue = c("dwfn","dwfn","dwfn","dwfn","dwfn","dwfn","dwfn","dwfn"),
                                 catch_scenario = c("0.2DWflat","0.2DWpower","0.2correct","0.2flat","0.2DWNoProc","0.2DWNoProcOldlKPrior","0.2DWNoProcOldlKPrior","0.2flat"),
                                 sigma_edev = c(rep(0.3,4),0.5,0.5,0.3,0.3),
                                 step_scenario = rep("5reg",8))

    model_config_df = rbind(model_config_df_base,model_config_df_sens)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

for(i in 1:nrow(model_config_df)){
        run_label_stem = paste0(model_config_df$cpue[i],"-exe",model_config_df$exec[i],"-sy",model_config_df$start_year[i],"-cs",model_config_df$catch_scenario[i],"-e",model_config_df$sigma_edev[i],"-ss",model_config_df$step_scenario[i],"_0")
        run_number = 35 + i
        run_number = sprintf("%04d", run_number)

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

    # Define default data
        tmp.data.priors = list(
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
            sigma_edev = model_config_df$sigma_edev[i],  # Prior SD for effort deviations
            
            # Observation error parameters
            sigmao_input = sigmao_input,
            obs_removals = catch_effort_annual$total_catch
        )
                        
        stan.data = c(tmp_data, tmp.data.priors)

        # Define scenario specific data parameters

        # start year
            if(model_config_df$start_year[i]!=1952){
                stan.data$T = stan.data$T - abs(diff(c(model_config_df$start_year[i],1952)))
                stan.data$index = stan.data$index[-seq(from=1,to=abs(diff(c(model_config_df$start_year[i],1952))),by=1),]
                stan.data$sigmao_mat = stan.data$sigmao_mat[-seq(from=1,to=abs(diff(c(model_config_df$start_year[i],1952))),by=1),]
                stan.data$effort = stan.data$effort[-seq(from=1,to=abs(diff(c(model_config_df$start_year[i],1952))),by=1)]
                stan.data$obs_removals = stan.data$obs_removals[-seq(from=1,to=abs(diff(c(model_config_df$start_year[i],1952))),by=1)]
                stan.data$t_dep = stan.data$t_dep - abs(diff(c(model_config_df$start_year[i],1952)))
            }

        # step scenario
            if(model_config_df$step_scenario[i]=="5reg" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$n_step = 5
                stan.data$n_periods = ceiling((stan.data$T-1) / stan.data$n_step)
            } else if(model_config_df$step_scenario[i]=="5reg" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$n_periods = ceiling((stan.data$T-1) / 5)
                stan.data$q_period_index = pmin(ceiling(1:stan.data$T / 5), stan.data$n_periods)
            }
        
        # multivariate priors
            if(model_config_df$exec[i] %in% c("BX","BFX")){
                stan.data$mv_prior_mean = mv_mean_catch_4d
                stan.data$mv_prior_sd = mv_prior_sd_4d
                stan.data$mv_prior_corr = mv_cor_catch_4d
            } else if(model_config_df$exec[i] %in% c("B","BF")){
                stan.data$mv_prior_mean = mv_mean_catch
                stan.data$mv_prior_sd = mv_prior_sd
                stan.data$mv_prior_corr = mv_cor_catch
            }

            if(model_config_df$start_year[i]==1988){
                stan.data$mv_prior_mean[4] = log(0.7)
                stan.data$mv_prior_sd[4] = 0.2
            }
        
        # catch scenario
            if(model_config_df$catch_scenario[i]=="0.2flat" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
            } else if(model_config_df$catch_scenario[i]=="0.2flat" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = rep(0.2,stan.data$T)
            }

            if(model_config_df$catch_scenario[i]=="0.2DWNoProc" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
                stan.data$PriorMean_logsigmap = log(0.001)
                stan.data$PriorSD_logsigmap = 0.05
            } else if(model_config_df$catch_scenario[i]=="0.2DWNoProc" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(2,rep(0.2,stan.data$T-1))
                stan.data$PriorMean_logsigmap = log(0.001)
                stan.data$PriorSD_logsigmap = 0.05
            }

            if(model_config_df$catch_scenario[i]=="0.2DWNoProcOldlKPrior" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
                stan.data$PriorMean_logsigmap = log(0.001)
                stan.data$PriorSD_logsigmap = 0.05
                stan.data$mv_prior_mean[1] = 13.8633826342304
                stan.data$mv_prior_sd[1] = 0.456389088374683
            } else if(model_config_df$catch_scenario[i]=="0.2DWNoProcOldlKPrior" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(2,rep(0.2,stan.data$T-1))
                stan.data$PriorMean_logsigmap = log(0.001)
                stan.data$PriorSD_logsigmap = 0.05
                stan.data$mv_prior_mean[1] = 13.8633826342304
                stan.data$mv_prior_sd[1] = 0.456389088374683
            }

            if(model_config_df$catch_scenario[i]=="0.2DWflat" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
            } else if(model_config_df$catch_scenario[i]=="0.2DWflat" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(10,rep(0.2,stan.data$T-1))
            }

            if(model_config_df$catch_scenario[i]=="0.2DWv2flat" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
            } else if(model_config_df$catch_scenario[i]=="0.2DWv2flat" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(100,rep(0.2,stan.data$T-1))
            }

            if(model_config_df$catch_scenario[i]=="0.2DWv3flat" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
            } else if(model_config_df$catch_scenario[i]=="0.2DWv3flat" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(5,rep(0.2,stan.data$T-1))
            }

            generate_power_decline = function(N, start_power = 0.5, end_power = 0.2) {
                years = 0:N
                decline_rate = (end_power / start_power)^(1/N)
                power = start_power * decline_rate^years
                return(power)
            }

            if(model_config_df$catch_scenario[i]=="0.2DWpower" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
            } else if(model_config_df$catch_scenario[i]=="0.2DWpower" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = c(10,generate_power_decline(stan.data$T-2))
            }

            if(model_config_df$catch_scenario[i]=="0.2correct" & model_config_df$exec[i] %in% c("B","BX")){
                stan.data$sigmac = 0.2
                stan.data$obs_removals[1] =  mean(stan.data$obs_removals[c(2:6)]/stan.data$effort[c(2:6)])*stan.data$effort[1]
            } else if(model_config_df$catch_scenario[i]=="0.2correct" & model_config_df$exec[i] %in% c("BF","BFX")){
                stan.data$sigmac = rep(0.2,stan.data$T)
                stan.data$obs_removals[1] =  mean(stan.data$obs_removals[c(2:6)]/stan.data$effort[c(2:6)])*stan.data$effort[1]
            }

        # exec scenario
            if(model_config_df$exec[i] =="B"){
                stan_c = stan_c.list[[1]]
            } else if(model_config_df$exec[i] =="BF"){
                stan_c = stan_c.list[[2]] 
            } else if(model_config_df$exec[i] =="BX"){
                stan_c = stan_c.list[[3]] 
            } else if(model_config_df$exec[i] =="BFX"){
                stan_c = stan_c.list[[4]] 
            }

    fit = fit_rstan(stan.data,
                    stan_c,
                    run_label = paste0(run_number,"-",run_label_stem),
                    exec_name = exec_name_vec[which(c("B","BF","BX","BFX") == model_config_df$exec[i])],
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

    print(fit, pars = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev"))
    
    quick_diagnostics(fit)
    # compare_marginals(ppc,fit,c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev", 
    #                    "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))

    t=as.data.table(summary(fit)$summary)
    t$name = rownames(summary(fit)$summary)
    t[n_eff<500|Rhat>1.01]
}
    
    
   