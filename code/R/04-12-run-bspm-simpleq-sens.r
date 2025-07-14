# Nicholas Ducharme-Barth
# 2025/07/10
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
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0_newlogK","updated_stan_data.RData"))
    newlogK_stan_data = updated_stan_data
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0_refined","updated_stan_data.RData"))
    refined_stan_data = updated_stan_data
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0","updated_stan_data.RData"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load additional cpue indices
    au_cpue_dt = fread(file.path(proj_dir,"data","input","AU_cpue.csv"))
    nz_cpue_dt = fread(file.path(proj_dir,"data","input","NZ_cpue.csv"))
    obs_cpue_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP.csv"))
    obs_cpue_no_PF_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_no_PF.csv"))
    obs_cpue_PF_only_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_PF_only.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# calc mean_se   
    mean_se = rep(NA,6)

    # dwfn cpue
    mean_se[1] = newlogK_stan_data$sigmao_input 
    # australia cpue
    mean_se[2] = mean(au_cpue_dt$se_log,na.rm=TRUE)
    # new zealand cpue
    mean_se[3] = mean(nz_cpue_dt$se_log,na.rm=TRUE)
    # Observer CPUE
    # calc SE from quantiles
    obs_cpue_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[4] = mean(obs_cpue_dt$se_from_quantiles,na.rm=TRUE)
    # calc SE from quantiles
    obs_cpue_no_PF_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[5] = mean(obs_cpue_no_PF_dt$se_from_quantiles,na.rm=TRUE)
    # calc SE from quantiles
    obs_cpue_PF_only_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[6] = mean(obs_cpue_PF_only_dt$se_from_quantiles,na.rm=TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name_vec = c("bspm_estqsimple_softdep_mvprior_x0_stt","bspm_estqsimple_softdep_mvprior_x0_sttgamma_flexsigmaC","bspm_estqsimple_softdep_fullmvprior_x0_sttgamma_flexsigmaC","bspm_estqsimpleV2_softdep_fullmvprior_x0_sttgamma_flexsigmaC")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
    model_config_df = rbind(expand.grid(exec=c("FSTTGF"),
                                  cpue=c("au","nz","obs","obsNoPF","obsPFonly"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b")), 
                            expand.grid(exec=c("FSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("p0.2","f0.1","p0.1"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b")),
                            expand.grid(exec=c("FSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b","ac","ae"),
                                  obs1954 = c("b","ac"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b")),
                            expand.grid(exec=c("FSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("b","r","newK"),
                                  shape=c("b","alt")),
                            expand.grid(exec=c("FSTTGF"),
                                  cpue=c("allnoObs","allobs","allnoPF","allPFonly"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b")),
                            expand.grid(exec=c("FSTTGF"),
                                  cpue=c("dwfn","au","obsNoPF"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newKX"),
                                  shape=c("b","alt")) 
                            
    )

    model_config_df = unique(model_config_df)
    # remove model 69
    model_config_df = subset(model_config_df,!(exec=="FSTTGF"&cpue=="dwfn"&sigma_catch=="f0.2"&obs1952=="b"&obs1954=="b"&n_step=="1"&qeff=="newK"&shape=="b"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

    for(i in 1:nrow(model_config_df)){
            run_label_stem = paste0(model_config_df$cpue[i],
                            "-exe",model_config_df$exec[i],
                            "-c",model_config_df$sigma_catch[i],
                            "-n",model_config_df$shape[i],
                            "-q",model_config_df$qeff[i],
                            "-s",model_config_df$n_step[i],
                            "-o52",model_config_df$obs1952[i],
                            "-o54",model_config_df$obs1954[i],"_0")
            run_number = 70 + i
            run_number = sprintf("%04d", run_number)
        if(model_config_df$qeff[i] == "b"){
            stan.data = updated_stan_data

            stan.data$full_mv_prior_mean = c(updated_stan_data$mv_prior_mean,updated_stan_data$mv_qdev_prior_mean)
            stan.data$full_mv_prior_sd = c(updated_stan_data$mv_prior_sd,updated_stan_data$mv_qdev_prior_sd)
            stan.data$full_mv_prior_corr = cbind(rbind(updated_stan_data$mv_prior_corr,matrix(rep(0,10),nrow=2)),
                                                        rbind(matrix(rep(0,10),nrow=5),updated_stan_data$mv_qdev_prior_corr))

        } else if(model_config_df$qeff[i] == "r"){
            stan.data = refined_stan_data

            stan.data$full_mv_prior_mean = c(refined_stan_data$mv_prior_mean,refined_stan_data$mv_qdev_prior_mean)
            stan.data$full_mv_prior_sd = c(refined_stan_data$mv_prior_sd,refined_stan_data$mv_qdev_prior_sd)
            stan.data$full_mv_prior_corr = cbind(rbind(refined_stan_data$mv_prior_corr,matrix(rep(0,10),nrow=2)),
                                                        rbind(matrix(rep(0,10),nrow=5),refined_stan_data$mv_qdev_prior_corr))
        } else if(model_config_df$qeff[i] == "newK"){
            stan.data = newlogK_stan_data
        } else if(model_config_df$qeff[i] == "newKX"){
            stan.data = newlogK_stan_data
            stan.data$full_mv_prior_mean[c(1,5)] = c(14.32,-1.32)
            stan.data$full_mv_prior_sd[c(1,5)] = c(0.2,0.2)
        }             
        
        stan.data$sigmac = model_config_df$sigma_catch[i]
        stan.data$n_step = model_config_df$n_step[i]  # years per period for catchability
        stan.data$n_periods = ceiling((stan.data$T-1) / stan.data$n_step)
        stan.data$fit_to_data = 1L
        if(model_config_df$exec[i] == "STT"){
            stan.data$nu_catch_rate = 0.1
        } else if (model_config_df$exec[i] == "STTGF"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
        } else if (model_config_df$exec[i] == "FSTTGF"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)
        } else if (model_config_df$exec[i] == "qV2FSTTGF"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean[-6])
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd[-6])
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr[-6,-6])
        }

                stan_c = switch(as.character(model_config_df$exec[i]),
                       "STT" = stan_c.list[[1]],
                       "STTGF" = stan_c.list[[2]],
                       "FSTTGF" = stan_c.list[[3]],
                       "qV2FSTTGF" = stan_c.list[[4]])    # bspm_estqsimple_softdep_mvprior_x0_stt
                
                exec_name = switch(as.character(model_config_df$exec[i]),
                       "STT" = exec_name_vec[1],
                       "STTGF" = exec_name_vec[2],
                       "FSTTGF" = exec_name_vec[3],
                       "qV2FSTTGF" = exec_name_vec[4])

        # Determine which CPUE index to fit and set sigmao_input accordingly

        if(model_config_df$cpue[i] %in% c("dwfn","au","nz","obs","obsNoPF","obsPFonly")){
            cpue_idx = switch(as.character(model_config_df$cpue[i]),
                            "dwfn" = 1,
                            "au" = 2, 
                            "nz" = 3,
                            "obs" = 4,
                            "obsNoPF" = 5,
                            "obsPFonly" = 6,
                            1)  # default to dwfn
                    
            # Set lambdas vector based on which index is being fit
            lambdas_vec = rep(0, 6)
            lambdas_vec[cpue_idx] = 1
            stan.data$lambdas = lambdas_vec
            stan.data$sigmao_input = mean_se[cpue_idx]
        } else {
            if(model_config_df$cpue[i] == "allnoObs"){
                stan.data$lambdas = c(1,1,1,0,0,0)
                stan.data$sigmao_input = 0.15
            } else if(model_config_df$cpue[i] == "allobs"){
                stan.data$lambdas = c(1,1,1,1,0,0)
                stan.data$sigmao_input = 0.15
            } else if(model_config_df$cpue[i] == "allnoPF"){
                stan.data$lambdas = c(1,1,1,0,1,0)
                stan.data$sigmao_input = 0.15
            } else if(model_config_df$cpue[i] == "allPFonly"){
                stan.data$lambdas = c(1,1,1,0,0,1)
                stan.data$sigmao_input = 0.15
            }
        }

        # catch scenarios
            generate_power_decline = function(N, start_power = 0.5, end_power = 0.2) {
                years = 0:N
                decline_rate = (end_power / start_power)^(1/N)
                power = start_power * decline_rate^years
                return(power)
            }
        
        # catch uncertainty scenario
            if(model_config_df$sigma_catch[i] == "f0.2"){
                stan.data$sigmac = rep(0.2,stan.data$T)
            } else if(model_config_df$sigma_catch[i] == "f0.1"){
                stan.data$sigmac = rep(0.1,stan.data$T)
            } else if(model_config_df$sigma_catch[i] == "p0.2"){
                stan.data$sigmac = generate_power_decline(stan.data$T-1, start_power = 0.5, end_power = 0.2)
            } else if(model_config_df$sigma_catch[i] == "p0.1"){
                stan.data$sigmac = generate_power_decline(stan.data$T-1, start_power = 0.5, end_power = 0.1)
            }

        # data uncertainties
            if(model_config_df$obs1952[i] == "ac"){
               stan.data$obs_removals[1] = stan.data$effort[1]*mean((stan.data$obs_removals/stan.data$effort)[2:6])
            } else if(model_config_df$obs1952[i] == "ae"){
                stan.data$effort[1] = stan.data$obs_removals[1]/mean((stan.data$obs_removals/stan.data$effort)[2:6])
            }

            if(model_config_df$obs1954[i] == "ac"){
                stan.data$obs_removals[3] = 0.5 * stan.data$obs_removals[3]
            }
        
        # alt shape
            if(model_config_df$shape[i] == "alt"){
                stan.data$mv_prior_mean[3] = log(2)
            }
        


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

        if(exec_name == "qV2FSTTGF"){
            print(fit, pars = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "sigma_qdev", 
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))
        } else {
           print(fit, pars = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff","rho", "sigma_qdev", 
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]")) 
        }
        print(fit, pars = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "sigma_qdev", 
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))
        print(stan.data$obs_removals[c(3,70)])
        
        quick_diagnostics(fit)
        # compare_marginals(ppc,fit,c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev", 
        #                    "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))

        t=as.data.table(summary(fit)$summary)
        t$name = rownames(summary(fit)$summary)
        na.omit(t)[n_eff<500|Rhat>1.01]
    }
    
    
   