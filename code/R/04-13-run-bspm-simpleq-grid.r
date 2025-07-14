# Nicholas Ducharme-Barth
# 2025/07/14
# Run BSPM with effort-based fishing mortality and updated prior structure
# Model grid

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
    exec_name_vec = c("bspm_estqsimple_softdep_fullmvprior_x0_sttgamma_flexsigmaC_OPT")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
    model_config_df = rbind(
                             expand.grid(exec=c("oFSTTGF"),
                                  cpue=c("dwfn","au","nz","obs","obsNoPF"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b","alt")) 
                            
    )

    model_config_df = unique(model_config_df)

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
            run_number = 99 + i
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
        } else if (model_config_df$exec[i] %in% c("FSTTGF","oFSTTGF")){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025
        } else if (model_config_df$exec[i] == "qV2FSTTGF"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean[-6])
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd[-6])
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr[-6,-6])

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025
        }

                stan_c = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = stan_c.list[[1]])    # bspm_estqsimple_softdep_mvprior_x0_stt
                
                exec_name = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = exec_name_vec[1])

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

        

        # calc realized prior
        stan.data.ppc = stan.data
        stan.data.ppc$fit_to_data = 0L
        stan.inits.ppc = replicate(5, stan_inits_func(Tm1 = (stan.data$T-1), 
                                                    n_periods = stan.data$n_periods,
                                                    exec_name = exec_name), 
                                simplify=FALSE)
    
        options(mc.cores = 5)
        ppc = sampling(object=stan_c,
                    data = stan.data.ppc,
                    init = stan.inits.ppc,
                    chains = 5,
                    warmup = 250,
                    iter = 7500,
                    thin = 1,
                    seed = 321,
                    control = list(adapt_delta = 0.99,max_treedepth=12))


        # summarize hmc samples
            hmc_samples_ppc = as.data.table(ppc) %>%
                        .[,iter:=1:.N] %>%
                        melt(.,id.vars="iter") %>%
                        .[,.(iter,variable,value)]

            names_dt = data.table(variable=unique(hmc_samples_ppc$variable)) %>%
                            .[,split:=sapply(variable,function(x)ifelse(length(grep("[",x,fixed=TRUE))==0,0,1))] %>%
                            .[,split_mat:=sapply(variable,function(x)ifelse(length(grep(",",x,fixed=TRUE))==0,0,1))] %>%
                            .[,row:=as.numeric(NA)] %>%
                            .[,col:=as.numeric(NA)] %>%
                            .[split_mat==1,row:=sapply(variable,function(x)as.numeric(strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][1]))] %>%
                            .[split_mat==1,col:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][2])))] %>%
                            .[split==1&split_mat==0,row:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(as.character(x),"[",fixed=TRUE)[[1]][2])))] %>%
                            .[,name:=variable] %>%
                            .[split==1,name:=sapply(variable,function(x)strsplit(as.character(x),"[",fixed=TRUE)[[1]][1])] %>%
                            .[,.(variable,name,row,col)]
                

            hmc_samples_ppc = merge(hmc_samples_ppc,names_dt,by="variable") %>%
                            .[!(name%in%c("C","sum1","sum2","p","sigmao2","dev","epsp"))] %>%
                            .[,.(iter,variable,name,row,col,value)]
            dir.create(file.path(proj_dir,"data","output","model_runs",paste0(run_number,"-",run_label_stem)),recursive=TRUE)
            fwrite(hmc_samples_ppc,file=file.path(proj_dir,"data","output","model_runs",paste0(run_number,"-",run_label_stem),"hmc_samples_ppc.csv"))        


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
            print(fit, pars = c("logK", "r", "shape","x0", "sigmap", "sigmao_add", "qeff", "sigma_qdev","nu_catch", 
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))
        } else {
           print(fit, pars = c("logK", "r", "shape","x0", "sigmap", "sigmao_add", "qeff","rho", "sigma_qdev", "nu_catch",
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]")) 
        }
        print(stan.data$obs_removals[c(3,70)])
        
        quick_diagnostics(fit)

        compare_marginals(ppc,fit,c("logK", "r", "shape","x0", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev","nu_catch",
                           "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]"))

        t=as.data.table(summary(fit)$summary)
        t$name = rownames(summary(fit)$summary)
        na.omit(t)[n_eff<500|Rhat>1.01]

            for(j in 1:5){
                retro_run_label_stem = paste0(model_config_df$cpue[i],
                            "-exe",model_config_df$exec[i],
                            "-c",model_config_df$sigma_catch[i],
                            "-n",model_config_df$shape[i],
                            "-q",model_config_df$qeff[i],
                            "-s",model_config_df$n_step[i],
                            "-o52",model_config_df$obs1952[i],
                            "-o54",model_config_df$obs1954[i],"_",j)
                fit = fit_rstan(stan.data,
                                    stan_c,
                                    run_label = paste0(run_number,"-",retro_run_label_stem),
                                    exec_name = exec_name,
                                    seed  = 321,
                                    chains = 5,
                                    n_thin = 10,
                                    iter_keep = 200,
                                    burnin.prop = 0.5,
                                    adapt_delta = 0.99,
                                    max_treedepth = 12,
                                    silent = FALSE,
                                    stan_save_dir=file.path(proj_dir,"data","output","model_runs"),
                                    n_cores=5)
            }
    }
    
    
   