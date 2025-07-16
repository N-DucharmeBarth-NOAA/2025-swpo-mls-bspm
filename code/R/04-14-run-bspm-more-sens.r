# Nicholas Ducharme-Barth
# 2025/07/15
# Run BSPM with effort-based fishing mortality and updated prior structure
# More one-offs

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
# calculate a reference F
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
    ref_0100_F_variability = calc_F("0100-dwfn-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
    model_config_df = rbind(
                             expand.grid(exec=c("oFLNF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "b",
                                  sy = "1952.10",
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(2,3,5,10),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "b",
                                  sy = "1952.10",
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult=c("2x","2xt","2xtr"),
                                  sigmaP = "b",
                                  sy = "1952.10",
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "0",
                                  sy = "1952.10",
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFSTTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "b",
                                  sy = c("1955.09","1988.05","1988.07","1988.09"),
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFLNFe","oFSTTGFe"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(2,5),
                                  qeff=c("newK"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "b",
                                  sy = "1952.10",
                                  sigma_e = "n"),
                            expand.grid(exec=c("oFSTTGFestF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c("f0.2"),
                                  obs1952 = c("b"),
                                  obs1954 = c("b"),
                                  n_step=c(1),
                                  qeff=c("newK","asF"),
                                  shape=c("b"),
                                  catch_mult="b",
                                  sigmaP = "b",
                                  sy = "1952.10",
                                  sigma_e = "n")                        
    )

    model_config_df = unique(model_config_df)
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executables
    exec_name_vec = c("bspm_estqsimple_softdep_fullmvprior_x0_sttgamma_flexsigmaC_OPT",
                      "bspm_estqsimple_softdep_fullmvprior_x0_LN_flexsigmaC_OPT",
                      "bspm_estq_softdep_fullmvprior_x0_LN_flexsigmaC_OPT",
                      "bspm_estq_softdep_fullmvprior_x0_sttgamma_flexsigmaC_OPT",
                      "bspm_estF_softdep_mvprior_x0_sttgamma_flexsigmaC_OPT")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

    for(i in 18:nrow(model_config_df)){
            run_label_stem = paste0(model_config_df$cpue[i],
                            "-exe",model_config_df$exec[i],
                            "-c",model_config_df$sigma_catch[i],
                            "-n",model_config_df$shape[i],
                            "-q",model_config_df$qeff[i],
                            "-s",model_config_df$n_step[i],
                            "-o52",model_config_df$obs1952[i],
                            "-o54",model_config_df$obs1954[i],
                            "-cm",model_config_df$catch_mult[i],
                            "-sP",model_config_df$sigmaP[i],
                            "-sy",model_config_df$sy[i],
                            "-sE",model_config_df$sigma_e[i],
                            "_0")
            run_number = 109 + i
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
        } else if (model_config_df$exec[i] == "oFLNF"){
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025
        } else if (model_config_df$exec[i] == "oFLNFe"){
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025

            if(model_config_df$sigma_e[i] != "n"){
                stan.data$sigma_edev = as.numeric(as.character(model_config_df$sigma_e[i]))
            } else {
                stan.data$sigma_edev = 0.3
            }
        } else if (model_config_df$exec[i] == "oFSTTGFe"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025

            if(model_config_df$sigma_e[i] != "n"){
                stan.data$sigma_edev = as.numeric(as.character(model_config_df$sigma_e[i]))
            } else {
                stan.data$sigma_edev = 0.3
            }
        } else if (model_config_df$exec[i] == "oFSTTGFestF"){
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)[1:4]
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)[1:4]
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)[1:4,1:4]

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025

            if(model_config_df$sigma_e[i] != "n"){
                stan.data$sigma_edev = as.numeric(as.character(model_config_df$sigma_e[i]))
            } else {
                stan.data$sigma_edev = 0.3
            }

            # sigmaF prior
            if(model_config_df$qeff[i] == "newK"){
                stan.data$PriorSD_sigmaf = ref_0100_F_variability$mean_F_mean*sqrt(pi/2)
            } else if(model_config_df$qeff[i] == "asF"){
                stan.data$PriorSD_sigmaf = 0.25*ref_0100_F_variability$mean_F_mean*sqrt(pi/2)
            }
            
        }

                stan_c = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = stan_c.list[[1]],
                       "oFLNF" = stan_c.list[[2]],
                       "oFLNFe" = stan_c.list[[3]],
                       "oFSTTGFe" = stan_c.list[[4]],
                       "oFSTTGFestF" = stan_c.list[[5]])    # bspm_estqsimple_softdep_mvprior_x0_stt
                
                exec_name = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = exec_name_vec[1],
                       "oFLNF" = exec_name_vec[2],
                       "oFLNFe" = exec_name_vec[3],
                       "oFSTTGFe" = exec_name_vec[4],
                       "oFSTTGFestF" = exec_name_vec[5])

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

        # catch mult
            if(model_config_df$catch_mult[i] == "2x"){
                stan.data$obs_removals = stan.data$obs_removals * 2
            } else if(model_config_df$catch_mult[i] == "2xt"){
                stan.data$obs_removals = stan.data$obs_removals * seq(from=1,to=2,length.out=length(stan.data$obs_removals))
            } else if(model_config_df$catch_mult[i] == "2xtr"){
                stan.data$obs_removals = stan.data$obs_removals * seq(from=2,to=1,length.out=length(stan.data$obs_removals))
            }

        # sigmaP
            if(model_config_df$sigmaP[i] == "0"){
                stan.data$PriorMean_logsigmap = log(0.001)
                stan.data$PriorSD_logsigmap = 0.05
            }

        # Start year adjustments
        start_year = as.numeric(strsplit(as.character(model_config_df$sy[i]),"[.]")[[1]][1])
        init_dep_prior = as.numeric(strsplit(as.character(model_config_df$sy[i]),"[.]")[[1]][2])/10
        stan.data$mv_prior_mean[4] = log(init_dep_prior)
        if(start_year != 1952){
            year_diff = abs(diff(c(start_year, 1952)))
            stan.data$T = stan.data$T - year_diff
            stan.data$index = stan.data$index[-seq(1, year_diff), , drop=FALSE]
            stan.data$sigmao_mat = stan.data$sigmao_mat[-seq(1, year_diff), , drop=FALSE]
            stan.data$effort = stan.data$effort[-seq(1, year_diff)]
            stan.data$sigmac = stan.data$sigmac[-seq(1, year_diff)]
            stan.data$obs_removals = stan.data$obs_removals[-seq(1, year_diff)]
            stan.data$t_dep = stan.data$t_dep - year_diff
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

        
        quick_diagnostics(fit)


        t=as.data.table(summary(fit)$summary)
        t$name = rownames(summary(fit)$summary)
        na.omit(t)[n_eff<500|Rhat>1.01]

    }
    
    
   