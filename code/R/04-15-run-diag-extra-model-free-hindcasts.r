# Nicholas Ducharme-Barth
# 2025/07/18
# Run BSPM with effort-based fishing mortality and updated prior structure
# hindcasts

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
                                  cpue="dwfn",
                                  peel=c(0,5,10,15,20))
                            
    )

    model_config_df = unique(model_config_df)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

    for(i in 1:nrow(model_config_df)){
            run_label_stem = paste0(model_config_df$cpue[i],
                            "-exe",model_config_df$exec[i],
                            "-peel",model_config_df$peel[i],"_0")
            run_number = 134 + i
            run_number = sprintf("%04d", run_number)
       
            stan.data = newlogK_stan_data
            
        stan.data$n_step = 1  # years per period for catchability
        stan.data$n_periods = ceiling((stan.data$T-1) / stan.data$n_step)
        stan.data$fit_to_data = 1L
            stan.data$nu_catch_gamma_shape = 2
            stan.data$nu_catch_gamma_rate = 0.1
            stan.data$sigmac = rep(0.2,stan.data$T)
            stan.data$mv_prior_mean = unname(stan.data$full_mv_prior_mean)
            stan.data$mv_prior_sd = unname(stan.data$full_mv_prior_sd)
            stan.data$mv_prior_corr = unname(stan.data$full_mv_prior_corr)

            # recenter x0 prior
            stan.data$mv_prior_mean[4] = 0
            stan.data$mv_prior_sd[4] = 0.025


                stan_c = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = stan_c.list[[1]])    # bspm_estqsimple_softdep_mvprior_x0_stt
                
                exec_name = switch(as.character(model_config_df$exec[i]),
                       "oFSTTGF" = exec_name_vec[1])

        # Determine which CPUE index to fit and set sigmao_input accordingly   
            # Set lambdas vector based on which index is being fit
            lambdas_vec = rep(0, 6)
            lambdas_vec[1] = 1
            stan.data$lambdas = lambdas_vec
            stan.data$sigmao_input = mean_se[1]

                stan.data$sigmac = rep(0.2,stan.data$T)

        # retro
        if(model_config_df$peel[i]>0){
                run_retro = model_config_df$peel[i]
                tmp_index = stan.data$index
                tmp_sigmao = stan.data$sigmao_mat
                tmp_T = stan.data$T

                tmp_idx = (tmp_T-(run_retro-1)):tmp_T
                tmp_index[tmp_idx,] = -999
                tmp_sigmao[tmp_idx,] = -999

                stan.data$index = tmp_index
                stan.data$sigmao_mat = tmp_sigmao
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
    
    
   