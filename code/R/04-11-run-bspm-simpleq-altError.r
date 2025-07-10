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
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0_refined","updated_stan_data.RData"))
    refined_stan_data = updated_stan_data
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0","updated_stan_data.RData"))
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name_vec = c("bspm_estqsimple_softdep_mvprior_x0_stt","bspm_estqsimple_softdep_mvprior_x0_sttgamma_flexsigmaC")
    stan_c.list = as.list(rep(NA,length(exec_name_vec)))
    for(i in 1:length(exec_name_vec)){
        stan_c.list[[i]] = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name_vec[i],".stan")), model_name = exec_name_vec[i])
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop model grid
    model_config_df = rbind( expand.grid(exec=c("STT","STTGF"),
                                  cpue=c("dwfn"),
                                  sigma_catch = c(0.2),
                                  n_step=c(1),
                                  qeff="b",
                                  shape="b")
    )

    model_config_df = unique(model_config_df)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

    for(i in 1:nrow(model_config_df)){
            run_label_stem = paste0(model_config_df$cpue[i],"-exe",model_config_df$exec[i],"-c",model_config_df$sigma_catch[i],"-n",model_config_df$shape[i],"-q",model_config_df$qeff[i],"-s",model_config_df$n_step[i],"_0")
            run_number = 66 + i
            run_number = sprintf("%04d", run_number)
        if(model_config_df$qeff[i] == "b"){
            stan.data = updated_stan_data
        } else if(model_config_df$qeff[i] == "r"){
            stan.data = refined_stan_data
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
            stan.data$sigmac = rep(stan.data$sigmac,stan.data$T)
        }

                stan_c = switch(model_config_df$exec[i],
                       "STT" = stan_c.list[[1]],
                       "STTGF" = stan_c.list[[1]])    # bspm_estqsimple_softdep_mvprior_x0_stt

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
        na.omit(t)[n_eff<500|Rhat>1.01]
    }
    
    
   