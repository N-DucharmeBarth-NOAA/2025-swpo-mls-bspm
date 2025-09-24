

# NOAA PIFSC
# 2025/09/24
# Run the diagnostic case model
# Model 0100: DWFN index; baseline priors

# Copyright (c) 2025 NOAA PIFSC
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(TAF)
    library(data.table)
    library(magrittr)
    library(rstan)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    source_dir_stem = file.path("boot","software")
    source_files = list.files(source_dir_stem)
    source_files = source_files[grep(".r",source_files,fixed=TRUE)]
    sapply(file.path(source_dir_stem,source_files),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for model outputs
    mkdir("model")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load inputs
    load(file.path("boot","data","updated_stan_data.RData"))
    ssp_scenarios = fread(file.path("data","ssp_scenarios.csv"))
    mean_se = fread(file.path("boot","data","mean_se.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name = "bspm_estqsimple_softdep_fullmvprior_x0_sttgamma_flexsigmaC_OPT"
    stan_c = stan_model(file=file.path("boot","software",paste0(exec_name,".stan")), model_name =exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs
    i = 1
    run_label_stem = paste0(c(ssp_scenarios$cpue[i],ssp_scenarios$shape[i]),collapse="-")    
    run_number = paste0("0",ssp_scenarios$run_number[i])

    stan.data = updated_stan_data

        stan.data$n_step =  1
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

        cpue_idx = switch(as.character(ssp_scenarios$cpue[i]),
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
        stan.data$sigmao_input = mean_se$mean_se[cpue_idx]

        # alt shape
            if(ssp_scenarios$shape[i] == "alt"){
                stan.data$mv_prior_mean[3] = log(2)
            }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run model
    cat("Note: Results may differ from original analysis due to parallel execution.\n")
    cat("Setting mc.cores = 1 for reproducibility.\n")
    options(mc.cores = 1)

    fit = fit_rstan(stan.data,
                        stan_c,
                        run_label = paste0(run_number,"-",run_label_stem,"_0"),
                        exec_name = exec_name,
                        seed = 321,
                        chains = 5,
                        n_thin = 10,
                        iter_keep = 200,
                        burnin.prop = 0.5,
                        adapt_delta = 0.99,
                        max_treedepth = 12,
                        silent = FALSE,
                        stan_save_dir = file.path("model"),
                        n_cores = 1)


print(fit, pars = c("logK", "r", "shape","x0", "sigmap", "sigmao_add", "qeff","rho", "sigma_qdev", "nu_catch",
                        "x[1]", "x[37]", "x[71]", "removals[3]", "removals[70]")) 
quick_diagnostics(fit)

