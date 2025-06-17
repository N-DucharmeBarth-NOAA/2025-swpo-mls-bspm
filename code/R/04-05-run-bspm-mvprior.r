# Nicholas Ducharme-Barth
# 2025/06/17
# Run BSPM with informative shape prior

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
    mv_mean_extreme = read.csv(file.path(proj_dir,"data","output","mv_mean_extreme.csv"))
    mv_cov_extreme = read.csv(file.path(proj_dir,"data","output","mv_cov_extreme.csv"))
    mv_cor_extreme = read.csv(file.path(proj_dir,"data","output","mv_cor_extreme.csv"))
    dep_pars = read.csv(file.path(proj_dir,"data","output","rel_dep_ssb_pars.csv"))
    sigmaF_par = read.csv(file.path(proj_dir,"data","output","sigmaF_pars_catch.csv"))
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare data matrices
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    time_years = catch_annual$year
    n_years = length(time_years)
    
    index_mat = matrix(-999, nrow = n_years, ncol = 2)
    se_mat = matrix(-999, nrow = n_years, ncol = 2)

    mean_se = mean(cpue_dt$SE)

    cpue_years = floor(cpue_dt$Time)
    for(i in 1:nrow(cpue_dt)) {
        year_idx = which(time_years == cpue_years[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, ] = cpue_dt$Obs[i]
            se_mat[year_idx, ] = cpue_dt$SE[i]/mean_se
        }
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare mvnorm prior
    # Convert to proper formats
    mv_prior_mean = as.numeric(mv_mean_extreme$x)
    mv_cov_matrix = as.matrix(mv_cov_extreme[,-1])  # Remove row names column
    mv_prior_corr = as.matrix(mv_cor_extreme[,-1])  # Remove row names column

    # Extract standard deviations from covariance matrix
    mv_prior_sd = sqrt(diag(mv_cov_matrix))

    # Optional: Validate that correlation matrix is positive definite
    if(min(eigen(mv_prior_corr)$values) < 1e-10) {
        warning("Correlation matrix has very small eigenvalues, adding small regularization")
        mv_prior_corr = mv_prior_corr + diag(1e-8, 3)
    }
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name = "bspm_estF_softdep_mvprior"
    stan_c = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name,".stan")), model_name =exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs
    i = 1
    run_label_stem = "2024cpueMVPrior_0"
    run_number = "0005"
            
    # specify data and priors
        tmp.data.priors = list(mv_prior_mean = mv_prior_mean,
        mv_prior_sd = mv_prior_sd,
        mv_prior_corr = mv_prior_corr,
        PriorMean_logsigmap = -2.9311037,
        PriorSD_logsigmap =  0.2661089,
        prior_depletion_meanlog = dep_pars$x[1],
        prior_depletion_sdlog = dep_pars$x[2]) 

        tmp_data = list(T=as.integer(n_years),
                        I = as.integer(ncol(index_mat)),
                        index=index_mat,
                        sigmao_mat=se_mat,
                        PriorSD_sigmao_add = 0.2,
                        lambdas=as.vector(c(1,0)),
                        t_dep = 37) # model starts t=1 in 1952, 36 is 1987, 37 is 1988

        tmp_data$sigmao_input = mean_se
        tmp_data$obs_removals = catch_annual$total_catch
        tmp_data$sigmac = 0.2
        tmp_data$PriorSD_sigmaf = sigmaF_par$x[1]                       
        stan.data = c(tmp_data,tmp.data.priors)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inits
            # set.seed(123)
            # stan.inits = replicate(4,stan_inits_func(Tm1 = (stan.data$T-1)),simplify=FALSE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run model
    fit = fit_rstan(stan.data,
                        stan_c,
                        run_label = paste0(run_number,"-",run_label_stem),
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

    print(fit, pars = c("logK", "r", "shape", "sigmap","sigmao_add","sigmaf","x[1]","x[37]","x[71]","removals[3]","removals[70]"))
    print(stan.data$obs_removals[c(3,70)])

#________________________________________________________________________________________________________________________________________________________________________________________________________
# save results
    saveRDS(fit, file.path(proj_dir,"data","output","model_runs",paste0(run_number,"-",run_label_stem,"-fit.rds")))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# run retros
    for(i in 1:5){
        run_label_stem = paste0("2024cpueMVPrior_",i)
        fit = fit_rstan(stan.data,
                            stan_c,
                            run_label = paste0(run_number,"-",run_label_stem),
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

