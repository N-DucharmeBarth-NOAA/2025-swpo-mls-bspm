# Nicholas Ducharme-Barth
# 2025/06/02
# Run BSPM

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
    logK_pars = read.csv(file.path(proj_dir,"data","output","logK_pars_extreme.csv"))
    rmax_pars = read.csv(file.path(proj_dir,"data","output","rmax_pars_extreme.csv"))
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare data matrices
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    time_years = catch_annual$year
    n_years = length(time_years)
    
    index_mat = matrix(-9999, nrow = n_years, ncol = 2)
    se_mat = matrix(-9999, nrow = n_years, ncol = 2)

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
# compile executable
    exec_name = "bspm_estF"
    stan_c = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name,".stan")), model_name =exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs
    i = 1
    run_label_stem = "2024cpueExPrior_0"
    run_number = "0001"
            
    # specify data and priors
        tmp.data.priors = list(PriorMean_logK = logK_pars$x[1],
        PriorSD_logK = logK_pars$x[2],
        PriorMean_logr = rmax_pars$x[1],
        PriorSD_logr = rmax_pars$x[2],
        PriorMean_logsigmap = -2.9311037,
        PriorSD_logsigmap =  0.2661089,
        PriorMean_logshape = log(2.0),
        PriorSD_logshape = 0.1)

        tmp_data = list(T=as.integer(n_years),
                        I = as.integer(ncol(index_mat)),
                        index=index_mat,
                        sigmao_mat=se_mat,
                        PriorSD_sigmao_add = 0.2,
                        lambdas=as.vector(c(1,0)))

        tmp_data$sigmao_input = mean_se
        tmp_data$obs_removals = catch_annual$total_catch
        tmp_data$sigmac = 1
        tmp_data$PriorSD_sigmaf = 0.025                        
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

    print(fit, pars = c("logK", "r", "shape", "sigmap","sigmao_add","sigmaf","x[1]","x[71]"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# save results
    saveRDS(fit, file.path(proj_dir,"data","output","model_runs",paste0(run_number,"-",run_label_stem,"-fit.rds")))



