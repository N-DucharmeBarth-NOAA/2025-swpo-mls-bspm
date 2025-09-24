

# NOAA PIFSC
# 2025/09/23
# Format the data for the Bayesian State-Space Surplus Production model (BSPM)

# Copyright (c) 2025 NOAA PIFSC
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(TAF)
    library(data.table)
    library(magrittr)
#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for processed data
    mkdir("data")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# bring in data
        load(file.path("boot","data","updated_stan_data.RData"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# process catch & effort data
# catch units are in numbers of individuals
# effort units are in 0.0001 hooks
        catch_effort_dt = data.table(year = 1952:2022,
                                     obs_catch = updated_stan_data$obs_removals,
                                     obs_effort = updated_stan_data$effort)

        write.taf(catch_effort_dt, dir="data")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# process index data
# index standard errors are standardized to a mean of 1 for each index

        fleet_names = c("dwfn","au","nz","obs","obsNoPF","obsPFonly")

        index_dt = as.data.table(updated_stan_data$index) %>%
                .[,year:=1952:2022] %>%
                melt(.,id.vars=c("year")) %>%
                setnames(.,c("variable","value"),c("fleet","idx")) %>%
                .[,fleet:=fleet_names[as.numeric(substr(fleet,2,2))]] %>%
                .[idx==-999,idx:=NA]

        se_dt = as.data.table(updated_stan_data$sigmao_mat) %>%
                .[,year:=1952:2022] %>%
                melt(.,id.vars=c("year")) %>%
                setnames(.,c("variable","value"),c("fleet","se")) %>%
                .[,fleet:=fleet_names[as.numeric(substr(fleet,2,2))]] %>%
                .[se==-999,se:=NA]
        index_dt = merge(index_dt,se_dt,by=c("year","fleet"))

        write.taf(index_dt, dir="data")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# write out priors
        # Multivariate prior parameters
        mv_params = c("logK", "log_r", "log_shape", "log_x0", "log_qeff", "atanh_rho", "log_sigma_qdev")
        mv_descriptions = c("Log carrying capacity", "Log maximum intrinsic rate of increase", 
                        "Log shape parameter", "Log initial depletion", "Log catchability",
                        "Arctanh catchability autocorrelation", "Log catchability variability")

        # Create multivariate prior entries (baseline)
        mv_dt = data.table(
        parameter = mv_params,
        distribution = "Multivariate Normal",
        mean = updated_stan_data$full_mv_prior_mean,
        sd = updated_stan_data$full_mv_prior_sd,
        description = mv_descriptions,
        scenario = "baseline"
        )

        # Alternative assumption for shape parameter
        mv_alt_dt = copy(mv_dt)
        mv_alt_dt[parameter == "log_shape", mean := log(2)]
        mv_alt_dt[, scenario := "alternative_shape"]

        # Individual priors (same for both scenarios)
        individual_dt = data.table(
        parameter = c("sigma_P", "sigma_O_add", "nu_catch"),
        distribution = c("Lognormal", "Half-Normal", "Gamma"),
        mean = c(updated_stan_data$PriorMean_logsigmap, 0, 2),
        sd = c(updated_stan_data$PriorSD_logsigmap, updated_stan_data$PriorSD_sigmao_add, 0.1),
        description = c("Process error", "Additional observation error", 
                        "Catch likelihood degrees of freedom"),
        scenario = "both"
        )

        # Combine and return
        prior_dt = rbindlist(list(mv_dt, mv_alt_dt, individual_dt))
        write.taf(prior_dt, dir="data")

        # Extract correlation matrix and create descriptive table
        mv_params = c("logK", "log_r", "log_shape", "log_x0", "log_qeff", "atanh_rho", "log_sigma_qdev")

        # Check actual dimensions first
        corr_data = updated_stan_data$full_mv_prior_corr
        n_params = sqrt(length(corr_data))

        # Get correlation matrix with correct dimensions
        corr_matrix = matrix(corr_data, nrow = n_params, ncol = n_params)
        rownames(corr_matrix) = mv_params[1:n_params]
        colnames(corr_matrix) = mv_params[1:n_params]

        # Convert to long format data.table with upper triangle only
        corr_dt = data.table()
        for(i in 1:(n_params-1)) {
                for(j in (i+1):n_params) {
                        corr_dt = rbind(corr_dt, data.table(
                        param1 = mv_params[i],
                        param2 = mv_params[j], 
                        correlation = corr_matrix[i, j]
                        ))
                }
        }
        write.taf(corr_dt, dir="data")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define scenarios
        ssp_scenarios =  expand.grid(  run_number = "1",
                                        cpue=c("dwfn","nz"),
                                        shape=c("b","alt")) 
                                
        ssp_scenarios$run_number = c("0100","0102","0105","0107")
        write.taf(ssp_scenarios, dir="data")
