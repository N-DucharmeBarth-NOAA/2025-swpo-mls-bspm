# Nicholas Ducharme-Barth
# 2025/06/01
# Prior pushforward check for BSPM parameters
# Fletcher-Schaefer production model

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in biological parameter results
    bio_params_dt = fread("./data/output/bspm_parameter_priors_filtered.csv")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read catch data
    catch_dt = fread("./data/input/catch.csv") %>%
               .[Time >= 1952, .(catch_n = Obs * 1000), by = Time] %>%
               setnames(., "Time", "time")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define production function
    set.seed(123)

    stochastic_surplus_production = function(seed, n, logK, r, sigmap, init_dep, T, harvest)
    {
        set.seed(seed)
        # fletcher-schaefer derived parameters following bdm (Edwards 2017; https://github.com/cttedwards/bdm/blob/master/R/bdm.R)
        dmsy = (1/n)^(1/(n-1))
        h = 2*dmsy
        m = r*h/4
        g = (n^((n/(n-1))))/(n-1)

        x = rep(NA,T)
        x[1] = init_dep
        dev_vector = mu_vector = H_vector = rep(NA,T)
        for(t in 2:T){
            H_vector[t] = H = min(c(exp(log(harvest[t-1]) - logK),x[t-1]))
            if(x[t-1]<=dmsy){
                mu_vector[t] = mu = x[t-1] + r * x[t-1] * (1 - x[t-1]/h) - H
            } else {
                mu_vector[t] = mu = x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1)) - H
            }
            if(mu > 0){
                mu = log(mu) - (sigmap^2)/2
            } else {
                mu = log(0.01) - (sigmap^2)/2
            }
            dev_vector[t] = rnorm(1)
            x[t] = exp(dev_vector[t]*sigmap + mu)
        }
        return(x)
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prior pushforward check
    set.seed(123)
    nsim = min(5000, nrow(bio_params_dt))
    sample_indices = sample(1:nrow(bio_params_dt), nsim, replace = FALSE)

    prior_dt = data.table(seed = 1:nsim,
                         n = rep(2, nsim),
                         logK = bio_params_dt$logK[sample_indices],
                         r = bio_params_dt$rmax[sample_indices],
                         sigmap = rep(0, nsim),
                         init_dep = runif(nsim,0.95,1.05))

    sim_dt.list = as.list(rep(NA, nsim))
    for(i in 1:nsim)
    {
        n = stochastic_surplus_production(seed = prior_dt$seed[i], 
                                         n = prior_dt$n[i], 
                                         logK = prior_dt$logK[i], 
                                         r = prior_dt$r[i], 
                                         sigmap = prior_dt$sigmap[i], 
                                         init_dep = prior_dt$init_dep[i], 
                                         T = nrow(catch_dt), 
                                         harvest = catch_dt$catch_n)
        roll_n = c(n[-length(n)], NA)
        sim_dt.list[[i]] = data.table(seed = rep(prior_dt$seed[i], nrow(catch_dt)),
                                     time = floor(catch_dt$time),
                                     dep = n,
                                     n = n*exp(prior_dt$logK[i]),
                                     u = catch_dt$catch_n/(roll_n*exp(prior_dt$logK[i])))
    }

    sim_dt = rbindlist(sim_dt.list)
