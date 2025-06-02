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
    library(ggplot2)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
    plot_dir = file.path(proj_dir, "plots", "pushforward")
    dir.create(plot_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in biological parameter results
    bio_params_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read catch data
    catch_dt = fread(file.path(proj_dir,"data","input","catch.csv")) %>%
               .[,Time := floor(Time)] %>%
               .[Time >= 1952, .(catch_n = sum(Obs * 1000)), by = Time] %>%
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
    nsim = min(10000, nrow(bio_params_dt))
    sample_indices = sample(1:nrow(bio_params_dt), nsim, replace = FALSE)

    prior_dt = data.table(seed = 1:nsim,
                         n = rep(2, nsim),
                         logK = log(rlnorm(nsim,log(10e5),0.5)),
                         r = bio_params_dt$rmax[sample_indices],
                         sigmap = rep(0, nsim),
                         init_dep = rep(1, nsim))

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
                                      rmax = rep(prior_dt$r[i], nrow(catch_dt)),
                                      logK = rep(prior_dt$logK[i], nrow(catch_dt)),
                                     time = floor(catch_dt$time),
                                     dep = n,
                                     n = n*exp(prior_dt$logK[i]),
                                     u = catch_dt$catch_n/(roll_n*exp(prior_dt$logK[i])))
    }

    sim_dt = rbindlist(sim_dt.list)

    # Add percentage change in n from 2010 to end of time series
    sim_dt[, pct_change_n := {
        n_2010 = n[time == 2010]
        n_end = n[time == max(time)]
        if(length(n_2010) > 0 & length(n_end) > 0) {
            ((n_end - n_2010) / n_2010) * 100
        } else {
            NA_real_
        }
    }, by = seed]

    # define filters
    # survival
        seed_surv = sim_dt[time==2022&dep>0.02&n>15000]$seed
        seed_flat = sim_dt[time==2022&dep>0.02&n>15000&pct_change_n > -5]$seed

#________________________________________________________________________________________________________________________________________________________________________________________________________
# plot

    set.seed(123)
    seed_plot = sample(unique(sim_dt$seed),500)
    p1 = sim_dt %>%
        .[seed%in%seed_plot] %>%
        .[,type:="Un-filtered"]
    p2 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_surv] %>%
        .[,type:="Filtered"]
    p3 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_flat] %>%
        .[,type:="Extreme filter"]


    p = rbind(p1,p2,p3) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Filtered","Extreme filter"))] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
		xlab("Year") +
        ylab("Depletion") +
        geom_path(aes(x=time,y=dep,color=rmax,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
    
    ggsave(filename=paste0("prior_push.sim_dep.png"), plot =p, device = "png", path = plot_dir,
  			scale = 1, width = 12, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    p = rbind(p1,p2,p3) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Filtered","Extreme filter"))] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
		xlab("Year") +
        ylab("Numbers") +
        geom_path(aes(x=time,y=n,color=rmax,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
    
    ggsave(filename=paste0("prior_push.sim_n.png"), plot =p, device = "png", path = plot_dir,
  			scale = 1, width = 12, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)
