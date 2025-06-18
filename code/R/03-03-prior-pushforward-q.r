# Nicholas Ducharme-Barth
# 2025/06/01
# Prior pushforward check for BSPM parameters with effort-based fishing mortality
# Fletcher-Schaefer production model with time-varying catchability

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
    plot_dir = file.path(proj_dir, "plots", "pushforward_effort")
    model_run_dir = file.path(proj_dir, "data", "output")
    dir.create(plot_dir, recursive = TRUE)
    dir.create(model_run_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read data
    bio_params_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))

    catch_dt = fread(file.path(proj_dir,"data","input","catch.csv")) %>%
               .[,Time := floor(Time)] %>%
               .[Time >= 1952, .(catch_n = sum(Obs * 1000)), by = Time] %>%
               setnames(., "Time", "time")

    effort_dt = fread(file.path(proj_dir, "data", "input", "WCPFC_L_PUBLIC_BY_FLAG_YR.CSV")) %>%
                .[, lat_short_d := ifelse(grepl("S$", lat_short), 
                                        -as.numeric(gsub("[NS]", "", lat_short)), 
                                        as.numeric(gsub("[NS]", "", lat_short))) + 2.5] %>%
                .[, lon_short_d := ifelse(grepl("W$", lon_short), 
                                        360 - as.numeric(gsub("[EW]", "", lon_short)), 
                                            as.numeric(gsub("[EW]", "", lon_short))) + 2.5] %>%
                .[lat_short_d > -60 & lat_short_d < 0 & lon_short_d < 230] %>%
                .[!(lat_short_d > -5 & lon_short_d > 210)] %>%
                .[,.(effort_scaled = sum(hhooks)/1e6, mls_n = sum(mls_n)),by=yy] %>%
                .[yy %in% 1952:2022] %>%
                setnames(., "yy", "time")

    # Merge data
    catch_effort_dt = merge(catch_dt, effort_dt, by = "time", all.x = TRUE)
    catch_effort_dt[is.na(effort_scaled), effort_scaled := mean(effort_scaled, na.rm = TRUE)]

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define effort-based production function
    set.seed(123)

    stochastic_surplus_production_effort = function(seed, n, logK, r, qeff, sigma_qdev, rho, n_step, sigmap, init_dep, T, effort)
    {
        set.seed(seed)
        # fletcher-schaefer parameters
        dmsy = (1/n)^(1/(n-1))
        h = 2*dmsy
        m = r*h/4
        g = (n^((n/(n-1))))/(n-1)
        
        # effort parameters (qeff already on natural scale)
        n_periods = ceiling(T / n_step)
        
        # simulate AR(1) qdev
        qdev_period = numeric(n_periods)
        qdev_period[1] = rnorm(1, 0, sigma_qdev)
        if(n_periods > 1) {
            for(p in 2:n_periods) {
                qdev_period[p] = rho * qdev_period[p-1] + rnorm(1, 0, sigma_qdev * sqrt(1 - rho^2))
            }
        }
        
        # calculate F time series
        F_mort = numeric(T-1)
        for(t in 1:(T-1)) {
            period = min(ceiling(t / n_step), n_periods)
            F_mort[t] = qeff * effort[t] * exp(qdev_period[period])
        }

        # population dynamics
        x = numeric(T)
        x[1] = init_dep
        removals = numeric(T-1)
        
        for(t in 2:T){
            # surplus production
            if(x[t-1] <= dmsy){
                biomass_prod = x[t-1] + r * x[t-1] * (1 - x[t-1]/h)
            } else {
                biomass_prod = x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1))
            }
            
            # fishing and process error
            biomass_prod = max(biomass_prod, 0.001)
            survival = exp(-F_mort[t-1])
            removals[t-1] = biomass_prod * (1 - survival) * exp(logK)
            
            mu = log(biomass_prod * survival) - (sigmap^2)/2
            x[t] = max(exp(rnorm(1, mu, sigmap)), 0.001)
        }
        
        return(list(x = x, F = F_mort, removals = removals))
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prior pushforward
    nsim = min(50000, nrow(bio_params_dt))
    sample_indices = sample(1:nrow(bio_params_dt), nsim, replace = FALSE)
    
    prior_dt = data.table(seed = 1:nsim,
                         n = bio_params_dt$shape[sample_indices],
                         logK = bio_params_dt$logK[sample_indices],
                         r = bio_params_dt$rmax[sample_indices],
                         qeff = runif(nsim, 1e-8,1e-1),
                         sigma_qdev = rlnorm(nsim,log(0.3),0.3),
                         rho = runif(nsim,0,1),
                         sigmap = rep(0, nsim),
                         init_dep = rep(1, nsim),
                         id = bio_params_dt$id[sample_indices])

    sim_dt.list = vector("list", nsim)
    
    for(i in 1:nsim)
    {
        result = stochastic_surplus_production_effort(
            seed = prior_dt$seed[i], 
            n = prior_dt$n[i], 
            logK = prior_dt$logK[i], 
            r = prior_dt$r[i], 
            qeff = prior_dt$qeff[i],
            sigma_qdev = prior_dt$sigma_qdev[i],
            rho = prior_dt$rho[i],
            n_step = 2,
            sigmap = prior_dt$sigmap[i], 
            init_dep = prior_dt$init_dep[i], 
            T = nrow(catch_effort_dt), 
            effort = catch_effort_dt$effort_scaled
        )
        
        # Calculate CPUE (marlin per 1000 hooks)
        predicted_cpue = c(result$removals, NA) / (catch_effort_dt$effort_scaled * 1000)
        observed_cpue = catch_effort_dt$catch_n / (catch_effort_dt$effort_scaled * 1000)
        
        sim_dt.list[[i]] = data.table(seed = rep(prior_dt$seed[i], nrow(catch_effort_dt)),
                                      rmax = rep(prior_dt$r[i], nrow(catch_effort_dt)),
                                      logK = rep(prior_dt$logK[i], nrow(catch_effort_dt)),
                                      qeff = rep(prior_dt$qeff[i], nrow(catch_effort_dt)),
                                      shape = rep(prior_dt$n[i], nrow(catch_effort_dt)),
                                      sigma_qdev = rep(prior_dt$sigma_qdev[i], nrow(catch_effort_dt)),
                                      rho = rep(prior_dt$rho[i], nrow(catch_effort_dt)),
                                      time = catch_effort_dt$time,
                                      dep = result$x,
                                      n = result$x * exp(prior_dt$logK[i]),
                                      u =  c(result$removals, NA)/(result$x * exp(prior_dt$logK[i])),
                                      effort = catch_effort_dt$effort_scaled,
                                      predicted_catch = c(result$removals, NA),
                                      observed_catch = catch_effort_dt$catch_n,
                                      predicted_cpue = predicted_cpue,
                                      observed_cpue = observed_cpue,
                                      id = rep(prior_dt$id[i], nrow(catch_effort_dt)))
    }

    sim_dt = rbindlist(sim_dt.list)

    # Add percentage change
    sim_dt[, pct_change_n := {
        n_2010 = n[time == 2010]
        n_end = n[time == max(time)]
        if(length(n_2010) > 0 & length(n_end) > 0) {
            ((n_end - n_2010) / n_2010) * 100
        } else {
            NA_real_
        }
    }, by = seed]

    # add max removals per scenario
    sim_dt[,maxCatch:=max(predicted_catch,na.rm=TRUE),by=seed]
    sim_dt[,avgCatch:=mean(predicted_catch,na.rm=TRUE),by=seed]

    # define filters
    # survival
        seed_surv = sim_dt[time==2022&dep>0.02&n>15000&rmax<1]$seed
        seed_flat = sim_dt[time==2022&dep<0.98&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10]$seed
        seed_catch = sim_dt[time==2022&dep<0.98&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10 & maxCatch < 350000 & avgCatch > 1000]$seed

#________________________________________________________________________________________________________________________________________________________________________________________________________
# plot

    set.seed(123)
    seed_plot = sample(unique(sim_dt$seed),500)
    p1 = sim_dt %>%
        .[seed%in%seed_plot] %>%
        .[,type:="Un-filtered"]
    p2 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_surv] %>%
        .[,type:="Survival filter"]
    p3 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_flat] %>%
        .[,type:="Trend filter"]
    p4 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_catch] %>%
        .[,type:="Catch filter"]

    p = rbind(p1,p2,p3,p4) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Survival filter","Trend filter","Catch filter"))] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
		xlab("Year") +
        ylab("Depletion") +
        geom_path(aes(x=time,y=dep,color=qeff,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("qeff",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("qeff",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))

    p

    pp = rbind(p1,p2,p3,p4) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Survival filter","Trend filter","Catch filter"))] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
		xlab("Year") +
        ylab("Removals") +
        geom_path(aes(x=time,y=predicted_catch,color=qeff,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("qeff",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("qeff",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
    pp

    ppp = rbind(p4) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Survival filter","Trend filter","Catch filter"))] %>%
        # Calculate median and 95% quantile by time, type, and sigmaF
        .[, .(
            median_removals = median(removals, na.rm = TRUE),
            q95_removals = quantile(removals, 0.95, na.rm = TRUE),
            q05_removals = quantile(removals, 0.05, na.rm = TRUE)
        ), by = .(time, type)] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
        xlab("Year") +
        ylab("Removals") +
        geom_ribbon(aes(x=time, ymin=q05_removals, ymax=q95_removals), alpha=0.3,fill="blue") +
        geom_line(aes(x=time, y=median_removals), size=1, color="blue") +
        geom_line(data=catch_dt,aes(x=time,y=catch_n),size=1) +
        geom_hline(yintercept=0) +
        theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                            panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                            panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                            strip.background =element_rect(fill="white"),
                            legend.key = element_rect(fill = "white"))

    ppp
