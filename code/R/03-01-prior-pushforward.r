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
    model_run_dir = file.path(proj_dir, "data", "output")
    dir.create(model_run_dir, recursive = TRUE)

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
                         n = bio_params_dt$shape[sample_indices],
                         logK = bio_params_dt$logK[sample_indices],
                         r = bio_params_dt$rmax[sample_indices],
                         sigmap = rep(0, nsim),
                         init_dep = rep(1, nsim),
                         id = bio_params_dt$id[sample_indices])

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
                                     u = catch_dt$catch_n/(roll_n*exp(prior_dt$logK[i])),
                                     id = rep(prior_dt$id[i], nrow(catch_dt)))
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
        seed_surv = sim_dt[time==2022&dep>0.02&n>15000&rmax<1]$seed
        seed_flat = sim_dt[time==2022&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10]$seed

        rmax_exprior_id = sim_dt[time==2022&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10]$id
        saveRDS(rmax_exprior_id, file.path(proj_dir,"data","output","rmax_exprior_id.rds"))
        fwrite(sim_dt,file.path(proj_dir,"data","output","sim_dt.csv"))

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

    # Create filtered datasets based on the three filter states
    sim_unfiltered_dt = sim_dt  # All simulations
    prior_unfiltered_dt = prior_dt

    sim_filtered_dt = sim_dt[seed %in% seed_surv]  # Survival filter
    prior_filtered_dt = prior_dt[seed %in% seed_surv]

    sim_extreme_dt = sim_dt[seed %in% seed_flat]  # Survival + stability filter  
    prior_extreme_dt = prior_dt[seed %in% seed_flat]

    # Define new filtered priors
    # logK
    logK_fn_filtered = function(par){-sum(dnorm(prior_filtered_dt$logK, mean = par[1], sd = par[2], log = TRUE))}
    logK_pars_filtered = nlminb(c(16.5, 0.4), logK_fn_filtered)$par
    write.csv(logK_pars_filtered, file = file.path(model_run_dir, "logK_pars_filtered.csv"))

    logK_fn_extreme = function(par){-sum(dnorm(prior_extreme_dt$logK, mean = par[1], sd = par[2], log = TRUE))}
    logK_pars_extreme = nlminb(c(16.5, 0.4), logK_fn_extreme)$par
    write.csv(logK_pars_extreme, file = file.path(model_run_dir, "logK_pars_extreme.csv"))

    # Plot logK priors
    png(filename = file.path(plot_dir, "prior.logK.png"), width = 6, height = 6, units = "in", bg = "white", res = 300)
    par(mfrow=c(3,1))

    # Filtered (survival only)
    hist(prior_filtered_dt$logK, freq=FALSE, breaks=50, xlab="logK", main="Prior: LogK - Filtered (Survival)")
    plot_x = seq(from=min(prior_filtered_dt$logK), to=max(prior_filtered_dt$logK), length.out=1000)
    plot_y = dnorm(plot_x, logK_pars_filtered[1], logK_pars_filtered[2])
    lines(plot_x, plot_y, col="red")
    legend("topright", c("Update prior: Filtered"), col=c("red"), lwd=3, bty="n")

    # Extreme (survival + stability)
    hist(prior_extreme_dt$logK, freq=FALSE, breaks=50, xlab="logK", main="Prior: LogK - Extreme (Survival + Stability)")
    plot_x2 = seq(from=min(prior_extreme_dt$logK), to=max(prior_extreme_dt$logK), length.out=1000)
    plot_y2 = dnorm(plot_x2, logK_pars_extreme[1], logK_pars_extreme[2])
    lines(plot_x2, plot_y2, col="red", lty=3)
    legend("topright", c("Update prior: Extreme"), col=c("red"), lwd=3, lty=3, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="logK", ylab="Density", 
        xlim=range(c(plot_x, plot_x2)), ylim=c(0, max(c(plot_y, plot_y2))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(density(prior_unfiltered_dt$logK), col="blue")
    legend("topright", c("Original prior", "Update prior: Filtered", "Update prior: Extreme"), 
        col=c("blue", "red", "red"), lty=c(1, 1, 3), lwd=3, bty="n")
    dev.off()

    # r (rmax)
    rmax_fn_filtered = function(par){-sum(dnorm(log(prior_filtered_dt$r), mean = par[1], sd = par[2], log = TRUE))}
    rmax_pars_filtered = nlminb(c(log(0.1), 0.4), rmax_fn_filtered)$par
    write.csv(rmax_pars_filtered, file = file.path(model_run_dir, "rmax_pars_filtered.csv"))

    rmax_fn_extreme = function(par){-sum(dnorm(log(prior_extreme_dt$r), mean = par[1], sd = par[2], log = TRUE))}
    rmax_pars_extreme = nlminb(c(log(0.1), 0.4), rmax_fn_extreme)$par
    write.csv(rmax_pars_extreme, file = file.path(model_run_dir, "rmax_pars_extreme.csv"))

    # Plot rmax priors
    png(filename = file.path(plot_dir, "prior.rmax.png"), width = 6, height = 6, units = "in", bg = "white", res = 300)
    par(mfrow=c(3,1))

    # Filtered (survival only)
    plot_x = seq(from=0, to=max(prior_filtered_dt$r), length.out=1000)
    hist(prior_filtered_dt$r, freq=FALSE, breaks=50, xlab="Rmax", main="Prior: Rmax - Filtered (Survival)", xlim=range(plot_x))
    plot_y = dlnorm(plot_x, rmax_pars_filtered[1], rmax_pars_filtered[2])
    lines(plot_x, plot_y, col="red")
    legend("topright", c("Update prior: Filtered"), col=c("red"), lwd=3, bty="n")

    # Extreme (survival + stability)
    plot_x2 = seq(from=0, to=max(prior_extreme_dt$r), length.out=1000)
    hist(prior_extreme_dt$r, freq=FALSE, breaks=50, xlab="Rmax", main="Prior: Rmax - Extreme (Survival + Stability)", xlim=range(plot_x))
    plot_y2 = dlnorm(plot_x2, rmax_pars_extreme[1], rmax_pars_extreme[2])
    lines(plot_x2, plot_y2, col="red", lty=3)
    legend("topright", c("Update prior: Extreme"), col=c("red"), lwd=3, lty=3, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="Rmax", ylab="Density", 
        xlim=range(c(plot_x, plot_x2)), ylim=c(0, max(c(plot_y, plot_y2))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(density(prior_unfiltered_dt$r), col="blue")
    legend("topright", c("Original prior", "Update prior: Filtered", "Update prior: Extreme"), 
        col=c("blue", "red", "red"), lty=c(1, 1, 3), lwd=3, bty="n")
    dev.off()

    # Shape parameter (n)
    shape_fn_filtered = function(par){-sum(dnorm(log(prior_filtered_dt$n), mean = par[1], sd = par[2], log = TRUE))}
    shape_pars_filtered = nlminb(c(log(2), 0.4), shape_fn_filtered)$par
    write.csv(shape_pars_filtered, file = file.path(model_run_dir, "shape_pars_filtered.csv"))

    shape_fn_extreme = function(par){-sum(dnorm(log(prior_extreme_dt$n), mean = par[1], sd = par[2], log = TRUE))}
    shape_pars_extreme = nlminb(c(log(2), 0.4), shape_fn_extreme)$par
    write.csv(shape_pars_extreme, file = file.path(model_run_dir, "shape_pars_extreme.csv"))

    # Plot shape parameter priors
    png(filename = file.path(plot_dir, "prior.shape.png"), width = 6, height = 6, units = "in", bg = "white", res = 300)
    par(mfrow=c(3,1))

    # Filtered (survival only)
    hist(prior_filtered_dt$n, freq=FALSE, breaks=50, xlab="Shape", main="Prior: Shape Parameter - Filtered (Survival)")
    plot_x = seq(from=min(prior_filtered_dt$n), to=max(prior_filtered_dt$n), length.out=1000)
    plot_y = dlnorm(plot_x, shape_pars_filtered[1], shape_pars_filtered[2])
    lines(plot_x, plot_y, col="red")
    legend("topright", c("Update prior: Filtered"), col=c("red"), lwd=3, bty="n")

    # Extreme (survival + stability)
    hist(prior_extreme_dt$n, freq=FALSE, breaks=50, xlab="Shape", main="Prior: Shape Parameter - Extreme (Survival + Stability)")
    plot_x2 = seq(from=min(prior_extreme_dt$n), to=max(prior_extreme_dt$n), length.out=1000)
    plot_y2 = dlnorm(plot_x2, shape_pars_extreme[1], shape_pars_extreme[2])
    lines(plot_x2, plot_y2, col="red", lty=3)
    legend("topright", c("Update prior: Extreme"), col=c("red"), lwd=3, lty=3, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="Shape", ylab="Density", 
        xlim=range(c(plot_x, plot_x2)), ylim=c(0, max(c(plot_y, plot_y2))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(density(prior_unfiltered_dt$n), col="blue")
    legend("topright", c("Original prior", "Update prior: Filtered", "Update prior: Extreme"), 
        col=c("blue", "red", "red"), lty=c(1, 1, 3), lwd=3, bty="n")
    dev.off()

    # Print summary statistics
    cat("Summary of filtered priors:\n")
    cat("=========================\n")
    cat("LogK parameters:\n")
    cat("Filtered (mean, sd):", round(logK_pars_filtered, 3), "\n")
    cat("Extreme (mean, sd):", round(logK_pars_extreme, 3), "\n")
    cat("\nRmax parameters:\n")
    cat("Filtered (log mean, log sd):", round(rmax_pars_filtered, 3), "\n")
    cat("Extreme (log mean, log sd):", round(rmax_pars_extreme, 3), "\n")
    cat("\nShape parameters:\n")
    cat("Filtered (mean, sd):", round(shape_pars_filtered, 3), "\n")
    cat("Extreme (mean, sd):", round(shape_pars_extreme, 3), "\n")
    cat("\nNumber of simulations retained:\n")
    cat("Original:", nrow(prior_unfiltered_dt), "\n")
    cat("Filtered:", nrow(prior_filtered_dt), "\n")
    cat("Extreme:", nrow(prior_extreme_dt), "\n")
