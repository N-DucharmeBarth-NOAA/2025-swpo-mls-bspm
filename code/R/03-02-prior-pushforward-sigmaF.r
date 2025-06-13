# Nicholas Ducharme-Barth
# 2025/06/13
# Prior pushforward check for BSPM parameters
# Fletcher-Schaefer production model
# add sigmaF

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

    stochastic_surplus_production_estF = function(seed, n, logK, r, sigmap, init_dep, T, sigmaF)
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
        F = abs(rnorm(T)) * sigmaF 
        for(t in 2:T){
            if(x[t-1]<=dmsy){
                mu_vector[t] = mu = (x[t-1] + r * x[t-1] * (1 - x[t-1]/h)) * (exp(-F[t-1]))
                H_vector[t-1] = (x[t-1] + r * x[t-1] * (1 - x[t-1]/h)) * (1-exp(-F[t-1]))
            } else {
                mu_vector[t] = mu = (x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1))) * (exp(-F[t-1]))
                H_vector[t-1] = (x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1))) * (1-exp(-F[t-1]))
            }
            if(mu > 0){
                mu = log(mu) - (sigmap^2)/2
            } else {
                mu = log(0.01) - (sigmap^2)/2
            }
            dev_vector[t] = rnorm(1)
            x[t] = exp(dev_vector[t]*sigmap + mu)
        }
        removals = exp(logK)*H_vector
        return(data.table(x=x,removals=removals))
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
                         sigmaF = abs(rnorm(nsim,0,0.1)),
                         init_dep = rep(1, nsim),
                         id = bio_params_dt$id[sample_indices])

    sim_dt.list = as.list(rep(NA, nsim))
    for(i in 1:nsim)
    {
        n = stochastic_surplus_production_estF(seed = prior_dt$seed[i], 
                                         n = prior_dt$n[i], 
                                         logK = prior_dt$logK[i], 
                                         r = prior_dt$r[i], 
                                         sigmap = prior_dt$sigmap[i], 
                                         init_dep = prior_dt$init_dep[i], 
                                         T = nrow(catch_dt), 
                                         sigmaF = prior_dt$sigmaF[i])
        roll_n = c(n$x[-length(n$x)], NA)
        sim_dt.list[[i]] = data.table(seed = rep(prior_dt$seed[i], nrow(catch_dt)),
                                      rmax = rep(prior_dt$r[i], nrow(catch_dt)),
                                      sigmaF = rep(prior_dt$sigmaF[i], nrow(catch_dt)),
                                      logK = rep(prior_dt$logK[i], nrow(catch_dt)),
                                     time = floor(catch_dt$time),
                                     dep = n$x,
                                     removals = n$removals,
                                     n = n$x*exp(prior_dt$logK[i]),
                                     u = n$removals/(roll_n*exp(prior_dt$logK[i])),
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

    # add max removals per scenario
    sim_dt[,maxCatch:=max(removals,na.rm=TRUE),by=seed]

    # define filters
    # survival
        seed_surv = sim_dt[time==2022&dep>0.02&n>15000&rmax<1]$seed
        seed_flat = sim_dt[time==2022&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10]$seed
        seed_catch = sim_dt[time==2022&dep>0.02&n>15000&rmax<1&pct_change_n > -5 & pct_change_n < 10 & maxCatch < 350000]$seed

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
        geom_path(aes(x=time,y=dep,color=rmax,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("Rmax",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
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
        geom_path(aes(x=time,y=removals,color=sigmaF,group=seed),alpha=0.5) +
        geom_hline(yintercept=0) +
        viridis::scale_color_viridis("SigmaF",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
		viridis::scale_fill_viridis("SigmaF",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
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

    # Create filtered datasets based on the three filter states
    sim_unfiltered_dt = sim_dt  # All simulations
    prior_unfiltered_dt = prior_dt

    sim_filtered_dt = sim_dt[seed %in% seed_surv]  # Survival filter
    prior_filtered_dt = prior_dt[seed %in% seed_surv]

    sim_extreme_dt = sim_dt[seed %in% seed_flat]  # Survival + stability filter  
    prior_extreme_dt = prior_dt[seed %in% seed_flat]

    sim_catch_dt = sim_dt[seed %in% seed_catch]  # Survival + stability + catch filter  
    prior_catch_dt = prior_dt[seed %in% seed_catch]

    # Define new filtered priors
    # logK
    logK_fn_filtered = function(par){-sum(dnorm(prior_filtered_dt$logK, mean = par[1], sd = par[2], log = TRUE))}
    logK_pars_filtered = nlminb(c(16.5, 0.4), logK_fn_filtered)$par

    logK_fn_extreme = function(par){-sum(dnorm(prior_extreme_dt$logK, mean = par[1], sd = par[2], log = TRUE))}
    logK_pars_extreme = nlminb(c(16.5, 0.4), logK_fn_extreme)$par

    logK_fn_catch = function(par){-sum(dnorm(prior_catch_dt$logK, mean = par[1], sd = par[2], log = TRUE))}
    logK_pars_catch = nlminb(c(16.5, 0.4), logK_fn_catch)$par
    write.csv(logK_pars_catch, file = file.path(model_run_dir, "logK_pars_catch.csv"))

    png(filename = file.path(plot_dir, "prior.logK.w_sigmaF.png"), width = 6, height = 8, units = "in", bg = "white", res = 300)
    par(mfrow=c(4,1))

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

    # Catch (survival + stability + catch)
    hist(prior_catch_dt$logK, freq=FALSE, breaks=50, xlab="logK", main="Prior: LogK - Catch (Survival + Stability + Catch)")
    plot_x3 = seq(from=min(prior_catch_dt$logK), to=max(prior_catch_dt$logK), length.out=1000)
    plot_y3 = dnorm(plot_x3, logK_pars_catch[1], logK_pars_catch[2])
    lines(plot_x3, plot_y3, col="red", lty=2)
    legend("topright", c("Update prior: Catch"), col=c("red"), lwd=3, lty=2, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="logK", ylab="Density", 
        xlim=range(c(plot_x, plot_x2, plot_x3)), ylim=c(0, max(c(plot_y, plot_y2, plot_y3))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(plot_x3, plot_y3, col="red", lty=2)
    lines(density(prior_unfiltered_dt$logK), col="blue")
    legend("topright", c("Original prior", "Update prior: Filtered", "Update prior: Extreme", "Update prior: Catch"), 
        col=c("blue", "red", "red", "red"), lty=c(1, 1, 3, 2), lwd=3, bty="n")
    dev.off()

    # r (rmax)
    rmax_fn_filtered = function(par){-sum(dnorm(log(prior_filtered_dt$r), mean = par[1], sd = par[2], log = TRUE))}
    rmax_pars_filtered = nlminb(c(log(0.1), 0.4), rmax_fn_filtered)$par

    rmax_fn_extreme = function(par){-sum(dnorm(log(prior_extreme_dt$r), mean = par[1], sd = par[2], log = TRUE))}
    rmax_pars_extreme = nlminb(c(log(0.1), 0.4), rmax_fn_extreme)$par

    rmax_fn_catch = function(par){-sum(dnorm(log(prior_catch_dt$r), mean = par[1], sd = par[2], log = TRUE))}
    rmax_pars_catch = nlminb(c(log(0.1), 0.4), rmax_fn_catch)$par
    write.csv(rmax_pars_catch, file = file.path(model_run_dir, "rmax_pars_catch.csv"))

    png(filename = file.path(plot_dir, "prior.rmax.w_sigmaF.png"), width = 6, height = 8, units = "in", bg = "white", res = 300)
    par(mfrow=c(4,1))

    # Filtered (survival only)
    hist(log(prior_filtered_dt$r), freq=FALSE, breaks=50, xlab="log(Rmax)", main="Prior: Rmax - Filtered (Survival)")
    plot_x = seq(from=min(log(prior_filtered_dt$r)), to=max(log(prior_filtered_dt$r)), length.out=1000)
    plot_y = dnorm(plot_x, rmax_pars_filtered[1], rmax_pars_filtered[2])
    lines(plot_x, plot_y, col="red")
    legend("topright", c("Update prior: Filtered"), col=c("red"), lwd=3, bty="n")

    # Extreme (survival + stability)
    hist(log(prior_extreme_dt$r), freq=FALSE, breaks=50, xlab="log(Rmax)", main="Prior: Rmax - Extreme (Survival + Stability)")
    plot_x2 = seq(from=min(log(prior_extreme_dt$r)), to=max(log(prior_extreme_dt$r)), length.out=1000)
    plot_y2 = dnorm(plot_x2, rmax_pars_extreme[1], rmax_pars_extreme[2])
    lines(plot_x2, plot_y2, col="red", lty=3)
    legend("topright", c("Update prior: Extreme"), col=c("red"), lwd=3, lty=3, bty="n")

    # Catch (survival + stability + catch)
    hist(log(prior_catch_dt$r), freq=FALSE, breaks=50, xlab="log(Rmax)", main="Prior: Rmax - Catch (Survival + Stability + Catch)")
    plot_x3 = seq(from=min(log(prior_catch_dt$r)), to=max(log(prior_catch_dt$r)), length.out=1000)
    plot_y3 = dnorm(plot_x3, rmax_pars_catch[1], rmax_pars_catch[2])
    lines(plot_x3, plot_y3, col="red", lty=2)
    legend("topright", c("Update prior: Catch"), col=c("red"), lwd=3, lty=2, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="log(Rmax)", ylab="Density", 
        xlim=range(c(plot_x, plot_x2, plot_x3)), ylim=c(0, max(c(plot_y, plot_y2, plot_y3))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(plot_x3, plot_y3, col="red", lty=2)
    lines(density(log(prior_unfiltered_dt$r)), col="blue")
    legend("topright", c("Original prior", "Update prior: Filtered", "Update prior: Extreme", "Update prior: Catch"), 
        col=c("blue", "red", "red", "red"), lty=c(1, 1, 3, 2), lwd=3, bty="n")
    dev.off()



    # Add sigmaF prior calculations (half-normal distribution)
    # Define half-normal density function
    dhalfnorm = function(x, sigma, log = FALSE) {
        if (log) {
            log(2) + dnorm(x, mean = 0, sd = sigma, log = TRUE) - log(sigma) - 0.5 * log(2 * pi)
        } else {
            (2 / (sigma * sqrt(2 * pi))) * exp(-x^2 / (2 * sigma^2))
        }
    }

    sigmaF_fn_filtered = function(par){-sum(dhalfnorm(prior_filtered_dt$sigmaF, sigma = par[1], log = TRUE))}
    sigmaF_pars_filtered = nlminb(c(0.05), sigmaF_fn_filtered)$par
    write.csv(sigmaF_pars_filtered, file = file.path(model_run_dir, "sigmaF_pars_filtered.csv"))

    sigmaF_fn_extreme = function(par){-sum(dhalfnorm(prior_extreme_dt$sigmaF, sigma = par[1], log = TRUE))}
    sigmaF_pars_extreme = nlminb(c(0.05), sigmaF_fn_extreme)$par
    write.csv(sigmaF_pars_extreme, file = file.path(model_run_dir, "sigmaF_pars_extreme.csv"))

    sigmaF_fn_catch = function(par){-sum(dhalfnorm(prior_catch_dt$sigmaF, sigma = par[1], log = TRUE))}
    sigmaF_pars_catch = nlminb(c(0.05), sigmaF_fn_catch)$par
    write.csv(sigmaF_pars_catch, file = file.path(model_run_dir, "sigmaF_pars_catch.csv"))

    # Plot sigmaF priors
    png(filename = file.path(plot_dir, "prior.sigmaF.w_sigmaF.png"), width = 6, height = 8, units = "in", bg = "white", res = 300)
    par(mfrow=c(4,1))

    # Filtered (survival only)
    hist(prior_filtered_dt$sigmaF, freq=FALSE, breaks=50, xlab="sigmaF", main="Prior: SigmaF - Filtered (Survival)")
    plot_x = seq(from=0, to=max(prior_filtered_dt$sigmaF), length.out=1000)
    plot_y = dhalfnorm(plot_x, sigma = sigmaF_pars_filtered[1])
    lines(plot_x, plot_y, col="red")
    legend("topright", c("Update prior: Filtered"), col=c("red"), lwd=3, bty="n")

    # Extreme (survival + stability)
    hist(prior_extreme_dt$sigmaF, freq=FALSE, breaks=50, xlab="sigmaF", main="Prior: SigmaF - Extreme (Survival + Stability)")
    plot_x2 = seq(from=0, to=max(prior_extreme_dt$sigmaF), length.out=1000)
    plot_y2 = dhalfnorm(plot_x2, sigma = sigmaF_pars_extreme[1])
    lines(plot_x2, plot_y2, col="red", lty=3)
    legend("topright", c("Update prior: Extreme"), col=c("red"), lwd=3, lty=3, bty="n")

    # Catch (survival + stability + catch)
    hist(prior_catch_dt$sigmaF, freq=FALSE, breaks=50, xlab="sigmaF", main="Prior: SigmaF - Catch (Survival + Stability + Catch)")
    plot_x3 = seq(from=0, to=max(prior_catch_dt$sigmaF), length.out=1000)
    plot_y3 = dhalfnorm(plot_x3, sigma = sigmaF_pars_catch[1])
    lines(plot_x3, plot_y3, col="red", lty=2)
    legend("topright", c("Update prior: Catch"), col=c("red"), lwd=3, lty=2, bty="n")

    # Comparison plot
    plot(plot_x, plot_y, col="red", type="l", xlab="sigmaF", ylab="Density", 
        xlim=range(c(plot_x, plot_x2, plot_x3)), ylim=c(0, max(c(plot_y, plot_y2, plot_y3))*1.2))
    lines(plot_x2, plot_y2, col="red", lty=3)
    lines(plot_x3, plot_y3, col="red", lty=2)
    lines(density(prior_unfiltered_dt$sigmaF), col="blue")
    # Add initial half-normal prior with sd = 0.025
    plot_x_initial = seq(from=0, to=max(c(plot_x, plot_x2, plot_x3)), length.out=1000)
    plot_y_initial = dhalfnorm(plot_x_initial, sigma = 0.025)
    lines(plot_x_initial, plot_y_initial, col="green", lwd=2)
    legend("topright", c("Initial half-normal (sd=0.025)", "Empirical prior", "Update prior: Filtered", "Update prior: Extreme", "Update prior: Catch"), 
        col=c("green", "blue", "red", "red", "red"), lty=c(1, 1, 1, 3, 2), lwd=c(2, 3, 3, 3, 3), bty="n")   
    dev.off()
