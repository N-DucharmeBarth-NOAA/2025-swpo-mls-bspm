# Nicholas Ducharme-Barth
# 2025/06/18
# Prior pushforward check for BSPM parameters with effort-based fishing mortality
# Fletcher-Schaefer production model with time-varying catchability and effort errors

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(GGally)
    library(MASS)
    library(mvtnorm)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    plot_dir = file.path(proj_dir, "plots", "pushforward","q")
    pushforward_output_dir = file.path(proj_dir, "data", "output","pushforward","q")
    dir.create(plot_dir, recursive = TRUE)
    dir.create(pushforward_output_dir, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read data
    bio_params_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))

    # replace logK values
    # original: log(rlnorm(n_generate,log(1e6),0.5))
    set.seed(777)
    old_logK = bio_params_dt$logK
    bio_params_dt$logK = log(rlnorm(nrow(bio_params_dt),log(2e6),1.25))
    
    # plot(density(bio_params_dt$logK))
    # lines(density(old_logK),col="red")

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
# define effort-based production function with dual error structure
    set.seed(123)

    generate_declining_qdev <- function(n_periods) {
        # Randomly choose decline type with broader variety
        decline_type <- sample(c("exponential", "linear", "power", "logistic"), 1, 
                                prob = c(0.4, 0.2, 0.2, 0.2))
        
        # Broader random parameters
        initial_level <- rlnorm(1, 0, 1)  
        noise_sd <- runif(1, 0.05, 0.5)       # Wide noise range
        periods <- 1:n_periods
        
        if(decline_type == "exponential") {
            decay_rate <- runif(1, 0.05, 2)   # Slow to very fast decay
            qdev_trend <- initial_level * exp(-decay_rate * (periods - 1))
            # Add extra boost to first 2-3 periods
            early_boost <- c(runif(min(3, n_periods), 0, 0.5), rep(0, max(0, n_periods-3)))
            qdev_trend <- qdev_trend + early_boost
        } else if(decline_type == "linear") {
            slope <- runif(1, -0.8, -0.02)      # Gentle to steep linear decline
            qdev_trend <- initial_level + slope * (periods - 1)
        } else if(decline_type == "power") {
            power <- runif(1, 0.3, 5.0)         # Wide range of power decay
            qdev_trend <- initial_level * (1 / periods)^(1/power)
        } else { # logistic decline
            midpoint <- runif(1, 2, n_periods-1)
            steepness <- runif(1, 0.5, 3.0)
            asymptote <- runif(1, -0.8, 0.2)
            qdev_trend <- asymptote + (initial_level - asymptote) / 
                        (1 + exp(steepness * (periods - midpoint)))
        }
        
        qdev_noise <- rnorm(n_periods, 0, noise_sd)
        return(qdev_trend + qdev_noise)
    }

    stochastic_surplus_production_effort = function(seed, n, logK, r, qeff, sigma_edev, n_step, sigmap, init_dep, T, effort)
    {
        set.seed(seed)
        # fletcher-schaefer parameters
        dmsy = (1/n)^(1/(n-1))
        h = 2*dmsy
        m = r*h/4
        g = (n^((n/(n-1))))/(n-1)
        
        # effort parameters (qeff already on natural scale)
        n_periods = ceiling((T-1) / n_step)  # Use T-1 for consistency with Stan
        
        # simulate AR(1) qdev (systematic catchability changes)
        qdev_period = generate_declining_qdev(n_periods)
        
        # calculate what rho and sigma_qdev would be based on the simulated declining qdev trend
            rho = cor(qdev_period[-1], qdev_period[-n_periods])
            rho = max(-0.99, min(0.99, rho))  # Keep in bounds

            innovations = qdev_period[-1] - rho * qdev_period[-n_periods]
            sigma_qdev = sd(innovations) / sqrt(1 - rho^2)
        
        # simulate annual effort deviations
        edev = rnorm(T-1, 0, sigma_edev)
        edev[3] = rnorm(1, abs(rnorm(1, 1, 0.3)), 0.2)   # Random large positive
        sigma_edev2 = sigma_edev^2
        
        # calculate F time series: F_t = qeff * exp(qdev_t) * effort_t * exp(edev_t - sigma_edev^2/2)
        F_mort = numeric(T-1)
        qdev = numeric(T-1)
        for(t in 1:(T-1)) {
            period = min(ceiling(t / n_step), n_periods)
            qdev[t] = qdev_period[period]
            F_mort[t] = qeff * exp(qdev[t]) * effort[t] * exp(edev[t] - sigma_edev2/2)
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
            biomass_prod = max(biomass_prod, 0.001)
            
            # Apply process error first, then fishing mortality
            # Process error is applied to biomass after production
            mu_process = log(biomass_prod) - (sigmap^2)/2
            biomass_after_process = max(exp(rnorm(1, mu_process, sigmap)), 0.001)
            
            # Apply fishing mortality
            survival = exp(-F_mort[t-1])
            removals[t-1] = biomass_after_process * (1 - survival) * exp(logK)
            
            # Population after fishing
            x[t] = max(biomass_after_process * survival, 0.001)
        }
        
        return(list(x = x, F = F_mort, qdev = qdev, edev = edev, removals = removals, rho=rep(rho,T), sigma_qdev=rep(sigma_qdev,T), sigma_edev = rep(sd(edev),T)))
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop initial catchability prior based on nominal scaled cpue and preliminary estimates of population scale
    scaled_cpue = (catch_effort_dt$catch_n/(catch_effort_dt$effort_scaled))
    
    min_pop = 250000 # from preliminary runs; with buffer
    max_pop = 3000000 # from preliminary runs; with buffer
    nsim_q = 1e4

    scaled_q = sample(scaled_cpue,nsim_q,replace=TRUE)/runif(nsim_q,min_pop,max_pop) 

    scaled_q_fn = function(par){-sum(dnorm(log(scaled_q), mean = par[1], sd = par[2], log = TRUE))}
    scaled_q_pars = nlminb(c(mean(log(scaled_q)), sd(log(scaled_q))), scaled_q_fn)$par

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prior pushforward with updated parameter structure
    nsim = min(50000, nrow(bio_params_dt))
    sample_indices = sample(1:nrow(bio_params_dt), nsim, replace = FALSE)
    
    prior_dt = data.table(seed = 1:nsim,
                         n = bio_params_dt$shape[sample_indices],
                         logK = bio_params_dt$logK[sample_indices],
                         r = bio_params_dt$rmax[sample_indices],
                         qeff = rlnorm(nsim,scaled_q_pars[1],scaled_q_pars[2]),
                         sigma_edev = rep(0.3, nsim),
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
            sigma_edev = prior_dt$sigma_edev[i],
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
                                      sigma_qdev = result$sigma_qdev,
                                      sigma_edev = result$sigma_edev,
                                      rho = result$rho,
                                      time = catch_effort_dt$time,
                                      dep = result$x,
                                      n = result$x * exp(prior_dt$logK[i]),
                                      F = c(result$F, NA),
                                      qdev = c(result$qdev, NA),
                                      edev = c(result$edev, NA),
                                      u = c(result$removals, NA)/(result$x * exp(prior_dt$logK[i])),
                                      effort = catch_effort_dt$effort_scaled,
                                      predicted_catch = c(result$removals, NA),
                                      observed_catch = catch_effort_dt$catch_n,
                                      predicted_cpue = predicted_cpue,
                                      observed_cpue = observed_cpue,
                                      id = rep(prior_dt$id[i], nrow(catch_effort_dt)))
    }

    sim_dt = rbindlist(sim_dt.list)

    # add metrics per scenario
    sim_dt[,maxF:=max(F,na.rm=TRUE),by=seed]
    sim_dt[,avgF:=mean(F,na.rm=TRUE),by=seed]
    sim_dt[,mindep:=min(dep,na.rm=TRUE),by=seed]
    sim_dt[,minN:=min(n,na.rm=TRUE),by=seed]
    sim_dt[,maxCatch:=max(predicted_catch,na.rm=TRUE),by=seed]
    sim_dt[,minCatch:=min(predicted_catch,na.rm=TRUE),by=seed]
    sim_dt[,avgCatch:=mean(predicted_catch,na.rm=TRUE),by=seed]

    # Period-specific metrics to capture your pattern
    sim_dt[, avg_catch_early := mean(predicted_catch[time <= 1960], na.rm=TRUE), by = seed]  # Early high period

    # define filters
    seed_surv = unique(sim_dt[mindep>0.01&minN>15000&shape<20]$seed)
    seed_f = unique(sim_dt[seed %in% seed_surv & maxF<2.5 & avgF<0.8]$seed)

    # Percentile-based catch filter 
    obs_catch_early_avg = mean(catch_effort_dt[time<1962]$catch_n)    # early period average
    obs_catch_overall_avg = mean(catch_effort_dt$catch_n)  # overall average  
    obs_max_catch = max(catch_effort_dt$catch_n)         # max
    obs_min_catch = min(catch_effort_dt$catch_n)          # min

    seed_catch = unique(sim_dt[
        seed %in% seed_f & 
        
        # Overall bounds
        maxCatch > 0.5 * obs_max_catch &
        maxCatch < 3.0 * obs_max_catch &
        minCatch > 0.5 * obs_min_catch &
        avgCatch > 0.5 * obs_catch_overall_avg &
        avgCatch < 2.0 * obs_catch_overall_avg &
        
        # Early period pattern  
        avg_catch_early > 0.5 * obs_catch_early_avg &
        avg_catch_early < 3 * obs_catch_early_avg
        
    ]$seed)

    rmax_exprior_id = unique(sim_dt[seed %in% seed_catch]$id)
    saveRDS(rmax_exprior_id, file.path(pushforward_output_dir,"rmax_exprior_id.rds"))
    fwrite(sim_dt,file.path(pushforward_output_dir,"sim_dt.csv"))
#________________________________________________________________________________________________________________________________________________________________________________________________________
# plot

    set.seed(123)
    seed_plot = sample(unique(sim_dt$seed),1000)
    plot_seed_catch =sample(unique(seed_catch),500)
    p1 = sim_dt %>%
        .[seed%in%seed_plot | seed %in% plot_seed_catch] %>%
        .[,type:="Un-filtered"]
    p2 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_surv | seed %in% plot_seed_catch] %>%
        .[,type:="Survival filter"]
    p3 = sim_dt %>%
        .[seed%in%seed_plot & seed %in% seed_f | seed %in% plot_seed_catch] %>%
        .[,type:="Fishing mortality filter"]
    p4 = sim_dt %>%
        .[seed %in% plot_seed_catch] %>%
        .[,type:="Catch filter"]

    p = rbind(p1,p2,p3,p4) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Survival filter","Fishing mortality filter","Catch filter"))] %>%
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

    p = rbind(p4) %>%
        .[,type:=factor(type,levels=c("Un-filtered","Survival filter","Fishing mortality filter","Catch filter"))] %>%
        # Calculate median and 95% quantile by time, type, and sigmaF
        .[, .(
            median_removals = median(predicted_catch, na.rm = TRUE),
            q95_removals = quantile(predicted_catch, 0.95, na.rm = TRUE),
            q05_removals = quantile(predicted_catch, 0.05, na.rm = TRUE)
        ), by = .(time, type)] %>%
        ggplot() +
        facet_wrap(~type) +
        ylim(0,NA) +
        xlab("Year") +
        ylab("Removals") +
        geom_ribbon(aes(x=time, ymin=q05_removals, ymax=q95_removals), alpha=0.3,fill="blue") +
        geom_line(aes(x=time, y=median_removals), linewidth=1, color="blue") +
        geom_line(data=catch_dt,aes(x=time,y=catch_n),linewidth=1) +
        geom_hline(yintercept=0) +
        theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                            panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                            panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                            strip.background =element_rect(fill="white"),
                            legend.key = element_rect(fill = "white"))

    ggsave(filename=paste0("prior_push.removals.png"), plot =p, device = "png", path = plot_dir,
  			scale = 1, width = 12, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)

    # Combine all three groups
        combined_df = rbind(p1,p2,p3,p4) %>%
                      .[,.(seed,type,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.) %>%
                       as.data.frame(.)
        
        p_pairs = ggpairs(combined_df,columns=c(3:8), aes(color=type, alpha=0.7),
                        upper=list(continuous="cor"),
                        lower=list(continuous="points"),
                        diag=list(continuous="densityDiag")) +
                viridis::scale_color_viridis("Filter",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
		        viridis::scale_fill_viridis("Filter",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                theme(
                        text = element_text(size = 20),
                        panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                        strip.background = element_rect(fill = "white"),
                        legend.key = element_rect(fill = "white")
                    )
        ggsave(filename="multivariate_pairs_comparison.png", plot=p_pairs, path=plot_dir,
            width=12, height=10, units="in", dpi=300)

    # Just for last filter
        combined_df = sim_dt %>%
                    .[seed %in% seed_catch] %>%
                    .[,type:="Catch filter"] %>%
                      .[,.(seed,type,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.) %>%
                       as.data.frame(.)
        
        p_pairs = ggpairs(combined_df,columns=c(3:8), aes(color=type, alpha=0.7),
                        upper=list(continuous="cor"),
                        lower=list(continuous="points"),
                        diag=list(continuous="densityDiag")) +
                viridis::scale_color_viridis("Filter",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
		        viridis::scale_fill_viridis("Filter",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                theme(
                        text = element_text(size = 20),
                        panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                        strip.background = element_rect(fill = "white"),
                        legend.key = element_rect(fill = "white")
                    )
        ggsave(filename="multivariate_catch_filter.png", plot=p_pairs, path=plot_dir,
            width=12, height=10, units="in", dpi=300)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define priors based on the pushforward
    # Create filtered datasets based on the three filter states
    sim_unfiltered_dt = sim_dt %>%  # All simulations
                     .[,.(seed,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.)

    sim_surv_dt = sim_dt[seed %in% seed_surv] %>%  # Survival filter
                    .[,.(seed,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.)
    sim_f_dt = sim_dt[seed %in% seed_f] %>%  # Fishing mortality filter  
                    .[,.(seed,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.)
    sim_catch_dt = sim_dt[seed %in% seed_catch] %>% # Fishing mortality filter  
                    .[,.(seed,rmax,shape,logK,qeff,rho,sigma_qdev,sigma_edev)] %>%
                       unique(.)

    # Create multivariate parameter matrices (log-transformed where appropriate)
    mv_params_catch = cbind(sim_catch_dt$logK, log(sim_catch_dt$rmax), log(sim_catch_dt$shape))
    colnames(mv_params_catch) = c("logK", "log_r", "log_shape")
    
    # MLE fitting for multivariate normal distributions with better error handling
    mv_mle_catch = function(par) {
        mu = par[1:3]
        # Use Cholesky decomposition parameterization for positive definiteness
        L_vec = par[4:9]  # 6 parameters for lower triangular Cholesky factor
        L = matrix(0, 3, 3)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)  # This ensures positive definiteness
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(mv_params_catch, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }

    # Better starting values using Cholesky decomposition
    init_cov_catch = cov(mv_params_catch)
    init_chol_catch = chol(init_cov_catch)  # Upper triangular
    init_chol_lower_catch = t(init_chol_catch)       # Convert to lower triangular

    init_catch = c(colMeans(mv_params_catch), 
                    init_chol_lower_catch[lower.tri(init_chol_lower_catch, diag=TRUE)])


    # Fit with more robust settings
    mv_fit_catch = nlminb(init_catch, mv_mle_catch,
                            control = list(eval.max = 1000, iter.max = 500))

    # Check convergence and extract parameters for FILTERED
    if(mv_fit_catch$convergence != 0 && mv_fit_catch$convergence != 4) {
        warning("Filtered MLE fitting may not have converged properly. Using sample moments instead.")
        mv_mean_catch = colMeans(mv_params_catch)
        mv_cov_catch = cov(mv_params_catch)
    } else {
        # Extract fitted parameters from Cholesky parameterization
        mv_mean_catch = mv_fit_catch$par[1:3]
        L_fitted_catch = matrix(0, 3, 3)
        L_fitted_catch[lower.tri(L_fitted_catch, diag=TRUE)] = mv_fit_catch$par[4:9]
        mv_cov_catch = L_fitted_catch %*% t(L_fitted_catch)
    }

    # Add names to means and covariance matrices
    names(mv_mean_catch) = c("logK", "log_r", "log_shape")
    dimnames(mv_cov_catch) = list(c("logK", "log_r", "log_shape"), c("logK", "log_r", "log_shape"))
    
    # Calculate correlations
    mv_cor_catch = cov2cor(mv_cov_catch)
    
    # Save multivariate parameters
    write.csv(mv_mean_catch, file.path(pushforward_output_dir, "mv_mean_catch.csv"))
    write.csv(mv_cov_catch, file.path(pushforward_output_dir, "mv_cov_catch.csv"))
    write.csv(mv_cor_catch, file.path(pushforward_output_dir, "mv_cor_catch.csv"))

    # Create ggpairs plots for effort-based pushforward analysis
        set.seed(123)
        n_samples = min(1000, nrow(sim_catch_dt))  # Sample size for comparison

        # Original unfiltered samples
        original_indices = sample(1:nrow(sim_unfiltered_dt), n_samples)
        original_df = data.frame(
            logK = sim_unfiltered_dt$logK[original_indices],
            log_r = log(sim_unfiltered_dt$rmax[original_indices]),
            log_shape = log(sim_unfiltered_dt$shape[original_indices]),
            type = "Original"
        )

        # Actual filtered samples (from your catch filter)
        filtered_indices = sample(1:nrow(sim_catch_dt), n_samples)
        actual_filtered_df = data.frame(
            logK = sim_catch_dt$logK[filtered_indices],
            log_r = log(sim_catch_dt$rmax[filtered_indices]),
            log_shape = log(sim_catch_dt$shape[filtered_indices]),
            type = "Actual Filtered"
        )

        # Samples from fitted filtered distribution
        fitted_filtered_samples = mvrnorm(n_samples, mv_mean_catch, mv_cov_catch)
        fitted_filtered_df = data.frame(
            logK = fitted_filtered_samples[,1],
            log_r = fitted_filtered_samples[,2], 
            log_shape = fitted_filtered_samples[,3],
            type = "Fitted Filtered"
        )

        # Combine all three groups
        combined_df = rbind(original_df, actual_filtered_df, fitted_filtered_df)
        combined_df$type = factor(combined_df$type, levels = c("Original", "Actual Filtered", "Fitted Filtered"))

        p_pairs = ggpairs(combined_df, columns=1:3, aes(color=type, alpha=0.7),
                        upper=list(continuous="cor"),
                        lower=list(continuous="points"),
                        diag=list(continuous="densityDiag")) +
                scale_color_manual(values = c("Original" = "gray", 
                                            "Actual Filtered" = "blue", 
                                            "Fitted Filtered" = "darkred")) +
                scale_fill_manual(values = c("Original" = "gray", 
                                            "Actual Filtered" = "blue", 
                                            "Fitted Filtered" = "darkred")) +
                theme(
                        text = element_text(size = 20),
                        panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                        strip.background = element_rect(fill = "white"),
                        legend.key = element_rect(fill = "white")
                    )

        ggsave(filename="multivariate_pairs_comparison.logk_rmax_shape.samples_fit.png", plot=p_pairs, path=plot_dir,
            width=12, height=10, units="in", dpi=300)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# Calculate joint prior for rho and sigma_qdev

    # Create bivariate parameter matrix on transformed scale
    qdev_params_transformed = cbind(atanh(sim_catch_dt$rho), log(sim_catch_dt$sigma_qdev))
    colnames(qdev_params_transformed) = c("atanh_rho", "log_sigma_qdev")
    
    # MLE fitting for bivariate normal distribution
    qdev_mle = function(par) {
        mu = par[1:2]
        L_vec = par[3:5]  # 3 parameters for 2x2 lower triangular Cholesky factor
        L = matrix(0, 2, 2)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)
        
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(qdev_params_transformed, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }

    # Starting values using Cholesky decomposition
    init_cov_qdev = cov(qdev_params_transformed)
    init_chol_qdev = chol(init_cov_qdev)
    init_chol_lower_qdev = t(init_chol_qdev)
    init_qdev = c(colMeans(qdev_params_transformed), 
                  init_chol_lower_qdev[lower.tri(init_chol_lower_qdev, diag=TRUE)])

    # Fit with robust settings
    qdev_fit = nlminb(init_qdev, qdev_mle, control = list(eval.max = 1000, iter.max = 500))

    # Extract parameters
    if(qdev_fit$convergence != 0 && qdev_fit$convergence != 4) {
        warning("qdev MLE fitting may not have converged properly. Using sample moments instead.")
        qdev_mean = colMeans(qdev_params_transformed)
        qdev_cov = cov(qdev_params_transformed)
    } else {
        qdev_mean = qdev_fit$par[1:2]
        L_fitted_qdev = matrix(0, 2, 2)
        L_fitted_qdev[lower.tri(L_fitted_qdev, diag=TRUE)] = qdev_fit$par[3:5]
        qdev_cov = L_fitted_qdev %*% t(L_fitted_qdev)
    }

    names(qdev_mean) = c("atanh_rho", "log_sigma_qdev")
    dimnames(qdev_cov) = list(c("atanh_rho", "log_sigma_qdev"), c("atanh_rho", "log_sigma_qdev"))
    qdev_cor = cov2cor(qdev_cov)
    
    # Save parameters
    write.csv(qdev_mean, file.path(pushforward_output_dir, "mv_qdev_mean.csv"))
    write.csv(qdev_cov, file.path(pushforward_output_dir, "mv_qdev_cov.csv"))
    write.csv(qdev_cor, file.path(pushforward_output_dir, "mv_qdev_cor.csv"))

    # Create ggpairs plots for rho and sigma_qdev
    set.seed(123)
    n_samples = min(1000, nrow(sim_catch_dt))

    # Original unfiltered samples
    original_indices = sample(1:nrow(sim_unfiltered_dt), n_samples)
    original_qdev_df = data.frame(
        rho = sim_unfiltered_dt$rho[original_indices],
        sigma_qdev = sim_unfiltered_dt$sigma_qdev[original_indices],
        type = "Original"
    )

    filtered_indices = sample(1:nrow(sim_catch_dt), n_samples)
    actual_qdev_df = data.frame(
        rho = sim_catch_dt$rho[filtered_indices],
        sigma_qdev = sim_catch_dt$sigma_qdev[filtered_indices],
        type = "Actual Filtered"
    )

    fitted_qdev_samples = mvrnorm(n_samples, qdev_mean, qdev_cov)
    fitted_qdev_df = data.frame(
        rho = tanh(fitted_qdev_samples[,1]),
        sigma_qdev = exp(fitted_qdev_samples[,2]),
        type = "Fitted"
    )

    combined_qdev_df = rbind(original_qdev_df, actual_qdev_df, fitted_qdev_df)
    combined_qdev_df$type = factor(combined_qdev_df$type, levels = c("Original", "Actual Filtered", "Fitted"))
    
    p_qdev_pairs = ggpairs(combined_qdev_df, columns=1:2, aes(color=type, alpha=0.7),
                    upper=list(continuous="cor"),
                    lower=list(continuous="points"),
                    diag=list(continuous="densityDiag")) +
            scale_color_manual(values = c("Original" = "gray", 
                                        "Actual Filtered" = "blue", 
                                        "Fitted" = "darkred")) +
            scale_fill_manual(values = c("Original" = "gray", 
                                        "Actual Filtered" = "blue", 
                                        "Fitted" = "darkred")) +
            theme(text = element_text(size = 20),
                  panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                  panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                  panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                  strip.background = element_rect(fill = "white"),
                  legend.key = element_rect(fill = "white"))

    ggsave(filename="multivariate_pairs_comparison.rho_sigma_qdev.png", plot=p_qdev_pairs, path=plot_dir,
        width=8, height=6, units="in", dpi=300)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# qeff parameter 
    qeff_fn_catch = function(par){-sum(dnorm(log(sim_catch_dt$qeff), mean = par[1], sd = par[2], log = TRUE))}
    qeff_pars_catch = nlminb(c(mean(log(sim_catch_dt$qeff)), sd(log(sim_catch_dt$qeff))), qeff_fn_catch)$par
    write.csv(qeff_pars_catch, file.path(pushforward_output_dir, "qeff_pars_catch.csv"))

    # Generate samples from fitted distribution
    set.seed(123)
    n_samples = min(1000, nrow(sim_catch_dt))
    fitted_qeff_samples = rlnorm(n_samples, qeff_pars_catch[1], qeff_pars_catch[2])

    # Plot qeff parameter priors - single panel
    png(filename = file.path(plot_dir, "prior.qeff.png"), width = 8, height = 6, units = "in", bg = "white", res = 300)
    
    # Plot density of original samples
    plot(density(sim_unfiltered_dt$qeff), col="gray", lwd=2, 
         main="qeff Prior", xlab="qeff", ylab="Density",
         xlim=c(0,0.4))
    
    # Add density of filtered samples
    lines(density(sim_catch_dt$qeff), col="blue", lwd=2)
    
    # Add density of fitted distribution samples
    lines(density(fitted_qeff_samples), col="darkred", lwd=2)
    
    legend("topright", c("Original", "Filtered", "Fitted"), 
           col=c("gray", "blue", "darkred"), lwd=2, bty="n")
    
    dev.off()

