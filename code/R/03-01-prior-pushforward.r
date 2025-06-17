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


#________________________________________________________________________________________________________________________________________________________________________________________________________
# multivariate priors
    library(MASS)
    library(GGally)
    library(corrplot)
    library(mvtnorm)
    
    # Create multivariate parameter matrices (log-transformed where appropriate)
    mv_params_filtered = cbind(prior_filtered_dt$logK, log(prior_filtered_dt$r), log(prior_filtered_dt$n))
    mv_params_extreme = cbind(prior_extreme_dt$logK, log(prior_extreme_dt$r), log(prior_extreme_dt$n))
    colnames(mv_params_filtered) = colnames(mv_params_extreme) = c("logK", "log_r", "log_n")
    
    # MLE fitting for multivariate normal distributions with better error handling
    mv_mle_filtered = function(par) {
        mu = par[1:3]
        # Use Cholesky decomposition parameterization for positive definiteness
        L_vec = par[4:9]  # 6 parameters for lower triangular Cholesky factor
        L = matrix(0, 3, 3)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)  # This ensures positive definiteness
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(mv_params_filtered, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }

    mv_mle_extreme = function(par) {
        mu = par[1:3]
        # Use Cholesky decomposition parameterization for positive definiteness
        L_vec = par[4:9]  # 6 parameters for lower triangular Cholesky factor
        L = matrix(0, 3, 3)
        L[lower.tri(L, diag=TRUE)] = L_vec
        Sigma = L %*% t(L)  # This ensures positive definiteness
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-10) return(1e10)
        
        tryCatch({
            -sum(dmvnorm(mv_params_extreme, mu, Sigma, log=TRUE))
        }, error = function(e) 1e10)
    }

    # Better starting values using Cholesky decomposition
    init_cov_filtered = cov(mv_params_filtered)
    init_chol_filtered = chol(init_cov_filtered)  # Upper triangular
    init_chol_lower_filtered = t(init_chol_filtered)       # Convert to lower triangular

    init_cov_extreme = cov(mv_params_extreme)
    init_chol_extreme = chol(init_cov_extreme)  # Upper triangular
    init_chol_lower_extreme = t(init_chol_extreme)       # Convert to lower triangular

    init_filtered = c(colMeans(mv_params_filtered), 
                    init_chol_lower_filtered[lower.tri(init_chol_lower_filtered, diag=TRUE)])

    init_extreme = c(colMeans(mv_params_extreme),
                    init_chol_lower_extreme[lower.tri(init_chol_lower_extreme, diag=TRUE)])

    # Fit with more robust settings
    mv_fit_filtered = nlminb(init_filtered, mv_mle_filtered,
                            control = list(eval.max = 1000, iter.max = 500))

    mv_fit_extreme = nlminb(init_extreme, mv_mle_extreme,
                        control = list(eval.max = 1000, iter.max = 500))

    # Check convergence and extract parameters for FILTERED
    if(mv_fit_filtered$convergence != 0 && mv_fit_filtered$convergence != 4) {
        warning("Filtered MLE fitting may not have converged properly. Using sample moments instead.")
        mv_mean_filtered = colMeans(mv_params_filtered)
        mv_cov_filtered = cov(mv_params_filtered)
    } else {
        # Extract fitted parameters from Cholesky parameterization
        mv_mean_filtered = mv_fit_filtered$par[1:3]
        L_fitted_filtered = matrix(0, 3, 3)
        L_fitted_filtered[lower.tri(L_fitted_filtered, diag=TRUE)] = mv_fit_filtered$par[4:9]
        mv_cov_filtered = L_fitted_filtered %*% t(L_fitted_filtered)
    }

    # Check convergence and extract parameters for EXTREME
    if(mv_fit_extreme$convergence != 0 && mv_fit_extreme$convergence != 4) {
        warning("Extreme MLE fitting may not have converged properly. Using sample moments instead.")
        mv_mean_extreme = colMeans(mv_params_extreme)
        mv_cov_extreme = cov(mv_params_extreme)
    } else {
        # Extract fitted parameters from Cholesky parameterization
        mv_mean_extreme = mv_fit_extreme$par[1:3]
        L_fitted_extreme = matrix(0, 3, 3)
        L_fitted_extreme[lower.tri(L_fitted_extreme, diag=TRUE)] = mv_fit_extreme$par[4:9]
        mv_cov_extreme = L_fitted_extreme %*% t(L_fitted_extreme)
    }

    # Add names to means and covariance matrices
    names(mv_mean_filtered) = names(mv_mean_extreme) = c("logK", "log_r", "log_n")
    dimnames(mv_cov_filtered) = dimnames(mv_cov_extreme) = list(c("logK", "log_r", "log_n"), c("logK", "log_r", "log_n"))
    
    # Calculate correlations
    mv_cor_filtered = cov2cor(mv_cov_filtered)
    mv_cor_extreme = cov2cor(mv_cov_extreme)
    
    # Save multivariate parameters
    write.csv(mv_mean_filtered, file.path(model_run_dir, "mv_mean_filtered.csv"))
    write.csv(mv_cov_filtered, file.path(model_run_dir, "mv_cov_filtered.csv"))
    write.csv(mv_cor_filtered, file.path(model_run_dir, "mv_cor_filtered.csv"))
    write.csv(mv_mean_extreme, file.path(model_run_dir, "mv_mean_extreme.csv"))
    write.csv(mv_cov_extreme, file.path(model_run_dir, "mv_cov_extreme.csv"))
    write.csv(mv_cor_extreme, file.path(model_run_dir, "mv_cor_extreme.csv"))
    
    
    # Create ggpairs plots
        set.seed(123)
        n_samples = min(1000, nrow(prior_filtered_dt))  # Sample size for comparison
        
        # Original samples (subsample for fair comparison)
        original_indices = sample(1:nrow(prior_unfiltered_dt), n_samples)
        original_df = data.frame(
            logK = prior_unfiltered_dt$logK[original_indices],
            log_r = log(prior_unfiltered_dt$r[original_indices]),
            log_n = log(prior_unfiltered_dt$n[original_indices]),
            type = "Original"
        )
        
        # Samples from fitted filtered distribution
        fitted_filtered_samples = mvrnorm(n_samples, mv_mean_filtered, mv_cov_filtered)
        fitted_filtered_df = data.frame(
            logK = fitted_filtered_samples[,1],
            log_r = fitted_filtered_samples[,2], 
            log_n = fitted_filtered_samples[,3],
            type = "Fitted Filtered"
        )
        
        # Samples from fitted extreme distribution
        fitted_extreme_samples = mvrnorm(n_samples, mv_mean_extreme, mv_cov_extreme)
        fitted_extreme_df = data.frame(
            logK = fitted_extreme_samples[,1],
            log_r = fitted_extreme_samples[,2],
            log_n = fitted_extreme_samples[,3], 
            type = "Fitted Extreme"
        )
        
        # Combine all three groups
        combined_df = rbind(original_df, fitted_filtered_df, fitted_extreme_df)
        combined_df$type = factor(combined_df$type, levels = c("Original", "Fitted Filtered", "Fitted Extreme"))
        
        p_pairs = ggpairs(combined_df, columns=1:3, aes(color=type, alpha=0.7),
                        upper=list(continuous="cor"),
                        lower=list(continuous="points"),
                        diag=list(continuous="densityDiag")) +
                scale_color_manual(values = c("Original" = "blue", 
                                            "Fitted Filtered" = "red", 
                                            "Fitted Extreme" = "darkred")) +
                scale_fill_manual(values = c("Original" = "blue", 
                                            "Fitted Filtered" = "red", 
                                            "Fitted Extreme" = "darkred")) +
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

    # Create additional comparison plot: Original vs Univariate vs Multivariate extreme
    set.seed(123)
    n_samples = min(1000, nrow(prior_extreme_dt))  # Sample size for comparison

    # Original samples (subsample for fair comparison)
    original_indices = sample(1:nrow(prior_extreme_dt), n_samples)
    original_df = data.frame(
        logK = prior_extreme_dt$logK[original_indices],
        log_r = log(prior_extreme_dt$r[original_indices]),
        log_n = log(prior_extreme_dt$n[original_indices]),
        type = "Original"
    )

    # Samples from univariate extreme distributions (independent sampling)
    univariate_extreme_df = data.frame(
        logK = rnorm(n_samples, logK_pars_extreme[1], logK_pars_extreme[2]),
        log_r = rnorm(n_samples, rmax_pars_extreme[1], rmax_pars_extreme[2]),
        log_n = rnorm(n_samples, shape_pars_extreme[1], shape_pars_extreme[2]),
        type = "Univariate Extreme"
    )

    # Samples from fitted multivariate extreme distribution
    fitted_mv_extreme_samples = mvrnorm(n_samples, mv_mean_extreme, mv_cov_extreme)
    fitted_mv_extreme_df = data.frame(
        logK = fitted_mv_extreme_samples[,1],
        log_r = fitted_mv_extreme_samples[,2],
        log_n = fitted_mv_extreme_samples[,3], 
        type = "Multivariate Extreme"
    )

    # Combine all three groups
    combined_extreme_df = rbind(original_df, univariate_extreme_df, fitted_mv_extreme_df)
    combined_extreme_df$type = factor(combined_extreme_df$type, 
                                    levels = c("Original", "Univariate Extreme", "Multivariate Extreme"))

    p_pairs_extreme = ggpairs(combined_extreme_df, columns=1:3, aes(color=type, alpha=0.7),
                    upper=list(continuous="cor"),
                    lower=list(continuous="points"),
                    diag=list(continuous="densityDiag")) +
            scale_color_manual(values = c("Original" = "blue", 
                                        "Univariate Extreme" = "orange", 
                                        "Multivariate Extreme" = "darkred")) +
            scale_fill_manual(values = c("Original" = "blue", 
                                    "Univariate Extreme" = "orange", 
                                    "Multivariate Extreme" = "darkred")) +
            theme(
                    text = element_text(size = 20),
                    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                    strip.background = element_rect(fill = "white"),
                    legend.key = element_rect(fill = "white")
                )

    ggsave(filename="univariate_vs_multivariate_extreme.png", plot=p_pairs_extreme, path=plot_dir,
        width=12, height=10, units="in", dpi=300)
