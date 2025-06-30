quick_diagnostics <- function(fit) {
  cat("=== STAN FIT DIAGNOSTICS ===\n")
  
  # Basic info
  cat("Chains:", fit@sim$chains, "\n")
  cat("Iterations:", fit@sim$iter, "\n")
  cat("Warmup:", fit@sim$warmup, "\n")
  
  # Convergence
  rhats <- summary(fit)$summary[,"Rhat"]
  cat("Max Rhat:", round(max(rhats, na.rm = TRUE), 4), "\n")
  cat("Parameters with Rhat > 1.01:", sum(rhats > 1.01, na.rm = TRUE), "\n")
  
  # Effective sample size
  n_eff <- summary(fit)$summary[,"n_eff"]
  cat("Min n_eff:", round(min(n_eff, na.rm = TRUE), 0), "\n")
  
  # Divergences
  div <- sum(get_divergent_iterations(fit))
  cat("Divergent transitions:", div, "\n")

  # Max tree depth
  tree <- sum(get_max_treedepth_iterations(fit))
  cat("Transitions exceeding max treedepth:", tree, "\n")
  
  # Energy diagnostics
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  energy <- sapply(sampler_params, function(x) x[,"energy__"])
  cat("Energy diagnostics check:\n")
  check_energy(fit)
  
  cat("========================\n")
}

compare_marginals <- function(fit1, fit2, params = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev")) {
  for(param in params) {
    samples1 <- extract(fit1, param)[[1]]
    samples2 <- extract(fit2, param)[[1]]
    
    cat("\n", param, ":\n")
    cat("  Model 1 (ppc) - Mean:", round(mean(samples1), 3), 
        "SD:", round(sd(samples1), 3), "\n")
    cat("  Model 2 (fit) - Mean:", round(mean(samples2), 3), 
        "SD:", round(sd(samples2), 3), "\n")
    cat("  Difference in means:", round(mean(samples2) - mean(samples1), 3), "\n")
  }
}

generate_prior_samples <- function(stan_data, n_samples = 1000) {
  
  # Extract prior parameters from stan_data
  mv_mean <- stan_data$mv_prior_mean
  mv_sd <- stan_data$mv_prior_sd  
  mv_corr <- stan_data$mv_prior_corr
  
  qeff_meanlog <- stan_data$prior_qeff_meanlog
  qeff_sdlog <- stan_data$prior_qeff_sdlog
  
  mv_qdev_mean <- stan_data$mv_qdev_prior_mean
  mv_qdev_sd <- stan_data$mv_qdev_prior_sd
  mv_qdev_corr <- stan_data$mv_qdev_prior_corr
  
  sigmap_meanlog <- stan_data$PriorMean_logsigmap
  sigmap_sdlog <- stan_data$PriorSD_logsigmap
  sigmao_add_sd <- stan_data$PriorSD_sigmao_add
  
  # Generate multivariate samples for logK, log_r, log_shape
  mv_cov <- diag(mv_sd) %*% mv_corr %*% diag(mv_sd)
  mv_samples <- MASS::mvrnorm(n_samples, mu = mv_mean, Sigma = mv_cov)
  
  # Transform to natural scale
  logK <- mv_samples[,1]  # already on log scale
  r <- exp(mv_samples[,2])  # transform from log scale
  shape <- exp(mv_samples[,3])  # transform from log scale
  
  # Generate qeff samples (independent lognormal)
  qeff <- rlnorm(n_samples, meanlog = qeff_meanlog, sdlog = qeff_sdlog)
  
  # Generate bivariate samples for rho and sigma_qdev (on transformed scale)
  qdev_cov <- diag(mv_qdev_sd) %*% mv_qdev_corr %*% diag(mv_qdev_sd)
  qdev_samples <- MASS::mvrnorm(n_samples, mu = mv_qdev_mean, Sigma = qdev_cov)
  
  # Transform to natural scale
  rho <- tanh(qdev_samples[,1])  # transform from atanh scale
  sigma_qdev <- exp(qdev_samples[,2])  # transform from log scale
  
  # Generate other parameter samples
  sigmap <- exp(rnorm(n_samples, mean = sigmap_meanlog, sd = sigmap_sdlog))
  sigmao_add <- abs(rnorm(n_samples, mean = 0, sd = sigmao_add_sd))
  
  # Return as a list mimicking rstan::extract() output
  prior_samples <- list(
    logK = logK,
    r = r, 
    shape = shape,
    sigmap = sigmap,
    sigmao_add = sigmao_add,
    qeff = qeff,
    rho = rho,
    sigma_qdev = sigma_qdev
  )
  
  return(prior_samples)
}

compare_marginals_pairs <- function(fit1, fit2, prior_samples = NULL, params = c("logK", "r", "shape", "sigmap", "sigmao_add", "qeff", "rho", "sigma_qdev"), type_levels = NULL) {
  
  # Extract samples from both fits
  samples_list <- list()
  for(param in params) {
    samples1 <- extract(fit1, param)[[1]]
    samples2 <- extract(fit2, param)[[1]]
    
    # Take same number of samples from each (use minimum)
    n_samples <- min(length(samples1), length(samples2))
    
    if (!is.null(prior_samples)) {
      # Include prior samples if provided
      prior_param <- prior_samples[[param]]
      n_prior <- min(length(prior_param), n_samples)
      
      samples_list[[param]] <- c(prior_param[1:n_prior], 
                                samples1[1:n_prior], 
                                samples2[1:n_prior])
      
      # Set default type levels if not provided
      if (is.null(type_levels)) {
        type_levels <- c("prior", "fit1", "fit2")
      } else if (length(type_levels) != 3) {
        stop("When prior_samples are provided, type_levels must have 3 elements")
      }
      
      type_labels <- c(rep(type_levels[1], n_prior), 
                      rep(type_levels[2], n_prior), 
                      rep(type_levels[3], n_prior))
      
    } else {
      # Original behavior without priors
      samples_list[[param]] <- c(samples1[1:n_samples], samples2[1:n_samples])
      
      # Set default type levels if not provided
      if (is.null(type_levels)) {
        type_levels <- c("fit1", "fit2")
      } else if (length(type_levels) != 2) {
        stop("When prior_samples are not provided, type_levels must have 2 elements")
      }
      
      type_labels <- c(rep(type_levels[1], n_samples), 
                      rep(type_levels[2], n_samples))
    }
  }
  
  # Create data frame
  df <- data.frame(samples_list)
  df$type <- factor(type_labels, levels = type_levels)
  
  # Create color mapping
  if (length(type_levels) == 3) {
    colors <- c("gray", "blue", "darkred")
  } else {
    colors <- c("blue", "darkred")
  }
  names(colors) <- type_levels
  
  # Create pairs plot using ggpairs' mapping approach
  p <- ggpairs(df, 
               columns = 1:length(params), 
               mapping = aes(color = type,fill=type, alpha = 0.7),
               upper = list(continuous = "cor"),
               lower = list(continuous = "points"), 
               diag = list(continuous = "densityDiag")) +
    theme(text = element_text(size = 12),
          panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # Apply color scale after creating the plot
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol) {
      if (!is.null(p[i,j])) {
        p[i,j] <- p[i,j] + 
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors)
      }
    }
  }
  
  return(p)
}
