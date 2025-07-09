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

fit_multivariate_prior <- function(param_matrix, param_names = NULL, use_mle = TRUE, 
                                   ridge_regularization = 1e-8, verbose = FALSE) {
    #' Fit multivariate normal distribution to parameter samples
    #' 
    #' @param param_matrix Matrix of N_samples x N_parameters 
    #' @param param_names Optional vector of parameter names for output labeling
    #' @param use_mle Logical, whether to use MLE fitting (TRUE) or sample moments (FALSE)
    #' @param ridge_regularization Numeric, regularization added to diagonal for numerical stability
    #' @param verbose Logical, whether to print convergence information
    #' 
    #' @return List containing:
    #'   - mean: Vector of fitted means
    #'   - cov: Fitted covariance matrix  
    #'   - cor: Fitted correlation matrix
    #'   - sd: Vector of standard deviations (diagonal of cov matrix)
    #'   - method: String indicating fitting method used
    #'   - convergence: Convergence information (if MLE used)
    #'   - ridge_applied: Logical, whether ridge regularization was applied
    
    # Input validation
    if(!is.matrix(param_matrix)) {
        param_matrix <- as.matrix(param_matrix)
    }
    
    n_samples <- nrow(param_matrix)
    n_params <- ncol(param_matrix)
    
    if(n_samples < n_params + 1) {
        stop("Need at least n_parameters + 1 samples for covariance estimation")
    }
    
    # Set parameter names if not provided
    if(is.null(param_names)) {
        param_names <- paste0("param_", 1:n_params)
    } else if(length(param_names) != n_params) {
        stop("Length of param_names must equal number of parameters")
    }
    
    colnames(param_matrix) <- param_names
    
    # Check for missing values
    if(any(is.na(param_matrix))) {
        warning("Missing values detected. Removing incomplete cases.")
        param_matrix <- param_matrix[complete.cases(param_matrix), , drop = FALSE]
        n_samples <- nrow(param_matrix)
    }
    
    # Sample moments (always computed as fallback)
    sample_mean <- colMeans(param_matrix)
    sample_cov <- cov(param_matrix)
    
    # Check for numerical issues and apply ridge regularization if needed
    ridge_applied <- FALSE
    cov_det <- det(sample_cov)
    min_eigenval <- min(eigen(sample_cov, only.values = TRUE)$values)
    
    if(cov_det < 1e-10 || min_eigenval < ridge_regularization) {
        if(verbose) {
            cat("Applying ridge regularization. Det =", cov_det, 
                ", Min eigenvalue =", min_eigenval, "\n")
        }
        warning("Sample covariance matrix is near-singular or ill-conditioned. Applying ridge regularization.")
        sample_cov <- sample_cov + diag(ridge_regularization, n_params)
        ridge_applied <- TRUE
    }
    
    sample_cor <- cov2cor(sample_cov)
    sample_sd <- sqrt(diag(sample_cov))
    
    # Add names
    names(sample_mean) <- param_names
    names(sample_sd) <- param_names
    dimnames(sample_cov) <- list(param_names, param_names)
    dimnames(sample_cor) <- list(param_names, param_names)
    
    # If not using MLE, return sample moments
    if(!use_mle) {
        return(list(
            mean = sample_mean,
            cov = sample_cov,
            cor = sample_cor,
            sd = sample_sd,
            method = "sample_moments",
            convergence = NA,
            ridge_applied = ridge_applied
        ))
    }
    
    # MLE fitting using Cholesky parameterization
    if(verbose) cat("Fitting multivariate normal distribution via MLE...\n")
    
    # Load required library
    if(!requireNamespace("mvtnorm", quietly = TRUE)) {
        stop("Package 'mvtnorm' is required for MLE fitting")
    }
    
    # Objective function for MLE
    mv_mle <- function(par) {
        mu <- par[1:n_params]
        
        # Cholesky parameterization for positive definiteness
        n_chol_params <- n_params * (n_params + 1) / 2
        L_vec <- par[(n_params + 1):(n_params + n_chol_params)]
        
        # Build lower triangular matrix
        L <- matrix(0, n_params, n_params)
        L[lower.tri(L, diag = TRUE)] <- L_vec
        
        # Reconstruct covariance matrix
        Sigma <- L %*% t(L)
        
        # Check for numerical issues
        if(any(!is.finite(Sigma)) || det(Sigma) <= 1e-12) {
            return(1e10)
        }
        
        # Compute negative log-likelihood
        tryCatch({
            -sum(mvtnorm::dmvnorm(param_matrix, mu, Sigma, log = TRUE))
        }, error = function(e) {
            return(1e10)
        })
    }
    
    # Initialize with sample moments (potentially regularized)
    init_chol <- tryCatch({
        chol(sample_cov)  # Upper triangular
    }, error = function(e) {
        # If Cholesky fails, add more ridge regularization
        if(verbose) cat("Initial Cholesky failed, adding more regularization\n")
        regularized_cov <- sample_cov + diag(ridge_regularization * 100, n_params)
        chol(regularized_cov)
    })
    
    init_chol_lower <- t(init_chol)  # Convert to lower triangular
    init_params <- c(sample_mean, 
                    init_chol_lower[lower.tri(init_chol_lower, diag = TRUE)])
    
    # Fit with robust settings
    fit_result <- tryCatch({
        nlminb(init_params, mv_mle,
               control = list(eval.max = 1000, iter.max = 1000, trace = ifelse(verbose, 1, 0)))
    }, error = function(e) {
        if(verbose) cat("MLE fitting failed:", e$message, "\n")
        list(convergence = -1)
    })
    
    # Check convergence
    converged <- (fit_result$convergence == 0 || fit_result$convergence == 4)
    
    if(!converged) {
        if(verbose) {
            cat("MLE fitting did not converge (code:", fit_result$convergence, 
                "). Using sample moments.\n")
        }
        
        return(list(
            mean = sample_mean,
            cov = sample_cov,
            cor = sample_cor,
            sd = sample_sd,
            method = "sample_moments_fallback",
            convergence = fit_result$convergence,
            ridge_applied = ridge_applied
        ))
    }
    
    # Extract fitted parameters
    fitted_mean <- fit_result$par[1:n_params]
    n_chol_params <- n_params * (n_params + 1) / 2
    L_vec <- fit_result$par[(n_params + 1):(n_params + n_chol_params)]
    
    # Reconstruct covariance matrix
    L_fitted <- matrix(0, n_params, n_params)
    L_fitted[lower.tri(L_fitted, diag = TRUE)] <- L_vec
    fitted_cov <- L_fitted %*% t(L_fitted)
    
    # Calculate derived quantities
    fitted_cor <- cov2cor(fitted_cov)
    fitted_sd <- sqrt(diag(fitted_cov))
    
    # Add names
    names(fitted_mean) <- param_names
    names(fitted_sd) <- param_names
    dimnames(fitted_cov) <- list(param_names, param_names)
    dimnames(fitted_cor) <- list(param_names, param_names)
    
    if(verbose) {
        cat("MLE fitting converged successfully.\n")
        cat("Log-likelihood:", -fit_result$objective, "\n")
    }
    
    return(list(
        mean = fitted_mean,
        cov = fitted_cov,
        cor = fitted_cor,
        sd = fitted_sd,
        method = "mle",
        convergence = fit_result$convergence,
        ridge_applied = ridge_applied
    ))
}

