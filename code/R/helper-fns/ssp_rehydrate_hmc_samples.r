# HMC Samples Expansion Function
# Reconstructs full hmc_samples from condensed parameters using Stan model equations

library(data.table)

#' Expand condensed HMC samples by reconstructing transformed parameters and generated quantities
#' @param model_dir Directory containing hmc_samples_compact.rda, fit_summary.csv, stan_data.csv
#' @return Full hmc_samples data.table
expand_hmc_samples <- function(model_dir) {
  
  # Load required files
  compact_path <- file.path(model_dir, "hmc_samples_compact.rda")
  summary_path <- file.path(model_dir, "fit_summary.csv")
  data_path <- file.path(model_dir, "stan_data.csv")
  
  if (!all(file.exists(compact_path, summary_path, data_path))) {
    stop("Missing required files in ", model_dir)
  }
  
  load(compact_path) # loads condensed_samples
  fit_summary <- fread(summary_path)
  stan_data <- fread(data_path)
  
  # Get model type and dispatch to appropriate expansion function
  exec_name <- fit_summary$exec[1]
  
  if (grepl("bspm_estq_softdep_mvprior", exec_name)) {
    return(expand_bspm_estq_softdep_mvprior(condensed_samples, stan_data))
  } else {
    stop("Expansion not implemented for model: ", exec_name)
  }
}

#' Expand bspm_estq_softdep_mvprior model
expand_bspm_estq_softdep_mvprior <- function(condensed_samples, stan_data) {
  
  # Extract data block parameters
  T_val <- stan_data[name == "T"]$value
  I_val <- stan_data[name == "I"]$value
  n_periods <- stan_data[name == "n_periods"]$value
  
  # Extract data arrays
  effort <- stan_data[name == "effort"]$value
  index <- get_stan_matrix(stan_data, "index", T_val, I_val)
  obs_removals <- stan_data[name == "obs_removals"]$value
  sigmao_mat <- get_stan_matrix(stan_data, "sigmao_mat", T_val, I_val)
  lambdas <- stan_data[name == "lambdas"]$value
  
  # Extract prior parameters for transformations
  mv_prior_mean <- get_stan_vector(stan_data, "mv_prior_mean", 3)
  mv_prior_sd <- get_stan_vector(stan_data, "mv_prior_sd", 3)
  mv_prior_corr <- get_stan_matrix(stan_data, "mv_prior_corr", 3, 3)
  
  mv_qdev_prior_mean <- get_stan_vector(stan_data, "mv_qdev_prior_mean", 2)
  mv_qdev_prior_sd <- get_stan_vector(stan_data, "mv_qdev_prior_sd", 2)
  mv_qdev_prior_corr <- get_stan_matrix(stan_data, "mv_qdev_prior_corr", 2, 2)
  
  prior_qeff_meanlog <- stan_data[name == "prior_qeff_meanlog"]$value
  prior_qeff_sdlog <- stan_data[name == "prior_qeff_sdlog"]$value
  PriorMean_logsigmap <- stan_data[name == "PriorMean_logsigmap"]$value
  PriorSD_logsigmap <- stan_data[name == "PriorSD_logsigmap"]$value
  PriorSD_sigmao_add <- stan_data[name == "PriorSD_sigmao_add"]$value
  sigma_edev <- stan_data[name == "sigma_edev"]$value
  
  # Get unique iterations for processing
  unique_iters <- unique(condensed_samples[, .(iter, chain, chain_iter)])
  
  # Process each iteration
  expanded_list <- list()
  
  for (i in 1:nrow(unique_iters)) {
    iter_info <- unique_iters[i]
    iter_samples <- condensed_samples[iter == iter_info$iter & chain == iter_info$chain]
    
    # Extract parameters for this iteration
    raw_mv_params <- get_parameter_vector(iter_samples, "raw_mv_params", 3)
    raw_logqeff <- get_parameter_scalar(iter_samples, "raw_logqeff")
    raw_qdev_params <- get_parameter_vector(iter_samples, "raw_qdev_params", 2)
    raw_logsigmap <- get_parameter_scalar(iter_samples, "raw_logsigmap")
    raw_sigmao_add <- get_parameter_scalar(iter_samples, "raw_sigmao_add")
    raw_epsp <- get_parameter_vector(iter_samples, "raw_epsp", T_val)
    raw_qdev_period <- get_parameter_vector(iter_samples, "raw_qdev_period", n_periods)
    raw_edev <- get_parameter_vector(iter_samples, "raw_edev", T_val - 1)
    
    # Compute transformed parameters following Stan model
    expanded_vars <- compute_transformed_parameters(
      raw_mv_params, raw_logqeff, raw_qdev_params, raw_logsigmap, raw_sigmao_add,
      raw_epsp, raw_qdev_period, raw_edev,
      mv_prior_mean, mv_prior_sd, mv_prior_corr,
      mv_qdev_prior_mean, mv_qdev_prior_sd, mv_qdev_prior_corr,
      prior_qeff_meanlog, prior_qeff_sdlog,
      PriorMean_logsigmap, PriorSD_logsigmap, PriorSD_sigmao_add,
      sigma_edev, effort, T_val, n_periods
    )
    
    # Compute generated quantities
    generated_vars <- compute_generated_quantities(
      expanded_vars, index, sigmao_mat, lambdas, I_val, T_val
    )
    
    # Combine all variables for this iteration
    all_vars <- c(expanded_vars, generated_vars)
    
    # Convert to data.table format
    iter_dt <- convert_to_hmc_format(all_vars, iter_info, condensed_samples)
    expanded_list[[i]] <- iter_dt
  }
  
  # Combine all iterations
  full_samples <- rbind(condensed_samples, rbindlist(expanded_list, fill = TRUE))
  
  return(full_samples)
}

#' Compute transformed parameters following bspm_estq_softdep_mvprior.stan
compute_transformed_parameters <- function(raw_mv_params, raw_logqeff, raw_qdev_params, 
                                         raw_logsigmap, raw_sigmao_add, raw_epsp, 
                                         raw_qdev_period, raw_edev,
                                         mv_prior_mean, mv_prior_sd, mv_prior_corr,
                                         mv_qdev_prior_mean, mv_qdev_prior_sd, mv_qdev_prior_corr,
                                         prior_qeff_meanlog, prior_qeff_sdlog,
                                         PriorMean_logsigmap, PriorSD_logsigmap, PriorSD_sigmao_add,
                                         sigma_edev, effort, T_val, n_periods) {
  
  # Multivariate transformation for logK, r, shape
  L_mv <- chol(mv_prior_corr)
  mv_params <- mv_prior_mean + diag(mv_prior_sd) %*% L_mv %*% raw_mv_params
  
  logK <- mv_params[1]
  r <- exp(mv_params[2])
  shape <- exp(mv_params[3])
  
  # Independent qeff parameter
  logqeff <- raw_logqeff * prior_qeff_sdlog + prior_qeff_meanlog
  qeff <- exp(logqeff)
  
  # Bivariate transformation for rho and sigma_qdev
  L_qdev <- chol(mv_qdev_prior_corr)
  qdev_params <- mv_qdev_prior_mean + diag(mv_qdev_prior_sd) %*% L_qdev %*% raw_qdev_params
  
  rho <- tanh(qdev_params[1])
  sigma_qdev <- exp(qdev_params[2])
  
  # Other parameters
  sigmap <- exp(raw_logsigmap * PriorSD_logsigmap + PriorMean_logsigmap)
  sigmao_add <- raw_sigmao_add * PriorSD_sigmao_add
  
  # Process error deviations
  dev <- numeric(T_val)
  epsp <- numeric(T_val)
  for (t in 1:T_val) {
    dev[t] <- raw_epsp[t] * sigmap
    epsp[t] <- exp(dev[t])
  }
  
  # Effort-based parameters with AR(1) structure
  qdev_period <- numeric(n_periods)
  qdev_period[1] <- raw_qdev_period[1] * sigma_qdev
  for (p in 2:n_periods) {
    qdev_period[p] <- rho * qdev_period[p-1] + raw_qdev_period[p] * sigma_qdev * sqrt(1 - rho^2)
  }
  
  # Map periods to time steps (assuming equal periods)
  n_step <- floor((T_val - 1) / n_periods)
  qdev <- numeric(T_val - 1)
  for (t in 1:(T_val - 1)) {
    period <- min(ceiling(t / n_step), n_periods)
    qdev[t] <- qdev_period[period]
  }
  
  # Annual effort deviations
  edev <- raw_edev * sigma_edev
  
  # Calculate fishing mortality
  F_vec <- numeric(T_val - 1)
  for (t in 1:(T_val - 1)) {
    F_vec[t] <- qeff * exp(qdev[t]) * effort[t] * exp(edev[t] - sigma_edev^2/2)
  }
  
  # Derived biological parameters
  n <- shape
  dmsy <- (1/n)^(1/(n-1))
  h <- 2 * dmsy
  m <- r * h / 4
  g <- n^(n/(n-1)) / (n-1)
  
  # Population dynamics
  x <- numeric(T_val)
  x[1] <- epsp[1]
  removals <- numeric(T_val - 1)
  
  for (t in 2:T_val) {
    if (x[t-1] <= dmsy) {
      x[t] <- (x[t-1] + r * x[t-1] * (1 - x[t-1]/h)) * exp(-F_vec[t-1]) * epsp[t]
      removals[t-1] <- (x[t-1] + r * x[t-1] * (1 - x[t-1]/h)) * epsp[t] * (1 - exp(-F_vec[t-1])) * exp(logK)
    } else {
      x[t] <- (x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1))) * exp(-F_vec[t-1]) * epsp[t]
      removals[t-1] <- (x[t-1] + g * m * x[t-1] * (1 - x[t-1]^(n-1))) * epsp[t] * (1 - exp(-F_vec[t-1])) * exp(logK)
    }
  }
  
  return(list(
    logK = logK, r = r, shape = shape, logqeff = logqeff, qeff = qeff,
    rho = rho, sigma_qdev = sigma_qdev, sigmap = sigmap, sigmao_add = sigmao_add,
    dev = dev, epsp = epsp, qdev_period = qdev_period, qdev = qdev, edev = edev,
    F = F_vec, n = n, dmsy = dmsy, h = h, m = m, g = g, x = x, removals = removals,
    # Backward compatibility raw parameters
    raw_logK = raw_mv_params[1], raw_logr = raw_mv_params[2], raw_logshape = raw_mv_params[3],
    raw_rho = raw_qdev_params[1], raw_sigma_qdev = raw_qdev_params[2]
  ))
}

#' Compute generated quantities (analytical catchability)
compute_generated_quantities <- function(params, index, sigmao_mat, lambdas, I_val, T_val) {
  
  # Analytical catchability calculation
  q <- numeric(I_val)
  sigmao <- sigmao_mat + params$sigmao_add
  sigmao2 <- sigmao^2
  
  for (i in 1:I_val) {
    sum1 <- 0
    sum2 <- 0
    p <- 0
    
    for (t in 1:T_val) {
      if (index[t, i] > 0 && params$x[t] > 0) {
        sum1 <- sum1 + log(index[t, i] / params$x[t]) / sigmao2[t, i]
        sum2 <- sum2 + 1 / sigmao2[t, i]
        p <- p + 1
      }
    }
    
    if (p > 2) {
      q[i] <- exp((0.5 * p + sum1) / sum2)
    } else {
      q[i] <- 0
    }
  }
  
  return(list(q = q, sigmao = sigmao, sigmao2 = sigmao2))
}

# Helper functions
get_stan_matrix <- function(stan_data, name, nrow, ncol) {
  data <- stan_data[name == name]
  matrix(data$value, nrow = nrow, ncol = ncol, byrow = TRUE)
}

get_stan_vector <- function(stan_data, name, length) {
  stan_data[name == name]$value[1:length]
}

get_parameter_vector <- function(iter_samples, name, length) {
  iter_samples[name == name][order(row)]$value[1:length]
}

get_parameter_scalar <- function(iter_samples, name) {
  iter_samples[name == name]$value[1]
}

convert_to_hmc_format <- function(vars, iter_info, template) {
  # Convert computed variables back to hmc_samples format
  # This creates the same structure as the original condensed_samples
  rows_list <- list()
  
  # Add scalar variables
  scalars <- c("logK", "r", "shape", "logqeff", "qeff", "rho", "sigma_qdev", 
               "sigmap", "sigmao_add", "n", "dmsy", "h", "m", "g",
               "raw_logK", "raw_logr", "raw_logshape", "raw_rho", "raw_sigma_qdev")
  
  for (var in scalars) {
    if (var %in% names(vars)) {
      rows_list[[length(rows_list) + 1]] <- data.table(
        run_id = template$run_id[1],
        iter = iter_info$iter,
        chain = iter_info$chain,
        chain_iter = iter_info$chain_iter,
        variable = var,
        name = var,
        row = as.numeric(NA),
        col = as.numeric(NA),
        value = vars[[var]],
        treedepth = template[iter == iter_info$iter & chain == iter_info$chain]$treedepth[1],
        divergent = template[iter == iter_info$iter & chain == iter_info$chain]$divergent[1],
        acceptance = template[iter == iter_info$iter & chain == iter_info$chain]$acceptance[1],
        stepsize = template[iter == iter_info$iter & chain == iter_info$chain]$stepsize[1],
        leapfrog = template[iter == iter_info$iter & chain == iter_info$chain]$leapfrog[1],
        energy = template[iter == iter_info$iter & chain == iter_info$chain]$energy[1]
      )
    }
  }
  
  # Add vector variables
  vectors <- list(
    list(name = "dev", var = "dev"),
    list(name = "epsp", var = "epsp"),
    list(name = "qdev_period", var = "qdev_period"),
    list(name = "qdev", var = "qdev"),
    list(name = "edev", var = "edev"),
    list(name = "F", var = "F"),
    list(name = "x", var = "x"),
    list(name = "removals", var = "removals"),
    list(name = "q", var = "q")
  )
  
  for (vec_info in vectors) {
    if (vec_info$var %in% names(vars)) {
      vec_data <- vars[[vec_info$var]]
      for (j in 1:length(vec_data)) {
        rows_list[[length(rows_list) + 1]] <- data.table(
          run_id = template$run_id[1],
          iter = iter_info$iter,
          chain = iter_info$chain,
          chain_iter = iter_info$chain_iter,
          variable = paste0(vec_info$name, "[", j, "]"),
          name = vec_info$name,
          row = j,
          col = as.numeric(NA),
          value = vec_data[j],
          treedepth = template[iter == iter_info$iter & chain == iter_info$chain]$treedepth[1],
          divergent = template[iter == iter_info$iter & chain == iter_info$chain]$divergent[1],
          acceptance = template[iter == iter_info$iter & chain == iter_info$chain]$acceptance[1],
          stepsize = template[iter == iter_info$iter & chain == iter_info$chain]$stepsize[1],
          leapfrog = template[iter == iter_info$iter & chain == iter_info$chain]$leapfrog[1],
          energy = template[iter == iter_info$iter & chain == iter_info$chain]$energy[1]
        )
      }
    }
  }
  
  return(rbindlist(rows_list))
}
