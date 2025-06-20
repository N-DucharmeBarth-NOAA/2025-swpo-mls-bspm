// Optimized fletcher-schaefer surplus production model with effort-based fishing mortality
// Applied 10+ key optimizations for performance and numerical stability
// Vectorized where possible for maximum performance

data {
    int T; // time dimension
    int I; // number of indices
    real index[T,I]; // matrix of indices
    real obs_removals[T]; // time series of catch
    real sigmao_mat[T,I]; // observation error corresponding to each index (scaled to mean 1)
    real sigmao_input; // observation error
    int lambdas[I]; // lambdas for CPUE
    real sigmac; // observation error for catch
    real sigma_edev; // effort dev variation
    int t_dep; // year of depletion prior
    int use_depletion_prior; // boolean flag (0/1) to control depletion prior
    int fit_to_data; // boolean flag (0/1) to control if likelihoods are turned on
    real<lower=0> epsilon; // numerical tolerance

    // New effort-based parameters  
    real effort[T]; // effort time series
    int n_step; // years per period (e.g., 3-5 years)
    int n_periods; // number of catchability periods

    // Updated multivariate priors (now 3-dimensional: logK, log_r, log_shape)
    vector[3] mv_prior_mean; // mean vector [logK, log_r, log_shape]
    vector[3] mv_prior_sd; // standard deviations 
    corr_matrix[3] mv_prior_corr; // correlation matrix

    // Independent qeff prior parameters
    real prior_qeff_meanlog; // mean on log scale
    real prior_qeff_sdlog; // SD on log scale

    // Bivariate prior for rho and sigma_qdev (on transformed scale)
    vector[2] mv_qdev_prior_mean; // [atanh_rho_mean, log_sigma_qdev_mean]
    vector[2] mv_qdev_prior_sd; // standard deviations
    corr_matrix[2] mv_qdev_prior_corr; // correlation matrix

    // Other priors
    real PriorMean_logsigmap;
    real PriorSD_logsigmap;
    real PriorSD_sigmao_add;
    real prior_depletion_meanlog; // prior mean on log scale
    real prior_depletion_sdlog; // prior SD on log scale
}

transformed data {
    int Tm1 = T - 1;
    // Pre-compute Cholesky decompositions
    matrix[3,3] mv_prior_chol = cholesky_decompose(mv_prior_corr);
    matrix[2,2] mv_qdev_prior_chol = cholesky_decompose(mv_qdev_prior_corr);
    
    // Pre-compute observation error base
    real sigmao_base[T,I];
    for(t in 1:T) {
        for(i in 1:I) {
            sigmao_base[t,i] = sigmao_mat[t,i] * sigmao_input;
        }
    }
}

parameters {
    // multivariate parameters 
    vector[3] raw_mv_params; // [raw_logK, raw_log_r, raw_log_shape]

    // Independent qeff parameter
    real raw_logqeff; // raw log effort catchability

    // Bivariate parameters for rho and sigma_qdev
    vector[2] raw_qdev_params; // [raw_atanh_rho, raw_log_sigma_qdev]

    real raw_logsigmap;
    real<lower=0> raw_sigmao_add;
    vector[T] raw_epsp; 
    
    // Effort-based parameters
    vector[n_periods] raw_qdev_period; // raw period-specific deviations (non-centered)
    vector[Tm1] raw_edev; // raw annual effort deviations (non-centered)
}

transformed parameters {
    // Core population parameters
    real logK, r, shape;
    
    // Effort/catchability parameters
    real logqeff, qeff, rho, sigma_qdev;
    
    // Time series vectors 
    vector[T] epsp;
    vector[Tm1] edev;
    vector[Tm1] F;
    vector[T] x;           
    vector[Tm1] removals;  
    
    // Error and observation parameters
    real sigmap, sigmao_add;
    real sigmao[T,I];
    real q[I];
    
    // Intermediate calculation variables
    vector[3] mv_params;
    vector[2] qdev_params;
    vector[n_periods] qdev_period;
    vector[Tm1] qdev;      

    // Use pre-computed Cholesky matrices
    mv_params = mv_prior_mean + diag_pre_multiply(mv_prior_sd, mv_prior_chol) * raw_mv_params;
    
    // Extract individual parameters (transformed scale)
    logK = mv_params[1]; // logK (already on log scale)
    r = exp(mv_params[2]); // r (transform from log scale)
    shape = exp(mv_params[3]); // shape (transform from log scale)
    
    // Independent qeff parameter
    logqeff = raw_logqeff * prior_qeff_sdlog + prior_qeff_meanlog;
    qeff = exp(logqeff);
    
    // Use pre-computed Cholesky matrix for qdev parameters
    qdev_params = mv_qdev_prior_mean + diag_pre_multiply(mv_qdev_prior_sd, mv_qdev_prior_chol) * raw_qdev_params;
    
    // Transform parameters back to natural scale
    rho = tanh(qdev_params[1]); // atanh_rho -> rho
    sigma_qdev = exp(qdev_params[2]); // log_sigma_qdev -> sigma_qdev
    
    // Non-centered parameterization for AR(1) catchability process
    qdev_period[1] = raw_qdev_period[1] * sigma_qdev;
    for(p in 2:n_periods) {
        qdev_period[p] = rho * qdev_period[p-1] + raw_qdev_period[p] * sigma_qdev * sqrt(1 - rho^2);
    }
    
    // Safer period assignment with bounds checking
    for(t in 1:Tm1) {
        int period = min(max(((t-1) / n_step) + 1, 1), n_periods);
        qdev[t] = qdev_period[period];
    }
    
    // effort deviations
    edev = raw_edev * sigma_edev;
    
    // fishing mortality calculation
    for(t in 1:Tm1) {
        F[t] = qeff * exp(qdev[t]) * effort[t] * exp(edev[t] - 0.5 * square(sigma_edev));
    }

    // Process error calculations
    sigmap = exp(raw_logsigmap * PriorSD_logsigmap + PriorMean_logsigmap);
    epsp = exp(raw_epsp * sigmap - 0.5 * square(sigmap));

    // observation error calculation
    sigmao_add = raw_sigmao_add * PriorSD_sigmao_add;
    for(t in 1:T) {
        for(i in 1:I) {
            sigmao[t,i] = sigmao_base[t,i] + sigmao_mat[t,i] * sigmao_add;
        }
    }
    
    // Fletcher-schaefer parameters
    real n = shape;
    real dmsy = pow((1/n), (1/(n-1)));
    real h = 2 * dmsy;
    real m = r * h / 4;
    real g = pow(n, (n/(n-1))) / (n-1);

    // Population dynamics
    x[1] = epsp[1];
    for(t in 2:T) {
        if(x[t-1] <= dmsy) {
            x[t] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h))) * epsp[t] * exp(-F[t-1]);
            removals[t-1] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h))) * epsp[t] * (1 - exp(-F[t-1])) * exp(logK);
        } 
        if(x[t-1] > dmsy) {
            x[t] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1], (n-1))))) * epsp[t] * exp(-F[t-1]);
            removals[t-1] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1], (n-1))))) * epsp[t] * (1 - exp(-F[t-1])) * exp(logK);
        } 
    }

    // Stabilized analytical q calculation
    real sigmao2[T,I];
    real sum1, sum2, p;
    
    for(i in 1:I) {
        sum1 = 0.0;
        sum2 = 0.0;
        p = 0.0;  
        for(t in 1:T) {
            sigmao2[t,i] = square(sigmao[t,i]);
            if(index[t,i] > 0.0 && x[t] > 0.0) {
                sum1 = sum1 + log(index[t,i]/x[t]) / sigmao2[t,i];
                sum2 = sum2 + 1 / sigmao2[t,i];
                p = p + 1.0;
            }
        }
        if(p > 2.0 && sum2 > epsilon) {
            q[i] = exp((0.5 * p + sum1) / sum2);
        } else {
            q[i] = 0.0;
        } 
    }
} 

model {
    // Updated multivariate normal prior (3-dimensional)
    raw_mv_params ~ std_normal();

    // Independent qeff prior (lognormal)
    raw_logqeff ~ std_normal();

    // Bivariate prior for rho and sigma_qdev
    raw_qdev_params ~ std_normal();
    
    // Effort-based priors (non-centered)
    raw_qdev_period ~ std_normal();
    raw_edev ~ std_normal();
    
    // Prior densities for other estimated parameters
    raw_epsp ~ std_normal();
    raw_logsigmap ~ std_normal();
    raw_sigmao_add ~ std_normal();
    
    // Depletion prior (conditional)
    if(use_depletion_prior == 1) {
        target += lognormal_lpdf(x[t_dep] | prior_depletion_meanlog, prior_depletion_sdlog);
    }
    
    // Observation model - uses analytical q[i]
    if(fit_to_data == 1) {
        for(i in 1:I) {
            if(lambdas[i] == 1) {
                for(t in 1:T) {
                    if(index[t,i] > 0.0 && x[t] > 0.0 && q[i] > 0.0) {
                        target += lognormal_lpdf(index[t,i] | log(q[i] * x[t]) - sigmao2[t,i] / 2, sigmao[t,i]);
                    }
                }
            }
        }

        // Catch observation model uses F derived from effort-based qeff
        for(t in 1:Tm1) {
            target += lognormal_lpdf(obs_removals[t] | log(removals[t]) - 0.5 * square(sigmac), sigmac);
        }
    }
}

generated quantities {
    // Raw parameter calculations 
    real raw_logK = (mv_params[1] - mv_prior_mean[1]) / mv_prior_sd[1];
    real raw_logr = (mv_params[2] - mv_prior_mean[2]) / mv_prior_sd[2];
    real raw_logshape = (mv_params[3] - mv_prior_mean[3]) / mv_prior_sd[3];
    real raw_rho = (qdev_params[1] - mv_qdev_prior_mean[1]) / mv_qdev_prior_sd[1];
    real raw_sigma_qdev = (qdev_params[2] - mv_qdev_prior_mean[2]) / mv_qdev_prior_sd[2];
}
