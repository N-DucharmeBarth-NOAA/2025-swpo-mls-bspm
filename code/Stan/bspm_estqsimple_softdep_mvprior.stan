// Modified fletcher-schaefer surplus production model with effort-based fishing mortality
// Updated prior structure based on pushforward analysis:
// - 3D multivariate prior for logK, log_r, log_shape  
// - Independent qeff parameter
// - Correlated rho and sigma_qdev parameters

data {
    int T; // time dimension
    int I; // number of indices
    real index[T,I]; // matrix of indices
    real obs_removals[T]; // time series of catch
    real sigmao_mat[T,I]; // observation error corresponding to each index (scaled to mean 1)
    real sigmao_input; // observation error
    int lambdas[I]; // lambdas for CPUE
    real sigmac; // observation error for catch
    int t_dep; // year of depletion prior
    int use_depletion_prior; // boolean flag (0/1) to control depletion prior
    int fit_to_data; // boolean flag (0/1) to control if likelihoods are turned on

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

transformed data{
    int Tm1;
    Tm1 = T-1;
    int Np1;
    Np1 = n_periods-1;
}

parameters {
    // Updated multivariate parameters (now 3-dimensional)
    vector[3] raw_mv_params; // [raw_logK, raw_log_r, raw_log_shape]

    // Independent qeff parameter
    real raw_logqeff; // raw log effort catchability

    // Bivariate parameters for rho and sigma_qdev
    vector[2] raw_qdev_params; // [raw_atanh_rho, raw_log_sigma_qdev]

    real raw_logsigmap;
    real<lower=0> raw_sigmao_add;
    real raw_epsp[T];
    
    // Effort-based parameters
    vector[Np1] raw_qdev_period; // raw period-specific deviations (non-centered)
}

transformed parameters {
    // leading parameters
    real logK;
    real r;
    real shape;
    real logqeff; // effort-based catchability (independent)
    real qeff; // transformed effort catchability
    real dev[T]; // recruitment deviates
    real epsp[T]; // process error (multiplicative)
    real sigmap;
    real rho; // AR correlation for catchability
    real sigma_qdev; // catchability variability

    // For backward compatibility - extract individual raw parameters
    real raw_logK;
    real raw_logr;
    real raw_logshape;
    real raw_rho;
    real raw_sigma_qdev;

    // Updated multivariate transformation (3-dimensional)
    vector[3] mv_params;
    mv_params = mv_prior_mean + diag_pre_multiply(mv_prior_sd, cholesky_decompose(mv_prior_corr)) * raw_mv_params;
    
    // Extract individual parameters (transformed scale)
    logK = mv_params[1]; // logK (already on log scale)
    r = exp(mv_params[2]); // r (transform from log scale)
    shape = exp(mv_params[3]); // shape (transform from log scale)
    
    // Extract individual raw parameters (for compatibility)
    raw_logK = (mv_params[1] - mv_prior_mean[1]) / mv_prior_sd[1];
    raw_logr = (mv_params[2] - mv_prior_mean[2]) / mv_prior_sd[2];
    raw_logshape = (mv_params[3] - mv_prior_mean[3]) / mv_prior_sd[3];
    
    // Independent qeff parameter
    logqeff = raw_logqeff * prior_qeff_sdlog + prior_qeff_meanlog;
    qeff = exp(logqeff);
    
    // Bivariate transformation for rho and sigma_qdev
    vector[2] qdev_params;
    qdev_params = mv_qdev_prior_mean + diag_pre_multiply(mv_qdev_prior_sd, cholesky_decompose(mv_qdev_prior_corr)) * raw_qdev_params;
    
    // Transform parameters back to natural scale
    rho = tanh(qdev_params[1]); // atanh_rho -> rho
    sigma_qdev = exp(qdev_params[2]); // log_sigma_qdev -> sigma_qdev
    
    // Extract individual raw parameters for rho and sigma_qdev (for compatibility)
    raw_rho = (qdev_params[1] - mv_qdev_prior_mean[1]) / mv_qdev_prior_sd[1];
    raw_sigma_qdev = (qdev_params[2] - mv_qdev_prior_mean[2]) / mv_qdev_prior_sd[2];
    
    // Effort-based fishing mortality calculation with dual error structure
    vector[n_periods] qdev_period; // transformed period deviations
    real qdev[Tm1]; // systematic catchability changes
    real F[Tm1];
    
    // Non-centered parameterization for AR(1) catchability process
    qdev_period[1] = 0;
    for(p in 2:n_periods) {
        qdev_period[p] = rho * qdev_period[p-1] + raw_qdev_period[p-1] * sigma_qdev * sqrt(1 - rho^2);
    }
    
    // Assign period-based qdev with temporal structure
    for(t in 1:Tm1) {
        int period = ((t-1) / n_step) + 1; // integer division
        period = min(period, n_periods); // safer bound checking
        qdev[t] = qdev_period[period];
    }
    
    // Calculate fishing mortality: F_t = qeff * exp(qdev_t) * effort_t
    for(t in 1:Tm1) {
        F[t] = qeff * exp(qdev[t]- sigma_qdev^2/2) * effort[t];
    }

    // process error
    sigmap = exp(raw_logsigmap*PriorSD_logsigmap + PriorMean_logsigmap); // lognormal prior
    real sigmap2;
    sigmap2 = pow(sigmap,2);

    // observation error
    real sigmao[T,I];
    real sigmao_sc;
    real sigmao_add;
    sigmao_add = raw_sigmao_add*PriorSD_sigmao_add;
    sigmao_sc = sigmao_input + sigmao_add;
    for(t in 1:T){
        dev[t] = raw_epsp[t]*sigmap;
        epsp[t] = exp(dev[t]-sigmap2/2);
        for(i in 1:I){
            sigmao[t,i] = sigmao_mat[t,i] * sigmao_sc;
        }
    }
    
    // fletcher-schaefer
    // parameters
    real n;
    real dmsy;
    real h;
    real m;
    real g;

    n = shape;
    dmsy = pow((1/n),(1/(n-1)));
    h = 2*dmsy;
    m = r*h/4;
    g = pow(n,(n/(n-1)))/(n-1);

    // population dynamics
    real removals[Tm1]; // catch (relative to logK) in each time step
    real x[T]; // population time series relative to logK
    x[1] = epsp[1];
    for(t in 2:T) {
        if(x[t-1]<=dmsy){
            x[t] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h)))*epsp[t]*(exp(-F[t-1]));
            removals[t-1] = ((x[t-1] + r * x[t-1] * (1 - x[t-1]/h)))*epsp[t]*(1-exp(-F[t-1]))*exp(logK);
        } 
        if(x[t-1]> dmsy){
            x[t] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1],(n-1)))))*epsp[t]*(exp(-F[t-1]));
            removals[t-1] = ((x[t-1] + g * m * x[t-1] * (1 - pow(x[t-1],(n-1)))))*epsp[t]*(1-exp(-F[t-1]))*exp(logK);
        } 
    }

    // Analytical q calculation for observation model (KEPT SEPARATE)
    real q[I]; // catchability for indices
    real sigmao2[T,I]; // observation error for each index
    real sum1;
    real sum2;
    real p;
    for(i in 1:I){
        sum1 = 0.0;
        sum2 = 0.0;
        p = 0.0;
        for(t in 1:T){
            sigmao2[t,i] = square(sigmao[t,i]);
            if(index[t,i]>0.0 && x[t]>0.0) {
                sum1 = sum1 + log(index[t,i]/x[t])/sigmao2[t,i];
                sum2 = sum2 + 1/sigmao2[t,i];
                p = p + 1.0;
            }
        }
        if(p>2.0){
            q[i] = exp((0.5 * p + sum1) / sum2);
        } else{
            q[i] = 0.0;
        } 
    }
} 

model {
    // Updated multivariate normal prior (3-dimensional)
    raw_mv_params ~ std_normal(); // Standard normal for raw parameters

    // Independent qeff prior (lognormal)
    raw_logqeff ~ std_normal(); // Standard normal for raw parameter

    // Bivariate prior for rho and sigma_qdev
    raw_qdev_params ~ std_normal(); // Standard normal for raw parameters
    
    // Effort-based priors (non-centered)
    raw_qdev_period ~ std_normal(); // All catchability deviations are standard normal
    
    // prior densities for other estimated parameters
    raw_epsp ~ std_normal();
    raw_logsigmap ~ std_normal();
    raw_sigmao_add ~ std_normal();
    
    // Lognormal prior depletion prior from Gedamke and Hoenig 2006 size based Z estimate
    if(use_depletion_prior == 1) {
        target += lognormal_lpdf(x[t_dep] | prior_depletion_meanlog, prior_depletion_sdlog);
    }
    
    // observation model - uses analytical q[i]
    real mu_index;
    for(i in 1:I){
        if(lambdas[i]==1){
            for(t in 1:T){
                if(index[t,i]>0.0 && x[t]>0.0 && q[i]>0.0) {
                    mu_index = log(q[i]*x[t]) - sigmao2[t,i]/2;
                    if(fit_to_data == 1) {
                        target += lognormal_lpdf(index[t,i] | mu_index,sigmao[t,i]);
                    }
                }
            }
        }
    }

    // Catch observation model uses F derived from effort-based qeff
    real mu_catch;
    for(t in 1:Tm1){
        mu_catch = log(removals[t]) - 0.5*sigmac^2;
        if(fit_to_data == 1) {
            target += lognormal_lpdf(obs_removals[t]|mu_catch,sigmac);
        }
    }
}
