// Modified fletcher-schaefer surplus production model with effort-based fishing mortality
// Original code from Bayesian biomass dynamic model (https://github.com/cttedwards/bdm/) by Charles Edwards
// This version uses effort-based F structure with dual error components:
// - qdev: systematic catchability changes (NOT bias corrected)
// - edev: effort measurement errors (bias corrected)

data {
    int T; // time dimension
    int I; // number of indices
    real index[T,I]; // matrix of indices
    real obs_removals[T]; // time series of catch
    real sigmao_mat[T,I]; // observation error corresponding to each index (scaled to mean 1)
    real sigmao_input; // observation error
    int lambdas[I]; // lambdas for CPUE
    real sigmac; // observation error for catch
    real sigma_edev; // effort deviate variance
    int t_dep; // year of depletion prior

    // New effort-based parameters  
    real effort[T]; // effort time series
    int n_step; // years per period (e.g., 3-5 years)
    int n_periods; // number of catchability periods

    // Updated multivariate priors (now 4-dimensional)
    vector[4] mv_prior_mean; // mean vector [logK, log_r, log_n, log_qeff]
    vector[4] mv_prior_sd; // standard deviations 
    corr_matrix[4] mv_prior_corr; // correlation matrix

    // Other priors
    real PriorMean_logsigmap;
    real PriorSD_logsigmap;
    real PriorSD_sigmao_add;
    real prior_sigma_qdev_sd; // prior SD for qdev variability
    real prior_rho_mean; // prior mean for rho
    real prior_rho_sd; // prior SD for rho
    real prior_depletion_meanlog; // prior mean on log scale
    real prior_depletion_sdlog; // prior SD on log scale
}

transformed data{
    int Tm1;
    Tm1 = T-1;
}

parameters {
    // Updated multivariate parameters (now 4-dimensional)
    vector[4] raw_mv_params; // [raw_logK, raw_log_r, raw_log_n, raw_log_qeff]

    real raw_logsigmap;
    real<lower=0> raw_sigmao_add;
    real raw_epsp[T];
    
    // Effort-based parameters
    real raw_sigma_qdev; // qdev variability
    real raw_sigma_edev; // edev variability
    real raw_rho; // AR correlation (to be transformed)
    vector[n_periods] raw_qdev_period; // raw period-specific deviations (non-centered)
    vector[Tm1] raw_edev; // raw annual effort deviations (non-centered)
}

transformed parameters {
    // leading parameters
    real logK;
    real r;
    real shape;
    real logqeff; // effort-based catchability (from multivariate prior)
    real qeff; // transformed effort catchability
    real dev[T]; // recruitment deviates
    real epsp[T]; // process error (multiplicative)
    real sigmap;
    real rho; // transformed AR correlation

    // For backward compatibility - extract individual raw parameters
    real raw_logK;
    real raw_logr;
    real raw_logshape;
    real raw_logqeff; // effort catchability

    // Updated multivariate transformation (4-dimensional)
    vector[4] mv_params;
    mv_params = mv_prior_mean + diag_pre_multiply(mv_prior_sd, cholesky_decompose(mv_prior_corr)) * raw_mv_params;
    
    // Extract individual parameters (transformed scale)
    logK = mv_params[1]; // logK (already on log scale)
    r = exp(mv_params[2]); // r (transform from log scale)
    shape = exp(mv_params[3]); // shape (transform from log scale)
    logqeff = mv_params[4]; // log effort catchability (already on log scale)
    qeff = exp(logqeff); // effort catchability (transform from log scale)
    
    // Extract individual raw parameters (for compatibility)
    raw_logK = (mv_params[1] - mv_prior_mean[1]) / mv_prior_sd[1];
    raw_logr = (mv_params[2] - mv_prior_mean[2]) / mv_prior_sd[2];
    raw_logshape = (mv_params[3] - mv_prior_mean[3]) / mv_prior_sd[3];
    raw_logqeff = (mv_params[4] - mv_prior_mean[4]) / mv_prior_sd[4];
    
    // Transform AR correlation parameter
    rho = tanh(raw_rho * prior_rho_sd + prior_rho_mean); // keeps rho in [-1,1]
    
    // Effort-based fishing mortality calculation with dual error structure
    real sigma_qdev; 
    vector[n_periods] qdev_period; // transformed period deviations
    real qdev[Tm1]; // systematic catchability changes
    real edev[Tm1]; // effort measurement errors
    real F[Tm1];
    
    // Transform error variabilities
    sigma_qdev = raw_sigma_qdev * prior_sigma_qdev_sd;
    sigma_edev = raw_sigma_edev * prior_sigma_edev_sd;
    
    // Non-centered parameterization for AR(1) catchability process
    qdev_period[1] = raw_qdev_period[1] * sigma_qdev;
    for(p in 2:n_periods) {
        qdev_period[p] = rho * qdev_period[p-1] + raw_qdev_period[p] * sigma_qdev * sqrt(1 - rho^2);
    }
    
    // Assign period-based qdev with temporal structure
    for(t in 1:Tm1) {
        int period = ((t-1) / n_step) + 1; // integer division
        period = min(period, n_periods); // safer bound checking
        qdev[t] = qdev_period[period];
    }
    
    // Independent annual effort deviations (non-centered)
    for(t in 1:Tm1) {
        edev[t] = raw_edev[t] * sigma_edev;
    }
    
    // Calculate fishing mortality: F_t = qeff * exp(qdev_t) * effort_t * exp(edev_t - sigma_edev^2/2)
    real sigma_edev2;
    sigma_edev2 = pow(sigma_edev, 2);
    for(t in 1:Tm1) {
        F[t] = qeff * exp(qdev[t]) * effort[t] * exp(edev[t] - sigma_edev2/2);
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
    // Updated multivariate normal prior (4-dimensional)
    raw_mv_params ~ std_normal(); // Standard normal for raw parameters

    // Effort-based priors (non-centered)
    raw_sigma_qdev ~ std_normal();
    raw_rho ~ std_normal(); // Prior on raw rho parameter
    raw_qdev_period ~ std_normal(); // All catchability deviations are standard normal
    raw_edev ~ std_normal(); // All effort deviations are standard normal
    
    // prior densities for other estimated parameters
    raw_epsp ~ std_normal();
    raw_logsigmap ~ std_normal();
    raw_sigmao_add ~ std_normal();
    
    // Lognormal prior depletion prior from Gedamke and Hoenig 2006 size based Z estimate
    target += lognormal_lpdf(x[t_dep] | prior_depletion_meanlog, prior_depletion_sdlog);
    
    // observation model - uses analytical q[i]
    real mu_index;
    for(i in 1:I){
        if(lambdas[i]==1){
            for(t in 1:T){
                if(index[t,i]>0.0 && x[t]>0.0 && q[i]>0.0) {
                    mu_index = log(q[i]*x[t]) - sigmao2[t,i]/2;
                    target += lognormal_lpdf(index[t,i] | mu_index,sigmao[t,i]);
                }
            }
        }
    }

    // Catch observation model uses F derived from effort-based qeff
    real mu_catch;
    for(t in 1:Tm1){
        mu_catch = log(removals[t]) - 0.5*sigmac^2;
        target += lognormal_lpdf(obs_removals[t]|mu_catch,sigmac);
    }
}
