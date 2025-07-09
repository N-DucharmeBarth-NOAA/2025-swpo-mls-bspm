# Flexible initialization function for all BSPM model variants
stan_inits_func = function(Tm1, n_periods = NULL, exec_name = "auto") {
    # Base initialization that all models share
    base_inits = list(
        raw_epsp = rnorm(Tm1 + 1, 0, 0.25),
        raw_logsigmap = rnorm(1, 0, 0.25),
        raw_sigmao_add = abs(rnorm(1, 0, 0.25))
    )
    
    if(exec_name %in% c("basic", "bspm_estF", "bspm_estF_softdep")) {
        # Basic F estimation models
        specific_inits = list(
            raw_logK = rnorm(1, 0, 0.25),
            raw_logr = rnorm(1, 0, 0.25),
            raw_logshape = rnorm(1, 0, 0.25),
            raw_sigmaf = abs(rnorm(1, 0, 0.25)),
            raw_F = abs(rnorm(Tm1, 0, 0.25))
        )
    } else if(exec_name %in% c("mvprior", "bspm_estF_softdep_mvprior","bspm_estF_mvprior")) {
        # Multivariate prior models (3D)
        specific_inits = list(
            raw_mv_params = rnorm(3, 0, 0.25),  # 3D multivariate
            raw_sigmaf = abs(rnorm(1, 0, 0.25)),
            raw_F = abs(rnorm(Tm1, 0, 0.25))
        )
    } else if(exec_name %in% c("effort", "bspm_estq_softdep_mvprior", "bspm_estq_flex")) {
        # Effort-based models with full parameter set
        if(is.null(n_periods)) {
            n_step = 3  # default
            n_periods = ceiling(Tm1 / n_step)
        }
        
        specific_inits = list(
            raw_mv_params = rnorm(3, 0, 0.25),  # 3D multivariate for logK, log_r, log_shape
            raw_logqeff = rnorm(1, 0, 0.25),  # Independent qeff parameter
            raw_qdev_params = rnorm(2, 0, 0.25),  # Bivariate rho/sigma_qdev
            raw_qdev_period = rnorm(n_periods, 0, 0.25),  # Period-specific catchability deviations
            raw_edev = rnorm(Tm1, 0, 0.25)  # Annual effort deviations
        )
    } else if (exec_name %in% c("bspm_estqsimple_softdep_mvprior")){ 
        # Effort-based models with full parameter set
        if(is.null(n_periods)) {
            n_step = 3  # default
            n_periods = ceiling(Tm1 / n_step)
        }
        
        specific_inits = list(
            raw_mv_params = rnorm(3, 0, 0.25),  # 3D multivariate for logK, log_r, log_shape
            raw_logqeff = rnorm(1, 0, 0.25),  # Independent qeff parameter
            raw_qdev_params = rnorm(2, 0, 0.25),  # Bivariate rho/sigma_qdev
            raw_qdev_period = rnorm(n_periods-1, 0, 0.25)  # Period-specific catchability deviations
        )
    } else if (exec_name %in% c("bspm_estqsimple_softdep_mvprior_x0")){ 
        # Effort-based models with full parameter set
        if(is.null(n_periods)) {
            n_step = 3  # default
            n_periods = ceiling(Tm1 / n_step)
        }
        
        specific_inits = list(
            raw_mv_params = rnorm(5, 0, 0.25),  # 5D multivariate for logK, log_r, log_shape, log_x0, log_qeff
            raw_qdev_params = rnorm(2, 0, 0.25),  # Bivariate rho/sigma_qdev
            raw_qdev_period = rnorm(n_periods-1, 0, 0.25)  # Period-specific catchability deviations
        )
    } else if(exec_name %in% c("bspm_estq_softdep_mvprior_x0", "bspm_estq_flex_x0")) {
        # Effort-based models with full parameter set
        if(is.null(n_periods)) {
            n_step = 3  # default
            n_periods = ceiling(Tm1 / n_step)
        }
        
        specific_inits = list(
            raw_mv_params = rnorm(4, 0, 0.25),  # 3D multivariate for logK, log_r, log_shape, log_x0
            raw_logqeff = rnorm(1, 0, 0.25),  # Independent qeff parameter
            raw_qdev_params = rnorm(2, 0, 0.25),  # Bivariate rho/sigma_qdev
            raw_qdev_period = rnorm(n_periods, 0, 0.25),  # Period-specific catchability deviations
            raw_edev = rnorm(Tm1, 0, 0.25)  # Annual effort deviations
        )
    } else {
        stop("Unknown model_type: ", exec_name)
    }
    
    # Combine base and specific initializations
    inits = c(base_inits, specific_inits)
    
    return(inits)
}
