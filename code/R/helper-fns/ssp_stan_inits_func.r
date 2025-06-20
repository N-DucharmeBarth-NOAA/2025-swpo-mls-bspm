# Flexible initialization function for all BSPM model variants
stan_inits_func = function(Tm1, n_periods = NULL, model_type = "auto") {
    
    # Base initialization that all models share
    base_inits = list(
        raw_epsp = rnorm(Tm1 + 1, 0, 0.25),
        raw_logsigmap = rnorm(1, 0, 0.25),
        raw_sigmao_add = abs(rnorm(1, 0, 0.25))
    )
    
    # Model-specific parameters
    if(model_type == "auto") {
        # Try to detect model type based on n_periods argument
        if(is.null(n_periods)) {
            model_type = "basic"  # Assume basic F model
        } else {
            model_type = "effort"  # Effort-based model
        }
    }
    
    if(model_type %in% c("basic", "bspm_estF", "bspm_estF_softdep")) {
        # Basic F estimation models
        specific_inits = list(
            raw_logK = rnorm(1, 0, 0.25),
            raw_logr = rnorm(1, 0, 0.25),
            raw_logshape = rnorm(1, 0, 0.25),
            raw_sigmaf = abs(rnorm(1, 0, 0.25)),
            raw_F = abs(rnorm(Tm1, 0, 0.25))
        )
    } else if(model_type %in% c("mvprior", "bspm_estF_softdep_mvprior")) {
        # Multivariate prior models (3D)
        specific_inits = list(
            raw_mv_params = rnorm(3, 0, 0.25),  # 3D multivariate
            raw_sigmaf = abs(rnorm(1, 0, 0.25)),
            raw_F = abs(rnorm(Tm1, 0, 0.25))
        )
    } else if(model_type %in% c("effort", "bspm_estq_softdep_mvprior", "bspm_estq_optimized")) {
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
    } else {
        stop("Unknown model_type: ", model_type)
    }
    
    # Combine base and specific initializations
    inits = c(base_inits, specific_inits)
    
    return(inits)
}

# Wrapper function for backward compatibility
stan_inits_func_auto = function(Tm1, n_periods = NULL, exec_name = NULL) {
    
    # Auto-detect model type from exec_name if provided
    if(!is.null(exec_name)) {
        if(grepl("estq.*mvprior", exec_name) || grepl("effort", exec_name)) {
            model_type = "effort"
        } else if(grepl("mvprior", exec_name)) {
            model_type = "mvprior"
        } else if(grepl("estF", exec_name)) {
            model_type = "basic"
        } else {
            model_type = "auto"
        }
    } else {
        model_type = "auto"
    }
    
    return(stan_inits_func(Tm1, n_periods, model_type))
}
