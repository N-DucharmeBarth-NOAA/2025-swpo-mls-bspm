# ISC SHARKWG
# 2024/06/11
# Robust prior pushforward for stochastic surplus production model
# Handles multiple Stan model structures including effort-based models

# Copyright (c) 2024 ISC SHARKWG
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_prior_pushforward = function(ssp_summary, stan_data, settings) {
    set.seed(settings$seed)
    chains = settings$chains
    n_samples = settings$iter_keep * chains
    
    # Helper function to safely extract data from stan_data
    get_stan_param = function(param_name, param_type = "Data", default_val = NULL) {
        result = stan_data[name == param_name & type == param_type]
        if (nrow(result) == 0) return(default_val)
        if (nrow(result) == 1) return(result$value)
        return(result[order(row)]$value)  # For vectors/matrices
    }
    
    # Helper function to check if parameter exists
    param_exists = function(param_name, param_type = "Data") {
        nrow(stan_data[name == param_name & type == param_type]) > 0
    }
    
    # Reconstruct correlation matrix helper
    reconstruct_corr_matrix = function(corr_name) {
        corr_data = stan_data[name == corr_name & type == "Data"]
        if (nrow(corr_data) == 0) return(NULL)
        
        max_row = max(corr_data$row)
        max_col = max(corr_data$col)
        corr_mat = matrix(NA, max_row, max_col)
        
        for(i in 1:nrow(corr_data)) {
            corr_mat[corr_data$row[i], corr_data$col[i]] = corr_data$value[i]
        }
        return(corr_mat)
    }
    
    # Check for different model configurations
    has_mv_prior = param_exists("mv_prior_mean")
    has_mv_qdev_prior = param_exists("mv_qdev_prior_mean")
    has_effort_params = param_exists("effort") || param_exists("n_periods")
    has_qeff = param_exists("prior_qeff_meanlog")
    
    # Initialize parameter storage
    logK = log_r = log_shape = r = shape = NULL
    raw_logK = raw_logr = raw_logshape = NULL
    logqeff = qeff = rho = sigma_qdev = NULL
    raw_logqeff = raw_rho = raw_sigma_qdev = NULL
    
    #--------------------------------------------------------------------
    # 1. Handle main multivariate prior (logK, log_r, log_shape)
    #--------------------------------------------------------------------
    if (has_mv_prior) {
        mv_mean = get_stan_param("mv_prior_mean")
        mv_sd = get_stan_param("mv_prior_sd")
        mv_corr_mat = reconstruct_corr_matrix("mv_prior_corr")
        
        if (is.null(mv_corr_mat)) {
            stop("Could not reconstruct mv_prior_corr matrix")
        }
        
        # Create covariance matrix
        mv_cov = diag(mv_sd) %*% mv_corr_mat %*% diag(mv_sd)
        
        # Sample from multivariate normal
        if (!requireNamespace("MASS", quietly = TRUE)) {
            stop("MASS package required for multivariate normal sampling")
        }
        mv_samples = MASS::mvrnorm(n_samples, mv_mean, mv_cov)
        
        # Extract parameters based on dimension
        logK = mv_samples[, 1]
        log_r = mv_samples[, 2]
        r = exp(log_r)
        
        # Calculate raw parameters
        raw_logK = (logK - mv_mean[1]) / mv_sd[1]
        raw_logr = (log_r - mv_mean[2]) / mv_sd[2]
        
        # Handle shape parameter (3rd dimension or separate)
        if (length(mv_mean) >= 3) {
            log_shape = mv_samples[, 3]
            shape = exp(log_shape)
            raw_logshape = (log_shape - mv_mean[3]) / mv_sd[3]
        } else {
            # Fallback for separate shape parameter
            if (param_exists("logshape", "PriorMean")) {
                shape_mean = get_stan_param("logshape", "PriorMean")
                shape_sd = get_stan_param("logshape", "PriorSD")
                log_shape = rnorm(n_samples, shape_mean, shape_sd)
                shape = exp(log_shape)
                raw_logshape = (log_shape - shape_mean) / shape_sd
            } else {
                # Default shape
                log_shape = rnorm(n_samples, 0, 0.5)
                shape = exp(log_shape)
                raw_logshape = log_shape / 0.5
            }
        }
        
    } else {
        # Fallback to univariate priors
        if (!is.null(ssp_summary$PriorMean_logK)) {
            logK = rnorm(n_samples, ssp_summary$PriorMean_logK, ssp_summary$PriorSD_logK)
            raw_logK = (logK - ssp_summary$PriorMean_logK) / ssp_summary$PriorSD_logK
        } else {
            stop("Cannot find logK prior parameters")
        }
        
        if (!is.null(ssp_summary$PriorMean_logr)) {
            log_r = rnorm(n_samples, ssp_summary$PriorMean_logr, ssp_summary$PriorSD_logr)
            raw_logr = (log_r - ssp_summary$PriorMean_logr) / ssp_summary$PriorSD_logr
            r = exp(log_r)
        } else {
            stop("Cannot find log_r prior parameters")
        }
        
        # Shape parameter
        if (param_exists("logshape", "PriorMean")) {
            shape_mean = get_stan_param("logshape", "PriorMean")
            shape_sd = get_stan_param("logshape", "PriorSD")
            log_shape = rnorm(n_samples, shape_mean, shape_sd)
            shape = exp(log_shape)
            raw_logshape = (log_shape - shape_mean) / shape_sd
        } else {
            log_shape = rnorm(n_samples, 0, 0.5)
            shape = exp(log_shape)
            raw_logshape = log_shape / 0.5
        }
    }
    
    #--------------------------------------------------------------------
    # 2. Handle qeff parameter (effort-based catchability)
    #--------------------------------------------------------------------
    if (has_qeff) {
        qeff_mean_log = get_stan_param("prior_qeff_meanlog")
        qeff_sd_log = get_stan_param("prior_qeff_sdlog")
        
        raw_logqeff = rnorm(n_samples, 0, 1)  # Standard normal
        logqeff = raw_logqeff * qeff_sd_log + qeff_mean_log
        qeff = exp(logqeff)
    }
    
    #--------------------------------------------------------------------
    # 3. Handle bivariate qdev prior (rho, sigma_qdev)
    #--------------------------------------------------------------------
    if (has_mv_qdev_prior) {
        qdev_mean = get_stan_param("mv_qdev_prior_mean")
        qdev_sd = get_stan_param("mv_qdev_prior_sd")
        qdev_corr_mat = reconstruct_corr_matrix("mv_qdev_prior_corr")
        
        if (!is.null(qdev_corr_mat)) {
            # Create covariance matrix
            qdev_cov = diag(qdev_sd) %*% qdev_corr_mat %*% diag(qdev_sd)
            
            # Sample from bivariate normal
            qdev_samples = MASS::mvrnorm(n_samples, qdev_mean, qdev_cov)
            
            # Transform back to natural scale
            rho = tanh(qdev_samples[, 1])  # atanh(rho) -> rho
            sigma_qdev = exp(qdev_samples[, 2])  # log(sigma_qdev) -> sigma_qdev
            
            # Calculate raw parameters
            raw_rho = (qdev_samples[, 1] - qdev_mean[1]) / qdev_sd[1]
            raw_sigma_qdev = (qdev_samples[, 2] - qdev_mean[2]) / qdev_sd[2]
        }
    }
    
    #--------------------------------------------------------------------
    # 4. Handle effort-based parameters
    #--------------------------------------------------------------------
    n_periods = get_stan_param("n_periods", default_val = NULL)
    n_step = get_stan_param("n_step", default_val = 3)
    sigma_edev = get_stan_param("sigma_edev", default_val = 0.2)
    effort = get_stan_param("effort", default_val = NULL)
    
    raw_qdev_period = raw_edev = qdev_period = edev = qdev = NULL
    
    if (has_effort_params && !is.null(n_periods)) {
        # Generate period-specific catchability deviations
        raw_qdev_period = matrix(rnorm(n_samples * n_periods, 0, 1), 
                                nrow = n_samples, ncol = n_periods)
        
        # Generate effort deviations
        T = get_stan_param("T", default_val = 30)
        Tm1 = T - 1
        raw_edev = matrix(rnorm(n_samples * Tm1, 0, 1), 
                         nrow = n_samples, ncol = Tm1)
        
        # Convert to actual deviations
        edev = raw_edev * sigma_edev
        
        # Calculate period-based catchability with AR(1) structure
        if (!is.null(rho) && !is.null(sigma_qdev)) {
            qdev_period = array(NA, dim = c(n_samples, n_periods))
            qdev = array(NA, dim = c(n_samples, Tm1))
            
            for (i in 1:n_samples) {
                # AR(1) process for periods
                qdev_period[i, 1] = raw_qdev_period[i, 1] * sigma_qdev[i]
                if (n_periods > 1) {
                    for (p in 2:n_periods) {
                        qdev_period[i, p] = rho[i] * qdev_period[i, p-1] + 
                                           raw_qdev_period[i, p] * sigma_qdev[i] * sqrt(1 - rho[i]^2)
                    }
                }
                
                # Assign period-based qdev to annual time series
                for (t in 1:Tm1) {
                    period = min(floor((t-1) / n_step) + 1, n_periods)
                    qdev[i, t] = qdev_period[i, period]
                }
            }
        }
    }
    
    #--------------------------------------------------------------------
    # 5. Handle standard parameters (process error, observation error)
    #--------------------------------------------------------------------
    
    # Process error
    sigmap_mean = get_stan_param("logsigmap", "PriorMean") %||%
                  get_stan_param("PriorMean_logsigmap") %||% -2.93
    sigmap_sd = get_stan_param("logsigmap", "PriorSD") %||%
                get_stan_param("PriorSD_logsigmap") %||% 0.27
    
    sigmap = exp(rnorm(n_samples, sigmap_mean, sigmap_sd))
    raw_logsigmap = (log(sigmap) - sigmap_mean) / sigmap_sd
    sigmap2 = sigmap^2
    
    # Observation error
    sigmao_input = get_stan_param("sigmao_input", default_val = 0.2)
    sigmao_add_sd = get_stan_param("sigmao_add", "PriorSD") %||%
                   get_stan_param("PriorSD_sigmao_add", default_val = 0.2)
    
    raw_sigmao_add = abs(rnorm(n_samples))
    sigmao_add = raw_sigmao_add * sigmao_add_sd
    sigmao_sc = sigmao_input + sigmao_add
    
    # Fishing mortality error (for F-based models)
    sigmaf_sd = get_stan_param("sigmaf", default_val = 0.3) %||%
               get_stan_param("PriorSD_sigmaf", default_val = 0.3)
    
    raw_sigmaf = abs(rnorm(n_samples))
    sigmaf = raw_sigmaf * sigmaf_sd
    
    #--------------------------------------------------------------------
    # 6. Generate time series and population dynamics
    #--------------------------------------------------------------------
    
    T = get_stan_param("T", default_val = 30)
    Tm1 = T - 1
    
    # Process error
    raw_epsp = epsp = matrix(NA, nrow = n_samples, ncol = T)
    for(j in 1:T) {
        epsp[, j] = rlnorm(n_samples, log(1.0) - (sigmap^2) / 2, sigmap)
    }
    raw_epsp = (log(epsp) - (sigmap^2) / 2) / sigmap
    
    # Fishing mortality
    F_matrix = matrix(NA, nrow = n_samples, ncol = Tm1)
    
    if (has_effort_params && !is.null(effort) && !is.null(qeff)) {
        # Effort-based F calculation
        for (i in 1:n_samples) {
            for (t in 1:Tm1) {
                qdev_val = if (!is.null(qdev)) qdev[i, t] else 0
                edev_val = if (!is.null(edev)) edev[i, t] else 0
                F_matrix[i, t] = qeff[i] * exp(qdev_val) * effort[t] * 
                               exp(edev_val - sigma_edev^2 / 2)
            }
        }
    } else {
        # Standard F estimation
        raw_F = matrix(abs(rnorm(n_samples * Tm1)), nrow = n_samples, ncol = Tm1)
        F_matrix = raw_F * sigmaf
    }
    
    # Derived parameters
    n = shape
    dmsy = (1 / n)^(1 / (n - 1))
    h = 2 * dmsy
    m = r * h / 4
    g = (n^(n / (n - 1))) / (n - 1)
    
    # Population dynamics
    x = matrix(NA, nrow = n_samples, ncol = T)
    x[, 1] = epsp[, 1]
    removals = matrix(NA, nrow = n_samples, ncol = Tm1)
    
    for(i in 1:n_samples) {
        for(t in 2:T) {
            if(x[i, t-1] <= dmsy[i]) {
                x[i, t] = (x[i, t-1] + r[i] * x[i, t-1] * (1 - x[i, t-1] / h[i])) * 
                         (exp(-F_matrix[i, t-1])) * epsp[i, t]
                removals[i, t-1] = ((x[i, t-1] + r[i] * x[i, t-1] * (1 - x[i, t-1] / h[i]))) * 
                                  epsp[i, t] * (1 - exp(-F_matrix[i, t-1])) * exp(logK[i])
            } else {
                x[i, t] = (x[i, t-1] + g[i] * m[i] * x[i, t-1] * (1 - x[i, t-1]^(n[i] - 1))) * 
                         (exp(-F_matrix[i, t-1])) * epsp[i, t]
                removals[i, t-1] = ((x[i, t-1] + g[i] * m[i] * x[i, t-1] * (1 - x[i, t-1]^(n[i] - 1)))) * 
                                  epsp[i, t] * (1 - exp(-F_matrix[i, t-1])) * exp(logK[i])
            }
        }
    }
    
    #--------------------------------------------------------------------
    # 7. Prepare output data tables
    #--------------------------------------------------------------------
    
    # Vector variables (only include those that exist)
    potential_vector_vars = c("raw_logK", "raw_logr", "raw_logsigmap", "raw_logshape", 
                             "logK", "r", "m", "sigmap", "sigmap2", "sigmao_sc", 
                             "shape", "n", "dmsy", "h", "g", "sigmao_add", "sigmaf",
                             "raw_logqeff", "logqeff", "qeff", "raw_rho", "raw_sigma_qdev", 
                             "rho", "sigma_qdev")
    
    vector_vars = potential_vector_vars[sapply(potential_vector_vars, function(v) {
        exists(v) && !is.null(get(v)) && is.vector(get(v))
    })]
    
    vector_dt.list = list()
    for(i in 1:length(vector_vars)) {
        var_val = get(vector_vars[i])
        vector_dt.list[[i]] = data.table(value = var_val) %>%
            .[, run_id := ssp_summary$run_id] %>%
            .[, iter := 1:n_samples] %>%
            .[, chain := sort(rep(1:chains, n_samples / chains))] %>%
            .[, chain_iter := rep(1:(n_samples / chains), chains)] %>%
            .[, variable := vector_vars[i]] %>%
            .[, name := vector_vars[i]] %>%
            .[, row := as.numeric(NA)] %>%
            .[, col := as.numeric(NA)] %>%
            .[, .(run_id, iter, chain, chain_iter, variable, name, row, col, value)] %>%
            .[order(chain, iter)]
    }
    vector_dt = rbindlist(vector_dt.list)
    
    # Matrix variables
    potential_matrix_vars = c("raw_epsp", "x", "removals", "F_matrix", "raw_qdev_period", 
                             "qdev_period", "raw_edev", "edev", "qdev")
    
    matrix_vars = potential_matrix_vars[sapply(potential_matrix_vars, function(v) {
        exists(v) && !is.null(get(v)) && is.matrix(get(v))
    })]
    
    matrix_dt.list = list()
    for(i in 1:length(matrix_vars)) {
        var_name = matrix_vars[i]
        # Handle F_matrix special case for naming consistency
        output_name = ifelse(var_name == "F_matrix", "F", var_name)
        
        matrix_dt.list[[i]] = as.data.table(get(var_name)) %>%
            .[, iter := 1:n_samples] %>%
            .[, chain := sort(rep(1:chains, n_samples / chains))] %>%
            .[, chain_iter := rep(1:(n_samples / chains), chains)] %>%
            melt(., id.vars = c("iter", "chain", "chain_iter")) %>%
            .[, row := as.numeric(gsub("V", "", variable))] %>%
            .[, run_id := ssp_summary$run_id] %>%
            .[, variable := paste0(output_name, "[", row, "]")] %>%
            .[, name := output_name] %>%
            .[, col := as.numeric(NA)] %>%
            .[, .(run_id, iter, chain, chain_iter, variable, name, row, col, value)] %>%
            .[order(chain, iter)]
    }
    matrix_dt = rbindlist(matrix_dt.list)
    out = rbind(vector_dt, matrix_dt)
    
    return(out)
}
