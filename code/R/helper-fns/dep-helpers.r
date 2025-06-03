#_____________________________________________________________________________________________________________________________
# define helper functions
transitional_mean_weight = function(W_inf, Wc, vbk, Z1, Z2, d) {
    # Equation 6 - modified to be in terms of weight
    # d time in years
    
    # Calculate each component separately for clarity
    term1 = W_inf
    
    numerator_fraction = Z1 * Z2 * (W_inf - Wc) * 
                        (Z1 + vbk + (Z2 - Z1) * exp(-(Z2 + vbk) * d))
    
    denominator_fraction = (Z1 + vbk) * (Z2 + vbk) * 
                            (Z1 + (Z2 - Z1) * exp(-Z2 * d))
    
    mean_weight = term1 - (numerator_fraction / denominator_fraction)
    
    return(mean_weight)
}

negative_log_likelihood = function(par, data, bio_params) {
    # Equation 8 implementation
    # define leading parameters
    Z1 = par[1]
    Z2 = par[2]

    # extract pars from bio_params
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    Lc = bio_params$Lc
    
    # define W_inf and Wc
    age_vector = seq(from=0, to=bio_params$max_age*2, length.out=50)
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    W_inf = weight_a*max(length_at_age)^weight_b
    Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

    # process data
    proc_dat = data[weight>Wc&year>=1952] %>%
               .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
               .[,d:=year-1952]

    # calculate d_vec
    d_vec = proc_dat$d

    # Calculate predicted lengths (using your transitional equation)
    predicted = transitional_mean_weight(W_inf, Wc, vbk, Z1, Z2, d_vec)
    
    log_likelihood = sum(dnorm(proc_dat$mean_wt, 
                         mean = predicted, 
                         sd = proc_dat$sd_wt/sqrt(proc_dat$N), 
                         log = TRUE))
    
    return(-log_likelihood)  # Return negative for minimization
}

rel_depletion = function(result, bio_params) {
    # Calculate relative depletion from Z1 to Z2
    # define parameters
    Z1 = result$par[1]
    Z2 = result$par[2]
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    cv_len = bio_params$cv_len
    maturity_a = bio_params$maturity_a
    l50 = bio_params$l50
    sex_ratio = bio_params$sex_ratio
    reproductive_cycle = bio_params$reproductive_cycle
    
    # age 
    age_vector = 1:max_age
    
    # length
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    len_lower = seq(from=0, to=ceiling(max(length_at_age*(1+cv_len))), by=1)
    len_upper = len_lower + 1
    length_vec = len_lower + 0.5

    # survival
    survival_at_age_Z1 = rep(NA, max_age)
    survival_at_age_Z1[1] = 1
    survival_at_age_Z2 = rep(NA, max_age)
    survival_at_age_Z2[1] = 1
    for(i in 2:length(age_vector)) {
        survival_at_age_Z1[i] = survival_at_age_Z1[i-1] * exp(-Z1)
        survival_at_age_Z2[i] = survival_at_age_Z2[i-1] * exp(-Z2)
    }
    survival_at_age_Z1[max_age] = survival_at_age_Z1[max_age-1] / (1 - exp(-Z1))
    survival_at_age_Z2[max_age] = survival_at_age_Z2[max_age-1] / (1 - exp(-Z2))

    # maturity
    maturity_b = -maturity_a/l50
    maturity_at_length = (exp(maturity_a+maturity_b*length_vec)) / (1+exp(maturity_a+maturity_b*length_vec))

    # calc maturity at age using PLA
    pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
    maturity_at_age = as.vector(matrix(maturity_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
    maturity_at_age = maturity_at_age/max(maturity_at_age)
    
    # weight
    weight_at_length = weight_a * length_vec ^ weight_b
    weight_at_age = as.vector(matrix(weight_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)

    # reproductive output per year (sex-ratio & reproductive cycle)
    reproduction_at_age = (maturity_at_age * weight_at_age * sex_ratio)/reproductive_cycle

    # spr calc
    # lifetime average eggs per recruit in fished and unfished conditions
    epr_Z1 = sum(reproduction_at_age*survival_at_age_Z1)
    epr_Z2 = sum(reproduction_at_age*survival_at_age_Z2)
    rel_dep_n = sum(survival_at_age_Z2)/sum(survival_at_age_Z1)
    rel_dep_ssb = epr_Z2/epr_Z1
    
    return(data.table(sample_id=bio_params$sample_id, rel_dep_n = rel_dep_n, rel_dep_ssb = rel_dep_ssb))
}

mean_wt_resid = function(result, bio_params, data) {
    # Calculate mean weight residuals
    # define parameters
    Z1 = result$par[1]
    Z2 = result$par[2]
    max_age = bio_params$max_age
    L1 = bio_params$L1
    L2 = bio_params$L2
    vbk = bio_params$vbk
    age1 = bio_params$age1
    age2 = bio_params$age2
    weight_a = bio_params$weight_a
    weight_b = bio_params$weight_b
    Lc = bio_params$Lc

    # define W_inf and Wc
    age_vector = seq(from=0, to=bio_params$max_age*2, length.out=50)
    length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    W_inf = weight_a*max(length_at_age)^weight_b
    Wc = weight_a*length_at_age[max(which(length_at_age<Lc))]^weight_b

    # process data
    proc_dat = data[weight>Wc&year>=1952] %>%
            .[,.(mean_wt=mean(weight),sd_wt=sd(weight),.N),by=year] %>%
            .[,d:=year-1952]

    # calculate d_vec
    d_vec = proc_dat$d

    # Calculate predicted lengths (using your transitional equation)
    predicted = transitional_mean_weight(W_inf, Wc, vbk, Z1, Z2, d_vec)
    
    proc_dat$pred_wt = predicted
    proc_dat$resid = proc_dat$mean_wt - proc_dat$pred_wt
    proc_dat$sample_id = bio_params$sample_id

    return(proc_dat[,.(sample_id, year, d, N, mean_wt, sd_wt, pred_wt, resid)])
}

#_____________________________________________________________________________________________________________________________
# Define processing function for a single parameter set
process_depletion_set = function(i, bio_params_dt, nz_wt) {
    
    # Extract parameters for current iteration
    params = bio_params_dt[i]
    
    # Initialize result list
    result = list(
        sample_id = params$sample_id,
        Z1 = NA_real_,
        Z2 = NA_real_,
        convergence = NA_integer_,
        log_likelihood = NA_real_,
        rel_dep_n = NA_real_,
        rel_dep_ssb = NA_real_,
        error_message = NA_character_
    )
    
    tryCatch({
        # Fit Z1 and Z2 using optimization
        opt_result = nlminb(start = c(params$M_ref + 0.1, params$M_ref + 0.2),
                           objective = negative_log_likelihood,
                           data = nz_wt,
                           bio_params = params,
                           lower = rep(params$M_ref, 2),
                           upper = rep(3, 2))
        
        # Store optimization results
        result$Z1 = opt_result$par[1]
        result$Z2 = opt_result$par[2]
        result$convergence = opt_result$convergence
        result$log_likelihood = opt_result$objective
        
        # Calculate relative depletion if optimization converged
        if(opt_result$convergence == 0) {
            depletion_result = rel_depletion(opt_result, params)
            result$rel_dep_n = depletion_result$rel_dep_n
            result$rel_dep_ssb = depletion_result$rel_dep_ssb
        }
        
    }, error = function(e) {
        result$error_message <<- as.character(e$message)
    })
    
    return(result)
}

# Define separate function to calculate residuals for successful runs
calculate_residuals_parallel = function(successful_dt, nz_wt) {
    
    cat("Calculating residuals for", nrow(successful_dt), "successful parameter sets...\n")
    
    # Setup parallel processing for residuals
    n_cores_resid = min(parallel::detectCores() - 1, nrow(successful_dt))
    cl_resid = makeCluster(n_cores_resid)
    
    # Export necessary objects to workers
    clusterEvalQ(cl_resid, {
        library(data.table)
        library(magrittr)
    })
    
    # Export functions and data to workers
    clusterExport(cl_resid, c("mean_wt_resid", "transitional_mean_weight", "nz_wt"))
    
    # Function to calculate residuals for one parameter set
    calc_residuals_single = function(i, successful_dt, nz_wt) {
        params = successful_dt[i]
        
        # Recreate opt_result structure from stored values
        opt_result = list(par = c(params$Z1, params$Z2))
        
        tryCatch({
            residuals = mean_wt_resid(opt_result, params, nz_wt)
            return(residuals)
        }, error = function(e) {
            return(data.table(sample_id = params$sample_id, error = as.character(e$message)))
        })
    }
    
    # Calculate residuals in parallel
    residuals_list = parLapply(cl_resid, 1:nrow(successful_dt), calc_residuals_single, 
                              successful_dt = successful_dt, nz_wt = nz_wt)
    
    # Stop the cluster
    stopCluster(cl_resid)
    
    # Combine residuals
    residuals_dt = rbindlist(residuals_list, fill = TRUE)
    
    return(residuals_dt)
}
