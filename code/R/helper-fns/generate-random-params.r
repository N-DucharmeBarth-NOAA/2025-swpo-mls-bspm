# Nicholas Ducharme-Barth
# 2025/05/31
# Generate random samples for BSPM parameters
# Based on biological parameter ranges and uncertainty
# Modified to include steepness parameter for sim_rmax function

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# function to generate random parameter samples
generate_random_params = function(n_samples = 1000, seed = NULL, min_weight_at_300cm = 150, max_weight_at_300cm = 350) {
    
    require(data.table)
    require(MASS) # for mvrnorm
    
    if(!is.null(seed)) {
        set.seed(seed)
    }
    
    # Function to check if weight at 300cm is within acceptable range
    check_weight_at_300cm = function(weight_a, weight_b, min_weight = 150, max_weight = 350) {
        weight_at_300cm = weight_a * (300 ^ weight_b)
        return(weight_at_300cm >= min_weight & weight_at_300cm <= max_weight)
    }
    
    # Generate more samples than needed to account for filtering
    n_generate = ceiling(n_samples * 1.5)  # Generate 50% extra to account for filtering
    
    # Generate correlated weight parameters using bivariate normal distribution
    # Based on fisheries literature, weight_a and weight_b are typically negatively correlated
    weight_a_mean = log(5.39942e-07)
    weight_b_mean = 3.58378
    weight_a_cv = 0.05  # CV for weight_a
    weight_b_cv = 0.05 # CV for weight_b (keeping this parameter tight)
    weight_correlation = -0.5  # negative correlation typical in L-W relationships
    
    # Convert CV to standard deviation for log-normal (weight_a) and normal (weight_b)
    weight_a_sd = weight_a_cv  # for log-normal, CV approximates sdlog when CV is small
    weight_b_sd = weight_b_mean * weight_b_cv
    
    # Covariance matrix for weight parameters
    weight_cov_matrix = matrix(c(weight_a_sd^2, 
                                weight_correlation * weight_a_sd * weight_b_sd,
                                weight_correlation * weight_a_sd * weight_b_sd, 
                                weight_b_sd^2), 
                              nrow = 2)
    
    # Generate correlated max_age and M_ref parameters
    # Negative correlation: higher M_ref associated with shorter max_age
    # M_ref will be log-normal, max_age will be normal
    max_age_target_mean = 15
    max_age_target_cv = 0.2
    max_age_meanlog = log(max_age_target_mean) - 0.5 * log(1 + max_age_target_cv^2)
    max_age_sdlog = sqrt(log(1 + max_age_target_cv^2))
    
    # Log-normal parameters for M_ref 
    M_ref_target_mean = 0.36 # 5.4/Amax Owen and Hamel
    M_ref_target_cv = 0.44  # Owen and Hamel
    M_ref_meanlog = log(M_ref_target_mean) - 0.5 * log(1 + M_ref_target_cv^2)
    M_ref_sdlog = sqrt(log(1 + M_ref_target_cv^2))
    
    age_M_correlation = -0.3  # loose negative correlation
        
    # Covariance matrix for log(max_age) and log(M_ref)
    age_logM_cov_matrix = matrix(c(max_age_sdlog^2, 
                                  age_M_correlation * max_age_sdlog * M_ref_sdlog,
                                  age_M_correlation * max_age_sdlog * M_ref_sdlog, 
                                  M_ref_sdlog^2), 
                                nrow = 2)
    
    # Keep generating until we have enough valid samples
    param_samples = data.table()
    attempts = 0
    max_attempts = 10
    
    while(nrow(param_samples) < n_samples && attempts < max_attempts) {
        attempts = attempts + 1
        
        # Generate correlated samples
        weight_params = mvrnorm(n = n_generate, 
                               mu = c(weight_a_mean, weight_b_mean), 
                               Sigma = weight_cov_matrix)
        
        # Generate correlated max_age and log(M_ref)
        age_logM_params = mvrnorm(n = n_generate,
                                 mu = c(max_age_meanlog, M_ref_meanlog),
                                 Sigma = age_logM_cov_matrix)
        
        # Initialize data.table to store samples
        temp_samples = data.table(
            # ID parameter for sim_rmax (renamed from sample_id for consistency)
            id = 1:n_generate,

            # Carrying capacity logK
            logK = log(rlnorm(n_generate,log(10e5),0.5)),
            
            # Correlated life history parameters (max_age and M_ref)
            max_age = round(exp(age_logM_params[,1])),  # Exponentiate and round
            M_ref = exp(age_logM_params[,2]),   # Exponentiate for log-normal
            
            # von Bertalanffy parameters (Francis parameterization - low correlation, treat as independent)
            L1 = rlnorm(n_generate, meanlog = log(60), sdlog = 0.2), # CV = 0.2, log-normal for positive values
            L2 = rlnorm(n_generate, meanlog = log(210), sdlog = 0.2), # CV = 0.2, log-normal for positive values
            vbk = rbeta(n_generate,shape1=6.5,shape2=3.5), # Beta distribution with mean ≈ 0.65 and reasonable variance
            
            cv_len = runif(n_generate, min = 0.05, max = 0.25),       # 0.1 to 0.25
            maturity_a = rnorm(n_generate, mean = -20, sd = abs(-20) * 0.2), # CV = 0.2 (keep normal for negative values)
            l50 = rlnorm(n_generate, meanlog = log(0.862069*214), sdlog = 0.2), # CV = 0.2, log-normal for positive values
            
            # Correlated weight parameters
            weight_a = exp(weight_params[,1]),                      # correlated log-normal
            weight_b = weight_params[,2],                           # correlated normal

            # sex-ratio
            sex_ratio = rnorm(n_generate, mean = 0.5, sd = 0.5 * 0.05), # CV = 0.05, small variability around 0.5
            
            # Steepness parameter (h) - required for sim_rmax
            # Using Beta distribution to constrain between 0.2 and 1.0, with mean around 0.7-0.8
            h = 0.2 + 0.8 * rbeta(n_generate, shape1 = 3, shape2 = 1.5), # Mean ≈ 0.73, constrained to [0.2, 1.0]

            # Length at first-capture; needed for F calc and depletion prior
            Lc = rlnorm(n_generate, meanlog = log(140), sdlog = 0.125),
            
            # Fixed parameters
            age1 = 0,
            age2 = 10,
            reproductive_cycle = 1
        )
        
        # Ensure L1 < L2 and biological constraints
        temp_samples[L1 >= L2, L1 := L2 - 10]  # Ensure L1 < L2
        temp_samples[L2 <= l50, L2 := l50 + 20]  # Ensure L2 > l50
        temp_samples[L2 <= 190, L2 := 190 + runif(.N, 5, 15)]  # Ensure L2 > 190
        
        # Ensure sex_ratio is bounded between 0 and 1
        temp_samples[sex_ratio < 0, sex_ratio := 0.01]
        temp_samples[sex_ratio > 1, sex_ratio := 0.99]
        
        # Apply weight check at 300cm
        temp_samples[, valid_weight := mapply(check_weight_at_300cm, weight_a, weight_b, 
                                            MoreArgs = list(min_weight = min_weight_at_300cm, max_weight = max_weight_at_300cm))]
        
        # Keep only valid samples
        valid_samples = temp_samples[valid_weight == TRUE]
        valid_samples[, valid_weight := NULL]  # Remove the check column
        
        # Append to our collection
        param_samples = rbind(param_samples, valid_samples)
        
        # Update sample IDs to be sequential
        param_samples[, id := 1:.N]
        
        cat("Attempt", attempts, ": Generated", nrow(valid_samples), "valid samples out of", n_generate, 
            "(", round(nrow(valid_samples)/n_generate*100, 1), "% passed weight check:", min_weight_at_300cm, "-", max_weight_at_300cm, "kg at 300cm)\n")
        cat("Total valid samples so far:", nrow(param_samples), "out of", n_samples, "needed\n")
    }
    
    # Take only the number of samples requested
    if(nrow(param_samples) >= n_samples) {
        param_samples = param_samples[1:n_samples]
        cat("Successfully generated", n_samples, "parameter sets that pass weight check (", min_weight_at_300cm, "-", max_weight_at_300cm, "kg at 300cm)\n")
    } else {
        cat("Warning: Only generated", nrow(param_samples), "valid samples out of", n_samples, "requested\n")
        cat("Consider relaxing the weight constraint or adjusting parameter distributions\n")
    }
    
    # Final check - report weight, steepness, and correlation statistics
    test_weights = param_samples[, weight_a * (300^weight_b)]
    cat("Weight at 300cm statistics:\n")
    cat("Range:", round(range(test_weights), 1), "kg\n")
    cat("Mean:", round(mean(test_weights), 1), "kg\n")
    cat("95% quantile:", round(quantile(test_weights, 0.95), 1), "kg\n")
    
    cat("Steepness (h) statistics:\n")
    cat("Range:", round(range(param_samples$h), 3), "\n")
    cat("Mean:", round(mean(param_samples$h), 3), "\n")
    cat("95% quantile:", round(quantile(param_samples$h, 0.95), 3), "\n")
    
    cat("Max_age and M_ref correlation:", round(cor(param_samples$max_age, param_samples$M_ref), 3), "\n")
    
    return(param_samples)
}
