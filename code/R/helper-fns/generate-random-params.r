# Nicholas Ducharme-Barth
# 2025/05/31
# Generate random samples for BSPM parameters
# Based on biological parameter ranges and uncertainty

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
    correlation = -0.5  # negative correlation typical in L-W relationships
    
    # Convert CV to standard deviation for log-normal (weight_a) and normal (weight_b)
    weight_a_sd = weight_a_cv  # for log-normal, CV approximates sdlog when CV is small
    weight_b_sd = weight_b_mean * weight_b_cv
    
    # Covariance matrix
    cov_matrix = matrix(c(weight_a_sd^2, 
                         correlation * weight_a_sd * weight_b_sd,
                         correlation * weight_a_sd * weight_b_sd, 
                         weight_b_sd^2), 
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
                               Sigma = cov_matrix)
        
        # Initialize data.table to store samples
        temp_samples = data.table(
            sample_id = 1:n_generate,

            # Carrying capacity logK
            logK = log(rlnorm(n_generate,log(10e5),0.5)),
            
            # Variable parameters with specified ranges
            max_age = round(runif(n_generate, min = 10, max = 20)),  # 15 +/- 5yrs
            M_ref = runif(n_generate, min = 0.2, max = 1.0),        # 0.2 - 1
            
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
            
            selex_l50 = runif(n_generate, min = 150, max = 180),    # 150 - 180
            selex_slope = runif(n_generate, min = 0.1, max = 1),    # 0.1 - 1
            selexNZ_l50 = runif(n_generate, min = 190, max = 210),  # 190 - 210
            selexNZ_slope = runif(n_generate, min = 0.1, max = 1),  # 0.1 - 1

            # sex-ratio
            sex_ratio = rnorm(n_generate, mean = 0.5, sd = 0.5 * 0.05), # CV = 0.05, small variability around 0.5
            
            # Fixed parameters
            age1 = 0,
            age2 = 10,
            reproductive_cycle = 1
        )
        
        # Ensure L1 < L2 and biological constraints
        temp_samples[L1 >= L2, L1 := L2 - 10]  # Ensure L1 < L2
        temp_samples[L2 <= l50, L2 := l50 + 20]  # Ensure L2 > l50
        temp_samples[L2 <= 190, L2 := 190 + runif(.N, 5, 15)]  # Ensure L2 > 190
        
        # Apply weight check at 300cm
        temp_samples[, valid_weight := mapply(check_weight_at_300cm, weight_a, weight_b, 
                                            MoreArgs = list(min_weight = min_weight_at_300cm, max_weight = max_weight_at_300cm))]
        
        # Keep only valid samples
        valid_samples = temp_samples[valid_weight == TRUE]
        valid_samples[, valid_weight := NULL]  # Remove the check column
        
        # Append to our collection
        param_samples = rbind(param_samples, valid_samples)
        
        # Update sample IDs to be sequential
        param_samples[, sample_id := 1:.N]
        
        cat("Attempt", attempts, ": Generated", nrow(valid_samples), "valid samples out of", n_generate, 
            "(", round(nrow(valid_samples)/n_generate*100, 1), "% passed weight check: 125-500kg at 300cm)\n")
        cat("Total valid samples so far:", nrow(param_samples), "out of", n_samples, "needed\n")
    }
    
    # Take only the number of samples requested
    if(nrow(param_samples) >= n_samples) {
        param_samples = param_samples[1:n_samples]
        cat("Successfully generated", n_samples, "parameter sets that pass weight check (125-500kg at 300cm)\n")
    } else {
        cat("Warning: Only generated", nrow(param_samples), "valid samples out of", n_samples, "requested\n")
        cat("Consider relaxing the weight constraint or adjusting parameter distributions\n")
    }
    
    # Final check - report weight statistics
    test_weights = param_samples[, weight_a * (300^weight_b)]
    cat("Weight at 300cm statistics:\n")
    cat("Range:", round(range(test_weights), 1), "kg\n")
    cat("Mean:", round(mean(test_weights), 1), "kg\n")
    cat("95% quantile:", round(quantile(test_weights, 0.95), 1), "kg\n")
    
    return(param_samples)
}
