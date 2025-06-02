# Nicholas Ducharme-Barth
# 2025/05/31
# Generate random samples for BSPM parameters
# Based on biological parameter ranges and uncertainty

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# function to generate random parameter samples
generate_random_params = function(n_samples = 1000, seed = NULL) {
    
    require(data.table)
    require(MASS) # for mvrnorm
    
    if(!is.null(seed)) {
        set.seed(seed)
    }
    
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
    
    # Generate correlated samples
    weight_params = mvrnorm(n = n_samples, 
                           mu = c(weight_a_mean, weight_b_mean), 
                           Sigma = cov_matrix)
    
    # Initialize data.table to store samples
    param_samples = data.table(
        sample_id = 1:n_samples,

        # Carrying capacity logK
        logK = log(rlnorm(n_samples,log(10e5),0.5)),
        
        # Variable parameters with specified ranges
        max_age = round(runif(n_samples, min = 10, max = 20)),  # 15 +/- 5yrs
        M_ref = runif(n_samples, min = 0.2, max = 1.0),        # 0.2 - 1
        
        # von Bertalanffy parameters (Francis parameterization - low correlation, treat as independent)
        L1 = rlnorm(n_samples, meanlog = log(60), sdlog = 0.2), # CV = 0.2, log-normal for positive values
        L2 = rlnorm(n_samples, meanlog = log(210), sdlog = 0.2), # CV = 0.2, log-normal for positive values
        vbk = rbeta(n_samples,shape1=6.5,shape2=3.5), # Beta distribution with mean ≈ 0.65 and reasonable variance
        
        cv_len = runif(n_samples, min = 0.05, max = 0.25),       # 0.1 to 0.25
        maturity_a = rnorm(n_samples, mean = -20, sd = abs(-20) * 0.2), # CV = 0.2 (keep normal for negative values)
        l50 = rlnorm(n_samples, meanlog = log(0.862069*214), sdlog = 0.2), # CV = 0.2, log-normal for positive values
        
        # Correlated weight parameters
        weight_a = exp(weight_params[,1]),                      # correlated log-normal
        weight_b = weight_params[,2],                           # correlated normal
        
        selex_l50 = runif(n_samples, min = 150, max = 180),    # 150 - 180
        selex_slope = runif(n_samples, min = 0.1, max = 1),    # 0.1 - 1
        selexNZ_l50 = runif(n_samples, min = 190, max = 210),  # 190 - 210
        selexNZ_slope = runif(n_samples, min = 0.1, max = 1),  # 0.1 - 1

        # sex-ratio
        sex_ratio = rnorm(n_samples, mean = 0.5, sd = 0.5 * 0.05), # CV = 0.05, small variability around 0.5
        
        # Fixed parameters
        age1 = 0,
        age2 = 10,
        reproductive_cycle = 1
    )
    
    # Ensure L1 < L2 (only constraint needed now since log-normal ensures positive values)
    param_samples[L1 >= L2, L1 := L2 - 10]  # Ensure L1 < L2
    
    return(param_samples)
}
