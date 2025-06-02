

# Nicholas Ducharme-Barth
# 2025/05/31
# Function to calculate average weight at age given F
# Assume knife-edge length-based selectivity

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Function to calculate average weight at age given F
calc_avg_capture_weight_at_age_fished = function(F_value, 
                                   max_age, L1, L2, vbk, age1, age2, cv_len,
                                   M_ref, selex_slope, selex_l50, 
                                   maturity_a, l50, weight_a, weight_b,
                                   selexNZ_l50,selexNZ_slope) {

  age_vector = 1:max_age
  
  # length
  length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
  len_lower = seq(from=0, to=ceiling(max(length_at_age*(1+cv_len))), by=1)
  len_upper = len_lower + 1
  length_vec = len_lower + 0.5
  
  # natural mortality
  mortality_at_age = M_ref * (max_age / age_vector)^(-1)
  
  # selectivity
  selex_length = 1 / (1 + exp(-selex_slope*(length_vec - selex_l50)))
  selexNZ_length = 1 / (1 + exp(-selexNZ_slope*(length_vec - selexNZ_l50)))
  
  # fishing mortality
  f_at_length = F_value * selex_length
  
  pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, 
                        len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
  f_at_age = as.vector(matrix(f_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
  selex_at_age = as.vector(matrix(selex_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
  selexNZ_at_age = as.vector(matrix(selexNZ_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)

  
  # survival
  survival_at_age_fished = rep(NA, max_age)
  survival_at_age_fished[1] = exp(-(mortality_at_age[1] + f_at_age[1]))
  for(i in 2:length(age_vector)){
    survival_at_age_fished[i] = survival_at_age_fished[i-1] * exp(-(mortality_at_age[i] + f_at_age[i]))
  }
  
  # weight
  weight_at_length = weight_a * length_vec ^ weight_b
  weight_at_age = as.vector(matrix(weight_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
  
  # Calculate numbers at age for fished population
  numbers_at_age_fished = c(1, survival_at_age_fished[-max_age])  # Assuming recruitment = 1
  
  # Calculate average weight weighted by numbers at age
  avg_weight_fished = sum(weight_at_age * numbers_at_age_fished * selexNZ_at_age) / sum(numbers_at_age_fished * selexNZ_at_age)
  
  return(avg_weight_fished)
}

# Objective function to minimize (squared difference from target)
objective_function = function(F_value, target_avg_weight, ...) {
  predicted_avg_weight = calc_avg_capture_weight_at_age_fished(F_value, ...)
  return((predicted_avg_weight - target_avg_weight)^2)
}


calc_avg_capture_weight_at_age_unfished = function( max_age, L1, L2, vbk, age1, age2, cv_len,
                                   M_ref, selex_slope, selex_l50, 
                                   maturity_a, l50, weight_a, weight_b,
                                   selexNZ_l50,selexNZ_slope) {

  age_vector = 1:max_age
  
  # length
  length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
  len_lower = seq(from=0, to=ceiling(max(length_at_age*(1+cv_len))), by=1)
  len_upper = len_lower + 1
  length_vec = len_lower + 0.5

  pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, 
                        len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
  
  # natural mortality
  mortality_at_age = M_ref * (max_age / age_vector)^(-1)
  
  # selectivity
  selex_length = 1 / (1 + exp(-selex_slope*(length_vec - selex_l50)))
  selex_at_age = as.vector(matrix(selex_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
  selexNZ_length = 1 / (1 + exp(-selexNZ_slope*(length_vec - selexNZ_l50)))
  selexNZ_at_age = as.vector(matrix(selexNZ_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)

  # survival
  survival_at_age = rep(NA, max_age)
  survival_at_age[1] = exp(-(mortality_at_age[1]))
  for(i in 2:length(age_vector)){
    survival_at_age[i] = survival_at_age[i-1] * exp(-(mortality_at_age[i]))
  }
  
  # weight
  weight_at_length = weight_a * length_vec ^ weight_b
  weight_at_age = as.vector(matrix(weight_at_length, nrow=1, ncol=length(length_vec)) %*% pla_LA)
  
  # Calculate numbers at age for fished population
  numbers_at_age = c(1, survival_at_age[-max_age])  # Assuming recruitment = 1
  
  # Calculate average weight weighted by numbers at age
  avg_weight_fished = sum(weight_at_age * numbers_at_age * selexNZ_at_age) / sum(numbers_at_age * selexNZ_at_age)
  
  return(avg_weight_fished)
}
