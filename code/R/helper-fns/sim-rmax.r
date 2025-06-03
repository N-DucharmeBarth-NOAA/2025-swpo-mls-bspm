# Nicholas Ducharme-Barth
# 2025/05/29
# Calc Rmax from biological parameters - Modified approach
# openMSE: https://openmse.com/features-assessment-models/3-sp/
# Stanley, R.D., M. McAllister, P. Starr and N. Olsen. 2009. Stock assessment for bocaccio (Sebastes paucispinis) in British Columbia waters. DFO Can. Sci. Advis. Sec. Res. Doc. 2009/055. xiv + 200 p. 
# https://www.dfo-mpo.gc.ca/csas-sccs/publications/resdocs-docrech/2009/2009_055-eng.htm

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

sim_rmax = function(id, max_age, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b, sex_ratio, reproductive_cycle, h) {
  require(data.table)
  set.seed(id)
  
  age_vector = 1:max_age
  
  # survival (assume constant M)
  p = exp(-M_ref) # adult survival
  survival_at_age = rep(NA, max_age)
  survival_at_age[1] = 1
  for(i in 2:length(age_vector)) {
    survival_at_age[i] = survival_at_age[i-1] * exp(-M_ref)
  }
  survival_at_age[max_age] = survival_at_age[max_age-1] / (1 - exp(-M_ref))
  
  # length
  length_at_age = L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vector - age1))) / (1.0 - exp(-vbk * (age2 - age1))) 
  len_lower = seq(from = 0, to = ceiling(max(length_at_age * (1 + cv_len))), by = 1)
  len_upper = len_lower + 1
  length_vec = len_lower + 0.5
  pla_LA = pla_function(length(length_vec), length(age_vector), age_vector, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len)
  
  # maturity
  maturity_b = -maturity_a / l50
  maturity_at_length = (exp(maturity_a + maturity_b * length_vec)) / (1 + exp(maturity_a + maturity_b * length_vec))
  maturity_at_age = as.vector(matrix(maturity_at_length, nrow = 1, ncol = length(length_vec)) %*% pla_LA)
  maturity_at_age = maturity_at_age / max(maturity_at_age)
  
  # Calculate age at maturity
  if(length(which(maturity_at_age < 0.5)) > 0 & length(which(maturity_at_age > 0.5)) > 0) {
    a_mat = 0.5 * max(which(maturity_at_age < 0.5)) + 0.5 * min(which(maturity_at_age > 0.5))
    amat50 = round(a_mat)
  } else if(length(which(maturity_at_age > 0.5)) > 0) {
    amat50 = min(which(maturity_at_age > 0.5))
    a_mat = amat50
  } else {
    amat50 = NA
    a_mat = NA
  }
  
  # weight
  weight_at_age = weight_a * length_at_age^weight_b
  
  # fecundity
  fecundity_at_age = (maturity_at_age * weight_at_age * sex_ratio) / reproductive_cycle
  
  # calculate spawners per recruit in the unfished
  phio = sum(survival_at_age * fecundity_at_age)
  
  # slope at origin of BH-SRR
  alpha = (4 * h) / ((1 - h) * phio)
  
  # solve for Rmax using the new approach
  solve_rmax = function(max_age, alpha, survival_at_age, fecundity_at_age, start_value = 0.1) {
    # Rmax equation openMSE
    objective_function = function(r_max) {
      sum_euler_lotka = alpha * sum(survival_at_age * fecundity_at_age * exp(-r_max * (1:max_age)))
      return((sum_euler_lotka - 1)^2)
    }
    
    # Use nlminb to minimize the objective function
    fit = nlminb(start = start_value, 
                 objective = objective_function,
                 lower = -10,  # reasonable bounds for r_max
                 upper = 10)
    
    return(fit$par)
  }
  
  rmax = solve_rmax(max_age, alpha, survival_at_age, fecundity_at_age)
  
  # Calculate derived quantities
  spr = sum(fecundity_at_age * survival_at_age)
  generation_time = sum(age_vector * fecundity_at_age * survival_at_age) / spr # Grant and Grant 1992
  
  if(generation_time * rmax < 0) {
    inflection_point = NA
  } else {
    inflection_point = 0.633 - 0.187 * log(generation_time * rmax) # Fowler 1988
  }
  
  return(data.table(
    id = id,
    rmax = round(rmax, digits = 4),
    epr_unfished = round(spr, digits = 4),
    alpha = round(alpha, digits = 4),
    h = round(h, digits = 4),
    generation_time = round(generation_time, digits = 4),
    inflection_point = round(inflection_point, digits = 4),
    amat50 = amat50,
    max_age = max_age
  ))
}
