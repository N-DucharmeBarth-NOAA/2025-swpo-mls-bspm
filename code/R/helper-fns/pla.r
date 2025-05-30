

# Nicholas Ducharme-Barth
# 2025/05/29
# Probability of length-at-age

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

pla_function = function(L, A, age_vec, len_lower, len_upper, L1, L2, vbk, age1, age2, cv_len) {
  # from Noel Cadigan
  pla <- matrix(0, nrow = L, ncol = A)
  
  for (a in 1:A) {
    ml <- L1 + (L2 - L1) * (1.0 - exp(-vbk * (age_vec[a] - age1))) / (1.0 - exp(-vbk * (age2 - age1)))
    sl <- cv_len * ml  # Convert coefficient of variation to standard deviation
    
    # Calculate probabilities
    for (z in 1:L) {
      len_lo_std <- (len_lower[z] - ml) / sl
      len_hi_std <- (len_upper[z] - ml) / sl
      pla[z, a] <- pnorm(len_hi_std) - pnorm(len_lo_std)
    }
    
    # Normalize the column
    pla[, a] <- pla[, a] / sum(pla[, a])
  }
  
  return(pla)
}
