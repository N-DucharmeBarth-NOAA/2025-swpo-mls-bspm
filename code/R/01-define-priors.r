

# Nicholas Ducharme-Barth
# 2025/05/29
# Define priors for BSPM
# rmax, K, dep_1986

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#_____________________________________________________________________________________________________________________________
# ...    

    id = 1
    max_age = 15 # Farley et al. 2021 +/- 5yrs
    M_ref = 0.3 # 2024 assessment; 0.2 - 1
    L1 = 60 # CAAL; needs variability
    L2 = 210 # CAAL; needs variability
    vbk = 0.8 # CAAL; needs variability
    age1 = 0 # fixed
    age2 = 10 # fixed
    cv_len = 0.15 # 0.1 to 0.2
    maturity_a = -20; # needs variability small
    l50 = (0.862069*214) # convert from LJFL to EOFL; L50 for females from Farley et al. 2021; needs variability
    weight_a = 5.39942e-07 # 2024 assessment; needs variability
    weight_b = 3.58378 # 2024 assessment; needs variability
    sex_ratio = 0.5 # fixed at 0.5
    reproductive_cycle = 1 # fixed at 1
    selex_l50 = 150 # 150 - 180
    selex_slope = 15 # 0.1 - 1

    target_average_weight = 100
    selexNZ_l50 = 205
    selexNZ_slope = 0.3

    # sim Rmax
    sim_rmax(id, max_age, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b,sex_ratio,reproductive_cycle)

    calc_avg_capture_weight_at_age_unfished( max_age, L1, L2, vbk, age1, age2, cv_len,
                                   M_ref, selex_slope, selex_l50, 
                                   maturity_a, l50, weight_a, weight_b,selexNZ_l50,selexNZ_slope)

    # calc F in 1986 given biology and target average weight
    result = optim(par = M_ref,
                fn = objective_function,
                target_avg_weight = target_average_weight,
                # Pass all your other parameters
                max_age = max_age,
                L1 = L1, L2 = L2, vbk = vbk,
                age1 = age1, age2 = age2, cv_len = cv_len,
                M_ref = M_ref, selex_slope = selex_slope, selex_l50 = selex_l50,
                maturity_a = maturity_a, l50 = l50,
                weight_a = weight_a, weight_b = weight_b,selexNZ_l50 = selexNZ_l50,selexNZ_slope = selexNZ_slope,
                method = "Brent", 
                lower = 0,         
                upper = 2)         
    F = result$par

    calc_spr(id, max_age, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b,sex_ratio,reproductive_cycle,F,selex_l50,selex_slope)
