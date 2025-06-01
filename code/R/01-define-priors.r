

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
    cv_age = 0.15
    A_max = 15 # Farley et al. 2021
    M_ref = 0.3 # 2024 assessment
    L1 = 60 # CAAL
    L2 = 210 # CAAL
    vbk = 0.8 # CAAL
    age1 = 0
    age2 = 10
    cv_len = 0.15
    maturity_a = -20
    l50 = (0.862069*214) # convert from LJFL to EOFL; L50 for females from Farley et al. 2021
    weight_a = 5.39942e-07 # 2024 assessment
    weight_b = 3.58378 # 2024 assessment
    sex_ratio = 0.5
    reproductive_cycle = 1

    sim_rmax(id, cv_age, A_max, M_ref, L1, L2, vbk, age1, age2, cv_len, maturity_a, l50, weight_a, weight_b,sex_ratio,reproductive_cycle)

