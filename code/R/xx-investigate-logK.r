# Nicholas Ducharme-Barth
# 2025/07/02
# Investigate parameter correlations with logK from models where F is free or driven by effort
# Investigate correlations with early recdevs
# For free F models, calculate realized effort devs and or catchability trend

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(rstan)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# read in some samples
    model01 = fread(file.path("data","output","model_runs","0001-2024cpueExPrior_0","hmc_samples.csv"))
    model03 = fread(file.path("data","output","model_runs","0005-2024cpueMVPrior_0","hmc_samples.csv"))
    model25 = fread(file.path("data","output","model_runs","0025-dwfn-c0.2-e0.3-s5_0","hmc_samples.csv"))

    leading01 = model01[name %in% c("raw_logK","raw_logr","raw_logsigmap","raw_sigmao_add","raw_epsp","raw_logshape","raw_sigmaf","raw_F")]
    unique01 = unique(leading01$variable)
