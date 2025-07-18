# Nicholas Ducharme-Barth
# 2025/07/02
# Post-hoc calculate leverage of individual observations on model outputs

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
# identify directories
    model_runs_path = file.path(proj_dir, "data", "output", "model_runs")
    all_dirs = list.files(model_runs_path, recursive = TRUE)
    all_dirs = all_dirs[grep("stan_data.csv", all_dirs, fixed = TRUE)]
    all_dirs = gsub("stan_data.csv", "", all_dirs, fixed = TRUE)
    all_dirs = gsub("/","",all_dirs)
    all_dirs = gsub("\\","",all_dirs,fixed=TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# update n_obs calc

for(i in 1:length(all_dirs)){
    stan_data = fread(file.path(model_runs_path,all_dirs[i],"stan_data.csv"))
    n_obs = ssp_count_total_observations(stan_data)
    fit_summary = fread(file.path(model_runs_path,all_dirs[i],"fit_summary.csv"))
    fit_summary$n_obs = n_obs
    fwrite(fit_summary,file=file.path(model_runs_path,all_dirs[i],"fit_summary.csv"))
}

