# Nicholas Ducharme-Barth
# 2025/07/15
# Use loo to compare models

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
    all_dirs = all_dirs[grep("fit_summary.csv", all_dirs, fixed = TRUE)]
    all_dirs = gsub("fit_summary.csv", "", all_dirs, fixed = TRUE)
    if (length(grep("-ppc", all_dirs, fixed = TRUE)) > 0) {
        all_dirs = all_dirs[-grep("-ppc", all_dirs, fixed = TRUE)]
    }
    all_dirs = all_dirs[grep("_0", all_dirs, fixed = TRUE)]
    all_dirs = gsub("/","",all_dirs)
    all_dirs = gsub("\\","",all_dirs,fixed=TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# loo compare
    model_directories = all_dirs[c(100,111,112,113,114,124,126)]
    results = ssp_compare_models_loo(model_directories,cache_loo=TRUE)

    # Access comparison results
    print(results$comparison)

