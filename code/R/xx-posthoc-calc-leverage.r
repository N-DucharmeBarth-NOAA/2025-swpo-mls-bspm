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
    all_dirs = all_dirs[grep("fit_summary.csv", all_dirs, fixed = TRUE)]
    all_dirs = gsub("fit_summary.csv", "", all_dirs, fixed = TRUE)
    if (length(grep("-ppc", all_dirs, fixed = TRUE)) > 0) {
        all_dirs = all_dirs[-grep("-ppc", all_dirs, fixed = TRUE)]
    }
    all_dirs = all_dirs[grep("_0", all_dirs, fixed = TRUE)]
    all_dirs = gsub("/","",all_dirs)
    all_dirs = gsub("\\","",all_dirs,fixed=TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# post-hoc loo_dt calc
    n_cores = parallelly::availableCores(omit = 1,logical=FALSE)

    for(i in 1:length(all_dirs)){
        dir = file.path("data","output","model_runs",all_dirs[i])
        if(!file.exists(file.path(dir,"loo_dt.csv"))){
            hmc_samples = fread(file.path(dir,"hmc_samples.csv"))
            stan_data = fread(file.path(dir,"stan_data.csv"))

            # likelihood
                tmp_likelihood = ssp_calc_likelihood(hmc_samples,stan_data)

                log_lik_mat = tmp_likelihood %>% 
                          .[lambda>0] %>%
                          na.omit(.) %>%
                          .[,value:=value*lambda] %>%
                          .[,obs_id:=paste0(T,"-",I)] %>%
                          .[,.(obs_id,iter,value)] %>%
                          dcast(.,iter~obs_id) %>%
                          .[,iter:=NULL] %>%
                          as.matrix(.)
                lik_mat = tmp_likelihood %>% 
                          .[lambda>0] %>%
                          na.omit(.) %>%
                          .[,value:=exp(value*lambda)] %>%
                          .[,obs_id:=paste0(T,"-",I)] %>%
                          .[,.(obs_id,iter,value)] %>%
                          dcast(.,iter~obs_id) %>%
                          .[,iter:=NULL] %>%
                          as.matrix(.)
                
                r_eff = loo::relative_eff(lik_mat, unique(hmc_samples[,.(iter,chain)])$chain, cores = n_cores)
                loo = loo::loo(log_lik_mat,r_eff=r_eff,cores=n_cores)
                
                loo_influence = loo::pareto_k_influence_values(loo)
                loo_dt = ssp_calc_loofluence(loo_influence,stan_data)
                fwrite(loo_dt,file.path(dir,"loo_dt.csv"))
        }
    }
    


                p = loo_dt %>%
                    ggplot() + 
                    xlab("Time") +
                    ylab("Value") +
                    facet_wrap(~I,scales="free_y",nrow=length(unique(loo_dt$I))) +
                    geom_path(aes(x=T,y=value)) +
                    geom_point(aes(x=T,y=value,fill=infl_cat),shape=21,col="black",size=3) +
                    viridis::scale_color_viridis("LOO\nPareto-k:\nInfluence", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                    viridis::scale_fill_viridis("LOO\nPareto-k:\nInfluence", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE)
