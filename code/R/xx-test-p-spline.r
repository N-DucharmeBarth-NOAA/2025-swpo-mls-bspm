# Nicholas Ducharme-Barth
# 2025/06/30
# Test splines for catchability

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(rstan)
    library(ggplot2)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
    dir_plot_fns = file.path(proj_dir,"code","R","plot-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    
    sapply(file.path(dir_plot_fns,(list.files(dir_plot_fns))),source)    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for model outputs
    dir.create(file.path(proj_dir,"data","output","model_runs"), showWarnings = FALSE, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load inputs
    catch_dt = fread(file.path(proj_dir,"data","input","catch.csv"))
    cpue_dt = fread(file.path(proj_dir,"data","input","cpue.csv"))
    
    # Load effort data 
    effort_dt = fread(file.path(proj_dir, "data", "input", "WCPFC_L_PUBLIC_BY_FLAG_YR.CSV")) %>%
                .[, lat_short_d := ifelse(grepl("S$", lat_short), 
                                        -as.numeric(gsub("[NS]", "", lat_short)), 
                                        as.numeric(gsub("[NS]", "", lat_short))) + 2.5] %>%
                .[, lon_short_d := ifelse(grepl("W$", lon_short), 
                                        360 - as.numeric(gsub("[EW]", "", lon_short)), 
                                            as.numeric(gsub("[EW]", "", lon_short))) + 2.5] %>%
                .[lat_short_d > -60 & lat_short_d < 0 & lon_short_d < 230] %>%
                .[!(lat_short_d > -5 & lon_short_d > 210)] %>%
                .[,.(effort_scaled = sum(hhooks)/1e6),by=yy] %>%
                .[yy %in% 1952:2022] %>%
                setnames(., "yy", "time")
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare catch data
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    # Merge with effort data
    catch_effort_annual = merge(catch_annual, effort_dt, by.x = "year", by.y = "time", all.x = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define stan data
    spline_options = expand.grid(num_knots=c(10,20,30,40,50,60),spline_degree=c(3))

    stan_data.list = as.list(rep(NA,nrow(spline_options)))

    for(i in seq_along(stan_data.list)){
        stan_data.list[[i]] = list(num_data = nrow(catch_effort_annual),
                    num_knots = spline_options$num_knots[i],
                    spline_degree = spline_options$spline_degree[i],
                    Y=catch_effort_annual[,.(nom_cpue=total_catch/effort_scaled)]$nom_cpue)
        stan_data.list[[i]]$X = 1:stan_data.list[[i]]$num_data
        stan_data.list[[i]]$knots = unname(quantile(stan_data.list[[i]]$X,probs=seq(from=0, to=1, length.out = stan_data.list[[i]]$num_knots)))
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile model
    exec_name = "p-spline"
    stan_c = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name,".stan")), model_name = exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# sample
    pred_dt.list = as.list(rep(NA,nrow(spline_options)))
    
    for(i in seq_along(pred_dt.list)){
        fit = sampling(stan_c, 
                    data = stan_data.list[[i]],
                    chains = 5,
                    warmup = 1000,
                    iter = 3000,
                    thin = 10,
                    seed = 123,
                    control = list(adapt_delta = 0.99,max_treedepth=12),
                    cores=5)
        quick_diagnostics(fit)
        
        pred_dt.list[[i]] = as.data.table(extract(fit,"Y_hat")) %>%
                        melt(.) %>%
                        .[,X:=as.numeric(as.factor(variable))] %>%
                        .[,.(med=quantile(value,probs=0.5),lower=quantile(value,probs=0.025),upper=quantile(value,probs=0.975)),by=X] %>%
                        .[,fit:=paste0("d",stan_data.list[[i]]$spline_degree,"-k",stan_data.list[[i]]$num_knots)]

        rm(list=c("fit"))
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# plot

    obs_dt = data.table(X=stan_data.list[[1]]$X,Y=stan_data.list[[1]]$Y)
    pred_dt = rbindlist(pred_dt.list) 

    p = obs_dt %>%
        ggplot() +
        geom_point(aes(x=X,y=log(Y)),shape=21,fill="black",color="white",size=3) +
        geom_path(data=pred_dt,aes(x=X,y=log(med),color=fit)) +
        # geom_hline(yintercept = 0) +
        get_ssp_theme()
    p

    p = merge(obs_dt,pred_dt,by="X") %>%
        .[,resid:=log(Y)-log(med)] %>%
        ggplot() +
        geom_path(aes(x=X,y=resid,color=fit)) +
        get_ssp_theme()
    p

    p = obs_dt %>%
        ggplot() +
        geom_point(aes(x=X,y=Y),shape=21,fill="black",color="white",size=3) +
        geom_path(data=pred_dt,aes(x=X,y=med,color=fit)) +
        # geom_hline(yintercept = 0) +
        get_ssp_theme()
    p
