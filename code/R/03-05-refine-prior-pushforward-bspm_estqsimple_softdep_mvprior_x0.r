# Nicholas Ducharme-Barth
# 2025/07/07
# Conduct prior pushforward for bspm_estqsimple_softdep_mvprior_x0.stan

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(rstan)
    library(ggplot2)
    library(GGally)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_output = file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0_refined")
    dir.create(dir_output,showWarnings = FALSE, recursive = TRUE)
    dir_plot = file.path(proj_dir,"plots","pushforward","bspm_estqsimple_softdep_mvprior_x0_refined")
    dir.create(dir_plot,showWarnings = FALSE, recursive = TRUE)
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load data
    hmc_samples = fread(file=file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0","hmc_samples.csv"))
    hmc_samples[,value:=as.numeric(value)]
    load(file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0","updated_stan_data.RData"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# calc obs catch and cpue (1000 hooks) reference points
    obs_max_catch = max(updated_stan_data$obs_removals,na.rm=TRUE)
    obs_min_catch = min(updated_stan_data$obs_removals,na.rm=TRUE)
    obs_catch_overall_avg = mean(updated_stan_data$obs_removals,na.rm=TRUE)

    obs_max_cpue_early = max((updated_stan_data$obs_removals/(1e5*updated_stan_data$effort))[c(1:10)],na.rm=TRUE)
    obs_min_cpue_early = min((updated_stan_data$obs_removals/(1e5*updated_stan_data$effort))[c(1:10)],na.rm=TRUE)
    obs_cpue_early_avg = mean((updated_stan_data$obs_removals/(1e5*updated_stan_data$effort))[c(1:10)],na.rm=TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# extract quantities and filter            

        # apply filtering
            shape_iter_samples = hmc_samples[name%in%c("shape")] %>%
                           setnames(.,"value","shape") %>%
                           .[,.(iter,shape)]
            logk_iter_samples = hmc_samples[name%in%c("logK")] %>%
                           setnames(.,"value","logK")
            q_iter_samples = hmc_samples[name%in%c("F")] %>%
                           setnames(.,"value","F") %>%
                           merge(.,data.table(effort=updated_stan_data$effort,row=1:length(updated_stan_data$effort)),by="row") %>%
                           .[,Fnum:=as.numeric(F)] %>%
                           .[,q:=Fnum/effort] %>%
                           .[,.(iter,row,q)]
            cpue_early_iter_samples = hmc_samples[name%in%c("x")] %>%
                        merge(.,logk_iter_samples[,.(iter,logK)],by="iter") %>%
                        merge(.,q_iter_samples[,.(iter,row,q)],by=c("iter","row")) %>%
                        .[,logK:=as.numeric(logK)] %>%
                        .[,value:=as.numeric(value)] %>%
                        .[,cpue:=value*exp(logK)*q/1e5] %>%
                        .[row<=10,.(avg_cpue_early=mean(cpue)),by=iter]
            x_iter_samples = hmc_samples[name%in%c("x")] %>%
                        merge(.,logk_iter_samples[,.(iter,logK)],by="iter") %>%
                        .[,logK:=as.numeric(logK)] %>%
                        .[,value:=as.numeric(value)] %>%
                        .[,n:=value*exp(logK)] %>%
                        .[,.(mindep=min(value),minN=min(n)),by=iter]
            f_iter_samples = hmc_samples[name%in%c("F")] %>%
                        .[,value:=as.numeric(value)] %>%
                        .[,.(maxF=max(value),avgF=mean(value)),by=iter]
            removals_iter_samples = hmc_samples[name%in%c("removals")] %>%
                        .[,value:=as.numeric(value)] %>%
                        .[,.(maxCatch=max(value),avgCatch=mean(value),minCatch=min(value)),by=iter]
            removals_early_iter_samples = hmc_samples[name%in%c("removals")&row<=10] %>%
                        .[,value:=as.numeric(value)] %>%
                        .[,.(avg_catch_early=mean(value)),by=iter]


            iter_samples = merge(x_iter_samples,f_iter_samples,by='iter') %>% 
                           merge(.,removals_iter_samples,by="iter") %>% 
                           merge(.,shape_iter_samples,by="iter") %>%
                           merge(.,removals_early_iter_samples,by="iter") %>%
                           merge(.,cpue_early_iter_samples,by="iter")

            iter_surve = iter_samples[mindep>0.01&minN>15000&shape<20]$iter
            iter_f = iter_samples[iter %in% iter_surve & maxF<2.5 & avgF<0.8]$iter
            iter_catch = iter_samples[iter %in% iter_surve & 
                                        # Overall bounds
                                    maxCatch > 0.5 * obs_max_catch &
                                    maxCatch < 5.0 * obs_max_catch &
                                    minCatch > 0.5 * obs_min_catch &
                                    avgCatch > 0.5 * obs_catch_overall_avg &
                                    avgCatch < 5.0 * obs_catch_overall_avg 
                                    ]$iter
            iter_cpue = iter_samples[iter %in% iter_catch & 
                                        # Overall bounds
                                    avg_cpue_early > 0.5 * obs_cpue_early_avg &
                                    avg_cpue_early < 5.0 * obs_cpue_early_avg 
                                    ]$iter

            leading_params = c("logK","r","shape","x0","qeff","rho","sigma_qdev","sigmap","sigmao_add")
            raw_leading_params = c("raw_logK","raw_logr","raw_logshape","raw_logx0","raw_logqeff","raw_rho","raw_sigma_qdev","raw_logsigmap","raw_sigmao_add")

            leading_params_initial_dt = dcast(hmc_samples[name %in% leading_params], 
                                iter ~ name, 
                                value.var = "value") %>%
                                .[,type:="initial"] %>%
                                setcolorder(.,c("iter","type",leading_params))
            leading_params_surve_dt = leading_params_initial_dt[iter %in% iter_surve] %>%
                               .[,type:="surve"] %>%
                                setcolorder(.,c("iter","type",leading_params))
            leading_params_f_dt = leading_params_initial_dt[iter %in% iter_f] %>%
                               .[,type:="f"] %>%
                                setcolorder(.,c("iter","type",leading_params))
            leading_params_catch_dt = leading_params_initial_dt[iter %in% iter_catch] %>%
                               .[,type:="catch"] %>%
                                setcolorder(.,c("iter","type",leading_params))
            leading_params_cpue_dt = leading_params_initial_dt[iter %in% iter_cpue] %>%
                               .[,type:="cpue"] %>%
                                setcolorder(.,c("iter","type",leading_params))
            leading_params_dt = rbind(leading_params_initial_dt,leading_params_surve_dt,leading_params_f_dt,leading_params_catch_dt,leading_params_cpue_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch","cpue"))] 

            raw_leading_params_initial_dt = dcast(hmc_samples[name %in% raw_leading_params], 
                                iter ~ name, 
                                value.var = "value") %>%
                                .[,type:="initial"] %>%
                                setcolorder(.,c("iter","type",raw_leading_params))
            raw_leading_params_surve_dt = raw_leading_params_initial_dt[iter %in% iter_surve] %>%
                               .[,type:="surve"] %>%
                                setcolorder(.,c("iter","type",raw_leading_params))
            raw_leading_params_f_dt = raw_leading_params_initial_dt[iter %in% iter_f] %>%
                               .[,type:="f"] %>%
                                setcolorder(.,c("iter","type",raw_leading_params))
            raw_leading_params_catch_dt = raw_leading_params_initial_dt[iter %in% iter_catch] %>%
                               .[,type:="catch"] %>%
                                setcolorder(.,c("iter","type",raw_leading_params))
            raw_leading_params_cpue_dt = raw_leading_params_initial_dt[iter %in% iter_cpue] %>%
                               .[,type:="cpue"] %>%
                                setcolorder(.,c("iter","type",raw_leading_params))
            raw_leading_params_dt = rbind(raw_leading_params_initial_dt,raw_leading_params_surve_dt,raw_leading_params_f_dt,raw_leading_params_catch_dt,raw_leading_params_cpue_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch","cpue"))]

            sample_ts_initial_dt = hmc_samples[name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="initial"]
            sample_ts_surve_dt = hmc_samples[iter %in% iter_surve&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="surve"]
            sample_ts_f_dt = hmc_samples[iter %in% iter_f&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="f"]                
            sample_ts_catch_dt = hmc_samples[iter %in% iter_catch&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="catch"]
            sample_ts_cpue_dt = hmc_samples[iter %in% iter_cpue&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="cpue"]
            sample_ts_dt = rbind(sample_ts_initial_dt,sample_ts_surve_dt,sample_ts_f_dt,sample_ts_catch_dt,sample_ts_cpue_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch","cpue"))]
            quantile_ts_dt = sample_ts_dt %>%
                             .[,.(lower=quantile(value,probs=0.025),mid=quantile(value,probs=0.5),upper=quantile(value,probs=0.975)),by=.(type,name,row)]

            # filtering results
                leading_params_dt[,.(total_removed=nrow(leading_params_dt[type=="initial"])-.N,initial=nrow(leading_params_dt[type=="initial"])),by=type]

            p_ts = quantile_ts_dt %>%
                    ggplot() +
                    facet_wrap(~name,scales="free_y") +
                    xlab("Time") +
                    ylab("Value") +
                    geom_hline(yintercept=0) +
                    geom_ribbon(aes(x=row,ymin=lower,ymax=upper,fill=type),alpha=0.3) +
                    geom_line(aes(x=row,y=mid,color=type),linewidth=1) +                    
                    viridis::scale_color_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            viridis::scale_fill_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            theme(
                                text = element_text(size = 20),
                                panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                                panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                                panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                                strip.background = element_rect(fill = "white"),
                                legend.key = element_rect(fill = "white")
                            )
            ggsave(filename=paste0("pushforward.ts.png"), plot =p_ts, device = "png", path = dir_plot,
  			scale = 1, width = 9, height = 6, units = c("in"),
  			dpi = 300, limitsize = TRUE)

            plot_samples = 500
            set.seed(123)
            p_leading = copy(leading_params_dt) %>%
                            .[, .SD[sample(.N, min(plot_samples, .N))], by = type] %>%
                            ggpairs(., columns = 3:ncol(leading_params_dt), aes(color = type, alpha = 0.4)) + 
                            viridis::scale_color_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            viridis::scale_fill_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            theme(
                                text = element_text(size = 20),
                                panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                                panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                                panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                                strip.background = element_rect(fill = "white"),
                                legend.key = element_rect(fill = "white")
                            )
            ggsave(filename=paste0("pairs.leading.png"), plot =p_leading, device = "png", path = dir_plot,
  			scale = 1, width = 12, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)

            set.seed(123)
            p_raw =  raw_leading_params_dt %>%
                            .[, .SD[sample(.N, min(plot_samples, .N))], by = type] %>%
                            ggpairs(., columns = 3:ncol(raw_leading_params_dt), aes(color = type, alpha = 0.4)) + 
                            viridis::scale_color_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            viridis::scale_fill_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            theme(
                                text = element_text(size = 20),
                                panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                                panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                                panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                                strip.background = element_rect(fill = "white"),
                                legend.key = element_rect(fill = "white")
                            )
            ggsave(filename=paste0("pairs.leading.raw.png"), plot =p_raw, device = "png", path = dir_plot,
  			scale = 1, width = 12, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE)        

#________________________________________________________________________________________________________________________________________________________________________________________________________
# develop new priors following pushforward
    # logK, r, shape, shape,qeff
        mv_params = leading_params_dt[type == "cpue", .(logK, r, shape,qeff)] %>%
                                .[, .(logK, log_r=log(r),log_shape=log(shape),log_qeff=log(qeff))] %>%
                                as.matrix(.)

        mv_prior_fitted = fit_multivariate_prior(
            param_matrix = mv_params,
            param_names = colnames(mv_params),
            use_mle = TRUE,
            verbose = TRUE
        )

        # splice in prior for x0
            pp_mv_prior_mean = c(mv_prior_fitted$mean[1:3],updated_stan_data$mv_prior_mean[4],mv_prior_fitted$mean[4])
            
            pp_mv_prior_corr = mv_prior_fitted$cor[1:3,1:3]
            pp_mv_prior_corr = rbind(pp_mv_prior_corr,rep(0,3),rep(0,3))
            pp_mv_prior_corr = cbind(pp_mv_prior_corr,rep(0,5),rep(0,5))
            pp_mv_prior_corr[4,4] = pp_mv_prior_corr[5,5] = 1
            pp_mv_prior_corr[1,5] = pp_mv_prior_corr[5,1] = mv_prior_fitted$cor[1,4]
            
            pp_mv_prior_sd = c(mv_prior_fitted$sd[1:3],updated_stan_data$mv_prior_sd[4],mv_prior_fitted$sd[4])
            names(pp_mv_prior_sd) = names(pp_mv_prior_mean) = colnames(pp_mv_prior_corr) = rownames(pp_mv_prior_corr) = c("logK","log_r","log_shape","log_x0","log_qeff")

        # reset cor to 0 for non logK-qeff and r-shape pairs
            pp_mv_prior_corr[which(abs(pp_mv_prior_corr)<0.1)] = 0
    
    # rho and sigma_qdev
        mv_qdev_params = leading_params_dt[type == "cpue", .(rho,sigma_qdev)] %>%
                                .[, .(atan_rho=atanh(rho), log_sigma_qdev=log(sigma_qdev))] %>%
                                as.matrix(.)

        mv_qdev_prior_fitted = fit_multivariate_prior(
            param_matrix = mv_qdev_params,
            param_names = colnames(mv_qdev_params),
            use_mle = TRUE,
            verbose = TRUE
        )

        pp_mv_qdev_prior_mean = mv_qdev_prior_fitted$mean
        pp_mv_qdev_prior_sd = mv_qdev_prior_fitted$sd
        pp_mv_qdev_prior_cor = mv_qdev_prior_fitted$cor

#________________________________________________________________________________________________________________________________________________________________________________________________________
# generate samples from new prior and compare
    # logK, r, shape, shape, x0, qeff
        new_samples_mv_cov = diag(pp_mv_prior_sd) %*% pp_mv_prior_corr %*% diag(pp_mv_prior_sd)
        new_samples_mv = MASS::mvrnorm(plot_samples, mu = pp_mv_prior_mean, Sigma = new_samples_mv_cov)

    # rho, sigma_qdev
        new_samples_qdev_mv_cov = diag(pp_mv_qdev_prior_sd) %*% pp_mv_qdev_prior_cor %*% diag(pp_mv_qdev_prior_sd)
        new_samples_qdev_mv = MASS::mvrnorm(plot_samples, mu = pp_mv_qdev_prior_mean, Sigma = new_samples_qdev_mv_cov)

        new_samples_dt = as.data.table(cbind(new_samples_mv,new_samples_qdev_mv)) %>%
        .[,type:="updated prior"] %>%
        .[,.(type=type,logK=logK,r=exp(log_r),shape=exp(log_shape),x0=exp(log_x0),log_qeff=log_qeff,atan_rho=atan_rho,sigma_qdev=exp(log_sigma_qdev))] %>%
                              .[,type:=factor(type,levels=c("initial","surve","f","catch","cpue","updated prior"))]

        plot_new_samples_dt = leading_params_dt[type%in%c("initial","catch","cpue"),.(type,logK, r, shape,x0,qeff,rho,sigma_qdev)] %>%
                              .[,.(type=type,logK=logK,r=r,shape=shape,x0=x0,log_qeff=log(qeff),atan_rho=atanh(rho),sigma_qdev=sigma_qdev)]%>%
                              rbind(.,new_samples_dt[,.(type,logK, r, shape,x0,log_qeff,atan_rho,sigma_qdev)]) %>%
                              .[,type:=factor(type,levels=c("initial","surve","f","catch","cpue","updated prior"))]
    
        set.seed(123)
        p_newp =  plot_new_samples_dt %>%
                            .[, .SD[sample(.N, min(plot_samples, .N))], by = type] %>%
                            ggpairs(., columns = 2:ncol(plot_new_samples_dt), aes(color = type, alpha = 0.4)) + 
                            viridis::scale_color_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            viridis::scale_fill_viridis("Filter", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                            theme(
                                text = element_text(size = 20),
                                panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                                panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
                                panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
                                strip.background = element_rect(fill = "white"),
                                legend.key = element_rect(fill = "white")
                            )
            ggsave(filename=paste0("pairs.leading.updated_prior.png"), plot =p_newp, device = "png", path = dir_plot,
  			scale = 1, width = 12, height = 9, units = c("in"),
  			dpi = 300, limitsize = TRUE) 
#________________________________________________________________________________________________________________________________________________________________________________________________________
# save new stan.data object
    old_updated_stan_data = updated_stan_data
    updated_stan_data = old_updated_stan_data
    
    updated_stan_data$mv_prior_mean = unname(pp_mv_prior_mean)
    updated_stan_data$mv_prior_sd = unname(pp_mv_prior_sd)
    updated_stan_data$mv_prior_corr = unname(pp_mv_prior_corr)
    updated_stan_data$mv_qdev_prior_mean = unname(pp_mv_qdev_prior_mean)
    updated_stan_data$mv_qdev_prior_sd = unname(pp_mv_qdev_prior_sd)
    updated_stan_data$mv_qdev_prior_corr = unname(pp_mv_qdev_prior_cor)
    updated_stan_data$fit_to_data = 1L

    save(updated_stan_data,file=file.path(dir_output,"updated_stan_data.RData"))

    run_id = "00-bspm_estqsimple_softdep_mvprior_x0-updatedprior"
    stan_data.list = as.list(rep(NA,length(updated_stan_data)))
                for(i in 1:length(updated_stan_data)){
                    tmp = updated_stan_data[[i]]
                    if(is.matrix(tmp)){
                        tmp = as.data.table(tmp) %>%
                        .[,row:=1:.N] %>%
                        melt(.,id.vars="row") %>%
                        .[,col:=as.numeric(factor(variable,levels=unique(variable)))] %>%
                        .[,name:=names(updated_stan_data[i])] %>%
                        .[,.(name,row,col,value)]
                    } else if(is.vector(tmp)&length(tmp)>1){
                        tmp = as.data.table(tmp) %>%
                        .[,row:=1:.N] %>%
                        melt(.,id.vars="row") %>%
                        .[,name:=names(updated_stan_data[i])] %>%
                        .[,col:=NA] %>%
                        .[,.(name,row,col,value)]
                    } else {
                        tmp=data.table(name=names(updated_stan_data[i]),row=NA,col=NA,value=updated_stan_data[[i]])
                    }
                    stan_data.list[[i]] = tmp

                    if(length(grep("PriorSD",names(updated_stan_data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorSD"] %>%
                            .[,name:=gsub("PriorSD_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMean",names(updated_stan_data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorMean"] %>%
                            .[,name:=gsub("PriorMean_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMin",names(updated_stan_data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorMin"] %>%
                            .[,name:=gsub("PriorMin_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMax",names(updated_stan_data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorMax"] %>%
                            .[,name:=gsub("PriorMax_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else {
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="Data"] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    }
                    stan_data.list[[i]] = tmp
                }

    updated_stan_data_dt = rbindlist(stan_data.list)
    fwrite(updated_stan_data_dt,file.path(dir_output,"updated_stan_data.csv"))
