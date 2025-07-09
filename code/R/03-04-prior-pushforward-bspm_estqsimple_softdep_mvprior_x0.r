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
    dir_output = file.path(proj_dir,"data","output","pushforward","bspm_estqsimple_softdep_mvprior_x0")
    dir.create(dir_output,showWarnings = FALSE, recursive = TRUE)
    dir_plot = file.path(proj_dir,"plots","pushforward","bspm_estqsimple_softdep_mvprior_x0")
    dir.create(dir_plot,showWarnings = FALSE, recursive = TRUE)
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    

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
    
    catch_annual = catch_dt[,.(total_catch = sum(Obs * 1000)), by = .(year = floor(Time))]
    setorder(catch_annual, year)
    
    # Merge with effort data
    catch_effort_annual = merge(catch_annual, effort_dt, by.x = "year", by.y = "time", all.x = TRUE)
    
    # Fill missing effort with mean (if any)
    if(any(is.na(catch_effort_annual$effort_scaled))) {
        catch_effort_annual[is.na(effort_scaled), effort_scaled := mean(effort_scaled, na.rm = TRUE)]
    }

    obs_max_catch = max(catch_effort_annual$total_catch,na.rm=TRUE)
    obs_min_catch = min(catch_effort_annual$total_catch,na.rm=TRUE)
    obs_catch_overall_avg = mean(catch_effort_annual$total_catch,na.rm=TRUE)
    obs_catch_early_avg = mean(catch_effort_annual[year<1962]$total_catch,na.rm=TRUE)

    # develop initial catchability prior based on nominal scaled cpue and same logK prior
        scaled_cpue = (catch_effort_annual$total_catch[1:10]/(catch_effort_annual$effort_scaled[1:10]))
        nsim_q = 1e4
        scaled_q = sample(scaled_cpue,nsim_q,replace=TRUE)/rlnorm(nsim_q,log(4e6),1.25) 

        scaled_q_fn = function(par){-sum(dnorm(log(scaled_q), mean = par[1], sd = par[2], log = TRUE))}
    qeff_pars_catch =  nlminb(c(mean(log(scaled_q)), sd(log(scaled_q))), scaled_q_fn)$par

    # Make input mv prior for logk, logr, logshape, logx0

    bio_params_dt = fread(file.path(proj_dir,"data","output","bspm_parameter_priors_filtered.csv"))
    lrs_corr = cor(bio_params_dt[,.(log(rmax),log(shape))])[1,2]

    mv_prior_mean = c(log(4e6),as.vector(apply(bio_params_dt[,.(log(rmax),log(shape))],2,mean)),0,qeff_pars_catch[1])
    mv_prior_sd = c(1.25,as.vector(apply(bio_params_dt[,.(log(rmax),log(shape))],2,sd)),0.025,qeff_pars_catch[2])
    mv_prior_corr = matrix(c(1,0,0,0,0,
                             0,1,lrs_corr,0,0,
                             0,lrs_corr,1,0,0,
                             0,0,0,1,0,
                             0,0,0,0,1),ncol=5)

    # Bivariate rho/sigma_qdev prior
    mv_qdev_mean = c(2,0.9)
    mv_qdev_sd = c(0.75,0.5)
    mv_qdev_cor = matrix(c(1,0,0,1),ncol=2)


#________________________________________________________________________________________________________________________________________________________________________________________________________
# load additional cpue indices
    au_cpue_dt = fread(file.path(proj_dir,"data","input","AU_cpue.csv"))
    nz_cpue_dt = fread(file.path(proj_dir,"data","input","NZ_cpue.csv"))
    obs_cpue_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP.csv"))
    obs_cpue_no_PF_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_no_PF.csv"))
    obs_cpue_PF_only_dt = fread(file.path(proj_dir,"data","input","obs-idx-with_OP_PF_only.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare data matrices    
    time_years = catch_effort_annual$year
    n_years = length(time_years)
    
    index_mat = matrix(-999, nrow = n_years, ncol = 6)
    se_mat = matrix(-999, nrow = n_years, ncol = 6)
    mean_se = rep(NA,6)

    # dwfn cpue
    cpue_years = floor(cpue_dt$Time)
    mean_se[1] = mean(cpue_dt$SE,na.rm=TRUE)
    for(i in 1:nrow(cpue_dt)) {
        year_idx = which(time_years == cpue_years[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 1] = cpue_dt$Obs[i]
            se_mat[year_idx, 1] = cpue_dt$SE[i]/mean_se[1]
        }
    }

    # australia cpue
    mean_se[2] = mean(au_cpue_dt$se_log,na.rm=TRUE)
    for(i in 1:nrow(au_cpue_dt)) {
        year_idx = which(time_years == au_cpue_dt$year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 2] = au_cpue_dt$obs[i]/mean(au_cpue_dt$obs)
            se_mat[year_idx, 2] = au_cpue_dt$se_log[i]/mean_se[2]
        }
    }

    # new zealand cpue
    mean_se[3] = mean(nz_cpue_dt$se_log,na.rm=TRUE)
    for(i in 1:nrow(nz_cpue_dt)) {
        year_idx = which(time_years == nz_cpue_dt$year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 3] = nz_cpue_dt$ob[i]/mean(nz_cpue_dt$ob)
            se_mat[year_idx, 3] = nz_cpue_dt$se_log[i]/mean_se[3]
        }
    }

    # Observer CPUE
    # calc SE from quantiles
    obs_cpue_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[4] = mean(obs_cpue_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_dt)) {
        year_idx = which(time_years == obs_cpue_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 4] = obs_cpue_dt$Estimate[i]/mean(obs_cpue_dt$Estimate)
            se_mat[year_idx, 4] = obs_cpue_dt$se_from_quantiles[i]/mean_se[4]
        }
    }

    # calc SE from quantiles
    obs_cpue_no_PF_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[5] = mean(obs_cpue_no_PF_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_no_PF_dt)) {
        year_idx = which(time_years == obs_cpue_no_PF_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 5] = obs_cpue_no_PF_dt$Estimate[i]/mean(obs_cpue_no_PF_dt$Estimate)
            se_mat[year_idx, 5] = obs_cpue_no_PF_dt$se_from_quantiles[i]/mean_se[5]
        }
    }

    # calc SE from quantiles
    obs_cpue_PF_only_dt[, se_from_quantiles := (Q97.5 - Q2.5) / (2 * 1.96)]
    mean_se[6] = mean(obs_cpue_PF_only_dt$se_from_quantiles,na.rm=TRUE)
    for(i in 1:nrow(obs_cpue_PF_only_dt)) {
        year_idx = which(time_years == obs_cpue_PF_only_dt$Year[i])
        if(length(year_idx) > 0) {
            index_mat[year_idx, 6] = obs_cpue_PF_only_dt$Estimate[i]/mean(obs_cpue_PF_only_dt$Estimate)
            se_mat[year_idx, 6] = obs_cpue_PF_only_dt$se_from_quantiles[i]/mean_se[6]
        }
    }
#________________________________________________________________________________________________________________________________________________________________________________________________________
# prepare multivariate prior parameters
    # Validate correlation matrices are positive definite
    if(min(eigen(mv_prior_corr)$values) < 1e-10) {
        warning("5D correlation matrix has very small eigenvalues, adding regularization")
        mv_prior_corr = mv_prior_corr + diag(1e-8, 5)
    }
    
    if(min(eigen(mv_qdev_cor)$values) < 1e-10) {
        warning("2D correlation matrix has very small eigenvalues, adding regularization")  
        mv_qdev_cor = mv_qdev_cor + diag(1e-8, 2)
    }
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# compile executable
    exec_name = "bspm_estqsimple_softdep_mvprior_x0" # bspm_estq_softdep_mvprior # bspm_estq_optimized
    stan_c = stan_model(file=file.path(proj_dir,"code","Stan",paste0(exec_name,".stan")), model_name = exec_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set-up model inputs

    lambda_vec = c(1,0,0,0,0,0)
    sigmao_input = mean_se[1]

 # Define effort parameters
    n_step = 1  # years per period for catchability
    n_periods = ceiling((n_years-1) / n_step)
            
    # specify data and priors
    tmp.data.priors = list(
        # 3D multivariate prior (logK, log_r, log_shape,log_x0,log_qeff)
        mv_prior_mean = mv_prior_mean,
        mv_prior_sd = mv_prior_sd,
        mv_prior_corr = mv_prior_corr,
        
        # Bivariate rho/sigma_qdev prior  
        mv_qdev_prior_mean = mv_qdev_mean,
        mv_qdev_prior_sd = mv_qdev_sd,
        mv_qdev_prior_corr = mv_qdev_cor,
        
        # Other priors
        PriorMean_logsigmap = -2.9311037,
        PriorSD_logsigmap = 0.2661089,
        PriorSD_sigmao_add = 0.2,
        prior_depletion_meanlog = 0,
        prior_depletion_sdlog = 1
    )

    tmp_data = list(
        T = as.integer(n_years),
        I = as.integer(ncol(index_mat)),
        index = index_mat,
        sigmao_mat = se_mat,
        lambdas = lambda_vec,
        t_dep = 37, # model starts t=1 in 1952, 37 is 1988
        use_depletion_prior = 0L, # Set to 1L to use prior, 0L to skip it
        fit_to_data = 0L, # Set to 1L to fit to data, 0L to skip likelihoods
        epsilon = 1e-10,
        
        # Effort-based parameters
        effort = catch_effort_annual$effort_scaled,
        n_step = as.integer(n_step),
        n_periods = as.integer(n_periods),
        
        # Observation error parameters
        sigmao_input = sigmao_input,
        obs_removals = catch_effort_annual$total_catch,
        sigmac = 0.5
    )
                     
    stan.data = c(tmp_data, tmp.data.priors)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# initial prior-pushforward

    chains = 5
    n_periods_val = if("n_periods" %in% names(stan.data)) stan.data$n_periods else NULL
    set.seed(333)
    stan.inits = replicate(chains, 
                                stan_inits_func(Tm1 = (stan.data$T-1), 
                                                    n_periods = n_periods_val,
                                                    exec_name = exec_name), 
                                simplify=FALSE)
    
    # options(mc.cores = min(c(chains,parallelly::availableCores(omit = 1,logical=FALSE))))
    options(mc.cores = 5)
    initial_prior_push = sampling(object=stan_c,
                data = stan.data,
                init = stan.inits,
                chains = chains,
                warmup = 250,
                iter = 3250,
                thin = 1,
                seed = 321,
                control = list(adapt_delta = 0.99,max_treedepth=12))


        # summarize hmc samples
            hmc_samples = as.data.table(initial_prior_push) %>%
                        .[,iter:=1:.N] %>%
                        melt(.,id.vars="iter") %>%
                        .[,.(iter,variable,value)]

            names_dt = data.table(variable=unique(hmc_samples$variable)) %>%
                            .[,split:=sapply(variable,function(x)ifelse(length(grep("[",x,fixed=TRUE))==0,0,1))] %>%
                            .[,split_mat:=sapply(variable,function(x)ifelse(length(grep(",",x,fixed=TRUE))==0,0,1))] %>%
                            .[,row:=as.numeric(NA)] %>%
                            .[,col:=as.numeric(NA)] %>%
                            .[split_mat==1,row:=sapply(variable,function(x)as.numeric(strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][1]))] %>%
                            .[split_mat==1,col:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][2])))] %>%
                            .[split==1&split_mat==0,row:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(as.character(x),"[",fixed=TRUE)[[1]][2])))] %>%
                            .[,name:=variable] %>%
                            .[split==1,name:=sapply(variable,function(x)strsplit(as.character(x),"[",fixed=TRUE)[[1]][1])] %>%
                            .[,.(variable,name,row,col)]
                

            hmc_samples = merge(hmc_samples,names_dt,by="variable") %>%
                            .[!(name%in%c("C","sum1","sum2","p","sigmao2","dev","epsp"))] %>%
                            .[,.(iter,variable,name,row,col,value)]
            fwrite(hmc_samples,file=file.path(dir_output,"hmc_samples.csv"))

#________________________________________________________________________________________________________________________________________________________________________________________________________
# extract quantities and filter            

        # apply filtering
            shape_iter_samples = hmc_samples[name%in%c("shape")] %>%
                           setnames(.,"value","shape") %>%
                           .[,.(iter,shape)]
            logk_iter_samples = hmc_samples[name%in%c("logK")] %>%
                           setnames(.,"value","logK") 
            x_iter_samples = hmc_samples[name%in%c("x")] %>%
                        merge(.,logk_iter_samples[,.(iter,logK)],by="iter") %>%
                        .[,n:=value*exp(logK)] %>%
                        .[,.(mindep=min(value),minN=min(n)),by=iter]
            f_iter_samples = hmc_samples[name%in%c("F")] %>%
                        .[,.(maxF=max(value),avgF=mean(value)),by=iter]
            removals_iter_samples = hmc_samples[name%in%c("removals")] %>%
                        .[,.(maxCatch=max(value),avgCatch=mean(value),minCatch=min(value)),by=iter]
            removals_early_iter_samples = hmc_samples[name%in%c("removals")&row<=10] %>%
                        .[,.(avg_catch_early=mean(value)),by=iter]


            iter_samples = merge(x_iter_samples,f_iter_samples,by='iter') %>% 
                           merge(.,removals_iter_samples,by="iter") %>% 
                           merge(.,shape_iter_samples,by="iter") %>%
                           merge(.,removals_early_iter_samples,by="iter")

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
            leading_params_dt = rbind(leading_params_initial_dt,leading_params_surve_dt,leading_params_f_dt,leading_params_catch_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch"))] 

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
            raw_leading_params_dt = rbind(raw_leading_params_initial_dt,raw_leading_params_surve_dt,raw_leading_params_f_dt,raw_leading_params_catch_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch"))]

            sample_ts_initial_dt = hmc_samples[name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="initial"]
            sample_ts_surve_dt = hmc_samples[iter %in% iter_surve&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="surve"]
            sample_ts_f_dt = hmc_samples[iter %in% iter_f&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="f"]                
            sample_ts_catch_dt = hmc_samples[iter %in% iter_catch&name%in%c("x","removals"),.(iter,name,row,value)] %>%
                                .[,type:="catch"]
            sample_ts_dt = rbind(sample_ts_initial_dt,sample_ts_surve_dt,sample_ts_f_dt,sample_ts_catch_dt) %>%
                                .[,type:=factor(type,levels=c("initial","surve","f","catch"))]
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
        mv_params = leading_params_dt[type == "catch", .(logK, r, shape,qeff)] %>%
                                .[, .(logK, log_r=log(r),log_shape=log(shape),log_qeff=log(qeff))] %>%
                                as.matrix(.)

        mv_prior_fitted = fit_multivariate_prior(
            param_matrix = mv_params,
            param_names = colnames(mv_params),
            use_mle = TRUE,
            verbose = TRUE
        )

        # splice in prior for x0
            pp_mv_prior_mean = c(mv_prior_fitted$mean[1:3],stan.data$mv_prior_mean[4],mv_prior_fitted$mean[4])
            
            pp_mv_prior_corr = mv_prior_fitted$cor[1:3,1:3]
            pp_mv_prior_corr = rbind(pp_mv_prior_corr,rep(0,3),rep(0,3))
            pp_mv_prior_corr = cbind(pp_mv_prior_corr,rep(0,5),rep(0,5))
            pp_mv_prior_corr[4,4] = pp_mv_prior_corr[5,5] = 1
            pp_mv_prior_corr[1,5] = pp_mv_prior_corr[5,1] = mv_prior_fitted$cor[1,4]
            
            pp_mv_prior_sd = c(mv_prior_fitted$sd[1:3],stan.data$mv_prior_sd[4],mv_prior_fitted$sd[4])
            names(pp_mv_prior_sd) = names(pp_mv_prior_mean) = colnames(pp_mv_prior_corr) = rownames(pp_mv_prior_corr) = c("logK","log_r","log_shape","log_x0","log_qeff")

        # reset cor to 0 for non logK-qeff and r-shape pairs
            pp_mv_prior_corr[which(abs(pp_mv_prior_corr)<0.1)] = 0
    
    # rho and sigma_qdev
        mv_qdev_params = leading_params_dt[type == "catch", .(rho,sigma_qdev)] %>%
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
                              .[,type:=factor(type,levels=c("initial","surve","f","catch","updated prior"))]

        plot_new_samples_dt = leading_params_dt[type%in%c("initial","catch"),.(type,logK, r, shape,x0,qeff,rho,sigma_qdev)] %>%
                              .[,.(type=type,logK=logK,r=r,shape=shape,x0=x0,log_qeff=log(qeff),atan_rho=atanh(rho),sigma_qdev=sigma_qdev)]%>%
                              rbind(.,new_samples_dt[,.(type,logK, r, shape,x0,log_qeff,atan_rho,sigma_qdev)]) %>%
                              .[,type:=factor(type,levels=c("initial","surve","f","catch","updated prior"))]
    
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
    updated_stan_data = stan.data
    
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

