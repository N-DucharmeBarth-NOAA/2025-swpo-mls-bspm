    

# Nicholas Ducharme-Barth
# 2025/06/02
# Wrapper function to run the stan model

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    # run_label = "0000_test_0"
    # exec_name = "01_stochastic_spm_vX"
    # seed  = 321
    # chains = 5
    # n_thin = 1
    # iter_keep = 200
    # burnin.prop = 0.5
    # adapt_delta = 0.8
    # max_treedepth = 10


    fit_rstan = function(stan.data,stan_c,run_label,exec_name,n_cores,seed,chains,n_thin,iter_keep,burnin.prop,adapt_delta,max_treedepth,
                        silent=TRUE,stan_code_dir = file.path(".","code","Stan"),stan_save_dir=file.path(".","data","output","model_runs")){
        # packages
        require(rstan)
        require(parallelly)
        require(data.table)
        require(magrittr)
        require(bayesplot)
        require(loo)

        # set-up
        if(missing(n_cores)){
            n_cores = parallelly::availableCores(omit = 1,logical=FALSE)
        }
        options(mc.cores = n_cores)
        # rstan_options(auto_write = TRUE)

        # define run characteristics
            run_id = paste0(run_label,"-",exec_name)
            run_number = strsplit(run_label,"_")[[1]][1]
            run_name = strsplit(run_label,"_")[[1]][2]
            run_retro = as.numeric(strsplit(run_label,"_")[[1]][2])
            iter_total = n_thin*iter_keep + burnin.prop*n_thin*iter_keep
            file.name = paste0(stan_code_dir,exec_name,".stan") # path to stan file

        # modify data if doing a retrospective run
            if(run_retro>0){
                tmp_index = stan.data$index
                tmp_sigmao = stan.data$sigmao_mat
                tmp_T = stan.data$T

                tmp_idx = (tmp_T-(run_retro-1)):tmp_T
                tmp_index[tmp_idx,] = -999
                tmp_sigmao[tmp_idx,] = -999

                stan.data$index = tmp_index
                stan.data$sigmao_mat = tmp_sigmao
            }

        # define the inits
            set.seed(seed)
            stan.inits = replicate(chains,stan_inits_func(Tm1 = (stan.data$T-1)),simplify=FALSE)
        
        # compile if necessary
        if(missing(stan_c)){
            rstan_options(auto_write = TRUE)
            stan_c = stan_model(file=file.name, model_name = exec_name)
        }

        # fit the model
            
            fit = sampling(object=stan_c,
                data = stan.data,
                init = stan.inits,
                chains = chains,
                warmup = burnin.prop*n_thin*iter_keep,
                iter = iter_total,
                thin = n_thin,
                seed = seed,
                control = list(adapt_delta = adapt_delta,max_treedepth=max_treedepth))

        # summarize hmc samples
            hmc_samples = as.data.table(fit) %>%
                        .[,iter:=1:.N] %>%
                        melt(.,id.vars="iter") %>%
                        .[,chain:=ceiling(iter/iter_keep)] %>%
                        .[,chain_iter:=iter-((chain-1)*iter_keep)] %>%
                        .[,.(iter,chain,chain_iter,variable,value)]

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
                
            np = as.data.table(nuts_params(fit)) %>%
                    setnames(.,c("Chain","Iteration","Value"),c("chain","chain_iter","value")) %>%
                    dcast(.,chain+chain_iter~Parameter) %>%
                    setnames(.,c("treedepth__", "divergent__","accept_stat__", "stepsize__","n_leapfrog__","energy__"),c("treedepth","divergent","acceptance","stepsize","leapfrog","energy")) %>%
                    .[,.(chain,chain_iter,treedepth,divergent,acceptance,stepsize,leapfrog,energy)]
                    
            hmc_samples = merge(hmc_samples,names_dt,by="variable") %>%
                            merge(.,np,by=c("chain","chain_iter")) %>%
                            .[!(name%in%c("C","sum1","sum2","p","sigmao2","dev","epsp"))] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value,treedepth,divergent,acceptance,stepsize,leapfrog,energy)]
            
            # add rhat and neff
            rhat_vec = rhat(fit)
            rhat_dt = data.table(variable=names(rhat_vec),rhat=rhat_vec) %>%
                            .[,split:=sapply(variable,function(x)ifelse(length(grep("[",x,fixed=TRUE))==0,0,1))] %>%
                            .[,split_mat:=sapply(variable,function(x)ifelse(length(grep(",",x,fixed=TRUE))==0,0,1))] %>%
                            .[,row:=as.numeric(NA)] %>%
                            .[,col:=as.numeric(NA)] %>%
                            .[split_mat==1,row:=sapply(variable,function(x)as.numeric(strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][1]))] %>%
                            .[split_mat==1,col:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][2])))] %>%
                            .[split==1&split_mat==0,row:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(as.character(x),"[",fixed=TRUE)[[1]][2])))] %>%
                            .[,name:=variable] %>%
                            .[split==1,name:=sapply(variable,function(x)strsplit(as.character(x),"[",fixed=TRUE)[[1]][1])] %>%
                            .[!(name%in%c("C","sum1","sum2","p","sigmao2","dev"))] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,variable,name,row,col,rhat)]
            
            neff_vec = neff_ratio(fit)
            neff_dt = data.table(variable=names(neff_vec),neff=neff_vec) %>%
                            .[,split:=sapply(variable,function(x)ifelse(length(grep("[",x,fixed=TRUE))==0,0,1))] %>%
                            .[,split_mat:=sapply(variable,function(x)ifelse(length(grep(",",x,fixed=TRUE))==0,0,1))] %>%
                            .[,row:=as.numeric(NA)] %>%
                            .[,col:=as.numeric(NA)] %>%
                            .[split_mat==1,row:=sapply(variable,function(x)as.numeric(strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][1]))] %>%
                            .[split_mat==1,col:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(strsplit(as.character(x),"[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][2])))] %>%
                            .[split==1&split_mat==0,row:=sapply(variable,function(x)as.numeric(gsub("]","",strsplit(as.character(x),"[",fixed=TRUE)[[1]][2])))] %>%
                            .[,name:=variable] %>%
                            .[split==1,name:=sapply(variable,function(x)strsplit(as.character(x),"[",fixed=TRUE)[[1]][1])] %>%
                            .[!(name%in%c("C","sum1","sum2","p","sigmao2","dev"))] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,variable,name,row,col,neff)]

        # summarize input data
                stan_data.list = as.list(rep(NA,length(stan.data)))
                for(i in 1:length(stan.data)){
                    tmp = stan.data[[i]]
                    if(is.matrix(tmp)){
                        tmp = as.data.table(tmp) %>%
                        .[,row:=1:.N] %>%
                        melt(.,id.vars="row") %>%
                        .[,col:=as.numeric(factor(variable,levels=unique(variable)))] %>%
                        .[,name:=names(stan.data[i])] %>%
                        .[,.(name,row,col,value)]
                    } else if(is.vector(tmp)&length(tmp)>1){
                        tmp = as.data.table(tmp) %>%
                        .[,row:=1:.N] %>%
                        melt(.,id.vars="row") %>%
                        .[,name:=names(stan.data[i])] %>%
                        .[,col:=NA] %>%
                        .[,.(name,row,col,value)]
                    } else {
                        tmp=data.table(name=names(stan.data[i]),row=NA,col=NA,value=stan.data[[i]])
                    }
                    stan_data.list[[i]] = tmp

                    if(length(grep("PriorSD",names(stan.data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorSD"] %>%
                            .[,name:=gsub("PriorSD_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMean",names(stan.data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorMean"] %>%
                            .[,name:=gsub("PriorMean_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMin",names(stan.data[i])))>0){
                        tmp = stan_data.list[[i]] %>%
                            .[,type:="PriorMin"] %>%
                            .[,name:=gsub("PriorMin_","",name,fixed=TRUE)] %>%
                            .[,run_id := run_id] %>%
                            .[,.(run_id,type,name,row,col,value)]
                    } else if(length(grep("PriorMax",names(stan.data[i])))>0){
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
            stan_data = rbindlist(stan_data.list)

        # summarize inits
                stan_inits.list = as.list(rep(NA,length(stan.inits)))

                strip_inits_name = function(x,target_names = names(stan.inits[[1]])){
                    new_x = x
                    for(i in 1:length(target_names)){
                        new_x = gsub(target_names[i],"",new_x,fixed=TRUE)
                    }
                    return(new_x)
                }

                for(i in 1:length(stan.inits)){
                    tmp = data.table(name=names(unlist(stan.inits[[i]])),value=unlist(stan.inits[[i]])) %>%
                          .[,row:=as.numeric(strip_inits_name(name))] %>%
                          .[,chain:=i] %>%
                          .[,run_id:=run_id] %>%
                          .[,.(run_id,name,chain,row,value)]
                    for(j in 1:nrow(tmp)){
                        if(!is.na(tmp$row[j])){
                            tmp$name[j] = gsub(as.character(tmp$row[j]),"",tmp$name[j],fixed=TRUE)
                        }
                    }
                    stan_inits.list[[i]] = tmp
                }
       
            stan_inits = rbindlist(stan_inits.list)

        # summarize fit
            n_par = length(unlist(stan.inits[[1]]))
    
            fit_summary = data.table(run_id = run_id,
                                exec = exec_name,
                                run_label = run_label,
                                run_number = run_number,
                                run_name = run_name,
                                run_retro = run_retro,
                                run_time_total=max(rowSums(get_elapsed_time(fit))),
                                run_time_avg=mean(rowSums(get_elapsed_time(fit))),
                                n_par = n_par,
                                low_bfmi = length(get_low_bfmi_chains(fit))/chains,
                                divergent = get_num_divergent(fit)/(chains*iter_keep),
                                treedepth = get_num_max_treedepth(fit)/(chains*iter_keep),
                                max_rhat = max(rhat(fit)[1:n_par]),
                                min_neff = min(neff_vec[1:n_par],na.rm=TRUE),
                                prop_rhat_105 = mean(rhat(fit)[1:n_par]<= 1.05),
                                prop_rhat_101 = mean(rhat(fit)[1:n_par]<= 1.01),
                                prop_neff_01 = mean(neff_ratio(fit)[1:n_par]<= 0.1),
                                prop_neff_400 = mean(neff_ratio(fit)[1:n_par]<= max(500/hmc_samples$iter)))
            
            # leading parameters
                leading_parameters = c("logK","r","sigmao_sc","sigmap","n")
                tmp_leading = hmc_samples[name%in%leading_parameters] %>%
                              .[name=="sigmao_sc",name:="sigmao"] %>%
                              .[,.(value=median(value)),by=.(run_id,name)] %>%
                              dcast(.,run_id~name)
                leading_parameters = c("logK","r","sigmao","sigmap","n")
                missing_leading = leading_parameters[!(leading_parameters %in% colnames(tmp_leading))]
                if(length(missing_leading)>0)
                {
                    # make dummy data
                    dummy = matrix(NA,ncol=length(missing_leading))
                    colnames(dummy) = missing_leading
                    tmp_leading = cbind(tmp_leading,as.data.table(dummy))
                }
                setcolorder(tmp_leading,c("run_id",leading_parameters))
                if(is.na(tmp_leading$sigmao)){
                    tmp_leading$sigmao = stan_data[name=="sigmao_sc"]$value
                }
                if(is.na(tmp_leading$sigmap)){
                    tmp_leading$sigmap = stan_data[name=="sigmap"]$value
                }
            # derived quants
                derived_quants = ssp_derived_quants(hmc_samples,stan_data,output="percentiles",percentile=0.5)

            # priors
                tmp_priors = stan_data[type!="Data"] %>%
                            .[,.(run_id,type,name,value)] %>%
                            dcast(.,run_id~type+name)

                possible_priors = c(paste0("PriorMean_",c("logK","logr","logshape","logsigmap","logsigmao")),
                                    paste0("PriorSD_",c("logK","logr","logshape","logsigmap","logsigmao")),
                                    "PriorMin_logk","PriorMax_logk")
                
                missing_priors = possible_priors[!(possible_priors %in% colnames(tmp_priors))]
                if(length(missing_priors)>0)
                {
                    # make dummy data
                    dummy = matrix(NA,ncol=length(missing_priors))
                    colnames(dummy) = missing_priors
                    tmp_priors = cbind(tmp_priors,as.data.table(dummy))
                }
                setcolorder(tmp_priors,c("run_id",possible_priors))

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
                # n_obs
                n_obs = tmp_likelihood %>% 
                          .[lambda>0] %>%
                          na.omit(.) %>%
                          .[iter==1] %>%
                          nrow(.)
                fit_summary$n_obs = n_obs
                # looic
                fit_summary$looic = loo$estimates["looic",1]
                # p_loo
                fit_summary$p_loo = loo$estimates["p_loo",1]
                # elpd_loo
                fit_summary$elpd_loo = loo$estimates["elpd_loo",1]
                # pareto-k
                fit_summary$prop_pareto_k = mean(loo$diagnostics$pareto_k<= 0.7)
                
                # total_ll
                total_ll = tmp_likelihood %>% 
                          .[lambda>0] %>%
                          na.omit(.) %>%
                          .[,value:=value*lambda] %>%
                          .[,.(tll=sum(value)),by=.(I,iter)]
                fit_summary$total_ll = median(total_ll$tll)
                # index_ll
                index_ll = tmp_likelihood %>% 
                          na.omit(.) %>%
                          .[,.(index_ll=sum(value)),by=.(I,iter)] %>%
                          .[,.(index_ll = median(index_ll)),by=I] %>%
                          .[order(I)]
                index_ll = as.data.table(matrix(index_ll$index_ll,ncol=nrow(index_ll)))
                colnames(index_ll) = paste0("index_ll_",1:ncol(index_ll))
                fit_summary = cbind(fit_summary,index_ll)



            # fits to data
                # index_rmse
                tmp_rmse = ssp_calc_rmse(hmc_samples,stan_data) %>% .[order(I)]
                
                fit_summary$median_rmse = median(tmp_rmse[lambdas>0,.(mean(rmse)),by=iter]$V1)

                index_rmse = as.data.table(matrix(tmp_rmse[,.(median(rmse)),by=I]$V1,ncol=max(tmp_rmse$I)))
                colnames(index_rmse) = paste0("index_rmse_",1:ncol(index_rmse))
                fit_summary = cbind(fit_summary,index_rmse)


            fit_summary = merge(fit_summary,tmp_leading,by="run_id") %>%
                          merge(.,derived_quants,by="run_id") %>%
                          merge(.,tmp_priors,by="run_id") %>%
                          setnames(.,"n","shape")

        
        # write-out summaries
            dir_output = file.path(stan_save_dir,run_label)
            dir.create(dir_output,recursive=TRUE)

            settings_dt = data.table(run_id=run_id,
                                    run_number=run_number,
                                    run_name=run_name,
                                    run_retro=run_retro,
                                    run_label=run_label,
                                    exec_name=exec_name,
                                    n_cores=n_cores,
                                    seed=seed,
                                    chains=chains,
                                    n_thin=n_thin,
                                    iter_keep=iter_keep,
                                    burnin.prop=burnin.prop,
                                    iter_total = iter_total,
                                    adapt_delta=adapt_delta,
                                    max_treedepth=max_treedepth,
                                    silent=silent,
                                    stan_code_dir = stan_code_dir)

            fwrite(fit_summary,file.path(dir_output,"fit_summary.csv"))
            fwrite(hmc_samples,file.path(dir_output,"hmc_samples.csv"))
            fwrite(settings_dt,file.path(dir_output,"settings.csv"))
            fwrite(neff_dt,file.path(dir_output,"neff.csv"))
            fwrite(rhat_dt,file.path(dir_output,"rhat.csv"))
            fwrite(stan_data,file.path(dir_output,"stan_data.csv"))
            fwrite(stan_inits,file.path(dir_output,"stan_inits.csv"))

        # return
            if(silent){
                return()
            } else {
                return(fit)
            }
    }
    
   