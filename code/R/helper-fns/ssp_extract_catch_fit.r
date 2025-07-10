ssp_extract_catch_fit = function(ssp_summary, samples_dt, stan_data, settings, sub_sample_prop=0.1, calc_std="TRUE"){
      require(data.table)
      require(magrittr)

      set.seed(settings$seed)
      tmp_samples = samples_dt[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)]
      sub_sample_n = round(max(tmp_samples$iter)*sub_sample_prop)
      if(sub_sample_n>max(tmp_samples$iter)){sub_sample_n=max(tmp_samples$iter)}
      if(sub_sample_n<1){sub_sample_n=1}
      samples_dt = tmp_samples[iter %in% sample(1:max(tmp_samples$iter),sub_sample_n,replace=FALSE)]
      iter_sampled = samples_dt[name=="removals"&row==1]$iter
      
      # Get observed catch data
      obs_catch_dt = stan_data[name=="obs_removals",.(run_id,row,value)] %>%
                    .[,iter:=0] %>%
                    .[,metric:="obs_catch"] %>%
                    .[,.(run_id,metric,iter,row,value)] %>%
                    .[value==-999,value:=NA]
      
      # Get observation error for catch
      sigmac_data = stan_data[name=="sigmac"]
      if(nrow(sigmac_data) == 1) {
      sigmac = rep(sigmac_data$value, stan_data[name=="T"]$value)
      } else {
      sigmac = sigmac_data[order(row)]$value
      }
      obs_se_dt = data.table(run_id=ssp_summary$run_id,
                      metric="sigmac",
                      iter=0,
                      row=obs_catch_dt$row,
                      value=sigmac[obs_catch_dt$row])
      
      # Get predicted catch
      pred_catch_dt = samples_dt[name=="removals",.(run_id,iter,row,value)] %>%
                     .[,metric:="pred_catch"] %>%
                     .[,.(run_id,metric,iter,row,value)]
      
      error_type = "LN"
      if(nrow(samples_dt[name=="nu_catch",.(run_id,iter,row,value)])>0){
            error_type = "ST"
            nu_catch = samples_dt[name=="nu_catch",.(run_id,iter,row,value)]$value
      }
      
      # Generate posterior predicted catch
      ppd_catch_list = list()
      for(i in 1:length(iter_sampled)){
            iter_val = iter_sampled[i]
            pred_vals = pred_catch_dt[iter==iter_val]$value
            if(error_type == "LN"){
                  ppd_vals = rlnorm(length(pred_vals), 
                        log(pred_vals) - 0.5*sigmac[1:length(pred_vals)]^2, 
                        sigmac[1:length(pred_vals)])
            } else if(error_type == "ST"){
                  ppd_vals = exp(rt(length(pred_vals), df = nu_catch) * sigmac[1:length(pred_vals)] + log(pred_vals))
            }

            
            ppd_catch_list[[i]] = data.table(run_id=ssp_summary$run_id,
                                           metric="ppd_catch",
                                           iter=iter_val,
                                           row=1:length(ppd_vals),
                                           value=ppd_vals)
      }
      ppd_catch_dt = rbindlist(ppd_catch_list)
      
      # Calculate PIT residuals
      y = as.vector(na.omit(obs_catch_dt$value))[1:(length(na.omit(obs_catch_dt$value))-1)] # no catch prediction for the terminal year
      yrep = dcast(ppd_catch_dt, run_id+metric+row~iter, value.var="value")
      pit_run_id = yrep$run_id
      pit_metric = rep("pit_residual", nrow(yrep))
      pit_iter = rep(0, nrow(yrep))
      pit_row = yrep$row
      yrep_mat = yrep %>%
                .[,run_id:=NULL] %>%
                .[,metric:=NULL] %>%
                .[,row:=NULL] %>%
                as.matrix(.)
      
      pit = rep(NA, length(y))
      for(i in 1:length(y)){
            if(!is.na(y[i])){
                  tmp_ecdf = ecdf(yrep_mat[i,])
                  pit[i] = tmp_ecdf(y[i])
            }
      }
      
      pit_dt = data.table(run_id=pit_run_id,
                         metric=pit_metric,
                         iter=pit_iter,
                         row=pit_row,
                         value=pit)
      
      # Calculate ordinary residuals and standardized residuals if requested
      y = as.vector(na.omit(obs_catch_dt$value))[1:(length(na.omit(obs_catch_dt$value))-1)] # no catch prediction for the terminal year
      keep = which(!is.na(y))
      
      if(calc_std=="TRUE"){
            leverage_mat = std_residual_mat = residual_mat = matrix(NA, nrow=length(y), ncol=length(iter_sampled))
            colnames(leverage_mat) = colnames(std_residual_mat) = colnames(residual_mat) = iter_sampled
      } else {
            residual_mat = matrix(NA, nrow=length(y), ncol=length(iter_sampled))
            colnames(residual_mat) = iter_sampled
      }
      
      y_scaled_norm = matrix(scale(y)/sqrt(sum(scale(y)^2)), nrow=length(y), ncol=1)
      
      for(i in 1:length(iter_sampled)){
            tmp_ypred = pred_catch_dt[iter==iter_sampled[i]][keep]$value
            residual_mat[,i] = y - tmp_ypred
            
            if(calc_std=="TRUE"){
                  tmp_ypred_scaled_norm = matrix(scale(tmp_ypred)/sqrt(sum(scale(tmp_ypred)^2)), nrow=length(tmp_ypred), ncol=1)
                  H = tmp_ypred_scaled_norm %*% matrix(scale(y)/sqrt(sum(scale(y)^2)), ncol=length(y), nrow=1)
                  
                  std_resid_init = as.vector(y_scaled_norm) - as.vector(tmp_ypred_scaled_norm)
                  mse = mean(sum(std_resid_init^2))
                  leverage_mat[,i] = diag(H)
                  std_residual_mat[,i] = std_resid_init/sqrt(mse*(1-diag(H)))
            }
      }
      
      # Convert matrices to data.tables
      if(calc_std=="TRUE"){
            leverage_dt = as.data.table(leverage_mat) %>%
                  .[,row:=pred_catch_dt[iter==iter_sampled[1]][keep]$row] %>%
                  melt(., id.vars="row") %>%
                  .[,variable:=sapply(variable, function(x)as.numeric(gsub("V","",x)))] %>%
                  setnames(., "variable", "iter") %>%
                  .[,run_id:=ssp_summary$run_id] %>%
                  .[,metric:="leverage"] %>%
                  .[,.(run_id,metric,iter,row,value)] %>%
                  .[,iter:=as.numeric(iter)]

            std_residual_dt = as.data.table(std_residual_mat) %>%
                  .[,row:=pred_catch_dt[iter==iter_sampled[1]][keep]$row] %>%
                  melt(., id.vars="row") %>%
                  .[,variable:=sapply(variable, function(x)as.numeric(gsub("V","",x)))] %>%
                  setnames(., "variable", "iter") %>%
                  .[,run_id:=ssp_summary$run_id] %>%
                  .[,metric:="std_residual"] %>%
                  .[,.(run_id,metric,iter,row,value)] %>%
                  .[,iter:=as.numeric(iter)]
      }

      residual_dt = as.data.table(residual_mat) %>%
            .[,row:=pred_catch_dt[iter==iter_sampled[1]][keep]$row] %>%
            melt(., id.vars="row") %>%
            .[,variable:=sapply(variable, function(x)as.numeric(gsub("V","",x)))] %>%
            setnames(., "variable", "iter") %>%
            .[,run_id:=ssp_summary$run_id] %>%
            .[,metric:="residual"] %>%
            .[,.(run_id,metric,iter,row,value)] %>%
            .[,iter:=as.numeric(iter)]
      
      # Combine all results
      if(calc_std=="TRUE"){
            out = rbind(obs_catch_dt,
                       obs_se_dt,
                       pred_catch_dt,
                       ppd_catch_dt,
                       pit_dt,
                       leverage_dt,
                       std_residual_dt,
                       residual_dt)
      } else {
            out = rbind(obs_catch_dt,
                       obs_se_dt,
                       pred_catch_dt,
                       ppd_catch_dt,
                       pit_dt,
                       residual_dt)
      }
      
      return(out)
}
