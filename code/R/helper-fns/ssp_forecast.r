# ISC SHARKWG
# 2024/06/11
# forecast - simplified version for when removals exist in samples_dt

# Copyright (c) 2024 ISC SHARKWG
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_forecast = function(ssp_summary,samples_dt,stan_data,settings,sub_sample_prop=0.1,forecast_years=5,forecast_type=c("Catch","U","MSY","Umsy")[1],avg_years=3,scalar=1,resample_raw_epsp="TRUE",exclude_last=0){
      set.seed(settings$seed)
      tmp_samples = samples_dt[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)]
      sub_sample_n = round(max(tmp_samples$iter)*sub_sample_prop)
      if(sub_sample_n>max(tmp_samples$iter)){sub_sample_n=max(tmp_samples$iter)}
      if(sub_sample_n<1){sub_sample_n=1}
      sub_samples_dt = tmp_samples[iter %in% sample(1:max(tmp_samples$iter),sub_sample_n,replace=FALSE)]
      iter_sampled = sub_samples_dt[name=="x"&row==1]$iter

      logK = sub_samples_dt[name=="logK"]$value
      r = sub_samples_dt[name=="r"]$value
      x = dcast(sub_samples_dt[name=="x",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
      raw_epsp = dcast(sub_samples_dt[name=="raw_epsp",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
      sigmap2 = sub_samples_dt[name=="sigmap2"]$value
      sigmap = sqrt(sub_samples_dt[name=="sigmap2"]$value)
      n = sub_samples_dt[name=="n"]$value
      dmsy = (1/n)^(1/(n-1))
      h = 2*dmsy
      m = sub_samples_dt[name=="m"]$value
      g = sub_samples_dt[name=="g"]$value
      T = ncol(x)

      # Handle removals matrix
      removals_mat_a = dcast(sub_samples_dt[name=="removals",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
      if(ncol(removals_mat_a)==T-1){
            sigmac = stan_data[name=="sigmac"]$value
            mu_catch = log(stan_data[name=="obs_removals"&row==T]$value) - 0.5*sigmac^2
            removals_mat_a = cbind(removals_mat_a,rlnorm(nrow(removals_mat_a),mu_catch,sigmac))
      }
      removals_mat = cbind(removals_mat_a,matrix(NA,nrow=length(iter_sampled),ncol=forecast_years))

      # Generate forecast process error
      raw_forecast_epsp = forecast_epsp = matrix(NA,nrow=length(iter_sampled),ncol=forecast_years)
      if(resample_raw_epsp == "TRUE"){
        for(i in 1:nrow(raw_forecast_epsp)){
                raw_forecast_epsp[i,] = sample(raw_epsp[i,], size = forecast_years, replace = TRUE)
                forecast_epsp[i,] = exp(raw_forecast_epsp[i,]*sigmap[i]-(sigmap[i]^2)/2)
        }
      } else {
        for(j in 1:ncol(forecast_epsp)){
                forecast_epsp[,j] = rlnorm(length(iter_sampled),log(1.0)-(sigmap^2)/2,sigmap)
        }
        raw_forecast_epsp = (log(forecast_epsp) - (sigmap^2)/2)/sigmap
      }

      forecast_x = matrix(NA,nrow=length(iter_sampled),ncol=forecast_years)
      
      # First year forecast
      for(i in 1:nrow(forecast_x)){
            C = min(c(exp(log(removals_mat[i,T]) - logK[i]),x[i,T]))
            
            if(x[i,T]<=dmsy[i]){
                  forecast_x[i,1] = (x[i,T] + r[i] * x[i,T] * (1 - x[i,T]/h[i]) - C)*forecast_epsp[i,1]
            } else {
                  forecast_x[i,1] = (x[i,T] + g[i] * m[i] * x[i,T] * (1 - x[i,T]^(n[i]-1)) - C)*forecast_epsp[i,1]
            }

            if(forecast_years == 1){
                removals_mat[i,(T+forecast_years)] = C*exp(logK[i])
            }
      }
      
      # Determine forecast catch/exploitation levels
      if(forecast_type == "Catch"){
            catch_vec = rep(NA,length(iter_sampled))
            for(i in 1:length(iter_sampled)){
                  catch_vec[i] = mean(removals_mat_a[i,][((length(removals_mat_a[i,])-(avg_years-1)):(length(removals_mat_a[i,])-exclude_last))])*scalar
            }
      } else if(forecast_type == "U"){
            u_vec = rep(NA,length(iter_sampled))
            for(i in 1:length(iter_sampled)){
                  u_vec[i] = mean(removals_mat_a[i,][((length(removals_mat_a[i,])-(avg_years-1)):(length(removals_mat_a[i,])-exclude_last))]/(x[i,][((length(x[i,])-(avg_years-1)):(length(x[i,])-exclude_last))]*exp(logK[i])))*scalar
            }
      } else if(forecast_type == "Umsy"){
            u_vec = (m/dmsy)*scalar
      } else {
            catch_vec = m * exp(logK) * scalar
      }

      # Forecast for years 2:forecast_years
      if(forecast_years > 1){
            for(i in 1:nrow(forecast_x)){
                  for(t in 2:forecast_years) {
                        if(forecast_type == "U"|forecast_type == "Umsy"){
                              C = forecast_x[i,t-1]*u_vec[i]
                        } else {
                              C = min(c(exp(log(catch_vec[i]) - logK[i]),forecast_x[i,t-1]))
                        }
                        removals_mat[i,(T+t-1)] = C*exp(logK[i])
                        if(forecast_x[i,t-1]<=dmsy[i]){
                              forecast_x[i,t] = (forecast_x[i,t-1] + r[i] * forecast_x[i,t-1] * (1 - forecast_x[i,t-1]/h[i]) - C)*forecast_epsp[i,t]
                        } else {
                              forecast_x[i,t] = (forecast_x[i,t-1] + g[i] * m[i] * forecast_x[i,t-1] * (1 - forecast_x[i,t-1]^(n[i]-1)) - C)*forecast_epsp[i,t]
                        }
                  }
                  if(forecast_type == "U"|forecast_type == "Umsy"){
                              removals_mat[i,(T+forecast_years)] = forecast_x[i,forecast_years]*u_vec[i]*exp(logK[i])
                  } else {
                              removals_mat[i,(T+forecast_years)] = min(c(exp(log(catch_vec[i]) - logK[i]),forecast_x[i,forecast_years]))*exp(logK[i])
                  }
            }
      }

      # Repackage results
      new_x = cbind(x,forecast_x)
      new_raw_epsp = cbind(raw_epsp,raw_forecast_epsp)

      new_x_dt = as.data.table(new_x) %>%
                  .[,iter:=iter_sampled] %>%
                            melt(.,id.vars=c("iter")) %>%
                            .[,row:=as.numeric(gsub("V","",variable))] %>%
              .[,run_id:=ssp_summary$run_id] %>%
              .[,variable:=paste0("x[",row,"]")] %>%
              .[,name:="x"] %>%
              .[,col:=as.numeric(NA)] %>%
              merge(.,unique(samples_dt[,.(iter,chain,chain_iter,treedepth,divergent,acceptance,stepsize,leapfrog,energy)]),by="iter") %>%
              .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value,treedepth,divergent,acceptance,stepsize,leapfrog,energy)] %>%
              .[order(chain,iter)]

      new_raw_epsp_dt = as.data.table(new_raw_epsp) %>%
                  .[,iter:=iter_sampled] %>%
                            melt(.,id.vars=c("iter")) %>%
                            .[,row:=as.numeric(gsub("V","",variable))] %>%
              .[,run_id:=ssp_summary$run_id] %>%
              .[,variable:=paste0("raw_epsp[",row,"]")] %>%
              .[,name:="raw_epsp"] %>%
              .[,col:=as.numeric(NA)] %>%
              merge(.,unique(samples_dt[,.(iter,chain,chain_iter,treedepth,divergent,acceptance,stepsize,leapfrog,energy)]),by="iter") %>%
              .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value,treedepth,divergent,acceptance,stepsize,leapfrog,energy)] %>%
              .[order(chain,iter)]

      removals_dt = as.data.table(removals_mat) %>%
                  .[,iter:=iter_sampled] %>%
                            melt(.,id.vars=c("iter")) %>%
                            .[,row:=as.numeric(gsub("V","",variable))] %>%
              .[,run_id:=ssp_summary$run_id] %>%
              .[,variable:=paste0("removals[",row,"]")] %>%
              .[,name:="removals"] %>%
              .[,col:=as.numeric(NA)] %>%
              merge(.,unique(samples_dt[,.(iter,chain,chain_iter,treedepth,divergent,acceptance,stepsize,leapfrog,energy)]),by="iter") %>%
              .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value,treedepth,divergent,acceptance,stepsize,leapfrog,energy)] %>%
              .[order(chain,iter)]
      
      # Combine results (removals case)
      out = samples_dt[!(name%in%c("raw_epsp","x","removals"))] %>%
            rbind(.,new_x_dt) %>%
            rbind(.,new_raw_epsp_dt) %>%
            rbind(.,removals_dt)

      return(out)
}
