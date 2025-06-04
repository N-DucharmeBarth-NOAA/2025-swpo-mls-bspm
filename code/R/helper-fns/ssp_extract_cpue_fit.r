    

# ISC SHARKWG
# 2024/06/11
# observed-expected, posterior-predictive checking (ppc) density, ecdf, pit residuals

# Copyright (c) 2024 ISC SHARKWG
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_extract_cpue_fit = function(ssp_summary,samples_dt,stan_data,settings,sub_sample_prop=0.1,active="TRUE",calc_std="TRUE"){
      set.seed(settings$seed)
      tmp_samples = samples_dt[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)]
      sub_sample_n = round(max(tmp_samples$iter)*sub_sample_prop)
      if(sub_sample_n>max(tmp_samples$iter)){sub_sample_n=max(tmp_samples$iter)}
      if(sub_sample_n<1){sub_sample_n=1}
      samples_dt = tmp_samples[iter %in% sample(1:max(tmp_samples$iter),sub_sample_n,replace=FALSE)]
      iter_sampled = samples_dt[name=="x"&row==1]$iter
      lambdas = stan_data[name=="lambdas"]$value
      if(active=="TRUE"){
            active_lambdas = which(lambdas==1)
      } else {
            active_lambdas = 1:length(lambdas)
      }

      # get obs
            obs_cpue_dt = stan_data[name=="index",.(run_id,row,col,value)] %>%
                          .[col%in%active_lambdas] %>%
                          setnames(.,"col","index") %>%
                          .[,iter:=0] %>%
                          .[,metric:="obs_cpue"] %>%
                          .[,.(run_id,metric,iter,row,index,value)] %>%
                          .[value==-999,value:=NA]            
    
    if(!("sigmao" %in% unique(samples_dt$name))){
            obs_se_dt = stan_data[name=="sigmao_mat",.(run_id,row,col,value)] %>%
                          .[col%in%active_lambdas] %>%
                          setnames(.,"col","index") %>%
                          .[,iter:=0] %>%
                          .[,metric:="sigmao"] %>%
                          .[,.(run_id,metric,iter,row,index,value)] %>%
                          .[value==-999,value:=NA] %>%
                          .[,value:=value*stan_data[name=="sigmao_sc"]$value]
            obs_se_mat = dcast(obs_se_dt[,.(row,index,value)],row~index) %>% .[,row:=NULL] %>% as.matrix(.)

      # get pred cpue & posterior predicted cpue
            x = dcast(samples_dt[name=="x",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            q = dcast(samples_dt[name=="q",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            q = q[,active_lambdas]
            if(is.null(dim(q))){
                  q = matrix(q,ncol=1,nrow=length(q))
            }

            ppd_cpue_array = pred_cpue_array = array(NA,dim=c(nrow(x),ncol(x),ncol(q)),dimnames=list(iter=iter_sampled,row=1:ncol(x),index=active_lambdas))
            for(i in 1:nrow(x)){
                  for(j in 1:ncol(x)){
                        for(k in 1:ncol(q)){
                              pred_cpue_array[i,j,k] = x[i,j] * q[i,k]
                              if(!is.na(obs_se_mat[j,k])){
                                    ppd_cpue_array[i,j,k] = rlnorm(1,log(x[i,j] * q[i,k])-0.5*(obs_se_mat[j,k]^2),obs_se_mat[j,k])
                              }
                        }
                  }
            }
    } else {
            obs_se_dt = samples_dt[name=="sigmao",.(run_id,iter,row,col,value)] %>%
                        .[col%in%active_lambdas] %>%
                        setnames(.,"col","index") %>%
                        .[,metric:="sigmao"] %>%
                        .[,.(run_id,metric,iter,row,index,value)] %>%
                        .[value<0,value:=NA] 
      
      # get pred cpue & posterior predicted cpue
            x = dcast(samples_dt[name=="x",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            q = dcast(samples_dt[name=="q",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            q = q[,active_lambdas]
            if(is.null(dim(q))){
                  q = matrix(q,ncol=1,nrow=length(q))
            }

            ppd_cpue_array = pred_cpue_array = array(NA,dim=c(nrow(x),ncol(x),ncol(q)),dimnames=list(iter=iter_sampled,row=1:ncol(x),index=active_lambdas))
            for(i in 1:nrow(x)){
                  tmp_obs_se_mat = dcast(obs_se_dt[iter==iter_sampled[i],.(row,index,value)],row~index) %>% .[,row:=NULL] %>% as.matrix(.)

                  for(j in 1:ncol(x)){
                        for(k in 1:ncol(q)){
                              pred_cpue_array[i,j,k] = x[i,j] * q[i,k]
                              if(!is.na(tmp_obs_se_mat[j,k])){
                                    ppd_cpue_array[i,j,k] = rlnorm(1,log(x[i,j] * q[i,k])-0.5*(tmp_obs_se_mat[j,k]^2),tmp_obs_se_mat[j,k])
                              }
                        }
                  }
            }
    } 

            pred_cpue_dt = as.data.table(pred_cpue_array) %>%
                        .[,iter:=as.numeric(as.character(iter))] %>%
                        .[,row:=as.numeric(as.character(row))] %>%
                        .[,index:=as.numeric(as.character(index))] %>%
                        .[order(iter,row,index)] %>%
                        .[,run_id:=ssp_summary$run_id] %>%
                        .[,metric:="pred_cpue"] %>%
                        .[,.(run_id,metric,iter,row,index,value)]
            ppd_cpue_dt = as.data.table(ppd_cpue_array) %>%
                        .[,iter:=as.numeric(as.character(iter))] %>%
                        .[,row:=as.numeric(as.character(row))] %>%
                        .[,index:=as.numeric(as.character(index))] %>%
                        .[order(iter,row,index)] %>%
                        .[,run_id:=ssp_summary$run_id] %>%
                        .[,metric:="ppd_cpue"] %>%
                        .[,.(run_id,metric,iter,row,index,value)]
      
      # ecdf & pit-residual
            y = as.vector(na.omit(obs_cpue_dt$value))
            yrep = dcast(ppd_cpue_dt,run_id+metric+row+index~iter)
            pit_run_id = yrep$run_id
            pit_metric = rep("pit_residual",nrow(yrep))
            pit_iter = rep(0,nrow(yrep))
            pit_row = yrep$row
            pit_index = yrep$index
            yrep = yrep %>%
                  .[,run_id:=NULL] %>%
                  .[,metric:=NULL] %>% 
                  .[,row:=NULL] %>%
                  .[,index:=NULL] %>%
                  as.matrix(.) 
            
            pit = rep(NA,length(y))
            for(i in 1:length(y)){
                  # y_ecdf = cumsum(table(floor(yrep[,i]/0.0001)*0.0001))
                  # x_ecdf = as.numeric(names(y_ecdf))
                  # y_ecdf = y_ecdf/tail(y_ecdf,n=1)
                  # if(length(which(x_ecdf<y[i]))==0){
                  #       pit[i]=0
                  # } else{
                  #       pit[i] = y_ecdf[max(which(x_ecdf<y[i]))]
                  # }
                  # tmp_group = pit_group[i]
                  tmp_ecdf = ecdf(yrep[i,])
                  pit[i] = tmp_ecdf(y[i])                
            }

            pit_dt = data.table(run_id=pit_run_id,
                              metric=pit_metric,
                              iter=pit_iter,
                              row=pit_row,
                              index=pit_index,
                              value=pit)
      
      # ordinary & standardized residuals; leverage
            y = as.vector(na.omit(obs_cpue_dt$value))
            # y_vec_norm = matrix(y/sqrt(sum(y^2)),nrow=length(y),ncol=1)
            y_scaled_norm = matrix(scale(y)/sqrt(sum(scale(y)^2)),nrow=length(y),ncol=1)
            keep = which(!is.na(obs_cpue_dt$value))
            if(calc_std=="TRUE"){
                leverage_mat = std_residual_mat = residual_mat = matrix(NA,nrow=length(y),ncol=length(iter_sampled))
                colnames(leverage_mat) = colnames(std_residual_mat) = colnames(residual_mat) = iter_sampled
            } else {
                residual_mat = matrix(NA,nrow=length(y),ncol=length(iter_sampled))
                colnames(residual_mat) = iter_sampled
            }
            
            for(i in 1:length(iter_sampled)){
                  tmp_ypred = pred_cpue_dt[iter==iter_sampled[i]][keep]$value
                  residual_mat[,i] = y-tmp_ypred
                  if(calc_std=="TRUE"){
                    tmp_ypred_scaled_norm = matrix(scale(tmp_ypred)/sqrt(sum(scale(tmp_ypred)^2)),nrow=length(tmp_ypred),ncol=1)
                    H = tmp_ypred_scaled_norm%*%matrix(scale(y)/sqrt(sum(scale(y)^2)),ncol=length(y),nrow=1)
                    # H = H*(1/(sqrt(sum(tmp_ypred^2)))^2)
                    # y_pred_init = as.vector(H%*%y_scaled_norm)

                    std_resid_init = as.vector(y_scaled_norm)-as.vector(tmp_ypred_scaled_norm)
                    mse = mean(sum(std_resid_init^2))
                    leverage_mat[,i] = diag(H)
                    std_residual_mat[,i] = std_resid_init/sqrt(mse*(1-diag(H)))
                  }
            }

            if(calc_std=="TRUE"){
                leverage_dt = as.data.table(leverage_mat) %>%
                    .[,row:=pred_cpue_dt[iter==iter_sampled[i]][keep]$row] %>%
                    .[,index:=pred_cpue_dt[iter==iter_sampled[i]][keep]$index] %>%
                    melt(.,id.vars=c("row","index")) %>%
                    .[,variable:=sapply(variable,function(x)as.numeric(gsub("V","",x)))] %>%
                    setnames(.,"variable","iter") %>%
                    .[,run_id:=ssp_summary$run_id] %>%
                    .[,metric:="leverage"] %>%
                    .[,.(run_id,metric,iter,row,index,value)] %>%
                    .[,iter:=as.numeric(iter)]

                std_residual_dt = as.data.table(std_residual_mat) %>%
                    .[,row:=pred_cpue_dt[iter==iter_sampled[i]][keep]$row] %>%
                    .[,index:=pred_cpue_dt[iter==iter_sampled[i]][keep]$index] %>%
                    melt(.,id.vars=c("row","index")) %>%
                    .[,variable:=sapply(variable,function(x)as.numeric(gsub("V","",x)))] %>%
                    setnames(.,"variable","iter") %>%
                    .[,run_id:=ssp_summary$run_id] %>%
                    .[,metric:="std_residual"] %>%
                    .[,.(run_id,metric,iter,row,index,value)] %>%
                    .[,iter:=as.numeric(iter)]
            }

            residual_dt = as.data.table(residual_mat) %>%
                  .[,row:=pred_cpue_dt[iter==iter_sampled[i]][keep]$row] %>%
                  .[,index:=pred_cpue_dt[iter==iter_sampled[i]][keep]$index] %>%
                  melt(.,id.vars=c("row","index")) %>%
                  .[,variable:=sapply(variable,function(x)as.numeric(gsub("V","",x)))] %>%
                  setnames(.,"variable","iter") %>%
                  .[,run_id:=ssp_summary$run_id] %>%
                  .[,metric:="residual"] %>%
                  .[,.(run_id,metric,iter,row,index,value)] %>%
                  .[,iter:=as.numeric(iter)]
    if(calc_std=="TRUE"){
        out = rbind(obs_cpue_dt,
                  obs_se_dt,
                  pred_cpue_dt,
                  ppd_cpue_dt,
                  pit_dt,
                  leverage_dt,
                  std_residual_dt,
                  residual_dt)
    } else {
        out = rbind(obs_cpue_dt,
                  obs_se_dt,
                  pred_cpue_dt,
                  ppd_cpue_dt,
                  pit_dt,
                  residual_dt)
    }
      return(out)
}    
