# ISC SHARKWG
# 2024/06/11
# Time series of derived quantities

# Copyright (c) 2024 ISC SHARKWG
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_derived_quants_ts = function(ssp_summary,samples_dt,stan_data,settings,sub_sample_prop=0.1){
      set.seed(settings$seed)
      tmp_samples = samples_dt[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)]
      sub_sample_n = round(max(tmp_samples$iter)*sub_sample_prop)
      if(sub_sample_n>max(tmp_samples$iter)){sub_sample_n=max(tmp_samples$iter)}
      if(sub_sample_n<1){sub_sample_n=1}
      samples_dt = tmp_samples[iter %in% sample(1:max(tmp_samples$iter),sub_sample_n,replace=FALSE)]
      iter_sampled = samples_dt[name=="x"&row==1]$iter

      logK = samples_dt[name=="logK"]$value
      r = samples_dt[name=="r"]$value
      raw_epsp = dcast(samples_dt[name=="raw_epsp",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
      D = dcast(samples_dt[name=="x",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
      sigmap2 = samples_dt[name=="sigmap2"]$value
      sigmap = sqrt(samples_dt[name=="sigmap2"]$value)
      n = samples_dt[name=="n"]$value
      Dmsy = samples_dt[name=="dmsy"]$value
      Pmsy = Dmsy * exp(logK)
      h = samples_dt[name=="h"]$value
      m = samples_dt[name=="m"]$value
      g = samples_dt[name=="g"]$value
      T = ncol(D)
      Umsy = as.numeric(m / Dmsy)
      Fmsy = as.numeric(-log(-Umsy+1))

      # Check if F is estimated in the model
      F_estimated <- "F" %in% unique(samples_dt$name)

      # Check if effort data is available
      effort_available <- "effort" %in% unique(stan_data$name)
      if(effort_available){
            effort_data = stan_data[name=="effort"]$value
      }

      if(!("removals" %in% unique(samples_dt$name))){
            removals = matrix(NA,nrow=nrow(D),ncol=ncol(D))
            for(i in 1:nrow(removals)){
                  removals[i,1:T] = stan_data[name=="removals"]$value
            }
      } else {
           removals = dcast(samples_dt[name=="removals",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.) 
           if(ncol(removals)==T-1){
                  sigmac_data = stan_data[name=="sigmac"]
                  if(nrow(sigmac_data) == 1) {
                  sigmac = rep(sigmac_data$value, stan_data[name=="T"]$value)
                  } else {
                  sigmac = sigmac_data[order(row)]$value
                  }
                  mu_catch = log(stan_data[name=="obs_removals"&row==T]$value) - 0.5*sigmac[T]^2
                  removals = cbind(removals,rlnorm(nrow(removals),mu_catch,sigmac[T]))
           }
      }

      # Extract F values if estimated
      if(F_estimated){
            F_estimated_matrix = dcast(samples_dt[name=="F",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            # Extract edev and qdev if they exist
            if("edev" %in% unique(samples_dt$name)){
                  edev_estimated_matrix = dcast(samples_dt[name=="edev",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            } else {
                  edev_estimated_matrix = NULL
            }
            if("qdev" %in% unique(samples_dt$name)){
                  qdev_estimated_matrix = dcast(samples_dt[name=="qdev",.(iter,row,value)],iter~row) %>% .[,iter:=NULL] %>% as.matrix(.)
            } else {
                  qdev_estimated_matrix = NULL
            }
      }

      epsp = dev = matrix(NA,nrow=nrow(raw_epsp),ncol=ncol(raw_epsp))
      surplus_production = epsilon_p = D_Dmsy=P_Pmsy=U_Umsy=F_Fmsy=F=U = P = edev = qdev = nominal_cpue = matrix(NA,nrow=nrow(D),ncol=ncol(D))
      
      for(i in 1:nrow(raw_epsp)){
         dev[i,] = raw_epsp[i,] * sigmap[i]
         
         P[i,] = exp(logK[i]) * D[i,]
         
         for(j in 1:ncol(U)){
            if(F_estimated && j <= ncol(F_estimated_matrix)){
                  F[i,j] = F_estimated_matrix[i,j]
                  U[i,j] = 1 - exp(-F[i,j])
                  # Extract edev and qdev if available
                  if(!is.null(edev_estimated_matrix) && j <= ncol(edev_estimated_matrix)){
                        edev[i,j] = edev_estimated_matrix[i,j]
                  }
                  if(!is.null(qdev_estimated_matrix) && j <= ncol(qdev_estimated_matrix)){
                        qdev[i,j] = qdev_estimated_matrix[i,j]
                  }
            } else {
                  U[i,j] = as.numeric(min(c(removals[i,j]/P[i,j],0.9999)))
                  F[i,j] = if(j < ncol(U)) as.numeric(-log(-U[i,j]+1)) else NA
                  # edev and qdev remain NA (already initialized as NA)
            }
            
            # Calculate nominal CPUE if effort is available and removals are estimated
            if(effort_available && j <= length(effort_data) && !is.na(removals[i,j])){
                  # Convert effort to thousand hooks and calculate CPUE (catch per 1000 hooks)
                  effort_thousand_hooks = effort_data[j] * 1e6 / 10
                  if(effort_thousand_hooks > 0){
                        nominal_cpue[i,j] = removals[i,j] / effort_thousand_hooks
                  } else {
                        nominal_cpue[i,j] = NA
                  }
            } else {
                  nominal_cpue[i,j] = NA
            }
            
            if(j<ncol(U)){
                  epsilon_p[i,j] = epsp[i,j] = exp(dev[i,j]-sigmap2[i]/2)
            } else {
                  epsilon_p[i,j] = rlnorm(1,log(1.0)-sigmap2[i]/2,sigmap[i])
            }

            if(D[i,j] <= Dmsy[i]){
                  surplus_production[i,j] = r[i] * D[i,j] * (1 - D[i,j]/h[i]) * epsilon_p[i,j]
            } else {
                  surplus_production[i,j] = g[i] * m[i] * D[i,j] * (1 - D[i,j]^(n[i]-1)) * epsilon_p[i,j]
            }
         }
         
         # Calculate ratios (handling NAs appropriately)
         F_Fmsy[i,] = as.numeric(F[i,]/Fmsy[i])
         U_Umsy[i,] = as.numeric(U[i,]/Umsy[i])
         P_Pmsy[i,] = P[i,]/Pmsy[i]
         D_Dmsy[i,] = D[i,]/Dmsy[i]    
      }

      # pad epsp & dev
      matrix_var = c("raw_epsp","dev","epsilon_p","surplus_production","D_Dmsy","P_Pmsy","U_Umsy","F_Fmsy","F","U","P","D","removals","edev","qdev","nominal_cpue")
      matrix_dt.list = as.list(rep(NA,length(matrix_var)))
      for(i in 1:length(matrix_dt.list)){
            matrix_dt.list[[i]] = as.data.table(get(matrix_var[i])) %>%
                                  .[,iter:=iter_sampled] %>%
                                  melt(.,id.vars=c("iter")) %>%
                                  .[,row:=as.numeric(gsub("V","",variable))] %>%
                    .[,run_id:=ssp_summary$run_id] %>%
                    .[,variable:=paste0(matrix_var[i],"[",row,"]")] %>%
                    .[,name:=matrix_var[i]] %>%
                    .[,col:=as.numeric(NA)] %>%
                    merge(.,unique(samples_dt[,.(iter,chain,chain_iter)]),by="iter") %>%
                    .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)] %>%
                    .[order(chain,iter)]
      }

      out = rbindlist(matrix_dt.list)

      return(out)
}
