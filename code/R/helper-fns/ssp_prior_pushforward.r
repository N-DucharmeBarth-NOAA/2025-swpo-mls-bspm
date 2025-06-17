    

# ISC SHARKWG
# 2024/06/11
# Prior pushforward for stochastic surplus production model

# Copyright (c) 2024 ISC SHARKWG
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_prior_pushforward = function(ssp_summary,stan_data,settings){
      set.seed(settings$seed)
      chains = settings$chains
      n_samples = settings$iter_keep*chains

      # Check if using multivariate priors
      has_mv_prior = "mv_prior_mean" %in% stan_data$name
      
      if (has_mv_prior) {
            # Extract multivariate prior parameters
            mv_mean = stan_data[name == "mv_prior_mean" & type == "Data"][order(row)]$value
            mv_sd = stan_data[name == "mv_prior_sd" & type == "Data"][order(row)]$value
            
            # Reconstruct correlation matrix
            mv_corr_data = stan_data[name == "mv_prior_corr" & type == "Data"]
            max_row = max(mv_corr_data$row)
            max_col = max(mv_corr_data$col)
            mv_corr_mat = matrix(NA, max_row, max_col)
            
            for(i in 1:nrow(mv_corr_data)) {
                  mv_corr_mat[mv_corr_data$row[i], mv_corr_data$col[i]] = mv_corr_data$value[i]
            }
            
            # Create covariance matrix
            mv_cov = diag(mv_sd) %*% mv_corr_mat %*% diag(mv_sd)
            
            # Sample from multivariate normal
            library(MASS)  # for mvrnorm
            mv_samples = mvrnorm(n_samples, mv_mean, mv_cov)
            
            # Extract individual parameters
            logK = mv_samples[, 1]
            log_r = mv_samples[, 2]
            log_shape = mv_samples[, 3]
            
            r = exp(log_r)
            shape = exp(log_shape)
            
            # Calculate raw parameters for compatibility
            raw_logK = (logK - mv_mean[1]) / mv_sd[1]
            raw_logr = (log_r - mv_mean[2]) / mv_sd[2]
            raw_logshape = (log_shape - mv_mean[3]) / mv_sd[3]
            
      } else {
            # Original univariate approach
            logK = rnorm(n_samples, ssp_summary$PriorMean_logK, ssp_summary$PriorSD_logK)
            raw_logK = (logK - ssp_summary$PriorMean_logK) / ssp_summary$PriorSD_logK
            r = exp(rnorm(n_samples, ssp_summary$PriorMean_logr, ssp_summary$PriorSD_logr))
            raw_logr = (log(r) - ssp_summary$PriorMean_logr) / ssp_summary$PriorSD_logr
            
            shape = exp(rnorm(n_samples, stan_data[name == "logshape" & type == "PriorMean"]$value, 
                              stan_data[name == "logshape" & type == "PriorSD"]$value))
            raw_logshape = (log(shape) - stan_data[name == "logshape" & type == "PriorMean"]$value) / 
                        stan_data[name == "logshape" & type == "PriorSD"]$value
      }

      sigmap = exp(rnorm(n_samples,stan_data[name=="logsigmap"&type=="PriorMean"]$value,stan_data[name=="logsigmap"&type=="PriorSD"]$value))
      raw_logsigmap = (log(sigmap) - stan_data[name=="logsigmap"&type=="PriorMean"]$value)/stan_data[name=="logsigmap"&type=="PriorSD"]$value

      raw_sigmao_add = abs(rnorm(n_samples))
      sigmao_add = raw_sigmao_add*stan_data[name=="sigmao_add"&type=="PriorSD"]$value
      sigmao_sc = stan_data[name=="sigmao_input"]$value + sigmao_add

      raw_sigmaf = abs(rnorm(n_samples))
      F = raw_F = matrix(NA,nrow=n_samples,ncol=stan_data[name=="T"]$value-1)
      sigmaf = raw_sigmaf*stan_data[name=="sigmaf"]$value


      raw_epsp = epsp = matrix(NA,nrow=n_samples,ncol=stan_data[name=="T"]$value)
      for(j in 1:ncol(epsp)){
            epsp[,j] = rlnorm(n_samples,log(1.0)-(sigmap^2)/2,sigmap)
            if(j<ncol(epsp)){
                  raw_F[,j] = abs(rnorm(n_samples))
                  F[,j] = raw_F[,j]*sigmaf
            }
      }

      raw_epsp = (log(epsp) - (sigmap^2)/2)/sigmap

      sigmap2 = sigmap^2
      n = shape
      dmsy = (1/n)^(1/(n-1))
      h = 2*dmsy
      m = r*h/4
      g = (n^(n/(n-1)))/(n-1)

      x = matrix(NA,nrow=n_samples,ncol=ncol(epsp))
      x[,1] = epsp[,1]
      T = ncol(x)
      removals = matrix(NA,nrow=n_samples,ncol=stan_data[name=="T"]$value-1)

      for(i in 1:nrow(x)){
            for(t in 2:T) {
                  if(x[i,t-1]<=dmsy[i]){
                  x[i,t] = (x[i,t-1] + r[i] * x[i,t-1] * (1 - x[i,t-1]/h[i]))*(exp(-F[i,t-1]))*epsp[i,t]
                  removals[i,t-1] = ((x[i,t-1] + r[i] * x[i,t-1] * (1 - x[i,t-1]/h[i])))*epsp[i,t]*(1-exp(-F[i,t-1]))*exp(logK[i])
                  } else {
                  x[i,t] = (x[i,t-1] + g[i] * m[i] * x[i,t-1] * (1 - x[i,t-1]^(n[i]-1)))*(exp(-F[i,t-1]))*epsp[i,t]
                  removals[i,t-1] = ((x[i,t-1] + g[i] * m[i] * x[i,t-1] * (1 - x[i,t-1]^(n[i]-1))))*epsp[i,t]*(1-exp(-F[i,t-1]))*exp(logK[i]);
                  }
            }
      }

      vector_var = c("raw_logK","raw_logr","raw_logsigmap","raw_logshape","logK","r","m","sigmap","sigmap2","sigmao_sc","shape","n","dmsy","h","g","sigmao_add","sigmaf")
      vector_dt.list = as.list(rep(NA,length(vector_var)))
      for(i in 1:length(vector_dt.list)){
            vector_dt.list[[i]] = data.table(value=get(vector_var[i])) %>%
                    .[,run_id:=ssp_summary$run_id] %>%
                    .[,iter:=1:n_samples] %>%
                    .[,chain:=sort(rep(1:chains,n_samples/chains))] %>%
                    .[,chain_iter:=rep(1:(n_samples/chains),chains)] %>%
                    .[,variable:=vector_var[i]] %>%
                    .[,name:=vector_var[i]] %>%
                    .[,row:=as.numeric(NA)] %>%
                    .[,col:=as.numeric(NA)] %>%
                    .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)] %>%
                    .[order(chain,iter)]
      }
      vector_dt = rbindlist(vector_dt.list)

      matrix_var = c("raw_epsp","x","removals")
      matrix_dt.list = as.list(rep(NA,length(matrix_var)))
      for(i in 1:length(matrix_dt.list)){
            matrix_dt.list[[i]] = as.data.table(get(matrix_var[i])) %>%
                                  .[,iter:=1:n_samples] %>%
                                  .[,chain:=sort(rep(1:chains,n_samples/chains))] %>%
                                  .[,chain_iter:=rep(1:(n_samples/chains),chains)] %>%
                                  melt(.,id.vars=c("iter","chain","chain_iter")) %>%
                                  .[,row:=as.numeric(gsub("V","",variable))] %>%
                    .[,run_id:=ssp_summary$run_id] %>%
                    .[,variable:=paste0(matrix_var[i],"[",row,"]")] %>%
                    .[,name:=matrix_var[i]] %>%
                    .[,col:=as.numeric(NA)] %>%
                    .[,.(run_id,iter,chain,chain_iter,variable,name,row,col,value)] %>%
                    .[order(chain,iter)]
      }
      matrix_dt = rbindlist(matrix_dt.list)

      out = rbind(vector_dt,matrix_dt)

      return(out)      
}    
