    

# Nicholas Ducharme-Barth
# 2024/03/05
# Calculate model likelihood

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_calc_likelihood = function(hmc_samples,stan_data){
    require(data.table)
    require(magrittr)

    # define dimensions and extract data
        T = stan_data[name=="T"]$value
        I = stan_data[name=="I"]$value
        n_ll = T*I
        max_iter = max(hmc_samples$iter)

        lambdas = stan_data[name=="lambdas"]$value
        index = dcast(stan_data[name=="index",.(row,col,value)],row~col) %>%
                .[,row:=NULL] %>%
                as.matrix(.)
        x = dcast(hmc_samples[name%in%c("x"),.(iter,row,value)],row~iter) %>%
                .[,row:=NULL] %>%
                as.matrix(.)
        q = dcast(hmc_samples[name%in%c("q"),.(iter,row,value)],row~iter) %>%
                .[,row:=NULL] %>%
                as.matrix(.)
        if("sigmao" %in% unique(hmc_samples$name)){
            sigmao_dt = hmc_samples[name%in%c("sigmao"),.(iter,row,col,value)] %>%
                        setnames(.,c("value"),c("sigmao")) %>%
                        .[,sigmao2:=sigmao^2]
        } else {
            sigmao_dt.list = as.list(rep(NA,max_iter))
            for(i in 1:max_iter){
                sigmao_dt.list[[i]] = stan_data[name=="sigmao_mat",.(row,col,value)] %>%
                                      .[,iter:=i] %>%
                                      .[,sigmao:=value*stan_data[name=="sigmao_sc"]$value]%>%
                                      .[,sigmao2:=sigmao^2] %>%
                                      .[,.(iter,row,col,sigmao,sigmao2)]
            }
            sigmao_dt = rbindlist(sigmao_dt.list)
        }
        # setkeyv(sigmao_dt,cols=c("iter","row","col"))

    # storage arrays
        ll_index = array(NA,dim=c(T,I,max_iter))
        log_lik = array(0,dim=c(n_ll,max_iter))
        ll = array(0,dim=c(max_iter))

        for(j in 1:max_iter){
            # iter = 1
            tmp_sigmao_mat = dcast(sigmao_dt[iter==j,.(row,col,sigmao)],row~col,value.var="sigmao") %>%
                            .[,row:=NULL] %>%
                            as.matrix(.)
            tmp_sigmao2_mat = tmp_sigmao_mat^2

            for(i in 1:I){
                for(t in 1:T){
                    if(index[t,i]>0 & x[t,j]>0 & q[i,j]>0) {
                        mu_index = log(q[i,j]*x[t,j]) - tmp_sigmao2_mat[t,i]/2;
                        ll_index[t,i,j] = dlnorm(index[t,i], mu_index,tmp_sigmao_mat[t,i],log=TRUE);
                        # mu_index = log(q[i,j]*x[t,j]) - sigmao_dt[.(j,t,i),]$sigmao2/2;
                        # ll_index[t,i,j] = dlnorm(index[t,i], mu_index,sigmao_dt[.(j,t,i),]$sigmao,log=TRUE);
                        # if(lambdas[i]==1){
                        #     ll[j] = ll[j] + ll_index[t,i,j];
                        #     log_lik[iter,j] = ll_index[t,i,j];
                        # } else {
                        #     log_lik[iter,j] = 0
                        # }
                    } # else {
                        # ll_index[t,i,j] = 0
                        # log_lik[iter,j] = 0
                    # }
                    # iter = iter + 1
                }
            }
            rm(list=c("tmp_sigmao_mat","tmp_sigmao2_mat"))
        }

        if(nrow(hmc_samples[name=="removals"])>0){
            if(max(hmc_samples[name=="removals"]$row,na.rm=TRUE)==(max(hmc_samples$row,na.rm=TRUE)-1)){
                # calc catch likelihood
                obs_removals = dcast(stan_data[name=="obs_removals",.(row,col,value)],row~col) %>%
                    .[,row:=NULL] %>%
                    as.matrix(.)
                ll_catch = array(NA,dim=c(max(hmc_samples[name=="removals"]$row,na.rm=TRUE),max_iter))
                removals = dcast(hmc_samples[name%in%c("removals"),.(iter,row,value)],row~iter) %>%
                    .[,row:=NULL] %>%
                    as.matrix(.)
                sigmac_data = stan_data[name=="sigmac"]
                if(nrow(sigmac_data) == 1) {
                    # Single value case: replicate to length T
                    sigmac = rep(sigmac_data$value, stan_data[name=="T"]$value)
                } else {
                    # Time-varying case: extract vector in row order
                    sigmac = sigmac_data[order(row)]$value
                }

                if(nrow(hmc_samples[name=="nu_catch"])>0){
                    nu_catch = hmc_samples[name=="nu_catch"][order(as.numeric(iter))]$value
                }

                for(j in 1:max_iter){
                    
                    for(t in 1:max(hmc_samples[name=="removals"]$row,na.rm=TRUE)){
                        if(nrow(hmc_samples[name=="nu_catch"])>0){
                            mu_catch = log(removals[t,j])
                            ll_catch[t,j] = dt((log(obs_removals[t,1]) - mu_catch) / sigmac[t], df = nu_catch[j], log = TRUE) - log(sigmac[t])
                        } else {
                            mu_catch = log(removals[t,j]) - 0.5*sigmac[t]^2
                            ll_catch[t,j] = dlnorm(obs_removals[t,1],mu_catch,sigmac[t],log=TRUE)
                        }
                    }
                }
                    ll_catch_dt = as.data.table(ll_catch) %>%
                    .[,T:=1:.N] %>%
                    melt(.,id.vars="T") %>%
                    setnames(.,c("variable","value"),c("iter","value")) %>%
                    .[,iter:=as.numeric(as.character(gsub("V","",iter)))] %>%
                    .[,name:="ll__"] %>%
                    merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                    .[,lambda:=1] %>%
                    .[,I:=0] %>%
                    .[,.(run_id,T,I,iter,name,lambda,value)]                

                    ll_dt = as.data.table(ll_index) %>%
                    setnames(.,c("V1","V2","V3","value"),c("T","I","iter","value")) %>%
                    .[,T:=as.numeric(as.character(T))] %>%
                    .[,I:=as.numeric(as.character(I))] %>%
                    .[,iter:=as.numeric(as.character(iter))] %>%
                    .[,name:="ll__"] %>%
                    .[,lambda:=lambdas[I]] %>%
                    merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                    .[,.(run_id,T,I,iter,name,lambda,value)]

                    out = rbind(ll_catch_dt,ll_dt)
            } else {
                out = as.data.table(ll_index) %>%
                setnames(.,c("V1","V2","V3","value"),c("T","I","iter","value")) %>%
                .[,T:=as.numeric(as.character(T))] %>%
                .[,I:=as.numeric(as.character(I))] %>%
                .[,iter:=as.numeric(as.character(iter))] %>%
                .[,name:="ll__"] %>%
                .[,lambda:=lambdas[I]] %>%
                merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                .[,.(run_id,T,I,iter,name,lambda,value)]
            }
        } else {
                out = as.data.table(ll_index) %>%
                setnames(.,c("V1","V2","V3","value"),c("T","I","iter","value")) %>%
                .[,T:=as.numeric(as.character(T))] %>%
                .[,I:=as.numeric(as.character(I))] %>%
                .[,iter:=as.numeric(as.character(iter))] %>%
                .[,name:="ll__"] %>%
                .[,lambda:=lambdas[I]] %>%
                merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                .[,.(run_id,T,I,iter,name,lambda,value)]
        }



    return(out)
}        
