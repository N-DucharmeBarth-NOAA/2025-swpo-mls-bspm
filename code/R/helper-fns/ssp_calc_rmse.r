    

# Nicholas Ducharme-Barth
# 2024/03/05
# Calculate index rmse

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_calc_rmse = function(hmc_samples,stan_data){
    require(data.table)
    require(magrittr)

    # define dimensions and extract data
        T = stan_data[name=="T"]$value
        I = stan_data[name=="I"]$value
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
        
        # catch
        obs_removals = stan_data[name=="obs_removals"]$value
        removals = dcast(hmc_samples[name%in%c("removals"),.(iter,row,value)],row~iter) %>%
                .[,row:=NULL] %>%
                as.matrix(.)

    # storage arrays
        nrmse_index = array(NA,dim=c(T,I,max_iter))
        nrmse_catch = array(NA,dim=c(nrow(removals),max_iter))

        for(j in 1:max_iter){
            # iter = 1
            for(i in 1:I){
                for(t in 1:T){
                    if(index[t,i]>0 & x[t,j]>0 & q[i,j]>0) {
                        nrmse_index[t,i,j] = (index[t,i] - q[i,j]*x[t,j])^2
                      
                    } 
                }
            }
        }

        for(j in 1:max_iter){
            for(t in 1:nrow(removals)){
                nrmse_catch[t,j] = (obs_removals[t] - removals[t,j])^2
            }
        }
        colnames(nrmse_catch) = 1:max_iter

        se_index = as.data.table(index) %>%
            melt(.,measure.vars=colnames(index)) %>%
            setnames(.,"variable","I") %>%
            .[,I:=as.numeric(as.character(I))] %>%
            .[value > -999] %>%
            .[,.(mu=mean(value)),by=I]

        index_rmse_dt = as.data.table(nrmse_index) %>%
                setnames(.,c("V1","V2","V3","value"),c("T","I","iter","value")) %>%
                .[,T:=as.numeric(as.character(T))] %>%
                .[,I:=as.numeric(as.character(I))] %>%
                .[,iter:=as.numeric(as.character(iter))] %>%
                .[,name:="se"] %>%
                merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                merge(.,se_index,by="I") %>%
                .[,.(rmse=sqrt(mean(value))/mu),by=.(run_id,I,iter)] %>%
                .[,lambda:=lambdas[I]] %>%
                .[,type:="index"] %>%
                .[,.(run_id,type,I,lambda,iter,rmse)] %>%
                unique(.)
        
        catch_rmse_dt = as.data.table(nrmse_catch) %>%
                        .[, T := 1:.N] %>%
                        melt(., id.vars = "T") %>%
                setnames(.,c("variable"),c("iter")) %>%
                .[,iter:=as.numeric(as.character(iter))] %>%
                .[,name:="se"] %>%
                merge(.,unique(hmc_samples[,.(run_id,iter)]),by="iter") %>%
                .[,mu:=mean(obs_removals[1:nrow(removals)])] %>%
                .[,.(rmse=sqrt(mean(value))/mu),by=.(run_id,iter)] %>%
                .[,type:="catch"] %>%
                .[,I:=NA] %>%
                .[,lambda:=NA] %>%
                .[,.(run_id,type,I,lambda,iter,rmse)] %>%
                unique(.)
        
        rmse_dt = rbind(index_rmse_dt,catch_rmse_dt)

    return(rmse_dt)
}        

