# Nicholas Ducharme-Barth
# 2024/03/04
# Calculate stochastic surplus production model derived quantities from hmc samples
# Modified to handle F estimation scenarios

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

ssp_derived_quants = function(hmc_samples,stan_data,output="percentiles",percentile=0.5){
        tmp1 = hmc_samples[,.(run_id,iter,name,value)] %>%
              .[name%in%c("m","logK","dmsy")] %>%
              dcast(.,run_id+iter~name) %>%
              .[,msy:=m*exp(logK)] %>%
              .[,Dmsy:=dmsy] %>%
              .[,Pmsy:=dmsy * exp(logK)] %>%
              .[,Umsy:=m/dmsy] %>%
              .[,Fmsy:=-log(-Umsy+1)] %>%
              .[,.(run_id,iter,msy,Dmsy,Pmsy,Umsy,Fmsy)] 
        
        if(!("removals" %in% unique(hmc_samples$name))){
            # Scenario 1: Neither removals nor F estimated - use observed catch
            tmp2 = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c("x[1]","logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,initial_population:=x*exp(logK)] %>%
                  .[,initial_depletion:=x] %>%
                  .[,initial_U:=ifelse(stan_data[name=="removals"&row==1]$value/initial_population>1,0.9999,stan_data[name=="removals"&row==1]$value/initial_population)] %>%
                  .[,initial_F:=-log(-initial_U+1)] %>%
                  .[,.(run_id,iter,initial_population,initial_depletion,initial_U,initial_F)]

            tmp3 = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c(paste0("x[",stan_data[name=="T"]$value,"]"),"logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,latest_population:=x*exp(logK)] %>%
                  .[,latest_depletion:=x] %>%
                  .[,latest_U:=ifelse(stan_data[name=="removals"&row==stan_data[name=="T"]$value]$value/latest_population>1,0.9999,stan_data[name=="removals"&row==stan_data[name=="T"]$value]$value/latest_population)] %>%
                  .[,latest_F:=-log(-latest_U+1)] %>%
                  .[,.(run_id,iter,latest_population,latest_depletion,latest_U,latest_F)]

            tmp_dt.list = as.list(rep(NA,4))
            for(i in 1:4){
                  tmp_dt.list[[i]] = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c(paste0("x[",stan_data[name=="T"]$value-i,"]"),"logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,population_Tm:=x*exp(logK)] %>%
                  .[,depletion_Tm:=x] %>%
                  .[,U_Tm:=ifelse(stan_data[name=="removals"&row==stan_data[name=="T"]$value-i]$value/population_Tm>1,0.9999,stan_data[name=="removals"&row==stan_data[name=="T"]$value-i]$value/population_Tm)] %>%
                  .[,F_Tm:=-log(-U_Tm+1)] %>%
                  .[,.(run_id,iter,population_Tm,depletion_Tm,U_Tm,F_Tm)] %>%
                  setnames(.,c("population_Tm","depletion_Tm","U_Tm","F_Tm"),paste0(c("population_T","depletion_T","U_T","F_T"),i))
            }
            recent_dt = merge(tmp_dt.list[[1]],tmp_dt.list[[2]],by=c("run_id","iter")) %>%
                        merge(tmp_dt.list[[3]],by=c("run_id","iter")) %>%
                        merge(tmp_dt.list[[4]],by=c("run_id","iter")) %>%
                        merge(.,tmp3[,.(run_id,iter,latest_population,latest_depletion)]) %>%
                        .[,recent_U:=(U_T1+U_T2+U_T3+U_T4)/4] %>%
                        .[,recent_F:=(F_T1+F_T2+F_T3+F_T4)/4] %>%
                        .[,recent_P:=(population_T3+population_T2+population_T1+latest_population)/4] %>%
                        .[,recent_D:=(depletion_T3+depletion_T2+depletion_T1+latest_depletion)/4] %>%
                        .[,.(run_id,iter,recent_U,recent_F,recent_P,recent_D)]
        } else {
            # Scenario 2: Both removals AND F estimated - use estimated F values
            initial_F = hmc_samples[name=="F"&row==1,.(run_id,iter,value)] %>%
                       setnames(.,"value","initial_F")
            latest_F = hmc_samples[name=="F"&row==(stan_data[name=="T"]$value-1),.(run_id,iter,value)] %>%
                      setnames(.,"value","latest_F")
            
            tmp2 = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c("x[1]","logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,initial_population:=x*exp(logK)] %>%
                  .[,initial_depletion:=x] %>%
                  merge(.,initial_F,by=c("run_id","iter")) %>%
                  .[,initial_U:=1-exp(-initial_F)] %>%
                  .[,.(run_id,iter,initial_population,initial_depletion,initial_U,initial_F)]

            # note that latest F only corresponds to T-1
            tmp3 = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c(paste0("x[",stan_data[name=="T"]$value,"]"),"logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,latest_population:=x*exp(logK)] %>%
                  .[,latest_depletion:=x] %>%
                  merge(.,latest_F,by=c("run_id","iter")) %>%
                  .[,latest_U:=1-exp(-latest_F)] %>%
                  .[,.(run_id,iter,latest_population,latest_depletion,latest_U,latest_F)]

            # For recent calculations, extract F values for the last 5 years excluding last year
            tmp_F.list = as.list(rep(NA,4))
            for(i in 1:4){
                tmp_F.list[[i]] = hmc_samples[name=="F"&row==(stan_data[name=="T"]$value-i),.(run_id,iter,value)] %>%
                                 setnames(.,"value",paste0("F_T",i))
            }
            
            # extract population for last 4 years including the last year
            tmp_dt.list = as.list(rep(NA,4))
            for(i in 1:4){
                  tmp_dt.list[[i]] = hmc_samples[,.(run_id,iter,variable,name,value)] %>%
                  .[variable%in%c(paste0("x[",stan_data[name=="T"]$value-i+1,"]"),"logK")] %>%
                  .[,.(run_id,iter,name,value)] %>%
                  dcast(.,run_id+iter~name) %>%
                  .[,population_Tm:=x*exp(logK)] %>%
                  .[,depletion_Tm:=x] %>%
                  merge(.,tmp_F.list[[i]],by=c("run_id","iter")) %>%
                  .[,F_Tm:=get(paste0("F_T",i))] %>%
                  .[,U_Tm:=1-exp(-F_Tm)] %>%
                  .[,.(run_id,iter,population_Tm,depletion_Tm,U_Tm,F_Tm)] %>%
                  setnames(.,c("population_Tm","depletion_Tm","U_Tm","F_Tm"),paste0(c("population_T","depletion_T","U_T","F_T"),i))
            }
            
            recent_dt = merge(tmp_dt.list[[1]],tmp_dt.list[[2]],by=c("run_id","iter")) %>%
                        merge(tmp_dt.list[[3]],by=c("run_id","iter")) %>%
                        merge(tmp_dt.list[[4]],by=c("run_id","iter")) %>%
                        merge(.,tmp3[,.(run_id,iter,latest_population,latest_depletion)]) %>%
                        .[,recent_U:=(U_T1+U_T2+U_T3+U_T4)/4] %>%
                        .[,recent_F:=(F_T1+F_T2+F_T3+F_T4)/4] %>%
                        .[,recent_P:=(population_T3+population_T2+population_T1+latest_population)/4] %>%
                        .[,recent_D:=(depletion_T3+depletion_T2+depletion_T1+latest_depletion)/4] %>%
                        .[,.(run_id,iter,recent_U,recent_F,recent_P,recent_D)]
        }

        tmp4 = merge(tmp1,tmp3,by=c("run_id","iter")) %>%
               .[,latest_P_Pmsy := latest_population/Pmsy] %>%
               .[,latest_D_Dmsy := latest_depletion/Dmsy] %>%
               .[,latest_U_Umsy := latest_U/Umsy] %>%
               .[,latest_F_Fmsy := latest_F/Fmsy] %>%
              .[,.(run_id,iter,latest_P_Pmsy,latest_D_Dmsy,latest_U_Umsy,latest_F_Fmsy)]
        
        tmp5 = merge(tmp1,recent_dt,by=c("run_id","iter")) %>%
               .[,recent_P_Pmsy := recent_P/Pmsy] %>%
               .[,recent_D_Dmsy := recent_D/Dmsy] %>%
               .[,recent_U_Umsy := recent_U/Umsy] %>%
               .[,recent_F_Fmsy := recent_F/Fmsy] %>%
              .[,.(run_id,iter,recent_P_Pmsy,recent_D_Dmsy,recent_U_Umsy,recent_F_Fmsy)]

        # merge...
        tmp = merge(tmp1,tmp2,by=c("run_id","iter")) %>%
              merge(.,tmp3,by=c("run_id","iter")) %>%
              merge(.,tmp4,by=c("run_id","iter")) %>%
              merge(.,tmp5,by=c("run_id","iter"))

        if(output=="percentiles"){
            tmp = tmp %>%
                    melt(.,id.vars = c("run_id","iter")) %>%
                    .[,.(value=quantile(value,probs=percentile,na.rm=TRUE)),by=.(run_id,variable)]
            tmp = tmp %>%
                    .[,variable:=factor(variable,levels=unique(tmp$variable),labels=paste0(unique(tmp$variable),"-p",percentile))] %>%
                    dcast(.,run_id~variable)          
        }
    return(tmp)
}
