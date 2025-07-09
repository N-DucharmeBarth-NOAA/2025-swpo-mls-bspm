                ssp_calc_loofluence = function(loo_influence,stan_data){
                    loo_dt = data.table(obs_id=names(loo_influence),influence = loo_influence) %>%
                            .[,T:=as.numeric(sapply(obs_id,function(x)strsplit(x,"-")[[1]][1]))] %>%
                            .[,I:=as.numeric(sapply(obs_id,function(x)strsplit(x,"-")[[1]][2]))] %>%
                            .[order(I,T)]

                            # match with obs
                            catch_obs = stan_data[name=="obs_removals"] %>%
                                        .[,.(row,col,value)] %>%
                                        .[,T:=row] %>%
                                        .[,I:=0]
                            index_obs = stan_data[name=="index"] %>%
                                        .[,.(row,col,value)] %>%
                                        .[,T:=row] %>%
                                        .[,I:=col]
                    loo_dt_catch = na.omit(merge(loo_dt,catch_obs[,.(I,T,value)],by=c("I","T"),all=TRUE))
                    loo_dt_index = na.omit(merge(loo_dt,index_obs[,.(I,T,value)],by=c("I","T"),all=TRUE))
                    loo_dt = rbind(loo_dt_catch,loo_dt_index) %>%
                            .[,infl_cat := ifelse(influence>1,"Very high",ifelse(influence>0.7,"High","Low"))] %>%
                            .[,infl_cat := factor(infl_cat,levels=c("Low","High","Very high"))]
                    return(loo_dt)
                }
