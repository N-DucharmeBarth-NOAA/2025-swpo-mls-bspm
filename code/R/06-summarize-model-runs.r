

# Nicholas Ducharme-Barth
# 2025/06/10
# summarize data from model runs

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(rstan)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)    

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for model outputs
    model_stem = file.path(proj_dir,"data","output","model_runs")

#_____________________________________________________________________________________________________________________________
# make summary file to seed the shiny app
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("fit_summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("fit_summary.csv","",all_dirs,fixed=TRUE)

    summary_df.list = lapply(all_dirs,function(x)as.data.frame(fread(file.path(model_stem,x,"fit_summary.csv"))))

    # # round numeric
    # for(i in 1:length(summary_df.list)){
    #     for(j in 1:ncol(summary_df.list[[i]])){
    #         if(class(summary_df.list[[i]][,j]) == "numeric"){
    #             # summary_df[,j] = round(summary_df[,j],digits=3)
    #             summary_df.list[[i]][,j] = trimws(format(round(summary_df.list[[i]][,j],digits=3),big.mark=",",scientific=FALSE))
    #         } else if(class(summary_df.list[[i]][,j]) == "integer64"){
    #             summary_df.list[[i]][,j] = trimws(format(round(as.numeric(summary_df.list[[i]][,j]),digits=3),big.mark=",",scientific=FALSE))
    #         }
    #     }
    #     summary_df.list[[i]] = as.data.table(summary_df.list[[i]])
    # }

    summary_dt = rbindlist(summary_df.list,fill=TRUE) %>% .[run_retro == 0]
    fwrite(summary_dt,file=file.path("data","output","summary.csv"))

#_____________________________________________________________________________________________________________________________
# calc retrospectives & mase
    summary_retro_vec = rep(NA,nrow(summary_dt))
    summary_hcxval_vec = rep(NA,nrow(summary_dt))
    summary_retro_pmsy_cover_vec = rep(NA,nrow(summary_dt))
    summary_retro_fmsy_cover_vec = rep(NA,nrow(summary_dt))
    retro_peels_length_vec = rep(NA,nrow(summary_dt))
    
    for(i in 1:nrow(summary_dt)){
        tmp_stem = paste0(strsplit(summary_dt$run_label[i],"_")[[1]][1])
        retro_peels = unname(sapply(all_dirs[grep(tmp_stem,all_dirs,fixed=TRUE)],function(x)as.numeric(gsub("/","",strsplit(x,"_")[[1]][2],fixed=TRUE))))
        retro_peels_length_vec[i] = length_retro_peels = length(retro_peels)-1
        if(length_retro_peels>0){
            retro_dt.list = as.list(rep(NA,length(retro_peels)))
            hcxval_dt.list = as.list(rep(NA,length(retro_peels)))
            med_pmsy = rep(NA,length(retro_peels))
            med_fmsy = rep(NA,length(retro_peels))
            cover_range_pmsy = matrix(NA,nrow=length(retro_peels),ncol=2)
            cover_range_fmsy = matrix(NA,nrow=length(retro_peels),ncol=2)

                    for(j in 1:length(retro_peels)){
                        # retro calcs
                            retro_dt.list[[j]] = ssp_derived_quants_ts(ssp_summary=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"fit_summary.csv")),
                                                                samples_dt = fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"hmc_samples.csv")),
                                                                stan_data=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"stan_data.csv")),
                                                                settings=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"settings.csv")),
                                                                sub_sample_prop=1)
                            if(retro_peels[j] == 0){
                                 pmsy_l95 = retro_dt.list[[j]] %>%
                                    .[name=="P_Pmsy"] %>%
                                    .[,.(l95=quantile(value,probs=0.025)),by=.(run_id,name,row)]
                                 pmsy_u95 = retro_dt.list[[j]] %>%
                                    .[name=="P_Pmsy"] %>%
                                    .[,.(u95=quantile(value,probs=0.975)),by=.(run_id,name,row)]
                                cover_range_pmsy[,1] = rev(tail(pmsy_l95$l95,n=max(retro_peels)+1))
                                cover_range_pmsy[,2] = rev(tail(pmsy_u95$u95,n=max(retro_peels)+1))
                                 fmsy_l95 = retro_dt.list[[j]] %>%
                                    .[name=="F_Fmsy"] %>%
                                    .[,.(l95=quantile(value,probs=0.025)),by=.(run_id,name,row)]
                                 fmsy_u95 = retro_dt.list[[j]] %>%
                                    .[name=="F_Fmsy"] %>%
                                    .[,.(u95=quantile(value,probs=0.975)),by=.(run_id,name,row)]
                                cover_range_fmsy[,1] = rev(tail(fmsy_l95$l95,n=max(retro_peels)+1))
                                cover_range_fmsy[,2] = rev(tail(fmsy_u95$u95,n=max(retro_peels)+1))
                            }

                            med = retro_dt.list[[j]] %>%
                                    .[name%in%c("P_Pmsy","F_Fmsy")] %>%
                                    .[,.(median=quantile(value,probs=0.5)),by=.(run_id,name,row)] %>%
                                    .[row==max(row)-retro_peels[j]]
                            med_pmsy[j] = med[name=="P_Pmsy"]$median
                            med_fmsy[j] = med[name=="F_Fmsy"]$median
                            
                            retro_dt.list[[j]] = retro_dt.list[[j]] %>%
                                                        .[name == "D"] %>%
                                                        .[,.(D=median(value)),by=.(run_id,name,row)] %>%
                                                        .[,retro:=retro_peels[j]] %>%
                                                        .[row%in%min(retro_dt.list[[j]]$row):(max(retro_dt.list[[j]]$row)-retro_peels[j])]
                        # hcxval
                            hcxval_dt.list[[j]] = ssp_extract_cpue_fit(ssp_summary=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"fit_summary.csv")),
                                                                samples_dt = fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"hmc_samples.csv")),
                                                                stan_data=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"stan_data.csv")),
                                                                settings=fread(file.path(model_stem,paste0(tmp_stem,"_",retro_peels[j]),"settings.csv")),
                                                                sub_sample_prop=1,
                                                                    active="TRUE",
                                                                    calc_std = "FALSE") %>%
                                                                .[metric %in% c("obs_cpue","pred_cpue")] %>%
                                                                .[,.(value=median(value)),by=.(run_id,metric,row)] %>%
                                                                .[,retro:=retro_peels[j]]

                    }

                    # retro calcs
                        retro_dt=rbindlist(retro_dt.list)

                        retro_vec = rep(NA,length(retro_peels)-1)
                        for(j in 1:length(retro_vec)){
                            tmp_peel = retro_dt[retro==retro_peels[j+1]] %>%
                                    .[row==max(row)]
                            tmp_base = retro_dt[retro==0&row==tmp_peel$row]
                            retro_vec[j] = (tmp_peel$D - tmp_base$D)/tmp_base$D
                        }
                        summary_retro_vec[i] = mean(retro_vec)
                    # hcxval calcs
                        hcxval_dt=rbindlist(hcxval_dt.list)
                        hcxval_maxrow = max(na.omit(hcxval_dt[metric=="obs_cpue"&retro==0])$row)
                        hcxval_maxrow_overall = max(hcxval_dt[metric=="obs_cpue"&retro==0]$row)
                        mase_top = rep(NA,length(retro_peels)-1)
                        mase_bottom = rep(NA,length(retro_peels)-1)
                        # for(j in 1:length(mase_top)){
                        #     tmp_base = hcxval_dt[metric=="obs_cpue"&retro==0&row==hcxval_maxrow-j+1]
                        #     tmp_pred = hcxval_dt[metric=="pred_cpue"&retro==retro_peels[j+1]&row==hcxval_maxrow-j+1]
                        #     tmp_naive = hcxval_dt[metric=="obs_cpue"&retro==retro_peels[j+1]&row==hcxval_maxrow-j]
                        #     mase_top[j] = tmp_base$value - tmp_pred$value
                        #     mase_bottom[j] = tmp_base$value - tmp_naive$value
                        # }
                        # summary_hcxval_vec[i] = mean(abs(mase_top))/mean(abs(mase_bottom)) 

                        for(j_prime in (1+(hcxval_maxrow_overall-hcxval_maxrow)):(length(retro_peels)-1)){
                            j = j_prime-(hcxval_maxrow_overall-hcxval_maxrow)
                            tmp_base = hcxval_dt[metric=="obs_cpue"&retro==0&row==hcxval_maxrow-j+1]
                            tmp_pred = hcxval_dt[metric=="pred_cpue"&retro==retro_peels[j_prime+1]&row==hcxval_maxrow-j+1]
                            tmp_naive = hcxval_dt[metric=="obs_cpue"&retro==retro_peels[j_prime+1]&row==hcxval_maxrow-j]
                            mase_top[j] = median(tmp_base$value) - median(tmp_pred$value)
                            mase_bottom[j] = median(tmp_base$value) - median(tmp_naive$value)
                        }
                        summary_hcxval_vec[i] = mean(abs(mase_top),na.rm=TRUE)/mean(abs(mase_bottom),na.rm=TRUE)
                    
                    # look at coverage for P/Pmsy and F/Fmsy  
                    cover_vec_pmsy = rep(NA,length(med_pmsy))
                    cover_vec_fmsy = rep(NA,length(med_fmsy))
                    for(j in 1:length(med_pmsy)){
                        if(med_pmsy[j]>cover_range_pmsy[j,1]&med_pmsy[j]<cover_range_pmsy[j,2]){
                            cover_vec_pmsy[j] = 1
                        } else {
                            cover_vec_pmsy[j] = 0
                        }
                        if(med_fmsy[j]>cover_range_fmsy[j,1]&med_fmsy[j]<cover_range_fmsy[j,2]){
                            cover_vec_fmsy[j] = 1
                        } else {
                            cover_vec_fmsy[j] = 0
                        }
                    }
                    summary_retro_pmsy_cover_vec[i] = mean(cover_vec_pmsy[-1],na.rm=TRUE)
                    summary_retro_fmsy_cover_vec[i] = mean(cover_vec_fmsy[-1],na.rm=TRUE)

                    # save relevant outputs for plotting later (coverage, retro, hcxval)
                    fwrite(retro_dt,file=file.path(model_stem,paste0(tmp_stem,"_0"),"retro_dt.csv"))
                    fwrite(hcxval_dt,file=file.path(model_stem,paste0(tmp_stem,"_0"),"hcxval_dt.csv"))
                    write.csv(med_pmsy,file=file.path(model_stem,paste0(tmp_stem,"_0"),"med_pmsy.csv"))
                    write.csv(med_fmsy,file=file.path(model_stem,paste0(tmp_stem,"_0"),"med_fmsy.csv"))
                    write.csv(cover_range_pmsy,file=file.path(model_stem,paste0(tmp_stem,"_0"),"cover_range_pmsy.csv"))
                    write.csv(cover_range_fmsy,file=file.path(model_stem,paste0(tmp_stem,"_0"),"cover_range_fmsy.csv"))
        } 
    }

    retro_hcxval_dt = data.table(run_id = summary_dt$run_id,retro=summary_retro_vec,hcxval=summary_hcxval_vec,coverage_pmsy=summary_retro_pmsy_cover_vec,coverage_fmsy=summary_retro_fmsy_cover_vec,n_retro=retro_peels_length_vec)
    fwrite(retro_hcxval_dt,file=file.path("data","output","retro_hcxval.csv"))
    fwrite(merge(summary_dt,retro_hcxval_dt,by="run_id"),file=file.path("data","output","summary_retro0.csv"))
