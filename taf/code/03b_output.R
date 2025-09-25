

# NOAA PIFSC
# 2025/09/24
# Extract results from the model ensemble

# Copyright (c) 2025 NOAA PIFSC
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(TAF)
    library(data.table)
    library(magrittr)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    source_dir_stem = file.path("boot","software")
    source_files = list.files(source_dir_stem)
    source_files = source_files[grep(".r",source_files,fixed=TRUE)]
    sapply(file.path(source_dir_stem,source_files),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# bring in data
    ssp_scenarios = fread(file.path("data","ssp_scenarios.csv"))
    
#________________________________________________________________________________________________________________________________________________________________________________________________________
# identify target model

    target_model_dirs = apply(ssp_scenarios,1,function(x)paste0("0",x[[1]],"-",x[[2]],"-",x[[3]],"_0"))

    all_data = lapply(file.path("model",target_model_dirs), load_model_data)
    posterior_dt_list = as.list(rep(NA,nrow(ssp_scenarios)))

    for(i in 1:nrow(ssp_scenarios)){
            run_label_stem = paste0(c(ssp_scenarios$cpue[i],ssp_scenarios$shape[i]),collapse="-")    
            run_number = paste0("0",ssp_scenarios$run_number[i])
            run_label = paste0(run_number,"-",run_label_stem,"_0")

            cpue_idx = switch(as.character(ssp_scenarios$cpue[i]),
                                    "dwfn" = 1,
                                    "au" = 2, 
                                    "nz" = 3,
                                    "obs" = 4,
                                    "obsNoPF" = 5,
                                    "obsPFonly" = 6,
                                    1)

        #________________________________________________________________________________________________________________________________________________________________________________________________________
        # make directory for model outputs
            mkdir("output")
            mkdir(file.path("output",run_label))

        #________________________________________________________________________________________________________________________________________________________________________________________________________
        # get relevant fit metrics  
            summary_dt = fread(file.path("model",run_label,"fit_summary.csv"))

            base_diags_dt = summary_dt[,c("run_label","divergent","max_rhat","min_neff")]
            base_diags_dt$divergent = base_diags_dt$divergent*1000
            base_diags_dt$min_neff = base_diags_dt$min_neff*1000
            base_diags_dt$index_rmse = as.data.frame(summary_dt)[,c("index_rmse_1","index_rmse_2","index_rmse_3","index_rmse_4","index_rmse_5","index_rmse_6")][1,cpue_idx]
            base_diags_dt$index_loglik = as.data.frame(summary_dt)[,c("index_ll_1","index_ll_2","index_ll_3","index_ll_4","index_ll_5","index_ll_6")][1,cpue_idx]
            base_diags_dt$index = ssp_scenarios$cpue[i]
            base_diags_dt$shape_prior = ssp_scenarios$shape[i]
            base_diags_dt = base_diags_dt[,.(run_label,index,shape_prior,divergent,max_rhat,min_neff,index_rmse,index_loglik)]

        #________________________________________________________________________________________________________________________________________________________________________________________________________
        # extract time series (medians): removals, N, D, F, D_Dmsy, F_Fmsy
                first_year = 1952
                posterior_dt_list[[i]] = ssp_derived_quants_ts(ssp_summary=fread(file.path("model",run_label,"fit_summary.csv")),
                            samples_dt=fread(file.path("model",run_label,"hmc_samples.csv")),
                            settings=fread(file.path("model",run_label,"settings.csv")),
                            stan_data=fread(file.path("model",run_label,"stan_data.csv")),
                            sub_sample_prop=1)
                derived_quants_ts = posterior_dt_list[[i]] %>%
                        .[name%in%c("removals","P","D","F","D_Dmsy","F_Fmsy")] %>%
                        .[,name:=factor(name,levels=c("removals","P","D","F","D_Dmsy","F_Fmsy"),labels=c("Removals","Numbers","Depletion","Fishing mortality","D_Dmsy","F_Fmsy"))] %>%
                        .[,row:=(row-1)+first_year] %>%
                        .[,.(value=median(value,na.rm=TRUE)),by=.(name,row)] %>%
                        .[,run_label:=run_label] %>%
                        setnames(.,"row","Year") %>%
                        dcast(.,run_label+Year~name)

        #________________________________________________________________________________________________________________________________________________________________________________________________________
        # extract observed cpue, predicted cpue and residual
                fit_dt = ssp_extract_cpue_fit(ssp_summary=fread(file.path("model",run_label,"fit_summary.csv")),
                            samples_dt=fread(file.path("model",run_label,"hmc_samples.csv")),
                            stan_data=fread(file.path("model",run_label,"stan_data.csv")),
                            settings=fread(file.path("model",run_label,"settings.csv")),
                            sub_sample_prop=1,
                            active="TRUE",
                            calc_std = "TRUE") %>%
                        .[metric %in% c("obs_cpue","pred_cpue","std_residual")] %>%
                        .[,metric:=factor(metric,levels=c("obs_cpue","pred_cpue","std_residual"),labels=c("Obs_CPUE","Pred_CPUE","Std_Residual"))] %>%
                        .[,row:=(row-1)+first_year] %>%
                        .[,.(value=median(value)),by=.(metric,row)] %>%
                        .[,run_label:=run_label] %>%
                        setnames(.,"row","Year") %>%
                        dcast(.,run_label+Year~metric)

        #________________________________________________________________________________________________________________________________________________________________________________________________________
        # write-out data files
            fwrite(base_diags_dt,file=file.path("output",run_label,"diagnostics.csv"))
            fwrite(derived_quants_ts,file=file.path("output",run_label,"derived_quantities.csv"))  
            fwrite(fit_dt,file=file.path("output",run_label,"fit.csv"))
            
                                        
    }

    tmp_summary = rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)

    yo_dt = tmp_summary[,.(run_label)] %>%
            unique(.) %>%
            .[,year_one:=extract_model_start_year(run_label)]

    posterior_dt = rbindlist(posterior_dt_list) %>% merge(., tmp_summary[, .(run_id, run_label)]) %>%
                merge(.,yo_dt,by="run_label") %>%  
                .[row >= 1, year := year_one + (row - 1)] %>%
                .[row < 1, year := year_one + (row - 1)] %>%
                .[name %in% c("Process error", "Process error (raw)") & row > 0, year := year + 1]


    mgmt_summary = posterior_dt[name %in% c("D","P","F","D_Dmsy","F_Fmsy","removals"), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name,year)] %>%
                    .[order(year, name)]
    
    fwrite(mgmt_summary,file=file.path("output","ens_mgmt_summary.csv"))
