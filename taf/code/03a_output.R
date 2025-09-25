

# NOAA PIFSC
# 2025/09/24
# Extract results from the diagnostic model
# Model 0100: DWFN index; baseline priors

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
    i = 1
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
        derived_quants_ts = ssp_derived_quants_ts(ssp_summary=fread(file.path("model",run_label,"fit_summary.csv")),
                            samples_dt=fread(file.path("model",run_label,"hmc_samples.csv")),
                            settings=fread(file.path("model",run_label,"settings.csv")),
                            stan_data=fread(file.path("model",run_label,"stan_data.csv")),
                            sub_sample_prop=1) %>%
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
    
                                 