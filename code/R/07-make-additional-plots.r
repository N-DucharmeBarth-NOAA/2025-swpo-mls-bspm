# Nicholas Ducharme-Barth
# 2025/06/10
# Make additional plots
# Retrospective analysis and hindcast cross-validation for BSPM

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(GGally)
    library(viridis)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
    dir_plot_fns = file.path(proj_dir,"code","R","plot-fns")
    model_stem = file.path(proj_dir,"data","output","model_runs")
    diag_model = "0100-dwfn-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0"

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)
    sapply(file.path(dir_plot_fns,(list.files(dir_plot_fns))),source)


    set_global_config(
        index_names = c("DWFN","AU","NZ","Obs (all)","Obs (NC,FJ&TO)","Obs (PF)"), 
        model_stem = file.path("data","output","model_runs"),
        height_per_panel = 350
    )
#________________________________________________________________________________________________________________________________________________________________________________________________________
# retrospective analysis - diagnostic
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("fit_summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("fit_summary.csv","",all_dirs,fixed=TRUE)

    summary_dt = fread(file.path(proj_dir,"data","output","summary.csv")) %>%
                  .[!is.na(run_label)] %>%
                  .[run_label==diag_model]

    model_retro_dt.list = list()

    for(i in 1:nrow(summary_dt)){
        tmp_stem = paste0(strsplit(summary_dt$run_label[i],"_")[[1]][1])
        retro_peels = unname(sapply(all_dirs[grep(tmp_stem,all_dirs,fixed=TRUE)],function(x)as.numeric(gsub("/","",strsplit(x,"_")[[1]][2],fixed=TRUE))))
        length_retro_peels = length(retro_peels)-1
        
        if(length_retro_peels>0){
            retro_dt.list = list()

            for(j in 1:length(retro_peels)){
                retro_dt.list[[j]] = ssp_derived_quants_ts(ssp_summary=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/fit_summary.csv")),
                                                    samples_dt = fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/hmc_samples.csv")),
                                                    stan_data=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/stan_data.csv")),
                                                    settings=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/settings.csv")),
                                                    sub_sample_prop=1)  %>%
                                                    .[,retro:=retro_peels[j]] %>%
                                                    .[row%in%1:(max(row)-retro_peels[j])]
            }

            year_one=1952
            retro_base = rbindlist(retro_dt.list) %>%
                .[row>=1,year:=year_one+(row-1)] %>%
                .[row<1,year:=year_one+(row-1)] %>%
                .[name%in%c("D","F","D_Dmsy","F_Fmsy")] %>%
                .[,name:=factor(name,levels=c("D","D_Dmsy","F","F_Fmsy"),labels=c("D",expression("D"/"D"["MSY"]),"F",expression("F"/"F"["MSY"])))] %>%
                .[,.(med=median(value,na.rm=TRUE),avg=mean(value,na.rm=TRUE),lp=quantile(value,probs=0.025,na.rm=TRUE),up=quantile(value,probs=0.975,na.rm=TRUE),lp50=quantile(value,probs=0.25,na.rm=TRUE),up50=quantile(value,probs=0.75,na.rm=TRUE)),by=.(retro,name,year)] %>%
                .[,retro:=factor(retro,levels=0:5,labels=2022:2017)] %>%
                .[,model_name:=summary_dt$run_number[i]]

            model_retro_dt.list[[i]] = retro_base
        }
    }
       
    retro_combined = rbindlist(model_retro_dt.list) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]
    
    p = retro_combined %>% ggplot() +
        ylab("Metric") +
        xlab("Year") +
        facet_grid(name~short_plot_names,scales="free_y",labeller = labeller(name = label_parsed, model_name = label_value))
    p = p + geom_ribbon(data=retro_combined[retro=="2022"],aes(x=year,ymin=lp,ymax=up),fill="gray80",alpha=0.6,linewidth=1.15)
    p = p + geom_ribbon(data=retro_combined[retro=="2022"],aes(x=year,ymin=lp50,ymax=up50),fill="gray40",alpha=0.6,linewidth=1.15)
    p = p + geom_path(data=retro_combined[retro=="2022"],aes(x=year,y=med),color="black",linewidth=1.15)
    p = p + geom_path(data=retro_combined[retro!="2022"],aes(x=year,y=med,color=retro),linewidth=0.8)
    p = p + geom_point(data=retro_combined[retro!="2022"&year==as.numeric(as.character(retro))],aes(x=year,y=med,fill=retro),shape=21,size=3)
    p = p + geom_hline(yintercept=0) +
                viridis::scale_color_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE,labels = scales::parse_format()) +
                viridis::scale_fill_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE,labels = scales::parse_format()) +
                theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                        strip.background =element_rect(fill="white"),
                        legend.key = element_rect(fill = "white"))       
    ggsave(filename="retro_analysis.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 12, units = c("in"),
            dpi = 300, limitsize = TRUE)   

#________________________________________________________________________________________________________________________________________________________________________________________________________
# hindcast cross-validation  - diagnostic
    model_hcxval_dt.list = list()

    for(i in 1:nrow(summary_dt)){
        tmp_stem = paste0(strsplit(summary_dt$run_label[i],"_")[[1]][1])
        retro_peels = unname(sapply(all_dirs[grep(tmp_stem,all_dirs,fixed=TRUE)],function(x)as.numeric(gsub("/","",strsplit(x,"_")[[1]][2],fixed=TRUE))))
        length_retro_peels = length(retro_peels)-1
        
        if(length_retro_peels>0){
            hcxval_dt.list = list()

            for(j in 1:length(retro_peels)){
                hcxval_dt.list[[j]] = ssp_extract_cpue_fit(ssp_summary=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/fit_summary.csv")),
                                                    samples_dt = fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/hmc_samples.csv")),
                                                    stan_data=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/stan_data.csv")),
                                                    settings=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/settings.csv")),
                                                    sub_sample_prop=1,
                                                        active="TRUE",
                                                        calc_std = "FALSE") %>%
                                                        .[metric %in% c("sigmao","obs_cpue","pred_cpue")] %>%
                                                        .[,retro:=retro_peels[j]]
            }

            year_one=1952
            hcxval_base = rbindlist(hcxval_dt.list) %>%
                .[row>=1,year:=year_one+(row-1)] %>%
                .[row<1,year:=year_one+(row-1)]

            hcxval_maxrow = max(na.omit(hcxval_base[metric=="obs_cpue"&retro==0])$row)
            hcxval_maxrow_overall = max(hcxval_base[metric=="obs_cpue"&retro==0]$row)

            hcxval_thin_dt.list = list()
            for(j_prime in (1+(hcxval_maxrow_overall-hcxval_maxrow)):(length(retro_peels)-1)){
                j = j_prime-(hcxval_maxrow_overall-hcxval_maxrow)
                hcxval_thin_dt.list[[j]] = rbind(hcxval_base[metric%in%c("obs_cpue","sigmao")&retro==0&row%in%1:(hcxval_maxrow-j+1),.(value=median(value)),by=.(run_id,metric,row,index,year,retro)],
                                                 hcxval_base[metric=="pred_cpue"&retro==retro_peels[j_prime+1]&row%in%1:(hcxval_maxrow-j+1),.(value=median(value)),by=.(run_id,metric,row,index,year,retro)])
            }
            hcxval_thin_dt = rbindlist(hcxval_thin_dt.list) %>%
                .[,retro:=factor(retro,levels=0:5,labels=2022:2017)] %>%
                .[,model_name:=summary_dt$run_number[i]]

            model_hcxval_dt.list[[i]] = hcxval_thin_dt
        }
    }
       
    hcxval_combined = rbindlist(model_hcxval_dt.list) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    obs_se_dt = hcxval_combined[metric=="sigmao"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        .[,se:=median(value),by=.(model_name,index,row,year)] %>%
        .[,se:=round(se,digits=3)] %>%
        .[,.(model_name,index,row,year,se)] %>%
        unique(.) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    obs_cpue_dt = hcxval_combined[metric=="obs_cpue"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        merge(.,obs_se_dt,by=c("model_name","index","row","year")) %>%
        .[,obs:=round(value,digits=3)] %>%
        .[,upper:=qlnorm(0.975,meanlog=log(obs),sdlog=se)] %>%
        .[,lower:=qlnorm(0.025,meanlog=log(obs),sdlog=se)] %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    p = hcxval_combined %>% ggplot() +
        ylab("Index") +
        xlab("Year") +
        geom_hline(yintercept=1,linetype="dashed") +
        facet_wrap(~short_plot_names,ncol=1)
    p = p + geom_segment(data=obs_cpue_dt,aes(x=year,xend=year,y=lower,yend=upper),linewidth=0.5)
    p = p + geom_point(data=obs_cpue_dt,aes(x=year,y=obs),color="black",fill="white",shape=21,size=3)
    p = p + geom_path(data=hcxval_combined[retro!="2022"&metric=="pred_cpue"],aes(x=year,y=value,color=retro),linewidth=0.8)
    p = p + geom_point(data=hcxval_combined[retro!="2022"&metric=="pred_cpue"&year==as.numeric(as.character(retro))+1],aes(x=year,y=value,fill=retro),shape=21,size=3)   
    p = p + geom_hline(yintercept=0) +
                    viridis::scale_color_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                    viridis::scale_fill_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                    theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                            panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                            panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                            strip.background =element_rect(fill="white"),
                            legend.key = element_rect(fill = "white"))     
    ggsave(filename="hcxval_analysis.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# model free hindcast
    target_model_dirs = sapply(c(diag_model,
                          "0136-dwfn-exeoFSTTGF-peel5_0",
                          "0137-dwfn-exeoFSTTGF-peel10_0",
                          "0138-dwfn-exeoFSTTGF-peel15_0",
                          "0139-dwfn-exeoFSTTGF-peel20_0"),function(x)file.path(model_stem,x))

    custom_params = get_default_params()
    # Model fits
        custom_params$fits$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
        custom_params$fits$active = TRUE  # TRUE | FALSE
        custom_params$fits$obs = TRUE  # TRUE | FALSE
        custom_params$fits$type = "Quantile"  # "Median" | "Spaghetti" | "Quantile"
        custom_params$fits$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
        custom_params$fits$resid = "PIT"  # "Ordinary" | "Standardized" | "PIT"
        custom_params$fits$ncol = 2
        custom_params$fits$resid_ncol = 1
        custom_params$fits$model_names = c("Diagnostic","-5yr (2017)", "-10yr (2012)", "-15yr (2007)", "-20yr (2002)")
    # Time series
        custom_params$ppts$var = c("Depletion (D)", "Population (P)", "F", "D_Dmsy", "F_Fmsy", "Removals", "Process error","Nominal CPUE","Catchability deviate")  # Any combination
        custom_params$ppts$show = "Posterior"  # "Prior" | "Posterior" | "Both"
        custom_params$ppts$combine = FALSE  # TRUE | FALSE
        custom_params$ppts$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
        custom_params$ppts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
        custom_params$ppts$ncol = 3
        custom_params$ppts$model_names = custom_params$fits$model_names
    
    
    # make plots
    p = generate_index_fit(target_model_dirs, params = custom_params$fits)
    p = p + geom_point(data=data.table(x=c(2017,2012,2007,2002),y=rep(0,4)),aes(x=x,y=y),
                        shape=21,color="black",size=3,
                        fill=viridis::viridis(5, begin = 0.1, end = 0.8, direction = -1, option = "H")[-1])
    ggsave(filename="model-free-hindcast.idx.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

    point_data = rbindlist(list(
                    data.table(
                    name = rep(custom_params$ppts$var[1],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[2],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[3],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[4],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[5],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[6],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(-0.2, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[7],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(-0.1, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[8],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    ),
                    data.table(
                    name = rep(custom_params$ppts$var[9],4),
                    x = c(2017, 2012, 2007, 2002),
                    y = rep(0, 4),
                    run_label = factor(custom_params$ppts$model_names[-1], levels = custom_params$ppts$model_names)  # Match existing levels
                    )
                    ))
    point_data[,name := factor(name,levels=custom_params$ppts$var)]
    
    p = generate_ppts(target_model_dirs, params = custom_params$ppts)
    p = p + geom_point(data=point_data,aes(x=x,y=y,fill=run_label),
                        shape=21,color="black",size=3)
    ggsave(filename="model-free-hindcast.ts.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width = 9, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# realized prior

    # Load data
    posterior_samples = fread(file.path(model_stem, diag_model, "hmc_samples.csv"))
    ppc_samples = fread(file.path(model_stem, diag_model, "hmc_samples_ppc.csv"))

    # Generate naive prior samples
    naive_prior_samples = ssp_prior_pushforward(
    ssp_summary = fread(file.path(model_stem, diag_model, "fit_summary.csv")),
    stan_data = fread(file.path(model_stem, diag_model, "stan_data.csv")),
    settings = fread(file.path(model_stem, diag_model, "settings.csv"))
    )

    # Variables of interest
    vars = c("logK", "r", "x0", "shape", "qeff", "rho", "sigma_qdev", "sigmap", "sigmao_add", "nu_catch")

    # Prepare data
    prep_data = function(dt, type_label, n_sample = 1000) {
    result = dt[variable %in% vars]
    if("divergent" %in% colnames(result)) result = result[divergent == 0]
    result[, value := as.numeric(as.character(value))]
    result = dcast(result, iter ~ variable, value.var = "value")
    result[, type := type_label]
    
    # Apply transformations
    if("qeff" %in% colnames(result)) result[, log_qeff := log(qeff)]
    if("rho" %in% colnames(result)) result[, atanh_rho := atanh(rho)]
    
    result[sample(.N, min(n_sample, .N))]
    }

    # Combine datasets
    plot_dt = rbind(
    prep_data(naive_prior_samples, "Naive Prior"),
    prep_data(ppc_samples, "Realized Prior"),
    prep_data(posterior_samples, "Posterior"),
    fill = TRUE
    )[, type := factor(type, levels = c("Naive Prior", "Realized Prior", "Posterior"))]

    # Update variables to include transformed versions
    plot_vars = c("logK", "r", "x0", "shape", "log_qeff", "atanh_rho", "sigma_qdev", "sigmap", "sigmao_add", "nu_catch")
    plot_cols = intersect(plot_vars, colnames(plot_dt))

    # Create pairs plot
    p_pairs = plot_dt %>%
    .[, .SD[sample(.N, min(1000, .N))], by = type] %>%
    ggpairs(columns = match(plot_cols, colnames(.)), aes(color = type, alpha = 0.4)) + 
    scale_color_viridis_d("Type", begin = 0.1, end = 0.8, option = "H") +
    scale_fill_viridis_d("Type", begin = 0.1, end = 0.8, option = "H") +
    theme(
        text = element_text(size = 20),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
        strip.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white")
    )

    ggsave(filename = "realized-prior-plot.png", plot = p_pairs, device = "png", 
        path = file.path(proj_dir, "plots"),
        scale = 1.25, width = 9, height = 9, units = "in",
        dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# ensemble data fits
 ensemble_models = c(diag_model,
                          "0102-nz-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0",
                          "0105-dwfn-exeoFSTTGF-cf0.2-nalt-qnewK-s1-o52b-o54b_0",
                          "0107-nz-exeoFSTTGF-cf0.2-nalt-qnewK-s1-o52b-o54b_0")

 target_model_dirs = sapply(ensemble_models,function(x)file.path(model_stem,x))

    custom_params = get_default_params()
    # Model fits
        custom_params$fits$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
        custom_params$fits$active = TRUE  # TRUE | FALSE
        custom_params$fits$obs = TRUE  # TRUE | FALSE
        custom_params$fits$type = "Quantile"  # "Median" | "Spaghetti" | "Quantile"
        custom_params$fits$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
        custom_params$fits$resid = "PIT"  # "Ordinary" | "Standardized" | "PIT"
        custom_params$fits$ncol = 1
        custom_params$fits$resid_ncol = 1
        custom_params$fits$model_names = c("Diagnostic","0102: NZ", "0105: DWFN & alt. shape", "0107: NZ & alt. shape")
    
    # make plots
    p = generate_index_fit(target_model_dirs, params = custom_params$fits)

    ggsave(filename="ensemble.idx_fit.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

    p = generate_catch_fit(target_model_dirs, params = custom_params$fits)

    ggsave(filename="ensemble.catch_fit.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# retrospective analysis - ensemble
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("fit_summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("fit_summary.csv","",all_dirs,fixed=TRUE)

    summary_dt = fread(file.path(proj_dir,"data","output","summary.csv")) %>%
                  .[!is.na(run_label)] %>%
                  .[run_label%in%ensemble_models]

    model_retro_dt.list = list()

    for(i in 1:nrow(summary_dt)){
        tmp_stem = paste0(strsplit(summary_dt$run_label[i],"_")[[1]][1])
        retro_peels = unname(sapply(all_dirs[grep(tmp_stem,all_dirs,fixed=TRUE)],function(x)as.numeric(gsub("/","",strsplit(x,"_")[[1]][2],fixed=TRUE))))
        length_retro_peels = length(retro_peels)-1
        
        if(length_retro_peels>0){
            retro_dt.list = list()

            for(j in 1:length(retro_peels)){
                retro_dt.list[[j]] = ssp_derived_quants_ts(ssp_summary=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/fit_summary.csv")),
                                                    samples_dt = fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/hmc_samples.csv")),
                                                    stan_data=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/stan_data.csv")),
                                                    settings=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/settings.csv")),
                                                    sub_sample_prop=1)  %>%
                                                    .[,retro:=retro_peels[j]] %>%
                                                    .[row%in%1:(max(row)-retro_peels[j])]
            }

            year_one=1952
            retro_base = rbindlist(retro_dt.list) %>%
                .[row>=1,year:=year_one+(row-1)] %>%
                .[row<1,year:=year_one+(row-1)] %>%
                .[name%in%c("D","F","D_Dmsy","F_Fmsy")] %>%
                .[,name:=factor(name,levels=c("D","D_Dmsy","F","F_Fmsy"),labels=c("D",expression("D"/"D"["MSY"]),"F",expression("F"/"F"["MSY"])))] %>%
                .[,.(med=median(value,na.rm=TRUE),avg=mean(value,na.rm=TRUE),lp=quantile(value,probs=0.025,na.rm=TRUE),up=quantile(value,probs=0.975,na.rm=TRUE),lp50=quantile(value,probs=0.25,na.rm=TRUE),up50=quantile(value,probs=0.75,na.rm=TRUE)),by=.(retro,name,year)] %>%
                .[,retro:=factor(retro,levels=0:5,labels=2022:2017)] %>%
                .[,model_name:=summary_dt$run_number[i]]

            model_retro_dt.list[[i]] = retro_base
        }
    }
       
    retro_combined = rbindlist(model_retro_dt.list) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]
    
    p = retro_combined %>% ggplot() +
        ylab("Metric") +
        xlab("Year") +
        facet_grid(name~short_plot_names,scales="free_y",labeller = labeller(name = label_parsed, model_name = label_value))
    p = p + geom_ribbon(data=retro_combined[retro=="2022"],aes(x=year,ymin=lp,ymax=up),fill="gray80",alpha=0.6,linewidth=1.15)
    p = p + geom_ribbon(data=retro_combined[retro=="2022"],aes(x=year,ymin=lp50,ymax=up50),fill="gray40",alpha=0.6,linewidth=1.15)
    p = p + geom_path(data=retro_combined[retro=="2022"],aes(x=year,y=med),color="black",linewidth=1.15)
    p = p + geom_path(data=retro_combined[retro!="2022"],aes(x=year,y=med,color=retro),linewidth=0.8)
    p = p + geom_point(data=retro_combined[retro!="2022"&year==as.numeric(as.character(retro))],aes(x=year,y=med,fill=retro),shape=21,size=3)
    p = p + geom_hline(yintercept=0) +
                viridis::scale_color_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE,labels = scales::parse_format()) +
                viridis::scale_fill_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE,labels = scales::parse_format()) +
                theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                        strip.background =element_rect(fill="white"),
                        legend.key = element_rect(fill = "white"))       
    ggsave(filename="ensemble.retro_analysis.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 12, units = c("in"),
            dpi = 300, limitsize = TRUE) 

#________________________________________________________________________________________________________________________________________________________________________________________________________
# hindcast cross-validation  - ensemble
    model_hcxval_dt.list = list()

    for(i in 1:nrow(summary_dt)){
        tmp_stem = paste0(strsplit(summary_dt$run_label[i],"_")[[1]][1])
        retro_peels = unname(sapply(all_dirs[grep(tmp_stem,all_dirs,fixed=TRUE)],function(x)as.numeric(gsub("/","",strsplit(x,"_")[[1]][2],fixed=TRUE))))
        length_retro_peels = length(retro_peels)-1
        
        if(length_retro_peels>0){
            hcxval_dt.list = list()

            for(j in 1:length(retro_peels)){
                hcxval_dt.list[[j]] = ssp_extract_cpue_fit(ssp_summary=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/fit_summary.csv")),
                                                    samples_dt = fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/hmc_samples.csv")),
                                                    stan_data=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/stan_data.csv")),
                                                    settings=fread(paste0(model_stem,"/",tmp_stem,"_",retro_peels[j],"/settings.csv")),
                                                    sub_sample_prop=1,
                                                        active="TRUE",
                                                        calc_std = "FALSE") %>%
                                                        .[metric %in% c("sigmao","obs_cpue","pred_cpue")] %>%
                                                        .[,retro:=retro_peels[j]]
            }

            year_one=1952
            hcxval_base = rbindlist(hcxval_dt.list) %>%
                .[row>=1,year:=year_one+(row-1)] %>%
                .[row<1,year:=year_one+(row-1)]

            hcxval_maxrow = max(na.omit(hcxval_base[metric=="obs_cpue"&retro==0])$row)
            hcxval_maxrow_overall = max(hcxval_base[metric=="obs_cpue"&retro==0]$row)

            hcxval_thin_dt.list = list()
            for(j_prime in (1+(hcxval_maxrow_overall-hcxval_maxrow)):(length(retro_peels)-1)){
                j = j_prime-(hcxval_maxrow_overall-hcxval_maxrow)
                hcxval_thin_dt.list[[j]] = rbind(hcxval_base[metric%in%c("obs_cpue","sigmao")&retro==0&row%in%1:(hcxval_maxrow-j+1),.(value=median(value)),by=.(run_id,metric,row,index,year,retro)],
                                                 hcxval_base[metric=="pred_cpue"&retro==retro_peels[j_prime+1]&row%in%1:(hcxval_maxrow-j+1),.(value=median(value)),by=.(run_id,metric,row,index,year,retro)])
            }
            hcxval_thin_dt = rbindlist(hcxval_thin_dt.list) %>%
                .[,retro:=factor(retro,levels=0:5,labels=2022:2017)] %>%
                .[,model_name:=summary_dt$run_number[i]]

            model_hcxval_dt.list[[i]] = hcxval_thin_dt
        }
    }
       
    hcxval_combined = rbindlist(model_hcxval_dt.list) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    obs_se_dt = hcxval_combined[metric=="sigmao"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        .[,se:=median(value),by=.(model_name,index,row,year)] %>%
        .[,se:=round(se,digits=3)] %>%
        .[,.(model_name,index,row,year,se)] %>%
        unique(.) %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    obs_cpue_dt = hcxval_combined[metric=="obs_cpue"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        merge(.,obs_se_dt,by=c("model_name","index","row","year")) %>%
        .[,obs:=round(value,digits=3)] %>%
        .[,upper:=qlnorm(0.975,meanlog=log(obs),sdlog=se)] %>%
        .[,lower:=qlnorm(0.025,meanlog=log(obs),sdlog=se)] %>%
                     .[,short_plot_names := sapply(model_name,function(x)strsplit(x,"-")[[1]][1])]

    p = hcxval_combined %>% ggplot() +
        ylab("Index") +
        xlab("Year") +
        geom_hline(yintercept=1,linetype="dashed") +
        facet_wrap(index_names[index]~short_plot_names,ncol=2)
    p = p + geom_segment(data=obs_cpue_dt,aes(x=year,xend=year,y=lower,yend=upper),linewidth=0.5)
    p = p + geom_point(data=obs_cpue_dt,aes(x=year,y=obs),color="black",fill="white",shape=21,size=3)
    p = p + geom_path(data=hcxval_combined[retro!="2022"&metric=="pred_cpue"],aes(x=year,y=value,color=retro),linewidth=0.8)
    p = p + geom_point(data=hcxval_combined[retro!="2022"&metric=="pred_cpue"&year==as.numeric(as.character(retro))+1],aes(x=year,y=value,fill=retro),shape=21,size=3)   
    p = p + geom_hline(yintercept=0) +
                    viridis::scale_color_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                    viridis::scale_fill_viridis("Peel",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                    theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                            panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                            panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                            strip.background =element_rect(fill="white"),
                            legend.key = element_rect(fill = "white"))     
    ggsave(filename="ensemble.hcxval_analysis.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
            scale = 1.25, width =9, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# ensemble data plots
 ensemble_models = c(diag_model,
                          "0102-nz-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0",
                          "0105-dwfn-exeoFSTTGF-cf0.2-nalt-qnewK-s1-o52b-o54b_0",
                          "0107-nz-exeoFSTTGF-cf0.2-nalt-qnewK-s1-o52b-o54b_0")

 target_model_dirs = sapply(ensemble_models,function(x)file.path(model_stem,x))

    custom_params = get_default_params()
    # Time series
        custom_params$ppts$var = c("Depletion (D)", "Population (P)","F", "D_Dmsy", "F_Fmsy", "Removals", "Process error","Nominal CPUE","Catchability deviate")  # Any combination
        custom_params$ppts$show = "Posterior"  # "Prior" | "Posterior" | "Both"
        custom_params$ppts$combine = TRUE  # TRUE | FALSE
        custom_params$ppts$prop = 1  # 0.01 to 1.00 (increments of 0.05)
        custom_params$ppts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
        custom_params$ppts$ncol = 3
        custom_params$ppts$model_names = c("Diagnostic","0102: NZ", "0105: DWFN & alt. shape", "0107: NZ & alt. shape")
        # Kobe & Majuro
        custom_params$kbmj$show = "Posterior"  # "Prior" | "Posterior" | "Both"
        custom_params$kbmj$combine = TRUE  # TRUE | FALSE
        custom_params$kbmj$prop = 1  # 0.01 to 1.00 (increments of 0.05)
        custom_params$kbmj$uncertainty = TRUE  # TRUE | FALSE
        custom_params$kbmj$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99
        custom_params$kbmj$resolution = 300  # 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500
        custom_params$kbmj$model_names = custom_params$ppts$model_names

        # Forecasts
        custom_params$forecasts$var = c("Depletion (D)","Population (P)", "F" ,"D_Dmsy", "F_Fmsy", "Removals", "Process error")  # Any combination
        custom_params$forecasts$combine = TRUE  # TRUE | FALSE
        custom_params$forecasts$prop = 1  # 0.01 to 1.00 (increments of 0.05)
        custom_params$forecasts$quants = 90  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
        custom_params$forecasts$nyears = 10  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
        custom_params$forecasts$resample_epsp = TRUE  # TRUE | FALSE
        custom_params$forecasts$type = "Catch"  # "Catch" | "U" | "MSY" | "Umsy"
        custom_params$forecasts$avg_year = 5  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        custom_params$forecasts$scalar = 1.0  # 0.01 to 5.0 (increments of 0.1)
        custom_params$forecasts$ncol = 3
        custom_params$forecasts$model_names = custom_params$ppts$model_names

    # make plots
    p = generate_ppts(target_model_dirs, params = custom_params$ppts)
    ggsave(filename="ensemble.ppts.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
     width = 12, height = 6, dpi = 300)

    p = generate_kb(target_model_dirs, params = custom_params$kbmj)
    ggsave(filename="ensemble.kb.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
     width = 9, height = 9, dpi = 300)
    
    p = generate_mj(target_model_dirs, params = custom_params$kbmj)
    ggsave(filename="ensemble.mj.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
     width = 9, height = 9, dpi = 300)
    
    p = generate_fcast(target_model_dirs, params = custom_params$forecasts)
    ggsave(filename="ensemble.forecast_5yr_meancatch.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
     width = 12, height = 6, dpi = 300)

    # extract specific forecast data
        all_data = lapply(target_model_dirs, load_model_data)
        # Generate forecasts for each model
        derived_dt_list = posterior_dt_list = forecast_posterior_dt_list = forecast_dt_list = list()
        for (i in seq_along(all_data)) {
            # Run forecasts
            forecast_dt_list[[i]] = ssp_forecast(
            ssp_summary = all_data[[i]]$summary,
            samples_dt = all_data[[i]]$samples,
            stan_data = all_data[[i]]$stan_data,
            settings = all_data[[i]]$settings,
            sub_sample_prop = custom_params$forecasts$prop,
            forecast_years = 30,
            forecast_type = custom_params$forecasts$type,
            avg_years = custom_params$forecasts$avg_year,
            scalar = custom_params$forecasts$scalar,
            resample_raw_epsp = ifelse(custom_params$forecasts$resample_epsp, "TRUE", "FALSE")
            )
            
            forecast_posterior_dt_list[[i]] = ssp_derived_quants_ts(
            ssp_summary = all_data[[i]]$summary,
            samples_dt = forecast_dt_list[[i]],
            stan_data = all_data[[i]]$stan_data,
            settings = all_data[[i]]$settings,
            sub_sample_prop = 1
            )

            posterior_dt_list[[i]] = ssp_derived_quants_ts(
            ssp_summary = all_data[[i]]$summary,
            samples_dt = all_data[[i]]$samples,
            stan_data = all_data[[i]]$stan_data,
            settings = all_data[[i]]$settings,
            sub_sample_prop = 1
            )

            derived_dt_list[[i]] = ssp_derived_quants(
            hmc_samples = all_data[[i]]$samples,
            stan_data = all_data[[i]]$stan_data,
            output="raw",
            percentile=0.5
            )
        }

        tmp_summary = rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)

        yo_dt = tmp_summary[,.(run_label)] %>%
            unique(.) %>%
            .[,year_one:=extract_model_start_year(run_label)]

        forecast_dt = rbindlist(forecast_dt_list) %>% merge(., tmp_summary[, .(run_id, run_label)]) %>%
                merge(.,yo_dt,by="run_label") %>%  
                .[row >= 1, year := year_one + (row - 1)] %>%
                .[row < 1, year := year_one + (row - 1)] %>%
                .[name %in% c("Process error", "Process error (raw)") & row > 0, year := year + 1]
        forecast_posterior_dt = rbindlist(forecast_posterior_dt_list) %>% merge(., tmp_summary[, .(run_id, run_label)]) %>%
                merge(.,yo_dt,by="run_label") %>%  
                .[row >= 1, year := year_one + (row - 1)] %>%
                .[row < 1, year := year_one + (row - 1)] %>%
                .[name %in% c("Process error", "Process error (raw)") & row > 0, year := year + 1]
        posterior_dt = rbindlist(posterior_dt_list) %>% merge(., tmp_summary[, .(run_id, run_label)]) %>%
                merge(.,yo_dt,by="run_label") %>%  
                .[row >= 1, year := year_one + (row - 1)] %>%
                .[row < 1, year := year_one + (row - 1)] %>%
                .[name %in% c("Process error", "Process error (raw)") & row > 0, year := year + 1]

        forecast_summary = forecast_posterior_dt[year>2022 & name %in% c("D","P","F","D_Dmsy","F_Fmsy","removals"), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name,year)] %>%
                    .[order(year, name)]

        mgmt_summary = posterior_dt[name %in% c("D","P","F","D_Dmsy","F_Fmsy","removals"), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name,year)] %>%
                    .[order(year, name)]

        # Define the conditions for each parameter
        forecast_probs = forecast_posterior_dt[
                        (name == "D_Dmsy" & value < 1) |
                        (name == "F_Fmsy" & value > 1) |
                        (name == "D" & value < 0.05)
                        ][, .(
                        count_meeting_condition = .N,
                        prop_meeting_condition = .N / forecast_posterior_dt[name == first(name) & year == first(year), .N]
                        ), by = .(name, year)] %>% .[order(year, name)] %>% .[year>2022]

        # calculate that D and F continue to increase/decrease
                D_gt_2022 = forecast_posterior_dt[name == "D" & year %in% c(2022, 2023:2052),
                        .(ref_val = value[year == 2022], 
                        curr_val = value[year > 2022],
                        curr_year = year[year > 2022]),
                        by = .(iter, chain, run_label)
                    ][, .(name = "D_gt_2022",
                        count_meeting_condition = sum(curr_val > ref_val, na.rm = TRUE),
                        prop_meeting_condition = mean(curr_val > ref_val, na.rm = TRUE)),
                    by = curr_year
                    ][order(curr_year)] %>%
                    setnames(.,"curr_year","year") %>%
                    .[,.(name,year,count_meeting_condition,prop_meeting_condition)]

                F_lt_2021 = forecast_posterior_dt[name == "F" & year %in% c(2021, 2021:2052),
                        .(ref_val = value[year == 2021], 
                        curr_val = value[year > 2021],
                        curr_year = year[year > 2021]),
                        by = .(iter, chain, run_label)
                    ][, .(name = "F_lt_2021",
                        count_meeting_condition = sum(curr_val < ref_val, na.rm = TRUE),
                        prop_meeting_condition = mean(curr_val < ref_val, na.rm = TRUE)),
                    by = curr_year
                    ][order(curr_year)] %>%
                    setnames(.,"curr_year","year") %>%
                    .[,.(name,year,count_meeting_condition,prop_meeting_condition)]
            
            forecast_probs = rbind(forecast_probs,D_gt_2022,F_lt_2021)

        mgmt_probs = posterior_dt[
                        (name == "D_Dmsy" & value < 1) |
                        (name == "F_Fmsy" & value > 1) |
                        (name == "D" & value < 0.05)
                        ][, .(
                        count_meeting_condition = .N,
                        prop_meeting_condition = .N / forecast_posterior_dt[name == first(name) & year == first(year), .N]
                        ), by = .(name, year)] %>% .[order(year, name)] %>% .[year>2022]

        # calc recent and latest removals
        latest_removals = posterior_dt[name %in% c("removals") & year %in% max(year), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name)] %>%
                    .[order(name)] %>%
                    .[,name:="latest_removals"]
        recent_removals = posterior_dt[name %in% c("removals") & year %in% seq(from=max(year),by=-1,length.out=4), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name)] %>%
                    .[order(name)] %>%
                    .[,name:="recent_removals"]        

        derived_dt = rbindlist(derived_dt_list) %>% melt(.,id.vars=c("run_id","iter")) %>%
                     setnames(.,"variable","name") %>%
                     rbind(., posterior_dt[name %in% c("removals") & year %in% max(year),.(name="latest_removals",value=mean(value)),by=.(run_id,iter)]) %>%
                     rbind(.,posterior_dt[name %in% c("removals") & year %in% seq(from=max(year),by=-1,length.out=4),.(name="recent_removals",value=mean(value)),by=.(run_id,iter)])

        derived_summary = derived_dt[name %in% c("msy","Dmsy","Fmsy","Pmsy","Plrp","latest_F","latest_depletion","latest_population",
                                                "latest_F_Fmsy","latest_D_Dmsy","recent_F_Fmsy","recent_D_Dmsy",
                                                "latest_D_Dlrp","recent_D_Dlrp","recent_F","recent_D","recent_P"), .(
                    mean = mean(value,na.rm=TRUE),
                    median = median(value,na.rm=TRUE),
                    q025 = quantile(value, 0.025,na.rm=TRUE),
                    q975 = quantile(value, 0.975,na.rm=TRUE),
                    sd = sd(value,na.rm=TRUE),
                    n = .N
                    ),by=.(name)] %>%
                    .[order(name)]
        
        derived_summary = rbind(derived_summary,latest_removals,recent_removals)
        
        derived_probs = derived_dt[
                        (name == "latest_D_Dmsy" & value < 1) |
                        (name == "recent_D_Dmsy" & value < 1) |
                        (name == "latest_F_Fmsy" & value > 1) |
                        (name == "recent_F_Fmsy" & value > 1) |
                        (name == "latest_D_Dlrp" & value < 1) |
                        (name == "recent_D_Dlrp" & value < 1)
                        ][, .(
                        count_meeting_condition = .N,
                        prop_meeting_condition = .N / derived_dt[name == first(name), .N]
                        ), by = .(name)] %>% .[order(name)]
        
        joint_probs = derived_dt[name %in% c("latest_D_Dmsy", "latest_F_Fmsy", "recent_D_Dmsy", "recent_F_Fmsy","latest_D_Dlrp","recent_D_Dlrp"),
                                .(latest_joint_kb = (value[name == "latest_D_Dmsy"] < 1) & (value[name == "latest_F_Fmsy"] > 1),
                                    recent_joint_kb = (value[name == "recent_D_Dmsy"] < 1) & (value[name == "recent_F_Fmsy"] > 1),
                                    latest_joint_mj = (value[name == "latest_D_Dlrp"] < 1) & (value[name == "latest_F_Fmsy"] > 1),
                                    recent_joint_mj = (value[name == "recent_D_Dlrp"] < 1) & (value[name == "recent_F_Fmsy"] > 1)),
                                by = .(iter, run_id)
                                ][, .(
                                rbind(
                                    .(name = "latest_joint_kb",
                                        count_meeting_condition = sum(latest_joint_kb),
                                        prop_meeting_condition = mean(latest_joint_kb)),
                                    .(name = "recent_joint_kb",
                                        count_meeting_condition = sum(recent_joint_kb),
                                        prop_meeting_condition = mean(recent_joint_kb)),
                                    .(name = "latest_joint_mj",
                                        count_meeting_condition = sum(latest_joint_mj),
                                        prop_meeting_condition = mean(latest_joint_mj)),
                                    .(name = "recent_joint_mj",
                                        count_meeting_condition = sum(recent_joint_mj),
                                        prop_meeting_condition = mean(recent_joint_mj))   
                                )
                                )]

        fwrite(forecast_summary,file=file.path(proj_dir,"data","output","ens_forecast_summary.csv"))
        fwrite(mgmt_summary,file=file.path(proj_dir,"data","output","ens_mgmt_summary.csv"))
        fwrite(forecast_probs,file=file.path(proj_dir,"data","output","ens_forecast_probs.csv"))
        fwrite(mgmt_probs,file=file.path(proj_dir,"data","output","ens_mgmt_probs.csv"))
        fwrite(derived_dt,file=file.path(proj_dir,"data","output","ens_derived_dt.csv"))
        fwrite(derived_summary,file=file.path(proj_dir,"data","output","ens_derived_summary.csv"))
        fwrite(derived_probs,file=file.path(proj_dir,"data","output","ens_derived_probs.csv"))
        fwrite(joint_probs,file=file.path(proj_dir,"data","output","ens_joint_probs.csv"))

