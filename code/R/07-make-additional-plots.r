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

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
    model_stem = file.path(proj_dir,"data","output","model_runs")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# retrospective analysis
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("fit_summary.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("fit_summary.csv","",all_dirs,fixed=TRUE)

    summary_dt = fread(file.path(proj_dir,"data","output","summary.csv")) %>%
                  .[!is.na(run_label)]

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
                .[name%in%c("D","U","D_Dmsy","U_Umsy")] %>%
                .[,name:=factor(name,levels=c("D","D_Dmsy","U","U_Umsy"),labels=c("D",expression("D"/"D"["MSY"]),"U",expression("U"/"U"["MSY"])))] %>%
                .[,.(med=median(value),avg=mean(value),lp=quantile(value,probs=0.025),up=quantile(value,probs=0.975),lp50=quantile(value,probs=0.25),up50=quantile(value,probs=0.75)),by=.(retro,name,year)] %>%
                .[,retro:=factor(retro,levels=0:5,labels=2022:2017)] %>%
                .[,model_name:=summary_dt$run_number[i]]

            model_retro_dt.list[[i]] = retro_base
        }
    }
       
    retro_combined = rbindlist(model_retro_dt.list)
    
    p = retro_combined %>% ggplot() +
        ylab("Metric") +
        xlab("Year") +
        facet_grid(name~model_name,scales="free_y",labeller = labeller(name = label_parsed, model_name = label_value))
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
# hindcast cross-validation
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
       
    hcxval_combined = rbindlist(model_hcxval_dt.list)

    obs_se_dt = hcxval_combined[metric=="sigmao"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        .[,se:=median(value),by=.(model_name,index,row,year)] %>%
        .[,se:=round(se,digits=3)] %>%
        .[,.(model_name,index,row,year,se)] %>%
        unique(.)

    obs_cpue_dt = hcxval_combined[metric=="obs_cpue"] %>%
        .[,.(model_name,index,row,year,value)] %>%
        merge(.,obs_se_dt,by=c("model_name","index","row","year")) %>%
        .[,obs:=round(value,digits=3)] %>%
        .[,upper:=qlnorm(0.975,meanlog=log(obs),sdlog=se)] %>%
        .[,lower:=qlnorm(0.025,meanlog=log(obs),sdlog=se)]

    p = hcxval_combined %>% ggplot() +
        ylab("Index") +
        xlab("Year") +
        geom_hline(yintercept=1,linetype="dashed") +
        facet_wrap(~model_name,ncol=1)
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
            scale = 1.25, width =6, height = 12, units = c("in"),
            dpi = 300, limitsize = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# fit to catch
    all_dirs = list.files(model_stem,recursive = TRUE)
    all_dirs = all_dirs[grep("hmc_samples.csv",all_dirs,fixed=TRUE)]
    all_dirs = gsub("hmc_samples.csv","",all_dirs,fixed=TRUE)
    all_dirs = all_dirs[grep("_0/",all_dirs,fixed=TRUE)]

    samples_ts.list = as.list(rep(NA,length(all_dirs))) 
    for(i in 1:length(all_dirs)){
        tmp_samples = fread(paste0(model_stem,"/",all_dirs[i],"hmc_samples.csv"))
        tmp_data = fread(paste0(model_stem,"/",all_dirs[i],"stan_data.csv"))

             
        tmp_ts = ssp_derived_quants_ts(ssp_summary=fread(paste0(model_stem,"/",all_dirs[i],"fit_summary.csv")),
                            samples_dt=tmp_samples,
                            settings=fread(paste0(model_stem,"/",all_dirs[i],"settings.csv")),
                            stan_data=tmp_data,
                            sub_sample_prop=1)
        tmp_catch = tmp_data[name=="obs_removals"] %>%
                    .[,.(row,value)] %>%
                    setnames(.,"value","obs_removals")
        
        # extract run_id from directory name
        run_id = paste0(strsplit(all_dirs[i],"/")[[1]][1])
        
        samples_ts.list[[i]] = merge(tmp_ts[name%in%c("removals"),.(run_id,iter,name,row,value)],tmp_catch,by="row") %>%
                               .[,run_id:=run_id] %>%
                               merge(.,summary_dt[,.(run_id,run_number)],by="run_id")

    }

    samples_ts = rbindlist(samples_ts.list) %>%
              .[,model_name:=run_number]

    year_one = 1952
    p = samples_ts %>%
        .[,.(obs=mean(obs_removals),med=median(value),avg=mean(value),lp=quantile(value,probs=0.025),up=quantile(value,probs=0.975)),by=.(model_name,name,row)] %>%
        .[row>=1,year:=year_one+(row-1)] %>%
        .[row<1,year:=year_one+(row-1)] %>%
        ggplot() +
        ylab("Catch") +
        xlab("Year") +
        facet_wrap(~model_name) +
        geom_ribbon(aes(x=year,ymin=lp,ymax=up,fill=model_name,group=model_name),alpha=0.2) +
        geom_point(aes(x=year,y=obs),size=3) +
        geom_line(aes(x=year,y=med,color=model_name),linewidth=1.15) +
        viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                                        panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                                        panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                                        strip.background =element_rect(fill="white"),
                                        legend.key = element_rect(fill = "white"))
            ggsave(filename="predicted_catch.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
                scale = 1.25, width =12, height = 9, units = c("in"),
                dpi = 300, limitsize = TRUE) 

#________________________________________________________________________________________________________________________________________________________________________________________________________
# index ppd
    samples_ts.list = as.list(rep(NA,length(all_dirs)))
    for(i in 1:length(all_dirs)){
        tmp_samples = fread(paste0(model_stem,"/",all_dirs[i],"hmc_samples.csv"))
        tmp_data = fread(paste0(model_stem,"/",all_dirs[i],"stan_data.csv"))

             
        ssp_summary=fread(paste0(model_stem,"/",all_dirs[i],"fit_summary.csv"))
        stan_data=fread(paste0(model_stem,"/",all_dirs[i],"stan_data.csv"))
        settings=fread(paste0(model_stem,"/",all_dirs[i],"settings.csv"))
        samples_dt=fread(paste0(model_stem,"/",all_dirs[i],"hmc_samples.csv"))
        fit_dt = ssp_extract_cpue_fit(ssp_summary=ssp_summary,
                            samples_dt=samples_dt,
                            stan_data=stan_data,
                            settings=settings,
                            sub_sample_prop=1,
                            active="TRUE",
                            calc_std = "FALSE")
        
        # extract run_id from directory name
        run_id = paste0(strsplit(all_dirs[i],"/")[[1]][1])
        
        samples_ts.list[[i]] = fit_dt %>%
                               .[,run_id:=run_id] %>%
                               merge(.,summary_dt[,.(run_id,run_number)],by="run_id") %>%
                               .[index==which(stan_data[name=="lambdas"]$value==1)]

    }
    year_one=1952
    samples_ts = rbindlist(samples_ts.list) %>%
              .[,model_name:=run_number] %>%
              .[row>=1,year:=year_one+(row-1)] %>%
              .[row<1,year:=year_one+(row-1)]
    
    obs_se_dt = samples_ts[metric=="sigmao"] %>%
                  .[,.(model_name,index,row,year,value)] %>%
                  .[,se:=median(value),by=.(model_name,index,row,year)] %>%
                  .[,se:=round(se,digits=3)] %>%
                  .[,.(model_name,index,row,year,se)] %>%
                  unique(.)

    obs_cpue_dt = samples_ts[metric=="obs_cpue"] %>%
                  .[,.(model_name,index,row,year,value)] %>%
                  merge(.,obs_se_dt[,.(model_name,index,row,year,se)],by=c("model_name","index","row","year")) %>%
                  .[,obs:=round(value,digits=3)] %>%
                  .[,upper:=qlnorm(0.975,meanlog=log(obs),sdlog=se)] %>%
                  .[,lower:=qlnorm(0.025,meanlog=log(obs),sdlog=se)]
        
    pred_cpue_dt = samples_ts[metric=="ppd_cpue",.(model_name,index,row,year,iter,value)] %>%
                    .[,.(median=median(value),upper=quantile(value,probs=0.975),lower=quantile(value,probs=0.025)),by=.(model_name,index,row,year)]

     p = pred_cpue_dt %>%
      ggplot() +
        ylab("Index") +
        xlab("Year") +
        geom_hline(yintercept=1,linetype="dashed") +
        facet_wrap(~model_name)
 
        p = p + geom_segment(data=obs_cpue_dt,aes(x=year,xend=year,y=lower,yend=upper),linewidth=0.5)
        p = p + geom_ribbon(aes(x=year,ymin=lower,ymax=upper,fill=as.character(index)),alpha=0.2) +
                geom_line(aes(x=year,y=median,color=as.character(index)),linewidth=1)
      
      p = p + geom_point(data=obs_cpue_dt,aes(x=year,y=obs),color="white",fill="black",shape=21,size=2)

      p = p + geom_hline(yintercept=0) +
                viridis::scale_color_viridis("Index",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                viridis::scale_fill_viridis("Index",begin = 0.1,end = 0.8,direction = -1,option = "H",discrete=TRUE,drop=FALSE) +
                theme(text = element_text(size = 20),panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                        panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                        panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                        strip.background =element_rect(fill="white"),
                        legend.key = element_rect(fill = "white"))
            ggsave(filename="ppd_index.png", plot = p, device = "png", path = file.path(proj_dir,"plots"),
                scale = 1.25, width =9, height = 12, units = c("in"),
                dpi = 300, limitsize = TRUE)
