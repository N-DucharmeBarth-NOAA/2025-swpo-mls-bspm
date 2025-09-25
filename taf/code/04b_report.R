

# NOAA PIFSC
# 2025/09/25
# Plot results from the model ensemble

# Copyright (c) 2025 NOAA PIFSC
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.



#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(viridis)
    library(bayesplot)
    library(GGally)
    library(MASS)
    library(randtests)
    library(loo)
    library(rstantools)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    source_dir_stem = "./boot/software/"
    source_files = list.files(source_dir_stem)
    source_files = source_files[grep(".r",source_files,fixed=TRUE)]
    sapply(paste0(source_dir_stem,source_files),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for plots
    dir_plots = "report"
    dir_model_runs = "model"
    mkdir(dir_plots)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# configure global settings
  set_global_config(
    index_names = c("DWFN","AU","NZ","Obs (all)","Obs (NC,FJ&TO)","Obs (PF)"), 
    model_stem = dir_model_runs,
    height_per_panel = 350
  )

#________________________________________________________________________________________________________________________________________________________________________________________________________
# ensemble data plots
    ssp_scenarios = fread(file.path("data","ssp_scenarios.csv"))
    ensemble_models = apply(ssp_scenarios,1,function(x)paste0("0",x[[1]],"-",x[[2]],"-",x[[3]],"_0"))

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
    ggsave(filename="ensemble.ppts.png", plot = p, device = "png", path = file.path(dir_plots),
     width = 12, height = 6, dpi = 300)

    p = generate_kb(target_model_dirs, params = custom_params$kbmj)
    ggsave(filename="ensemble.kb.png", plot = p, device = "png", path = file.path(dir_plots),
     width = 9, height = 9, dpi = 300)
    
    p = generate_mj(target_model_dirs, params = custom_params$kbmj)
    ggsave(filename="ensemble.mj.png", plot = p, device = "png", path = file.path(dir_plots),
     width = 9, height = 9, dpi = 300)
    
    p = generate_fcast(target_model_dirs, params = custom_params$forecasts)
    ggsave(filename="ensemble.forecast_5yr_meancatch.png", plot = p, device = "png", path = file.path(dir_plots),
     width = 12, height = 6, dpi = 300)
