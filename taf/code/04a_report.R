

# NOAA PIFSC
# 2025/09/25
# Plot results from the diagnostic model
# Model 0100: DWFN index; baseline priors

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
# define model directories
    ssp_scenarios = fread(file.path("data","ssp_scenarios.csv"))
    model_list = apply(ssp_scenarios,1,function(x)paste0("0",x[[1]],"-",x[[2]],"-",x[[3]],"_0"))[1]
    
    model_dirs = file.path(dir_model_runs, model_list)
    
    # verify model directories exist
    missing_models = model_dirs[!dir.exists(model_dirs)]
    if(length(missing_models) > 0) {
        stop("Missing model directories: ", paste(missing_models, collapse = ", "))
    }

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set custom parameters for analysis
    custom_params = get_default_params()
    
    # HMC diagnostics
    custom_params$hmc$leading_params = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","x0","nu_catch")  # Any combination
    custom_params$hmc$raw = FALSE  # TRUE (transformed) | FALSE (raw)
    custom_params$hmc$diag = "Divergences"  # "None" | "Divergences" | "Max. treedepth"
    custom_params$hmc$eps = FALSE  # TRUE | FALSE
    custom_params$hmc$lags = 30  # 5, 10, 15, 20, 25, 30, 35, 40, 45, 50
    custom_params$hmc$scheme = "brightblue"  # "blue" | "brightblue" | "gray" | "darkgray" | "green" | "pink" | "purple" | "red" | "teal" | "yellow" | "viridis" | "viridisA" | "viridisB" | "viridisC" | "viridisD" | "viridisE"

    # PPC settings
    custom_params$ppc$scheme = "brightblue"  # "blue" | "brightblue" | "gray" | "darkgray" | "green" | "pink" | "purple" | "red" | "teal" | "yellow" | "viridis" | "viridisA" | "viridisB" | "viridisC" | "viridisD" | "viridisE"
    custom_params$ppc$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$ppc$active = TRUE  # TRUE | FALSE
    custom_params$ppc$group = FALSE  # TRUE (aggregate) | FALSE (group by index)
    custom_params$ppc$stat = "median"  # "mean" | "median" | "sd" | "mad" | c("mean", "sd") | c("median", "mad") | etc.
    custom_params$ppc$qqdist = "uniform"  # "uniform" | "normal"

    # Model fits
    custom_params$fits$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$fits$active = TRUE  # TRUE | FALSE
    custom_params$fits$obs = TRUE  # TRUE | FALSE
    custom_params$fits$type = "Quantile"  # "Median" | "Spaghetti" | "Quantile"
    custom_params$fits$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$fits$resid = "PIT"  # "Ordinary" | "Standardized" | "PIT"
    custom_params$fits$ncol = 2
    custom_params$fits$resid_ncol = 1
    custom_params$fits$model_names = c("Diagnostic")

    # Prior-posterior parameters
    custom_params$ppp$leading_params = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","x0","nu_catch")  # Any combination
    custom_params$ppp$raw = TRUE  # TRUE (transformed) | FALSE (raw)
    custom_params$ppp$show = "Both"  # "Prior" | "Posterior" | "Both"
    custom_params$ppp$combine = FALSE  # TRUE | FALSE
    custom_params$ppp$ncol = 3
    custom_params$ppp$model_names = c("Diagnostic")

    # Time series
    custom_params$ppts$var = c("Depletion (D)", "Population (P)", "D_Dmsy", "F_Fmsy", "Removals", "Process error","Nominal CPUE","Catchability deviate")  # Any combination
    custom_params$ppts$show = "Posterior"  # "Prior" | "Posterior" | "Both"
    custom_params$ppts$combine = FALSE  # TRUE | FALSE
    custom_params$ppts$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$ppts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$ppts$ncol = 3
    custom_params$ppts$model_names = c("Diagnostic")

    # Kobe & Majuro
    custom_params$kbmj$show = "Posterior"  # "Prior" | "Posterior" | "Both"
    custom_params$kbmj$combine = FALSE  # TRUE | FALSE
    custom_params$kbmj$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$kbmj$uncertainty = TRUE  # TRUE | FALSE
    custom_params$kbmj$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99
    custom_params$kbmj$resolution = 300  # 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500
    custom_params$kbmj$model_names = c("Diagnostic")

    # Forecasts
    custom_params$forecasts$var = c("Depletion (D)","Population (P)", "D_Dmsy", "F_Fmsy", "Removals", "Process error")  # Any combination
    custom_params$forecasts$combine = FALSE  # TRUE | FALSE
    custom_params$forecasts$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$forecasts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$forecasts$nyears = 10  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    custom_params$forecasts$resample_epsp = TRUE  # TRUE | FALSE
    custom_params$forecasts$type = "Catch"  # "Catch" | "U" | "MSY" | "Umsy"
    custom_params$forecasts$avg_year = 5  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    custom_params$forecasts$scalar = 1.0  # 0.01 to 5.0 (increments of 0.1)
    custom_params$forecasts$ncol = 3
    custom_params$forecasts$model_names = c("Diagnostic")


#________________________________________________________________________________________________________________________________________________________________________________________________________
# generate all plots with batch processing

    # generate complete analysis
    all_plots = generate_all_plots(
        model_dirs = model_dirs,
        output_dir = dir_plots,
        params = custom_params,
        save_plots = TRUE,
        plot_format = "png",
        width = 14,
        height = 10,
        dpi = 300,
        parallel = TRUE,
        n_cores = 1,
        comparison_only = FALSE
    )
