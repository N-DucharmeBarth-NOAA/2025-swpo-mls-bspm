# Nicholas Ducharme-Barth
# 2025/06/04
# Generate All SSP Model Analysis Plots

# Copyright (c) 2025 Nicholas Ducharme-Barth
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
# define directories
    proj_dir = this.path::this.proj()
    dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
    dir_plot_fns = file.path(proj_dir,"code","R","plot-fns")
    dir_model_runs = file.path(proj_dir,"data","output","model_runs")
    dir_plots = file.path(proj_dir,"plots","sens-0006")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source plot functions
    sapply(file.path(dir_plot_fns,(list.files(dir_plot_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# make directory for plots
    dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# configure global settings
    set_global_config(
        year_one = 1952,
        index_names = c("DWFN CPUE", "DWFN CPUE"),
        model_stem = dir_model_runs,
        height_per_panel = 350
    )

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define model directories
    model_list = c(
        "0005-2024cpueMVPrior_0",
        "0006-2024cpueEffortQeff_0"
    )
    
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
    custom_params$hmc$leading_params = c("logK", "r", "sigmao_add", "sigmap", "shape", "sigmaf","qeff","rho","sigma_qdev")  # Any combination
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

    # Prior-posterior parameters
    custom_params$ppp$leading_params = c("logK", "r", "sigmao_add", "sigmap", "shape", "sigmaf","qeff","rho","sigma_qdev")  # Any combination
    custom_params$ppp$raw = TRUE  # TRUE (transformed) | FALSE (raw)
    custom_params$ppp$show = "Both"  # "Prior" | "Posterior" | "Both"
    custom_params$ppp$combine = FALSE  # TRUE | FALSE
    custom_params$ppp$ncol = 3
    custom_params$ppp$nrow = 4

    # Time series
    custom_params$ppts$var = c("Depletion (D)", "Population (P)", "D_Dmsy", "F_Fmsy", "Removals", "Process error","Nominal CPUE","Effort deviate","Catchability deviate")  # Any combination
    custom_params$ppts$show = "Posterior"  # "Prior" | "Posterior" | "Both"
    custom_params$ppts$combine = FALSE  # TRUE | FALSE
    custom_params$ppts$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$ppts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$ppts$ncol = 3
    custom_params$ppts$nrow = 3  # Let ggplot calculate rows

    # Kobe & Majuro
    custom_params$kbmj$show = "Posterior"  # "Prior" | "Posterior" | "Both"
    custom_params$kbmj$combine = FALSE  # TRUE | FALSE
    custom_params$kbmj$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$kbmj$uncertainty = TRUE  # TRUE | FALSE
    custom_params$kbmj$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99
    custom_params$kbmj$resolution = 300  # 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500

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
    custom_params$forecasts$nrow = 2  # Let ggplot calculate rows

#________________________________________________________________________________________________________________________________________________________________________________________________________
# generate all plots with batch processing
    cat("Starting SSP Model Analysis Plot Generation\n")
    cat("==========================================\n")
    cat("Models to analyze:", length(model_dirs), "\n")
    cat("Base output directory:", dir_plots, "\n\n")
    
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
        n_cores = 2,
        comparison_only = FALSE
    )

#________________________________________________________________________________________________________________________________________________________________________________________________________
# save analysis metadata
    cat("\nSaving analysis metadata...\n")
    
    analysis_metadata = list(
        analysis_date = Sys.time(),
        models_analyzed = model_list,
        model_directories = model_dirs,
        parameters_used = custom_params,
        output_directory = dir_plots,
        plot_format = "png",
        n_models = length(model_dirs),
        functions_used = c(
            "HMC diagnostics: 6 functions",
            "PPC diagnostics: 7 functions", 
            "Model fits: 3 functions",
            "Management plots: 2 functions",
            "Forecasts: 1 function"
        )
    )
    
    saveRDS(analysis_metadata, file = file.path(dir_plots, "analysis_metadata.rds"))
    
    # save session info for reproducibility
    writeLines(
        capture.output(sessionInfo()),
        con = file.path(dir_plots, "session_info.txt")
    )

#________________________________________________________________________________________________________________________________________________________________________________________________________
# print completion summary
    cat("\n" , rep("=", 60), "\n", sep = "")
    cat("SSP MODEL ANALYSIS COMPLETED SUCCESSFULLY\n")
    cat(rep("=", 60), "\n", sep = "")
    cat("Models analyzed:", length(model_dirs), "\n")
    cat("Total plots generated:", 
        length(list.files(dir_plots, pattern = "\\.png$", recursive = TRUE)), "\n")
    cat("Analysis metadata saved to: analysis_metadata.rds\n")
    cat("\nKey output files:\n")
    cat("- Individual model diagnostics: ./", basename(dir_plots), "/[model_name]/\n", sep = "")
    cat("- Multi-model comparisons: ./", basename(dir_plots), "/comparisons/\n", sep = "")
    cat("- Key comparison plots: ./", basename(dir_plots), "/key_comparison_*.png\n", sep = "")
    cat("- Presentation plots: ./", basename(dir_plots), "/presentation_*.png\n", sep = "")
    cat(rep("=", 60), "\n", sep = "")
