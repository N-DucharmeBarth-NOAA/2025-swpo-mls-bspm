# Nicholas Ducharme-Barth
# 2025/07/20
# Generate plots for one-off sensitivity models by group

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(GGally)
library(viridis)
library(parallel)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
proj_dir = this.path::this.proj()
dir_helper_fns = file.path(proj_dir,"code","R","helper-fns")
dir_plot_fns = file.path(proj_dir,"code","R","plot-fns")
model_stem = file.path(proj_dir,"data","output","model_runs")
dir_output = file.path(proj_dir,"data","output")

set_global_config(
    index_names = c("DWFN","AU","NZ","Obs (all)","Obs (NC,FJ&TO)","Obs (PF)"), 
    model_stem = file.path("data","output","model_runs"),
    height_per_panel = 350
)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
sapply(file.path(dir_helper_fns,(list.files(dir_helper_fns))),source)
sapply(file.path(dir_plot_fns,(list.files(dir_plot_fns))),source)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load model configuration
model_desc_config = fread(file.path(dir_output,"model_desc_config.csv"))
model_desc_config[,short_name := sprintf("%04d", short_name)]

# diagnostic model info
diag_model = "0100-dwfn-exeoFSTTGF-cf0.2-nb-qnewK-s1-o52b-o54b_0"
diag_dir = file.path(model_stem, diag_model)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# create and save short labels CSV
# create_short_labels <- function(descriptions) {
#     labels = sapply(descriptions, function(desc) {
#         # Extract key terms and abbreviate
#         desc = gsub("Fits to the ", "", desc)
#         desc = gsub(" index", "", desc)
#         desc = gsub(" observer data", " obs", desc)
#         desc = gsub("New Caledonia, Fiji, Tonga", "NC/FJ/TO", desc)
#         desc = gsub("French Polynesia", "FP", desc)
#         desc = gsub("Australia", "AU", desc)
#         desc = gsub("New Zealand", "NZ", desc)
#         desc = gsub("longline ", "", desc)
#         desc = gsub("Combined alternative ", "", desc)
#         desc = gsub("Recalculates ", "", desc)
#         desc = gsub("Time-varying ", "Time-var ", desc)
#         desc = gsub("Constant ", "", desc)
#         desc = gsub("Applied ", "", desc)
#         desc = gsub("Model starts in ", "Start ", desc)
#         desc = gsub("Direct F estimation with ", "Direct F ", desc)
#         desc = gsub("Replaced.*with ", "", desc)
#         desc = gsub("Modified priors.*estimates", "Higher K prior", desc)
        
#         # Truncate if still too long
#         if(nchar(desc) > 25) {
#             words = strsplit(desc, " ")[[1]]
#             desc = paste(words[1:min(3, length(words))], collapse=" ")
#         }
#         return(desc)
#     })
#     return(labels)
# }

# create short labels CSV for manual editing
# oneoff_models = model_desc_config[group != "Final ensemble"]
# oneoff_models[, short_label := create_short_labels(description)]

label_file = file.path(dir_output, "oneoff_short_labels_hand_edit.csv")
# if(!file.exists(label_file)) {
#     fwrite(oneoff_models[, .(short_name, group, description, short_label)], label_file)
#     cat("Created short labels file:", label_file, "\n")
#     cat("Edit manually if desired, then re-run script\n")
# }

# read back the labels
label_lookup = fread(label_file)
setkey(label_lookup, short_name)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# process groups function
process_group <- function(group_name, model_desc_config, label_lookup, model_stem, diag_dir, proj_dir) {
    cat("Processing group:", group_name, "\n")
    
    # get models in group
    group_models = model_desc_config[group == group_name]
    
    # create output directory
    safe_group_name = gsub("[^A-Za-z0-9_-]", "_", group_name)
    plot_dir = file.path(proj_dir, "plots", "one-off", safe_group_name)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # prepare model directories and labels
    model_dirs = c(diag_dir, file.path(model_stem, paste0(group_models$short_name, "*")))
    
    # expand wildcards and filter existing directories
    all_dirs = unlist(lapply(model_dirs, function(pattern) {
        if(grepl("\\*", pattern)) {
            matches = list.files(dirname(pattern), pattern = basename(gsub("\\*", ".*", pattern)), full.names = TRUE)
            return(matches[dir.exists(matches)])
        } else {
            return(ifelse(dir.exists(pattern), pattern, character(0)))
        }
    }))
    
    if(length(all_dirs) < 2) {
        cat("  Insufficient models found, skipping group\n")
        return(NULL)
    }
    
    # create labels
    group_lookup = label_lookup[short_name %in% group_models$short_name]
    group_labels = paste0(group_lookup$short_name, ": ", group_lookup$short_label)
    labels = c("Diagnostic", group_labels)
    labels = labels[1:length(all_dirs)]  # trim to match found directories
    
    # setup parameters
    custom_params = get_default_params()
    custom_params$fits$model_names = labels
    custom_params$fits$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$fits$active = TRUE  # TRUE | FALSE
    custom_params$fits$obs = TRUE  # TRUE | FALSE
    custom_params$fits$type = "Quantile"  # "Median" | "Spaghetti" | "Quantile"
    custom_params$fits$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$fits$resid = "PIT"  # "Ordinary" | "Standardized" | "PIT"
    custom_params$fits$ncol = 2
    custom_params$fits$resid_ncol = 1
    
    custom_params$ppp$model_names = labels
    custom_params$ppp$leading_params = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","x0","nu_catch")  # Any combination
    custom_params$ppp$raw = TRUE  # TRUE (transformed) | FALSE (raw)
    custom_params$ppp$show = "Both"  # "Prior" | "Posterior" | "Both"
    custom_params$ppp$combine = FALSE  # TRUE | FALSE
    custom_params$ppp$ncol = 3
    
    custom_params$ppts$model_names = labels
    custom_params$ppts$var = c("Depletion (D)", "Population (P)", "D_Dmsy", "F_Fmsy", "Removals", "Process error","Nominal CPUE","Catchability deviate")  # Any combination
    custom_params$ppts$show = "Posterior"  # "Prior" | "Posterior" | "Both"
    custom_params$ppts$combine = FALSE  # TRUE | FALSE
    custom_params$ppts$prop = 0.25  # 0.01 to 1.00 (increments of 0.05)
    custom_params$ppts$quants = 95  # 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    custom_params$ppts$ncol = 3
    
    # generate plots
    tryCatch({
        p1 = generate_index_fit(all_dirs, custom_params$fits)
        ggsave(file.path(plot_dir, "index_fit.png"), p1, width = 12, height = 8, dpi = 300)
        
        p2 = generate_catch_fit(all_dirs, custom_params$fits)
        ggsave(file.path(plot_dir, "catch_fit.png"), p2, width = 12, height = 6, dpi = 300)
        
        p3 = generate_ppp(all_dirs, custom_params$ppp)
        ggsave(file.path(plot_dir, "prior_posterior_params.png"), p3, width = 12, height = 10, dpi = 300)
        
        p4 = generate_ppts(all_dirs, custom_params$ppts)
        ggsave(file.path(plot_dir, "prior_posterior_timeseries.png"), p4, width = 12, height = 8, dpi = 300)
        
        cat("  Generated plots for", length(all_dirs), "models\n")
        return(group_name)
        
    }, error = function(e) {
        cat("  Error generating plots:", e$message, "\n")
        return(NULL)
    })
}

#________________________________________________________________________________________________________________________________________________________________________________________________________
# process each group (parallel or sequential)
groups = unique(model_desc_config[group != "Final ensemble"]$group)

# set number of cores (adjust as needed)
n_cores = min(parallel::detectCores() - 1, length(groups))

if(n_cores > 1) {
    cat("Processing", length(groups), "groups in parallel using", n_cores, "cores\n")
    
    # export necessary objects for parallel processing
    cluster_vars = list(
        model_desc_config = model_desc_config,
        label_lookup = label_lookup,
        model_stem = model_stem,
        diag_dir = diag_dir,
        proj_dir = proj_dir
    )
    
    results = parallel::mclapply(groups, function(group_name) {
        process_group(group_name, cluster_vars$model_desc_config, cluster_vars$label_lookup, 
                     cluster_vars$model_stem, cluster_vars$diag_dir, cluster_vars$proj_dir)
    }, mc.cores = n_cores)
    
} else {
    cat("Processing", length(groups), "groups sequentially\n")
    results = lapply(groups, process_group, model_desc_config, label_lookup, model_stem, diag_dir, proj_dir)
}

cat("Completed one-off sensitivity plots\n")
