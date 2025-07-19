# Nicholas Ducharme-Barth
# 2025/07/18
# Generate median and 95% credible intervals for leading parameters across all models

# Copyright (c) 2025 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define directories
proj_dir = this.path::this.proj()
dir_model_runs = file.path(proj_dir,"data","output","model_runs")
dir_output = file.path(proj_dir,"data","output")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load model configuration
model_desc_config = fread(file.path(dir_output,"model_desc_config.csv"))
model_desc_config[,short_name := sprintf("%04d", short_name)]

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define leading parameters from #tbl-priors (plus sigmaf for some models)
leading_params = c("logK", "r", "shape", "x0", "qeff", "rho", "sigma_qdev", "sigmap", "sigmao_add", "nu_catch", "sigmaf")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# function to calculate parameter summary statistics
calculate_parameter_summary = function(hmc_data, param_name) {
  if (!param_name %in% hmc_data$variable) {
    return(list(median = NA, lower = NA, upper = NA))
  }
  param_values = hmc_data[variable == param_name]$value
  list(
    median = median(param_values, na.rm = TRUE),
    lower = quantile(param_values, 0.025, na.rm = TRUE),
    upper = quantile(param_values, 0.975, na.rm = TRUE)
  )
}

#________________________________________________________________________________________________________________________________________________________________________________________________________
# process all models
results_list = list()

for (i in 1:nrow(model_desc_config)) {
  model_name = model_desc_config$short_name[i]
  
  # find model directory
  model_runs = list.files(dir_model_runs, pattern = paste0("^", model_name))
  if (length(model_runs) == 0) next
  
  # load hmc samples
  hmc_file = file.path(dir_model_runs, model_runs[1], "hmc_samples.csv")
  if (!file.exists(hmc_file)) next
  
  hmc_samples = fread(hmc_file)
  
  # calculate summaries for each parameter (initialize with NAs)
  model_results = data.table(
    model = model_name,
    parameter = leading_params,
    median = NA_real_,
    lower_ci = NA_real_,
    upper_ci = NA_real_
  )
  
  # only calculate stats for parameters that exist in this model
  available_params = intersect(leading_params, unique(hmc_samples$variable))
  
  for (param in available_params) {
    stats = calculate_parameter_summary(hmc_samples, param)
    model_results[parameter == param, `:=`(
      median = stats$median,
      lower_ci = stats$lower,
      upper_ci = stats$upper
    )]
  }
  
  results_list[[model_name]] = model_results
}

#________________________________________________________________________________________________________________________________________________________________________________________________________
# combine results and save
all_results = rbindlist(results_list)

# wide format for table display
results_wide = dcast(all_results, parameter ~ model, 
                    value.var = c("median", "lower_ci", "upper_ci"))

# long format for detailed analysis  
results_long = all_results

# save outputs with parameter availability summary
fwrite(results_wide, file.path(dir_output, "parameter_credible_intervals_wide.csv"))
fwrite(results_long, file.path(dir_output, "parameter_credible_intervals_long.csv"))

# create parameter availability matrix
param_availability = all_results[, .(available = !is.na(median)), by = .(model, parameter)]
availability_wide = dcast(param_availability, model ~ parameter, value.var = "available")
setcolorder(availability_wide, c("model", leading_params))
fwrite(availability_wide, file.path(dir_output, "parameter_availability_matrix.csv"))

# summary
cat("Processed", length(results_list), "models\n")
cat("Results saved to data/output/\n")
cat("Parameter availability matrix also saved\n")
