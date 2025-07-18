# Function to calculate total number of data observations from stan_data.csv
# stan_data <- fread(stan_data_path)
ssp_count_total_observations <- function(stan_data) {
  
  # Load required packages
  if (!require(data.table)) {
    stop("data.table package is required")
  }
  
  # Filter for Data entries only
  data_entries <- stan_data[type == "Data"]
  
  # Get lambda values to identify active indices
  lambdas <- data_entries[name == "lambdas"]
  active_indices <- lambdas[value == 1, row]
  
  # Count index observations (non-missing values, -999 indicates missing)
  # Filter for active indices only
  active_index_obs <- data_entries[name == "index" & value != -999 & col %in% active_indices]
  n_index_obs <- nrow(active_index_obs)
  
  # Count removals/catch observations (subtract 1 because last obs is not fit)
  removals_obs <- data_entries[name == "obs_removals"]
  n_removals_obs <- nrow(removals_obs) - 1
  
  # Total observations
  total_obs <- n_index_obs + n_removals_obs
  
  return(total_obs)
}

# Example usage:
# total <- count_total_observations("path/to/stan_data.csv")
