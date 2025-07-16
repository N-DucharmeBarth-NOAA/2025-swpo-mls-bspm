ssp_compare_models_loo = function(model_dirs, 
                              model_runs_path = file.path(proj_dir, "data", "output", "model_runs"),
                              n_cores = parallelly::availableCores(omit = 1, logical = FALSE),
                              force_recalc = FALSE,
                              cache_loo = FALSE) {
    
    # Load required packages
    require(data.table)
    require(magrittr)
    require(rstan)
    require(loo)
    
    # Initialize list to store LOO objects
    loo_objects = list()
    
    # Process each directory
    for(i in 1:length(model_dirs)) {
        dir_name = model_dirs[i]
        dir_path = file.path(model_runs_path, dir_name)
        
        # Check if directory exists
        if(!dir.exists(dir_path)) {
            warning(paste("Directory does not exist:", dir_path))
            next
        }
        
        # Check if LOO has already been calculated (unless forcing recalculation)
        loo_obj_file = file.path(dir_path, "loo_object.rds")
        
        if(file.exists(loo_obj_file) && !force_recalc && cache_loo) {
            # Load existing LOO object
            cat("Loading existing LOO object for:", dir_name, "\n")
            loo_objects[[dir_name]] = readRDS(loo_obj_file)
            
        } else {
            # Calculate LOO
            cat("Calculating LOO for:", dir_name, "\n")
            
            # Check required files exist
            hmc_file = file.path(dir_path, "hmc_samples.csv")
            stan_data_file = file.path(dir_path, "stan_data.csv")
            
            if(!file.exists(hmc_file) || !file.exists(stan_data_file)) {
                warning(paste("Required files missing for:", dir_name))
                next
            }
            
            # Load data
            hmc_samples = fread(hmc_file)
            stan_data = fread(stan_data_file)
            
            # Calculate likelihood
            tmp_likelihood = ssp_calc_likelihood(hmc_samples, stan_data)
            
            # Prepare log-likelihood matrix
            log_lik_mat = tmp_likelihood %>% 
                          .[lambda > 0] %>%
                          na.omit(.) %>%
                          .[, value := value * lambda] %>%
                          .[, obs_id := paste0(T, "-", I)] %>%
                          .[, .(obs_id, iter, value)] %>%
                          dcast(., iter ~ obs_id) %>%
                          .[, iter := NULL] %>%
                          as.matrix(.)
            
            # Prepare likelihood matrix for relative efficiency
            lik_mat = tmp_likelihood %>% 
                          .[lambda > 0] %>%
                          na.omit(.) %>%
                          .[, value := exp(value * lambda)] %>%
                          .[, obs_id := paste0(T, "-", I)] %>%
                          .[, .(obs_id, iter, value)] %>%
                          dcast(., iter ~ obs_id) %>%
                          .[, iter := NULL] %>%
                          as.matrix(.)
            
            # Calculate relative efficiency
            r_eff = loo::relative_eff(lik_mat, 
                                     unique(hmc_samples[, .(iter, chain)])$chain, 
                                     cores = n_cores)
            
            # Calculate LOO
            loo_obj = loo::loo(log_lik_mat, r_eff = r_eff, cores = n_cores)
            
            # Store LOO object
            loo_objects[[dir_name]] = loo_obj
            
            # Save LOO object for future use (only if caching is enabled)
            if(cache_loo) {
                saveRDS(loo_obj, loo_obj_file)
            }
        }
    }
    
    # Remove any NULL entries (from failed calculations)
    loo_objects = loo_objects[!sapply(loo_objects, is.null)]
    
    # Check if we have at least 2 models to compare
    if(length(loo_objects) < 2) {
        stop("Need at least 2 successfully processed models for comparison")
    }
    
    # Perform LOO comparison
    cat("Comparing", length(loo_objects), "models...\n")
    comparison_result = loo::loo_compare(loo_objects)
    
    # Print summary
    cat("\nModel comparison results:\n")
    print(comparison_result)
    
    # Return both the comparison and individual LOO objects
    return(list(
        comparison = comparison_result,
        loo_objects = loo_objects,
        model_names = names(loo_objects)
    ))
}
