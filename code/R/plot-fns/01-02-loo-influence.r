#' Generate loo influence plot
generate_loo_influence <- function(model_dir, params = NULL) {
  if (is.null(params)) params <- get_default_params()$hmc
  
  # Load data
  data <- load_model_data(model_dir)
  
  # Check if rhat data is available
  if (is.null(data$loo)) {
    stop("Loo influence data not available. Ensure loo_dt.csv exists in the model directory.")
  }
  
  # Use loaded loo data
  plot_dt <- data$loo %>%
            .[,infl_cat := factor(infl_cat,levels=c("Low","High","Very high"))]
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Create loo influence plot
  p <- plot_dt %>% ggplot() + 
                    xlab("Time") +
                    ylab("Value") +
                    facet_wrap(~I,scales="free_y",nrow=length(unique(plot_dt$I))) +
                    geom_path(aes(x=T,y=value)) +
                    geom_point(aes(x=T,y=value,fill=infl_cat),shape=21,col="black",size=3) +
                    viridis::scale_color_viridis("LOO\nPareto-k:\nInfluence", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                    viridis::scale_fill_viridis("LOO\nPareto-k:\nInfluence", begin = 0.1, end = 0.8, direction = 1, option = "H", discrete = TRUE) +
                    get_ssp_theme()
  
  return(p)
}
