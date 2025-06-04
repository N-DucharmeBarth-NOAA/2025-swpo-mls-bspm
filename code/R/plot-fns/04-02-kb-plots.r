# Standalone Plot Functions for SSP Model Analysis
# Part 4B: Kobe Plot Function

# Note: This requires Parts 1-3 and 4A to be loaded first

# =============================================================================
# KOBE PLOT FUNCTION
# =============================================================================

#' Generate Kobe Plot (P/P_MSY vs F/F_MSY)
#' @param model_dirs Vector of model directories
#' @param params List of parameters (uses defaults if NULL)
generate_kb <- function(model_dirs, params = NULL) {
  if (is.null(params)) params <- get_default_params()$kbmj
  
  # Load data for all models
  all_data <- lapply(model_dirs, load_model_data)
  tmp_summary <- rbindlist(lapply(all_data, function(x) x$summary), fill = TRUE)
  
  # Parameter mapping for Kobe plot
  parameter_map <- cbind(
    c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Surplus production"),
    c("D", "P", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "removals", "dev", "raw_epsp", "surplus_production")
  )
  colnames(parameter_map) <- c("input", "grab")
  target_par <- c("P_Pmsy", "F_Fmsy")
  
  # Generate posterior and prior data
  posterior_dt_list <- list()
  prior_dt_list <- list()
  
  for (i in seq_along(all_data)) {
    if (params$show %in% c("Both", "Prior")) {
      tmp_samples <- ssp_prior_pushforward(
        ssp_summary = all_data[[i]]$summary,
        stan_data = all_data[[i]]$stan_data,
        settings = all_data[[i]]$settings
      )
      prior_dt_list[[i]] <- ssp_derived_quants_ts(
        ssp_summary = all_data[[i]]$summary,
        samples_dt = tmp_samples,
        stan_data = all_data[[i]]$stan_data,
        settings = all_data[[i]]$settings,
        sub_sample_prop = params$prop
      )
    }
    
    if (params$show %in% c("Both", "Posterior")) {
      posterior_dt_list[[i]] <- ssp_derived_quants_ts(
        ssp_summary = all_data[[i]]$summary,
        samples_dt = all_data[[i]]$samples,
        stan_data = all_data[[i]]$stan_data,
        settings = all_data[[i]]$settings,
        sub_sample_prop = params$prop
      )
    }
  }
  
  # Process prior data
  if (params$show %in% c("Both", "Prior")) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[, "input"]), "grab"]
    prior_dt <- rbindlist(prior_dt_list) %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[match(target_par, parameter_map[, "input"]), "grab"], 
                        labels = parameter_map[match(target_par, parameter_map[, "input"]), "input"])] %>%
      .[, type := "Prior"] %>%
      .[, type := factor(type, levels = c("Prior", "Posterior"))]
  }
  
  # Process posterior data
  if (params$show %in% c("Both", "Posterior")) {
    grab_parameters <- parameter_map[match(target_par, parameter_map[, "input"]), "grab"]
    posterior_dt <- rbindlist(posterior_dt_list) %>%
      .[name %in% grab_parameters] %>%
      .[, name := factor(name, levels = parameter_map[match(target_par, parameter_map[, "input"]), "grab"], 
                        labels = parameter_map[match(target_par, parameter_map[, "input"]), "input"])] %>%
      .[, type := "Posterior"] %>%
      .[, type := factor(type, levels = c("Prior", "Posterior"))]
  }
  
  # Combine data
  if (params$show == "Both") {
    plot_dt <- rbind(prior_dt, posterior_dt) %>%
      merge(., tmp_summary[, .(run_id, run_label)])
  } else if (params$show == "Posterior") {
    plot_dt <- posterior_dt %>%
      merge(., tmp_summary[, .(run_id, run_label)])
  } else {
    plot_dt <- prior_dt %>%
      merge(., tmp_summary[, .(run_id, run_label)])
  }
  
  # Handle combine option
  if (params$combine) {
    plot_dt <- plot_dt %>%
      .[type == "Posterior", run_label := "Combined posterior"] %>%
      .[type == "Prior", run_label := "Combined prior"] %>%
      .[, run_label := factor(run_label, levels = c("Combined prior", "Combined posterior"))]
  }
  
  if (nrow(plot_dt) == 0) {
    stop("No data available for plotting")
  }
  
  # Generate uncertainty contours if requested
  contour_dt <- NULL
  if (params$uncertainty) {
    contour_points_dt <- plot_dt[row == max(plot_dt$row)] %>%
      .[, .(type, run_label, name, iter, value)] %>%
      dcast(., type + run_label + iter ~ name) %>%
      .[, group_id := paste0(type, "-", run_label)]
    
    unique_id <- unique(contour_points_dt$group_id)
    contour_dt_list <- list()
    
    for (i in seq_along(unique_id)) {
      tmp_contour <- contour_points_dt[group_id == unique_id[i]]
      
      if (nrow(tmp_contour) > 10) {  # Need enough points for contours
        tryCatch({
          mv.kde <- MASS::kde2d(tmp_contour$P_Pmsy, tmp_contour$F_Fmsy, 
                               n = as.numeric(params$resolution),
                               lims = c(range(tmp_contour$P_Pmsy) * c(0.75, 1.25), 
                                       range(tmp_contour$F_Fmsy) * c(0.75, 1.25)))
          
          dx <- diff(mv.kde$x[1:2])
          dy <- diff(mv.kde$y[1:2])
          sz <- sort(mv.kde$z)
          c1 <- cumsum(sz) * dx * dy
          
          dimnames(mv.kde$z) <- list(mv.kde$x, mv.kde$y)
          dc <- reshape2::melt(mv.kde$z)
          dc$prob <- approx(sz, 1 - c1, dc$value)$y
          dc <- as.data.table(dc[, c("Var1", "Var2", "prob")]) %>%
            dcast(., Var1 ~ Var2)
          
          kd <- MASS::kde2d(tmp_contour$P_Pmsy, tmp_contour$F_Fmsy, n = as.numeric(params$resolution))
          kd$x <- dc$Var1
          kd$y <- as.numeric(colnames(dc)[-1])
          kd$z <- as.matrix(dc)[, -1]
          
          cntr <- contourLines(kd, levels = as.numeric(params$quants) / 100)
          
          if (length(cntr) > 0) {
            cntr_dt <- rbindlist(lapply(seq_along(cntr), function(j) {
              dt <- as.data.table(cntr[[j]]) %>%
                .[, cntr_id := j] %>%
                .[, .(cntr_id, x, y)] %>%
                setnames(., c("x", "y"), target_par) %>%
                .[, group_id := unique_id[i]]
              rbind(dt, dt[1])  # Close the contour
            }))
            
            contour_dt_list[[i]] <- merge(unique(tmp_contour[, .(type, run_label, group_id)]), cntr_dt, by = "group_id")
          }
        }, error = function(e) {
          warning(paste("Could not generate contours for", unique_id[i], ":", e$message))
        })
      }
    }
    
    if (length(contour_dt_list) > 0) {
      contour_dt <- rbindlist(contour_dt_list, fill = TRUE) %>%
        .[, plot_id := paste0(group_id, "-", cntr_id)] %>%
        .[, type := sapply(group_id, function(x) strsplit(x, "-")[[1]][1])] %>%
        .[, run_label := sapply(group_id, function(x) strsplit(x, "-")[[1]][2])] %>%
        .[, type := factor(type, levels = c("Prior", "Posterior"))]
      
      if (params$combine) {
        contour_dt <- contour_dt %>%
          .[, run_label := factor(run_label, levels = c("Combined prior", "Combined posterior"))]
      }
    }
  }
  
  # Summarize trajectory data
  plot_dt <- plot_dt %>%
    .[!is.na(value)] %>%
    .[, .(med = median(value)), by = .(run_label, type, name, row)] %>%
    dcast(., run_label + type + row ~ name)
  
  # Create Kobe plot
  p <- plot_dt %>%
    ggplot() +
    xlab(expression(P/P["MSY"])) +
    ylab(expression(F/F["MSY"])) +
    coord_fixed(ylim = c(0, 2.25), xlim = c(0, 2.25)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    # Kobe plot quadrants
    geom_polygon(data = data.table(P_Pmsy = c(0, 1, 1, 0, 0), F_Fmsy = c(0, 0, 1, 1, 0)), 
                aes(x = P_Pmsy, y = F_Fmsy), fill = "yellow", alpha = 0.2) +
    geom_polygon(data = data.table(P_Pmsy = c(0, 1, 1, 0, 0), F_Fmsy = c(1, 1, 5, 5, 1)), 
                aes(x = P_Pmsy, y = F_Fmsy), fill = "red", alpha = 0.2) +
    geom_polygon(data = data.table(P_Pmsy = c(1, 5, 5, 1, 1), F_Fmsy = c(1, 1, 5, 5, 1)), 
                aes(x = P_Pmsy, y = F_Fmsy), fill = "orange", alpha = 0.2) +
    geom_polygon(data = data.table(P_Pmsy = c(1, 5, 5, 1, 1), F_Fmsy = c(0, 0, 1, 1, 0)), 
                aes(x = P_Pmsy, y = F_Fmsy), fill = "green", alpha = 0.2) +
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 1, linewidth = 1.15, color = "black") +
    geom_vline(xintercept = 1, linewidth = 1.15, color = "black") +
    geom_path(aes(x = P_Pmsy, y = F_Fmsy, linetype = type, color = run_label), linewidth = 1.25)
  
  # Add uncertainty contours if available
  if (!is.null(contour_dt)) {
    p <- p + geom_polygon(data = contour_dt, 
                         aes(x = P_Pmsy, y = F_Fmsy, color = run_label, fill = run_label, 
                             group = plot_id, linetype = type), alpha = 0.1)
  }
  
  # Add start and end points
  p <- p +
    geom_point(data = plot_dt[row == 1], aes(x = P_Pmsy, y = F_Fmsy, color = run_label), 
              shape = 21, fill = "white", size = 3) +
    geom_point(data = plot_dt[row == max(plot_dt$row)], aes(x = P_Pmsy, y = F_Fmsy, fill = run_label), 
              color = "black", shape = 21, size = 3) +
    scale_linetype_manual("Distribution type", values = c("dotted", "solid"), drop = FALSE) +
    viridis::scale_color_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    viridis::scale_fill_viridis("Model run", begin = 0.1, end = 0.8, direction = -1, option = "H", discrete = TRUE, drop = FALSE) +
    get_ssp_theme()
  
  return(p)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Generate Kobe plot
# model_dirs <- c("./path/to/model1/", "./path/to/model2/")
# 
# # Use defaults
# p1 <- generate_kb(model_dirs)
# 
# # Custom parameters
# custom_kobe_params <- list(
#   show = "Posterior",
#   combine = FALSE,
#   prop = 0.5,
#   uncertainty = TRUE,
#   quants = 90,
#   resolution = 200
# )
# p2 <- generate_kb(model_dirs, custom_kobe_params)
