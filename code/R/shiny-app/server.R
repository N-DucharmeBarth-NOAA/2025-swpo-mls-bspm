# server.R - Complete server logic matching the styled UI

server <- function(input, output, session) {
  
  # Helper function for null operator
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  
  ref_table_reduced = summary_dt %>% as.data.frame(.)

  # One reactive that does exactly what's needed
  filtered_id = reactive({
    req(input$summary_table_rows_selected)  # Correct table name
    keep_models = c(ref_table_reduced[input$summary_table_rows_selected, ]$run_label)
    return(keep_models)  
  })

  # Get model directories for selected models
  selected_model_dirs <- reactive({
    req(filtered_id())
    file.path("data","output","model_runs", filtered_id())
  })
  
  # Parameter builders - convert UI inputs to function parameters
  get_hmc_params <- reactive({
    list(
      leading_params = input$hmc.leading_params %||% c("logK", "x0", "r"),
      raw = input$hmc.raw %||% TRUE,
      diag = input$hmc.diag %||% "None",
      eps = input$hmc.eps %||% TRUE,
      lags = as.numeric(input$hmc.lags %||% 30),
      scheme = input$hmc.scheme %||% "brightblue"
    )
  })
  
  get_ppc_params <- reactive({
    list(
      scheme = input$ppc.scheme %||% "brightblue",
      prop = as.numeric(input$ppc.prop %||% 0.25),
      active = input$ppc.active %||% TRUE,
      group = input$ppc.group %||% TRUE,
      stat = input$ppc.stat %||% "median",
      qqdist = input$ppc.qqdist %||% "uniform"
    )
  })
  
  get_fits_params <- reactive({
    list(
      prop = as.numeric(input$fits.prop %||% 0.25),
      active = input$fits.active %||% TRUE,
      obs = input$fits.obs %||% TRUE,
      type = input$fits.type %||% "Median",
      quants = as.numeric(input$fits.quants %||% 95),
      resid = input$fits.resid %||% "PIT"
    )
  })
  
  get_ppp_params <- reactive({
    list(
      leading_params = input$ppp.leading_params %||% c("logK", "x0", "r"),
      raw = input$ppp.raw %||% TRUE,
      show = input$ppp.show %||% "Both",
      combine = input$ppp.combine %||% FALSE
    )
  })
  
  get_ppts_params <- reactive({
    list(
      var = input$ppts.var %||% c("Depletion (D)", "F_Fmsy", "Removals", "Process error (mult.)"),
      show = input$ppts.show %||% "Both",
      combine = input$ppts.combine %||% FALSE,
      prop = as.numeric(input$ppts.prop %||% 0.25),
      quants = as.numeric(input$ppts.quants %||% 95)
    )
  })
  
  get_kbmj_params <- reactive({
    list(
      show = input$kbmj.show %||% "Both",
      combine = input$kbmj.combine %||% FALSE,
      prop = as.numeric(input$kbmj.prop %||% 0.25),
      uncertainty = input$kbmj.uncertainty %||% TRUE,
      quants = as.numeric(input$kbmj.quants %||% 95),
      resolution = as.numeric(input$kbmj.resolution %||% 100)
    )
  })
  
  get_forecast_params <- reactive({
    list(
      var = input$forecasts.var %||% c("Depletion (D)", "F_Fmsy", "Removals", "Process error"),
      combine = input$forecasts.combine %||% FALSE,
      prop = as.numeric(input$forecasts.prop %||% 0.25),
      quants = as.numeric(input$forecasts.quants %||% 95),
      nyears = as.numeric(input$forecasts.nyears %||% 5),
      resample_epsp = input$forecasts.resample_epsp %||% TRUE,
      type = input$forecasts.type %||% "Catch",
      avg_year = as.numeric(input$forecasts.avg_year %||% 3),
      scalar = as.numeric(input$forecasts.scalar %||% 1)
    )
  })
  
  # Generic plot renderer - handles all common plotting logic
  render_plot <- function(plot_func, params_func = NULL, single_only = FALSE, 
                         default_params = NULL, output_name = "Plot") {
    renderPlot({
      # Check if models are selected
      if (is.null(filtered_id()) || length(filtered_id()) == 0) {
        return(ggplot() + 
              annotate("text", x = 0.5, y = 0.5, 
                      label = "Please select one or more models from the summary table",
                      hjust = 0.5, vjust = 0.5, size = 5) +
              theme_void())
      }
      
      # Handle single-only plots
      if(single_only && length(filtered_id()) > 1) {
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = paste("Select only one model for", output_name),
                       hjust = 0.5, vjust = 0.5, size = 5) +
               theme_void() +
               labs(title = paste("Multiple Models Selected -", output_name, "Unavailable")))
      }
      
      # Get parameters
      if (!is.null(params_func)) {
        params <- params_func()
      } else if (!is.null(default_params)) {
        params <- default_params
      } else {
        params <- list()
      }
      
      # Generate plot
      tryCatch({
        plot_func(selected_model_dirs(), params)
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                  label = paste("Error generating", output_name, ":\n", e$message),
                  hjust = 0.5, vjust = 0.5, size = 4, color = "red") +
          theme_void() +
          labs(title = paste("Error:", output_name))
      })
    }, height = 700, width = 900)  # Default dimensions
  }
  
  # =============================================================================
  # MODEL SUMMARY TABLE 
  # =============================================================================
  
  output$summary_table <- DT::renderDataTable({
    summary_df <- summary_dt %>%
                 as.data.frame(., stringsAsFactors = FALSE)
    summary_DT <- DT::datatable(summary_df, filter = 'top', rownames = FALSE,
    options = list(scrollX = TRUE, search = list(regex = TRUE, caseInsensitive = FALSE), pageLength = 25))
    return(summary_DT)
  })
  outputOptions(output, "summary_table", suspendWhenHidden = FALSE)
  
  # =============================================================================
  # HMC DIAGNOSTICS
  # =============================================================================
  
  output$plots_hmc_parcoord <- render_plot(
    generate_hmc_parcoord, get_hmc_params, single_only = TRUE, output_name = "Parallel coordinates"
  )
  
  output$plots_hmc_pairs <- render_plot(
    generate_hmc_pairs, get_hmc_params, single_only = TRUE, output_name = "Pairs plots"
  )
  
  output$plots_hmc_trace <- render_plot(
    generate_hmc_trace, get_hmc_params, single_only = TRUE, output_name = "Trace plots"
  )
  
  output$plots_hmc_rhat <- render_plot(
    generate_hmc_rhat, get_hmc_params, single_only = TRUE, output_name = "R-hat plots"
  )
  
  output$plots_hmc_neff <- render_plot(
    generate_hmc_neff, get_hmc_params, single_only = TRUE, output_name = "N_eff plots"
  )
  
  output$plots_hmc_acf <- render_plot(
    generate_hmc_acf, get_hmc_params, single_only = TRUE, output_name = "Autocorrelation plots"
  )
  
  # =============================================================================
  # PPC PLOTS
  # =============================================================================
  
  output$plots_ppc_dens <- render_plot(
    generate_ppc_dens, get_ppc_params, single_only = TRUE, output_name = "PPC density overlay"
  )
  
  output$plots_ppc_ecdf <- render_plot(
    generate_ppc_ecdf, get_ppc_params, single_only = TRUE, output_name = "PPC ECDF"
  )
  
  output$plots_ppc_pit_ecdf <- render_plot(
    generate_ppc_pit_ecdf, get_ppc_params, single_only = TRUE, output_name = "PPC PIT ECDF"
  )
  
  output$plots_ppc_stat <- render_plot(
    generate_ppc_stat, get_ppc_params, single_only = TRUE, output_name = "PPC test statistics"
  )
  
  output$plots_ppc_loo_pit <- render_plot(
    generate_ppc_loo_pit, get_ppc_params, single_only = TRUE, output_name = "LOO PIT"
  )
  
  output$plots_ppc_loo_qq <- render_plot(
    generate_ppc_loo_qq, get_ppc_params, single_only = TRUE, output_name = "LOO Q-Q"
  )
  
  output$plots_ppc_loo_interval <- render_plot(
    generate_ppc_loo_interval, get_ppc_params, single_only = TRUE, output_name = "LOO intervals"
  )
  
  # =============================================================================
  # MODEL FITS
  # =============================================================================
  
  output$plots_index_fit <- render_plot(
    generate_index_fit, get_fits_params, output_name = "Index fit"
  )
  
  output$plots_index_fit_ppd <- render_plot(
    generate_index_fit_ppd, get_fits_params, output_name = "Index fit with PPD"
  )
  
  output$plots_index_fit_residuals <- render_plot(
    generate_index_fit_residuals, get_fits_params, output_name = "Index residuals"
  )
  
  # =============================================================================
  # PRIOR-POSTERIOR PLOTS
  # =============================================================================
  
  output$plots_ppp <- render_plot(
    generate_ppp, get_ppp_params, output_name = "Prior-posterior parameters"
  )
  
  output$plots_ppts <- render_plot(
    generate_ppts, get_ppts_params, output_name = "Prior-posterior time series"
  )
  
  # =============================================================================
  # KOBE & MAJURO PLOTS
  # =============================================================================
  
  output$plots_kb <- render_plot(
    generate_kb, get_kbmj_params, output_name = "Kobe plot"
  )
  
  output$plots_mj <- render_plot(
    generate_mj, get_kbmj_params, output_name = "Majuro plot"
  )
  
  # =============================================================================
  # FORECASTS
  # =============================================================================
  
  output$plots_fcast <- render_plot(
    generate_fcast, get_forecast_params, output_name = "Forecasts"
  )
  
  # =============================================================================
  # SESSION MANAGEMENT
  # =============================================================================
  
  # Track current tab for debugging
  observeEvent(input$sidebarmenu, {
    cat("Current tab:", input$sidebarmenu, "\n")
  })
  
}
