# server.R - Complete server logic with full plotting functionality

server <- function(input, output, session) {
  
  # Reactive values for plot dimensions
  plot_width <- reactive({ input$plot_width })
  plot_height <- reactive({ input$plot_height })
  
  # Parameter builders - convert UI inputs to function parameters
  get_hmc_params <- reactive({
    list(
      leading_params = input$hmc_params %||% c("logK", "x0", "r"),
      raw = input$hmc_raw %||% TRUE,
      diag = input$hmc_diag %||% "None",
      eps = input$hmc_eps %||% TRUE,
      lags = input$hmc_lags %||% 30,
      scheme = input$hmc_scheme %||% "brightblue"
    )
  })
  
  get_ppc_params <- reactive({
    list(
      scheme = input$ppc_scheme %||% "brightblue",
      prop = (input$ppc_prop %||% 25) / 100,  # Convert percentage to decimal
      active = input$ppc_active %||% TRUE,
      group = input$ppc_group %||% TRUE,
      stat = input$ppc_stat %||% "median",
      qqdist = input$ppc_qqdist %||% "uniform"
    )
  })
  
  get_fits_params <- reactive({
    list(
      prop = (input$fits_prop %||% 25) / 100,  # Convert percentage to decimal
      active = input$fits_active %||% TRUE,
      obs = input$fits_obs %||% TRUE,
      type = input$fits_type %||% "Median",
      quants = input$fits_quants %||% 95,
      resid = input$fits_resid %||% "PIT",
      ncol = NULL, nrow = NULL, resid_ncol = NULL, resid_nrow = NULL
    )
  })
  
  get_management_params <- reactive({
    list(
      show = "Both",
      combine = input$mgmt_combine %||% FALSE,
      prop = (input$mgmt_prop %||% 25) / 100,  # Convert percentage to decimal
      uncertainty = input$mgmt_uncertainty %||% TRUE,
      quants = input$mgmt_quants %||% 95,
      resolution = 100
    )
  })
  
  get_forecast_params <- reactive({
    list(
      var = c("Depletion (D)", "F_Fmsy", "Removals", "Process error"),
      combine = input$mgmt_combine %||% FALSE,
      prop = (input$mgmt_prop %||% 25) / 100,
      quants = input$mgmt_quants %||% 95,
      nyears = input$forecast_years %||% 5,
      resample_epsp = TRUE,
      type = input$forecast_type %||% "Catch",
      avg_year = 3,
      scalar = 1,
      ncol = NULL, nrow = NULL
    )
  })
  
  # Generic plot renderer - handles all common plotting logic
  render_plot <- function(plot_func, params_func = NULL, single_only = FALSE, 
                         default_params = NULL, output_name = "Plot") {
    renderPlot({
      req(selected_model_ids(), length(selected_model_ids()) > 0)
      
      # Handle single-only plots
      if(single_only && length(selected_model_ids()) > 1) {
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = paste(output_name, "only available for single models.\nSelect only one row in the summary table."), 
                       hjust = 0.5, vjust = 0.5, size = 6) + 
               theme_void())
      }
      
      # Check if any models are selected
      if(length(selected_model_ids()) == 0) {
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                       label = "No models selected.\nClick on rows in the Model Summary table to select models.", 
                       hjust = 0.5, vjust = 0.5, size = 6) + 
               theme_void())
      }
      
      # Get parameters
      params <- if(!is.null(params_func)) params_func() else default_params
      
      # Get model path(s)
      if(length(selected_model_ids()) == 1) {
        model_path <- get_model_path(selected_model_ids()[1])
      } else {
        model_path <- selected_model_paths()
      }
      
      # Safe plot generation with error handling
      tryCatch({
        plot_func(model_path, params)
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message), 
                   hjust = 0.5, vjust = 0.5, size = 6, color = "red") +
          theme_void()
      })
    }, width = plot_width, height = plot_height)
  }
  
  # =============================================================================
  # MODEL SELECTION AND SUMMARY TABLE
  # =============================================================================
  
  # Create reference table for row selection
  ref_table_reduced <- summary_dt %>% as.data.frame(.)
  
  # Render summary table with row selection
  output$summary_table <- DT::renderDataTable({
    summary_df <- summary_dt %>% as.data.frame(., stringsAsFactors = FALSE)
    summary_DT <- DT::datatable(summary_df, 
                                filter = 'top',
                                rownames = FALSE,
                                selection = 'multiple',  # Allow multiple row selection
                                options = list(scrollX = TRUE, 
                                              pageLength = 25,
                                              search = list(regex = TRUE, caseInsensitive = FALSE)))
    return(summary_DT)
  })
  
  # Make sure table rendering isn't suspended when hidden
  outputOptions(output, "summary_table", suspendWhenHidden = FALSE)
  
  # Reactive to get selected model IDs from table row selection
  selected_model_ids <- reactive({
    req(input$summary_table_rows_selected)
    selected_rows <- input$summary_table_rows_selected
    model_ids <- ref_table_reduced[selected_rows, ]$model_id  # Adjust column name as needed
    return(model_ids)
  })
  
  # Get selected model paths
  selected_model_paths <- reactive({
    req(selected_model_ids())
    sapply(selected_model_ids(), get_model_path)
  })
  
  # =============================================================================
  # HMC DIAGNOSTICS - All single model only
  # =============================================================================
  
  output$hmc_trace_plot <- render_plot(
    generate_hmc_trace, get_hmc_params, single_only = TRUE, output_name = "Trace plots"
  )
  
  output$hmc_parcoord_plot <- render_plot(
    generate_hmc_parcoord, get_hmc_params, single_only = TRUE, output_name = "Parallel coordinates plots"
  )
  
  output$hmc_pairs_plot <- render_plot(
    generate_hmc_pairs, get_hmc_params, single_only = TRUE, output_name = "Pairs plots"
  )
  
  output$hmc_acf_plot <- render_plot(
    generate_hmc_acf, get_hmc_params, single_only = TRUE, output_name = "Autocorrelation plots"
  )
  
  output$hmc_rhat_plot <- render_plot(
    generate_hmc_rhat, get_hmc_params, single_only = TRUE, output_name = "R-hat plots"
  )
  
  output$hmc_neff_plot <- render_plot(
    generate_hmc_neff, get_hmc_params, single_only = TRUE, output_name = "N_eff plots"
  )
  
  # =============================================================================
  # CPUE POSTERIOR PREDICTIVE CHECKS
  # =============================================================================
  
  output$cpue_ppc_dens_plot <- render_plot(
    generate_ppc_dens, get_ppc_params, output_name = "CPUE PPC density"
  )
  
  output$cpue_ppc_ecdf_plot <- render_plot(
    generate_ppc_ecdf, get_ppc_params, output_name = "CPUE PPC ECDF"
  )
  
  output$cpue_ppc_pit_ecdf_plot <- render_plot(
    generate_ppc_pit_ecdf, get_ppc_params, output_name = "CPUE PPC PIT ECDF"
  )
  
  output$cpue_ppc_stat_plot <- render_plot(
    generate_ppc_stat, get_ppc_params, output_name = "CPUE PPC statistics"
  )
  
  output$cpue_ppc_loo_pit_plot <- render_plot(
    generate_ppc_loo_pit, get_ppc_params, single_only = TRUE, output_name = "CPUE LOO PIT"
  )
  
  output$cpue_ppc_loo_qq_plot <- render_plot(
    generate_ppc_loo_qq, get_ppc_params, single_only = TRUE, output_name = "CPUE LOO Q-Q"
  )
  
  output$cpue_ppc_loo_interval_plot <- render_plot(
    generate_ppc_loo_interval, get_ppc_params, single_only = TRUE, output_name = "CPUE LOO intervals"
  )
  
  # =============================================================================
  # CATCH POSTERIOR PREDICTIVE CHECKS
  # =============================================================================
  
  output$catch_ppc_dens_plot <- render_plot(
    generate_catch_ppc_dens, get_ppc_params, output_name = "Catch PPC density"
  )
  
  output$catch_ppc_ecdf_plot <- render_plot(
    generate_catch_ppc_ecdf, get_ppc_params, output_name = "Catch PPC ECDF"
  )
  
  output$catch_ppc_pit_ecdf_plot <- render_plot(
    generate_catch_ppc_pit_ecdf, get_ppc_params, output_name = "Catch PPC PIT ECDF"
  )
  
  output$catch_ppc_stat_plot <- render_plot(
    generate_catch_ppc_stat, get_ppc_params, output_name = "Catch PPC statistics"
  )
  
  output$catch_ppc_loo_pit_plot <- render_plot(
    generate_catch_ppc_loo_pit, get_ppc_params, single_only = TRUE, output_name = "Catch LOO PIT"
  )
  
  output$catch_ppc_loo_qq_plot <- render_plot(
    generate_catch_ppc_loo_qq, get_ppc_params, single_only = TRUE, output_name = "Catch LOO Q-Q"
  )
  
  output$catch_ppc_loo_interval_plot <- render_plot(
    generate_catch_ppc_loo_interval, get_ppc_params, single_only = TRUE, output_name = "Catch LOO intervals"
  )
  
  # =============================================================================
  # MODEL FITS
  # =============================================================================
  
  output$index_fit_plot <- render_plot(
    generate_index_fit, get_fits_params, output_name = "Index fit"
  )
  
  output$index_fit_ppd_plot <- render_plot(
    generate_index_fit_ppd, get_fits_params, output_name = "Index fit with PPD"
  )
  
  output$index_residuals_plot <- render_plot(
    generate_index_fit_residuals, get_fits_params, output_name = "Index residuals"
  )
  
  output$catch_fit_plot <- render_plot(
    generate_catch_fit, get_fits_params, output_name = "Catch fit"
  )
  
  output$catch_fit_ppd_plot <- render_plot(
    generate_catch_fit_ppd, get_fits_params, output_name = "Catch fit with PPD"
  )
  
  output$catch_residuals_plot <- render_plot(
    generate_catch_fit_residuals, get_fits_params, output_name = "Catch residuals"
  )
  
  # =============================================================================
  # MANAGEMENT PLOTS
  # =============================================================================
  
  output$kobe_plot <- render_plot(
    generate_kb, get_management_params, output_name = "Kobe plot"
  )
  
  output$majuro_plot <- render_plot(
    generate_mj, get_management_params, output_name = "Majuro plot"
  )
  
  output$forecasts_plot <- render_plot(
    generate_fcast, get_forecast_params, output_name = "Forecasts"
  )
  
  output$prior_posterior_params_plot <- render_plot(
    generate_ppp, default_params = get_default_params()$ppp, output_name = "Prior-posterior parameters"
  )
  
  output$prior_posterior_ts_plot <- render_plot(
    generate_ppts, default_params = get_default_params()$ppts, output_name = "Prior-posterior time series"
  )
  
  # =============================================================================
  # SESSION MANAGEMENT
  # =============================================================================
  
  # Track current tab for conditional panels
  observeEvent(input$sidebarmenu, {
    # This reactive helps conditional panels work correctly
    # The input$sidebarmenu value is used in conditionalPanel conditions
  })
  
  # Clean up resources on session end
  session$onSessionEnded(function() {
    # Cleanup any temporary files or resources if needed
  })
  
}
