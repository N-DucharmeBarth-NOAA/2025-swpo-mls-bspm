server <- function(input, output, session) {

  # Helper function for null operator
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  
  # pixel height for each panel. i.e row height when plotting by species
  height_per_panel = 350

  ref_table_reduced = summary_dt %>%
                as.data.frame(.)

  output$summarytable = DT::renderDataTable({
    summary_df = summary_dt %>%
                 as.data.frame(.,stringsAsFactors=FALSE)
    summary_DT = DT::datatable(summary_df, filter = 'top',rownames=FALSE,
    options = list(scrollX = TRUE, search = list(regex = TRUE, caseInsensitive = FALSE),pageLength = 25))
    return(summary_DT)
  })
  outputOptions(output, "summarytable", suspendWhenHidden = FALSE)

  filtered_id = reactive({
    req(input$summarytable_rows_selected)
    keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$run_label)
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
        leading_params = input$hmc.leading_params %||% c("logK", "shape", "r"),
        raw = input$hmc.raw %||% TRUE,
        diag = input$hmc.diag %||% "None",
        eps = input$hmc.eps %||% TRUE,
        lags = as.numeric(input$hmc.lags %||% 30),
        scheme = input$hmc.scheme %||% "brightblue"
      )
    })

    get_ppc_index_params <- reactive({
      list(
        scheme = input$ppc_index.scheme %||% "brightblue",
        prop = as.numeric(input$ppc_index.prop %||% 0.25),
        active = input$ppc_index.active %||% TRUE,
        group = input$ppc_index.group %||% TRUE,
        stat = input$ppc_index.stat %||% "median",
        qqdist = input$ppc_index.qqdist %||% "uniform"
      )
    })

    get_ppc_catch_params <- reactive({
      list(
        scheme = input$ppc_catch.scheme %||% "brightblue",
        prop = as.numeric(input$ppc_catch.prop %||% 0.25),
        stat = input$ppc_catch.stat %||% "median",
        qqdist = input$ppc_catch.qqdist %||% "uniform"
      )
    })

    get_index_fit_params <- reactive({
      list(
        prop = as.numeric(input$index_fit.prop %||% 0.25),
        active = input$index_fit.active %||% TRUE,
        obs = input$index_fit.obs %||% TRUE,
        type = input$index_fit.type %||% "Median",
        quants = as.numeric(input$index_fit.quants %||% 95),
        resid = input$index_fit.resid %||% "PIT",
        ncol = as.numeric(input$index_fit.ncol %||% 2),
        resid_ncol = as.numeric(input$index_fit.resid_ncol %||% 1)
      )
    })

    get_catch_fit_params <- reactive({
      list(
        prop = as.numeric(input$catch_fit.prop %||% 0.25),
        obs = input$catch_fit.obs %||% TRUE,
        type = input$catch_fit.type %||% "Median",
        quants = as.numeric(input$catch_fit.quants %||% 95),
        resid = input$catch_fit.resid %||% "PIT"
      )
    })

  # Generic plot renderer - handles all common plotting logic
  render_plot <- function(plot_func, params_func = NULL, single_only = FALSE, 
                         default_params = NULL, output_name = "Plot") {
    renderPlot({
      req(filtered_id(), length(filtered_id()) > 0)
      
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
  # PPC INDEX PLOTS
  # =============================================================================
  
  output$plots_ppc_dens <- render_plot(
    generate_ppc_dens, get_ppc_index_params, single_only = TRUE, output_name = "PPC density overlay"
  )
  
  output$plots_ppc_ecdf <- render_plot(
    generate_ppc_ecdf, get_ppc_index_params, single_only = TRUE, output_name = "PPC ECDF"
  )
  
  output$plots_ppc_pit_ecdf <- render_plot(
    generate_ppc_pit_ecdf, get_ppc_index_params, single_only = TRUE, output_name = "PPC PIT ECDF"
  )
  
  output$plots_ppc_stat <- render_plot(
    generate_ppc_stat, get_ppc_index_params, single_only = TRUE, output_name = "PPC test statistics"
  )
  
  output$plots_ppc_loo_pit <- render_plot(
    generate_ppc_loo_pit, get_ppc_index_params, single_only = TRUE, output_name = "LOO PIT"
  )
  
  output$plots_ppc_loo_qq <- render_plot(
    generate_ppc_loo_qq, get_ppc_index_params, single_only = TRUE, output_name = "LOO Q-Q"
  )
  
  output$plots_ppc_loo_interval <- render_plot(
    generate_ppc_loo_interval, get_ppc_index_params, single_only = TRUE, output_name = "LOO intervals"
  )

  # =============================================================================
  # PPC CATCH PLOTS
  # =============================================================================
  
  output$plots_ppc_catch_dens <- render_plot(
    generate_catch_ppc_dens, get_ppc_catch_params, single_only = TRUE, output_name = "Catch PPC density overlay"
  )
  
  output$plots_ppc_catch_ecdf <- render_plot(
    generate_catch_ppc_ecdf, get_ppc_catch_params, single_only = TRUE, output_name = "Catch PPC ECDF"
  )
  
  output$plots_ppc_catch_pit_ecdf <- render_plot(
    generate_catch_ppc_pit_ecdf, get_ppc_catch_params, single_only = TRUE, output_name = "Catch PPC PIT ECDF"
  )
  
  output$plots_ppc_catch_stat <- render_plot(
    generate_catch_ppc_stat, get_ppc_catch_params, single_only = TRUE, output_name = "Catch PPC test statistics"
  )
  
  output$plots_ppc_catch_loo_pit <- render_plot(
    generate_catch_ppc_loo_pit, get_ppc_catch_params, single_only = TRUE, output_name = "Catch LOO PIT"
  )
  
  output$plots_ppc_catch_loo_qq <- render_plot(
    generate_catch_ppc_loo_qq, get_ppc_catch_params, single_only = TRUE, output_name = "Catch LOO Q-Q"
  )
  
  output$plots_ppc_catch_loo_interval <- render_plot(
    generate_catch_ppc_loo_interval, get_ppc_catch_params, single_only = TRUE, output_name = "Catch LOO intervals"
  )

  # =============================================================================
  # INDEX FIT PLOTS
  # =============================================================================
  
  output$plots_index_fit <- render_plot(
    generate_index_fit, get_index_fit_params, output_name = "Index fit"
  )
  
  output$plots_index_fit_ppd <- render_plot(
    generate_index_fit_ppd, get_index_fit_params, output_name = "Index fit with PPD"
  )
  
  output$plots_index_fit_residuals <- render_plot(
    generate_index_fit_residuals, get_index_fit_params, output_name = "Index fit residuals"
  )

  # =============================================================================
  # CATCH FIT PLOTS
  # =============================================================================
  
  output$plots_catch_fit <- render_plot(
    generate_catch_fit, get_catch_fit_params, output_name = "Catch fit"
  )
  
  output$plots_catch_fit_ppd <- render_plot(
    generate_catch_fit_ppd, get_catch_fit_params, output_name = "Catch fit with PPD"
  )
  
  output$plots_catch_fit_residuals <- render_plot(
    generate_catch_fit_residuals, get_catch_fit_params, output_name = "Catch fit residuals"
  )
}
