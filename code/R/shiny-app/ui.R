# CSS styling 
css <- htmltools::HTML(
    "#summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
        transform:rotateX(180deg);
    }
    #summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
        transform:rotateX(180deg);
    }"
)

# Wrapper function for plot boxes to reduce repetition
plot_box <- function(title, plot_output_id, collapsed = TRUE, help_text = "Select only one model.") {
  box(title = title, 
      solidHeader = TRUE, 
      collapsible = TRUE, 
      collapsed = collapsed, 
      status = "primary", 
      width = 12,
      p(help_text),
      plotOutput(plot_output_id, height = "auto")
  )
}

ui = dashboardPage(
  header = dashboardHeader(title="BSPM Model Analysis"),
  sidebar = dashboardSidebar(
    br(),
    br(),
    sidebarMenu(id="sidebarmenu",
      menuItem("Introduction", tabName="introduction"),
      menuItem("Summary table", tabName="table"),
      menuItem("Bayesian diags: Convergence", tabName="plots_hmc"),
      menuItem("Bayesian diags: PPC Index", tabName = "plots_ppc_index"),
      menuItem("Bayesian diags: PPC Catch", tabName = "plots_ppc_catch"),
      menuItem("Model fits: Index", tabName = "plots_index_fit"),
      menuItem("Model fits: Catch", tabName = "plots_catch_fit"),
      menuItem("Pr. & Post: Params", tabName = "plots_ppp"),
      menuItem("Pr. & Post: Timeseries", tabName = "plots_ppts"),
      menuItem("Kobe & Majuro", tabName = "plots_kbmj"),
      menuItem("Forecasts", tabName = "plots_forecasts"),
      selected = "introduction"
    ),

    # Prior-Posterior Parameter Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_ppp'",
      awesomeCheckboxGroup(
        inputId = "ppp.leading_params",
        label = "Parameter", 
        choices = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","x0","sigmaf"),
        selected = c("logK", "r")
      ),
      switchInput(
        inputId = "ppp.raw",  
        label = "Transform parameter?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      radioButtons(
        inputId = "ppp.show",
        label = "Show distributions", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "ppp.combine",  
        label = "Combine models?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      numericInput(
        inputId = "ppp.ncol",
        label = "Plot columns",
        value = 3,
        min = 1,
        max = 6,
        step = 1
      )
    ),

    # Prior-Posterior Time Series Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_ppts'",
      awesomeCheckboxGroup(
        inputId = "ppts.var",
        label = "Time series variable", 
        choices = c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Process error (mult.)", "Surplus production", "Effort deviate", "Catchability deviate", "Nominal CPUE"),
        selected = c("Depletion (D)", "Population (P)")
      ),
      radioButtons(
        inputId = "ppts.show",
        label = "Show distributions", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "ppts.combine",  
        label = "Combine models?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "ppts.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "ppts.quants",  
        label = "Quantiles (%)",
        choices = c(1, 5, seq(from = 10, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      numericInput(
        inputId = "ppts.ncol",
        label = "Plot columns",
        value = 3,
        min = 1,
        max = 6,
        step = 1
      )
    ),

    # Kobe & Majuro Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_kbmj'",
      radioButtons(
        inputId = "kbmj.show",
        label = "Show distributions", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "kbmj.combine",  
        label = "Combine models?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "kbmj.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "kbmj.uncertainty",  
        label = "Show uncertainty contours?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "kbmj.quants",  
        label = "Quantiles (%)",
        choices = c(1, 5, seq(from = 10, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "kbmj.resolution",  
        label = "Contour resolution",
        choices = seq(from = 50, to = 500, by = 25),
        selected = "300",
        grid = TRUE
      )
    ),

    # Forecast Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_forecasts'",
      awesomeCheckboxGroup(
        inputId = "forecast.var",
        label = "Forecast variables", 
        choices = c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Surplus production"),
        selected = c("Depletion (D)", "Removals")
      ),
      switchInput(
        inputId = "forecast.combine",  
        label = "Combine models?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "forecast.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecast.quants",  
        label = "Quantiles (%)",
        choices = c(1, 5, seq(from = 10, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecast.nyears",  
        label = "Forecast years",
        choices = seq(from = 1, to = 20, by = 1),
        selected = "10",
        grid = TRUE
      ),
      switchInput(
        inputId = "forecast.resample_epsp",  
        label = "Resample process error?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "forecast.avg_year",  
        label = "Average years for forecast",
        choices = seq(from = 1, to = 10, by = 1),
        selected = "5",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecast.scalar",  
        label = "Forecast scalar",
        choices = seq(from = 0.1, to = 5.0, by = 0.1),
        selected = "1",
        grid = TRUE
      ),
      numericInput(
        inputId = "forecast.ncol",
        label = "Plot columns",
        value = 3,
        min = 1,
        max = 6,
        step = 1
      )
    ),

    # Only show these on the plotting tabs - not Introduction and Summary table tabs
    conditionalPanel(condition="input.sidebarmenu == 'plots_hmc'",
     awesomeCheckboxGroup(
        inputId = "hmc.leading_params",
        label = "Parameter", 
        choices = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","x0","sigmaf","lp__"),
        selected = c("logK", "r")
      ),
      switchInput(
        inputId = "hmc.raw",  
        label = "Transform parameter?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      radioButtons(
        inputId = "hmc.diag",
        label = "Error type", 
        choices = c("None", "Divergences", "Max. treedepth"),
        selected = "None"
      ),
      switchInput(
        inputId = "hmc.eps",  
        label = "Include devs?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "hmc.lags",  
        label = "ACF lags (choose a value)",
        choices = seq(from = 5, to = 50, by = 5),
        selected = "30",
        grid = TRUE
      ),
      selectInput(
        inputId = "hmc.scheme",
        label = "Select bayesplot color scheme", 
        choices = c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE"),
        selected = "brightblue",
        multiple = FALSE
      )
    ),

    # PPC Index Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_ppc_index'",
      selectInput(
        inputId = "ppc_index.scheme",
        label = "Select bayesplot color scheme", 
        choices = c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE"),
        selected = "brightblue",
        multiple = FALSE
      ),
      sliderTextInput(
        inputId = "ppc_index.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "ppc_index.active",  
        label = "Only fitted indices?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      switchInput(
        inputId = "ppc_index.group",  
        label = "Aggregate observations for PPC?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      awesomeCheckboxGroup(
        inputId = "ppc_index.stat",
        label = "PPC statistic\n(choose 1 or 2)", 
        choices = c("mean", "median", "sd", "mad"),
        selected = "median"
      ),
      radioButtons(
        inputId = "ppc_index.qqdist",
        label = "QQ distribution", 
        choices = c("uniform", "normal"),
        selected = "uniform"
      )
    ),

    # PPC Catch Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_ppc_catch'",
      selectInput(
        inputId = "ppc_catch.scheme",
        label = "Select bayesplot color scheme", 
        choices = c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE"),
        selected = "brightblue",
        multiple = FALSE
      ),
      sliderTextInput(
        inputId = "ppc_catch.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      awesomeCheckboxGroup(
        inputId = "ppc_catch.stat",
        label = "PPC statistic\n(choose 1 or 2)", 
        choices = c("mean", "median", "sd", "mad"),
        selected = "median"
      ),
      radioButtons(
        inputId = "ppc_catch.qqdist",
        label = "QQ distribution", 
        choices = c("uniform", "normal"),
        selected = "uniform"
      )
    ),

    # Index Fit Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_index_fit'",
      sliderTextInput(
        inputId = "index_fit.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "index_fit.active",  
        label = "Only fitted indices?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      switchInput(
        inputId = "index_fit.obs",  
        label = "Show observation error?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      radioButtons(
        inputId = "index_fit.type",
        label = "Plot type", 
        choices = c("Median", "Spaghetti", "Quantile"),
        selected = "Median"
      ),
      sliderTextInput(
        inputId = "index_fit.quants",  
        label = "Quantiles (%)",
        choices = c(1, 5, seq(from = 10, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      radioButtons(
        inputId = "index_fit.resid",
        label = "Residual type", 
        choices = c("Ordinary", "Standardized", "PIT"),
        selected = "PIT"
      ),
      numericInput(
        inputId = "index_fit.ncol",
        label = "Plot columns",
        value = 2,
        min = 1,
        max = 6,
        step = 1
      ),
      numericInput(
        inputId = "index_fit.resid_ncol",
        label = "Residual plot columns",
        value = 1,
        min = 1,
        max = 6,
        step = 1
      )
    ),

    # Catch Fit Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_catch_fit'",
      sliderTextInput(
        inputId = "catch_fit.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "catch_fit.obs",  
        label = "Show observation error?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      radioButtons(
        inputId = "catch_fit.type",
        label = "Plot type", 
        choices = c("Median", "Spaghetti", "Quantile"),
        selected = "Median"
      ),
      sliderTextInput(
        inputId = "catch_fit.quants",  
        label = "Quantiles (%)",
        choices = c(1, 5, seq(from = 10, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      radioButtons(
        inputId = "catch_fit.resid",
        label = "Residual type", 
        choices = c("Ordinary", "Standardized", "PIT"),
        selected = "PIT"
      )
    ),

    br(),
    br(),
    tags$footer(
      div(style="text-align:center",
        tags$p("version 0.0.1"),
        tags$p(paste("Copyright", format(Sys.time(),"%Y"), "NOAA Fisheries - PIFSC"))
      )
    )
  ), # End of sidebar

  body = dashboardBody(
    # Force sidebar-tab connection with JavaScript
    tags$script(HTML("
      $(document).ready(function(){
        $('.sidebar-menu a').on('click', function(){
          var tabName = $(this).attr('data-value');
          if(tabName){
            Shiny.setInputValue('sidebarmenu', tabName);
          }
        });
      });
    ")),
    tags$head(tags$style(HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}') )),
    tags$head(tags$style(css)),
    # Start of main tab stuff
    tabItems(
      # **** Introduction ****
      tabItem(tabName="introduction", h2("Introduction"),
        fluidRow(column(12, includeMarkdown(paste0("./introduction.md"))))
      ), # End of introduction tab

      # **** Summary table ****
      tabItem(tabName="table", h2("Summary table"),
        fluidRow(box(title="Model metrics", collapsed=FALSE, solidHeader=TRUE, collapsible=TRUE, status="primary", width=12,
         DT::dataTableOutput("summarytable")))
      ), # End of table tab

      # **** Bayesian diagnostics plots ****
      tabItem(tabName="plots_hmc", h2("Bayesian diagnostics: Convergence"),
        fluidRow(
          plot_box("Parcoord", "plots_hmc_parcoord", collapsed = FALSE),
          plot_box("Pairs", "plots_hmc_pairs"),
          plot_box("Trace", "plots_hmc_trace"),
          plot_box("Rhat", "plots_hmc_rhat"),
          plot_box("Effective sample size", "plots_hmc_neff"),
          plot_box("Autocorrelation plots", "plots_hmc_acf")
        )
      ), # End of plots_hmc tab

      # PPC Index Tab
      tabItem(tabName = "plots_ppc_index", 
        h2("Bayesian diagnostics: Posterior Predictive Checking (PPC) - Index"),
        fluidRow(
          plot_box("Density overlay", "plots_ppc_dens", collapsed = FALSE),
          plot_box("ECDF", "plots_ppc_ecdf"),
          plot_box("ECDF (PIT)", "plots_ppc_pit_ecdf"),
          plot_box("Test statistics", "plots_ppc_stat"),
          plot_box("LOO-PIT", "plots_ppc_loo_pit"),
          plot_box("LOO-PIT QQ", "plots_ppc_loo_qq"),
          plot_box("LOO Posterior Predicted Interval", "plots_ppc_loo_interval")
        )
      ), # End of plots_ppc_index tab

      # PPC Catch Tab
      tabItem(tabName = "plots_ppc_catch", 
        h2("Bayesian diagnostics: Posterior Predictive Checking (PPC) - Catch"),
        fluidRow(
          plot_box("Density overlay", "plots_ppc_catch_dens", collapsed = FALSE),
          plot_box("ECDF", "plots_ppc_catch_ecdf"),
          plot_box("ECDF (PIT)", "plots_ppc_catch_pit_ecdf"),
          plot_box("Test statistics", "plots_ppc_catch_stat"),
          plot_box("LOO-PIT", "plots_ppc_catch_loo_pit"),
          plot_box("LOO-PIT QQ", "plots_ppc_catch_loo_qq"),
          plot_box("LOO Posterior Predicted Interval", "plots_ppc_catch_loo_interval")
        )
      ), # End of plots_ppc_catch tab

      # Index Fit Tab
      tabItem(tabName = "plots_index_fit", 
        h2("Model Fits: Index Fit"),
        fluidRow(
          plot_box("Index fit", "plots_index_fit", collapsed = FALSE, help_text = "Select one or more models.\nObserved vs predicted index data."),
          plot_box("Index fit with PPD", "plots_index_fit_ppd", help_text = "Select one or more models.\nIndex fits using posterior predictive distributions."),
          plot_box("Index fit residuals", "plots_index_fit_residuals", help_text = "Select one or more models.\nResidual analysis for index data.")
        )
      ), # End of plots_index_fit tab

      # Catch Fit Tab
      tabItem(tabName = "plots_catch_fit", 
        h2("Model Fits: Catch Fit"),
        fluidRow(
          plot_box("Catch fit", "plots_catch_fit", collapsed = FALSE, help_text = "Select one or more models.\nObserved vs predicted catch data."),
          plot_box("Catch fit with PPD", "plots_catch_fit_ppd", help_text = "Select one or more models.\nCatch fits using posterior predictive distributions."),
          plot_box("Catch fit residuals", "plots_catch_fit_residuals", help_text = "Select one or more models.\nResidual analysis for catch data.")
        )
      ), # End of plots_catch_fit tab

      # Prior-Posterior Parameter Tab
      tabItem(tabName = "plots_ppp", 
        h2("Prior & Posterior: Leading parameters"),
        fluidRow(
          plot_box("Parameter prior-posterior distributions", "plots_ppp", collapsed = FALSE, help_text = "Select one or more models.\nCompare prior and posterior distributions for model parameters.")
        )
      ), # End of plots_ppp tab

      # Prior-Posterior Time Series Tab
      tabItem(tabName = "plots_ppts", 
        h2("Prior & Posterior: Time series"),
        fluidRow(
          plot_box("Time series prior-posterior comparison", "plots_ppts", collapsed = FALSE, help_text = "Select one or more models.\nCompare prior and posterior distributions for time series variables.")
        )
      ), # End of plots_ppts tab

      # Kobe & Majuro Tab
      tabItem(tabName = "plots_kbmj", 
        h2("Kobe & Majuro Plots"),
        fluidRow(
          plot_box("Kobe plot", "plots_kobe", collapsed = FALSE, help_text = "Select one or more models.\nKobe plot: P/P_MSY vs F/F_MSY for stock status assessment."),
          plot_box("Majuro plot", "plots_majuro", help_text = "Select one or more models.\nMajuro plot: Depletion vs F/F_MSY for alternative stock status visualization.")
        )
      ), # End of plots_kbmj tab

      # Forecasts Tab
      tabItem(tabName = "plots_forecasts", 
        h2("Catch-based forecasts"),
        fluidRow(
          plot_box("Forecast projections", "plots_forecast", collapsed = FALSE, help_text = "Select one or more models.\nModel forecasts showing projected future population dynamics and management metrics.")
        )
      ) # End of plots_forecasts tab

    ) # End of tabItems
  ) # End of dashboardBody
)
