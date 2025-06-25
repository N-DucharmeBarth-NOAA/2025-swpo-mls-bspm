# ui.R - Complete user interface styled like the example

# CSS styling (source this from css.r in your app.R)
css <- htmltools::HTML(
    "#summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
        transform:rotateX(180deg);
    }
    #summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
        transform:rotateX(180deg);
    }"
)

# Define UI
ui <- dashboardPage(
  header = dashboardHeader(title = "BSPM Model Analysis"),
  
  sidebar = dashboardSidebar(
    br(),
    br(),
    sidebarMenu(id = "sidebarmenu",
      menuItem("Introduction", tabName = "introduction"),
      menuItem("Summary table", tabName = "table"),
      menuItem("Bayesian diags: Convergence", tabName = "plots_hmc"),
      menuItem("Bayesian diags: PPC", tabName = "plots_tab_ppc"),
      menuItem("Model fits", tabName = "plots_model_fits"),
      menuItem("Pr. & Post: params", tabName = "plots_tab_ppp"),
      menuItem("Pr. & Post: time-series", tabName = "plots_tab_ppts"),
      menuItem("Kobe & Majuro", tabName = "plots_tab_kbmj"),
      menuItem("Forecasts", tabName = "plots_tab_forecasts"),
      selected = "introduction"  # Add this to set default selection
    ),
    
    # HMC Diagnostics Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_hmc'",
      awesomeCheckboxGroup(
        inputId = "hmc.leading_params",
        label = "Parameter", 
        choices = c("logK", "x0", "r", "sigmao_add", "sigmap", "shape", "sigmaf", "ll_q", "lp__"),
        selected = c("logK", "x0", "r")
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
      awesomeRadio(
        inputId = "hmc.diag",
        label = "Error type", 
        choices = c("None", "Divergences", "Max. treedepth"),
        selected = "None"
      ),
      switchInput(
        inputId = "hmc.eps",  
        label = "Include eps?",
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
      pickerInput(
        inputId = "hmc.scheme",
        label = "Select bayesplot color scheme", 
        choices = c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE"),
        selected = "brightblue",
        multiple = FALSE
      )
    ),
    
    # PPC Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_tab_ppc'",
      pickerInput(
        inputId = "ppc.scheme",
        label = "Select bayesplot color scheme", 
        choices = c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE"),
        selected = "brightblue",
        multiple = FALSE
      ),
      sliderTextInput(
        inputId = "ppc.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "ppc.active",  
        label = "Only fitted indices?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      switchInput(
        inputId = "ppc.group",  
        label = "Aggregate observations for PPC?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      awesomeCheckboxGroup(
        inputId = "ppc.stat",
        label = "PPC statistic\n(choose 1 or 2)", 
        choices = c("mean", "median", "sd", "mad"),
        selected = "median"
      ),
      awesomeRadio(
        inputId = "ppc.qqdist",
        label = "QQ distribution", 
        choices = c("uniform", "normal"),
        selected = "uniform"
      )
    ),
    
    # Model Fits Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_model_fits'",
      sliderTextInput(
        inputId = "fits.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      switchInput(
        inputId = "fits.active",  
        label = "Only fitted indices?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      switchInput(
        inputId = "fits.obs",  
        label = "Show obs. error?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      awesomeRadio(
        inputId = "fits.type",
        label = "Show", 
        choices = c("Median", "Spaghetti", "Quantile"),
        selected = "Median"
      ),
      sliderTextInput(
        inputId = "fits.quants",  
        label = "Credible interval (%)",
        choices = c(1, seq(from = 5, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      awesomeRadio(
        inputId = "fits.resid",
        label = "Residual type", 
        choices = c("Ordinary", "Standardized", "PIT"),
        selected = "PIT"
      )
    ),
    
    # Prior-Posterior Parameters Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_tab_ppp'",
      awesomeCheckboxGroup(
        inputId = "ppp.leading_params",
        label = "Parameter", 
        choices = c("logK", "x0", "r", "sigmao_add", "sigmap", "shape", "sigmaf", "ll_q"),
        selected = c("logK", "x0", "r")
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
      awesomeRadio(
        inputId = "ppp.show",
        label = "Show", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "ppp.combine",  
        label = "Combine posterior?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      )
    ),
    
    # Prior-Posterior Time Series Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_tab_ppts'",
      pickerInput(
        inputId = "ppts.var",
        label = "Select metric(s)", 
        choices = c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Process error (mult.)", "Surplus production"),
        selected = c("Depletion (D)", "F_Fmsy", "Removals", "Process error (mult.)"),
        options = list(`actions-box` = TRUE), 
        multiple = TRUE
      ),
      awesomeRadio(
        inputId = "ppts.show",
        label = "Show", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "ppts.combine",  
        label = "Combine posterior?",
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
        label = "Credible interval (%)",
        choices = c(1, seq(from = 5, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      )
    ),
    
    # Kobe & Majuro Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_tab_kbmj'",
      awesomeRadio(
        inputId = "kbmj.show",
        label = "Show", 
        choices = c("Prior", "Posterior", "Both"),
        selected = "Both"
      ),
      switchInput(
        inputId = "kbmj.combine",  
        label = "Combine posterior?",
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
        label = "Show terminal year uncertainty?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "kbmj.quants",  
        label = "Credible interval (%)",
        choices = c(1, seq(from = 5, to = 95, by = 5), 99),
        selected = "95",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "kbmj.resolution",  
        label = "Contour resolution",
        choices = seq(from = 50, to = 500, by = 25),
        selected = "100",
        grid = TRUE
      )
    ),
    
    # Forecasts Controls
    conditionalPanel(condition = "input.sidebarmenu == 'plots_tab_forecasts'",
      pickerInput(
        inputId = "forecasts.var",
        label = "Select metric(s)", 
        choices = c("Depletion (D)", "Population (P)", "U", "F", "D_Dmsy", "P_Pmsy", "U_Umsy", "F_Fmsy", "Removals", "Process error", "Process error (raw)", "Surplus production"),
        selected = c("Depletion (D)", "F_Fmsy", "Removals", "Process error"),
        options = list(`actions-box` = TRUE), 
        multiple = TRUE
      ),
      switchInput(
        inputId = "forecasts.combine",  
        label = "Combine posterior?",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      sliderTextInput(
        inputId = "forecasts.prop",  
        label = "Sub-sample proportion",
        choices = c(0.01, seq(from = 0.05, to = 1, by = 0.05)),
        selected = "0.25",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecasts.quants",  
        label = "Credible interval (%)",
        choices = c(1, seq(from = 5, to = 100, by = 5)),
        selected = "95",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecasts.nyears",  
        label = "Years in forecast period",
        choices = 1:20,
        selected = "5",
        grid = TRUE
      ),
      switchInput(
        inputId = "forecasts.resample_epsp",  
        label = "Re-sample historical process?",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      awesomeRadio(
        inputId = "forecasts.type",
        label = "Forecast type", 
        choices = c("Catch", "U", "MSY", "Umsy"),
        selected = "Catch"
      ),
      sliderTextInput(
        inputId = "forecasts.avg_year",  
        label = "Average (Catch/U) for final n yrs.",
        choices = 1:10,
        selected = "3",
        grid = TRUE
      ),
      sliderTextInput(
        inputId = "forecasts.scalar",  
        label = "Catch/U multiplier",
        choices = c(0.01, seq(from = 0.1, to = 5, by = 0.1)),
        selected = "1",
        grid = TRUE
      )
    ),
    
    br(),
    br(),
    tags$footer(
      div(style = "text-align:center",
        tags$p("version 0.0.1"),
        tags$p(paste("Copyright", format(Sys.time(), "%Y"), "BSPM Analysis"))
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
    tags$head(tags$style(HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}'))),
    tags$head(tags$style(css)),
    
    # Start of main tab content
    tabItems(
      # Introduction Tab
      tabItem(tabName = "introduction", 
        h2("Introduction"),
        fluidRow(
          column(12, includeMarkdown("introduction.md"))
        )
      ),

      # Summary Table Tab
      tabItem(tabName = "table", 
        h2("Summary table"),
        fluidRow(
          box(title = "Model metrics", collapsed = FALSE, solidHeader = TRUE, collapsible = TRUE, status = "primary", width = 12,
            DT::dataTableOutput("summary_table")
          )
        )
      ),

      # HMC Diagnostics Tab
      tabItem(tabName = "plots_hmc", 
        h2("Bayesian diagnostics: Convergence"),
        fluidRow(
          box(title = "Parcoord", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_parcoord", height = "auto")
          ),
          box(title = "Pairs", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_pairs", height = "auto")
          ),
          box(title = "Trace", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_trace", height = "auto")
          ),
          box(title = "Rhat", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_rhat", height = "auto")
          ),
          box(title = "Effective sample size", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_neff", height = "auto")
          ),
          box(title = "Autocorrelation plots", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_hmc_acf", height = "auto")
          )
        )
      ),

      # PPC Tab
      tabItem(tabName = "plots_tab_ppc", 
        h2("Bayesian diagnostics: Posterior Predictive Checking (PPC)"),
        fluidRow(
          box(title = "Density overlay", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_dens", height = "auto")
          ),
          box(title = "ECDF", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_ecdf", height = "auto")
          ),
          box(title = "ECDF (PIT)", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_pit_ecdf", height = "auto")
          ),
          box(title = "Test statistics", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_stat", height = "auto")
          ),
          box(title = "LOO-PIT", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_loo_pit", height = "auto")
          ),
          box(title = "LOO-PIT QQ", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_loo_qq", height = "auto")
          ),
          box(title = "LOO Posterior Predicted Interval", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select only one model."),
            plotOutput("plots_ppc_loo_interval", height = "auto")
          )
        )
      ),

      # Model Fits Tab
      tabItem(tabName = "plots_model_fits", 
        h2("Model fits"),
        fluidRow(
          box(title = "Index fit", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_index_fit", height = "auto")
          ),
          box(title = "Index fit: posterior predicted", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_index_fit_ppd", height = "auto")
          ),
          box(title = "Index fit: residuals", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_index_fit_residuals", height = "auto")
          )
        )
      ),

      # Prior-Posterior Parameters Tab
      tabItem(tabName = "plots_tab_ppp", 
        h2("Prior & Posterior: Leading parameters"),
        fluidRow(
          box(title = "Distributions", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_ppp", height = "auto")
          )
        )
      ),

      # Prior-Posterior Time Series Tab
      tabItem(tabName = "plots_tab_ppts", 
        h2("Prior & Posterior: Time-series quantities"),
        fluidRow(
          box(title = "Time-series plots", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_ppts", height = "auto")
          )
        )
      ),

      # Kobe & Majuro Tab
      tabItem(tabName = "plots_tab_kbmj", 
        h2("Kobe & Majuro plots"),
        fluidRow(
          box(title = "Kobe plot", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_kb", height = "auto")
          ),
          box(title = "Majuro plot", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_mj", height = "auto")
          )
        )
      ),

      # Forecasts Tab
      tabItem(tabName = "plots_tab_forecasts", 
        h2("Forecasts"),
        fluidRow(
          box(title = "Forecast plots", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            p("Select at least one model."),
            plotOutput("plots_fcast", height = "auto")
          )
        )
      )
    ) # End of tabItems
  ) # End of dashboardBody
)
