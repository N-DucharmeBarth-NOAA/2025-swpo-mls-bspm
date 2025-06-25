# ui.R - Complete user interface with full plotting functionality

# Create custom theme
my_theme <- create_theme(
  adminlte_color(
    light_blue = "#434C5E"
  ),
  adminlte_sidebar(
    width = "400px",
    dark_bg = "#D8DEE9",
    dark_hover_bg = "#81A1C1",
    dark_color = "#2E3440"
  ),
  adminlte_global(
    content_bg = "#FFF",
    box_bg = "#D8DEE9", 
    info_box_bg = "#D8DEE9"
  )
)

# Constants for UI elements
COLOR_SCHEMES <- c("blue", "brightblue", "gray", "darkgray", "green", "pink", "purple", "red", "teal", "yellow", "viridis", "viridisA", "viridisB", "viridisC", "viridisD", "viridisE")
HMC_PARAMETERS <- c("logK", "x0", "r", "sigmao_add", "sigmap", "shape", "qeff", "rho", "sigma_qdev", "sigmaf", "ll_q")
STAT_CHOICES <- c("mean", "median", "sd", "mad", "min", "max", "q25", "q75")
FIT_TYPES <- c("Median", "Mean", "Quantile")
RESID_TYPES <- c("PIT", "Standardized", "Pearson")

# Define UI
ui <- dashboardPage(
  
  # Header
  dashboardHeader(title = "BSPM Model Analysis - Complete"),
  
  # Sidebar
  dashboardSidebar(
    width = 400,
    sidebarMenu(
      id = "sidebarmenu",  # Add ID for conditional panels
      menuItem("Introduction", tabName = "introduction", icon = icon("info-circle")),
      menuItem("Model Summary", tabName = "model_summary", icon = icon("table")),
      
      # HMC Diagnostics
      menuItem("HMC Diagnostics", tabName = "hmc_plots", icon = icon("chart-line"),
               menuSubItem("Trace & Pairs", tabName = "hmc_main"),
               menuSubItem("Convergence", tabName = "hmc_convergence")),
      
      # Posterior Predictive Checks
      menuItem("Posterior Predictive", icon = icon("chart-bar"),
               menuSubItem("CPUE PPC", tabName = "cpue_ppc"),
               menuSubItem("Catch PPC", tabName = "catch_ppc")),
      
      # Model Fits
      menuItem("Model Fits", icon = icon("chart-area"),
               menuSubItem("Index Fits", tabName = "index_fits"),
               menuSubItem("Catch Fits", tabName = "catch_fits")),
      
      # Management & Forecasts
      menuItem("Management", tabName = "management_plots", icon = icon("fish")),
      

    ),
    
    br(),
    
    # Model selection
    h4("Model Selection"),
    checkboxGroupInput("selected_models",
                      "Select models to compare:",
                      choices = model_choices,
                      selected = if(length(model_choices) > 0) model_choices[1] else NULL),
    
    br(),
    
    # Global plot options
    h4("Plot Options"),
    numericInput("plot_width", "Plot Width (px):", value = 900, min = 400, max = 1400, step = 50),
    numericInput("plot_height", "Plot Height (px):", value = 700, min = 300, max = 1200, step = 50),
    
    # Context-specific parameters
    conditionalPanel(
      condition = "input.tabs == 'hmc_main' || input.tabs == 'hmc_convergence'",
      h4("HMC Parameters"),
      selectInput("hmc_params", "Parameters:",
                  choices = HMC_PARAMETERS,
                  selected = c("logK", "x0", "r"),
                  multiple = TRUE),
      selectInput("hmc_diag", "Diagnostics:",
                  choices = c("None", "Divergences", "Max. treedepth"),
                  selected = "None"),
      checkboxInput("hmc_raw", "Use transformed parameters", value = TRUE),
      checkboxInput("hmc_eps", "Include eps parameter", value = TRUE),
      numericInput("hmc_lags", "Autocorrelation lags:", value = 30, min = 5, max = 50, step = 5),
      selectInput("hmc_scheme", "Color scheme:", choices = COLOR_SCHEMES, selected = "brightblue")
    ),
    
    conditionalPanel(
      condition = "input.sidebarmenu == 'cpue_ppc' || input.sidebarmenu == 'catch_ppc'",
      h4("PPC Parameters"),
      numericInput("ppc_prop", "Posterior draws (%):", value = 25, min = 1, max = 100, step = 5),
      checkboxInput("ppc_active", "Show active draws", value = TRUE),
      checkboxInput("ppc_group", "Group by index", value = TRUE),
      selectInput("ppc_stat", "Test statistics:",
                  choices = STAT_CHOICES,
                  selected = "median",
                  multiple = TRUE),
      selectInput("ppc_qqdist", "QQ distribution:",
                  choices = c("uniform", "normal"),
                  selected = "uniform"),
      selectInput("ppc_scheme", "Color scheme:", choices = COLOR_SCHEMES, selected = "brightblue")
    ),
    
    conditionalPanel(
      condition = "input.sidebarmenu == 'index_fits' || input.sidebarmenu == 'catch_fits'",
      h4("Model Fit Parameters"),
      numericInput("fits_prop", "Posterior draws (%):", value = 25, min = 1, max = 100, step = 5),
      checkboxInput("fits_active", "Show active draws", value = TRUE),
      checkboxInput("fits_obs", "Show observations", value = TRUE),
      selectInput("fits_type", "Fit type:", choices = FIT_TYPES, selected = "Median"),
      numericInput("fits_quants", "Credible interval (%):", value = 95, min = 50, max = 99, step = 5),
      selectInput("fits_resid", "Residual type:", choices = RESID_TYPES, selected = "PIT")
    ),
    
    conditionalPanel(
      condition = "input.sidebarmenu == 'management_plots'",
      h4("Management Plot Parameters"),
      numericInput("mgmt_prop", "Posterior draws (%):", value = 25, min = 1, max = 100, step = 5),
      checkboxInput("mgmt_uncertainty", "Show uncertainty", value = TRUE),
      numericInput("mgmt_quants", "Credible interval (%):", value = 95, min = 50, max = 99, step = 5),
      numericInput("forecast_years", "Forecast years:", value = 5, min = 1, max = 20, step = 1),
      selectInput("forecast_type", "Forecast type:", 
                  choices = c("Catch", "MSY", "F"), selected = "Catch"),
      checkboxInput("mgmt_combine", "Combine plots", value = FALSE)
    ),

  ),
  
  # Body
  dashboardBody(
    use_theme(my_theme),
    
    tabItems(
      
      # Introduction
      tabItem(tabName = "introduction", 
        fluidRow(
          box(title = "BSPM Model Analysis - Complete Suite", status = "primary", solidHeader = TRUE, width = 12,
            includeMarkdown("introduction.md")
          )
        )
      ),
      
      # Model summary
      tabItem(tabName = "model_summary", 
        h2("Model Summary"),
        p("Click on table rows to select models for comparison. Selected models will be used in all plot tabs."),
        fluidRow(
          box(title = "Model Comparison Table", status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("summary_table")
          )
        )
      ),
      
      # HMC Main Diagnostics
      tabItem(tabName = "hmc_main", 
        h2("HMC Diagnostics - Main"),
        fluidRow(
          box(title = "Trace Plots", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 12,
            plotOutput("hmc_trace_plot", height = "auto")
          ),
          box(title = "Parallel Coordinates", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("hmc_parcoord_plot", height = "auto")
          ),
          box(title = "Pairs Plot", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("hmc_pairs_plot", height = "auto")
          ),
          box(title = "Autocorrelation", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            plotOutput("hmc_acf_plot", height = "auto")
          )
        )
      ),
      
      # HMC Convergence
      tabItem(tabName = "hmc_convergence", 
        h2("HMC Diagnostics - Convergence"),
        fluidRow(
          box(title = "R-hat Diagnostics", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("hmc_rhat_plot", height = "auto")
          ),
          box(title = "Effective Sample Size", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("hmc_neff_plot", height = "auto")
          )
        )
      ),
      
      # CPUE PPC
      tabItem(tabName = "cpue_ppc", 
        h2("CPUE Posterior Predictive Checks"),
        fluidRow(
          box(title = "Density Overlays", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("cpue_ppc_dens_plot", height = "auto")
          ),
          box(title = "Empirical CDF", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("cpue_ppc_ecdf_plot", height = "auto")
          ),
          box(title = "PIT ECDF", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("cpue_ppc_pit_ecdf_plot", height = "auto")
          ),
          box(title = "Test Statistics", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("cpue_ppc_stat_plot", height = "auto")
          ),
          box(title = "LOO PIT", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("cpue_ppc_loo_pit_plot", height = "auto")
          ),
          box(title = "LOO Q-Q", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("cpue_ppc_loo_qq_plot", height = "auto")
          ),
          box(title = "LOO Intervals", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("cpue_ppc_loo_interval_plot", height = "auto")
          )
        )
      ),
      
      # Catch PPC
      tabItem(tabName = "catch_ppc", 
        h2("Catch Posterior Predictive Checks"),
        fluidRow(
          box(title = "Density Overlays", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("catch_ppc_dens_plot", height = "auto")
          ),
          box(title = "Empirical CDF", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("catch_ppc_ecdf_plot", height = "auto")
          ),
          box(title = "PIT ECDF", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("catch_ppc_pit_ecdf_plot", height = "auto")
          ),
          box(title = "Test Statistics", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("catch_ppc_stat_plot", height = "auto")
          ),
          box(title = "LOO PIT", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("catch_ppc_loo_pit_plot", height = "auto")
          ),
          box(title = "LOO Q-Q", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("catch_ppc_loo_qq_plot", height = "auto")
          ),
          box(title = "LOO Intervals", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 4,
            plotOutput("catch_ppc_loo_interval_plot", height = "auto")
          )
        )
      ),
      
      # Index Fits
      tabItem(tabName = "index_fits", 
        h2("Index Model Fits"),
        fluidRow(
          box(title = "Index Fit", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("index_fit_plot", height = "auto")
          ),
          box(title = "Index Fit with PPD", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("index_fit_ppd_plot", height = "auto")
          ),
          box(title = "Index Residuals", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            plotOutput("index_residuals_plot", height = "auto")
          )
        )
      ),
      
      # Catch Fits
      tabItem(tabName = "catch_fits", 
        h2("Catch Model Fits"),
        fluidRow(
          box(title = "Catch Fit", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("catch_fit_plot", height = "auto")
          ),
          box(title = "Catch Fit with PPD", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("catch_fit_ppd_plot", height = "auto")
          ),
          box(title = "Catch Residuals", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            plotOutput("catch_residuals_plot", height = "auto")
          )
        )
      ),
      
      # Management plots
      tabItem(tabName = "management_plots", 
        h2("Management Plots"),
        fluidRow(
          box(title = "Kobe Plot", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("kobe_plot", height = "auto")
          ),
          box(title = "Majuro Plot", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary", width = 6,
            plotOutput("majuro_plot", height = "auto")
          ),
          box(title = "Forecasts", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 12,
            plotOutput("forecasts_plot", height = "auto")
          ),
          box(title = "Prior-Posterior Parameters", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("prior_posterior_params_plot", height = "auto")
          ),
          box(title = "Prior-Posterior Time Series", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary", width = 6,
            plotOutput("prior_posterior_ts_plot", height = "auto")
          )
        )
      ),

    )
  )
)
