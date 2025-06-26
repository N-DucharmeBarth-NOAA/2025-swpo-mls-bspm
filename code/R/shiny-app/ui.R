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
      selected = "introduction"
    ),
    # Only show these on the plotting tabs - not Introduction and Summary table tabs
    conditionalPanel(condition="input.sidebarmenu == 'plots_hmc'",
     awesomeCheckboxGroup(
        inputId = "hmc.leading_params",
        label = "Parameter", 
        choices = c("logK","r","sigmao_add","sigmap","shape","qeff","rho","sigma_qdev","lp__"),
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
      ) # End of plots_ppc_catch tab

    ) # End of tabItems
  ) # End of dashboardBody
)
