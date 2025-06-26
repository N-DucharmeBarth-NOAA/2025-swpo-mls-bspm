# CSS styling 
css <- htmltools::HTML(
    "#summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
        transform:rotateX(180deg);
    }
    #summary_table > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
        transform:rotateX(180deg);
    }"
)

ui = dashboardPage(
  header = dashboardHeader(title="BSPM Model Analysis"),
  sidebar = dashboardSidebar(
    br(),
    br(),
    sidebarMenu(id="sidebarmenu",
      menuItem("Introduction", tabName="introduction"),
      menuItem("Summary table", tabName="table"),
      menuItem("Bayesian diags: Convergence", tabName="plots_hmc"),
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
      ) # End of plots.hmc tab
    ) # End of tabItems
  ) # End of dashboardBody
)
