# app.R - Main application entry point
# Source global.R first (for data and common functions)
source("./global.R")

# Source ui/server
source("./ui.R")
source("./server.R")

# Run the app
app <- shinyApp(ui = ui, server = server)
runApp(app, port = 8888, host = "127.0.0.1")
