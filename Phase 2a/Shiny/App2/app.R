#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput(inputId ='runs','Choose the number of runs for the simulation', value = 10000, min=5000, max=20000,step=1000)
        ),
        mainPanel(
            span('This panel provides a small scale online simulation for the probability estimates of CSL346 bronze study.', style = 'font-size:16px; color:grey'),
            tabsetPanel(
                type = 'tab',
                tabPanel('Probability Plot',
                         0
                ),
                tabPanel('Probability Table',
                         0
                )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
}

# Run the application 
shinyApp(ui = ui, server = server)
