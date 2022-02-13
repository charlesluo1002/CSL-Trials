library(shiny)
library(ggplot2)
library(gridExtra)


ui <- fluidPage(
  titlePanel(title = 'Test App'),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = 'n',label = 'Sample Size',min = 1,max = 100,value = 25),
      textInput(inputId = 'title', label = 'Title of the plot', 'Customized plot in Rshiny'),
      radioButtons(inputId = 'color', label = 'Choose color', choices = list('Blue','Green')),
      submitButton('Click to apply changes')
    ),
    mainPanel(
      plotOutput(outputId = 'box')
    )
  )
)


server <- function(input, output){
  
  output$box <- renderPlot({
    boxplot(rnorm(input$n), col = input$color, main = input$title, xlab = 'sample data')
  })
  
}


shinyApp(ui = ui, server = server)





