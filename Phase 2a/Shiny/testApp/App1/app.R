library(shiny)

ui <- fluidPage(
  titlePanel(title = 'Test App 3'),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'var',label = '1. Select the variable', choices = c('Sepal.Length'=1, 'Sepdal.Width'=2,'Petal.Length'=3,'Petal.Width'=4),selected=3, selectize = T),
      sliderInput(inputId = 'bin', label = '2. Select the number of bins:', min = 5, max = 30, value = 15),
      radioButtons(inputId = 'color', label = '3. Choose the color of the historgram', choices = c('Blue','Green', 'black'), selected = 'black')
    ),
    mainPanel(
      tabsetPanel(
        type = 'tab',
        tabPanel('Histogram',
                 textOutput(outputId = 'text1'),
                 plotOutput(outputId = 'myhist')),
        tabPanel('Summary',
                 tableOutput(outputId = 'summary')),
        tabPanel('Data',
                 dataTableOutput(outputId = 'data')),
        tabPanel('Extra',
                 img(src='cat-pet-animal-domestic-104827.jpeg')),
        tabPanel('Extra2',
                 verbatimTextOutput(outputId = 'prints'))
      )
    )
  )
)


server <- function(input, output){
  output$text1 <- renderText({
    colm = as.numeric(input$var)
    paste('The variable that is analyzed is ', names(iris[colm]))
  })
  output$myhist <- renderPlot({
    colm = as.numeric(input$var)
    hist(iris[,colm],
         col = input$color,
         xlim = c(0, max(iris[,colm])),
         main = paste('Histogram of', names(iris[colm])),
         breaks = seq(0, max(iris[,colm]), l = input$bin + 1),
         xlab = names(iris[colm])
    )
  })
  output$data <- renderDataTable({
    iris
  })
  output$prints <- renderPrint({
    summary(iris)
  })
  output$summary <- renderTable({
    as.data.frame.matrix(summary(iris))
  }, na = '')
}


shinyApp(ui = ui, server = server)





