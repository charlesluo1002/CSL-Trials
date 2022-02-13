library(shiny)
library(ggplot2)
library(gridExtra)


ui <- fluidPage(
  titlePanel(title = 'test app 2'),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = 'sepallength',label = 'sepal length',min = 4.3,max = 7.9,value = 5.8, step = 0.1),
      sliderInput(inputId = 'sepalwidth',label = 'sepal width',min = 2,max = 4.4,value = 3, step = 0.1),
      sliderInput(inputId = 'petallength',label = 'petal length',min = 1,max = 6.9,value = c(1.6,5.1), step = 0.1),
      sliderInput(inputId = 'petalwidth',label = 'petal width',min = 0.1,max = 2.5,value = 1.3, step = 0.3, animate = animationOptions(interval = 3000, loop = TRUE))
    ),
    mainPanel(
      plotOutput(outputId = 'box')
    )
  )
)


server <- function(input, output){
  output$box <- renderPlot({
  s1 <- subset(iris, Sepal.Length >= 4.3 & Sepal.Length < input$sepallength, select = c(Sepal.Length, Species))
  s2 <- subset(iris, Sepal.Width >= 2 & Sepal.Width < input$sepalwidth, select = c(Sepal.Width, Species))
  s3 <- subset(iris, Petal.Length >= input$petallength[1] & Petal.Length < input$petallength[2], select = c(Petal.Length, Species))
  s4 <- subset(iris, Petal.Width >= 0.1 & Petal.Width < input$petalwidth, select = c(Petal.Width, Species))
  t1 <- ggplot(s1, aes(x = Species, y = Sepal.Length)) + geom_boxplot() + theme_bw()
  t2 <- ggplot(s2, aes(x = Species, y = Sepal.Width)) + geom_boxplot() + theme_bw() 
  t3 <- ggplot(s3, aes(x = Species, y = Petal.Length)) + geom_boxplot() + theme_bw() 
  t4 <- ggplot(s4, aes(x = Species, y = Petal.Width)) + geom_boxplot() + theme_bw()
  grid.arrange(t1,t2,t3,t4, nrow = 1)
  })
}

shinyApp(ui = ui, server = server)





