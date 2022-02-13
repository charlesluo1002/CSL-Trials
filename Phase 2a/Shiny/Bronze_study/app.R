library(shiny)
library(ggplot2)
library(reshape2)

# Define UI for application
ui <- fluidPage(
    
    # Application title
    titlePanel("CSL346 Bronze study simulation"),
    tabsetPanel(type = 'tab',
                # Sidebar
                tabPanel('Default simulation outputs with 10M runs',
                         sidebarLayout(
                             sidebarPanel(
                                 checkboxInput(inputId = 'p1',
                                               label = 'Check if you want P(HR<=k | reject H0) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p2',
                                               label = 'Check if you want P(HR>k | not reject H0) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p3',
                                               label = 'Check if you want P(reject H0 | HR<=k) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p4',
                                               label = 'Check if you want P(not reject H0 | HR>k) plotted',
                                               value = TRUE),
                                 submitButton('Apply Changes')
                             ),
                             
                             # main panel with several tabs
                             mainPanel(
                                 span('This panel provides the 10m-runs pre-run simulation outputs for the probability estimates of CSL346 bronze study.', style = 'font-size:16px; color:grey'),
                                 tabsetPanel(
                                     type = 'tab',
                                     tabPanel('Probability Estimation Plot',
                                              verbatimTextOutput(outputId = 'text1'),
                                              plotOutput(outputId = 'plot1'),
                                              span(textOutput(outputId = 'text2'), style = 'color:red')
                                     ),
                                     tabPanel('Probability Estimation table',
                                              tableOutput(outputId = 'table1')),
                                     tabPanel('Full report',
                                              tags$iframe(style = "height:560px;width:100%",
                                                          src = "Bronze study simulations 2.pdf"))
                                 )
                             )
                         )
                ),
                tabPanel('Customized Simulation with small runs',
                         sidebarLayout(
                             sidebarPanel(
                                 sliderInput(inputId ='runs','Choose the number of runs for the simulation', value = 2000, min=2000, max=20000,step=100),
                                 numericInput(inputId = 'sampleSize', label='Enter the sample size:', value=75),
                                 numericInput(inputId = 'SD', label = 'Enter the assumed SD: ', value = 0.75),
                                 numericInput(inputId = 'alpha', label = 'Enter the assumed alpha ', value = 0.1),
                                 checkboxInput(inputId = 'p5',
                                               label = 'Check if you want P(HR<=k | reject H0) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p6',
                                               label = 'Check if you want P(HR>k | not reject H0) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p7',
                                               label = 'Check if you want P(reject H0 | HR<=k) plotted',
                                               value = TRUE),
                                 checkboxInput(inputId = 'p8',
                                               label = 'Check if you want P(not reject H0 | HR>k) plotted',
                                               value = TRUE),
                                 submitButton('Apply Changes')
                             ),
                             mainPanel(
                                 span('A sample simulation of user chosen sample size, one-tail alpha and SD of the treatment effects on log ACR is presented below for the probability estimates of CSL346 bronze study.', style = 'font-size:16px; color:grey'),
                                 tabsetPanel(
                                     type = 'tab',
                                     tabPanel('Probability Estimation Plot',
                                              verbatimTextOutput(outputId = 'text3'),
                                              plotOutput(outputId = 'plot2'),
                                              span(textOutput(outputId = 'text4'), style = 'color:red')
                                     ),
                                     tabPanel('Probability Estimation table',
                                              tableOutput(outputId = 'table2'))
                                 )
                             )
                         )
                )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    if (!exists('quick_estimation_table')) {load('bronze_study_simulation_outputs.RData')}

    if (!exists('df_post')) {df_post <- read.csv('PosteriorSample.csv')}
    
    output$text1 <- renderText({
        paste(c('The following probabilities are plotted:',names(quick_estimation_table)[2:5][c(input$p1,input$p2,input$p3,input$p4)]),collapse = '\n')
    })
    
    output$plot1 <- renderPlot({
        table_long2 <- melt(quick_estimation_table[,c(T,input$p1,input$p2,input$p3,input$p4,F,F,F,F,F,F)], id='k')  # convert to long format
        ggplot(data=table_long2,
               aes(x=k, y=value, colour = variable)) +
            geom_line(size=1.5) +
            geom_vline(xintercept=c(0.7,0.75,0.8, 0.85), linetype="dotted") + 
            #geom_line(size=2, linetype = rep(1:4,each=length(k_sequence)))
            labs(title = 'probability estimates', x = 'Hazard ratio cutoff k', y = 'Probability', color = 'Legend Title\n') +
            scale_color_manual(values = c("blue", "green", 'red', 'yellow')[c(input$p1,input$p2,input$p3,input$p4)]) +
            theme(
                plot.title = element_text(size=12, face="bold.italic", hjust = 0.5),
                axis.title.x = element_text(size=10, face="bold"),
                axis.title.y = element_text(size=10, face="bold")
            )
    })
    
    output$text2 <- renderText({
        paste('Assurance of the study =',assurance)
    })
    
    output$table1 <- renderTable({
        quick_estimation_table[(quick_estimation_table$k*1000)%%50 ==0,1:5]
    }, digits=3)
    
    # panel 2
    output$text3 <- renderText({
        paste(c('The following probabilities are plotted:',names(quick_estimation_table)[2:5][c(input$p5,input$p6,input$p7,input$p8)]),collapse = '\n')
    })
    
    quick_simulation <- function(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=100, n=1, l=1000, k_sequence = (500:1000)/1000){
        if (is.null(S)) S = ceiling(2*(qnorm(alpha_as)+qnorm(1-power_as))^2*sd_as^2/(mu_as)^2)
        sigma0 <- sd_as*sqrt(2/S)
        regionA_cutoff = qt(alpha_as, 2*S-2)*sigma0
        estimation_table = data.frame(matrix(0, length(k_sequence), 11), stringsAsFactors=F)
        colnames(estimation_table) <-  colnames(estimation_table) <- c('k','P(HR<=k | reject H0)', 'P(HR>k | not reject H0)', 'P(reject H0 | HR<=k)', 'P(not reject H0 | HR>k)', 't1o1', 't0o0','t1','t0','o1','o0')
        estimation_table$k = k_sequence
        assurance = 0
        # simulate data
        for (i in 1:l){
            sample_row <- as.numeric(df_post[i,])
            for (j in 1:n){
                gamma0 <- rnorm(1, sample_row[5], sqrt(sample_row[4]))
                gamma0hat <- rnorm(m,gamma0, sigma0)
                theta0 <- sample_row[11] + sample_row[6] * gamma0 + rnorm(m, 0, sqrt(sample_row[3]))
                for (kk in 1:length(k_sequence)){
                    k = k_sequence[kk]
                    estimation_table$t1o1[kk] = estimation_table$t1o1[kk] + sum(gamma0hat <= regionA_cutoff & theta0 <= log(k))
                    estimation_table$t0o0[[kk]] = estimation_table$t0o0[kk] + sum(theta0 > log(k) & gamma0hat > regionA_cutoff)
                    o1 = sum(gamma0hat <= regionA_cutoff)
                    t1 =sum(theta0 <= log(k))
                    estimation_table$o1[kk] = estimation_table$o1[kk] + o1
                    estimation_table$o0[kk] = estimation_table$o0[kk] + m - o1
                    estimation_table$t1[kk] = estimation_table$t1[kk] + t1
                    estimation_table$t0[kk] = estimation_table$t0[kk] + m - t1
                }
            }
            #if (i%%50 == 0) cat(i,'rounds done.\n')
        }
        
        estimation_table$`P(HR<=k | reject H0)` = with(estimation_table, t1o1/o1)
        estimation_table$`P(HR>k | not reject H0)` = with(estimation_table,t0o0/o0)
        estimation_table$`P(reject H0 | HR<=k)` = with(estimation_table,t1o1/t1)
        estimation_table$`P(not reject H0 | HR>k)` = with(estimation_table,t0o0/t0)
        return(estimation_table)
    }
    
    output$plot2 <- renderPlot({
        set.seed(1)
        simple_estimation_table <- quick_simulation(df_post, mu_as=-0.315, sd_as=input$SD, alpha_as=input$alpha, power_as=0.9, S=input$sampleSize, m=100, n=1, l=input$runs/100, k_sequence = (500:1000)/1000)
        
        table_long <- melt(simple_estimation_table[,c(T,input$p5,input$p6,input$p7,input$p8,F,F,F,F,F,F)], id='k')  # convert to long format
        
        ggplot(data=table_long,
               aes(x=k, y=value, colour = variable)) +
            geom_line(size=1.5) +
            geom_vline(xintercept=c(0.7,0.75,0.8, 0.85), linetype="dotted") + 
            #geom_line(size=2, linetype = rep(1:4,each=length(k_sequence)))
            labs(title = 'probability estimates', x = 'Hazard ratio cutoff k', y = 'Probability', color = 'Legend Title\n') +
            scale_color_manual(values = c("blue", "green", 'red', 'yellow')[c(input$p5,input$p6,input$p7,input$p8)]) +
            theme(
                plot.title = element_text(size=12, face="bold.italic", hjust = 0.5),
                axis.title.x = element_text(size=10, face="bold"),
                axis.title.y = element_text(size=10, face="bold")
            )
    })
    
    output$text4 <- renderText({
        set.seed(1)
        paste('Assurance =', quick_simulation(df_post, mu_as=-0.315, sd_as=input$SD, alpha_as=input$alpha, power_as=0.9, S=input$sampleSize, m=100, n=1, l=input$runs/100, k_sequence = (500:1000)/1000)$o1[1]/input$runs)
    })
    
    output$table2 <- renderTable({
        set.seed(1)
        simple_estimation_table <- quick_simulation(df_post, mu_as=-0.315, sd_as=input$SD, alpha_as=0.1, power_as=0.9, S=input$sampleSize, m=100, n=1, l=input$runs/100, k_sequence = (500:1000)/1000)
        
        simple_estimation_table[(quick_estimation_table$k*1000)%%50 ==0,1:5]
    }, digits=3)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
