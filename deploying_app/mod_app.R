library(shiny)
library(deSolve)
library(tidyverse)
library(gridExtra)

source('expose_shiny.R')

ui <- fluidPage(
  titlePanel("Treg dynamics: Incumbent model"),
  
  sidebarLayout(
    
    sidebarPanel(
      position = "left",
      h4("Rates to explore"),
      fluidRow(
        withMathJax(),
        column(4,
               numericInput(
                 "rho_D", label = helpText('$$\\alpha_D$$ Division of displaceable'), value = 0.0004, min = 0, max = 0.01, step = 0.0001
                 )),
        column(4,
               numericInput(
                 "rho_I", label = helpText('$$\\alpha_I$$ Division of incumbent'), value = 0.04, min = 0, max = 1, step = 0.01
                 )),
        column(4,
               numericInput(
                 "delta_D", label = helpText('$$\\delta_D$$ Loss of displaceable'), value = 0.027, min = 0, max = 0.5, step = 0.002
               )),
        column(4,
               numericInput(
                 "psi", label = helpText('$$\\psi$$ Influx'), value = 0.01, min = 0, max = 1, step = 0.01
               )),
        column(4,
               numericInput(
                 "alpha", label = helpText('$$\\mu_f$$ Thymic exit'), value = 0.87, min = 0, max = 1, step = 0.01
               )),
        column(4,
               numericInput(
                 "beta", label = helpText('$$\\mu_b$$ Backcirculation'), value = 0.01, min = 0, max = 1, step = 0.01
               ))
      ),
      
        fluidRow(
        column(3,
               numericInput(
                 "f1", label = helpText('fraction Khi in y1'), value = 0.5, min = 0, max = 1, step = 0.1
               )),
        column(3,
               numericInput(
                 "f2", label = helpText('fraction Khi in y2'), value = 0.5, min = 0, max = 1, step = 0.1
               )),
        column(3,
               numericInput(
                 "f3", label = helpText('fraction Khi in y3'), value = 0.5, min = 0, max = 1, step = 0.1
               )),
        column(3,
               numericInput(
                 "f4", label = helpText('fraction Khi in y4'), value = 0.5, min = 0, max = 1, step = 0.1
               ))
        )),
    
    mainPanel(
      # Output: plots ----
      fluidRow(
        column(6,plotOutput(outputId="distPlot1", width="600px",height="500px")),
        column(6,plotOutput(outputId="distPlot2", width="600px",height="800px"))
      )
    )
    
    
    )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  simulate <- reactive({
    init_cond <- c(y1=input$f1 * exp(9.7), y2= (1 - input$f1) * exp(9.7), y3=input$f2 * exp(14.5), y4= (1 - input$f2) * exp(14.5), 
                   y5=input$f3 * exp(9.25), y6= (1 - input$f3) * exp(9.25), y7=input$f4 * exp(11.45), y8= (1 - input$f4) * exp(1.45),
                   y9=0, y10=0, y11=0, y12=0)
    params <-  c(psi=input$psi, rho_D=input$rho_D, delta_D=input$delta_D + input$rho_D,
                 alpha=input$alpha, rho_I=input$rho_I, delta_I=input$rho_I, beta=input$beta)
    
    
    ## Predictions
    init_pred1 <- ode(y=init_cond, times=c(40, 45), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
    init_cond1 <- c(init_pred1[1] + init_pred1[9], init_pred1[2] + init_pred1[10], init_pred1[3] + init_pred1[11],
                    init_pred1[4] + init_pred1[12], init_pred1[5], init_pred1[6], init_pred1[7], init_pred1[8],
                    y9=0,y10=0,y11=0,y12=0)
    R_ode_pred1 <- data.frame(ode(y=init_cond1, times=c(45, ts_pred1), func=shm_chi, parms=params, ageatBMT=45)) %>%
      filter(time != 45) %>%
      mutate(time_seq = time,
             counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
             counts_per = y3 + y4 + y5 + y6 + y11 + y12,
             total_counts = counts_thy + counts_per,
             Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
             Nfd_per = (y11 + y12)/(counts_per * chivec4),
             donor_ki_thy = (y9)/(y9 + y10),
             donor_ki_per = (y11)/(y11 + y12),
             host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
             host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
             ageBMT_bin = 'agebin1') %>%
      select(time_seq, ageBMT_bin, counts_thy, counts_per, Nfd_thy, Nfd_per, donor_ki_thy, donor_ki_per, host_ki_thy, host_ki_per)
    
    init_pred2 <- ode(y=init_cond, times=c(40, 66), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
    init_cond2 <- c(init_pred2[1] + init_pred2[9], init_pred2[2] + init_pred2[10], init_pred2[3] + init_pred2[11],
                    init_pred2[4] + init_pred2[12], init_pred2[5], init_pred2[6], init_pred2[7], init_pred2[8],
                    y9=0,y10=0,y11=0,y12=0)
    
    R_ode_pred2 <- data.frame(ode(y=init_cond2,  times=c(66, ts_pred2), func=shm_chi, parms=params, ageatBMT=66)) %>%
      filter(time != 66) %>%
      mutate(time_seq = time,
             counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
             counts_per = y3 + y4 + y5 + y6 + y11 + y12,
             total_counts = counts_thy + counts_per,
             Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
             Nfd_per = (y11 + y12)/(counts_per * chivec4),
             donor_ki_thy = (y9)/(y9 + y10),
             donor_ki_per = (y11)/(y11 + y12),
             host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
             host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
             ageBMT_bin = 'agebin2') %>%
      select(time_seq, ageBMT_bin, counts_thy, counts_per, Nfd_thy, Nfd_per, donor_ki_thy, donor_ki_per, host_ki_thy, host_ki_per)
    
    
    init_pred3 <-  ode(y=init_cond, times=c(40, 76), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
    init_cond3 <- c(init_pred3[1] + init_pred3[9], init_pred3[2] + init_pred3[10], init_pred3[3] + init_pred3[11],
                    init_pred3[4] + init_pred3[12], init_pred3[5], init_pred3[6], init_pred3[7], init_pred3[8],
                    y9=0,y10=0,y11=0,y12=0)
    
    R_ode_pred3 <-data.frame(ode(y=init_cond3,  times=c(76, ts_pred3), func=shm_chi, parms=params, ageatBMT=76)) %>%
      filter(time != 76) %>%
      mutate(time_seq = time,
             counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
             counts_per = y3 + y4 + y5 + y6 + y11 + y12,
             total_counts = counts_thy + counts_per,
             Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
             Nfd_per = (y11 + y12)/(counts_per * chivec4),
             donor_ki_thy = (y9)/(y9 + y10),
             donor_ki_per = (y11)/(y11 + y12),
             host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
             host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
             ageBMT_bin = 'agebin3') %>%
      select(time_seq, ageBMT_bin, counts_thy, counts_per, Nfd_thy, Nfd_per, donor_ki_thy, donor_ki_per, host_ki_thy, host_ki_per)
    
    
    init_pred4 <-  ode(y=init_cond, times=c(40, 118), func=shm_chi, parms=params, ageatBMT=40)[2,2:13]
    init_cond4 <- c(init_pred4[1] + init_pred4[9], init_pred4[2] + init_pred4[10], init_pred4[3] + init_pred4[11],
                    init_pred4[4] + init_pred4[12], init_pred4[5], init_pred4[6], init_pred4[7], init_pred4[8],
                    y9=0,y10=0,y11=0,y12=0)
    
    R_ode_pred4 <- data.frame(ode(y=init_cond4,  times=c(118, ts_pred4), func=shm_chi, parms=params, ageatBMT=118)) %>%
      filter(time != 118) %>%
      mutate(time_seq = time,
             counts_thy = y1 + y2 + y7 + y8 + y9 + y10,
             counts_per = y3 + y4 + y5 + y6 + y11 + y12,
             total_counts = counts_thy + counts_per,
             Nfd_thy = (y9 + y10)/(counts_thy * chivec4),
             Nfd_per = (y11 + y12)/(counts_per * chivec4),
             donor_ki_thy = (y9)/(y9 + y10),
             donor_ki_per = (y11)/(y11 + y12),
             host_ki_thy = (y1 + y7)/(y1 + y2 + y7 + y8),
             host_ki_per = (y3 + y5)/(y3 + y4 + y5 + y6),
             ageBMT_bin = 'agebin4') %>%
      select(time_seq, ageBMT_bin, counts_thy, counts_per, Nfd_thy, Nfd_per, donor_ki_thy, donor_ki_per, host_ki_thy, host_ki_per)
    
    
    rbind(R_ode_pred1, R_ode_pred2, R_ode_pred3, R_ode_pred4)
    
  })
  
  output$distPlot1 <- renderPlot({
    R_pred <- simulate() 
    
    ## Total counts
    Counts_pred <- R_pred %>%
      select(time_seq, ageBMT_bin, contains('counts')) %>%
      rename(Thymus = counts_thy,
             Periphery = counts_per) %>%
      gather(c(Thymus, Periphery), key='location', value = 'total_counts')
    
    p1 <- ggplot() +
      geom_line(data = Counts_pred, aes(x = time_seq, y = total_counts, color = ageBMT_bin)) +
      geom_point(data = counts_data, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
      labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
      scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
      scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
      scale_y_continuous(limits = c(5e3, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
      facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
      guides(fill = 'none') + myTheme + theme(legend.position = c(0.85, 0.2))
    
    # Normalised donor fractions
    Nfd_pred <- R_pred %>%
      select(time_seq, ageBMT_bin, contains('Nfd')) %>%
      rename(Thymus = Nfd_thy,
             Periphery = Nfd_per) %>%
      gather(c(Thymus, Periphery), key='location', value = 'Nfd')
    
    p2 <- ggplot() +
      geom_line(data = Nfd_pred, aes(x = time_seq, y = Nfd, color = ageBMT_bin)) +
      geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Nfd, color = ageBMT_bin), size=2) +
      labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
      scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
      scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
      scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
      facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
      guides(fill='none', col = 'none')+ myTheme
    
    ptlist <- list(p1, p2)
    grid.arrange(grobs=ptlist, heights=,nrow=length(ptlist))
    
  })
  
  
  output$distPlot2 <- renderPlot({
    R_pred <- simulate()
    # Thymic Ki67 fractions
    ki_thy_pred <- R_pred %>%
      select(time_seq, ageBMT_bin, contains('ki_thy')) %>%
      rename(Donor = donor_ki_thy,
             Host = host_ki_thy) %>%
      gather(c(Donor, Host), key='subcomp', value = 'prop_ki')
    
    
    p3 <- ggplot() +
      geom_line(data = ki_thy_pred, aes(x = time_seq, y = prop_ki*100, color = subcomp)) +
      geom_point(data = filter(ki_data, location == "Thymus"), aes(x = age.at.S1K, y = prop_ki*100, color = subcomp), size=1.5) +
      labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in thymic naive Tregs") +
      scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
      scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
      facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
      guides(fill='none') + myTheme + theme(legend.title = element_blank(), legend.position = c(0.875, 0.925))
    
    # Peripheral Ki67 fractions
    ki_per_pred <- R_pred %>%
      select(time_seq, ageBMT_bin, contains('ki_per')) %>%
      rename(Donor = donor_ki_per,
             Host = host_ki_per) %>%
      gather(c(Donor, Host), key='subcomp', value = 'prop_ki')
    
    
    p4 <- ggplot() +
      geom_line(data = ki_per_pred, aes(x = time_seq, y = prop_ki*100, color = subcomp)) +
      geom_point(data = filter(ki_data, location == "Periphery"), aes(x = age.at.S1K, y = prop_ki*100, color = subcomp), size=1.5) +
      labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in peripheral naive Tregs") +
      scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
      scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
      facet_wrap(~ ageBMT_bin, scales = 'free', labeller = as_labeller(fac_labels))+
      guides(fill='none', col='none') + myTheme + theme(legend.title = element_blank(), legend.position = c(0.875, 0.925))
    
    ptlist <- list(p3, p4)
    grid.arrange(grobs=ptlist, heights=,nrow=length(ptlist))
  })
}

shinyApp(ui = ui, server = server)