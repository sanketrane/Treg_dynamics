library(shiny)

ui <- fluidPage(
  titlePanel("Understanding naive Treg dynamics using the Incumbent model"),
  
  sidebarLayout(
    
    mainPanel(
      img(src = "Inc_p1.jpg", height = 352, width = 672),
      h4('Assumptions: Incumbent model'),
      p("
      1. Host compartment is heterogeneous – Displaceable and incumbent subsets.
      Donor compartment is homogeneous."),
      p("2. Displaceable in host and donor is maintained by:"),
      p("(i) influx from FoxP3- SP4, (ii) loss, and (iii) homeostatic division."),
      p("3. Host incumbent compartmental size is maintained by loss ≤ self-renewal."),
      p("4. All naïve Tregs circulate between thymus and periphery.")
    ),
    
    sidebarPanel(
      position = "left",
      h4("Parameters to explore"),
      fluidRow(
        column(6,
      sliderInput(
        "rho_D", label = "rho_D: rate of division of Displaceable cells", value = 0.0, min = 0, max = 1, step = 0.01
      )),
      column(6,
      sliderInput(
        "rho_I", label = "rho_I: rate of division of Incumbent cells", value = 0.0, min = 0, max = 1, step = 0.01
        ))
      ))
    )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#007bc2", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
}

shinyApp(ui = ui, server = server)