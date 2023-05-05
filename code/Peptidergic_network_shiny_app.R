library(shiny)
library(visNetwork)
library(webshot2)

ui <- fluidPage(
  sliderInput("centralGravity", "centralGravity", min = 0.1, max = 10, value = 1),
  sliderInput("springConstant", "springConstant", min = 0, max = 0.1, value = 0.01),
  sliderInput("nodeDistance", "nodeDistance", min = 1, max = 500, value = 250),
  visNetworkOutput("network"),
  actionButton("captureBtn", "Capture Screenshot")
)

server <- function(input, output) {
  
  output$network <- renderVisNetwork({
    
    visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
      visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 23, type = "full") %>%
      visPhysics(
        solver = "repulsion", 
        repulsion = list(
          nodeDistance = input$nodeDistance, 
          centralGravity = input$centralGravity,
          springConstant = input$springConstant)
      ) %>%
      visEdges(
        smooth = list(style = "curve-style: haystack;"),
        scaling = list(min = 1, max = 15),
        color = list(inherit = TRUE, opacity = 0.7),
        arrows = list(to = list(
          enabled = TRUE,
          scaleFactor = 1.2, type = "arrow"
        ))
      ) %>%
      visNodes(
        borderWidth = 0.3,
        color = list(border = "black"),
        opacity = 0.9,
        shape = "dot",
        scaling = list(min = 15, max = 55),
        font = list(color = "black", size = 25),
      ) %>%
      visOptions(highlightNearest = TRUE, width = 2000, height = 1600) %>%
      visInteraction(
        navigationButtons = FALSE,
        dragNodes = TRUE, dragView = FALSE,
        zoomView = TRUE)
    
      
  })
  
  observeEvent(input$captureBtn, {
    webshot2::webshot(
      url = "http://127.0.0.1:4990/",  # Replace with your network visualization URL
      file = "pictures/peptidergic_networks_full_modules.png",  # Replace with your desired filename
      vwidth = 2000,  # Replace with your desired width
      vheight = 1600,  # Replace with your desired height
      cliprect = c(160, 120, 1750, 1480),  # Replace with the dimensions of your network visualization
      zoom = 2,  # Replace with your desired zoom level
      delay = 1  # Replace with your desired delay time (in seconds)
    )
  })
  
}

shinyApp(ui, server)

