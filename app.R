
library(shiny)
library(golubEsets)
library(RColorBrewer)
library(pheatmap)
library(grid)

# -------------------------
# Daten laden
# -------------------------
data(Golub_Train)

x <- exprs(Golub_Train)

# Samples benennen
sample_labels <- paste(
  pData(Golub_Train)$Samples,
  pData(Golub_Train)$ALL.AML,
  sep = "_"
)
colnames(x) <- sample_labels

# Log2-Transformation (robust)
xLog <- log2(pmax(x, 1))

# Varianzberechnung einmal
geneVariance <- apply(xLog, 1, var)
sortedGenes <- names(sort(geneVariance, decreasing = TRUE))

# Patientenzahlen
num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  
  titlePanel("Heatmap of Patients and Genes"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Anzahl Patienten"),
      tags$ul(
        tags$li(paste("ALL:", num_ALL)),
        tags$li(paste("AML:", num_AML))
      ),
      
      sliderInput("numberOfGenes",
        "Number of Genes",
        min = 10,
        max = 100,
        value = 50
      ),
      
      selectInput("distMea",
        "Distance Measure",
        choices = c("euclidean", "maximum",
                    "manhattan", "canberra",
                    "binary", "minkowski")
      ),
      
      selectInput("clustMeth",
        "Clustering Method",
        choices = c("ward.D", "ward.D2",
                    "single", "complete",
                    "average", "mcquitty",
                    "median", "centroid")
      )
    ),
    
    mainPanel(
      plotOutput("heatmap", height = "900px")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output) {
  
  selectedGenes <- reactive({
    xLog[sortedGenes[1:input$numberOfGenes], , drop = FALSE]
  })
  
  output$heatmap <- renderPlot({
    
    mat <- selectedGenes()
    
    # pheatmap MUSS mit print() gerendert werden
    print(
      pheatmap(
        mat,
        clustering_distance_rows = input$distMea,
        clustering_distance_cols = input$distMea,
        clustering_method = input$clustMeth,
        color = colorRampPalette(brewer.pal(8, "Blues"))(25),
        fontsize = 8
      )
    )
  })
}

shinyApp(ui, server)
