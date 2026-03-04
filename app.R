
library(shiny)
library(golubEsets)
library(RColorBrewer)
library(pheatmap)

options(shiny.error = browser)   # <--- DEBUG-MODUS: zeigt alle Fehler

# -------------------------
# Daten laden und vorbereiten
# -------------------------
data(Golub_Train)

x <- exprs(Golub_Train)

sample_labels <- paste(
  pData(Golub_Train)$Samples,
  pData(Golub_Train)$ALL.AML,
  sep = "_"
)

colnames(x) <- sample_labels

xLog <- log2(pmax(x, 1))

num_ALL <- sum(pData(Golub_Train)$ALL.AML == "ALL")
num_AML <- sum(pData(Golub_Train)$ALL.AML == "AML")

geneVariance <- apply(xLog, 1, var)
sortedGenes <- names(sort(geneVariance, decreasing = TRUE))

annotation <- data.frame(
  Leukemia = pData(Golub_Train)$ALL.AML
)
rownames(annotation) <- colnames(xLog)

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
                  value = 50),
      
      selectInput("distMea",
                  "Distance Measure",
                  choices = c("euclidean", "maximum",
                              "manhattan", "canberra",
                              "binary", "minkowski")),
      
      selectInput("clustMeth",
                  "Clustering Method",
                  choices = c("ward.D", "ward.D2",
                              "single", "complete",
                              "average", "mcquitty",
                              "median", "centroid"))
    ),
    
    mainPanel(
      plotOutput("heatmap", height = "900px"),
      verbatimTextOutput("debug")   # <-- DEBUG-FENSTER
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
  
  output$debug <- renderPrint({
    list(
      selected_genes_dim = dim(selectedGenes()),
      colnames_match = identical(colnames(selectedGenes()), rownames(annotation)),
      annotation_rows = rownames(annotation)[1:5],
      matrix_cols = colnames(selectedGenes())[1:5]
    )
  })
  
  output$heatmap <- renderPlot({
    
    mat <- selectedGenes()
    
    print(
      pheatmap(
        mat,
        clustering_distance_rows = input$distMea,
        clustering_distance_cols = input$distMea,
        clustering_method = input$clustMeth,
        annotation_col = annotation,
        color = colorRampPalette(brewer.pal(8, "Blues"))(25),
        fontsize = 9,
        main = "Heatmap der Top-variablen Gene"
      )
    )
  })
}

shinyApp(ui, server)
