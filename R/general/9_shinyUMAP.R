#' ShinyUMAP
#'
#' Interactive shiny app to explore single-cell results
#' 
#' @param seurat.rds Seurat Object containing UMAP info
#' @param mymarkers, optional - FindAllMarkers results
#' @param prefix, optional - sample ID
#' @param gene.description, optional - if available will add a short gene description to the markers table. DF V1 = gene_name V2 = description
#' @param OUTPUT
#'
#' @author Vincent Hahaut
#' 
#' @return Open a shiny app
#' 
#' @export

shinyUMAP <- function(seurat.rds = NULL, mymarkers = NULL, prefix = "mySample", OUTPUT = "~/Desktop/screenshots/", gene.description = NULL){
  
  print("0. Load Packages")
  
  list.of.packages <- c("ggplot2", "tidyverse", "DT", "shiny", "Seurat", "jpg")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(tidyverse))
  #suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(shiny))

  mycolors <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#E41A1C","#377EB8","#4DAF4A",
                "#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
                "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                "#B3B3B3","#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2")
  
  dir.create("~/Desktop/screenshots/")
  
  print("1. Load Dataset - It may take some time")
  
  # Get UMAP Info
  if(is.character(seurat.rds)){
    seurat <- read_rds(seurat.rds)
  } else if(is(seurat.rds) == "Seurat"){
    seurat <- seurat.rds
  }
  
  embedding <- Seurat::Embeddings(Seurat::Reductions(seurat, "umap"))
  UMAP <- data_frame(UMAP_1 = embedding[,1], UMAP_2 = embedding[,2])
  markers <- row.names(seurat)
  
  if(is.character(mymarkers)){
    mymarkersTable <- read_tsv(mymarkers)
  } else {
    mymarkersTable <- mymarkers
  }
  
  print("2. Initialize the ShinyApp")
  
  ui <- shiny::fluidPage(
    
    # Application title
    shiny::titlePanel("Explore Markers"),
    
    # SideBar with:
    # Input for markers to select
    # 
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("marker", label = "Select your markers", choices = markers, multiple = F),
        shiny::actionButton("button", label = "Go!"),
        shiny::p("Save current UMAP graphic to Desktop"),
        shiny::plotOutput("clusterPlot"),
        shiny::plotOutput("idPlot")
      ),
      
      # Display colored UMAP
      shiny::mainPanel(align="center",
        shiny::titlePanel("Dimension Reduction"),
        shiny::plotOutput("umapPlot", width = "80%"),
        shiny::br(),
        shiny::titlePanel("Cluster Markers"),
        DT::dataTableOutput("markerTable")
      )
    )
  )
  
  # Draw interractively the UMAP
  server <- function(input, output) {
    
    output$umapPlot <- shiny::renderPlot({
      
      out_plot <- ggplot2::ggplot(UMAP, ggplot2::aes_string(x = "UMAP_1", y = "UMAP_2", color = Seurat::FetchData(seurat, vars = input$marker)[,1])) +
        ggplot2::geom_point(size = 0.8, alpha = 0.75) +
        scale_color_gradientn(colours = c("grey", RColorBrewer::brewer.pal(10, "Spectral"))) +
        theme_bw() +
        ggtitle(label = prefix, subtitle = paste0("UMAP - ", input$marker)) +
        theme(plot.title = element_text(color="black", size=14, face="bold"))
      
      return(out_plot)
      
    })
    
    # Only need to create clusterPlot / idPlot plots once:
    if(!exists("idPlot")){
      
      print("Create Cluster and ID plots")
      
      output$clusterPlot <- shiny::renderPlot({
      
        # Get UMAP - seurat_clusters
        clusterCenter <- UMAP %>% 
          mutate(seurat_clusters = as.character(seurat$seurat_clusters)) %>% 
          group_by(seurat_clusters) %>%
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        
        clusterPlot <- ggplot2::ggplot() +
          ggplot2::geom_point(data = UMAP, aes(x = UMAP_1, y = UMAP_2, color = seurat$seurat_clusters), size = 0.8, alpha = 0.75) +
          scale_color_manual(values = mycolors) +
          theme_bw() +
          theme(legend.position = "none",
                axis.title = element_blank()) +
          geom_text(data = clusterCenter, aes(x = UMAP_1, y = UMAP_2, label = seurat_clusters), color = "black", fontface = "bold")
        
        return(clusterPlot)
      })
      
      output$idPlot <- shiny::renderPlot({
        
        # Get UMAP - Samples
        idPlot <- ggplot2::ggplot() +
          ggplot2::geom_point(data = UMAP, aes(x = UMAP_1, y = UMAP_2, color = seurat$ID, fill = seurat$ID), shape = 21, size = 0.5, alpha = 0.1) +
          scale_color_manual(values = mycolors) +
          scale_fill_manual(values = mycolors) +
          theme_bw() +
          theme(legend.position = "bottom",
                axis.title = element_blank(),
                legend.title = element_blank()) +
          guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2, alpha = 1)))
      
        return(idPlot)
      })
   
    # Cluster Marker Table
    if(!is.null(gene.description)){
      output$markerTable <- DT::renderDT(mymarkersTable %>%
                                           left_join(gene.description, by = c("gene" =  "V1")) %>%
                                           dplyr::select(gene, avg_log2FC, cluster, p_val, p_val_adj, V2) %>% 
                                           mutate(cluster = as.character(cluster),
                                                  p_val = as.numeric(round(p_val, 6)),
                                                  p_val_adj = as.numeric(round(p_val, 6)),
                                                  avg_log2FC = as.numeric(round(avg_log2FC, 3))),
                                         class = "display nowrap compact",
                                         filter = "top",
                                         options = list(scrollX = TRUE,
                                                        search = list(regex = TRUE, caseInsensitive = FALSE)
                                         ))
      
    } else {
      output$markerTable <- DT::renderDT(mymarkersTable %>%
                                           select(gene, avg_log2FC, cluster, p_val, p_val_adj) %>% 
                                           mutate(cluster = as.character(cluster),
                                                  p_val = as.numeric(round(p_val, 6)),
                                                  p_val_adj = as.numeric(round(p_val, 6)),
                                                  avg_log2FC = as.numeric(round(avg_log2FC, 3))),
                                         class = "display nowrap compact",
                                         filter = "top",
                                         options = list(scrollX = TRUE,
                                                        search = list(regex = TRUE, caseInsensitive = FALSE)
                                         ))
    }
    
    print("Good to Go!")
    }
  
    # Save graphic if button pressed
    observeEvent(input$button, { 

      myplot <- ggplot2::ggplot(UMAP, ggplot2::aes_string(x = "UMAP_1", y = "UMAP_2", color = Seurat::FetchData(seurat, vars = input$marker)[,1])) +
        ggplot2::geom_point(size = 0.8, alpha = 0.75) +
        scale_color_viridis_c() +
        theme_bw() +
        ggtitle(label = prefix, subtitle = paste0("UMAP - ", input$marker)) +
        theme(plot.title = element_text(color="black", size=18), 
              axis.title = element_text(color="black", size=16),
              axis.text = element_text(color="black", size=14),
              legend.text = element_text(color="black", size=14),
              legend.title = element_blank())
      
      ggsave(plot = myplot, 
             dpi = 120,
             width = 8,
             height = 8,
             filename = paste0(OUTPUT, prefix, "-UMAP_", as.character(input$marker), ".jpg"))
      
      print(paste0("Saved Graphic '", input$marker, "' to ", OUTPUT))
    })
    
  }
  
  
  # Run the application
  shiny::shinyApp(ui = ui, server = server)
  
}
