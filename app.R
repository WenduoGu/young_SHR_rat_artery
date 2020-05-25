# Load packages ----
library(shiny)
library(ggplot2)
library(Seurat)
library(cowplot)
library(imager)
library(esquisse)
library(shinyWidgets)
library(DT)
library(plotly)

# Load data
load("data/Shiny_DFs.RData")


load("data/aggr.combined_right_celltype.Rdata")
load("data/aggr.combined_SMC.RData")
load("data/aggr.combined_EC.RData")
load("data/aggr.combined_MSC.RData")
load("data/aggr.combined_Immune.RData")
empty1 <- load.image("data/empty.png")
Figure_2A <- load.image("data/Fig2A.png")
Figure_3A <- load.image("data/Figure 3A SMC DEGs.png")
Figure_S7A <- load.image("data/Figure S7A EC DEGs.png")
Figure_5A <- load.image("data/Figure 5A MSC DEGs.png")


### Get plots
plot1 <- DimPlot(aggr.combined, label = TRUE, reduction = "tsne", pt.size = 0.4, label.size = 5, repel = TRUE)+
  scale_color_manual(values = c("#35978f", "#33a02c", "#f46d43", "#d53e4f", "#5e4fa2"))

Idents(aggr.combined_SMC) <- aggr.combined_SMC$clusters2
plot_SMC_dim <- DimPlot(aggr.combined_SMC, label = TRUE, reduction = "tsne", pt.size = 0.7, label.size = 5, repel = TRUE)+
  scale_color_manual(values = c("#c364c5" , "#cb4154" ,  "#3bb08f", "#ff6e4a","#1974d2" ,"#a5694f","#ffcf48" , "#f78fa7"  ,"#80daeb"  ))

Idents(aggr.combined_EC) <- aggr.combined_EC$clusters2
plot_EC_dim <- DimPlot(aggr.combined_EC, label = TRUE, reduction = "tsne", pt.size = 2, label.size = 5, repel = TRUE)+
  scale_color_manual(values = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d"  ))

Idents(aggr.combined_MSC) <- aggr.combined_MSC$clusters2
plot_MSC_dim <- DimPlot(aggr.combined_MSC, label = TRUE, reduction = "tsne", pt.size = 2, label.size = 5, repel = TRUE) #+
#scale_color_manual(values = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d"  ))

Idents(aggr.combined_Immune) <- aggr.combined_Immune$clusters2
plot_Immune_dim <- DimPlot(aggr.combined_Immune, label = TRUE, reduction = "tsne", pt.size = 2.3, label.size = 5, repel = TRUE)+
  scale_color_manual(values = c(rgb(138/255,168/255,218/255) , rgb(226/255,170/255,70/255), rgb(159/255,165/255,131/255) ,rgb(135/255,219/255,221/255) , 
                                rgb(242/255,142/255,144/255), rgb(186/255,112/255,110/255),
                                rgb(230/255,173/255,241/255) ,   rgb(138/255,238/255,178/255)))

genes <- rownames(aggr.combined)
Idents(aggr.combined_SMC) <- aggr.combined_SMC$stim2
Idents(aggr.combined_EC) <- aggr.combined_EC$stim2
Idents(aggr.combined_MSC) <- aggr.combined_MSC$stim2
Idents(aggr.combined_Immune) <- aggr.combined_Immune$stim2


library(shiny)
ui <- navbarPage(title = "Hypertension in young SHR rats", 
                 
                 #### page 1
                 tabPanel("FeaturePlot with major cell types",
                          # titlePanel("FeaturePlot"),

                          sidebarLayout(
                            sidebarPanel(
                              helpText("Select genes"),
                              selectInput("gene",
                                          label = "Choose a gene (symbol)",
                                          choices = genes),

                              #  filterDF_UI("filtering"),

                              width = 3,
                              plotOutput("empty1"),
                              plotOutput("Figure_2A"),
                              filterDF_UI("filtering")
                            ),

                            mainPanel(
                              tags$h2("Plot with all cells from major cell types"),

                              plotOutput("featureplot_overview"),
                              plotOutput("featureplot_splitview"),
                              radioButtons(
                                inputId = "dataset",
                                label = "Cell type markers and their intersection with gene classifications:",
                                choices =  c("All_markers_logFC_0.1",
                                             "Cell_type_markers_logFC_0.41",
                                             "Cell_type_markers_logFC_0.41_genes",
                                             "Cell_type_markers_logFC_0.41_major_class",
                                             "Cell_type_markers_logFC_0.41_minor_class",
                                             "Cell_type_markers_logFC_0.41_intersected",
                                             "Cell_type_markers_logFC_0.41_intersected_major_class",
                                             "Cell_type_markers_logFC_0.41_intersected_minor_class"),
                                inline = TRUE),

                              DT::dataTableOutput("table_all")
                            )
                          )),
                 ### page 2
                 tabPanel("SMC DEGs",
                          titlePanel("Differentially expressed genes in SMCs in young SHR rats"),
                          
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Select genes"),
                              selectInput("genes1", 
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Myh11", "Pcp4"),multiple = TRUE),
                              
                              #  filterDF_UI("filtering"),
                              
                              width = 3,
                              plotOutput("empty_SMC_page2"),
                              plotOutput("Figure_3A"),
                              
                              filterDF_UI("efiltering_SMC_DEGs_page2")
                              
                            ),
                            
                            mainPanel(
                              tags$h2("SMCs DEGs in SHR MA and AA"),
                              
                              plotOutput("vlnplot_dotplot_SMC_page2"),
                              
                              #  plotlyOutput("Volcano_DEGs_in_MA_SMC"),
                              plotlyOutput("Volcano_DEGs_in_SMC_page2"),
                              
                              radioButtons(
                                inputId = "dataset_SMC_DEGs_page2",
                                label = "SMC differentially expressed genes and their intersection with gene classifications:",
                                choices =  c("SMC_DEGs_logFC_0.10",
                                             "SMC_DEGs_logFC_0.41",
                                             "SMC_DEGs_logFC_0.41_genes",
                                             "SMC_DEGs_logFC_0.41_major_class",
                                             "SMC_DEGs_logFC_0.41_minor_class",
                                             "SMC_DEGs_logFC_0.41_intersected",
                                             "SMC_DEGs_logFC_0.41_intersected_major_class",
                                             "SMC_DEGs_logFC_0.41_intersected_minor_class"),
                                inline = TRUE),
                              
                              DT::dataTableOutput("table_SMC_DEGs_page2")
                            )
                          ))
                 ,
                 ### page 3
                 tabPanel("EC DEGs",
                          titlePanel("Differentially expressed genes in ECs in young SHR rats"),

                          sidebarLayout(
                            sidebarPanel(
                              helpText("Select genes"),
                              selectInput("genes2",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Cdh5","LOC306079"),multiple = TRUE),

                              #  filterDF_UI("filtering"),

                              width = 3,
                              plotOutput("empty_EC_page3"),
                              plotOutput("Figure_S7A"),

                              filterDF_UI("filtering_EC_DEGs_page3")

                            ),

                            mainPanel(
                              tags$h2("ECs DEGs in SHR MA and AA"),

                              plotOutput("vlnplot_dotplot_EC_page3"),

                              #  plotlyOutput("Volcano_DEGs_in_MA_EC"),
                              plotlyOutput("Volcano_DEGs_in_EC_page3"),

                              radioButtons(
                                inputId = "dataset_EC_DEGs_page3",
                                label = "EC differentially expressed genes and their intersection with gene classifications:",
                                choices =  c("EC_DEGs_logFC_0.10",
                                             "EC_DEGs_logFC_0.41",
                                             "EC_DEGs_logFC_0.41_genes",
                                             "EC_DEGs_logFC_0.41_major_class",
                                             "EC_DEGs_logFC_0.41_minor_class",
                                             "EC_DEGs_logFC_0.41_intersected",
                                             "EC_DEGs_logFC_0.41_intersected_major_class",
                                             "EC_DEGs_logFC_0.41_intersected_minor_class"),
                                inline = TRUE),

                              DT::dataTableOutput("table_EC_DEGs_page3")
                            )
                          ))
                 ,
                 #### page 4
                 tabPanel("MSC DEGs",
                          titlePanel("Differentially expressed genes in MSCs in young SHR rats"),

                          sidebarLayout(
                            sidebarPanel(
                              helpText("Select genes"),
                              selectInput("genes3",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Col8a1", "Alpl"),multiple = TRUE),

                              #  filterDF_UI("filtering"),

                              width = 3,
                              plotOutput("empty_MSC_page4"),
                              plotOutput("Figure_5A"),

                              filterDF_UI("filtering_MSC_DEGs_page4")

                            ),

                            mainPanel(
                              tags$h2("MSCs DEGs in SHR MA and AA"),

                              plotOutput("vlnplot_dotplot_MSC_page4"),

                              #  plotlyOutput("Volcano_DEGs_in_MA_MSC"),
                              plotlyOutput("Volcano_DEGs_in_MSC_page4"),

                              radioButtons(
                                inputId = "dataset_MSC_DEGs_page4",
                                label = "MSC differentially expressed genes and their intersection with gene classifications:",
                                choices =  c("MSC_DEGs_logFC_0.10",
                                             "MSC_DEGs_logFC_0.41",
                                             "MSC_DEGs_logFC_0.41_genes",
                                             "MSC_DEGs_logFC_0.41_major_class",
                                             "MSC_DEGs_logFC_0.41_minor_class",
                                             "MSC_DEGs_logFC_0.41_intersected",
                                             "MSC_DEGs_logFC_0.41_intersected_major_class",
                                             "MSC_DEGs_logFC_0.41_intersected_minor_class"),
                                inline = TRUE),

                              DT::dataTableOutput("table_MSC_DEGs_page4")
                            )
                          ))
                 ,
                 ### page 5
                 tabPanel("SMC sub-clusters",
                          titlePanel("SMC sub-clusters in young SHR rats and differentially expressed genes"),
                          
                          sidebarLayout(
                            sidebarPanel(
                              
                              tags$h4("Sub-cluster markers"),
                              
                              helpText("Select genes"),
                              selectInput("genes_SMC_sub_markers", 
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Myh11")),
                              
                              plotOutput("empty_SMC"),
                              
                              plotOutput("empty_SMC3"),
                              
                              
                              filterDF_UI("filtering_SMC_markers"),
                              
                              width = 3,
                              # plotOutput("empty_SMC2"),
                              
                              tags$h4("DEGs analysis"),
                              helpText("Select genes"),
                              selectInput("genes_SMC_sub_DEGs", 
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Myh11","Acta2"), multiple = TRUE),
                              
                              plotOutput("empty_SMC4"),
                              filterDF_UI("filtering_SMC_DEGs")
                              
                            ),
                            
                            mainPanel(
                              
                              tags$h2("SMC sub-cluster markers"),
                              
                              plotOutput("featureplot_vlnplot_SMC"),
                              
                              #  plotlyOutput("Volcano_DEGs_in_MA_SMC"),
                              plotOutput("featureplot_split_SMC"),
                              
                              
                              radioButtons(
                                inputId = "dataset_SMC_markers",
                                label = "SMC subcluster markers:",
                                choices =  c("SMC_subcluster_markers_logFC_0.10", "SMC_subcluster_markers_logFC_0.41"),
                                inline = TRUE),
                              
                              DT::dataTableOutput("table_SMC_markers"),
                              tags$h2("SMC sub-cluster differentially expressed genes"),
                              
                              plotOutput("vlnplot_dotplot_SMC"),
                              
                              radioButtons(
                                inputId = "dataset_SMC_DEGs",
                                label = "SMC subcluster DEGs:",
                                choices = c("SMC_subcluster_DEGs_logFC_0.10", "SMC_subcluster_DEGs_logFC_0.41"),
                                inline = TRUE),
                              
                              DT::dataTableOutput("table_SMC_DEGs")
                            )
                          ))
                 ,
                 ### page 6
                 tabPanel("EC sub-clusters",
                          titlePanel("EC sub-clusters in young SHR rats and differentially expressed genes"),

                          sidebarLayout(
                            sidebarPanel(

                              tags$h4("Sub-cluster markers"),

                              helpText("Select genes"),
                              selectInput("genes_EC_sub_markers",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Cdh5")),

                              plotOutput("empty_EC"),

                              plotOutput("empty_EC3"),


                              filterDF_UI("filtering_EC_markers"),

                              width = 3,
                              # plotOutput("empty_EC2"),

                              tags$h4("DEGs analysis"),
                              helpText("Select genes"),
                              selectInput("genes_EC_sub_DEGs",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Cdh5","LOC306079"), multiple = TRUE),

                              plotOutput("empty_EC4"),
                              filterDF_UI("filtering_EC_DEGs")

                            ),

                            mainPanel(

                              tags$h2("EC sub-cluster markers"),

                              plotOutput("featureplot_vlnplot_EC"),

                              #  plotlyOutput("Volcano_DEGs_in_MA_EC"),
                              plotOutput("featureplot_split_EC"),


                              radioButtons(
                                inputId = "dataset_EC_markers",
                                label = "EC subcluster markers:",
                                choices =  c("EC_subcluster_markers_logFC_0.10", "EC_subcluster_markers_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_EC_markers"),
                              tags$h2("EC sub-cluster differentially expressed genes"),

                              plotOutput("vlnplot_dotplot_EC"),

                              radioButtons(
                                inputId = "dataset_EC_DEGs",
                                label = "EC subcluster DEGs:",
                                choices =  c("EC_subcluster_DEGs_logFC_0.10", "EC_subcluster_DEGs_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_EC_DEGs")
                            )
                          ))
                 ,
                 ### page 7
                 tabPanel("MSC sub-clusters",
                          titlePanel("MSC sub-clusters in young SHR rats and differentially expressed genes"),

                          sidebarLayout(
                            sidebarPanel(

                              tags$h4("Sub-cluster markers"),

                              helpText("Select genes"),
                              selectInput("genes_MSC_sub_markers",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Col8a1")),

                              plotOutput("empty_MSC"),

                              plotOutput("empty_MSC3"),


                              filterDF_UI("filtering_MSC_markers"),

                              width = 3,
                              # plotOutput("empty_MSC2"),

                              tags$h4("DEGs analysis"),
                              helpText("Select genes"),
                              selectInput("genes_MSC_sub_DEGs",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Col8a1","Alpl"), multiple = TRUE),

                              plotOutput("empty_MSC4"),
                              filterDF_UI("filtering_MSC_DEGs")

                            ),

                            mainPanel(

                              tags$h2("MSC sub-cluster markers"),

                              plotOutput("featureplot_vlnplot_MSC"),

                              #  plotlyOutput("Volcano_DEGs_in_MA_MSC"),
                              plotOutput("featureplot_split_MSC"),


                              radioButtons(
                                inputId = "dataset_MSC_markers",
                                label = "MSC subcluster markers:",
                                choices =  c("MSC_subcluster_markers_logFC_0.10", "MSC_subcluster_markers_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_MSC_markers"),
                              tags$h2("MSC sub-cluster differentially expressed genes"),

                              plotOutput("vlnplot_dotplot_MSC"),

                              radioButtons(
                                inputId = "dataset_MSC_DEGs",
                                label = "MSC subcluster DEGs:",
                                choices =  c("MSC_subcluster_DEGs_logFC_0.10", "MSC_subcluster_DEGs_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_MSC_DEGs")
                            )
                          )) 
                 ,
                 ### page 8
                 tabPanel("Immune sub-clusters",
                          titlePanel("Immune sub-clusters in young SHR rats and differentially expressed genes"),

                          sidebarLayout(
                            sidebarPanel(

                              tags$h4("Sub-cluster markers"),

                              helpText("Select genes"),
                              selectInput("genes_Immune_sub_markers",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Cd14")),

                              plotOutput("empty_Immune"),

                              plotOutput("empty_Immune3"),


                              filterDF_UI("filtering_Immune_markers"),

                              width = 3,
                              # plotOutput("empty_Immune2"),

                              tags$h4("DEGs analysis"),
                              helpText("Select genes"),
                              selectInput("genes_Immune_sub_DEGs",
                                          label = "Choose genes (symbols)",
                                          choices = genes, selected = c("Cd14","Cxcl1"), multiple = TRUE),

                              plotOutput("empty_Immune4"),
                              filterDF_UI("filtering_Immune_DEGs")

                            ),

                            mainPanel(

                              tags$h2("Immune sub-cluster markers"),

                              plotOutput("featureplot_vlnplot_Immune"),

                              #  plotlyOutput("Volcano_DEGs_in_MA_Immune"),
                              plotOutput("featureplot_split_Immune"),


                              radioButtons(
                                inputId = "dataset_Immune_markers",
                                label = "Immune subcluster markers:",
                                choices =  c("Immune_subcluster_markers_logFC_0.10", "Immune_subcluster_markers_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_Immune_markers"),
                              tags$h2("Immune sub-cluster differentially expressed genes"),

                              plotOutput("vlnplot_dotplot_Immune"),

                              radioButtons(
                                inputId = "dataset_Immune_DEGs",
                                label = "Immune subcluster DEGs:",
                                choices =  c("Immune_subcluster_DEGs_logFC_0.10", "Immune_subcluster_DEGs_logFC_0.41"),
                                inline = TRUE),

                              DT::dataTableOutput("table_Immune_DEGs")
                            )
                          ))
                 
)


library(shiny)
server <- function(input, output, session) {
  
  ##########for first tabpanel - aggr.combined
  output$featureplot_overview <- renderPlot({
    plot_grid(plot1,
              FeaturePlot(aggr.combined, features = input$gene, split.by = NULL, reduction = "tsne", pt.size = 1.3,
                          cols = c("lightgrey", "red"), order = TRUE)+
                theme(legend.position = "none", axis.line = element_line(size = 0.5)) + ggtitle(paste(input$gene, "overview", sep = " ")) +
                theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0),
                       axis.text = element_text(size = 18), axis.ticks = element_line(size = 0.5), axis.title = element_text(size = 18)),

              VlnPlot(aggr.combined, features = input$gene, pt.size = 0.4,  cols = c("#35978f", "#33a02c", "#f46d43", "#d53e4f", "#5e4fa2")
              )+NoLegend()+
                theme(axis.title.y = element_text(size = 20, face  = "italic", color = "black", angle = 0, vjust = 0.5),
                      axis.title.x = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.background = element_blank(),
                      plot.title = element_text(face = 4, size = 18))+
                #geom_violin(colour= "transparent", scale = "width")+
                ylab("")+
                ggtitle(input$gene),

              ncol = 3)},
    width = 1400,height = 350)

  output$featureplot_splitview <- renderPlot({
    FeaturePlot(aggr.combined, features = input$gene, split.by = "stim1", reduction = "tsne", pt.size = 1.3,
                cols = c("lightgrey", "red"), order = TRUE)+
      theme(legend.position = "none", axis.line = element_blank()) + ggtitle(paste(input$gene, "splitview", sep = " ")) +
      theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0.03),
             axis.text = element_text(size = 18), axis.ticks = element_line(size = 4), axis.title = element_text(size = 18))



  },
  width = 1400,height = 350 )

  #######for first tabpanel - aggr.combined
  output$empty1 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$Figure_2A <- renderPlot({

    par(cex.main = 2)
    plot(Figure_2A, axes = FALSE, main = "Figure_2A")
  }
  #, width = 1600,height = 380
  )

  ##### reactive data for data for first tabpanel - aggr.combined

  res_filter <- callModule(
    module = filterDF,
    id = "filtering",
    data_table = reactive({
      CELLTYPE_dfs[[input$dataset]]
    }),
    data_vars = shiny::reactive({
      intersect( c("Type", colnames(CELLTYPE_dfs[[1]])[c(1:2,6)]),
                 colnames(CELLTYPE_dfs[[input$dataset]]) )
    }),
    drop_ids = FALSE
  )

  output$table_all <- DT::renderDataTable({
    DT::datatable( if (length(intersect( c("Type", colnames(CELLTYPE_dfs[[1]])[c(1:2,6)]),
                                         colnames(CELLTYPE_dfs[[input$dataset]]))) == 0 ){CELLTYPE_dfs[[input$dataset]]}
                   else{res_filter$data_filtered()},  options = list(pageLength = 10))
  })


  ##### page 2
  ##########for SMC DEGs tabpanel
  output$vlnplot_dotplot_SMC_page2 <- renderPlot({
    plot_grid(
      
      VlnPlot(aggr.combined_SMC, features = input$genes1[length(input$genes1)], pt.size =0,  cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8),
                                                                                                      rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(), 
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes1[length(input$genes1)]),
      
      DotPlot(aggr.combined_SMC, features = input$genes1, cols = c("white", "red")) + 
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") + 
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1), 
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),
      
      ncol = 2)}, 
    width = 1400,height = 350)
  
  #######for SMC DEGs tabpanel
  output$empty_SMC_page2 <- renderPlot({
    plot(empty1, axes = FALSE)
  })
  
  output$Figure_3A <- renderPlot({
    
    par(cex.main = 1.5)
    plot(Figure_3A, axes = FALSE, main = "Figure_3A")
  })
  
  ##########################plot input and plot
  i = 1; j = 1 ## SHR-MA and WKY-MA
  dataframe_SMC <- as.data.frame(DEGs_list[[1]][[1]])
  dataframe_SMC$gene <- rownames(dataframe_SMC)
  temp_name_SMC <- names(DEGs_list[[1]])[1]
  top_peaks_SMC <- dataframe_SMC[with(dataframe_SMC, order(avg_logFC, p_val)),][1:10,]
  top_peaks_SMC <- rbind(top_peaks_SMC, dataframe_SMC[with(dataframe_SMC, order(-avg_logFC, p_val)),][1:10,])
  a_SMC <- list()
  for (s in seq_len(nrow(top_peaks_SMC))) {
    m <- top_peaks_SMC[s, ]
    a_SMC[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }
  
  a_SMC[[nrow(top_peaks_SMC)+1]] <- list(x = 0.4, y = 1.05, text = temp_name_SMC, showarrow = F, 
                                         xref='paper', yref='paper',font = list(size = 17))
  
  # output$Volcano_DEGs_in_MA_SMC <- renderPlotly({
  #   plot_ly(data = dataframe, x = dataframe$avg_logFC, y = -log10(dataframe$p_val), 
  #           text = dataframe$gene, 
  #           mode = "markers", color = dataframe$group_color) %>% 
  #     layout(title = temp_name)  %>% layout(annotations = a)
  # })
  
  ## SHR-AA and WKY-AA
  dataframe_SMC_AA <- as.data.frame(DEGs_list[[2]][[1]])
  dataframe_SMC_AA$gene <- rownames(dataframe_SMC_AA)
  temp_name_SMC_AA <- names(DEGs_list[[2]])[1]
  top_peaks_SMC_AA <- dataframe_SMC_AA[with(dataframe_SMC_AA, order(avg_logFC, p_val)),][1:10,]
  top_peaks_SMC_AA <- rbind(top_peaks_SMC_AA, dataframe_SMC_AA[with(dataframe_SMC_AA, order(-avg_logFC, p_val)),][1:10,])
  a_SMC_AA <- list()
  for (s in seq_len(nrow(top_peaks_SMC_AA))) {
    m <- top_peaks_SMC_AA[s, ]
    a_SMC_AA[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }
  a_SMC_AA[[nrow(top_peaks_SMC_AA)+1]] <- list(x = 0.7 , y = 1.05, text = temp_name_SMC_AA, showarrow = F, 
                                               xref='paper', yref='paper', font = list(size = 17))
  # output$Volcano_DEGs_in_AA_SMC <- renderPlotly({
  #   plot_ly(data = dataframe_AA, x = dataframe_AA$avg_logFC, y = -log10(dataframe_AA$p_val), 
  #           text = dataframe_AA$gene, 
  #           mode = "markers", color = dataframe_AA$group_color) %>% 
  #     layout(title = temp_name_AA)  %>% layout(annotations = a_AA)
  # })
  
  
  output$Volcano_DEGs_in_SMC_page2 <- renderPlotly({
    
    subplot(
      style(plot_ly(data = dataframe_SMC, x = dataframe_SMC$avg_logFC, y = -log10(dataframe_SMC$p_val), 
                    text = dataframe_SMC$gene, 
                    mode = "markers", color = dataframe_SMC$group_color) # %>% layout(title = temp_name)  
            %>% layout(annotations = a_SMC), showlegend = FALSE),
      
      
      plot_ly(data = dataframe_SMC_AA, x = dataframe_SMC_AA$avg_logFC, y = -log10(dataframe_SMC_AA$p_val), 
              text = dataframe_SMC_AA$gene, 
              mode = "markers", color = dataframe_SMC_AA$group_color) # %>% layout(title = temp_name_AA) 
      %>% layout(annotations = a_SMC_AA) #%>% add_text(size = 8)
      
      , nrows = 1, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE) 
  })
  
  
  
  res_filter_SMC <- callModule(
    module = filterDF,
    id = "efiltering_SMC_DEGs_page2",
    data_table = reactive({
      SMC_DEGs_dfs[[input$dataset_SMC_DEGs_page2]]
    }),
    data_vars = shiny::reactive({
      intersect( c("Type", colnames(SMC_DEGs_dfs[[1]])[c(1:2,7)]),
                 colnames(SMC_DEGs_dfs[[input$dataset_SMC_DEGs_page2]]) )
    }),
    drop_ids = FALSE
  )
  
  output$table_SMC_DEGs_page2 <- DT::renderDataTable({
    DT::datatable( if (length(intersect( c("Type",colnames(SMC_DEGs_dfs[[1]])[c(1:2,7)]),
                                         colnames(SMC_DEGs_dfs[[input$dataset_SMC_DEGs_page2]]))) == 0 ){SMC_DEGs_dfs[[input$dataset_SMC_DEGs_page2]]}
                   else{res_filter_SMC$data_filtered()},  options = list(pageLength = 10))
  })

  ##### page 3
  ##########for EC DEGs tabpanel
  output$vlnplot_dotplot_EC_page3 <- renderPlot({
    plot_grid(

      VlnPlot(aggr.combined_EC, features = input$genes2[length(input$genes2)], pt.size =0,  cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8),
                                                                                                     rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(),
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes2[length(input$genes2)]),

      DotPlot(aggr.combined_EC, features = input$genes2, cols = c("white", "red")) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") +
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1),
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),

      ncol = 2)},
    width = 1400,height = 350)

  #######for EC DEGs tabpanel
  output$empty_EC_page3 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$Figure_S7A <- renderPlot({

    par(cex.main = 1.5)
    plot(Figure_S7A, axes = FALSE, main = "Figure_S7A")
  })

  ##########################plot input and plot
  i = 1; j = 2 ## SHR-MA and WKY-MA
  dataframe_EC <- as.data.frame(DEGs_list[[1]][[2]])
  dataframe_EC$gene <- rownames(dataframe_EC)
  temp_name_EC <- names(DEGs_list[[1]])[2]
  top_peaks_EC <- dataframe_EC[with(dataframe_EC, order(avg_logFC, p_val)),][1:10,]
  top_peaks_EC <- rbind(top_peaks_EC, dataframe_EC[with(dataframe_EC, order(-avg_logFC, p_val)),][1:10,])
  a_EC <- list()
  for (s in seq_len(nrow(top_peaks_EC))) {
    m <- top_peaks_EC[s, ]
    a_EC[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }

  a_EC[[nrow(top_peaks_EC)+1]] <- list(x = 0.4, y = 1.05, text = temp_name_EC, showarrow = F,
                                       xref='paper', yref='paper',font = list(size = 17))

  # output$Volcano_DEGs_in_MA_EC <- renderPlotly({
  #   plot_ly(data = dataframe, x = dataframe$avg_logFC, y = -log10(dataframe$p_val),
  #           text = dataframe$gene,
  #           mode = "markers", color = dataframe$group_color) %>%
  #     layout(title = temp_name)  %>% layout(annotations = a)
  # })

  ## SHR-AA and WKY-AA
  dataframe_EC_AA <- as.data.frame(DEGs_list[[2]][[2]])
  dataframe_EC_AA$gene <- rownames(dataframe_EC_AA)
  temp_name_EC_AA <- names(DEGs_list[[2]])[2]
  top_peaks_EC_AA <- dataframe_EC_AA[with(dataframe_EC_AA, order(avg_logFC, p_val)),][1:10,]
  top_peaks_EC_AA <- rbind(top_peaks_EC_AA, dataframe_EC_AA[with(dataframe_EC_AA, order(-avg_logFC, p_val)),][1:10,])
  a_EC_AA <- list()
  for (s in seq_len(nrow(top_peaks_EC_AA))) {
    m <- top_peaks_EC_AA[s, ]
    a_EC_AA[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }
  a_EC_AA[[nrow(top_peaks_EC_AA)+1]] <- list(x = 0.7 , y = 1.05, text = temp_name_EC_AA, showarrow = F,
                                             xref='paper', yref='paper', font = list(size = 17))
  # output$Volcano_DEGs_in_AA_EC <- renderPlotly({
  #   plot_ly(data = dataframe_AA, x = dataframe_AA$avg_logFC, y = -log10(dataframe_AA$p_val),
  #           text = dataframe_AA$gene,
  #           mode = "markers", color = dataframe_AA$group_color) %>%
  #     layout(title = temp_name_AA)  %>% layout(annotations = a_AA)
  # })


  output$Volcano_DEGs_in_EC_page3 <- renderPlotly({

    subplot(
      style(plot_ly(data = dataframe_EC, x = dataframe_EC$avg_logFC, y = -log10(dataframe_EC$p_val),
                    text = dataframe_EC$gene,
                    mode = "markers", color = dataframe_EC$group_color) # %>% layout(title = temp_name)
            %>% layout(annotations = a_EC), showlegend = FALSE),


      plot_ly(data = dataframe_EC_AA, x = dataframe_EC_AA$avg_logFC, y = -log10(dataframe_EC_AA$p_val),
              text = dataframe_EC_AA$gene,
              mode = "markers", color = dataframe_EC_AA$group_color) # %>% layout(title = temp_name_AA)
      %>% layout(annotations = a_EC_AA) #%>% add_text(size = 8)

      , nrows = 1, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE)
  })



  res_filter_EC <- callModule(
    module = filterDF,
    id = "filtering_EC_DEGs_page3",
    data_table = reactive({
      EC_DEGs_dfs[[input$dataset_EC_DEGs_page3]]
    }),
    data_vars = shiny::reactive({
      intersect( c("Type", colnames(EC_DEGs_dfs[[1]])[c(1:2,7)]),
                 colnames(EC_DEGs_dfs[[input$dataset_EC_DEGs_page3]]) )
    }),
    drop_ids = FALSE
  )

  output$table_EC_DEGs_page3 <- DT::renderDataTable({
    DT::datatable( if (length(intersect( c("Type",colnames(EC_DEGs_dfs[[1]])[c(1:2,7)]),
                                         colnames(EC_DEGs_dfs[[input$dataset_EC_DEGs_page3]]))) == 0 ){EC_DEGs_dfs[[input$dataset_EC_DEGs_page3]]}
                   else{res_filter_EC$data_filtered()},  options = list(pageLength = 10))
  })

  ### page 4
  ##########for MSC DEGs tabpanel
  output$vlnplot_dotplot_MSC_page4 <- renderPlot({
    plot_grid(

      VlnPlot(aggr.combined_MSC, features = input$genes3[length(input$genes3)], pt.size =0,  cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8),
                                                                                                      rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(),
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes3[length(input$genes3)]),

      DotPlot(aggr.combined_MSC, features = input$genes3, cols = c("white", "red")) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") +
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1),
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),

      ncol = 2)},
    width = 1400,height = 350)

  #######for MSC DEGs tabpanel
  output$empty_MSC_page4 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$Figure_5A <- renderPlot({

    par(cex.main = 1.5)
    plot(Figure_5A, axes = FALSE, main = "Figure_5A")

  })

  ##########################plot input and plot
  i = 1; j = 1 ## SHR-MA and WKY-MA
  dataframe_MSC <- as.data.frame(DEGs_list[[1]][[3]])
  dataframe_MSC$gene <- rownames(dataframe_MSC)
  temp_name_MSC <- names(DEGs_list[[1]])[3]
  top_peaks_MSC <- dataframe_MSC[with(dataframe_MSC, order(avg_logFC, p_val)),][1:10,]
  top_peaks_MSC <- rbind(top_peaks_MSC, dataframe_MSC[with(dataframe_MSC, order(-avg_logFC, p_val)),][1:10,])
  a_MSC <- list()
  for (s in seq_len(nrow(top_peaks_MSC))) {
    m <- top_peaks_MSC[s, ]
    a_MSC[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }

  a_MSC[[nrow(top_peaks_MSC)+1]] <- list(x = 0.4, y = 1.05, text = temp_name_MSC, showarrow = F,
                                         xref='paper', yref='paper',font = list(size = 17))

  # output$Volcano_DEGs_in_MA_MSC <- renderPlotly({
  #   plot_ly(data = dataframe, x = dataframe$avg_logFC, y = -log10(dataframe$p_val),
  #           text = dataframe$gene,
  #           mode = "markers", color = dataframe$group_color) %>%
  #     layout(title = temp_name)  %>% layout(annotations = a)
  # })

  ## SHR-AA and WKY-AA
  dataframe_MSC_AA <- as.data.frame(DEGs_list[[2]][[3]])
  dataframe_MSC_AA$gene <- rownames(dataframe_MSC_AA)
  temp_name_MSC_AA <- names(DEGs_list[[2]])[3]
  top_peaks_MSC_AA <- dataframe_MSC_AA[with(dataframe_MSC_AA, order(avg_logFC, p_val)),][1:10,]
  top_peaks_MSC_AA <- rbind(top_peaks_MSC_AA, dataframe_MSC_AA[with(dataframe_MSC_AA, order(-avg_logFC, p_val)),][1:10,])
  a_MSC_AA <- list()
  for (s in seq_len(nrow(top_peaks_MSC_AA))) {
    m <- top_peaks_MSC_AA[s, ]
    a_MSC_AA[[s]] <- list(
      x = m[["avg_logFC"]],
      y = -log10(m[["p_val"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40,
      font = list(size = 14)
    )
  }
  a_MSC_AA[[nrow(top_peaks_MSC_AA)+1]] <- list(x = 0.7 , y = 1.05, text = temp_name_MSC_AA, showarrow = F,
                                               xref='paper', yref='paper', font = list(size = 17))
  # output$Volcano_DEGs_in_AA_MSC <- renderPlotly({
  #   plot_ly(data = dataframe_AA, x = dataframe_AA$avg_logFC, y = -log10(dataframe_AA$p_val),
  #           text = dataframe_AA$gene,
  #           mode = "markers", color = dataframe_AA$group_color) %>%
  #     layout(title = temp_name_AA)  %>% layout(annotations = a_AA)
  # })


  output$Volcano_DEGs_in_MSC_page4 <- renderPlotly({

    subplot(
      style(plot_ly(data = dataframe_MSC, x = dataframe_MSC$avg_logFC, y = -log10(dataframe_MSC$p_val),
                    text = dataframe_MSC$gene,
                    mode = "markers", color = dataframe_MSC$group_color) # %>% layout(title = temp_name)
            %>% layout(annotations = a_MSC), showlegend = FALSE),


      plot_ly(data = dataframe_MSC_AA, x = dataframe_MSC_AA$avg_logFC, y = -log10(dataframe_MSC_AA$p_val),
              text = dataframe_MSC_AA$gene,
              mode = "markers", color = dataframe_MSC_AA$group_color) # %>% layout(title = temp_name_AA)
      %>% layout(annotations = a_MSC_AA) #%>% add_text(size = 8)

      , nrows = 1, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE)
  })



  res_filter_MSC <- callModule(
    module = filterDF,
    id = "filtering_MSC_DEGs_page4",
    data_table = reactive({
      MSC_DEGs_dfs[[input$dataset_MSC_DEGs_page4]]
    }),
    data_vars = shiny::reactive({
      intersect( c("Type", colnames(MSC_DEGs_dfs[[1]])[c(1:2,7)]),
                 colnames(MSC_DEGs_dfs[[input$dataset_MSC_DEGs_page4]]) )
    }),
    drop_ids = FALSE
  )

  output$table_MSC_DEGs_page4 <- DT::renderDataTable({
    DT::datatable( if (length(intersect( c("Type",colnames(MSC_DEGs_dfs[[1]])[c(1:2,7)]),
                                         colnames(MSC_DEGs_dfs[[input$dataset_MSC_DEGs_page4]]))) == 0 ){MSC_DEGs_dfs[[input$dataset_MSC_DEGs_page4]]}
                   else{res_filter_MSC$data_filtered()},  options = list(pageLength = 10))
  })

  
  ###page 5
  output$featureplot_vlnplot_SMC <- renderPlot({
    plot_grid(plot_SMC_dim,
              FeaturePlot(aggr.combined_SMC, features = input$genes_SMC_sub_markers, split.by = NULL, reduction = "tsne", pt.size = 1, 
                          cols = c("lightgrey", "red"), order = TRUE)+
                theme(legend.position = "none", axis.line = element_line(size = 0.5)) + ggtitle(paste(input$genes_SMC_sub_markers, "overview", sep = " ")) + 
                theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0),
                       axis.text = element_text(size = 18), axis.ticks = element_line(size = 0.5), axis.title = element_text(size = 18)),
              
              # VlnPlot(aggr.combined_SMC1, features = input$genes_SMC_sub_markers, pt.size = 0.4,  cols = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d" )
              # )+NoLegend()+
              #   theme(axis.title.y = element_text(size = 20, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              #         axis.title.x = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
              #         panel.border = element_rect(colour = "black", fill=NA, size=1),
              #         panel.background = element_blank(),
              #         plot.title = element_text(face = 4, size = 18))+
              #   #geom_violin(colour= "transparent", scale = "width")+
              #   ylab("")+
              #   ggtitle(input$genes_SMC_sub_markers),
              
              ncol = 3)}, 
    width = 1400,height = 350)
  
  output$featureplot_split_SMC <- renderPlot({
    FeaturePlot(aggr.combined_SMC, features = input$genes_SMC_sub_markers, split.by = "stim2", reduction = "tsne", pt.size = 1, 
                cols = c("lightgrey", "red"), order = TRUE)+
      theme(legend.position = "none", axis.line = element_blank()) + ggtitle(paste(input$genes_SMC_sub_markers, "splitview", sep = " ")) + 
      theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0.03),
             axis.text = element_text(size = 18), axis.ticks = element_line(size = 4), axis.title = element_text(size = 18))
    
    
    
  }, 
  width = 1400,height = 350 )
  
  # output$empty_SMC2 <- renderPlot({
  #   plot(empty1, axes = FALSE)
  # })
  output$empty_SMC3 <- renderPlot({
    plot(empty1, axes = FALSE)
  })
  
  res_filter_SMC_markers <- callModule(
    module = filterDF,
    id = "filtering_SMC_markers",
    data_table = reactive({
      SMC_subcluster_markers[[input$dataset_SMC_markers]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )
  
  output$empty_SMC4 <- renderPlot({
    plot(empty1, axes = FALSE)
  })
  
  output$table_SMC_markers <- DT::renderDataTable({
    DT::datatable( res_filter_SMC_markers$data_filtered(),  options = list(pageLength = 8))
  })
  
  
  ##########for SMC DEGs tabpanel
  Idents(aggr.combined_SMC) <- aggr.combined_SMC$stim2
  output$vlnplot_dotplot_SMC <- renderPlot({
    plot_grid(
      
      VlnPlot(aggr.combined_SMC, features = input$genes_SMC_sub_DEGs[length(input$genes_SMC_sub_DEGs)], pt.size =0,  
              cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8), rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(), 
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes_SMC_sub_DEGs[length(input$genes_SMC_sub_DEGs)]),
      
      DotPlot(aggr.combined_SMC, features = input$genes_SMC_sub_DEGs, cols = c("white", "red")) + 
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") + 
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1), 
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),
      
      ncol = 2)}, 
    width = 1400,height = 350)
  
  #######for SMC DEGs tabpanel
  output$empty_SMC <- renderPlot({
    plot(empty1, axes = FALSE)
  })
  
  
  res_filter_SMC_DEGs <- callModule(
    module = filterDF,
    id = "filtering_SMC_DEGs",
    data_table = reactive({
      SMC_subcluster_DEGs[[input$dataset_SMC_DEGs]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )
  
  output$table_SMC_DEGs <- DT::renderDataTable({
    DT::datatable( res_filter_SMC_DEGs$data_filtered(),  options = list(pageLength = 8))
  })
  
  ##### page 6
  output$featureplot_vlnplot_EC <- renderPlot({
    plot_grid(plot_EC_dim,
              FeaturePlot(aggr.combined_EC, features = input$genes_EC_sub_markers, split.by = NULL, reduction = "tsne", pt.size = 1.7,
                          cols = c("lightgrey", "red"), order = TRUE)+
                theme(legend.position = "none", axis.line = element_line(size = 0.5)) + ggtitle(paste(input$genes_EC_sub_markers, "overview", sep = " ")) +
                theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0),
                       axis.text = element_text(size = 18), axis.ticks = element_line(size = 0.5), axis.title = element_text(size = 18)),

              # VlnPlot(aggr.combined_EC1, features = input$genes_EC_sub_markers, pt.size = 0.4,  cols = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d" )
              # )+NoLegend()+
              #   theme(axis.title.y = element_text(size = 20, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              #         axis.title.x = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
              #         panel.border = element_rect(colour = "black", fill=NA, size=1),
              #         panel.background = element_blank(),
              #         plot.title = element_text(face = 4, size = 18))+
              #   #geom_violin(colour= "transparent", scale = "width")+
              #   ylab("")+
              #   ggtitle(input$genes_EC_sub_markers),

              ncol = 3)},
    width = 1400,height = 350)

  output$featureplot_split_EC <- renderPlot({
    FeaturePlot(aggr.combined_EC, features = input$genes_EC_sub_markers, split.by = "stim2", reduction = "tsne", pt.size = 1.7,
                cols = c("lightgrey", "red"), order = TRUE)+
      theme(legend.position = "none", axis.line = element_blank()) + ggtitle(paste(input$genes_EC_sub_markers, "splitview", sep = " ")) +
      theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0.03),
             axis.text = element_text(size = 18), axis.ticks = element_line(size = 4), axis.title = element_text(size = 18))



  },
  width = 1400,height = 350 )

  # output$empty_EC2 <- renderPlot({
  #   plot(empty1, axes = FALSE)
  # })
  output$empty_EC3 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  res_filter_EC_markers <- callModule(
    module = filterDF,
    id = "filtering_EC_markers",
    data_table = reactive({
      EC_subcluster_markers[[input$dataset_EC_markers]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$empty_EC4 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$table_EC_markers <- DT::renderDataTable({
    DT::datatable( res_filter_EC_markers$data_filtered(),  options = list(pageLength = 8))
  })


  ##########for EC DEGs tabpanel
  Idents(aggr.combined_EC) <- aggr.combined_EC$stim2
  output$vlnplot_dotplot_EC <- renderPlot({
    plot_grid(

      VlnPlot(aggr.combined_EC, features = input$genes_EC_sub_DEGs[length(input$genes_EC_sub_DEGs)], pt.size =0,
              cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8), rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(),
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes_EC_sub_DEGs[length(input$genes_EC_sub_DEGs)]),

      DotPlot(aggr.combined_EC, features = input$genes_EC_sub_DEGs, cols = c("white", "red")) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") +
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1),
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),

      ncol = 2)},
    width = 1400,height = 350)

  #######for EC DEGs tabpanel
  output$empty_EC <- renderPlot({
    plot(empty1, axes = FALSE)
  })


  res_filter_EC_DEGs <- callModule(
    module = filterDF,
    id = "filtering_EC_DEGs",
    data_table = reactive({
      EC_subcluster_DEGs[[input$dataset_EC_DEGs]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$table_EC_DEGs <- DT::renderDataTable({
    DT::datatable( res_filter_EC_DEGs$data_filtered(),  options = list(pageLength = 8))
  })

  ##### page 7
  output$featureplot_vlnplot_MSC <- renderPlot({
    plot_grid(plot_MSC_dim,
              FeaturePlot(aggr.combined_MSC, features = input$genes_MSC_sub_markers, split.by = NULL, reduction = "tsne", pt.size = 1.7,
                          cols = c("lightgrey", "red"), order = TRUE)+
                theme(legend.position = "none", axis.line = element_line(size = 0.5)) + ggtitle(paste(input$genes_MSC_sub_markers, "overview", sep = " ")) +
                theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0),
                       axis.text = element_text(size = 18), axis.ticks = element_line(size = 0.5), axis.title = element_text(size = 18)),

              # VlnPlot(aggr.combined_MSC1, features = input$genes_MSC_sub_markers, pt.size = 0.4,  cols = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d" )
              # )+NoLegend()+
              #   theme(axis.title.y = element_text(size = 20, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              #         axis.title.x = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
              #         panel.border = element_rect(colour = "black", fill=NA, size=1),
              #         panel.background = element_blank(),
              #         plot.title = element_text(face = 4, size = 18))+
              #   #geom_violin(colour= "transparent", scale = "width")+
              #   ylab("")+
              #   ggtitle(input$genes_MSC_sub_markers),

              ncol = 3)},
    width = 1400,height = 350)

  output$featureplot_split_MSC <- renderPlot({
    FeaturePlot(aggr.combined_MSC, features = input$genes_MSC_sub_markers, split.by = "stim2", reduction = "tsne", pt.size = 1.7,
                cols = c("lightgrey", "red"), order = TRUE)+
      theme(legend.position = "none", axis.line = element_blank()) + ggtitle(paste(input$genes_MSC_sub_markers, "splitview", sep = " ")) +
      theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0.03),
             axis.text = element_text(size = 18), axis.ticks = element_line(size = 4), axis.title = element_text(size = 18))



  },
  width = 1400,height = 350 )

  # output$empty_MSC2 <- renderPlot({
  #   plot(empty1, axes = FALSE)
  # })
  output$empty_MSC3 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  res_filter_MSC_markers <- callModule(
    module = filterDF,
    id = "filtering_MSC_markers",
    data_table = reactive({
      MSC_subcluster_markers[[input$dataset_MSC_markers]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$empty_MSC4 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$table_MSC_markers <- DT::renderDataTable({
    DT::datatable( res_filter_MSC_markers$data_filtered(),  options = list(pageLength = 8))
  })


  ##########for MSC DEGs tabpanel
  Idents(aggr.combined_MSC) <- aggr.combined_MSC$stim2
  output$vlnplot_dotplot_MSC <- renderPlot({
    plot_grid(

      VlnPlot(aggr.combined_MSC, features = input$genes_MSC_sub_DEGs[length(input$genes_MSC_sub_DEGs)], pt.size =0,
              cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8), rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(),
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes_MSC_sub_DEGs[length(input$genes_MSC_sub_DEGs)]),

      DotPlot(aggr.combined_MSC, features = input$genes_MSC_sub_DEGs, cols = c("white", "red")) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") +
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1),
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),

      ncol = 2)},
    width = 1400,height = 350)

  #######for MSC DEGs tabpanel
  output$empty_MSC <- renderPlot({
    plot(empty1, axes = FALSE)
  })


  res_filter_MSC_DEGs <- callModule(
    module = filterDF,
    id = "filtering_MSC_DEGs",
    data_table = reactive({
      MSC_subcluster_DEGs[[input$dataset_MSC_DEGs]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$table_MSC_DEGs <- DT::renderDataTable({
    DT::datatable( res_filter_MSC_DEGs$data_filtered(),  options = list(pageLength = 8))
  })
  ##### page 8
  output$featureplot_vlnplot_Immune <- renderPlot({
    plot_grid(plot_Immune_dim,
              FeaturePlot(aggr.combined_Immune, features = input$genes_Immune_sub_markers, split.by = NULL, reduction = "tsne", pt.size = 1.7,
                          cols = c("lightgrey", "red"), order = TRUE)+
                theme(legend.position = "none", axis.line = element_line(size = 0.5)) + ggtitle(paste(input$genes_Immune_sub_markers, "overview", sep = " ")) +
                theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0),
                       axis.text = element_text(size = 18), axis.ticks = element_line(size = 0.5), axis.title = element_text(size = 18)),

              # VlnPlot(aggr.combined_Immune1, features = input$genes_Immune_sub_markers, pt.size = 0.4,  cols = c("#dd3497","#fd8d3c", "#6a51a3","#CB1313","#1d91c0",  "#41ab5d" )
              # )+NoLegend()+
              #   theme(axis.title.y = element_text(size = 20, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              #         axis.title.x = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(),
              #         panel.border = element_rect(colour = "black", fill=NA, size=1),
              #         panel.background = element_blank(),
              #         plot.title = element_text(face = 4, size = 18))+
              #   #geom_violin(colour= "transparent", scale = "width")+
              #   ylab("")+
              #   ggtitle(input$genes_Immune_sub_markers),

              ncol = 3)},
    width = 1400,height = 350)

  output$featureplot_split_Immune <- renderPlot({
    FeaturePlot(aggr.combined_Immune, features = input$genes_Immune_sub_markers, split.by = "stim2", reduction = "tsne", pt.size = 1.7,
                cols = c("lightgrey", "red"), order = TRUE)+
      theme(legend.position = "none", axis.line = element_blank()) + ggtitle(paste(input$genes_Immune_sub_markers, "splitview", sep = " ")) +
      theme( plot.title = element_text(color="black", size=18, face=4, hjust = 0.03),
             axis.text = element_text(size = 18), axis.ticks = element_line(size = 4), axis.title = element_text(size = 18))



  },
  width = 800,height = 350 )

  # output$empty_Immune2 <- renderPlot({
  #   plot(empty1, axes = FALSE)
  # })
  output$empty_Immune3 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  res_filter_Immune_markers <- callModule(
    module = filterDF,
    id = "filtering_Immune_markers",
    data_table = reactive({
      Immune_subcluster_markers[[input$dataset_Immune_markers]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$empty_Immune4 <- renderPlot({
    plot(empty1, axes = FALSE)
  })

  output$table_Immune_markers <- DT::renderDataTable({
    DT::datatable( res_filter_Immune_markers$data_filtered(),  options = list(pageLength = 8))
  })


  ##########for Immune DEGs tabpanel
  Idents(aggr.combined_Immune) <- aggr.combined_Immune$stim2
  output$vlnplot_dotplot_Immune <- renderPlot({
    plot_grid(

      VlnPlot(aggr.combined_Immune, features = input$genes_Immune_sub_DEGs[length(input$genes_Immune_sub_DEGs)], pt.size =0,
              cols = c(rgb(35/255,128/255,157/255, 0.8), rgb(218/255,123/255,70/255, 0.8), rgb(38/255,184/255,177/255, 0.8), rgb(246/255,61/255,24/255, 0.6)))+NoLegend()+
        theme(axis.title.y = element_text(size = 17, face  = "italic", color = "black", angle = 0, vjust = 0.5),
              axis.title.x = element_blank(),
              #  axis.text = element_blank(),
              axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              # axis.ticks = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
              panel.background = element_blank(),
              plot.title = element_text(face = "italic", size = 20))+
        geom_violin(colour= "transparent", scale = "width")+
        ylab("")+
        ggtitle(input$genes_Immune_sub_DEGs[length(input$genes_Immune_sub_DEGs)]),

      DotPlot(aggr.combined_Immune, features = input$genes_Immune_sub_DEGs, cols = c("white", "red")) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) + xlab("") +
        theme(axis.text.y = element_text(angle=0, hjust=1, face = "italic"))+
        theme(axis.line = element_line(size = 1),
              axis.ticks = element_line(size = 1),
              axis.text.y = element_text(size = 20
              ),
              axis.text.x = element_text(size = 17),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 17),
              axis.title =  element_text(size = 0))+
        ggtitle("")+
        coord_flip(),

      ncol = 2)},
    width = 800,height = 350)

  #######for Immune DEGs tabpanel
  output$empty_Immune <- renderPlot({
    plot(empty1, axes = FALSE)
  })


  res_filter_Immune_DEGs <- callModule(
    module = filterDF,
    id = "filtering_Immune_DEGs",
    data_table = reactive({
      Immune_subcluster_DEGs[[input$dataset_Immune_DEGs]]
    }),
    data_vars = shiny::reactive({
      c("Dysregulation", "Sub-cluster", "p_val", "avg_logFC")
    }),
    drop_ids = FALSE
  )

  output$table_Immune_DEGs <- DT::renderDataTable({
    DT::datatable( res_filter_Immune_DEGs$data_filtered(),  options = list(pageLength = 8))
  })

  
}

shinyApp(ui, server)



