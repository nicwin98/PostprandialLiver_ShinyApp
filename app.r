library(EnhancedVolcano)
library(tidyverse)
library(DT)
library(shiny)
library(shinyjs)
library(shinythemes) 
library(writexl)
library(plotly)

resH_shiny <- read.csv("resHealthy.csv")
resN_shiny <- read.csv("resNAFLD.csv")
resC_shiny <- read.csv("resCirrhosis.csv")
t_box_shiny_h <- read.csv("shiny_box_healthy.csv")
t_box_shiny_n <- read.csv("shiny_box_NAFLD.csv")
t_box_shiny_c <- read.csv("shiny_box_cirrhosis.csv")
GOBP_H <- read.csv("GOBP_Genes_PF_H.csv")
GOBP_N <- read.csv("GOBP_Genes_PF_N.csv")
GOBP_C <- read.csv("GOBP_Genes_PF_C.csv")
GOBP_H$Intervention <- rep("Postprandial vs Fasting - Healthy", length(GOBP_H$Description))
GOBP_H$GeneFraction <- sapply(GOBP_H$GeneRatio, function(x) eval(parse(text=x)))
GOBP_N$Intervention <- rep("Postprandial vs Fasting - NAFLD", length(GOBP_N$Description))
GOBP_N$GeneFraction <- sapply(GOBP_N$GeneRatio, function(x) eval(parse(text=x)))
GOBP_C$Intervention <- rep("Postprandial vs Fasting - Cirrhosis", length(GOBP_C$Description))
GOBP_C$GeneFraction <- sapply(GOBP_C$GeneRatio, function(x) eval(parse(text=x)))
t_box_shiny_h$Intervention <- factor(t_box_shiny_h$Intervention, levels = c("Fasting", "Postprandial"))
t_box_shiny_n$Intervention <- factor(t_box_shiny_n$Intervention, levels = c("Fasting", "Postprandial"))
t_box_shiny_c$Intervention <- factor(t_box_shiny_c$Intervention, levels = c("Fasting", "Postprandial"))

ui <- fluidPage(theme = shinytheme("paper"),
            navbarPage(title = "Differential Expression Analysis - Postprandial Liver Study",
                       tabPanel(title = "Preface",
                                fluidRow(column(12, wellPanel(tags$h4("Transcriptomic analysis of human liver samples assessing the effect of a standardized meal in healthy individuals, NAFLD, and cirrhosis patients as presented in:"),
                                                              tags$br(),
                                                              tags$i(tags$h5("A web-based browsable resource of hepatic gene expression in health and liver disease")),
                                                              tags$h6("Josephine Grandt, Christian D. Johansen, Anne Sofie H. Jensen, Mikkel P. Werge, Elias B. Rashu1, Andreas Møller, 
                                                                      Anders E. Junker, Lise Hobolth, Christian Mortensen, Mogens Vyberg, Reza Rafiolsadat Serizawa, Søren Møller, Lise Lotte Gluud, Nicolai J. Wewer Albrechtsen"),
                                                              tags$h6("Paper: X, DOI: X"),
                                                              tags$br(),
                                                              tags$h6("This app contains information on RNA sequencing data from liver samples taken from healthy individuals, NAFLD patients, and cirrhosis patients. The subjects were
                                                                      randomized into two subgroups: fasting and postprandial (intervention). Individuals in the postprandial subgroup ingested a standardized meal prior to sampling, while 
                                                                      the fasted individuals were fasted overnight. When comparing postprandial vs fasting within each group (healthy, NAFLD, and cirrhosis), differences in sequencing 
                                                                      runs were accounted for in the bioinformatic analysis (Design: ~ Sequencing + GroupIntervention). For detailed information on study design and methods 
                                                                      please refer to the abovementioned publication and the R codes available at", 
                                                                      tags$a("https://github.com/nicwin98/PostprandialLiverStudy"), ".")))),
                                fluidRow(column(12, wellPanel(
                                  tags$h5("The app contains three types of tabs where you can find the following information:"),
                                  tags$ul(
                                    tags$li(tags$strong("Differential Expression Analysis:")),
                                    tags$ul(
                                      tags$li("Results from analyses of postprandial vs fasting in: Healthy | NAFLD | Cirrhosis"),
                                      tags$ul(
                                        tags$li("An adjustable volcano plot of all analyzed genes."),
                                        tags$li("A table of all genes displaying the results of the differential expression analysis."),
                                          tags$ul(
                                            tags$li("The information in the table is available for download."))))),
                                  tags$ul(
                                    tags$li(tags$strong("Gene Ontology - Biological Pathways")),
                                    tags$ul(
                                      tags$li("All significant GOBPs are displayed in a table."),
                                      tags$li("Create a dotplot by clicking on your GOBPs of interest."),
                                      tags$ul(
                                        tags$li("The table information and the dotplot are available for download.")))),
                                  tags$ul(
                                    tags$li(tags$strong("Single Gene Expression:")),
                                    tags$ul(
                                      tags$li("Boxplot showing the sample variation for individual genes - separated for each comparison."),
                                      tags$ul(
                                        tags$li("The boxplot and datapoints are available for download.")))),
                                  tags$br(),
                                  tags$h6("Any inquires related to the application and publication should be directed 
                                                              to the corresponding author, Nicolai J. Wewer Albrechtsen: nicolai.albrechtsen@sund.ku.dk."),
                                  tags$h6("An additional app has been created in relation to the abovementioned publication:", tags$a("https://weweralbrechtsenlab.shinyapps.io/PLS_Groups/"), 
                                          "."))))),  
                    navbarMenu(title = "Differential Expression Analysis", 
                       tabPanel(title = "Healthy - Postprandial vs Fasting",
                                fluidRow(column(4,
                                                wellPanel(
                                                  sliderInput(
                                                    inputId = "num1",
                                                    label = "Choose a log2FoldChange cutoff",
                                                    value = 1,
                                                    min = 0.0,
                                                    max = 5.0,
                                                    step = 0.5),
                                                  sliderInput(
                                                    inputId = "xlim1",
                                                    label = "Adjust x-axis",
                                                    value = c(-6,7),
                                                    min = -11,
                                                    max = 7,
                                                    step = 1),
                                                  selectInput(
                                                    inputId = "radio1",
                                                    label = "Choose FDR adjusted p-value:",
                                                    selected = 0.05,
                                                    choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                         column(8,plotOutput("volcanoH"))),
                                fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                              tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                              tags$h6("The table displays information on:"),
                                                              tags$ul(
                                                                tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                                tags$li("Gene Symbol"),
                                                                tags$li("Gene Name"),
                                                                tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                                tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                                tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                                tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                                tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                                tags$li("P-value: Not corrected for multiple testing"),
                                                                tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                              tags$h6("Sometimes adjusted p-values, p-values, and stat values are missing. This is likely due to:"),
                                                              tags$ul(
                                                                tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                                tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab)."),
                                                                tags$li("Missing stat, pvalue, and padj: A global filtering step were applied post-hoc to avoid extreme fold changes
                                                                        for very lowly expressed genes (in the samples compared).")),
                                                              tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in the livers of those who ingested a 
                                                                  standardized meal (Postprandial) compared to the group of fasted individuals (Fasting). This tab only contains information on healthy individuals.")))),
                                fluidRow(column(12, wellPanel(DTOutput("listH")))),
                                fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE1",
                                                                           label = "Clear selected rows"),
                                                              downloadButton("printDEselect1", 'Download info on all selected genes'))))),
                  
                  tabPanel(title = "NAFLD - Postprandial vs Fasting",
                           fluidRow(column(4,
                                           wellPanel(
                                             sliderInput(
                                               inputId = "num2",
                                               label = "Choose a log2FoldChange cutoff",
                                               value = 1,
                                               min = 0.0,
                                               max = 5.0,
                                               step = 0.5),
                                             sliderInput(
                                               inputId = "xlim2",
                                               label = "Adjust x-axis",
                                               value = c(-9,8),
                                               min = -9,
                                               max = 8,
                                               step = 1),
                                             selectInput(
                                               inputId = "radio2",
                                               label = "Chosse FDR adjusted p-value:",
                                               selected = 0.05,
                                               choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                    column(8,plotOutput("volcanoN"))),
                           fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                         tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                         tags$h6("The table displays information on:"),
                                                         tags$ul(
                                                           tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                           tags$li("Gene Symbol"),
                                                           tags$li("Gene Name"),
                                                           tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                           tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                           tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                           tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                           tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                           tags$li("P-value: Not corrected for multiple testing"),
                                                           tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                         tags$h6("Sometimes adjusted p-values, p-values, and stat values are missing. This is likely due to:"),
                                                         tags$ul(
                                                           tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                           tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab)."),
                                                           tags$li("Missing stat, pvalue, and padj: A global filtering step were applied post-hoc to avoid extreme fold changes
                                                                        for very lowly expressed genes (in the samples compared).")),
                                                         tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in the livers of those who ingested a 
                                                                  standardized meal (Postprandial) compared to the group of fasted individuals (Fasting). This tab only contains information on NAFLD patients.")))),
                           fluidRow(column(12, wellPanel(DTOutput("listN")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE2",
                                                                      label = "Clear selected rows"),
                                                         downloadButton("printDEselect2", 'Download info on all selected genes'))))),
                  tabPanel(title = "Cirrhosis - Postprandial vs Fasting",
                           fluidRow(column(4,
                                           wellPanel(
                                             sliderInput(
                                               inputId = "num3",
                                               label = "Choose a log2FoldChange cutoff",
                                               value = 1,
                                               min = 0.0,
                                               max = 5.0,
                                               step = 0.5),
                                             sliderInput(
                                               inputId = "xlim3",
                                               label = "Adjust x-axis",
                                               value = c(-4,4),
                                               min = -10,
                                               max = 5,
                                               step = 1),
                                             selectInput(
                                               inputId = "radio3",
                                               label = "Chosse FDR adjusted p-value:",
                                               selected = 0.05,
                                               choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                    column(8,plotOutput("volcanoC"))),
                           fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                         tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                         tags$h6("The table displays information on:"),
                                                         tags$ul(
                                                           tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                           tags$li("Gene Symbol"),
                                                           tags$li("Gene Name"),
                                                           tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                           tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                           tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                           tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                           tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                           tags$li("P-value: Not corrected for multiple testing"),
                                                           tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                         tags$h6("Sometimes adjusted p-values, p-values, and stat values are missing. This is likely due to:"),
                                                         tags$ul(
                                                           tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                           tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab)."),
                                                           tags$li("Missing stat, pvalue, and padj: A global filtering step were applied post-hoc to avoid extreme fold changes
                                                                        for very lowly expressed genes (in the samples compared).")),
                                                         tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in the livers of those who ingested a 
                                                                  standardized meal (Postprandial) compared to the group of fasted individuals (Fasting). This tab only contains information on cirrhosis patients.")))),
                           fluidRow(column(12, wellPanel(DTOutput("listC")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE3",
                                                                      label = "Clear selected rows"),
                                                         downloadButton("printDEselect3", 'Download info on all selected genes')))))),
              navbarMenu(title = "Gene Ontology - Biological Pathways",    
                  tabPanel(title = "GOBP - Healthy - Postprandial vs Fasting",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched for those who received a standardized meal (Postprandial) compared
                                                            to those who fasted (Fasting). This tab only contains information for healthy individuals.")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_H", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_H", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_H1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_H2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_H")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_H",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_H")))),
                  
                  tabPanel(title = "GOBP - NAFLD - Postprandial vs Fasting",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched for those who received a standardized meal (Postprandial) compared
                                                            to those who fasted (Fasting). This tab only contains information for NAFLD patients.")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_N", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_N", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_N1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_N2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_N")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_N",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_N")))),
                  
                  tabPanel(title = "GOBP - Cirrhosis - Postprandial vs Fasting",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched for those who received a standardized meal (Postprandial) compared
                                                            to those who fasted (Fasting). This tab only contains information for cirrhosis patients.")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect_C", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot_C", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_C1",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radioGO_C2",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable_C")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP_C",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot_C"))))),
              
              navbarMenu(title = "Single Gene Expression - Boxplot",
                        tabPanel(title = "Healthy - Postprandial vs Fasting",
                                 fluidRow(column(12, wellPanel(tags$h5("The boxplot below displays expression counts normalized with DESeq2 (plotCount function) across treatment groups"),
                                                               tags$h6("The counts are NOT comparable across genes, as they are not normalized to gene length. The boxplot shows the median, 25", 
                                                                       tags$sup("th"),", and 75", tags$sup("th"),"percentiles. Points are displayed as outliers if they are above or below 
                                                                 1.5 times the interquartile range.")))),
                                 fluidRow(column(12,
                                                 wellPanel(
                                                   selectizeInput(
                                                     inputId = 'gene_h',
                                                     label = "Select your gene of interest",
                                                     choices = NULL),
                                                   downloadButton("printboxplot_H", 'Download plot as a PDF file'),
                                                   downloadButton("printdatatable_H", 'Download expression values from plot')))),
                                 fluidRow(column(12, plotlyOutput("expressionboxplot_H")))),
                        
                        tabPanel(title = "NAFLD - Postprandial vs Fasting",
                                 fluidRow(column(12, wellPanel(tags$h5("The boxplot below displays expression counts normalized with DESeq2 (plotCount function) across treatment groups"),
                                                               tags$h6("The counts are NOT comparable across genes, as they are not normalized to gene length. The boxplot shows the median, 25", 
                                                                       tags$sup("th"),", and 75", tags$sup("th"),"percentiles. Points are displayed as outliers if they are above or below 
                                                                 1.5 times the interquartile range.")))),
                                 fluidRow(column(12,
                                                 wellPanel(
                                                   selectizeInput(
                                                     inputId = 'gene_n',
                                                     label = "Select your gene of interest",
                                                     choices = NULL),
                                                   downloadButton("printboxplot_N", 'Download plot as a PDF file'),
                                                   downloadButton("printdatatable_N", 'Download expression values from plot')))),
                                 fluidRow(column(12, plotlyOutput("expressionboxplot_N")))),
                        
                        tabPanel(title = "Cirrhosis - Postprandial vs Fasting",
                                 fluidRow(column(12, wellPanel(tags$h5("The boxplot below displays expression counts normalized with DESeq2 (plotCount function) across treatment groups"),
                                                               tags$h6("The counts are NOT comparable across genes, as they are not normalized to gene length. The boxplot shows the median, 25", 
                                                                       tags$sup("th"),", and 75", tags$sup("th"),"percentiles. Points are displayed as outliers if they are above or below 
                                                                 1.5 times the interquartile range.")))),
                                 fluidRow(column(12,
                                                 wellPanel(
                                                   selectizeInput(
                                                     inputId = 'gene_c',
                                                     label = "Select your gene of interest",
                                                     choices = NULL),
                                                   downloadButton("printboxplot_C", 'Download plot as a PDF file'),
                                                   downloadButton("printdatatable_C", 'Download expression values from plot')))),
                                 fluidRow(column(12, plotlyOutput("expressionboxplot_C")))))
                  
                  
            )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'gene_h', choices = as.vector(unique(t_box_shiny_h$GeneSymbol)), selected = character(0), server = TRUE)
  data_box_H <- reactive({filter(t_box_shiny_h, GeneSymbol == input$gene_h, .preserve = TRUE)})
  updateSelectizeInput(session, 'gene_n', choices = as.vector(unique(t_box_shiny_n$GeneSymbol)), selected = character(0), server = TRUE)
  data_box_N <- reactive({filter(t_box_shiny_n, GeneSymbol == input$gene_n, .preserve = TRUE)})
  updateSelectizeInput(session, 'gene_c', choices = as.vector(unique(t_box_shiny_c$GeneSymbol)), selected = character(0), server = TRUE)
  data_box_C <- reactive({filter(t_box_shiny_c, GeneSymbol == input$gene_c, .preserve = TRUE)})
  
  reactive_H <- reactive({resH_shiny})
  reactive_N <- reactive({resN_shiny})
  reactive_C <- reactive({resC_shiny})
  
  output$volcanoH <- renderPlot({
    EnhancedVolcano(reactive_H(),
                    lab = resH_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Postprandial vs Fasting - Healthy',
                    pCutoff = as.numeric(input$radio1),
                    FCcutoff = input$num1,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim1),
                    ylim = c(0,6.2))
  })
  
  output$listH <- renderDT({datatable(resH_shiny, 
                      options = list(
                        scrollX = TRUE,
                        pageLength = 10,
                        lengthMenu = c(5,10,25,50,200),
                        filter = "bottom"
                      )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset1 <- DT::dataTableProxy("listH")
  shiny::observeEvent(input$clearRowsDE1, {
    DT::selectRows(DTreset1, NULL)
  })
  
  output$printDEselect1 = downloadHandler('postprandialVSfasting_Healthy_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listH_rows_selected
  writexl::write_xlsx(resH_shiny[srows_data, , drop = FALSE], path = file)
  })
  
  output$volcanoN <- renderPlot({
    EnhancedVolcano(reactive_N(),
                    lab = resN_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Postprandial vs Fasting - NAFLD',
                    pCutoff = as.numeric(input$radio2),
                    FCcutoff = input$num2,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim2),
                    ylim = c(0,6))
  })
  
  output$listN <- renderDT({datatable(resN_shiny, 
                                       options = list(
                                         scrollX = TRUE,
                                         pageLength = 10,
                                         lengthMenu = c(5,10,25,50,200),
                                         filter = "bottom"
                                       )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset2 <- DT::dataTableProxy("listN")
  shiny::observeEvent(input$clearRowsDE2, {
    DT::selectRows(DTreset2, NULL)
  })
  
  output$printDEselect2 = downloadHandler('postprandialVSfasting_NAFLD_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listN_rows_selected
  writexl::write_xlsx(resN_shiny[srows_data, , drop = FALSE], path = file)
  })

  output$volcanoC <- renderPlot({
    EnhancedVolcano(reactive_C(),
                    lab = resC_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Postprandial vs Fasting - Cirrhosis',
                    pCutoff = as.numeric(input$radio3),
                    FCcutoff = input$num3,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim3),
                    ylim = c(0,4.2))
  })
  
  output$listC <- renderDT({datatable(resC_shiny, 
                                       options = list(
                                         scrollX = TRUE,
                                         pageLength = 10,
                                         lengthMenu = c(5,10,25,50,200),
                                         filter = "bottom"
                                       )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset3 <- DT::dataTableProxy("listC")
  shiny::observeEvent(input$clearRowsDE3, {
    DT::selectRows(DTreset3, NULL)
  })
  
  output$printDEselect3 = downloadHandler('postprandialVSfasting_Cirrhosis_GeneList_Selected.xlsx', content = function(file) 
  {srows_data <- input$listC_rows_selected
  writexl::write_xlsx(resC_shiny[srows_data, , drop = FALSE], path = file)
  })
  
  output$expressionboxplot_H <- renderPlotly({
    validate(
      need(input$gene_h, "Select a gene above to generate the boxplot.")
    )
    ebplot <- ggplotly(ggplot(data_box_H(), aes(x=Intervention, y=count, fill = Intervention)) +
                         geom_boxplot(outlier.shape = NA) +
                         geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
                         ylab("Count") +
                         xlab("Intervention") +
                         ggtitle(paste0("Healthy individuals - ", input$gene_h," human liver RNA expression")) +
                         scale_color_manual(values = c("blue","red")) +
                         theme_bw())
    print(ebplot)
  })
  
  output$printboxplot_H <- downloadHandler(
    filename = function() {
      paste0('ExpressionBoxPlot_PLS_Healthy_',input$gene_h ,"_", Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, 
             ggplot(data_box_H(), aes(x=Intervention, y=count, fill = Intervention)) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
               ylab("Count") +
               xlab("Intervention") +
               ggtitle(paste0("Healthy individuals - ", input$gene_h," human liver RNA expression")) +
               scale_color_manual(values = c("blue","red")) +
               theme_bw(), 
             dpi = 700, height = 7, width = 7)
    })
  
  output$printdatatable_H <- downloadHandler(
    filename = function() {paste('NormCounts_PLS_Healthy_', input$gene_h ,"_", Sys.Date(), '.xlsx', sep='')},
    content = function(file) 
    {writexl::write_xlsx(dplyr::filter(data_box_H()), path = file)
    })
 
  output$expressionboxplot_N <- renderPlotly({
    validate(
      need(input$gene_n, "Select a gene above to generate the boxplot.")
    )
    ebplot <- ggplotly(ggplot(data_box_N(), aes(x=Intervention, y=count, fill = Intervention)) +
                         geom_boxplot(outlier.shape = NA) +
                         geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
                         ylab("Count") +
                         xlab("Intervention") +
                         ggtitle(paste0("NAFLD patients - ", input$gene_n," human liver RNA expression")) +
                         scale_color_manual(values = c("blue","red")) +
                         theme_bw())
    print(ebplot)
  })
  
  output$printboxplot_N <- downloadHandler(
    filename = function() {
      paste0('ExpressionBoxPlot_PLS_NAFLD_',input$gene_n ,"_", Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, 
             ggplot(data_box_N(), aes(x=Intervention, y=count, fill = Intervention)) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
               ylab("Count") +
               xlab("Intervention") +
               ggtitle(paste0("NAFLD patients - ", input$gene_n," human liver RNA expression")) +
               scale_color_manual(values = c("blue","red")) +
               theme_bw(), 
             dpi = 700, height = 7, width = 7)
    })
  
  output$printdatatable_N <- downloadHandler(
    filename = function() {paste('NormCounts_PLS_NAFLD_', input$gene_n ,"_", Sys.Date(), '.xlsx', sep='')},
    content = function(file) 
    {writexl::write_xlsx(dplyr::filter(data_box_N()), path = file)
    })
  
  output$expressionboxplot_C <- renderPlotly({
    validate(
      need(input$gene_c, "Select a gene above to generate the boxplot.")
    )
    ebplot <- ggplotly(ggplot(data_box_C(), aes(x=Intervention, y=count, fill = Intervention)) +
                         geom_boxplot(outlier.shape = NA) +
                         geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
                         ylab("Count") +
                         xlab("Intervention") +
                         ggtitle(paste0("Cirrhosis patients - ", input$gene_c," human liver RNA expression")) +
                         scale_color_manual(values = c("blue","red")) +
                         theme_bw())
    print(ebplot)
  })
  
  output$printboxplot_C <- downloadHandler(
    filename = function() {
      paste0('ExpressionBoxPlot_PLS_Cirrhosis_',input$gene_c ,"_", Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, 
             ggplot(data_box_C(), aes(x=Intervention, y=count, fill = Intervention)) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter(aes(col = Intervention), alpha = 0.6, color = "black", width = 0.15) +
               ylab("Count") +
               xlab("Intervention") +
               ggtitle(paste0("Cirrhosis patients - ", input$gene_c," human liver RNA expression")) +
               scale_color_manual(values = c("blue","red")) +
               theme_bw(), 
             dpi = 700, height = 7, width = 7)
    })
  
  output$printdatatable_C <- downloadHandler(
    filename = function() {paste('NormCounts_PLS_Cirrhosis_', input$gene_c ,"_", Sys.Date(), '.xlsx', sep='')},
    content = function(file) 
    {writexl::write_xlsx(dplyr::filter(data_box_C()), path = file)
    })
  
  output$GOBPtable_H <- renderDT(GOBP_H,
                               filter = "bottom",
                               options = list(
                                 autoWidth = TRUE,
                                 scrollX = TRUE,
                                 rowCallback = JS(
                                   "function(row, data) {",
                                   "for (i = 1; i < data.length; i++) {",
                                   "if (data[i]<0.01){",
                                   "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                   "}",
                                   "}",
                                   "}"),
                                 columnDefs = list(list(width = '250px', targets = c(3)),
                                                   list(width = '100px', targets = c(1)),
                                                   list(width = '55px', targets = c(6,7)),
                                                   list(visible=FALSE, targets=c(9,10))),
                                 pageLength = 10,
                                 lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_H <- DT::dataTableProxy("GOBPtable_H")
  shiny::observeEvent(input$clearRowsGOBP_H, {
    DT::selectRows(DTreset_GO_H, NULL)
  })
  
  output$GOBPlot_H <- renderPlot({
    validate(
      need(input$GOBPtable_H_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_H_rows_selected
    GOBP_H <- GOBP_H[srows_data, , drop = FALSE]
    ggplot(GOBP_H, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Intervention, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_H = downloadHandler('Postprandial_Healthy_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_H_rows_selected
  writexl::write_xlsx(GOBP_H[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_H <- downloadHandler(
    filename = function() {
      paste('DotPlot_Postprandial_Healthy_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_H_rows_selected
      GOBP_H <- GOBP_H[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_H, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Intervention, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_H1), height = as.numeric(input$radioGO_H2))
    })
  
  output$GOBPtable_N <- renderDT(GOBP_N,
                                  filter = "bottom",
                                  options = list(
                                    autoWidth = TRUE,
                                    scrollX = TRUE,
                                    rowCallback = JS(
                                      "function(row, data) {",
                                      "for (i = 1; i < data.length; i++) {",
                                      "if (data[i]<0.01){",
                                      "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                      "}",
                                      "}",
                                      "}"),
                                    columnDefs = list(list(width = '250px', targets = c(3)),
                                                      list(width = '100px', targets = c(1)),
                                                      list(width = '55px', targets = c(6,7)),
                                                      list(visible=FALSE, targets=c(9,10))),
                                    pageLength = 10,
                                    lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_N <- DT::dataTableProxy("GOBPtable_N")
  shiny::observeEvent(input$clearRowsGOBP_N, {
    DT::selectRows(DTreset_GO_N, NULL)
  })
  
  output$GOBPlot_N <- renderPlot({
    validate(
      need(input$GOBPtable_N_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_N_rows_selected
    GOBP_N <- GOBP_N[srows_data, , drop = FALSE]
    ggplot(GOBP_N, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Intervention, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_N = downloadHandler('Postprandial_NAFLD_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_N_rows_selected
  writexl::write_xlsx(GOBP_N[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_N <- downloadHandler(
    filename = function() {
      paste('DotPlot_Postprandial_NAFLD_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_N_rows_selected
      GOBP_N <- GOBP_N[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_N, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Intervention, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_N1), height = as.numeric(input$radioGO_N2))
    })
  
  output$GOBPtable_C <- renderDT(GOBP_C,
                                  filter = "bottom",
                                  options = list(
                                    autoWidth = TRUE,
                                    scrollX = TRUE,
                                    rowCallback = JS(
                                      "function(row, data) {",
                                      "for (i = 1; i < data.length; i++) {",
                                      "if (data[i]<0.01){",
                                      "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                      "}",
                                      "}",
                                      "}"),
                                    columnDefs = list(list(width = '250px', targets = c(3)),
                                                      list(width = '100px', targets = c(1)),
                                                      list(width = '55px', targets = c(6,7)),
                                                      list(visible=FALSE, targets=c(9,10))),
                                    pageLength = 10,
                                    lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset_GO_C <- DT::dataTableProxy("GOBPtable_C")
  shiny::observeEvent(input$clearRowsGOBP_C, {
    DT::selectRows(DTreset_GO_C, NULL)
  })
  
  output$GOBPlot_C <- renderPlot({
    validate(
      need(input$GOBPtable_C_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_C_rows_selected
    GOBP_C <- GOBP_C[srows_data, , drop = FALSE]
    ggplot(GOBP_C, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Intervention, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect_C = downloadHandler('Postprandial_Cirrhosis_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_C_rows_selected
  writexl::write_xlsx(GOBP_C[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot_C <- downloadHandler(
    filename = function() {
      paste('DotPlot_Postprandial_Cirrhosis_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_C_rows_selected
      GOBP_C <- GOBP_C[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_C, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Intervention, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radioGO_C1), height = as.numeric(input$radioGO_C2))
    })
  
}

shinyApp(ui, server)