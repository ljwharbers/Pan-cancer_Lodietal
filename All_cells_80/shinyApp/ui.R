library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(shinyBS)
library(shinyLP)
library(shinythemes)
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
msigdbr_list <- readRDS("hallmarks.rds")
msigdbr_list

tcels <- a("T-cells", href="http://apps.lambrechtslab.sites.vib.be/T-cells/")
bcels <- a("B-cells", href="http://apps.lambrechtslab.sites.vib.be/B-cells")
mmcels <- a("Myeloid cells", href="http://apps.lambrechtslab.sites.vib.be/Myeloid-cells")
dc <- a("Dendritic cells (DCs)",href="http://apps.lambrechtslab.sites.vib.be/Dendritic-cells")
ecels <- a("Endothelial cells (ECs)", href="http://apps.lambrechtslab.sites.vib.be/Endothelial-cells/")

jumbotron_modified <- function (header, content, button = TRUE, ...)
{
  button_label = c(...)
  if (button) {
    div(class = "jumbotron", h1(header), p(content), p(a(class = "btn btn-primary btn-lg button",
                                                         id = "tabBut", button_label)))
  }
  else {
    div(class = "jumbotron",
        style="background-image: url(https://media.istockphoto.com/id/1184711342/photo/fluorescent-imaging-immunofluorescence-of-cancer-cells-growing-in-2d-with-nuclei-in-blue.jpg?s=612x612&w=0&k=20&c=f_OSML-2BbiKAHYI02p2ONAyXwT7FQO7iKHh8_wPdek=);width: 100%;",
        h1(header,style="color:white"), p(content))
  }
}


### Start server code
shinyUI(fluidPage(
  ### HTML formatting of error messages

  tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))),
  list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 14px; }"))),


  ### Page title
  titlePanel(title=div(img(src="vib_logo.png", height = '50', width = '130'), HTML("Pan-Cancer atlas scRNA-Seq")),windowTitle ="Pan-Cancer atlas" ),
  navbarPage(
    NULL,
        tabPanel(
      HTML("Home Page"),
      jumbotron_modified("PanCancer Atlas",content="",button = FALSE),
      hr(),
      # h5("Welcome to our PanCancer shinyApp. The aim of the app is to allow the easy exploration of gene and gene signature expressions throughout different cancer types."),
      # h5("Due to the large size of our PanCancer atlas and the limitations in processing power of ShinyApps in this app we had to randomly sample 40% of cells from our Atlas. However, we have other dedicated shinyApps with all cells for the following cell types:"),
      # h5("Note! Here I changed the files so all cells are present (instead of 40%). Problem: to much RAM needed. Fixed?"),
      tags$style("#project-grid {
                      display: grid;
                      grid-template-columns: 600px 1fr;
                      grid-gap: 10px;
                      }"),
      # h4("Exploring the Tumor Microenvironment with Single-Cell RNA Sequencing"),
      # div(id = "project-grid",
      #   div(img(src = 'Fig1A_cropped.png', style = 'border-radius: 50%', width = '600')),
      #   div(
      #   h5("The tumor microenvironment (TME) plays a crucial role in shaping tumor evolution and determining therapy effectiveness.
      #   Single-cell RNA sequencing (scRNA-seq) has revolutionized our understanding of the TME by enabling detailed insights at the cellular level.
      #   While previous studies often focused on a specific cell type and cancer type or combined data from various sources, our approach is distinct.
      #   We present a comprehensive, uniform pan-cancer scRNA-seq dataset, comprising 683,184 high-quality single-cell transcriptomes from 234 treatment-naïve samples across", tags$b("9 different cancer types"),", all generated in-house to minimize technological variability.
      #   Within this dataset, we have identified ", tags$b("70 shared subclusters"), " that provide deeper insights into TME heterogeneity.
      #   ")
      #   )
      #   ),


      h4("Exploring the Tumor Microenvironment with Single-Cell RNA Sequencing"),
      h5("
      The tumor microenvironment (TME) plays a crucial role in shaping tumor evolution and determining therapy effectiveness. Single-cell RNA sequencing (scRNA-seq) has revolutionized our understanding of the TME by enabling detailed insights at the cellular level.
      While previous studies often focused on a specific cell type and cancer type or combined data from various sources, our approach is distinct.
      We present a comprehensive, uniform pan-cancer scRNA-seq dataset, comprising 683,184 high-quality single-cell transcriptomes from 234 treatment-naïve samples across",
      tags$b("9 different cancer types,"), " all generated in-house to minimize technological variability.
      Within this dataset, we have identified ", tags$b("70 shared subclusters"), " that provide deeper insights into TME heterogeneity.
      "),
      img(src = 'Fig1A_cropped.png', width = '500'),
      br(),
      br(),
      h4("Welcome to our Pan-Cancer ShinyApps"),
      h5("These ShinyApps are designed to facilitate the exploration of gene and gene signature expression across 9 cancer types.
      It allows users to examine individual gene expression, analyze customizable gene signatures, rank tumors,
      and explore gene expression changes within specific subclusters.
      Besides this major app where all 70 shared subclusters can be explored simultaneously,
      we have individual apps for individual cell types.
      "),
      h5("Below, you can find links to the ShinyApps for the individual cell types:"),
      # h5("We have other dedicated shinyApps with all cells for the following cell types:"),
      h5("•" ,tcels),
      h5("•" ,bcels),
      h5("•" ,mmcels),
      h5("•" ,dc),
      h5("•" ,ecels),
      br(),
      img(src="UMAPs_cropped.png", width = '1200'),
      br(),
      br(),
      br(),
      br(),
      h5("Reference: Lodi  et al. 2024, in preparation."),
    ),

    ### Tab1.c1: violinplot / boxplot
    tabPanel(
      HTML("Gene Expression (vlnplot/boxplot)"),
      h4("Cell information / gene expression violin plot / box plot"),
      fluidRow(
        column(
          3, style="border-right: 2px solid black",
          selectInput("sc1c1inp1", "Cell information (X-axis):",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Cell information to group cells by",
                   content = c("Select categorical cell information to group cells by",
                               "- Single cells are grouped by this categorical covariate",
                               "- Plotted as the X-axis of the violin plot / box plot")),
          selectInput("sc1c1inp2", "Gene name (Y-axis):", choices=NULL,multiple = TRUE) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Gene to plot",
                   content = c("Gene to plot on Y-axis",
                               "- Can be continuous cell information (e.g. nUMIs / scores)",
                               "- Can also be gene expression")),
          radioButtons("sc1c1typ", "Plot type:",
                       choices = c("violin", "boxplot"),
                       selected = "violin", inline = TRUE),
          checkboxInput("sc1c1pts", "Show data points", value = FALSE),
          actionButton("goButton_vio", "Generate / Update plot", class = "btn-success"),
          br(),br(),
          actionButton("sc1c1togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1c1togL % 2 == 1",
            selectInput("sc1c1sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1c1sub1.ui"),
            actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary")
          ),
          conditionalPanel(
            condition = "input.sc1c1togL % 2 == 1",
            selectInput("sc1c1sub3", "",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp2),
            uiOutput("sc1c1sub3.ui"),
            actionButton("sc1c1sub3all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1c1sub3non", "Deselect all groups", class = "btn btn-primary")
          ), br(), br(),
          actionButton("sc1c1tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1c1tog % 2 == 1",
            sliderInput("sc1c1siz", "Data point size:",
                        min = 0, max = 4, value = 0.25, step = 0.25),
            radioButtons("sc1c1psz", "Plot size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Medium", inline = TRUE),
            radioButtons("sc1c1fsz", "Font size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Small", inline = TRUE))
        ), # End of column (6 space)
        column(9, uiOutput("sc1c1oup.ui"),
               downloadButton("sc1c1oup.pdf", "Download PDF"),
               downloadButton("sc1c1oup.png", "Download PNG"), br(),
               div(style="display:inline-block",
                   numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px",
                                min = 4, max = 20, value = 6, step = 0.5)),
               div(style="display:inline-block",
                   numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px",
                                min = 4, max = 20, value = 10, step = 0.5))
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    ),     # End of tab (2 space)

  tabPanel(
    HTML("Gene Expression (bubbleplot/heatmap)"),
    h4("Gene expression bubbleplot / heatmap"),
    "In this tab, users can visualise the gene expression patterns of ",
    "multiple genes grouped by categorical cell information (e.g. cancer type / cell subtype).", br(),
    "The normalised expression are averaged, log-transformed and then plotted.",
    br(),br(),
    fluidRow(
      column(
        3, style="border-right: 2px solid black",
        textAreaInput("sc1bb1inp", HTML("List of gene names <br />
                                            (Max 50 genes, separated <br />
                                             by , or ; or newline):"),
                      height = "200px",
                      value = paste0(sc1def$genes, collapse = ", ")) %>%
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "List of genes to plot on bubbleplot / heatmap",
                 content = c("Input genes to plot",
                             "- Maximum 50 genes (due to ploting space limitations)",
                             "- Genes should be separated by comma, semicolon or newline")),
        selectInput("sc1bb1grp", "Group by:",
                    choices = sc1conf[grp == TRUE]$UI,
                    selected = sc1conf[grp == TRUE]$UI[9]) %>%
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "Cell information to group cells by",
                 content = c("Select categorical cell information to group cells by",
                             "- Single cells are grouped by this categorical covariate",
                             "- Plotted as the X-axis of the bubbleplot / heatmap")),
        radioButtons("sc1bb1plt", "Plot type:",
                     choices = c("Bubbleplot", "Heatmap"),
                     selected = "Bubbleplot", inline = TRUE),
        checkboxInput("sc1bb1scl", "Scale gene expression", value = TRUE),
        checkboxInput("sc1bb1row", "Cluster rows (genes)", value = TRUE),
        checkboxInput("sc1bb1col", "Cluster columns (samples)", value = FALSE),
        br(),

        #### added:
        actionButton("goButton_bb1", "Generate / Update plot", class = "btn-success"),
        br(),br(),
        #### END

        actionButton("sc1bb1togL", "Toggle to subset cells"),
        conditionalPanel(
          condition = "input.sc1bb1togL % 2 == 1",
          selectInput("sc1bb1sub1", "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput("sc1bb1sub1.ui"),
          actionButton("sc1bb1sub1all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc1bb1sub1non", "Deselect all groups", class = "btn btn-primary")
        ),


        #### added:
        conditionalPanel(
          condition = "input.sc1bb1togL % 2 == 1",
          selectInput("sc1bb1sub3", "",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp2),
          uiOutput("sc1bb1sub3.ui"),
          actionButton("sc1bb1sub3all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc1bb1sub3non", "Deselect all groups", class = "btn btn-primary")
        ),
        br(), br(),
        #### END

        actionButton("sc1bb1tog", "Toggle graphics controls"),
        conditionalPanel(
          condition = "input.sc1bb1tog % 2 == 1",
          radioButtons("sc1bb1cols", "Colour scheme:",
                       choices = c("White-Red", "Blue-Yellow-Red",
                                   "Yellow-Green-Purple"),
                       selected = "Blue-Yellow-Red"),
          radioButtons("sc1bb1psz", "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE),
          radioButtons("sc1bb1fsz", "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Medium", inline = TRUE))
      ), # End of column (6 space)
      column(9, h4(htmlOutput("sc1bb1oupTxt")),
             uiOutput("sc1bb1oup.ui"),
             downloadButton("sc1bb1oup.pdf", "Download PDF"),
             downloadButton("sc1bb1oup.png", "Download PNG"), br(),
             div(style="display:inline-block",
                 numericInput("sc1bb1oup.h", "PDF / PNG height:", width = "138px",
                              min = 4, max = 20, value = 10, step = 0.5)),
             div(style="display:inline-block",
                 numericInput("sc1bb1oup.w", "PDF / PNG width:", width = "138px",
                              min = 4, max = 20, value = 10, step = 0.5))
      )  # End of column (6 space)
    )    # End of fluidRow (4 space)
  ),      # End of tab (2 space)





  ### Tab1.bb2: Gene expr two variables
  tabPanel(
    HTML("Gene Expression (bubbleplot/heatmap) 2 groups"),
    h4("Gene expression bubbleplot / heatmap"),
    "In this tab, users can visualise the gene expression patterns of ",
    "one gene grouped by two types of categorical cell information (e.g. cancer type / cell subtype).", br(),
    "The normalised expression are averaged, log-transformed and then plotted.",
    br(),br(),


    # Input gene (new)
    fluidRow(
      column(
        3, style="border-right: 2px solid black",
        # 6, style="border-right: 2px solid black", h4("Gene expression 1"),
        fluidRow(
          column(
            6, selectInput("sc1bb2inp", "Gene name:", choices=NULL) %>%
              helper(type = "inline", size = "m", fade = TRUE,
                     title = "Select gene you want to plot" # ,
                     # content = c("Select gene to colour cells by gene expression",
                     #             paste0("- Gene expression are coloured in a ",
                     #                    "White-Red colour scheme which can be ",
                     #                    "changed in the plot controls")))
              )
          ) # end of column
        ) # end of fluidRow
        ,

        #### Group by 1 (row)
        selectInput("sc1bb2grp1", "Group rows by:",
                    choices = sc1conf[grp == TRUE]$UI,
                    selected = sc1conf[grp == TRUE]$UI[9]) %>%
          #### ? button
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "Cell information to group cells by",
                 content = c("Select categorical cell information to group cells by",
                             "- Single cells are grouped by this categorical covariate",
                             "- Plotted as the X-axis of the bubbleplot / heatmap")),


        #### Group by 2 (column)
        selectInput("sc1bb2grp2", "Group columns by:",
                    choices = sc1conf[grp == TRUE]$UI,
                    selected = sc1conf[grp == TRUE]$UI[6]) %>%
          #### ? button
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "Cell information to group cells by",
                 content = c("Select categorical cell information to group cells by",
                             "- Single cells are grouped by this categorical covariate",
                             "- Plotted as the Y-axis of the bubbleplot / heatmap")),

        #### Choose plot type
        radioButtons("sc1bb2plt", "Plot type:",
                     choices = c("Bubbleplot", "Heatmap"),
                     selected = "Bubbleplot", inline = TRUE),

        #### Select options
        checkboxInput("sc1bb2scl", "Scale gene expression (by rows)", value = FALSE),
        checkboxInput("sc1bb2row", "Cluster rows", value = FALSE),
        checkboxInput("sc1bb2col", "Cluster columns", value = FALSE),

        br(),

        #### added:
        actionButton("goButton_bb2", "Generate / Update plot", class = "btn-success"),
        br(),br(),

        #### Toggle to subset cells
        actionButton("sc1bb2togL", "Toggle to subset cells"),

        #### show if clicked on button
        conditionalPanel(
          condition = "input.sc1bb2togL % 2 == 1",

          #### Select which metadata column (variable) to subset from
          selectInput("sc1bb2sub1", "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1conf[grp == TRUE]$UI[9]),

          #### gives all the levels (options) for chosen metadata column (variable)
          uiOutput("sc1bb2sub1.ui"),

          actionButton("sc1bb2sub1all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc1bb2sub1non", "Deselect all groups", class = "btn btn-primary")

        ),

        #### added:
        conditionalPanel(
          condition = "input.sc1bb2togL % 2 == 1",
          selectInput("sc1bb2sub3", "",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1),
          uiOutput("sc1bb2sub3.ui"),
          actionButton("sc1bb2sub3all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc1bb2sub3non", "Deselect all groups", class = "btn btn-primary")
        ),
        ####

        br(), br(),

        actionButton("sc1bb2tog", "Toggle graphics controls"),
        conditionalPanel(
          condition = "input.sc1bb2tog % 2 == 1",
          radioButtons("sc1bb2cols", "Colour scheme:",
                       choices = c("White-Red", "Blue-Yellow-Red",
                                   "Yellow-Green-Purple"),
                       selected = "Blue-Yellow-Red"),
          radioButtons("sc1bb2psz", "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE),
          radioButtons("sc1bb2fsz", "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE))
      ), # End of column (6 space)

      column(9, h4(htmlOutput("sc1bb2oupTxt")),
             uiOutput("sc1bb2oup.ui"),
             downloadButton("sc1bb2oup.pdf", "Download PDF"),
             downloadButton("sc1bb2oup.png", "Download PNG"), br(),
             div(style="display:inline-block",
                 numericInput("sc1bb2oup.h", "PDF / PNG height:", width = "138px",
                              min = 4, max = 20, value = 10, step = 0.5)),
             div(style="display:inline-block",
                 numericInput("sc1bb2oup.w", "PDF / PNG width:", width = "138px",
                              min = 4, max = 20, value = 20, step = 0.5))
      )  # End of column (6 space)
    )    # End of fluidRow (4 space)
  ),      # End of tab (2 space)

    ### Tab1.c2: Proportion plot
    tabPanel(
      HTML("Cells proportion"),
      h4("Proportion / cell numbers across different cell information"),
      fluidRow(
        column(
          3, style="border-right: 2px solid black",
          selectInput("sc1c2inp1", "Cell information to plot (X-axis):",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp2) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Cell information to plot cells by",
                   content = c("Select categorical cell information to plot cells by",
                               "- Plotted as the X-axis of the proportion plot")),
          selectInput("sc1c2inp2", "Cell information to group / colour by:",
                      choices = sc1conf[grp == TRUE]$UI,
                      selected = sc1def$grp1) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Cell information to group / colour cells by",
                   content = c("Select categorical cell information to group / colour cells by",
                               "- Proportion / cell numbers are shown in different colours")),
          radioButtons("sc1c2typ", "Plot value:",
                       choices = c("Proportion", "CellNumbers"),
                       selected = "Proportion", inline = TRUE),
          checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE),
          actionButton("sc1c2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1c2togL % 2 == 1",
            selectInput("sc1c2sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1c2sub1.ui"),
            actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary")
          ),
          conditionalPanel(
            condition = "input.sc1c2togL % 2 == 1",
            selectInput("sc1c2sub3", "",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp2),
            uiOutput("sc1c2sub3.ui"),
            actionButton("sc1c2sub3all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1c2sub3non", "Deselect all groups", class = "btn btn-primary")
          ), br(), br(),
          actionButton("sc1c2tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1c2tog % 2 == 1",
            radioButtons("sc1c2psz", "Plot size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Medium", inline = TRUE),
            radioButtons("sc1c2fsz", "Font size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Medium", inline = TRUE))
        ), # End of column (6 space)
        column(9, uiOutput("sc1c2oup.ui"),
               downloadButton("sc1c2oup.pdf", "Download PDF"),
               downloadButton("sc1c2oup.png", "Download PNG"), br(),
               div(style="display:inline-block",
                   numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px",
                                min = 4, max = 20, value = 8, step = 0.5)),
               div(style="display:inline-block",
                   numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px",
                                min = 4, max = 20, value = 10, step = 0.5))
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    ),     # End of tab (2 space)

    ### Tab gene signature # Written by Lucas

    tabPanel(
      HTML("Gene signature (boxplot)"),
      fluidRow(
        column(
          3, style="border-right: 2px solid black",
          selectInput("sc1d1_set", "Select gene signature",multiple = F,selectize=F,
                      choices = names(msigdbr_list), selected = "HALLMARK_INFLAMMATORY_RESPONSE") ,
          uiOutput("sc1d1_list.ui"),
          selectInput("sc1d1grp", "Group by:",
                      choices = c("Cancer_type_early_adv","Cell_subtype"),
                      selected = "Cancer_type_early_adv") %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Cell information to group cells by",
                   content = c("- If Cancer_type_stage is selected, the average is calculated for each sample and the plot is grouped based on their Cancer type.",
                               "- If Cell_subtype is selected, the average will be calculated for each subtype in each sample.")),
          actionButton("goButton_exp", "Generate / Update plot", class = "btn-success"),
          br(),
          br(),
          actionButton("sc1d1togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1d1togL % 2 == 1",
            selectInput("sc1d1sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1d1sub1.ui"),
            actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary")
          ),
          ### New toggle
          conditionalPanel(
            condition = "input.sc1d1togL % 2 == 1",
            selectInput("sc1d1sub3", "",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp2),
            uiOutput("sc1d1sub3.ui"),
            actionButton("sc1d1sub3all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1d1sub3non", "Deselect all groups", class = "btn btn-primary")
          ), br(), br(),
          actionButton("sc1d1tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1d1tog % 2 == 1",
            radioButtons("sc1d1psz", "Plot size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Small", inline = TRUE),
            radioButtons("sc1d1fsz", "Font size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Small", inline = TRUE),
            sliderInput("sc1d1siz", "Point size:",
                        min = 0, max = 4, value = 0.25, step = 0.25),
            checkboxInput("sc1d1lab1", "Show cell info labels", value = TRUE),
            radioButtons("sc1d1col1", "Colour:",
                         choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
                                     "Yellow-Green-Purple","Black-Violet-Yellow"),
                         selected = "Black-Violet-Yellow"))
        ), # End of column (6 space)
        column(9, h4(htmlOutput("sc1d1oupTxt")),
               uiOutput("sc1d1oup.ui"),
               downloadButton("sc1d1oup.pdf", "Download PDF"),
               downloadButton("sc1d1oup.png", "Download PNG"), br(),
               div(style="display:inline-block",
                   numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px",
                                min = 4, max = 20, value = 6, step = 0.5)),
               div(style="display:inline-block",
                   numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px",
                                min = 4, max = 20, value = 14, step = 0.5))
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    )      # End of tab (2 space)
    ,
    ### New tab gene signature with correlations # Written by Lucas

    tabPanel(
      HTML("Gene signature (correlation)"),
      h5("Due to the size of our datasets the correlation analysis can take up to 3 minutes."),
      fluidRow(
        column(
          3, style="border-right: 2px solid black",
          selectInput("sc1d2_set", "Select gene signature",multiple = F,selectize=F,
                      choices = names(msigdbr_list), selected = "HALLMARK_INFLAMMATORY_RESPONSE") ,
          uiOutput("sc1d2_list.ui"),
          # selectInput("sc1d2grp", "Group by:",
          #             choices = c("Cancer_type_stage","Cell_subtype"),
          #             selected = "Cancer_type_stage") %>%
          #   helper(type = "inline", size = "m", fade = TRUE,
          #          title = "Cell information to group cells by",
          #          content = c("- If Cancer_type_stage is selected, the average is calculated for each sample and the plot is grouped based on their Cancer type.",
          #                      "- If Cell_subtype is selected, the average will be calculated for each subtype in each sample.")),
          actionButton("goButton_correlation", "Generate / Update plot", class = "btn-success"),
          br(), br(),
          actionButton("sc1d2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1d2togL % 2 == 1",
            selectInput("sc1d2sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1d2sub1.ui"),
            actionButton("sc1d2sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1d2sub1non", "Deselect all groups", class = "btn btn-primary")
          ),
          #New toggle
          conditionalPanel(
            condition = "input.sc1d2togL % 2 == 1",
            selectInput("sc1d2sub3", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp2),
            uiOutput("sc1d2sub3.ui"),
            actionButton("sc1d2sub3all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1d2sub3non", "Deselect all groups", class = "btn btn-primary")
          ), br(), br(),
          actionButton("sc1d2tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1d2tog % 2 == 1",
            radioButtons("sc1d2psz", "Plot size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Small", inline = TRUE),
            radioButtons("sc1d2fsz", "Font size:",
                         choices = c("Small", "Medium", "Large"),
                         selected = "Small", inline = TRUE),
            sliderInput("sc1d2siz", "Point size:",
                        min = 0, max = 4, value = 0.25, step = 0.25),
            checkboxInput("sc1d2lab1", "Show cell info labels", value = TRUE),
            radioButtons("sc1d2col1", "Colour:",
                         choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
                                     "Yellow-Green-Purple","Black-Violet-Yellow"),
                         selected = "Black-Violet-Yellow"))
        ), # End of column (6 space)
        column(9, h4(htmlOutput("sc1d2oupTxt")),
               uiOutput("sc1d2oup.ui"),
               downloadButton("sc1d2oup.pdf", "Download PDF"),
               downloadButton("sc1d2oup.png", "Download PNG"), br(),
               div(style="display:inline-block",
                   numericInput("sc1d2oup.h", "PDF / PNG height:", width = "138px",
                                min = 4, max = 20, value = 6, step = 0.5)),
               div(style="display:inline-block",
                   numericInput("sc1d2oup.w", "PDF / PNG width:", width = "138px",
                                min = 4, max = 20, value = 14, step = 0.5))
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    )      # End of tab (2 space)
    ,

    ##################################
    br(),
    p("", style = "font-size: 125%;"),
    p(em("This webpage was made by Lucas Maciel using "), a("ShinyCell",
                                                            href = "https://github.com/SGDDNB/ShinyCell",target="_blank"), " as base code."),
    br(),br(),br(),br(),br()
  )))



