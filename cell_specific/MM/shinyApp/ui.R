library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
msigdbr_list <- readRDS("hallmarks.rds")
msigdbr_list

# url <- a("https://doi.org/10.1038/s41586-022-05242-7", href="https://doi.org/10.1038/s41586-022-05242-7")
# rdr <- a("Google Drive",href="https://drive.google.com/drive/folders/1poq4Lo5AxVp0WpG1EMgIjIeDR4q98zcA")
# Mesenchymal_like_genes <- c('AA',"BB")
### Start server code
shinyUI(fluidPage(
  ### HTML formatting of error messages

  tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))),
  list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 14px; }"))),


  ### Page title
  titlePanel(title=div(img(src="vib_logo.png", height = '50', width = '130'), HTML("Pan-Cancer Myeloid cells scRNA-Seq")),windowTitle ="Myeloid cells" ),
  navbarPage(
    NULL,
    # fluidRow(column(12,align="left",
    #                 h5("Reference: Karras  et al. 2022 doi:",url,". Seurat object available at ",rdr))),
    ### Tab1.a1: cellInfo vs geneExpr on dimRed
    tabPanel(
      HTML("Metadata vs Gene Expression"),
      h4("Cell information vs gene expression on reduced dimensions"),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              6, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                             selected = sc1def$dimred[1])),
            column(
              6,selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                            selected = sc1def$dimred[2]))
          )
        ), # End of column (6 space)
        column(
          3,br(),
          actionButton("sc1a1togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a1togL % 2 == 1",
            selectInput("sc1a1sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1a1sub1.ui"),
            actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary")
          )
        ), # End of column (6 space)
        column(
          6, br(),
          actionButton("sc1a1tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a1tog0 % 2 == 1",
            fluidRow(
              column(
                6, sliderInput("sc1a1siz", "Point size:",
                               min = 0, max = 4, value = 0.25, step = 0.25),
                radioButtons("sc1a1psz", "Plot size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                radioButtons("sc1a1fsz", "Font size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE)
              ),
              column(
                6, radioButtons("sc1a1asp", "Aspect ratio:",
                                choices = c("Square", "Fixed", "Free"),
                                selected = "Square", inline = TRUE),
                checkboxInput("sc1a1txt", "Show axis text", value = FALSE)
              )
            )
          )
        )  # End of column (6 space)
      ),   # End of fluidRow (4 space)
      fluidRow(
        column(
          6, style="border-right: 2px solid black",
          fluidRow(
            column(
              4, selectInput("sc1a1inp1", "Cell information:",
                             choices = sc1conf$UI,
                             selected = sc1def$meta1)
            ),
            column(
              6, br(),
              actionButton("sc1a1tog1", "Left-panel plot controls"),
              conditionalPanel(
                condition = "input.sc1a1tog1 % 2 == 1",
                radioButtons("sc1a1col1", "Colour (Continuous data):",
                             choices = c("Blue-White-Red","Grey-Blue","White-Red","Blue-Yellow-Red","Yellow-Green-Purple","Black-Violet-Yellow"),
                             selected = "Blue-Yellow-Red"),
                radioButtons("sc1a1ord1", "Plot order:",
                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                             selected = "Original", inline = TRUE),
                checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE)
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1a1oup1.ui"))),
          downloadButton("sc1a1oup1.pdf", "Download PDF"),
          downloadButton("sc1a1oup1.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 6, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 11, step = 0.5)),
        ), # End of column (6 space)
        column(
          6,
          fluidRow(
            column(
              4, selectInput("sc1a1inp2", "Type gene name:", choices=NULL) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c("Select gene to colour cells by gene expression",
                                   "- Type your gene of interest",
                                   paste0("- Gene expression are coloured in a ",
                                          "White-Red colour scheme which can be ",
                                          "changed in the plot controls\n")))
            ),
            column(
              6, br(),
              actionButton("sc1a1tog2", "Right-panel plot controls"),
              conditionalPanel(
                condition = "input.sc1a1tog2 % 2 == 1",
                radioButtons("sc1a1col2", "Colour:",
                             choices = c("Blue-White-Red","Grey-Blue","White-Red","Blue-Yellow-Red","Yellow-Green-Purple","Black-Violet-Yellow"),
                             selected = "White-Red"),
                radioButtons("sc1a1ord2", "Plot order:",
                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                             selected = "Max-1st", inline = TRUE)
              )
            )
          ) ,
          fluidRow(column(12, uiOutput("sc1a1oup2.ui"))),
          downloadButton("sc1a1oup2.pdf", "Download PDF"),
          downloadButton("sc1a1oup2.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 6, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 11, step = 0.5))
        )  # End of column (6 space)
      ),    # End of fluidRow (4 space)
      br(),
      #Create new fluidrow just to display the numbers
      fluidRow(actionButton("sc1a1tog9", "Show cell numbers / statistics",class = "btn-success"),align = "center",
               conditionalPanel(
                 condition = "input.sc1a1tog9 % 2 == 1",
                 h4("Cell numbers / statistics"),
                 dataTableOutput("sc1a1.dt")) )

    ),     # End of tab (2 space)

    ### Tab1.a2: cellInfo vs cellInfo on dimRed
    tabPanel(
      HTML("Metadata vs Metadata"),
      h4("Cell information vs cell information on dimension reduction"),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              6, selectInput("sc1a2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                             selected = sc1def$dimred[1])),
            column(
              6, selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                             selected = sc1def$dimred[2]))
          )
        ), # End of column (6 space)
        column(
          3, br(),
          actionButton("sc1a2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a2togL % 2 == 1",
            selectInput("sc1a2sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1a2sub1.ui"),
            actionButton("sc1a2sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1a2sub1non", "Deselect all groups", class = "btn btn-primary")
          )
        ), # End of column (6 space)
        column(
          6, br(),
          actionButton("sc1a2tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a2tog0 % 2 == 1",
            fluidRow(
              column(
                6, sliderInput("sc1a2siz", "Point size:",
                               min = 0, max = 4, value = 0.25, step = 0.25),
                radioButtons("sc1a2psz", "Plot size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                radioButtons("sc1a2fsz", "Font size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE)
              ),
              column(
                6, radioButtons("sc1a2asp", "Aspect ratio:",
                                choices = c("Square", "Fixed", "Free"),
                                selected = "Square", inline = TRUE),
                checkboxInput("sc1a2txt", "Show axis text", value = FALSE)
              )
            )
          )
        )  # End of column (6 space)
      ),   # End of fluidRow (4 space)
      fluidRow(
        column(
          6, style="border-right: 2px solid black",
          fluidRow(
            column(
              4, selectInput("sc1a2inp1", "Cell information 1:",
                             choices = sc1conf$UI,
                             selected = sc1def$meta1) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                       title = "Cell information to colour cells by",
                       content = c("Select cell information to colour cells",
                                   "- Categorical covariates have a fixed colour palette",
                                   paste0("- Continuous covariates are coloured in a ",
                                          "Blue-Yellow-Red colour scheme, which can be ",
                                          "changed in the plot controls")))
            ),
            column(
              6, br(),
              actionButton("sc1a2tog1", "Left-panel plot controls"),
              conditionalPanel(
                condition = "input.sc1a2tog1 % 2 == 1",
                radioButtons("sc1a2col1", "Colour (Continuous data):",
                             choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
                                         "Yellow-Green-Purple","Black-Violet-Yellow"),
                             selected = "Blue-Yellow-Red"),
                radioButtons("sc1a2ord1", "Plot order:",
                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                             selected = "Original", inline = TRUE),
                checkboxInput("sc1a2lab1", "Show cell info labels", value = TRUE)
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1a2oup1.ui"))),
          downloadButton("sc1a2oup1.pdf", "Download PDF"),
          downloadButton("sc1a2oup1.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1a2oup1.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 6, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1a2oup1.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 11, step = 0.5))
        ), # End of column (6 space)
        column(
          6,
          fluidRow(
            column(
              4, selectInput("sc1a2inp2", "Cell information 2:",
                             choices = sc1conf$UI,
                             selected = sc1def$meta2) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                       title = "Cell information to colour cells by",
                       content = c("Select cell information to colour cells",
                                   "- Categorical covariates have a fixed colour palette",
                                   paste0("- Continuous covariates are coloured in a ",
                                          "Blue-Yellow-Red colour scheme, which can be ",
                                          "changed in the plot controls")))
            ),
            column(
              6, br(),
              actionButton("sc1a2tog2", "Right-panel plot controls"),
              conditionalPanel(
                condition = "input.sc1a2tog2 % 2 == 1",
                radioButtons("sc1a2col2", "Colour (Continuous data):",
                             choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
                                         "Yellow-Green-Purple","Black-Violet-Yellow"),
                             selected = "Blue-Yellow-Red"),
                radioButtons("sc1a2ord2", "Plot order:",
                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                             selected = "Original", inline = TRUE),
                checkboxInput("sc1a2lab2", "Show cell info labels", value = TRUE)
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1a2oup2.ui"))),
          downloadButton("sc1a2oup2.pdf", "Download PDF"),
          downloadButton("sc1a2oup2.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 6, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 11, step = 0.5))
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    ),     # End of tab (2 space)

    ### Tab1.a3: geneExpr vs geneExpr on dimRed
    tabPanel(
      HTML("Gene Expression (umap)"),
      h4("Gene expression on dimension reduction"),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              6, selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                             selected = sc1def$dimred[1])),
            column(
              6,selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                            selected = sc1def$dimred[2]))
          )
        ), # End of column (6 space)
        column(
          3, br(),
          actionButton("sc1a3togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a3togL % 2 == 1",
            selectInput("sc1a3sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1a3sub1.ui"),
            actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary")
          )
        ), # End of column (6 space)
        column(
          6, br(),
          actionButton("sc1a3tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a3tog0 % 2 == 1",
            fluidRow(
              column(
                6, sliderInput("sc1a3siz", "Point size:",
                               min = 0, max = 4, value = 0.25, step = 0.25),
                radioButtons("sc1a3psz", "Plot size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                radioButtons("sc1a3fsz", "Font size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Small", inline = TRUE)
              ),
              column(
                6, radioButtons("sc1a3asp", "Aspect ratio:",
                                choices = c("Square", "Fixed", "Free"),
                                selected = "Square", inline = TRUE),
                checkboxInput("sc1a3txt", "Show axis text", value = FALSE) ,
                radioButtons("sc1a3col1", "Colour:",
                             choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
                                         "Yellow-Green-Purple","Black-Violet-Yellow"),
                             selected = "White-Red"),
                radioButtons("sc1a3ord1", "Plot order:",
                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                             selected = "Max-1st", inline = TRUE)
              )
            )
          )
        )  # End of column (6 space)
      ),   # End of fluidRow (4 space)
      fluidRow(
        column(
          12, style="border-right: 2px solid black",
          fluidRow(
            column(
              4, selectInput("sc1a3inp1", "Type gene names:", choices=NULL,multiple = TRUE) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                       title = "Gene expression to colour cells by",
                       content = c("Select gene to colour cells by gene expression",
                                   paste0("- Gene expression are coloured in a ",
                                          "White-Red colour scheme which can be ",
                                          "changed in the plot controls")))
            ),
            column(3,
                   br(),
                   actionButton("goButton_umap", "Generate / Update plot", class = "btn-success"),
            )
          ),
          fluidRow(column(12, uiOutput("sc1a3oup1.ui"))),
          downloadButton("sc1a3oup1.pdf", "Download PDF"),
          downloadButton("sc1a3oup1.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 6, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 15, step = 0.5))
        ) # End of column (6 space)
        # column(
        #   6,
        #   fluidRow(
        #     column(
        #       4, selectInput("sc1a3inp2", "Type gene name 2:", choices=NULL) %>%
        #         helper(type = "inline", size = "m", fade = TRUE,
        #                title = "Gene expression to colour cells by",
        #                content = c("Select gene to colour cells by gene expression",
        #                            paste0("- Gene expression are coloured in a ",
        #                                   "Grey-Blue colour scheme which can be ",
        #                                   "changed in the plot controls")))
        #     ),
        #     column(
        #       6, br(),
        #       actionButton("sc1a3tog2", "Right-panel plot controls"),
        #       conditionalPanel(
        #         condition = "input.sc1a3tog2 % 2 == 1",
        #         radioButtons("sc1a3col2", "Colour:",
        #                      choices = c("Blue-White-Red","Grey-Blue","White-Red", "Blue-Yellow-Red",
        #                                  "Yellow-Green-Purple","Black-Violet-Yellow"),
        #                      selected = "White-Red"),
        #         radioButtons("sc1a3ord2", "Plot order:",
        #                      choices = c("Max-1st", "Min-1st", "Original", "Random"),
        #                      selected = "Max-1st", inline = TRUE)
        #       )
        #     )
        #   ),
        #   fluidRow(column(12, uiOutput("sc1a3oup2.ui"))),
        #   downloadButton("sc1a3oup2.pdf", "Download PDF"),
        #   downloadButton("sc1a3oup2.png", "Download PNG"), br(),
        #   div(style="display:inline-block",
        #       numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "138px",
        #                    min = 4, max = 20, value = 6, step = 0.5)),
        #   div(style="display:inline-block",
        #       numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "138px",
        #                    min = 4, max = 20, value = 6, step = 0.5))
        # )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    ),     # End of tab (2 space)

    ### Tab1.b2: Gene coexpression plot
    tabPanel(
      HTML("Gene coexpression (umap)"),
      h4("Coexpression of two genes on reduced dimensions"),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              5, selectInput("sc1b2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                             selected = sc1def$dimred[1])),
            column(
              5,selectInput("sc1b2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                            selected = sc1def$dimred[2]))
          )
        ), # End of column (6 space)
        column(
          3, br(),
          actionButton("sc1b2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1b2togL % 2 == 1",
            selectInput("sc1b2sub1", "Cell information to subset:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1),
            uiOutput("sc1b2sub1.ui"),
            actionButton("sc1b2sub1all", "Select all groups", class = "btn btn-primary"),
            actionButton("sc1b2sub1non", "Deselect all groups", class = "btn btn-primary")
          )
        ), # End of column (6 space)
        column(
          6, br(),
          actionButton("sc1b2tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1b2tog0 % 2 == 1",
            fluidRow(
              column(
                6, sliderInput("sc1b2siz", "Point size:",
                               min = 0, max = 4, value = 0.25, step = 0.25),
                radioButtons("sc1b2psz", "Plot size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                radioButtons("sc1b2fsz", "Font size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Small", inline = TRUE)
              ),
              column(
                6, radioButtons("sc1b2asp", "Aspect ratio:",
                                choices = c("Square", "Fixed", "Free"),
                                selected = "Square", inline = TRUE),
                checkboxInput("sc1b2txt", "Show axis text", value = FALSE)
              )
            )
          )
        )  # End of column (6 space)
      ),   # End of fluidRow (4 space)
      fluidRow(
        column(
          3,
          selectInput("sc1b2inp1", "Gene 1:", choices=NULL) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Gene expression to colour cells by",
                   content = c("Select genes to colour cells by gene expression",
                               paste0("- Gene expression are coloured in a ",
                                      "White-Red colour scheme which can be ",
                                      "changed in the plot controls"))),
          selectInput("sc1b2inp2", "Gene 2:", choices=NULL) %>%
            helper(type = "inline", size = "m", fade = TRUE,
                   title = "Gene expression to colour cells by",
                   content = c("Select gene to colour cells by gene expression",
                               paste0("- Gene expression are coloured in a ",
                                      "White-Red colour scheme which can be ",
                                      "changed in the plot controls"))),
          actionButton("sc1b2tog1", "Toggle plot controls"),
          conditionalPanel(
            condition = "input.sc1b2tog1 % 2 == 1",
            radioButtons("sc1b2col1", "Colour:",
                         choices = c("Red (Gene1); Blue (Gene2)",
                                     "Orange (Gene1); Blue (Gene2)",
                                     "Red (Gene1); Green (Gene2)",
                                     "Green (Gene1); Blue (Gene2)"),
                         selected = "Red (Gene1); Blue (Gene2)"),
            radioButtons("sc1b2ord1", "Plot order:",
                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                         selected = "Max-1st", inline = TRUE)
          )
        ), # End of column (6 space)
        column(
          6, style="border-right: 2px solid black",
          uiOutput("sc1b2oup1.ui"),
          downloadButton("sc1b2oup1.pdf", "Download PDF"),
          downloadButton("sc1b2oup1.png", "Download PNG"), br(),
          div(style="display:inline-block",
              numericInput("sc1b2oup1.h", "PDF / PNG height:", width = "138px",
                           min = 4, max = 20, value = 7, step = 0.5)),
          div(style="display:inline-block",
              numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px",
                           min = 4, max = 20, value = 7, step = 0.5))
        ), # End of column (6 space)
        column(
          3, uiOutput("sc1b2oup2.ui"),
          downloadButton("sc1b2oup2.pdf", "Download PDF"),
          downloadButton("sc1b2oup2.png", "Download PNG"),
          br(), h4("Cell numbers"),
          dataTableOutput("sc1b2.dt")
        )  # End of column (6 space)
      )    # End of fluidRow (4 space)
    ),     # End of tab (2 space)

    ### Tab gene signature # Written by Lucas

    tabPanel(
      HTML("Gene signature (boxplot/umap)"),
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

    ##################################
    br(),
    p("", style = "font-size: 125%;"),
    p(em("This webpage was made by Lucas Maciel using "), a("ShinyCell",
                                                            href = "https://github.com/SGDDNB/ShinyCell",target="_blank"), " as base code."),
    br(),br(),br(),br(),br()
  )))



