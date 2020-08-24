ui <-  fluidPage(
  
  #tags$head(
  #  tags$style(HTML("
  #
  #    .selectize-input {
  #      height: 250px;
  #      font-size: 24pt;
  #      padding-top: 5px;
  #    }
  #
  #  "))),
  
  tags$head(tags$style(HTML('
    #selectedGenes{
        //color: #fc6600;
    }
    
    #selectedGenes_postFilter{                        
      //color: #fc6600;
    }
    
    #genesToTrack {
      //color: #fc6600;
    }
    
    #selectedGenes_preFilter{                        
      //color: #03ac13;
    }
    
    .selectize-input {
      height: 215px; 
    }
    
    .selectize-input > div {
      //color: #03ac13 !important;
    }
  
    [for="fdr"] {
      color : #777777;
    }
    
    [for="fdr_txt"] {
      color : #777777;
    }
    
    #label_pvalue {
      color : #777777;
    }
    
    '))),
  
  fluidRow(
    useShinyjs(),
    column(12, h2(id='title', "Enhanded MA Plot -- maDataClean.RData")),
    column(12,
           column(2,
                  h4("Load..."),
                  actionButton(inputId="loadData", label="Browse files..."),
                  br(),
                  br(),
                  h4("Select P-value"),
                  sliderTextInput(inputId="fdr", 
                                  label="Pre-defined selection",
                                  grid = TRUE,
                                  force_edges = TRUE,
                                  choices = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                                  selected = handle$fdrVal),
                  numericInput(inputId="fdr_txt", 
                               label="Manual input (0-1):", 
                               value=handle$fdrVal, 
                               min = 0, max = 1, step = 0.01),
                  h5(id="label_pvalue")
           ),
           column(1,
                  br(),
                  br(),
                  br(),
                  br(),
                  plotOutput("legendPlot", height = 200, width = 100),
                  #plotOutput("PFoldlegendPlot", height = 200),
                  br(),
                  br(),
                  downloadButton(outputId="buttonSaveMAPlot", label="Figure"),
                  hidden(actionButton(inputId="buttonMAplotClearHighlight", label="Clear Selected"))
           ),
           column(7,
                  plotlyOutput("maPlot", height = 400),
                  br(),
            ),
           column(2,
                  h4(id="label_track_genes", "Tracked genes"),
                  disabled(textAreaInput(inputId="genesToTrack", label=NULL,
                                         value="", placeholder="Genes to track in all plots here...",
                                         width="100%", rows="13")),
                  actionButton(inputId="buttonEmptySelectedGenes", label="Reset", width = '100%'))
    ),
    column(12,
          column(4,
                 column(6,
                        h4("Search Genes"),
                        # Ref for selection event (https://stackoverflow.com/questions/50168069/r-shiny-selectize-selected-event)
                        selectizeInput(inputId='selectGenesByName', label=NULL, choices=c('A','B','C','D'), multiple=TRUE,
                                       options = list(
                                         render = I("{item: function(item, escape) {return '<div class=\"item\" onclick=\"Shiny.onInputChange(\\\'selectGenesByName_click\\\', \\\'' + escape(item.value) + '\\\')\">' + escape(item.value) + '</div>';}}")
                                       ))
                 ),
                 column(6,
                        h4(id="label_selected_genes", "Selected genes..."),
                        disabled(textAreaInput(inputId="selectedGenes", label=NULL,
                                               value="",placeholder="Genes selected in the tab will appear here...",
                                               width="100%", rows="10"))
                 ),
                 br(),
                 column(12,
                        radioGroupButtons(
                          inputId='filterOccurence', choices = c('Keep all', 'Keep singles', 'Keep multiples'), selected = c('Keep all'), status = 'default',
                          justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"), individual=T) 
                        )
                 )
          ),
          column(4,
                 h4("Filters"),
                 fluidRow(
                        column(4, checkboxInput(inputId = "filter_chk_cutOffX", label = div("Cut-off at", br(), "(X-axis)"))),
                        column(5, sliderInput("filter_slider_cutOffX", label=NULL, step=0.1, min=0, max=100, value=c(0,100))),
                        column(2, checkboxInput(inputId="filter_chk_cutOffX_reverse", label="Reverse interval"))
                  ),
                 fluidRow(
                        column(4, checkboxInput(inputId = "filter_chk_cutOffY", label = div("Cut-off at", br(), "(Y-axis)"))),
                        column(5, sliderInput("filter_slider_cutOffY", label=NULL, step=0.1, min=0, max=100, value=c(0,100))),
                        column(2, checkboxInput(inputId="filter_chk_cutOffY_reverse", label="Reverse interval"))
                  ),
                 fluidRow(
                        column(6, checkboxInput(inputId="filter_chk_topK", label="Top Genes By P-value")),
                        column(6, numericInput(inputId='filter_val_topK', label=NULL, value = 5))
                        ),
                 fluidRow(
                        column(4, checkboxInput(inputId="filter_chk_ctg", label="Categories")),
                        column(8, 
                          checkboxGroupButtons(
                             inputId = "filter_val_ctg", label = NULL,
                             choices = c("Up", "Down"),
                             selected = c("Up", "Down"), 
                             justified = TRUE, status = "default",
                             checkIcon = list(yes = icon("ok",lib = "glyphicon"),individual=T))
                           )
                  )
                   #column(12, 
                  #        column(4, checkboxInput(inputId = "filter_chk_topK", label = "Top K")),
                  #        column(8, numericInput(inputId='filter_val_topK', label=NULL, value = 5))
                  #        ),
                  # column(6,
                  #        checkboxInput(inputId = "filter_chk_ctg", label = "Categories"),
                  #        checkboxGroupButtons(
                  #          inputId = "filter_val_ctg", label = NULL,
                  #          choices = c("Up", "Down"),
                  #          selected = c("Up", "Down"), 
                  #          justified = TRUE, status = "default",
                  #          checkIcon = list(yes = icon("ok",lib = "glyphicon"),individual=T)
                  #        )
                  # ),
                  # column(6,
                  #        checkboxInput(inputId = "filter_chk_cutOffX", label = "Cut-off (X)"),
                  #        #textInput(inputId='filter_val_cutOffX_A', label=NULL, placeholder ="A"),
                  #        #textInput(inputId='filter_val_cutOffX_B', label=NULL, placeholder ="B"),
                  #        
                  #        checkboxInput(inputId = "filter_chk_cutOffY", label = "Cut-off at y-axis"),
                  #        textInput(inputId='filter_val_cutOffY_A', label=NULL, placeholder ="A"),
                  #        textInput(inputId='filter_val_cutOffY_B', label=NULL, placeholder ="B"),
                  # )
          ),
          column(2,
                   h4(id="label_selectedGenes_postFilter", "Selection after filters (0)"),
                   disabled(textAreaInput(inputId="selectedGenes_postFilter", label=NULL,
                                          value="",placeholder="Genes after the filters will appear here...",
                                          rows="12"))
          ),
          column(2,
                 actionButton(inputId="buttonTrackSelectedGenes", label="Track Selected Genes", width = '100%'),
                 h4("Notes"),
                 textAreaInput(inputId="notes", label=NULL,
                               value=handle$notes,placeholder="Write your notes and observations here...",
                               width="100%", rows="8"),
                 fluidRow(
                    column(6, actionButton(inputId="buttonTableView", label="Show Table"),),
                    column(6, offset=0, downloadButton(outputId="buttonDownloadGenes", label="Data"))
                    )
          )
    ),
    column(2,

           #### RIGHT MENU
           # checkboxInput(inputId="checkboxFavorite", label="Add to important tabs", value = FALSE, width = "100%"),
           # textAreaInput(inputId="userNotes", label=NULL,
           #               value="",placeholder="Take your notes here...",
           #               width="100%", rows="10"),
           # actionButton(inputId="buttonEmptyUserNotes", label="Reset"),
           hidden(
             h6(id="label_selectedGenes_preFilter", "Count = 0"),
             textAreaInput(inputId="selectedGenes_preFilter", label=NULL,
                           value="",placeholder="Select the genes through lasso or by name and add here",
                           rows="15")
           )
    ) #### END RIGHT MENU
  ) # END FLUID ROW
)

