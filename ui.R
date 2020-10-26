# Copyright: Ali Sheharyar (Texas AM University at Qatar), Michael Aupetit (Qatar Computing Research Institute)
# October 25, 2020
# Code Version 2
# This file is part of "Enhanced MA plot"
# 
# "Enhanced MA plot" is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version. (GPL-3 or later)
# 
# "Enhanced MA plot" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with "Enhanced MA plot" in the "COPYING" file  If not, see <https://www.gnu.org/licenses/>.
#
# Please cite the Github the code as: 
# https://github.com/alisheharyar/Enhanced_MA_Plot

### SHINY UI
ui <-  fluidPage(

  
  tags$head(tags$style(HTML('
    #selectedGenesByLasso{
        //color: #fc6600;
    }
    
    #selectedGenes_postFilter{                        
      color: orange;
    }
    
    #genesToTrack {
      color: LimeGreen;
    }
    
    #selectedGenes_preFilter{                        
      //color: #03ac13;
    }
    
    .selectize-input {
      height: 180px; 
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
    
    #filter_val_down {
    background-color:rgb(91,156,213)
    }
    #filter_val_notsig {
    background-color:darkgray
    }
    #filter_val_up {
    background-color:rgb(255,127,127)
    }

    
    '))),
  
  fluidRow(
    useShinyjs(),
    column(12, h2(id='title', "Enhanded MA Plot")),
    column(12,
           column(2,
                  h4("Open MA data")%>% 
                    helper(type = "inline",
                           title = "User manual",
                           content = c("Open an MA data .csv or .RData file clicking on \"Browse...\" ",
                                       "The .RData file must contain a data frame named \"MAdata\".",
                                       "The \"MAdata\" data frame or the table in the .csv file must contain columns",
                                       "named \"geneName\", \"baseMean\", \"log2FoldChange\", \"pAdj\" typical output from DESeq2 pipeline.",
                                       "Other columns will be ignored, but will be copied in the exported data on \"Save filtered MA data\".",
                                       " ",
                                       "You can also use example test data by clicking the \"Test Data\" button",
                                       "When a checked icon appears on the \"Test Data\" button, that means the interface displays the example test data",
                                       " ",
                                       "The MA plot displays each gene as a point, with the following:",
                                       "x-axis: log2(basemean+1), with basemean=R*G",
                                       "y-axis: log2FoldChange, with  log2FoldChange=log2(Fc)=log2(R)-log(G)", 
                                       "color: Bluish for negative and Reddish for positive, grey for non-significant difference",
                                       "where R and G are the signal intensities or read counts of a given gene",
                                       "in the two conditions under study.",
                                       " ",
                                       "Other functions are given as tooltips on hover/? of the corresponding widgets.",
                                       "",
                                       "Credits:",
                                       "Copyright: Dr. Ali Sheharyar (Texas AM University at Qatar), Dr. Michael Aupetit (Qatar Computing Research Institute)",
                                       "October 25, 2020. Code Version 2. Released under license GPL-3 or later",
                                       "Design study conducted by Mrs. Talar Boghos Yacoubian (MSc HBKU)",
                                       "Supervisors: Dr. Dena Al Thani (HBKU CSE), Dr. Michael Aupetit (QCRI)",
                                       "Expert users: Ms. Dina Aljogol (HBKU CHLS) and Pr. Borbala Mifsud (HBKU CHLS)",
                                       "Please refer/cite the code as: ",
                                       "https://github.com/alisheharyar/Enhanced_MA_Plot"
                                       ),
                           size = "l"),
                  fileInput("loadData", label=NULL,
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv",
                              ".RData")
                  ),
                  fluidRow(
                    column(6,actionButton(inputId="buttonLoadTestData", label="Test Data...", width = '100%')),
                    column(6,actionButton(inputId="buttonResetUI", label="Reset...", width = '100%'))),
                  bsTooltip("buttonLoadTestData", "Checked means test data are displayed.",
                            "right", trigger = "hover"),
                  bsTooltip("buttonResetUI", "Reset all selections and tracked genes.",
                            "right", trigger = "hover"),
                  br(),
                  br(),
                  br(),
                  h4(id="label_pvalue","P-value Cut-off (FDR)"),
                  sliderTextInput(inputId="fdr", 
                                  label="Pre-defined",
                                  grid = TRUE,
                                  force_edges = TRUE,
                                  choices = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                                  selected = handle$fdrVal),
                  numericInput(inputId="fdr_txt", 
                               label="Manual ]0,1[", 
                               value=handle$fdrVal, 
                               min = 0, max = 1, step = 0.001),
           ),
           column(1,
                  br(),
                  br(),
                  br(),
                  plotOutput("legendPlot", height = 200, width = 100),
                  br(),
                  h4("Figure"),
                  downloadButton(outputId="buttonSaveMAPlotPNG", label=".png"),
                  bsTooltip("buttonSaveMAPlotPNG", "Save the plot as a .png",
                            "right", trigger = "hover"),
                  
                  downloadButton(outputId="buttonSaveMAPlotRDATA", label=".RData"),
                  bsTooltip("buttonSaveMAPlotRDATA", "Save the plot as \"MAplot\" ggplot plot into a .RData file. Load the plot MAplot_zzz.RData file then try \"> MAplot\" then \"> MAplot+coord_fixed()\" in the RStudio console.",
                            "right", trigger = "hover")
           ),
           column(7,
                  plotlyOutput("maPlot", height = 400),
                  br()
            ),
           column(2,
                  h4(id="label_track_genes", "Tracked genes (0)"),
                  disabled(textAreaInput(inputId="genesToTrack", label=NULL,
                                         value="", placeholder="Search/Lasso->Filter->Track to see Genes here...",
                                         width="100%", rows="13")),
                  actionButton(inputId="buttonClearTrackedGenes", label="Clear Tracked", width = '100%')
                  )
    ),
    column(12,
          column(4,
                 column(6,
                        h4("Search Genes"),
                        # Ref for selection event: https://stackoverflow.com/questions/50168069/r-shiny-selectize-selected-event)
                        # Ref for copy/pasting into a selectize: https://github.com/rstudio/shiny/issues/1663
                        selectizeInput(inputId='selectGenesByName', label=NULL, choices=NULL, multiple=TRUE,
                                       options = list(
                                         splitOn = I("(function() { return /[,; ]/; })()"), #allow for copy-paste
                                         create = I("function(input, callback){return {value: input,text: input};}"),
                                         render = I("{item: function(item, escape) {return '<div class=\"item\" onclick=\"Shiny.onInputChange(\\\'selectGenesByName_click\\\', \\\'' + escape(item.value) + '\\\')\">' + escape(item.value) + '</div>';}}")
                                       ))%>% 
                          helper(type = "inline",
                                 title = "Gene search by names",
                                 content = c("Copy-paste gene names here.",
                                             "Those not available in the loaded MA data will be indicated and can be copied.",
                                             "Those available can be filtered and tracked..."),
                                 size = "m")
                 ),
                 column(6,
                        h4(id="label_selected_genes", "Lassoed Genes (0)"),
                        disabled(textAreaInput(inputId="selectedGenesByLasso", label=NULL,
                                               value="",placeholder="Genes selected by lasso/box in the plot will appear here... (Double-click the plot to clear)",
                                               width="100%", rows="9"))
                 ),
                 br(),
                 column(12,
                        radioGroupButtons(
                          inputId='filterKeep', choices = c('Keep all', 'Keep singles', 'Keep multiples'), selected = c('Keep all'), status = 'default',
                          justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"), individual=T) 
                        )
                 )
          ),
          column(4,
                 h4("Set Filters")%>% 
                   helper(type = "inline",
                          title = "Set Filters",
                          content = c("Apply filters by checking boxes and set values.",
                                      "All filters are applied together (logical AND).", 
                                      "Filtered genes are highlighted in orange in the MA plot."),
                          size = "m"),
                 fluidRow(
                        column(4, h5(div("Log2(Mean Expr.)"))),
                        column(5, sliderInput("filter_slider_cutOffX", label=NULL, step=0.1, min=0, max=100, value=c(0,100))),
                        column(2, checkboxInput(inputId="filter_chk_cutOffX_reverse", label="Reverse interval"))
                  ),
                 fluidRow(
                        column(4, h5(div("Log2(Fold Change)"))),
                        column(5, sliderInput("filter_slider_cutOffY", label=NULL, step=1, min=0, max=16, value=0)),
                        column(2, checkboxInput(inputId="filter_chk_cutOffY_reverse", label="Reverse interval"))
                  ),
                 fluidRow(
                        column(8, h5(id="TopK_genes","No Top/Bottom rank filter by P-value if 0")
                               %>% 
                                 helper(type = "inline",
                                        title = "Top/Bottom-K Genes by P-value",
                                        content = c("No Top/Bottom rank filter on P-value if K=0",
                                                    "Select the |K| most significant genes (lowest P-value) if K>0.",
                                                    "Select the |K| least significant genes (highest P-value) if K<0."),
                                        size = "m")),
                        column(4, numericInput(inputId='filter_val_topK', label=NULL, value = 0))
                        ),
                 fluidRow(
                        column(4, 
                               bsButton(
                                 inputId = "filter_val_down", 
                                 label = "Down -",
                                 type = "toggle",
                                 block=TRUE,
                                 value = TRUE, 
                                 icon =  icon("ok",lib = "glyphicon"))
                        ),
                        column(4, 
                               bsButton(
                                 inputId = "filter_val_notsig", 
                                 label = "Not Sig.",
                                 type = "toggle",
                                 block=TRUE,
                                 value = TRUE, 
                                 icon =  character(0))
                        ),
                        column(4, 
                               bsButton(
                                 inputId = "filter_val_up", 
                                 label = "Up +",
                                 type = "toggle",
                                 block=TRUE,
                                 value = TRUE, 
                                 icon =  icon("ok",lib = "glyphicon"))
                        )
                  )
          ),
          column(2,
                   h4(id="label_selectedGenes_postFilter", "Filtered Genes (0)"),
                   disabled(textAreaInput(inputId="selectedGenes_postFilter", label=NULL,
                                          value="",placeholder="Genes after the filters will appear here...",
                                          rows="9")),
                 actionButton(inputId="buttonClearSelectedGenes", label="Clear Selections", width = '100%')
          ),
          column(2,
                 h4("Track Filtered Genes"),
                 actionButton(inputId="buttonTrackSelectedGenes", label="Append to Tracked", width = '100%'),
                 h4("Notes"),
                 textAreaInput(inputId="notes", label=NULL,
                               value=handle$notes,placeholder="Write yobservations here... (Saved only with .RData)",
                               width="100%", rows="4"),
                 h4("Save Filtered MA data"),
                 fluidRow(
                    column(4, actionButton(inputId="buttonTableView", label="Details...")),
                    column(4, offset=0, downloadButton(outputId="buttonDownloadGenesCSV", label=".csv")),
                    column(4, offset=0, downloadButton(outputId="buttonDownloadGenesRDATA", label=".RData")),
                    bsTooltip("buttonTableView", "Show a Table view of the filtered data.",
                              "right", trigger = "hover"),
                    bsTooltip("buttonDownloadGenesCSV", "Save the filtered data as a CSV file.",
                              "left", trigger = "hover"),
                    bsTooltip("buttonDownloadGenesRDATA", "Save a .RData file containing the filtered data in MAdata dataframe, the plot in MAplot ggplot, and the notes as a text in MAnotes.",
                              "left", trigger = "hover")
                 )
          )
    )
  ) # END FLUID ROW
)

