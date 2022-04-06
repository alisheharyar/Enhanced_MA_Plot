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

server <- function(input, output,session) {
  
  
  
  ###########################
  ###########################
  ######## TRIGGERS #########
  ###########################
  ###########################
  
  observe_helpers(withMathJax = TRUE) # enable help buttons
  
  #### FLush out all reactive
  
  #output triggers - triggers that will appear when we want to plot/print
  trigMAplot<-makeReactiveTrigger() #to refresh ma plot
  
  #input triggers - triggers that will appead when we want to save information from updated input
  trigMAselected<-makeReactiveTrigger() #to save the gene info brushed in the MA Plot
  
  trigMAcore<- makeReactiveTrigger()
  
   
  ############################################
  ############ LEFT SIDE PANEL ###############
  ############################################
  
  # LOADING DATA
  # MA data must come as table .CSV or as a "MAdata" dataframe saved in .RData file. 
  # Columns must be named: geneName, baseMean, log2FoldChange, pAdj
  observeEvent(input$loadData,{
    inFile <- input$loadData
    
    if (is.null(inFile))
    {
      print("NO DATA LOADED")
      return(NULL)
    }else{
      print("DATA LOADED")
        fileExt<-unlist(strsplit(inFile$datapath,"[.]"))[2]
        if (fileExt=="csv"){
          maData<-as.data.table(read.csv(inFile$datapath, header = TRUE))
        }else{
          load(inFile$datapath)
          maData<-as.data.table(MAdata)
        }
        # initialize the interface with the new data
        handle<<-initVar(handle)
        handle<<-initData(session,handle,maData)
        
        if (is.null(handle$MAdataCur)){
          resetAll()
          # replot
          trigMAplot$trigger()
        }else{
          initUI(session)
          
          updateActionButton(session,inputId="buttonLoadTestData", icon = character(0))
          # replot
          trigMAplot$trigger()
        }
      }
  })

  ### LOAD TEST DATA
  observeEvent(input$buttonLoadTestData,{
      showModal(modalDialog(
        tagList(), 
        title="This will reset all selections and replace currently loaded data, do you want to continue?",
        footer = tagList(actionButton("confirmLoadTest", "Yes, load test data"),
                         modalButton("No, cancel")
        )
      ))
  })
  observeEvent(input$confirmLoadTest,{  
    inFile <- input$confirmLoadTest
    if (is.null(inFile))
    {
    }else{
      print("DATA LOADED")
      load("MAdata.RData")
      maData<-as.data.table(MAdata)
      
      # initialize the interface with the new data
      handle<<-initVar(handle)
      handle<<-initData(session,handle,maData)
      initUI(session)
      
      updateActionButton(session,inputId="buttonLoadTestData", icon =  icon("ok",lib = "glyphicon"))
      
      # replot
      trigMAplot$trigger()
    }
    # Remove modal dialog
    removeModal()
  })  
  
  
  ## RESET UI
  observeEvent(input$buttonResetUI,{
    showModal(modalDialog(
      tagList(), 
      title="This will reset all selected and tracked genes, do you want to continue?",
      footer = tagList(actionButton("confirmResetUI", "Yes, reset"),
                       modalButton("No, cancel")
      )
    ))
  })
observeEvent(input$confirmResetUI,{ # empty the list
    
  if (!is.null(input$confirmResetUI))
  {
    resetAll()
    # replot
    trigMAplot$trigger()
  }
  # Remove modal dialog
  removeModal()
})
  
resetAll = function(){
  handle$selectedMAdata<<-NULL
  handle$selectedGenes<<-NULL
  
  handle<<-initVar(handle)
  handle<<-initData(session,handle,NULL)
  initUI(session)
  updateActionButton(session,inputId="buttonLoadTestData", icon = character(0))
  
}


observeEvent(input$genesToTrack,{
   
    handle$genesToTrack<<-cleanStrGenesToTrack(input$genesToTrack)

    # trigger plots
    trigMAcore$trigger()
    
    html(id="label_track_genes", paste("Tracked Genes (",length(handle$genesToTrack),")"))
    
  })

  
  observeEvent(input$buttonTrackSelectedGenes,{ # refresh plots
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, 
                        value = implode(sort(unique(c(cleanStrGenesToTrack(input$genesToTrack),
                                                      cleanStrGenesToTrack(input$selectedGenes_postFilter)))),sep=" "))
  })
  
  ########################################
  ############ MA PLOT TAB ###############
  ########################################

  #to include the legend in the plots
  output$legendPlot <- renderPlot({
    plotLegend(handle)
  })
  
  output$PFoldlegendPlot <- renderPlot({
    plotPFoldLegend(handle)
  })

  ###save the slider data in fdr when it is updated and trigger the update of the maplot
  

  ### Slider Input of P-value cut-off
  observeEvent(input$fdr, {

    #updateLog(handle, "FDR Slider Moved")

    handle$fdrVal <<- isolate(input$fdr)
    #update the color gene column in colorGene and update handle$geneColor
    handle$geneColor <<- assignColorToGene(handle, handle$MAdataCur)
    
    # update the manual input
    updateTextInput(session, "fdr_txt", value=handle$fdrVal)
    
    html(id="label_pvalue", paste("P-value Cut-off (FDR) =",handle$fdrVal))
    
    trigMAcore$trigger()
  })
  
  ### Manual Text Input of P-value cut-off
  observeEvent(input$fdr_txt, {
    handle$fdrVal <<- isolate(min(max(PvalLIM,input$fdr_txt),1-PvalLIM))
    handle$geneColor <<- assignColorToGene(handle, handle$MAdataCur)
    
    # update the manual input (needed when going out of range)
    updateTextInput(session, "fdr_txt", value=handle$fdrVal)
    
    html(id="label_pvalue", paste("P-value Cut-off (FDR) =",handle$fdrVal))
    
    trigMAcore$trigger()
  })

  
  #### MA CORE MANAGER
  # all input trigger the core
  # all output are triggered by the core
  observe({
    trigMAcore$depend()

    # update genes to highlight remanent in current MA data
    handle$MAindToTrack<<-which(is.element(handle$MAdataCur$geneName,handle$genesToTrack))
    # update genes to highlight transient in current MA data (green dots)
    handle$MAindToHighlight<<-which(is.element(handle$MAdataCur$geneName,handle$selectedGenes))
    
    trigMAplot$trigger()
  })

  output$maPlot <- renderPlotly({
    #the rendering of the maplot will be triggered by the trigMaPlot triggers
    trigMAplot$depend()
    
    if (!is.null(handle$MAdataCur)){
      
      filterX <- NULL
      filterY <- NULL
      
        filterX1 <-  c(handle$ranges$X[1], input$filter_slider_cutOffX[1])
        filterX2 <- c(input$filter_slider_cutOffX[2], handle$ranges$X[2])
        filterX <- list(R1=filterX1, R2=filterX2, internal=!input$filter_chk_cutOffX_reverse)
     
        filterY1 = c(handle$ranges$Y[1], -input$filter_slider_cutOffY[1])
        filterY2 = c(input$filter_slider_cutOffY[1], handle$ranges$Y[2])
        filterY <- list(R1=filterY1, R2=filterY2, internal=input$filter_chk_cutOffY_reverse)
      
       g <- plotMA(handle, showHighlight=TRUE, title=NULL, discrete=F, 
                  filterX=filterX, filterY=filterY)
      
      handle$curPlot <<- g
      
      p <- ggplotly(g, source="maPlot", tooltip=c('text'))
      
      p <- layout(p, dragmode = "select")
      event_register(p, "plotly_selected")
      
    }else{
      p<-ggplot()+annotate("text", x = 0, y = 0, label = "Please load some MA data to start...")+theme_void()
    }
  })
  
  # Save the plot as PNG format
  output$buttonSaveMAPlotPNG<-downloadHandler(
    filename=function(){
      paste("MAplot_",Sys.Date(),".png",sep="")
    },
    content=function(file) {
      p <- ggplot2::last_plot()
      #save(p, file=file)
      ggsave(file=file, device="png")
    }
  )
  
  # Save the plot as RData format
  output$buttonSaveMAPlotRDATA<-downloadHandler(
    filename=function(){
      paste("MAplot_",Sys.Date(),".RData",sep="")
    },
    content=function(file) {
      MAplot <- ggplot2::last_plot()
      save(MAplot, file=file)
    }
  )
  
  
  
  ## LASSO SELECTION: capture the brush event on maplot and copy data in handle 
observe({
  # updateLog(handle, "Data Brushed on MA Plot")

  brushData <- event_data("plotly_selected", source = "maPlot")

  if(length(brushData)>0){

    handle$maGeneBrushSelInfo <<- brushData

    ind = which(is.element(round(log2(handle$MAdataCur$baseMean + 1),6), round(brushData$x,6)) & is.element(round(handle$MAdataCur$log2FoldChange,6), round(brushData$y,6)))

    handle$maGeneBrushSelInfo <<- handle$MAdataCur[ind,]
  }else{
    handle$maGeneBrushSelInfo <<- NULL
  }

    trigMAselected$trigger()
  })


  #update Selected genes textArea
  observe({
    trigMAselected$depend()

    if(is.null(handle$maGeneBrushSelInfo)) {
      value="Genes selected by lasso/box in the plot will appear here... (double-click the plot to clear)"
      # if(!is.null(handle$MAindToHighlight)) {
      #   handle$MAindToHighlight <<- NULL
      #   trigMAplot$trigger()
      # }
    }else{
      value <- handle$maGeneBrushSelInfo$geneName
    }

    updateTextAreaInput(session, inputId="selectedGenesByLasso", label = NULL, value = implode(sort(value), "  "))
    html(id="label_selected_genes", paste0("Lassoed Genes (",length(value),")"))
  })

  ## SAVE NOTES
  observeEvent(input$notes, {
    handle$notes <<- input$notes
  })
  
  ## Empty list of selected genes if button Clear Selections is pressed
  observe({ # empty the list
    input$buttonClearSelectedGenes
    
    handle$selectedMAdata<<-NULL
    handle$selectedGenes<<-NULL
    initSELECT(session)
  })
  ## Empty list of tracked genes if button Clear Tracked is pressed
  observe({ # empty the list
    input$buttonClearTrackedGenes
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = "")
  })

  
  ### POPUP TABLE VIEW OF MA DATA OF CURRENT SELECTED GENES
  observeEvent(input$buttonTableView,{
    showModal(modalDialog(DT::renderDataTable(as.data.frame(handle$selectedMAdata),options=list(scrollX=400)),
                            title="Table View"))
  })
  
  ### DOWNLOAD BUTTON OF MA DATA OF CURRENT SELECTED GENES - csv format
  output$buttonDownloadGenesCSV<-downloadHandler(
    filename=function(){
      paste0("MAdataSelected_",Sys.Date(),".csv")
    },
    content=function(file){
      write.csv(handle$selectedMAdata, file=file, row.names = FALSE) 
    })

  ### DOWNLOAD BUTTON OF MA TEST DATA - csv format
  output$buttonDownloadTestDataCSV<-downloadHandler(
    filename=function(){
      paste0("MAdataTEST.csv")
    },
    content=function(file){
      load("MAdata.RData")
      maDataTEST<-as.data.table(MAdata)
      
      write.csv(maDataTEST, file=file, row.names = FALSE) 
    })
  
  
  ### DOWNLOAD BUTTON OF MA DATA OF CURRENT SELECTED GENES, PLOTS AND NOTES - RData format
  output$buttonDownloadGenesRDATA<-downloadHandler(
    filename=function(){
      paste0("MAdataSelected_",Sys.Date(),".RData")
    },
    content=function(file){
      MAdata<-handle$selectedMAdata
      MAplot <- ggplot2::last_plot()
      MAnotes <- handle$notes
      save(MAdata, MAplot, MAnotes, MAnotes, file=file) 
    })
  

  # COPY-PASTING GENES IN SEARCH BOX BY GENE NAMES
  observe({
    ## Put all gene names in capital letters
    allSelectedGenes<-toupper(input$selectGenesByName)
    
    if(!is.null(allSelectedGenes)){
      ## check if gene is in available MA data
      isin<-is.element(allSelectedGenes,handle$MAdataCur$geneName)
      
      ## Display list of genes NOT in the MA data
      if (any(!isin))
      {
        showModal(modalDialog(
          title = "Warning!",
          paste0("These genes are not available in the current MA data, they will be removed: ",toString(allSelectedGenes[!isin]))
        ))
      }
      
      ## Update the list with only available genes
      if (any(isin))
      {
        updateSelectInput(session, inputId='selectGenesByName', choices=handle$MAdataCur[,"geneName"],
                        selected=allSelectedGenes[isin])
      }else{
        updateSelectInput(session, inputId='selectGenesByName', choices=handle$MAdataCur[,"geneName"],
                          selected=character(0))
      }
    }
    
    })
  
  # FILTERING
  observe({
    input$filter_val_topK # trigger on change of checkbox of topK genes
    input$filterKeep  # trigger on change of boolean mixing of genes-by-names and genes-by-lasso sets
    input$selectedGenesByLasso # trigger on lasso selection
    input$selectGenesByName # trigger on change of search genes by names
    input$fdr  #trigger on change of fdr P-value cut-off
    input$filter_chk_cutOffX_reverse
    input$filter_chk_cutOffY_reverse
    
    
    # Update displayed text of Filter top/bottom K genes
    if (input$filter_val_topK!=0)
    {
      if (input$filter_val_topK>0)
        html(id="TopK_genes", paste("Keep ",input$filter_val_topK," most significant (lowest P-value)"))
      if (input$filter_val_topK<0)
        html(id="TopK_genes", paste("Keep ",abs(input$filter_val_topK)," least significant (highest P-value)"))
    }else{
      html(id="TopK_genes", paste("No Top/Bottom filter by P-value"))
    }
    
    # update toggle buttons
    if(input$filter_val_up)
      updateButton(session, inputId="filter_val_up", icon = icon("ok",lib = "glyphicon"))
    else
      updateButton(session, inputId="filter_val_up", icon = character(0))
    
    if(input$filter_val_notsig)
      updateButton(session, inputId="filter_val_notsig", icon = icon("ok",lib = "glyphicon"))
    else
      updateButton(session, inputId="filter_val_notsig", icon = character(0))
    
    if(input$filter_val_down)
      updateButton(session, inputId="filter_val_down", icon = icon("ok",lib = "glyphicon"))
    else
      updateButton(session, inputId="filter_val_down", icon = character(0))
    
    
    #### IF DATA LOADED
    if (!is.null(handle$MAdataCur)){
      
    
    # make a data frame from the selected genes
    gnames <- c(input$selectGenesByName, cleanStrGenesToTrack(input$selectedGenesByLasso))
    
    
    # filter based on occurence
    if(input$filterKeep == 'Keep all') {
      # don't do anything
    } 
    else if(input$filterKeep == 'Keep singles') {
      gnames <- unique(gnames[!is.element(gnames,gnames[duplicated(gnames)])])
    } 
    else if(input$filterKeep == 'Keep multiples') {
      gnames <- unique(gnames[duplicated(gnames)])
    }
    
    maDataSelected <- handle$MAdataCur %>% filter(geneName %in% gnames)
    
    
    
    # Filter on x-cutoff
      A <- input$filter_slider_cutOffX[1]
      B <- input$filter_slider_cutOffX[2]
      
      maDataSelected$log2mean <- log2(maDataSelected$baseMean + 1)
      
      # A<B: internal interval
      if(!input$filter_chk_cutOffX_reverse) { 
        maDataSelected <- maDataSelected %>% filter(A <= log2mean & log2mean <= B)
      }
      # A>B: external interval
      else {
        maDataSelected <- maDataSelected %>% filter(log2mean < A | B < log2mean)
      }

    
    # Filter on y-cutoff
      A <- -(input$filter_slider_cutOffY)
      B <- (input$filter_slider_cutOffY)
      
      #  internal interval
      if(input$filter_chk_cutOffY_reverse) { 
        maDataSelected <- maDataSelected %>% filter(A <= log2FoldChange & log2FoldChange <= B)
      }
      #  external interval
      else {
        maDataSelected <- maDataSelected %>% filter(log2FoldChange < A | B < log2FoldChange)
      }


      
    # Filter on red/grey/blue colors
      maDataSelected$sig = setMAColor(handle$fdrVal, maDataSelected)
      maDataSelected1 <- NULL
      maDataSelected2 <- NULL
      maDataSelected3 <- NULL
      
      if(input$filter_val_up)
        maDataSelected1 <- maDataSelected %>% filter(sig == 1)
      if(input$filter_val_down)
        maDataSelected2 <- maDataSelected %>% filter(sig == 2)
      if(input$filter_val_notsig)
        maDataSelected3 <- maDataSelected %>% filter(sig == 3)
      
      maDataSelected <- rbind(maDataSelected1, maDataSelected2, maDataSelected3) 
    
      # Filter top/bottom K genes
      if (input$filter_val_topK!=0)
      {
        maDataSelected <- maDataSelected %>% top_n(-input$filter_val_topK, wt=pAdj)
      }
      
    
    sel <- maDataSelected$geneName
    
    updateTextAreaInput(session, inputId="selectedGenes_postFilter", label = NULL, value = implode(sort(sel), "  "))
    html(id="label_selectedGenes_postFilter", paste("Filtered Genes (",length(sel), ")"))
    
    handle$selectedMAdata<<-maDataSelected
    handle$selectedGenes<<-sel
    
    trigMAcore$trigger()
    }
  })
}