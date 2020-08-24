server <- function(input, output,session) {
  
  
  
  ###########################
  ###########################
  ######## TRIGGERS #########
  ###########################
  ###########################
  
  #### FLush out all reactive
  
  #output triggers - triggers that will appear when we want to plot/print
  trigMAplot<-makeReactiveTrigger() #to refresh ma plot
  
  #input triggers - triggers that will appead when we want to save information from updated input
  trigMAselected<-makeReactiveTrigger() #to save the gene info brushed in the MA Plot
  
  trigMAcore<- makeReactiveTrigger()
  
  updateSliderInput(session, "filter_slider_cutOffX", 
                    min=handle$ranges$X[1], 
                    max=handle$ranges$X[2],
                    value=handle$ranges$X)
  
  updateSliderInput(session, "filter_slider_cutOffY", 
                    min=handle$ranges$Y[1], 
                    max=handle$ranges$Y[2],
                    value=handle$ranges$Y)
  
  
  ############################################
  ############ LEFT SIDE PANEL ###############
  ############################################
  observeEvent(input$loadData, {
    showModal(modalDialog(fileInput(inputId="loadMA", label="MA data"),
      title = "Load data...", "Here to load data"
    ))
  })

  observeEvent(input$genesToTrack,{
   
    handle$genesToTrack<<-cleanStrGenesToTrack(input$genesToTrack)

    # trigger plots
    trigMAcore$trigger()
    
    html(id="label_track_genes", paste("Tracked Genes (",length(handle$genesToTrack),")"))
    
  })

  observeEvent(input$buttonUniqueGeneTrack,{ # return all genes without duplicates
    strGeneToTrack<-implode(unique(cleanStrGenesToTrack(input$genesToTrack)),sep=",")
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = strGeneToTrack)
  })
  observeEvent(input$buttonMultiGeneTrack,{ # return only the genes who are there multiple times
    cln<-cleanStrGenesToTrack(input$genesToTrack)
    strGeneToTrack<-implode(unique(cln[duplicated(cln)]),sep=",")
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = strGeneToTrack )
  })
  observeEvent(input$buttonSingleGeneTrack,{ # return only the genes who are there a single time
    cln<-cleanStrGenesToTrack(input$genesToTrack)
    strGeneToTrack<-implode(unique(cln[!is.element(cln,cln[duplicated(cln)])]),sep=",")
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = strGeneToTrack )
  })
  observeEvent(input$buttonEmptyGeneTrack,{ # empty the list
    handle$genesToTrack<<-NULL
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = "")
  })
  
  observeEvent(input$buttonTrackSelectedGenes,{ # refresh plots
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, 
                        value = implode(sort(c(cleanStrGenesToTrack(input$genesToTrack),cleanStrGenesToTrack(input$selectedGenes_postFilter))),sep="  "))
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
  

  observeEvent(input$fdr, {

    #updateLog(handle, "FDR Slider Moved")

    handle$fdrVal <<- isolate(input$fdr)
    #update the color gene column in colorGene and update handle$geneColor
    handle$geneColor <<- assignColorToGene(handle, maData)
    
    # update the manual input
    updateTextInput(session, "fdr_txt", value=handle$fdrVal)
    
    html(id="label_pvalue", paste("Current P-value = ",handle$fdrVal,""))
    
    trigMAcore$trigger()
  })
  
  observeEvent(input$fdr_txt, {
    handle$fdrVal <<- isolate(input$fdr_txt)
    handle$geneColor <<- assignColorToGene(handle, maData)
    
    html(id="label_pvalue", paste("Current P-value = ",handle$fdrVal,""))
    
    trigMAcore$trigger()
  })

  
  #### MA CORE MANAGER
  # all input trigger the core
  # all output are triggered by the core
  observe({
    trigMAcore$depend()

    # update genes to highlight in current MA data
    handle$MAindToTrack<<-which(is.element(handle$MAdataCur$name,handle$genesToTrack))

    trigMAplot$trigger()
  })

  output$maPlot <- renderPlotly({
    #the rendering of the maplot will be triggered by the trigMaPlot triggers
    trigMAplot$depend()
    input$filter_chk_cutOffX
    
    filterX <- NULL
    filterY <- NULL
    
    if(input$filter_chk_cutOffX) {
      filterX1 <-  c(handle$ranges$X[1], input$filter_slider_cutOffX[1])
      filterX2 <- c(input$filter_slider_cutOffX[2], handle$ranges$X[2])
      filterX <- list(R1=filterX1, R2=filterX2, internal=!input$filter_chk_cutOffX_reverse)
    }
    
    if(input$filter_chk_cutOffY) {
      filterY1 = c(handle$ranges$Y[1], input$filter_slider_cutOffY[1])
      filterY2 = c(input$filter_slider_cutOffY[2], handle$ranges$Y[2])
      filterY <- list(R1=filterY1, R2=filterY2, internal=!input$filter_chk_cutOffY_reverse)
    }
    
    #g <- plotMA(handle, showHighlight=FALSE, title=NULL, discrete=(input$fdr_scale=='Discrete'))
    g <- plotMA(handle, showHighlight=FALSE, title=NULL, discrete=F, 
                filterX=filterX, filterY=filterY)
    
    handle$curPlot <<- g

    p <- ggplotly(g, source="maPlot", tooltip=c('text'))
    
    p <- layout(p, dragmode = "select")
    event_register(p, "plotly_selected")
    
    #print(p)
  })
  
  # Save the plot
  output$buttonSaveMAPlot<-downloadHandler(
    filename=function(){
      paste("enhanced_MA_plot_",Sys.Date(),".png",sep="")
    },
    content=function(file) {
      p <- ggplot2::last_plot()
      #save(p, file=file)
      ggsave(file=file, device="png")
    }
  )
  
  #capture the brush event on maplot and copy data in handle 
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

  observe({
    
    d <- event_data("plotly_click", source = "maPlot")
  
    handle$MAindToHighlight <<- NULL
    trigMAplot$trigger()
  })

  #update Selected genes textArea
  observe({
    trigMAselected$depend()

    if(is.null(handle$maGeneBrushSelInfo)) {
      value="Click and drag events in MA Plot (double click to clear)"
      if(!is.null(handle$MAindToHighlight)) {
        handle$MAindToHighlight <<- NULL
        trigMAplot$trigger()
      }
    }else{
      #value=c(handle$selectedGenes,handle$maGeneBrushSelInfo$name)
      #handle$selectedMAdata<<-rbind(handle$selectedMAdata,handle$maGeneBrushSelInfo)
      #handle$selectedGenes<<-value
      value <- handle$maGeneBrushSelInfo$name
    }

    updateTextAreaInput(session, inputId="selectedGenes", label = NULL, value = implode(sort(value), "  "))
    html(id="label_selected_genes", paste0("Lasso Selection (",length(value),")"))
  })

  observeEvent(input$notes, {
    handle$notes <<- input$notes
  })
  
  ## Empty list of selected genes if tab changes or if button reset is pressed
  observe({ # empty the list
    input$tabs
    input$buttonEmptySelectedGenes
    
    handle$selectedHIdata<<-NULL
    handle$selectedMAdata<<-NULL
    handle$selectedGenes<<-NULL
    updateTextAreaInput(session, inputId="selectedGenes", label = NULL, value = "")
    updateTextAreaInput(session, inputId="selectedGenes_preFilter", label = NULL, value = "")
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = "")
    html(id="label_selected_genes", paste("Lasso Selection ( 0 )"))
    html(id="label_selectedGenes_preFilter", "Count = 0")
  })
  
  ### RESET USER NOTES
  observeEvent(input$buttonEmptyUserNotes,{
    updateTextAreaInput(session, inputId="userNotes", label = NULL, value = "")
  })
  
  ### POPUP TABLE VIEW OF MA DATA OF CURRENT SELECTED GENES
  observeEvent(input$buttonTableView,{
    showModal(modalDialog(DT::renderDataTable(as.data.frame(handle$selectedMAdata),options=list(scrollX=400)),
                            title="Table View"))
  })
  
  ### DOWNLOAD BUTTON OF MA DATA OF CURRENT SELECTED GENES
  output$buttonDownloadGenes<-downloadHandler(
    filename=function(){
      paste("MA_data_",Sys.Date(),".R",sep="")
      },
    content=function(file){
      #write.csv(handle$selectedMAdata,file)
      save(handle, file=file)
      }
  )

  ### Populate the list of all genes
  updateSelectInput(session, inputId='selectGenesByName', choices=maDataClean5[,"external_gene_name"])
  
  ### Handler for the 'Track Genes' button for the 'Selection By Name' list. 
  observeEvent(input$buttonTrackSelectedGenesByName, {
    updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = implode(sort(c(cleanStrGenesToTrack(input$genesToTrack),cleanStrGenesToTrack(implode(input$selectGenesByName)))),sep=","))
  }, ignoreNULL=FALSE)
  
  
  observeEvent(input$buttonMAplotClearHighlight, {
    handle$MAindToHighlight <<- NULL
    trigMAplot$trigger()
  })
  
  observe({
    #trigMAcore$depend()
    
    # update genes to highlight in current MA data
    handle$MAindToHighlight<<-which(is.element(handle$MAdataCur$name,input$selectGenesByName_click))
    print(handle$MAindToHighlight)
    
    trigMAplot$trigger()
  })
 
  # Filtering
  observe({
    input$filter_chk_topK
    input$filter_chk_ctg
    input$filter_chk_cutOffX
    input$filter_chk_cutOffY
    #input$input$selectedGenes_preFilter
    input$filterOccurence
    input$selectedGenes
    input$selectGenesByName
    
    # make a data frame from the selected genes
    gnames <- c(input$selectGenesByName, cleanStrGenesToTrack(input$selectedGenes))
    
    
    # filter based on occurence
    if(input$filterOccurence == 'Keep all') {
      # don't do anything
    } 
    else if(input$filterOccurence == 'Keep singles') {
      gnames <- unique(gnames[!is.element(gnames,gnames[duplicated(gnames)])])
    } 
    else if(input$filterOccurence == 'Keep multiples') {
      gnames <- unique(gnames[duplicated(gnames)])
    }
    
    maDataSelected <- handle$MAdataCur %>% filter(name %in% gnames)
    
    # Filter on x-cutoff
    if(input$filter_chk_cutOffX) {
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
    }
    
    # Filter on y-cutoff
    if(input$filter_chk_cutOffY) {
      A <- input$filter_slider_cutOffY[1]
      B <- input$filter_slider_cutOffY[2]
      
      if(input$filter_chk_cutOffY_reverse) {
        B <- input$filter_slider_cutOffY[1]
        A <- input$filter_slider_cutOffY[2]
      }
      
      # A<B: internal interval
      if(!input$filter_chk_cutOffY_reverse) { 
        maDataSelected <- maDataSelected %>% filter(A <= log2FoldChange & log2FoldChange <= B)
      }
      # A>B: external interval
      else {
        maDataSelected <- maDataSelected %>% filter(log2FoldChange < A | B < log2FoldChange)
      }
    }

    # Filter top K genes
    if(input$filter_chk_topK)
      maDataSelected <- maDataSelected %>% top_n(-input$filter_val_topK, wt=padj)
    
    # Filter on red/blue colors
    if(input$filter_chk_ctg) {
      maDataSelected$sig = setMAColor(handle$fdrVal, maDataSelected)
      
      maDataSelected1 <- NULL
      maDataSelected2 <- NULL
      if('Up' %in% input$filter_val_ctg)
        maDataSelected1 <- maDataSelected %>% filter(sig == 1)
      if('Down' %in% input$filter_val_ctg)
        maDataSelected2 <- maDataSelected %>% filter(sig == 2)
      
      maDataSelected <- rbind(maDataSelected1, maDataSelected2) 
    }
    
    sel <- maDataSelected$name
    
    updateTextAreaInput(session, inputId="selectedGenes_postFilter", label = NULL, value = implode(sort(sel), "  "))
    html(id="label_selectedGenes_postFilter", paste("Selected Genes (",length(sel), ")"))
    
    handle$selectedMAdata<<-maDataSelected
    handle$selectedGenes<<-sel
  })
}