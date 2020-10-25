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

library(DT) # R interface to js data table library for displaying an interative table
library(data.table) # faster than data.frame data structure
library(reshape2) 
library(plotly) 
library(dplyr)
library(shinyjs)
library(shinyWidgets)
library(shinyBS) # tooltip
library(shinyhelper) # help button

source ("utils.R") #triggers for the interface
source ("ggmaplot2.R") #new maplot with yellow color for NA data
source ("plots.R")


#### LOADING MA DATA 


# EXPECTED COLUMN NAMES: 
#geneName, baseMean, log2FoldChange, pAdj, (detectionCall)


handle<-NULL

#### HANDLE GLOBAL VARIABLES
PvalLIM=10^-8

initVar<-function(handle){
  #list of global state variables 
  handle$genesToTrack<-NULL # name sof genes to track in all plots
  
  handle$MAdataCur<-NULL ## store current MA data to use as input in MA interface
  handle$MAindToTrack<-NULL ## store index to highlight in MA plot of current MA data
  handle$MAindToHighlight<-NULL ## This only highlights one gene (Ali)
  
  
  
  handle$userID <- NULL
  handle$fdrVal <- 0.05 #value from the fdr slider
  handle$maGeneBrushSelInfo <- NULL #the info (rows of df displayed - geneName baseMean log2FoldChange pAdj detection_call) from the brush selected in MA plot
  handle$selectedMAdata<- NULL # rows of MA data whose genes are currently selected
  handle$selectedGenes<- NULL  # names of the currently selected genes
  
 
  
  handle$maColor <- c("#B31B21", "#1465AC", "darkgray", "yellow","green") # for the cases (up, down, non-significant, not available,tracked) of MA Data
  handle$defaultBarColor <- "#000000"
  handle$maColorLabels <- c("Up-regulated", "Down-regulated", "Not significant", "N/A", "Tracked")
  handle$defaultBarColorLabel <- "Gene not in MA Data"
  
  my_rgb <- function(r,g,b) { rgb(r,g,b,maxColorValue = 255) }
  handle$legend <- data.frame(x=0, y=9:1, 
                              colour=c(handle$maColor[5], handle$maColor[4], 
                                       my_rgb(192,0,0), my_rgb(255,51,51), my_rgb(255,127,127), "darkgray",
                                       my_rgb(91,156,213), my_rgb(46,117,182), my_rgb(32,78,121)),
                              label=c("Tracked", "N/A", "Up-high", "Up-mid", "Up-low", "Not significant", "Down-low", "Down-mid", "Down-high"), 
                              stringsAsFactors = FALSE)
  
  handle$maPlotHighlight <- NULL #gene name to highlight on MA Plot based on the selected gene in hi-C Bar Chart
  
  handle$ranges <-NULL
  handle$geneColor<-NULL
  
  return(handle)
}

## RESET ALL UI
initUI<-function(session)
{
  updateTextAreaInput(session, inputId="selectedGenesByLasso", label = NULL, value = "")
  updateTextAreaInput(session, inputId="selectedGenes_postFilter", label = NULL, value = "")
  updateSelectInput(session, inputId='selectGenesByName', choices=character(0),selected=character(0))
  updateTextAreaInput(session, inputId="genesToTrack", label = NULL, value = "")
  
  html(id="label_selected_genes", paste("Lasso Selection ( 0 )"))
  html(id="label_selectedGenes_preFilter", "Count = 0")
}

## RESET FILTER SELECTIONS AREA
initSELECT<-function(session)
{
  updateTextAreaInput(session, inputId="selectedGenesByLasso", label = NULL, value = "")
  updateTextAreaInput(session, inputId="selectedGenes_postFilter", label = NULL, value = "")
  updateSelectInput(session, inputId='selectGenesByName', choices=character(0),selected=character(0))
  
  html(id="label_selected_genes", paste("Lasso Selection ( 0 )"))
  html(id="label_selectedGenes_preFilter", "Count = 0")
}





# INIT DATA and data-dependent variables (called when loading data)
initData<-function(session,handle,maData){
  if(is.null(maData)){
    
  }else{
    handle$MAdataCur<-maData ## store current MA data to use as input in MA interface
    if (is.null(handle$MAdataCur$detectionCall)){
      handle$MAdataCur$detectionCall=1
    }
    
    # convert to capital letters all gene names
    handle$MAdataCur$geneName<-toupper(handle$MAdataCur$geneName)
    
    # pre-computing for mapping the up/down to -1 to 1 range. (-ve for down/blue and +ve for red/up genes)
    handle$MAdataCur$pFold <- handle$MAdataCur$pAdj*(handle$MAdataCur$log2FoldChange/abs(handle$MAdataCur$log2FoldChange))
    
    
    # axes ranges to set the cut-off filter slider ranges
    handle$ranges <- list(X=c(floor(min(log2(handle$MAdataCur$baseMean),2)), ceiling(max(log2(handle$MAdataCur$baseMean),2))), 
                          Y=c(floor(min(handle$MAdataCur$log2FoldChange, 2)), ceiling(max(handle$MAdataCur$log2FoldChange, 2))))
    
    
    handle$geneColor<-assignColorToGene(handle, maData) #assign color to gene names based on fdr significance value in MA Data
  
    ### update ranges of sliders in the Filter panel
    updateSliderInput(session, "filter_slider_cutOffX", 
                      min=handle$ranges$X[1], 
                      max=handle$ranges$X[2],
                      value=handle$ranges$X)
    
    updateSliderInput(session, "filter_slider_cutOffY", 
                      min=0,
                      max=ceiling(max(abs(handle$ranges$Y[1]),abs(handle$ranges$Y[2]))),
                      value=0)
    
    ### Populate the list of all genes in Search list
    updateSelectInput(session, inputId='selectGenesByName', choices=handle$MAdataCur[,"geneName"])
    
  }
  return(handle)
}

# call init variables
handle<-initVar(handle)


#data log function
handle$log <- data.frame(timestamp = "", event = "", stringsAsFactors = FALSE)

updateLog <- function(handle, event) {
  handle$log <<- rbind(handle$log, c(as.character(Sys.time()), event))
}

cat("Setting up the application...\n")


notesFile <<- 'notes.R'
if(file.exists(notesFile)){
  local({ 
    load(file=notesFile)
    handle$notes <<- notes
  })
}

onStop(function() {
  print(handle$notes)
  notes <- handle$notes
  save(notes, file=notesFile)
  cat("Writing notes to the file\n")
})