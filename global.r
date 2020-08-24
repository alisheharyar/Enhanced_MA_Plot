#library(shiny)
#library(ggplot2) 
#library(HiCcompare) #for the hic plot; to download it, we used setRepositories
library(DT) # R interface to js data table library for displaying an interative table
library(data.table) # faster than data.frame data structure
#library(ggpubr) #ggmaplot
#library(hicutils) #to change it to a matrix but i think we're not using it now
library(reshape2) 
library(plotly) 
#library(RColorBrewer)
#library(heatmaply) #alternative way to draw the heatmap
library(dplyr)
#library(shinyBS) # for tooltip and popup
#library(cetcolor) # color palette

#library(enrichplot)
#library(DOSE)
#BiocManager::install("enrichplot")

library(shinyjs)
library(shinyWidgets)
#library(grid)
#library(ggnewscale)
#library(circlize)

source ("utils.R") #triggers for the interface
source ("ggmaplot2.R") #new maplot with yellow color for NA data
#source ("hiCDataProcessing.R") #data processing to plot hic data
source ("plots.R")



### FAKE DATA FOR HISTORY SANKEY
datSK <- data.frame(From=c(rep("A",3), rep("B", 3)),
                    To=c(rep(c("X", "Y", "Z"),2)),
                    Weight=c(5,7,6,2,9,4))


#### LOADING MA DATA (DEBUG)

load("../Data/maDataClean5Dina.RData") # maDataClean5
maData<-as.data.table(maDataClean5[1:500,])
# EXPECTED COLUMN NAMES: baseMean log2FoldChange padj detection_call 
# COLUMN NAMES AS LOADED:
# c("GeneID","Base mean","log2(FC)","StdErr","Wald-Stats","P-value","P-adj",
#    "GeneIDshort","refseq_mrna","ensembl_gene_id",
#    "chromosome_name","start_position","end_position","external_gene_name","shortchr")

colnames(maData)<-c("GeneID","baseMean","log2FoldChange","StdErr","Wald-Stats","P-value","padj",
                    "GeneIDshort","refseq_mrna","ensembl_gene_id",
                    "chromosome_name","start_position","end_position","name","shortchr")
# add $detection_call 
maData$detection_call=1

#### LOADING HIC DATA (DEBUG)


handle<-NULL


handle$CHR1cur<-1
handle$CHR2cur<-1

handle$pathData1<-"../Data/"
handle$prefixData1<-"Dina"
handle$pathData2<-"../Data/"
handle$prefixData2<-"GM12878"

#### HANDLE GLOBAL VARIABLES
#list of global constant variables 


#list of global state variables 

handle$genesToTrack<-NULL # name sof genes to track in all plots

handle$MAdataCur<-maData ## store current MA data to use as input in MA tab
handle$MAindToTrack<-NULL ## store index to highlight in MA plot of current MA data
handle$MAindToHighlight<-NULL ## This only highlights one gene (Ali)

# pre-computing for mapping the up/down to -1 to 1 range. (-ve for down/blue and +ve for red/up genes)
handle$MAdataCur$pFold <- handle$MAdataCur$padj*(handle$MAdataCur$log2FoldChange/abs(handle$MAdataCur$log2FoldChange))

names(handle$MAContColors) <- levels(handle$MAdataCur$pFold_factor)



handle$userID <- NULL
handle$fdrVal <- 0.05 #value from the fdr slider
handle$maGeneBrushSelInfo <- NULL #the info (rows of df displayed - name baseMean log2FoldChange padj detection_call) from the brush selected in MA plot
handle$selectedMAdata<- NULL # rows of MA data whose genes are currently selected
handle$selectedGenes<- NULL  # names of the currently selected genes


handle$maColor <- c("#B31B21", "#1465AC", "darkgray", "yellow","green") # for the cases (up, down, non-significant, not available,tracked) of MA Data
handle$defaultBarColor <- "#000000"
handle$maColorLabels <- c("Up-regulated", "Down-regulated", "Not significant", "N/A", "Tracked")
handle$defaultBarColorLabel <- "Gene not in MA Data"

my_rgb <- function(r,g,b) { rgb(r,g,b,maxColorValue = 255) }
handle$legend <- data.frame(x=0, y=9:1, 
                            colour=c(handle$maColor[5], handle$maColor[4], my_rgb(192,0,0), my_rgb(255,51,51), my_rgb(255,127,127), "darkgray",
                                     my_rgb(91,156,213), my_rgb(46,117,182), my_rgb(32,78,121)),
                            label=c("Tracked", "N/A", "Up-high", "Up-mid", "Up-low", "Not significant", "Down-low", "Down-mid", "Down-high"), 
                            #colour=c(handle$maColor, handle$defaultBarColor), 
                            #label=c(handle$maColorLabels, handle$defaultBarColorLabel), 
                            stringsAsFactors = FALSE)

handle$maPlotHighlight <- NULL #gene name to highlight on MA Plot based on the selected gene in hi-C Bar Chart

# axes ranges to set the cut-off filter slider ranges
handle$ranges <- list(X=c(round(min(log2(handle$MAdataCur$baseMean)),2), round(max(log2(handle$MAdataCur$baseMean)),2)), 
                      Y=c(round(min(handle$MAdataCur$log2FoldChange), 2), round(max(handle$MAdataCur$log2FoldChange), 2))
                      )

#initialization
handle$geneColor=assignColorToGene(handle, maData) #assign color to gene names based on fdr significance value in MA Data

#data log function
handle$log <- data.frame(timestamp = "", event = "", stringsAsFactors = FALSE)

updateLog <- function(handle, event) {
  handle$log <<- rbind(handle$log, c(as.character(Sys.time()), event))
}

cat("Doing application setup\n")
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