


###########################
###########################
##### Tooltip PLOT #########
###########################
###########################

#### try out this code for improved interactivity
#### https://stackoverflow.com/questions/48597530/how-to-change-a-plot-when-hovering-over-elements-in-shiny

showTooltipPanel<-function(hover,myTooltip)
{
  # calculate point position INSIDE the image as percent of total dimensions
  # from left (horizontal) and from top (vertical)
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  # calculate distance from left and bottom side of the picture in pixels
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  # create style property fot tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                  "left:", (1.2*left_px) , "px; top:", (1.2*top_px), "px;")
  
    
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(myTooltip))
  )
}


###########################
###########################
##### Legend PLOT #########
###########################
###########################

plotLegend = function (handle) {
  
  ggplot() + annotate("rect", xmin=0, xmax=1, ymin=1, ymax=9, fill="white")+
    annotate("rect", ymin=handle$legend$y, ymax=handle$legend$y+1, xmin=0, xmax=0.15, fill=handle$legend$colour)+
    annotate("text", x= handle$legend$x+0.2, y= handle$legend$y+0.5, label=handle$legend$label, size=4, colour="black", hjust=0) +
    theme_void()
}

# to just extract the legend, used the cowplot package
# ref: https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
plotPFoldLegend = function(handle) {
  data <- handle$MAdataCur
  fdr <- handle$fdrVal
  T <- fdr
  
  data <- data.frame(log2mean=log2(data$baseMean + 1),
                     log2FC = data$log2FoldChange, 
                     padj = data$padj, 
                     pFold = data$padj*(data$log2FoldChange/abs(data$log2FoldChange)),
                     pFoldDummy = cut(handle$MAdataCur$pFold, 
                                      breaks=c(-1, -T-2*(1-T)/3*2, -T-(1-T)/3 ,-T, 0, T, T+(1-T)/3,T+2*(1-T)/3, 1), 
                                      labels = c("Down-gray", "Down-low", "Down-mid", "Down-high", 
                                                 "Up-high", "Up-mid", "Up-low", "Up-gray")),
                     stringsAsFactors = FALSE)
  
  # reorder the factor levels
  # ref: http://www.cookbook-r.com/Manipulating_data/Changing_the_order_of_levels_of_a_factor/
  data$pFoldDummy <- factor(data$pFoldDummy, levels=c("Up-high", "Up-mid", "Up-low", "Up-gray", 
                                                      "Down-gray", "Down-low", "Down-mid", "Down-high"))
  
  my_rgb <- function(r,g,b) { rgb(r,g,b,maxColorValue = 255) }
  
  cols <- c( "Up-gray"=my_rgb(166,166,166), "Up-low"=my_rgb(255,127,127), "Up-mid"=my_rgb(255,51,51), "Up-high"=my_rgb(192,0,0), 
             "Down-gray"=my_rgb(166,166,166), "Down-high"=my_rgb(32,78,121), "Down-mid"=my_rgb(46,117,182), "Down-low"=my_rgb(91,156,213))
  
  p <- ggplot()
  p <- p + geom_point(data=data, 
                      aes(x = log2mean, y = log2FC, fill=pFoldDummy, colour=pFoldDummy), 
                      size = 7, pch = 22) + 
    labs(colour='P-value', fill='P-value') + 
    scale_color_manual(values=cols, guide=F) +
    scale_fill_manual(values=cols) +
    theme(
      legend.key.size = unit(0, 'lines'),
      legend.spacing = unit(0, 'cm')
      )
    
  legend <- cowplot::get_legend(p)
  
  grid.newpage()
  grid.draw(legend)
}

###########################
###########################
######## MA PLOT ##########
###########################
###########################

plotMA = function(handle, showHighlight=FALSE, title=NULL, discrete=T, filterX, filterY) {
  
  p <- ggplot()
  
  
  if(!is.null(filterX)) {
    if(filterX$internal) {
      p <- p + annotate("rect", fill="gray",
               xmin=filterX$R1[1]-500, xmax=filterX$R1[2],  
               ymin=-500, ymax=500)
      
      p <- p + annotate("rect", fill="gray",
               xmin=filterX$R2[1], xmax=filterX$R2[2]+500,  
               ymin=-500, ymax=500)
    } else {
      p <- p + annotate("rect", fill="gray",
                        xmin=filterX$R1[2], xmax=filterX$R2[1],  
                        ymin=-500, ymax=500)
    }
  }
  
  if(!is.null(filterY)) {
    if(filterY$internal) {
      p <- p + annotate("rect", fill="gray",
                        xmin=-500, xmax=500,  
                        ymin=filterY$R1[1]-500, ymax=filterY$R1[2])
      
      p <- p + annotate("rect", fill="gray",
                        xmin=-500, xmax=500,  
                        ymin=filterY$R2[1], ymax=filterY$R2[2]+500)
    } else {
      p <- p + annotate("rect", fill="gray",
                        xmin=-500, xmax=500,  
                        ymin=filterY$R1[2], ymax=filterY$R2[1])
    }
  }
  
  p <- p + coord_cartesian(xlim=handle$ranges$X, ylim = handle$ranges$Y)
  p <- ggmaplot2 ( data=handle$MAdataCur, indToTrack=handle$MAindToTrack, indToHighlight=handle$MAindToHighlight  ,fdr = handle$fdrVal, fc = 1.5, 
                  genenames = handle$MAdataCur$name, detection_call = NULL, sizepts = 0.4,                       
                  font.label = c(12, "plain", "black"), 
                  font.legend = "bold", font.main = "bold", 
                  label.rectangle = FALSE, 
                  palette = handle$maColor, 
                  top = 15, select.top.method = "padj", 
                  main = title, 
                  xlab = "Log2 mean expression", ylab = "Log2 fold change", 
                  ggtheme = theme_classic(),legend = "top", 
                  discrete=discrete, p)
  
  if (showHighlight) 
    p <- p +  
      geom_point(data=maData[which(is.element(maData$name, unique(handle$maPlotHighlight$gene))),], 
                 aes(x=log2(baseMean+1), y=log2FoldChange), color="#000000", size=3)
  
  return( p )
}

###########################
###########################
######## COLORING #########
######## FUNCTIONS ########
###########################
###########################

##assign color to gene names based on fdr significance value in MA Data
assignColorToGene = function(handle, maData){
  geneColor <- NULL
  geneColor$gene = maData$name
  sig=setMAColor(handle$fdrVal, maData) #, handle$indToTrack)
  geneColor$color = handle$maColor[sig]
  return(geneColor) 
}

#setMAColor gives codes 1, 2, 3, 4, 5 for the cases (up, down, non-significant, not available, tracked) of MA Data
setMAColor = function(fdrVal, maData){ #}, indToTrack) {
  sig <- rep(3, nrow(maData))
  sig[which(maData$padj <= fdrVal & maData$log2FoldChange < 0 & abs(maData$log2FoldChange) >= 
              maData$log2FoldChange & maData$detection_call == 1)] = 2
  sig[which(maData$padj <= fdrVal & maData$log2FoldChange > 0 & abs(maData$log2FoldChange) >= 
              maData$log2FoldChange & maData$detection_call == 1)] = 1
  sig[which(is.na(maData$padj))] = 4
  #sig[indToTrack] = 5
  return(sig)}