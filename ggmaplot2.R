# ggmaplot2 is modified from https://github.com/kassambara/ggpubr/blob/master/R/ggmaplot.R
# ggmaplot from library("ggpubr")  
# it serves at plotting the MAplot in Enhanced MA plot V2
# library("ggpubr")

parse_font <- function(font){
  if(is.null(font)) res <- NULL
  else if(inherits(font, "list")) res <- font
  else{
    # matching size and face
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
    if(length(face) == 0) face <- NULL else face <- font[face]
    color <- setdiff(font, c(size, face))
    if(length(color) == 0) color <- NULL
    res <- list(size=size, face = face, color = color)
  }
  res
}

ggmaplot2<-function (data, indToTrack=NULL,indToHighlight=NULL, fdr = 0.05, fc = 1.5, genenames = NULL, detectionCall = NULL,
                     sizepts = 5, 
                     font.label = c(12, "plain", "black"), font.legend = "bold", font.main = "bold", 
                     label.rectangle = FALSE,
                     palette = c("#B31B21", "#1465AC", "darkgray", "yellow","green"), top = 15,
                     select.top.method = c("pAdj", "fc"), main = NULL, 
                     xlab = "Log2(1 + Mean expression)",
                     ylab = "Log2(Fold change)", 
                     ggtheme = theme_classic(), 
                     legend = "top", 
                     discrete=T, p,
                     ...)
{
  #'MA-plot from means and log fold changes
  #'@description Make MA-plot which is a scatter plot of log2 fold changes (M, on
  #'  the y-axis) versus the average expression signal (A, on the x-axis). \code{M
  #'  = log2(x/y)} and \code{A = (log2(x) + log2(y))/2 = log2(xy)*1/2}, where x
  #'  and y are respectively the mean of the two groups being compared.
  #'@inheritParams ggboxplot
  #'@inheritParams ggpar
  #'@param data an object of class DESeqResults, get_diff, DE_Results, matrix or
  #'  data frame containing the columns baseMean (or baseMeanLog2),
  #'  log2FoldChange, and padj. Rows are genes.
  #'
  #'  Format accepted for the input data: \
  #'  \code{baseMean | log2FoldChange | pAdj | (detectionCall)}. This is a typical output from
  #'  DESeq2 pipeline. we use log2(baseMean+!) as the x-axis variable.
  #'  
  #'  Terminology:
  #'
  #'  \itemize{ \item baseMean: the mean expression of genes in the two groups.
  #'  \item log2FoldChange: the log2 fold changes of group 2 compared to group 1
  #'  \item pAdj: the adjusted p-value of the used statiscal test. }
  #'@param fdr Accepted false discovery rate for considering genes as
  #'  differentially expressed.
  #'@param fc the fold change threshold. Only genes with a fold change >= fc and
  #'  pAdj <= fdr are considered as significantly differentially expressed.
  #'@param genenames a character vector of length nrow(data) specifying gene names
  #'  corresponding to each row. Used for point labels.
  #'@param detectionCall a numeric vector with length = nrow(data), specifying if
  #'  the genes is expressed (value = 1) or not (value = 0). For example
  #'  detectionCall = c(1, 1, 0, 1, 0, 1). Default is NULL. If detectionCall
  #'  column is available in data, it will be used.
  #'@param size points size.
  #'@param alpha numeric value betwenn 0 an 1 specifying point alpha for
  #'  controlling transparency. For example, use alpha = 0.5.
  #'@param font.label a vector of length 3 indicating respectively the size (e.g.:
  #'  14), the style (e.g.: "plain", "bold", "italic", "bold.italic") and the
  #'  color (e.g.: "red") of point labels. For example \emph{font.label = c(14,
  #'  "bold", "red")}.
  #'@param label.rectangle logical value. If TRUE, add rectangle underneath the
  #'  text, making it easier to read.
  #'@param top the number of top genes to be shown on the plot. Use top = 0 to
  #'  hide to gene labels.
  #'@param select.top.method methods to be used for selecting top genes. Allowed
  #'  values include "padj" and "fc" for selecting by adjusted p values or fold
  #'  changes, respectively.
  #'@param label.select character vector specifying some labels to show.
  #'@param ... other arguments to be passed to \code{\link{ggpar}}.
  #'@return returns a ggplot.
  # 
  
  data_orig <- data
  # ggmaplot2 is modified from https://github.com/kassambara/ggpubr/blob/master/R/ggmaplot.R
  # ggmaplot from library("ggpubr") except for the colorscale: 
  # we added a color for the NA values in the pAdj column of data.
  # palette is a list of colors for (up, down, non-significant, not available) values of pAdj respectively
  # by Default: Palette=c("#B31B21", "#1465AC", "darkgray", "yellow")
  #
  # PARAMETERS FOR TESTING
  # data=maData
  # main = expression("Group 1" %->% "Group 2")
  # xlab = "Log2 mean expression"
  # ylab = "Log2 fold change"
  # fdr = 0.05
  # fc = 2
  # palette = c("#B31B21", "#1465AC", "darkgray", "yellow")
  # genenames = as.vector(maData$name)
  # legend = "top"
  # top = 20
  # ggtheme = ggplot2::theme_minimal()
  # select.top.method = "pAdj" #c("pAdj", "fc")
  # label.rectangle=TRUE
  # 
  # sizepts=0.4
  # font.label = c("bold", 11)
  # font.legend = "bold"
  # font.main = "bold"
  
  if (!base::inherits(data, c("matrix", "data.frame", "DataFrame",
                              "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  
  detectionCall<-data$detectionCall
  
  if (!is.null(detectionCall)) {
    if (nrow(data) != length(detectionCall))
      stop("detectionCall must be a numeric vector of length = nrow(data)")
  }  else if ("detectionCall" %in% colnames(data)) {
    detectionCall <- as.vector(data$detectionCall)
  }  else detectionCall = rep(1, nrow(data))
  
  if (is.null(legend)) 
    legend <- c(0.12, 0.9)
  ss <- base::setdiff(c("baseMean", "log2FoldChange", "pAdj"), 
                      colnames(data))
  if (length(ss) > 0) 
    stop("The colnames of data must contain: ", paste(ss, 
                                                      collapse = ", "))
  if (is.null(genenames)) 
  {
    genenames <- rownames(data)
  }else if (length(genenames) != nrow(data)) 
    stop("genenames should be of length nrow(data).")
  sig = setMAColor(fdr, data) #,indToTrack)
    # sig <- rep(3, nrow(data))
  # sig[which(data$pAdj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= 
  #             log2(fc) & detectionCall == 1)] = 2
  # sig[which(data$pAdj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= 
  #             log2(fc) & detectionCall == 1)] = 1
  # sig[which(is.na(data$pAdj))] = 4
  
 
  T <- fdr
  data <- data.frame(name = genenames, 
                     #mean = data$baseMean,
                     log2mean=log2(data$baseMean + 1),
                     log2FC = data$log2FoldChange, 
                     pAdj = data$pAdj, 
                     sig = sig,
                     pFold = data$pFold,
                     pFoldDummy = cut(handle$MAdataCur$pFold, 
                                      #breaks=c(-1, -T-2*(1-T)/3, -T-(1-T)/3 ,-T, 0, T, T+(1-T)/3,T+(1-T)/3*2, 1), 
                                      breaks=c(-1, -T, -2*T/3 ,-T/3, 0, T/3, 2*T/3, T,  1), 
                                      labels = c("Down-gray", "Down-low", "Down-mid", "Down-high", 
                                                 "Up-high", "Up-mid", "Up-low", "Up-gray")),
                     tracked=rep(FALSE,nrow(data)), 
                     highlighted=rep(FALSE,nrow(data)), 
                     stringsAsFactors = FALSE)
  
  
  
  if (!is.null(indToTrack)) data$tracked[indToTrack]<-TRUE
  
  if (!is.null(indToHighlight) && length(indToHighlight)>0) data$highlighted[indToHighlight]<-TRUE
  
  
   
  . <- NULL
  data$sig <- as.factor(data$sig)
  DOTlev <- levels(data$sig) %>% as.numeric()
  palette <- palette[DOTlev]
  new.levels <- c(paste0("Up "), #, sum(sig == 1)), 
                  paste0("Down "), #,sum(sig == 2)), 
                  "Non-Significant","Not Available","Tracked") %>% .[DOTlev]
  data$sig <- factor(data$sig, labels = new.levels)
  
   
  # select.top.method <- match.arg(select.top.method)
  if (select.top.method == "pAdj") 
  {
    data <- data[order(data$pAdj), ]
    
  }else if (select.top.method == "fc") 
    data <- data[order(abs(data$log2FC), decreasing = TRUE), 
                 ]
  labs_data <- stats::na.omit(data)
  labs_data <- subset(labs_data, pAdj <= fdr & name != "" & 
                        abs(log2FC) >= log2(fc))
  labs_data <- utils::head(labs_data, top)
  
  
  font.label <- parse_font(font.label)
  font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color <- ifelse(is.null(font.label$color), "black", 
                             font.label$color)
  font.label$face <- ifelse(is.null(font.label$face), "plain", 
                            font.label$face)
  set.seed(42)
  mean <- log2FC <- sig <- name <- pAdj <- NULL
  # p <- ggplot(data, aes(x =log2mean, y = log2FC, label=  name)) + 
  #   geom_point(aes(color = sig), size = sizepts) +

  
  if(is.null(p))
    p <- ggplot()
  
  if(discrete) {
    p <- p + geom_point(data=data, 
              aes(x = log2mean, y = log2FC, label=name, fill = sig, color=sig), 
              size = 1, pch = 21, show.legend = F)
    
    p <- p + geom_point(data=data[data$tracked,], 
                        aes(x = log2mean, y = log2FC, label=name, fill = sig), 
                        size = 2, color = "green", pch = 21, show.legend = F)

    p <- p + geom_point(data=data[data$highlighted,], 
                 aes(x = log2mean, y = log2FC, label=  name,fill = sig), 
                 size = 4, color = "orange", pch = 21, show.legend = F)
    
  } else {

    
    my_rgb <- function(r,g,b) { rgb(r,g,b,maxColorValue = 255) }
    
    # pFoldDummy -1 to 1:
    cols <- c( "Up-gray"=my_rgb(166,166,166), "Up-low"=my_rgb(255,127,127), "Up-mid"=my_rgb(255,51,51), "Up-high"=my_rgb(192,0,0), 
               "Down-gray"=my_rgb(166,166,166), "Down-high"=my_rgb(32,78,121), "Down-mid"=my_rgb(46,117,182), "Down-low"=my_rgb(91,156,213))
    
    specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

    # highlight
    s <- data[data$highlighted,]
    
    if(dim(s)[1] > 0) {
      df_line <- data.frame(x=c(handle$ranges$X[1]-100, s[1,]$log2mean, s[1,]$log2mean),
                            y=c(s[1,]$log2FC, s[1,]$log2FC, handle$ranges$Y[1]-100))
      
      
      p <- p + geom_path(data=df_line, aes(x, y, group=1), linetype = "dashed", color="orange")
    }
    
    p <- p+ geom_point(data=data[data$highlighted,], 
                       aes(x = log2mean, y = log2FC), 
                       size = 3, color = "orange", fill=NA, stroke = 0.5, pch = 21, show.legend = F)
                       #size = 4, fill='orange', color = "orange", pch = 21, show.legend = F)
    
    p <- p + geom_point(data=data[data$tracked,], 
                        aes(x = log2mean, y = log2FC), 
                        size = 3, color = handle$maColor[5], fill=NA, stroke = 0.5, pch = 21, show.legend = F)
    
    p <- p + geom_point(data=data,
                        aes(x = log2mean, y = log2FC, label=name, 
                            fill=pFoldDummy, colour=pFoldDummy, 
                            text=paste0('</br>Gene name: ', name, '</br>',
                                       'P-value: ', format(pAdj, digits=5), '</br>',
                                       'Log2(1 + Mean Expression): ', round(log2mean, digits = 2), '</br>',
                                       'Log2(Fold Change): ', round(log2FC, digits = 2), '</br>')), 
                        size = 1.5, stroke=0,
                        pch = 21) + 
      
      labs(colour='P-value', fill='P-value') + 
      scale_color_manual(values = cols, guide=F) + 
      scale_fill_manual(values = cols, guide=F)
    
    p <- p + ggtheme #ggpar(p, ggtheme = ggtheme)
    p <- p + 
      theme(
        legend.position = c(.95, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
      )
    
  }
  
  
    
  
  # if (label.rectangle) {
  #   p <- p + ggrepel::geom_label_repel(data = labs_data,
  #                                      mapping = aes(label = name), box.padding = unit(0.35,
  #                                                                                      "lines"), point.padding = unit(0.3, "lines"),
  #                                      force = 1, fontface = font.label$face, size = font.label$size/3,
  #                                      color = font.label$color)
  # } else {
  #   p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name), 
  #                                     box.padding = unit(0.35, "lines"), point.padding = unit(0.3, 
  #                                                                                             "lines"), force = 1, fontface = font.label$face, 
  #                                     size = font.label$size/3, color = font.label$color)
  # }
  p <- p + scale_x_continuous(breaks = seq(0, max(data$log2mean), 2)) 
  
  p <- p + labs(x = xlab, y = ylab, title = main, color = "") 
  #p <- p + geom_hline(yintercept = c(0, -log2(fc), log2(fc)), # dashed lines across the plot
  #                    linetype = c(1, 2, 2), 
  #                    color = c("black", "black", "black"))
  
  if(discrete) {
    p <- ggpar(p, palette = palette)
    p <- ggpar(p, ggtheme = ggtheme)
  }
  
  p <- p + theme(legend.position = "none")
  #p <- p + theme(legend.position = "right")
  p
}