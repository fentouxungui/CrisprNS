#' Draw a paired sample
#' Draw sgRNAs logFC-logFC scatter plot for each gene type
#' @param object
#' @param mode paired or single. paired: using two samples - logFC:logFC; single: one sample - logFC:MeanExpression
#' @param Group a vector. length 2 for paired mode, when more than two groups, or you want change x and y samples. length 1 for single mode.
#' default using first level as xSample, second level as ySample.
#' @param Replicate paired mode, samples from which replicate, used together with groups parameter when you have more than one replicate.
#' @param xSample the sample for x coordinate, if no ySample defined, using single mode!
#' @param ySample the sample for y coordinate
#' @param limitation x,y coordinate limit
#' @param breaks x,y coordinate breaks
#'
#' @return a plot
#' @export
#'
#' @examples
#' DrawPaired(test,Replicate = "rep2")
#' DrawPaired(test,xSample = compares[2],ySample = compares[4])
DrawPaired <- function(object,
                   Group = NULL,
                   Replicate = NULL,
                   xSample = NULL,
                   ySample = NULL,
                   limitation = c(-4,4),
                   breaks = c(-4,-3,-2,-1,0,1,2,3,4)){
  # if xSample and ySample both not defined
  if (is.null(xSample) & is.null(ySample)) {
    if (is.null(Group)) { # group == NULL
      if (length(levels(object@Samples.meta.data$Group)) == 1) { # Only found one group
        stop("Error! Only one group found, Please use DrawSingle!")
      } else { # using level1 as xsample.group, level2 as ysamle.group
        xSample.group <- levels(object@Samples.meta.data$Group)[1]
        ySample.group <- levels(object@Samples.meta.data$Group)[2]
      }
    }else if(length(Group) == 2 & sum(Group %in% levels(object@Samples.meta.data$Group)) == 2){ # Group are set,and right
      xSample.group <- Group[1]; ySample.group <- Group[2]
    }else { # error!
      stop("Error! Please check the Group parameter!")
    }
    # set replicate
    if(is.null(Replicate)){
      if (length(levels(object@Samples.meta.data$Replicates)) == 1) {
        Sample.replicate <- levels(object@Samples.meta.data$Replicates)
      }else{
        stop("Error! Please input Replicate parameter!")
        }
    }else if(length(Replicate) == 1 & Replicate %in% object@Samples.meta.data$Replicates){
      Sample.replicate <- Replicate
    }else {
      stop("Error! Please input correct Replicate parameter!")
    }
    xsample <- getXYsample(object = object, slot="logFC", Group = xSample.group, Replicate = Sample.replicate)
    ysample <-getXYsample(object = object, slot="logFC", Group = ySample.group, Replicate = Sample.replicate)
  } else if(!is.null(xSample) & !is.null(ySample)){ # if xSample and ySample both defined
    if (xSample %in% colnames(object@logfc) & ySample %in% colnames(object@logfc)) { # if xSample and ySamle both exist
      xsample <- xSample
      ysample <- ySample
    }else{ # if xSample and ySamle are not both exist
      stop("Error! Please check if the samples exist in colnames(object@logfc)!")
    }
  }else { # if xSample and ySample are not both defined
    stop("Error, Please define xSample and ySample at the same time!")
  }
 tmp.data <- as.data.frame(object@logfc[,c(xsample,ysample)])
 DrawCore(object = object,tmp.data = tmp.data,xsample = xsample, ysample = ysample,limitation = limitation,breaks = breaks)
}



#' Draw a single Sample (final stage vs initial stage)
#' Draw sgRNAs logFC-logMean scatter plot for each gene type
#' @param object
#' @param Group
#' @param Replicate
#' @param xSample
#' @param x.limitation
#' @param y.limitation
#' @param x.breaks
#' @param y.breaks
#'
#' @return a plot
#' @export
#'
#' @examples
#' DrawSingle(test,Group = "WT",Replicate = "rep1")
DrawSingle <- function(object,
                       Group = NULL,
                       Replicate = NULL,
                       xSample = NULL,
                       x.limitation = c(-3,1),
                       y.limitation = c(0,9),
                       x.breaks = c(-3,-2,-1,0,1,2,3),
                       y.breaks = c(0,2,3,4,5,6,7,8)
                       ){
  # if xSample and ySample both not defined
  if (is.null(xSample)) {
    if (is.null(Group)) { # group == NULL
      stop("Error! Please specify a Group!")
    }else if(length(Group) == 1 & Group %in% levels(object@Samples.meta.data$Group)){ # Group are set,and right
      xSample.group <- Group
    }else { # error!
      stop("Error! Please check the Group parameter!")
    }
    # set replicate
    if(is.null(Replicate)){
      if (length(levels(object@Samples.meta.data$Replicates)) == 1) {
        Sample.replicate <- levels(object@Samples.meta.data$Replicates)
      }else{
        stop("Error! Please input Replicate parameter!")
      }
    }else if(length(Replicate) == 1 & Replicate %in% object@Samples.meta.data$Replicates){
      Sample.replicate <- Replicate
    }else {
      stop("Error! Please input correct Replicate parameter!")
    }
    xsample <- getXYsample(object = object, slot="logFC", Group = xSample.group, Replicate = Sample.replicate)
    ysample <- getXYsample(object = object, slot="logMean", Group = xSample.group, Replicate = Sample.replicate)

  } else if(!is.null(xSample)){ # if xSample and ySample both defined
    if (xSample %in% colnames(object@logfc)) { # if xSample and ySamle both exist
      xsample <- xSample
      ysample <- gsub("logFC","logMean",xsample)
    }else{ # if xSample and ySamle are not both exist
      stop("Error! Please check if the samples exist in colnames(object@logfc)!")
    }
  }
  tmp.data <- as.data.frame(cbind(object@logfc[,xsample],object@logMean[,ysample]))
  colnames(tmp.data) <- c(xsample,ysample)
  DrawCoreSingle(object=object,tmp.data=tmp.data,xsample=xsample,ysample=ysample,x.limitation=x.limitation, y.limitation=y.limitation,x.breaks=x.breaks,y.breaks=y.breaks)
}


#' Core draw code used for DrawPaired function
#'
#' @param object
#' @param tmp.data
#' @param xsample
#' @param ysample
#' @param limitation
#' @param breaks
#'
#' @return a plot
#' @export
#'
#' @examples
DrawCore <- function(object,tmp.data,xsample,ysample,limitation,breaks){
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggpubr))
  p0 <- ggplot2::ggplot(tmp.data, aes(x=get(xsample), y=get(ysample))) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    geom_point(size = 0.3,color = "#bed2ed") +
    #lims(x = limitation, y = limitation) +
    theme_minimal() +
    coord_fixed() +
    scale_x_continuous( breaks = breaks, limits = limitation) +
    scale_y_continuous( breaks = breaks, limits = limitation) +
    labs(x = xsample, y = ysample, title="Total sgRNAs")

  Genetypes <- unique(object@sgRNAs.meta.data$Type)

  for (i in 1:length(Genetypes)) {
    p <- ggplot2::ggplot(tmp.data, aes(x=get(xsample), y=get(ysample))) +
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
      geom_point(size = 0.3,color = "#bed2ed") +
      #lims(x=c(-6,6),y=c(-6,6)) +
      geom_point(data = tmp.data[object@sgRNAs.meta.data$Type == Genetypes[i],],size = 0.3,color = "red") +
      theme_minimal() +
      coord_fixed() +
      scale_x_continuous(breaks = breaks, limits = limitation) +
      scale_y_continuous(breaks = breaks, limits = limitation) +
      labs(x = xsample, y = ysample, title= paste("Highlight",Genetypes[i], "sgRNAs",sep = " "))
    assign(paste("p",i,sep = ""),p)
  }

  p.res <- c()

  for (i in 0:length(Genetypes)) {
    p.res <- append(p.res, paste("p",i,sep = ""))
  }
  ggpubr::ggarrange(plotlist =  mget(p.res),labels = LETTERS[1:(1+length(Genetypes))],nrow = ceiling((1+length(Genetypes))/2),ncol = 2)
}

#' Core draw code used for DrawSingle function
#'
#' @param object
#' @param tmp.data
#' @param xsample
#' @param ysample
#' @param x.limitation
#' @param y.limitation
#' @param x.breaks
#' @param y.breaks
#'
#' @return a plot
#' @export
#'
#' @examples
DrawCoreSingle <- function(object,tmp.data,xsample,ysample,x.limitation, y.limitation,x.breaks,y.breaks){
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggpubr))
  # Total sgRNAs
  p0 <- ggplot2::ggplot(tmp.data, aes(x=get(xsample), y=get(ysample))) +
    geom_vline(xintercept = 0,colour = "gray") +
    geom_point(cex = 1.5,colour = "black",alpha=0.2) +
    theme_bw() +
    scale_x_continuous(breaks = x.breaks,limits = x.limitation) +
    scale_y_continuous(breaks = y.breaks,limits = y.limitation) +
    labs(x=xsample,y=ysample,title="Total sgRNAs")

  Genetypes <- unique(object@sgRNAs.meta.data$Type)

  for (i in 1:length(Genetypes)) {
    p <- ggplot2::ggplot(tmp.data, aes(x=get(xsample), y=get(ysample))) +
      geom_vline(xintercept = 0,colour = "gray") +
      geom_point(cex = 1.5,colour = "black",alpha=0.2) +
      theme_bw() +
      #lims(x=c(-6,6),y=c(-6,6)) +
      geom_point(data = tmp.data[object@sgRNAs.meta.data$Type == Genetypes[i],],color = "red",alpha = 0.4) +
      scale_x_continuous(breaks = x.breaks,limits = x.limitation) +
      scale_y_continuous(breaks = y.breaks,limits = y.limitation) +
      labs(x = xsample, y = ysample, title= paste("Highlight",Genetypes[i], "sgRNAs",sep = " "))
    assign(paste("p",i,sep = ""),p)
  }

  p.res <- c()

  for (i in 0:length(Genetypes)) {
    p.res <- append(p.res, paste("p",i,sep = ""))
  }
  ggpubr::ggarrange(plotlist =  mget(p.res),labels = LETTERS[1:(1+length(Genetypes))],nrow = ceiling((1+length(Genetypes))/2),ncol = 2)
}


