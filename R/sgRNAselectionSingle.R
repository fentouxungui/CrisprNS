#' cut logFC-logMean plot to select candidate sgRNAs/Genes
#' using nontargeting-sgRNAs as control for logFC cutoff, also a minimum log mean expression value cutoff!
#' @param object 
#' @param Replicate 
#' @param Group 
#' @param controlsgRNA 
#' @param Quantile.control 
#' @param min.logMean 
#' @param xSample 
#' @param x.limitation 
#' @param y.limitation 
#' @param x.breaks 
#' @param y.breaks 
#'
#' @return sgRNAs list - a data frame
#' @export
#'
#' @examples
#' CutSingle(test,Replicate = "rep1",Group = "WT",Quantile.control = 1)
CutSingle <- function(object,
                      Replicate = NULL,
                      Group = NULL,
                      controlsgRNA = "nontargeting",
                      Quantile.control = 0.95,
                      min.logMean = 2.5,
                      xSample = NULL,
                      x.limitation = c(-3,1),
                      y.limitation = c(0,9),
                      x.breaks = c(-3,-2,-1,0,1,2,3),
                      y.breaks = c(0,2,3,4,5,6,7,8)){
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
  
  controlsgRNAIndex <- object@sgRNAs.meta.data$Type == controlsgRNA
  controlsgRNAdata <- tmp.data[controlsgRNAIndex,]
  logFC.cut <- -(quantile(controlsgRNAdata[,1],Quantile.control))
  
  message("Cutoffs:")
  message(paste(xsample,logFC.cut,sep = ": "))
  message(paste(ysample,min.logMean,sep = ": "))
  
  sgRNAIndex <- which(tmp.data[,1] < logFC.cut & tmp.data[,2] >  min.logMean)
  tmp.data <- cbind(object@sgRNAs.meta.data[,1:2],tmp.data)
  # check sgRNA position!
  p <- ggplot2::ggplot(tmp.data, aes(x=get(xsample), y=get(ysample))) +
    geom_vline(xintercept = 0,colour = "gray") +
    geom_vline(xintercept = logFC.cut,colour = "black",linetype = "dotted") +
    geom_hline(yintercept = min.logMean,colour = "black",linetype = "dotted") +
    geom_point(cex = 1.5,colour = "black",alpha=0.2) +
    theme_bw() +
    #lims(x=c(-6,6),y=c(-6,6)) +
    geom_point(data = tmp.data[sgRNAIndex,],color = "red",alpha = 0.4) +
    scale_x_continuous(breaks = x.breaks,limits = x.limitation) +
    scale_y_continuous(breaks = y.breaks,limits = y.limitation) +
    labs(x = xsample, y = ysample, title= paste("Highlight selected sgRNAs",sep = " "))
  
  print(p)
  pdf(file = paste("Highlight selected sgRNAs for",xsample,"and",ysample,".pdf",sep = " "),width = 8,height = 8)
  print(p)
  dev.off()
  selected.data <- tmp.data[sgRNAIndex,]
  annotation.lib <- object@sgRNAs.meta.data[,!(colnames(object@sgRNAs.meta.data) == "Gene")]
  selected.data <- merge(selected.data, annotation.lib, all.x = TRUE, by = "sgRNA")
  return(selected.data)
}


#' cut logFC-logMean plot to select candidate sgRNAs/Genes in batch mode
#' see function - CutSingle
#'
#' @param object 
#' @param Replicates 
#' @param Groups 
#' @param controlsgRNA 
#' @param Quantile.control 
#' @param min.logMean 
#' @param x.limitation 
#' @param y.limitation 
#' @param x.breaks 
#' @param y.breaks 
#'
#' @return sgRNAs list - a data frame
#' @export
#'
#' @examples
#'  BatchCutSingle(test)
BatchCutSingle <- function(object,
                           Replicates = NULL,
                           Groups = NULL,
                           controlsgRNA = "nontargeting",
                           Quantile.control = 0.95,
                           min.logMean = 2.5,
                           x.limitation = c(-3,1),
                           y.limitation = c(0,9),
                           x.breaks = c(-3,-2,-1,0,1,2,3),
                           y.breaks = c(0,2,3,4,5,6,7,8)){
  Groups.levels <- levels(object@Samples.meta.data$Group)
  Replicate.levels <- levels(object@Samples.meta.data$Replicates)
  GeneTypes <- unique(object@sgRNAs.meta.data$Type)
  # check the parameter
  if (is.null(Groups)) {groups.used <- Groups.levels[1:2]}
  else{
    groups.used <- Groups
    if (length(groups.used) != sum(groups.used %in% Groups.levels)) {
    stop("Error! Please check your Groups parameter!")
  }}
  if (is.null(Replicates)) { Replicates <- Replicate.levels} else {
    if (length(Replicates) != sum(Replicates %in% Replicate.levels)) {
      stop("Error! Please check your Replicates parameter!")
    }
  }
  # calculating for each sample
  for (Agroup in Groups.levels) {
      for (replicates in Replicates) {
        message(paste("> Handling",Agroup, "-",replicates,":",sep = " "))
        eachRep <- CutSingle(object,
                             Replicate = replicates,
                             Group = Agroup,
                             controlsgRNA = controlsgRNA,
                             Quantile.control = Quantile.control,
                             min.logMean = min.logMean,
                             x.limitation = x.limitation,
                             y.limitation = y.limitation,
                             x.breaks = x.breaks,
                             y.breaks = y.breaks)
        eachRep <- eachRep[,1:5]
        
        colnames(eachRep) <- c("sgRNA","Gene","logFC","logMean","Type")
        eachRep$Replicate <- replicates
        eachRep$Group <- Agroup
        if (exists("merged")) {
          merged <- rbind(merged,eachRep)
        }else{
          merged <- eachRep
        }
    }
  }

  return(merged)
}
