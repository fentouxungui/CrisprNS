#' select sgRNAs
#' select sgRNAs from two logFC samples using control sgRNAs to set cutoff for nessary/unnessary sgRNAs. the specific nessary sgRNAs for ySample will be selected!
#'
#' @param object
#' @param Replicate samples from which replicate(if you have more than two groups, please use group pararmeter to define which two group)
#' @param Groups define using which two groups, default is the first two groups,and the first group is used for xSample, second for ySample.
#' @param xSample
#' @param ySample
#' @param controlsgRNA
#' @param quantile.cut a length 2 vector,to define the range of unnessary sgRNAs.
#' @param limitation
#' @param breaks
#'
#' @return a sgRNA list - a data frame
#' @export
#'
#' @examples
#' QuantileCutPair(test,Replicate = "rep1")
#' QuantileCutPair(test,Replicate = "rep1",Groups = c("Metabolism","WT"))
QuantileCutPair <- function(object,
                        Replicate = NULL,
                        Groups = NULL,
                        xSample = "xsample",
                        ySample = "ysample",
                        controlsgRNA = "nontargeting",
                        quantile.cut = c(0.95, 0.95),
                        limitation = c(-1,1),
                        breaks = c(-2,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,1.5,2)){
  Groups.levels <- levels(object@Samples.meta.data$Group)
  Replicate.levels <- levels(object@Samples.meta.data$Replicates)
  GeneTypes <- unique(object@sgRNAs.meta.data$Type)

  if (!controlsgRNA %in% GeneTypes) {
    stop("Error, Please define correct controlsgRNA-Gene types!")
  }
  if (is.null(Replicate)) {
    if (!xSample %in% colnames(object@logfc) | !ySample %in% colnames(object@logfc)) {
      stop("Error! Please input correct samples names defined in fc slot!")
    }
    tmp.data <- object@logfc[,c(xSample,ySample)]
  } else {
    if (!Replicate %in% Replicate.levels) {
      stop("Error! Please input correct replicate names defined in Replicate slot!")
    }
    if (is.null(Groups)) {
      xSample <- getXYsample(object = object, slot="logFC", Group = Groups.levels[1], Replicate = Replicate)
      ySample <- getXYsample(object = object, slot="logFC", Group = Groups.levels[2], Replicate = Replicate)
      # ySample <- paste(compareName,"logFC",sep = ".")
    } else if(length(Groups) != 2 | sum(Groups %in% Groups.levels) != 2){
      stop("Please input correct Groups!")
    } else {
      xSample <- getXYsample(object = object, slot="logFC", Group = Groups[1], Replicate = Replicate)
      ySample <- getXYsample(object = object, slot="logFC", Group = Groups[2], Replicate = Replicate)
    }
    tmp.data <-object@logfc[,c(xSample,ySample)]
  }
  controlsgRNAIndex <- object@sgRNAs.meta.data$Type == controlsgRNA
  controlsgRNAdata <- tmp.data[controlsgRNAIndex,]
  distance.cut.table <- sqrt(apply(apply(controlsgRNAdata,2,function(x)(x*x)),2,function(x)quantile(x, quantile.cut)))
  distance.cut <- c(distance.cut.table[1,1],distance.cut.table[2,2])

  message("Distance to origin Cutoff:")
  message(paste(xSample,distance.cut[1],sep = ": "))
  message(paste(ySample,distance.cut[2],sep = ": "))
  sgRNAIndex <- which(abs(tmp.data[,1]) < distance.cut[1] & tmp.data[,2] <  distance.cut[2]*(-1))
  tmp.data <- cbind(object@sgRNAs.meta.data[,1:2],tmp.data)
  # check sgRNA position!
  p <- ggplot2::ggplot(tmp.data, aes(x=get(xSample), y=get(ySample))) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    geom_hline(yintercept = -distance.cut[2],colour = "black",linetype = "dotted") +
    geom_vline(xintercept = c(-distance.cut[1],distance.cut[1]),colour = "black",linetype = "dotted") +
    geom_point(size = 0.3,color = "#bed2ed") +
    #lims(x=c(-6,6),y=c(-6,6)) +
    geom_point(data = tmp.data[sgRNAIndex,],size = 0.3,color = "red") +
    theme_minimal() +
    coord_fixed() +
    scale_x_continuous(breaks = breaks, limits = limitation) +
    scale_y_continuous(breaks = breaks, limits = limitation) +
    labs(x = xSample, y = ySample, title= paste("Highlight selected sgRNAs",sep = " "))
  print(p)
  pdf(file = paste("Highlight selected sgRNAs for",xSample,"and",ySample,".pdf",sep = " "),width = 8,height = 8)
  print(p)
  dev.off()
  selected.data <- tmp.data[sgRNAIndex,]
  annotation.lib <- object@sgRNAs.meta.data[,!(colnames(object@sgRNAs.meta.data) == "Gene")]
  selected.data <- merge(selected.data, annotation.lib, all.x = TRUE, by = "sgRNA")
  return(selected.data)
}



#' select sgRNAs in batch mode
#' see function - QuantileCutPair
#'
#' @param object
#' @param Replicates
#' @param Groups
#' @param controlsgRNA
#' @param quantile.cut
#' @param limitation
#' @param breaks
#'
#' @return selected sgRNA list | a data frame
#' @export
#'
#' @examples
#' selected <- BatchQuantileCutPair(test)
BatchQuantileCutPair <- function(object,
                             Replicates = NULL,
                             Groups = NULL,
                             controlsgRNA = "nontargeting",
                             quantile.cut = c(0.95, 0.95),
                             limitation = c(-1,1),
                             breaks = c(-2,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,1.5,2)){
  Groups.levels <- levels(object@Samples.meta.data$Group)
  Replicate.levels <- levels(object@Samples.meta.data$Replicates)
  GeneTypes <- unique(object@sgRNAs.meta.data$Type)

  if (is.null(Replicates)) { Replicates <- Replicate.levels}
  for (replicates in Replicates) {
    message(paste("> Handling",replicates,":",sep = " "))
    eachRep <- QuantileCutPair(object,
                           Replicate = replicates,
                           Groups = Groups,
                           controlsgRNA = controlsgRNA,
                           quantile.cut = quantile.cut,
                           limitation = limitation,
                           breaks = breaks)
    eachRep <- eachRep[,1:5]
    if (is.null(Groups)) {groups.used <- Groups.levels[1:2]}else{groups.used <- Groups}
    colnames(eachRep) <- c("sgRNA","Gene",groups.used,"Type")
    eachRep$Replicate <- replicates
    if (exists("merged")) {
      merged <- rbind(merged,eachRep)
    }else{
      merged <- eachRep
    }
  }
  return(merged)
}
