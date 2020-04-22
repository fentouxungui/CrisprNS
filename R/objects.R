#' Crispr Negative Selection Object
#'
#' @slot project character.
#' @slot raw.data data.frame.
#' @slot data data.frame.
#' @slot fc data.frame.
#' @slot compare character.
#' @slot meta.data data.frame.
#' @slot replicate character.
#' @slot version ANY.
#'
#' @return
#' @export
#'
#' @examples
CrisprNS <- setClass(Class = "CrisprNS",
                     slots = list(project = "character",
                                  raw.data = "data.frame",
                                  data = "data.frame",
                                  fc = "data.frame",
                                  mean = "data.frame",
                                  meta.data = "data.frame",
                                  Group = "character", #c("Control","Case"),
                                  Replicate = "character", #c("rep1","rep2","rep3"),
                                  Stage = "character", #c("initial","final"))
                                  library.metadata = "data.frame",
                                  Gene.Types = "character",
                                  sgRNA.Index = "data.frame",
                                  ControlMedian = "numeric",
                                  QuantileCut = "list",
                                  version = "ANY"),
                     #package = "",
                     S3methods = FALSE)


#' Create a Crispr Negative Selection Object
#'
#' @param counts
#' @param project
#' @param meta.data
#' @param Group
#' @param Replicate
#' @param Stage
#'
#' @return
#' @export
#'
#' @examples
CreateCrisprNSObject <- function(counts,
                                 project = "CrisprNegativeSelection",
                                 library.metadata,
                                 meta.data, # colnames = c("Samples","Group","Replicates","Stage")
                                 Group, # c("WT","Case"),
                                 Replicate, # c("rep1","rep2","rep3"),
                                 Stage # c("day0","day14")
                                 )
  {
  # check full parameters!
  if (is.null(x = meta.data) | is.null(y = Group) | is.null(z = Replicate) | is.null(w = Stage) | is.null(u = library.metadata)) {
    stop("Please input full parameters!")
  }
  # read in counts data - mageck output
  raw.data <- counts
  if (!identical(colnames(raw.data)[1:2],c("sgRNA","Gene"))) {
    stop("please input mageck-output with first colume names: sgRNA and Gene!")
  }
  # check meta.data
  if (!identical(colnames(meta.data),c("Samples","Group","Replicates","Stage"))) {
    stop("Please check the meta data file, which should include (Samples,Group,Replicates,Stage)")
  }
  if (sum(meta.data$Samples %in% colnames(raw.data)) != length(meta.data$Samples) ) {
    stop("Some samples in metadata did not exist in counts data!")
  }
  # check library meta data
  if (!identical(colnames(library.metadata)[1:3],c("ID","Type","Gene"))) {
     stop("The library meta data should start with ID,TYPE,Gene!")
   }
  Gene.Types <- unique(library.metadata$Type)
  message(paste(length(Gene.Types),"Types of Genes are found in library:"))
  CountStats <- table(library.metadata$Type)
  message(paste(names(CountStats),collapse = ", "))
  message(paste(as.vector(CountStats),collapse = ", "))
  sgRNA.Index <- data.frame(row.names = raw.data$sgRNA)
  for (type in Gene.Types) {
    choosen <- dplyr::filter(library.metadata,Type == type)
    sgRNA.Index$tmp <- raw.data$Gene %in% choosen$Gene
    colnames(sgRNA.Index)[which(colnames(sgRNA.Index) == "tmp")] <- type
  }
  # Add some function to tell how many sgRNAs/Genes are found in the raw.data



  # create object
  object <- new(
    Class = "CrisprNS",
    raw.data = raw.data,
    meta.data = meta.data,
    Group = Group,
    Replicate = Replicate,
    Stage = Stage,
    library.metadata = library.metadata,
    Gene.Types = Gene.Types,
    sgRNA.Index = sgRNA.Index
  )
  return(object)
}

#' Normalize raw data counts
#'
#' @param object
#' @param initialStage
#' @param controlsgRNA
#'
#' @return
#' @export
#'
#' @examples
NormalizeData <- function(object,
                          initialStage = "Day0",
                          controlsgRNA = "nontargeting"){
  if (class(object) != "CrisprNS") {
    stop("Error, Please input a CrisprNS Object created by CreateCrisprNSObject function")
  }
  Groups <- object@Group
  Replicates <- object@Replicate
  Genetypes <- object@Gene.Types
  if (!controlsgRNA %in% Genetypes) {
    stop("Error, Please define correct genetypes form Gene.Types slot!")
  }
  if (!initialStage %in% object@Stage) {
    stop("Error! the initial stage not exist in predefined Stages!")
  }
  # Normalize to same library size
  message("1. Normalizing to same library size for each sample!")
  samples <- object@meta.data$Samples
  rawdata <- as.data.frame(object@raw.data[,samples])
  tmpdata <- as.data.frame(apply(rawdata,2,function(x)((x/sum(x))*1e6 + 1)))
  object@data <- tmpdata
  # Calculate Foldchange for each pair (final stage vs initial stage)
  message("2. Calculate Foldchange for each pair (final stage vs initial stage)")
  finalStage <- setdiff(object@Stage,initialStage)
  res.tmp <- data.frame(row.names = 1:length(rownames(tmpdata)))
  for (groupinfo in Groups) {
    for (replicateinfo in Replicates) {
      initialsample <- extracname(samples,x=groupinfo,y = replicateinfo,z = initialStage)
      finalsample <- extracname(samples,x=groupinfo,y = replicateinfo,z = finalStage)
      logfc <- log(tmpdata[,finalsample]/tmpdata[,initialsample])
      res.tmp$tmp <- logfc
      compareName <- gsub(initialStage,paste(finalStage,"vs",initialStage,sep = "."),initialsample)
      colnames(res.tmp)[which(colnames(res.tmp) == "tmp")] <- compareName
    }
  }
  # Normalizing logFC to the median control sgRNAs
  message("3. Normalizing logFC to the median control sgRNAs!")
  ControlMedian <- as.numeric(apply(res.tmp[object@sgRNA.Index[,controlsgRNA],],2,median))
  object@ControlMedian <- ControlMedian
  for (i in 1:length(ControlMedian)) {
    message(paste(colnames(res.tmp)[i],":",ControlMedian[i],sep = " "))
  }
  res.tmp <- sweep(res.tmp,2,ControlMedian,"-")
  object@fc <- res.tmp

  return(object)
}


extracname <- function(vector,x,y,z="."){
  tmp <- vector
  for (i in c(x, y, z)) {
    tmp <- grep(i, tmp, value = TRUE)
  }
  return(tmp)
}

#' Draw scatter plot
#'
#' @param object
#' @param xSample
#' @param ySample
#' @param limitation
#' @param breaks
#'
#' @return
#' @export
#'
#' @examples
DrawFC <- function(object,
                   Replicate = NULL,
                   xSample = "xsample",
                   ySample = "ysample",
                   limitation = c(-4,4),
                   breaks = c(-4,-3,-2,-1,0,1,2,3,4)){
  if (is.null(Replicate)) {
    if (!xSample %in% colnames(object@fc) | !ySample %in% colnames(object@fc)) {
      stop("Error! Please input correct samples names defined in fc slot!")
    }
    tmp.data <- object@fc[,c(xSample,ySample)]
  } else {
    if (!Replicate %in% object@Replicate) {
      stop("Error! Please input correct replicate names defined in Replicate slot!")
    }
    xSample = extracname(colnames(object@fc),Replicate,object@Group[1])
    ySample = extracname(colnames(object@fc),Replicate,object@Group[2])
    tmp.data <-object@fc[,c(xSample,ySample)]
  }

  p0 <- ggplot2::ggplot(tmp.data, aes(x=get(xSample), y=get(ySample))) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    geom_point(size = 0.3,color = "#bed2ed") +
    #lims(x = limitation, y = limitation) +
    theme_minimal() +
    coord_fixed() +
    scale_x_continuous( breaks = breaks, limits = limitation) +
    scale_y_continuous( breaks = breaks, limits = limitation) +
    labs(x = xSample, y = ySample, title="Total sgRNAs")

  Genetypes <- object@Gene.Types
  GeneIndex <- object@sgRNA.Index

  for (i in 1:length(Genetypes)) {
    p <- ggplot2::ggplot(tmp.data, aes(x=get(xSample), y=get(ySample))) +
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
      geom_point(size = 0.3,color = "#bed2ed") +
      #lims(x=c(-6,6),y=c(-6,6)) +
      geom_point(tmp.data = tmp.data[GeneIndex[,Genetypes[i]],],size = 0.3,color = "red") +
      theme_minimal() +
      coord_fixed() +
      scale_x_continuous(breaks = breaks, limits = limitation) +
      scale_y_continuous(breaks = breaks, limits = limitation) +
      labs(x = xSample, y = ySample, title= paste("Highlight",Genetypes[i], "sgRNAs",sep = " "))
    assign(paste("p",i,sep = ""),p)
  }

  p.res <- c()

  for (i in 0:length(Genetypes)) {
    p.res <- append(p.res, paste("p",i,sep = ""))
  }
  ggpubr::ggarrange(plotlist =  mget(p.res),labels = LETTERS[1:(1+length(Genetypes))],nrow = ceiling((1+length(Genetypes))/2),ncol = 2)
}

QuantileCut <- function(object,
                        Replicate = NULL,
                        xSample = "xsample",
                        ySample = "ysample",
                        controlsgRNA = "nontargeting",
                        quantile.cut = c(0.95, 0.95),
                        limitation = c(-1,1),
                        breaks = c(-2,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,1.5,2)){
  if (!controlsgRNA %in% object@Gene.Types) {
    stop("Error, Please define correct genetypes form Gene.Types slot!")
  }
  if (is.null(Replicate)) {
    if (!xSample %in% colnames(object@fc) | !ySample %in% colnames(object@fc)) {
      stop("Error! Please input correct samples names defined in fc slot!")
    }
    tmp.data <- object@fc[,c(xSample,ySample)]
  } else {
    if (!Replicate %in% object@Replicate) {
      stop("Error! Please input correct replicate names defined in Replicate slot!")
    }
    xSample = extracname(colnames(object@fc),Replicate,object@Group[1])
    ySample = extracname(colnames(object@fc),Replicate,object@Group[2])
    tmp.data <-object@fc[,c(xSample,ySample)]
  }
  controlsgRNAdata <- tmp.data[object@sgRNA.Index[,controlsgRNA],]
  distance.cut.table <- sqrt(apply(apply(controlsgRNAdata,2,function(x)(x*x)),2,function(x)quantile(x, quantile.cut)))
  distance.cut <- c(distance.cut.table[1,1],distance.cut.table[2,2])
  # cut.info <- data.frame(row.names = colnames(tmp.data))
  # cut.info <- cbind(cut.info,quantile.cut,distance.cut)
  # object@QuantileCut <- append(object@QuantileCut,list(tmp = cut.info))
  # names(object@QuantileCut)[which(names(object@QuantileCut) == "tmp")] <- paste(xSample,"and",ySample,"using",controlsgRNA,"As-reference",sep = "-")
  message("Distance to origin Cutoff:")
  message(paste(xSample,distance.cut[1],sep = ": "))
  message(paste(ySample,distance.cut[2],sep = ": "))
  sgRNAIndex <- which(abs(tmp.data[,1]) < distance.cut[1] & tmp.data[,2] <  distance.cut[2]*(-1))
  tmp.data <- cbind(object@raw.data[,1:2],tmp.data)
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
  selected.data <- merge(selected.data,object@library.metadata, all.x = TRUE, by.x = "sgRNA",by.y = "ID")
  return(selected.data)
}


BatchQuantileCut <- function(object,
                             Replicates = "ALL",
                             controlsgRNA = "nontargeting",
                             quantile.cut = c(0.95, 0.95),
                             limitation = c(-1,1),
                             breaks = c(-2,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,1.5,2)){
  if (Replicates == "ALL") { Replicates <- object@Replicate}
  for (replicates in Replicates) {
    message(paste("> Handling",replicates,":",sep = " "))
    eachRep <- QuantileCut(object,
                Replicate = replicates,
                controlsgRNA = controlsgRNA,
                quantile.cut = quantile.cut,
                limitation = limitation,
                breaks = breaks)
    eachRep <- eachRep[,1:5]
    colnames(eachRep) <- c("sgRNA","Gene",object@Group,"Type")
    eachRep$Replicate <- replicates
    if (exists("merged")) {
      merged <- rbind(merged,eachRep)
    }else{
      merged <- eachRep
    }
  }
  return(merged)
}

#' Static Genes in all replicates
#'
#' @param AnDataFrame
#'
#' @return
#' @export
#'
#' @examples
staticsgRNAs <- function(AnDataFrame){
  tmp <- reshape2::dcast(AnDataFrame,Gene~Replicate)
  tmp$sum <- apply(tmp[2:length(colnames(tmp))],1, sum)
  tmp <- dplyr::arrange(tmp,desc(sum))
  return(tmp)
}
