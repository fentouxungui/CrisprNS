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
#'
#' @examples
CrisprNS <- setClass(Class = "CrisprNS",
                     slots = list(project = "character",
                                  sgRNAs.meta.data = "data.frame", # sgRNAs id,gene,sgRNA.Index,and other info from library.metadata
                                  Samples.meta.data = "data.frame", # Do not change this! colnames: Samples,Group,Replicate,Stage,ControlsgRNAsMedian
                                  counts = "matrix",
                                  normalized = "matrix",
                                  logfc = "matrix",
                                  logMean = "matrix",
                                  library.meta.data = "data.frame",
                                  QuantileCutLog = "list",
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
#' sgRNAs.meta.data <- read.delim("../../../../HCT116cellline/results/mageck/Metabolism.count.txt",stringsAsFactors = FALSE)
#' library.meta.data <- read.csv("../../../../library/Library2-APC_SL_metabolism+essential+nontargeting.csv",stringsAsFactors = FALSE)
#' sample.meta.data <- data.frame(Samples = grep("Day",colnames(sgRNAs.meta.data),value = TRUE),
#'                               Group = factor(c(rep("WT",6),rep("Metabolism",6)),levels = c("WT","Metabolism")),
#'                               Replicates = factor(rep(paste("rep",1:3,sep = ""),4),levels = c(paste("rep",1:3,sep = ""))),
#'                               Stage = factor(rep(rep(paste("Day",c(0,14),sep = ""),each = 3),2),levels = c(paste("Day",c(0,14),sep = ""))))
#' test <- CreateCrisprNSObject(project = "CrisprNegativeSelection",
#'                             mageck.outputs = sgRNAs.meta.data,
#'                             library.meta.data = library.meta.data,
#'                             sample.meta.data = sample.meta.data)
CreateCrisprNSObject <- function(project = "CrisprNegativeSelection",
                                 mageck.outputs,
                                 library.meta.data,
                                 sample.meta.data # colnames = c("Samples","Group","Replicates","Stage") the last three should be factor!
                                 )
  {
  # load library
  suppressMessages(library(dplyr))

  # check full parameters!
  if (is.null(x = library.meta.data) | is.null(y = mageck.outputs) | is.null(z = sample.meta.data)) {
    stop("Please input full parameters!")
  }
  # check mageck outputs
  if (!identical(colnames(mageck.outputs)[1:2],c("sgRNA","Gene"))) {
    stop("please input mageck-output with first two column names should be sgRNA and Gene!")
  }
  # check sample meta data
  BasicSampleinfo <- c("Samples","Group","Replicates","Stage")
  if (sum(BasicSampleinfo %in% colnames(sample.meta.data)) != length(BasicSampleinfo)) {
    stop(paste("Please check the sample meta data info, which should include:", paste(BasicSampleinfo,collapse = ","),sep = " "))
  }
  if (sum(sample.meta.data$Samples %in% colnames(mageck.outputs)) != length(sample.meta.data$Samples) ) {
    stop("Some samples in sample meta data did not exist in the mageck outputs!")
  }
  # check library meta data
  keys <- c("ID","Type","Gene")
  if (sum(keys %in% colnames(library.meta.data)) != length(keys)) {
     stop(paste("The library meta data should include:",paste(keys,collapse = ","),sep = " "))
  }
  # check library meta data - Gene types
  message("Checking Gene types in libraray meta data!")
  GeneTypes <- unique(library.meta.data$Type)
  message(paste(length(GeneTypes),"types of Genes are found in library:"))
  CountStats <- table(library.meta.data$Type)
  message(paste(names(CountStats),collapse = ", "))
  message(paste(as.vector(CountStats),collapse = ", "))
  message(paste("Total",length(library.meta.data$ID),"sgRNAs are found, which corresponding to",length(unique(library.meta.data$Gene)),"genes!",sep = " "))
  message(paste(rep("-",50),sep = ""))

  # check if all the sgRNAs detected are found in library meta data
  message("Checking sgRNAs found in Mageck-outputs!")
  sgRNAs.meta.data <- mageck.outputs[,c("sgRNA","Gene")]
  if (sum(sgRNAs.meta.data$sgRNA %in% library.meta.data$ID) != length(sgRNAs.meta.data$sgRNA)) {
    sgRNAs.lost <- sgRNAs.meta.data$sgRNA[!sgRNAs.meta.data$sgRNA %in% library.meta.data$ID]
    message(paste("Attention! some sgRNAs are not found in library:",paste(sgRNAs.lost,collapse = ","),sep = " "))
    message("This may be a bug caused by R or Linux ~")
    message(paste("From the Mageck-output, Total",length(sgRNAs.meta.data$sgRNA),"are found, which corresponding to",length(unique(sgRNAs.meta.data$Gene)),"genes!",sep = " "))
  } else{
    message(paste("From the Mageck-output, Total",length(sgRNAs.meta.data$sgRNA),"are found, which corresponding to",length(unique(sgRNAs.meta.data$Gene)),"genes!",sep = " "))
  }
  message(paste(rep("-",50),sep = ""))

  # filter library meta data and add types and other infos to sgRNAs meta data
  ## remove shared 'Gene' column from library meta data
  library.meta.data.part <- library.meta.data[,!(colnames(library.meta.data) == "Gene")]
  ## add infos to sgRNAs meta data
  sgRNAs.ID <- sgRNAs.meta.data$sgRNA
  sgRNAs.meta.data <- merge(sgRNAs.meta.data,library.meta.data.part,all.x = TRUE,by.x = "sgRNA", by.y = "ID",sort = FALSE)
  rownames(sgRNAs.meta.data) <- sgRNAs.meta.data$sgRNA
  sgRNAs.meta.data <- sgRNAs.meta.data[sgRNAs.ID,]

  ## in case of NA in sgRNAs.meta.data - this may caused by a bug in R or Linux!
  if (NA %in% sgRNAs.meta.data$Type) {
    NAindex <- which(is.na(sgRNAs.meta.data$Type))
    NAGene <- sgRNAs.meta.data$Gene[NAindex]
    columns <- which(colnames(sgRNAs.meta.data) != "sgRNA")
    for (i in 1:length(NAindex)) {
      if (!NAGene %in% sgRNAs.meta.data$Gene) {
        stop("Error! Sorry! the Gene infos of the lost sgRNAs can not be filled with sgRNAs from same Gene!")
      }
      sgRNAs.meta.data[NAindex[i],columns] <- sgRNAs.meta.data[sgRNAs.meta.data$Gene == NAGene,columns][1,]
      message(paste("Attention! infos of sgRNA",sgRNAs.meta.data$sgRNA[NAindex[i]],"is filled using other sgRNAs from same gene for the NA value!",sep = " "))
    }
    message(paste(rep("-",50),sep = ""))
  }


  # preprocess the Sample meta data
  rownames(sample.meta.data) <- sample.meta.data$Samples
  sample.meta.data <- sample.meta.data[,(!colnames(sample.meta.data) == "Samples")]
  # set first level as control or initial stage!
  ## Group set
  if (length(unique(sample.meta.data$Group)) == 1) {
    message(paste("Only on Group is found:",levels(sample.meta.data$Group),"!",sep = " "))
  } else {
    message(paste(length(unique(sample.meta.data$Group)),"Groups are found:",paste(levels(sample.meta.data$Group),collapse = ", "),"!",sep = " "))
    message(paste("Default using",levels(sample.meta.data$Group)[1],"as the control sample!",sep =" "))
    message("You can use relvel(objcet@sample.meta.data$Group) to reset Group, the first level will be set as the control!")
  }
  message(paste(rep("-",50),sep = ""))

  ## Replicate set
  if (length(unique(sample.meta.data$Replicates)) == 1) {
    message(paste("Only on Replicate is found:",levels(sample.meta.data$Replicates),"!",sep = " "))
  } else {
    message(paste(length(unique(sample.meta.data$Replicates)),"Replicates are found:",paste(levels(sample.meta.data$Replicates),collapse = ", "),"!",sep = " "))
  }
  message(paste(rep("-",50),sep = ""))

  ## Stage set
  if (length(unique(sample.meta.data$Stage)) != 2) {
    stop("Stage set is wrong, Please define two stages: initital stage and final stage!")
  } else {
    message(paste(length(unique(sample.meta.data$Stage)),"Stages are found:",paste(levels(sample.meta.data$Stage),collapse = ", "),"!",sep = " "))
    message(paste(levels(sample.meta.data$Stage)[1],"is set as the initial stage, and ",levels(sample.meta.data$Stage)[2],"is set as the final stage!",sep = " "))
    message("You can use relvel(objcet@sample.meta.data$Stage) to reset stage, the first level will be set as the initial stage!")
  }
  message(paste(rep("-",50),sep = ""))

  # create object
  object <- new(
    Class = "CrisprNS",
    sgRNAs.meta.data = sgRNAs.meta.data,
    Samples.meta.data = sample.meta.data,
    counts = as.matrix(mageck.outputs[,rownames(sample.meta.data)]),
    library.meta.data = library.meta.data)

  return(object)
}
