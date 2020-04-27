#' Normalize raw data counts
#'
#' @param object
#' @param initialStage
#' @param controlsgRNA
#'
#' @return CrisprNS object with normalized, logFC, logMean slots
#' @export
#'
#' @examples
#' test <- NormalizeData(test,initialStage = "Day0",controlsgRNA = "nontargeting")
NormalizeData <- function(object,
                          initialStage = NULL,
                          controlsgRNA = "nontargeting"){
  # check if normalization has been taken!
  if (dim(object@normalized)[1] == 0) {
    message("The normalization has been done before, the old value will be removed!")
  }
  groups <- levels(object@Samples.meta.data$Group)
  replicateS <- levels(object@Samples.meta.data$Replicates)
  Genetypes <- unique(object@sgRNAs.meta.data$Type)
  # Check the parameters!
  if (class(object) != "CrisprNS") {
    stop("Error, Please input a CrisprNS Object created by CreateCrisprNSObject function")
  }
  if (is.null(initialStage)) {
    initialStage <- levels(object@Samples.meta.data$Stage)[1]
  } else {
    if(!initialStage %in% levels(object@Samples.meta.data$Stage)){
      stop("Error! the initial stage not exist in predefined Stages!")
    }
  }
  if (!controlsgRNA %in% Genetypes) {
    stop(paste("Error, Please define correct gene types from:", paste(Genetypes,collapse = ", "),sep = " "))
  }
  # change the initial stage if the initialStage is not the default
  if (initialStage != levels(object@Samples.meta.data$Stage)[1]) {
    object@Samples.meta.data$Stage <- relevel(object@Samples.meta.data$Stage,initialStage)
    message(paste("waring: The initial stage is reset to",initialStage,sep = " "))
  }

  # Normalize to same library size
  message(">> 1. Normalizing to same library size for each sample!")
  object@normalized <- apply(object@counts,2,function(x)((x/sum(x))*1e6 + 1))
  message(paste(rep("-",50),sep = ""))

  # Calculate Foldchange and log mean expression value for each pair (final stage vs initial stage)
  message(">> 2. Calculating log(Foldchange) and log mean expression value for each pair: final stage vs initial stage.")
  finalStage <- setdiff(levels(object@Samples.meta.data$Stage),initialStage)
  ## define a data.frame to store logFC and mean
  logfc.res <- data.frame(row.names = 1:nrow(object@counts))
  mean.res <- data.frame(row.names = 1:nrow(object@counts))

  for (group in groups) {
    for (replicates in replicateS) {
      initialsample <- extractSample(object@Samples.meta.data, Group = group, Replicate = replicates, Stage = initialStage)
      finalsample <- extractSample(object@Samples.meta.data, Group = group, Replicate = replicates, Stage = finalStage)
      # if no samples are found under this parameter, then break
      if (length(initialsample) == 0 | length(finalsample) == 0) {
        message(paste("Attention! No samples are found using this parameter: group -",group,"replicate -",replicates, sep = " "))
        break()
      }
      # if found different numbers of samples, then stop!
      if (length(initialsample) != 1 | length(finalsample) != 1) {
        stop("Error! found duplicate samples with same annotation, Please check the sample meta data!")
      }

      message(paste("@Using",initialsample,"as initial stage sample! And",finalsample,"as final stage sample!"sep = " "))
      # save results for each compare
      logfc <- log(object@normalized[,finalsample]/object@normalized[,initialsample])
      logfc.res$tmp <- logfc
      # remove common part in the name
      compareName <- MergeName(finalsample, initialsample, split = "\\.",combineBy = ".vs.")
      colnames(logfc.res)[which(colnames(logfc.res) == "tmp")] <- paste(compareName,"logFC",sep=".")

      mean.expression <- vector_mean(object@normalized[,finalsample],object@normalized[,initialsample])
      mean.res$tmp <- log(mean.expression)
      colnames(mean.res)[which(colnames(mean.res) == "tmp")] <- paste(compareName,"logMean",sep=".")
    }
  }
  object@logMean <- as.matrix(mean.res)
  message(paste(rep("-",50),sep = ""))


  # Normalizing logFC to the median control sgRNAs
  message(">> 3. Normalizing logFC to the median of the control sgRNAs!")
  ControlsgRNAindex <- which(object@sgRNAs.meta.data$Type == controlsgRNA)
  ControlMedian <- as.numeric(apply(logfc.res[ControlsgRNAindex,],2,median))
  message("Details of the control sgRNA's median value for each compare:")
  for (i in 1:length(ControlMedian)) {
    message(paste(colnames(logfc.res)[i],":",ControlMedian[i],sep = " "))
  }
  logfc.res <- sweep(logfc.res,2,ControlMedian,"-")
  object@logfc <- as.matrix(logfc.res)

  return(object)
}

#
