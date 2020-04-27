#' Calculate the mean value in each position for two vectors
#'
#' @param vector1 a numeric vector
#' @param vector2 a numeric vector
#'
#' @return a vector
#' @export
#'
#' @examples
#' vector_mean(c(1,2,3),c(4,5,6))
vector_mean <- function(vector1,vector2){
  tmp <- c()
  if (length(vector1) == length(vector2)) {
    for (i in 1:length(vector1)) {
      tmp <- append(tmp,mean(c(vector1[i],vector2[i])))
    }
  }
  return(tmp)
}



#' Merge two character with removing the common part in the second character
#'
#' @param characterA
#' @param characterB
#' @param split how characterA and characterB are separated by
#' @param combineBy join character
#'
#' @return a character without common part
#' @export
#'
#' @examples
#' >>characterA <- "HCT116.Metabolism.Day14.rep3"
#' >>characterB <- "HCT116.Metabolism.Day0.rep3"
#' >>MergeName(characterA = characterA, characterB = characterB)
#' "HCT116.Metabolism.Day14.rep3-vs-Day0"
MergeName <- function(characterA,characterB,split = "\\.",combineBy = "-vs-"){
  commonPart <- intersect(unlist(strsplit(characterA,split=split)),unlist(strsplit(characterB,split=split)))
  for (commonpart in commonPart) {
    characterB <- gsub(commonpart,"",characterB)
  }
  # remove split character in head
  characterB <- gsub(paste("^",split,"*",sep = ""),"",characterB)
  # remove split character in tail
  characterB <- gsub(paste(split,"*","$",sep = ""),"",characterB)
  tmp <- paste(characterA,characterB,sep = combineBy)
  return(tmp)
}

#' extract rownames only for samples.meta.data slot from CrisprNS object
#' the data frame should have colnames: Group, Replicates, Stage
#'
#' @param ADataFrame
#' @param Group
#' @param Replicates
#' @param Stage
#'
#' @return a character or a character vector
#' @export
#'
#' @examples
#' extractSample(object@Samples.meta.data, Group = "WT", Replicate = "rep1", Stage = "Day0")
extractSample <- function(ADataFrame,Group, Replicates, Stage){
  tmp <- ADataFrame[ADataFrame$Group == Group & ADataFrame$Replicates == Replicates & ADataFrame$Stage == Stage,]
  return(rownames(tmp))
}


# extracname <- function(vector,x,y,z="."){
#   tmp <- vector
#   for (i in c(x, y, z)) {
#     tmp <- grep(i, tmp, value = TRUE)
#   }
#   return(tmp)
# }

#' Get xSample/ySample Name in logFC slot or logMean slot
#' 
#' @param object 
#' @param slot 
#' @param Group 
#' @param Replicate 
#'
#' @return a character | a name in logFC slot or logMean slot
#' @export
#'
#' @examples
#' getXYsample(object,slot="logFC",Group = "WT", Replicate = "rep1")
getXYsample <- function(object,
                        slot = "logFC",
                        Group = "WT",
                        Replicate = "rep1"){
  initialStage <- levels(object@Samples.meta.data$Stage)[1]
  finalStage <- levels(object@Samples.meta.data$Stage)[2]
  initialsample <- extractSample(object@Samples.meta.data, Group = Group, Replicate = Replicate, Stage = initialStage)
  finalsample <- extractSample(object@Samples.meta.data, Group = Group, Replicate = Replicate, Stage = finalStage)
  compareName <- MergeName(finalsample, initialsample, split = "\\.",combineBy = ".vs.")
  XYsample <-  paste(compareName,slot,sep=".")
  return(XYsample)
}