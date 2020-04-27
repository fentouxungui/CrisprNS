#' Static Genes/sgRNAs in all replicates /and groups combination
#' return a sgRNAs' or Genes' distribution in each replicate /and groups combination
#'
#' @param AnDataFrame result from BatchQuantileCutPair
#'
#' @return a data frame
#' @export
#'
#' @examples
#' static(b,formula.y = "Type",formula.x = "Replicate")
#' static(b,formula.y = "Type",formula.x = c("Replicate","Group"))
static <- function(AnDataFrame,formula.y= "Gene", formula.x = "Replicate"){
  suppressMessages(library(dplyr))
  suppressMessages(library(reshape2))
  dcast.formula <- paste(formula.y,"~",paste(formula.x,collapse = " + "),sep = " ")
  tmp <- reshape2::dcast(AnDataFrame,formula = formula(dcast.formula))
  tmp$sum <- apply(tmp[2:length(colnames(tmp))],1, sum)
  tmp <- dplyr::arrange(tmp,desc(sum))
  return(tmp)
}

