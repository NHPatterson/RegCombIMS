#' @export
#' @title Get correlation matrix from Cardinal dataset
#' @description Makes a correlation matrix using base R's 'cor' function, names columns and rows mz
#' @param combined_cardinaldata an object of class "MSImageSet", first dataset
#' @return a correlation matrix with row and columns names as mzs

getCorrelationMat <- function(combined_cardinaldata) {
  
  correlation_matrix = cor(t(iData(combined_cardinaldata)))
  
  row.names(correlation_matrix) <- Cardinal::mz(combined_cardinaldata)
  colnames(correlation_matrix) <- Cardinal::mz(combined_cardinaldata)
  
  return(correlation_matrix)
  
}