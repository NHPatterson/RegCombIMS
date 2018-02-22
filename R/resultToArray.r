#' @export
#' @title Cardinal resultSet to R array
#' @description converts a Cardinal ResultSet data (stored as 1D vector) to 2D array, for now, only works for SpatialShrunkenCentroids and PCA
#' @param cardinaldata an object of class "MSImageSet"
#' @param resultdata an object of class "ResultSet"
#' @return image_matrix

resultToArray<- function(cardinaldata, resultdata, idx){
  
  image_matrix <- matrix(0, nrow=max(cardinaldata$x),ncol=max(cardinaldata$y))
  
  xy <- coord(cardinaldata)[,1:2]
  
  if(class(resultdata)[1] == "SpatialShrunkenCentroids"){
    resultdf <- data.frame(x= pData(resultdata)$x, y = pData(resultdata)$y, resultdata@resultData[[1]]$probabilities)
  }
  
  if(class(resultdata)[1] == "PCA"){
    resultdf <- data.frame(x= pData(resultdata)$x, y = pData(resultdata)$y, resultdata@resultData[[1]]$scores)
  }
  
  
  for(i in 1:nrow(xy)){
    image_matrix[xy[i,1],xy[i,2]] <- resultdf[i,idx+2]
  }
  
  return(image_matrix)
}