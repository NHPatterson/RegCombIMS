#' @export
#' @title Reduce Cardinal dataset's coordinates to minimum
#' @description reduceCardinalCoord is used to set the minimum y and x coordinate values in a Cardinal dataset to zero. 
#' @param cardinaldata an object of class "MSImageSet"
#' @return cardinaldata

reduceCardinalCoord <- function(cardinaldata){
  
  #grab coordinate data
  df <- coord(cardinaldata)[1:2]
  #placeholder if there is no samplename
  df$sample <- "placeholder"
  #split by sample
  dfs <- split(df, df$sample)
  
  #loop through samples and reduce each's coordinates
  for(i in 1:length(dfs)){
    #store original x & y coordinates
    dfs[[i]]$origx <- dfs[[i]]$x
    dfs[[i]]$origy <- dfs[[i]]$y
    
    #reduced coordinates
    dfs[[i]]$x <- dfs[[i]]$x - (min(dfs[[i]]$x) - 1)
    dfs[[i]]$y <- dfs[[i]]$y - (min(dfs[[i]]$y) - 1)
    
  }
  
  #recombine sample split dataframe
  df <-do.call(rbind, dfs)
  
  #store old coordinate values and add new to cardinal dataset
  cardinaldata$origx <- df$origx
  cardinaldata$origy <- df$origy
  cardinaldata$x     <- df$x
  cardinaldata$y     <- df$y
  
  #regenerate the necessary cardinal metadata
  
  cardinaldata <- regeneratePositions(cardinaldata)
  
  return(cardinaldata)
  
}