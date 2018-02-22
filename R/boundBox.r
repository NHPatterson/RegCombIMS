#' @export
#' @title Bounding box of zeros around Cardinaldata
#' @description Adds zeros around "MSImageSet" object to so that datasets can be combined after registration
#' @param cardinaldata an object of class "MSImageSet"
#' @param xmax maximum of bounding box in x
#' @param ymax maximum of bounding box in y
#' @param boolean, whether to center the imaging data in the bounding box
#' @return msset an object of class "MSImageSet"

boundBox <- function(cardinaldata, xmax=300, ymax=250, center=F){
  
  #size of zero box
  xseq <- c(1:xmax)
  yseq <- c(1:ymax)
  
  #grab pixel and intensity data
  intensity_data <- data.frame(pData(cardinaldata)[1:2],t(iData(cardinaldata)))
  
  #make empty data.frame expanded to size of bounding box
  empty_df<-data.frame(expand.grid(xseq,yseq),sampleNames(cardinaldata))
  empty_df[4:(nrow(cardinaldata)+3)] <- 0
  colnames(empty_df) <- c("x","y","sample")
  
  #sort pixels by x then y
  empty_df <- arrange(empty_df,x,y)
  intensity_data <- arrange(intensity_data,x,y)
  
  #make character vector of pixel names combining x & y
  zero_pixels  <- paste0("X",empty_df$x,"Y",empty_df$y) 
  
  #center data in bounding box
  if(center == T){
    mpx <- round(max(pData(cardinaldata)$x) / 2, 0)
    mpy <- round(max(pData(cardinaldata)$y) / 2, 0)
    
    mpx <- (xmax / 2) - mpx
    mpy <- (xmax / 2) - mpy
    
    cdata_pixels <- paste0("X",intensity_data$x + mpx,"Y",intensity_data$y + mpy)
    
  } else {
    
    #bounding box extends after the data
    
    cdata_pixels <- paste0("X",intensity_data$x,"Y",intensity_data$y)
    
  }
  
  #fill empty pixels with intensity data
  non_zero_idx = which(zero_pixels %in% cdata_pixels)
  empty_df[non_zero_idx, 4:length(empty_df)] <- intensity_data[3:length(intensity_data)]
  
  #make cardinal data
  msset <- MSImageSet(spectra=t(empty_df[,4:length(empty_df)]), coord=empty_df[c("x","y","sample")], 
                      mz=Cardinal::mz(cardinaldata),
                      processingData = processingData(cardinaldata),
                      protocolData = protocolData(cardinaldata),
                      experimentData = experimentData(cardinaldata))
  
  msset$binary_mask = 0
  msset$binary_mask[non_zero_idx] = 1
  
  return(msset)
}