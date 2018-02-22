#' @export
#' @title Apply NiftyReg transform found through registration on entire "MSImageSet"
#' @description Applies the NiftyReg transform on the entire MSIimageSet iteratively
#' @param cardinaldata an object of class "MSImageSet"
#' @param niftyregtform the transformation developed during the RNiftyReg registration step
#' @param padding zero padding used around the borders of the image to avoid out of frame pixels
#' @param interpolation, passed to RNiftyReg. Defaults to nearest neighbor
#' @import RNiftyReg
#' @return cardinaldata_registered a registered "MSImageSet"


applyNifty <- function(cardinaldata, niftyregtform, padding=0, interpolation = 0){
  
  ##get padded images as list, apply transformation to image
  image_matrix <- lapply(1:nrow(cardinaldata), function(signal){
    image_matrix <- matrix(0, nrow=max(cardinaldata$x)+padding,ncol=max(cardinaldata$y)+padding)
    
    if(ndim(imageData(cardinaldata)@positionArray) > 3){
      stop("Registration only works on 2-D images and cardinal datasets containing only one sample right now")
    }
    
    if(ndim(imageData(cardinaldata)@positionArray) == 2){
      image_matrix[1:max(cardinaldata$x), 1:max(cardinaldata$y)] <- Cardinal::imageData(cardinaldata)[signal,,]
    }else if(ndim(imageData(cardinaldata)@positionArray) == 3){
      image_matrix[1:max(cardinaldata$x), 1:max(cardinaldata$y)] <- Cardinal::imageData(cardinaldata)[signal,,,]
    }
    
    image_matrix[is.na(image_matrix)] <- 0
    image_matrix<-applyTransform(forward(niftyregtform), image_matrix, nearest=T, interpolation=interpolation)
    return(image_matrix)
  })
  
  ##get transformed images into a data.frame and build MSImageSet from data.frame:
  aa <- data.frame(x=as.vector(row(image_matrix[[1]])), y=as.vector(col(image_matrix[[1]])), sample=sampleNames(cardinaldata))
  
  #matrix to df column
  for(i in 1:length(image_matrix)){
    aa[3+i] <- as.vector(image_matrix[[i]])
  }
  
  #remove empty pixels
  aa <-aa[which(rowSums(aa[4:length(aa)]) != 0),]
  
  #generate new cardinal dataset
  cardinaldata_registered <- MSImageSet(spectra=t(aa[,4:length(aa)]), coord=aa[c("x","y","sample")], 
                                        mz=Cardinal::mz(cardinaldata),
                                        processingData = processingData(cardinaldata),
                                        protocolData = protocolData(cardinaldata),
                                        experimentData = experimentData(cardinaldata))
  
  return(cardinaldata_registered)
  
}