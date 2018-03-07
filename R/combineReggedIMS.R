#' @export
#' @title Combine registered IMS datasets into a single format
#' @description Combines datasets with an mz_offset (soon:or with a featureData label)
#' @param cardinaldata1 an object of class "MSImageSet", first dataset
#' @param cardinaldata2 an object of class "MSImageSet", second dataset
#' @param mz_offset offset of the m/z axis between datasets
#' @param sample_name passed as 'sample' in resultant "MSImageSet"
#' @return a combined offset "MSImageSet"


combineReggedIMS <- function(cardinaldata1,cardinaldata2,mz_offset=13000,sample_name="combined_ds"){
  
  #create bounding boxes around datasets
  cardinaldata1 <-boundBox(cardinaldata1, xmax=max(c(max(cardinaldata1$x),max(cardinaldata2$x))),ymax=max(c(max(cardinaldata1$y),max(cardinaldata2$y))))
  cardinaldata2 <-boundBox(cardinaldata2, xmax=max(c(max(cardinaldata1$x),max(cardinaldata2$x))),ymax=max(c(max(cardinaldata1$y),max(cardinaldata2$y))))

  mz_cardinaldata1 <- Cardinal::mz(cardinaldata1)
  mz_cardinaldata2 <- Cardinal::mz(cardinaldata2)
  
  binary_mask_cardinaldata1 <- cardinaldata1$binary_mask
  binary_mask_cardinaldata2 <- cardinaldata2$binary_mask
  
  #sort coordinates, and data
  cardinaldata1_df <- arrange(cbind(pData(cardinaldata1)[,1:3],t(iData(cardinaldata1)),binary_mask_cardinaldata1),x,y)
  cardinaldata2_df <- arrange(cbind(pData(cardinaldata2)[,1:3],t(iData(cardinaldata2)),binary_mask_cardinaldata2),x,y)
  
  keep_idx = which(cardinaldata1_df$binary_mask + cardinaldata2_df$binary_mask == 2)
  
  
  #organize 
  cardinaldata1_2_df <- cbind(cardinaldata1_df[1:(ncol(cardinaldata1_df)-1)][keep_idx,], cardinaldata2_df[4:(ncol(cardinaldata2_df) -1)][keep_idx,])
  levels(cardinaldata1_2_df$sample)[1] <- sample_name
  
  ##make new MSIageSet:
  combined <- MSImageSet(spectra=t(cardinaldata1_2_df[,!names(cardinaldata1_2_df) %in% c("x","y","sample")]), 
                          coord=cardinaldata1_2_df[c("x","y","sample")], 
                          mz=c(mz_cardinaldata1,mz_cardinaldata2+mz_offset),
                          processingData = processingData(cardinaldata1),
                          protocolData = protocolData(cardinaldata1),
                          experimentData = experimentData(cardinaldata1))
  
  
  return(combined)
}