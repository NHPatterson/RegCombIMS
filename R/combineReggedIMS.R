#' @export
#' @title Combine registered IMS datasets into a single format
#' @description Combines datasets with an mz_offset (soon:or with a featureData label)
#' @param cardinaldata1 an object of class "MSImageSet", first dataset
#' @param cardinaldata2 an object of class "MSImageSet", second dataset
#' @param ds1_name name added to the data under @featureData$ds_origin for first dataset
#' @param ds1_name name added to the data under @featureData$ds_origin for second dataset
#' @param combined_name passed as 'sample' in resultant "MSImageSet"
#' @return a combined  "MSImageSet"


combineReggedIMS <- function(cardinaldata1,cardinaldata2,ds1_name = "dataset1",ds2_name="dataset2" ,combined_name="combined_ds"){
  
  ##determine which pixels are retained during combination (must be present in both datasets)
  #use 'X..Y..' for unique indexing
  ds1_xy <- paste0('X',cardinaldata1$x,'Y',cardinaldata1$y)
  ds2_xy <- paste0('X',cardinaldata2$x,'Y',cardinaldata2$y)
  
  ds_intersection <- intersect(ds1_xy,ds2_xy)
  
  ds1_retained_pixels <- which(ds1_xy %in% ds_intersection)
  ds2_retained_pixels <- which(ds2_xy %in% ds_intersection)
  
  #subset retained pixels
  cardinaldata1 <- cardinaldata1[,ds1_retained_pixels]
  cardinaldata2 <- cardinaldata2[,ds2_retained_pixels]
  
  if ('ds_origin' %in% varLabels(cardinaldata1@featureData)){
  } else {
    cardinaldata1@featureData$ds_origin <- ds1_name
  }
  
  cardinaldata2@featureData$ds_origin <- ds2_name

  mz_sortvector <- sort(c(Cardinal::mz(cardinaldata1),Cardinal::mz(cardinaldata2)),index.return=T)

  ##make new MSIageSet:
  combined <- MSImageSet(spectra=rbind(iData(cardinaldata1), 
                                         iData(cardinaldata2))[mz_sortvector$ix,], 
                          coord=data.frame(coord(cardinaldata1)[c("x","y")],sample=combined_name,row.names=NULL), 
                          mz=mz_sortvector$x,
                          processingData = processingData(cardinaldata1),
                          protocolData = protocolData(cardinaldata1),
                          experimentData = experimentData(cardinaldata1))
  
  
  combined@featureData$ds_origin = c(cardinaldata1@featureData$ds_origin, cardinaldata2@featureData$ds_origin)[mz_sortvector$ix]
  
  return(combined)
}
