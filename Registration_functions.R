
##Functions to do registration and combination of IMS datasets: 

#reduces the coordinates of a Cardinal dataset to minimum, only works for one sample at a time in this implementation
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
  
  #regenerate the necessary cardinal metdata
  
  cardinaldata <- regeneratePositions(cardinaldata)
  
  return(cardinaldata)
  
}

#converts a Cardinal ResultSet image to array, ensures that the dimensions are the same between the datasets
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

#zero fills pixels outside the the imaging area for combination later
bound_box <- function(cardinaldata, xmax=300, ymax=250, center=F){
  
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

#applies a transformation to each ion image iteratively and rebuilds the final dataset
#interpolation = 0 is nearest neighbors, see RNiftyReg for other options
ApplyNifty <- function(cardinaldata,niftyregtform,padding=0, interpolation = 0){
  
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

#combines two registered datasets
combineReggedIMS <- function(cardinaldata1,cardinaldata2,mz_offset=13000,sample_name="combined_ds"){
  
  #create bounding boxes around datasets
  cardinaldata1 <-bound_box(cardinaldata1, xmax=max(c(max(cardinaldata1$x),max(cardinaldata2$x))),ymax=max(c(max(cardinaldata1$y),max(cardinaldata2$y))))
  cardinaldata2 <-bound_box(cardinaldata2, xmax=max(c(max(cardinaldata1$x),max(cardinaldata2$x))),ymax=max(c(max(cardinaldata1$y),max(cardinaldata2$y))))

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




