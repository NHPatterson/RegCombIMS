#' @export
#' @title Query correlation matrix for top markers accross combined data sets
#' @description 
#' @param combined_cardinaldata an object of class "MSImageSet", combined dataset
#' @param correlation_matrix object generated from getCorrelationMatrix
#' @param query_dataset name of the dataset in 'ds_origin' that we wish to query for
#' @param query_against name of the dataset in 'ds_origin' that we will query against, can be multiple, passing 'all' will query all datasets
#' @param query_mz mz in the query_dataset one wishes to find correlations to in query_against dataset(s)
#' @param plot_ions boolean, whether to make a plot of top ions
#' @param top_n number of ions to plot
#' @param ... arguments passed to image, useful for passing layout=c() of image with multiple queries
#' @return data.frame of mz and r, sorted in descending order by r for the query_mz, first value is always mz_query

correlationQuery <- function(combined_cardinaldata, correlation_matrix,
                              query_dataset = 'mydata', query_against = c('otherdata','all'), 
                              query_mz = 885.55, plot_ions=T, top_n = 9,...){
  
  if (query_dataset %in% unique(combined_cardinaldata@featureData$ds_origin) == F){
    print('Error : queried dataset not found in @featureData$ds_origin')
    stop()
  }
  
  if (query_against != 'all' & query_against %in% unique(combined_cardinaldata@featureData$ds_origin) == F){
    print('Error: dataset to query for not found in @featureData$ds_origin')
    stop()
  }
  
  if (query_against == 'all'){
    
    peak_indices = !combined_cardinaldata@featureData$ds_origin %in% query_dataset
    
  } else{
    
    peak_indices = combined_cardinaldata@featureData$ds_origin %in% query_against
    
  }
  
  query_return = correlation_matrix[features(combined_cardinaldata, mz= query_mz),peak_indices]
  
  top_mzs = as.numeric(names(sort(query_return, decreasing = T)))
  top_rs = as.numeric(sort(query_return, decreasing = T))
  
  if (plot_ions == T){
    image(combined_cardinaldata, mz=query_mz ,strip=F, main = paste('Queried ion', 'm/z :', query_mz ),...)
    
    
    for (i in 1:top_n){
      image(combined_cardinaldata, mz=top_mzs[i], main=paste0(i, ",m/z:",round(top_mzs[i],3), ',r=',  round(top_rs[i],3)), strip=F)
      
    }
  }
  
  return(data.frame(mz = c(query_mz,top_mzs), r = c(1,top_rs)))
}


  