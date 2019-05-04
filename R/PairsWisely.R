PairsWisely <- function(Rasterstack, Species = "All"){
  
  library(raster)
  library(tidyverse)
  library(Matrix)
   
  print("Making data frame!")
  
  if(class(Rasterstack)=="RasterBrick"){
    
    t1 <- Sys.time()
    print("Getting the grid values")
    Valuedf <- data.frame(raster::getValues(Rasterstack)) %>% as.matrix
    
  } else{
    
    Valuedf <- lapply(1:length(Rasterstack), function(a){
      
      getValues(Rasterstack[[a]])
      
    }) %>% bind_cols %>% as.data.frame()
    
  }
  
  Valuedf[is.na(Valuedf)] <- 0
  Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")
  
  print(paste0("Data frame size = ", dim(Valuedf)))
  
  if (Species != "All"){
    Valuedf <- Valuedf[, Species]
  }
  
  print(paste0("Data frame size = ", dim(Valuedf)))
  
  if (any(Matrix::colSums(Valuedf) == 0)) {
    print("Removing some species with no ranging data :(")
    Valuedf <- Valuedf[, -which(Matrix::colSums(Valuedf) == 0)]
  }
  
  RangeOverlap <- matrix(NA, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
  dimnames(RangeOverlap) <- list(colnames(Valuedf), colnames(Valuedf))
  print("Calculating Overlap")
  for (x in 1:(ncol(Valuedf) - 1)) {
    print(colnames(Valuedf)[x])
    TrainGrids <- Valuedf[, x]
    SubRangedf <- Valuedf[which(TrainGrids == 1), x:ncol(Valuedf)]
    if (!is.null(dim(SubRangedf))) {
      RangeOverlap[x, x:ncol(Valuedf)] <- apply(SubRangedf,
                                                2, function(a) length(which(a == 1)))
    }
    else RangeOverlap[x, x:ncol(Valuedf)] <- SubRangedf
  }
  
  x = ncol(Valuedf)
  print(colnames(Valuedf)[x])
  TrainGrids <- Valuedf[, x]
  SubRangedf <- Valuedf[which(TrainGrids == 1), x:ncol(Valuedf)]
  RangeOverlap[x, x:ncol(Valuedf)] <- length(which(SubRangedf ==
                                                     1))
  FullRangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)),
                      nrow(RangeOverlap))
  FullRangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)),
                      nrow(RangeOverlap))
  RangeAdj <- RangeOverlap/(FullRangeA + FullRangeB - RangeOverlap)
  TimeTaken = Sys.time() - t1
  print(TimeTaken)
  diag(RangeAdj) <- NA
  RangeAdj[lower.tri(RangeAdj)] <- t(RangeAdj)[!is.na(t(RangeAdj))]
  diag(RangeAdj) <- 1
  
  RangeAdj
}
