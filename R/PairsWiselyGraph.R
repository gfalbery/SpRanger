
# Greg's range extraction function ####

PairsWiselyGraph <- function(Rasterstack, Species = "All"){

  library(raster); library(tidyverse); library(Matrix); library(igraph)

  t1 <- Sys.time()

  print("Getting the grid values")

  Valuedf <- data.frame(getValues(Rasterstack)) %>% as.matrix
  Valuedf[is.na(Valuedf)] <- 0

  Valuedf <- Valuedf %>% as("dgCMatrix")

  bipgraph <- graph.incidence(Valuedf, weighted = NULL)

  Spgraph <- bipartite.projection(bipgraph)$proj1

  if(Species != "All") Valuedf <- Valuedf[,Species]

  if(any(colSums(Valuedf)==0)){ print("Removing some species with no ranging data :(")

    Valuedf <- Valuedf[,-which(colSums(Valuedf)==0)]

  }

  RangeOverlap <- matrix(NA, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
  dimnames(RangeOverlap) <- list(colnames(Valuedf),colnames(Valuedf))

  print("Calculating Overlap")

  for(x in 1:(ncol(Valuedf)-1)){

    print(colnames(Valuedf)[x])

    TrainGrids <- Valuedf[,x]

    SubRangedf <- Valuedf[which(TrainGrids==1),x:ncol(Valuedf)]

    if(!is.null(dim(SubRangedf))){

      RangeOverlap[x,x:ncol(Valuedf)] <- apply(SubRangedf, 2, function(a) length(which(a==1)))

    } else   RangeOverlap[x,x:ncol(Valuedf)] <- SubRangedf

  }

  x = ncol(Valuedf)

  print(colnames(Valuedf)[x])

  TrainGrids <- Valuedf[,x]

  SubRangedf <- Valuedf[which(TrainGrids==1),x:ncol(Valuedf)]

  RangeOverlap[x,x:ncol(Valuedf)] <- length(which(SubRangedf==1))

  FullRangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
  FullRangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))

  RangeAdj <- RangeOverlap/(FullRangeA + FullRangeB - RangeOverlap) # Weighted evenly

  TimeTaken = Sys.time() - t1

  print(TimeTaken)

  diag(RangeAdj) <- NA

  RangeAdj[lower.tri(RangeAdj)] <- t(RangeAdj)[!is.na(t(RangeAdj))]
  diag(RangeAdj) <- 1

  return(RangeAdj)

  detach(package:raster)

}
