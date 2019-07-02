
# PairsWisening ####

PairsWisening <- function(RasterStack1, RasterStack2, Species = "All", Area = F){

  library(raster)
  library(tidyverse)
  library(Matrix)

  t1 <- Sys.time()
  print("Getting started")

  Shared <- intersect(names(RasterStack1), names(RasterStack2))
  Shared <- sort(Shared)
  Total <- length(RasterStack1)+length(RasterStack2)

  if(length(Shared)<(Total)/2){
    print(paste("Losing", Total-length(Shared)))
  }

  if (Species != "All"){
    Shared <- intersect(Species, Shared)
  }

  if(Area){AreaFun <- function(a){sum(a, na.rm = T)} } else {AreaFun <- function(a){length(which(a>0))}}

  RangeOverlap <-
    FullRangeA <-
    FullRangeB <- as.list(rep(NA, length(Shared)))

  for(i in 1:length(Shared)){

    print(Shared[i])

    Values1 <- values(RasterStack1[[Shared[i]]])
    Values2 <- values(RasterStack2[[Shared[i]]])

    TrainGrids <- Values1
    SubRangedf <- Values2[which(TrainGrids > 0)]
    RangeOverlap[[Shared[i]]] <- AreaFun(SubRangedf)
    FullRangeA[[Shared[i]]] <- AreaFun(Values1)
    FullRangeB[[Shared[i]]] <- AreaFun(Values2)

  }

  RangeAdj <-
    unlist(RangeOverlap)/(unlist(FullRangeA) + unlist(FullRangeB) - unlist(RangeOverlap))

  names(RangeAdj) <- Shared

  t2 <- Sys.time()

  print(t2 - t1)

  DF1 <- data.frame(
    Sp = Shared,
    Stack1 = FullRangeA %>% unlist,
    Stack2 = FullRangeB %>% unlist,
    Overlap = RangeOverlap %>% unlist,
    POverlap = RangeAdj

  )

  return(DF1)

}
