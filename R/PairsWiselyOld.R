

PairsWiselyOld <- function(Rasterstack){

  Valuedf <- data.frame(getValues(Rasterstack))
  Valuedf2 <- reshape2::melt(Valuedf)
  Valuedf2$x <- rep(1:Rasterstack[[1]]@ncols, Rasterstack[[1]]@nrows)
  Valuedf2$y <- rep(Rasterstack[[1]]@nrows:1, each = Rasterstack[[1]]@ncols)

  Rangedf <- Valuedf2 %>%
    dplyr::rename(Host = variable, Presence = value)

  Rangedf$GridID <- with(Rangedf, paste(x, y))

  Range0 <- levels(Rangedf$Host)[which(table(Rangedf$Host)==0)] # Hosts that have no spatial records??
  Rangedf <- droplevels(Rangedf)
  Rangedf <- Rangedf[order(Rangedf$Host),]

  # Using igraph to project it

  RangeOverlap <- matrix(0, nrow = nlevels(Rangedf$Host), ncol = nlevels(Rangedf$Host))
  dimnames(RangeOverlap) <- list(levels(Rangedf$Host),levels(Rangedf$Host))

  for(x in levels(Rangedf$Host)){

    if(x == first(levels(Rangedf$Host))) t1 <- Sys.time()

    Grids <- Rangedf[Rangedf$Host==x,"GridID"]
    SubRangedf <- Rangedf[Rangedf$GridID %in% Grids,]

    RangeOverlap[x,] <- table(SubRangedf$Host)

    print(x)

    if(x == last(levels(Rangedf$Host))) t2 <- Sys.time()

  }

  RangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
  RangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))

  RangeAdj <- RangeOverlap/(RangeA + RangeB - RangeOverlap)

  TimeTaken = Sys.time() - t1

  print(TimeTaken)

  return(RangeAdj)

}
