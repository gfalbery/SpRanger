
IntersectGet <-
  function(Rasterstack, PairList, Names = "All"){

    library(raster); library(tidyverse); library(Matrix)

    t1 <- Sys.time()
    print("Getting the grid values")

    if(class(Rasterstack)=="RasterBrick"){

      Valuedf <- data.frame(raster::getValues(Rasterstack)) %>% as.matrix

    }else{

      Valuedf <- lapply(1:length(Rasterstack), function(a){

        getValues(Rasterstack[[a]])

      }) %>% bind_cols %>% as.data.frame()

    }

    colnames(Valuedf) <- names(Rasterstack)

    Valuedf[is.na(Valuedf)] <- 0
    Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")

    print(paste0("Data frame size = ", dim(Valuedf)))

    if (Species != "All"){
      Valuedf <- Valuedf[, Species]
    }

    if (any(Matrix::colSums(Valuedf) == 0)) {
      print("Removing some species with no ranging data :(")
      Valuedf <- Valuedf[, -which(Matrix::colSums(Valuedf) == 0)]
    }

    print(paste0("Data frame size = ", dim(Valuedf)))

    lapply(1:nrow(PairList), function(b){

      Sp1 <- PairList[1,"Sp.1"]
      Sp2 <- PairList[1,"Sp.2"]

      Vector1 <- Valuedf[,Sp1]
      Vector1[Valuedf[,Sp2]==0] <- 0

      Vector1 %>% return

    }) -> ValuesVectors#

    names(ValuesVectors) <- paste(PairList[,1],
                                  PairList[,2],
                                  sep = ".")

    return(ValuesVectors)

  }
