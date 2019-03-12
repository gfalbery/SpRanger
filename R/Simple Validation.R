# Simple Validation #####

Validate <- function(VirusAssocs, Network){

  a = 1

  ValidList <- list()

  for(a in a:length(VirusAssocs)){

    print(names(VirusAssocs)[a])

    pHosts <- VirusAssocs[[a]]

    pHosts2 <- intersect(pHosts, rownames(Network))

    if(length(pHosts2)>1){

      FocalNet <- Network[pHosts2,]

      ValidEst <- list()

      for(b in pHosts2){

        pHosts4 <- setdiff(pHosts2, b)

        pHosts3 <- setdiff(colnames(FocalNet), pHosts4)

        Estimates <- FocalNet[pHosts4, pHosts3]

        if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)

        Ests <- data.frame(Sp = names(sort(colSums(Estimates), decreasing = T)),
                           Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates)) %>%
          mutate(Focal = ifelse(Sp==b, 1, 0),
                 Iteration = b)

        rownames(Ests) <- Ests$Sp

        ValidEst[[b]] <- Ests

      }

      ValidList[[names(VirusAssocs)[a]]] <- ValidEst

    } else {

      ValidList[[names(VirusAssocs)[a]]] <- NA

      print("Hosts Not Found!")

    }
  }

  return(ValidList)

}
