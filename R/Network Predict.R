
# Network Prediction only ####

NetworkPredict <- function(HostList, Network, Fun = colSums, IncludeObserved = F){

  require(Matrix)

  pHosts <- intersect(HostList, rownames(Network))

  if(length(pHosts)>0){

    ValidEst <- list()

    if(!IncludeObserved){

      pHosts2 <- setdiff(colnames(Network), pHosts)

    }else{

      pHosts2 <- colnames(Network)

    }

    Estimates <- Network[pHosts, pHosts2] #%>% as.matrix

    if(is.null(dim(Estimates))){
      Estimates <- t(data.frame(Estimates))
    }

    PredictFunction <- Fun

    ValidEst <- tibble(Sp = names(sort(PredictFunction(Estimates), decreasing = T)),
                       Count = sort(PredictFunction(Estimates), decreasing = T)/nrow(Estimates))

    ValidEst$Rank <- nrow(ValidEst) - rank(ValidEst$Count, ties.method = "average") + 1

  } else {

    ValidEst <- NA

    print("Hosts Not Found!")

  }

  if(IncludeObserved){

    ValidEst$Observed <- as.numeric(ValidEst$Sp %in% HostList)

  }

  return(ValidEst)

}
