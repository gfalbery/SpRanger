
# Network Prediction only ####

NetworkPredict <- function(HostList, Network){

  HostList
  pHosts <- intersect(HostList, rownames(Network))

  if(length(pHosts)>1){

    FocalNet <- Network[pHosts,]

    ValidEst <- list()

    pHosts2 <- setdiff(colnames(FocalNet), pHosts)

    Estimates <- FocalNet[pHosts, pHosts2]

    if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)

    ValidEst <- tibble(Sp = names(sort(colSums(Estimates), decreasing = T)),
                       Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates))

    ValidEst$Rank <- nrow(ValidEst) - rank(ValidEst$Count, ties.method = "average") + 1

  } else {

    ValidEst <- NA

    print("Hosts Not Found!")

  }

  return(ValidEst)

}
