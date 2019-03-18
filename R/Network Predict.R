
# Network Prediction only ####

NetworkPredict <- function(HostList, Network){

  pHosts <- intersect(HostList, rownames(Network))

  if(length(pHosts)>0){

    ValidEst <- list()

    pHosts2 <- setdiff(colnames(Network), pHosts)

    Estimates <- Network[pHosts, pHosts2]

    if(is.null(dim(Estimates))){
      Estimates <- rbind(Estimates, Estimates)
    }

    ValidEst <- tibble(Sp = names(sort(colSums(Estimates), decreasing = T)),
                       Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates))

    ValidEst$Rank <- nrow(ValidEst) - rank(ValidEst$Count, ties.method = "average") + 1

  } else {

    ValidEst <- NA

    print("Hosts Not Found!")

  }

  return(ValidEst)

}
