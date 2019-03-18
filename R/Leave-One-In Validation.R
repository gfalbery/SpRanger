
NetworkValidate2 <- function(HostList, Network){

  pHosts <- HostList
  pHosts2 <- intersect(pHosts, rownames(Network))

  if(length(pHosts2) > 1) {

    FocalNet <- Network[pHosts2,]

    ValidEst <- list()

    for(b in pHosts2) {

      pHosts3 <- setdiff(colnames(FocalNet), b)

      Estimates <- FocalNet[b,pHosts3]

      Ests <- data.frame(Sp = pHosts3,
                         Count = Estimates) %>%
        mutate(Focal = ifelse(Sp %in% pHosts2, "Observed", "Predicted"), Iteration = b) %>%
        slice(order(Count))

      ValidEst[[b]] <- Ests

    }

    ValidEst <- ValidEst %>% bind_rows %>% group_by(Sp) %>%
      dplyr::summarise(Count = mean(Count)) %>%
      slice(order(Count, decreasing = T)) %>%
      mutate(Focal = ifelse(Sp %in% pHosts2, "Observed", "Predicted"))

    ValidEst$Rank <- nrow(ValidEst) - rank(ValidEst$Count, ties.method = "average") + 1

  } else {

    ValidEst <- NA
    print("Hosts Not Found!")

  }

  return(ValidEst)

}
