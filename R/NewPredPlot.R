

PredPlot <- function(Virus = NULL,
                     HostList = NULL,
                     Threshold = 10,
                     Validate = TRUE,
                     FocalDisplay = c(1,0),
                     Facet = FALSE,
                     Map = TRUE,
                     Tree = TRUE,
                     Summarise = TRUE,
                     Legend = "right"){

  require(dplyr); require(stringr); library(tidyverse); require(ggtree); require(cowplot)

  if(!is.null(Virus)&!is.null(HostList)) stop("Only specify one of Virus or Host please :)")

  if(!is.null(Virus)){

    Virus <- Virus %>% str_replace_all(" ", "_")

    if(!Virus%in%names(GAMValid)) stop("Virus not found :( try entering hosts manually?") else print("We have this virus! Wahey.")

    Df <- GAMValid[[Virus]]

    Df$Rank <- nrow(Df) - rank(Df$Count)

    PredHosts <- Df %>% filter(Rank < Threshold) #%>% select(Sp) %>% unlist

  } else if(!is.null(HostList)){

    HostList <- HostList %>% str_replace_all(" ", "_")

    if(length(HostList)<2) print("This is only one host! Predictions might be rubbish.", immediate = T)

    if(length(intersect(HostList, rownames(AllSums)))==0){

      stop("None of your hosts are in our dataset! :( are they marine, or are there synonyms?")

    } else if(length(intersect(HostList, rownames(AllSums)))<length(HostList)){

      print(paste("Warning: Some hosts not found:",
                  paste(setdiff(HostList, rownames(AllSums)), collapse = ", "),
                  "; tread carefully"))
    }

    pHosts2 <- intersect(HostList, rownames(AllSums))

    if(length(pHosts2)>0){

      FocalNet <- AllSums[pHosts2,]

      ValidEst <- list()

      if(Validate){

        for(b in pHosts2){

          pHosts4 <- setdiff(pHosts2, b)

          pHosts3 <- setdiff(colnames(FocalNet), pHosts4)

          Estimates <- FocalNet[pHosts4, pHosts3]/length(AllSims)

          if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)

          Ests <- data.frame(Sp = names(sort(colSums(Estimates), decreasing = T)),
                             Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates)) %>%
            mutate(Focal = ifelse(Sp==b, 1, 0),
                   Iteration = b)

          rownames(Ests) <- Ests$Sp

          ValidEst[[b]] <- Ests

        }

        ValidDF <- ValidEst %>%
          bind_rows() %>%
          group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>%
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = as.factor(Focal))

      } else {

        Estimates <- FocalNet[pHosts2,]/length(AllSims)

        if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)

        Ests <- data.frame(Sp = names(sort(colSums(Estimates), decreasing = T)),
                           Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates)) %>%
          mutate(Focal = ifelse(Sp%in%HostList,1,0))

        rownames(Ests) <- Ests$Sp

        ValidDF <- Ests %>%
          slice(order(Focal)) %>%
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = as.factor(Focal))

      }

    } else print("Hosts Not Found!")

    Df <- ValidDF %>% as.data.frame()
    Df$Rank <- nrow(Df) - rank(Df$Count)
    if(1 %in% FocalDisplay) Df[Df$Focal==1,"Include"] <- 1
    if(0 %in% FocalDisplay) Df <- Df %>% mutate(Include = ifelse(Rank<Threshold&Focal==0,1,0))
    PredHosts <- Df %>% filter(Include == 1) #%>% select(Sp) %>% unlist

  }

  Df <- Df %>% mutate(Focal = factor(c("Predicted","Observed")[(3-as.numeric(Focal))]))

  VirusName <- str_replace_all(Virus, "_", " ")

  if(Facet){

    Plots <- list()

    for(x in FocalDisplay){

      Plots[[x]] <- Df %>% filter(Focal = c("Observed","Predicted")[FocalDisplay[x+1]]) %>%
        select(Sp) %>% RangePlot(Map = Map, Tree = Tree) +
        ggtitle(c("Observed","Predicted")[FocalDisplay[x+1]])

    }

    MapPlot <- Plots %>% arrange_ggplot2(nrow = 2)

  } else {

    MapPlot <- Df %>% select(Sp) %>% RangePlot(Map = Map,
                                               Tree = Tree)

  }

  ReturnList <- list()

  ReturnList[["MapPlot"]] <- MapPlot

  if(Validate){

    ValidPlot <- ggplot(Df, aes(Focal, Count, colour = Focal, alpha = Focal)) +
      ggforce::geom_sina() +
      scale_alpha_manual(values = c(0.3, 1)) +
      geom_text(data = Df[1,], inherit.aes = F, aes(x = 1.5, y = max(Df$Count)*1.1,
                                                    label = paste("Mean Focal Rank =",
                                                                  mean(Df[Df$Focal=="Observed","Rank"])))) +
      ggtitle("Model Success Rate")

    ReturnList[["ValidPlot"]] <- ValidPlot
  }

  if(Summarise){

    SummariseDF <- Panth1 %>% select(Sp, hOrder, MSW05_Family) %>%
      dplyr::rename(Order = hOrder, Family = MSW05_Family) %>%
      right_join(PredHosts, by = "Sp") %>%
      dplyr::rename(Species = Sp)

    ReturnList[["Summarise"]] <- SummariseDF

  }

  return(ReturnList)

}
