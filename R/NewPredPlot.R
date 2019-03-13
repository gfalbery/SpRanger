

PredPlot <- function(Virus = NULL,
                     HostList = NULL,
                     Threshold = 10,
                     Validate = TRUE,
                     FocalDisplay = c("Observed","Predicted"),
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

      if(Validate){

        ValidDF <- NetworkValidate(HostList = HostList, Network = AllSums)

        FullDF <- as.data.frame(ValidDF)

        SubValidDF <- dplyr::filter(ValidDF, Focal == "Predicted")
        SubValidDF <- SubValidDF %>% mutate(SubRank = nrow(SubValidDF) - rank(Count) + 1)

        ValidDF <- left_join(ValidDF, SubValidDF)
        ValidDF <- ValidDF %>% mutate(Include = ifelse(Focal == "Observed"|SubRank<Threshold, 1, 0)) %>%
          filter(Include == 1)

      } else {

        ValidDF <- NetworkPredict(HostList = HostList, Network = AllSums)

        FullDF <- as.data.frame(ValidDF)

        ValidDF <- ValidDF %>% filter(Rank<=Threshold)

      }

    } else stop("Hosts Not Found!")

    SpList <- list()

    if("Observed"%in%FocalDisplay) SpList$Observed <- HostList
    if("Predicted"%in%FocalDisplay) SpList$Predicted <- setdiff(ValidDF$Sp, HostList)

  }

  VirusName <- str_replace_all(Virus, "_", " ")

  if(Facet&length(FocalDisplay)==2){

    Plots <- list()

    for(x in FocalDisplay){

      Title = x

      PlotSp <- SpList[[x]]

      Plots[[x]] <-

        RangePlot(Mammals = PlotSp, Map = Map, Tree = Tree) +
        ggtitle(as.character(paste(Title, "Hosts"))) +
        theme(legend.position = Legend)

    }

    MapPlot <- Plots %>% arrange_ggplot2(nrow = 2)

  } else {

    MapPlot <- RangePlot(Mammals = unlist(SpList),
                         Map = Map,
                         Tree = Tree) +
      theme(legend.position = Legend)

  }

  ReturnList <- list()

  ReturnList[["MapPlot"]] <- MapPlot

  if(Validate){

    FullDF$Focal <- FullDF$Focal %>% factor(levels = c("Predicted", "Observed"))

    ValidPlot <- ggplot(FullDF, aes(Focal, Count, colour = Focal, alpha = Focal)) +
      ggforce::geom_sina() +
      scale_alpha_manual(values = c(0.3, 1)) +
      geom_text(data = FullDF[1,], inherit.aes = F, aes(x = 1.5, y = max(FullDF$Count)*1.1,
                                                        label = paste("Mean Focal Rank =",
                                                                      mean(FullDF[FullDF$Focal=="Observed","Rank"])))) +
      ggtitle("Model Success Rate")

    ReturnList[["ValidPlot"]] <- ValidPlot
  }

  if(Summarise){

    SummariseDF <- Panth1 %>%
      dplyr::select(Sp, hOrder, MSW05_Family) %>%
      dplyr::rename(Order = hOrder,
                    Family = MSW05_Family) %>%
      filter(Sp%in%unlist(SpList)) %>%
      dplyr::rename(Species = Sp)

    ReturnList[["Summarise"]] <- SummariseDF

  }

  return(ReturnList)

}
