

PredPlot <- function(Virus = NULL,
                     Hosts = NULL,
                     Threshold = 10,
                     Focal = c(1,0),
                     Facet = FALSE){

  require(dplyr); require(stringr); library(tidyverse); require(ggtree)

  if(!is.null(Virus)&!is.null(Hosts)) stop("Only specify one of Virus or Host please :)")

  if(!is.null(Virus)){

    Virus <- Virus %>% str_replace_all(" ", "_")

    if(!Virus%in%names(GAMValid)) stop("Virus not found :( try entering hosts manually?") else print("We have this virus! Wahey.")

    Df <- GAMValid[[Virus]]

    Df$Rank = nrow(Df) - rank(Df$Count)

    PredHosts <- Df %>% filter(Rank < Threshold) #%>% select(Sp) %>% unlist

  } else if(!is.null(Hosts)){

    Hosts <- Hosts %>% str_replace_all(" ", "_")

    if(length(Hosts)<2) warning("This is only one host! Predictions might be rubbish.")

    if(length(intersect(Hosts, rownames(AllSims[[1]])))==0){
      stop("None of your hosts are in our dataset! :( are they marine, or are there synonyms?")

    } else if(length(intersect(Hosts, rownames(AllSims[[1]])))<length(Hosts)){

      warning(paste("Some hosts not found:",
                    paste(setdiff(Hosts, rownames(AllSims[[1]]), collapse = ", ")),
                    "; tread carefully"))
    }

    if(length(pHosts2)>0){

      FocalNet <- AllSums[pHosts2,]

      ValidEst <- list()

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

    } else print("Hosts Not Found!")

    GAMValid <- ValidEst %>%
      bind_rows() %>%
      group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>%
      slice(order(Count, decreasing = T))

    Df <- GAMValid
    Df$Rank = nrow(Df) - rank(Df$Count)
    PredHosts <- Df %>% filter(Rank < Threshold) #%>% select(Sp) %>% unlist

  }

  PredHostPolygons <- FullPolygons %>% filter(Host%in%PredHosts$Sp) %>%
    left_join(PredHosts, by = c("Host" = "Sp")) %>%
    mutate(Host = factor(Host, levels = Df[order(Df$Rank, decreasing = TRUE),"Sp"] %>% unlist))

  VirusName <- str_replace_all(Virus, "_", " ")

  Groups <- ifelse(STFull$tip.label%in%Hosts,"Known",ifelse(STFull$tip.label%in%PredHosts,"Predicted",""))

  groupInfo <- split(STFull$tip.label, Groups)
  chiroptera <- groupOTU(STFull, groupInfo)

  plot2 <- ggtree(chiroptera, aes(color = group, alpha = group)) +
    scale_colour_manual(values = c("black", "red", "blue")) +
    scale_alpha_manual(values = c(0.01,0.5,1)) +
    theme(legend.position = "top")

  g2 = ggplotGrob(plot2)

  xmin = -1.9*(10^7)
  xmax = -1.2*(10^7)
  ymin = -9*(10^6)
  ymax = 9*(10^6)

  rectborder <- 60000
  rect = data.frame(long = c(xmin-rectborder, xmin-rectborder, xmax+rectborder, xmax+rectborder),
                    lat = c(ymin-rectborder, ymax+rectborder, ymax+rectborder, ymin-rectborder))

  if(Facet == FALSE){

    ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) +
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) +
      labs(
        title = paste("Predicted", VirusName, "Hosts")) +
      coord_fixed() +
      scale_x_continuous(breaks = -10:10*2000000) +
      scale_y_continuous(breaks = -5:5*2000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat),
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) %>% return #Facet_wrap(~Host)

  } else {

    ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) +
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) + # alpha = max(Rank)-Rank)) +
      labs(#alpha = "Inverse Rank",
        title = paste("Predicted", VirusName, "Hosts")) +
      coord_fixed() +
      #theme(legend.position = "none") +
      scale_x_continuous(breaks = -10:10*2000000) +
      scale_y_continuous(breaks = -5:5*2000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat),
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) +
      facet_wrap(~Focal, ncol = 1,
                 labeller = labeller(Focal = c("0" = "Predicted", "1" ="Known"))) %>% return #Facet_wrap(~Host)
  }

}
