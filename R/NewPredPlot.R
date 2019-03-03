

PredPlot <- function(Virus = NULL,
                     Hosts = NULL,
                     Threshold = 10,
                     Validate = TRUE,
                     Focal = c(1,0),
                     Facet = FALSE,
                     Legend = "right"){

  require(dplyr); require(stringr); library(tidyverse); require(ggtree)

  if(!is.null(Virus)&!is.null(Hosts)) stop("Only specify one of Virus or Host please :)")

  if(!is.null(Virus)){

    Virus <- Virus %>% str_replace_all(" ", "_")

    if(!Virus%in%names(GAMValid)) stop("Virus not found :( try entering hosts manually?") else print("We have this virus! Wahey.")

    Df <- GAMValid[[Virus]]

    Df$Rank <- nrow(Df) - rank(Df$Count)

    PredHosts <- Df %>% filter(Rank < Threshold) #%>% select(Sp) %>% unlist

  } else if(!is.null(Hosts)){

    Hosts <- Hosts %>% str_replace_all(" ", "_")

    if(length(Hosts)<2) print("This is only one host! Predictions might be rubbish.", immediate = T)

    if(length(intersect(Hosts, rownames(AllSums)))==0){

      stop("None of your hosts are in our dataset! :( are they marine, or are there synonyms?")

    } else if(length(intersect(Hosts, rownames(AllSums)))<length(Hosts)){

      print(paste("Warning: Some hosts not found:",
                  paste(setdiff(Hosts, rownames(AllSums)), collapse = ", "),
                  "; tread carefully"))
    }

    pHosts2 <- intersect(Hosts, rownames(AllSums))

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

        GAMValid <- ValidEst %>%
          bind_rows() %>%
          group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>%
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = as.factor(Focal))

      } else {

        #pHosts3 <- setdiff(colnames(FocalNet), pHosts2)
        #Estimates <- FocalNet[pHosts2, pHosts3]/length(AllSims)

        Estimates <- FocalNet[pHosts2,]/length(AllSims)

        if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)

        Ests <- data.frame(Sp = names(sort(colSums(Estimates), decreasing = T)),
                           Count = sort(colSums(Estimates), decreasing = T)/nrow(Estimates)) %>%
          mutate(Focal = ifelse(Sp%in%Hosts,1,0))

        rownames(Ests) <- Ests$Sp

        GAMValid <- Ests %>%
          slice(order(Focal)) %>%
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = as.factor(Focal))

      }

    } else print("Hosts Not Found!")

    Df <- GAMValid %>% as.data.frame()
    Df$Rank <- nrow(Df) - rank(Df$Count)
    Df <- Df %>% mutate(Include = ifelse(Rank<Threshold&Focal==0,1,0))
    if(1 %in% Focal) Df[Df$Focal==1,"Include"] <- 1
    PredHosts <- Df %>% filter(Include == 1) #%>% select(Sp) %>% unlist

  }

  PredHostPolygons <- FullPolygons %>% filter(Host%in%PredHosts$Sp) %>%
    left_join(PredHosts, by = c("Host" = "Sp")) %>%
    mutate(Host = factor(Host, levels = Df[order(Df$Rank, decreasing = TRUE),"Sp"] %>% unlist))

  VirusName <- str_replace_all(Virus, "_", " ")

  Groups <- ifelse(STFull$tip.label%in%Hosts,"Known",ifelse(STFull$tip.label%in%PredHosts$Sp,"Predicted",""))

  groupInfo <- split(STFull$tip.label, Groups)
  chiroptera <- groupOTU(STFull, groupInfo)

  plot2 <- ggtree(chiroptera, aes(color = group, alpha = group)) +
    scale_colour_manual(values = c("black", "red", "blue")) +
    scale_alpha_manual(values = c(0.01,0.5,1))

  g2 <- ggplotGrob(plot2)

  xmin <- -2.1*(10^7)
  xmax <- -1.45*(10^7)
  ymin <- -9*(10^6)
  ymax <- 9*(10^6)

  rectborder <- 60000
  rect <- data.frame(long = c(xmin-rectborder, xmin-rectborder, xmax+rectborder, xmax+rectborder),
                     lat = c(ymin-rectborder, ymax+rectborder, ymax+rectborder, ymin-rectborder))

  if(Facet == FALSE){

    MapPlot <- ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) +
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) +
      labs(
        title = paste("Predicted", VirusName, "Hosts")) +
      coord_fixed() +
      theme(legend.position = Legend) +
      scale_x_continuous(breaks = -10:10*5000000) +
      scale_y_continuous(breaks = -5:5*5000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat),
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) %>% return #Facet_wrap(~Host)

  } else {

    MapPlot <- ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) +
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) + # alpha = max(Rank)-Rank)) +
      labs(#alpha = "Inverse Rank",
        title = paste0("Predicted", VirusName, " Hosts")) +
      coord_fixed() +
      theme(legend.position = Legend) +
      scale_x_continuous(breaks = -10:10*5000000) +
      scale_y_continuous(breaks = -5:5*5000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat),
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) +
      facet_wrap(~Focal, ncol = 1,
                 labeller = labeller(Focal = c("0" = "Predicted", "1" ="Known"))) %>% return #Facet_wrap(~Host)

  }

  if(Validate){

    ValidPlot <- ggplot(Df, aes(Focal, Count, colour = Focal)) +
      ggforce::geom_sina() +
      geom_text(data = Df[1,], inherit.aes = F, aes(x = 1.5, y = max(Df$Count)*1.1, label = paste("Mean Focal Rank =", mean(Df[Df$Focal==1,"Rank"])))) +
      ggtitle("Model Success Rate")

    return(list(MapPlot,
                ValidPlot))
  }else{

    return(MapPlot)

  }

}
