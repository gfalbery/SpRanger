# Host Range Plot ####

RangePlot <- function(Virus = NULL,
                     HostList = NULL,
                     Map = TRUE,
                     Tree = TRUE,
                     Legend = "right"){

  require(dplyr); require(stringr); library(tidyverse); require(ggtree); require(cowplot)

  if(is.null(Virus)&is.null(HostList)) stop("One of Virus or HostList must be specified!")
  if(!is.null(Virus)&!is.null(HostList)) stop("Only specify one of Virus or Host please :)")

  if(!is.null(Virus)){

    Virus <- Virus %>% str_replace_all(" ", "_")

    if(!Virus%in%names(VirusAssocs)) stop("Virus not found :( try entering hosts manually?") else print("We have this virus! Wahey.")

    HostList <- VirusAssocs[[Virus]]

  }

  if(any(!HostList%in%FullPolygons$Host)) print(paste0("Some hosts not in our data:", paste(setdiff(HostList, FullPolygons$Host), colllapse = ", ")))

  RangePolygons <- FullPolygons %>% filter(Host%in%HostList)

  VirusName <- str_replace_all(Virus, "_", " ")

  MapPlot <- ggplot(RangePolygons, aes(long, lat, group = paste(Host, group))) +
    geom_polygon(fill = NA, aes(colour = Host)) +
    coord_fixed() +
    theme(legend.position = Legend) +
    scale_x_continuous(breaks = -10:10*5000000) +
    scale_y_continuous(breaks = -5:5*5000000)+
    theme(plot.title = element_text(hjust=0.5))

  if(Map){

    MapPlot <- ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) +
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) +
      labs(
        title = paste("Predicted", VirusName, "Hosts")) +
      coord_fixed() +
      theme(legend.position = Legend) +
      scale_x_continuous(breaks = -10:10*5000000) +
      scale_y_continuous(breaks = -5:5*5000000)
  }

  if(!is.null(Virus)) MapPlot <- MapPlot + ggtitle(paste(VirusName, "Hosts"))

  if(Tree){

    Groups <- ifelse(STFull$tip.label%in%HostList,"Known","")

    groupInfo <- split(STFull$tip.label, Groups)
    chiroptera <- groupOTU(STFull, groupInfo)

    plot2 <- ggtree(chiroptera, aes(color = group, alpha = group)) +
      scale_colour_manual(values = c("black", "red")) +
      scale_alpha_manual(values = c(0.05,1))

    g2 <- ggplotGrob(plot2)

    xmin <- -2.1*(10^7)
    xmax <- -1.45*(10^7)
    ymin <- -9*(10^6)
    ymax <- 9*(10^6)

    rectborder <- 60000
    rect <- data.frame(long = c(xmin-rectborder, xmin-rectborder, xmax+rectborder, xmax+rectborder),
                       lat = c(ymin-rectborder, ymax+rectborder, ymax+rectborder, ymin-rectborder))

    MapPlot <- MapPlot +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat),
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax)
  }

  return(MapPlot)

}
