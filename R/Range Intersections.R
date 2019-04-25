
IntersectGet <- function(Rasterstack, Names = "All", Predicate = F){

  require(Matrix)

  if(all(Names == "All")) Names = names(Rasterstack) else{
    if(length(setdiff(Names, names(Rasterstack)))>0){
      Names = intersect(Names, names(Rasterstack))
      print("Warning: some names not in the rasterstack!")
    }
  }

  N = length(Names)

  if(length(Predicate)==0){

    RasterList <- lapply(1:N, function(a){

      print(Names[a])

      lapply(1:N, function(b){

        if(a<b){

          raster::intersect(Rasterstack[[Names[a]]], Rasterstack[[Names[b]]])

        }

      })
    })
  } else {

    if(length(setdiff(Names, rownames(Predicate))>0)){
      Names <- intersect(Names, colnames(Predicate))
      print("Warning: Some species not in predicate matrix!!")
    }

    if(length(setdiff(rownames(Predicate), Names)>0)){
      Names <- intersect(Names, colnames(Predicate))
      print("Warning: Some predicate species not in rasterstack!!")
    }

    Predicate <- Predicate[Names, Names]
    Predicate[upper.tri(Predicate)] <- NA

    Predicate2 <- Predicate %>% reshape2::melt() %>% rename(Sp = Var1, Sp2 = Var2) %>%
      mutate(Sp = as.character(Sp), Sp2 = as.character(Sp2)) %>%
      na.omit() %>% filter(!as.numeric(value)==0)

    RasterList <- lapply(1:nrow(Predicate2), function(a){

      print(paste0(a, ": ", paste(Predicate2[a,c("Sp","Sp2")], collapse = ".")))

      raster::intersect(Rasterstack[[Predicate2[a,"Sp"]]], Rasterstack[[Predicate2[a,"Sp2"]]])

    })

    names(RasterList) <- paste(Predicate2$Sp, Predicate2$Sp2, sep = ".")

  }

}
