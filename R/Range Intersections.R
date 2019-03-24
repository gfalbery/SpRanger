
IntersectGet <- function(Rasterstack, Names = "All", Predicate = F){

  if(Names == "All") Names = names(Rasterstack)

  N = length(Names)

  if(!Predicate){

    lapply(1:N, function(a){

      print(Names[a])

      lapply(1:N, function(b){

        if(a<b){

          raster::intersect(MammalRanges[[Names[a]]], MammalRanges[[Names[b]]])

        }

      })
    })
  } else {

    if(length(setdiff(Names, rownames(Predicate))>0)){
      Names <- intersect(Names, colnames(Predicate))
      print("Warning: Some species not in predicate matrix!!")
    }

    lapply(1:N, function(a){

      print(Names[a])

      lapply(1:N, function(b){

        if(a<b){

          if(Predicate[Names[a],Names[b]]>0){

            raster::intersect(MammalRanges[[Names[a]]], MammalRanges[[Names[b]]])

          }
        }
      })
    })
  }
}
