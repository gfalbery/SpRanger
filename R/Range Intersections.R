
IntersectGet <- function(Rasterstack, Names = "All", Predicate = F){

  if(Names == "All") Names = names(Rasterstack)

  N = length(Names)

  if(Predicate){

    lapply(1:N, function(a){

      print(Names[a])

      lapply(1:N, function(b){

        if(a<b){

          if(RangeAdj[Names[a],Names[b]]>0){

            raster::intersect(MammalRanges[[Names[a]]], MammalRanges[[Names[b]]])

          }
        }
      })
    })
  } else {

    lapply(1:N, function(a){

      print(Names[a])

      lapply(1:N, function(b){

        if(a<b){

          raster::intersect(MammalRanges[[Names[a]]], MammalRanges[[Names[b]]])

        }
      })
    })
  }
}
