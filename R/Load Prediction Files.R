
# Load files needed for host prediction ####

LoadPredFiles <- function(DataRoot){

  load(paste0(DataRoot,"/FullPolygons.Rdata"))
  load(paste0(DataRoot,"/GAMValidation.Rdata"))
  load(paste0(DataRoot,"/AllSums.Rdata"))
  load(paste0(DataRoot,"/WorldMap.Rdata"))


}
