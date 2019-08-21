# Buffer 2 ####

buffer2 <- function(r, dist) {

  if("raster"%in%class(r)|"RasterLayer"%in%class(r)){
    projR <- projectRaster(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- raster::buffer(projR, dist)
    projRb <- projectRaster(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
    projRb[!is.na(projRb)] <- 1
  }

  if("sf"%in%class(r)){
    projR <- st_transform(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- st_buffer(projR, dist)
    projRb <- st_transform(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
  }

  return(projRb)

}
