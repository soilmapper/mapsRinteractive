#' @title spdf_from_extent
#' @description Create a spatialPolygonsDataFrame from extent of a spatial object.
#'
#' @param x A spatial object.
#'
#' @details If x is projected, the spatialpolygonsdataframe will also be projected
#'
#' @examples
#' require(raster) #load required package.
#' r1<-raster::raster(ext=extent(c(0,10,0,10)), res=1, vals=1:100) #create example raster.
#' spdf<-spdf_from_extent(r1) #convert the raster extent to SpatialPolygonsdataFrame.
#' plot(spdf) #Plot results.
#'
#' @return SpatialPolygonsDataFrame.
#'
#' @import 'raster'
#' @import 'sp'
#'
#' @export
spdf_from_extent<-function(x){
  x1<-as.vector(extent(x))[1]
  x2<-as.vector(extent(x))[2]
  y1<-as.vector(extent(x))[3]
  y2<-as.vector(extent(x))[4]
  coords<-rbind(c(x1,y1), c(x1, y2), c(x2,y2), c(x2, y1), c(x1,y1))
  p = Polygon(coords)
  ps = Polygons(list(p),1)
  sps = SpatialPolygons(list(ps))
  if(!is.na(crs(x)))crs(sps)<-crs(x)
  sps$id<-1
  return(sps)
}
