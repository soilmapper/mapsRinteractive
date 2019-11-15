#' @title check
#'
#' @description Checks attributes, geometries and projections of spatial datasets.
#'
#' @inheritParams mri
#'
#' @return A list with checked and corrected datasets together with a vector of logged feedback.
#'
#' @details Intended for checking data in functions of mapsRinteractive.
#'
#' @keywords Internal
#'
check<-function(x = NULL, y = NULL, z = NULL, field = NULL, edge = 0, filter = 1,
                resolution = NULL){

  #are required spatial data provided?
  if(is.null(x)) {stop(call. =F, 'Raster data (argument: x) is missing.')}
  if(is.null(z)) {stop(call. =F, 'Point location data (argument: z) is missing.')}
  if(is.null(field)) {stop(call. =F, 'Attribute column name (argument: field) is missing.')}
  if(!(field %in% names(z))) {stop(call. =F, 'Argument field is not a column name in z.')}

  #are other arguemnts appropriate

  if(!is.null(filter)) if(is.na(filter) | filter<1 | 20<filter){stop(call. =F, 'Argument filter is out of range.')}
  if(is.null(filter)) {stop(call. =F, 'Argument filter is out of range.')}
  if(!is.null(resolution)) if(is.na(resolution) | resolution <0.1 | 10000<resolution){stop(call. =F, 'Argument resolution is out of range.')}
  if(!is.null(edge)) if(is.na(edge) | edge<0|10000<edge){stop(call. =F, 'Argument edge is out of range.')}

  #create boundary polygon (y), if not provided
  if(is.null(y)) {
    y<-spdf_from_extent(x)
    t<-'No SpatialPolygonsDataFrame (argument: y) provided. The analyses will be performed within the intersection of the raster and a buffer zone around the soil observations. See details.'
    feedback<-fback(t)
  }

  #are data projected?
  if(is.na(crs(x))){stop(call.=F, 'Raster (argument: x) not projected.')}
  if(is.na(crs(y))){stop(call.=F, 'SpatialPolygonsDataFrame (argument: Y) not projected.')}
  if(is.na(crs(z))){stop(call.=F, 'SpatialPointsDataFrame (argument: Z) not projected.')}

  #are data projected onto the same coordinate system?
  if(compareCRS(x, y)==F){
    t<-'The coordinate systems of the raster data and the polygon data are not the same. The polygon data will be projected to the coordinate system of the raster.'
    feedback<-fback(t)
    y<-spTransform(x=y, CRSobj=crs(x))
  }
  if(compareCRS(x, z)==F){
    t<-'The coordinate systems of the raster data and the point data are not the same. The point data will be projected to the coordinate system of the raster.'
    feedback<-fback(t)
    z<-spTransform(x=z, CRSobj=crs(x))
  }

  #is the raster data in a polar or cartesian coordinate system? Reproject to web mercator if the latter
  prj<-crs(x) #store projection
  web_mercator<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs'
  if (is.projected(crs(x))==F) x<-spTransform(x=x, CRSobj=CRS(web_mercator))
  if (is.projected(crs(y))==F) y<-spTransform(x=y, CRSobj=CRS(web_mercator))
  if (is.projected(crs(z))==F) x<-spTransform(x=z, CRSobj=CRS(web_mercator))

  #If the raster dataset has more then one layer, use only the first one.
  if(nlayers(x)>1){
    x<-x[[1]]
    t<-'The raster dataset has more than one layer. Only the first layer will be used.'
    feedback<-fback(t)
  }

  #is the polygon overlapping the raster (check again after clipping of polygon)?
  if(!is.null(y)) if(length(y[spdf_from_extent(x),])==0) {
    stop(call. =F, 'the raster (argument x) and the boundary polygon (argument y) are not overlapping. Consider skipping argument y')
  }

  #clip the polygon, if the raster does not fill it
  y<-rgeos::gIntersection(y, spdf_from_extent(x), byid=F)

  #shrink polygon, if a buffer (argument edge) is specified
  y_large<-y #save a copy
  if(edge>0)if(!is.null(buffer(y,-edge))) {
    y <- buffer(y,-edge)} else {
     stop(call. =F, 'The buffer is too large, please chose a smaller value of argument: edge.')
    }

  #Prepare and check point location data
  ##rename attribute column
  names(z)[names(z)==field]<-'obs'
  ##drop other columns
  z<-z[,'obs']
  ##are data numeric
  if(!is.numeric(z$obs)){stop(call. =F, 'Point observation data (specified by arguments:  z and field) not numeric.')}
  ##are there enough point observation data left after NA values and infinite values have been excluded
  z<-z[,'obs'] #drop other columns
  z<-z[!is.na(z$obs),] #omit points with NA values
  z<-z[is.finite(z@data$obs),] #omit points with infinte values
  if(nrow(z)<10) {stop(call. =F, 'There are fewer than ten non-NA values in the point observation dataset.')}

  #clip the polygon, if the buffered points does not fill it
  d<-as.matrix(stats::dist(coordinates(z)))#distance matrix
  f<-function(d){return(sort(unique(d))[2])}
  d<-apply (d, 1, f) #rowmins of distance matrix
  width<-1.5*kth(d, 2) #largest minimum distance to
  pt_buff<-buffer(z, width, dissolve=T)
  if(is.null(rgeos::gIntersection(y, pt_buff, byid=F))){
    stop(call. =F, 'The prepared datasets with polygon and point data do not overlap.')
  } else{
    y<-rgeos::gIntersection(y, pt_buff, byid=F)
    }

  #is the polygon overlapping the raster (check again after clipping of polygon)?
  if(!is.null(y)) if(length(y[spdf_from_extent(x),])==0) {
    stop(call. =F, 'the raster and the boundary polygon are not overlapping. Consider skipping argument y')
  }

  #make sure y is a SpatialPolygonsDataFrame and not just of class SpatialPolygons
  y$id<-1

  #mean filter raster data
  if(filter>1){
    aa<-is.na(x)
    if(even(filter)) filter<-filter+1 #if filter is an even number, add one cell
    w<-matrix(1, filter, filter)
    x<-raster::focal(x, w=w, fun=mean, na.rm=T, NAonly = F, pad=T)
    x[aa]<-NA
    t<-'Your raster have been mean filtered. See argument filter.'
    feedback<-fback(t)
  }

  #crop and mask the raster by the polygon
  x<-crop(x, y)
  x<-mask(x, y)

  #is the raster resolution too coarse in relation to area of analysis.
  resmax<-signif((area(y)/10),1)
  if(!is.null(resolution)) if(resolution>resmax){
    stop(call. =F, paste0('The specified raster resolution (argument: resolution) is too coarse in relation to the delineated area of analyses. In this case a resolution smaller than ', resmax, ' m is allowed.'))
    }
  if(is.null(resolution)) if(res(x)[1]>resmax){
    stop(call. =F, paste0('The raster resolution (argument: resolution) is too coarse in relation to the delineated area of analyses. In this case a resolution samller than ', resmax, ' m is required.'))
  }

  #is the raster resolution too fine in relation to area of analysis.
  resmin<-signif((area(y)/(10000*10000)),1)
  if(!is.null(resolution)) if(resolution<resmin){
    stop(call. =F, paste0('The specified raster resolution (argument: resolution) is too fine in relation to the delineated area of analyses. In this case a resolution larger than ', resmin, ' m is allowed.'))
  }
  if(is.null(resolution)) if(res(x)[1]<resmin){
    t<-'The raster resolution is very fine in relation to the area. The raster will be resampled. Consider specifying a coarser resolution (argument: resolution).'
    feedback<-fback(t)
    template<-raster(res=resmin, ext=extent(x), crs=prj)
    x<-resample (x=x, y=template, method='bilinear')
    t<-paste0('Your raster will be resampled to a cell size of ', resmin, ' m to speed up the processing. If you do not want this to happen. Please specify the desired resolution (argument: resolution).')
    feedback<-fback(t)
  }

  #rename raster field
  names(x)<-'map'

  #are some point location data outside y
  if(nrow(z)>nrow(z[y,])){
    z<-z[y,]
    t<-'Some point locations are outside the area specified for analysis. These will not be used. If you have specified an edge (buffer inside the polygon), consider a smaller value for this argument.'
    feedback<-fback(t)
  }

  #extract raster data for point locations and omit points where map data is NA
  z$map<-extract(x=x, y=z, df=F)
  z<-z[complete.cases(z@data),]

  #are there enough point data, wihin the intersection of x and y?
  if (nrow(z)<10) {
  t<-'There are fewer than ten point locations with both raster and observation data. The performance of the map adaptation/evaluation may not be very good.'
  feedback<-fback(t)
  }

  #check wheter there is any variation in the point location or and or the raster data
  obs<-z$obs[!is.na(z$obs)]
  if(var(obs)==0) {stop(call. =F, 'After data preparation, the variance in point observation data is zero.')}
  map<-z$map[!is.na(z$map)]
  if(var(map)==0) {stop(call. =F, 'After data preparation, the variance in raster data is zero.')}

  #check number of point locations with fulfilled data requirements
  t<-paste0('In total ', nrow(z), ' point locations fulfilled the data and geometry requirements.')
  if(nrow(z)==0) stop(call. =F)
  feedback<-fback(t)

  #final cropping and masking of the raster by the polygon
  x<-crop(x, y)
  x<-mask(x, y)

  #return results
  return(list(x=x, y=y, z=z, feedback=feedback))
}
