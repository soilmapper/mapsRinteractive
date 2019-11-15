#' @title ordkrige
#'
#' @description Regression kriging using a standardized variogram.
#'
#' @param check_data Locigal. Shall attributes, geometries and projections of
#' the input data (arguments x, y and z) be checked. Default = TRUE.
#'
#' @inheritParams mri
#'
#' @return A list with 1) a raster layer with predicted values and 2) a
#' SpatialPolygonsDataFrame with cross-validation data. For details, see
#' mri function.
#'
#' @details This is the ordinary kriging function called by the mri function.
#' It uses a standardized variogram and requires a raster template for which
#' predictions are made. For details, see documentation of the mri function.
#'
#' @export
ordkrige<-function(x=NULL, y=NULL, z=NULL, field = NULL, edge = 0, filter = 1, resolution = NULL,
                   md = 'Sph', rg = NULL, ng = 0.1, check_data=T){

  #check input data
  if(check_data){
    a<-check(x=x, y=y, z=z, field=field, edge = edge, filter=filter, resolution=resolution)
    x<-a[[1]]; y<-a[[2]]; z<-a[[3]]; feedback<-a[[4]]
  }

  #compute range (argument rg) if not specified by user
  if(is.null(rg)) rg<-0.5*sqrt(area(y))

  ##parameterize standardized semivariogram model
  sill<-var(z$obs, na.rm=T)
  mod<- vgm(psill = (1-ng)*sill, model= md, range = rg,nugget= ng*sill)

  #cross validate
  for (i in 1:nrow(z)){
    obs_pred<-krige(obs~1, locations=z[-i,], newdata=z[i,], model = mod, debug.level=0)
    z[i,'ordkrig_cv']<-obs_pred@data['var1.pred']
  }

  #ordinary kriging to raster
  crs(x)<-crs(z) #to fix a bug
  gsmod<-gstat(formula=obs~1, locations=z, model=mod)
  ordkrig<-interpolate(x, gsmod,  fun=predict, debug.level=0)
  ordkrig<-mask(ordkrig, x)

  #return objects
  return(list(ordkrig, z))
}

