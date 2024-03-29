#' Function to create stack for integration points
#'
#' @param mesh INLA mesh (e.g., \code{$mesh} of object generated by \code{\link{MakeSpatialRegion}}).
#' @param data \code{data.frame} with at least columns "X" and "Y" (or \code{coordnames}) giving locations, plus optional other covariate. Alternatively, a \code{SpatialPointsDataFrame} object.
#' @param area Area around each integration point (e.g., the \code{$w} element of the object generated by \code{\link{MakeSpatialRegion}}).
#' @param tag name for the stack.
#' @param coordnames Names for coordinates in data (if \code{data} is not a \code{SpatialPointsDataFrame}).
#' @param InclCoords should coordinates be included in data? (defaults to \code{FALSE}).
#' @return An INLA stack for integration
#'
#' @export
#' @import INLA
MakeIntegrationStackPts <- function(mesh, data, area, tag='mesh',
                                 coordnames=c("X","Y"), InclCoords=FALSE) {
  
  #check for coordinates in covariate data
  if(class(data)!="SpatialPointsDataFrame" & !all(coordnames%in%names(data))) stop("Coordinates not in the data")

  #get coordindate and covariate column names
  if(class(data)=="SpatialPointsDataFrame") {
    coordnames <- colnames(data@coords)
    Names <- colnames(data@data)
  } else {
    Names  <- names(data)[!names(data)%in%coordnames]
  }
  
  #if no coordinate names given default to x and y and join covariate names
  if(InclCoords) Names  <- c(coordnames, Names)

  #get mesh triangle centroids
  Points <- cbind(c(mesh$loc[,1]), c(mesh$loc[,2]))
  
  #set column names of centroids to cordinate names from covariate data
  colnames(Points) <- coordnames
  
  #get value of nearest covariate to each mesh centroids
  NearestCovs <- GetNearestCovariate(points=Points, covs=data)
  
  #set intercept to 1 for each mesh element
  NearestCovs$Intercept <- rep(1,nrow(NearestCovs))
  
  #add coordinate to data part of spatialpoints object if they are given
  if(InclCoords) {
    NearestCovs@data[,colnames(NearestCovs@coords)] <- NearestCovs@coords
  }

  # Projector matrix for integration points.
  projmat.ip <- Matrix::Diagonal(mesh$n, rep(1, mesh$n))  # from mesh to integration points

  #create inla stack
  stk.ip <- inla.stack(data=list(resp=rep(0,mesh$n), e=area),
                       A=list(1,projmat.ip), tag=tag,
                       effects=list(NearestCovs@data, list(i=1:mesh$n)))
  
   return(stk.ip)
}
