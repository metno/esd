#' Function to plot contours on stereographic projection of a lon/lat grid 
#'
#' @param field A field object or a matrix with longitude and latitude attributes which is the output of esd::map.field
#' @param lat A numeric vector of latitudes
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param nx Length of output grid along x-direction. Default: nx=100
#' @param ny Length of output grid along y-direction. Default: ny=100
#' @param nlevels Number of levels in contour plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param projection map projection "lonlat" (default), "sphere", "np", "sp" 
#' @param breaks Contour levels. If the argument 'breaks' is not defined, breaks = pretty(zlim, n=nlevels)
#' @param zlim Range of contours. If the arguments 'breaks' and 'zlim' are not defined, zlim = range(field, na.rm=TRUE) 
#' @param nlevels Number of contour levels. If the arguments 'breaks' and 'nlevels' are not defined, nlevels=5
#' @param main main title
#' @param add a boolean; If TRUE, add contour lines to an existing plot. If FALSE, create new plot with contour lines
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export spherical_contour
spherical_contour <- function(field, FUN="mean", lonR = NULL, latR = NULL, 
                              nx=100, ny=100, projection="sphere",
                              breaks=NULL, zlim=NULL, nlevels=5,
                              xlab=NULL, ylab=NULL, main=NULL, add=FALSE, 
                              col="grey70", drawlabels=TRUE, verbose = FALSE) {
  if(verbose) print("contour_projection")
  # Create two-dimensional map
  if(inherits(field, "field")) field <- map(field, FUN=FUN, plot=FALSE)
  # Extract longitude and latitude from the field object
  lonxy <- rep(lon(field), length(lat(field)))
  latxy <- sort(rep(lat(field), length(lon(field))))
  
  # Transform longitude/latitude coordinates to spherical
  transformed_coords <- cartesian2sphere(lonxy, latxy, lonR=lonR, latR=latR, 
                                         verbose=verbose)
  
  if (!requireNamespace("MBA", quietly = TRUE)) {
    stop("Package 'MBA' needed to use spherical_contour. Please install it.")
  } else {
    
    # Prepare data for MBA interpolation
    data <- data.frame(X = transformed_coords$X[transformed_coords$visible], 
                       Y = transformed_coords$Y[transformed_coords$visible], 
                       field = as.vector(field)[transformed_coords$visible])
    # Interpolation using MBA package
    mba_out <- MBA::mba.surf(data, no.X = nx, no.Y = ny)
    # Plot the contours on the transformed spherical coordinates
    Z_interp <- matrix(mba_out$xyz.est$z, nrow = nx, ncol = ny)
    if(is.null(breaks)) {
      if(is.null(zlim)) zlim <- range(field, na.rm=TRUE)
      if(is.null(nlevels)) nlevels <- 5
      breaks <- pretty(zlim, n=nlevels)
    }
    contour(mba_out$xyz.est$x, mba_out$xyz.est$y, Z_interp, levels=breaks, 
            drawlabels = drawlabels, xlab = xlab, ylab = ylab, main = main, 
            col=col, add = add)
  }
}