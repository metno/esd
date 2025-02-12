#' Function to plot contours on spherical grid 
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
#' @param main main title
#' @param add a boolean; If TRUE, add contour lines to an existing plot. If FALSE, create new plot with contour lines. 
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export spherical_contour
spherical_contour <- function(field, FUN="mean", lonR = NULL, latR = NULL, 
                              nx=100, ny=100, nlevels=5, #method="akima",
                              xlab=NULL, ylab=NULL, main=NULL, add=FALSE, 
                              col="grey25", drawlabels=TRUE, verbose = FALSE) {
  if(verbose) print("spherical_contour")
  # Create two-dimensional map
  if(inherits(field, "field")) field <- map(field, FUN=FUN, plot=FALSE)
  # Extract longitude and latitude from the field object
  lonxy <- rep(lon(field), length(lat(field)))
  latxy <- sort(rep(lat(field), length(lon(field))))
  
  # Transform longitude/latitude coordinates to spherical
  transformed_coords <- esd::cartesian2sphere(lonxy, latxy, lonR=lonR, latR=latR, 
                                              verbose=verbose)
  ## KMP 2025-02-10: I tried 2 methods for interpolating to the spherical grid, 
  ##   but will remove one of them to reduce reliance on external libraries
  ##   I'm going with the method inclded in the MBA package, as akima was not as apt at dealing with very large/dense grids
  #if(method=="akima") {
  #  if (!requireNamespace("akima", quietly = TRUE)) {
  #    stop("Package 'akima' needed to use spherical_contour. Please install it or use type='MBA' to interpolate with the MBA package.")
  #  } else {
  #    ## Interpolated using akima package
  #    interpolation <- akima::interp(transformed_coords$X, transformed_coords$Z, 
  #                                   field, duplicate = "mean")
  #    contour(interpolation$x, interpolation$y, interpolation$z, nlevels = nlevels, 
  #            drawlabels = drawlabels, xlab = xlab, ylab = ylab, main = main, 
  #            col=col, add = add)
  #  }
  #} else if(method=="MBA") {
    if (!requireNamespace("MBA", quietly = TRUE)) {
      stop("Package 'MBA' needed to use spherical_contour. Please install it.")
    } else {
      # Prepare data for MBA interpolation
      data <- data.frame(X = transformed_coords$X, Y = transformed_coords$Z, 
                         field = as.vector(field))
      # Interpolation using MBA package
      mba_out <- MBA::mba.surf(data, no.X = nx, no.Y = ny)
      # Plot the contours on the transformed spherical coordinates
      Z_interp <- matrix(mba_out$xyz.est$z, nrow = nx, ncol = ny)
      contour(mba_out$xyz.est$x, mba_out$xyz.est$y, Z_interp, nlevels = nlevels, 
              drawlabels = drawlabels, xlab = xlab, ylab = ylab, main = main, 
              col=col, add = add)
    }
  #}
}