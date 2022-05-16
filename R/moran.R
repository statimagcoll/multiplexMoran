#' Compute Moran's I matrix
#'
#' This code was modified from the original in the spdep package.
#'
#' @param x A matrix of variables to compute the bivariate Moran's I
#' @param listw A listw object specifying the weighting matrix.
#' @param comp A vector of integers, strings, or factors indicating broad tissue compartments in which to compute Moran's I.
#' @param zero.policy How to deal with zones without neighbors. If TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param NAOK If TRUE, NAs are removed in computations.
#' @param demean If TRUE demean columns of x. Might not want to if x is dummy variables.
#'
#' @return Returns a list of matrices, where each matrix is the multivariate Moran's I for each tissue compartment.
#' @importFrom spdep lag.listw
#' @export
moran = function(x, listw, comp=NULL, zero.policy = FALSE, NAOK = FALSE, demean=TRUE){
  x = as.matrix(x)
  C = nrow(x)
  M = ncol(x)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  if (n1 != C)
    stop("objects of different length")
  if(demean){
    xx <- colMeans(x, na.rm = NAOK)
    z <- sweep(x, 2, xx, '-')
  } else {
    z <- x
  }
  zz <- sqrt(colMeans(z^2, na.rm = NAOK))
  #K <- (length(x) * sum(z^4, na.rm = NAOK))/(zz^2)
  S0 = sum(lag.listw(listw, rep(1, C), zero.policy = zero.policy, NAOK = NAOK))
  # apply to all columns of z
  lz <- apply(z, 2, function(w) lag.listw(listw, w, zero.policy = zero.policy, NAOK = NAOK) )
  I <- crossprod(z, lz)/outer(zz, zz)/S0
  res <- I
  res
}

#' Compute Moran's Local Index of Spatial Association (LISA) matrix
#'
#' This code was modified from the moran function above. LISA is directional and is the association centered around marker y.
#' The column means of LISA are equal to Moran's I.
#'
#' @param y A vector to compute the index relative to.
#' @param x A matrix of variables to compute the lisa with y.
#' @param listw A listw object specifying the weighting matrix.
#' @param comp A vector of integers, strings, or factors indicating broad tissue compartments in which to compute Moran's I.
#' @param zero.policy How to deal with zones without neighbors. If TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param NAOK If TRUE, NAs are removed in computations.
#' @param demean If TRUE demean columns of x. Might not want to if x is dummy variables.
#'
#' @return Returns a list of matrices, where each matrix is the multivariate Moran's I for each tissue compartment.
#' @importFrom spdep lag.listw
#' @export
lisa = function(y, x, listw, comp=NULL, zero.policy = FALSE, NAOK = FALSE, demean=TRUE){
  x = as.matrix(cbind(y, x) )
  C = nrow(x)
  M = ncol(x)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  if (n1 != C)
    stop("objects of different length")
  if(demean){
    xx <- colMeans(x, na.rm = NAOK)
    z <- sweep(x, 2, xx, '-')
  } else {
    z <- x
  }
  zz <- sqrt(colMeans(z^2, na.rm = NAOK))
  #K <- (length(x) * sum(z^4, na.rm = NAOK))/(zz^2)
  S0 = sum(lag.listw(listw, rep(1, C), zero.policy = zero.policy, NAOK = NAOK))
  # apply to all columns of z
  lz <- apply(z, 2, function(w) lag.listw(listw, w, zero.policy = zero.policy, NAOK = NAOK) )
  # multiply each column of lz with z, scale out the variance, scale out the weights and sample, scale backup the sample
  LISA <- sweep(sweep(lz, 1, z, '*'), 2, zz[1] * zz, '/' )/S0*C
  LISA
}


#' Convert coordinates to a listw object
#'
#' Parameters for this should be based on biological knowledge. The algorithm does
#' a k nearest neighbor search and then selects all cells that are less than maxDist
#' distance from the center as neighbors with nonzero weights. If k is chosen too
#' small then it's possible for the weight matrix to be nonsymmetric because the k
#' nearest neighbors are not symmetric. Choosing k large seems to be pretty quick still.
#'
#' @param x A matrix of coordinates (of cell centroids)
#' @param k The maximum number of neighbors each cell can have
#' @param maxDist The maximum distance a neighbor can be from the central cell
#' @param comp A vector of integers, strings, or factors, indicating what tissue compartment each cell belongs to
#' @param W The weights used for the surrounding cells. defaults to 1
#'
#' @return Returns a list of matrices, where each matrix is the multivariate Moran's I for each tissue compartment.
#' @importFrom Matrix sparseMatrix
#' @importFrom FNN get.knn
#' @importFrom spdep mat2listw
#' @export
coords2listw = function(x, k=200, maxDist=50, comp=NULL, W=1){
  nd = get.knn(x, k=k)
  # if(is.null(W)){
  #   W = 1/c(nd$nn.dist)
  # }
  if(!is.null(maxDist)){
    W = c(as.numeric(nd$nn.dist< maxDist))
  }
  listw = spdep::mat2listw(sparseMatrix(i = rep(1:nrow(x), k), j = c(nd$nn.index), x=W))
}

