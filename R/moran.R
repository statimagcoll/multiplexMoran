#' Compute Moran's I matrix
#'
#' This code was modified from the original in the spdep package.
#'
#' @param x A matrix of variables to compute the bivariate Moran's I
#' @param listw A listw object specifying the weighting matrix.
#' @param comp A vector of integers, strings, or factors indicating broad tissue compartments in which to compute Moran's I.
#' @param zero.policy How to deal with zones without neighbors. If TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param NAOK If TRUE, NAs are removed in computations.
#'
#' @return Returns a list of matrices, where each matrix is the multivariate Moran's I for each tissue compartment.
#' @importFrom spdep lag.listw
#' @export
moran = function(x, listw, comp=NULL, zero.policy = FALSE, NAOK = FALSE){
  x = as.matrix(x)
  C = nrow(x)
  M = ncol(x)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  if (n1 != C)
    stop("objects of different length")
  xx <- colMeans(x, na.rm = NAOK)
  z <- sweep(x, 2, xx, '-')
  zz <- sqrt(colMeans(z^2, na.rm = NAOK))
  #K <- (length(x) * sum(z^4, na.rm = NAOK))/(zz^2)
  S0 = sum(lag.listw(listw, rep(1, C), zero.policy = zero.policy, NAOK = NAOK))
  # apply to all columns of z
  lz <- apply(z, 2, function(w) lag.listw(listw, w, zero.policy = zero.policy, NAOK = NAOK) )
  I <- crossprod(z, lz)/outer(zz, zz)/S0
  res <- I
  res
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

