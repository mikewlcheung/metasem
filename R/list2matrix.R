#' Convert a List of Symmetric Matrices into a Stacked Matrix
#' 
#' It converts a list of symmetric matrices into a stacked matrix. Dimensions
#' of the symmetric matrices have to be the same. It tries to preserve the
#' dimension names if possible. Dimension names will be created if there are no
#' dimension names in the first symmetric matrix.
#' 
#' 
#' @param x A list of \eqn{k}{k} \eqn{p}{p} x \eqn{p}{p} symmetric matrices.
#' @param diag Logical. If it is \code{TRUE}, \code{\link[OpenMx]{vech}} is
#' used to vectorize the (covariance) matrices. If it is \code{FALSE},
#' \code{\link[OpenMx]{vechs}} is used to vectorize the (correlation) matrices.
#' @return A \eqn{k}{k} x \eqn{p*}{p*} stacked matrix where \eqn{p* = p(p-1)/2
#' }{p* = p(p-1)/2} for \code{diag}=\code{FALSE} or \eqn{p* = p(p+1)/2 }{p* =
#' p(p+1)/2} for \code{diag}=\code{TRUE}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' C1 <- matrix(c(1,0.5,0.4,0.5,1,0.2,0.4,0.2,1), ncol=3)  
#' C2 <- matrix(c(1,0.4,NA,0.4,1,NA,NA,NA,NA), ncol=3)  
#' 
#' ## A list without dimension names 
#' list2matrix(list(C1, C2))
#' #      x2_x1 x3_x1 x3_x2
#' # [1,]   0.5   0.4   0.2
#' # [2,]   0.4    NA    NA
#' 
#' dimnames(C1) <- list( c("x","y","z"), c("x","y","z") )
#' dimnames(C2) <- list( c("x","y","z"), c("x","y","z") )
#' 
#' ## A list with dimension names
#' list2matrix(list(C1, C2))
#' #      y_x z_x z_y
#' # [1,] 0.5 0.4 0.2
#' # [2,] 0.4  NA  NA
#' 
list2matrix <- function(x, diag=FALSE) {
  if (length(x) == 0) {
    stop("\"x\" cannot be an empty list.")
  }

  if (!all(sapply(x, is.matrix))) {
    stop("All elements of \"x\" must be matrices.")
  }

  if (!is.list(x))
    stop("Input \"x\" must be a list of matrices.")

  ## if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
  if (length(unique(sapply(x, function(x) dim(x)[[1]]))) != 1)
    stop("Dimensions of matrices in \"x\" have to be the same in order to stack them together.")
  
  if (is.null(dimnames(x[[1]]))) {
    oldNames <- paste("x", 1:dim(x[[1]])[[1]], sep = "")
  } else {
    oldNames <- dimnames(x[[1]])[[1]]
  }
   
  if (diag) {
    psNames <- vech(outer(oldNames, oldNames, paste, sep = "_"))
    ## out <- t(sapply(x, function(x) {(vech(x))}))
    ## out is a list
    out <- lapply(x, vech)    
  } else {
    psNames <- vechs(outer(oldNames, oldNames, paste, sep = "_"))
    ## Fix a bug found by Steffen Zitzmann when x is a 2x2 matrix with diag=FALSE
    ## It returns 1xn vector rather than nx1 because of sapply().
    ## out <- t(sapply(x, function(x) {(vechs(x))}))
    out <- lapply(x, vechs)
  }
  
  ## convert the list into a matrix
  ## out <- matrix(unlist(out), nrow=length(out), byrow=TRUE)
  out <- do.call(rbind, out)
  
  dimnames(out) <- list(names(x), psNames)
  out
}
