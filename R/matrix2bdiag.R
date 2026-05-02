#' Convert a Matrix into a Block Diagonal Matrix
#' 
#' It converts a matrix into a block diagonal matrix.
#' 
#' Each row of \code{x} is converted into a symmetric matrix via
#' \code{\link[metaSEM]{vec2symMat}}. Then the list of the symmetric matrices
#' is converted into a block diagonal matrix via a function written by Scott
#' Chasalow posted at
#' http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html.
#' 
#' @param x A \eqn{k}{k} x \eqn{p}{p} matrix of numerics or characters.
#' @param \dots Further arguments to be passed to
#' \code{\link[metaSEM]{vec2symMat}}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{vec2symMat}}
#' @keywords utilities
#' @examples
#' 
#' (m1 <- matrix(1:12, ncol=6, byrow=TRUE))
#' #      [,1] [,2] [,3] [,4] [,5] [,6]
#' # [1,]    1    2    3    4    5    6
#' # [2,]    7    8    9   10   11   12
#' 
#' matrix2bdiag(m1)
#' #      [,1] [,2] [,3] [,4] [,5] [,6]
#' # [1,]    1    2    3    0    0    0
#' # [2,]    2    4    5    0    0    0
#' # [3,]    3    5    6    0    0    0
#' # [4,]    0    0    0    7    8    9
#' # [5,]    0    0    0    8   10   11
#' # [6,]    0    0    0    9   11   12
#' 
matrix2bdiag <- function(x, ...) {
  tmp <- split(as.matrix(x), row(x))
  # Use bdiagMat() to handle string matrices
  out <- bdiagMat(lapply(tmp, vec2symMat, ...))
  as.matrix(out)
}
