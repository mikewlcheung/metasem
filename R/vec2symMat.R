#' Convert a Vector into a Symmetric Matrix
#'
#' It converts a vector into a symmetric matrix by filling up the elements into
#' the lower triangle of the matrix.
#'
#'
#' @param x A vector of numerics or characters
#' @param diag Logical. If it is \code{TRUE} (the default), the diagonals of
#' the created matrix are replaced by elements of \code{x}; otherwise, the
#' diagonals of the created matrix are replaced by "1".
#' @param byrow Logical. If it is \code{FALSE} (the default), the created
#' matrix is filled by columns; otherwise, the matrix is filled by rows.
#' @return A symmetric square matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{matrix2bdiag}}
#' @keywords utilities
#' @examples
#'
#' vec2symMat(1:6)
#' #      [,1] [,2] [,3]
#' # [1,]    1    2    3
#' # [2,]    2    4    5
#' # [3,]    3    5    6
#'
#' vec2symMat(1:6, diag=FALSE)
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    1    2    3
#' # [2,]    1    1    4    5
#' # [3,]    2    4    1    6
#' # [4,]    3    5    6    1
#'
#' vec2symMat(letters[1:6])
#' #      [,1] [,2] [,3]
#' # [1,] "a"  "b"  "c" 
#' # [2,] "b"  "d"  "e" 
#' # [3,] "c"  "e"  "f" 
#'
vec2symMat <- function (x, diag=TRUE, byrow=FALSE) {
  m <- length(x)
  d <- if (diag) 1 else -1
  n <- floor((sqrt(1 + 8*m) - d)/2)
  if (m != n*(n + d)/2) 
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- Diag(n)

  ## Row major
  if (byrow) {
    mat[upper.tri(mat, diag=diag)] <- x
    index <- lower.tri(mat)
    mat[index] <- t(mat)[index]  
  } else {
  ## Column major: default behavior
    mat[lower.tri(mat, diag=diag)] <- x
    # Just mirroring the matrix, exclude the diagonals
    ## mat[upper.tri(mat, diag=FALSE)] <- mat[lower.tri(mat, diag=FALSE)]
    ## Corrected a bug
    index <- upper.tri(mat)
    mat[index] <- t(mat)[index]  
  }
  mat
}
