#' Create a Block Diagonal Matrix by Repeating the Input
#'
#' It creates a block diagonal matrix by repeating the input matrix several
#' times.
#'
#'
#' @param x A numeric or character matrix (or values)
#' @param times Number of times of \code{x} to be repeated.
#' @return A numeric or character block diagonal matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{bdiagMat}},
#' \code{\link[metaSEM]{matrix2bdiag}}
#' @keywords utilities
#' @examples
#'
#' ## Block diagonal matrix of numerics
#' bdiagRep( matrix(1:4,nrow=2,ncol=2), 2 )
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    3    0    0
#' # [2,]    2    4    0    0
#' # [3,]    0    0    1    3
#' # [4,]    0    0    2    4
#'
#' ## Block diagonal matrix of characters
#' bdiagRep( matrix(letters[1:4],nrow=2,ncol=2), 2 )
#' #      [,1] [,2] [,3] [,4]
#' # [1,] "a"  "c"  "0"  "0" 
#' # [2,] "b"  "d"  "0"  "0" 
#' # [3,] "0"  "0"  "a"  "c" 
#' # [4,] "0"  "0"  "b"  "d" 
#'
bdiagRep <- function(x, times) {
  bdiagMat( replicate(times, x, simplify=FALSE) )
}
