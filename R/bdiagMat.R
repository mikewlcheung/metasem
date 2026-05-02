#' Create a Block Diagonal Matrix
#'
#' It creates a block diagonal matrix from a list of numeric or character
#' matrices.
#'
#'
#' @param x A list of numeric or character matrices (or values)
#' @return A numeric or character block diagonal matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{bdiagRep}},
#' \code{\link[metaSEM]{matrix2bdiag}}
#' @references It was based on a function posted by Scott Chasalow at
#' http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html.
#' @keywords utilities
#' @examples
#'
#' ## Block diagonal matrix of numbers
#' bdiagMat( list(matrix(1:4,nrow=2,ncol=2),
#'                matrix(5:6,nrow=1,ncol=2)) )
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    3    0    0
#' # [2,]    2    4    0    0
#' # [3,]    0    0    5    6
#'
#' ## Block diagonal matrix of characters
#' bdiagMat( list(matrix(letters[1:4],nrow=2,ncol=2),
#'                matrix(letters[5:6],nrow=1,ncol=2)) )
#' #      [,1] [,2] [,3] [,4]
#' # [1,] "a"  "c"  "0"  "0" 
#' # [2,] "b"  "d"  "0"  "0" 
#' # [3,] "0"  "0"  "e"  "f" 
#'
bdiagMat <- function(x){
  if(!is.list(x)) stop("\"x\" must be a list.")
  n <- length(x)
  if(n==0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
  (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
}
