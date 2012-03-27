vec2symMat <- function (x, diag=TRUE) {
  m <- length(x)
  d <- if (diag) 1 else -1
  n <- floor((sqrt(1 + 8*m) - d)/2)
  if (m != n*(n + d)/2) 
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- diag(n)

  mat[lower.tri(mat, diag=diag)] <- x
  # Just mirroring the matrix, exclude the diagonals
  ## mat[upper.tri(mat, diag=FALSE)] <- mat[lower.tri(mat, diag=FALSE)]
  ## Corrected a bug
  index <- upper.tri(mat)
  mat[index] <- t(mat)[index]  
  mat
}
