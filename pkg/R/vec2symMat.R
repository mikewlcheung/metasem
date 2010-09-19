vec2symMat <- function (x, diag=TRUE, byrow=FALSE) {
  m <- length(x)
  d <- if (diag) 1 else -1
  n <- floor((sqrt(1 + 8*m) - d)/2)
  if (m != n*(n + d)/2) 
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- diag(n)

  if (byrow) {
    mat[upper.tri(mat, diag=diag)] <- x
  } else {
    mat[lower.tri(mat, diag=diag)] <- x
  }
  mat+t(mat)-diag(diag(mat))
}
