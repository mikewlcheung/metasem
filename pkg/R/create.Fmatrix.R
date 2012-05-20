create.Fmatrix <- function(x, name) {
  x <- as.logical(x)
  Fmatrix <- diag(as.numeric(x))[x, , drop=FALSE]
  if (missing(name)) {
    as.mxMatrix(Fmatrix)
  } else {
    as.mxMatrix(Fmatrix, name=name)
  }
}
