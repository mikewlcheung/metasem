matrix2bdiag <- function(x, ...) {
  tmp <- split(as.matrix(x), row(x))
  # Use .bdiag() as it handles string matrices
  out <- .bdiag(lapply(tmp, vec2symMat, ...))
  as.matrix(out)
}
