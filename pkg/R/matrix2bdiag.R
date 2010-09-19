matrix2bdiag <- function(x, ...) {
  tmp <- split(as.matrix(x), row(x))
  out <- bdiag(lapply(tmp, vec2symMat, ...))
  as.matrix(out)
}
