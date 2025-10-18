#Diag <- function(x, nrow, ncol) {
Diag <- function(x, ...) {
  ## if (inherits(x, "character") & missing(...)) {
  if (inherits(x, "character") && length(list(...)) == 0) {
    p <- length(x)
    out <- matrix(0, nrow=p, ncol=p)
    diag(out) <- x
  } else {
    out <- diag(x, ...)
    #out <- diag(x, nrow=nrow, ncol=ncol)
  }
  out
}

`Diag<-` <- function(x, value) {
  diag(x) <- value
  x
}
