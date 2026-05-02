#Diag <- function(x, nrow, ncol) {


#' Matrix Diagonals
#'
#' Extract or replace the diagonal of a matrix, or construct a diagonal matrix
#' with the same behaviors as \code{diag} prior to R-3.0.0.
#'
#' Starting from R-3.0.0, \code{diag(x)} returns a numeric matrix with NA in the
#' diagonals when x is a character vector. Although this follows what the
#' manual says, this breaks metaSEM. The \code{Diag} has the same functions
#' as \code{diag} except that \code{Diag(x)} works for a character vector of x
#' by returning a square matrix of character "0" with \code{x} as the
#' diagonals.
#'
#' @aliases Diag Diag<-
#' @param x A matrix, vector or 1D array, or missing.
#' @param ... Optional dimensions (\code{nrow} and \code{ncol}) for the result
#' when \code{x} is not a matrix.
#' @param value Either a single value or a vector of length equal to that of
#' the current diagonal. Should be of a mode which can be coerced to that of
#' \code{x}.
#' @note See
#' http://r.789695.n4.nabble.com/Behaviors-of-diag-with-character-vector-in-R-3-0-0-td4663735.html
#' for the discussion.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link{diag}}
#' @keywords utilities
#' @examples
#'
#' v <- c("a", "b")
#' Diag(v)
#'
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

#' @rdname Diag
`Diag<-` <- function(x, value) {
  diag(x) <- value
  x
}
