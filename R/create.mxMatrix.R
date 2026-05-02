#' Create a Vector into MxMatrix-class
#'
#' It converts a vector into \code{MxMatrix-class} via \code{mxMatrix}.
#'
#' If there are non-numeric values in \code{x}, they are treated as the labels
#' of the free parameters. If an "*" is present, the numeric value on the
#' left-hand side will be treated as the starting value for a free parameter or
#' a fixed value for a fixed parameter. If it is a matrix of numeric values,
#' there are no free parameters in the output matrix. \code{nrow} and
#' \code{ncol} will be calculated from the length of \code{x} unless
#' \code{type="Full"} is specified.
#'
#' @param x A character or numeric vector
#' @param type Matrix type similar to those listed in
#' \code{\link[OpenMx]{mxMatrix}}
#' @param ncol The number of columns. It is necessary when \code{type="Full"}.
#' It is ignored and determined by the length of \code{x} for the other types
#' of matrices.
#' @param nrow The number of rows. It is necessary when \code{type="Full"}. It
#' is ignored and determined by the length of \code{x} for the other types of
#' matrices.
#' @param as.mxMatrix Logical. If it is \code{TRUE}, the output is a matrix of
#' \code{MxMatrix-class}. If it is \code{FALSE}, it is a numeric matrix.
#' @param byrow Logical. If \code{FALSE} (the default) the matrix is filled by
#' columns, otherwise the matrix is filled by rows.
#' @param \dots Further arguments to be passed to
#' \code{\link[OpenMx]{mxMatrix}}. Please note that \code{type}, \code{nrow},
#' \code{ncol}, \code{values}, \code{free} and \code{labels} will be created
#' automatically. Thus, these arguments except labels should be avoided in
#' \dots
#' @return A \code{\link[OpenMx]{MxMatrix-class}} object with the same
#' dimensions as \code{x}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[OpenMx]{mxMatrix}},
#' \code{\link[metaSEM]{create.mxMatrix}},
#' \code{\link[metaSEM]{create.Fmatrix}}
#' @keywords utilities
#' @examples
#'
#' ## a and b are free parameters with starting values and labels
#' (a1 <- c(1:4, "5*a", 6, "7*b", 8, 9))
#'
#' (mat1 <- create.mxMatrix(a1, ncol=3, nrow=3, name="mat1"))
#'
#' ## Arrange the elements by row
#' (mat2 <- create.mxMatrix(a1, ncol=3, nrow=3, as.mxMatrix=FALSE, byrow=TRUE))
#'
#' (a3 <- c(1:3, "4*f4", "5*f5", "6*f6"))
#'
#' (mat3 <- create.mxMatrix(a3, type="Symm", name="mat3"))
#'
#' ## Create character matrix
#' (mat4 <- create.mxMatrix(a3, type="Symm", as.mxMatrix=FALSE))
#'
#' ## Arrange the elements by row
#' (mat5 <- create.mxMatrix(a3, type="Symm", as.mxMatrix=FALSE, byrow=TRUE))
#'
#' (mat6 <- create.mxMatrix(a3, type="Diag", lbound=6:1, name="mat6"))
#'
create.mxMatrix <- function(x, type=c("Full","Symm","Diag","Stand"), ncol=NA, nrow=NA, as.mxMatrix=TRUE, byrow=FALSE, ...) {
  my.mx <- function(y) {
    # suppress warnings
    ## warn <- options()$warn
    ## options(warn=-1)
    values <- suppressWarnings(as.numeric(y))  # They are NA for characters 
    free <- is.na(values)    # They are TRUE for parameters with labels 
    freePara1 <- y[free]     # Extract free parameters

    # check if there are any free parameters
    if (length(freePara1)>0) {
      freePara2 <- strsplit(freePara1, "*", fixed=TRUE)
      # replace NA with starting values 0.5 before "0.5*a"
      values[free] <- sapply(freePara2, function(y){ as.numeric(y[1]) })
      labels <- rep(NA, length(y))
      labels[free] <- sapply(freePara2, function(y){ y[2] })

      ## Replace TRUE by FALSE in free when there are definition variables
      free[grep("data.", labels)] <- FALSE
      
      out <- list(values=values, free=free, labels=labels)
    } else {
      out <- list(values=values, free=free, labels=NA)
    }
    ## options(warn=warn)
    out
  }

  type <- match.arg(type)
  
  if (as.mxMatrix) {
    
    switch(type,
           Full = { if (length(x) != ncol*nrow)
                      stop("Length of \"x\" does not match the dimensions of \"ncol\" and \"nrow\".\n")
                    untitled1 <- matrix(x, ncol=ncol, nrow=nrow, byrow=byrow)
                    untitled1 <- as.mxMatrix(untitled1, ...)
                  },
           Symm = { no.var <- (sqrt(1 + 8 * length(x)) - 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a symmetric matrix.\n")
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Symm", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, byrow=byrow, ...)   
                  },
           Diag = { no.var <- length(x)
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Diag", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, ...)                                         
                   },
           Stand = { no.var <- (sqrt(1 + 8 * length(x)) + 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a standardized matrix.\n")
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, byrow=byrow, ...)   
                  })
  } else {
        switch(type,
           Full = { if (length(x) != ncol*nrow)
                      stop("Length of \"x\" does not match the dimensions of \"ncol\" and \"nrow\".\n")
                    untitled1 <- matrix(x, ncol=ncol, nrow=nrow, byrow=byrow)                    
                  },
           Symm = { no.var <- (sqrt(1 + 8 * length(x)) - 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a symmetric matrix.\n")
                    untitled1 <- matrix(0, ncol=no.var, nrow=no.var)
                    if (byrow) {
                        untitled1[upper.tri(untitled1, diag = TRUE)] <- x
                        untitled1[lower.tri(untitled1)] <- t(untitled1)[lower.tri(untitled1)]
                    } else {
                        untitled1[lower.tri(untitled1, diag = TRUE)] <- x
                        untitled1[upper.tri(untitled1)] <- t(untitled1)[upper.tri(untitled1)]
                    }
                  },
           Diag = { untitled1 <- Diag(x)                                         
                   },
           Stand = { no.var <- (sqrt(1 + 8 * length(x)) + 1)/2
                     if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                       stop("Length of \"x\" does not match the no. of elements for a standardized matrix.\n")
                     untitled1 <- Diag(rep(1,no.var))
                     if (byrow) {
                         untitled1[upper.tri(untitled1, diag = FALSE)] <- x
                         untitled1[lower.tri(untitled1)] <- t(untitled1)[lower.tri(untitled1)]
                     } else {
                         untitled1[lower.tri(untitled1, diag = FALSE)] <- x
                         untitled1[upper.tri(untitled1)] <- t(untitled1)[upper.tri(untitled1)]  
                     }
                  })    
  }
  untitled1
}

