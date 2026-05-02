#' Test Positive Definiteness of a List of Square Matrices
#' 
#' It tests the positive definiteness of a square matrix or a list of square
#' matrices. It returns \code{TRUE} if the matrix is positive definite. It
#' returns \code{FALSE} if the matrix is either non-positive definite or not
#' symmetric. Variables with \code{NA} in the diagonals will be removed before
#' testing. It returns \code{NA} when there are missing correlations even after
#' deleting the missing variables.
#' 
#' 
#' @param x A square matrix or a list of square matrices
#' @param check.aCov If it is \code{TRUE}, it mirrors the checking in
#' \code{\link[metaSEM]{asyCov}}.
#' @param cor.analysis Whether the input matrix is a correlation or a
#' covariance matrix. It is ignored when \code{check.aCov=FALSE}.
#' @param tol Tolerance (relative to largest variance) for numerical lack of
#' positive-definiteness in \code{x}. It is adopted from
#' \code{\link[MASS]{mvrnorm}}.
#' @return If the input is a matrix, it returns \code{TRUE}, \code{FALSE} or
#' \code{NA}. If the input is a list of matrices, it returns a list of
#' \code{TRUE}, \code{FALSE} or \code{NA}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' A <- diag(1,3)
#' is.pd(A)
#' # TRUE
#' 
#' B <- matrix(c(1,2,2,1), ncol=2)
#' is.pd(B)
#' # FALSE
#' 
#' is.pd(list(A, B))
#' # TRUE FALSE
#' 
#' C <- A
#' C[2,1] <- C[1,2] <- NA
#' is.pd(C)
#' # NA
#' 
is.pd <- function(x, check.aCov=FALSE, cor.analysis=TRUE, tol=1e-06) {
    if (is.list(x)) {
        return(sapply(x, is.pd, check.aCov=check.aCov, cor.analysis=cor.analysis, tol=tol))
    }
    else {
        ## Criteria based on asyCov()
        if (check.aCov) {
            if (cor.analysis) Diag(x)[is.na(Diag(x))] <- 1 else Diag(x)[is.na(Diag(x))] <- mean(Diag(x), na.rm=TRUE)
            x[is.na(x)] <- 0
        } else {
        ## Normal definition of pd    
            miss.index <- is.na(Diag(x)) 
            x <- x[!miss.index, !miss.index]
        }
        
        ## Catch the error when there are NA in the matrix
        lambda <- tryCatch(eigen(x, only.values = TRUE)$values, error=function(e) e)
        ## Return NA when there are NA in the matrix
        if (inherits(lambda, "error")) {
            out <- NA
        } else {
            # lambda_k/lambda_1 > tol
            ## if (lambda[length(lambda)]/lambda[1] > tol) {
          
            ## Use the definition in MASS::mvrnorm
            if (all(lambda >= -tol*abs(lambda[1L]))) {
                out <- TRUE
            } else {
                out <- FALSE
            }
        }
    }
    return(out)
}
