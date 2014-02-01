is.pd <- function(x, tol = 1e-06) {
    if (is.list(x)) {
        return(sapply(x, is.pd, tol = tol))
    }
    else {
        miss.index <- is.na(Diag(x)) 
        x <- x[!miss.index, !miss.index]
        ## Catch the error when there are NA in the matrix
        lambda <- tryCatch(eigen(x, only.values = TRUE)$values, error=function(e) e)
        ## Return NA when there are NA in the matrix
        if (inherits(lambda, "error")) {
            out <- NA
        } else {
            # lambda_k/lambda_1 > tol
            if (lambda[length(lambda)]/lambda[1] > tol) {
                out <- TRUE
            } else {
                out <- FALSE
            }
        }
    }
    return(out)
}
