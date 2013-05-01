is.pd <- function(x, tol = 1e-06) {
    if (is.list(x)) {
        return(sapply(x, is.pd, tol = tol))
    }
    else {
        miss.index <- is.na(Diag(x)) 
        x <- x[!miss.index, !miss.index]
        lambda <- eigen(x, only.values = TRUE)$values
        # lambda_k/lambda_1 > tol
        if (lambda[length(lambda)]/lambda[1] > tol) {
            return(TRUE)
        }
        else {
            return(FALSE)
        }
    }
}
