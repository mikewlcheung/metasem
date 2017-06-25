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
