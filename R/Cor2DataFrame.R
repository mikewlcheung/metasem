Cor2DataFrame <- function(x, n, v.na.replace=TRUE, row.names.unique=FALSE,
                          cor.analysis=TRUE, acov="weighted", ...) {
    if (length(x) != length(n)) stop("Lengths of 'x' and 'n' are different.\n")
    
    if (cor.analysis) {
        my.df <- list2matrix(x=suppressWarnings(lapply(x, cov2cor)), diag=FALSE)
    } else {
        my.df <- list2matrix(x=x, diag=TRUE)
    }

    acovR <- asyCov(x=x, n=n, cor.analysis=cor.analysis, acov=acov, ...)

    ## NA is not allowed in definition variables
    ## They are replaced by 1e10
    if (v.na.replace) acovR[is.na(acovR)] <- 1e10

    data <- suppressWarnings(data.frame(my.df, acovR, check.names=FALSE))
    
    ## Use unique row names if the row names are duplicated.
    if (row.names.unique) rownames(data) <- make.names(names(x), unique=TRUE)    

    list(data=data, n=n, obslabels=colnames(x[[1]]),
         ylabels=dimnames(my.df)[[2]], vlabels=dimnames(acovR)[[2]])
}

tssem2DataFrame <- function(X, v.na.replace = TRUE, row.names.unique = FALSE,
                            cor.analysis = TRUE, acov="weighted", ...) {
    
    out <- Cor2DataFrame(x=X$data, n=X$n, v.na.replace=v.na.replace,
                         row.names.unique=row.names.unique,
                         cor.analysis=cor.analysis, acov=acov, ...)

    ## Variable names in the inputs > 2
    if (length(names(X)) > 2) {
        out$data <- data.frame(out$data, data.frame(X[-c(1,2)]), check.names=FALSE)
    }

    out
}
