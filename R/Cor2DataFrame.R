Cor2DataFrame <- function(x, n, v.na.replace=TRUE, row.names.unique=FALSE,
                          cor.analysis=TRUE, acov="weighted",
                          append.vars=TRUE, ...) {

    ## x is a list of "data", "n", ...
    if (all(c("data", "n") %in% names(x))) {
        my.cor <- x$data
        n <- x$n
        obslabels <- colnames(x$data[[1]])
    } else {
        ## x is just a list of correlation matrices. "n" is provided as an argument.
        my.cor <- x
        obslabels <- colnames(x[[1]])
    }

    if (length(my.cor) != length(n)) stop("Lengths of 'x' and 'n' are different.\n")
    
    if (cor.analysis) {
        ## Standardize and then vechs()
        my.df <- list2matrix(x=suppressWarnings(lapply(my.cor, cov2cor)), diag=FALSE)
    } else {
        ## vech()
        my.df <- list2matrix(x=my.cor, diag=TRUE)
    }

    acovR <- asyCov(x=my.cor, n=n, cor.analysis=cor.analysis, acov=acov, ...)

    ## NA is not allowed in definition variables
    ## They are replaced by 1e10
    if (v.na.replace) acovR[is.na(acovR)] <- 1e10

    ## x is a list of "data", "n", and moderators, and append
    ## Append the moderators x[-c(1,2)] into data
    if ( all(c(c("data", "n") %in% names(x),
               length(names(x))>2,
               append.vars)) )  {
        data <- suppressWarnings(data.frame(my.df, acovR, x[-c(1,2)], check.names=FALSE))        
    } else {
        data <- suppressWarnings(data.frame(my.df, acovR, check.names=FALSE))
    }
        
    ## Use unique row names if the row names are duplicated.
    if (row.names.unique) rownames(data) <- make.names(names(x), unique=TRUE)    

    list(data=data, n=n, obslabels=obslabels, ylabels=dimnames(my.df)[[2]],
         vlabels=dimnames(acovR)[[2]])
}
