## http://tolstoy.newcastle.edu.au/R/e6/help/09/05/15051.html
## Starting values based on simple average of all available cases
## http://tolstoy.newcastle.edu.au/R/e10/help/10/04/2866.html
## Not fully tested for missing data yet
# a <- array(unlist(X), c(9,9,11))
# apply(a, c(1,2),mean, na.rm=TRUE)
startValues <- function(x, cor.analysis = TRUE) {
    no.var <- max(sapply(x, ncol))
    my.start <- matrix(rowMeans(sapply(x, as.vector), na.rm = TRUE), nrow = no.var)
    out <- as.matrix(nearPD(my.start, corr = cor.analysis)$mat)
    dimnames(out) <- dimnames(x[[1]])
    out
}

