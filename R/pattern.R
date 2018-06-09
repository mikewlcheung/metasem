pattern.n <- function(x, n) {    
    if (!is.list(x)) stop("\"x\" must be a list of matrices.\n")
    if (length(x)!=length(n)) stop("The lengths of \"x\" and \"n\" must be the same.\n")
    
    fun <- function(x1, n1) {
        ## x2: a copy of x1
        x2 <- x1
        ## replace NA with 0
        x2[is.na(x1)] <-0
        ## replace not NA with the sample size
        x2[!is.na(x1)] <- n1
        x2}
    my.df <- mapply(fun, x, n, SIMPLIFY=FALSE)
    Reduce('+', my.df)
}

pattern.na <- function(x, show.na=TRUE,
                       type=c("tssem", "osmasem")) {

    type <- match.arg(type)

    if (type=="tssem") {
        out <- Reduce("+", lapply(x, is.na))
        if (show.na==FALSE) {
            out <- length(x)-out
        }    
    } else {
        out <- x$data[, x$ylabels]
        out <- split(out, seq_len(nrow(out)))
        out <- lapply(out, function(x) { x <- is.na(x)
                        matrix(x, ncol=1, nrow=length(x)) %*% 
                        matrix(x, ncol=length(x), nrow=1) })
        out <- Reduce("+", out)
        dimnames(out) <- list(x$ylabels, x$ylabels)
        if (show.na==FALSE) {
            out <- nrow(x$data)-out
        }
    }    
    out
} 

