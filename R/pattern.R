#' Display the Accumulative Sample Sizes for the Covariance Matrix
#' 
#' It displays the accumulative sample sizes for the covariance matrix.
#' 
#' 
#' @param x A list of square matrices
#' @param n A vector of sample sizes.
#' @return A square matrix of the accumulative sample sizes of the input
#' matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' ## Show the pattern of missing data
#' pattern.n(Hunter83$data, Hunter83$n)
#' 
#' #             Ability Knowledge Work sample Supervisor
#' # Ability        3815      3372        3281       3605
#' # Knowledge      3372      3532        2998       3322
#' # Work sample    3281      2998        3441       3231
#' # Supervisor     3605      3322        3231       3765
#' 
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



#' Display the Pattern of Missing Data of a List of Square Matrices
#' 
#' It displays the pattern of missing data (or pattern of data that are
#' present) of a list of square matrices with the same dimensions.
#' 
#' 
#' @param x A list of square matrices
#' @param show.na If it is \code{TRUE}, it shows the pattern of missing data.
#' If it is \code{FALSE}, it shows the pattern of data that are present.
#' @param type If it is \code{tssem}, it reports the pattern of missing
#' correlations for the tssem approach. If it is \code{osmasem}, it reports the
#' pattern of missing correlations for the data created by
#' \code{\link[metaSEM]{Cor2DataFrame}}.
#' @return A square matrix of numerical values with the same dimensions of the
#' input matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' ## Show the pattern of missing data
#' pattern.na(Hunter83$data, show.na=TRUE)
#' 
#' #             Ability Knowledge Work sample Supervisor
#' # Ability           1         3           3          2
#' # Knowledge         3         2           4          3
#' # Work sample       3         4           2          3
#' # Supervisor        2         3           3          1
#' 
#' ## Show the pattern of data that are present
#' pattern.na(Hunter83$data, show.na=FALSE)
#' 
#' #             Ability Knowledge Work sample Supervisor
#' # Ability          13        11          11         12
#' # Knowledge        11        12          10         11
#' # Work sample      11        10          12         11
#' # Supervisor       12        11          11         13
#' 
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

