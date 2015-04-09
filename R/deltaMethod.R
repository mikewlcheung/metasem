deltaMethod <- function (fn, variables, Mean, Cov) {   
    Cov <- as.matrix(Cov)
    n <- length(Mean)
    if ((dim(Cov)[1] != n) || (dim(Cov)[2] != n)) 
        stop("Dimensions of \"Mean\" and \"Cov\" are different.\n")          
    if (!is.list(fn))
      fn <- list(fn)
    ## Convert into formulas
    fn <- lapply(fn, function(x) as.formula(paste("~", x, sep="")))
    for (i in 1:n) assign(variables[i], Mean[i])
        grad <- t(sapply(fn, function(x) as.numeric(attr(eval(deriv(x, variables)), "gradient"))))
    grad %*% Cov %*% t(grad)
}
