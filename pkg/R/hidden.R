## http://tolstoy.newcastle.edu.au/R/e6/help/09/05/15051.html
## Starting values based on simple average of all available cases
## http://tolstoy.newcastle.edu.au/R/e10/help/10/04/2866.html
# a <- array(unlist(X), c(9,9,11))
# apply(a, c(1,2),mean, na.rm=TRUE)
.startValues <- function(x, cor.analysis = TRUE) {
    no.var <- max(sapply(x, ncol))
    my.start <- matrix(rowMeans(sapply(x, as.vector), na.rm = TRUE), nrow = no.var)
    out <- as.matrix(nearPD(my.start, corr = cor.analysis)$mat)
    dimnames(out) <- dimnames(x[[1]])
    out
}

.indepwlsChisq <- function(S, acovS, cor.analysis = TRUE) {
    no.var <- ncol(S)
    sampleS <- mxMatrix("Full", ncol = no.var, nrow = no.var, values = c(S), free = FALSE, 
        name = "sampleS")
    if (cor.analysis) {
        impliedS <- mxMatrix("Iden", ncol = no.var, nrow = no.var, free = FALSE, 
            name = "impliedS")          # a bit redundant, but it simplies the later codes
        ps <- no.var * (no.var - 1)/2
        vecS <- mxAlgebra(vechs(sampleS - impliedS), name = "vecS")        
    } else {
        impliedS <- mxMatrix("Diag", ncol = no.var, nrow = no.var, values = diag(S), 
            free = TRUE, name = "impliedS")
        ps <- no.var * (no.var + 1)/2
        vecS <- mxAlgebra(vech(sampleS - impliedS), name = "vecS")        
    }
    if (ncol(acovS) != ps) 
        stop("No. of dimension of \"S\" does not match the dimension of \"acovS\"\n")
    
    # Inverse of asymptotic covariance matrix
    invacovS <- tryCatch(solve(acovS), error = function(e) e)
    if (inherits(invacovS, "error")) 
        stop(print(invacovS))    
    invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invacovS), free = FALSE, 
        name = "invAcov")
    obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")
    objective <- mxAlgebraObjective("obj")
    
    indep.out <- tryCatch(mxRun(mxModel(model = "Independent model", impliedS, sampleS, 
        vecS, invAcov, obj, objective), silent=TRUE, suppressWarnings=TRUE), error = function(e) e)
    
    if (inherits(indep.out, "error")) 
        stop(print(indep.out))
    else return(indep.out@output$Minus2LogLikelihood)
}

.minus2LL <- function(x, n, model=c("saturated", "independent")) {
  if (is.list(x)) {
    if (length(x) != length(n))
      stop("Lengths of \"x\" and \"n\" are not the same.\n")
    return(mapply(.minus2LL, x=x, n=n, model=model))
  } else {
    miss.index <- is.na(diag(x))
    x <- x[!miss.index, !miss.index]
    if (!is.pd(x))
      stop("\"x\" is not positive definite.\n")
    no.var <- ncol(x)
    vars <- paste("v", 1:no.var, sep="")
    dimnames(x) <- list(vars, vars)
    obsCov <- mxData(observed=x, type='cov', numObs=n)
    ## if (missing(model))
    ##   stop("\"model\" was not specified.\n")
    model <- match.arg(model)
    switch(model,
           saturated = expCov <- mxMatrix("Symm", nrow=no.var, ncol=no.var, free=TRUE, 
                                          value=vech(x), name="expCov"),
           independent = expCov <- mxMatrix("Diag", nrow=no.var, ncol=no.var, free=TRUE, 
                                            value=diag(x), name="expCov")
     )
    objective <- mxMLObjective(covariance = "expCov", dimnames=vars)
    fit <- tryCatch(mxRun(mxModel("model", expCov, obsCov, objective), silent=TRUE, 
                    suppressWarnings=TRUE), error = function(e) e)
    if (inherits(fit, "error")) 
          stop(print(fit))
    fit@output$Minus2LogLikelihood
  }
}

## http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html
## http://tolstoy.newcastle.edu.au/R/help/04/05/1322.html
## It works for character matrices.
