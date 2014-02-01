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
  sampleS <- mxMatrix("Full", ncol = no.var, nrow = no.var, values = c(S),
                      free = FALSE, name = "sampleS")
  if (cor.analysis) {
    impliedS <- mxMatrix("Iden", ncol = no.var, nrow = no.var, free = FALSE,
                         name = "impliedS")     # a bit redundant, but it simplies the codes later
    ps <- no.var * (no.var - 1)/2
    vecS <- mxAlgebra(vechs(sampleS - impliedS), name = "vecS")
  } else {
    impliedS <- mxMatrix("Diag", ncol = no.var, nrow = no.var, values = Diag(S),
                         free = TRUE, name = "impliedS")
    ps <- no.var * (no.var + 1)/2
    vecS <- mxAlgebra(vech(sampleS - impliedS), name = "vecS")
  }
  if (ncol(acovS) != ps)
    stop("Dimensions of \"S\" do not match the dimensions of \"acovS\"\n")
    
  # Inverse of asymptotic covariance matrix
  invacovS <- tryCatch(solve(acovS), error = function(e) e)
  if (inherits(invacovS, "error"))
    stop(print(invacovS))    
  invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invacovS),
                      free = FALSE, name = "invAcov")
  obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")
  objective <- mxAlgebraObjective("obj")
  mx.model <- mxModel(model = "Independent model", impliedS, sampleS,
                      vecS, invAcov, obj, objective)
  mx.model <- mxOption(mx.model, "Calculate Hessian", "No") 
  mx.model <- mxOption(mx.model, "Standard Errors"  , "No")     
    
  indep.out <- tryCatch(mxRun(mx.model, silent=TRUE,
                              suppressWarnings=TRUE), error = function(e) e)
    
  if (inherits(indep.out, "error"))
    stop(print(indep.out))
  else return(indep.out@output$Minus2LogLikelihood)
}

.minus2LL <- function(x, n) {
## model=c("saturated", "independent")
  if (is.list(x)) {
    if (length(x) != length(n))
      stop("Lengths of \"x\" and \"n\" are not the same.\n")
    fit.list <- mapply(.minus2LL, x=x, n=n)
    out <- apply(matrix(unlist(fit.list), nrow=3), 1, sum)
    names(out) <- c("SaturatedLikelihood", "IndependenceLikelihood", "independenceDoF")
    out
  } else {
    miss.index <- is.na(Diag(x))
    x <- x[!miss.index, !miss.index]
    if (!is.pd(x))
      stop("\"x\" is not positive definite.\n")
    no.var <- ncol(x)
    vars <- paste("v", 1:no.var, sep="")
    dimnames(x) <- list(vars, vars)
    obsCov <- mxData(observed=x, type='cov', numObs=n)
    ## if (missing(model))
    ##   stop("\"model\" was not specified.\n")
    ## model <- match.arg(model)
    ## switch(model,
    ##        saturated = expCov <- mxMatrix("Symm", nrow=no.var, ncol=no.var, free=TRUE, 
    ##                                       values=vech(x), name="expCov"),
    ##        independent = expCov <- mxMatrix("Diag", nrow=no.var, ncol=no.var, free=TRUE, 
    ##                                         values=diag(x), name="expCov")
    ##  )
    expCov <- mxMatrix("Diag", nrow=no.var, ncol=no.var, free=TRUE, 
                       values=Diag(x), name="expCov")    
    objective <- mxMLObjective(covariance = "expCov", dimnames=vars)
    mx.model <- mxModel("model", expCov, obsCov, objective)
    mx.model <- mxOption(mx.model, "Calculate Hessian", "No") 
    mx.model <- mxOption(mx.model, "Standard Errors"  , "No")
    
    fit <- tryCatch(mxRun(mx.model, silent=TRUE, suppressWarnings=TRUE), error = function(e) e)
    if (inherits(fit, "error")) 
          stop(print(fit))
    #fit@output$Minus2LogLikelihood
    ## c(unlist(summary(fit)[c("SaturatedLikelihood", "IndependenceLikelihood", "independenceDoF")]))
    c(unlist(fit@output[c("SaturatedLikelihood", "IndependenceLikelihood")]), independenceDoF=no.var*(no.var-1)/2)
  }
}

## http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html
## http://tolstoy.newcastle.edu.au/R/help/04/05/1322.html
## It works for character matrices.

## Calculate R2 for meta and meta3 objects
.R2 <- function(object) {
  no.y <- object$no.y
    ## meta3 or meta3ML class
  if ( any(c("meta3","meta3X") %in% class(object)) ) {
    ## Tau2 with predictors
    Tau2_2model <- tryCatch( eval(parse(text="mxEval(Tau2_2, object$mx.fit)")), error = function(e) NA )
    Tau2_3model <- tryCatch( eval(parse(text="mxEval(Tau2_3, object$mx.fit)")), error = function(e) NA )
    Tau2_2base <- tryCatch( eval(parse(text="mxEval(Tau2_2, object$mx0.fit$mx.fit)")), error = function(e) NA )
    Tau2_3base <- tryCatch( eval(parse(text="mxEval(Tau2_3, object$mx0.fit$mx.fit)")), error = function(e) NA )           

    R2_2 <- max((1-Tau2_2model/Tau2_2base), 0)
    R2_3 <- max((1-Tau2_3model/Tau2_3base), 0)
    R2.values <- matrix(c(Tau2_2base, Tau2_2model, R2_2, Tau2_3base, Tau2_3model, R2_3), ncol=2)
    dimnames(R2.values) <- list(c("Tau2 (no predictor)", "Tau2 (with predictors)", "R2"),
                                c("Level 2", "Level 3"))
    } else {
    ## Return NA when Tau2 values cannot be retrieved, e.g., constrained at lbound or no name
    Tau2model <- sapply(1:no.y, function(x) { temp <- paste("mxEval(Tau2_",x,"_",x,", object$mx.fit)", sep="")
                        tryCatch( eval(parse(text=temp)), error = function(e) NA ) })
    Tau2base <- sapply(1:no.y, function(x) { temp <- paste("mxEval(Tau2_",x,"_",x,", object$mx0.fit$mx.fit)", sep="")
                       tryCatch( eval(parse(text=temp)), error = function(e) NA ) })
    R2_2 <- pmax((1-Tau2model/Tau2base), 0)
    R2.values <- rbind(Tau2base, Tau2model, R2_2)
    dimnames(R2.values) <- list( c("Tau2 (no predictor)", "Tau2 (with predictors)", "R2"),
                                 paste("y", 1:no.y, sep="") )       
    ## dimnames(R2.values) <- list( paste("y", c(t(outer(1:no.y, c(": Tau2 (no predictor)", ": Tau2 (with predictors)", ": R2"),
    ##                                                   paste, sep=""))), sep=""), c("Estimate")) 
    }
  
  R2.values
} 

## Calculate I2 for meta and meta3 objects
.I2 <- function(object, my.mx) {
  I2 <- object$I2
  no.y <- object$no.y
  if (missing(my.mx)) {
    my.mx <- summary(object$mx.fit)
  }
  my.ci <- my.mx$CI
   
  ## meta3 class
  if ( "meta3" %in% class(object) ) {
    I2.names <- c("I2q","I2hm","I2am")
    ## Rearrange different I2 into the standard format
    I2.names <- I2.names[ c("I2q","I2hm","I2am") %in% I2 ]
    I2.names <- c(outer(I2, c("_2","_3"), paste, sep=""))
   
    ## Wald test, no CI
    if (is.null(dimnames(my.ci))) {
      I2.values <- matrix(NA, nrow=length(I2.names), ncol=3)
      I2.values[,2] <- eval(parse(text = paste("mxEval(c(", paste(I2.names, collapse=","), "), object$mx.fit)", sep="")))
    } else {
      dimnames(my.ci)[[1]] <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                                    {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      I2.values <- my.ci[paste(I2.names, "[1,1]", sep=""), , drop=FALSE]

    }
    ## Truncate within 0 and 1
    ## I2.values <- apply(I2.values, c(1,2), function(x) max(x,0))
    I2.values <- ifelse(I2.values<0, 0, I2.values)
    I2.values <- ifelse(I2.values>1, 1, I2.values)        
    I2.names <- sub("I2q_2", "I2_2 (Typical v: Q statistic)", I2.names)
    I2.names <- sub("I2q_3", "I2_3 (Typical v: Q statistic)", I2.names)
    I2.names <- sub("I2hm_2", "I2_2 (Typical v: harmonic mean)", I2.names)
    I2.names <- sub("I2hm_3", "I2_3 (Typical v: harmonic mean)", I2.names)
    I2.names <- sub("I2am_2", "I2_2 (Typical v: arithmetic mean)", I2.names)
    I2.names <- sub("I2am_3", "I2_3 (Typical v: arithmetic mean)", I2.names)
    I2.names <- sub("ICC_2", "ICC_2 (tau^2/(tau^2+tau^3))", I2.names)
    I2.names <- sub("ICC_3", "ICC_3 (tau^3/(tau^2+tau^3))", I2.names)    
    dimnames(I2.values) <- list(I2.names, c("lbound", "Estimate", "ubound"))
    
    } else {
    ## meta class
               
    I2.names <- c("I2q","I2hm","I2am")
    ## Rearrange different I2 into the standard format
    I2.names <- I2.names[ c("I2q","I2hm","I2am") %in% I2 ]

    ## Wald test, no CI
    if (is.null(dimnames(my.ci))) {
      I2.values <- matrix(NA, nrow=length(I2.names)*no.y, ncol=3)
      I2.values[,2] <- eval(parse(text = paste("mxEval(I2_values, object$mx.fit)", sep="")))
    } else {## LB CI  model.name <- "Meta analysis with ML."
      dimnames(my.ci)[[1]] <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                                    {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      I2.values <- my.ci[paste("I2_values[", 1:(length(I2.names)*no.y), ",1]", sep=""), ,drop=FALSE]
    }    
    ## Truncate into 0
    ## I2.values <- apply(I2.values, c(1,2), function(x) max(x,0))
    I2.values <- ifelse(I2.values<0, 0, I2.values)
    I2.values <- ifelse(I2.values>1, 1, I2.values) 
    I2.names <- paste(": ", I2.names, sep="")
    I2.names <- paste("Intercept", c( outer(1:no.y, I2.names, paste, sep="")), sep="")
    I2.names <- sub("I2q", "I2 (Q statistic)", I2.names)
    I2.names <- sub("I2hm", "I2 (harmonic mean)", I2.names)
    I2.names <- sub("I2am", "I2 (arithmetic mean)", I2.names)
    dimnames(I2.values) <- list(I2.names, c("lbound", "Estimate", "ubound"))
    }
  
  I2.values
}

## Calculate asymptotic sampling covariance matrix
## fn=list("x1*x2/x3", "x2^x1/x3")
## variables=c("x1", "x2", "x3")
## Mean and Cov based on variables
.delta <- function (fn, variables, Mean, Cov) { 
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

.asyCov <- function(x, n, cor.analysis=TRUE, type=c("individual", "fixed", "random")) {
    type <- match.arg(type)
    switch( type,
           individual = {
               ## Diagonals of starting values
               my.start <- diag(.startValues(x, cor.analysis=TRUE))
               ## Replace diagonals with 1.0
               my.cor <- lapply(x, function (x) { my.na <- is.na(Diag(x))
                                                  Diag(x)[my.na] <- my.start[my.na]
                                                  x })
               ## Replace missing variables with 0.0
               my.cor <- lapply(my.cor, function (x) { x[is.na(x)] <- 0; x }) },
           fixed = {
               my.fixed <- tssem1REM(my.df=x, n=n, cor.analysis=cor.analysis, RE.type="Zero", silent=TRUE)
               my.cor <- vec2symMat(coef(my.fixed, select="fixed"), diag=!cor.analysis)
               my.cor <- rep(list(my.cor), length(n)) },
           random = {
               my.random <- tssem1REM(my.df=x, n=n, cor.analysis=cor.analysis, RE.type="Diag", silent=TRUE)
               my.cor <- vec2symMat(coef(my.random, select="fixed"), diag=!cor.analysis)
               my.cor <- rep(list(my.cor), length(n)) })

    ## Calculate the asymptotic sampling covariance matrix of the correlation matrix
    asyCov(x=my.cor, n=n, cor.analysis=cor.analysis)
}
