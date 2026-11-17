## http://tolstoy.newcastle.edu.au/R/e6/help/09/05/15051.html
## Starting values based on simple average of all available cases
## http://tolstoy.newcastle.edu.au/R/e10/help/10/04/2866.html
# a <- array(unlist(X), c(9,9,11))
# apply(a, c(1,2),mean, na.rm=TRUE)
.startValues <- function(x, cor.analysis=TRUE, n=NULL) {
  no.var <- max(sapply(x, ncol))
  my.data <- do.call(cbind, lapply(x, as.vector))
  if (is.null(n)) {
      ## Unweighted mean
      my.start <- apply(my.data, 1, weighted.mean, na.rm=TRUE)
  } else {
      ## Weigted mean by sample sizes
      my.start <- apply(my.data, 1, weighted.mean, w=n, na.rm=TRUE)
  }
  my.start <- matrix(my.start, nrow = no.var)
  out <- as.matrix(nearPD(my.start, corr = cor.analysis)$mat)
  dimnames(out) <- dimnames(x[[1]])
  out
}

.indepwlsChisq <- function(S, aCov, cor.analysis = TRUE) {
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
  if (ncol(aCov) != ps)
    stop("Dimensions of \"S\" do not match the dimensions of \"aCov\"\n")
    
  # Inverse of asymptotic covariance matrix
  invaCov <- tryCatch(chol2inv(chol(aCov)), error = function(e) e)
  if (inherits(invaCov, "error"))
    stop(print(invaCov))    
  invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invaCov),
                      free = FALSE, name = "invAcov")
  obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")

  mx.model <- mxModel(model = "Independent model", impliedS, sampleS,
                      vecS, invAcov, obj, mxFitFunctionAlgebra("obj"))
  mx.model <- mxOption(mx.model, "Calculate Hessian", "No") 
  mx.model <- mxOption(mx.model, "Standard Errors"  , "No")     
    
  indep.out <- tryCatch(mxRun(mx.model, silent=TRUE,
                              suppressWarnings=TRUE), error = function(e) e)
    
  if (inherits(indep.out, "error")) {
      return(NA)
      ## stop(print(indep.out))
  }  else {
      return(indep.out@output$Minus2LogLikelihood)
  }
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
    x <- x[!miss.index, !miss.index, drop=FALSE]
    if (!isTRUE(is.pd(x)))
      stop("\"x\" is not positive definite.\n")
    no.var <- ncol(x)
    vars <- paste("v", 1:no.var, sep="")
    dimnames(x) <- list(vars, vars)
    obsCov <- mxData(observed=x, type='cov', numObs=n)
    objective <- mxExpectationNormal(covariance = "expCov", dimnames=vars)
    fit.baseline <- function(type) {
      expCov <- mxMatrix(type, nrow=no.var, ncol=no.var, free=TRUE,
                         values=if (type == "Symm") vech(x) else diag(x), name="expCov")
      model <- mxModel("model", expCov, obsCov, objective, mxFitFunctionML())
      model <- mxOption(model, "Calculate Hessian", "No")
      model <- mxOption(model, "Standard Errors", "No")
      fit <- mxRun(model, silent=TRUE, suppressWarnings=TRUE)
      if (!is.finite(fit@output$fit) || !fit@output$status$code %in% c(0L, 1L))
        stop("The ", type, " baseline model did not converge; fit statistics are unavailable.")
      fit@output$fit
    }
    c(SaturatedLikelihood=fit.baseline("Symm"),
      IndependenceLikelihood=fit.baseline("Diag"),
      independenceDoF=no.var*(no.var-1)/2)
  }
}

## http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html
## http://tolstoy.newcastle.edu.au/R/help/04/05/1322.html
## It works for character matrices.

## Calculate R2 for meta and meta3L objects
.R2 <- function(object) {
  no.y <- object$no.y
    ## meta3L or meta3LML class
  if ( any(c("meta3L","meta3LFIML") %in% class(object)) ) {
    ## Tau2 with predictors
    Tau2_2model <- tryCatch( mxEvalByName("Tau2_2", object$mx.fit), error = function(e) NA )
    Tau2_3model <- tryCatch( mxEvalByName("Tau2_3", object$mx.fit), error = function(e) NA )
    Tau2_2base <- tryCatch( mxEvalByName("Tau2_2", object$mx0.fit$mx.fit), error = function(e) NA )
    Tau2_3base <- tryCatch( mxEvalByName("Tau2_3", object$mx0.fit$mx.fit), error = function(e) NA )           

    R2_2 <- max((1-Tau2_2model/Tau2_2base), 0)
    R2_3 <- max((1-Tau2_3model/Tau2_3base), 0)
    R2.values <- matrix(c(Tau2_2base, Tau2_2model, R2_2, Tau2_3base, Tau2_3model, R2_3), ncol=2)
    dimnames(R2.values) <- list(c("Tau2 (no predictor)", "Tau2 (with predictors)", "R2"),
                                c("Level 2", "Level 3"))
    } else {
    ## Return NA when Tau2 values cannot be retrieved, e.g., constrained at lbound or no name
    Tau2model <- sapply(1:no.y, function(x)
                        tryCatch( mxEvalByName(paste0("Tau2_",x,"_",x), object$mx.fit), error = function(e) NA ))
    Tau2base <- sapply(1:no.y, function(x)
                       tryCatch( mxEvalByName(paste0("Tau2_",x,"_",x), object$mx0.fit$mx.fit), error = function(e) NA ))
    R2_2 <- pmax((1-Tau2model/Tau2base), 0)
    R2.values <- rbind(Tau2base, Tau2model, R2_2)
    dimnames(R2.values) <- list( c("Tau2 (no predictor)", "Tau2 (with predictors)", "R2"),
                                 paste("y", 1:no.y, sep="") )       
    ## dimnames(R2.values) <- list( paste("y", c(t(outer(1:no.y, c(": Tau2 (no predictor)", ": Tau2 (with predictors)", ": R2"),
    ##                                                   paste, sep=""))), sep=""), c("Estimate")) 
    }
  
  R2.values
} 

## Calculate I2 for meta and meta3L objects
.I2 <- function(object, my.mx) {
  I2 <- object$I2
  no.y <- object$no.y
  if (missing(my.mx)) {
    my.mx <- summary(object$mx.fit)
  }

  ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
  my.ci <- my.mx$CI
  if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[,1:3, drop=FALSE]
      
  ## meta3L class
  if ( "meta3L" %in% class(object) ) {
    I2.names <- c(outer(I2, c("_2","_3"), paste, sep=""))
   
    ## Wald test, no CI
    if (is.null(dimnames(my.ci))) {
      I2.values <- matrix(NA, nrow=length(I2.names), ncol=3)
      I2.values[,2] <- unlist(lapply(I2.names, function(x) c(mxEvalByName(x, object$mx.fit))))
    } else {
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                           {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- my.ci[!is.na(name), , drop=FALSE]
      dimnames(my.ci)[[1]] <- name[!is.na(name)]
      I2.values <- my.ci
      ## I2.values <- my.ci[paste(I2.names, "[1,1]", sep=""), , drop=FALSE]
    }
    ## Truncate within 0 and 1
    I2.values[I2.values<0] <- 0
    I2.values[I2.values>1] <- 1    
    ## I2.values <- apply(I2.values, c(1,2), function(x) max(x,0))
    ## I2.values <- ifelse(I2.values<0, 0, I2.values)
    ## I2.values <- ifelse(I2.values>1, 1, I2.values) 
  
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
      I2.values[,2] <- mxEvalByName("I2_values", object$mx.fit)
    } else {## LB CI  model.name <- "Meta analysis with ML."
        ## modified to remove NA in row names
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                           {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- my.ci[!is.na(name), , drop=FALSE]
      dimnames(my.ci)[[1]] <- name[!is.na(name)]
      I2.values <- my.ci
      #I2.values <- my.ci[paste("I2_values[", 1:(length(I2.names)*no.y), ",1]", sep=""), ,drop=FALSE]
    }    
    ## Truncate into 0
    I2.values[I2.values<0] <- 0
    I2.values[I2.values>1] <- 1

    ## I2.values <- ifelse(I2.values<0, 0, I2.values)
    ## I2.values <- ifelse(I2.values>1, 1, I2.values)
    
    I2.names <- paste(": ", I2.names, sep="")
    I2.names <- paste("Intercept", c( outer(1:no.y, I2.names, paste, sep="")), sep="")
    I2.names <- sub("I2q", "I2 (Q statistic)", I2.names)
    I2.names <- sub("I2hm", "I2 (harmonic mean)", I2.names)
    I2.names <- sub("I2am", "I2 (arithmetic mean)", I2.names)
    dimnames(I2.values) <- list(I2.names, c("lbound", "Estimate", "ubound"))
    }
  
  I2.values
}

## Modified from rmseaConfidenceIntervalHelper() in OpenMx
.rmseaCI <- function(chi.squared, df, N, G=1, lower=.05, upper=.95){

    if (is.na(chi.squared) || is.na(df) || df < 0) {
        return(c(lower=NA_real_, upper=NA_real_))
    } else if (df==0) {
        return(c(lower=0, upper=0))
    } else {
        pChiSqFun <- function(x, val, degf, goal){ goal - pchisq(val, degf, ncp=x) }  

        # Lower confidence interval
        if (pchisq(chi.squared, df=df, ncp=0) >= upper) { #sic
            lower.lam <- tryCatch( uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=upper)$root,
                                   error=function(e) e )
            if (inherits(lower.lam, "error")) lower.lam <- NA
            # solve pchisq(ch, df=df, ncp=x) == upper for x            
        } else{
            lower.lam <- 0
        }

        # Upper confidence interval
        if (pchisq(chi.squared, df=df, ncp=0) >= lower) { #sic
            upper.lam <- tryCatch( uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=lower)$root,
                                   error=function(e) e )
            if (inherits(upper.lam, "error")) upper.lam <- NA
            # solve pchisq(ch, df=df, ncp=x) == lower for x
        } else{
            upper.lam <- 0
        }

        ## Normal-theory maximum likelihood convention, matching the point estimate.
        lower.rmsea <- sqrt(G)*sqrt(max((lower.lam)/N/df, 0))
        upper.rmsea <- sqrt(G)*sqrt(max((upper.lam)/N/df, 0))

        return(c(lower=lower.rmsea, upper=upper.rmsea))
    }
}

## Fit indices of the target model shared by summary.wls() and summary.tssem1FEM()
## df=0 follows lavaan: p=NA, RMSEA=0 with CI (0,0), TLI=1, CFI by its formula.
## df<0 (more free parameters than observed statistics): all indices are NA.
.fitIndices <- function(tT, dfT, tB, dfB, N, G=1) {
    na.out <- list(p=NA_real_, RMSEA=NA_real_, RMSEA.CI=c(lower=NA_real_, upper=NA_real_),
                   TLI=NA_real_, CFI=NA_real_)
    if (is.na(tT) || is.na(dfT)) return(na.out)
    if (dfT < 0) {
        warning("The model has more free parameters than observed statistics (df = ", dfT,
                "). It is not identified, and the fit indices are not computed.\n")
        return(na.out)
    }
    if (dfT == 0) {
        p <- NA_real_
        RMSEA <- 0
        TLI <- 1
    } else {
        p <- pchisq(tT, df=dfT, lower.tail=FALSE)
        RMSEA <- sqrt(G)*sqrt(max((tT-dfT)/N,0)/dfT)
        ## Hu and Bentler (1998)
        TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    }
    RMSEA.CI <- .rmseaCI(chi.squared=tT, df=dfT, N=N, G=G, lower=.05, upper=.95)
    numerator <- max(tT-dfT, 0)
    denominator <- max(tT-dfT, tB-dfB, 0)
    CFI <- if (isTRUE(denominator == 0)) 1 else 1 - numerator/denominator
    list(p=p, RMSEA=RMSEA, RMSEA.CI=RMSEA.CI, TLI=TLI, CFI=CFI)
}

## Use chol2inv() first and ginv() if there are problems
.solve <- function(x, parameters) {
    var.names <- colnames(x)
    
    acov <- tryCatch( 2*chol2inv(chol(x)), error = function(e) e )
    # Issue a warning instead of error message
    if (inherits(acov, "error")) {
      warning("Error in solving the Hessian matrix. Generalized inverse is used. The standard errors may not be trustworthy.\n")
      acov <- 2*MASS::ginv(x)
    }    
    dimnames(acov) <- list(var.names, var.names)
    
    if (!missing(parameters)) {
        acov <- acov[parameters, parameters, drop=FALSE]
    }    
    
    acov    
}

## Generate names for OSMASEM
## x: a vector of observed variables
## output: a vector of correlation coefficients and sampling covariance matrix
## This function is used to select data in data created by Cor2DataFrame()
.genCorNames <- function(x, cor.analysis=TRUE) {
    ## Names for the correlation elements
    if (cor.analysis) {
        psNames <- vechs(outer(x, x, paste, sep = "_"))
    } else {
        psNames <- vech(outer(x, x, paste, sep = "_"))
    }
    ## Distinct pairs can have the same concatenation when variable names
    ## themselves contain underscores. Keep historical names where unique.
    psNames <- make.unique(as.character(psNames), sep="_dup")
    
    ## Names for the sampling covariance matrix of the correlation vector
    psCovNames <- paste("C(", outer(psNames, psNames, paste, sep = " "), ")", sep="")
    psCovNames <- make.unique(as.character(vech(matrix(psCovNames, nrow=length(psNames), ncol=length(psNames)))), sep="_dup")

    list(ylabels=psNames, vlabels=psCovNames)
}
## OpenMx's default "Feasibility tolerance" (5e-2) lets constraints be
## violated by up to 0.05 at the solution, which can change the estimates of
## models with active inequality constraints. Models with mxConstraints are
## fitted with a tight tolerance unless the user has set one, either globally
## with mxOption(NULL, ...) or on the model.
.feasibility.tolerance <- 1e-6

.has.constraints <- function(model) {
  length(model@constraints) > 0 ||
    any(vapply(model@submodels, .has.constraints, logical(1)))
}

.tighten.feasibility <- function(model) {
  if (!.has.constraints(model) ||
      !is.null(model@options[["Feasibility tolerance"]]) ||
      as.numeric(mxOption(NULL, "Feasibility tolerance")) != 0.05)
    return(model)
  mxOption(model, "Feasibility tolerance", .feasibility.tolerance)
}

## After fitting a model with mxConstraints: rerun once from the solution if
## the status code is neither 0 nor 1, and warn if a constraint is still
## violated by more than 10 times the feasibility tolerance in use.
.check.constraints <- function(fit, rerun=TRUE, ...) {
  if (!inherits(fit, "MxModel") || !.has.constraints(fit)) return(fit)
  if (rerun && !(fit@output$status$code %in% c(0, 1))) {
    refit <- tryCatch(mxRun(fit, silent=TRUE, suppressWarnings=TRUE, ...),
                      error=function(e) e)
    if (!inherits(refit, "error")) fit <- refit
  }
  tolerance <- fit@options[["Feasibility tolerance"]]
  if (is.null(tolerance)) tolerance <- mxOption(NULL, "Feasibility tolerance")
  violation <- max(abs(as.numeric(fit@output$constraintFunctionValues)), 0)
  if (is.finite(violation) && violation > 10*as.numeric(tolerance))
    warning("A constraint is violated by ", signif(violation, 3),
            " at the solution (OpenMx status code ", fit@output$status$code,
            "). The estimates may be wrong. You are advised to rerun the ",
            "model, e.g., with rerun().",
            call.=FALSE)
  fit
}

## mxTryHard() prints its progress to the console and its notes as messages,
## even with silent=TRUE. Internal fallbacks run it quietly; warnings and
## errors are kept.
.quiet <- function(expr) {
  value <- NULL
  utils::capture.output(value <- suppressMessages(expr))
  value
}
