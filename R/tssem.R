tssem1FEM <- function(my.df, n, cor.analysis=TRUE, model.name=NULL,
                     cluster=NULL, suppressWarnings=TRUE, silent=TRUE,
                     run=TRUE, ...) {
  if (!is.null(cluster)) {
    data.cluster <- tapply(my.df, cluster, function(x) {x})
    n.cluster <- tapply(n, cluster, function(x) {x})
    out <- list()
    for (i in 1:length(data.cluster)) {
      ## Need to correct it to tssem1()
      out[[i]] <- tssem1FEM(data.cluster[[i]], n.cluster[[i]], 
                           cor.analysis=cor.analysis, model.name=model.name,
                           suppressWarnings=suppressWarnings, ...)
    }
    names(out) <- names(data.cluster) 
    class(out) <- "tssem1FEM.cluster"
    out
  } else {
    ## Check whether all studies have the same dimensions  
    my.range <- range(sapply(my.df, function(x) {ncol(x)}))
    if ( !all.equal(my.range[1], my.range[2]) )
      stop("Dimensions of groups are not the same!\n")
    
    no.groups <- length(my.df)
    no.var <- ncol(my.df[[1]])

    ## Get the original variable names
    original.names <- colnames(my.df[[1]])
    
    var.names <- paste("x", 1:no.var, sep = "")
    ## Convert variable labels to x1, x2, ...
    my.df <- lapply(my.df, function(x) {dimnames(x) <- list(var.names, var.names); x})
    total.n <- sum(n)   
    
    ## Check positive definiteness of data
    isPD <- is.pd(my.df)
    if (!all(isPD)) 
        stop(paste("Group ", (1:no.groups)[!isPD], " is not positive definite.", sep = ""))
    
    ## Prepare starting values based on covariance matrices
    sv <- .startValues(my.df, cor.analysis = FALSE)
    # matrix of labels; only use the lower triangle
    ps.labels <- outer(1:no.var, 1:no.var, function(x, y) paste("s", x, y, sep = ""))
    if (cor.analysis==TRUE) {
      S <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=TRUE, values=vechs(cov2cor(sv)),
                    labels=vechs(ps.labels), name="S")
    } else {
      ## Set lower bound on variances
      lbound <- matrix(NA, nrow=no.var, ncol=no.var)
      Diag(lbound) <- 0.00001
      S <- mxMatrix(type="Symm", nrow=no.var, ncol=no.var, free=TRUE, values=vech(sv),
                    labels=vech(ps.labels), lbound=lbound, name="S")      
    }
       
    ## Index for missing variables: only check the diagonals only!!!
    miss.index <- lapply(my.df, function(x) { is.na(Diag(x)) })
    
    ## complete.index <- NULL
    ## for (i in no.groups:1) {
    ##     if (sum(!miss.index[[i]], na.rm = TRUE) == no.var) 
    ##         complete.index = i
    ## }
    ## if (is.null(complete.index)) 
    ##     stop("It is expected that at least one study has complete data!")
     
    ## if ( sum(!miss.index[[1]], na.rm = TRUE) != no.var )
    ##     stop("There should be no missing data in the first group.")
    
    for (i in 1:no.groups) {
        no.var.i <- sum(!miss.index[[i]], na.rm = TRUE)
        miss.index.i <- miss.index[[i]]
        my.df.i <- my.df[[i]][!miss.index.i, !miss.index.i]
        ## dimnames(my.df.i) <- list(var.names[!miss.index.i], var.names[!miss.index.i])
        
        # Prepare matrices for calculations
        if (cor.analysis) {
            if (is.null(model.name)) model.name <- "TSSEM1 Analysis of Correlation Matrix"
            SD <- Diag(paste(sqrt(Diag(sv)), "*sd", i, "_", 1:no.var, sep=""))
            # Fixed a bug when SD is a scalar
            SD <- SD[!miss.index.i, , drop=FALSE]
            SD <- as.mxMatrix(SD)   
        } else {
            if (is.null(model.name)) model.name <- "TSSEM1 Analysis of Covariance Matrix"
            SD <- Diag(rep(1, no.var))
            SD <- SD[!miss.index.i, , drop=FALSE]
            SD <- as.mxMatrix(SD)
        }
        # Expected covariance matrix
        expC <- mxAlgebra(SD %&% S, name="expC", dimnames=list(var.names[!miss.index.i],
                          var.names[!miss.index.i]))
        # Create mxModel
        ## model <- mxModel(paste("group",i,sep=""), S, Dmatrix, expC,
        ##                      mxData(observed=my.df.i, type="cov", numObs=n[i]),
        ##                      mxMLObjective(expC, dimnames=var.names[!miss.index[[i]]]))

        fitFunction <- mxFitFunctionML()
        
        ## Fixed a bug when my.df[[i]][, , ] is a scalar
        g.model <- paste("g", i, " <- mxModel(\"g", i, "\", S, SD, expC, mxData(observed=my.df[[",i,"]][!miss.index[[",i,"]],!miss.index[[",i,"]], drop=FALSE], type=\"cov\", numObs=n[", i,
                "]), fitFunction, mxExpectationNormal(covariance=\"expC\", means=NA, dimnames=var.names[!miss.index[[", i, "]]]))", sep = "")
        eval(parse(text = g.model))       
    }
    
    ## if (cor.analysis==TRUE) {
    ##     tssem1.model <- paste("tssem1 <- mxModel(\"", model.name, "\", ", paste("S", 
    ##         1:no.groups, sep = "", collapse = ","), ", ", paste("D", 1:no.groups, 
    ##         sep = "", collapse = ","), ", ", paste("expC", 1:no.groups, sep = "", 
    ##         collapse = ","), ", ", paste("g", 1:no.groups, sep = "", collapse = ","), 
    ##         ", mxAlgebra(", paste("g", 1:no.groups, ".objective", sep = "", collapse = "+"), 
    ##         ", name=\"obj\"), mxAlgebraObjective(\"obj\"))", sep = "")
    ## } else {
    ##     tssem1.modelb <- paste("tssem1 <- mxModel(\"", model.name, "\", ", paste("S", 
    ##         1:no.groups, sep = "", collapse = ","), ", ", paste("g", 1:no.groups, 
    ##         sep = "", collapse = ","), ", mxAlgebra(", paste("g", 1:no.groups, ".objective", 
    ##         sep = "", collapse = "+"), ", name=\"obj\"), mxAlgebraObjective(\"obj\"))", 
    ##         sep = "")
    ## }
    tssem1.model <- paste("tssem1 <- mxModel(\"", model.name, "\", S, ",
                          paste("g", 1:no.groups, sep = "", collapse = ","),
                          ", mxAlgebra(", paste("g", 1:no.groups, ".objective", sep = "", collapse = "+"),
                          ", name=\"obj\"), mxFitFunctionAlgebra(algebra =\"obj\"))", sep = "")
    eval(parse(text = tssem1.model))

    ## Return mx model without running the analysis
    if (run==FALSE) return(tssem1)
        
    # try to run it with error message as output
    ## mx.fit <- mxRun(tssem1)
    mx.fit <- tryCatch(mxRun(tssem1, suppressWarnings = suppressWarnings, silent=silent, ...),
                           error = function(e) e)
    if (inherits(mx.fit, "error")) 
        warning(print(mx.fit))
    
    pooledS <- eval(parse(text = "mxEval(S, mx.fit)"))

    # Fixed a bug that all elements have to be inverted before selecting some of them
    if (cor.analysis) {
        #Hessian_S <- 0.5*mx.fit@output$calculatedHessian[vechs(ps.labels), vechs(ps.labels)]
        acovS <- tryCatch( 2*solve(mx.fit@output$calculatedHessian)[vechs(ps.labels),
                           vechs(ps.labels)], error = function(e) e) 
    } else {
        #Hessian_S <- 0.5*mx.fit@output$calculatedHessian[vech(ps.labels), vech(ps.labels)]
        acovS <-  tryCatch( 2*solve(mx.fit@output$calculatedHessian)[vech(ps.labels), 
                            vech(ps.labels)], error = function(e) e)
    }
    # Issue a warning instead of error message
    if (inherits(acovS, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acovS))
    } else {
      # Fixed a bug in a few lines later in dimnames(acovS) when acovS is a scalar
      acovS <- as.matrix(acovS)
    }
    
    # check dimnames
    if (is.null(dimnames(my.df[[1]]))) {
        dimnames(pooledS) <- list(var.names, var.names)
    } else {
        df.dim <- dimnames(my.df[[1]])
        dimnames(pooledS) <- df.dim
        # create matrix of labels for ps
        if (cor.analysis) {
            psMatnames <- vechs(outer(unlist(df.dim[1]), unlist(df.dim[1]), paste, 
                sep = " "))
        } else {
            psMatnames <- vech(outer(unlist(df.dim[1]), unlist(df.dim[1]), paste, 
                sep = " "))
        }
        #dimnames(Hessian_S) <- list(psMatnames, psMatnames)
        dimnames(acovS) <- list(psMatnames, psMatnames)
    }

    ## Calculate 2LL of the saturated and independence models and the DF of independence model
    baseMinus2LL <- tryCatch(.minus2LL(x=my.df, n=n), error = function(e) e)   

    out <- list(call = match.call(), cor.analysis=cor.analysis, data=my.df, pooledS = pooledS,
                acovS = acovS, n = n, 
                modelMinus2LL = mx.fit@output$Minus2LogLikelihood,
                baseMinus2LL = baseMinus2LL, mx.model=tssem1, mx.fit = mx.fit,
                original.names=original.names)
    class(out) <- "tssem1FEM"
    return(out)
  }
}

tssem1REM <- function(my.df, n, cor.analysis=TRUE, RE.type=c("Symm", "Diag", "Zero"), RE.startvalues=0.1, RE.lbound = 1e-10,
                      I2="I2q", model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  ## It handles missing effect sizes rather than missing correlations. Thus, it is more flexible than tssem1FEM().
  ## ACOV is calculated without missing data by assuming 1 and 0 for the missing variances and covariances.
  ## Missing values are indicated by the missing effect sizes.

  ## Get the original variable names
  original.names <- colnames(my.df[[1]])

  if (cor.analysis) {
      ## Replace diagonals with 1.0
      my.complete <- lapply(my.df, function (x) { Diag(x)[is.na(Diag(x))] <- 1; x })
  } else {
      ## Replace diagonals with the mean of diagonals
      my.complete <- lapply(my.df, function (x) { Diag(x)[is.na(Diag(x))] <- mean(Diag(x), na.rm=TRUE); x })
  }      
  ## Replace missing variables with 0.0 regardless of cor.analysis
  my.complete <- lapply(my.complete, function (x) { x[is.na(x)] <- 0; x })
  
  ## Calculate the asymptotic sampling covariance matrix of the correlation matrix
  acovR <- asyCov(x=my.complete, n=n, cor.analysis=cor.analysis)

  ## Fixed a bug that my.df is covariance matrix while cor.analysis is TRUE
  ## When cor.analysis=TRUE, the old version just takes the lower triangle without converting covariance into correlation.
  if (cor.analysis) {
    ## Convert possible covariance matrices into correlation matrices
    ## When there are NA in diagonas, they become 1 after cov2cor()
    ## It is fine as the diagonals are not used in cor.analysis=TRUE
    ES <- list2matrix(x=suppressWarnings(lapply(my.df, cov2cor)), diag=FALSE)
  } else {
    ES <- list2matrix(x=my.df, diag=TRUE)
  } 
  ## no. of effect sizes
  no.es <- ncol(ES)

  if (is.null(model.name)) {
    if (cor.analysis) {
      model.name <- "TSSEM1 (Random Effects Model) Analysis of Correlation Matrix"
    } else {
      model.name <- "TSSEM1 (Random Effects Model) Analysis of Covariance Matrix"
    }
  }

  RE.type <- match.arg(RE.type)
  switch( RE.type,
         Symm = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.startvalues=RE.startvalues,
                               RE.lbound=RE.lbound, suppressWarnings=TRUE, silent=silent, run=run, ...),
## Prior to R-3.0.0
##       Diag = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2,
##                             RE.constraints=Diag(x=paste(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es, sep=""),
##                                            nrow=no.es, ncol=no.es), RE.lbound=RE.lbound),
         Diag = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2,
                               RE.constraints=Diag(x=paste(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es, sep="")),
                               RE.lbound=RE.lbound, suppressWarnings=TRUE, silent=silent, run=run, ...),
         Zero = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.constraints=matrix(0, ncol=no.es, nrow=no.es),
                               suppressWarnings=TRUE, silent=silent, run=run, ...) )

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.fit)
  
  ## if (RE.diag.only==TRUE) {
  ##   ## No covariance between random effects
  ##   mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2,
  ##                    RE.constraints=diag(x=paste(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es, sep=""),
  ##                                        nrow=no.es, ncol=no.es), RE.lbound=RE.lbound)    
  ## } else {
  ##   mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.startvalues=RE.startvalues, RE.lbound=RE.lbound)
  ## }

  out <- list(total.n=sum(n), cor.analysis=cor.analysis, RE.type=RE.type, no.es=no.es, original.names=original.names)
  out <- c(out, mx.fit)
  class(out) <- c("tssem1REM", "meta")
  return(out)
}

tssem1 <- function(my.df, n, method=c("FEM", "REM"), cor.analysis=TRUE, cluster=NULL,
                   RE.type=c("Symm", "Diag", "Zero"), RE.startvalues=0.1, RE.lbound=1e-10, I2="I2q",
                   model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  method <- match.arg(method)
  switch(method,
    FEM = out <- tssem1FEM(my.df=my.df, n=n, cor.analysis=cor.analysis, model.name=model.name,
                          cluster=cluster, suppressWarnings=suppressWarnings, silent=silent, run=run, ...),
    REM = out <- tssem1REM(my.df=my.df, n=n, cor.analysis=cor.analysis, RE.type=RE.type,
                          RE.startvalues=RE.startvalues, RE.lbound=RE.lbound, I2=I2,
                          model.name=model.name, suppressWarnings=suppressWarnings, silent=silent, run=run, ...) )
  out  
}
  
## wls <- function(S, acovS, n, impliedS, matrices, cor.analysis = TRUE,
##                 intervals.type =c("z", "LB"), model.name, suppressWarnings = TRUE, ...) {
##     impliedS@name <- "impliedS"
##     no.var <- ncol(S)
##     sampleS <- mxMatrix("Full", ncol = no.var, nrow = no.var, values = c(S), free = FALSE, 
##         name = "sampleS")
    
##     intervals.type <- match.arg(intervals.type)
##     # Default is z
##     switch(intervals.type,
##            z = intervals <- FALSE,
##            LB = intervals <- TRUE)
    
##     if (cor.analysis) {
##         if (missing(model.name)) model.name <- "WLS Analysis of Correlation Structure"
##         ps <- no.var * (no.var - 1)/2
##         vecS <- mxAlgebra(vechs(sampleS - impliedS), name = "vecS")
##     } else {
##         if (missing(model.name)) model.name <- "WLS Analysis of Covariance Structure"
##         ps <- no.var * (no.var + 1)/2
##         vecS <- mxAlgebra(vech(sampleS - impliedS), name = "vecS")
##     }
##     if (ncol(acovS) != ps) 
##         stop("No. of dimension of \"S\" does not match the dimension of \"acovS\"\n")
    
##     # Inverse of asymptotic covariance matrix
##     invacovS <- tryCatch(solve(acovS), error = function(e) e)
##     if (inherits(invacovS, "error")) {
##         cat("Error in inverting \"acovS\":\n")
##         stop(print(invacovS))
##     }
    
##     invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invacovS), free = FALSE, 
##         name = "invAcov")
##     obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")
##     objective <- mxAlgebraObjective("obj")
    
##     if (missing(matrices)) {
##         text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", ", "impliedS", 
##             ", sampleS, vecS, invAcov, obj, objective, mxCI(\"impliedS\")), intervals=", 
##             intervals, ", suppressWarnings = ", suppressWarnings, ", ...)", sep = "")
##     } else {
##         matName1 <- sapply(matrices, function(x) {
##             x@name
##         })
##         matName2 <- sapply(matName1, function(x) {
##             paste("\"", x, "\"", sep = "")
##         })
        
##         text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", ", "impliedS", 
##             ", sampleS, vecS, invAcov, obj, objective, ", paste(matName1, collapse = ", "), 
##             ", mxCI(c(", paste(matName2, collapse = ", "), "))", "), intervals=", 
##             intervals, ", suppressWarnings = ", suppressWarnings, ", ... )", sep = "")
##     }
    
##     mx.fit <- tryCatch(eval(parse(text = text1)), error = function(e) e)
    
##     # try to run it with error message as output
##     if (inherits(mx.fit, "error")) {
##         cat("Error in running the mxModel:\n")
##         stop(print(mx.fit))
##     } else {
##         out <- list(call = match.call(), noObservedStat=ps, n=n, cor.analysis=cor.analysis,
##                     indepModelChisq=.indepwlsChisq(S=S, acovS=acovS, cor.analysis=cor.analysis),
##                     indepModelDf=no.var*(no.var-1)/2, mx.fit=mx.fit)
##         class(out) <- 'wls'
##     }
##     out
## }

wls <- function(Cov, asyCov, n, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, diag.constraints=FALSE,
                cor.analysis=TRUE, intervals.type=c("z", "LB"), mx.algebras=NULL, model.name=NULL,
                suppressWarnings=TRUE, silent=TRUE, run=TRUE,...) {
  if (is.null(Smatrix)) {
    stop("\"Smatrix\" matrix is not specified.\n")
  } else {
    if (is.matrix(Smatrix)) Smatrix <- as.mxMatrix(Smatrix)
    p <- nrow(Smatrix@values)
    Smatrix@name <- "Smatrix"  
  }

  if (is.null(Amatrix)) {
    Amatrix <- as.mxMatrix(matrix(0, nrow=p, ncol=p), name="Amatrix")
  } else {
    if (is.matrix(Amatrix)) Amatrix <- as.mxMatrix(Amatrix)  
    Amatrix@name <- "Amatrix"
  }

  if (is.null(Fmatrix)) {
    Fmatrix <- as.mxMatrix(Diag(rep(p,1)), name="Fmatrix")
  } else {
    if (is.matrix(Fmatrix)) Fmatrix <- as.mxMatrix(Fmatrix)  
    Fmatrix@name <- "Fmatrix"
  }

  Id <- as.mxMatrix(Diag(rep(p,1)), name="Id")

  ## No. of observed variables
  no.var <- ncol(Cov)
  if (is.pd(Cov)) {
    sampleS <- as.mxMatrix(Cov, name="sampleS")
  } else {
    stop("\"Cov\" is not positive definite.\n")
  }
    
  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
         z = intervals <- FALSE,
         LB = intervals <- TRUE)

  # Inverse of asymptotic covariance matrix
  if (is.pd(asyCov)) {
    invacovS <- tryCatch(solve(asyCov), error = function(e) e)
  } else {
    stop("\"asyCov\" is not positive definite.\n")
  }

  ## It appears that solve() does not fail
  if (inherits(invacovS, "error")) {
    cat("Error in inverting \"asyCov\":\n")
    stop(print(invacovS))
  }
  invAcov <- as.mxMatrix(invacovS, name="invAcov")
  impliedS1 <- mxAlgebra( (Fmatrix%*%solve(Id-Amatrix))%&%Smatrix, name="impliedS1" )

  ## Assuming no constraint
  Constraints <- 0
  if (cor.analysis) {
    if (is.null(model.name)) model.name <- "WLS Analysis of Correlation Structure"
    ps <- no.var * (no.var - 1)/2
    
    ## Count no. of dependent variables including both observed and latent variables
    ## Since it is correlation structure, Smatrix@values=1 and Smatrix@free=FALSE on the diagonals.
    Constraints <- Diag(Smatrix@free)
    
    ## Setup nonlinear constraints on diagonals
    if (diag.constraints & (sum(Constraints)>0)) {      
      One <- mxMatrix("Full", values=1, ncol=1, nrow=sum(Constraints), free=FALSE, name="One")
      select <- create.Fmatrix(Constraints, name="select")
      constraint <- mxConstraint( select%*%diag2vec(solve(Id-Amatrix)%&%Smatrix)==One, name="constraint" )
      impliedS <- mxAlgebra(impliedS1, name="impliedS")
    } else {
    ## Use 1-impliedS for error variances  
      Diag(Smatrix@free) <- FALSE
      Diag(Smatrix@values)[Constraints] <- 0
      Diag(Smatrix@labels)[Constraints] <- NA
       
    ## Error variances are computed rather than estimated
    ## Ematrix = 1 - diag(impliedS)
    ## If diag(Smatrix) are fixed at 1, Ematrix are 0.
      Ematrix <- mxAlgebra(Fmatrix%&%Id - vec2diag(diag2vec(impliedS1)), name="Ematrix")
      impliedS <- mxAlgebra( impliedS1 + Ematrix, name="impliedS")
    }
    vecS <- mxAlgebra(vechs(sampleS - impliedS), name="vecS")
    
  } else {
    if (is.null(model.name)) model.name <- "WLS Analysis of Covariance Structure"
    ps <- no.var * (no.var + 1)/2
    impliedS <- mxAlgebra(impliedS1, name="impliedS")
    vecS <- mxAlgebra(vech(sampleS - impliedS), name="vecS")
  }  
  
  if (ncol(asyCov) != ps) 
    stop("No. of dimension of \"Cov\" does not match the multiplier of the dimension of \"asyCov\"\n")
    
  obj <- mxAlgebra( t(vecS) %&% invAcov, name = "obj" )
  objective <- mxFitFunctionAlgebra(algebra="obj")

  mx.model <- mxModel(model=model.name, Fmatrix, Amatrix, Smatrix, Id, impliedS1, impliedS,
                      vecS, invAcov, obj, objective, sampleS, mxCI(c("Amatrix", "Smatrix")))

  ## Add additional mxAlgebras
  if (!is.null(mx.algebras)) {
    for (i in 1:length(mx.algebras)) {
      mx.model <- mxModel(mx.model, mx.algebras[[i]])
    }
    mx.model <- mxModel(mx.model, mxCI(names(mx.algebras)))
  }
  ## mx.model <- eval(parse(text = text1))

  ## Add constraints on diagonals
  if (cor.analysis) {
    if (diag.constraints) {
      mx.model <- mxModel(mx.model, One, select, constraint)
    } else {
      mx.model <- mxModel(mx.model, Ematrix, mxCI(c("Amatrix")))
    }
  }
  
  ## if (diag.constraints) {
  ##   text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", Fmatrix, Amatrix, Smatrix, Id, impliedS, vecS, invAcov, ",
  ##                  "obj, objective, One, select, constraint, sampleS, mxCI(c(\"Amatrix\", \"Smatrix\"))), intervals=", 
  ##                   intervals, ", suppressWarnings = ", suppressWarnings, ", ...)", sep="")
  ## } else {
  ##   text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", Fmatrix, Amatrix, Smatrix, Id, impliedS, vecS, invAcov, ",
  ##                  "obj, objective, sampleS, mxCI(c(\"Amatrix\", \"Smatrix\"))), intervals=", 
  ##                   intervals, ", suppressWarnings = ", suppressWarnings, ", ...)", sep="")
  ## }

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.model)
  
  mx.fit <- tryCatch( mxRun(mx.model, intervals=intervals, suppressWarnings=suppressWarnings, silent=silent, ...), error = function(e) e)

  # try to run it with error message as output
  if (inherits(mx.fit, "error")) {
      cat("Error in running the mxModel:\n")
      warning(print(mx.fit))
  } else {
      out <- list(call=match.call(), Cov=Cov, asyCov=asyCov, noObservedStat=ps, n=n, cor.analysis=cor.analysis, diag.constraints=diag.constraints, Constraints=Constraints,
                  indepModelChisq=.indepwlsChisq(S=Cov, acovS=asyCov, cor.analysis=cor.analysis),
                  indepModelDf=no.var*(no.var-1)/2, mx.model=mx.model, mx.fit=mx.fit, mx.algebras=names(mx.algebras), intervals.type=intervals.type)
      class(out) <- 'wls'
  }
  out
}


tssem2 <- function(tssem1.obj, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, diag.constraints=FALSE,
                   intervals.type = c("z", "LB"), mx.algebras=NULL, model.name=NULL, suppressWarnings=TRUE,
                   silent=TRUE, run=TRUE, ...) {
  if ( !is.element( class(tssem1.obj)[1], c("tssem1FEM.cluster", "tssem1FEM", "tssem1REM")) )
      stop("\"tssem1.obj\" must be of neither class \"tssem1FEM.cluster\", class \"tssem1FEM\" or \"tssem1REM\".")
      
  switch(class(tssem1.obj)[1],
         tssem1FEM.cluster = { out <- lapply(tssem1.obj, tssem2, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
                                             diag.constraints=diag.constraints, intervals.type=intervals.type,
                                             mx.algebras=mx.algebras, 
                                             model.name=model.name, suppressWarnings=suppressWarnings, silent=silent, run=run, ...)
                              class(out) <- "wls.cluster" },
         tssem1FEM = { cor.analysis <- tssem1.obj$cor.analysis
                      # check the call to determine whether it is a correlation or covariance analysis
                      ## cor.analysis <- tssem1.obj$call[[match("cor.analysis", names(tssem1.obj$call))]]
                      ## # if not specified, the default in tssem1() is cor.analysis=TRUE
                      if (cor.analysis==TRUE) {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Fixed Effects Model) Analysis of Correlation Structure"
                      } else {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Fixed Effects Model) Analysis of Covariance Structure"
                      # to handle symbolic F(T) vs. logical FALSE(TRUE)
                      ## cor.analysis <- as.logical(as.character(cor.analysis))
                      }
                      ## Use the original varible names for the observed covariance matrix 
                      pooledS <- tssem1.obj$pooledS
                      dimnames(pooledS) <- list(tssem1.obj$original.names, tssem1.obj$original.names)
                      out <- wls(Cov=pooledS, asyCov=tssem1.obj$acovS, n=sum(tssem1.obj$n),
                                 Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
                                 diag.constraints=diag.constraints, cor.analysis=cor.analysis,
                                 intervals.type=intervals.type, mx.algebras=mx.algebras,
                                 model.name=model.name, suppressWarnings=suppressWarnings,
                                 silent=silent, run=run, ...) },
         tssem1REM = { cor.analysis <- tssem1.obj$cor.analysis
                      ## Extract the pooled correlation matrix
                      pooledS <- vec2symMat( coef(tssem1.obj, select="fixed"), diag=!cor.analysis)
                      ## Extract the asymptotic covariance matrix of the pooled correlations
                      asyCov <- vcov(tssem1.obj, select="fixed")
                      
                      if (cor.analysis==TRUE) {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Random Effects Model) Analysis of Correlation Structure"
                      } else {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Random Effects Model) Analysis of Covariance Structure"
                      }
                      ## Use the original varible names for the observed covariance matrix
                      dimnames(pooledS) <- list(tssem1.obj$original.names, tssem1.obj$original.names)                       
                      out <- wls(Cov=pooledS, asyCov=asyCov, n=tssem1.obj$total.n,
                                 Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix, diag.constraints=diag.constraints,
                                 cor.analysis=cor.analysis, intervals.type=intervals.type, mx.algebras=mx.algebras,
                                 model.name=model.name, suppressWarnings = suppressWarnings, silent=silent, run=run, ...) })
  out
}
    