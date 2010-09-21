tssem1 <- function(my.df, n, start.values, cor.analysis = TRUE,
                   suppressWarnings = TRUE, ...) {
    no.groups <- length(my.df)
    no.var <- max(sapply(my.df, ncol))
    var.names <- paste("x", 1:no.var, sep = "")
    # matrix of labels; only use the lower triangle
    ps.labels <- outer(1:no.var, 1:no.var, function(x, y) paste("C", x, y, sep = ""))
    total.n <- sum(n)
    
    ## Check positive definiteness of data
    isPD <- is.pd(my.df)
    if (!all(isPD)) 
        stop(paste("Group(s) ", (1:no.groups)[isPD], " are not positive definite."), sep = "")
    
    ## Prepare starting values
    if (missing(start.values)) 
        sv <- .startValues(my.df, cor.analysis = cor.analysis)
    
    ## Index for missing variables: only check the diagonals only!!!
    miss.index <- lapply(my.df, function(x) { is.na(diag(x)) })
    
    ## complete.index <- NULL
    ## for (i in no.groups:1) {
    ##     if (sum(!miss.index[[i]], na.rm = TRUE) == no.var) 
    ##         complete.index = i
    ## }
    ## if (is.null(complete.index)) 
    ##     stop("It is expected that at least one study has complete data!")
     
    if ( sum(!miss.index[[1]], na.rm = TRUE) != no.var )
        stop("There should be no missing data in the first group.")
    
    for (i in 1:no.groups) {
        no.var.i <- sum(!miss.index[[i]], na.rm = TRUE)
        
        # Prepare matrices for calculations
        if (cor.analysis) {
            model.name <- "TSSEM1 Analysis of Correlation Matrix"
            S.matrix <- paste("S", i, " <- mxMatrix('Stand', nrow=", no.var.i, ", ncol=", 
                no.var.i, ", free=TRUE, values=vechs(sv[!miss.index[[", i, "]],!miss.index[[", 
                i, "]]]), name=\"S", i, "\", labels=vechs(ps.labels[!miss.index[[", 
                i, "]],!miss.index[[", i, "]]]))", sep = "")
            D.matrix <- paste("D", i, " <- mxMatrix('Diag', nrow=", no.var.i, ", ncol=", 
                no.var.i, ", free=TRUE, values=sqrt(diag(my.df[[", i, "]])[!miss.index[[", 
                i, "]]]), labels=c(", paste("\"D", i, var.names[!miss.index[[i]]], 
                  "\"", sep = "", collapse = ","), "), name=\"D", i, "\")", sep = "")
            # Expected covariance matrix
            expC.algebra <- paste("expC", i, " <- mxAlgebra(D", i, " %&% S", i, ", name=\"expC", 
                i, "\", dimnames=list(var.names[!miss.index[[", i, "]]], var.names[!miss.index[[", 
                i, "]]]))", sep = "")
            # Create mxModel
            g.model <- paste("g", i, " <- mxModel(\"g", i, "\", S", i, ", D", i, 
                ", expC", i, ", mxData(observed=my.df[[", i, "]][!miss.index[[", 
                i, "]],!miss.index[[", i, "]]], type=\"cov\", numObs=n[", i, "]), mxMLObjective(\"expC", 
                i, "\", dimnames=var.names[!miss.index[[", i, "]]]))", sep = "")
            
            eval(parse(text = S.matrix))
            eval(parse(text = D.matrix))
            eval(parse(text = expC.algebra))
            eval(parse(text = g.model))
        } else {
            model.name <- "TSSEM1 Analysis of Covariance Matrix"
            S.matrix <- paste("S", i, " <- mxMatrix('Symm', nrow=", no.var.i, ", ncol=", 
                no.var.i, ", free=TRUE, values=vech(sv[!miss.index[[", i, "]],!miss.index[[", 
                i, "]]]), name=\"S", i, "\", labels=vech(ps.labels[!miss.index[[", 
                i, "]],!miss.index[[", i, "]]]))", sep = "")
            g.model <- paste("g", i, " <- mxModel(\"g", i, "\", S", i, ", mxData(observed=my.df[[", 
                i, "]][!miss.index[[", i, "]],!miss.index[[", i, "]]], type=\"cov\", numObs=n[", 
                i, "]), mxMLObjective(\"S", i, "\", dimnames=var.names[!miss.index[[", 
                i, "]]]))", sep = "")
            
            eval(parse(text = S.matrix))
            eval(parse(text = g.model))
        }
    }
    
    if (cor.analysis) {
        tssem1.model <- paste("tssem1 <- mxModel(\"", model.name, "\", ", paste("S", 
            1:no.groups, sep = "", collapse = ","), ", ", paste("D", 1:no.groups, 
            sep = "", collapse = ","), ", ", paste("expC", 1:no.groups, sep = "", 
            collapse = ","), ", ", paste("g", 1:no.groups, sep = "", collapse = ","), 
            ", mxAlgebra(", paste("g", 1:no.groups, ".objective", sep = "", collapse = "+"), 
            ", name=\"obj\"), mxAlgebraObjective(\"obj\"))", sep = "")
    } else {
        tssem1.model <- paste("tssem1 <- mxModel(\"", model.name, "\", ", paste("S", 
            1:no.groups, sep = "", collapse = ","), ", ", paste("g", 1:no.groups, 
            sep = "", collapse = ","), ", mxAlgebra(", paste("g", 1:no.groups, ".objective", 
            sep = "", collapse = "+"), ", name=\"obj\"), mxAlgebraObjective(\"obj\"))", 
            sep = "")
    }
    eval(parse(text = tssem1.model))
    
    # try to run it with error message as output
    tssem1.fit <- tryCatch(mxRun(tssem1, suppressWarnings = suppressWarnings, ...),
                           error = function(e) e)
    if (inherits(tssem1.fit, "error")) 
        stop(print(tssem1.fit))
    
    pooledS <- eval(parse(text = "mxEval(S1, tssem1.fit)"))
    
    if (cor.analysis) {
        #Hessian_S <- 0.5*tssem1.fit@output$calculatedHessian[vechs(ps.labels), vechs(ps.labels)]
        acovS <- tryCatch( 2*solve(tssem1.fit@output$calculatedHessian[vechs(ps.labels),
                           vechs(ps.labels)]), error = function(e) e) 
    } else {
        #Hessian_S <- 0.5*tssem1.fit@output$calculatedHessian[vech(ps.labels), vech(ps.labels)]
        acovS <-  tryCatch( 2*solve(tssem1.fit@output$calculatedHessian[vech(ps.labels), 
                            vech(ps.labels)]), error = function(e) e)
    }
    # Issue a warning instead of error message
    if (inherits(acovS, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acovS))
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

    independentMinus2LL <- tryCatch(sum(.minus2LL(x=my.df, n=n, model="independent")), error = function(e) e)
    saturatedMinus2LL <- tryCatch(sum(.minus2LL(x=my.df, n=n, model="saturated")), error = function(e) e)
    
    out <- list(call = match.call(), data=my.df, pooledS = pooledS, acovS = acovS, total.n = total.n, 
                modelMinus2LL = tssem1.fit@output$Minus2LogLikelihood,
                independentMinus2LL = independentMinus2LL, saturatedMinus2LL = saturatedMinus2LL,
                tssem1.fit = tssem1.fit)
    class(out) <- "tssem1"
    return(out)
}



wls <- function(S, acovS, n, impliedS, matrices, cor.analysis = TRUE,
                intervals.type =c("z", "LB"), suppressWarnings = TRUE, ...) {
    impliedS@name <- "impliedS"
    no.var <- ncol(S)
    sampleS <- mxMatrix("Full", ncol = no.var, nrow = no.var, values = c(S), free = FALSE, 
        name = "sampleS")
    
    intervals.type <- match.arg(intervals.type)
    # Default is z
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)
    
    if (cor.analysis) {
        model.name <- "Correlation structure"
        ps <- no.var * (no.var - 1)/2
        vecS <- mxAlgebra(vechs(sampleS - impliedS), name = "vecS")
    } else {
        model.name <- "Covariance structure"
        ps <- no.var * (no.var + 1)/2
        vecS <- mxAlgebra(vech(sampleS - impliedS), name = "vecS")
    }
    if (ncol(acovS) != ps) 
        stop("No. of dimension of \"S\" does not match the dimension of \"acovS\"\n")
    
    # Inverse of asymptotic covariance matrix
    invacovS <- tryCatch(solve(acovS), error = function(e) e)
    if (inherits(invacovS, "error")) {
        cat("Error in inverting \"acovS\":\n")
        stop(print(invacovS))
    }
    
    invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invacovS), free = FALSE, 
        name = "invAcov")
    obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")
    objective <- mxAlgebraObjective("obj")
    
    if (missing(matrices)) {
        text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", ", "impliedS", 
            ", sampleS, vecS, invAcov, obj, objective, mxCI(\"impliedS\")), intervals=", 
            intervals, ", suppressWarnings = ", suppressWarnings, ", ...)", sep = "")
    } else {
        matName1 <- sapply(matrices, function(x) {
            x@name
        })
        matName2 <- sapply(matName1, function(x) {
            paste("\"", x, "\"", sep = "")
        })
        
        text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", ", "impliedS", 
            ", sampleS, vecS, invAcov, obj, objective, ", paste(matName1, collapse = ", "), 
            ", mxCI(c(", paste(matName2, collapse = ", "), "))", "), intervals=", 
            intervals, ", suppressWarnings = ", suppressWarnings, ", ... )", sep = "")
    }
    
    wls.fit <- tryCatch(eval(parse(text = text1)), error = function(e) e)
    
    # try to run it with error message as output
    if (inherits(wls.fit, "error")) {
        cat("Error in running the mxModel:\n")
        stop(print(wls.fit))
    } else {
        out <- list(call = match.call(), noObservedStat=ps, n=n, 
                    indepModelChisq=.indepwlsChisq(S=S, acovS=acovS, cor.analysis=cor.analysis),
                    indepModelDf=no.var*(no.var-1)/2, wls.fit=wls.fit)
        class(out) <- 'wls'
    }
    out
}


tssem2 <- function(tssem1.obj, impliedS, matrices, intervals.type = c("z", "LB"),
                   suppressWarnings = TRUE, ...) {
  if (!is.element("tssem1", class(tssem1.obj)))
    stop("\"tssem1.obj\" must be an object of class \"tssem1\".")
  # check the call to determine whether it is a correlation or covariance analysis
  cor.analysis <- tssem1.obj$call[[match("cor.analysis", names(tssem1.obj$call))]]
  # if not specified, the default in tssem1() is cor.analysis=TRUE
  if (is.null(cor.analysis)) {
     cor.analysis <- TRUE
  } else {
     # to handle symbolic F(T) vs. logical FALSE(TRUE)
     cor.analysis <- as.logical(as.character(cor.analysis))
  }
  wls(S=tssem1.obj$pooledS, acovS=tssem1.obj$acovS, n=tssem1.obj$total.n, impliedS=impliedS,
      matrices=matrices, cor.analysis = cor.analysis, intervals.type = intervals.type,
      suppressWarnings = suppressWarnings, ...)
}
