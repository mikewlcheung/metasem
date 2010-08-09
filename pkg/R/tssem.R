Minus2LL <- function(x, n, type=c("saturated", "independent")) {
  if (is.list(x)) {
    if (length(x) != length(n))
      stop("Lengths of \"x\" and \"n\" are not the same.\n")
    return(mapply(Minus2LL, x=x, n=n, type=type))
  } else {
    miss.index <- is.na(diag(x))
    x <- x[!miss.index, !miss.index]
    if (!is.pd(x))
      stop("\"x\" is not positive definite.\n")
    no.var <- ncol(x)
    vars <- paste("v", 1:no.var, sep="")
    dimnames(x) <- list(vars, vars)
    obsCov <- mxData(observed=x, type='cov', numObs=n)
    if (missing(type))
      stop("\"type\" was not specified.\n")
    type <- match.arg(type)
    switch(type,
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


tssem1 <- function(my.df, n, start.values, cor.analysis = TRUE, ...) {
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
        sv <- startValues(my.df, cor.analysis = cor.analysis)
    
    ## Index for missing variables: only check the diagonals only!!!
    miss.index <- lapply(my.df, function(x) { is.na(diag(x)) })
    
    complete.index <- NULL
    for (i in no.groups:1) {
        if (sum(!miss.index[[i]], na.rm = TRUE) == no.var) 
            complete.index = i
    }
    if (is.null(complete.index)) 
        stop("It is expected that at least one study has complete data!")
    
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
    tssem1.fit <- tryCatch(mxRun(tssem1, ...), error = function(e) e)
    if (inherits(tssem1.fit, "error")) 
        stop(print(tssem1.fit))
    
    pooledS <- eval(parse(text = paste("mxEval(S", complete.index, ",tssem1.fit)", 
        sep = "")))
    
    if (cor.analysis) {
        #Hessian_S <- 0.5*tssem1.fit@output$calculatedHessian[vechs(ps.labels), vechs(ps.labels)]
        acovS <- 2 * solve(tssem1.fit@output$calculatedHessian[vechs(ps.labels), 
            vechs(ps.labels)])
    } else {
        #Hessian_S <- 0.5*tssem1.fit@output$calculatedHessian[vech(ps.labels), vech(ps.labels)]
        acovS <- 2 * solve(tssem1.fit@output$calculatedHessian[vech(ps.labels), 
            vech(ps.labels)])
    }
    
    # check dimnames
    if (is.null(dimnames(my.df[[complete.index]]))) {
        dimnames(pooledS) <- list(var.names, var.names)
    } else {
        df.dim <- dimnames(my.df[[complete.index]])
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

    independentMinus2LL <- tryCatch(sum(Minus2LL(x=my.df, n=n, type="independent")), error = function(e) e)
    saturatedMinus2LL <- tryCatch(sum(Minus2LL(x=my.df, n=n, type="saturated")), error = function(e) e)
    
    out <- list(pooledS = pooledS, acovS = acovS, total.n = total.n, cor.analysis = cor.analysis,
                modelMinus2LL = tssem1.fit@output$Minus2LogLikelihood,
                independentMinus2LL = independentMinus2LL, saturatedMinus2LL = saturatedMinus2LL,
                tssem1.fit = tssem1.fit)
    class(out) <- "tssem1"
    return(out)
}



wls <- function(S, acovS, n, impliedS, matrices, cor.analysis = TRUE, intervals = FALSE, ...) {
    impliedS@name <- "impliedS"
    no.var <- ncol(S)
    sampleS <- mxMatrix("Full", ncol = no.var, nrow = no.var, values = c(S), free = FALSE, 
        name = "sampleS")
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
    if (inherits(invacovS, "error")) 
        stop(print(invacovS))
    
    invAcov <- mxMatrix("Full", ncol = ps, nrow = ps, values = c(invacovS), free = FALSE, 
        name = "invAcov")
    obj <- mxAlgebra(t(vecS) %&% invAcov, name = "obj")
    objective <- mxAlgebraObjective("obj")
    
    if (missing(matrices)) {
        text1 <- paste("mxRun(mxModel(model=\"", model.name, "\", ", "impliedS", 
            ", sampleS, vecS, invAcov, obj, objective, mxCI(\"impliedS\")), intervals=", 
            intervals, ",...)", sep = "")
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
            intervals, ",... )", sep = "")
    }
    
    out <- tryCatch(eval(parse(text = text1)), error = function(e) e)
    
    # try to run it with error message as output
    if (inherits(out, "error")) {
        stop(print(out))
    }
    #else {
    #    class(wls.fit) <- 'wls'
    #}
    out
}


tssem2 <- function(tssem1.obj, impliedS, matrices, intervals = FALSE, ...) {
  if (!is.element("tssem1", class(tssem1.obj)))
    stop("\"tssem1.obj\" must be an object of class \"tssem1\".")
  wls(S=tssem1.obj$pooledS, acovS=tssem1.obj$acovS, n=tssem1.obj$total.n, impliedS=impliedS,
      matrices=matrices, cor.analysis = tssem1.obj$cor.analysis, intervals = intervals, ...)
}
