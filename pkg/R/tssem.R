tssem1FE <- function(my.df, n, cor.analysis=TRUE, model.name=NULL,
                     cluster=NULL, suppressWarnings=TRUE, ...) {
  if (!is.null(cluster)) {
    data.cluster <- tapply(my.df, cluster, function(x) {x})
    n.cluster <- tapply(n, cluster, function(x) {x})
    out <- list()
    for (i in 1:length(data.cluster)) {
      ## Need to correct it to tssem1()
      out[[i]] <- tssem1FE(data.cluster[[i]], n.cluster[[i]], 
                           cor.analysis=cor.analysis, model.name=model.name,
                           suppressWarnings=suppressWarnings, ...)
    }
    names(out) <- names(data.cluster) 
    class(out) <- "tssem1FE.cluster"
    out
  } else {
    my.range <- range(sapply(my.df, function(x) {ncol(x)}))
    if ( !all.equal(my.range[1], my.range[2]) )
      stop("Dimensions of all groups are not the same!\n")
    
    no.groups <- length(my.df)
    no.var <- ncol(my.df[[1]])
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
      diag(lbound) <- 0.00001
      S <- mxMatrix(type="Symm", nrow=no.var, ncol=no.var, free=TRUE, values=vech(sv),
                    labels=vech(ps.labels), lbound=lbound, name="S")      
    }
       
    ## Index for missing variables: only check the diagonals only!!!
    miss.index <- lapply(my.df, function(x) { is.na(diag(x)) })
    
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
            SD <- diag(paste(sqrt(diag(sv)), "*sd", i, "_", 1:no.var, sep=""))
            SD <- SD[!miss.index.i, ]
            SD <- as.mxMatrix(SD)   
        } else {
            if (is.null(model.name)) model.name <- "TSSEM1 Analysis of Covariance Matrix"
            SD <- diag(rep(1, no.var))
            SD <- SD[!miss.index.i, ]
            SD <- as.mxMatrix(SD)
        }
        # Expected covariance matrix
        expC <- mxAlgebra(SD %&% S, name="expC", dimnames=list(var.names[!miss.index.i],
                          var.names[!miss.index.i]))
        # Create mxModel
        ## model <- mxModel(paste("group",i,sep=""), S, Dmatrix, expC,
        ##                      mxData(observed=my.df.i, type="cov", numObs=n[i]),
        ##                      mxMLObjective(expC, dimnames=var.names[!miss.index[[i]]]))
        g.model <- paste("g", i, " <- mxModel(\"g", i, "\", S, SD, expC, mxData(observed=my.df[[",i,"]][!miss.index[[",i,"]],!miss.index[[",i,"]]], type=\"cov\", numObs=n[", i,
                "]), mxMLObjective(\"expC\", dimnames=var.names[!miss.index[[", i, "]]]))", sep = "")
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
                          ", name=\"obj\"), mxAlgebraObjective(\"obj\"))", sep = "")
    eval(parse(text = tssem1.model))
    
    # try to run it with error message as output
    ## tssem1.fit <- mxRun(tssem1)
    tssem1.fit <- tryCatch(mxRun(tssem1, suppressWarnings = suppressWarnings, ...),
                           error = function(e) e)
    if (inherits(tssem1.fit, "error")) 
        stop(print(tssem1.fit))
    
    pooledS <- eval(parse(text = "mxEval(S, tssem1.fit)"))
    
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
    
    out <- list(call = match.call(), cor.analysis=cor.analysis, data=my.df, pooledS = pooledS,
                acovS = acovS, total.n = total.n, 
                modelMinus2LL = tssem1.fit@output$Minus2LogLikelihood,
                independentMinus2LL = independentMinus2LL, saturatedMinus2LL = saturatedMinus2LL,
                tssem1.fit = tssem1.fit)
    class(out) <- "tssem1FE"
    return(out)
  }
}

tssem1RE <- function(my.df, n, cor.analysis=TRUE, RE_diag=FALSE, RE.startvalues=0.1, RE.lbound = 1e-10,
                     model.name=NULL, suppressWarnings=TRUE, ...) {
  ## Replace diagonals with 1.0
  my.complete <- lapply(my.df, function (x) { diag(x)[is.na(diag(x))] <- 1; x })
  ## Replace missing variables with 0.0
  my.complete <- lapply(my.complete, function (x) { x[is.na(x)] <- 0; x })
  
  ## Calculate the asymptotic sampling covariance matrix of the correlation matrix
  acovR <- asyCov(x=my.complete, n=n, cor.analysis=cor.analysis)
  
  ## Convert the correlation matrices into effect sizes; default of diag=FALSE for cor matrix
  ES <- list2matrix(x=my.df, diag=!cor.analysis)
  ## no. of effect sizes
  no.es <- ncol(ES)

  if (is.null(model.name)) {
    if (cor.analysis) {
      model.name <- "TSSEM1 (Random Effects Model) Analysis of Correlation Matrix"
    } else {
      model.name <- "TSSEM1 (Random Effects Model) Analysis of Covariance Matrix"
    }
  }
  
  if (RE_diag==TRUE) {
    ## No covariance between random effects
    meta.fit <- meta(y=ES, v=acovR, model.name=model.name,
                     RE.constraints=diag(x=paste(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es, sep=""),
                                         nrow=no.es, ncol=no.es))    
  } else {
    meta.fit <- meta(y=ES, v=acovR, model.name=model.name, RE.startvalues=RE.startvalues, RE.lbound = RE.lbound)
  }

  out <- list(total.n=sum(n), cor.analysis=cor.analysis, RE_diag=RE_diag, no.es=no.es)
  out <- c(out, meta.fit)
  class(out) <- c("tssem1RE", "meta")
  return(out)
}

tssem1 <- function(my.df, n, method=c("FE", "RE"), cor.analysis=TRUE, cluster=NULL,
                   RE_diag=FALSE, RE.startvalues=0.1, RE.lbound = 1e-10,
                   model.name=NULL, suppressWarnings=TRUE, ...) {
  method <- match.arg(method)
  switch(method,
    FE = out <- tssem1FE(my.df=my.df, n=n, cor.analysis=cor.analysis, model.name=model.name,
                          cluster=cluster, suppressWarnings=suppressWarnings, ...),
    RE = out <- tssem1RE(my.df=my.df, n=n, cor.analysis=cor.analysis, RE_diag=RE_diag,
                          RE.startvalues=RE.startvalues, RE.lbound=RE.lbound,
                          model.name=model.name, suppressWarnings=suppressWarnings, ...) )
  out  
}
  
wls <- function(S, acovS, n, impliedS, matrices, cor.analysis = TRUE,
                intervals.type =c("z", "LB"), model.name, suppressWarnings = TRUE, ...) {
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
        if (missing(model.name)) model.name <- "WLS Analysis of Correlation Structure"
        ps <- no.var * (no.var - 1)/2
        vecS <- mxAlgebra(vechs(sampleS - impliedS), name = "vecS")
    } else {
        if (missing(model.name)) model.name <- "WLS Analysis of Covariance Structure"
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
        out <- list(call = match.call(), noObservedStat=ps, n=n, cor.analysis=cor.analysis,
                    indepModelChisq=.indepwlsChisq(S=S, acovS=acovS, cor.analysis=cor.analysis),
                    indepModelDf=no.var*(no.var-1)/2, wls.fit=wls.fit)
        class(out) <- 'wls'
    }
    out
}


tssem2 <- function(tssem1.obj, impliedS, matrices, intervals.type = c("z", "LB"),
                   model.name=NULL, suppressWarnings = TRUE, ...) {
  if ( !is.element( class(tssem1.obj)[1], c("tssem1FE.cluster", "tssem1FE", "tssem1RE")) )
      stop("\"tssem1.obj\" must be of neither class \"tssem1FE.cluster\", class \"tssem1FE\" or \"tssem1RE\".")
      
  switch(class(tssem1.obj)[1],
         tssem1FE.cluster = { out <- lapply(tssem1.obj, tssem2, impliedS=impliedS, matrices=matrices,
                                            intervals.type=intervals.type, model.name=model.name,
                                            suppressWarnings=suppressWarnings, ...)
                              class(out) <- "wls.cluster" },
         tssem1FE = { cor.analysis <- tssem1.obj$cor.analysis
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
                      out <- wls(S=tssem1.obj$pooledS, acovS=tssem1.obj$acovS, n=tssem1.obj$total.n,
                                 impliedS=impliedS, matrices=matrices, cor.analysis=cor.analysis,
                                 intervals.type=intervals.type, model.name=model.name,
                                 suppressWarnings = suppressWarnings, ...) },
         tssem1RE = { cor.analysis <- tssem1.obj$cor.analysis
                     ## Extract the pooled correlation matrix
                      pooledS <- vec2symMat( coef(tssem1.obj, select="fixed"), diag=!cor.analysis)
                     ## Extract the asymptotic covariance matrix of the pooled correlations
                      acovS <- vcov(tssem1.obj, select="fixed")
                      
                      if (cor.analysis==TRUE) {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Random Effects Model) Analysis of Correlation Structure"
                      } else {
                        if (is.null(model.name)) model.name <- "TSSEM2 (Random Effects Model) Analysis of Covariance Structure"
                      }
                      out <- wls(S=pooledS, acovS=acovS, n=tssem1.obj$total.n, impliedS=impliedS,
                                 matrices=matrices, cor.analysis=cor.analysis, intervals.type=intervals.type,
                                 model.name=model.name, suppressWarnings = suppressWarnings, ...) })
  out
}
    
