summary.wls <- function(object, ...) {
    if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")

    n <- object$n
    tT <- object$wls.fit@output$Minus2LogLikelihood
    ## Adjust for the constraints on the diagonals
    dfT <- object$noObservedStat - summary(object$wls.fit)$estimatedParameters + object$noConstraints
    tB <- object$indepModelChisq
    dfB <- object$indepModelDf
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)
    sampleS <- mxEval(sampleS, object$wls.fit)
    impliedS <- mxEval(impliedS, object$wls.fit)
    # it relies on that "Correlation structure" is used as the model name
    ## if (!is.na(match("Correlation structure", object$wls.fit@name))) {
    ##    cor.analysis <- TRUE
    ## } else {
    ##    cor.analysis <- FALSE
    ## }
    cor.analysis <- object$cor.analysis

    ## Hu and Bentler (1998) Psychological Methods
    ## Protect RMSEA divided by 0
    if (min((tT-dfT),0)==0) {
      RMSEA <- 0
    } else {
      RMSEA <- sqrt(max((tT-dfT)/(n-1),0)/dfT)
    }
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
	## FIXME: better to use -2LL+2r where r is the no. of free parameters (Mplus, p. 22; R ?AIC)
	## This will be more consistent with logLik(). However, it seems that df not npar is used in logLik().
    AIC <- tT-2*dfT
    BIC <- tT-log(n)*dfT
    if (cor.analysis) {
      # diag(impliedS)!=1
      SRMR <- sqrt(mean(vechs(sampleS-impliedS)^2))
    } else {
      # standardize it according to Hu & Bentler (1998)
      stand <- diag(1/sqrt(diag(sampleS)))
      SRMR <- sqrt(mean(vech(stand %*% (sampleS-impliedS) %*% stand)^2))
    }
    
    stat <- matrix(c(n, tT, dfT, p, tB, dfB, object$noConstraints, RMSEA, SRMR, TLI, CFI, AIC, BIC), ncol=1)
    rownames(stat) <- c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Chi-square of independent model",
                        "DF of independent model", "No. of constraints imposed on \"S\"", "RMSEA", "SRMR", "TLI", "CFI", "AIC", "BIC")
    colnames(stat) <- "Value"
  
    # calculate coefficients
    my.mx <- summary(object$wls.fit)
    ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
    my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1 
    # For example, P[1,2], L[1,2], ...
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name

    my.ci <- my.mx$CI
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      coefficients <- my.para[, -c(1:4)]
	  intervals.type="z"
    } else {
      ## # model.name: may vary in diff models
      ## model.name <- object$call[[match("model.name", names(object$call))]]
      ## # if not specified, the default is "Structure" as it can be either "Correlation Structure" or "Covariance Structure" 
      ## if (is.null(model.name)) {
      ##   model.name <- "Structure."
      ## } else {
      ##   model.name <- paste(model.name, ".", sep="")
      ## }          
      ## name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
      ##                  {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)

      ## Simply remove the part before "."
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, ".", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)      

      my.ci <- data.frame(name, my.ci)
      coefficients <- merge(my.para, my.ci, by=c("name"))
      dimnames(coefficients)[[1]] <- coefficients$name
      coefficients <- coefficients[order(coefficients$matrix, coefficients$row, coefficients$col), ]
      coefficients <- coefficients[, -c(1:4, 8)]
	  intervals.type="LB"
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    ## Mx.status1 <- object$wls.fit@output$status[[1]]    
    ## libMatrix <- installed.packages()
    ## out <- list(call=object$call, coefficients=coefficients, stat=stat, intervals.type=intervals.type,
    ##             Mx.status1=Mx.status1, R.version=as.character(getRversion()),
    ##             OpenMx.version=libMatrix["OpenMx", "Version"],
    ##             metaSEM.version=libMatrix["metaSEM", "Version"], date=date())
    out <- list(call=object$call, coefficients=coefficients, stat=stat, intervals.type=intervals.type,
                Mx.status1=object$wls.fit@output$status[[1]])
    class(out) <- "summary.wls"
    out
}

print.summary.wls <- function(x, ...) {
    if (!is.element("summary.wls", class(x)))
    stop("\"x\" must be an object of class \"summary.wls\".")
    ## call.text <- deparse(x$call)
    ## cat("Call:\n")
    ## for (i in 1:length(call.text)) {
    ##     cat(call.text[i], "\n")
    ## }	

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
    
    ## cat("\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")
}

summary.tssem1FEM <- function(object, ...) {
    if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")
    
    # it relies on the model name
    ## if (!is.na(match("TSSEM1 Analysis of Correlation Matrix", object$tssem1.fit@name))) {
    ##    cor.analysis <- TRUE
    ## } else {
    ##    cor.analysis <- FALSE
    ## }   

    ## if ("Correlation" %in% unlist(strsplit(object$tssem1.fit@name, " "))) {
    ##   cor.analysis <- TRUE
    ## } else {
    ##   cor.analysis <- FALSE
    ## }
    cor.analysis <- object$cor.analysis
    
    # Calculate the no. of variables based on the implied S
    mx.fit <- summary(object$tssem1.fit)

    ## Fixed a warning in R CMD check
    ## summary.tssem1: no visible binding for global variable ‘S1’
    no.var <- ncol(object$pooledS)

    # Fixed if there are incomplete data in the first group
    ## no.var <- ncol(mxEval(S1, object$tssem1.fit))
    if (cor.analysis)
      ps <- no.var*(no.var-1)/2
    else ps <- no.var*(no.var+1)/2
    n <- object$total.n
    tT <- object$modelMinus2LL - object$saturatedMinus2LL 
    dfT <- mx.fit$degreesOfFreedom
    tB <- object$independentMinus2LL - object$saturatedMinus2LL
    dfB <- mx.fit$degreesOfFreedom + ps
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)

    no.groups <- length(object$tssem1.fit@submodels)
    RMSEA <- sqrt(no.groups)*sqrt(max((tT-dfT)/(n-1),0)/dfT)
    ## Hu and Bentler (1998)
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
    ## FIXME: definitions of AIC
    AIC <- tT-2*dfT
    BIC <- tT-log(n)*dfT
    
    ## Index for missing variables: only check the diagonals only!!!
    #miss.index <- lapply(x$data, function(x) { is.na(diag(x)) })
    srmr <- function(sampleS, pooledS, cor.analysis) {
      index <- is.na(diag(sampleS))
      # Selection matrix to select complete data
      Sel <- diag(rep(1, ncol(pooledS)))
      Sel <- Sel[!index, ]
      sampleS <- sampleS[!index, !index]
      if (cor.analysis) {
        vechs( cov2cor(sampleS)- Sel %*% pooledS %*% t(Sel) )^2
      } else {
        stand <- diag(1/sqrt(diag(sampleS)))
        vech( stand %*% (sampleS - Sel %*% pooledS %*% t(Sel) ) %*% stand )^2
      }
    }
    SRMR <- sqrt(mean(unlist(sapply(object$data, srmr, pooledS=object$pooledS, cor.analysis=cor.analysis))))
    #SRMR <- sqrt(mean(mapply(srmr, x$data, miss.index,
    #                         MoreArgs=list(pooledS=x$pooledS, cor.analysis=cor.analysis))))
    
    stat <- matrix(c(n, tT, dfT, p, tB, dfB, RMSEA, SRMR, TLI, CFI, AIC, BIC), ncol=1)
    rownames(stat) <- c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Chi-square of independent model",
                        "DF of independent model", "RMSEA", "SRMR", "TLI", "CFI", "AIC", "BIC")
    colnames(stat) <- "Value"

    # calculate coefficients    
    my.para <- summary(object$tssem1.fit)$parameters
    ## my.para <- my.para[my.para$matrix=="S1", ]
    my.para <- my.para[my.para$matrix=="S", ]
    #Sel <- grep("^S", my.para$matrix, value=TRUE)
    #my.para <- subset(my.para, my.para$matrix==Sel)
    my.para <- my.para[order(my.para$row, my.para$col), ]
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name
    coefficients <- my.para[, c(5,6)]
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    ## Mx.status1 <- object$tssem1.fit@output$status[[1]]   
    ## libMatrix <- installed.packages()    
    ## out <- list(call=object$call, coefficients=coefficients, stat=stat, Mx.status1=Mx.status1,
    ##             R.version=as.character(getRversion()), OpenMx.version=libMatrix["OpenMx", "Version"],
    ##             metaSEM.version=libMatrix["metaSEM", "Version"], date=date())
    out <- list(call=object$call, coefficients=coefficients, stat=stat, Mx.status1=object$tssem1.fit@output$status[[1]])
    class(out) <- "summary.tssem1FEM"
    out
}
   
print.summary.tssem1FEM <- function(x, ...) {
    if (!is.element("summary.tssem1FEM", class(x)))
    stop("\"x\" must be an object of class \"summary.tssem1FEM\".")
    ## call.text <- deparse(x$call)
    
    ## cat("Call:\n")
    ## # Ad-hoc solution to remove very long call text
	## if ( length(call.text)>30 ) {
    ##     cat(call.text[1], "\n")
    ## } else { 		
    ##     for (i in 1:length(call.text)) {
    ##         cat(call.text[i], "\n")
    ##     }	
    ## }
    
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Coefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)

    ## cat("\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("OpenMx status:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")    
}

print.tssem1FEM <- function(x, ...) {
    if (!is.element("tssem1FEM", class(x)))
      stop("\"x\" must be an object of class \"tssem1FEM\".")
    ## call.text <- deparse(x$call)

    ## cat("Call:\n")
    ## # Ad-hoc solution to remove very long call text
	## if ( length(call.text)>30 ) {
    ##     cat(call.text[1], "\n")
    ## } else { 	
    ##     for (i in 1:length(call.text)) {
    ##         cat(call.text[i], "\n")
    ##     }	
    ## }
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

print.tssem1FEM.cluster <- function(x, ...) {
    if (!is.element("tssem1FEM.cluster", class(x)))
      stop("\"x\" must be an object of class \"tssem1FEM.cluster\".")

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

print.tssem1REM <- function(x, ...) {
    if (!is.element("tssem1REM", class(x)))
      stop("\"x\" must be an object of class \"tssem1REM\".")
    ## call.text <- deparse(x$call)

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")
    print(summary.default(x), ...)
}

summary.tssem1REM <- function(object, ...) {
    if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")
    summary.meta(object)
}
   
print.wls <- function(x, ...) {
    if (!is.element("wls", class(x)))
      stop("\"x\" must be an object of class \"wls\".")
    ## call.text <- deparse(x$call)
    ## cat("Call:\n")
    ## for (i in 1:length(call.text)) {
    ##     cat(call.text[i], "\n")
    ## }
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")    
    cat("Structure:\n")	
    print(summary.default(x), ...)
}

print.meta <- function(x, ...) {
    if (!is.element("meta", class(x)))
      stop("\"x\" must be an object of class \"meta\".")
    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")     
    cat("Structure:\n")
    print(summary.default(x), ...)
}

summary.meta <- function(object, homoStat=TRUE, ...) {
    if (!is.element("meta", class(object)))
    ## if (!("meta" %in% class(object)))
    stop("\"object\" must be an object of class \"meta\".")

    # calculate coefficients    
    my.mx <- summary(object$meta.fit)
    ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
    my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1
    ## OpenMx1.1: y1, y2 and x1 appear in col
    my.para$col <- sub("[a-z]", "", my.para$col)  # Fixed for OpenMx1.1
    # For example, P[1,2], L[1,2], ...
    my.para$label <- my.para$name
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    
    my.ci <- my.mx$CI
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      # remove rows with missing labels
      my.para <- my.para[!is.na(my.para$label), ]
      coefficients <- my.para[, -c(1:4,7)]
      dimnames(coefficients)[[1]] <- my.para$label 
    } else {
      # model.name: may vary in diff models
      model.name <- object$call[[match("model.name", names(object$call))]]
      # if not specified, the default in meta() is "Meta analysis with ML"
      if (is.null(model.name)) {
        model.name <- "Meta analysis with ML."
      } else {
        model.name <- paste(model.name, ".", sep="")
      }
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- data.frame(name, my.ci)
      my.para <- merge(my.para, my.ci, by=c("name"))
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      # remove rows with missing labels
      my.para <- my.para[!is.na(my.para$label), ]
      coefficients <- my.para[, -c(1:4,7,9)]
      dimnames(coefficients)[[1]] <- my.para$label 
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    intervals.type <- object$call[[match("intervals.type", names(object$call))]]
    # default
    if (is.null(intervals.type))
      intervals.type <- "z"

    ## Homogeneity statistic
    no.y <- object$no.y
    no.v <- no.y*(no.y+1)/2
	# Remove studies that have missing x. Make sure that studies are the same in calculating Q.stat and meta()
    if (homoStat) {
      Q.stat <- homoStat(y=object$data[!object$miss.x, 1:no.y], v=object$data[!object$miss.x, (no.y+1):(no.y+no.v)])
    } else {
      Q.stat <- list(Q=NA, Q.df=NA, pval=NA)
    }

    ## Calculate I2 only if no.x=0
    if (object$no.x==0) {
      ## Test whether it is meta2 or meta3
      if (object$type=="meta2") {
        ## Need to do sth for meta2
        heter.values <- NULL
      } else {
        heter.indices <- object$call[[match("heter.indices", names(object$call))]]
        if (is.null(heter.indices)) {
          heter.indices <- "I2hm"
        } else {
          ## remove the first "c" character
          heter.indices <- as.character(heter.indices)[-1]
        }
        heter.names <-c(outer(heter.indices, c("_2","_3"), paste, sep=""))
        
        ## Wald test, no CI
        if (is.null(dimnames(my.ci))) {
          my.heter <- paste("mxEval(c(", paste(heter.names, collapse=","), "), object$meta.fit)", sep="")
          heter.values <- matrix(NA, nrow=length(heter.names), ncol=3)
          heter.values[,2] <- eval(parse(text = my.heter))
        } else {
        ## LB CI  
        # model.name <- "Meta analysis with ML."
          heter.values <- my.mx$CI[paste(model.name, heter.names, "[1,1]", sep=""), ]
        }

        heter.names <- sub("I2q_2", "I2_2 (Q statistic)", heter.names)
        heter.names <- sub("I2q_3", "I2_3 (Q statistic)", heter.names)
        heter.names <- sub("I2hm_2", "I2_2 (harmonic mean)", heter.names)
        heter.names <- sub("I2hm_3", "I2_3 (harmonic mean)", heter.names)
        heter.names <- sub("I2am_2", "I2_2 (arithmetic mean)", heter.names)
        heter.names <- sub("I2am_3", "I2_3 (arithmetic mean)", heter.names)
        dimnames(heter.values) <- list(heter.names, c("lbound", "Estimate", "ubound"))
      }
    } else {
      ## no.x != 0
      heter.values <- NULL
    }
    
    ## Mx.status1 <- object$meta.fit@output$status[[1]]   
    ## libMatrix <- installed.packages()    
    ## out <- list(call=object$call, type=object$type, Q.stat=Q.stat, intervals.type=intervals.type,
    ##             heter.values=heter.values, no.studies=my.mx$numObs,
    ##             obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
    ##             df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
    ##             coefficients=coefficients, Mx.status1=Mx.status1, R.version=as.character(getRversion()),
    ##             OpenMx.version=libMatrix["OpenMx", "Version"], metaSEM.version=libMatrix["metaSEM", "Version"],
    ##             date=date())
    out <- list(call=object$call, type=object$type, Q.stat=Q.stat, intervals.type=intervals.type,
                heter.values=heter.values, no.studies=my.mx$numObs,
                obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
                df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
                coefficients=coefficients, Mx.status1=object$meta.fit@output$status[[1]])
    class(out) <- "summary.meta"
    out
}

print.summary.meta <- function(x, ...) {
    if (!is.element("summary.meta", class(x)))
    stop("\"x\" must be an object of class \"summary.meta\".")

    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nQ statistic on homogeneity of effect sizes:", x$Q.stat[["Q"]])
    cat("\nDegrees of freedom of the Q statistic:", x$Q.stat[["Q.df"]])
    cat("\nP value of the Q statistic:", x$Q.stat[["pval"]])
    cat("\n\nNumber of studies (or clusters):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL)        

    if (x$type=="meta3") {
    ## Print heterogeneity indices if no x in the call
    if ( is.null(x$call[[match("x", names(x$call))]]) ) {
      switch(x$intervals.type,
             z =  { cat("\n\nHeterogeneity indices:\n")
                    printCoefmat(x$heter.values[,2, drop=FALSE], ...) },
             LB = { cat("\n\nHeterogeneity indices (and 95% likelihood-based CIs):\n")
                    printCoefmat(x$heter.values, ...) } )
    }
    }
    
    ## cat("\n\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("\nOpenMx status1:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")      
}

summary.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")

    # calculate coefficients    
    my.mx <- summary(object$reml.fit)
    ## my.para <- my.mx$parameters       # Worked up to OpenMx1.0.6
    my.para <- my.mx$parameters[, 1:6]   # Fixed for OpenMx1.1 
    # For example, P[1,2], L[1,2], ...
    my.para$label <- my.para$name
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    
    my.ci <- my.mx$CI
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      # remove rows with missing labels
      my.para <- my.para[!is.na(my.para$label), ]
      coefficients <- my.para[, -c(1:4,7)]
      dimnames(coefficients)[[1]] <- my.para$label 
    } else {
      # model.name: may vary in diff models
      model.name <- object$call[[match("model.name", names(object$call))]]
      # if not specified, the default in meta() is "Variance component with REML"
      if (is.null(model.name)) {
        model.name <- "Variance component with REML."
      } else {
        model.name <- paste(model.name, ".", sep="")
      }      
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, model.name, fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      # remove duplicate elements in my.ci from my.para$name
      name.sel <- name %in% my.para$name
      my.ci <- data.frame(name=my.para$name, my.ci[name.sel, ,drop=FALSE])
      my.para <- merge(my.para, my.ci, by=c("name"))
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      # remove rows with missing labels
      my.para <- my.para[!is.na(my.para$label), ]
      coefficients <- my.para[, -c(1:4,7,9)]
      dimnames(coefficients)[[1]] <- my.para$label 
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

    intervals.type <- object$call[[match("intervals.type", names(object$call))]]
    # default
    if (is.null(intervals.type))
      intervals.type <- "z"

    ## Mx.status1 <- object$reml.fit@output$status[[1]]   
    ## libMatrix <- installed.packages()    
    ## out <- list(call=object$call, intervals.type=intervals.type, no.studies=my.mx$numObs,
    ##             obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
    ##             df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
    ##             coefficients=coefficients, Mx.status1=Mx.status1, R.version=as.character(getRversion()),
    ##             OpenMx.version=libMatrix["OpenMx", "Version"], metaSEM.version=libMatrix["metaSEM", "Version"],
    ##             date=date())
    out <- list(call=object$call, intervals.type=intervals.type, no.studies=my.mx$numObs,
                obsStat=my.mx$observedStatistics, estPara=my.mx$estimatedParameters,
                df=my.mx$degreesOfFreedom, Minus2LL=my.mx$Minus2LogLikelihood,
                coefficients=coefficients, Mx.status1=object$reml.fit@output$status[[1]])
    class(out) <- "summary.reml"
    out
}

print.summary.reml <- function(x, ...) {
    if (!is.element("summary.reml", class(x)))
    stop("\"x\" must be an object of class \"summary.reml\".")
    
    ## cat("Call:\n")
    ## cat(deparse(x$call))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation"),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    cat("\nNumber of studies (or clusters):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL)        

    ## cat("\n\nR version:", x$R.version)
    ## cat("\nOpenMx version:", x$OpenMx.version)
    ## cat("\nmetaSEM version:", x$metaSEM.version)
    ## cat("\nDate of analysis:", x$date)
    cat("OpenMx status:", x$Mx.status1, "(\"0\" and \"1\": considered fine; other values indicate problems)\n")
    ## cat("\nSee http://openmx.psyc.virginia.edu/wiki/errors for the details.\n\n")      
}

print.reml <- function(x, ...) {
    if (!is.element("reml", class(x)))
      stop("\"x\" must be an object of class \"reml\".")
    ## cat("Call:\n")
    ## cat(deparse(x$call))

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")     
    cat("Structure:\n")
    print(summary.default(x), ...)
}

## FIXME: no name when there is only one element
vcov.meta <- function(object, select=c("all", "fixed", "random"), ...) {
    if (!is.element("meta", class(object)))
    stop("\"object\" must be an object of class \"meta\".")

    # labels of the parameters    
    ## my.name <- summary(object$meta.fit)$parameters$name
    my.name <- names( omxGetParameters(object$meta.fit) )
    my.name <- my.name[!is.na(my.name)]
    
    select <- match.arg(select)
    switch( select,
         ## all = my.name <- my.name,
         fixed =  my.name <- my.name[ grep("Intercept|Slope", my.name) ],
         random = my.name <- my.name[ grep("Tau2", my.name) ]
         )    
    ## switch( select,
    ##        ## all = my.name <- my.name,
    ##        fixed = if (no.x==0) {
    ##          my.name <- paste("Intercept", 1:no.y, sep="")
    ##        } else {
    ##          my.name <- c( paste("Intercept", 1:no.y, sep=""),
    ##                        outer(1:no.y, 1:no.x, function(y, x) paste("Slope", y,"_", x, sep = "")) )
    ##        },
    ##        random = if ("tssem1REM" %in% class(object)) {
    ##                     my.name <- paste("Tau2_", 1:no.y,"_", 1:no.y,sep="")
    ##                 } else {
    ##                     my.name <- vech(outer(1:no.y, 1:no.y, function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    ##                 }
    ##        )

    acov <- tryCatch( 2*solve(object$meta@output$calculatedHessian[my.name, my.name]), error = function(e) e)
    # Issue a warning instead of error message
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      warning(print(acov))
    } else {
      return(acov)
    }
}

vcov.tssem1FEM <- function(object, ...) {
  if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")
  object$acovS
}

vcov.tssem1FEM.cluster <- function(object, ...) {
  if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
  lapply(object, vcov.tssem1FEM)
}  

vcov.tssem1REM <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")
  vcov.meta(object, select, ...)
}

vcov.wls <- function(object, ...) {
  if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")
  acovS <- tryCatch( 2*solve(object$wls.fit@output$calculatedHessian), error = function(e) e ) 
  # Issue a warning instead of error message
  if (inherits(acovS, "error")) {
    cat("Error in solving the Hessian matrix.\n")
    warning(print(acovS))
  } else {
    my.mx <- summary(object$wls.fit)
    # For example, P[1,2], L[1,2], ...
    my.names <- with(my.mx$parameters[, 2:4], paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(acovS) <- list(my.names, my.names)
    acovS
  }
}

vcov.wls.cluster <- function(object, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, vcov.wls)
}
  
vcov.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")

    ## # labels of the parameters    
    ## my.name <- summary(object$reml.fit)$parameters$name
    my.name <- names(omxGetParameters(object$reml.fit))
    my.name <- my.name[!is.na(my.name)]
    acov <- tryCatch( 2*solve(object$reml@output$calculatedHessian[my.name, my.name]), error = function(e) e) 
    
    if (inherits(acov, "error")) {
      cat("Error in solving the Hessian matrix.\n")
      stop(print(acov))
    } else {
      return(acov)
    }
}
  

coef.meta <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("meta", class(object)))
    stop("\"object\" must be an object of class \"meta\".")

  ## # labels of the parameters    
  ## my.name <- summary(object$meta.fit)$parameters$name
  ## my.name <- my.name[!is.na(my.name)]

  my.para <- omxGetParameters(object$meta.fit)
  select <- match.arg(select)
  switch( select,
         fixed =  my.para <- my.para[ grep("Intercept|Slope", names(my.para)) ],
         random = my.para <- my.para[ grep("Tau2", names(my.para)) ]
         )
  ## switch( select,
  ##        ## all = my.name <- my.name,
  ##        fixed =  my.name <- my.name[ grep("Intercept|Slope", my.name) ],
  ##        random = my.name <- my.name[ grep("Tau2", my.name) ]
  ##        )
  ## object$meta.fit@output$estimate[my.name]
  my.para
}

coef.tssem1FEM <- function(object, ...) {  
    if (!is.element("tssem1FEM", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM\".")
    object$pooledS
}

coef.tssem1REM <- function(object, select=c("all", "fixed", "random"), ...) {
  if (!is.element("tssem1REM", class(object)))
    stop("\"object\" must be an object of class \"tssem1REM\".")
  coef.meta(object, select, ...)
}

coef.tssem1FEM.cluster <- function(object, ...) {
    if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
    lapply(object, coef.tssem1FEM)
}
  
coef.wls <- function(object, ...) {
    if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")
    ## object$wls.fit@output$estimate
    my.mx <- summary(object$wls.fit)
    my.coef <- my.mx$parameters$Estimate
    # For example, P[1,2], L[1,2], ...
    names(my.coef) <- with(my.mx$parameters[, 2:4], paste(matrix,"[",row,",",col,"]",sep=""))
    my.coef
}

coef.wls.cluster <- function(object, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, coef.wls)
}

coef.reml <- function(object, ...) {
    if (!is.element("reml", class(object)))
    stop("\"object\" must be an object of class \"reml\".")
    # labels of the parameters    
    ## my.name <- summary(object$reml.fit)$parameters$name
    my.name <- names(omxGetParameters(object$reml.fit))
    my.name <- my.name[!is.na(my.name)]
    object$reml.fit@output$estimate[my.name]
}

anova.meta <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$meta.fit)
  comparison <- lapply(list(...), function(x) x$meta.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

anova.wls <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$wls.fit)
  comparison <- lapply(list(...), function(x) x$wls.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

anova.reml <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$reml.fit)
  comparison <- lapply(list(...), function(x) x$reml.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.tssem1FEM.cluster <- function(object, ...) {
    if (!is.element("tssem1FEM.cluster", class(object)))
    stop("\"object\" must be an object of class \"tssem1FEM.cluster\".")
    lapply(object, summary.tssem1FEM)
}

summary.wls.cluster <- function(object, ...) {
    if (!is.element("wls.cluster", class(object)))
    stop("\"object\" must be an object of class \"wls.cluster\".")
    lapply(object, summary.wls)
}
