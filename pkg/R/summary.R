summary.wls <- function(x, ...) {
    if (!is.element("wls", class(x)))
    stop("\"x\" must be an object of class \"wls\".")

    n <- x$n
    tT <- x$wls.fit@output$Minus2LogLikelihood 
    dfT <- x$noObservedStat - summary(x$wls.fit)$estimatedParameters
    tB <- x$indepModelChisq
    dfB <- x$indepModelDf
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)
    sampleS <- mxEval(sampleS, x$wls.fit)
    impliedS <- mxEval(impliedS, x$wls.fit)
    # it relies on that "Correlation structure" is used as the model name
    if (!is.na(match("Correlation structure", x$wls.fit@name))) {
       cor.analysis <- TRUE
    } else {
       cor.analysis <- FALSE
    }   
 
    RMSEA <- sqrt(max((tT-dfT)/(n-1),0)/dfT)
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
    AIC <- tT-2*dfT
    BIC <- tT-log(n)*dfT
    if (cor.analysis) {
      # diag(impliedS)!=1
      SRMR <- sqrt(mean(vechs(sampleS-impliedS)^2))
    } else {
      # standardize it according to Hu & Benterl (1998)
      stand <- diag(1/sqrt(diag(sampleS)))
      SRMR <- sqrt(mean(vech(stand %*% (sampleS-impliedS) %*% stand)^2))
    }
    
    stat <- matrix(c(n, tT, dfT, p, tB, dfB, RMSEA, SRMR, TLI, CFI, AIC, BIC), ncol=1)
    rownames(stat) <- c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Chi-square of independent model",
                        "DF of independent model", "RMSEA", "SRMR", "TLI", "CFI", "AIC", "BIC")
    colnames(stat) <- "Value"
  
    # calculate coefficients    
    my.para <- summary(x$wls.fit)$parameters
    # For example, P[1,2], L[1,2], ...
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name

    my.ci <- summary(x$wls.fit)$CI
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      my.para <- my.para[order(my.para$matrix, my.para$row, my.para$col), ]
      coefficients <- my.para[, -c(1:4)]
    } else {
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, "structure.", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- data.frame(name, my.ci)
      coefficients <- merge(my.para, my.ci, by=c("name"))
      dimnames(coefficients)[[1]] <- coefficients$name
      coefficients <- coefficients[order(coefficients$matrix, coefficients$row, coefficients$col), ]
      coefficients <- coefficients[, -c(1:4, 8)]
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    out <- list(call=x$call, coefficients=coefficients, stat=stat)
    class(out) <- "summary.wls"
    out
}

print.summary.wls <- function(x, ...) {
    if (!is.element("summary.wls", class(x)))
    stop("\"x\" must be an object of class \"summary.wls\".")
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)
    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
}

summary.tssem1 <- function(x, ...) {
    if (!is.element("tssem1", class(x)))
    stop("\"x\" must be an object of class \"tssem1\".")
    
    # it relies on the model name
    if (!is.na(match("TSSEM1 Analysis of Correlation Matrix", x$tssem1.fit@name))) {
       cor.analysis <- TRUE
    } else {
       cor.analysis <- FALSE
    }   
 
    # Calculate the no. of variables based on the implied S
    mx.fit <- summary(x$tssem1.fit)
    # FIXME: what if there are incomplete data in S1
    no.var <- ncol(mxEval(S1, x$tssem1.fit))
    if (cor.analysis)
      ps <- no.var*(no.var-1)/2
    else ps <- no.var*(no.var+1)/2
    n <- x$total.n
    tT <- x$modelMinus2LL - x$saturatedMinus2LL 
    dfT <- mx.fit$degreesOfFreedom
    tB <- x$independentMinus2LL - x$saturatedMinus2LL
    dfB <- mx.fit$degreesOfFreedom + ps
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)

    no.groups <- length(x$tssem1.fit@submodels)
    RMSEA <- sqrt(no.groups)*sqrt(max((tT-dfT)/(n-1),0)/dfT)
    TLI <- (tB/dfB - tT/dfT)/(tB/dfB-1)
    CFI <- 1 - max((tT-dfT),0)/max((tT-dfT),(tB-dfB),0)
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
    SRMR <- sqrt(mean(unlist(sapply(x$data, srmr, pooledS=x$pooledS, cor.analysis=cor.analysis))))
    #SRMR <- sqrt(mean(mapply(srmr, x$data, miss.index,
    #                         MoreArgs=list(pooledS=x$pooledS, cor.analysis=cor.analysis))))
    
    stat <- matrix(c(n, tT, dfT, p, tB, dfB, RMSEA, SRMR, TLI, CFI, AIC, BIC), ncol=1)
    rownames(stat) <- c("Sample size", "Chi-square of target model", "DF of target model",
                        "p value of target model", "Chi-square of independent model",
                        "DF of independent model", "RMSEA", "SRMR", "TLI", "CFI", "AIC", "BIC")
    colnames(stat) <- "Value"

    # calculate coefficients    
    my.para <- summary(x$tssem1.fit)$parameters
    my.para <- my.para[my.para$matrix=="S1", ]
    #Sel <- grep("^S", my.para$matrix, value=TRUE)
    #my.para <- subset(my.para, my.para$matrix==Sel)
    my.para <- my.para[order(my.para$row, my.para$col), ]
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name
    coefficients <- my.para[, c(5,6)]
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    out <- list(call=x$call, coefficients=coefficients, stat=stat)
    class(out) <- "summary.tssem1"
    out
}

   
print.summary.tssem1 <- function(x, ...) {
    if (!is.element("summary.tssem1", class(x)))
    stop("\"x\" must be an object of class \"summary.tssem1\".")
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)
    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
}

print.tssem1 <- function(x, ...) {
    if (!is.element("tssem1", class(x)))
      stop("\"x\" must be an object of class \"tssem1\".")
    cat("Call:\n")
    print(x$call)
    cat("\nStructure:\n")
    print(summary.default(x), ...)
}

print.wls <- function(x, ...) {
    if (!is.element("wls", class(x)))
      stop("\"x\" must be an object of class \"wls\".")
    cat("Call:\n")
    print(x$call)
    cat("\nStructure:\n")
    print(summary.default(x), ...)
}
