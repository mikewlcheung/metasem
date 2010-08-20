summary.wls <- function(obj, ...) {
    if (!is.element("wls", class(obj)))
    stop("\"obj\" must be an object of class \"wls\".")

    n <- obj$n
    tT <- obj$wls.fit@output$Minus2LogLikelihood 
    dfT <- obj$noObservedStat - summary(obj$wls.fit)$estimatedParameters
    tB <- obj$indepModelChisq
    dfB <- obj$indepModelDf
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)
    sampleS <- mxEval(sampleS, obj$wls.fit)
    impliedS <- mxEval(impliedS, obj$wls.fit)
    # it relies on that "Correlation structure" is used as the model name
    if (!is.na(match("Correlation structure", obj$wls.fit@name))) {
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
    my.para <- summary(obj$wls.fit)$parameters
    # For example, P[1,2], L[1,2], ...
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name

    my.ci <- summary(obj$wls.fit)$CI
    # Determine if CIs on parameter estimates are present
    if (is.null(dimnames(my.ci))) {
      my.para$lbound <- my.para$Estimate - qnorm(.975)*my.para$Std.Error
      my.para$ubound <- my.para$Estimate + qnorm(.975)*my.para$Std.Error
      coefficients <- my.para[, -c(1:4)]
    } else {
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, "structure.", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
      my.ci <- data.frame(name, my.ci)
      coefficients <- merge(my.para, my.ci, by=c("name"))
      dimnames(coefficients)[[1]] <- coefficients$name
      coefficients <- coefficients[, -c(1:4, 8)]
    }
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    out <- list(call=obj$call, coefficients=coefficients, stat=stat)
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

