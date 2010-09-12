summary.wls <- function(object, ...) {
    if (!is.element("wls", class(object)))
    stop("\"object\" must be an object of class \"wls\".")

    n <- object$n
    tT <- object$wls.fit@output$Minus2LogLikelihood 
    dfT <- object$noObservedStat - summary(object$wls.fit)$estimatedParameters
    tB <- object$indepModelChisq
    dfB <- object$indepModelDf
    p <- pchisq(tT, df=dfT, lower.tail=FALSE)
    sampleS <- mxEval(sampleS, object$wls.fit)
    impliedS <- mxEval(impliedS, object$wls.fit)
    # it relies on that "Correlation structure" is used as the model name
    if (!is.na(match("Correlation structure", object$wls.fit@name))) {
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
    my.para <- summary(object$wls.fit)$parameters
    # For example, P[1,2], L[1,2], ...
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name

    my.ci <- summary(object$wls.fit)$CI
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
    
    out <- list(call=object$call, coefficients=coefficients, stat=stat)
    class(out) <- "summary.wls"
    out
}

print.summary.wls <- function(x, ...) {
    if (!is.element("summary.wls", class(x)))
    stop("\"x\" must be an object of class \"summary.wls\".")
    call.text <- deparse(x$call)
    cat("Call:\n")
    if (is.na(call.text[2]))
       cat(call.text[1])
    else cat(call.text[c(1,2)])
    cat("\n\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)
    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
}

summary.tssem1 <- function(object, ...) {
    if (!is.element("tssem1", class(object)))
    stop("\"object\" must be an object of class \"tssem1\".")
    
    # it relies on the model name
    if (!is.na(match("TSSEM1 Analysis of Correlation Matrix", object$tssem1.fit@name))) {
       cor.analysis <- TRUE
    } else {
       cor.analysis <- FALSE
    }   
 
    # Calculate the no. of variables based on the implied S
    mx.fit <- summary(object$tssem1.fit)
    # FIXME: what if there are incomplete data in S1
    no.var <- ncol(mxEval(S1, object$tssem1.fit))
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
    my.para <- my.para[my.para$matrix=="S1", ]
    #Sel <- grep("^S", my.para$matrix, value=TRUE)
    #my.para <- subset(my.para, my.para$matrix==Sel)
    my.para <- my.para[order(my.para$row, my.para$col), ]
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    dimnames(my.para)[[1]] <- my.para$name
    coefficients <- my.para[, c(5,6)]
    coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
    coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))
    
    out <- list(call=object$call, coefficients=coefficients, stat=stat)
    class(out) <- "summary.tssem1"
    out
}

   
print.summary.tssem1 <- function(x, ...) {
    if (!is.element("summary.tssem1", class(x)))
    stop("\"x\" must be an object of class \"summary.tssem1\".")
    call.text <- deparse(x$call)
    cat("Call:\n")
    if (is.na(call.text[2]))
       cat(call.text[1])
    else cat(call.text[c(1,2)])
    cat("\n\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)
    cat("\nGoodness-of-fit indices:\n")
    printCoefmat(x$stat, ...)
}

print.tssem1 <- function(x, ...) {
    if (!is.element("tssem1", class(x)))
      stop("\"x\" must be an object of class \"tssem1\".")
    call.text <- deparse(x$call)
    cat("Call:\n")
    if (is.na(call.text[2]))
       cat(call.text[1])
    else cat(call.text[c(1,2)])
    cat("\n\nStructure:\n")
    print(summary.default(x), ...)
}

print.wls <- function(x, ...) {
    if (!is.element("wls", class(x)))
      stop("\"x\" must be an object of class \"wls\".")
    call.text <- deparse(x$call)
    cat("Call:\n")
    if (is.na(call.text[2]))
       cat(call.text[1])
    else cat(call.text[c(1,2)])
    cat("\n\nStructure:\n")
    print(summary.default(x), ...)
}

print.meta <- function(x, ...) {
    if (!is.element("meta", class(x)))
      stop("\"x\" must be an object of class \"meta\".")
    cat("Call:\n")
    cat(deparse(x$call))
    cat("\n\nStructure:\n")
    print(summary.default(x), ...)
}

summary.meta <- function(object, ...) {
    if (!is.element("meta", class(object)))
    stop("\"object\" must be an object of class \"meta\".")

    # calculate coefficients    
    my.para <- summary(object$meta.fit)$parameters
    # For example, P[1,2], L[1,2], ...
    my.para$label <- my.para$name
    my.para$name <- with(my.para, paste(matrix,"[",row,",",col,"]",sep=""))
    
    my.ci <- summary(object$meta.fit)$CI
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
      name <- sapply(unlist(dimnames(my.ci)[1]), function(x)
                       {strsplit(x, "Meta analysis.", fixed=TRUE)[[1]][2]}, USE.NAMES=FALSE)
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
    
    out <- list(call=object$call, coefficients=coefficients)
    class(out) <- "summary.meta"
    out
}

print.summary.meta <- function(x, ...) {
    if (!is.element("summary.meta", class(x)))
    stop("\"x\" must be an object of class \"summary.meta\".")
    cat("Call:\n")
    cat(deparse(x$call))
    cat("\n\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)
    ## cat("\nGoodness-of-fit indices:\n")
    ## printCoefmat(x$stat, ...)
}
