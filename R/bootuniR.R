## Generate population correlation matrices from uniR1 object
bootuniR1 <- function(x, Rep, nonPD.pop=c("replace", "nearPD", "accept")) {
  if (!is.element("uniR1", class(x)))
      stop("'x' must be an object of class 'uniR1'.")
  
  rCorPop(Sigma=x$r.mean, V=diag(vechs(x$r.SD^2)),
          corr=TRUE, k=Rep, nonPD.pop=nonPD.pop)
}

bootuniR2 <- function(model, data, n, ...) {
  if (!requireNamespace("lavaan", quietly=TRUE))
    stop("\"lavaan\" package is required for this function.")

  out <- lapply(data, function(x) {try( lavaan::sem(model=model, sample.cov=x,
                                                    sample.nobs=n, ...))})
  ## ## remove replications with errors
  ## out <- out[sapply(out, function(x) !inherits(x, "try-error"))]
  class(out) <- "bootuniR2"
  out
}

summary.bootuniR2 <- function(object, probs=c(0, 0.1, 0.5, 0.9, 1),
                              cutoff.chisq.pvalue=0.05, cutoff.CFI=0.9,
                              cutoff.SRMR=0.1, cutoff.RMSEA=0.05, ...) {
  if (!is.element("bootuniR2", class(object)))
    stop("'object' must be an object of class 'bootuniR2'.")

  Rep.total <- length(object)

  ## remove replications with errors
  object <- object[sapply(object, function(x) !inherits(x, "try-error"))]
  ## remove not converged replications
  object <- object[sapply(object, lavaan::lavInspect, what="converged")]

  para <- sapply(object, lavaan::coef)

  Quantile1 <- t(apply(para, 1, quantile, probs=probs))
  colnames(Quantile1) <- paste0("Quantile: ", colnames(Quantile1))
  coefficients <- cbind(Mean=apply(para, 1, mean, na.rm=TRUE),
                        SD=apply(para, 1, sd, na.rm=TRUE),
                        Quantile=Quantile1)

  fit <- t(sapply(object, lavaan::fitMeasures,
                  fit.measures=c("chisq", "df", "pvalue", "cfi", "srmr", "rmsea")))

  Quantile2 <- t(apply(fit[, c("chisq", "cfi", "srmr", "rmsea")], 2, quantile, probs=probs))
  colnames(Quantile2) <- paste0("Quantile: ", colnames(Quantile2))

  fitindices <- cbind(Mean=apply(fit[, c("chisq", "cfi", "srmr", "rmsea")], 2, mean, na.rm=TRUE),
                      SD=apply(fit[, c("chisq", "cfi", "srmr", "rmsea")], 2, sd, na.rm=TRUE),
                      Quantile=Quantile2)

  chisq.df <- mean(fit[, "df"], na.rm=TRUE)
  pvalue.reject <- mean(ifelse(fit[, "pvalue"] < cutoff.chisq.pvalue, yes=1, no=0), na.rm=TRUE)*100
  cfi <- mean(ifelse(fit[, "cfi"] > cutoff.CFI, yes=1, no=0), na.rm=TRUE)*100
  srmr <- mean(ifelse(fit[, "srmr"] < cutoff.SRMR, yes=1, no=0), na.rm=TRUE)*100
  rmsea <- mean(ifelse(fit[, "rmsea"] < cutoff.RMSEA, yes=1, no=0), na.rm=TRUE)*100
  fitsummary <- list(chisq.df=chisq.df, pvalue.reject=pvalue.reject, cfi=cfi, srmr=srmr,
                     rmsea=rmsea, cutoff.chisq.pvalue=cutoff.chisq.pvalue, cutoff.CFI=cutoff.CFI,
                     cutoff.SRMR=cutoff.SRMR, cutoff.RMSEA=cutoff.RMSEA, Rep.total=Rep.total,
                     Rep.success=length(object))

  out <- list(coefficients=coefficients, fitindices=fitindices, fitsummary=fitsummary)
  class(out) <- "summary.bootuniR2"
  out
}

print.summary.bootuniR2 <- function(x, ...) {
  if (!is.element("summary.bootuniR2", class(x)))
    stop("\"x\" must be an object of class \"summary.bootuniR2\".")

  cat("Summary of the coefficients:\n")
  printCoefmat(x$coefficients)

  cat("\nSummary of the goodness-of-fit indices:\n")
  printCoefmat(x$fitindices)

  cat("\nChisq test (df): ", x$fitsummary$chisq.df)
  cat("\nPercentage of rejection (", x$fitsummary$cutoff.chisq.pvalue, "): ", 
                                     x$fitsummary$pvalue.reject)
  cat("\nPercentage of CFI >", x$fitsummary$cutoff.CFI, ": ", x$fitsummary$cfi)
  cat("\nPercentage of SRMR <", x$fitsummary$cutoff.SRMR, ": ", x$fitsummary$srmr)
  cat("\nPercentage of RMSEA <", x$fitsummary$cutoff.RMSEA, ": ", x$fitsummary$rmsea)
  cat("\nNumber of total replications:", x$fitsummary$Rep.total)
  cat("\nNumber of successful replications:", x$fitsummary$Rep.success, "\n")
}
