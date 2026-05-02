## Generate population correlation matrices from uniR1 object


#' Parametric bootstrap on the univariate R (uniR) object
#' 
#' It generates correlation matrices with the parametric bootstrap on the
#' univariate R (uniR) object.
#' 
#' This function implements the parametric bootstrap approach suggested by Yu
#' et al. (2016). It is included in this package for research interests. Please
#' refer to Cheung (2018) for the issues associated with this parametric
#' bootstrap approach.
#' 
#' @param x An object of class 'uniR1'
#' @param Rep Number of replications of the parametric bootstrap
#' @param nonPD.pop If it is \code{replace}, generated non-positive definite
#' matrices are replaced by generated new ones which are positive definite. If
#' it is \code{nearPD}, they are replaced by nearly positive definite matrices
#' by calling \code{Matrix::nearPD()}. If it is \code{accept}, they are
#' accepted.
#' @return An object of the generated correlation matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{rCor}}, \code{\link[metaSEM]{bootuniR2}},
#' \code{\link[metaSEM]{Nohe15}}
#' @references Cheung, M. W.-L. (2018). Issues in solving the problem of effect
#' size heterogeneity in meta-analytic structural equation modeling: A
#' commentary and simulation study on Yu, Downes, Carter, and O'Boyle (2016).
#' \emph{Journal of Applied Psychology}, \bold{103}, 787-803.
#' 
#' Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle, E. H. (2016). The
#' problem of effect size heterogeneity in meta-analytic structural equation
#' modeling.  \emph{Journal of Applied Psychology}, \bold{101}, 1457-1473.
#' @keywords bootuniR
bootuniR1 <- function(x, Rep, nonPD.pop=c("replace", "nearPD", "accept")) {
  if (!is.element("uniR1", class(x)))
      stop("'x' must be an object of class 'uniR1'.")
  
  rCorPop(Sigma=x$r.mean, V=diag(vechs(x$r.SD^2)),
          corr=TRUE, k=Rep, nonPD.pop=nonPD.pop)
}



#' Fit Models on the bootstrapped correlation matrices
#' 
#' It fits structural equation models on the bootstrapped correlation matrices.
#' 
#' This function fits the lavaan model with the bootstrapped correlation
#' matrices. It implements the parametric bootstrap approach suggested by Yu et
#' al. (2016). It is included in this package for research interests. Please
#' refer to Cheung (2018) for the issues associated with this parametric
#' bootstrap approach.
#' 
#' @param model A model in \code{\link[lavaan]{sem}} syntax.
#' @param data A list of correlation matrices.
#' @param n Sample size in fitting the structural equation models
#' @param \dots Further arguments to be passed to \code{\link[lavaan]{sem}}.
#' @return A list of the fitted object from \code{\link[lavaan]{sem}}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{bootuniR2}},
#' \code{\link[metaSEM]{tssemParaVar}}, \code{\link[metaSEM]{Nohe15}}
#' @references Cheung, M. W.-L. (2018). Issues in solving the problem of effect
#' size heterogeneity in meta-analytic structural equation modeling: A
#' commentary and simulation study on Yu, Downes, Carter, and O'Boyle (2016).
#' \emph{Journal of Applied Psychology}, \bold{103}, 787-803.
#' 
#' Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle, E. H. (2016). The
#' problem of effect size heterogeneity in meta-analytic structural equation
#' modeling.  \emph{Journal of Applied Psychology}, \bold{101}, 1457-1473.
#' @keywords bootuniR
bootuniR2 <- function(model, data, n, ...) {
  ## if (!requireNamespace("lavaan", quietly=TRUE))
  ##   stop("\"lavaan\" package is required for this function.")

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
