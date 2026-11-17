#' Rerun models via mxTryHard()
#'
#' It reruns models via mxTryHard(). Subclasses of the supported classes
#' are accepted; subclass metadata is preserved for clustered lists.
#'
#'
#' @param object An object of either class \code{tssem1FEM}, class
#' \code{tssem1REM}, class \code{wls}, class \code{meta}, class \code{reml},
#' class \code{osmasem}, class \code{osmasem3L}, and class \code{MxModel}.
#' @param autofixtau2 Logical. Whether automatically fixes elements of tau2
#' with NaN standard errors after a retry. Plain NA standard errors produce a
#' warning because their cause is unclear and are not automatically fixed.
#' It only works for objects of class
#' \code{tssem1REM}, class \code{meta}, and class \code{osmasem}.
#' @param extraTries The number of attempts to run the model in addition to the
#' first.
#' @param \dots Further arguments to be passed to
#' \code{\link[OpenMx]{mxTryHard}}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords tssem meta osmasem osmasem3L wls
#' @examples
#'
#' \donttest{
#' random1 <- tssem1(Digman97$data, Digman97$n, method="REM", RE.type="Diag")
#' random1_rerun <- rerun(random1)
#' summary(random1_rerun)
#' }
#'
.rerun.tau.zero <- function(formula) {
  if ("exp" %in% all.names(formula, functions=TRUE)) -Inf else 0
}

.rerun.invalid.tau.se <- function(se) {
  nan <- is.nan(se)
  if (any(is.na(se) & !nan))
    warning("Some random-effect standard errors are NA; their cause is unclear, so those parameters were not automatically fixed.", call.=FALSE)
  nan
}

rerun <- function(object, autofixtau2=FALSE, extraTries=10, ...) {
  if (!inherits(object, c("wls", "tssem1FEM", "tssem1REM", "meta",
                                      "meta3LFIML", "reml",
                                      "tssem1FEM.cluster", "wls.cluster",
                                      "osmasem", "osmasem2", "osmasem3L",
                                      "MxModel", "mxsem")))
    stop("'object' must be an object of neither class 'meta', 'meta3LFIML',
'wls', 'reml', 'tssem1FEM', 'tssem1REM', 'tssem1FEM.cluster', 'wls.cluster',
'osmasem', 'osmasem2', 'osmasem3L', 'MxModel', or 'mxsem'.")

  ## Run a rerun without autofixtau2 to minimize over-fixing
  ## Many of the NA in SEs may disappear after rerunning it.
  if (autofixtau2 & inherits(object, c("tssem1REM", "meta",
                                                   "osmasem", "osmasem2"))) {
    object <- rerun(object, autofixtau2=FALSE, extraTries=extraTries, ...)
  }

  ## Automatically fix the problematic Tau2 into 0 for tssem1REM and meta objects
  ## osmasem2 object is similar to meta object as it uses Tau2
  if (autofixtau2 & inherits(object, c("tssem1REM", "meta",
                                                   "osmasem2"))) {
    se_tau2 <- suppressWarnings(sqrt(diag(vcov(object, select="random"))))
    nan_idx <- .rerun.invalid.tau.se(se_tau2)
    if (any(nan_idx)) {
      object$mx.fit <- omxSetParameters(object$mx.fit, names(se_tau2)[nan_idx],
                                        free=FALSE, values=0)
    }
  }

  ## Automatically fix the problematic Tau2 into 0 for osmasem objects
  if (autofixtau2 & inherits(object, "osmasem")) {
    se_tau2 <- suppressWarnings(sqrt(diag(vcov(object, select="random"))))
    nan_idx <- .rerun.invalid.tau.se(se_tau2)

    if (any(nan_idx)) {
      ## exp(-Inf) = 0; an untransformed or squared SD uses zero.
      values <- .rerun.tau.zero(object$mx.fit$Tau2$formula)
      object$mx.fit <- omxSetParameters(object$mx.fit, names(se_tau2)[nan_idx],
                                        free=FALSE, values=values)
    }
  }

  ## Cluster related objects
  if (inherits(object, c("tssem1FEM.cluster", "wls.cluster"))) {
    out <- object
    out[] <- lapply(object, rerun, extraTries=extraTries, ...)
    ## Pure MxModel objects
  } else if (inherits(object, "MxModel")) {
    out <- mxTryHard(object, extraTries=extraTries, greenOK=TRUE, paste=FALSE,
                     bestInitsOutput=FALSE, ...)
    ## Other metaSEM objects
  } else {
    out <- object
    object$mx.fit <- mxOption(object$mx.fit, "Calculate Hessian", "Yes")
    object$mx.fit <- mxOption(object$mx.fit, "Standard Errors", "Yes")
    out$mx.fit <- mxTryHard(object$mx.fit, extraTries=extraTries,
                            greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE,
                            intervals = isTRUE(object$intervals.type == "LB"),
                            ...)
  }
  out
}
