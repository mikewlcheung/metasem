rerun <- function(object, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta", "meta3X", "reml", "tssem1FEM.cluster", "wls.cluster", "osmasem", "MxModel")))
    stop("\"object\" must be an object of neither class \"meta\", \"meta3X\", \"wls\", \"reml\", \"tssem1FEM\", \"tssem1REM\", \"tssem1FEM.cluster\", \"wls.cluster\", \"osmasem\", or \"MxModel\".")

  ## Cluster related objects
  if (is.element(class(object)[1], c("tssem1FEM.cluster", "wls.cluster"))) {
      out <- lapply(object, rerun, ...)
      class(out) <- class(object)[1]
  ## Pure MxModel objects    
  } else if (is.element(class(object)[1], "MxModel")) {
      out <- mxTryHard(object, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, ...)
  ## Other metaSEM objects    
  } else {
    out <- object   
    # No LB option
    if (is.null(object$intervals.type)) {
      out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, ...)
    } else {
      switch(object$intervals.type,
             z = out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, ...),
             LB =out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, intervals=TRUE, ...))
    }
  }
  out
}
