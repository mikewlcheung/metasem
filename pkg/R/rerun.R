rerun <- function(object, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta", "meta3X", "reml", "tssem1FEM.cluster", "wls.cluster")))
    stop("\"object\" must be an object of neither class \"meta\", \"meta3X\", \"wls\", \"reml\", \"tssem1FEM\", \"tssem1REM\", \"tssem1FEM.cluster\" or \"wls.cluster\".")

  if (is.element(class(object)[1], c("tssem1FEM.cluster", "wls.cluster"))) {
    out <- lapply(object, rerun, ...)
    class(out) <- class(object)[1]
  } else {
    out <- object   
    # No LB option
    if (is.null(object$intervals.type)) {
      out$mx.fit <- mxRun(object$mx.fit, ...)
    } else {
      switch(object$intervals.type,
             z = out$mx.fit <- mxRun(object$mx.fit, ...),
             LB =out$mx.fit <- mxRun(object$mx.fit, intervals=TRUE, ...))
    }
  }
  out
}
