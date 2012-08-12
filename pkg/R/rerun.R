rerun <- function(object, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta", "reml")))
      stop("\"object\" must be an object of neither class \"meta\", \"wls\", \"reml\", \"tssem1FEM\" or \"tssem1REM\".")
  
  out <- object
  out$mx.fit <- mxRun(object$mx.fit, ...)
  out
}
