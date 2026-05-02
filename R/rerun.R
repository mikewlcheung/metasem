#' Rerun models via mxTryHard()
#' 
#' It reruns models via mxTryHard().
#' 
#' 
#' @param object An object of either class \code{tssem1FEM}, class
#' \code{tssem1REM}, class \code{wls}, class \code{meta}, class \code{reml},
#' class \code{osmasem}, class \code{osmasem3L}, and class \code{MxModel}.
#' @param autofixtau2 Logical. Whether automatically fixes elements of tau2
#' with NA of standard errors. It only works for objects of class
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
rerun <- function(object, autofixtau2=FALSE, extraTries=10, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta",
                                      "meta3LFIML", "reml",
                                      "tssem1FEM.cluster", "wls.cluster",
                                      "osmasem", "osmasem2", "osmasem3L",
                                      "MxModel", "sem")))
    stop("'object' must be an object of neither class 'meta', 'meta3LFIML',
'wls', 'reml', 'tssem1FEM', 'tssem1REM', 'tssem1FEM.cluster', 'wls.cluster',
'osmasem', 'osmasem2', 'osmasem3L', 'MxModel', or 'mxRAMModel'.")

  ## Run a rerun without autofixtau2 to minimize over-fixing
  ## Many of the NA in SEs may disappear after rerunning it.
  if (autofixtau2 & is.element(class(object)[1], c("tssem1REM", "meta",
                                                   "osmasem", "osmasem2"))) {
    object <- rerun(object, autofixtau2=FALSE, extraTries=extraTries, ...)
  }
  
  ## Automatically fix the problematic Tau2 into 0 for tssem1REM and meta objects
  ## osmasem2 object is similar to meta object as it uses Tau2
  if (autofixtau2 & is.element(class(object)[1], c("tssem1REM", "meta",
                                                   "osmasem2"))) {
    ## Get the Tau2 with NA is SE
    tau2nan <- suppressWarnings(sqrt(diag(vcov(object, select="random"))))
    
    ## Check if tau2nan exists; otherwise names(tau2nan) returns error        
    if (any(tau2nan <- is.nan(tau2nan))) {
      ## Keep elements with NAN (TRUE)
      tau2nan <- tau2nan[tau2nan]
      ## Fix the Tau2 at 0
      object$mx.fit <- omxSetParameters(object$mx.fit, names(tau2nan),
                                        free=FALSE, values=0)
    }
  }

  ## Automatically fix the problematic Tau2 into 0 for osmasem objects
  if (autofixtau2 & is.element(class(object)[1], "osmasem")) {
    ## Check if there are negative variances
    tau2nan <- suppressWarnings(sqrt(diag(vcov(object, select="random"))))

    if (any(tau2nan <- is.nan(tau2nan))) {
      tau2nan <- tau2nan[tau2nan]
      
      ## FIXME: Need to check how robust it is for "Symm" and "User"
      ## Dirty fix to check whether the transformation is expLog or sqSD
      ## Remove all white spaces
      my.transform <- gsub("\\s", "", deparse(object$mx.fit$Tau2$formula)) 
      if (grep("exp", my.transform)) {
        ## exp(-Inf)=0
        values = -Inf
        ## Assuming sdSD
      } else {
        ## 0^2=0
        values = 0
      }
      object$mx.fit <- omxSetParameters(object$mx.fit, names(tau2nan),
                                        free=FALSE, values=values)
    }
  }    
   
  ## Cluster related objects
  if (is.element(class(object)[1], c("tssem1FEM.cluster", "wls.cluster"))) {
    out <- lapply(object, rerun, extraTries=extraTries, ...)
    class(out) <- class(object)[1]
    ## Pure MxModel objects    
  } else if (is.element(class(object)[1], "MxModel")) {
    out <- mxTryHard(object, extraTries=extraTries, greenOK=TRUE, paste=FALSE,
                     bestInitsOutput=FALSE, ...)
    ## Other metaSEM objects    
  } else {

    out <- object
    object$mx.fit <- mxOption(object$mx.fit, "Calculate Hessian", "Yes")
    object$mx.fit <- mxOption(object$mx.fit, "Standard Errors", "Yes")
    
    ## No LB option
    if (is.null(object$intervals.type)) {
      out$mx.fit <- mxTryHard(object$mx.fit, extraTries=extraTries,
                              greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE,
                              ...)
    } else {
      switch(object$intervals.type,
             z = out$mx.fit <- mxTryHard(object$mx.fit, extraTries=extraTries,
                                         greenOK=TRUE, paste=FALSE,
                                         bestInitsOutput=FALSE, ...),
             LB =out$mx.fit <- mxTryHard(object$mx.fit, extraTries=extraTries,
                                         greenOK=TRUE, paste=FALSE,
                                         bestInitsOutput=FALSE,
                                         intervals=TRUE, ...))
    }
  }
  ## Run it again as the SEs sometimes disappear
  ## fit <- out$mx.fit
  ## fit <- mxOption(fit, "Calculate Hessian", "Yes")
  ## fit <- mxOption(fit, "Standard Errors", "Yes")
  ## out$mx.fit <- mxRun(fit, silent=TRUE, intervals=(intervals.type=="LB"),
  ##                     suppressWarnings=TRUE)
  out
}
