rerun <- function(object, autofixtau2=FALSE, extraTries=10, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta",
                                      "meta3LFIML", "reml",
                                      "tssem1FEM.cluster", "wls.cluster",
                                      "osmasem", "osmasem3L", "MxModel",
                                      "mxRAMmodel")))
    stop("\"object\" must be an object of neither class \"meta\", \"meta3LFIML\", \"wls\", \"reml\", \"tssem1FEM\", \"tssem1REM\", \"tssem1FEM.cluster\", \"wls.cluster\", \"osmasem\", \"osmasem3L\", \"MxModel\", or \"mxRAMModel\".")

  ## Run a rerun without autofixtau2 to minimize over-fixing
  ## Many of the NA in SEs may disappear after rerunning it.
  if (autofixtau2 & is.element(class(object)[1], c("tssem1REM", "meta", "osmasem"))) {
    object <- rerun(object, autofixtau2=FALSE, extraTries=extraTries, ...)
  }
  
  ## Automatically fix the problematic Tau2 into 0 for tssem1REM and meta objects
  if (autofixtau2 & is.element(class(object)[1], c("tssem1REM", "meta"))) {
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
  out
}
