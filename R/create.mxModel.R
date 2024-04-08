## FIXME: when there are constraints with replace.constraints=TRUE and
## intervals.type="LB", it returns an error because some parameters are
## replaced with the new parameters in the constraints. However, the
## names of these new parameters are not in the CI object.
create.mxModel <- function(model.name="mxModel", RAM=NULL, data=NULL,
                           Cov=NULL, means=NULL, numObs,
                           intervals.type=c("z", "LB"), startvalues=NULL,
                           replace.constraints=FALSE, mxModel.Args=NULL,
                           run=TRUE, silent=TRUE, ...) {

  intervals.type <- match.arg(intervals.type)
    
  Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
  Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
  Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
  Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
    
  ## Some basic checking in RAM
  checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

  ## Extract all observed and latent variable names
  var.names <- colnames(Fmatrix$values)

  ## Without raw data
  if (is.null(data))  {
    ## Without means
    if (is.null(means)) {
      mx.data <- mxData(observed=Cov, type="cov", numObs=numObs)
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix",
                                 F="Fmatrix", dimnames=var.names)
    } else {
      ## With means
      mx.data <- mxData(observed=Cov, type="cov", means=means, numObs=numObs)
      expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                                 M="Mmatrix", dimnames=var.names)
    }
  } else {
    ## With raw data
    mx.data <- mxData(observed=data, type="raw")
    expFun <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix",
                               M="Mmatrix", dimnames=var.names)
  }    

  ## Create an incomplete model, which will be used to store other mx objects.
  mx.model <- mxModel(model.name, mx.data, expFun, mxFitFunctionML())
  
  ## Collate the starting values from RAM and add them to startvalues
  para.labels <- c(Amatrix$labels[Amatrix$free], Smatrix$labels[Smatrix$free],
                   Mmatrix$labels[Mmatrix$free])
  para.values <- c(Amatrix$values[Amatrix$free], Smatrix$values[Smatrix$free],
                   Mmatrix$values[Mmatrix$free])
  ## Name the starting values with names, which is consistent with the startvalues
  names(para.values) <- para.labels
  para.values <- as.list(para.values)
  ## Remove starting values from para.values if they are overlapped with startvalues    
  para.values[names(para.values) %in% names(startvalues)] <- NULL
  startvalues <- c(startvalues, para.values)
  
  ## Extract a local copy for ease of reference
  ## Remove starting values for ease of matching
  A <- as.symMatrix(RAM$A)
  S <- as.symMatrix(RAM$S)
  M <- as.symMatrix(RAM$M)
  mxalgebras <- RAM$mxalgebras
  
  ## Any names of the constraints == parameters?
  ## If yes, these parameters are replaced by the constraints
  index <- sapply(mxalgebras, function(x) {
    ## Convert R language to a vector string
    ## form[1]: "=="
    ## form[2]: "m"
    ## form[3]: "p1 * cos(p2 * data.x) + p2 * sin(p1 * data.x)"
    form <- as.character(x$formula)
    if (form[1]=="==" & form[2]%in% para.labels) TRUE else FALSE
  })
    
  ############################################# 
  ## Need to replace parameters with mxalgebras,
  ## if any TRUE and replace.constraints==TRUE
  if (any(index) & replace.constraints) {

    ## Extract constraints that needed to be replaced
    mxalgebras.const <- mxalgebras[index]

    for (i in seq_along(mxalgebras.const)) {
      form <- as.character(mxalgebras.const[[i]]$formula)

      ## Replace the A matrix
      if (any(grep(form[2], A))) {
        A[which(form[2]==A)] <- form[3]
      }

      ## Replace the S matrix
      if (any(grep(form[2], S))) {
        S[which(form[2]==S)] <- form[3]
      }

      ## Replace the M matrix
      if (any(grep(form[2], M))) {
        M[which(form[2]==M)] <- form[3]
      }
    }

    ## Remove the constraints so they won't be added again
    mxalgebras[index] <- NULL
  }

  ## Check whether there are replacements
  ## Remove the starting values before comparisons
  if (all(A==as.symMatrix(RAM$A))) {
    mx.model <- mxModel(mx.model, Amatrix)
  } else {
    A <- as.mxAlgebra(A, startvalues=startvalues, name="Amatrix")
    mx.model <- mxModel(mx.model, A$mxalgebra, A$parameters, A$list)
  }

  if (all(S==as.symMatrix(RAM$S))) {
    mx.model <- mxModel(mx.model, Smatrix)
  } else {
    S <- as.mxAlgebra(S, startvalues=startvalues, name="Smatrix")
    mx.model <- mxModel(mx.model, S$mxalgebra, S$parameters, S$list)
  }
         
  ## Create an identity matrix from the no. of columens of Fmatrix,
  ## including all latent and observed variables
  Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")

  ## Expected covariance matrix and means of the observed and latent variables
  Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
  expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")   
    
  ## Add the mean structure only if there are means
  if (!is.null(data) | !is.null(means)) {

    if (all(M==as.symMatrix(RAM$M))) {
      mx.model <- mxModel(mx.model, Mmatrix)
    } else {
      M <- as.mxAlgebra(M, startvalues=startvalues, name="Mmatrix")
      mx.model <- mxModel(mx.model, M$mxalgebra, M$parameters, M$list)
    }
        
    expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, expMean,
                        mxCI(c("Amatrix", "Smatrix", "Mmatrix")))
  } else {
    ## No mean structure
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expCov, 
                        mxCI(c("Amatrix", "Smatrix")))
  }
  
  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  ## A list of mxalgebras required SE or CI
  mxalgebras.ci <- NULL
  
  ## Add mxAlgebra and mxConstraint from RAM$mxalgebra
  if (!is.null(mxalgebras)) {
    for (i in seq_along(mxalgebras)) {
      mx.model <- mxModel(mx.model, mxalgebras[[i]])
      ## Name of the mxalgebra
      name.mxalgebra <- names(mxalgebras)[i]
      ## Check if the name constains constraint1, constraint2, ...,
      ## If no, they are mxalgebra, not mxconstraints. Include them in mxCI. 
      if (!grepl("^constraint[0-9]", name.mxalgebra)) {
        mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
        mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
      }
    }
  }
    
  if (run==FALSE) return(mx.model)
  
  ## Default is z
  mx.fit <- tryCatch(mxRun(mx.model, intervals=(intervals.type=="LB"),
                           suppressWarnings=TRUE, silent=TRUE, ...),
                     error=function(e) e)

  ## Check if any errors
  if (inherits(mx.fit, "error")) {
    mx.fit <- mxTryHard(mx.model, extraTries=50, intervals=FALSE, silent=TRUE)                        
    mx.fit <- tryCatch(mxRun(mx.fit, intervals=(intervals.type=="LB"),
                             suppressWarnings=TRUE, silent=TRUE, ...),
                       error=function(e) e)
    if (inherits(mx.fit, "error")) {   
      warning("Error in running mxModel.\n")
    }
  }

  out <- list(mx.fit=mx.fit, RAM=RAM, data=data, mxalgebras=mxalgebras.ci,
              intervals.type=intervals.type)
  class(out) <- "mxRAMmodel"
  out
}

summary.mxRAMmodel <- function(object, robust=FALSE, ...) {
  if (!is.element("mxRAMmodel", class(object)))
    stop("\"object\" must be an object of class \"mxRAMmodel\".")

  # calculate coefficients    
  my.mx <- summary(object$mx.fit)
  ## Exclude lbound ubound etc
  my.para <- my.mx$parameters[, 1:6, drop=FALSE]   

  # Determine if CIs on parameter estimates are present
  if (object$intervals.type=="z") {

    ## Replace the SEs with robust SEs
    if (robust) {
      my.robust <- suppressMessages(imxRobustSE(object$mx.fit))
      my.para[, "Std.Error"] <- my.robust[my.para$name]
    }
        
    my.para$lbound <- with(my.para, Estimate - qnorm(.975)*Std.Error)
    my.para$ubound <- with(my.para, Estimate + qnorm(.975)*Std.Error)
    coefficients <- my.para[, -c(1:4), drop=FALSE]
    dimnames(coefficients)[[1]] <- my.para$name
      
  } else {
 
    ## Convert a data frame with length of 0 in my.mx$CI and remove the last column "note"
    my.ci <- my.mx$CI
    if (length(my.ci)==0) my.ci <- NULL else my.ci <- my.ci[, 1:3, drop=FALSE]        
    
    ## Select the elements matched my.para (excluded I2)  
    my.ci <- my.ci[row.names(my.ci) %in% my.para$name, ]
      
    my.ci <- data.frame(name=row.names(my.ci), my.ci)
    my.para <- merge(my.para, my.ci, by=c("name"))      
    coefficients <- my.para[, -c(1:4,8)]
    dimnames(coefficients)[[1]] <- my.para$name
    # NA for LBCI
    coefficients$Std.Error <- NA
  }
  
  coefficients$"z value" <- coefficients$Estimate/coefficients$Std.Error
  coefficients$"Pr(>|z|)" <- 2*(1-pnorm(abs(coefficients$"z value")))

  informationCriteria <- my.mx$informationCriteria
  ## Better column names
  colnames(informationCriteria) <- c("df Penalty", "Parameters Penalty",
                                     "Sample-Size Adjusted")

  ## Get the mxalgebras
  if (!is.null(object$mxalgebras)) {
    if (object$intervals.type=="z") {
      estimate <- eval(parse(text=paste0("mxEval(rbind(",
                                         paste(object$mxalgebras, collapse=","),
                                         "), object$mx.fit)")))
      SE <- eval(parse(text=paste0("mxSE(rbind(",
                                   paste(object$mxalgebras, collapse=","),
                                   "), model=object$mx.fit, silent=TRUE)")))
      mxalgebras <- cbind(lbound=estimate - 1.96*SE,
                          estimate=estimate,
                          ubound=estimate + 1.96*SE)
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    } else {
      my.ci <- my.mx$CI
      index <- NULL
      for (i in seq_along(object$mxalgebras)) {
        ## Get the names of the mxalgebras combined with the model name
        index <- c(index, grep(paste(object$mx.fit$name, object$mxalgebras[i], sep="."),
                               rownames(my.mx$CI)))
      }
      mxalgebras <- my.mx$CI[index, c("lbound", "estimate", "ubound")]
      dimnames(mxalgebras) <- list(object$mxalgebras,
                                   c("lbound", "estimate", "ubound"))
    }
  } else {
    mxalgebras <- NULL
  }

  out <- list(coefficients=coefficients, mxalgebras=mxalgebras,
              intervals.type=object$intervals.type,
              robust=robust, no.studies=my.mx$numObs,
              obsStat=my.mx$observedStatistics,
              estPara=my.mx$estimatedParameters, df=my.mx$degreesOfFreedom,
              Minus2LL=my.mx$Minus2LogLikelihood,
              Mx.status1=object$mx.fit@output$status[[1]],
              informationCriteria=informationCriteria)
    class(out) <- "summary.mxRAMmodel"
    out
}

print.summary.mxRAMmodel <- function(x, ...) {
    if (!is.element("summary.mxRAMmodel", class(x))) {
      stop("\"x\" must be an object of class \"summary.mxRAMmodel\".")
    }
    
    cat("95% confidence intervals: ")
    switch(x$intervals.type,
           z = cat("z statistic approximation (robust=", x$robust, ")", sep=""),
           LB = cat("Likelihood-based statistic") )

    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, ...)

    if (!is.null(x$mxalgebras)) {
      cat("\nMxalgebras:\n")
      print(x$mxalgebras)
    }
    
    cat("\nInformation Criteria:\n")
    print(x$informationCriteria)
   
    cat("\nNumber of subjects (or studies):", x$no.studies)
    cat("\nNumber of observed statistics:", x$obsStat)
    cat("\nNumber of estimated parameters:", x$estPara)
    cat("\nDegrees of freedom:", x$df)
    cat("\n-2 log likelihood:", x$Minus2LL, "\n")        
    cat("OpenMx status1:", x$Mx.status1, "(\"0\" or \"1\": The optimization is considered fine.\nOther values may indicate problems.)\n")

    if (!(x$Mx.status1 %in% c(0,1))) warning("OpenMx status1 is neither 0 or 1. You are advised to 'rerun' it again.\n")
}

coef.mxRAMmodel <- function(object, ...) {
  if (!is.element("mxRAMmodel", class(object)))
    stop("\"object\" must be an object of class \"mxRAMmodel\".")

  coef(object$mx.fit)
}

vcov.mxRAMmodel <- function(object, robust=FALSE, ...) {
  if (!is.element("mxRAMmodel", class(object)))
    stop("\"object\" must be an object of class \"mxRAMmodel\".")
  
  if (robust) {
    suppressMessages(imxRobustSE(object$mx.fit, details=TRUE)$cov)
  } else {
    vcov(object$mx.fit)
  } 
}

anova.mxRAMmodel <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

plot.mxRAMmodel <- function(x, manNames=NULL, latNames=NULL,
                            labels=c("labels", "RAM"), what="est", nCharNodes=0,
                            nCharEdges=0, layout=c("tree", "circle", "spring",
                                                   "tree2", "circle2"),
                            sizeMan=8, sizeLat=8, edge.label.cex=1.3,
                            color="white", weighted=FALSE, ...) {

  if (!requireNamespace("semPlot", quietly=TRUE))    
    stop("\"semPlot\" package is required for this function.")
    
  if (!inherits(x, "mxRAMmodel"))
    stop("'mxRAMmodel' object is required.\n")
  
  A <- x$mx.fit@matrices$Amatrix$values
  S <- x$mx.fit@matrices$Smatrix$values
  F <- x$mx.fit@matrices$Fmatrix$values
  M <- x$mx.fit@matrices$Mmatrix$values
  RAM <- x$RAM

  ## When there are definition variables, data in the first role are used in
  ## the output. Better to replace it with their means.
  for (i in seq_len(nrow(S))) 
    for (j in seq_len(ncol(S))) {
      if (grepl("data.", RAM$S[i, j])) {
        tmp <- strsplit(RAM$S[i, j], "data.", fixed=TRUE)[[1]][2]
        S[i, j] <- eval(parse(text=paste0("mean(x$data$", tmp, ", na.rm=TRUE)")))
      }
    }     

  for (i in seq_len(nrow(A))) 
    for (j in seq_len(ncol(A))) {
      if (grepl("data.", RAM$A[i, j])) {
        tmp <- strsplit(RAM$A[i, j], "data.", fixed=TRUE)[[1]][2]
        A[i, j] <- eval(parse(text=paste0("mean(x$data$", tmp, ", na.rm=TRUE)")))
      }
    }     
 
  for (j in seq_len(ncol(M))) {
    if (grepl("data.", RAM$M[1, j])) {
      tmp <- strsplit(RAM$M[1, j], "data.", fixed=TRUE)[[1]][2]
      M[1, j] <- eval(parse(text=paste0("mean(x$data$", tmp, ", na.rm=TRUE)")))
    }
  }    
  
  ## index of observed variables
  index_obs <- (apply(F, 2, sum)==1)
  allNames <- colnames(A)
  manNames <- allNames[index_obs]
  latNames <- allNames[!index_obs]
        
  sem.plot <- semPlot::ramModel(A=A, S=S, F=F, M=M, manNames=manNames, latNames=latNames)

  invisible( semPlot::semPaths(sem.plot, what=what, nCharNodes=nCharNodes,
                               nCharEdges=nCharEdges, layout=match.arg(layout),
                               sizeMan=sizeMan, sizeLat=sizeLat,
                               edge.label.cex=edge.label.cex, color=color,
                               weighted=weighted, ...) )
}
