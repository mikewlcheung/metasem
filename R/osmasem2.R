osmasem2 <- function(model.name="osmasem2", RAM, data, cor.analysis=TRUE,
                     RE.type.Sigma=c("Diag", "Symm", "Zero"),
                     RE.type.Mu=c("Diag", "Symm", "Zero"),
                     RE.type.SigmaMu=c("Zero", "Full"),
                     mean.analysis=FALSE, intervals.type=c("z", "LB"),
                     startvalues=NULL, replace.constraints=FALSE,
                     mxModel.Args=NULL, run=TRUE, ...) {

  intervals.type <- match.arg(intervals.type, c("z", "LB"))
  RE.type.Sigma <- match.arg(RE.type.Sigma, c("Diag", "Symm", "Zero"))
  RE.type.Mu <- match.arg(RE.type.Mu, c("Diag", "Symm", "Zero"))
  RE.type.SigmaMu <- match.arg(RE.type.SigmaMu, c("Zero", "Full"))
  
  ## Names of the observed variables
  obslabels <- data$obslabels
  no.p <- length(obslabels)
  ## No. of correlation/covariance elements
  no.ps <- length(data$ylabels)

  ## Check if it is a correlation or covariance analysis
  if ((cor.analysis) & (no.ps != no.p*(no.p-1)/2)) {
      stop("Lengths of 'data$obslabels' and 'data$ylabels' do not match for 'cor.analysis=TRUE'.\n")   
  } else if ((cor.analysis==FALSE) & (no.ps != no.p*(no.p+1)/2)) {
     stop("Lengths of 'data$obslabels' and 'data$ylabels' do not match for 'cor.analysis=FALSE'.\n")
  }
  
  ## No sampling covariance matrix of the means in the dataset
  if (is.null(data$vMlabels) & mean.analysis) {
    mean.analysis <- FALSE
    warning("'mean.analysis' is 'TRUE' but the data do not contain the means.
'mean.analysis' is ignored.\n")
  } 

  if (mean.analysis & (length(data$vMlabels) != no.p*(no.p+1)/2)) {
    stop("The length of 'data$vMlabels' does not match that of the number of variables.\n")
  }

  
  #### Create heterogeneity matrix on the covariance structure (no.ps x no.ps)
  labels <- outer(seq_len(no.ps), seq_len(no.ps),
                  function(x, y) {paste0("tau2Cov", x, "_", y)})
  lbound <- matrix(NA, ncol=no.ps, nrow=no.ps)
  ## lbound only applies on the diagonals
  diag(lbound) <- 1e-10
  ## starting values
  values <- matrix(0, ncol=no.ps, nrow=no.ps)
  if (cor.analysis) {
    diag(values) <- 0.2
  } else {
    diag(values) <- 1
  }
  
  switch(RE.type.Sigma,
         Symm = TauCov <- mxMatrix("Symm", ncol=no.ps, nrow=no.ps, free=TRUE,
                                   labels=vech(labels), lbound=vech(lbound),
                                   values=vech(values), name="TauCov"),
         Diag = TauCov <- mxMatrix("Diag", ncol=no.ps, nrow=no.ps, free=TRUE,
                                   labels=Diag(labels), lbound=diag(lbound),
                                   values=diag(values), name="TauCov"),
         Zero = TauCov <- mxMatrix("Zero", ncol=no.ps, nrow=no.ps, free=FALSE,
                                   name="TauCov"))
   
  ## Known sampling covariance matrix
  VCov <- mxMatrix("Symm", ncol=no.ps, nrow=no.ps, free=FALSE,
                   labels=paste0("data.", data$vlabels), name="VCov")
  
  mx.model <- mxModel(model=model.name,
                      mxData(observed=data$data, type="raw"),
                      mxFitFunctionML(), TauCov, VCov)
  
  #### Create heterogeneity matrix on the mean structure
  if (mean.analysis) {

    ## Heterogeneity of the means (no.p x no.p)
    labels <- outer(seq_len(no.p), seq_len(no.p),
                    function(x, y) { paste0("tau2Mean", x, "_", y)})
    lbound <- matrix(NA, ncol=no.p, nrow=no.p)
    diag(lbound) <- 1e-10
    values <- matrix(0, ncol=no.p, nrow=no.p)
    diag(values) <- 0.5
  
    switch(RE.type.Mu,
           Symm = TauMean <- mxMatrix("Symm", ncol=no.p, nrow=no.p, free=TRUE,
                                     labels=vech(labels), lbound=vech(lbound),
                                     values=vech(values), name="TauMean"),
           Diag = TauMean <- mxMatrix("Diag", ncol=no.p, nrow=no.p, free=TRUE,
                                     labels=Diag(labels), lbound=diag(lbound),
                                     values=diag(values), name="TauMean"),
           Zero = TauMean <- mxMatrix("Zero", ncol=no.p, nrow=no.p, free=FALSE,
                                     name="TauMean"))

    ## Known sampling covariance matrix of the means
    VMean <- mxMatrix("Symm", ncol=no.p, nrow=no.p, free=FALSE,
                      labels=paste0("data.", data$vMlabels), name="VMean")

    ## Heterogeneity covariance between Cov and Means: A no.p x no.ps rectangle matrix 
    switch(RE.type.SigmaMu,
           Zero = TauCovMean <- mxMatrix("Full", nrow=no.p, ncol=no.ps,
                                         free=FALSE, name="TauCovMean"),
           Full = TauCovMean <- mxMatrix("Full", nrow=no.p, ncol=no.ps,
                                         free=TRUE, values=0,
                                         labels=
                                           outer(data$obslabels,
                                                 data$ylabels,
                                                 function(y, z)
                                                 {paste0("tau2CM",y,"_",z)}),
                                         name="TauCovMean"))
    
    ## Check if the sampling covariances (a no.p x no.ps rectangle) is present in the data
    ## If not, create one with zeros.
    if (is.null(data$VyMlabels)) {
      VCovMean <- mxMatrix("Zero", nrow=no.p, ncol=no.ps, name="VCovMean")
    } else {
      VCovMean <- mxMatrix("Full", nrow=no.p, ncol=no.ps, free=FALSE, values=0,
                           labels=paste0("data.", data$VyMlabels),
                           name="VCovMean")
    }
    
    ## Heterogeneity matrix of everything
    TauTotal <- mxAlgebra(rbind(cbind(TauCov, t(TauCovMean)),
                                cbind(TauCovMean, TauMean)), name="TauTotal")

    ## Known sampling covariance matrix of everything
    VTotal <- mxAlgebra(rbind(cbind(VCov, t(VCovMean)),
                              cbind(VCovMean, VMean)), name="VTotal")

    ## Expected covariance matrix=Tau + V
    expCov <- mxAlgebra(TauTotal + VTotal, name="expCov")
    
    ## expMean: the final mean of both Cov and Means or just Cov
    mx.model <- mxModel(mx.model,
                        mxExpectationNormal(covariance='expCov',
                                           means='expMean',
                                           dimnames=c(data$ylabels,
                                                      data$obslabels)),
                        TauMean, TauCovMean, TauTotal,
                        VMean, VCovMean, VTotal, expCov)       
  } else {
    ## No mean structure, only cov
    TauTotal <- mxAlgebra(TauCov, name="TauTotal")
    VTotal <- mxAlgebra(VCov, name="VTotal")
    
    ## Expected covariance matrix=Tau + V
    expCov <- mxAlgebra(TauTotal + VTotal, name="expCov")

    mx.model <- mxModel(mx.model, TauTotal, VTotal, expCov,
                        mxExpectationNormal(covariance='expCov',
                                            means='expMean',
                                            dimnames=data$ylabels))
  }
 
    
  #### Prepare saturated and independent models
  if (cor.analysis) {
    ## Independent model
    SigmaInd <- mxMatrix("Iden", ncol=no.p, nrow=no.p, name="SigmaInd")
    vecSigmaInd <- mxAlgebra(t(vechs(SigmaInd)), name="vecSigmaInd")
    ## Saturated model
    SigmaSat <- mxMatrix("Stan", ncol=no.p, nrow=no.p, free=TRUE,
                         name="SigmaSat")
    vecSigmaSat <- mxAlgebra(t(vechs(SigmaSat)), name="vecSigmaSat") 
  } else {
    ## Analysis of covariance matrices
    ## Independent model
    SigmaInd <- mxMatrix("Diag", ncol=no.p, nrow=no.p, values=1, free=TRUE,
                         name="SigmaInd")
    vecSigmaInd <- mxAlgebra(t(vech(SigmaInd)), name="vecSigmaInd")
    ## Saturated model
    SigmaSat <- mxMatrix("Symm", ncol=no.p, nrow=no.p, free=TRUE,
                         name="SigmaSat")
    vecSigmaSat <- mxAlgebra(t(vech(SigmaSat)), name="vecSigmaSat") 
  }

  ## The mean structure is always Full in the saturated and independent models.
  ## TODO Goodness-of-fit indices may be over-fitted when the mean
  ## structure is included.
  if (mean.analysis) {
    ## Saturated model
    MuFull <- mxMatrix("Full", ncol=no.p, nrow=1, values=0, free=TRUE,
                       name="MuFull")
    mx.sat <- mxModel(mx.model, SigmaSat, vecSigmaSat, MuFull,
                      mxAlgebra(cbind(vecSigmaSat, MuFull), name="expMean"))
    mx.ind <- mxModel(mx.model, SigmaInd, vecSigmaInd, MuFull,
                      mxAlgebra(cbind(vecSigmaInd, MuFull), name="expMean"))
  } else {
    ## No mean structure
    mx.sat <- mxModel(mx.model, SigmaSat, vecSigmaSat,
                      mxAlgebra(vecSigmaSat, name="expMean"))
    mx.ind <- mxModel(mx.model, SigmaInd, vecSigmaInd,
                      mxAlgebra(vecSigmaInd, name="expMean"))
  }

    
  #### Prepare model speciification
  Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
  Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
  Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
  Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
   
  ## Some basic checking in RAM
  checkRAM(Amatrix, Smatrix, cor.analysis=cor.analysis)

  ## Extract all observed and latent variable names
  all.names <- colnames(Fmatrix$values)

  ## Note. startvalues may overwrite the starting values in RAM
  ## Collate the starting values from RAM and add them to startvalues
  para.labels <- c(Amatrix$labels[Amatrix$free], Smatrix$labels[Smatrix$free],
                   Mmatrix$labels[Mmatrix$free])
  para.values <- c(Amatrix$values[Amatrix$free], Smatrix$values[Smatrix$free],
                   Mmatrix$values[Mmatrix$free])
  ## Name the starting values with names, which is consistent with the startvalues
  names(para.values) <- para.labels
  para.values <- as.list(para.values)

  ## Prepare startvalues
  if (is.null(startvalues)) {
    ## Note. startvalues may overwrite the starting values in RAM
    startvalues <- para.values
  } else {
    ## Remove para.values if they are overlapped with startvalues    
    para.values[names(para.values) %in% names(startvalues)] <- NULL
    startvalues <- c(startvalues, para.values)

    ## Replace the startvalues in Amatrix, Smatrix, and Mmatrix
    for (i in seq_along(startvalues)) {
      Amatrix$values[Amatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
      Smatrix$values[Smatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
      Mmatrix$values[Mmatrix$labels==names(startvalues)[i]] <- startvalues[[i]]
    }
  }
  
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
      
      ## Remove the parameters from para.labels so they won't appear in mxCI
      para.labels <- para.labels[!para.labels==form[2]]
    }

    ## Remove the constraints so they won't be added again
    mxalgebras[index] <- NULL
  }

  ## New parameter labels including those in constraints
  new.para.labels <- unique(c(A, S, M))
  ## Get the variable names
  new.para.labels <- all.vars(parse(text=new.para.labels))
  ## Drop the definition variables
  new.para.labels <- new.para.labels[!grepl("data.", new.para.labels)]
 
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

  Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")
  Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
  
  ## Expected covariance matrix of the observed variables
  expSigma <- mxAlgebra(Fmatrix %&% (Id_A %&% Smatrix), name="expSigma")

  ## Apply constraints on the diagonals to ensure unity
  ## if there are free parameters on Smatrix
  if (cor.analysis & sum(Diag(Smatrix@free))) {
    Constraints <- Diag(Smatrix@free)
    ## Use nonlinear constraints to impose diagonals as 1
    One <- mxMatrix("Full", values=1, ncol=1, nrow=sum(Constraints),
                    free=FALSE, name="One")
    select <- create.Fmatrix(Constraints, name="select")
    constraint <- mxConstraint(select %*% diag2vec(Id_A %&% Smatrix) == One, 
                               name="constraint")
    vecSigma <- mxAlgebra(t(vechs(expSigma)), name="vecSigma")

    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expSigma, vecSigma,
                        One, select, constraint)
  } else {
    vecSigma <- mxAlgebra(t(vech(expSigma)), name="vecSigma")
    
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expSigma, vecSigma)
  }


  ## Combine everything to create a single model
  ## Mean structure
  if (mean.analysis) {   
    if (all(M==as.symMatrix(RAM$M))) {
      mx.model <- mxModel(mx.model, Mmatrix)
    } else {
      M <- as.mxAlgebra(M, startvalues=startvalues, name="Mmatrix")
      mx.model <- mxModel(mx.model, M$mxalgebra, M$parameters, M$list)
    }        
    expMu <- mxAlgebra(Mmatrix %*% t(Id_A) %*% t(Fmatrix), name="expMu")
    ## Expected Cov and means combined
    expMean <- mxAlgebra(cbind(vecSigma, expMu), name="expMean")
    mx.model <- mxModel(mx.model, expMu, expMean,
                        mxCI(new.para.labels))
  } else {
    ## No mean structure; only cov or cor
    expMean <- mxAlgebra(vecSigma, name="expMean")
    mx.model <- mxModel(mx.model, expMean,
                        mxCI(new.para.labels))
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

  ## Ensure to provide SEs
  mx.model <- mxOption(mx.model, "Calculate Hessian", "Yes")
  mx.model <- mxOption(mx.model, "Standard Errors", "Yes")
  
  if (run) {
    ## Default is z
    mx.fit <- tryCatch(mxRun(mx.model, intervals=FALSE,
                             suppressWarnings=TRUE, silent=TRUE),
                       error=function(e) e)

    ## Check if any errors
    if (inherits(mx.fit, "error") | !(mx.fit$output$status$code %in% c(0,1))) {
      mx.fit <- mxTryHard(mx.model, extraTries=10, intervals=FALSE, silent=TRUE)                        
      mx.fit <- tryCatch(mxRun(mx.fit, intervals=(intervals.type=="LB"),
                               suppressWarnings=TRUE, silent=TRUE, ...),
                         error=function(e) e)
      if (inherits(mx.fit, "error")) {   
        warning("Error in running mxModel.\n")
      }
    }
  } else {
    mx.fit <- mx.model
  }

  out <- list(mx.fit=mx.fit, mx.ind=mx.ind, mx.sat=mx.sat,
              cor.analysis=cor.analysis, mean.analysis=mean.analysis,
              n=sum(data$n), intervals.type=intervals.type)
  class(out) <- "osmasem2"
  out
}  

coef.osmasem2 <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem2", class(object)))
        stop("\"object\" must be an object of class \"osmasem2\".")
    
    mx.coef <- coef(object$mx.fit)
    ## index for untransformed random effects (not the correct ones!) 
    index <- grep("tau2", names(mx.coef))

    select <- match.arg(select)
    switch(select,
           fixed =  mx.coef <- mx.coef[-index],
           random = mx.coef <- mx.coef[index])
    mx.coef    
}

vcov.osmasem2 <- function(object, select=c("fixed", "all", "random"),
                         robust=FALSE, ...) {
    if (!is.element("osmasem2", class(object)))
        stop("\"object\" must be an object of class \"osmasem2\".")

    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names(omxGetParameters(object$mx.fit))
    my.name <- my.name[!is.na(my.name)]

    ## index for untransformed random effects (not the correct ones!) 
    index.random <- grep("tau2", my.name) 

    select <- match.arg(select)   
    if (length(index.random) != 0) {             
        switch(select,
               fixed =  my.name <- my.name[-index.random],
               random = my.name <- my.name[index.random])
    }
    
    if (robust) {
      out <- suppressMessages(imxRobustSE(object$mx.fit, details=TRUE))$cov
    } else {
      out <- vcov(object$mx.fit)
    }
  
    out[my.name, my.name]
}


anova.osmasem2 <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.osmasem2 <- function(object, fitIndices=FALSE, numObs,
                             robust=FALSE, mxTryHard=FALSE, ...) {
  if (!is.element("osmasem2", class(object)))
    stop("\"object\" must be an object of class \"osmasem2\".")

  if (missing(numObs)) numObs <- object$n
  
  ## Calculate chi-square statistic and other fit indices
  if (fitIndices) {

    ## Turn off some functions to improve the speed
    mx.ind <- object$mx.ind
    mx.sat <- object$mx.sat
    mx.ind <- mxOption(mx.ind, "Calculate Hessian", "No")
    mx.ind <- mxOption(mx.ind, "Standard Errors", "No")
    mx.sat <- mxOption(mx.sat, "Calculate Hessian", "No")
    mx.sat <- mxOption(mx.sat, "Standard Errors", "No")

    if (mxTryHard) {
      Sat.stat <- tryCatch(mxTryHard(mx.sat, silent=TRUE),
                           error = function(e) e)
      Ind.stat <- tryCatch(mxTryHard(mx.ind, silent=TRUE),
                           error = function(e) e)
    } else {
      Sat.stat <- tryCatch(mxRun(mx.sat, silent=TRUE, suppressWarnings=TRUE),
                           error = function(e) e)
      Ind.stat <- tryCatch(mxRun(mx.ind, silent=TRUE, suppressWarnings=TRUE),
                           error = function(e) e)
    }
    
    ## Sat.model=FALSE for saturated model if there are either errors or nonconvergent.
    if (inherits(Sat.stat, "error") | !(Sat.stat$output$status$code %in% c(0,1))) {
      Sat.model <- FALSE
    } else {
      Sat.stat <- summary(Sat.stat)
      if (Sat.stat$degreesOfFreedom <= 0) {
        Sat.model <- FALSE
      } else {
        Sat.model <- TRUE
      }
    }

    ## Ind.model=FALSE for independence model if there are either errors or nonconvergent.
    ## NA is acceptable for independence model as it means that optimization was not attempted.
    if (inherits(Ind.stat, "error") | !(Ind.stat$output$status$code %in% c(NA, 0, 1))) {
      Ind.model <- FALSE
    } else {
      Ind.stat <- summary(Ind.stat)
      if (Ind.stat$degreesOfFreedom <= 0) { 
        Ind.model <- FALSE
      } else {
        Ind.model <- TRUE
      }
    }

    ## If sat.model is defined
    if (Sat.model) {
      if (Ind.model) {
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                       SaturatedDoF=Sat.stat$degreesOfFreedom,
                       IndependenceLikelihood=Ind.stat$Minus2LogLikelihood,
                       IndependenceDoF=Ind.stat$degreesOfFreedom,
                       numObs=numObs)
      } else {
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,numObs,
                       SaturatedDoF=Sat.stat$degreesOfFreedom,
                       numObs=numObs, ...)
        warning("There are errors in either fitting the independence model or its degree of freedom being non-positive.\n")
      }
      ## if Sat.model is not defined, no chi-square and fit indices    
    } else {
      out <- summary(object$mx.fit, numObs=numObs, ...)
      warning("There are errors in either fitting the saturated model or its degree of freedom being non-positive.\n")
    }
    ## fitIndices=FALSE 
  } else {
    out <- summary(object$mx.fit, numObs=numObs, ...)
  }
  
  ## Modified the SE with robust SE
  if (robust) {
    robust.SE <- suppressMessages(imxRobustSE(object$mx.fit))
    out$parameters[, "Std.Error"] <- robust.SE[out$parameters$name]
  }

  ## Additional output to the mx summary
  out$parameters$`z value` <- with(out$parameters, Estimate/Std.Error)
  out$parameters$`Pr(>|z|)` <- 2*(1-pnorm(abs(out$parameters$`z value`)))
  out
}

