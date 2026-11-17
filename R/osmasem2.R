## Restrict data created by Cor2DataFrame() to 'subset.variables', kept in the
## order of the variables in the data. The columns of the correlations or
## covariances, the means, and their sampling covariances are selected by the
## names that Cor2DataFrame() gives them. VyMlabels, which has no naming
## rule, is selected by position: osmasem2() fills it into a no.p x no.ps
## matrix by column (rows: obslabels; columns: ylabels).
.subset.Cor2DataFrame <- function(data, subset.variables, cor.analysis) {
  .check.subset.variables(subset.variables, data$obslabels)
  if (anyDuplicated(data$ylabels) || anyDuplicated(data$vlabels) || anyDuplicated(colnames(data$data)))
    stop("Ambiguous correlation or sampling covariance names; recreate the data with Cor2DataFrame().", call.=FALSE)
  obslabels <- data$obslabels[data$obslabels %in% subset.variables]
  if (cor.analysis && length(obslabels) < 2)
    stop("'subset.variables' must contain at least two variables.\n")
  p <- length(data$obslabels)
  index <- matrix(NA_integer_, p, p)
  index[lower.tri(index, diag=!cor.analysis)] <- seq_len(length(data$ylabels))
  index[upper.tri(index)] <- t(index)[upper.tri(index)]
  rows <- match(obslabels, data$obslabels)
  selected <- index[rows, rows, drop=FALSE]
  cols <- if (cor.analysis) vechs(selected) else vech(selected)
  if (anyNA(cols) || length(data$vlabels) != length(data$ylabels)*(length(data$ylabels)+1L)/2L)
    stop("Correlation and sampling covariance labels do not match the data variables.", call.=FALSE)
  vindex <- matrix(NA_integer_, length(data$ylabels), length(data$ylabels))
  vindex[lower.tri(vindex, diag=TRUE)] <- seq_along(data$vlabels)
  vindex[upper.tri(vindex)] <- t(vindex)[upper.tri(vindex)]
  selected.v <- vech(vindex[cols, cols, drop=FALSE])
  if (!is.null(data$VyMlabels)) {
    if (length(data$VyMlabels) != length(data$obslabels)*length(data$ylabels) ||
        anyNA(cols))
      stop("'VyMlabels' does not match the variables in the data.\n")
    VyM <- matrix(data$VyMlabels, nrow=length(data$obslabels))
    data$VyMlabels <- c(VyM[rows, cols, drop=FALSE])
  }
  if (!is.null(data$vMlabels)) {
    midx <- matrix(NA_integer_, p, p)
    midx[lower.tri(midx, diag=TRUE)] <- seq_along(data$vMlabels)
    midx[upper.tri(midx)] <- t(midx)[upper.tri(midx)]
    data$vMlabels <- data$vMlabels[vech(midx[rows, rows, drop=FALSE])]
  }
  data$obslabels <- obslabels
  data$ylabels <- data$ylabels[cols]
  data$vlabels <- data$vlabels[selected.v]
  absent <- setdiff(c(data$ylabels, data$vlabels, data$vMlabels, data$VyMlabels),
                    colnames(data$data))
  if (length(absent))
    stop("The data do not have the columns ", paste(absent[seq_len(min(3, length(absent)))], collapse=", "),
         if (length(absent) > 3) ", ...",
         ". Use data created by Cor2DataFrame() with 'subset.variables'.\n")
  data
}

#' @rdname osmasem
osmasem2 <- function(model.name="osmasem2", RAM, data, cor.analysis=TRUE,
                     RE.type.Sigma=c("Diag", "Symm", "Zero"),
                     RE.type.Mu=c("Symm", "Diag", "Zero"),
                     RE.type.SigmaMu=c("Zero", "Full"),
                     mean.analysis=FALSE, startvalues=NULL,
                     replace.constraints=FALSE, mxModel.Args=NULL,
                     run=TRUE, subset.variables=NULL, ...) {

  dots <- as.list(substitute(list(...)))[-1L]
  if ("intervals.type" %in% names(dots)) {
    interval <- eval(dots[[which(names(dots) == "intervals.type")[1L]]], parent.frame())
    if (!identical(interval, "z"))
      stop("osmasem2() no longer supports likelihood-based intervals; remove 'intervals.type'.", call.=FALSE)
    warning("'intervals.type' is deprecated in osmasem2(); z intervals are used.", call.=FALSE)
  }
  if (is.character(startvalues) && length(startvalues) == 1L && startvalues %in% c("z", "LB"))
    stop("An older positional 'intervals.type' appears in 'startvalues'; name later arguments. osmasem2() no longer supports likelihood-based intervals.", call.=FALSE)

  .check.ram.names(RAM)
  if (anyDuplicated(data$ylabels) || anyDuplicated(colnames(data$data)))
    stop("Data contain duplicate statistic or column names; regenerate them with Cor2DataFrame().", call.=FALSE)
  .check.no.period(data$obslabels, "Observed variable names")
  .check.no.period(colnames(data$data), "Data column names")
  .check.no.period(subset.variables, "subset.variables")
  .check.parameter.names(names(startvalues))

  ## Define a local copy of the function to solve issues in JASP
  `%&%` <- OpenMx::`%&%`
  `vechs` <- OpenMx::`vechs`
  `diag2vec` <- OpenMx::`diag2vec`
  
  ## If the above does not work, try to do it explicitly:
  # `%&%` <- function (x, y) {
  #   return(x %*% y %*% t(x))
  # }
  # `vechs` <- function (x) {
  #   return(x[lower.tri(x, diag = FALSE)])
  # }
  # `diag2vec` <- function (x) {
  #   return(as.matrix(diag(as.matrix(x))))
  # }
  
  
  RE.type.Sigma <- match.arg(RE.type.Sigma, c("Diag", "Symm", "Zero"))
  RE.type.Mu <- match.arg(RE.type.Mu, c("Symm", "Diag", "Zero"))
  RE.type.SigmaMu <- match.arg(RE.type.SigmaMu, c("Zero", "Full"))

  if (!is.null(subset.variables))
    data <- .subset.Cor2DataFrame(data, subset.variables, cor.analysis)
  
  ## Names of the observed variables
  obslabels <- data$obslabels
  no.p <- length(obslabels)
  ## No. of correlation/covariance elements
  no.ps <- length(data$ylabels)

  ## Check if it is a correlation or covariance analysis
  if ((cor.analysis) & (no.ps != no.p*(no.p-1)/2)) {
      stop("Lengths of 'data$obslabels' and 'data$ylabels' do not match that in 'cor.analysis=TRUE'.\n")   
  } else if ((cor.analysis==FALSE) & (no.ps != no.p*(no.p+1)/2)) {
     stop("Lengths of 'data$obslabels' and 'data$ylabels' do not match that in 'cor.analysis=FALSE'.\n")
  }
  
  ## No sampling covariance matrix of the means in the dataset
  if (is.null(data$vMlabels) & mean.analysis) {
    mean.analysis <- FALSE
    warning("'mean.analysis=TRUE' but the data do not contain the means.
'mean.analysis=FALSE' is used.\n")
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
   
  ## Known sampling covariance matrix of correlation/covariance coefficients
  VCov <- mxMatrix("Symm", ncol=no.ps, nrow=no.ps, free=FALSE,
                   labels=paste0("data.", data$vlabels), name="VCov")
  
  mx.model <- mxModel(model=model.name,
                      mxData(observed=data$data, type="raw"),
                      mxFitFunctionML(), TauCov, VCov)
  
  #### Create heterogeneity matrix on the mean structure
  if (mean.analysis) {

    ## Heterogeneity of the means (no.p x no.p)
    labels <- outer(seq_len(no.p), seq_len(no.p),
                    function(x, y) {paste0("tau2Mean", x, "_", y)})
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
    
    ## Check if the sampling covariances between Cov and Means
    ## (a no.p x no.ps rectangle matrix ) is present in the data
    ## If not, create one with zeros.
    if (is.null(data$VyMlabels)) {
      VCovMean <- mxMatrix("Zero", nrow=no.p, ncol=no.ps, name="VCovMean")
    } else {
      VCovMean <- mxMatrix("Full", nrow=no.p, ncol=no.ps, free=FALSE, values=0,
                           labels=paste0("data.", data$VyMlabels),
                           name="VCovMean")
    }
    
    ## Heterogeneity matrix of everything: cor/cov coefficients + means
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
  RAM <- .match.RAM.obs(RAM, obslabels)
  .check.replacement.parameters(RAM, replace.constraints, mxModel.Args)
  Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
  Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
  Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
  Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
  ## Bounds from inequalities in the lavaan syntax (see lavaan2RAM())
  Amatrix <- .addbound(Amatrix, attr(RAM, "lbound"), attr(RAM, "ubound"))
  Smatrix <- .addbound(Smatrix, attr(RAM, "lbound"), attr(RAM, "ubound"))
  Mmatrix <- .addbound(Mmatrix, attr(RAM, "lbound"), attr(RAM, "ubound"))
   
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
    ## Use the starting values in RAM as the starting values
    startvalues <- para.values
  } else {
    ## Note. startvalues from users overwrite the starting values in RAM
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

  if (isTRUE(replace.constraints)) {
    selection <- .select.substitutions(mxalgebras, para.labels, attr(RAM, "lbound"), attr(RAM, "ubound"))
    A <- .apply.replacements(A, selection$replacements)
    S <- .apply.replacements(S, selection$replacements)
    M <- .apply.replacements(M, selection$replacements)
    mxalgebras <- selection$remaining
  }

  ## Preserve variance restrictions when an originally free diagonal is
  ## replaced, including replacements by constants.
  Constraints <- diag(Smatrix$free) | diag(S) != diag(as.symMatrix(RAM$S))

  ## A fixed residual variance does not guarantee that an observed total
  ## variance is one. Include every selected observed variable; identities
  ## and redundant restrictions are removed by the constraint builder.
  if (cor.analysis) {
    filter <- Fmatrix$values
    if (any(Fmatrix$free) || any(rowSums(filter == 1) != 1L) ||
        any(rowSums(filter != 0) != 1L))
      stop("Correlation analysis requires Fmatrix to select individual observed variables.", call.=FALSE)
    Constraints <- Constraints | colSums(filter == 1) > 0L
  }

  ## Build changed matrices with the same constant/algebra handling as sem().
  add.matrix <- function(model, original, symbolic, source) {
    if (all(symbolic == as.symMatrix(source)))
      return(mxModel(model, original))
    if (!length(all.vars(parse(text=symbolic)))) {
      values <- matrix(vapply(symbolic, function(z) eval(str2lang(z)), numeric(1)),
                       nrow=nrow(symbolic), ncol=ncol(symbolic),
                       dimnames=dimnames(symbolic))
      return(mxModel(model, as.mxMatrix(values, name=original$name)))
    }
    replacement <- as.mxAlgebra(symbolic, startvalues=startvalues,
                                lbound=attr(RAM, "lbound"),
                                ubound=attr(RAM, "ubound"), name=original$name)
    mxModel(model, replacement$mxalgebra, replacement$parameters, replacement$list)
  }

  mx.model <- add.matrix(mx.model, Amatrix, A, RAM$A)
  mx.model <- add.matrix(mx.model, Smatrix, S, RAM$S)

  Id <- as.mxMatrix(diag(ncol(Fmatrix$values)), name="Id")
  Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
  
  ## Expected covariance matrix of the observed variables
  expSigma <- mxAlgebra(Fmatrix %&% (Id_A %&% Smatrix), name="expSigma")

  ## Correlation analysis uses the off-diagonal elements only, whether or
  ## not the variances need constraints.
  if (cor.analysis) {
    vecSigma <- mxAlgebra(t(vechs(expSigma)), name="vecSigma")
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expSigma, vecSigma)

    mx.model <- .osmasem2.unit.constraints(mx.model, A, S, Constraints,
                                           data$data, RAM$mxalgebras)
  } else {
    vecSigma <- mxAlgebra(t(vech(expSigma)), name="vecSigma")
    
    mx.model <- mxModel(mx.model, Fmatrix, Id, Id_A, expSigma, vecSigma)
  }


  ## Combine everything to create a single model
  ## Mean structure
  if (mean.analysis) {   
    mx.model <- add.matrix(mx.model, Mmatrix, M, RAM$M)
    expMu <- mxAlgebra(Mmatrix %*% t(Id_A) %*% t(Fmatrix), name="expMu")
    ## Expected Cov and means combined
    expMean <- mxAlgebra(cbind(vecSigma, expMu), name="expMean")
    mx.model <- mxModel(mx.model, expMu, expMean)
  } else {
    ## No mean structure; only cov or cor
    expMean <- mxAlgebra(vecSigma, name="expMean")
    mx.model <- mxModel(mx.model, expMean)
  }

  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  ## Add mxAlgebra and mxConstraint from RAM$mxalgebra
  if (!is.null(mxalgebras)) {
    for (i in seq_along(mxalgebras)) {
      mx.model <- mxModel(mx.model, mxalgebras[[i]])
    }
  }

  #### Dropped mxCI in this version
  
  ## ## A list of mxalgebras required SE or CI
  ## mxalgebras.ci <- NULL
  
  ## ## Add mxAlgebra and mxConstraint from RAM$mxalgebra
  ## if (!is.null(mxalgebras)) {
  ##   for (i in seq_along(mxalgebras)) {
  ##     mx.model <- mxModel(mx.model, mxalgebras[[i]])
  ##     ## Name of the mxalgebra
  ##     name.mxalgebra <- names(mxalgebras)[i]
  ##     ## Check if the name constains constraint1, constraint2, ...,
  ##     ## If no, they are mxalgebra, not mxconstraints. Include them in mxCI. 
  ##     if (!grepl("^constraint[0-9]", name.mxalgebra)) {
  ##       mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
  ##       mxalgebras.ci <- c(mxalgebras.ci, name.mxalgebra)
  ##     }
  ##   }
  ## }

  ## Ensure to provide SEs
  mx.model <- mxOption(mx.model, "Calculate Hessian", "Yes")
  mx.model <- mxOption(mx.model, "Standard Errors", "Yes")
  mx.model <- .tighten.feasibility(mx.model)
  
  if (run) {
    ## Default is z
    mx.fit <- tryCatch(mxRun(mx.model, intervals=FALSE,
                             suppressWarnings=TRUE, silent=TRUE),
                       error=function(e) e)

    ## Check if any errors
    if (inherits(mx.fit, "error") || !(mx.fit$output$status$code %in% c(0,1))) {
      first.failure <- if (inherits(mx.fit, "error")) conditionMessage(mx.fit) else
        paste0("OpenMx status ", mx.fit$output$status$code)
      mx.fit <- tryCatch(.quiet(mxTryHard(mx.model, extraTries=30, intervals=FALSE, silent=TRUE)),
                         error=function(e) e)
      if (inherits(mx.fit, "error"))
        stop("osmasem2() failed: ", first.failure, "; mxTryHard() also failed: ",
             conditionMessage(mx.fit), call.=FALSE)
      extra <- list(...)
      extra$intervals.type <- NULL
      mx.fit <- tryCatch(do.call(mxRun, c(list(mx.fit, intervals=FALSE,
                               suppressWarnings=TRUE, silent=TRUE), extra)),
                         error=function(e) e)
      if (inherits(mx.fit, "error")) {   
        warning("Error in running mxModel.\n")
      }
    }
    mx.fit <- .check.constraints(mx.fit)
    if (inherits(mx.fit, "MxModel"))
      .warn.Tau.not.psd(mxEvalByName("TauTotal", mx.fit, compute=TRUE), "TauTotal")
  } else {
    mx.fit <- mx.model
  }
  
  out <- list(mx.fit=mx.fit, mx.ind=mx.ind, mx.sat=mx.sat, RAM=RAM,
              cor.analysis=cor.analysis, mean.analysis=mean.analysis,
              n=sum(data$n), data=data)
  class(out) <- "osmasem2"
  out
}  

coef.osmasem2 <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem2", class(object)))
        stop("\"object\" must be an object of class \"osmasem2\".")
    
    mx.coef <- coef(object$mx.fit)
    index <- names(mx.coef) %in% .osmasem2.random.labels(object$mx.fit)

    select <- match.arg(select)
    switch(select,
           fixed =  mx.coef <- mx.coef[!index],
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

    index.random <- my.name %in% .osmasem2.random.labels(object$mx.fit)

    select <- match.arg(select)
    switch(select,
           fixed =  my.name <- my.name[!index.random],
           random = my.name <- my.name[index.random])
    
    out <- .osmasem.vcov(object$mx.fit, robust=robust)
  
    out[my.name, my.name, drop=FALSE]
}


## Labels of the free parameters in the heterogeneity matrices
.osmasem2.random.labels <- function(mx.fit) {
  unlist(lapply(c("TauCov", "TauMean", "TauCovMean"), function(x) {
    m <- mx.fit[[x]]
    if (methods::is(m, "MxMatrix")) m$labels[m$free] else NULL
  }))
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
      Sat.stat <- tryCatch(.quiet(mxTryHard(mx.sat, silent=TRUE)),
                           error = function(e) e)
      Ind.stat <- tryCatch(.quiet(mxTryHard(mx.ind, silent=TRUE)),
                           error = function(e) e)
    } else {
      Sat.stat <- tryCatch(mxRun(mx.sat, silent=TRUE, suppressWarnings=TRUE),
                           error = function(e) e)
      Ind.stat <- tryCatch(mxRun(mx.ind, silent=TRUE, suppressWarnings=TRUE),
                           error = function(e) e)
    }
    
    ## Sat.model=FALSE for saturated model if there are either errors or nonconvergent.
    if (inherits(Sat.stat, "error") || !(Sat.stat$output$status$code %in% c(0,1))) {
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
    if (inherits(Ind.stat, "error") || !(Ind.stat$output$status$code %in% c(NA, 0, 1))) {
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
                       numObs=numObs, ...)
      } else {
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
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
  
  ## Modified the SE with robust SE, or with model-based SE when there are
  ## mxConstraints
  if (robust || .has.constraints(object$mx.fit)) {
    SE <- sqrt(diag(.osmasem.vcov(object$mx.fit, robust=robust)))
    out$parameters[, "Std.Error"] <- SE[out$parameters$name]
  }

  ## Additional output to the mx summary
  out$parameters$`z value` <- out$parameters$Estimate / out$parameters$Std.Error
  out$parameters$`Pr(>|z|)` <- 2*pnorm(-abs(out$parameters$`z value`))
  out
}


## Sparse polynomials in scalar expressions. Non-polynomial calls are opaque
## atoms: cancellation is algebraic, never inferred from a Jacobian at one
## parameter value (which could discard a real, locally flat restriction).
.osmasem.polynomials <- function() {
  clean <- function(x) x[x != 0]
  add <- function(x, y) {
    keys <- union(names(x), names(y))
    z <- setNames(numeric(length(keys)), keys)
    z[names(x)] <- x
    before <- z[names(y)]
    z[names(y)] <- before + y
    ## Cancel roundoff relative to the two summands, never relative to an
    ## unrelated large coefficient (a small new monomial can be meaningful).
    cancelled <- abs(z[names(y)]) <=
      16 * .Machine$double.eps * (abs(before) + abs(y))
    z[names(y)[cancelled]] <- 0
    clean(z)
  }
  multiply <- function(x, y) {
    z <- numeric(0)
    for (i in seq_along(x)) for (j in seq_along(y)) {
      factors <- c(strsplit(names(x)[i], "\034", fixed=TRUE)[[1]],
                   strsplit(names(y)[j], "\034", fixed=TRUE)[[1]])
      key <- paste(sort(factors[nzchar(factors) & factors != "1"]), collapse="\034")
      if (!nzchar(key)) key <- "1"
      term <- setNames(x[i]*y[j], key)
      z <- add(z, term)
      if (length(z) > 2000L) stop("Polynomial expansion too large.")
    }
    z
  }
  constant <- function(x) {
    if (length(x) != 1L || !is.finite(x))
      stop("Polynomial constants must be finite scalars.")
    if (x == 0) numeric(0) else setNames(x, "1")
  }
  product <- multiply
  parse.one <- function(e) {
    if (is.numeric(e)) return(constant(e))
    if (is.call(e)) {
      op <- as.character(e[[1]])
      if (op == "(") return(parse.one(e[[2]]))
      if (op %in% c("+","-")) {
        x <- parse.one(e[[2]])
        if (length(e) == 2L) return(if (op == "-") -x else x)
        y <- parse.one(e[[3]])
        return(add(x, if (op == "-") -y else y))
      }
      if (op == "*") return(product(parse.one(e[[2]]),parse.one(e[[3]])))
      if (op == "/") {
        denominator <- parse.one(e[[3]])
        if (!length(denominator)) stop("Polynomial division by zero.")
        if (length(denominator) == 1L && identical(names(denominator), "1"))
          return(product(parse.one(e[[2]]), constant(1 / denominator[[1]])))
      }
      if (op == "^" && is.numeric(e[[3]]) && e[[3]] >= 0 &&
          e[[3]] <= 10 && e[[3]] == as.integer(e[[3]])) {
        z <- constant(1)
        for (i in seq_len(e[[3]])) z <- product(z,parse.one(e[[2]]))
        return(z)
      }
      if (!length(all.vars(e))) return(constant(eval(e)))
    }
    setNames(1,deparse1(e))
  }
  to.text <- function(x) {
    if (!length(x)) return("0")
    paste(vapply(seq_along(x), function(i) {
      factors <- strsplit(names(x)[i],"\034",fixed=TRUE)[[1]]
      paste(c(sprintf("%.17g",x[i]), paste0("(",factors,")")),collapse="*")
    },character(1)),collapse=" + ")
  }
  list(add=add, multiply=product, constant=constant, parse=parse.one, text=to.text)
}

## For an acyclic RAM, enumerate the finite directed paths to get exact
## polynomial total variances. Cyclic RAMs use the unexpanded inverse below.
.osmasem.variance.polynomials <- function(A, S, p) {
  n <- nrow(A)
  Ap <- lapply(A,function(z) p$parse(str2lang(z)))
  Sp <- lapply(S,function(z) p$parse(str2lang(z)))
  dim(Ap) <- dim(Sp) <- c(n,n)
  T <- matrix(vector("list",n*n),n,n)
  state <- integer(n)
  visit <- function(i) {
    if (state[i] == 1L) stop("Cyclic RAM.")
    if (state[i] == 2L) return(invisible(NULL))
    state[i] <<- 1L
    T[[i,i]] <<- p$constant(1)
    for (k in seq_len(n)) if (length(Ap[[i,k]])) {
      visit(k)
      for (j in seq_len(n)) T[[i,j]] <<-
        p$add(T[[i,j]],p$multiply(Ap[[i,k]],T[[k,j]]))
    }
    state[i] <<- 2L
  }
  for (i in seq_len(n)) visit(i)
  lapply(seq_len(n),function(i) {
    z <- numeric(0)
    for (j in seq_len(n)) for (k in seq_len(n))
      if (length(T[[i,j]]) && length(Sp[[j,k]]) && length(T[[i,k]]))
        z <- p$add(z,p$multiply(p$multiply(T[[i,j]],Sp[[j,k]]),T[[i,k]]))
    p$add(z,p$constant(-1))
  })
}

## Build unit-variance restrictions at every distinct observed moderator
## pattern. Simplify through the user's equalities solely to detect identities
## and duplicate restrictions; the fitted parameterization is left intact.
.osmasem2.unit.constraints <- function(model, A, S, candidates, data, equalities) {
  if (!any(candidates)) return(model)
  selection <- tryCatch(.select.substitutions(equalities, unique(c(A,S))),
                        error=function(e) list(replacements=NULL))
  diagnostic.A <- .apply.replacements(A,selection$replacements)
  diagnostic.S <- .apply.replacements(S,selection$replacements)
  defvars <- grep("^data\\.",all.vars(parse(text=c(A,S,diagnostic.A,diagnostic.S))),value=TRUE)
  columns <- sub("^data\\.","",defvars)
  patterns <- if (length(columns)) unique(data[,columns,drop=FALSE]) else data.frame(.row=1)
  if (anyNA(patterns)) stop("Definition variables in correlation constraints must not be missing.")
  polynomial <- .osmasem.polynomials()
  basis <- list()
  independent <- function(z) {
    if (!length(z)) return(FALSE)
    z <- z / max(abs(z))
    ## Track a separate cancellation scale for every monomial. A global
    ## row tolerance would erase small, genuinely independent restrictions.
    scale <- abs(z)
    for (b in basis) {
      pivot <- names(b)[which.max(abs(b))]
      if (pivot %in% names(z)) {
        adjustment <- -b * z[[pivot]] / b[[pivot]]
        keys <- union(names(scale), names(adjustment))
        accumulated <- setNames(numeric(length(keys)), keys)
        accumulated[names(scale)] <- scale
        accumulated[names(adjustment)] <- accumulated[names(adjustment)] + abs(adjustment)
        scale <- accumulated
        z <- polynomial$add(z, adjustment)
        z <- z[abs(z) > 16 * .Machine$double.eps * scale[names(z)]]
      }
    }
    if (!length(z)) return(FALSE)
    if (length(z) == 1L && identical(names(z),"1"))
      stop("The unit-variance restrictions are incompatible across studies; specify residual variances consistent with the moderator effects.",call.=FALSE)
    ## Normalize new pivots so later rows never divide by roundoff-sized
    ## coefficients. This rank test is independent of the overall row scale.
    basis[[length(basis)+1L]] <<- z / max(abs(z))
    TRUE
  }
  substitute.data <- function(x,row) {
    ## Parentheses preserve negative constants under powers when the AST is
    ## serialized. Keep full double precision, including scale() columns.
    subs <- setNames(lapply(columns,function(v) call("(",as.numeric(row[[v]]))),defvars)
    x[] <- vapply(x,function(z) deparse1(.substitute.symbols(str2lang(z),subs),
      control=c("keepNA","keepInteger","niceNames","digits17")),character(1))
    x
  }
  matrix.text <- function(x) paste0("rbind(",paste(apply(x,1,function(r)
    paste0("cbind(",paste(r,collapse=","),")")),collapse=","),")")
  count <- 0L
  for (row in seq_len(nrow(patterns))) {
    a <- substitute.data(A,patterns[row,,drop=FALSE])
    s <- substitute.data(S,patterns[row,,drop=FALSE])
    da <- substitute.data(diagnostic.A,patterns[row,,drop=FALSE])
    ds <- substitute.data(diagnostic.S,patterns[row,,drop=FALSE])
    check <- tryCatch(.osmasem.variance.polynomials(da,ds,polynomial),error=function(e) NULL)
    actual <- tryCatch(.osmasem.variance.polynomials(a,s,polynomial),error=function(e) NULL)
    for (i in which(candidates)) {
      if (!is.null(check) && !independent(check[[i]])) next
      count <- count+1L
      if (!is.null(actual)) {
        coefficients <- actual[[i]]
        ## Scale up small restrictions so OpenMx does not treat them as
        ## satisfied merely because all coefficients are below its tolerance.
        ## Never scale down: retain the tolerance in variance units.
        if (length(coefficients))
          coefficients <- coefficients / min(1,max(abs(coefficients)))
        expression <- paste0("(",polynomial$text(coefficients),") == 0")
      } else {
        selector <- matrix(as.numeric(seq_len(nrow(A))==i),nrow=1)
        covariance <- paste0("solve(Id - ",matrix.text(a),") %&% ",matrix.text(s))
        expression <- paste0(matrix.text(selector)," %*% (",covariance,") %*% t(",matrix.text(selector),") == 1")
      }
      constraint <- eval(call("mxConstraint",str2lang(expression),
                              name=paste0("unitVariance",count)))
      model <- mxModel(model,constraint)
    }
  }
  model
}
