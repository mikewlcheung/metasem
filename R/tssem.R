tssem1FEM <- function(Cov, n, cor.analysis=TRUE, model.name=NULL,
                      cluster=NULL, suppressWarnings=TRUE, silent=TRUE,
                      run=TRUE, ...) {
    
  ## Function to mark diagonals as NA when all correlation coefficients of the variable are NA
  ## assuming matrices are symmetric
  addDiagNA <- function(x) {
      add_NA <- function(y) {
          for(i in 1:nrow(y)) {
              ## All elements in the row except the diagonal are NA
              if (sum(is.na(y[, i]))==(nrow(y)-1)) y[i,i] <- NA
          }
          y
      }
      if (is.list(x)) lapply(x, add_NA) else add_NA(x)
  }
  ## Mark the diagonals as NA when all correlation coefficients are NA
  Cov <- addDiagNA(Cov)

  ## Smooth non-symmetric matrices probably due to rounding
  ## assuming no input error in the matrices
  if (is.list(Cov)) {
      Cov <- lapply(Cov, function(x) {x <- (x+t(x))/2; x})
  } else {
      Cov <- (Cov+t(Cov))/2
  }

  ## Split the data by cluster
  if (!is.null(cluster)) {
    data.cluster <- tapply(Cov, cluster, function(x) {x})
    n.cluster <- tapply(n, cluster, function(x) {x})
    out <- list()
    for (i in seq_along(data.cluster)) {
      ## Need to correct it to tssem1()
      out[[i]] <- tssem1FEM(data.cluster[[i]], n.cluster[[i]],
                           cor.analysis=cor.analysis, model.name=model.name,
                           suppressWarnings=suppressWarnings, ...)
    }
    names(out) <- names(data.cluster)
    class(out) <- "tssem1FEM.cluster"
    out
  } else {
    ## Throw an error when df and n are in different lengths.
    if (length(Cov) != length(n)) stop("Lengths of \"df\" and \"n\" are different.\n") 
      
    ## Check whether all studies have the same dimensions
    my.range <- range(sapply(Cov, function(x) {ncol(x)}))
    if ( !all.equal(my.range[1], my.range[2]) )
      stop("Dimensions of groups are not the same!\n")

    no.groups <- length(Cov)
    no.var <- ncol(Cov[[1]])

    ## Get the original variable names
    original.names <- colnames(Cov[[1]])

    var.names <- paste("x", 1:no.var, sep = "")
    ## Convert variable labels to x1, x2, ...
    Cov <- lapply(Cov, function(x) {dimnames(x) <- list(var.names, var.names); x})
    total.n <- sum(n)

    ## Check positive definiteness of data
    ## Print warning rather than stop
    isPD <- is.pd(Cov)
    if (!all(isPD))
        warning(paste("Group ", (1:no.groups)[!isPD], " is not positive definite.", sep = ""))

    ## Prepare starting values based on covariance matrices
    sv <- .startValues(Cov, cor.analysis = FALSE)
      
    ## matrix of labels; only use the lower triangle
    ## Fixed a bug reported by John Ma when there are more than 120 variables by replacing " " with "_"  
    ps.labels <- outer(1:no.var, 1:no.var, function(x, y) paste0("s", x, "_", y))
    if (cor.analysis==TRUE) {
      S <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=TRUE, values=vechs(cov2cor(sv)),
                    labels=vechs(ps.labels), name="S")
    } else {
      ## Set lower bound on variances
      lbound <- matrix(NA, nrow=no.var, ncol=no.var)
      Diag(lbound) <- 0.00001
      S <- mxMatrix(type="Symm", nrow=no.var, ncol=no.var, free=TRUE, values=vech(sv),
                    labels=vech(ps.labels), lbound=lbound, name="S")
    }
      
    ## Index for missing variables: only check the diagonals only!!!
    miss.index <- lapply(Cov, function(x) { is.na(Diag(x)) })

    ## complete.index <- NULL
    ## for (i in no.groups:1) {
    ##     if (sum(!miss.index[[i]], na.rm = TRUE) == no.var)
    ##         complete.index = i
    ## }
    ## if (is.null(complete.index))
    ##     stop("It is expected that at least one study has complete data!")

    ## if ( sum(!miss.index[[1]], na.rm = TRUE) != no.var )
    ##     stop("There should be no missing data in the first group.")

    for (i in 1:no.groups) {
        no.var.i <- sum(!miss.index[[i]], na.rm = TRUE)
        miss.index.i <- miss.index[[i]]
        Cov.i <- Cov[[i]][!miss.index.i, !miss.index.i]
        ## dimnames(Cov.i) <- list(var.names[!miss.index.i], var.names[!miss.index.i])

        # Prepare matrices for calculations
        if (cor.analysis) {
            if (is.null(model.name)) model.name <- "TSSEM1 Correlation"
            SD <- Diag(paste(sqrt(Diag(sv)), "*sd", i, "_", 1:no.var, sep=""))
            # Fixed a bug when SD is a scalar
            SD <- SD[!miss.index.i, , drop=FALSE]
            SD <- as.mxMatrix(SD)
        } else {
            if (is.null(model.name)) model.name <- "TSSEM1 Covariance"
            SD <- Diag(rep(1, no.var))
            SD <- SD[!miss.index.i, , drop=FALSE]
            SD <- as.mxMatrix(SD)
        }
        # Expected covariance matrix
        expC <- mxAlgebra(SD %&% S, name="expC", dimnames=list(var.names[!miss.index.i],
                          var.names[!miss.index.i]))

        fitFunction <- mxFitFunctionML()

        ## Fixed a bug when Cov[[i]][, , ] is a scalar
        g.model <- paste("g", i, " <- mxModel(\"g", i, "\", S, SD, expC, mxData(observed=Cov[[",i,"]][!miss.index[[",i,"]],!miss.index[[",i,"]], drop=FALSE], type=\"cov\", numObs=n[", i,
                "]), fitFunction, mxExpectationNormal(covariance=\"expC\", means=NA, dimnames=var.names[!miss.index[[", i, "]]]))", sep = "")
        eval(parse(text = g.model))
    }

    ## Replaced mxFitFunctionAlgebra() with mxFitFunctionMultigroup()  
    tssem1.model <- paste0("tssem1 <- mxModel(\"", model.name, "\", S, ",
                          paste0("g", 1:no.groups, collapse = ","),
                          ", mxFitFunctionMultigroup(c(", paste0("\"g", 1:no.groups, "\"", collapse=","), ")))")
    ## mgFitFn <- paste0("mxFitFunctionMultigroup(c(", paste0("\"g", 1:no.groups, "\"", collapse=","), ")")
      
    eval(parse(text = tssem1.model))

    ## Return mx model without running the analysis
    if (run==FALSE) return(tssem1)

    # try to run it with error message as output
    ## mx.fit <- mxRun(tssem1)
    mx.fit <- tryCatch(mxRun(tssem1, suppressWarnings = suppressWarnings, silent=silent, ...),
                       error = function(e) e)

    ## It is useful to return error for computer simulation.  
    if (inherits(mx.fit, "error")) {
        warning(print(mx.fit))
        out <- mx.fit
    } else {
        ## Calculate 2LL of the saturated and independence models and the DF of independence model
        baseMinus2LL <- tryCatch(.minus2LL(x=Cov, n=n), error = function(e) e)

        out <- list(call = match.call(), cor.analysis=cor.analysis, data=Cov, n = n,
                    baseMinus2LL = baseMinus2LL, mx.model=tssem1, mx.fit = mx.fit,
                    original.names=original.names)
        class(out) <- "tssem1FEM"
    }
    return(out)
      
  }
}

tssem1REM <- function(Cov, n, cor.analysis=TRUE, RE.type=c("Diag", "Symm", "Zero", "User"),
                      RE.startvalues=0.1, RE.lbound = 1e-10, RE.constraints=NULL,
                      I2="I2q", acov=c("weighted", "individual", "unweighted"),
                      asyCovOld=FALSE, model.name=NULL, suppressWarnings=TRUE, silent=TRUE,
                      run=TRUE, ...) {
  ## It handles missing effect sizes rather than missing correlations. Thus, it is more flexible than tssem1FEM().
  ## ACOV is calculated without missing data by assuming 1 and 0 for the missing variances and covariances.
  ## Missing values are indicated by the missing effect sizes.

  acov <- match.arg(acov, c("weighted", "individual", "unweighted"))
    
  ## Throw an error when df and n are in different lengths.
  if (length(Cov) != length(n)) stop("Lengths of \"df\" and \"n\" are different.\n")   
    
  ## Get the original variable names
  original.names <- colnames(Cov[[1]])

  ## Calculate the asymptotic sampling covariance matrix of the correlation matrix
  if (asyCovOld) {
      ## SEM approach is sensitive to NA. Additional care is required.      
        
      if (cor.analysis) {
        ## Replace diagonals with 1.0
        my.complete <- lapply(Cov, function (x) { Diag(x)[is.na(Diag(x))] <- 1; x })
      } else {
        ## Replace diagonals with the mean of diagonals
        my.complete <- lapply(Cov, function (x) { Diag(x)[is.na(Diag(x))] <- mean(Diag(x), na.rm=TRUE); x })
      }
      
      ## Replace missing variables with 0.0 regardless of cor.analysis
      my.complete <- lapply(my.complete, function (x) { x[is.na(x)] <- 0; x })
          
      acovR <- asyCovOld(x=my.complete, n=n, cor.analysis=cor.analysis, dropNA=FALSE, as.matrix=TRUE, acov=acov)
  } else {
      ## asyCov() will break down when there are NA
      acovR <- asyCov(x=Cov, n=n, cor.analysis=cor.analysis, as.matrix=TRUE, acov=acov)
  }
    
  ## Fixed a bug that Cov is covariance matrix while cor.analysis is TRUE
  ## When cor.analysis=TRUE, the old version just takes the lower triangle without converting covariance into correlation.
  if (cor.analysis) {
    ## Convert possible covariance matrices into correlation matrices
    ## When there are NA in diagonas, they become 1 after cov2cor()
    ## It is fine as the diagonals are not used in cor.analysis=TRUE
    ES <- list2matrix(x=suppressWarnings(lapply(Cov, cov2cor)), diag=FALSE)
  } else {
    ES <- list2matrix(x=Cov, diag=TRUE)
  }
  ## no. of effect sizes
  no.es <- ncol(ES)

  if (is.null(model.name)) {
    if (cor.analysis) {
      model.name <- "TSSEM1 Correlation"
    } else {
      model.name <- "TSSEM1 Covariance"
    }
  }

  RE.type <- match.arg(RE.type, c("Diag", "Symm", "Zero", "User"))
  switch( RE.type,
         Symm = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.startvalues=RE.startvalues,
                               RE.lbound=RE.lbound, suppressWarnings=TRUE, silent=silent, run=run, ...),
         Diag = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2,
                               RE.constraints=Diag(paste0(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es)),
                               RE.lbound=RE.lbound, suppressWarnings=TRUE, silent=silent, run=run, ...),
         Zero = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.constraints=matrix(0, ncol=no.es, nrow=no.es),
                               suppressWarnings=TRUE, silent=silent, run=run, ...),
         User = mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.constraints=RE.constraints,
                               suppressWarnings=TRUE, silent=silent, run=run, ...) )

  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.fit)

  ## if (RE.diag.only==TRUE) {
  ##   ## No covariance between random effects
  ##   mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2,
  ##                    RE.constraints=diag(x=paste(RE.startvalues, "*Tau2_", 1:no.es, "_", 1:no.es, sep=""),
  ##                                        nrow=no.es, ncol=no.es), RE.lbound=RE.lbound)
  ## } else {
  ##   mx.fit <- meta(y=ES, v=acovR, model.name=model.name, I2=I2, RE.startvalues=RE.startvalues, RE.lbound=RE.lbound)
  ## }

  out <- list(total.n=sum(n), cor.analysis=cor.analysis, RE.type=RE.type, no.es=no.es, original.names=original.names)
  out <- c(out, mx.fit)
  class(out) <- c("tssem1REM", "meta")
  return(out)
}



#' First Stage of the Two-Stage Structural Equation Modeling (TSSEM)
#' 
#' It conducts the first stage analysis of TSSEM by pooling
#' correlation/covariance matrices. \code{tssem1FEM()} and \code{tssem1REM()}
#' use fixed- and random-effects models, respectively. \code{tssem1()} is a
#' wrapper of these functions.
#' 
#' 
#' @aliases tssem1 tssem1FEM tssem1REM
#' @param Cov A list of correlation/covariance matrices
#' @param n A vector of sample sizes
#' @param method Either \code{"REM"} (default if missing) or \code{"FEM"}.  If
#' it is "REM",a random-effects meta-analysis will be applied. If it is "FEM",
#' a fixed-effects meta-analysis will be applied.
#' @param cor.analysis Logical. The output is either a pooled correlation or a
#' covariance matrix.
#' @param cluster A character vector in \code{tssem3L1} and \code{tssemRobust1}
#' or a vector of characters or numbers indicating the clusters in
#' \code{tssem1}. Analyses will be conducted for each cluster. It will be
#' ignored when \code{method="REM"}.
#' @param RE.type Either \code{"Diag"}, \code{"Symm"}, \code{"Zero"} or
#' \code{"User"}. If it is \code{"Diag"} (default if missing), a diagonal
#' matrix is used for the random effects meaning that the random effects are
#' independent. If it is \code{"Symm"}, a symmetric matrix is used for the
#' random effects on the covariances among the correlation (or covariance)
#' vectors. If it is \code{"Zero"}, there is no random effects which is similar
#' to the conventional Generalized Least Squares (GLS) approach to
#' fixed-effects analysis.  \code{"User"}, the user has to specify the variance
#' component via the \code{RE.constraints} argument. This argument will be
#' ignored when \code{method="FEM"}.
#' @param RE.startvalues Starting values on the diagonals of the variance
#' component of the random effects. It will be ignored when
#' \code{method="FEM"}.
#' @param RE.lbound Lower bounds on the diagonals of the variance component of
#' the random effects. It will be ignored when \code{method="FEM"}.
#' @param RE.constraints A \eqn{p*}{p*} x \eqn{p*}{p*} matrix specifying the
#' variance components of the random effects, where \eqn{p*}{p*} is the number
#' of effect sizes. If the input is not a matrix, it is converted into a matrix
#' by \code{as.matrix()}. The default is that all covariance/variance
#' components are free. The format of this matrix follows
#' \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance components can
#' be constrained equally by using the same labels. If a zero matrix is
#' specified, it becomes a fixed-effects meta-analysis.
#' @param I2 Possible options are \code{"I2q"}, \code{"I2hm"} and
#' \code{"I2am"}. They represent the \code{I2} calculated by using a typical
#' within-study sampling variance from the Q statistic, the harmonic mean and
#' the arithmetic mean of the within-study sampling variances (Xiong, Miller, &
#' Morris, 2010). More than one options are possible. If
#' \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity
#' indices will be constructed.
#' @param acov If it is \code{individual}, the sampling variance-covariance
#' matrices are calculated based on individual correlation/covariance matrix.
#' If it is either \code{unweighted} or \code{weighted} (the default), the
#' average correlation/covariance matrix is calculated based on the unweighted
#' or weighted mean with the sample sizes. The average correlation/covariance
#' matrix is used to calculate the sampling variance-covariance matrices. This
#' argument is ignored with the \code{method="FEM"} argument.
#' @param asyCovOld Whether the old \code{asyCov} is used. See
#' \code{\link[metaSEM]{asyCov}}.
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed. It
#' is passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
#' @return Either an object of class \code{tssem1FEM} for fixed-effects TSSEM,
#' an object of class \code{tssem1FEM.cluster} for fixed-effects TSSEM with
#' \code{cluster} argument, or an object of class \code{tssem1REM} for
#' random-effects TSSEM.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{wls}}, \code{\link[metaSEM]{Cheung09}},
#' \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Digman97}},
#' \code{\link[metaSEM]{issp89}}, \code{\link[metaSEM]{issp05}}
#' @references Cheung, M. W.-L. (2014). Fixed- and random-effects meta-analytic
#' structural equation modeling: Examples and analyses in R. \emph{Behavior
#' Research Methods}, \bold{46}, 29-40.
#' 
#' Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation
#' models. \emph{Structural Equation Modeling}, \bold{20}, 429-454.
#' 
#' Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
#' modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10},
#' 40-64.
#' 
#' Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing
#' covariance matrices in meta-analytic structural equation modeling.
#' \emph{Structural Equation Modeling}, \bold{16}, 28-53.
#' @keywords tssem
tssem1 <- function(Cov, n, method=c("REM", "FEM"), cor.analysis=TRUE, cluster=NULL,
                   RE.type=c("Diag", "Symm", "Zero", "User"), RE.startvalues=0.1, RE.lbound=1e-10,
                   RE.constraints=NULL, I2="I2q", acov=c("weighted", "individual", "unweighted"),
                   asyCovOld=FALSE, model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  method <- match.arg(method)
  switch(method,
    FEM = out <- tssem1FEM(Cov=Cov, n=n, cor.analysis=cor.analysis, model.name=model.name,
                          cluster=cluster, suppressWarnings=suppressWarnings, silent=silent, run=run, ...),
    REM = out <- tssem1REM(Cov=Cov, n=n, cor.analysis=cor.analysis, RE.type=RE.type,
                          RE.startvalues=RE.startvalues, RE.lbound=RE.lbound, RE.constraints=RE.constraints,
                          I2=I2, acov=acov, asyCovOld=asyCovOld, model.name=model.name,
                          suppressWarnings=suppressWarnings, silent=silent, run=run, ...) )
  out
}

## Known bug: wls() will fall into loop when the Amatrix is zero


#' Conduct a Correlation/Covariance Structure Analysis with WLS
#' 
#' It fits a correlation or covariance structure with weighted least squares
#' (WLS) estimation method where the inverse of the asymptotic covariance
#' matrix is used as the weight matrix. \code{tssem2} conducts the second stage
#' analysis of the two-stage structural equation modeling (TSSEM).
#' \code{tssem2} is a wrapper of \code{wls}.
#' 
#' 
#' @aliases wls tssem2
#' @param tssem1.obj An object of either class \code{tssem1FEM}, class
#' \code{tssem1FEM.cluster} or class \code{tssem1REM} returned from
#' \code{tssem1()}
#' @param Cov A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
#' where \eqn{p}{p} is the number of variables.
#' @param aCov A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
#' matrix of either \code{\link[OpenMx]{vechs}} \code{(Cov)} or
#' \code{\link[OpenMx]{vech}} \code{(Cov)} where \eqn{p* = p(p-1)/2 }{p* =
#' p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for
#' covariance matrix.
#' @param n Sample size.
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param Amatrix If \code{RAM} is not specified, an \code{Amatrix} is
#' required. An asymmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL}, a matrix of zero
#' will be created. If it is a matrix, it will be converted into
#' \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
#' @param Smatrix If \code{RAM} is not specified, an \code{Smatrix} is
#' required. A symmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Fmatrix A filter matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
#' identity matrix with the same dimensions of \code{Cov} will be created. If
#' it is a matrix, it will be converted into
#' \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function. It
#' is not required when there is no latent variable.
#' @param diag.constraints Logical. This argument is ignored when
#' \code{cor.analysis=FALSE}. If \code{diag.constraints=TRUE}, the diagonals of
#' the model implied matrix would be constrained at 1 by nonlinear constraints.
#' The drawback is that standard error will not be generated. Parametric
#' bootstrap is used to estimate the standard error by drawing samples from
#' \eqn{\mathcal{N}(vech(Cov), asyCov)}{N(vech(Cov), asyCov)} for covariance
#' analysis and \eqn{\mathcal{N}(vechs(Cov), asyCov)}{N(vechs(Cov), asyCov)}
#' for correlation analysis while asyCov is treated as fixed. This process is
#' computationally intensive. A better approach is to request likelihood-based
#' confidence intervals (CIs) by specifying \code{intervals.type="LB"}.  If
#' \code{diag.constraints=FALSE} and \code{cor.analysis=TRUE}, the diagonals
#' are automatically constrained as ones by treating the error variances as
#' computed values rather than as parameters. Since the error variances are not
#' parameters, they are not reported.
#' @param cor.analysis Logical. Analysis of correlation or covariance
#' structure. If \code{cor.analysis=TRUE}, \code{\link[OpenMx]{vechs}} is used
#' to vectorize \code{S}; otherwise, \code{\link[OpenMx]{vech}} is used to
#' vectorize \code{S}.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% Wald CIs based on the z statistic. If
#' it is \code{LB}, it calculates the 95\% likelihood-based CIs on the
#' parameter estimates. Please note that the z values and their associated p
#' values are based on the z statistic. They are not related to the
#' likelihood-based CIs.
#' @param mx.algebras A list of \code{\link[OpenMx]{mxMatrix}} or
#' \code{\link[OpenMx]{mxAlgebra}} objects on the \code{Amatrix},
#' \code{Smatrix}, and \code{Fmatrx}. It can be used to define new functions of
#' parameters and their LBCIs. For example, if the regression coefficients to
#' calculate an indirect effect are stored in A[1,2] and A[1,3], we may define
#' \code{list(ind=mxAlgebra(Amatrix[1,2]*Amatrix[1,3], name="ind"))} See the
#' examples in \code{\link[metaSEM]{Becker92}} and
#' \code{\link[metaSEM]{Hunter83}}. It should be noted that Fmatrix, Amatrix,
#' Smatrix, Iden (a \eqn{p}{p} x \eqn{p}{p} identity matrix), sampleS (sample
#' correlation or covariance matrix), impliedS1, impliedS (model implied
#' correlation or covariance matrix), vecS, invAcov, obj, One, select and
#' constraint and Ematrix (computed error variances when
#' \code{diag.constraints=FALSE}) have been defined internally. You should not
#' create new matrices using these names.
#' @param mxModel.Args A list of arguments passed to
#' \code{\link[OpenMx]{mxModel}}. These include, for example, additional
#' \code{\link[OpenMx]{mxMatrix}} and \code{\link[OpenMx]{mxConstraint}}.
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}. If it is missing, the default is "TSSEM2 (or
#' WLS) Analysis of Correlation Structure" for \code{cor.analysis=TRUE} and
#' "TSSEM2 (or WLS) Analysis of Covariance Structure" for
#' \code{cor.analysis=FALSE}.
#' @param subset.variables An optional character vector of variable names to
#' select variables in the analysis. For example, there are 10 variables in
#' \code{Cov}, say, x1 to x10. We may use \code{c("x1", "x2", "x3")} to select
#' three variables in the analysis. Please note that this argument does not
#' reorder the data. That is, \code{c("x3", "x2", "x1")} is the same as
#' \code{c("x1", "x2", "x3")}.
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed.
#' The argument to be passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Further arguments to be passed to \code{\link[OpenMx]{mxRun}}.
#' @return An object of class \code{wls} with a list of \item{call}{The matched
#' call} \item{Cov}{Input data of either a covariance or correlation matrix}
#' \item{asyCov}{The asymptotic covariance matrix of the input data}
#' \item{noObservedStat}{Number of observed statistics} \item{n}{Sample size}
#' \item{cor.analysis}{logical} \item{noConstraints}{Number of constraints
#' imposed on S} \item{indepModelChisq}{Chi-square statistic of the independent
#' model returned by \code{.indepwlsChisq} } \item{indepModelDf}{Degrees of
#' freedom of the independent model returned by \code{.indepwlsChisq}}
#' \item{mx.fit}{A fitted object returned from \code{\link[OpenMx]{mxRun}}}
#' @note If the input is a list of \code{tssem1.obj}, it returns a list of
#' results for each cluster.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{Becker92}},
#' \code{\link[metaSEM]{Digman97}}, \code{\link[metaSEM]{Hunter83}},
#' \code{\link[metaSEM]{issp89}}, \code{\link[metaSEM]{issp05}}
#' @references Bentler, P.M., & Savalei, V. (2010). Analysis of correlation
#' structures: current status and open problems. In Kolenikov, S., Thombs, L.,
#' & Steinley, D. (Eds.). \emph{Recent Methodological Developments in Social
#' Science Statistics} (pp. 1-36). Hoboken, NJ: Wiley.
#' 
#' Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as multiple-group
#' structural equation models. \emph{Structural Equation Modeling}, \bold{17},
#' 481-509.
#' 
#' Cheung, M. W.-L. (2014). Fixed- and random-effects meta-analytic structural
#' equation modeling: Examples and analyses in R. \emph{Behavior Research
#' Methods}, \bold{46}, 29-40.
#' 
#' Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
#' modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10},
#' 40-64.
#' 
#' Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing
#' covariance matrices in meta-analytic structural equation modeling.
#' \emph{Structural Equation Modeling}, \bold{16}, 28-53.
#' 
#' Joreskog, K. G., Sorbom, D., Du Toit, S., & Du Toit, M. (1999). \emph{LISREL
#' 8: New Statistical Features.} Chicago: Scientific Software International.
#' 
#' McArdle, J. J., & MacDonald, R. P. (1984). Some algebraic properties of the
#' Reticular Action Model for moment structures. \emph{British Journal of
#' Mathematical and Statistical Psychology}, \bold{37}, 234-251.
#' @keywords tssem
#' @examples
#' 
#' \donttest{
#' #### Analysis of correlation structure
#' R1.labels <- c("a1", "a2", "a3", "a4")
#' 
#' R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
#'                0.22, 1.00, 0.30, 0.22,
#'                0.24, 0.30, 1.00, 0.24,
#'                0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4,
#'                dimnames=list(R1.labels, R1.labels))
#' n <- 1000
#' acovR1 <- asyCov(R1, n)
#' 
#' #### One-factor CFA model using lavaan specification
#' model1 <- "f =~ a1 + a2 + a3 + a4"
#' 
#' RAM1 <- lavaan2RAM(model1, obs.variables=R1.labels)
#' 
#' wls.fit1a <- wls(Cov=R1, aCov=acovR1, n=n, RAM=RAM1,
#'                  cor.analysis=TRUE, intervals="LB")
#' summary(wls.fit1a)
#' 
#' ## One-factor CFA model using RAM specification
#' (A1 <- cbind(matrix(0, nrow=5, ncol=4),
#'              matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),
#'              ncol=1)))
#' 
#' (S1 <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))
#' 
#' ## The first 4 variables are observed while the last one is latent.
#' (F1 <- create.Fmatrix(c(1,1,1,1,0), name="F1"))
#' 
#' wls.fit1b <- wls(Cov=R1, aCov=acovR1, n=n, Fmatrix=F1, Smatrix=S1, Amatrix=A1,
#'                  cor.analysis=TRUE, intervals="LB")
#' summary(wls.fit1b)
#' 
#' ## Select 3 variables to analyze
#' model2 <- "f =~ a1 + a2 + a3"
#' 
#' RAM2 <- lavaan2RAM(model2, obs.variables=R1.labels[-4])
#' 
#' wls.fit1c <- wls(Cov=R1, aCov=acovR1, n=n, RAM=RAM2,
#'                  cor.analysis=TRUE, subset.variables=c("a1", "a2", "a3"))
#' summary(wls.fit1c)
#' 
#' #### Multiple regression analysis using lavaan specification
#' R2.labels <- c("y", "x1", "x2")
#' 
#' R2 <- matrix(c(1.00, 0.22, 0.24, 
#'                0.22, 1.00, 0.30, 
#'                0.24, 0.30, 1.00), ncol=3, nrow=3,
#'                dimnames=list(R2.labels, R2.labels))
#' acovR2 <- asyCov(R2, n)
#' 
#' model3 <- "y ~ x1 + x2
#'            ## Variances of x1 and x2 are 1
#'            x1 ~~ 1*x1
#'            x2 ~~ 1*x2
#'            ## x1 and x2 are correlated
#'            x1 ~~ x2"
#' 
#' RAM3 <- lavaan2RAM(model3, obs.variables=R2.labels)
#' 
#' wls.fit2a <- wls(Cov=R2, aCov=acovR2, n=n, RAM=RAM3,
#'                  cor.analysis=TRUE, intervals="z")
#' summary(wls.fit2a)
#' 
#' 
#' #### Multiple regression analysis using RAM specification
#' 
#' ## A2: Regression coefficents
#' #    y x1 x2
#' # y  F T  T 
#' # x1 F F  F 
#' # x2 F F  F 
#' (A2 <- mxMatrix("Full", ncol=3, nrow=3, byrow=TRUE,
#'                free=c(FALSE, rep(TRUE, 2), rep(FALSE, 6)), name="A2"))
#' 
#' ## S2: Covariance matrix of free parameters
#' #    y x1 x2
#' # y  T F  F 
#' # x1 F F  F 
#' # x2 F T  F
#' (S2 <- mxMatrix("Symm", ncol=3, nrow=3, values=c(0.2,0,0,1,0.2,1),
#'                 labels=c("Var_y", NA, NA, NA, "Cov_x1_x2", NA),
#'                 free=c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE), name="S2"))
#' 
#' ## F may be ignored as there is no latent variable.
#' wls.fit2b <- wls(Cov=R2, aCov=acovR2, n=n, Amatrix=A2, Smatrix=S2,
#'                  cor.analysis=TRUE, intervals="LB")
#' summary(wls.fit2b)
#' 
#' 
#' #### Analysis of covariance structure using lavaan specification
#' R3.labels=c("a1", "a2", "a3", "a4")
#' 
#' R3 <- matrix(c(1.50, 0.22, 0.24, 0.18,
#'                0.22, 1.60, 0.30, 0.22,
#'                0.24, 0.30, 1.80, 0.24,
#'                0.18, 0.22, 0.24, 1.30), ncol=4, nrow=4,
#'                dimnames=list(R3.labels, R3.labels))
#' n <- 1000
#' acovS3 <- asyCov(R3, n, cor.analysis=FALSE)
#' 
#' model3 <- "f =~ a1 + a2 + a3 + a4"
#' 
#' RAM3 <- lavaan2RAM(model3, obs.variables=R3.labels)
#' 
#' wls.fit3a <- wls(Cov=R3, aCov=acovS3, n=n, RAM=RAM3,
#'                  cor.analysis=FALSE)
#' summary(wls.fit3a)
#' 
#' #### Analysis of covariance structure using RAM specification
#' (A3 <- cbind(matrix(0, nrow=5, ncol=4),
#'              matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),ncol=1)))
#' 
#' (S3 <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))
#' 
#' F3 <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
#' (F3 <- create.Fmatrix(F3, name="F3", as.mxMatrix=FALSE))
#' 
#' wls.fit3b <- wls(Cov=R3, aCov=acovS3, n=n, Amatrix=A3, Smatrix=S3,
#'                 Fmatrix=F3, cor.analysis=FALSE)
#' summary(wls.fit3b)
#' }
#' 
wls <- function(Cov, aCov, n, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, 
                diag.constraints=FALSE, cor.analysis=TRUE, intervals.type=c("z", "LB"), 
                mx.algebras=NULL, mxModel.Args=NULL, subset.variables=NULL,
                model.name=NULL, suppressWarnings=TRUE, 
                silent=TRUE, run=TRUE, ...) {

    ## Filter out variables not used in the analysis
    if (!is.null(subset.variables)) {
        obslabels <- rownames(Cov)
        index  <- !(subset.variables %in% obslabels)
        if (any(index)) {
            stop(paste0(paste0(subset.variables[index], collapse = ", "), " are not in the \"Cov\".\n"))
        }
        index <- (obslabels %in% subset.variables)
        ## A square matrix of dimensions subset.variables
        index.mat <- matrix(TRUE, nrow=length(index), ncol=length(index))
        index.mat[!index, ] <- index.mat[, !index] <- FALSE

        if (cor.analysis) {
            subset.acov <- vechs(index.mat)
        } else {
            subset.acov <- vech(index.mat)
        }

        ## Filter the variables
        Cov <- Cov[subset.variables, subset.variables]
        aCov <- aCov[subset.acov, subset.acov]
    }
    
    ## Read RAM first. If it is not specified, read individual matrices
    if (!is.null(RAM)) {
        Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
        Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
        Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
    } else {
        if (is.null(Smatrix)) {
            stop("\"Smatrix\" matrix is not specified.\n")
        } else if (is.matrix(Smatrix)) {
            Smatrix <- as.mxMatrix(Smatrix, name="Smatrix")
        } else {
            ## Change the name of the input mxMatrix
            Smatrix@name <- "Smatrix"
        }
        
        ## No. of observed and latent variables
        p <- nrow(Smatrix@values)  

        if (is.null(Amatrix)) {
            ## If Amatrix is not specified, use a zero matrix.
            Amatrix <- as.mxMatrix(matrix(0, nrow=p, ncol=p), name="Amatrix")
        } else if (is.matrix(Amatrix)) {
            Amatrix <- as.mxMatrix(Amatrix, name="Amatrix")
        } else {
            ## Change the name of the input mxMatrix
            Amatrix@name <- "Amatrix"
        }

        if (is.null(Fmatrix)) {
            ## If Fmatrix is not specified, use an identity matrix.
            Fmatrix <- as.mxMatrix(Diag(rep(p,1)), name="Fmatrix")
        } else if (is.matrix(Fmatrix)) {
            Fmatrix <- as.mxMatrix(Fmatrix, name="Fmatrix")
        } else {
            ## Change the name of the input mxMatrix
            Fmatrix@name <- "Fmatrix"
        }    
    }      
    
  ## CheckRAM
  checkRAM(Amatrix=Amatrix, Smatrix=Smatrix, cor.analysis=cor.analysis)

  ## No. of observed and latent variables
  p <- nrow(Smatrix@values)      
  ## A pxp identity matrix
  Id <- as.mxMatrix(Diag(rep(p, 1)), name="Id")

  ## No. of observed variables
  no.var <- ncol(Cov)
  if (is.pd(Cov)) {
    sampleS <- as.mxMatrix(Cov, name="sampleS")
  } else {
    stop("\"Cov\" is not positive definite.\n")
  }

  intervals.type <- match.arg(intervals.type)
  # Default is z
  switch(intervals.type,
         z = intervals <- FALSE,
         LB = intervals <- TRUE)

  # Inverse of asymptotic covariance matrix
  if (is.pd(aCov)) {
    invaCov <- tryCatch(chol2inv(chol(aCov)), error = function(e) e)
    ## It appears that solve() does not fail
    if (inherits(invaCov, "error")) {
      cat("Error in inverting \"aCov\":\n")
      stop(print(invaCov))
    }
    invAcov <- as.mxMatrix(invaCov, name="invAcov")
  } else {
    stop("\"aCov\" is not positive definite.\n")
  }

  if (cor.analysis) {
    if (is.null(model.name)) model.name <- "WLS Correlation"
    ps <- no.var * (no.var - 1)/2
  } else {
    if (is.null(model.name)) model.name <- "WLS Covariance"
    ps <- no.var * (no.var + 1)/2
  }
  
  if (ncol(aCov) != ps)
    stop("No. of dimension of \"Cov\" does not match the multiplier of the dimension of \"aCov\"\n")
  
  ## Assuming no constraint
  Constraints <- 0
  
  if (cor.analysis) {
    
    ## Count no. of dependent variables including both observed and latent variables
    ## Since it is correlation structure, Smatrix@values=1 and Smatrix@free=FALSE on the diagonals.
    Constraints <- Diag(Smatrix@free)
    ## Use nonlinear constraints to impose diagonals as 1
    if (diag.constraints & (sum(Constraints)>0)) {
      
      One <- mxMatrix("Full", values=1, ncol=1, nrow=sum(Constraints), free=FALSE, name="One")
      select <- create.Fmatrix(Constraints, name="select")
      constraint <- mxConstraint( select%*%diag2vec(solve(Id-Amatrix)%&%Smatrix)==One, 
                                  name="constraint" )
      impliedS <- mxAlgebra( (Fmatrix%*%solve(Id-Amatrix))%&%Smatrix, name="impliedS" )
      
      vecS <- mxAlgebra(vechs(sampleS - impliedS), name="vecS")
      obj <- mxAlgebra( t(vecS) %&% invAcov, name = "obj" )
      objective <- mxFitFunctionAlgebra(algebra="obj")
      
      mx.model <- mxModel(model=model.name, Fmatrix, Amatrix, Smatrix, Id, impliedS,
                          vecS, invAcov, obj, objective, sampleS, One, select, 
                          constraint, mxCI(c("Amatrix", "Smatrix")))
    } else {
    ## Consider error variances as functions of parameters  
      ## check types of variables
      ## dv: variables pointed by some variables
      dv <- apply(Amatrix$free, 1, any)
#       ## iv: variables pointed towards other variables
#       iv <- apply(Amatrix$free, 2, any)
#       ## med: mediators
#       med <- iv & dv
      
      ## S1: Smatrix without error variances on the dependent variables
      ## Fix the diagonal elements associated with dv to 0
      S1 <- Smatrix
      S1$name <- "S1"
      diag(S1$labels)[dv] <- NA
      diag(S1$values)[dv] <- 0
      diag(S1$free)[dv] <- FALSE

      ## Use it rather than Smatrix in mxCI()
      S1_labels <- S1$labels
      diag(S1_labels) <- NA
      S1_labels <- c(na.omit(unique(c(S1_labels))))

      ## A function to extract levels of path directions started from IVs to DVs
      path <- function(Amatrix) {
        ## assuming non-recursive models
        if (is.matrix(Amatrix)) Amatrix <- as.mxMatrix(Amatrix)
        A <- Amatrix$free
        
        if (any(Diag(A))) stop("Diagonals on 'A' must be zeros\n")
        ## no. of variables
        p <- ncol(A)
        dv <- apply(A, 1, any)
        ## iv: variables pointed towards other variables
        iv <- apply(A, 2, any)
        
        Level <- list()
        ## IVs
        Level[[1]] <- iv==TRUE&dv==FALSE
        ## no. of variables counted
        count <- length((1:p)[Level[[1]]])
        
        ## do until all p variables are countered 
        while(p > count) {
          level <- length(Level)
          ## predicated by previous level
          cand1 <- A[, Level[[level]], drop=FALSE]
          cand1 <- apply(cand1, 1, any)
          ## predicted by cand1
          cand2 <- A[, cand1, drop=FALSE]
          cand2 <- apply(cand2, 1, any)
          Level[[level+1]] <- cand1==TRUE&cand2==FALSE
          ## no. of elements counted
          count <- count + length((1:p)[Level[[level+1]]])
        }
        Level
      }
      
      ## Levels of directions
      Level <- path(Amatrix)

      ## no error variance involved for ALL variables
      impliedS1 <- mxAlgebra( solve(Id-Amatrix)%&%S1, name="impliedS1")
      
      ## setup for error variances for mediators and dvs
      ## Excluded level1 as they are ivs (started with i=2)
      for (i in 2:length(Level)) {
        ## filter the diagonals of mediators
        text1 <- paste("sel",i, " <- as.mxMatrix(diag(Level[[",i,"]]), name='sel",i,"')", sep="" )
        eval(parse(text=text1))
        ## extract the error variances of mediators based on the previous model implied covariance matrix (i-1)
        text2 <- paste("E",i, " <- mxAlgebra(vec2diag(diag2vec(Id-impliedS",i-1,"))*sel",i,", name='E",i,"')", sep="" )
        eval(parse(text=text2))
        ## Implied covariance matrix with estimated error variances on mediators for (i)
        text3 <- paste("E", 2:i, sep="", collapse="+")
        text4 <- paste("impliedS",i, " <- mxAlgebra(solve(Id-Amatrix)%&%(S1+",text3,"), name='impliedS",i,"')", sep="")
        eval(parse(text=text4))
      }
      
      ## Final Smatrix including all error variances
      text5 <- paste("E", 2:length(Level), sep="", collapse="+")
      ## Smatrix=S1+E2+E3...
      text6 <- paste("Smatrix <- mxAlgebra(S1+", text5, ", name='Smatrix')", sep="")
      eval(parse(text=text6))
      
      impliedS <- mxAlgebra((Fmatrix%*%solve(Id-Amatrix))%&%Smatrix, name="impliedS")
      
      vecS <- mxAlgebra(vechs(sampleS - impliedS), name="vecS")
      obj <- mxAlgebra( t(vecS) %&% invAcov, name = "obj" )
      objective <- mxFitFunctionAlgebra(algebra="obj")
      
      mx.model <- mxModel(model=model.name, Fmatrix, Amatrix, Smatrix, Id, impliedS1, 
                          impliedS, vecS, invAcov, obj, objective, sampleS, 
                          S1, mxCI(c("Amatrix", S1_labels)))
      
      ## Add the matrices into mx.model
      text6a <- paste(",sel",2:length(Level), sep="", collapse="")
      text6b <- paste(",E",2:length(Level), sep="", collapse="")
      text6c <- paste(",impliedS",2:length(Level), sep="", collapse="")
      text7 <- paste("mx.model <- mxModel(mx.model",text6a, text6b, text6c, ")", sep="")
      eval(parse(text=text7))
    }  ## if (diag.constraints & (sum(Constraints)>0))
  } else {
    ## analysis of covariance rather than correlation matrix
    impliedS <- mxAlgebra( (Fmatrix%*%solve(Id-Amatrix))%&%Smatrix, name="impliedS" )
    vecS <- mxAlgebra(vech(sampleS - impliedS), name="vecS")
    
    obj <- mxAlgebra( t(vecS) %&% invAcov, name = "obj" )
    objective <- mxFitFunctionAlgebra(algebra="obj")
    
    mx.model <- mxModel(model=model.name, Fmatrix, Amatrix, Smatrix, Id, impliedS,
                        vecS, invAcov, obj, objective, sampleS, mxCI(c("Amatrix", "Smatrix")))  
  }

  ## Add additiona mxalgebras from RAM
  ## Initialize as NULL  
  algebra.names <- NULL  
  if (!is.null(RAM$mxalgebras)) {
    for (i in seq_along(RAM$mxalgebras)) {
      mx.model <- mxModel(mx.model, RAM$mxalgebras[[i]])
    }
    ## check if they are mxalgebra, not mxconstraint
    algebra.names <- names(RAM$mxalgebras)
    isalgebra <- !grepl("^constraint[0-9]+", algebra.names)
    if (any(isalgebra)) {
        ## Remove names for mxconstraints
        algebra.names <- algebra.names[isalgebra]
        mx.model <- mxModel(mx.model, mxCI(algebra.names))
    }
  }  
    
  ## Add additional mxAlgebras
  if (!is.null(mx.algebras)) {
    for (i in seq_along(mx.algebras)) {
      mx.model <- mxModel(mx.model, mx.algebras[[i]])
    }
    mx.model <- mxModel(mx.model, mxCI(names(mx.algebras)))
  }
    
  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
      for (i in seq_along(mxModel.Args)) {
          mx.model <- mxModel(mx.model, mxModel.Args[[i]])
      }
  }
    
  ## Return mx model without running the analysis
  if (run==FALSE) return(mx.model)

  mx.fit <- tryCatch( mxRun(mx.model, intervals=intervals, suppressWarnings=suppressWarnings, silent=silent, ...), 
                      error = function(e) e)

  # try to run it with error message as output
  if (inherits(mx.fit, "error")) {
      cat("Error in running the mxModel:\n")
      warning(print(mx.fit))
      out <- mx.fit
  } else {
      out <- list(call=match.call(), Cov=Cov, aCov=aCov, noObservedStat=ps, n=n, cor.analysis=cor.analysis, 
                  diag.constraints=diag.constraints, Constraints=Constraints,
                  indepModelChisq=.indepwlsChisq(S=Cov, aCov=aCov, cor.analysis=cor.analysis),
                  indepModelDf=no.var*(no.var-1)/2, mx.model=mx.model, mx.fit=mx.fit,
                  mx.algebras=c(names(mx.algebras), algebra.names), 
                  intervals.type=intervals.type)
      class(out) <- 'wls'
  }
  out
}


#' @rdname wls
tssem2 <- function(tssem1.obj, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, diag.constraints=FALSE,
                   intervals.type = c("z", "LB"), mx.algebras=NULL, mxModel.Args=NULL, subset.variables=NULL,
                   model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  if ( !is.element( class(tssem1.obj)[1], c("tssem1FEM.cluster", "tssem1FEM", "tssem1REM")) )
      stop("\"tssem1.obj\" must be of neither class \"tssem1FEM.cluster\", class \"tssem1FEM\" or \"tssem1REM\".")

  switch(class(tssem1.obj)[1],
         tssem1FEM.cluster = { out <- lapply(tssem1.obj, tssem2, RAM=RAM, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
                                             diag.constraints=diag.constraints, intervals.type=intervals.type,
                                             mx.algebras=mx.algebras, mxModel.Args=mxModel.Args,
                                             subset.variables=subset.variables,
                                             model.name=model.name, suppressWarnings=suppressWarnings, silent=silent,
                                             run=run, ...)
                              class(out) <- "wls.cluster" },
         tssem1FEM = { if (tssem1.obj$cor.analysis==TRUE) {
                        if (is.null(model.name)) model.name <- "TSSEM2 Correlation"
                      } else {
                        if (is.null(model.name)) model.name <- "TSSEM2 Covariance"
                      # to handle symbolic F(T) vs. logical FALSE(TRUE)
                      ## cor.analysis <- as.logical(as.character(cor.analysis))
                      }
                      ## Use the original varible names for the observed covariance matrix
                      pooledS <- coef.tssem1FEM(tssem1.obj)                      
                      ## dimnames(pooledS) <- list(tssem1.obj$original.names, tssem1.obj$original.names)
                      out <- wls(Cov=coef.tssem1FEM(tssem1.obj), aCov=vcov.tssem1FEM(tssem1.obj),
                                 n=sum(tssem1.obj$n), RAM=RAM, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
                                 diag.constraints=diag.constraints, cor.analysis=tssem1.obj$cor.analysis,
                                 intervals.type=intervals.type, mx.algebras=mx.algebras, mxModel.Args=mxModel.Args,
                                 subset.variables=subset.variables,
                                 model.name=model.name, suppressWarnings=suppressWarnings,
                                 silent=silent, run=run, ...) },
         tssem1REM = { cor.analysis <- tssem1.obj$cor.analysis
                      ## Extract the pooled correlation matrix
                      pooledS <- vec2symMat( coef(tssem1.obj, select="fixed"), diag=!cor.analysis)
                      dimnames(pooledS) <- list(tssem1.obj$original.names, tssem1.obj$original.names) 
                      ## Extract the asymptotic covariance matrix of the pooled correlations                      
                      aCov <- vcov(tssem1.obj, select="fixed")
                      aCov.names <- .genCorNames(x=tssem1.obj$original.names, cor.analysis=cor.analysis)
                      dimnames(aCov) <- list(aCov.names$ylabels, aCov.names$ylabels)
      
                      if (cor.analysis==TRUE) {
                        if (is.null(model.name)) model.name <- "TSSEM2 Correlation"
                      } else {
                        if (is.null(model.name)) model.name <- "TSSEM2 Covariance"
                      }
                      ## Use the original varible names for the observed covariance matrix
                      dimnames(pooledS) <- list(tssem1.obj$original.names, tssem1.obj$original.names)
                      out <- wls(Cov=pooledS, aCov=aCov, n=tssem1.obj$total.n, RAM=RAM,
                                 Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix, diag.constraints=diag.constraints,
                                 cor.analysis=cor.analysis, intervals.type=intervals.type, mx.algebras=mx.algebras,
                                 subset.variables=subset.variables,
                                 mxModel.Args=mxModel.Args, model.name=model.name, suppressWarnings = suppressWarnings,
                                 silent=silent, run=run, ...) })
  out
}




#' Estimate the heterogeneity (SD) of the parameter estimates of the TSSEM
#' object
#' 
#' It estimates the heterogeneity of the parameter estimates of the TSSEM
#' objects using either the bootstrap or the delta methods.
#' 
#' The bootstrap method is based on the discussion in Cheung (2018) and Yu et
#' al. (2016). The delta method is an alternative method to obtain the
#' heterogeneity.
#' 
#' @param tssem1.obj An object of class \code{tssem1REM} returned from
#' \code{tssem1()}
#' @param tssem2.obj An object of class \code{wls} returned from
#' \code{tssem2()} or \code{wls()}
#' @param method If it is \code{bootstrap}, random correlation matrices are
#' sampled from the \code{tssem1.obj} by the parametric bootstrap. If it is
#' \code{delta}, the delta method is used to estimate the heterogeneity of the
#' parameter estimates.
#' @param interval The desired interval, e.g., .8 or .95.
#' @param Rep The number of parametric bootstrap. It is ignored when the method
#' is \code{delta}.
#' @param output Either a \code{data.frame} or \code{matrices} of the output.
#' @param nonPD.pop If it is \code{replace}, generated non-positive definite
#' matrices are replaced by generated new ones which are positive definite. If
#' it is \code{nearPD}, they are replaced by nearly positive definite matrices
#' by calling \code{Matrix::nearPD()}. If it is \code{accept}, they are
#' accepted.
#' @return Either a \code{data.frame} or \code{matrices} of the output.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{bootuniR1}}, \code{\link[metaSEM]{bootuniR2}},
#' \code{\link[metaSEM]{Nohe15}}
#' @references Cheung, M. W.-L. (2018). Issues in solving the problem of effect
#' size heterogeneity in meta-analytic structural equation modeling: A
#' commentary and simulation study on Yu, Downes, Carter, and O'Boyle (2016).
#' \emph{Journal of Applied Psychology}, \bold{103}, 787-803.
#' 
#' Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle, E. H. (2016). The
#' problem of effect size heterogeneity in meta-analytic structural equation
#' modeling.  \emph{Journal of Applied Psychology}, \emph{101}, 1457-1473.
#' @keywords tssem
tssemParaVar <- function(tssem1.obj, tssem2.obj, method=c("bootstrap", "delta"),
                         interval=0.8, Rep=50, output=c("data.frame", "matrices"),
                         nonPD.pop=c("replace", "nearPD", "accept")) {

    if (interval <=0 | interval >=1) stop("'interval' must be within 0 and 1.\n")
    
    ## Means of stage 1
    Means <- eval(parse(text="mxEval(Inter, tssem1.obj$mx.fit)"))
    ## Variance componet of the correlation coefficients
    V <- eval(parse(text="mxEval(Tau, tssem1.obj$mx.fit)"))

    mx.model <- tssem2.obj$mx.fit
    ## Maximize the speed
    mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
    mx.model <- mxOption(mx.model, "Standard Errors", "No")
  
    method <- match.arg(method)
    output <- match.arg(output)

    ## delta method
    if (method=="delta") {

        if (!requireNamespace("numDeriv", quietly = TRUE))
            stop("\"numDeriv\" package is required for this function.")

        ## Gradiant function
        ## x: Correlation coefficients
        fun1 <- function(x) {
            sampleS <- as.mxMatrix(vec2symMat(x, diag=FALSE), name="sampleS")
            model <- mxModel(mx.model, sampleS=sampleS)
            fit <- mxRun(model, silent=TRUE)
            coef(fit)
        }
        
        ## Jacobian matrix
        grad <- numDeriv::jacobian(fun1, x=Means)
        Tau2 <- grad %*% V %*% t(grad)
    } else {
        ## bootstrap method
        boot.cor <- rCorPop(Sigma=vec2symMat(Means, diag=FALSE),
                            V=V, k=Rep, nonPD.pop=nonPD.pop)

        boot.coef <- t(sapply(boot.cor, function(x) {
                               sampleS <- as.mxMatrix(x, name="sampleS")
                               model <- mxModel(mx.model, sampleS=sampleS)
                               fit <- mxRun(model, silent=TRUE)
                               coef(fit)
        }))
        Tau2 <- cov(boot.coef, use="complete.obs") 
    }

    ## Parameters for checking
    para <- coef(mxRun(mx.model, silent=TRUE))
    if (is.null(dimnames(Tau2))) {
        dimnames(Tau2) <- list(names(para), names(para))
    }
    SD <- sqrt(diag(Tau2))

    if (output=="data.frame") {
        z <- qnorm((1-interval)/2, lower.tail=FALSE)
        out <- data.frame(Parameter=para, SD=SD, lbound=para-z*SD, ubound=para+z*SD)
        lbound <- paste0(interval*100, "% lbound")
        ubound <- paste0(interval*100, "% ubound")
        colnames(out) <- c("Parameter", "SD", "lbound", "ubound")
    } else {
        out <- list(para=para, SD=SD, Tau2=Tau2)
    }
    out
}



## deltaCV <- function(tssem1.obj, tssem2.obj) {
    
##     if (!requireNamespace("numDeriv", quietly = TRUE)) 
##         stop("\"numDeriv\" package is required for this function.")

##     ## Means
##     Means <- mxEval(Inter, tssem1.obj$mx.fit)
##     ## Variance componet of the correlation coefficients
##     V <- mxEval(Tau, tssem1.obj$mx.fit)

##     mx.model <- tssem2.obj$mx.fit
##     ## Maximize the speed
##     mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
##     mx.model <- mxOption(mx.model, "Standard Errors", "No")
    
##     ## Gradiant function
##     ## x: Correlation coefficients
##     fun <- function(x) {
##         sampleS <- as.mxMatrix(vec2symMat(x, diag=FALSE), name="sampleS")
##         model <- mxModel(mx.model, sampleS=sampleS)
##         fit <- mxRun(model, silent=TRUE)
##         coef(fit)
##     }

##     grad <- numDeriv::jacobian(fun, x=Means)
##     Tau2 <- grad %*% V %*% t(grad)

##     para <- coef(mxRun(mx.model, silent=TRUE))
##     dimnames(Tau2) <- list(names(para), names(para))
    
##     SD <- sqrt(diag(Tau2))
    
##     list(para=para, SD=SD, Tau2=Tau2)
## }



## delta2 <- function(tssem1.obj, Amatrix = NULL, Smatrix = NULL, Fmatrix = NULL, 
##                    diag.constraints = FALSE) {
    
##     if (!requireNamespace("numDeriv", quietly = TRUE)) 
##         stop("\"numDeriv\" package is required for this function.")

##     ## Means
##     Means <- mxEval(Inter, tssem1.obj$mx.fit)
##     ## Variance componet of the correlation coefficients
##     V <- mxEval(Tau, tssem1.obj$mx.fit)

##     ## Setup the model without running it
##     mx.model <- tssem2(tssem1.obj=tssem1.obj, Amatrix=Amatrix, Smatrix=Smatrix,
##                        Fmatrix=Fmatrix, diag.constraints=diag.constraints,
##                        run=FALSE)

##     ## Maximize the speed
##     mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
##     mx.model <- mxOption(mx.model, "Standard Errors"  , "No")
    
##     ## Gradiant function
##     ## x: Correlation coefficients
##     fun <- function(x) {
##         sampleS <- as.mxMatrix(vec2symMat(x, diag=FALSE), name="sampleS")
##         model <- mxModel(mx.model, sampleS=sampleS)
##         fit <- mxRun(model, silent=TRUE)
##         coef(fit)
##     }

##     grad <- numDeriv::jacobian(fun, x=Means)
##     Tau2 <- grad %*% V %*% t(grad)

##     para <- coef(mxRun(mx.model, silent=TRUE))
##     dimnames(Tau2) <- list(names(para), names(para))

##     list(para=para, Tau2=Tau2)
## }
