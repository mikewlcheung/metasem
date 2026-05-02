#' Create or Generate the Model Implied Correlation or Covariance Matrices
#' 
#' It creates or generates the model implied correlation or covariance matrices
#' based on the RAM model specification.
#' 
#' This function can be used to generate the model implied correlation matrix
#' for the standardized parameters with the \code{corr=TRUE} argument. Suppose
#' we want to calculate the population correlation matrix for a mediation model
#' with x, m, and y. We only need to specify the population path coefficients
#' among x, m, and y in the \code{Amatrix}. We do not need to specify the
#' population error variances of m and y. We treat the error variances as
#' unknown parameters by giving them starting values in the \code{Smatrix}
#' matrix. When the covariance matrix is requested by specifying
#' \code{corr=FALSE}, it simply calculates the population model covariance
#' matrix by treating the values in \code{Smatrix} as the population values.
#' 
#' @aliases impliedR rimpliedR
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param Amatrix If \code{RAM} is not specified, an \code{Amatrix} is
#' required. An asymmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Smatrix If \code{RAM} is not specified, an \code{Smatrix} is
#' required. A symmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Fmatrix A filter matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is missing, an identity matrix
#' with the same dimensions of \code{Smatrix} will be created, which means that
#' all variables are observed. If it is a matrix, it will be converted into
#' \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function. It
#' is not required when there is no latent variable.
#' @param Mmatrix An optional matrix of the mean vector. It is assumed zeros if
#' missing.
#' @param AmatrixSD Standard deviations (SD) of the elements in the
#' \code{Amatrix}. If it is missing, a matrix of zero is created.
#' @param SmatrixSD Standard deviations (SD) of the elements in the
#' \code{Smatrix}. If it is missing, a matrix of zero is created.
#' @param k Number of studies.
#' @param corr Logical. The output is either the model implied correlation
#' matrix or the covariance matrix.
#' @param labels A character vector of the observed and latent variables with
#' the same dimensions as that in the \code{Amatrix} and \code{Smatrix}.
#' @param nonPD.pop If it is \code{replace}, generated non-positive definite
#' matrices are replaced by generated new ones which are positive definite. If
#' it is \code{nearPD}, they are replaced by nearly positive definite matrices
#' by calling \code{Matrix::nearPD()}. If it is \code{accept}, they are
#' accepted.
#' @param \dots Not used.
#' @return A list of RAM matrices, the model implied correlation or covariance
#' matrix of the observed variables (\code{SigmaObs}), of both observed and
#' latent variables (\code{SigmaAll}), the minimum fit (\code{minFit}) which
#' should be zero, and the status code of the optimization (\code{status})
#' which should also be zero when the optimization is fine. The last object is
#' \code{mx.fit} which is the output after running the model. It can be used in
#' the diagnosis.
#' @note It is important to ensure that all the population values in
#' \code{Amatrix} must be set as fixed parameters; otherwise, these values may
#' be altered with the \code{corr=TRUE} argument. When there is an error or
#' warning message about the status code, there is a high chance that some of
#' the values in \code{Amatrix} are incorrectly set as free parameters.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' set.seed(100)
#' 
#' ## A one-factor CFA model
#' model <- "f =~ 0.3*x1 + 0.4*x2 + 0.5*x3
#'           f ~~ 1*f"
#' 
#' RAM <- lavaan2RAM(model)
#' 
#' impliedR(RAM, corr=TRUE)
#' 
#' ## A simple mediation model
#' ## All are population parameters in the A matrix
#' A1 <- matrix(c(0, 0, 0,
#'                0.3, 0, 0,
#'                0.2, 0.4, 0), nrow=3, ncol=3, byrow=TRUE,
#'              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
#' A1             
#' 
#' ## Variance of x is fixed at 1 while the other variances are free.
#' S1 <- matrix(c(1, 0, 0,
#'                0, "0.1*ErrVarM",0,
#'                0, 0, "0.1*ErrVarY"), nrow=3, ncol=3,
#'              dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
#' S1
#' 
#' impliedR(Amatrix=A1, Smatrix=S1)
#' 
#' ## SD of A1
#' A1SD <- matrix(c(0, 0, 0,
#'                  0.1, 0, 0,
#'                  0.1, 0.1, 0), nrow=3, ncol=3, byrow=TRUE,
#'                dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
#' A1SD
#' 
#' rimpliedR(Amatrix=A1, Smatrix=S1, AmatrixSD=A1SD, k=2)
#' 
#' ## A CFA model
#' A2 <- matrix(c(0, 0, 0, 0.3,
#'                0, 0, 0, 0.4,
#'                0, 0, 0, 0.5,
#'                0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
#'              dimnames=list(c("x1", "x2", "x3", "f"),
#'                            c("x1", "x2", "x3", "f")))
#' A2
#' 
#' ## Variance of f is fixed at 1 while the other variances are free.
#' S2 <- matrix(c("0.7*Err1", 0, 0, 0,
#'                 0, "0.7*Err2", 0, 0,
#'                 0, 0, "0.7*Err3", 0,
#'                 0, 0, 0, 1), nrow=4, ncol=4,
#'             dimnames=list(c("x1", "x2", "x3", "f"), c("x1", "x2", "x3", "f")))
#' S2
#' 
#' F2 <- create.Fmatrix(c(1,1,1,0), as.mxMatrix=FALSE)
#' F2
#' 
#' ## Model implied correlation matrix
#' impliedR(Amatrix=A2, Smatrix=S2, Fmatrix=F2, corr=TRUE)
#' 
#' ## Model implied covariance matrix
#' impliedR(Amatrix=A2, Smatrix=S2, Fmatrix=F2, corr=FALSE)
#' 
#' ## SD of A2
#' A2SD <- matrix(c(0, 0, 0, 0.1,
#'                  0, 0, 0, 0.1,
#'                  0, 0, 0, 0.1,
#'                  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
#'                dimnames=list(c("x1", "x2", "x3", "f"),
#'                              c("x1", "x2", "x3", "f")))               
#' A2SD
#' 
#' ## SD of S2: correlated between x1 and x2
#' S2SD <- matrix(c(0, 0.1, 0, 0,
#'                  0.1, 0, 0, 0,
#'                  0, 0, 0, 0.1,
#'                  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
#'                dimnames=list(c("x1", "x2", "x3", "f"),
#'                              c("x1", "x2", "x3", "f")))               
#' S2SD
#' 
#' rimpliedR(Amatrix=A2, Smatrix=S2, Fmatrix=F2, AmatrixSD=A2SD,
#'           SmatrixSD=S2SD, k=2)
#' 
impliedR <- function(RAM, Amatrix, Smatrix, Fmatrix, Mmatrix, corr=TRUE,
                     labels, ...) {

  if (!missing(RAM)) {
    Amatrix <- RAM$A
    Smatrix <- RAM$S
    Fmatrix <- RAM$F
    Mmatrix <- RAM$M
  }
    
  if (missing(Smatrix)) {
    stop("\"Smatrix\" matrix is not specified.\n")
  } else {
    if (is.matrix(Smatrix)) Smatrix <- as.mxMatrix(Smatrix)
    ## No. of observed and latent variables
    p <- nrow(Smatrix@values)
    Smatrix@name <- "Smatrix"
  }

  if (missing(Amatrix)) {
    stop("\"Amatrix\" matrix is not specified.\n")
  } else {
    if (is.matrix(Amatrix)) Amatrix <- as.mxMatrix(Amatrix)
    Amatrix@name <- "Amatrix"
  }

  if (missing(Fmatrix)) {
    Fmatrix <- as.mxMatrix(Diag(p), name="Fmatrix")
  } else {
    if (is.matrix(Fmatrix)) Fmatrix <- as.mxMatrix(Fmatrix)
    Fmatrix@name <- "Fmatrix"
  }
  
  if (corr | missing(Mmatrix)) {
    Mmatrix <- as.mxMatrix(matrix(0, nrow=1, ncol=p), name="Mmatrix")
  } else {
    if (is.matrix(Mmatrix)) Mmatrix <- as.mxMatrix(Mmatrix)
    Mmatrix@name <- "Mmatrix"
  }
  
  ## A pxp identity matrix
  Id <- as.mxMatrix(Diag(p), name="Id")

  ## Model implied correlation/covariance matrix including latent variables
  SigmaAll <- mxAlgebra( solve(Id-Amatrix) %&% Smatrix, name="SigmaAll" )

  ## Model implied correlation/covariance matrix of the observed variables
  SigmaObs <- mxAlgebra( Fmatrix %&% SigmaAll, name="SigmaObs" )

  ## Model implied mean
  Mu <- mxAlgebra( Mmatrix %*% t(Fmatrix %*% solve(Id-Amatrix)), name="Mu")
  
  if (corr) {
    ## Create One vector for the diagonal constraint
    One <- create.mxMatrix(rep(1,p), type="Full", ncol=1, nrow=p, name="One")

    ## Ensure observed and latent are standardized
    minFit <- mxAlgebra( sum((One-diag2vec(SigmaAll))^2), name="minFit" )

    model <- mxModel(model="impliedR", Amatrix, Smatrix, Fmatrix, Mmatrix,
                     Id, One, SigmaAll, SigmaObs, Mu, minFit,
                     mxFitFunctionAlgebra("minFit"))
  } else {
    ## Covariance matrix, no need for the constraint
    model <- mxModel(model="impliedSigma", Amatrix, Smatrix, Fmatrix, Mmatrix,
                     Id, SigmaAll, SigmaObs, Mu)
  }

  mx.fit <- mxRun(model, silent=TRUE)

  A <- eval(parse(text = "mxEval(Amatrix, mx.fit)"))
  S <- eval(parse(text = "mxEval(Smatrix, mx.fit)"))
  F <- eval(parse(text = "mxEval(Fmatrix, mx.fit)"))
  M <- eval(parse(text = "mxEval(Mmatrix, mx.fit)"))
  SigmaObs <- eval(parse(text = "mxEval(SigmaObs, mx.fit)"))
  SigmaAll <- eval(parse(text = "mxEval(SigmaAll, mx.fit)"))
  Mu <- eval(parse(text = "mxEval(Mu, mx.fit)"))

  ## Create the labels for the matrices
  ## Index for the observed variables
  index <- apply(Fmatrix@values, 1, function(x) which(x==1))
  
  if (missing(labels)) {
    if (!is.null(dimnames(Smatrix@values))) {
      labels <- colnames(Smatrix@values)
    } else if (!is.null(dimnames(Amatrix@values))) {
      labels <- colnames(Amatrix@values)
    } else if (!is.null(dimnames(Fmatrix@values))) {
      labels <- colnames(Fmatrix@values)
    } else {
      labels <- NULL
    }
  } else if (length(labels)!=p) {
    warning("Length of \"labels\" is different from the number of variables.\n")
  }
      
  if (!is.null(labels)) {
    labels.obs <- labels[index]
    dimnames(A) <- dimnames(S) <- dimnames(SigmaAll) <- list(labels, labels)
    dimnames(SigmaObs) <- list(labels.obs, labels.obs)
    dimnames(F) <- list(labels.obs, labels)
    dimnames(M) <- list("1", labels)
    dimnames(Mu) <- list("1", labels.obs)
  }
  
  if (corr) {
    ## minFit is the amount of misfit on the constraints
    ## It should be close to zero.
    minFit <- c(eval(parse(text = "mxEval(minFit, mx.fit)")))
    status <- c(mx.fit$output$status[[1]])
  } else {
    ## It is zero by definition for covariance matrix.
    minFit <- 0
    status <- 0
  }

  if (!isTRUE(all.equal(minFit, 0))) {
    warning("The diagonals of the correlation matrix are not zero! ",
            "You should not trust the results.\n")
  }
  
  if (status!=0) {
    warning("The status code of optimization is non-zero. ",
            "Please check if there are too many free parameters in your population model.\n")
  }

  out <- list(A=A, S=S, F=F, M=M, SigmaObs=SigmaObs, SigmaAll=SigmaAll, Mu=Mu,
              corr=corr, minFit=minFit, status=status, mx.fit=mx.fit)
  class(out) <- "impliedR"
  out
}

print.impliedR <- function(x, ...) {
    if (!is.element("impliedR", class(x)))
        stop("\"x\" must be an object of class \"impliedR\".")
    cat("Amatrix:\n")
    print(x$A)
    cat("\nSmatrix:\n")
    print(x$S)
    cat("\nFmatrix:\n")
    print(x$F)
    cat("\nMmatrix:\n")
    print(x$M)    
    cat("\nModel implied matrix of the observed variables:\n")
    print(x$SigmaObs)
    cat("\nModel implied matrix of the observed and latent variables:\n")
    print(x$SigmaAll)
    cat("\nModel implied vector of the observed means:\n")
    print(x$Mu)
    cat("\nCorrelation matrix:", x$corr)
    cat("\nSigma of the observed variables is positive definite:", is.pd(x$SigmaObs))
    cat("\nSigma of both the observed and latent variables is positive definite:", is.pd(x$SigmaAll))
    if (x$corr) {
        cat("\nMinimum value of the fit function (it should be close to 0 for correlation analysis): ", x$minFit)
        cat("\nStatus code of the optimization (it should be 0 for correlation analysis): ", x$status)
    }
    cat("\n")
}

## Generate model implied matrices from random parameters
## Random effects are independent.
#' @rdname impliedR
rimpliedR <- function(RAM, Amatrix, Smatrix, Fmatrix, AmatrixSD, SmatrixSD, k=1,
                      corr=TRUE, nonPD.pop=c("replace", "nearPD", "accept")) {

  ## Only values are used in matrices
  if (!missing(RAM)) {
    Amatrix <- RAM$A
    Smatrix <- RAM$S
    Fmatrix <- RAM$F
  }
 
  ## No. of observed variables
  p <- ncol(Amatrix)
  ## No. of elements in Amatrix
  # n <- p*p
  
  ## All variables are observed.
  if (missing(Fmatrix)) Fmatrix <- Diag(p)

  ## If missing SD matrices, use a zero matrix  
  if (missing(AmatrixSD)) AmatrixSD <- matrix(0, ncol=p, nrow=p)
  if (missing(SmatrixSD)) SmatrixSD <- matrix(0, ncol=p, nrow=p)
  
  if (!all(sapply(list(dim(Amatrix), dim(Smatrix), dim(SmatrixSD)),
                  FUN=identical, dim(AmatrixSD))))
    stop("Dimensions of \"Amatrix\", \"Smatrix\", \"AmatrixSD\", and \"SmatrixSD\" must be the same.")   
    
  ## Try to get the labels of all variables from A and then S
  labels <- colnames(Amatrix)
  if (is.null(labels)) labels <- colnames(Smatrix)
  
  ## Select the labels of the observed variables
  if (!is.null(labels)) labels <- labels[as.logical(colSums(Fmatrix))]
    
  ## A vector of means of Amatrix by column major
  A.mean <- as.numeric(Amatrix)  
  ## A diagonal matrix of variances of Amatrix by column major
  A.var <- diag(c(AmatrixSD^2))

  ## Model implied R or S
  impR1 <- impliedR(Amatrix=Amatrix, Smatrix=Smatrix, corr=corr)
  ## A vector of means of Smatrix
  if (corr) {
    S.mean <- vechs(impR1$S)
    S.var <- diag(vechs(SmatrixSD^2))
  } else {
    S.mean <- vech(impR1$S)
    S.var <- diag(vech(SmatrixSD^2))
  }
    
  nonPD.pop <- match.arg(nonPD.pop)  
  ## Count for nonPD matrices
  nonPD.count <- 0
  
  ## Generate a model implied R and return if it is PD
  genImpR <- function() {
    ## Generate random A matrix
    A <- matrix(mvtnorm::rmvnorm(n=1, mean=A.mean, sigma=A.var), ncol=p, nrow=p)
    ## Generate random S matrix
    S <- mvtnorm::rmvnorm(n=1, mean=S.mean, sigma=S.var)
    ## Convert S back to a pxp matrix
    S <- vec2symMat(S, diag=!corr)
    ## Replace the diagonals from the model implied R
    ## Elements are either 1 or starting values
    if (corr) diag(S) <- diag(Smatrix)
   
    impR2 <- impliedR(Amatrix=A, Smatrix=S, Fmatrix=Fmatrix, corr=corr)
    
    ## isPD includes: status=0 and PD
    list(R=impR2$SigmaObs, isPD=(impR2$status==0 & is.pd(impR2$SigmaObs)))
  }

  ## Generate random correlation matrices
  genCor <- function() {
    impR3 <- genImpR()
    R <- impR3$R
    isPD <- impR3$isPD
    
    ## R is nonPD
    if (!isPD) {
      ## global rather than local assignment
      nonPD.count <<- nonPD.count+1
      switch(nonPD.pop,
             replace = while (!isPD) {
               impR4 <- genImpR()
               R <- impR4$R
               ## isPD includes: status=0 and PD
               isPD <- impR4$isPD
               nonPD.count <<- nonPD.count+1
             },
             nearPD = {R <- as.matrix(Matrix::nearPD(R, corr=corr, 
                                                     keepDiag=corr)$mat)},
             accept = {} )
    }
    ## Ad hoc, R may not be symmetric due to the precision
    R[lower.tri(R)] <- t(R)[lower.tri(t(R))]
    if (!is.null(labels)) dimnames(R) <- list(labels, labels)
    R
  }  
  
  ## Repeat it k times
  ## Simplify it when AmatrixSD=0 and SmatrixSD=0
  if (all(c(AmatrixSD, SmatrixSD)==0)) {
    tmp <- genCor()
    out <- replicate(n=k, tmp, simplify=FALSE)
  } else {
    out <- replicate(n=k, genCor(), simplify=FALSE) 
  }
  
  attr(out, "k") <- k
  attr(out, "nonPD.count") <- nonPD.count
  out
}
