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
