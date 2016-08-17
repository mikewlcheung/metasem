impliedSigma <- function(Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, labels=NULL,
                         cor.analysis=TRUE, ...) {
  if (is.null(Smatrix)) {
    stop("\"Smatrix\" matrix is not specified.\n")
  } else {
    if (is.matrix(Smatrix)) Smatrix <- as.mxMatrix(Smatrix)
    ## No. of observed and latent variables
    p <- nrow(Smatrix@values)
    Smatrix@name <- "Smatrix"
  }

  if (is.null(Amatrix)) {
    stop("\"Amatrix\" matrix is not specified.\n")
  } else {
    if (is.matrix(Amatrix)) Amatrix <- as.mxMatrix(Amatrix)
    Amatrix@name <- "Amatrix"
  }

  if (is.null(Fmatrix)) {
    Fmatrix <- as.mxMatrix(Diag(p), name="Fmatrix")
  } else {
    if (is.matrix(Fmatrix)) Fmatrix <- as.mxMatrix(Fmatrix)
    Fmatrix@name <- "Fmatrix"
  }

  ## A pxp identity matrix
  Id <- as.mxMatrix(Diag(p), name="Id")

  ## Model implied correlation/covariance matrix including latent variables
  SigmaAll <- mxAlgebra( solve(Id-Amatrix)%&%Smatrix, name="SigmaAll" )
  
  ## Model implied correlation/covariance matrix of the observed variables
  SigmaObs <- mxAlgebra( Fmatrix%&%SigmaAll, name="SigmaObs" )
  
  if (cor.analysis) {
     ## Create One vector for the diagonal constraint   
     One <- create.mxMatrix(rep(1,p), type="Full", ncol=1, nrow=p, name="One")

     ## Ensure observed and latent are standardized
     minFit <- mxAlgebra( sum((One-diag2vec(SigmaAll))^2), name="minFit" )

     model <- mxModel(model="impliedSigma", Amatrix, Smatrix, Fmatrix, Id, One,
                      SigmaAll, SigmaObs, minFit,
                      mxFitFunctionAlgebra("minFit"))
  } else {
     ## Covariance matrix, no need for the constraint
     model <- mxModel(model="impliedSigma", Amatrix, Smatrix, Fmatrix, Id, 
                      SigmaAll, SigmaObs)    
  }

  mx.fit <- mxRun(model, silent=TRUE)

  A <- eval(parse(text = "mxEval(Amatrix, mx.fit)"))
  S <- eval(parse(text = "mxEval(Smatrix, mx.fit)"))
  F <- eval(parse(text = "mxEval(Fmatrix, mx.fit)"))
  SigmaObs <- eval(parse(text = "mxEval(SigmaObs, mx.fit)"))
  SigmaAll <- eval(parse(text = "mxEval(SigmaAll, mx.fit)"))

  ## Create the labels for the matrices
  if (!is.null(labels)) {
      if (length(labels==p)) {
          fmatrix <- Fmatrix@values
          index <- apply(fmatrix, 1, function(x) which(x==1))
          labels.obs <- labels[index]
          dimnames(A) <- dimnames(S) <- dimnames(SigmaAll) <- list(labels, labels)
          dimnames(SigmaObs) <- list(labels.obs, labels.obs)
          dimnames(F) <- list(labels.obs, labels)
      } else {
          warning("The length of \"labels\" is different from that in \"Smatrix\".\n")
      }
  }
          
  if (cor.analysis) {
      minFit <- c(eval(parse(text = "mxEval(minFit, mx.fit)")))
      status <- c(mx.fit$output$status[[1]])
  } else {
      minFit <- 0
      status <- 0
  }

  if (!isTRUE(all.equal(minFit, 0))) warning("The diagonals of the correlation matrix are not zero! ",
                                     "You should not trust the results.\n")

  if (status!=0) warning("The status code of optimization is non-zero. ",
                         "Please check if there are too many free parameters in your population model.\n")  

  out <- list(A=A, S=S, F=F, SigmaObs=SigmaObs, SigmaAll=SigmaAll, minFit=minFit,
              status=status, mx.fit=mx.fit)
  class(out) <- "impliedSigma"
  out
}

print.impliedSigma <- function(x, ...) {
  if (!is.element("impliedSigma", class(x)))
    stop("\"x\" must be an object of class \"uniR1\".")
  print(x[-8])
}

