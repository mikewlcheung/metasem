## Generate sample correlation matrices
rCor <- function(Sigma, V, n, corr=TRUE, raw.data=FALSE,
                 nonPD.pop=c("replace", "nearPD", "accept"),
                 nonPD.sam=c("stop", "nearPD")) {
    
    nonPD.pop <- match.arg(nonPD.pop)
    nonPD.sam <- match.arg(nonPD.sam)
    
    ## Generate population matrices   
    P <- rCorPop(Sigma=Sigma, V=V, k=length(n), corr=corr,
                 nonPD.pop=nonPD.pop)
    
    ## Generate sample matrices
    R <- rCorSam(Sigma=P, n=n, corr=corr, raw.data=raw.data,
                 nonPD.sam=nonPD.sam)
    list(P=P, R=R)
}

## Generate population correlation matrices
rCorPop <- function(Sigma, V, k, corr=TRUE, 
                    nonPD.pop=c("replace", "nearPD", "accept")) {

  ## Convert them to matrices
  if (!is.matrix(Sigma)) Sigma <- as.matrix(Sigma)
  if (!is.matrix(V)) V <- as.matrix(V)    
    
  if (!is.pd(Sigma)) stop("'Sigma' is not positive definite.\n")
  if (!is.pd(V)) stop("'V' is not positive definite.\n")
  
  if (corr) {
    Sigma_v <- vechs(cov2cor(Sigma))
  } else {
    Sigma_v <- vech(Sigma)
  }
  
  if (length(Sigma_v) != ncol(V)) stop("Dimensions of 'Sigma' and 'V' are different.\n")
  
  nonPD.pop <- match.arg(nonPD.pop)
  
  ## Count for nonPD matrices
  nonPD.count <- 0
  
  ## Generate 1 sample
  genCor <- function() {
    dat <- mvtnorm::rmvnorm(n=1, mean=Sigma_v, sigma=V)
    R <- vec2symMat(dat, diag=!corr)
    
    isPD <- is.pd(R)
    
    ## Generated R is nonPD
    if (!isPD) {
      ## global rather than local assignment
      nonPD.count <<- nonPD.count+1
      switch(nonPD.pop,
             replace = while (!isPD) {
               dat <- mvtnorm::rmvnorm(n=1, mean=Sigma_v, sigma=V)
               R <- vec2symMat(dat, diag=!corr)
               isPD <- is.pd(R)
               nonPD.count <<- nonPD.count+1
             },
             nearPD = {R <- as.matrix(Matrix::nearPD(R, corr=corr, keepDiag=corr)$mat)},
             accept = {} )
    }
    ## Ad hoc, R may not be symmetric due to the precision
    R[lower.tri(R)] <- t(R)[lower.tri(t(R))]
    #if (!is.null(dimnames(Sigma))) dimnames(R) <- dimnames(Sigma)
    R
  }
  ## Repeat it k times
  out <- replicate(n=k, genCor(), simplify=FALSE)
  if (!is.null(dimnames(Sigma))) {
    out <- lapply(out, function(x) {dimnames(x) <- dimnames(Sigma); x})
  }
  
  attr(out, "Sigma") <- Sigma
  attr(out, "V") <- V
  attr(out, "k") <- k
  attr(out, "corr") <- corr
  attr(out, "nonPD.count") <- nonPD.count
  attr(out, "nonPD.pop") <- nonPD.pop
  class(out) <- "CorPop"
  out
}

summary.CorPop <- function(object, ...) {
  if (!is.element("CorPop", class(object)))
      stop("\"object\" must be an object of class \"CorPop\".")
  
  corr <- attr(object, "corr")
  
  ## Adhoc fix when object is not a list of matrices.
  ## One variable or two variables with corr=TRUE
  if ( ncol(object[[1]])==1 | (ncol(object[[1]])==2 & corr==TRUE) ) {
      fix <- TRUE
  } else {
      fix <- FALSE
  }

  ## R: average correlation matrix based on the generated data
  if (corr) {
      my.df <- t(sapply(object, vechs))
      if (fix) my.df <- matrix(my.df, ncol=1)
      R <- vec2symMat(colMeans(my.df), diag=FALSE)
  } else {
      my.df <- t(sapply(object, vech))
      if (fix) my.df <- matrix(my.df, ncol=1)
      R <- vec2symMat(colMeans(my.df), diag=TRUE)
  }
  
  Sigma <- attr(object, "Sigma")
  dimnames(R) <- dimnames(Sigma)
  V_Samp <- cov(my.df) 

  ## Add dimnames to V_Samp
  V.names <- outer(colnames(Sigma), colnames(Sigma), paste0)
  if (corr) {
    V.names <- OpenMx::vechs(V.names)
  } else {
    V.names <- OpenMx::vech(V.names)
  }
  dimnames(V_Samp) <- list(V.names, V.names)

  V_Pop <- attr(object, "V") 
  dimnames(V_Pop) <- list(V.names, V.names)
  
  out <- list(Sigma=Sigma, V_Pop=V_Pop, R=R, V_Samp=V_Samp,
              k=attr(object, "k"), corr=corr,
              nonPD.count=attr(object, "nonPD.count"),
              nonPD.pop=attr(object, "nonPD.pop"))
  class(out) <- "summary.CorPop"
  out
}

print.summary.CorPop <- function(x, ...) {
  if (!is.element("summary.CorPop", class(x)))
    stop("\"x\" must be an object of class \"summary.CorPop\".")
  
  cat("Population Sigma:\n")
  print(x$Sigma)
  cat("\nEmpirical Sigma:\n")
  print(x$R)
  cat("\nPopulation V:\n")
  print(x$V_Pop)
  cat("\nEmpirical V:\n")
  print(x$V_Samp)
  cat("\nMethod to handle non-positive definite matrices:", x$nonPD.pop)
  cat("\nNumber of samples requested:", x$k, "\n")
  cat("\nCount of total samples generated:", x$k+x$nonPD.count, "\n")
  cat("\nCount of non-positive definite matrices:", x$nonPD.count, "\n")
}  
                  
## Generate sample correlation matrices
rCorSam <- function(Sigma, n, corr=TRUE, raw.data=FALSE, 
                    nonPD.sam=c("stop", "nearPD")) {
  ## Convert Sigma into a list
  if (!is.list(Sigma)) Sigma <- list(Sigma)
  
  ## Convert Sigma into matrices
  Sigma <- lapply(Sigma, as.matrix)

  nonPD.sam <- match.arg(nonPD.sam)
    
  ## Return more than 1 matrices when either Sigma is a list or n is a vector.
  ## Not checked if the lengths are not matched.  
  if ( length(Sigma)>1 | length(n)>1 ) {
    mapply(rCorSam, Sigma=Sigma, n=n, corr=corr, raw.data=raw.data, 
           nonPD.sam=nonPD.sam, SIMPLIFY=FALSE)
  } else {
    ## Output 1 matrix
    
    ## The input population matrix must be positive definite;
    ## otherwise, sample matrices cannot be generated.
    sigma <- Sigma[[1]]
    if (!is.pd(sigma)) {
      switch(nonPD.sam,
             stop = stop("Sigma is not positive definite.\n"),
             nearPD = {warning("Sigma is not positive definite. 'NearPD' is applied.\n")
                       sigma <- as.matrix(Matrix::nearPD(sigma, corr=corr, 
                                                         keepDiag=corr)$mat)})
    }
    
    ## Generate raw data
    if (raw.data) {
      ## In principle, the following trunk is not necessary.
      ## It may break and return an error when Sigma is PD but it does not

      # dat <- try(mvtnorm::rmvnorm(n=n, mean=rep(0, ncol(sigma)), sigma=sigma), silent=TRUE)
      # if (inherits(dat, "try-error")) {
      #     switch( nonPD.sam,
      #        stop = stop("Sigma is not positive definite.\n"),
      #        nearPD = {warning("Sigma is not positive definite. 'NearPD' is applied.\n")
      #            sigma <- as.matrix(Matrix::nearPD(sigma, corr=corr, keepDiag=corr)$mat)
      #            ## Assuming it works fine after applying nearPD.
      #            dat <- mvtnorm::rmvnorm(n=n, mean=rep(0, ncol(sigma)), sigma=sigma) })
      # }     
      dat <- mvtnorm::rmvnorm(n=n, mean=rep(0, ncol(sigma)), sigma=sigma)
      if (corr) out <- cor(dat) else out <- cov(dat)
    } else {
      ## Generate coviances matrices directly from a Wishart distribution
      ## n=1: one set of generated data
      # out <- try(stats::rWishart(n=1, df=(n-1), Sigma=sigma)/(n-1), silent=TRUE)
      # if (inherits(out, "try-error")) {
      #     switch( nonPD.sam,
      #        stop = stop("Sigma is not positive definite.\n"),
      #        nearPD = {warning("Sigma is not positive definite. 'NearPD' is applied.\n")
      #            sigma <- as.matrix(Matrix::nearPD(sigma, corr=corr, keepDiag=corr)$mat)
      #            out <- stats::rWishart(n=1, df=(n-1), Sigma=sigma)/(n-1) })
      # }  
      out <- stats::rWishart(n=1, df=(n-1), Sigma=sigma)/(n-1)
      out <- out[, , 1]
      if (corr) out <- cov2cor(out)
    }
      ## ## Ad hoc, out is not symmetric due to the precision issue
      ## Don't do it for scalar correlations/variances
      if (is.matrix(out)) out[lower.tri(out)] <- t(out)[lower.tri(t(out))]
      if (!is.null(dimnames(sigma))) dimnames(out) <- dimnames(sigma)
      out
  }
}

