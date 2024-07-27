Cor2DataFrame <- function(x, n, v.na.replace=TRUE, cor.analysis=TRUE,
                          acov=c("weighted", "individual", "unweighted"),
                          Means, row.names.unique=FALSE, append.vars=TRUE,
                          asyCovOld=FALSE, ...) {
  
  acov <- match.arg(acov, c("weighted", "individual", "unweighted"))
  
  ## x is a list of "data", "n", ...
  if (all(c("data", "n") %in% names(x))) {
    my.cov <- x$data
    n <- x$n
    obslabels <- colnames(x$data[[1]])
  } else {
    ## x is just a list of correlation matrices. "n" is provided as an argument.
    my.cov <- x
    obslabels <- colnames(x[[1]])
  }

  if (length(my.cov) != length(n)) stop("Lengths of 'x' and 'n' are different.\n")
  
  if (cor.analysis) {
    ## Standardize and then vechs()
    my.df <- list2matrix(x=suppressWarnings(lapply(my.cov, cov2cor)), diag=FALSE)
  } else {
    ## vech()
    my.df <- list2matrix(x=my.cov, diag=TRUE)
  }

  ## NA is not allowed in definition variables
  ## Check if there are NAs in my.cov and replace them with the weighted or unweighted averages
  if (any(sapply(my.cov, is.na)) & v.na.replace) {
    ## Replace NA with 0 before calculations
    my.x <- lapply(my.cov, function(z) {z[is.na(z)] <- 0; z} )
    if (acov=="unweighted") {
      ## x: original covariance matrices in the input
      ## Unweighted means = sum of r/(no. of studies)
      cov.mean <- Reduce("+", my.x)/pattern.na(x, show.na = FALSE)
    } else if (acov=="weighted") {
      my.x <- mapply("*", my.x, n, SIMPLIFY = FALSE)
      ## Weighted means = Cummulative sum of r*n/(sum of n)
      cov.mean <- Reduce("+", my.x)/pattern.n(x, n)
    }
    
    my.cov <- lapply(my.cov,
                     function(z) {
                       na.index <- is.na(z)
                       z[na.index] <- cov.mean[na.index]
                       z})
  }
   
  if (asyCovOld) {
    acovR <- asyCovOld(x=my.cov, n=n, cor.analysis=cor.analysis, acov=acov, ...)
  } else {
    acovR <- asyCov(x=my.cov, n=n, cor.analysis=cor.analysis, acov=acov, ...)
  }

  ## x is a list of "data", "n", and moderators, and append
  ## Append the moderators x[-c(1,2)] into data
  if (all(c(c("data", "n") %in% names(x), length(names(x))>2, append.vars))) {
    data <- suppressWarnings(data.frame(my.df, acovR, x[-c(1,2)], check.names=FALSE))        
  } else {
    data <- suppressWarnings(data.frame(my.df, acovR, check.names=FALSE))
  }
  
  ## Use unique row names if the row names are duplicated.
  if (row.names.unique) rownames(data) <- make.names(names(x), unique=TRUE)    

  #### Additional means
  if (!missing(Means)) {

    ## Some basic checks
    if (nrow(Means) != length(n)) {
      stop("Number of rows of 'Means' and length of 'n' are different.\n")
    }
    if (ncol(Means) != length(obslabels)) {
      stop("Number of columns of 'Means' and covariance matrices are different.\n")
    }
    if (!identical(colnames(my.cov[[1]]), colnames(Means))) {
      stop("The variable names are not in the same order in 'x' and 'Means'.
The results are likely incorrect unless this is what you want.\n")
    }

    ## Sampling covariance matrices of the means: covariance matrices/n
    acov_mean <- mapply(function(x, y) {x/y}, my.cov, n, SIMPLIFY=FALSE)
    acov_mean <- t(sapply(acov_mean, function(x) vech(x)))

    ## Variable names of p (sampling covariance matrix of the means
    pCovNames <- matrix(paste("C(",
                              outer(obslabels, obslabels, paste, sep = " "), 
                              ")", sep=""),
                        nrow=length(obslabels), ncol=length(obslabels))
    pCovNames <- vech(pCovNames)
    colnames(acov_mean) <- pCovNames
    
    if (row.names.unique) {
      rownames(acov_mean) <- make.names(rownames(acov_mean), unique=TRUE)    
    }
    
    out <- list(data=cbind(data, Means, acov_mean), n=n, obslabels=obslabels,
                ylabels=colnames(my.df), vlabels=colnames(acovR),
                vMlabels=pCovNames, VyMlabels=NULL)
  } else {
    ## Without the means
    out <- list(data=data, n=n, obslabels=obslabels, ylabels=colnames(my.df),
                vlabels=colnames(acovR), VyMlabels=NULL)
  }

  ## obslabels: labels of the means (Means) and dimnames in the
  ## correlation/covariance matrices (Cov)
  ## ylabels: vech(Cov) generated from list2matrix()
  ## vlabels: Acov of vech(Cov) generated from asyCov()
  ## Mlabels: labels of the means (not included as they are identical to obslabels)
  ## vMlabels: Acov of Means
  ## VyMlabels: Acov of Cov and Means (not included as Cov and Means are independent)

  out
}
