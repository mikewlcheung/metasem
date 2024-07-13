Cor2DataFrame <- function(x, n, v.na.replace=TRUE, cor.analysis=TRUE,
                          acov=c("weighted", "individual", "unweighted"),
                          Means, row.names.unique=FALSE, append.vars=TRUE,
                          asyCovOld=FALSE, ...) {

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

  if (asyCovOld) {
    acovR <- asyCovOld(x=my.cov, n=n, cor.analysis=cor.analysis, acov=acov, ...)
  } else {
    acovR <- asyCov(x=my.cov, n=n, cor.analysis=cor.analysis, acov=acov, ...)
  }

  ## NA is not allowed in definition variables
  ## They are replaced by 1e10
  if (v.na.replace) acovR[is.na(acovR)] <- 1e10

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
      stop("The variable names are not in the same order in 'x' and 'Means'. The results are likely incorrect unless this is what you want.\n")
    }
    
    ## Sampling covariance matrices of the means: covariance matrices/n
    ## NA are replaced with 10^5
    acov_mean <- mapply(function(x, y) {
      out <- x/y
      out[is.na(out)] <- 1e10
      out},
      my.cov, n, SIMPLIFY=FALSE)
    acov_mean <- t(sapply(acov_mean, function(x) vech(x)))

    ## Variable names of p (sampling covariance matrix of the means
    pCovNames <- matrix(paste("M(",
                              outer(obslabels, obslabels, paste, sep = " "), 
                              ")", sep=""),
                        nrow=length(obslabels), ncol=length(obslabels))
    pCovNames <- vech(pCovNames)
    colnames(acov_mean) <- pCovNames
    
    if (row.names.unique) {
      rownames(acov_mean) <- make.names(rownames(acov_mean), unique=TRUE)    
    }
    
    list(data=cbind(data, Means, acov_mean), n=n, obslabels=obslabels,
         ylabels=colnames(my.df), vlabels=colnames(acovR),
         Meanvlabels=pCovNames)
  } else {
    ## Without the means
    list(data=data, n=n, obslabels=obslabels, ylabels=colnames(my.df),
         vlabels=colnames(acovR))
  }
}
