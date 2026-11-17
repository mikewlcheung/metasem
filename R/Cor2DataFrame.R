#' Convert correlation or covariance matrices into a dataframe of correlations
#' or covariances with their sampling covariance matrices
#' 
#' It converts the correlation or covariance matrices into a dataframe of
#' correlations or covariances with their asymptotic sampling covariance
#' matrices. It uses the \code{asyCov} at the backend.
#' Variable and appended moderator column names must not contain periods;
#' rename them consistently, e.g., \code{x.1} to \code{x_1}, and check for
#' duplicates. Multi-column moderator names generated with periods are
#' converted to underscores; collisions are rejected. Study names may contain periods. Generated pair labels that would collide are made unique with a \code{_dup} suffix; existing objects with ambiguous duplicate labels must be regenerated.
#' 
#' 
#' @param x A list of data with correlation/covariance matrix in \code{x$data}
#' and sample sizes \code{x$n}. Additional variables in \code{x}, in any
#' order relative to \code{data} and \code{n}, can be attached. Named matrices are aligned to the first matrix by variable name.
#' @param n If \code{x} is a list of correlation matrices without \code{x$data}
#' and \code{x$n}, a vector of sample sizes \code{n} must be provided.
#' @param v.na.replace Logical. For weighted or unweighted average sampling
#' covariance matrices, replace missing covariance entries by their available
#' weighted or unweighted mean before calculating sampling covariances. An
#' entry missing from every study cannot be replaced and raises an error.
#' @param cor.analysis Logical. The output is either a correlation or
#' covariance matrix.
#' @param acov If it is \code{weighted}, the average correlation/covariance
#' matrix is calculated based on the weighted mean with the sample sizes. The
#' average correlation/covariance matrix is used to calculate the sampling
#' variance-covariance matrices.
#' @param Means An optional matrix of means. The number of rows must be the
#' same as the length of \code{n}. The sampling covariance matrices of the
#' means are calculated by the covariance matrices divided by the sample sizes.
#' Therefore, it is important to make sure that covariance matrices (not
#' correlation matrices) are used in \code{x} when \code{Means} are included;
#' otherwise, the calculated sampling covariance matrices of the means are
#' incorrect.
#' @param row.names.unique Logical. If it is \code{FALSE} (the default), unique
#' row names are not created. If \code{TRUE}, the rows are named by the names
#' of the studies, made unique.
#' @param append.vars Whether to append the additional variables to the output
#' dataframe.
#' @param asyCovOld Whether to use the legacy Hessian-based method. It requires
#' complete covariances among retained variables in its working matrix.
#' \code{asyCov()} is preferred. See \code{\link[metaSEM]{asyCov}}.
#' @param \dots Further arguments to be passed to
#' \code{\link[metaSEM]{asyCov}}.
#' @return A list of components: (1) a data frame of correlations or
#' covariances with their sampling covariance matrices; (2) a vector of sample
#' sizes; (3) labels of the correlations; and (4) labels of their sampling
#' covariance matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{asyCov}}, \code{\link[metaSEM]{osmasem}},
#' \code{\link[metaSEM]{create.vechsR}}, \code{\link[metaSEM]{create.Tau2}},
#' \code{\link[metaSEM]{create.V}}
#' @keywords osmasem
#' @examples
#' 
#' \donttest{
#' ## Provide a list of correlation matrices and a vector of sample sizes as the inputs
#' my.df1 <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)
#' 
#' ## Add Lag time as a variable
#' my.df1$data <- data.frame(my.df1$data, Lag=Nohe15A1$Lag, check.names=FALSE)
#' 
#' ## Data
#' my.df1$data
#' 
#' ## Sample sizes
#' my.df1$n
#' 
#' ## ylabels
#' my.df1$ylabels
#' 
#' ## vlabels
#' my.df1$vlabels
#' 
#' #### Simplified version to do it
#' my.df2 <- Cor2DataFrame(Nohe15A1)
#' }
#' 
Cor2DataFrame <- function(x, n, v.na.replace=TRUE, cor.analysis=TRUE,
                          acov=c("weighted", "individual", "unweighted"),
                          Means, row.names.unique=FALSE, append.vars=TRUE,
                          asyCovOld=FALSE, ...) {
  
  acov <- match.arg(acov, c("weighted", "individual", "unweighted"))
  
  ## x is a list of "data", "n", ...
  if (all(c("data", "n") %in% names(x))) {
    my.cov <- x$data
    n <- x$n
  } else {
    ## x is just a list of correlation matrices. "n" is provided as an argument.
    my.cov <- x
  }

  .check.cov.names(my.cov)
  if (!missing(Means)) .check.no.period(colnames(Means), "Mean variable names")
  if (append.vars && all(c("data", "n") %in% names(x)))
    .check.no.period(setdiff(names(x), c("data", "n")), "Moderator column names")

  if (length(my.cov) != length(n)) stop("Lengths of 'x' and 'n' are different.\n")
  my.cov <- .align.study.matrices(my.cov)
  obslabels <- colnames(my.cov[[1]])
  
  if (cor.analysis) {
    ## Standardize and then vechs()
    my.df <- list2matrix(x=suppressWarnings(lapply(my.cov, cov2cor)), diag=FALSE)
  } else {
    ## vech()
    my.df <- list2matrix(x=my.cov, diag=TRUE)
  }

  ## NA is not allowed in definition variables
  ## Check if there are NAs in my.cov and replace them with the weighted or unweighted averages
  if (any(sapply(my.cov, is.na)) && v.na.replace && acov != "individual") {
    ## Replace NA with 0 before calculations
    my.x <- lapply(my.cov, function(z) {z[is.na(z)] <- 0; z} )
    if (acov=="unweighted") {
      ## x: original covariance matrices in the input
      ## Unweighted means = sum of r/(no. of studies)
      counts <- pattern.na(my.cov, show.na = FALSE)
      if (any(counts == 0 & Reduce("|", lapply(my.cov, is.na))))
        stop("Cannot replace an entry missing from every covariance matrix.")
      cov.mean <- Reduce("+", my.x)/counts
    } else if (acov=="weighted") {
      my.x <- mapply("*", my.x, n, SIMPLIFY = FALSE)
      ## Weighted means = Cummulative sum of r*n/(sum of n)
      counts <- pattern.n(my.cov, n)
      if (any(counts == 0 & Reduce("|", lapply(my.cov, is.na))))
        stop("Cannot replace an entry missing from every covariance matrix.")
      cov.mean <- Reduce("+", my.x)/counts
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

  ## x is a list of "data", "n", and moderators. Append the moderators,
  ## selected by name rather than by position, into data.
  moderators <- if (all(c("data", "n") %in% names(x)))
    x[setdiff(names(x), c("data", "n"))] else NULL
  if (length(moderators) && append.vars) {
    for (m in moderators) if (is.matrix(m) || is.data.frame(m))
      .check.no.period(colnames(m), "Moderator column names")
    data <- suppressWarnings(data.frame(my.df, acovR, moderators, check.names=FALSE))
    base.cols <- ncol(my.df) + ncol(acovR)
    if (ncol(data) > base.cols) {
      moderator.cols <- seq.int(base.cols + 1L, ncol(data))
      colnames(data)[moderator.cols] <- gsub(".", "_", colnames(data)[moderator.cols], fixed=TRUE)
    }
  } else {
    data <- suppressWarnings(data.frame(my.df, acovR, check.names=FALSE))
  }

  ## Use unique row names (the names of the studies) if they are duplicated.
  if (row.names.unique) {
    study.names <- names(my.cov)
    if (is.null(study.names)) study.names <- as.character(seq_along(my.cov))
    rownames(data) <- make.names(study.names, unique=TRUE)
  }

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
  if (anyDuplicated(colnames(out$data)))
    stop("Generated correlation, sampling covariance, mean, or moderator column names collide; rename the conflicting variables.", call.=FALSE)

  ## obslabels: labels of the means (Means) and dimnames in the
  ## correlation/covariance matrices (Cov)
  ## ylabels: vech(Cov) generated from list2matrix()
  ## vlabels: Acov of vech(Cov) generated from asyCov()
  ## Mlabels: labels of the means (not included as they are identical to obslabels)
  ## vMlabels: Acov of Means
  ## VyMlabels: Acov of Cov and Means (not included as Cov and Means are independent)

  out
}
