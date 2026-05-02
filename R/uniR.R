#' First Stage analysis of the univariate R (uniR) approach
#' 
#' It conducts the first stage analysis of the uniR analysis by pooling
#' elements of the correlation coefficients individually.
#' 
#' This function implements the univariate r approach proposed by Viswesvaran
#' and Ones (1995) to conduct meta-analytic structural equation modeling
#' (MASEM). It uses Schmidt and Hunter's approach to combine correlation
#' coefficients. It is included in this package for research interests. The
#' two-stage structural equation modeling (TSSEM) approach is preferred (e.g.,
#' Cheung, 2015; Cheung & Chan, 2005).
#' 
#' @param Cor A list of correlation matrices
#' @param n A vector of sample sizes
#' @param \dots Further arguments which are currently ignored
#' @return An object of class \code{uniR1} of the original data, the sample
#' sizes, the harmonic mean of sample sizes, the average correlation matrix,
#' the standard errors of the correlation matrix, and the standard deviations
#' (heterogeneity) of the correlation matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{uniR2}}, \code{\link[metaSEM]{Becker09}}
#' @references Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural
#' equation modeling approach}. Chichester, West Sussex: John Wiley & Sons,
#' Inc.
#' 
#' Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
#' modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10},
#' 40-64.
#' 
#' Schmidt, F. L., & Hunter, J. E. (2015). \emph{Methods of meta-analysis:
#' Correcting error and bias in research findings (3rd ed.)}. Thousand Oaks,
#' CA: Sage.
#' 
#' Viswesvaran, C., & Ones, D. S. (1995). Theory testing: Combining
#' psychometric meta-analysis and structural equations modeling.
#' \emph{Personnel Psychology}, \bold{48}, 865-885.
#' @keywords uniR
uniR1 <- function(Cor, n, ...) {
  ## Average correlation: Schmidt and Hunter (2015, Eq. 3.1)
  ## muliplied cor matrices with n
  r.sum <- mapply("*", Cor, n, SIMPLIFY = FALSE)
  ## replace NA with 0
  r.sum <- lapply(r.sum, function(x) {x[is.na(x)] <- 0; x})
  ## cummulative sum of r*n
  r.sum <- Reduce("+", r.sum)
  ## cumulative n
  n.sum <- pattern.n(Cor, n)
  ## average r weighted by n
  r.mean <- r.sum/n.sum

  ## Average squared error: Schmidt and Hunter (2015, p. 101, Eq. 3.2)
  r.diff2 <- lapply(Cor, function(x) (x-r.mean)^2)
  ## muliplied it by n
  r.diff2 <- mapply("*", r.diff2, n, SIMPLIFY = FALSE)
  ## replace NA with 0
  r.diff2 <- lapply(r.diff2, function(x) {x[is.na(x)] <- 0; x})
  ## cummulative sum of r*n
  r.diff2 <- Reduce("+", r.diff2)
  ## frequency-weighted average squared error
  r.S2 <- r.diff2/n.sum

  ## Average sampling error: Second improved approximation of Schmidt and Hunter (2015, p. 101, Eq. 3.7)
  ## No. of studies
  K <- pattern.na(Cor, show.na = FALSE)
  ## N.mean = sum of sample sizes / no. of studies
  N.mean <- n.sum/K
  r.SE2 <- (1-r.mean^2)^2 / (N.mean-1)

  ## Heterogeneity variance: Schmidt and Hunter (2015)
  r.SD2 <- r.S2 - r.SE2
  ## if negative, truncates to 0
  r.SD2[r.SD2<0] <- 0

  ## set the diagonals at NA
  diag(r.SE2) <- diag(r.SD2) <- NA

  ## cumulative n of the lower triangle
  n.harmonic <- n.sum[lower.tri(n.sum)]
  n.harmonic <- round( length(n.harmonic)/sum(1/n.harmonic), 0 )
  out <- list(data=Cor, n=n, r.mean=r.mean, r.SE=sqrt(r.SE2),
              r.SD=sqrt(r.SD2), n.harmonic=n.harmonic)
  class(out) <- "uniR1"
  out
}

print.uniR1 <- function(x, ...) {
  if (!is.element("uniR1", class(x)))
    stop("\"x\" must be an object of class \"uniR1\".")
  cat("\nTotal sample sizes: ", sum(x$n))
  cat("\nHarmonic mean of the sample sizes: ", x$n.harmonic, "\n")
  cat("\nAverage correlation matrix: ", "\n")
  print(x$r.mean)
  cat("\nSampling error (SE) of the average correlation matrix: ", "\n")
  print(x$r.SE)
  cat("\nPopulation heterogeneity (SD) of the average correlation matrix: ", "\n")
  print(x$r.SD)
}


## Assuming that there are dimnames in the pooled correlation matrix
## it may require some safety check with try() for simulation studies
uniR2mx <- function(x, RAM=NULL, Amatrix = NULL, Smatrix = NULL, Fmatrix = NULL,
                    model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    if (!is.element("uniR1", class(x)))
        stop("\"x\" must be an object of class \"uniR1\".")

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
            stop("\"Amatrix\" matrix is not specified.\n")
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

    ## No. of observed and latent variables
    p <- nrow(Smatrix@values) 
    
    r.mean <- x$r.mean

    ## If no variable name is given, use x1, x2, ...
    if (is.null(colnames(r.mean))) {
        var.labels <- paste0("x", seq_len(ncol(r.mean)))
        dimnames(r.mean) <- list(var.labels, var.labels)
    } else {
        var.labels <- colnames(r.mean)
    }

    ## Added labels for latent variables
    if (length(var.labels) < p) {
        var.labels <- c(var.labels, paste0("f", seq_len(p-length(var.labels))))
    }

    if (is.null(model.name)) model.name <- "UniR2"
    exp <- mxExpectationRAM(A="Amatrix", S="Smatrix", F="Fmatrix", dimnames=var.labels)

    mx.model <- mxModel(model.name,
                        mxData(observed=r.mean, type ="cov", numObs=x$n.harmonic),
                        Amatrix, Smatrix, Fmatrix, exp, mxFitFunctionML())

    ## Return mx model without running the analysis
    if (run==FALSE) return(mx.model)
    mxRun(mx.model, suppressWarnings=suppressWarnings, silent=silent, ...)
}

uniR2lavaan <- function(x, model, ...) {
  ## if (!requireNamespace("lavaan", quietly=TRUE))
  ##   stop("\"lavaan\" package is required for this function.")

  if (!is.element("uniR1", class(x)))
    stop("\"x\" must be an object of class \"uniR1\".")

  lavaan::sem(model=model, sample.nobs=x$n.harmonic, sample.cov=x$r.mean, fixed.x=FALSE,
              ...)
}
