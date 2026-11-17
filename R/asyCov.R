## Normalize intentional scalar reuse and align comparable named studies.
.prepare.asycov.input <- function(x, n) {
  if (!is.numeric(n) || is.complex(n) || !length(n) || anyNA(n) ||
      any(!is.finite(n)) || any(n <= 0))
    stop("'n' must contain positive finite sample sizes.", call.=FALSE)
  if (!is.list(x) || is.data.frame(x)) return(list(x=x, n=n))
  if (!length(x)) stop("'x' must contain at least one matrix.", call.=FALSE)
  if (length(x) != length(n)) {
    if (length(n) == 1L) n <- rep(n, length(x))
    else if (length(x) == 1L) x <- rep(x, length(n))
    else stop("Lengths of 'x' and 'n' must match, unless one has length one.", call.=FALSE)
  }
  ## Individual covariance matrices of different sizes remain supported.
  if (all(vapply(x, is.matrix, logical(1))) &&
      length(unique(vapply(x, nrow, integer(1)))) == 1L &&
      length(unique(vapply(x, ncol, integer(1)))) == 1L)
    x <- .align.study.matrices(x)
  list(x=x, n=n)
}

## Check only matrices whose entries are known. Missing entries are not zeros:
## a partial matrix may still have a positive-definite completion. Every fully
## specified 1- or 2-variable principal submatrix must itself be positive definite.
.check.asycov.matrix <- function(x, context="input", alternatives=FALSE,
                                  work.limit=100000L) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x))
    stop("The ", context, " matrix must be a square numeric matrix.", call.=FALSE)
  describe <- function(position) {
    rows <- rownames(x); cols <- colnames(x)
    row <- if (is.null(rows)) as.character(position[1L]) else rows[position[1L]]
    col <- if (is.null(cols)) as.character(position[2L]) else cols[position[2L]]
    paste0("[", row, ", ", col, "]")
  }
  nonfinite <- which(!is.na(x) & !is.finite(x), arr.ind=TRUE)
  if (nrow(nonfinite))
    stop("The ", context, " matrix contains non-finite observed entries, including ",
         describe(nonfinite[1L, ]), ".", call.=FALSE)
  asymmetric <- which(is.na(x) != is.na(t(x)), arr.ind=TRUE)
  if (nrow(asymmetric))
    stop("The ", context, " matrix has an asymmetric missing-value pattern at ",
         describe(asymmetric[1L, ]), "; provide both entries of each observed pair.", call.=FALSE)
  keep <- !is.na(diag(x))
  x <- x[keep, keep, drop=FALSE]
  if (!nrow(x)) return(invisible(NULL))
  other <- if (context == "weighted average") "acov='individual' or acov='unweighted'"
           else if (context == "unweighted average") "acov='individual' or acov='weighted'"
           else "acov='weighted' or acov='unweighted'"
  message <- paste0("The ", context, " matrix is not positive definite. Check the study matrices",
      if (alternatives) paste0("; consider ", other, " if appropriate. Other acov methods have different assumptions and may also fail") else " or provide a suitable working matrix", ".")
  check <- function(m) {
    if (any(!is.finite(m)) || !isSymmetric(unname(m), tol=1e-8) ||
        inherits(tryCatch(chol(m), error=function(e) e), "error"))
      stop(message, call.=FALSE)
  }
  for (i in seq_len(nrow(x))) check(x[i, i, drop=FALSE])
  if (!anyNA(x)) {
    check(x)
  } else if (nrow(x) > 1L) {
    for (i in seq_len(nrow(x)-1L)) for (j in (i+1L):nrow(x)) {
      pair <- x[c(i,j), c(i,j), drop=FALSE]
      if (!anyNA(pair)) check(pair)
    }
    ## Maximal fully observed principal submatrices are maximal cliques of
    ## the graph of known covariances. Bron-Kerbosch with pivoting checks each
    ## maximal clique; every other fully observed principal submatrix is
    ## contained in one of them. Searching all cliques can be exponential.
    known <- !is.na(x) & !is.na(t(x))
    diag(known) <- FALSE
    work <- 0L
    incomplete <- FALSE
    visit <- function(selected, candidates, excluded) {
      if (incomplete) return(invisible(NULL))
      work <<- work + 1L
      if (work > work.limit) {
        incomplete <<- TRUE
        return(invisible(NULL))
      }
      if (!length(candidates) && !length(excluded)) {
        if (length(selected) > 2L) check(x[selected, selected, drop=FALSE])
        return(invisible(NULL))
      }
      possible <- c(candidates, excluded)
      if (!length(possible)) return(invisible(NULL))
      degree <- vapply(possible, function(u) sum(known[u, candidates]), integer(1))
      pivot <- possible[which.max(degree)]
      branch <- candidates[!known[pivot, candidates]]
      for (v in branch) {
        neighbors <- which(known[v, ])
        visit(c(selected, v), intersect(candidates, neighbors),
              intersect(excluded, neighbors))
        if (incomplete) break
        candidates <- candidates[candidates != v]
        excluded <- c(excluded, v)
      }
      invisible(NULL)
    }
    visit(integer(), seq_len(nrow(x)), integer())
    if (incomplete)
      warning("Positive-definiteness validation of the ", context,
              " matrix was incomplete because the search limit was reached; ",
              "some fully observed submatrices were not checked.", call.=FALSE)
  }
  invisible(NULL)
}

## The Hessian approach needs all covariances among retained variables.
.check.asycovold.missing <- function(x) {
  present <- !is.na(diag(x))
  retained <- x[present, present, drop=FALSE]
  if (!nrow(retained))
    stop("asyCovOld() requires at least one observed variable.", call.=FALSE)
  if (anyNA(retained))
    stop("asyCovOld() requires all covariances among retained variables; use asyCov() for partial matrices or supply an explicit working matrix.", call.=FALSE)
  invisible(NULL)
}

#### Generate asymptotic covariance matrix of correlation/covariance matrix by *column major*
#' @rdname asyCov
asyCovOld <- function(x, n, cor.analysis=TRUE, dropNA=FALSE, as.matrix=TRUE,
                      acov=c("individual", "unweighted", "weighted"),
                      suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    .check.cov.names(x)
    input <- .prepare.asycov.input(x, n)
    x <- input$x
    n <- input$n

    makeSymmetric <- function(x) {x[lower.tri(x)] = t(x)[lower.tri(t(x))]; x}
    
    if (is.list(x)) {

        ## whether to use average correlation matrix
        acov <- match.arg(acov, c("individual", "unweighted", "weighted"))
        if (acov == "individual") {
            for (i in seq_along(x)) {
                .check.asycovold.missing(x[[i]])
                .check.asycov.matrix(x[[i]], paste0("study ", i), alternatives=TRUE)
            }
        }
        if (acov != "individual") {
            ## Replace NA with 0 before calculations
            my.x <- lapply(x, function(x) {x[is.na(x)] <- 0; x} )

            if (acov=="unweighted") {
                ## Unweighted means = sum of correlations/no. of studies
                my.x <- Reduce("+", my.x)/pattern.na(x, show.na = FALSE)
            } else {
                my.x <- mapply("*", my.x, n, SIMPLIFY = FALSE)
                ## Weighted means = Cumulative sum of r*n/sum of n
                my.x <- Reduce("+", my.x)/pattern.n(x, n)
            }

            .check.asycovold.missing(my.x)
            .check.asycov.matrix(my.x, paste(acov, "average"), alternatives=TRUE)
            ## cov2cor() replaces even missing diagonals by 1. Preserve
            ## entirely unreported variables so the legacy path can drop them.
            if (cor.analysis) {
                missing.diag <- is.na(diag(my.x))
                my.x <- suppressWarnings(cov2cor(my.x))
                diag(my.x)[missing.diag] <- NA_real_
            }

            ## Make sure that it is symmetric
            ## Fix the error message from OpenMx::verifyCovarianceMatrix
            my.x <- makeSymmetric(my.x)
            
            ## Repeat it to k studies
            x <- replicate(length(x), my.x, simplify = FALSE)    
        }
        ## whether to use average correlation matrix

        
        # Check if it returns a matrix or a list
        if (as.matrix) {
          # No. of variables
          if (length(unique(vapply(x, nrow, integer(1)))) != 1L)
            stop("The dimensions of matrices should be the same in order to stack them together!")
          
          if (is.null(dimnames(x[[1]]))) {
            oldNames <- paste0("x", 1:dim(x[[1]])[[1]])
          } else {
            oldNames <- dimnames(x[[1]])[[1]]
          }
          psCovNames <- .genCorNames(oldNames, cor.analysis=cor.analysis)$vlabels

          ## Fixed a bug before v.0.7-0 that uses only the first n 
          ## out.list <- lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
          ##                    suppressWarnings = suppressWarnings, dropNA = FALSE, ...)
          out.list <- mapply(asyCovOld, x, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = FALSE, ..., SIMPLIFY=FALSE)          
          #output
          out <- t(sapply(out.list, function(x) {(vech(x))}))
          
          ## http://stackoverflow.com/questions/17772916/using-correlation-matrices-for-meta-analytic-sem
          # Fixed a bug when 2x2 matrices of correlation
          # The output are incorrectly arranged in a row rather than in a column.
          if (dim(out)[1]==1) out <- t(out)
          ## A list(NULL, psCovNames) is required
          dimnames(out) <- list(NULL, psCovNames)
          
          out
        } else {
          #output
          ## lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
          ##                    suppressWarnings = suppressWarnings, dropNA = dropNA, ...)
          mapply(asyCovOld, x, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = dropNA, ..., SIMPLIFY=FALSE)          
        }
        
    } else {
      
        # Assumption: check the diagonals for missing data only
        miss.index <- is.na(Diag(x))
        .check.asycovold.missing(x)
        x.new <- x[!miss.index, !miss.index, drop=FALSE]
        .check.asycov.matrix(x.new, "input")
        p <- nrow(x.new)
        if (is.null(dimnames(x))) {
            oldNames <- paste("x", 1:nrow(x), sep = "")
            cNames <- oldNames[!miss.index]
            dimnames(x.new) <- list(cNames, cNames)
        } else {
            oldNames <- dimnames(x)[[1]]
            cNames <- dimnames(x.new)[[1]]
        }
        # create matrix of labels for ps
        psOldNames <- .genCorNames(oldNames, cor.analysis=cor.analysis)$ylabels
        full.index <- matrix(NA_integer_, length(oldNames), length(oldNames))
        full.index[lower.tri(full.index, diag=!cor.analysis)] <- seq_along(psOldNames)
        full.index[upper.tri(full.index)] <- t(full.index)[upper.tri(full.index)]
        kept.index <- full.index[!miss.index, !miss.index, drop=FALSE]
        acovName <- psOldNames[if (cor.analysis) vechs(kept.index) else vech(kept.index)]
        
        if (cor.analysis) {
            S <- mxMatrix("Stand", nrow = p, ncol = p, free = TRUE, values = jitter(vechs(cov2cor(x.new))), 
                name = "S", labels = acovName)
            D <- mxMatrix("Diag", nrow = p, ncol = p, free = TRUE, values = sqrt(Diag(x.new)), 
                name = "D")
            modelName <- "Asymptotic covariance matrix of correlation matrix"
        } else {
            S <- mxMatrix("Symm", nrow = p, ncol = p, free = TRUE, values = jitter(vech(x.new)), 
                name = "S", labels = acovName)
            D <- mxMatrix("Iden", nrow = p, ncol = p, name = "D")
            modelName <- "Asymptotic covariance matrix of covariance matrix"
        }
        expCov <- mxAlgebra(D %&% S, name = "expCov", dimnames = list(cNames, cNames))

        cModel <- mxModel(model = modelName, mxData(x.new, "cov", numObs = n), S, 
                          D, expCov, mxFitFunctionML(),
                          mxExpectationNormal("expCov", means=NA, dimnames = cNames))

        ## Return mx model without running the analysis
        if (run==FALSE) return(cModel)
        
        # try to run it with error message as output
        mxFit <- tryCatch(mxRun(cModel, silent=silent, suppressWarnings=suppressWarnings, ...), 
                          error = function(e) e)
        if (inherits(mxFit, "error")) {
            stop(print(mxFit))
        }
        # Need to multiply 2 to the inverse of Hessian matrix
        # http://openmx.psyc.virginia.edu/thread/360
        # Fixed a bug that all elements have to be inverted before selecting some of them
        acovS <- tryCatch(2 * solve(mxFit@output$calculatedHessian)[acovName, acovName, drop=FALSE], 
                              error = function(e) e)
        if (inherits(acovS, "error")) {
            stop(print(acovS))
        }

        ## No need to do it as [, drop=FALSE] has been added
        ## # When the dimensions are 1x1, dimnames are removed. Added them explicitly
        ## dimnames(acovS) <- list(acovName, acovName)

        if (dropNA) {          
          out <- acovS
        } else {          
          # oldNames include data for NA
          p <- length(psOldNames)
          out <- matrix(NA, nrow=p, ncol=p, dimnames=list(psOldNames, psOldNames))
          out[acovName, acovName] <- acovS
        }
        return(out)
    }
    
}

## Faster version without dropNA, suppressWarnings, silent, and run arguments.


#' Compute Asymptotic Covariance Matrix of a Correlation/Covariance Matrix
#'
#' Variable names and data column names must not contain periods. Rename them consistently, for example \code{x.1} to \code{x_1}, and check that the new names are unique. Study, group, and cluster values may contain periods.
#' 
#' It computes the asymptotic sampling covariance matrix of a
#' correlation/covariance matrix under the assumption of multivariate
#' normality.
#' 
#' 
#' @aliases asyCov asyCovOld
#' @param x A correlation/covariance matrix or a list of correlation/covariance
#' matrices. Missing-value patterns must be symmetric, and observed entries
#' must be finite. In \code{asyCov()}, missing input elements can propagate to sampling
#' covariance entries. In covariance analysis, entries that can be determined from
#' the available values are retained. A known zero times an unknown covariance
#' is treated as zero. A missing statistic has a missing row and column in its
#' sampling covariance matrix. For lists, \code{acov="weighted"} and
#' \code{acov="unweighted"} use available-study averages as a common working
#' matrix before calculating sampling covariances. A pair missing from every
#' study remains unknown. The resulting missing sampling covariances may need
#' an explicit approximation before fitting a meta-analysis. \code{asyCovOld()}
#' identifies missing variables using the diagonal. A one-variable covariance
#' matrix is supported; a one-variable correlation matrix has no correlations
#' to analyze.
#' @param n Positive finite sample size or vector of sample sizes. For a list
#' of studies, lengths must match unless one input has length one; named
#' study matrices of equal size are aligned by variable names.
#' @param cor.analysis Logical. The output is either a correlation or
#' covariance matrix.
#' @param dropNA Logical. If it is \code{TRUE}, the resultant dimensions will
#' be reduced by dropping the missing variables. If it is \code{FALSE}, the
#' resultant dimensions are the same as the input by keeping the missing
#' variables.
#' @param as.matrix Logical. If it is \code{TRUE} and \code{x} is a list of
#' correlation/covariance matrices with the same dimensions, the asymptotic
#' covariance matrices will be column vectorized and stacked together. If it is
#' \code{FALSE}, the output will be a list of asymptotic covariance matrices.
#' Note that if it is \code{TRUE}, \code{dropNA} will be \code{FALSE}
#' automatically. This option is useful when passing the asymptotic covariance
#' matrices to \code{\link[metaSEM]{meta}}
#' @param acov If it is \code{individual}, the sampling
#' variance-covariance matrices are calculated based on the individual
#' correlation/covariance matrix. If it is either \code{unweighted} or
#' \code{weighted}, the average correlation/covariance matrix is calculated
#' based on the unweighted or weighted mean with the sample sizes. The average
#' correlation/covariance matrix is used to calculate the sampling
#' variance-covariance matrices. Fully specified non-positive-definite input
#' or pooled working matrices are rejected. With missing entries, fully
#' specified principal submatrices are checked using a bounded search. A
#' warning is issued if the search limit is reached. Even a complete search
#' does not guarantee that the partial matrix has a positive-definite completion.
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed. It
#' is passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots It is ignored in \code{asyCov}. The additional arguments will
#' be passed to \code{\link[OpenMx]{mxRun}} in \code{asyCovOld}.
#' @return An asymptotic covariance matrix of the vectorized
#' correlation/covariance matrix or a list of these matrices. If
#' \code{as.matrix}=\code{TRUE} and \code{x} is a list of matrices, the output
#' is a stacked matrix.
#' @note Before 1.2.6, \code{asyCov} used an SEM approach based on Cheung and
#' Chan (2004). After 1.2.6, \code{asyCov} was rewritten based on Olkin and
#' Siotani (1976) for correlation matrix and Yuan and Bentler (2007, p. 371)
#' for covariance matrix. Arguments such as \code{dropNA},
#' \code{suppressWarnings}, \code{silent}, and \code{run} were dropped. The
#' original version was renamed to \code{asyCovOld} for compatibility. It is
#' retained as a legacy Hessian-based method; \code{asyCov} is preferred.
#' \code{asyCovOld} requires all covariances among variables with observed
#' diagonals. Supply an explicit working matrix or use \code{asyCov} for
#' partially observed pairs.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @references Cheung, M. W.-L., & Chan, W. (2004). Testing dependent
#' correlation coefficients via structural equation modeling.
#' \emph{Organizational Research Methods}, \bold{7}, 206-223.
#' 
#' Olkin, I., & Siotani, M. (1976). Asymptotic distribution of functions of a
#' correlation matrix. In S. Ideka (Ed.), \emph{Essays in probability and
#' statistics} (pp. 235-251). Shinko Tsusho.
#' 
#' Yuan, K.-H., & Bentler, P. M. (2007). Robust procedures in structural
#' equation modeling. In S.-Y. Lee (Ed.), \emph{Handbook of Latent Variable and
#' Related Models} (pp. 367-397). Elsevier/North-Holland.
#' @keywords meta-analysis
#' @examples
#' 
#' \donttest{
#' C1 <- matrix(c(1,0.5,0.4,0.5,1,0.2,0.4,0.2,1), ncol=3)  
#' asyCov(C1, n=100)
#' 
#' ## Data with missing values
#' C2 <- matrix(c(1,0.4,NA,0.4,1,NA,NA,NA,NA), ncol=3)  
#' C3 <- matrix(c(1,0.2,0.2,1), ncol=2)
#' 
#' ## Output is a stacked matrix of asymptotic covariance matrices
#' asyCov(list(C1,C2), n=c(100,50), as.matrix=TRUE)
#' 
#' ## Output is a stacked matrix of asymptotic covariance matrices
#' asyCov(list(C3,C3), n=c(100,50), as.matrix=TRUE)
#' 
#' ## Output is a list of asymptotic covariance matrices using the old version
#' asyCovOld(list(C1,C2,C3), n=c(100,50,50), dropNA=TRUE, as.matrix=FALSE)
#' }
#' 
asyCov <- function(x, n, cor.analysis=TRUE, as.matrix=TRUE,
                   acov=c("weighted", "individual", "unweighted"), ...) {

    .check.cov.names(x)
    input <- .prepare.asycov.input(x, n)
    x <- input$x
    n <- input$n

    ## Copy lower tri to upper tri
    makeSymmetric <- function(x) {x[upper.tri(x)] = t(x)[upper.tri(t(x))]; x}
    
    ## Main function to calculate acov for R or S multiplied by N
    N_Acov <- function(x, cor.analysis, context="input", alternatives=FALSE) {

        .check.asycov.matrix(x, context, alternatives)
        ## No. of variables
        p <- ncol(x)
        if (cor.analysis && p < 2L)
            stop("At least two variables are required for correlation analysis.")
        ## No. of correlations
        ps <- p*(p-1)/2
        
        ## Formula in cor handles NA automatically
        if (cor.analysis) {

            ## Make sure that it is a correlation matrix
            x <- suppressWarnings(cov2cor(x))  ## Suppress warnings in case of NA
  
            ## n*acov
            n_acov <- function(R, s, t, u, v) {
                0.5*R[s,t]*R[u,v]*(R[s,u]^2 + R[s,v]^2 + R[t,u]^2 + R[t,v]^2) + 
                    R[s,u]*R[t,v] + R[s,v]*R[t,u] - 
                    (R[s,t]*R[s,u]*R[s,v] + R[t,s]*R[t,u]*R[t,v] + R[u,s]*R[u,t]*R[u,v] + 
                     R[v,s]*R[v,t]*R[v,u])
            }
    
            ## index of variables  [1:p, 1:p] for correlation 1:p_s
            index <- matrix(NA, ncol=2, nrow=ps)
            k <- 1
        
            for (j in 1:(p-1))
                for (i in (j+1):p) {
                    ## cat(i, j, "\n")
                    index[k, ] <- c(i, j)
                    k <- k+1    
                }
    
            out <- matrix(NA, ncol=ps, nrow=ps)
            
            for (j in 1:ps) 
                for (i in j:ps) {
                    out[i, j] <- n_acov(R=x, s=index[i,1], t=index[i,2], u=index[j,1], v=index[j,2])
                }

            ## Copy the lower tri to upper tri
            out <- makeSymmetric(out)

            ## Analysis of covariance matrix
        } else {

            ## Normal-theory covariance of sample covariance elements:
            ## n Cov(s_ij, s_kl) = x_ik*x_jl + x_il*x_jk.
            ## Work in vech order and retain unknown inputs instead of
            ## substituting zero. A known zero times an unknown finite
            ## covariance is nevertheless known to be zero.
            pairs <- which(lower.tri(x, diag=TRUE), arr.ind=TRUE)
            i <- pairs[, 1L]
            j <- pairs[, 2L]
            multiply <- function(a, b) {
                product <- a*b
                known.zero <- (!is.na(a) & a == 0 & is.na(b)) |
                              (!is.na(b) & b == 0 & is.na(a))
                product[known.zero] <- 0
                product
            }
            out <- multiply(x[i, i, drop=FALSE], x[j, j, drop=FALSE]) +
                   multiply(x[i, j, drop=FALSE], x[j, i, drop=FALSE])

            ## A missing sample covariance itself has no sampling variance
            ## or covariance in the returned matrix, even if the formula is
            ## otherwise determined by zeros.
            missing <- is.na(x[cbind(i, j)])
            out[missing, ] <- out[, missing] <- NA_real_
        }
            
        ## Check dimnames is defined. If not, use x1, x2, etc
        if (is.null(colnames(x))) {
            obsvars <- paste0("x", seq_len(p))
        } else {
            obsvars <- colnames(x)
        }
            
        psvars <- .genCorNames(obsvars, cor.analysis=cor.analysis)$ylabels

        ## out is n*acov, not acov
        dimnames(out) <- list(psvars, psvars)

        ## Avoid problems in is.pd() or OpenMx.
        ## makeSymmetric(out)
        out
    }
    
    ## Only 1 matrix input
    ## Note. n_acov is n*acov!
    ## The output uses n.
    if (!is.list(x)) {

        out <- N_Acov(x, cor.analysis)
        
        if (length(n)==1) {
            out <- out/n
        } else {            
            ## n is a vecotr, return a list
            out <- lapply(n, function(z) out/z)
        }
        
        ## from if (!is.list(x))    
    } else {
            
        ## A list of matrix input
        acov <- match.arg(acov, c("weighted", "individual", "unweighted"))

        ## Calculate acov using individual correlation matrices
        if (acov=="individual") {

            ## A list of n*acov
            out <- lapply(seq_along(x), function(i)
                N_Acov(x[[i]], cor.analysis, paste0("study ", i), alternatives=TRUE))
            ## Divide them with n
            out <- mapply("/", out, n, SIMPLIFY=FALSE)

            ## Calculate acov using an average correlation matrix
        } else {
            
            ## Replace NA with 0 before calculations
            my.x <- lapply(x, function(z) {z[is.na(z)] <- 0; z} )
            
            if (acov=="unweighted") {
                ## Unweighted means = sum of correlations/no. of studies w/o NA
                my.x <- Reduce("+", my.x)/pattern.na(x, show.na=FALSE)
            } else {
                ## Weighted means = cummulative sum of r*n/sum of n w/o NA
                my.x <- mapply("*", my.x, n, SIMPLIFY=FALSE)                
                my.x <- Reduce("+", my.x)/pattern.n(x, n)
            }

            ## my.x is one correlation matrix
            out <- N_Acov(my.x, cor.analysis, paste(acov, "average"), alternatives=TRUE)
            out <- lapply(n, function(z) out/z)        
            
        } # if (acov=="individual")
    } ## from if (!is.list(x))


    ## Check if it is a list of matrices
    ## If 1 matrix, return it.
    ## If a list of matrices and as.matrix, stack them as a matrix
    if (is.list(out) & as.matrix) { 
 
        ## Assume the first one is with variable names, p*
        ## FIXME: psvars is a vector in ess but a column vector in RStudio. Different environments?
        psvars <- c(colnames(out[[1]]))
  
        ## Variable name p**, e.g., "C(x1_x2 x2_x3)"
        ## Name crash with matrixcalc::vech
        pssvars <- OpenMx::vech(outer(psvars, psvars,
                                      function(y, z) paste0("C(", y, " ", z, ")")))
        pssvars <- make.unique(as.character(pssvars), sep="_dup")
        
        ## Return a packed matrix selecting the lower triangles including the diagonals
        out <- t(sapply(out, vech))

        ## Ad-hoc to handle when there is only a column
        if (length(pssvars)==1) out <- matrix(out, ncol=1)

        colnames(out) <- pssvars

        ## Return a list of matrices
    } else {
        names(out) <- names(x)
    }

    out
}
