#### Generate asymptotic covariance matrix of correlation/covariance matrix by *column major*
asyCovOld <- function(x, n, cor.analysis=TRUE, dropNA=FALSE, as.matrix=TRUE,
                      acov=c("individual", "unweighted", "weighted"),
                      suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    makeSymmetric <- function(x) {x[lower.tri(x)] = t(x)[lower.tri(t(x))]; x}
    
    if (is.list(x)) {

        ## whether to use average correlation matrix
        acov <- match.arg(acov, c("individual", "unweighted", "weighted"))
        if (acov != "individual") {
            ## Replace NA with 0 before calculations
            my.x <- lapply(x, function(x) {x[is.na(x)] <- 0; x} )

            if (acov=="unweighted") {
                ## Unweighted means = sum of correlations/no. of studies
                my.x <- Reduce("+", my.x)/pattern.na(x, show.na = FALSE)
            } else {
                my.x <- mapply("*", my.x, n, SIMPLIFY = FALSE)
                ## Weighted means = Cummulative sum of r*n/sum of n
                my.x <- Reduce("+", my.x)/pattern.n(x, n)
            }

            ## Make sure that the diagonals are 1 for correlation analysis
            if (cor.analysis) my.x <- cov2cor(my.x)

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
          if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
            stop("The dimensions of matrices should be the same in order to stack them together!")
          
          if (is.null(dimnames(x[[1]]))) {
            oldNames <- paste0("x", 1:dim(x[[1]])[[1]])
          } else {
            oldNames <- dimnames(x[[1]])[[1]]
          }
          psOldNames <- outer(oldNames, oldNames, paste, sep = "_")          
          if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
          } else {
            psOldNames <- vech(psOldNames)
          }
          # cov/var of psOldNames
          psCovNames <- paste("C(", outer(psOldNames, psOldNames, paste, sep = " "), ")", sep="")
          psCovNames <- vech(matrix(psCovNames, nrow=length(psOldNames), ncol=length(psOldNames)))

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
        x.new <- x[!miss.index, !miss.index]
        if (!is.pd(x.new)) 
            stop("x is not positive definite!\n")
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
        psOldNames <- outer(oldNames, oldNames, paste, sep = "_")
        psMatnames <- outer(cNames, cNames, paste, sep = "_")
        
        if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
            acovName <- vechs(psMatnames)
            S <- mxMatrix("Stand", nrow = p, ncol = p, free = TRUE, values = jitter(vechs(cov2cor(x.new))), 
                name = "S", labels = acovName)
            D <- mxMatrix("Diag", nrow = p, ncol = p, free = TRUE, values = sqrt(Diag(x.new)), 
                name = "D")
            modelName <- "Asymptotic covariance matrix of correlation matrix"
        } else {
            psOldNames <- vech(psOldNames)
            acovName <- vech(psMatnames)
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
asyCov <- function(x, n, cor.analysis=TRUE, as.matrix=TRUE,
                   acov=c("weighted", "individual", "unweighted"), ...) {

    ## Copy lower tri to upper tri
    makeSymmetric <- function(x) {x[upper.tri(x)] = t(x)[upper.tri(t(x))]; x}
    
    ## Main function to calculate acov for R or S multiplied by N
    N_Acov <- function(x, cor.analysis) {

        ## No. of variables
        p <- ncol(x)
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

            if (!requireNamespace("matrixcalc", quietly=TRUE))    
                stop("\"matrixcalc\" package is required for this function.")
           
            ## index of NA in ps x px acov
            index_NA <- vech(is.na(x))

            ## Replace NA with 0 before calculations
            x[is.na(x)] <- 0
            
            ## Generalized inverse of the D matrix
            Dm <- matrixcalc::D.matrix(p)
            Dg <- MASS::ginv(Dm)
            ## Yuan & Bentler (Handbook of Latent Variable and Related Models, p. 371)
            out <- 2*Dg %*% (x %x% x) %*% t(Dg)

            ## Restore NA
            out[index_NA, ] <- out[, index_NA] <- NA
        }
            
        ## Check dimnames is defined. If not, use x1, x2, etc
        if (is.null(colnames(x))) {
            obsvars <- paste0("x", seq_len(p))
        } else {
            obsvars <- colnames(x)
        }
            
        if (cor.analysis) {
            psvars <- vechs(outer(obsvars, obsvars,
                                  function(x, y) paste(x, y, sep="_")))
        } else {
            psvars <- vech(outer(obsvars, obsvars,
                                 function(x, y) paste(x, y, sep="_")))
        }

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
            out <- lapply(x, N_Acov, cor.analysis=cor.analysis)
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
            out <- N_Acov(my.x, cor.analysis)
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
