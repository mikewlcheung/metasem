#### Generate asymptotic covariance matrix of correlation/covariance matrix by *column major*
asyCov <- function(x, n, cor.analysis = TRUE, dropNA = FALSE, as.matrix = TRUE,
                   silent = TRUE, suppressWarnings = TRUE, ...) {
    if (is.list(x)) {
      
        # Check if it returns a matrix or a list
        if (as.matrix) {
          # No. of variables
          if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
            stop("The dimensions of matrices should be the same in order to stack them together!")
          
          if (is.null(dimnames(x[[1]]))) {
            oldNames <- paste("x", 1:dim(x[[1]])[[1]], sep = "")
          } else {
            oldNames <- dimnames(x[[1]])[[1]]
          }
          psOldNames <- outer(oldNames, oldNames, paste, sep = "")          
          if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
          } else {
            psOldNames <- vech(psOldNames)
          }
          # cov/var of psOldNames
          psCovNames <- paste("cov(", outer(psOldNames, psOldNames, paste, sep = "_"), ")", sep="")
          psCovNames <- vech(matrix(psCovNames, nrow=length(psOldNames), ncol=length(psOldNames)))
          out.list <- lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = FALSE, ...)
          #output
          out <- t(sapply(out.list, function(x) {(vech(x))}))
          dimnames(out)[[2]] <- psCovNames
          out
        } else {
          #output
          lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = dropNA, ...)
        }
        
    } else {
      
        # Assumption: check the diagonals for missing data only
        miss.index <- is.na(diag(x))
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
        psOldNames <- outer(oldNames, oldNames, paste, sep = "")
        psMatnames <- outer(cNames, cNames, paste, sep = "")
        
        if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
            acovName <- vechs(psMatnames)
            S <- mxMatrix("Stand", nrow = p, ncol = p, free = TRUE, values = jitter(vechs(cov2cor(x.new))), 
                name = "S", labels = acovName)
            D <- mxMatrix("Diag", nrow = p, ncol = p, free = TRUE, values = sqrt(diag(x.new)), 
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
        objective <- mxMLObjective("expCov", dimnames = cNames)
        cModel <- mxModel(model = modelName, mxData(x.new, "cov", numObs = n), S, 
            D, expCov, objective)
        
        # try to run it with error message as output
        mxFit <- tryCatch(mxRun(cModel, silent = silent, suppressWarnings = suppressWarnings, ...), 
                          error = function(e) e)
        if (inherits(mxFit, "error")) {
            stop(print(mxFit))
        }
        # Need to multiply 2 to the inverse of Hessian matrix
        # http://openmx.psyc.virginia.edu/thread/360
        acovS <- tryCatch(2 * solve(mxFit@output$calculatedHessian[acovName, acovName] ), 
                              error = function(e) e)
        if (inherits(acovS, "error")) {
            stop(print(acovS))
        }

        # When the dimensions are 1x1, dimnames are removed. Added them explicitly
        dimnames(acovS) <- list(acovName, acovName)

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

