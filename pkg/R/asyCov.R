#### Generate asymptotic covariance matrix of correlation/covariance matrix by *column major*
asyCov <- function(x, n, cor.analysis = TRUE, silent = TRUE, suppressWarnings = TRUE, ...) {
    if (is.list(x)) {
        lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent, suppressWarnings = suppressWarnings, ...)
    } else {
        # Assumption: check the diagonals for missing data only
        miss.index <- is.na(diag(x))
        x.new <- x[!miss.index, !miss.index]
        if (!is.pd(x.new)) 
            stop("x is not positive definite!\n")
        p <- nrow(x.new)
        if (is.null(dimnames(x))) {
            cNames <- paste("x", 1:nrow(x), sep = "")
            cNames <- cNames[!miss.index]
            dimnames(x.new) <- list(cNames, cNames)
        } else {
            cNames <- dimnames(x.new)[[1]]
        }
        # create matrix of labels for ps
        psMatnames <- outer(cNames, cNames, paste, sep = "")
        
        if (cor.analysis) {
            acovName <- vechs(psMatnames)
            S <- mxMatrix("Stand", nrow = p, ncol = p, free = TRUE, values = jitter(vechs(cov2cor(x.new))), 
                name = "S", labels = acovName)
            D <- mxMatrix("Diag", nrow = p, ncol = p, free = TRUE, values = sqrt(diag(x.new)), 
                name = "D")
            modelName <- "Asymptotic covariance matrix of correlation matrix"
        } else {
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

        return(acovS)
    }
}

